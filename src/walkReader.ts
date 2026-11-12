import { RemoteFile } from 'generic-filehandle2'

import { refNameLookup, resolveRefName } from './refNames.ts'
import { checkRange } from './tabixSource.ts'
import {
  blocksOver,
  canonicalLink,
  chunkRef,
  decodeChunk,
  parseBlocks,
  pieceEnds,
  pieceLinks,
  referenceNodes,
  unloadedRow,
} from './walkChunk.ts'
import {
  ATTACH_BP,
  departureOwnTiles,
  parseDepartureRow,
  squeezer,
} from './walkDeparture.ts'
import { WalkFile, WalkFormatError, once, raceAbort } from './walkFile.ts'
import {
  WalkGraph,
  byteBudgetError,
  edgeKey,
  isPlaceholder,
  joinPieces,
  keptNodesError,
  keptStepsError,
  megabytes,
  stepBudgetError,
  surelyKept,
  walkCut,
  walkHeader,
  walksKept,
} from './walkRows.ts'

import type {
  AbridgedPiece,
  ChunkRef,
  DecodedChunk,
  PieceBlocks,
  PlacedNode,
  RefPath,
} from './walkChunk.ts'
import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkFilehandle, WalkUnit } from './walkFile.ts'
import type { BudgetNames, WalkHeader, WalkRow } from './walkRows.ts'

// A window reads whole tiles and keeps the steps on it, so a narrow window
// reads several times what it keeps: 464 haplotypes across the 1.5 kb ABCA7
// repeat read 5.0 M steps and keep 1.16 M, and the cut takes 0.4 s from warm
// files. The step budget bounds what the cut keeps, and rows are read to this
// multiple of it.
const READ_STEPS_PER_KEPT = 8

// bp of reference a cut keeps past each end of the window
const WALK_CONTEXT = 1000

// detail records an expanding cut reads at once
const DETAIL_READS = 16

// the runs off R a reader keeps the placeholder ids of, and the most ids it
// gives before it starts again, which stay above walkChunk's SKIPPED
const RUN_IDS = 1 << 20
const MAX_RUN_COUNT = 2 ** 30

// The budgets a cut keeps to unless its caller passes others: what a browser
// window can hold. Infinity lifts one
export const WALK_BUDGETS = {
  stepBudget: 4_000_000,
  byteBudget: 8_000_000,
  departureByteBudget: 2_000_000,
  attachBp: ATTACH_BP,
}

/**
 * A window whose departure read would fetch more than its budget. Zooming in
 * does not lower the read, since a record is filed over its whole stretch,
 * so this is not a RegionTooLargeError and shows no zoom-in notice.
 */
export class DepartureReadError extends Error {
  override name = 'DepartureReadError'
  readonly bytes: number
  readonly limit: number

  constructor(bytes: number, limit: number) {
    super(
      `Too many departures from the reference over this window to read (${megabytes(bytes)} against a budget of ${megabytes(limit)}). A collapsed repeat usually causes this: haplotypes hop between its copies, and each hop is a departure`,
    )
    this.bytes = bytes
    this.limit = limit
  }
}

// Reads `items`, at most `limit` at once and none once `signal` aborts, and
// hands each to `use` in order
export async function readInOrder<T, R>(
  items: T[],
  limit: number,
  signal: AbortSignal,
  read: (item: T) => Promise<R>,
  use: (item: T, result: R) => void,
) {
  const reads: Promise<R>[] = []
  const next = () => {
    if (reads.length < items.length && !signal.aborted) {
      const started = read(items[reads.length]!)
      reads.push(started)
      started.then(next, () => {})
    }
  }
  for (let k = 0; k < limit; k++) {
    next()
  }
  for (let i = 0; i < items.length; i++) {
    signal.throwIfAborted()
    use(items[i]!, await reads[i]!)
  }
}

export interface WalkCutOptions {
  signal?: AbortSignal
  // which walks to decode, by `sample#haplotype#contig`; undefined is all
  keep?: (name: string) => boolean
  // the most steps a cut keeps, WALK_BUDGETS' by default; Infinity has no
  // limit. The rows of the window's tiles and expanded detail records may
  // hold READ_STEPS_PER_KEPT times as many; those of tiles fetched for own
  // steps are not counted
  stepBudget?: number
  // the most compressed bytes a cut fetches from the tiles under its window
  // and context and the tiles fetched for records' own steps, WALK_BUDGETS'
  // by default; Infinity has no limit. The detail records `expand` reads
  // count only against stepBudget
  byteBudget?: number
  // the most compressed bytes the departure read fetches, WALK_BUDGETS' by
  // default; Infinity has no limit
  departureByteBudget?: number
  // the most bp of own sequence a departure gets its steps attached for,
  // ATTACH_BP by default. The cut fetches every tile holding own steps of
  // such a record that its read lacks, which gives those records their steps
  // and nothing else
  attachBp?: number
  /**
   * @deprecated Ignored: every eligible record's tiles are fetched, their
   * bytes counted against byteBudget
   */
  attachFetches?: number
  context?: number
  // false lists no links, for a caller that draws no edges
  links?: boolean
  // keep the bases of the nodes, from a file written with --sequences
  sequences?: boolean
  // how a zoom-in notice names the budgets
  budgetNames?: BudgetNames
  // the most nodes a cut keeps; Infinity, the default, has no limit
  nodeLimit?: number
  // keep the reference over the window and its context alone, and split each
  // walk where it steps onto the reference off them
  contained?: boolean
  // the abridged walks to read from their detail records, by
  // `sample#haplotype#contig`; none by default. Of such a walk the cut reads
  // the records of the pieces whose steps it keeps or attaches
  expand?: (name: string) => boolean
}

export interface RegionCutOptions extends Omit<WalkCutOptions, 'keep'> {
  // the `sample#haplotype` walks to keep beside the reference; every walk when
  // undefined or empty
  haplotypes?: string[]
  // the sample whose sequence the region names, the file's first reference
  // by default
  referenceSample?: string
}

export type WalkReaderHeader = WalkHeader & { tile: number }

const overlaps = (unit: WalkUnit, span: { start: number; end: number }) =>
  unit.length > 0 && unit.start < span.end && unit.end > span.start

const bytesOf = (units: WalkUnit[]) =>
  units.reduce((sum, unit) => sum + unit.length, 0)

// The stretch a cut draws: the window and its context
const drawnSpan = (
  region: { start: number; end: number },
  context: number,
) => ({
  start: Math.max(0, region.start - context),
  end: region.end + context,
})

// The tiles an attach fetch reads for a record: those its `ownTiles` names
export function attachTiles(units: WalkUnit[], ownTiles: number[]) {
  return units.filter(
    u => u.kind === 'tile' && u.length > 0 && ownTiles.includes(Number(u.name)),
  )
}

const pieceKey = (
  name: string,
  fragStart: number | string,
  piece: number | string,
) => `${name}\t${fragStart}\t${piece}`

// A unit to read, and the base a window must reach to need it
interface Needed {
  unit: WalkUnit
  from: number
}

// An abridged row a tile's decode met, with its tile, the tile's R and rows,
// and the rows its piece joins
interface Met extends Needed {
  piece: AbridgedPiece
  path: RefPath
  tileRows: string[][]
  into: WalkRow[]
}

// An abridged row with its blocks row
interface Read extends Met {
  blocks: PieceBlocks
  on: ReturnType<typeof blocksOver>
}

// Where each piece of `rows` starts on its contig
const pieceStarts = (rows: WalkRow[]) => {
  const starts = new Map<string, number>()
  for (const row of rows) {
    starts.set(pieceKey(row.name, row.fragStart, row.piece), row.hapOffset)
  }
  return starts
}

const keptNodes = (kept: Set<number>) => {
  let count = 0
  for (const id of kept) {
    count += isPlaceholder(id) ? 0 : 1
  }
  return count
}

// The directory beside a data file's URL, before any query or fragment
const directoryOf = (url: string) => url.replace(/(?=$|[?#])/, '.idx')

// The records of departure lines over `span`, each line once, of the walks
// kept, with the tiles holding each one's own steps
function departureRecords(
  lines: string[],
  span: { start: number; end: number },
  keep: ((name: string) => boolean) | undefined,
) {
  const records: WalkDeparture[] = []
  const ownTiles = new Map<WalkDeparture, number[]>()
  for (const line of new Set(lines)) {
    const cols = line.split('\t')
    if (
      cols[3] === 'departure' &&
      +cols[1]! < span.end &&
      +cols[2]! > span.start
    ) {
      const record = parseDepartureRow(cols)
      if (keep === undefined || keep(record.name)) {
        records.push(record)
        ownTiles.set(record, departureOwnTiles(cols))
      }
    }
  }
  return { records, ownTiles }
}

/**
 * A `gfa-to-tabix --walks` file of walk format 8, `<prefix>.<sample>.walks.gz`
 * and its directory `<prefix>.<sample>.walks.gz.idx`, read a window at a
 * time. At open the directory, the header and the sequence's unit table; for
 * a window, the tiles and departure blocks overlapping the window and its
 * context, in one read wherever their bytes touch, each tile decoded once;
 * then every tile holding the own steps of records within `attachBp` the
 * read lacks, which give those steps alone, the `blocks` unit of
 * each tile that abridges a walk kept, and the detail records of the walks
 * `expand` names. The
 * table gives each read's bytes before it is made;
 * a window whose tiles are past the byte budget reads no row, and one whose
 * attach tiles take it past reads none of them. Only the walks `keep`
 * accepts are decoded. The step budget bounds the steps the cut keeps, and
 * rows of the window holding READ_STEPS_PER_KEPT times as many stop the
 * read.
 */
export class WalkReader {
  private readonly file: WalkFile

  constructor(file: WalkFile) {
    this.file = file
  }

  /**
   * The walk file at `url` and its directory at `url` plus `.idx`, before any
   * query string, each opened by `open`, a RemoteFile by default
   */
  static open(
    url: string,
    open: (url: string) => WalkFilehandle = location =>
      new RemoteFile(location),
  ) {
    return WalkReader.fromFilehandles({
      data: open(url),
      directory: open(directoryOf(url)),
    })
  }

  // The walk file `data` with its directory, the `.idx` file
  static fromFilehandles({
    data,
    directory,
  }: {
    data: WalkFilehandle
    directory: WalkFilehandle
  }) {
    return new WalkReader(new WalkFile(data, directory))
  }

  private readonly readHeader = once(async () => {
    const header = walkHeader(await this.file.header())
    const { tile } = header
    if (tile === undefined) {
      throw new WalkFormatError("The walk file's #walks line states no tile:i:")
    }
    return { ...header, tile }
  })

  // the id of each run of steps off R a cut left unloaded, by its walk and
  // its stretch of the contig: the same run has the same id in every cut
  private runIds = new Map<string, Map<string, number>>()
  private runsHeld = 0
  private runCount = 0

  private runId(name: string, start: number, end: number) {
    const ids = this.runIds.get(name) ?? new Map<string, number>()
    this.runIds.set(name, ids)
    const key = `${start}-${end}`
    const held = ids.get(key)
    const id = held ?? -++this.runCount
    ids.set(key, id)
    this.runsHeld += held === undefined ? 1 : 0
    return id
  }

  // Past RUN_IDS runs the ids are let go, between cuts: a run met again takes
  // a new id, and none is given twice
  private forgetRuns() {
    if (this.runsHeld > RUN_IDS || this.runCount > MAX_RUN_COUNT) {
      this.runIds = new Map()
      this.runsHeld = 0
      this.runCount = this.runCount > MAX_RUN_COUNT ? 0 : this.runCount
    }
  }

  private readonly readLookup = once(async () =>
    refNameLookup(await this.file.sequences()),
  )

  // The header, read once for every caller: one caller giving up stops only
  // its own wait, and a read that fails is tried again by the next caller
  header(opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.readHeader(), opts.signal)
  }

  // The reference sequences the walk file holds windows on
  referenceSequences(opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.file.sequences(), opts.signal)
  }

  /**
   * The walk file's sequence for `refName` on `sample`, or for `refName` as
   * the file spells it: `chr6` on `GRCh38` is `GRCh38#0#chr6`. Undefined when
   * the file holds neither.
   */
  async sequenceFor(
    sample: string,
    refName: string,
    opts: { signal?: AbortSignal } = {},
  ) {
    const lookup = await raceAbort(this.readLookup(), opts.signal)
    return resolveRefName(lookup, sample, refName)
  }

  // The units the file holds for a reference sequence, in file order
  units(refName: string, opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.file.units(refName), opts.signal)
  }

  // Where the reference's last node on `refName` ends, from the extents of
  // its tiles; 0 for a sequence the file lacks
  async sequenceLength(refName: string, opts: { signal?: AbortSignal } = {}) {
    const units = await this.units(refName, opts)
    return units.reduce(
      (end, u) => (u.kind === 'tile' ? Math.max(end, u.end) : end),
      0,
    )
  }

  // The compressed bytes a cut of `region` reads before any attach fetch or
  // expansion: its tiles, their blocks units and its departure blocks
  async bytesFor(
    refName: string,
    region: { start: number; end: number },
    opts: { signal?: AbortSignal; context?: number } = {},
  ) {
    checkRange(refName, region.start, region.end)
    const drawn = drawnSpan(region, opts.context ?? WALK_CONTEXT)
    const units = await this.units(refName, opts)
    const tiles = new Set(
      units
        .filter(u => u.kind === 'tile' && overlaps(u, drawn))
        .map(u => u.name),
    )
    return bytesOf(
      units.filter(u =>
        u.kind === 'blocks'
          ? tiles.has(u.name)
          : (u.kind === 'tile' || u.kind === 'departures') &&
            overlaps(u, drawn),
      ),
    )
  }

  // Every pass of the reference over `region`, by start, from the ref rows of
  // the tiles there: a node the reference passes twice there comes twice
  async referenceNodes(
    refName: string,
    region: { start: number; end: number },
    opts: { signal?: AbortSignal } = {},
  ) {
    checkRange(refName, region.start, region.end)
    const tiles = (await this.units(refName, opts)).filter(
      u => u.kind === 'tile' && overlaps(u, region),
    )
    const read = await this.file.read(tiles, opts.signal)
    return referenceNodes(read.flat().map(line => line.split('\t')))
      .filter(
        node =>
          node.refName === refName &&
          node.end > region.start &&
          node.start < region.end,
      )
      .sort((a, b) => a.start - b.start || a.end - b.end)
  }

  async cut(
    refName: string,
    region: { start: number; end: number },
    opts: WalkCutOptions = {},
  ) {
    const {
      signal,
      keep,
      context = WALK_CONTEXT,
      stepBudget = WALK_BUDGETS.stepBudget,
      byteBudget = WALK_BUDGETS.byteBudget,
      departureByteBudget = WALK_BUDGETS.departureByteBudget,
      attachBp = ATTACH_BP,
      links = true,
      contained = false,
      nodeLimit = Infinity,
      expand,
    } = opts
    checkRange(refName, region.start, region.end)
    this.forgetRuns()
    const [header, units] = await Promise.all([
      this.header({ signal }),
      this.units(refName, { signal }),
    ])
    const { tile } = header
    const drawn = drawnSpan(region, context)
    const tiles = units.filter(u => u.kind === 'tile' && overlaps(u, drawn))
    const departureBlocks = units.filter(
      u => u.kind === 'departures' && overlaps(u, drawn),
    )
    const tileBytesTo = (end: number) =>
      bytesOf(tiles.filter(u => u.start < end + context))
    const checkBytes = (
      bytesTo: (end: number) => number,
      attaching = false,
    ) => {
      if (Number.isFinite(byteBudget)) {
        const tooLarge = byteBudgetError(
          bytesTo,
          byteBudget,
          region,
          tile,
          opts.budgetNames,
          attaching,
        )
        if (tooLarge) {
          throw tooLarge
        }
      }
    }
    checkBytes(tileBytesTo)
    const departureBytes = bytesOf(departureBlocks)
    if (departureBytes > departureByteBudget) {
      throw new DepartureReadError(departureBytes, departureByteBudget)
    }
    const reads = new AbortController()
    const stop = () => {
      reads.abort(signal?.reason)
    }
    signal?.throwIfAborted()
    signal?.addEventListener('abort', stop, { once: true })
    const rows: WalkRow[] = []
    // the rows of tiles fetched to attach, which give records their own steps
    // and nothing else
    const attachRows: WalkRow[] = []
    const graph = new WalkGraph({
      refName,
      start: region.start,
      end: region.end,
      context,
      sequences: opts.sequences,
    })
    // the name filters each tile was decoded under, by its offset
    const decodedBy = new Map<number, ((name: string) => boolean)[]>()
    const decodedFor = (unit: WalkUnit, name: string) =>
      decodedBy.get(unit.offset)?.some(took => took(name)) === true
    const ends: DecodedChunk['ends'] = []
    const surely = new Set<number>()
    const bounds = {
      contained,
      start: region.start - context,
      end: region.end + context,
    }
    const filtered = keep !== undefined
    const unitsOf = (kind: string) =>
      new Map(units.filter(u => u.kind === kind).map(u => [u.name, u]))
    const detailUnits = unitsOf('detail')
    const blocksUnits = unitsOf('blocks')
    // the units read past the window's tiles, each with the base a window
    // must reach to need it
    const extra: Needed[] = []
    const bytesWith = (more: Needed[]) => (end: number) =>
      tileBytesTo(end) +
      bytesOf(extra.filter(x => x.from < end + context).map(x => x.unit)) +
      bytesOf(more.filter(x => x.from < end + context).map(x => x.unit))
    // the abridged rows decoded and not yet read from a blocks row or a
    // detail record
    let met: Met[] = []
    const detailsRead = new Set<string>()
    const refs = new Map<WalkUnit, ChunkRef>()
    // the drawn stretch, as far as the reference's passes over it reach
    const over = { ...drawn }
    // the ends of the records whose own steps the cut attaches
    const attachLeaves = new Set<string>()
    const attachRejoins = new Set<string>()
    const checkNodes = (
      chunk: { nodes: Map<number, unknown>; rows: WalkRow[] },
      readTo: number,
    ) => {
      if (Number.isFinite(nodeLimit)) {
        surelyKept(chunk, graph.nodes, bounds, surely)
        if (surely.size > nodeLimit) {
          surely.clear()
          surelyKept({ nodes: graph.nodes, rows }, graph.nodes, bounds, surely)
        }
        const tooMany = keptNodesError(
          surely.size,
          nodeLimit,
          region,
          filtered,
          readTo,
        )
        if (tooMany) {
          throw tooMany
        }
      }
    }
    // links as pairs of handles, each written the one way round
    const addLinks = (pairs: number[]) => {
      if (links) {
        for (let k = 0; k < pairs.length; k += 2) {
          const [a, b] = canonicalLink(pairs[k]!, pairs[k + 1]!)
          graph.addLinkEnds(a, b)
        }
      }
    }
    const placeNodes = (nodes: Map<number, PlacedNode>) => {
      for (const [id, node] of nodes) {
        graph.placeNode(id, node)
      }
    }
    const place = (tileCut: DecodedChunk, readTo: number) => {
      placeNodes(tileCut.nodes)
      for (const { id, ...node } of tileCut.revisits) {
        graph.placeNode(id, node)
      }
      graph.addSequences(tileCut.sequences)
      graph.addChunk(tileCut.chunk)
      for (const end of tileCut.ends) {
        ends.push(end)
      }
      checkNodes(tileCut, readTo)
    }
    // a tile is decoded into `into`, with its nodes and the links between
    // pieces, for the names `only` takes that no earlier decode of it took,
    // and a cut over its node limit stops there
    const decode = (
      unit: WalkUnit,
      lines: string[],
      into: WalkRow[],
      from: number,
      only?: (name: string) => boolean,
    ) => {
      const wanted = only ?? keep ?? (() => true)
      const earlier = decodedBy.get(unit.offset) ?? []
      decodedBy.set(unit.offset, [...earlier, wanted])
      const tileRows = lines.map(line => line.split('\t'))
      const tileCut = decodeChunk(
        tileRows,
        Number(unit.name),
        name => wanted(name) && !earlier.some(took => took(name)),
        opts.sequences,
      )
      const { path } = tileCut.chunk
      const hull = into === rows ? path.hull(drawn) : undefined
      over.start = Math.min(over.start, hull?.start ?? over.start)
      over.end = Math.max(over.end, hull?.end ?? over.end)
      for (const row of tileCut.rows) {
        into.push(row)
      }
      for (const piece of tileCut.abridged) {
        met.push({ unit, piece, path, tileRows, into, from })
      }
      place(tileCut, unit.end)
    }
    const readDetail = async ({ unit }: { unit: WalkUnit }) => {
      const [lines = []] = await this.file.read([unit], reads.signal)
      return lines.map(line => line.split('\t'))
    }
    const decodeDetail = (
      {
        ref,
        piece,
        into,
      }: { ref: ChunkRef; piece: AbridgedPiece; into: WalkRow[] },
      held: string[][],
    ) => {
      const tileCut = decodeChunk(
        held,
        piece.chunkStart,
        name => name === piece.name,
        opts.sequences,
        ref,
      )
      for (const decodedRow of tileCut.rows) {
        if (decodedRow.name === piece.name) {
          into.push(decodedRow)
        }
      }
      place(tileCut, piece.chunkStart + tile)
    }
    const onDrawn = (span: { start: number; end: number }) =>
      span.start < over.end && span.end > over.start
    // each window piece's blocks to build whole: those from the first block
    // over the drawn stretch to the last, of the pieces that join it
    const coresOf = (pieces: Read[], others: Met[]) => {
      const lines = new Map<
        string,
        { piece: number; on: boolean; read?: Read }[]
      >()
      for (const read of pieces) {
        const { name, fragStart, piece } = read.piece
        const key = `${name}\t${fragStart}`
        const line = lines.get(key) ?? []
        lines.set(key, line)
        line.push({ piece, on: read.on.first >= 0, read })
      }
      for (const { piece } of others) {
        lines.get(`${piece.name}\t${piece.fragStart}`)?.push({
          piece: piece.piece,
          on: onDrawn(piece.span),
        })
      }
      for (const row of rows) {
        lines.get(`${row.name}\t${row.fragStart}`)?.push({
          piece: row.piece,
          on:
            row.span === undefined ||
            onDrawn({ start: row.span[0], end: row.span[1] }),
        })
      }
      const cores = new Map<Read, [number, number]>()
      for (const line of lines.values()) {
        line.sort((a, b) => a.piece - b.piece)
        let from = 0
        for (let i = 1; i <= line.length; i++) {
          if (i === line.length || line[i]!.piece !== line[i - 1]!.piece + 1) {
            let first = -1
            let last = -1
            for (let j = from; j < i; j++) {
              if (line[j]!.on) {
                first = first < 0 ? j : first
                last = j
              }
            }
            for (let j = Math.max(first, 0); j <= last; j++) {
              const read = line[j]!.read
              if (read) {
                cores.set(read, [
                  j === first ? read.on.first : 0,
                  j === last ? read.on.last : read.blocks.begins.length - 1,
                ])
              }
            }
            from = i
          }
        }
      }
      return cores
    }
    const attachNeeds = ({ piece, blocks }: Read) => {
      const { leave, rejoin } = pieceEnds(piece, blocks)
      return (
        (blocks.after > 0 &&
          attachLeaves.has(
            edgeKey(piece.name, leave.pos, leave.reversed, leave.hap),
          )) ||
        (blocks.before > 0 &&
          attachRejoins.has(
            edgeKey(piece.name, rejoin.pos, rejoin.reversed, rejoin.hap),
          ))
      )
    }
    const loaded = (item: Met) =>
      detailsRead.has(item.piece.abridged.detail.name)
    const detailOf = (item: Met) => {
      const { name } = item.piece.abridged.detail
      const unit = detailUnits.get(name)
      if (!unit) {
        throw new WalkFormatError(
          `${refName}'s unit table lists no detail unit ${name}, which an abridged row of ${item.piece.name} names`,
        )
      }
      const ref = refs.get(item.unit) ?? chunkRef(item.tileRows)
      refs.set(item.unit, ref)
      return { unit, ref, piece: item.piece, into: item.into }
    }
    // every abridged row met is read. A walk the cut expands is read from
    // its detail records: those of its pieces over the drawn stretch and of
    // the tiles fetched for own steps, by its abridged rows alone, and those
    // of its other pieces a record attaches steps of. Every other piece is
    // built from its blocks row
    const resolve = async (attaching: boolean) => {
      const batch = met
      met = []
      const details: ReturnType<typeof detailOf>[] = []
      const wantDetail = (item: Met) => {
        const { name } = item.piece.abridged.detail
        if (!detailsRead.has(name)) {
          detailsRead.add(name)
          details.push(detailOf(item))
        }
      }
      for (const item of batch) {
        if (
          expand?.(item.piece.name) === true &&
          (attaching || onDrawn(item.piece.span))
        ) {
          wantDetail(item)
        }
      }
      const blocked = batch.filter(item => !loaded(item))
      const tilesMet = new Map<WalkUnit, number>()
      for (const { unit, from } of blocked) {
        tilesMet.set(unit, Math.min(tilesMet.get(unit) ?? from, from))
      }
      const wanted = [...tilesMet].map(([tileUnit, from]) => {
        const unit = blocksUnits.get(tileUnit.name)
        if (!unit) {
          throw new WalkFormatError(
            `${refName}'s unit table lists no blocks unit ${tileUnit.name}, which the abridged rows of that tile need`,
          )
        }
        return { unit, from }
      })
      if (wanted.length > 0) {
        checkBytes(bytesWith(wanted), attaching)
      }
      const lines = await this.file.read(
        wanted.map(w => w.unit),
        reads.signal,
      )
      const listed = new Map<WalkUnit, Map<string, string[]>>()
      for (const [k, tileUnit] of [...tilesMet.keys()].entries()) {
        extra.push(wanted[k]!)
        const byPiece = new Map<string, string[]>()
        for (const line of lines[k]!) {
          const cols = line.split('\t')
          byPiece.set(
            pieceKey(cols[4] ?? '', cols[5] ?? '', cols[6] ?? ''),
            cols,
          )
        }
        listed.set(tileUnit, byPiece)
      }
      const pieces = blocked.map(item => {
        const { name, fragStart, piece } = item.piece
        const cols = listed
          .get(item.unit)
          ?.get(pieceKey(name, fragStart, piece))
        if (!cols) {
          throw new WalkFormatError(
            `${refName}'s blocks unit ${item.unit.name} holds no row for piece ${piece} of ${name}, which its tile abridges`,
          )
        }
        const blocks = parseBlocks(cols)
        return { ...item, blocks, on: blocksOver(blocks, over) }
      })
      const cores = attaching
        ? new Map<Read, [number, number]>()
        : coresOf(pieces, batch.filter(loaded))
      for (const read of pieces) {
        if (expand?.(read.piece.name) === true && attachNeeds(read)) {
          wantDetail(read)
        }
      }
      await readInOrder(
        details,
        DETAIL_READS,
        reads.signal,
        readDetail,
        decodeDetail,
      )
      const starts = attaching ? pieceStarts(attachRows) : pieceStarts(rows)
      for (const { piece } of batch) {
        starts.set(
          pieceKey(piece.name, piece.fragStart, piece.piece),
          piece.hapOffset,
        )
      }
      for (const read of pieces) {
        const { piece, blocks } = read
        const next = starts.get(
          pieceKey(piece.name, piece.fragStart, piece.piece + 1),
        )
        const pieceEnd =
          piece.hapOffset + blocks.before + blocks.bp + blocks.after
        if (next !== undefined && next !== pieceEnd) {
          throw new WalkFormatError(
            `${refName}'s blocks row of ${piece.name}, piece ${piece.piece} in tile ${read.unit.name}, ends at ${pieceEnd} on its contig, and the line's next piece starts at ${next}`,
          )
        }
        if (!loaded(read)) {
          const built = unloadedRow(
            piece,
            read.blocks,
            read.path,
            over,
            cores.get(read),
            (start, end) => this.runId(piece.name, start, end),
          )
          read.into.push(built.row)
          placeNodes(built.nodes)
          addLinks(built.links)
          ends.push({
            name: piece.name,
            fragStart: piece.fragStart,
            piece: piece.piece,
            first: built.first,
            last: built.last,
          })
          checkNodes({ nodes: built.nodes, rows: [built.row] }, read.unit.end)
        }
      }
    }
    // the blocks units and detail records the tiles decoded need, then the
    // links between pieces
    const finish = async (attaching = false) => {
      await resolve(attaching)
      if (links) {
        for (const [a, b] of pieceLinks(ends)) {
          graph.addLinkEnds(a, b)
        }
      }
    }
    let cut: ReturnType<typeof walkCut>
    try {
      const read = await this.file.read(
        [...tiles, ...departureBlocks],
        reads.signal,
      )
      const { records, ownTiles } = departureRecords(
        read.slice(tiles.length).flat(),
        drawn,
        keep,
      )
      for (const d of records) {
        if (d.ownBp > 0 && d.ownBp <= attachBp) {
          attachLeaves.add(
            edgeKey(d.name, d.leavePos, d.leaveReversed, d.hapStart),
          )
          attachRejoins.add(
            edgeKey(
              d.name,
              d.rejoinPos,
              d.rejoinReversed,
              d.hapStart + d.ownBp,
            ),
          )
        }
      }
      tiles.forEach((unit, k) => {
        decode(unit, read[k]!, rows, unit.start)
      })
      await finish()
      const tooManyRead = stepBudgetError(
        rows,
        stepBudget,
        READ_STEPS_PER_KEPT,
        region,
        tile,
        filtered,
        opts.budgetNames,
      )
      if (tooManyRead) {
        throw tooManyRead
      }
      const cutRows = () =>
        walkCut(joinPieces(rows), graph.nodes, refName, records, {
          contained,
          context,
          attachBp,
          squeeze: squeezer(header.holes, refName),
          drawn: bounds,
          attach: joinPieces(attachRows),
          placeholder: (name, start, end) => this.runId(name, start, end),
        })
      cut = cutRows()
      // the tiles holding own steps the read lacks, the walks wanting each,
      // and where on R the first record wanting it starts
      const fetches = new Map<WalkUnit, { names: Set<string>; from: number }>()
      for (const { departure } of cut.unattached) {
        const { name } = departure
        const from = Math.min(
          departure.leavePos,
          departure.rejoinPos,
          region.end + context - 1,
        )
        for (const unit of attachTiles(units, ownTiles.get(departure) ?? [])) {
          if (!decodedFor(unit, name)) {
            const wanting = fetches.get(unit)
            fetches.set(unit, {
              names: (wanting?.names ?? new Set()).add(name),
              from: Math.min(wanting?.from ?? from, from),
            })
          }
        }
      }
      if (fetches.size > 0) {
        const fetching = [...fetches].map(([unit, { from }]) => ({
          unit,
          from,
        }))
        checkBytes(bytesWith(fetching), true)
        const fetched = await this.file.read([...fetches.keys()], reads.signal)
        for (const [k, [unit, { names, from }]] of [...fetches].entries()) {
          extra.push(fetching[k]!)
          decode(unit, fetched[k]!, attachRows, from, name => names.has(name))
        }
        await finish(true)
        cut = cutRows()
      }
    } catch (error) {
      reads.abort()
      throw error
    } finally {
      signal?.removeEventListener('abort', stop)
    }
    const tooMany =
      keptStepsError(
        cut.fragments,
        cut.departures,
        stepBudget,
        region,
        filtered,
      ) ?? keptNodesError(keptNodes(cut.kept), nodeLimit, region, filtered)
    if (tooMany) {
      throw tooMany
    }
    addLinks(cut.links)
    const { kept, fragments, departures } = cut
    return {
      graph,
      kept,
      fragments,
      departures,
      gaps: cut.gapsOf(),
      placeholders: kept.size - keptNodes(kept),
    }
  }

  /**
   * The graph under a region as tables: the reference sample's sequence over
   * `region`, with `context` bp either side, and the walks of the haplotypes
   * asked for. `region.refName` may be a contig (`chr6`) or the file's whole
   * name; the result names the sequence as the file does.
   */
  async cutRegion(
    region: { refName: string; start: number; end: number },
    opts: RegionCutOptions = {},
  ) {
    const { haplotypes, referenceSample, ...rest } = opts
    const { signal } = opts
    const header = await this.header({ signal })
    const sample = referenceSample ?? header.references[0] ?? ''
    const refName = await this.sequenceFor(sample, region.refName, { signal })
    if (refName === undefined) {
      const names = await this.referenceSequences({ signal })
      throw new Error(
        `The walk file has no sequence ${region.refName} on ${sample}; it holds ${names.slice(0, 5).join(', ')}${names.length > 5 ? ', …' : ''}`,
      )
    }
    const { graph, kept, fragments, departures, gaps } = await this.cut(
      refName,
      region,
      { ...rest, keep: walksKept(haplotypes, refName) },
    )
    return {
      tables: graph.tables(kept, fragments, departures, gaps),
      refName,
      sample,
    }
  }
}
