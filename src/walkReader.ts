import { RemoteFile } from 'generic-filehandle2'

import { refNameLookup, resolveRefName } from './refNames.ts'
import { checkRange } from './tabixSource.ts'
import {
  chunkRef,
  decodeChunk,
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
import {
  WalkFile,
  WalkFormatError,
  once,
  raceAbort,
  unitName,
} from './walkFile.ts'
import {
  WalkGraph,
  byteBudgetError,
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
  // the abridged walks to read in full from their detail records, by
  // `sample#haplotype#contig`; none by default. A detail record is read
  // where its extent overlaps the window and its context, and in a tile
  // fetched for a record's own steps; any other abridged piece of the walk
  // is left out, as the cut has no step of it to keep
  expand?: (name: string) => boolean
  // for checks that hold a cut expanding nothing equal to one expanding
  // every walk: read each expanded piece as far as a cut that leaves it
  // unloaded knows it, without its own steps after its last step on R where
  // the window's tiles lack the next piece, and expand none in a tile
  // fetched for own steps
  expandStated?: boolean
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

const pieceKey = (name: string, fragStart: number, piece: number) =>
  `${name}\t${fragStart}\t${piece}`

// A decoded piece without its steps after its last step on R
function toLastStepOnR(row: WalkRow): WalkRow {
  const { ids, rev, refStarts, against } = row.steps
  let n = row.n
  while (n > 1 && Number.isNaN(refStarts?.[n - 1])) {
    n--
  }
  return {
    ...row,
    n,
    interior: Math.max(0, (row.interior ?? 0) - (row.n - n)),
    steps: {
      ids: ids.subarray(0, n),
      rev: rev.subarray(0, n),
      refStarts: refStarts?.subarray(0, n),
      against: against?.subarray(0, n),
    },
  }
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
 * A `gfa-to-tabix --walks` file of walk format 7, `<prefix>.<sample>.walks.gz`
 * and its directory `<prefix>.<sample>.walks.gz.idx`, read a window at a
 * time. At open the directory, the header and the sequence's unit table; for
 * a window, the tiles and departure blocks overlapping the window and its
 * context, in one read wherever their bytes touch, each tile decoded once;
 * then every tile holding the own steps of records within `attachBp` the
 * read lacks, which give those steps alone, the detail records of the
 * walks `expand` names, and the `blocks` unit of each tile with an abridged
 * piece over the window that the cut does not expand, which gives the
 * piece's steps on the reference and the bp of its own between them. The
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

  // The compressed bytes a cut of `region` reads before any attach fetch,
  // expansion or blocks unit: its tiles and departure blocks
  async bytesFor(
    refName: string,
    region: { start: number; end: number },
    opts: { signal?: AbortSignal; context?: number } = {},
  ) {
    checkRange(refName, region.start, region.end)
    const drawn = drawnSpan(region, opts.context ?? WALK_CONTEXT)
    const units = await this.units(refName, opts)
    return bytesOf(
      units.filter(
        u =>
          (u.kind === 'tile' || u.kind === 'departures') && overlaps(u, drawn),
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
      expandStated = false,
    } = opts
    checkRange(refName, region.start, region.end)
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
    // the detail records of the expanded walks, each once, with the R of
    // their tile and the rows they join; a path's abridged pieces in a tile
    // share one
    let details: { ref: ChunkRef; piece: AbridgedPiece; into: WalkRow[] }[] = []
    const detailsRead = new Set<string>()
    // the abridged pieces over the window that the cut does not expand, by
    // tile, each with its tile's R
    let unloaded = new Map<
      WalkUnit,
      { piece: AbridgedPiece; path: RefPath }[]
    >()
    // where each piece of the window's tiles starts on its contig
    const pieceStarts = new Map<string, number>()
    let placeholders = 0
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
        const over = keptNodesError(
          surely.size,
          nodeLimit,
          region,
          filtered,
          readTo,
        )
        if (over) {
          throw over
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
    // and a cut over its node limit stops there. An abridged piece of a
    // window tile is read from its detail record where the cut expands its
    // walk and built from its blocks row where it does not, in either case
    // only where the record's extent overlaps the drawn stretch; one of a
    // tile fetched for own steps is read where the cut expands its walk
    const decode = (
      unit: WalkUnit,
      lines: string[],
      into: WalkRow[],
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
      const windowed = into === rows
      let ref: ChunkRef | undefined
      for (const row of tileCut.rows) {
        into.push(row)
        if (windowed) {
          pieceStarts.set(
            pieceKey(row.name, row.fragStart, row.piece),
            row.hapOffset,
          )
        }
      }
      for (const piece of tileCut.abridged) {
        const { name } = piece.abridged.detail
        const detail = detailUnits.get(name)
        if (!detail) {
          throw new WalkFormatError(
            `${refName}'s unit table lists no detail unit ${name}, which an abridged row of ${piece.name} names`,
          )
        }
        if (windowed) {
          pieceStarts.set(
            pieceKey(piece.name, piece.fragStart, piece.piece),
            piece.hapOffset,
          )
        }
        if (windowed ? overlaps(detail, drawn) : !expandStated) {
          if (expand?.(piece.name) === true) {
            if (!detailsRead.has(name)) {
              detailsRead.add(name)
              ref ??= chunkRef(tileRows)
              details.push({ ref, piece, into })
            }
          } else if (windowed) {
            unloaded.set(unit, [
              ...(unloaded.get(unit) ?? []),
              { piece, path: tileCut.chunk.path },
            ])
          }
        }
      }
      place(tileCut, unit.end)
    }
    const nextPieceRead = (row: WalkRow) =>
      pieceStarts.has(pieceKey(row.name, row.fragStart, row.piece + 1))
    const readDetail = async ({ piece }: { piece: AbridgedPiece }) => {
      const unit = detailUnits.get(piece.abridged.detail.name)
      const [lines = []] = await this.file.read(
        unit ? [unit] : [],
        reads.signal,
      )
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
          into.push(
            expandStated && !nextPieceRead(decodedRow)
              ? toLastStepOnR(decodedRow)
              : decodedRow,
          )
        }
      }
      place(tileCut, piece.chunkStart + tile)
    }
    // each unloaded piece's row from its tile's blocks unit, a run of its own
    // after its last step on R as long as the next piece's start says
    const buildUnloaded = async () => {
      const pending = [...unloaded]
      unloaded = new Map()
      const blocksOf = pending.map(([tileUnit]) => {
        const unit = blocksUnits.get(tileUnit.name)
        if (!unit) {
          throw new WalkFormatError(
            `${refName}'s unit table lists no blocks unit ${tileUnit.name}, which the abridged rows of that tile need`,
          )
        }
        return unit
      })
      if (pending.length > 0) {
        checkBytes(
          end =>
            tileBytesTo(end) +
            bytesOf(
              blocksOf.filter((_, k) => pending[k]![0].start < end + context),
            ),
        )
      }
      const read = await this.file.read(blocksOf, reads.signal)
      pending.forEach(([tileUnit, pieces], k) => {
        const lists = new Map<string, string>()
        for (const line of read[k]!) {
          const [, , , , name = '', fragStart, piece, list = ''] =
            line.split('\t')
          lists.set(pieceKey(name, Number(fragStart), Number(piece)), list)
        }
        for (const { piece, path } of pieces) {
          const list = lists.get(
            pieceKey(piece.name, piece.fragStart, piece.piece),
          )
          if (list === undefined) {
            throw new WalkFormatError(
              `${unitName(blocksOf[k]!)} holds no row for piece ${piece.piece} of ${piece.name}, which tile ${tileUnit.name} abridges`,
            )
          }
          const next = pieceStarts.get(
            pieceKey(piece.name, piece.fragStart, piece.piece + 1),
          )
          const { row, nodes } = unloadedRow(
            piece,
            list,
            path,
            next === undefined ? 0 : next - piece.hapOffset - piece.abridged.bp,
            () => --placeholders,
          )
          rows.push(row)
          placeNodes(nodes)
          checkNodes({ nodes, rows: [row] }, tileUnit.end)
        }
      })
    }
    // the detail records and blocks units the tiles decoded asked for, then
    // the links between pieces
    const finish = async () => {
      const pending = details
      details = []
      await readInOrder(
        pending,
        DETAIL_READS,
        reads.signal,
        readDetail,
        decodeDetail,
      )
      await buildUnloaded()
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
      tiles.forEach((unit, k) => {
        decode(unit, read[k]!, rows)
      })
      const { records, ownTiles } = departureRecords(
        read.slice(tiles.length).flat(),
        drawn,
        keep,
      )
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
        checkBytes(
          end =>
            tileBytesTo(end) +
            bytesOf(
              [...fetches]
                .filter(([, w]) => w.from < end + context)
                .map(([unit]) => unit),
            ),
          true,
        )
        const fetched = await this.file.read([...fetches.keys()], reads.signal)
        for (const [k, [unit, { names }]] of [...fetches].entries()) {
          decode(unit, fetched[k]!, attachRows, name => names.has(name))
        }
        await finish()
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
      ) ?? keptNodesError(cut.kept.size, nodeLimit, region, filtered)
    if (tooMany) {
      throw tooMany
    }
    const { kept, fragments, departures } = cut
    return { graph, kept, fragments, departures, gaps: cut.gapsOf() }
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
