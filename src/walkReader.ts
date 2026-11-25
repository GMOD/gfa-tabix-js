import { RemoteFile } from 'generic-filehandle2'

import { refNameLookup, resolveRefName } from './refNames.ts'
import { checkRange } from './tabixSource.ts'
import { panSNHaplotypeOrSample } from './util.ts'
import {
  blockLinks,
  blocksOver,
  bpRow,
  canonicalLink,
  chunkRef,
  decodeChunk,
  ownRunsRow,
  ownStretch,
  tileLinks,
  parseBlocks,
  referenceNodes,
  tileNodeIds,
  tileSteps,
  unloadedRow,
} from './walkChunk.ts'
import { WalkTally } from './walkCounts.ts'
import {
  ATTACH_BP,
  attachesOwnSteps,
  departureOwnTiles,
  parseDepartureRow,
  squeezer,
} from './walkDeparture.ts'
import { WalkFile, WalkFormatError, once, raceAbort } from './walkFile.ts'
import {
  WalkGraph,
  byteBudgetError,
  countBytesError,
  countStepsError,
  edgeKey,
  isPlaceholder,
  joinPieces,
  keptNodesError,
  keptStepsError,
  megabytes,
  readSteps,
  stepBudgetError,
  surelyKept,
  surelyKeptError,
  unattachedRecords,
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
  PieceEnds,
  RefPath,
} from './walkChunk.ts'
import type { CountUnit, CountWalk } from './walkCounts.ts'
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
export const WALK_CONTEXT = 1000

// detail records an expanding cut reads at once
const DETAIL_READS = 16

// The budgets a cut keeps to unless its caller passes others: what a browser
// window can hold. Infinity lifts one
export const WALK_BUDGETS = Object.freeze({
  stepBudget: 4_000_000,
  byteBudget: 8_000_000,
  departureByteBudget: 2_000_000,
  attachBp: ATTACH_BP,
})

// The same cut rules in both; a figure has no window to protect
export const CUT_PRESETS = Object.freeze({
  browser: Object.freeze({
    context: WALK_CONTEXT,
    contained: false,
    attachBp: ATTACH_BP,
    stepBudget: WALK_BUDGETS.stepBudget,
    byteBudget: WALK_BUDGETS.byteBudget,
    departureByteBudget: WALK_BUDGETS.departureByteBudget,
  }),
  figure: Object.freeze({
    context: WALK_CONTEXT,
    contained: false,
    attachBp: ATTACH_BP,
    stepBudget: Infinity,
    byteBudget: Infinity,
    departureByteBudget: Infinity,
  }),
}) satisfies Record<string, WalkCutOptions>

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

// The links from each piece's last step to the first of its line's next
const piecesJoined = (ends: PieceEnds[]) => {
  const lines = new Map<string, Map<number, PieceEnds>>()
  for (const piece of ends) {
    const key = `${piece.name}\t${piece.fragStart}`
    const line = lines.get(key) ?? new Map<number, PieceEnds>()
    lines.set(key, line)
    line.set(piece.piece, piece)
  }
  const pairs: number[] = []
  for (const line of lines.values()) {
    for (const [k, piece] of line) {
      const next = line.get(k + 1)
      if (next && piece.last >= 0 && next.first >= 0) {
        pairs.push(piece.last, next.first)
      }
    }
  }
  return pairs
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
  // ATTACH_BP by default, which bounds the nodes one departure brings into a
  // cut. A departure gets its steps where it leaves or rejoins on the window
  // or its context, and the cut fetches every tile holding own steps of such
  // a record that its read lacks, which gives those records their steps and
  // nothing else
  attachBp?: number
  /**
   * @deprecated Ignored: every eligible record's tiles are fetched, their
   * bytes counted against byteBudget
   */
  attachFetches?: number
  context?: number
  // which links the cut lists: 'walked', the default, those its steps take;
  // 'between', also every link a reference row or route of the units it
  // reads takes between two of its nodes; 'none', for a caller that draws no
  // edges
  links?: 'walked' | 'between' | 'none'
  // keep the bases of the nodes, from a file written with --sequences
  sequences?: boolean
  // how a zoom-in notice names the budgets
  budgetNames?: BudgetNames
  // the most nodes a cut keeps; Infinity, the default, has no limit. A
  // placeholder is no node, and is not counted
  nodeLimit?: number
  // the most placeholders a cut keeps, each a node in its tables that a
  // layout places; Infinity, the default, has no limit
  placeholderLimit?: number
  // keep the reference over the window and its context alone, and split each
  // walk where it steps onto the reference off them
  contained?: boolean
  // the abridged walks to read from their detail records, by
  // `sample#haplotype#contig`; none by default. Of such a walk the cut reads
  // the records of the pieces whose steps it keeps or attaches
  expand?: (name: string) => boolean
}

export type WalkCountOptions = Pick<
  WalkCutOptions,
  | 'signal'
  | 'context'
  | 'contained'
  | 'attachBp'
  | 'stepBudget'
  | 'byteBudget'
  | 'departureByteBudget'
  | 'expand'
  | 'budgetNames'
>

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

// The tiles and departure blocks a cut over `drawn` reads first
const windowUnits = (
  units: WalkUnit[],
  drawn: { start: number; end: number },
) => ({
  tiles: units.filter(u => u.kind === 'tile' && overlaps(u, drawn)),
  departureBlocks: units.filter(
    u => u.kind === 'departures' && overlaps(u, drawn),
  ),
})

function checkDepartureBytes(departureBlocks: WalkUnit[], budget: number) {
  const bytes = bytesOf(departureBlocks)
  if (bytes > budget) {
    throw new DepartureReadError(bytes, budget)
  }
}

// The tiles an attach fetch reads for a record: those its `ownTiles` names,
// in file order
export function attachTiles(units: WalkUnit[]) {
  const byStart = new Map<number, number[]>()
  units.forEach((u, k) => {
    const start = Number(u.name)
    if (u.kind === 'tile' && u.length > 0 && !Number.isNaN(start)) {
      const held = byStart.get(start)
      if (held) {
        held.push(k)
      } else {
        byStart.set(start, [k])
      }
    }
  })
  return (ownTiles: number[]) =>
    [...new Set(ownTiles)]
      .flatMap(start => byStart.get(start) ?? [])
      .sort((a, b) => a - b)
      .map(k => units[k]!)
}

const pieceKey = (
  name: string,
  fragStart: number | string,
  piece: number | string,
) => `${name}\t${fragStart}\t${piece}`

const keyOf = (row: { name: string; fragStart: number; piece: number }) =>
  pieceKey(row.name, row.fragStart, row.piece)

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

// A detail unit to read, and the pieces a cut takes from it, each with the
// bp its abridged row states
interface Detail {
  unit: WalkUnit
  ref: ChunkRef
  piece: AbridgedPiece
  into: WalkRow[]
  pieces: Map<string, number>
}

// Where each piece of `rows` starts on its contig
const pieceStarts = (rows: WalkRow[]) => {
  const starts = new Map<string, number>()
  for (const row of rows) {
    starts.set(pieceKey(row.name, row.fragStart, row.piece), row.hapOffset)
  }
  return starts
}

// A piece of a line a cut holds, whether it lies over the drawn stretch, and
// what the cut reads it by
interface LinePiece<T> {
  piece: number
  on: boolean
  item?: T
}

// Of each line's pieces, in piece order, by runs of consecutive piece numbers,
// those from a run's first piece over the drawn stretch to its last, as the
// items they are read by, each with whether it is that first or that last
function betweenDrawn<T>(lines: LinePiece<T>[][]) {
  const between: { item: T; first: boolean; last: boolean }[] = []
  for (const line of lines) {
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
          const { item } = line[j]!
          if (item !== undefined) {
            between.push({ item, first: j === first, last: j === last })
          }
        }
        from = i
      }
    }
  }
  return between
}

function beforeDrawn<T>(lines: LinePiece<T>[][]) {
  const before: T[] = []
  for (const line of lines) {
    for (let i = 0; i + 1 < line.length; i++) {
      const { item, piece } = line[i]!
      const next = line[i + 1]!
      if (item !== undefined && next.on && next.piece === piece + 1) {
        before.push(item)
      }
    }
  }
  return before
}

// A piece of a window's tile the counts may read: a walk row, or an abridged
// row
type CountPiece = { unit: CountUnit; tile: WalkUnit } & (
  | { walk: CountWalk; moved?: undefined }
  | { moved: AbridgedPiece; walk?: undefined }
)

const split = (lines: string[]) => lines.map(line => line.split('\t'))

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
 * the `blocks` unit of each of those tiles that abridges a walk kept and not
 * expanded; the detail records of the pieces of the walks `expand` names
 * that the cut builds; then every tile holding own steps the read lacks of
 * a record within `attachBp` that leaves or rejoins on the window or its
 * context, which give those steps alone, from the tile's rows and the
 * records, with no `blocks` unit read. A record that only jumps over them
 * costs no read past its line. The table gives
 * each read's bytes before it is made; a window whose tiles are past the
 * byte budget reads no row, and one whose attach tiles take it past reads
 * none of them. Only the walks `keep` accepts are decoded. The step budget
 * bounds the steps the cut keeps, and rows of the window holding
 * READ_STEPS_PER_KEPT times as many stop the read.
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

  /**
   * The compressed bytes a cut of `region` reads when it keeps every walk
   * and expands none: its tiles, their blocks units, its departure blocks,
   * which this reads as the cut does, and the tiles named for their own
   * steps by its records within `attachBp` that leave or rejoin on the
   * window or its context. A cut that keeps fewer walks reads no more, and
   * one that expands walks no more but for their detail units
   */
  async bytesFor(
    refName: string,
    region: { start: number; end: number },
    opts: {
      signal?: AbortSignal
      context?: number
      attachBp?: number
      departureByteBudget?: number
    } = {},
  ) {
    const {
      signal,
      context = WALK_CONTEXT,
      attachBp = ATTACH_BP,
      departureByteBudget = WALK_BUDGETS.departureByteBudget,
    } = opts
    checkRange(refName, region.start, region.end)
    const drawn = drawnSpan(region, context)
    const units = await this.units(refName, { signal })
    const { tiles, departureBlocks } = windowUnits(units, drawn)
    checkDepartureBytes(departureBlocks, departureByteBudget)
    const tileNames = new Set(tiles.map(u => u.name))
    const { records, ownTiles } = departureRecords(
      (await this.file.read(departureBlocks, signal)).flat(),
      drawn,
      undefined,
    )
    const attachNames = new Set(
      records
        .filter(d => attachesOwnSteps(d, drawn, attachBp))
        .flatMap(d => ownTiles.get(d) ?? []),
    )
    const blocks = units.filter(
      u => u.kind === 'blocks' && tileNames.has(u.name),
    )
    const attached = units.filter(
      u =>
        u.kind === 'tile' &&
        u.length > 0 &&
        attachNames.has(Number(u.name)) &&
        !tileNames.has(u.name),
    )
    return (
      bytesOf(tiles) +
      bytesOf(blocks) +
      bytesOf(departureBlocks) +
      bytesOf(attached)
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
      links = 'walked',
      contained = false,
      nodeLimit = Infinity,
      placeholderLimit = Infinity,
      expand,
    } = opts
    checkRange(refName, region.start, region.end)
    const [header, units] = await Promise.all([
      this.header({ signal }),
      this.units(refName, { signal }),
    ])
    const { tile } = header
    const drawn = drawnSpan(region, context)
    const { tiles, departureBlocks } = windowUnits(units, drawn)
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
    checkDepartureBytes(departureBlocks, departureByteBudget)
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
    const graph = new WalkGraph(
      {
        refName,
        start: region.start,
        end: region.end,
        context,
        sequences: opts.sequences,
      },
      header.revisited,
    )
    // the name filters each tile was decoded under, by its offset
    const decodedBy = new Map<number, ((name: string) => boolean)[]>()
    const decodedFor = (unit: WalkUnit, name: string) =>
      decodedBy.get(unit.offset)?.some(took => took(name)) === true
    const surely = new Set<number>()
    const surelyPlaceholders = new Set<number>()
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
    // R of each tile decoded, by its name, and the blocks units read
    const tilePaths = new Map<string, { path: RefPath; rowStarts: number[] }>()
    const pieceEnds: PieceEnds[] = []
    const blocksRead = new Map<WalkUnit, string[]>()
    // the id of each run of a contig a placeholder stands for, by the walk's
    // name, the run's start and its end: -1 for the first run, then -2
    const placeholderIds = new Map<string, Map<number, Map<number, number>>>()
    let minted = 0
    const placeholderId = (name: string, start: number, end: number) => {
      const byStart =
        placeholderIds.get(name) ?? new Map<number, Map<number, number>>()
      placeholderIds.set(name, byStart)
      const byEnd = byStart.get(start) ?? new Map<number, number>()
      byStart.set(start, byEnd)
      const id = byEnd.get(end) ?? -(minted + 1)
      minted = Math.max(minted, -id)
      byEnd.set(end, id)
      return id
    }
    // the abridged rows decoded and not yet read from a blocks row or a
    // detail record
    let met: Met[] = []
    const refs = new Map<WalkUnit, ChunkRef>()
    // the drawn stretch, as far as the reference's passes over it reach
    const over = { ...drawn }
    const ownRecords = new Map<string, WalkDeparture[]>()
    const checkNodes = (
      chunk: { nodes: Map<number, unknown>; rows: WalkRow[] },
      readTo: number,
    ) => {
      if (Number.isFinite(nodeLimit) || Number.isFinite(placeholderLimit)) {
        surelyKept(chunk, graph.nodes, bounds, surely, surelyPlaceholders)
        if (
          surely.size > nodeLimit ||
          surelyPlaceholders.size > placeholderLimit
        ) {
          surely.clear()
          surelyPlaceholders.clear()
          surelyKept(
            { nodes: graph.nodes, rows },
            graph.nodes,
            bounds,
            surely,
            surelyPlaceholders,
          )
        }
        const tooMany =
          keptNodesError(surely.size, nodeLimit, region, filtered, readTo) ??
          keptNodesError(
            surelyPlaceholders.size,
            placeholderLimit,
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
      if (links !== 'none') {
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
      const { path } = tileCut
      if (links === 'between' && !tilePaths.has(unit.name)) {
        const taken = tileLinks(tileRows, path)
        tilePaths.set(unit.name, { path, rowStarts: taken.rowStarts })
        addLinks(taken.links)
        for (const ends of taken.ends) {
          pieceEnds.push(ends)
        }
      }
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
    const readDetail = async ({ unit }: Detail) => {
      const [lines = []] = await this.file.read([unit], reads.signal)
      return lines.map(line => line.split('\t'))
    }
    const decodeDetail = (
      { unit, ref, piece, into, pieces }: Detail,
      held: string[][],
    ) => {
      const tileCut = decodeChunk(
        held,
        piece.chunkStart,
        name => name === piece.name,
        opts.sequences,
        ref,
      )
      const taken = tileCut.rows.filter(row => pieces.has(keyOf(row)))
      if (taken.length !== pieces.size) {
        throw new WalkFormatError(
          `${refName}'s detail unit ${unit.name} holds walk rows for ${taken.length} of the ${pieces.size} pieces of ${piece.name} that the cut reads from it`,
        )
      }
      for (const row of taken) {
        const { bp } = row.steps
        const rowBp = bp?.reduce((sum, stepBp) => sum + stepBp, 0)
        const hapBp = pieces.get(keyOf(row))
        if (bp?.every(stepBp => stepBp >= 0) === true && rowBp !== hapBp) {
          throw new WalkFormatError(
            `${refName}'s detail unit ${unit.name} holds ${rowBp} bp of steps for piece ${row.piece} of ${row.name}, and its abridged row states a hapBp of ${hapBp}`,
          )
        }
        into.push(row)
      }
      place({ ...tileCut, rows: taken }, piece.chunkStart + tile)
    }
    const onDrawn = (span: { start: number; end: number }) =>
      span.start < over.end && span.end > over.start
    // the lines `held` names a piece of, each with its pieces among `held`
    // and the rows decoded, in piece order
    const linesOf = <T>(
      held: (LinePiece<T> & { name: string; fragStart: number })[],
    ) => {
      const lines = new Map<string, LinePiece<T>[]>()
      for (const { name, fragStart, ...piece } of held) {
        const key = `${name}\t${fragStart}`
        const line = lines.get(key) ?? []
        lines.set(key, line)
        line.push(piece)
      }
      for (const row of rows) {
        lines.get(`${row.name}\t${row.fragStart}`)?.push({
          piece: row.piece,
          on:
            row.span === undefined ||
            onDrawn({ start: row.span[0], end: row.span[1] }),
        })
      }
      return [...lines.values()].map(line =>
        line.sort((a, b) => a.piece - b.piece),
      )
    }
    const bySpan = <T extends Met>(item: T) => {
      const { name, fragStart, piece, span } = item.piece
      return { name, fragStart, piece, on: onDrawn(span), item }
    }
    // each built piece's blocks to build whole: those from the first block
    // over the drawn stretch to the last, of the pieces that join it
    const coresOf = (pieces: Read[], others: Met[]) => {
      const cores = new Map<Read, [number, number]>()
      for (const { item, first, last } of betweenDrawn(
        linesOf([
          ...pieces.map(read => ({ ...bySpan(read), on: read.on.first >= 0 })),
          ...others.map(other => ({ ...bySpan(other), item: undefined })),
        ]),
      )) {
        cores.set(item, [
          first ? item.on.first : 0,
          last ? item.on.last : item.blocks.begins.length - 1,
        ])
      }
      return cores
    }
    const expanding = (item: Met) => expand?.(item.piece.name) === true
    const ownRunsOf = (piece: AbridgedPiece) =>
      (ownRecords.get(piece.name) ?? []).flatMap(d => {
        const stretch = ownStretch(piece, d)
        return stretch
          ? [
              {
                ...stretch,
                attached: attachesOwnSteps(d, drawn, attachBp),
                pos: stretch.trailing ? d.leavePos : d.rejoinPos,
                reversed: stretch.trailing ? d.leaveReversed : d.rejoinReversed,
                key: stretch.trailing
                  ? edgeKey(d.name, d.leavePos, d.leaveReversed, d.hapStart)
                  : edgeKey(
                      d.name,
                      d.rejoinPos,
                      d.rejoinReversed,
                      d.hapStart + d.ownBp,
                    ),
              },
            ]
          : []
      })
    const blocksOf = async (items: Met[], attaching: boolean) => {
      const tilesMet = new Map<WalkUnit, number>()
      for (const { unit, from } of items) {
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
      const unread = wanted.filter(w => !blocksRead.has(w.unit))
      if (unread.length > 0) {
        checkBytes(bytesWith(unread), attaching)
      }
      const fetched = await this.file.read(
        unread.map(w => w.unit),
        reads.signal,
      )
      unread.forEach((w, k) => {
        extra.push(w)
        blocksRead.set(w.unit, fetched[k]!)
      })
      const listed = new Map<WalkUnit, Map<string, string[]>>()
      for (const [k, tileUnit] of [...tilesMet.keys()].entries()) {
        const byPiece = new Map<string, string[]>()
        for (const line of blocksRead.get(wanted[k]!.unit)!) {
          const cols = line.split('\t')
          byPiece.set(
            pieceKey(cols[4] ?? '', cols[5] ?? '', cols[6] ?? ''),
            cols,
          )
        }
        listed.set(tileUnit, byPiece)
      }
      return new Map(
        items.map(item => {
          const { name, piece, hapBp } = item.piece
          const cols = listed.get(item.unit)?.get(keyOf(item.piece))
          if (!cols) {
            throw new WalkFormatError(
              `${refName}'s blocks unit ${item.unit.name} holds no row for piece ${piece} of ${name}, which its tile abridges`,
            )
          }
          const blocks = parseBlocks(cols)
          const listBp = blocks.before + blocks.bp + blocks.after
          if (listBp !== hapBp) {
            throw new WalkFormatError(
              `${refName}'s blocks row of ${name}, piece ${piece} in tile ${item.unit.name}, adds up to ${listBp} bp, and its abridged row states a hapBp of ${hapBp}`,
            )
          }
          return [item, { ...item, blocks, on: blocksOver(blocks, over) }]
        }),
      )
    }
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
      return {
        unit,
        ref,
        piece: item.piece,
        into: item.into,
        pieces: new Map<string, number>(),
      }
    }
    // The cut builds the pieces whose span lies over the drawn stretch, those
    // between two such pieces of a line, those holding the own steps of a
    // record it attaches and, contained with context, those whose own
    // steps the context reaches before a piece over the drawn stretch: an
    // expanded walk's from its detail unit, any other from its tile's blocks
    // unit. Of a tile fetched to attach it builds the pieces holding those
    // own steps, an unexpanded one from its abridged row and the record. Any
    // other piece of a window's tile is its bp
    const resolve = async (attaching: boolean) => {
      const batch = met
      met = []
      const lines = attaching ? [] : linesOf(batch.map(bySpan))
      const drawnPieces = new Set(betweenDrawn(lines).map(held => held.item))
      const reached = new Set(
        contained && context > 0 && !attaching ? beforeDrawn(lines) : [],
      )
      const ownRuns = new Map(batch.map(item => [item, ownRunsOf(item.piece)]))
      const built = batch.filter(item => {
        const runs = ownRuns.get(item) ?? []
        return (
          drawnPieces.has(item) ||
          runs.some(run => run.attached) ||
          (reached.has(item) && runs.some(run => run.trailing))
        )
      })
      const details = new Map<string, Detail>()
      for (const item of built.filter(expanding)) {
        const { name } = item.piece.abridged.detail
        const detail = details.get(name) ?? detailOf(item)
        details.set(name, detail)
        detail.pieces.set(keyOf(item.piece), item.piece.hapBp)
      }
      const blocks = await blocksOf(
        attaching ? [] : batch.filter(item => !expanding(item)),
        attaching,
      )
      const fromBlocks = built.flatMap(item => {
        const read = blocks.get(item)
        return read ? [{ read, runs: ownRuns.get(item) ?? [] }] : []
      })
      const cores = coresOf(
        fromBlocks.map(({ read }) => read),
        built.filter(expanding),
      )
      await readInOrder(
        [...details.values()],
        DETAIL_READS,
        reads.signal,
        readDetail,
        decodeDetail,
      )
      const starts = attaching ? pieceStarts(attachRows) : pieceStarts(rows)
      for (const { piece } of batch) {
        starts.set(keyOf(piece), piece.hapOffset)
      }
      for (const { piece, unit } of batch) {
        const next = starts.get(
          pieceKey(piece.name, piece.fragStart, piece.piece + 1),
        )
        const pieceEnd = piece.hapOffset + piece.hapBp
        if (next !== undefined && next !== pieceEnd) {
          throw new WalkFormatError(
            `${refName}'s abridged row of ${piece.name}, piece ${piece.piece} in tile ${unit.name}, ends at ${pieceEnd} on its contig, and the line's next piece starts at ${next}`,
          )
        }
      }
      for (const { read, runs } of fromBlocks) {
        const { piece } = read
        for (const run of runs) {
          const listed = run.trailing ? read.blocks.after : read.blocks.before
          if (listed !== run.end - run.start) {
            throw new WalkFormatError(
              `${refName}'s blocks row of ${piece.name}, piece ${piece.piece} in tile ${read.unit.name}, ${run.trailing ? 'ends' : 'starts'} with ${listed} bp off the reference, and a departure record's own steps there are ${run.end - run.start} bp by its abridged row`,
            )
          }
        }
        const made = unloadedRow(
          piece,
          read.blocks,
          read.path,
          over,
          cores.get(read),
          (start, end) => placeholderId(piece.name, start, end),
        )
        read.into.push(made.row)
        placeNodes(made.nodes)
        checkNodes({ nodes: made.nodes, rows: [made.row] }, read.unit.end)
      }
      const unbuilt = new Set(batch)
      for (const item of built) {
        unbuilt.delete(item)
        if (attaching && !expanding(item)) {
          const { piece } = item
          const made = ownRunsRow(
            piece,
            (ownRuns.get(item) ?? []).filter(run => run.attached),
            item.path,
            (start, end) => placeholderId(piece.name, start, end),
          )
          item.into.push(made.row)
          placeNodes(made.nodes)
        }
      }
      for (const item of attaching ? [] : unbuilt) {
        item.into.push(bpRow(item.piece))
      }
    }
    // whether `more` read with the units read so far passes the byte budget
    const overBytes = (more: Needed[], attach: boolean) =>
      Number.isFinite(byteBudget) &&
      byteBudgetError(
        bytesWith(more),
        byteBudget,
        region,
        tile,
        opts.budgetNames,
        attach,
      ) !== undefined
    // Whether none of the refusals a cut makes before keptStepsError can
    // come, from its tiles' lines and records before any is decoded: a row
    // its decode refuses; its node limit, which no more nodes than the
    // tiles' ref rows and routes name can reach where no detail record or
    // attach fetch brings more; its placeholder limit, which nothing reaches
    // where no tile abridges a walk it keeps unexpanded and no attach fetch
    // comes; the byte budget, over the blocks units it reads and every tile
    // it can fetch to attach; and the read budget, over the most its rows
    // and pieces can hold
    const noEarlierRefusal = (
      read: string[][],
      counted: (ReturnType<typeof tileSteps> | 'edge')[],
      records: WalkDeparture[],
      ownTiles: Map<WalkDeparture, number[]>,
    ) => {
      const all = tiles.map((_, k) => {
        const tally = counted[k]
        return tally === 'edge' ? tileSteps(read[k]!, keep, expand) : tally
      })
      const blocks = tiles
        .filter((_, k) => all[k]?.blocks === true)
        .map(unit => ({ unit: blocksUnits.get(unit.name), from: unit.start }))
      const unread = blocks.flatMap(({ unit, from }) =>
        unit && !blocksRead.has(unit) ? [{ unit, from }] : [],
      )
      const attaching = new Map<WalkUnit, number>()
      const decoded = new Set(tiles.map(unit => unit.offset))
      const tilesNamed = new Map<number, WalkUnit[]>()
      for (const unit of units) {
        if (unit.kind === 'tile' && unit.length > 0) {
          const start = Number(unit.name)
          tilesNamed.set(start, [...(tilesNamed.get(start) ?? []), unit])
        }
      }
      for (const d of records) {
        if (attachesOwnSteps(d, bounds, attachBp)) {
          const from = Math.min(
            d.leavePos,
            d.rejoinPos,
            region.end + context - 1,
          )
          for (const start of ownTiles.get(d) ?? []) {
            for (const unit of tilesNamed.get(start) ?? []) {
              if (!decoded.has(unit.offset)) {
                attaching.set(unit, Math.min(attaching.get(unit) ?? from, from))
              }
            }
          }
        }
      }
      const fetching = [...attaching].map(([unit, from]) => ({ unit, from }))
      const nodesFit = () =>
        !Number.isFinite(nodeLimit) ||
        (fetching.length === 0 &&
          all.every(t => t?.expanded === false) &&
          read
            .slice(0, tiles.length)
            .reduce((ids, lines) => tileNodeIds(lines, ids), new Set<number>())
            .size <= nodeLimit)
      const placeholdersFit =
        !Number.isFinite(placeholderLimit) ||
        (fetching.length === 0 && all.every(t => t?.blocks === false))
      const held = all.reduce(
        (sum, t) => ({
          n: sum.n + (t?.n ?? 0),
          interior: sum.interior + (t?.interior ?? 0),
        }),
        { n: 0, interior: 0 },
      )
      return (
        all.every(t => t !== undefined) &&
        blocks.every(({ unit }) => unit !== undefined) &&
        readSteps(held) <= stepBudget * READ_STEPS_PER_KEPT &&
        !(unread.length > 0 && overBytes(unread, false)) &&
        !(fetching.length > 0 && overBytes([...unread, ...fetching], true)) &&
        placeholdersFit &&
        nodesFit()
      )
    }
    // keptStepsError's notice, given before any tile is decoded where the
    // steps the tiles strictly inside the drawn stretch surely keep
    // (tileSteps in walkChunk.ts) pass the step budget and no earlier
    // refusal can come, so the finished cut would throw it
    const earlyRefusal = (
      read: string[][],
      records: WalkDeparture[],
      ownTiles: Map<WalkDeparture, number[]>,
    ) => {
      const counted = tiles.map((unit, k) =>
        unit.start > drawn.start && unit.end < drawn.end
          ? tileSteps(read[k]!, keep, expand)
          : 'edge',
      )
      const kept = counted.reduce(
        (sum, tally) => sum + (tally === 'edge' ? 0 : (tally?.kept ?? 0)),
        0,
      )
      return kept > stepBudget &&
        noEarlierRefusal(read, counted, records, ownTiles)
        ? surelyKeptError(
            tiles.map((unit, k) => {
              const tally = counted[k]
              return tally === 'edge' ? unit : { ...unit, steps: tally?.kept }
            }),
            stepBudget,
            region,
            context,
            filtered,
          )
        : undefined
    }
    let cut: ReturnType<typeof walkCut>
    try {
      const jumps =
        links === 'between'
          ? tiles.flatMap(windowTile => {
              const unit = blocksUnits.get(windowTile.name)
              return unit ? [{ unit, from: windowTile.start, windowTile }] : []
            })
          : []
      checkBytes(bytesWith(jumps))
      const read = await this.file.read(
        [...tiles, ...departureBlocks, ...jumps.map(x => x.unit)],
        reads.signal,
      )
      const tilesAndRecords = tiles.length + departureBlocks.length
      jumps.forEach((x, k) => {
        extra.push(x)
        blocksRead.set(x.unit, read[tilesAndRecords + k]!)
      })
      const { records, ownTiles } = departureRecords(
        read.slice(tiles.length, tilesAndRecords).flat(),
        drawn,
        keep,
      )
      for (const d of records) {
        if (d.ownBp > 0) {
          const held = ownRecords.get(d.name) ?? []
          ownRecords.set(d.name, held)
          held.push(d)
        }
      }
      const early = Number.isFinite(stepBudget)
        ? earlyRefusal(read, records, ownTiles)
        : undefined
      if (early) {
        throw early
      }
      tiles.forEach((unit, k) => {
        decode(unit, read[k]!, rows, unit.start)
      })
      for (const { unit, windowTile } of jumps) {
        const { path, rowStarts } = tilePaths.get(windowTile.name)!
        for (const line of blocksRead.get(unit)!) {
          const piece = blockLinks(line.split('\t'), path, rowStarts)
          addLinks(piece.links)
          pieceEnds.push(piece.ends)
        }
      }
      await resolve(false)
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
      const joined = joinPieces(rows)
      const cutOptions = {
        contained,
        context,
        attachBp,
        squeeze: squeezer(header.holes, refName),
        drawn: bounds,
        placeholder: placeholderId,
        links: links !== 'none',
      }
      // the tiles holding own steps the read lacks, the walks wanting each,
      // and where on R the first record wanting it starts
      const fetches = new Map<WalkUnit, { names: Set<string>; from: number }>()
      const tilesOwning = attachTiles(units)
      for (const departure of unattachedRecords(
        joined,
        graph.nodes,
        records,
        cutOptions,
      )) {
        const { name } = departure
        const from = Math.min(
          departure.leavePos,
          departure.rejoinPos,
          region.end + context - 1,
        )
        for (const unit of tilesOwning(ownTiles.get(departure) ?? [])) {
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
        await resolve(true)
      }
      cut = walkCut(joined, graph.nodes, refName, records, {
        ...cutOptions,
        attach: joinPieces(attachRows),
      })
    } catch (error) {
      reads.abort()
      throw error
    } finally {
      signal?.removeEventListener('abort', stop)
    }
    const nodes = keptNodes(cut.kept)
    const placeholders = cut.kept.size - nodes
    const tooMany =
      keptStepsError(
        cut.fragments,
        cut.departures,
        stepBudget,
        region,
        filtered,
      ) ??
      keptNodesError(nodes, nodeLimit, region, filtered) ??
      keptNodesError(placeholders, placeholderLimit, region, filtered)
    if (tooMany) {
      throw tooMany
    }
    addLinks(cut.links)
    addLinks(piecesJoined(pieceEnds))
    const { kept, fragments, departures } = cut
    return {
      graph,
      kept,
      fragments,
      departures,
      get gaps() {
        return cut.gapsOf()
      },
      placeholders,
    }
  }

  /**
   * How many of the file's haplotypes step each node and link over `region`,
   * and how often, read as a cut that keeps every walk reads, with the same
   * `context`, `contained`, `attachBp`, budgets and `expand`. It counts each
   * piece the cut builds whole, and each piece holding the own steps of a
   * record the cut attaches, while it decodes the rows, keeping no step and
   * building no graph. Every count states why it may fall short, or that it
   * is exact (README, Counts)
   */
  async counts(
    refName: string,
    region: { start: number; end: number },
    opts: WalkCountOptions = {},
  ) {
    const {
      signal,
      context = WALK_CONTEXT,
      contained = false,
      stepBudget = WALK_BUDGETS.stepBudget,
      byteBudget = WALK_BUDGETS.byteBudget,
      departureByteBudget = WALK_BUDGETS.departureByteBudget,
      attachBp = ATTACH_BP,
      expand,
    } = opts
    checkRange(refName, region.start, region.end)
    const [header, units, sequences] = await Promise.all([
      this.header({ signal }),
      this.units(refName, { signal }),
      this.referenceSequences({ signal }),
    ])
    const drawn = drawnSpan(region, context)
    const { tiles, departureBlocks } = windowUnits(units, drawn)
    const unitsOf = (kind: string) =>
      new Map(units.filter(u => u.kind === kind).map(u => [u.name, u]))
    const blocksUnits = unitsOf('blocks')
    const detailUnits = unitsOf('detail')
    const extra: Needed[] = []
    const checkBytes = (more: Needed[], attaching = false) => {
      const bytesTo = (end: number) =>
        bytesOf(tiles.filter(u => u.start < end + context)) +
        bytesOf(extra.filter(x => x.from < end + context).map(x => x.unit)) +
        bytesOf(more.filter(x => x.from < end + context).map(x => x.unit))
      const tooLarge = Number.isFinite(byteBudget)
        ? countBytesError(
            bytesTo,
            byteBudget,
            region,
            header.tile,
            opts.budgetNames,
            attaching,
          )
        : undefined
      if (tooLarge) {
        throw tooLarge
      }
      for (const x of more) {
        extra.push(x)
      }
    }
    const blocksOf = (tileUnits: WalkUnit[]) =>
      tileUnits.flatMap(tileUnit => {
        const unit = blocksUnits.get(tileUnit.name)
        return unit ? [{ unit, from: tileUnit.start }] : []
      })
    const prefetched = expand === undefined ? blocksOf(tiles) : []
    checkBytes(prefetched)
    checkDepartureBytes(departureBlocks, departureByteBudget)
    const read = await this.file.read(
      [...tiles, ...departureBlocks, ...prefetched.map(x => x.unit)],
      signal,
    )
    const { records, ownTiles } = departureRecords(
      read.slice(tiles.length, tiles.length + departureBlocks.length).flat(),
      drawn,
      undefined,
    )
    const listed = new Map(
      prefetched.map((x, k) => [
        x.unit,
        read[tiles.length + departureBlocks.length + k]!,
      ]),
    )
    const tally = new WalkTally(
      [
        ...new Set([
          ...sequences.map(panSNHaplotypeOrSample),
          ...header.haplotypes,
        ]),
      ],
      header.revisited === 0,
    )
    const parsed = tiles.map((unit, k) =>
      tally.tile(split(read[k]!), Number(unit.name), true),
    )
    const over = { ...drawn }
    for (const unit of parsed) {
      const hull = unit.r.path.hull(drawn)
      over.start = Math.min(over.start, hull?.start ?? over.start)
      over.end = Math.max(over.end, hull?.end ?? over.end)
    }
    for (const unit of parsed) {
      tally.placeReference(unit, over)
    }
    const onOver = (span: { start: number; end: number }) =>
      span.start < over.end && span.end > over.start
    const attached = records.filter(d => attachesOwnSteps(d, drawn, attachBp))
    const ownOf = new Map<string, WalkDeparture[]>()
    for (const d of attached) {
      const held = ownOf.get(d.name) ?? []
      ownOf.set(d.name, held)
      held.push(d)
    }
    const holdsOwn = (name: string, start: number, end: number) =>
      (ownOf.get(name) ?? []).some(
        d => start < d.hapStart + d.ownBp && end > d.hapStart,
      )
    const pieces: CountPiece[] = parsed.flatMap((unit, k) => [
      ...unit.walks.map(walk => ({ unit, tile: tiles[k]!, walk })),
      ...unit.abridged.map(moved => ({ unit, tile: tiles[k]!, moved })),
    ])
    const lines = new Map<string, LinePiece<CountPiece>[]>()
    for (const item of pieces) {
      const { name, fragStart, piece, span } = item.walk ?? item.moved
      const key = pieceKey(name, fragStart, '')
      const line = lines.get(key) ?? []
      lines.set(key, line)
      line.push({ piece, on: onOver(span), item })
    }
    const inOrder = [...lines.values()].map(line =>
      line.sort((a, b) => a.piece - b.piece),
    )
    const chosen = new Set(betweenDrawn(inOrder).map(held => held.item))
    for (const item of contained && context > 0 ? beforeDrawn(inOrder) : []) {
      chosen.add(item)
    }
    for (const item of pieces) {
      const { name, hapOffset } = item.walk ?? item.moved
      if (ownOf.has(name) && !chosen.has(item)) {
        const bp = item.walk
          ? tally.walk(item.unit, item.walk, false)
          : item.moved.hapBp
        if (holdsOwn(name, hapOffset, hapOffset + bp)) {
          chosen.add(item)
        }
      }
    }
    const rowSteps = parsed.flatMap(unit =>
      unit.walks.map(({ n, interior, name }) => ({
        n,
        interior,
        name,
        chunkStart: unit.start,
      })),
    )
    const tooManyRead = countStepsError(
      rowSteps,
      stepBudget,
      READ_STEPS_PER_KEPT,
      region,
      header.tile,
      opts.budgetNames,
    )
    if (tooManyRead) {
      throw tooManyRead
    }
    let stepsRead = rowSteps.reduce((sum, row) => sum + readSteps(row), 0)
    let stopped = false
    // the detail units of the pieces `expand` names that fit the step budget,
    // each with the pieces it gives; the rest stay unloaded
    const detailsOf = (wanted: { unit: CountUnit; moved: AbridgedPiece }[]) => {
      const details = new Map<
        string,
        { unit: WalkUnit; tile: CountUnit; pieces: Map<string, AbridgedPiece> }
      >()
      for (const { unit, moved: abridged } of wanted) {
        const { name } = abridged.abridged.detail
        const detail = detailUnits.get(name)
        if (!detail) {
          throw new WalkFormatError(
            `${refName}'s unit table lists no detail unit ${name}, which an abridged row of ${abridged.name} names`,
          )
        }
        const held = details.get(name) ?? {
          unit: detail,
          tile: unit,
          pieces: new Map<string, AbridgedPiece>(),
        }
        details.set(name, held)
        held.pieces.set(keyOf(abridged), abridged)
      }
      return [...details.values()].filter(({ pieces: given }) => {
        const steps = [...given.values()].reduce(
          (sum, piece) => sum + piece.abridged.steps,
          0,
        )
        stopped ||= stepsRead + steps > stepBudget * READ_STEPS_PER_KEPT
        stepsRead += stopped ? 0 : steps
        return !stopped
      })
    }
    const countDetails = (
      details: ReturnType<typeof detailsOf>,
      unitLines: string[][],
    ) => {
      details.forEach(({ unit, tile, pieces: given }, k) => {
        const held = tally.detail(split(unitLines[k]!), tile)
        const taken = held.walks.filter(walk => given.has(keyOf(walk)))
        if (taken.length !== given.size) {
          throw new WalkFormatError(
            `${refName}'s detail unit ${unit.name} holds walk rows for ${taken.length} of the ${given.size} pieces that the counts read from it`,
          )
        }
        for (const walk of taken) {
          const bp = tally.walk(held, walk, true)
          const hapBp = given.get(keyOf(walk))!.hapBp
          if (bp !== hapBp) {
            throw new WalkFormatError(
              `${refName}'s detail unit ${unit.name} holds ${bp} bp of steps for piece ${walk.piece} of ${walk.name}, and its abridged row states a hapBp of ${hapBp}`,
            )
          }
        }
      })
      return new Set(
        details.flatMap(({ pieces: given }) => [...given.values()]),
      )
    }
    const abridgedChosen = pieces.flatMap(item =>
      item.moved && chosen.has(item)
        ? [{ unit: item.unit, tile: item.tile, moved: item.moved }]
        : [],
    )
    const details = detailsOf(
      abridgedChosen.filter(item => expand?.(item.moved.name) === true),
    )
    const expanded = new Set(
      details.flatMap(({ pieces: given }) => [...given.values()]),
    )
    const fromBlocks = abridgedChosen.filter(item => !expanded.has(item.moved))
    const blocksWanted =
      expand === undefined
        ? []
        : blocksOf([...new Set(fromBlocks.map(item => item.tile))])
    checkBytes(blocksWanted)
    signal?.throwIfAborted()
    const more = await this.file.read(
      [...details.map(d => d.unit), ...blocksWanted.map(x => x.unit)],
      signal,
    )
    blocksWanted.forEach((x, k) => {
      listed.set(x.unit, more[details.length + k]!)
    })
    const starts = new Map<string, number>()
    for (const item of pieces) {
      const held = item.walk ?? item.moved
      starts.set(keyOf(held), held.hapOffset)
    }
    for (const { moved: abridged, tile } of pieces.flatMap(item =>
      item.moved ? [{ moved: item.moved, tile: item.tile }] : [],
    )) {
      const next = starts.get(
        pieceKey(abridged.name, abridged.fragStart, abridged.piece + 1),
      )
      const pieceEnd = abridged.hapOffset + abridged.hapBp
      if (next !== undefined && next !== pieceEnd) {
        throw new WalkFormatError(
          `${refName}'s abridged row of ${abridged.name}, piece ${abridged.piece} in tile ${tile.name}, ends at ${pieceEnd} on its contig, and the line's next piece starts at ${next}`,
        )
      }
    }
    for (const unit of parsed) {
      for (const reference of unit.references) {
        if (
          reference.to > reference.from &&
          onOver({
            start: unit.r.path.starts[reference.from]!,
            end: unit.r.path.end(reference.to - 1),
          })
        ) {
          tally.countReference(unit, reference)
        }
      }
    }
    for (const item of pieces) {
      if (item.walk && chosen.has(item)) {
        tally.walk(item.unit, item.walk, true)
      }
    }
    countDetails(details, more)
    const blocksRows = new Map<WalkUnit, Map<string, string[]>>()
    for (const { moved: abridged, tile, unit } of fromBlocks) {
      const blocksUnit = blocksUnits.get(tile.name)
      const byPiece =
        blocksRows.get(tile) ??
        new Map(
          split(blocksUnit ? (listed.get(blocksUnit) ?? []) : []).map(cols => [
            pieceKey(cols[4] ?? '', cols[5] ?? '', cols[6] ?? ''),
            cols,
          ]),
        )
      blocksRows.set(tile, byPiece)
      const cols = byPiece.get(keyOf(abridged))
      if (!cols) {
        throw new WalkFormatError(
          `${refName}'s blocks unit ${tile.name} holds no row for piece ${abridged.piece} of ${abridged.name}, which its tile abridges`,
        )
      }
      const blocks = parseBlocks(cols)
      const listBp = blocks.before + blocks.bp + blocks.after
      if (listBp !== abridged.hapBp) {
        throw new WalkFormatError(
          `${refName}'s blocks row of ${abridged.name}, piece ${abridged.piece} in tile ${tile.name}, adds up to ${listBp} bp, and its abridged row states a hapBp of ${abridged.hapBp}`,
        )
      }
      tally.blocks(unit, abridged, blocks)
    }
    const windowTiles = new Set(tiles.map(u => u.name))
    const tileUnits = new Map(
      units
        .filter(u => u.kind === 'tile' && u.length > 0)
        .map(u => [Number(u.name), u]),
    )
    const fetches = new Map<WalkUnit, { names: Set<string>; from: number }>()
    for (const d of attached) {
      const from = Math.min(d.leavePos, d.rejoinPos, region.end + context - 1)
      for (const unit of (ownTiles.get(d) ?? []).flatMap(start => {
        const held = tileUnits.get(start)
        return held ? [held] : []
      })) {
        if (!windowTiles.has(unit.name)) {
          const wanting = fetches.get(unit)
          fetches.set(unit, {
            names: (wanting?.names ?? new Set()).add(d.name),
            from: Math.min(wanting?.from ?? from, from),
          })
        }
      }
    }
    const fetching = [...fetches].map(([unit, { from }]) => ({ unit, from }))
    checkBytes(fetching, true)
    signal?.throwIfAborted()
    const fetched = await this.file.read(
      fetching.map(x => x.unit),
      signal,
    )
    const owned = [...fetches.values()].flatMap(({ names }, k) => {
      const unit = tally.tile(
        split(fetched[k]!),
        Number(fetching[k]!.unit.name),
        false,
      )
      for (const walk of unit.walks) {
        if (
          names.has(walk.name) &&
          holdsOwn(
            walk.name,
            walk.hapOffset,
            walk.hapOffset + tally.walk(unit, walk, false),
          )
        ) {
          tally.walk(unit, walk, true)
        }
      }
      return unit.abridged
        .filter(
          abridged =>
            names.has(abridged.name) &&
            holdsOwn(
              abridged.name,
              abridged.hapOffset,
              abridged.hapOffset + abridged.hapBp,
            ),
        )
        .map(moved => ({ unit, moved }))
    })
    const ownDetails = detailsOf(
      owned.filter(item => expand?.(item.moved.name) === true),
    )
    signal?.throwIfAborted()
    const loaded = countDetails(
      ownDetails,
      await this.file.read(
        ownDetails.map(d => d.unit),
        signal,
      ),
    )
    for (const { moved } of owned) {
      if (!loaded.has(moved)) {
        tally.unread(moved)
      }
    }
    tally.departures(records)
    return tally.result()
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
