import { LocalFile, RemoteFile } from 'generic-filehandle2'

import {
  MissingSequenceError,
  refNameLookup,
  resolveRefName,
} from './refNames.ts'
import { checkRange } from './tabixSource.ts'
import { panSNHaplotypeOrSample, panSNMatchesPrefix } from './util.ts'
import {
  blockLinks,
  blocksOver,
  blocksSteps,
  bpRow,
  canonicalLink,
  chunkRef,
  decodeChunk,
  ownRunsRow,
  ownStretch,
  tileLinks,
  parseBlocks,
  referenceNodes,
  statedAbridgedSteps,
  tileNodeIds,
  tileSteps,
  unloadedRow,
} from './walkChunk.ts'
import {
  WalkContainer,
  WalkFormatError,
  once,
  raceAbort,
} from './walkContainer.ts'
import { WalkTally } from './walkCounts.ts'
import {
  ATTACH_BP,
  attachesOwnSteps,
  departureOwnTiles,
  parseDepartureRow,
  squeezer,
} from './walkDeparture.ts'
import { FORMS_STEP_BUDGET, walkDepartures } from './walkDepartures.ts'
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
import {
  blocksPiece,
  refuseRevisited,
  runsTile,
  skippedPiece,
  windowRuns,
} from './walkRuns.ts'

import type {
  AbridgedPiece,
  ChunkRef,
  DecodedChunk,
  PieceBlocks,
  PlacedNode,
  PieceEnds,
  RefPath,
} from './walkChunk.ts'
import type { WalkFilehandle, WalkUnit } from './walkContainer.ts'
import type { CountUnit, CountWalk } from './walkCounts.ts'
import type { WalkDeparture } from './walkDeparture.ts'
import type { BudgetNames, WalkHeader, WalkRow } from './walkRows.ts'
import type { RunsPiece } from './walkRuns.ts'

// A window reads whole tiles and keeps the steps on it, so a narrow window
// reads several times what it keeps: 464 haplotypes across the 1.5 kb ABCA7
// repeat read 5.0 M steps and keep 1.16 M, and the cut takes 0.4 s from warm
// files. The step budget bounds what the cut keeps, and rows are read to this
// multiple of it.
const READ_STEPS_PER_KEPT = 8

// bp of reference a cut keeps past each end of the window
export const WALK_CONTEXT = 1000

// reads of detail records an expanding cut makes at once
const DETAIL_READS = 16

// the most compressed bytes of detail records one of those reads fetches,
// but for a record larger than that alone
const DETAIL_READ_BYTES = 2_000_000

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

// `items` as runs to read at once: those whose units' bytes touch in the
// file, in file order, up to `maxBytes` a run, the runs in the order of
// their first item among `items`
export function touchingRuns<T extends { unit: WalkUnit }>(
  items: T[],
  maxBytes: number,
) {
  const runs: { items: T[]; bytes: number; first: number }[] = []
  const place = new Map(items.map((item, k) => [item, k]))
  for (const item of [...items].sort((a, b) => a.unit.offset - b.unit.offset)) {
    const run = runs.at(-1)
    const last = run?.items.at(-1)
    const at = place.get(item) ?? 0
    if (
      run &&
      last &&
      last.unit.offset + last.unit.length === item.unit.offset &&
      run.bytes + item.unit.length <= maxBytes
    ) {
      run.items.push(item)
      run.bytes += item.unit.length
      run.first = Math.min(run.first, at)
    } else {
      runs.push({ items: [item], bytes: item.unit.length, first: at })
    }
  }
  return runs.sort((a, b) => a.first - b.first).map(run => run.items)
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
  // the abridged walks to read from their detail records: those a function
  // accepts by `sample#haplotype#contig`, or those under any of a list's
  // PanSN prefixes (`HG00097`, `HG00097#1`); none by default. Of such a walk
  // the cut reads the records of the pieces whose steps it keeps or attaches
  expand?: string[] | ((name: string) => boolean)
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

// the sample whose contig a range's `refName` names, the file's first
// reference by default
interface OnReference {
  referenceSample?: string
}

export interface CutForRangeOptions extends WalkCutOptions, OnReference {}

export interface GraphForRangeOptions
  extends Omit<WalkCutOptions, 'keep'>, OnReference {
  // the `sample#haplotype` walks to keep beside the reference; every walk when
  // undefined or empty
  haplotypes?: string[]
}

export interface CountsForRangeOptions extends WalkCountOptions, OnReference {}

export type WalkRunsOptions = Pick<
  WalkCutOptions,
  | 'signal'
  | 'keep'
  | 'context'
  | 'contained'
  | 'byteBudget'
  | 'departureByteBudget'
  | 'budgetNames'
>

export interface RunsForRangeOptions extends WalkRunsOptions, OnReference {
  // the `sample#haplotype` walks to keep beside the reference, where `keep`
  // is not given; every walk when undefined or empty
  haplotypes?: string[]
}

export interface BytesForRangeOptions extends OnReference {
  signal?: AbortSignal
  context?: number
  attachBp?: number
  departureByteBudget?: number
}

export interface FormsForRangeOptions extends OnReference {
  // the most steps the cut keeps, FORMS_STEP_BUDGET at least
  stepBudget?: number
  // the most compressed bytes the cut fetches from the range's tiles and
  // their blocks units; WALK_BUDGETS' by default
  byteBudget?: number
  // the most compressed bytes the cut's departure read fetches;
  // WALK_BUDGETS' by default
  departureByteBudget?: number
  signal?: AbortSignal
}

export type WalkFileHeader = WalkHeader & { tile: number }

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

type OneOf<T> = {
  [K in keyof T]: Pick<T, K> & Partial<Record<Exclude<keyof T, K>, never>>
}[keyof T]

type WalksSource = OneOf<{
  walksUrl: string
  walksPath: string
  walksFilehandle: WalkFilehandle
}>

type IdxSource = OneOf<{
  idxUrl: string
  idxPath: string
  idxFilehandle: WalkFilehandle
}>

type NoIdx = Partial<Record<keyof IdxSource, never>>

/**
 * Where the walk file and its directory are: exactly one of `walksUrl`,
 * `walksPath` or `walksFilehandle`, and one of `idxUrl`, `idxPath` or
 * `idxFilehandle`, which may be left out beside a URL or a path: the
 * directory is then `<walks>.idx`, before any query string or fragment of
 * a URL
 */
export type WalkFileOptions =
  | (Exclude<WalksSource, { walksFilehandle: WalkFilehandle }> &
      (IdxSource | NoIdx))
  | (Extract<WalksSource, { walksFilehandle: WalkFilehandle }> & IdxSource)

const directoryOf = (url: string) => url.replace(/(?=$|[?#])/, '.idx')

const sourcesGiven = (...values: unknown[]) =>
  values.filter(value => value !== undefined).length

function openWalkFile(opts: WalkFileOptions) {
  const idxGiven = sourcesGiven(opts.idxUrl, opts.idxPath, opts.idxFilehandle)
  if (
    sourcesGiven(opts.walksUrl, opts.walksPath, opts.walksFilehandle) !== 1 ||
    idxGiven > 1 ||
    (idxGiven === 0 && sourcesGiven(opts.walksUrl, opts.walksPath) === 0)
  ) {
    throw new TypeError(
      'A WalkFile takes one of walksUrl, walksPath or walksFilehandle, and at most one of idxUrl, idxPath or idxFilehandle, which a walksFilehandle needs',
    )
  }
  const walks =
    opts.walksFilehandle ??
    (opts.walksPath !== undefined
      ? new LocalFile(opts.walksPath)
      : new RemoteFile(opts.walksUrl))
  const idx =
    opts.idxFilehandle ??
    (opts.idxPath !== undefined
      ? new LocalFile(opts.idxPath)
      : opts.idxUrl !== undefined
        ? new RemoteFile(opts.idxUrl)
        : opts.walksPath !== undefined
          ? new LocalFile(`${opts.walksPath}.idx`)
          : new RemoteFile(directoryOf(opts.walksUrl)))
  return new WalkContainer(walks, idx)
}

// The abridged walks a cut reads from their detail records: those a
// function accepts, or those under any of a list's PanSN prefixes
const expandFilter = (expand: WalkCutOptions['expand']) =>
  Array.isArray(expand)
    ? (name: string) => expand.some(prefix => panSNMatchesPrefix(name, prefix))
    : expand

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
 *
 * Every range method takes 0-based half-open coordinates and a `refName`
 * that is a contig of the reference sample (`chr6`, on `referenceSample`,
 * the file's first reference by default) or the file's own name for it
 * (`GRCh38#0#chr6`), and throws MissingSequenceError for a sequence the
 * file lacks.
 */
export class WalkFile {
  private readonly container: WalkContainer

  constructor(opts: WalkFileOptions) {
    this.container = openWalkFile(opts)
  }

  private readonly readHeader = once(async () => {
    const header = walkHeader(await this.container.header())
    const { tile } = header
    if (tile === undefined) {
      throw new WalkFormatError("The walk file's #walks line states no tile:i:")
    }
    return { ...header, tile }
  })

  private readonly readLookup = once(async () =>
    refNameLookup(await this.container.sequences()),
  )

  // The header, read once for every caller: one caller giving up stops only
  // its own wait, and a read that fails is tried again by the next caller
  getHeader(opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.readHeader(), opts.signal)
  }

  // The reference sequences the walk file holds windows on
  getReferenceSequences(opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.container.sequences(), opts.signal)
  }

  /**
   * The walk file's sequence for `contig` on `sample`, or for `contig` as
   * the file spells it: `chr6` on `GRCh38` is `GRCh38#0#chr6`. Undefined when
   * the file holds neither.
   */
  async getSequenceName(
    sample: string,
    contig: string,
    opts: { signal?: AbortSignal } = {},
  ) {
    const lookup = await raceAbort(this.readLookup(), opts.signal)
    return resolveRefName(lookup, sample, contig)
  }

  // The units the file holds for a reference sequence, by the file's name
  // for it, in file order; none for a sequence the file lacks
  getUnits(refName: string, opts: { signal?: AbortSignal } = {}) {
    return raceAbort(this.container.units(refName), opts.signal)
  }

  // Where the reference's last node on `refName`, the file's name for it,
  // ends, from the extents of its tiles; 0 for a sequence the file lacks
  async getSequenceLength(
    refName: string,
    opts: { signal?: AbortSignal } = {},
  ) {
    const units = await this.getUnits(refName, opts)
    return units.reduce(
      (end, u) => (u.kind === 'tile' ? Math.max(end, u.end) : end),
      0,
    )
  }

  // The file's name for `refName` and the reference sample it was read on
  private async resolve(
    refName: string,
    opts: { referenceSample?: string; signal?: AbortSignal },
  ) {
    const { signal } = opts
    const header = await this.getHeader({ signal })
    const sample = opts.referenceSample ?? header.references[0] ?? ''
    const name = await this.getSequenceName(sample, refName, { signal })
    if (name === undefined) {
      throw new MissingSequenceError(
        refName,
        sample,
        await this.getReferenceSequences({ signal }),
      )
    }
    return { refName: name, sample }
  }

  /**
   * The compressed bytes a cut of the range reads when it keeps every walk
   * and expands none: its tiles, their blocks units, its departure blocks,
   * which this reads as the cut does, and the tiles named for their own
   * steps by its records within `attachBp` that leave or rejoin on the
   * window or its context. A cut that keeps fewer walks reads no more, and
   * one that expands walks no more but for their detail units
   */
  async getBytesForRange(
    refName: string,
    start: number,
    end: number,
    opts: BytesForRangeOptions = {},
  ) {
    const sequence = await this.resolve(refName, opts)
    return this.bytesFor(sequence.refName, { start, end }, opts)
  }

  private async bytesFor(
    refName: string,
    region: { start: number; end: number },
    opts: BytesForRangeOptions,
  ) {
    const {
      signal,
      context = WALK_CONTEXT,
      attachBp = ATTACH_BP,
      departureByteBudget = WALK_BUDGETS.departureByteBudget,
    } = opts
    checkRange(refName, region.start, region.end)
    const drawn = drawnSpan(region, context)
    const units = await this.getUnits(refName, { signal })
    const { tiles, departureBlocks } = windowUnits(units, drawn)
    checkDepartureBytes(departureBlocks, departureByteBudget)
    const tileNames = new Set(tiles.map(u => u.name))
    const { records, ownTiles } = departureRecords(
      (await this.container.read(departureBlocks, signal)).flat(),
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

  // Every pass of the reference over the range, by start, from the ref rows
  // of the tiles there: a node the reference passes twice there comes twice
  async getReferenceNodes(
    refName: string,
    start: number,
    end: number,
    opts: { referenceSample?: string; signal?: AbortSignal } = {},
  ) {
    const sequence = await this.resolve(refName, opts)
    return this.referenceNodes(sequence.refName, { start, end }, opts)
  }

  private async referenceNodes(
    refName: string,
    region: { start: number; end: number },
    opts: { signal?: AbortSignal },
  ) {
    checkRange(refName, region.start, region.end)
    const tiles = (await this.getUnits(refName, opts)).filter(
      u => u.kind === 'tile' && overlaps(u, region),
    )
    const read = await this.container.read(tiles, opts.signal)
    return referenceNodes(read.flat().map(line => line.split('\t')))
      .filter(
        node =>
          node.refName === refName &&
          node.end > region.start &&
          node.start < region.end,
      )
      .sort((a, b) => a.start - b.start || a.end - b.end)
  }

  /**
   * The cut of the range: its graph, the nodes it keeps, each walk's
   * fragments and departures over it and the gaps between them, and
   * `tables()`, which makes the graph's tables of all of them
   */
  async getCutForRange(
    refName: string,
    start: number,
    end: number,
    opts: CutForRangeOptions = {},
  ) {
    const sequence = await this.resolve(refName, opts)
    return this.cut(sequence.refName, { start, end }, opts)
  }

  private async cut(
    refName: string,
    region: { start: number; end: number },
    opts: WalkCutOptions,
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
    } = opts
    const expand = expandFilter(opts.expand)
    checkRange(refName, region.start, region.end)
    const [header, units] = await Promise.all([
      this.getHeader({ signal }),
      this.getUnits(refName, { signal }),
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
    const surelyPlaceholders = Number.isFinite(placeholderLimit)
      ? new Set<number>()
      : undefined
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
          (surelyPlaceholders?.size ?? 0) > placeholderLimit
        ) {
          surely.clear()
          surelyPlaceholders?.clear()
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
            surelyPlaceholders?.size ?? 0,
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
    const readDetails = (run: Detail[]) =>
      this.container.read(
        run.map(detail => detail.unit),
        reads.signal,
      )
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
    const readBlocks = async (wanted: Needed[]) => {
      const unread = wanted.filter(w => !blocksRead.has(w.unit))
      const fetched = await this.container.read(
        unread.map(w => w.unit),
        reads.signal,
      )
      unread.forEach((w, k) => {
        extra.push(w)
        blocksRead.set(w.unit, fetched[k]!)
      })
    }
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
      await readBlocks(unread)
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
      const wanted = [...details.values()]
      const held = new Map<Detail, string[]>()
      let decoded = 0
      await readInOrder(
        touchingRuns(wanted, DETAIL_READ_BYTES),
        DETAIL_READS,
        reads.signal,
        readDetails,
        (run, fetched) => {
          run.forEach((detail, k) => held.set(detail, fetched[k] ?? []))
          while (decoded < wanted.length && held.has(wanted[decoded]!)) {
            const detail = wanted[decoded++]!
            decodeDetail(detail, split(held.get(detail) ?? []))
            held.delete(detail)
          }
        },
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
      all: ReturnType<typeof tileSteps>[],
      records: WalkDeparture[],
      ownTiles: Map<WalkDeparture, number[]>,
    ) => {
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
    // steps the cut surely keeps pass the step budget and no earlier refusal
    // can come, so the finished cut would throw it. The count is the rows'
    // (tileSteps in walkChunk.ts) of the tiles strictly inside the drawn
    // stretch; where that fits and the steps the tiles' abridged rows state
    // could take it over, the cut reads the blocks units it would read
    // anyway and adds what it builds of each abridged piece on the drawn
    // stretch (blocksSteps), in any tile of the window
    const earlyRefusal = async (
      read: string[][],
      records: WalkDeparture[],
      ownTiles: Map<WalkDeparture, number[]>,
    ) => {
      const counted = tiles.map((unit, k) =>
        unit.start > drawn.start && unit.end < drawn.end
          ? tileSteps(read[k]!, keep, expand)
          : 'edge',
      )
      const rowSteps = counted.map(tally =>
        tally === 'edge' ? undefined : tally?.kept,
      )
      const kept = rowSteps.reduce<number>((sum, n) => sum + (n ?? 0), 0)
      const all = () =>
        tiles.map((_, k) => {
          const tally = counted[k]
          return tally === 'edge' ? tileSteps(read[k]!, keep, expand) : tally
        })
      // an edge tile's `steps` are those of its abridged pieces on the drawn
      // stretch, and `whole` what its rows keep where the cut draws all of it
      const refusal = (
        steps: (number | undefined)[],
        whole: (number | undefined)[] = [],
      ) =>
        surelyKeptError(
          tiles.map((unit, k) => {
            const start = Math.max(unit.start, drawn.start)
            const end = Math.min(unit.end, drawn.end)
            const drawnShare =
              unit.end > unit.start
                ? (end - start) / (unit.end - unit.start)
                : 1
            const count = steps[k]
            return count === undefined
              ? unit
              : counted[k] === 'edge'
                ? {
                    start,
                    end,
                    steps: count,
                    likely: count + (whole[k] ?? 0) * drawnShare,
                  }
                : { ...unit, steps: count }
          }),
          stepBudget,
          region,
          context,
          filtered,
        )
      const fromBlocks = async () => {
        const tallies = all()
        const abridging = tiles.flatMap((windowTile, k) => {
          const unit = blocksUnits.get(windowTile.name)
          return tallies[k]?.blocks === true && unit
            ? [{ unit, from: windowTile.start, k }]
            : []
        })
        const builtOf = async () => {
          await readBlocks(abridging)
          return abridging.map(({ unit, k }) =>
            blocksSteps(
              read[k]!,
              blocksRead.get(unit) ?? [],
              bounds,
              keep,
              expand,
            ),
          )
        }
        const built = noEarlierRefusal(read, tallies, records, ownTiles)
          ? await builtOf()
          : undefined
        const steps = [...rowSteps]
        abridging.forEach(({ k }, i) => {
          steps[k] = (steps[k] ?? 0) + (built?.[i] ?? 0)
        })
        return built?.every(n => n !== undefined) === true
          ? refusal(
              steps,
              tallies.map(tally => tally?.kept),
            )
          : undefined
      }
      return kept > stepBudget
        ? noEarlierRefusal(read, all(), records, ownTiles)
          ? refusal(rowSteps)
          : undefined
        : kept +
              tiles.reduce(
                (sum, _, k) =>
                  sum + statedAbridgedSteps(read[k]!, keep, expand),
                0,
              ) >
            stepBudget
          ? fromBlocks()
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
      const read = await this.container.read(
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
        ? await earlyRefusal(read, records, ownTiles)
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
        const fetched = await this.container.read(
          [...fetches.keys()],
          reads.signal,
        )
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
      refName,
      graph,
      kept,
      fragments,
      departures,
      get gaps() {
        return cut.gapsOf()
      },
      placeholders,
      tables: () => graph.tables(kept, fragments, departures, cut.gapsOf()),
    }
  }

  /**
   * How many of the file's haplotypes step each node and link over the range,
   * and how often, read as a cut that keeps every walk reads, with the same
   * `context`, `contained`, `attachBp`, budgets and `expand`. It counts each
   * piece the cut builds whole, and each piece holding the own steps of a
   * record the cut attaches, while it decodes the rows, keeping no step and
   * building no graph. Every count states why it may fall short, or that it
   * is exact (docs/counts.md)
   */
  async getCountsForRange(
    refName: string,
    start: number,
    end: number,
    opts: CountsForRangeOptions = {},
  ) {
    const sequence = await this.resolve(refName, opts)
    return this.counts(sequence.refName, { start, end }, opts)
  }

  private async counts(
    refName: string,
    region: { start: number; end: number },
    opts: WalkCountOptions,
  ) {
    const {
      signal,
      context = WALK_CONTEXT,
      contained = false,
      stepBudget = WALK_BUDGETS.stepBudget,
      byteBudget = WALK_BUDGETS.byteBudget,
      departureByteBudget = WALK_BUDGETS.departureByteBudget,
      attachBp = ATTACH_BP,
    } = opts
    const expand = expandFilter(opts.expand)
    checkRange(refName, region.start, region.end)
    const [header, units, sequences] = await Promise.all([
      this.getHeader({ signal }),
      this.getUnits(refName, { signal }),
      this.getReferenceSequences({ signal }),
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
    const read = await this.container.read(
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
    const more = await this.container.read(
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
    const fetched = await this.container.read(
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
      await this.container.read(
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
   * The graph under the range as tables: the reference sample's sequence
   * over it, with `context` bp either side, and the walks of the haplotypes
   * asked for, every one by default; `refName` names the sequence as the
   * file does, and `sequences` holds the nodes' bases when asked for.
   */
  async getGraphForRange(
    refName: string,
    start: number,
    end: number,
    opts: GraphForRangeOptions = {},
  ) {
    const { haplotypes, ...rest } = opts
    const sequence = await this.resolve(refName, opts)
    const cut = await this.cut(
      sequence.refName,
      { start, end },
      { ...rest, keep: walksKept(haplotypes, sequence.refName) },
    )
    return {
      tables: cut.tables(),
      refName: sequence.refName,
      sample: sequence.sample,
      sequences: opts.sequences === true ? cut.graph.sequences : undefined,
    }
  }

  /**
   * The walks getCutForRange returns for the range with the same options, as
   * runs (walkRuns.ts), read from the rows without decoding a step: each
   * walk row's route choices and each abridged piece's blocks row, from the
   * tiles and departure blocks over the window and its context and those
   * tiles' blocks units, with no detail unit and no tile past them. Walks
   * are those `keep` accepts, or of the `haplotypes` asked for beside the
   * reference's, every one by default. It refuses a file whose header does
   * not state revisited:i:0 with a RunsRefusedError, and a window past the
   * byte budgets as a cut does.
   */
  async getRunsForRange(
    refName: string,
    start: number,
    end: number,
    opts: RunsForRangeOptions = {},
  ) {
    const { haplotypes, ...rest } = opts
    const sequence = await this.resolve(refName, opts)
    return this.runs(
      sequence.refName,
      { start, end },
      {
        ...rest,
        keep: rest.keep ?? walksKept(haplotypes, sequence.refName),
      },
    )
  }

  private async runs(
    refName: string,
    region: { start: number; end: number },
    opts: WalkRunsOptions = {},
  ) {
    const {
      signal,
      keep,
      context = WALK_CONTEXT,
      contained = false,
      byteBudget = WALK_BUDGETS.byteBudget,
      departureByteBudget = WALK_BUDGETS.departureByteBudget,
    } = opts
    checkRange(refName, region.start, region.end)
    const [header, units] = await Promise.all([
      this.getHeader({ signal }),
      this.getUnits(refName, { signal }),
    ])
    const refused = refuseRevisited(header.revisited)
    if (refused) {
      throw refused
    }
    const drawn = drawnSpan(region, context)
    const { tiles, departureBlocks } = windowUnits(units, drawn)
    const tileBytesTo = (end: number) =>
      bytesOf(tiles.filter(u => u.start < end + context))
    const checkBytes = (bytesTo: (end: number) => number) => {
      const tooLarge = Number.isFinite(byteBudget)
        ? byteBudgetError(
            bytesTo,
            byteBudget,
            region,
            header.tile,
            opts.budgetNames,
          )
        : undefined
      if (tooLarge) {
        throw tooLarge
      }
    }
    checkBytes(tileBytesTo)
    checkDepartureBytes(departureBlocks, departureByteBudget)
    const read = await this.container.read(
      [...tiles, ...departureBlocks],
      signal,
    )
    const { records } = departureRecords(
      read.slice(tiles.length).flat(),
      drawn,
      keep,
    )
    const decoded = tiles.map((unit, k) =>
      runsTile(split(read[k]!), Number(unit.name), keep),
    )
    const blocksUnits = new Map(
      units.filter(u => u.kind === 'blocks').map(u => [u.name, u]),
    )
    const wanted = tiles.flatMap((tileUnit, k) => {
      const unit = blocksUnits.get(tileUnit.name)
      if (decoded[k]!.abridged.length > 0 && !unit) {
        throw new WalkFormatError(
          `${refName}'s unit table lists no blocks unit ${tileUnit.name}, which the abridged rows of that tile need`,
        )
      }
      return unit && decoded[k]!.abridged.length > 0
        ? [{ unit, from: tileUnit.start, k }]
        : []
    })
    checkBytes(
      end =>
        tileBytesTo(end) +
        bytesOf(wanted.filter(x => x.from < end + context).map(x => x.unit)),
    )
    const blocksRead = await this.container.read(
      wanted.map(x => x.unit),
      signal,
    )
    const blocksRows = new Map<number, Map<string, string[]>>()
    wanted.forEach(({ k }, w) => {
      blocksRows.set(
        k,
        new Map(
          split(blocksRead[w]!).map(cols => [
            pieceKey(cols[4] ?? '', cols[5] ?? '', cols[6] ?? ''),
            cols,
          ]),
        ),
      )
    })
    // the drawn stretch as far as the reference's passes over it reach, and
    // the window and its context as walkCut widens them
    const over = { ...drawn }
    const bounds = {
      start: region.start - context,
      end: region.end + context,
    }
    const window = { ...bounds }
    for (const { path } of decoded) {
      const hull = path.hull(drawn)
      over.start = Math.min(over.start, hull?.start ?? over.start)
      over.end = Math.max(over.end, hull?.end ?? over.end)
      const reach = path.hull(bounds)
      window.start = Math.min(window.start, reach?.start ?? window.start)
      window.end = Math.max(window.end, reach?.end ?? window.end)
    }
    const onOver = (span: { start: number; end: number }) =>
      span.start < over.end && span.end > over.start
    const abridged = decoded.flatMap((tile, k) =>
      tile.abridged.map(piece => ({ piece, k })),
    )
    const lines = new Map<string, LinePiece<AbridgedPiece>[]>()
    for (const { piece } of abridged) {
      const key = `${piece.name}\t${piece.fragStart}`
      const line = lines.get(key) ?? []
      lines.set(key, line)
      line.push({ piece: piece.piece, on: onOver(piece.span), item: piece })
    }
    for (const tile of decoded) {
      for (const piece of tile.pieces) {
        lines.get(`${piece.name}\t${piece.fragStart}`)?.push({
          piece: piece.piece,
          on: piece.span === undefined || onOver(piece.span),
        })
      }
    }
    const inOrder = [...lines.values()].map(line =>
      line.sort((a, b) => a.piece - b.piece),
    )
    const drawnPieces = new Set(betweenDrawn(inOrder).map(held => held.item))
    const reached = new Set(
      contained && context > 0 ? beforeDrawn(inOrder) : [],
    )
    const ownRecords = new Map<string, WalkDeparture[]>()
    for (const d of records) {
      if (d.ownBp > 0) {
        ownRecords.set(d.name, [...(ownRecords.get(d.name) ?? []), d])
      }
    }
    const built = (piece: AbridgedPiece) => {
      const runs = (ownRecords.get(piece.name) ?? []).flatMap(d => {
        const stretch = ownStretch(piece, d)
        return stretch
          ? [{ ...stretch, attached: attachesOwnSteps(d, drawn, ATTACH_BP) }]
          : []
      })
      return (
        drawnPieces.has(piece) ||
        runs.some(run => run.attached) ||
        (reached.has(piece) && runs.some(run => run.trailing))
      )
    }
    const pieces: RunsPiece[] = decoded.flatMap(tile => tile.pieces)
    for (const { piece, k } of abridged) {
      const cols = blocksRows.get(k)?.get(keyOf(piece))
      if (built(piece) && !cols) {
        throw new WalkFormatError(
          `${refName}'s blocks unit ${tiles[k]!.name} holds no row for piece ${piece.piece} of ${piece.name}, which its tile abridges`,
        )
      }
      pieces.push(
        cols && built(piece)
          ? blocksPiece(piece, cols, decoded[k]!.path)
          : skippedPiece(piece),
      )
    }
    return windowRuns({
      refName,
      pieces,
      records,
      window,
      contained,
      context,
      holes: header.holes,
      revisited: header.revisited ?? Number.NaN,
    })
  }

  /**
   * Every haplotype's departures over the range, for structuralForms, from
   * one cut without links that expands no walk. A haplotype the header names
   * with no walk in the cut is listed with nothing covered.
   */
  async getFormsForRange(
    refName: string,
    start: number,
    end: number,
    opts: FormsForRangeOptions = {},
  ) {
    const {
      stepBudget = FORMS_STEP_BUDGET,
      byteBudget,
      departureByteBudget,
      signal,
    } = opts
    const [header, sequence] = await Promise.all([
      this.getHeader({ signal }),
      this.resolve(refName, opts),
    ])
    const cut = await this.cut(
      sequence.refName,
      { start, end },
      {
        links: 'none',
        stepBudget: Math.max(stepBudget, FORMS_STEP_BUDGET),
        byteBudget,
        departureByteBudget,
        signal,
      },
    )
    const haplotypes = walkDepartures(
      {
        fragments: cut.fragments,
        departures: cut.departures,
        nodes: cut.graph.nodes,
      },
      sequence.refName,
    )
    const seen = new Set(haplotypes.map(h => h.haplotype))
    for (const haplotype of header.haplotypes) {
      if (!seen.has(haplotype)) {
        haplotypes.push({ haplotype, covered: [], departures: [] })
      }
    }
    return haplotypes
  }
}

export type CutForRange = Awaited<ReturnType<WalkFile['getCutForRange']>>

export type GraphForRange = Awaited<ReturnType<WalkFile['getGraphForRange']>>

export type RunsForRange = Awaited<ReturnType<WalkFile['getRunsForRange']>>
