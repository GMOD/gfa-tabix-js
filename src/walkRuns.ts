import {
  RefPath,
  abridgedPiece,
  blockPositions,
  checkedId,
  idOf,
  interiorAt,
  numbers,
  parseBlocks,
  parseChoices,
  parseCopy,
  refRows,
  tokens,
  traverse,
  NONE,
  WalkRowError,
} from './walkChunk.ts'
import { adjacentSteps, departureStretch, squeezer } from './walkDeparture.ts'
import { chainCrossings } from './walkDepartures.ts'
import { consecutivePieces, edgeKey } from './walkRows.ts'

import type { GraphTables } from './util.ts'
import type { AbridgedPiece, Route } from './walkChunk.ts'
import type { Hole, WalkDeparture } from './walkDeparture.ts'
import type { ChainEvent } from './walkDepartures.ts'

/**
 * What a walk does past one end of a fragment, read from the departure that
 * leaves its last step on the reference, or rejoins its first: `beyond`, it
 * is next on the reference past the drawn stretch on that side, the same way
 * round; `elsewhere`, it is next on the reference reversed, or back on the
 * other side; `edge`, no departure, and the end lies on the drawn stretch's
 * edge, so whatever the walk does next lies off it; `ends`, no departure and
 * the end lies inside the drawn stretch: the walk has no step on the
 * reference sequence past it.
 */
export const END_KINDS = ['beyond', 'elsewhere', 'edge', 'ends'] as const

export type EndKind = (typeof END_KINDS)[number]

/**
 * One end of each fragment, by END_KINDS index in `kinds`. For `beyond` and
 * `elsewhere`: `pos`, the departure's breakpoint off the fragment, its leave
 * for a `before` end and its rejoin for an `after` one, `reversed` 1 where
 * that step runs against the reference, `ownBp` the walk's bp between the
 * fragment's stretch and that step, a lead or tail included, and
 * `missingBp` the part between the contig's W lines. For `ends`, `ownBp` is
 * the walk's bp the read holds past the stretch, all of it off the
 * reference; for `edge` it is NaN. `pos` is NaN for both.
 */
export interface WalkEnds {
  kinds: Uint8Array
  pos: Float64Array
  reversed: Uint8Array
  ownBp: Float64Array
  missingBp: Float64Array
}

/**
 * A window's walks as runs, a fragment per walk-table walk of the cut, in
 * its order: the reference's fragments first, then by name and contig
 * position. A stretch is consecutive steps on the reference, each starting on
 * the base the last ended on, all one way; it never spans a hole. A run is
 * what lies between two stretches of a fragment. A fragment's stretches lie
 * at `offsets[p]` up to `offsets[p + 1]`, each over `[starts, starts + bp)`
 * of the reference, `against` 1 where the walk takes it against the
 * reference, `offBefore` the bp of the walk's own run before it, 0 for a
 * fragment's first, and `loaded` 0 where that run holds a placeholder, a
 * stretch of a piece the read did not load. `leadBp` and `tailBp` are the
 * bp off the reference before a fragment's first stretch and after its last,
 * a contained cut's context. The reference's nodes are those of its
 * fragments, and the drawn stretch spans them.
 */
export interface WalkRuns {
  refName: string
  drawn: { refName: string; start: number; end: number }
  revisited: number
  holes: { start: number; end: number }[]
  reference: { names: string[]; starts: Float64Array; lengths: Int32Array }
  walks: {
    names: string[]
    starts: Float64Array
    ends: Float64Array
    trimmedBefore: Float64Array
    trimmedAfter: Float64Array
    offsets: Int32Array
    leadBp: Float64Array
    leadLoaded: Uint8Array
    tailBp: Float64Array
    tailLoaded: Uint8Array
    before: WalkEnds
    after: WalkEnds
  }
  stretches: {
    starts: Float64Array
    bp: Float64Array
    against: Uint8Array
    offBefore: Float64Array
    loaded: Uint8Array
  }
  // the departures over the drawn stretch of the walks with no fragment,
  // such as one that jumps the window, as the cut's departure table lists
  // them, by walk, then contig position and leave
  departing: {
    names: string[]
    hapStarts: Float64Array
    ownBp: Float64Array
    missingBp: Float64Array
    leavePos: Float64Array
    leaveReversed: Uint8Array
    rejoinPos: Float64Array
    rejoinReversed: Uint8Array
  }
}

// A WalkRuns' departing departures, of walks with no fragment
export type WalkDeparting = WalkRuns['departing']

/**
 * A file whose header does not state revisited:i:0: a run off the reference
 * may hold a node the reference passes elsewhere, so the runs would not say
 * what the steps do. A host cuts the window instead.
 */
export class RunsRefusedError extends Error {
  override name = 'RunsRefusedError'
}

export const refuseRevisited = (revisited: number | undefined) =>
  revisited === 0
    ? undefined
    : new RunsRefusedError(
        `The walk file's header states ${revisited === undefined ? 'no revisited:i:' : `revisited:i:${revisited}`}, and runs are read only from a file stating revisited:i:0; cut the window instead`,
      )

// A stretch, as one piece's steps along its tile's R from position `a` to
// `b` in the walk's order
interface OnItem {
  kind: 'on'
  piece: number
  path: RefPath
  a: number
  b: number
  against: boolean
}

// Steps off R, each step's bp; one placeholder where not `loaded`
interface OffItem {
  kind: 'off'
  piece: number
  bp: number[]
  loaded: boolean
}

// A piece the read does not build, as a cut holds it: its bp alone
interface SkipItem {
  kind: 'skip'
  piece: number
  bp: number
}

type Item = OnItem | OffItem | SkipItem

export interface RunsPiece {
  name: string
  fragStart: number
  piece: number
  hapOffset: number
  items: Item[]
}

const breaksHeld = new WeakMap<RefPath, { n: number; breaks: number[] }>()

// The positions of R that start past a hole, where the one before ends
function breaksOf(path: RefPath) {
  const n = path.starts.length
  let held = breaksHeld.get(path)
  if (held?.n !== n) {
    const breaks: number[] = []
    for (let i = 1; i < n; i++) {
      if (path.starts[i] !== path.end(i - 1)) {
        breaks.push(i)
      }
    }
    held = { n, breaks }
    breaksHeld.set(path, held)
  }
  return held.breaks
}

function firstAbove(sorted: number[], value: number) {
  let lo = 0
  let hi = sorted.length
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    if (sorted[mid]! <= value) {
      lo = mid + 1
    } else {
      hi = mid
    }
  }
  return lo
}

// Whether R holds a hole between positions `a` and `b`, one apart
const holeBetween = (breaks: number[], a: number, b: number) =>
  breaks[firstAbove(breaks, Math.max(a, b) - 1)] === Math.max(a, b)

// A piece's items as its steps come, each stretch of R split at holes and
// joined to the stretch before where it carries on from it
class PieceItems {
  readonly items: Item[] = []
  private readonly path: RefPath
  private readonly piece: number

  constructor(path: RefPath, piece: number) {
    this.path = path
    this.piece = piece
  }

  // R's positions `from` to `to`, in the walk's order
  along(from: number, to: number, against: boolean) {
    const breaks = breaksOf(this.path)
    const lo = Math.min(from, to)
    const hi = Math.max(from, to)
    const cuts: number[] = []
    for (
      let k = firstAbove(breaks, lo);
      k < breaks.length && breaks[k]! <= hi;
      k++
    ) {
      cuts.push(breaks[k]!)
    }
    let at = from
    for (const cut of against ? cuts.toReversed() : cuts) {
      const end = against ? cut : cut - 1
      this.add(at, end, against)
      at = against ? end - 1 : end + 1
    }
    this.add(at, to, against)
  }

  private add(from: number, to: number, against: boolean) {
    const last = this.items.at(-1)
    if (
      last?.kind === 'on' &&
      last.path === this.path &&
      last.against === against &&
      last.b + (against ? -1 : 1) === from &&
      !holeBetween(breaksOf(this.path), last.b, from)
    ) {
      last.b = to
    } else {
      this.items.push({
        kind: 'on',
        piece: this.piece,
        path: this.path,
        a: from,
        b: to,
        against,
      })
    }
  }

  off(bp: number, loaded: boolean) {
    const last = this.items.at(-1)
    if (last?.kind === 'off' && last.loaded === loaded) {
      last.bp.push(bp)
    } else {
      this.items.push({ kind: 'off', piece: this.piece, bp: [bp], loaded })
    }
  }
}

/**
 * A tile's R, its routes and its node lengths, its kept pieces as runs
 * items, the `ref` and `walk` rows decoded and the `abridged` rows listed
 * for blocksPiece, each walk row's span
 */
export function runsTile(
  rows: string[][],
  tileStart: number,
  keep: ((name: string) => boolean) | undefined,
) {
  const path = new RefPath()
  const lengths = new Map<number, number>()
  const routes: Route[] = []
  const pieces: (RunsPiece & { span?: { start: number; end: number } })[] = []
  const abridged: AbridgedPiece[] = []
  const kept = (name: string) => keep === undefined || keep(name)
  const references = refRows(rows).map(cols => {
    const handles = tokens(cols[9]!, { id: 0 }, [])
    const first = path.handles.length
    path.add(handles, numbers(cols[10]!), +cols[6]!)
    return { cols, first, n: handles.length }
  })
  for (const { cols, first, n } of references) {
    if (kept(cols[4]!) && n > 0) {
      const items = new PieceItems(path, +cols[7]!)
      items.along(first, first + n - 1, false)
      pieces.push({
        name: cols[4]!,
        fragStart: +cols[5]!,
        piece: +cols[7]!,
        hapOffset: +cols[6]!,
        items: items.items,
      })
    }
  }
  const previous = { id: 0 }
  for (const cols of rows) {
    if (cols[3] === 'nodes') {
      const ids = numbers(cols[4]!)
      const nodeLengths = numbers(cols[5]!)
      let id = 0
      for (let k = 0; k < ids.length; k++) {
        id = checkedId(id + ids[k]!)
        lengths.set(id, nodeLengths[k]!)
      }
    } else if (cols[3] === 'route') {
      const entry = tokens(cols[4]!, previous, [])
      const interior = tokens(cols[5]!, previous, [])
      const exit = tokens(cols[6]!, previous, [])
      routes.push({ entry: entry[0] ?? NONE, interior, exit: exit[0] ?? NONE })
    }
  }
  for (const cols of rows) {
    const name = cols[4]!
    if (cols[3] === 'abridged' && kept(name)) {
      abridged.push(abridgedPiece(cols, tileStart))
    } else if (cols[3] === 'walk' && kept(name)) {
      const piece = +cols[7]!
      const items = new PieceItems(path, piece)
      const [start, startCopy] = parseCopy(cols[9]!)
      const n = +cols[8]!
      let steps = 0
      traverse(path, routes, start, startCopy, parseChoices(cols[10]!), n, {
        step(h, i) {
          items.along(i, i, path.handles[i] !== h)
          steps++
        },
        stretch(from, to, forward) {
          items.along(forward ? from + 1 : from - 1, to, !forward)
          steps += Math.abs(to - from)
        },
        route(index, flip) {
          const route = routes[index]!
          for (let k = 0; k < route.interior.length; k++) {
            items.off(lengths.get(idOf(interiorAt(route, flip, k))) ?? 0, true)
          }
          steps += route.interior.length
          return true
        },
      })
      if (steps !== n) {
        throw new WalkRowError(`a walk row decodes to ${steps} steps, not ${n}`)
      }
      pieces.push({
        name,
        fragStart: +cols[5]!,
        piece,
        hapOffset: +cols[6]!,
        items: items.items,
        span: { start: +cols[1]!, end: +cols[2]! },
      })
    }
  }
  return { path, pieces, abridged }
}

/**
 * An abridged piece as runs items from its `blocks` row: each block a
 * stretch of the tile's R and each run off R one placeholder
 */
export function blocksPiece(
  piece: AbridgedPiece,
  cols: string[],
  path: RefPath,
): RunsPiece {
  const blocks = parseBlocks(cols)
  if (blocks.before + blocks.bp + blocks.after !== piece.hapBp) {
    throw new WalkRowError(
      `the blocks row of ${piece.name}, piece ${piece.piece}, adds up to ${blocks.before + blocks.bp + blocks.after} bp, and its abridged row states a hapBp of ${piece.hapBp}`,
    )
  }
  const { from, to } = blockPositions(piece, blocks, path)
  const items = new PieceItems(path, piece.piece)
  if (blocks.before > 0) {
    items.off(blocks.before, false)
  }
  for (let k = 0; k < from.length; k++) {
    if (blocks.off[k]! > 0) {
      items.off(blocks.off[k]!, false)
    }
    const against = blocks.against[k]!
    items.along(
      against ? to[k]! - 1 : from[k]!,
      against ? from[k]! : to[k]! - 1,
      against,
    )
  }
  if (blocks.after > 0) {
    items.off(blocks.after, false)
  }
  return {
    name: piece.name,
    fragStart: piece.fragStart,
    piece: piece.piece,
    hapOffset: piece.hapOffset,
    items: items.items,
  }
}

// An abridged piece a read does not build: its bp, as a cut holds it
export const skippedPiece = (piece: AbridgedPiece): RunsPiece => ({
  name: piece.name,
  fragStart: piece.fragStart,
  piece: piece.piece,
  hapOffset: piece.hapOffset,
  items: [{ kind: 'skip', piece: piece.piece, bp: piece.hapBp }],
})

const stepsIn = (item: Item) =>
  item.kind === 'on'
    ? Math.abs(item.b - item.a) + 1
    : item.kind === 'off'
      ? item.bp.length
      : 1

// The position of R step k of a stretch takes
const positionAt = (item: OnItem, k: number) =>
  item.against ? item.a - k : item.a + k

const stepBp = (item: Item, k: number) =>
  item.kind === 'on'
    ? item.path.end(positionAt(item, k)) -
      item.path.starts[positionAt(item, k)]!
    : item.kind === 'off'
      ? item.bp[k]!
      : item.bp

// bp of an item's first k steps
const bpBefore = (item: Item, k: number) =>
  item.kind === 'on'
    ? item.against
      ? item.path.before[item.a + 1]! - item.path.before[item.a + 1 - k]!
      : item.path.before[item.a + k]! - item.path.before[item.a]!
    : item.kind === 'off'
      ? item.bp.slice(0, k).reduce((sum, bp) => sum + bp, 0)
      : k > 0
        ? item.bp
        : 0

const itemBp = (item: Item) => bpBefore(item, stepsIn(item))

// The steps of a stretch whose pass lies over [start, end), as the first and
// last of them in the walk's order; none where `first` > `last`
function stepsOver(item: OnItem, start: number, end: number) {
  const { path } = item
  const lo = Math.min(item.a, item.b)
  const hi = Math.max(item.a, item.b)
  let from = lo
  let past = hi + 1
  while (from < past) {
    const mid = (from + past) >> 1
    if (path.end(mid) > start) {
      past = mid
    } else {
      from = mid + 1
    }
  }
  let low = from
  let high = hi + 1
  while (low < high) {
    const mid = (low + high) >> 1
    if (path.starts[mid]! < end) {
      low = mid + 1
    } else {
      high = mid
    }
  }
  const p1 = from
  const p2 = low - 1
  return item.against
    ? { first: item.a - p2, last: item.a - p1 }
    : { first: p1 - item.a, last: p2 - item.a }
}

interface Step {
  i: number
  k: number
}

interface Span {
  first: Step
  last: Step
  firstItem: OnItem
  lastItem: OnItem
  leadBp: number
  leadLoaded: boolean
  tailBp: number
  tailLoaded: boolean
}

// The step after `at` in a run, or before it
function stepPast(items: Item[], at: Step, backward: boolean) {
  let { i, k } = at
  if (backward) {
    k--
    while (k < 0 && i > 0) {
      i--
      k = stepsIn(items[i]!) - 1
    }
  } else {
    k++
    while (i < items.length && k >= stepsIn(items[i]!)) {
      i++
      k = 0
    }
  }
  return { i, k }
}

// A run's off-R steps past `from` back to `to`, or on to it, kept up to
// `context` bp: whole steps, but of a placeholder only the bp left
function contextBp(
  items: Item[],
  from: Step,
  to: Step,
  context: number,
  backward: boolean,
) {
  let left = context
  let bp = 0
  let loaded = true
  const before = (a: Step, b: Step) => a.i < b.i || (a.i === b.i && a.k < b.k)
  let at = stepPast(items, from, backward)
  while (
    left > 0 &&
    at.i >= 0 &&
    at.i < items.length &&
    at.k >= 0 &&
    (backward ? !before(at, to) : before(at, to))
  ) {
    const item = items[at.i]!
    const size = stepBp(item, at.k)
    const placeholder = item.kind === 'off' && !item.loaded
    bp += placeholder && size > left ? left : size
    loaded &&= item.kind === 'off' && item.loaded
    left -= size
    at = stepPast(items, at, backward)
  }
  return { bp, loaded }
}

// The spans of a run a cut keeps, as walkSpans in walkRows.ts finds them
function runSpans(
  items: Item[],
  window: { start: number; end: number },
  contained: boolean,
  context: number,
) {
  const over = items.map(item =>
    item.kind === 'on'
      ? { item, ...stepsOver(item, window.start, window.end) }
      : undefined,
  )
  const spans: Span[] = []
  let from: Step = { i: 0, k: 0 }
  let first: (Step & { item: OnItem }) | undefined
  let last: (Step & { item: OnItem }) | undefined
  const close = (at: Step) => {
    if (first && last) {
      const lead = contained
        ? contextBp(items, first, from, context, true)
        : { bp: 0, loaded: true }
      const tail = contained
        ? contextBp(items, last, at, context, false)
        : { bp: 0, loaded: true }
      spans.push({
        first,
        last,
        firstItem: first.item,
        lastItem: last.item,
        leadBp: lead.bp,
        leadLoaded: lead.loaded,
        tailBp: tail.bp,
        tailLoaded: tail.loaded,
      })
    }
    first = undefined
    last = undefined
  }
  items.forEach((item, i) => {
    const on = over[i]
    const n = stepsIn(item)
    if (!contained) {
      if (on && on.first <= on.last) {
        first ??= { i, k: on.first, item: on.item }
        last = { i, k: on.last, item: on.item }
      }
    } else if (item.kind === 'skip') {
      close({ i, k: 0 })
      from = { i: i + 1, k: 0 }
    } else if (on) {
      if (on.first > on.last) {
        close({ i, k: 0 })
        from = { i: i + 1, k: 0 }
      } else {
        if (on.first > 0) {
          close({ i, k: 0 })
          from = { i, k: on.first }
        }
        first ??= { i, k: on.first, item: on.item }
        last = { i, k: on.last, item: on.item }
        if (on.last < n - 1) {
          close({ i, k: on.last + 1 })
          from = { i: i + 1, k: 0 }
        }
      }
    }
  })
  close({ i: items.length, k: 0 })
  return spans
}

interface Departed {
  pos: number
  reversed: boolean
  ownBp: number
  missingBp: number
}

interface Fragment {
  name: string
  start: number
  end: number
  trimmedBefore: number
  trimmedAfter: number
  leadBp: number
  leadLoaded: boolean
  tailBp: number
  tailLoaded: boolean
  stretches: Stretch[]
  // the breakpoints of the first stretch's first step and the last's last
  firstBreak: number
  lastBreak: number
  before?: Departed
  after?: Departed
  // the walk's bp the read holds before the first stretch and after the last
  readBefore: number
  readAfter: number
  reference?: { names: string[]; starts: number[]; lengths: number[] }
}

interface Stretch {
  lo: number
  hi: number
  against: boolean
  offBefore: number
  loaded: boolean
}

// A departure of a walk with no fragment, computed or a `record`
interface Departing {
  name: string
  hapStart: number
  ownBp: number
  missingBp: number
  leavePos: number
  leaveReversed: boolean
  rejoinPos: number
  rejoinReversed: boolean
  record: boolean
}

// Departing departures as the cut sorts its departure table: by walk name,
// then contig position, then leave, computed ones before records
function departingTable(departing: Departing[]): WalkRuns['departing'] {
  const sorted = departing.toSorted(
    (a, b) =>
      (a.name < b.name ? -1 : a.name > b.name ? 1 : 0) ||
      a.hapStart - b.hapStart ||
      a.leavePos - b.leavePos ||
      Number(a.record) - Number(b.record),
  )
  return {
    names: sorted.map(d => d.name),
    hapStarts: Float64Array.from(sorted, d => d.hapStart),
    ownBp: Float64Array.from(sorted, d => d.ownBp),
    missingBp: Float64Array.from(sorted, d => d.missingBp),
    leavePos: Float64Array.from(sorted, d => d.leavePos),
    leaveReversed: Uint8Array.from(sorted, d => (d.leaveReversed ? 1 : 0)),
    rejoinPos: Float64Array.from(sorted, d => d.rejoinPos),
    rejoinReversed: Uint8Array.from(sorted, d => (d.rejoinReversed ? 1 : 0)),
  }
}

// A record as a fragment's end: its rejoin past a `leave`, else its leave
const departed = (
  d: WalkDeparture | undefined,
  leave: boolean,
): Departed | undefined =>
  d
    ? {
        pos: leave ? d.rejoinPos : d.leavePos,
        reversed: leave ? d.rejoinReversed : d.leaveReversed,
        ownBp: d.ownBp,
        missingBp: d.missingBp,
      }
    : undefined

// The breakpoint of a step on R where the walk leaves it, or rejoins it
const breakpoint = (lo: number, hi: number, against: boolean, leave: boolean) =>
  against === leave ? lo : hi

export interface RunsInput {
  refName: string
  pieces: RunsPiece[]
  records: WalkDeparture[]
  window: { start: number; end: number }
  contained: boolean
  context: number
  holes: Hole[]
  revisited: number
}

/**
 * The runs of a window's pieces, as a cut over the same window keeps them:
 * each line's pieces in a row of piece numbers joined, each kept from its
 * first step on `window` to its last, split where a contained cut splits
 * it, and each fragment's ends from the departure at its first and last
 * steps on the reference, computed in their piece or read from `records`
 */
export function windowRuns({
  refName,
  pieces,
  records,
  window,
  contained,
  context,
  holes,
  revisited,
}: RunsInput): WalkRuns {
  const squeeze = squeezer(holes, refName)
  const leaves = new Map<string, WalkDeparture>()
  const rejoins = new Map<string, WalkDeparture>()
  for (const d of records) {
    leaves.set(edgeKey(d.name, d.leavePos, d.leaveReversed, d.hapStart), d)
    rejoins.set(
      edgeKey(d.name, d.rejoinPos, d.rejoinReversed, d.hapStart + d.ownBp),
      d,
    )
  }
  const fragments: Fragment[] = []
  const jumping: Departing[] = []
  for (const run of consecutivePieces(pieces)) {
    const items = run.flatMap(piece => piece.items)
    const { name } = run[0]!
    const runStart = run[0]!.hapOffset
    const haps = new Float64Array(items.length + 1)
    haps[0] = runStart
    items.forEach((item, i) => {
      haps[i + 1] = haps[i]! + itemBp(item)
    })
    const runEnd = haps[items.length]!
    const hapAt = ({ i, k }: Step) => haps[i]! + bpBefore(items[i]!, k)
    const coordinates = (item: OnItem, k: number) => {
      const p = positionAt(item, k)
      return { lo: item.path.starts[p]!, hi: item.path.end(p) }
    }
    const spans = runSpans(items, window, contained, context)
    if (spans.length === 0 && run[0]!.name !== refName) {
      for (const d of pieceDepartures(items, haps, window, squeeze)) {
        jumping.push({ ...d, name: run[0]!.name })
      }
    }
    for (const span of spans) {
      const { first, last, firstItem, lastItem } = span
      const start = hapAt(first) - span.leadBp
      const lastEnd = hapAt({ i: last.i, k: last.k + 1 })
      const end = lastEnd + span.tailBp
      const stretches: Stretch[] = []
      let off = 0
      let offLoaded = true
      let offSteps = false
      for (let i = first.i; i <= last.i; i++) {
        const item = items[i]!
        if (item.kind === 'on') {
          const from = i === first.i ? first.k : 0
          const to = i === last.i ? last.k : stepsIn(item) - 1
          const a = coordinates(item, from)
          const b = coordinates(item, to)
          const lo = Math.min(a.lo, b.lo)
          const hi = Math.max(a.hi, b.hi)
          const held = stretches.at(-1)
          if (
            held &&
            !offSteps &&
            held.against === item.against &&
            (item.against ? hi === held.lo : lo === held.hi)
          ) {
            held.lo = Math.min(held.lo, lo)
            held.hi = Math.max(held.hi, hi)
          } else {
            stretches.push({
              lo,
              hi,
              against: item.against,
              offBefore: off,
              loaded: offLoaded,
            })
          }
          off = 0
          offLoaded = true
          offSteps = false
        } else {
          off += itemBp(item)
          offLoaded &&= item.kind === 'off' && item.loaded
          offSteps = true
        }
      }
      const head = coordinates(firstItem, first.k)
      const tail = coordinates(lastItem, last.k)
      const firstBreak = breakpoint(head.lo, head.hi, firstItem.against, false)
      const lastBreak = breakpoint(tail.lo, tail.hi, lastItem.against, true)
      const computed = (
        own: number,
        item: OnItem,
        k: number,
        leave: boolean,
        at: number,
        against: boolean,
      ): Departed | undefined => {
        const c = coordinates(item, k)
        const pos = breakpoint(c.lo, c.hi, item.against, !leave)
        return own === 0 &&
          item.against === against &&
          squeeze(pos) === squeeze(at)
          ? undefined
          : { pos, reversed: item.against, ownBp: own, missingBp: 0 }
      }
      // the next step on R past item i in its piece, `by` 1 or -1, and the
      // bp off R before it
      const nearestOn = (i: number, by: number) => {
        const { piece } = items[i]!
        let j = i + by
        let own = 0
        while (items[j]?.piece === piece && items[j]?.kind !== 'on') {
          own += itemBp(items[j]!)
          j += by
        }
        const found = items[j]
        return found?.kind === 'on' && found.piece === piece
          ? { item: found, own }
          : undefined
      }
      const next = nearestOn(last.i, 1)
      const previous = nearestOn(first.i, -1)
      const after =
        last.k < stepsIn(lastItem) - 1
          ? undefined
          : next
            ? computed(
                next.own,
                next.item,
                0,
                true,
                lastBreak,
                lastItem.against,
              )
            : departed(
                leaves.get(edgeKey(name, lastBreak, lastItem.against, lastEnd)),
                true,
              )
      const before =
        first.k > 0
          ? undefined
          : previous
            ? computed(
                previous.own,
                previous.item,
                stepsIn(previous.item) - 1,
                false,
                firstBreak,
                firstItem.against,
              )
            : departed(
                rejoins.get(
                  edgeKey(name, firstBreak, firstItem.against, hapAt(first)),
                ),
                false,
              )
      const reference =
        name === refName ? referenceSteps(items, first, last) : undefined
      fragments.push({
        name,
        start,
        end,
        trimmedBefore: start - runStart,
        trimmedAfter: runEnd - end,
        leadBp: span.leadBp,
        leadLoaded: span.leadLoaded,
        tailBp: span.tailBp,
        tailLoaded: span.tailLoaded,
        stretches,
        firstBreak,
        lastBreak,
        before,
        after,
        readBefore: hapAt(first) - runStart,
        readAfter: runEnd - lastEnd,
        reference,
      })
    }
  }
  const held = new Set(fragments.map(f => f.name))
  const recorded = records.map(d => ({ ...d, record: true }))
  return runsTables(
    fragments,
    refName,
    revisited,
    holes,
    [...jumping, ...recorded].filter(d => !held.has(d.name)),
  )
}

// The departures between two steps on R of one piece of a run with no
// fragment that lie over `window`, as walkCut computes them
function pieceDepartures(
  items: Item[],
  haps: Float64Array,
  window: { start: number; end: number },
  squeeze: (bp: number) => number,
) {
  const out: Omit<Departing, 'name'>[] = []
  let leave: { item: OnItem; i: number } | undefined
  let own = 0
  items.forEach((item, i) => {
    if (item.kind === 'on') {
      if (leave?.item.piece === item.piece) {
        const from = leave.item
        const p = positionAt(from, stepsIn(from) - 1)
        const leavePos = from.against ? from.path.starts[p]! : from.path.end(p)
        const q = positionAt(item, 0)
        const rejoinPos = item.against ? item.path.end(q) : item.path.starts[q]!
        const ends = {
          ownBp: own,
          leavePos,
          leaveReversed: from.against,
          rejoinPos,
          rejoinReversed: item.against,
        }
        const [lo, hi] = departureStretch(ends)
        if (
          !adjacentSteps(ends, squeeze) &&
          lo < window.end &&
          hi > window.start
        ) {
          out.push({
            ...ends,
            hapStart: haps[leave.i + 1]!,
            missingBp: 0,
            record: false,
          })
        }
      }
      leave = { item, i }
      own = 0
    } else {
      own += itemBp(item)
    }
  })
  return out
}

// The reference's nodes a fragment of its own steps on, in order
function referenceSteps(items: Item[], first: Step, last: Step) {
  const names: string[] = []
  const starts: number[] = []
  const lengths: number[] = []
  for (let i = first.i; i <= last.i; i++) {
    const item = items[i]!
    if (item.kind === 'on') {
      const from = i === first.i ? first.k : 0
      const to = i === last.i ? last.k : stepsIn(item) - 1
      for (let k = from; k <= to; k++) {
        const p = positionAt(item, k)
        names.push(String(idOf(item.path.handles[p]!)))
        starts.push(item.path.starts[p]!)
        lengths.push(item.path.end(p) - item.path.starts[p]!)
      }
    }
  }
  return { names, starts, lengths }
}

const byWalk = (refName: string) => (a: Fragment, b: Fragment) =>
  Number(b.name === refName) - Number(a.name === refName) ||
  (a.name < b.name ? -1 : a.name > b.name ? 1 : a.start - b.start)

function endsTable(n: number): WalkEnds {
  return {
    kinds: new Uint8Array(n),
    pos: new Float64Array(n),
    reversed: new Uint8Array(n),
    ownBp: new Float64Array(n),
    missingBp: new Float64Array(n),
  }
}

// A fragment's tables, its ends told apart by the drawn stretch
function runsTables(
  unsorted: Fragment[],
  refName: string,
  revisited: number,
  holes: Hole[],
  departing: Departing[],
): WalkRuns {
  const fragments = unsorted.toSorted(byWalk(refName))
  const names: string[] = []
  const starts: number[] = []
  const lengths: number[] = []
  let drawnStart = Infinity
  let drawnEnd = -Infinity
  for (const f of fragments) {
    if (f.reference) {
      for (let k = 0; k < f.reference.names.length; k++) {
        names.push(f.reference.names[k]!)
        starts.push(f.reference.starts[k]!)
        lengths.push(f.reference.lengths[k]!)
      }
      for (const s of f.stretches) {
        drawnStart = Math.min(drawnStart, s.lo)
        drawnEnd = Math.max(drawnEnd, s.hi)
      }
    }
  }
  const drawn = {
    refName,
    start: Number.isFinite(drawnStart) ? drawnStart : Number.NaN,
    end: Number.isFinite(drawnEnd) ? drawnEnd : Number.NaN,
  }
  const n = fragments.length
  const offsets = new Int32Array(n + 1)
  fragments.forEach((f, p) => {
    offsets[p + 1] = offsets[p]! + f.stretches.length
  })
  const total = offsets[n]!
  const stretches = {
    starts: new Float64Array(total),
    bp: new Float64Array(total),
    against: new Uint8Array(total),
    offBefore: new Float64Array(total),
    loaded: new Uint8Array(total),
  }
  const before = endsTable(n)
  const after = endsTable(n)
  const setEnd = (
    ends: WalkEnds,
    p: number,
    d: Departed | undefined,
    at: number,
    against: boolean,
    leading: boolean,
    readBp: number,
  ) => {
    const highSide = against === leading
    const pastDrawn = (pos: number) =>
      highSide ? pos >= drawn.end : pos <= drawn.start
    const kind: EndKind = d
      ? d.reversed === against && pastDrawn(d.pos)
        ? 'beyond'
        : 'elsewhere'
      : pastDrawn(at)
        ? 'edge'
        : 'ends'
    ends.kinds[p] = END_KINDS.indexOf(kind)
    ends.pos[p] = d ? d.pos : Number.NaN
    ends.reversed[p] = d?.reversed === true ? 1 : 0
    ends.ownBp[p] = d ? d.ownBp : kind === 'ends' ? readBp : Number.NaN
    ends.missingBp[p] = d ? d.missingBp : 0
  }
  fragments.forEach((f, p) => {
    f.stretches.forEach((s, k) => {
      const at = offsets[p]! + k
      stretches.starts[at] = s.lo
      stretches.bp[at] = s.hi - s.lo
      stretches.against[at] = s.against ? 1 : 0
      stretches.offBefore[at] = s.offBefore
      stretches.loaded[at] = s.loaded ? 1 : 0
    })
    const head = f.stretches[0]!
    const tail = f.stretches.at(-1)!
    setEnd(before, p, f.before, f.firstBreak, head.against, true, f.readBefore)
    setEnd(after, p, f.after, f.lastBreak, tail.against, false, f.readAfter)
  })
  return {
    refName,
    drawn,
    revisited,
    holes: holes
      .filter(h => h.refName === refName)
      .map(({ start, end }) => ({ start, end }))
      .sort((a, b) => a.start - b.start),
    reference: {
      names,
      starts: Float64Array.from(starts),
      lengths: Int32Array.from(lengths),
    },
    walks: {
      names: fragments.map(f => f.name),
      starts: Float64Array.from(fragments, f => f.start),
      ends: Float64Array.from(fragments, f => f.end),
      trimmedBefore: Float64Array.from(fragments, f => f.trimmedBefore),
      trimmedAfter: Float64Array.from(fragments, f => f.trimmedAfter),
      offsets,
      leadBp: Float64Array.from(fragments, f => f.leadBp),
      leadLoaded: Uint8Array.from(fragments, f => (f.leadLoaded ? 1 : 0)),
      tailBp: Float64Array.from(fragments, f => f.tailBp),
      tailLoaded: Uint8Array.from(fragments, f => (f.tailLoaded ? 1 : 0)),
      before,
      after,
    },
    stretches,
    departing: departingTable(departing),
  }
}

/**
 * A cut's tables read back as runs: each walk-table walk a fragment, its
 * steps on the reference grouped into stretches and those off it into runs,
 * and its ends from the departure table's departure leaving its last step on
 * the reference, or rejoining its first. `holes` are the walk file's on
 * `refName`, which the tables do not state. This is what WalkFile.getRunsForRange
 * must return for the same window and options.
 */
export function fragmentRuns(
  { nodes, walks, departures, revisited }: GraphTables,
  refName: string,
  holes: Hole[] = [],
): WalkRuns {
  const refStarts = walks.refStarts ?? new Float64Array(walks.steps.length)
  const leaving = new Map<string, number>()
  const rejoining = new Map<string, number>()
  departures.names.forEach((name, k) => {
    leaving.set(
      edgeKey(
        name,
        departures.leavePos[k]!,
        departures.leaveReversed[k] === 1,
        departures.hapStarts[k]!,
      ),
      k,
    )
    rejoining.set(
      edgeKey(
        name,
        departures.rejoinPos[k]!,
        departures.rejoinReversed[k] === 1,
        departures.hapStarts[k]! + departures.ownBp[k]!,
      ),
      k,
    )
  })
  const recorded = (k: number | undefined, leave: boolean) =>
    k === undefined
      ? undefined
      : {
          pos: leave ? departures.rejoinPos[k]! : departures.leavePos[k]!,
          reversed:
            (leave
              ? departures.rejoinReversed[k]
              : departures.leaveReversed[k]) === 1,
          ownBp: departures.ownBp[k]!,
          missingBp: departures.missingBp[k]!,
        }
  const fragments: Fragment[] = walks.names.map((name, p) => {
    const from = walks.offsets[p]!
    const to = walks.offsets[p + 1]!
    const stretches: Stretch[] = []
    let off = 0
    let offLoaded = true
    let offSteps = false
    let lead = 0
    let leadLoaded = true
    let firstOn = -1
    let lastOn = -1
    for (let s = from; s < to; s++) {
      const node = walks.steps[s]!
      const bp = nodes.lengths[node]!
      const placeholder = nodes.placeholders?.[node] === 1
      const at = refStarts[s]!
      if (Number.isNaN(at)) {
        if (firstOn < 0) {
          lead += bp
          leadLoaded &&= !placeholder
        } else {
          off += bp
          offLoaded &&= !placeholder
          offSteps = true
        }
      } else {
        const against =
          (walks.reversed[s]! ^ (nodes.reversed?.[node] ?? 0)) === 1
        const lo = at
        const hi = at + bp
        const held = stretches.at(-1)
        if (
          held &&
          !offSteps &&
          held.against === against &&
          (against ? hi === held.lo : lo === held.hi)
        ) {
          held.lo = Math.min(held.lo, lo)
          held.hi = Math.max(held.hi, hi)
        } else {
          stretches.push({
            lo,
            hi,
            against,
            offBefore: firstOn < 0 ? 0 : off,
            loaded: firstOn < 0 ? true : offLoaded,
          })
        }
        firstOn = firstOn < 0 ? s : firstOn
        lastOn = s
        off = 0
        offLoaded = true
        offSteps = false
      }
    }
    const stepAt = (s: number) => {
      const node = walks.steps[s]!
      const lo = refStarts[s]!
      return {
        lo,
        hi: lo + nodes.lengths[node]!,
        against: (walks.reversed[s]! ^ (nodes.reversed?.[node] ?? 0)) === 1,
      }
    }
    const head = stepAt(firstOn)
    const tail = stepAt(lastOn)
    const firstBreak = breakpoint(head.lo, head.hi, head.against, false)
    const lastBreak = breakpoint(tail.lo, tail.hi, tail.against, true)
    const start = walks.starts[p]!
    const end = walks.ends[p]!
    const trimmedBefore = walks.trimmedBefore?.[p] ?? 0
    const trimmedAfter = walks.trimmedAfter?.[p] ?? 0
    const tailBp = off
    const reference =
      name === refName
        ? {
            names: Array.from(
              walks.steps.subarray(from, to),
              node => nodes.names[node]!,
            ),
            starts: Array.from(refStarts.subarray(from, to)),
            lengths: Array.from(
              walks.steps.subarray(from, to),
              node => nodes.lengths[node]!,
            ),
          }
        : undefined
    return {
      name,
      start,
      end,
      trimmedBefore,
      trimmedAfter,
      leadBp: lead,
      leadLoaded,
      tailBp,
      tailLoaded: offLoaded,
      stretches,
      firstBreak,
      lastBreak,
      before: recorded(
        rejoining.get(edgeKey(name, firstBreak, head.against, start + lead)),
        false,
      ),
      after: recorded(
        leaving.get(edgeKey(name, lastBreak, tail.against, end - tailBp)),
        true,
      ),
      readBefore: trimmedBefore + lead,
      readAfter: trimmedAfter + tailBp,
      reference,
    }
  })
  const held = new Set(walks.names)
  const departing = departures.names.flatMap((name, k) =>
    held.has(name)
      ? []
      : [
          {
            name,
            hapStart: departures.hapStarts[k]!,
            ownBp: departures.ownBp[k]!,
            missingBp: departures.missingBp[k]!,
            leavePos: departures.leavePos[k]!,
            leaveReversed: departures.leaveReversed[k] === 1,
            rejoinPos: departures.rejoinPos[k]!,
            rejoinReversed: departures.rejoinReversed[k] === 1,
            record: departures.fromSteps[k] === 0,
          },
        ],
  )
  return runsTables(fragments, refName, revisited, holes, departing)
}

// Whether an end of END_KINDS index `kind` rejoins the reference
const endRejoins = (kind: number) => kind <= END_KINDS.indexOf('elsewhere')

/**
 * A fragment's stretches and runs as the events walkDepartures reads, two
 * stretches with no own bp between them that follow on, holes removed, as
 * one; with `ends`, a departure at either end that rejoins the reference
 * too, as a cut lists it
 */
export function runEvents(runs: WalkRuns, p: number, ends = false) {
  const { stretches, walks } = runs
  const squeeze = squeezer(
    runs.holes.map(h => ({ ...h, refName: runs.refName })),
    runs.refName,
  )
  const events: ChainEvent[] = []
  if (ends && endRejoins(walks.before.kinds[p]!)) {
    const pos = walks.before.pos[p]!
    events.push(
      {
        kind: 'on',
        lo: pos,
        hi: pos,
        reversed: walks.before.reversed[p] === 1,
      },
      { kind: 'own', bp: walks.before.ownBp[p]!, loaded: false },
    )
  }
  let held: Extract<ChainEvent, { kind: 'on' }> | undefined
  for (let k = walks.offsets[p]!; k < walks.offsets[p + 1]!; k++) {
    const lo = stretches.starts[k]!
    const hi = lo + stretches.bp[k]!
    const reversed = stretches.against[k] === 1
    const own = stretches.offBefore[k]!
    if (
      held &&
      own === 0 &&
      held.reversed === reversed &&
      squeeze(reversed ? held.lo : held.hi) === squeeze(reversed ? hi : lo)
    ) {
      held.lo = Math.min(held.lo, lo)
      held.hi = Math.max(held.hi, hi)
    } else {
      if (held) {
        events.push({
          kind: 'own',
          bp: own,
          loaded: stretches.loaded[k] === 1,
        })
      }
      held = { kind: 'on', lo, hi, reversed }
      events.push(held)
    }
  }
  if (ends && endRejoins(walks.after.kinds[p]!)) {
    const pos = walks.after.pos[p]!
    events.push(
      { kind: 'own', bp: walks.after.ownBp[p]!, loaded: false },
      { kind: 'on', lo: pos, hi: pos, reversed: walks.after.reversed[p] === 1 },
    )
  }
  return events
}

export interface RunDeparture {
  leavePos: number
  rejoinPos: number
  ownBp: number
  invertedBp: number
  loaded: boolean
}

/**
 * A fragment's departures in its own order, as walkDepartures reads a chain:
 * where it leaves a stretch along its way and rejoins the next, with the bp
 * of its own runs between and of the stretches it takes the other way,
 * `loaded` false where a run between holds a placeholder. Its way is the
 * one walkDepartures reads it in.
 */
export function runDepartures(runs: WalkRuns, p: number): RunDeparture[] {
  const { flip, crossings } = chainCrossings(runEvents(runs, p))
  return (flip ? crossings.toReversed() : crossings).map(c => ({
    leavePos: flip ? c.to : c.from,
    rejoinPos: flip ? c.from : c.to,
    ownBp: c.ownBp,
    invertedBp: c.invertedBp,
    loaded: c.loaded,
  }))
}
