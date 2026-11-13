import type { Abridged, WalkRow } from './walkRows.ts'

// A tile of gfa-to-tabix's walk format 8, whose rows from `kind` on are
// format 6's (that repo's docs/walks-format.md). A tile lists its `ref` rows,
// which joined in hapOffset order are the tile's reference path R; a `nodes` row with the
// lengths of the nodes its routes visit; `route` rows, each a way off R; a
// `walk` row per piece of every other path, naming the routes it takes in
// order and, where R holds a node more than once, the occurrence a step
// takes; and an `abridged` row per piece the writer moved to a detail
// record. A `ref` row's `passes` say how many times the reference passes each
// of its steps' nodes. A handle is 2 * id + 1 when reversed.

const NONE = -1

export interface Route {
  entry: number
  interior: number[]
  exit: number
}

export interface PlacedNode {
  refName: string
  start: number
  end: number
  rank: number
  // how many times the reference's paths pass the node, from the chunk's ref
  // rows; 0 for a node off its R, whether off the reference or on it only
  // where these rows do not hold it
  passes: number
  // a run of a piece's steps off R the cut did not load, placed on its contig
  placeholder?: true
}

const idOf = (handle: number) => (handle - (handle & 1)) / 2

// so an id fits an Int32Array, clear of the negative ids placeholders take
const ID_LIMIT = 2 ** 31

function checkedId(id: number) {
  if (id < 0 || id >= ID_LIMIT) {
    throw new WalkRowError(
      `a row holds node id ${id}, where walk format 8 allows ids from 0 to below 2^31`,
    )
  }
  return id
}

// the other strand's handle; arithmetic, which holds past 2^31
const otherStrand = (handle: number) => (handle & 1 ? handle - 1 : handle + 1)

// The handles a column of delta tokens names, each 2 * (id - previous id) +
// reverse, from `previous.id`, which ends on the last id read
function tokens(field: string, previous: { id: number }, out: number[]) {
  if (field === '.') {
    return out
  }
  let v = 0
  let neg = false
  for (let p = 0; p <= field.length; p++) {
    const ch = p < field.length ? field.charCodeAt(p) : 44
    if (ch === 44) {
      const d = neg ? -v : v
      const r = d & 1
      previous.id = checkedId(previous.id + (d - r) / 2)
      out.push(previous.id * 2 + r)
      v = 0
      neg = false
    } else if (ch === 45) {
      neg = true
    } else {
      v = v * 10 + (ch - 48)
    }
  }
  return out
}

function numbers(field: string) {
  return field === '.' ? [] : field.split(',').map(Number)
}

// A ref row's `passes`: how many times the reference's paths pass each
// step's node, `.` when every one is 1
function passCounts(cols: string[], n: number) {
  const field = cols[11]
  const counts =
    field === '.' ? new Array<number>(n).fill(1) : numbers(field ?? '')
  if (counts.length !== n || counts.some(c => !Number.isInteger(c) || c < 1)) {
    throw new WalkRowError(`a ref row's passes do not count its ${n} steps`)
  }
  return counts
}

const forwardEnds = (a: number, b: number) => 2 - (a & 1) - (b & 1)

// A link one way round, so a link and its reverse are one: the way with more
// forward ends, as a walk along the reference takes it, else the smaller of a
// and the flipped b first
export function canonicalLink(a: number, b: number): [number, number] {
  const forward = forwardEnds(a, b)
  return forward > 1 || (forward === 1 && a <= otherStrand(b))
    ? [a, b]
    : [otherStrand(b), otherStrand(a)]
}

// A handle token and the occurrence in R its `@k` names, NONE for the default
function parseCopy(token: string): [number, number] {
  const at = token.indexOf('@')
  return at < 0 ? [+token, NONE] : [+token.slice(0, at), +token.slice(at + 1)]
}

// The chunk's reference path R: its handles, the bp of its steps before each
// position, where each position starts on the reference, and each node's
// occurrences in R
export class RefPath {
  readonly handles: number[] = []
  readonly before: number[] = [0]
  readonly starts: number[] = []
  readonly occurrences = new Map<number, number[]>()

  add(handles: number[], lengths: number[], start: number) {
    for (let k = 0; k < handles.length; k++) {
      const h = handles[k]!
      const i = this.handles.length
      this.handles.push(h)
      this.before.push(this.before[i]! + lengths[k]!)
      this.starts.push(start)
      start += lengths[k]!
      const id = idOf(h)
      const list = this.occurrences.get(id)
      if (list) {
        list.push(i)
      } else {
        this.occurrences.set(id, [i])
      }
    }
  }

  // bp of R's steps before position i
  bp(i: number) {
    return this.before[i]!
  }

  first(handle: number) {
    return this.occurrences.get(idOf(handle))?.[0]
  }

  // occurrence k of the handle's node, counting R's positions holding it from 0
  copy(handle: number, k: number) {
    return this.occurrences.get(idOf(handle))?.[k]
  }

  nearest(handle: number, i: number) {
    const at = this.occurrences.get(idOf(handle))
    if (!at) {
      return undefined
    }
    let lo = 0
    let hi = at.length
    while (lo < hi) {
      const mid = (lo + hi) >> 1
      if (at[mid]! <= i) {
        lo = mid + 1
      } else {
        hi = mid
      }
    }
    const before = at[lo - 1]
    const after = at[lo]
    if (before === undefined || after === undefined) {
      return before ?? after
    }
    return after - i <= i - before ? after : before
  }

  // The position a stretch along R from position i, reading R the way `h`
  // does there, first holds `target`, i itself included
  reach(i: number, h: number, target: number) {
    const forward = this.handles[i] === h
    let found: number | undefined
    for (const j of this.occurrences.get(idOf(target)) ?? []) {
      if (forward && j >= i && this.handles[j] === target) {
        if (found === undefined || j < found) {
          found = j
        }
      } else if (
        !forward &&
        j <= i &&
        otherStrand(this.handles[j]!) === target
      ) {
        if (found === undefined || j > found) {
          found = j
        }
      }
    }
    return found
  }

  handleAt(i: number, forward: boolean) {
    return forward ? this.handles[i]! : otherStrand(this.handles[i]!)
  }

  // the position that starts at `start` on the reference
  startingAt(start: number) {
    const i = firstAtLeast(this.starts, start)
    return this.starts[i] === start ? i : undefined
  }

  endingAt(end: number) {
    const i = firstAtLeast(this.starts, end) - 1
    return i >= 0 && this.end(i) === end ? i : undefined
  }

  // the position `bp` of R's steps past position i
  past(i: number, bp: number) {
    const target = this.before[i]! + bp
    const j = firstAtLeast(this.before, target)
    return this.before[j] === target ? j : undefined
  }

  end(i: number) {
    return this.starts[i]! + this.before[i + 1]! - this.before[i]!
  }

  // from the first base to the last of the positions over `span`
  hull(span: { start: number; end: number }) {
    const at = firstAtLeast(this.starts, span.start)
    const first = at > 0 && this.end(at - 1) > span.start ? at - 1 : at
    const past = firstAtLeast(this.starts, span.end)
    return past > first
      ? { start: this.starts[first]!, end: this.end(past - 1) }
      : undefined
  }
}

function firstAtLeast(sorted: number[], value: number) {
  let lo = 0
  let hi = sorted.length
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    if (sorted[mid]! < value) {
      lo = mid + 1
    } else {
      hi = mid
    }
  }
  return lo
}

// What a walk row's traversal reports, in step order
interface Visitor {
  // one step: the first, or a route's landing, at position i of R
  step(handle: number, i: number): void
  // the steps along R after position `from` up to and including `to`
  stretch(from: number, to: number, forward: boolean): void
  // a route's steps off R, reversed and flipped when `flip`; false stops
  // the traversal
  route(index: number, flip: boolean): boolean
}

export class WalkRowError extends Error {}

const orientedEntry = (route: Route, flip: boolean) =>
  flip ? (route.exit === NONE ? NONE : otherStrand(route.exit)) : route.entry

const orientedExit = (route: Route, flip: boolean) =>
  flip ? (route.entry === NONE ? NONE : otherStrand(route.entry)) : route.exit

// A route's interior handle k, taken the way round `flip` says
const interiorAt = (route: Route, flip: boolean, k: number) =>
  flip
    ? otherStrand(route.interior[route.interior.length - 1 - k]!)
    : route.interior[k]!

// A walk row's choices, four numbers each: route index, flip, and the
// occurrences its `:k` and `@k` name, NONE for the defaults
function parseChoices(field: string) {
  const choices: number[] = []
  if (field !== '.') {
    let previous = 0
    const parts = [0, NONE, NONE]
    let part = 0
    let v = 0
    let neg = false
    let any = false
    for (let p = 0; p <= field.length; p++) {
      const ch = p < field.length ? field.charCodeAt(p) : 44
      if (ch === 44 || ch === 58 || ch === 64) {
        if (any) {
          parts[part] = neg ? -v : v
        }
        v = 0
        neg = false
        any = false
        if (ch === 44) {
          const d = parts[0]!
          const r = d & 1
          previous += (d - r) / 2
          choices.push(previous, r, parts[1]!, parts[2]!)
          parts[0] = 0
          parts[1] = NONE
          parts[2] = NONE
          part = 0
        } else {
          part = ch === 58 ? 1 : 2
        }
      } else if (ch === 45) {
        neg = true
      } else {
        v = v * 10 + (ch - 48)
        any = true
      }
    }
  }
  return choices
}

// Follows a walk row of `n` steps as docs/walks-format.md decodes it
function traverse(
  path: RefPath,
  routes: Route[],
  start: number,
  startCopy: number,
  choices: number[],
  n: number,
  visit: Visitor,
) {
  const routeOf = (index: number) => {
    const route = routes[index]
    if (!route) {
      throw new WalkRowError(`route ${index} is not in the chunk`)
    }
    return route
  }
  let emitted = 0
  let k = 0
  if (
    choices.length > 0 &&
    orientedEntry(routeOf(choices[0]!), choices[1] === 1) === NONE
  ) {
    if (!visit.route(choices[0]!, choices[1] === 1)) {
      return
    }
    emitted += routeOf(choices[0]!).interior.length
    k = 4
  }
  let h = start
  let i = startCopy === NONE ? path.first(h) : path.copy(h, startCopy)
  if (i === undefined) {
    throw new WalkRowError('a walk row starts off its chunk')
  }
  visit.step(h, i)
  emitted++
  // along R to the occurrence `copy` names, or to the first position holding
  // `target`
  const follow = (target: number, copy: number) => {
    const forward = path.handles[i!] === h
    const j =
      copy === NONE ? path.reach(i!, h, target) : path.copy(target, copy)
    if (
      j === undefined ||
      (forward ? j < i! : j > i!) ||
      path.handleAt(j, forward) !== target
    ) {
      throw new WalkRowError('a walk row runs off its chunk')
    }
    if (j !== i) {
      visit.stretch(i!, j, forward)
      emitted += Math.abs(j - i!)
    }
    i = j
    h = target
  }
  for (; k < choices.length; k += 4) {
    const index = choices[k]!
    const flip = choices[k + 1] === 1
    const route = routeOf(index)
    const entry = orientedEntry(route, flip)
    if (entry === NONE) {
      throw new WalkRowError('an open route after the first choice')
    }
    follow(entry, choices[k + 2]!)
    if (!visit.route(index, flip)) {
      return
    }
    emitted += route.interior.length
    const exit = orientedExit(route, flip)
    if (exit === NONE) {
      return
    }
    h = exit
    const copy = choices[k + 3]!
    i = copy === NONE ? path.nearest(h, i) : path.copy(h, copy)
    if (i === undefined) {
      throw new WalkRowError('a route lands off its chunk')
    }
    visit.step(h, i)
    emitted++
  }
  const tail = n - emitted
  if (tail < 0) {
    throw new WalkRowError(`a walk row decodes to more than its ${n} steps`)
  }
  if (tail > 0) {
    const forward = path.handles[i] === h
    const j = i + (forward ? tail : -tail)
    if (j < 0 || j >= path.handles.length) {
      throw new WalkRowError('a walk row runs off its chunk')
    }
    visit.stretch(i, j, forward)
  }
}

// A walk row's steps as node ids, reverse flags, where each lies on the
// reference, NaN off R, whether each runs against R at its pass, and its bp,
// -1 where `lengthOf` has none
function decodeWalk(
  path: RefPath,
  routes: Route[],
  cols: string[],
  choices: number[],
  lengthOf: (id: number) => number,
) {
  const n = +cols[8]!
  const ids = new Int32Array(n)
  const rev = new Uint8Array(n)
  const refStarts = new Float64Array(n).fill(Number.NaN)
  const against = new Uint8Array(n)
  const bp = new Int32Array(n)
  let at = 0
  const push = (h: number, i?: number) => {
    if (at >= n) {
      throw new WalkRowError('a walk row decodes to more steps than it states')
    }
    if (i === undefined) {
      bp[at] = lengthOf(idOf(h))
    } else {
      refStarts[at] = path.starts[i]!
      against[at] = path.handles[i] === h ? 0 : 1
      bp[at] = path.before[i + 1]! - path.before[i]!
    }
    ids[at] = idOf(h)
    rev[at++] = h & 1
  }
  const [start, startCopy] = parseCopy(cols[9]!)
  traverse(path, routes, start, startCopy, choices, n, {
    step: push,
    stretch(from, to, forward) {
      if (forward) {
        for (let j = from + 1; j <= to; j++) {
          push(path.handles[j]!, j)
        }
      } else {
        for (let j = from - 1; j >= to; j--) {
          push(otherStrand(path.handles[j]!), j)
        }
      }
    },
    route(index, flip) {
      const route = routes[index]!
      for (let k = 0; k < route.interior.length; k++) {
        push(interiorAt(route, flip, k))
      }
      return true
    },
  })
  if (at !== n) {
    throw new WalkRowError(`a walk row decodes to ${at} steps, not ${n}`)
  }
  return { ids, rev, refStarts, against, bp }
}

// A chunk as decoded, kept until the cut knows its nodes: tables() then
// places the cut's nodes off R
export interface Chunk {
  path: RefPath
  routes: Route[]
  // walk rows in path-name order, choices parsed
  walks: { cols: string[]; choices: number[] }[]
  lengths: Map<number, number>
  // for a detail record, the node lengths of its chunk's R, which it shares
  refLengths?: Map<number, number>
}

export interface DecodedChunk {
  // node bases from the rows' SQ:Z: tags, when asked for
  sequences: Map<number, string>
  // the reference's rows and the kept walk rows, steps decoded
  rows: WalkRow[]
  // the kept `abridged` rows, whose steps a detail record or a `blocks` row
  // gives
  abridged: AbridgedPiece[]
  // R's nodes with their positions, and the kept walks' nodes off R with
  // their lengths, placed by placeOffR
  nodes: Map<number, PlacedNode>
  // R's later passes over a node in `nodes`
  revisits: (PlacedNode & { id: number })[]
  chunk: Chunk
}

// The bases an SQ:Z: tag at or after column `from` joins, one node's per
// length, as `ids` names them
function sliceBases(
  cols: string[],
  from: number,
  ids: number[],
  lengths: number[],
  out: Map<number, string>,
) {
  const tag = cols.slice(from).find(c => c.startsWith('SQ:Z:'))
  if (tag === undefined) {
    return
  }
  let at = 5
  for (let k = 0; k < ids.length; k++) {
    const length = lengths[k]!
    out.set(ids[k]!, tag.slice(at, at + length))
    at += length
  }
}

// An `abridged` line, `path fragStart hapOffset piece nsteps hapBp detail`: a
// piece the writer moved to a detail record, with its span on R and its bp
// on its contig
function abridgedPiece(cols: string[], tileStart: number) {
  const count = (k: number, least = 0) => {
    const n = Number(cols[k])
    if (
      cols[k] === undefined ||
      cols[k] === '' ||
      !Number.isSafeInteger(n) ||
      n < least
    ) {
      throw new WalkRowError(`not an abridged row: ${cols.join('\t')}`)
    }
    return n
  }
  const abridged: Abridged = {
    steps: count(8),
    detail: {
      refName: cols[0]!,
      start: tileStart,
      name: `${tileStart}:${count(10)}`,
    },
  }
  return {
    chunkStart: tileStart,
    name: cols[4]!,
    fragStart: count(5),
    hapOffset: count(6),
    piece: count(7),
    hapBp: count(9, 1),
    span: { start: count(1), end: count(2) },
    abridged,
  }
}

export type AbridgedPiece = ReturnType<typeof abridgedPiece>

// A step that stands for steps of a piece the cut did not build, on R or
// off it: it has their bp and is in no fragment
export const SKIPPED = -(2 ** 31)

const DIGIT_0 = 48

// A `blocks` row's list, the piece in path order: each block as where it
// begins on R, its bp and whether it runs against R, the bp off R before
// each, the bp off R before the first and after the last, and the bp from
// the first block to the last
export function parseBlocks(cols: string[]) {
  const list = cols[7] ?? ''
  const begins: number[] = []
  const lengths: number[] = []
  const against: boolean[] = []
  const off: number[] = []
  let at = 0
  let bp = 0
  let before = 0
  let run = 0
  let value = 0
  let digits = 0
  let negative = false
  let inBlock = false
  let begin = 0
  let reversed = false
  for (let p = 0; p <= list.length; p++) {
    const ch = p < list.length ? list.charCodeAt(p) : 44
    const digit = ch - DIGIT_0
    if (digit >= 0 && digit <= 9) {
      value = value * 10 + digit
      digits++
    } else if (ch === 45 && !inBlock && digits === 0 && !negative) {
      negative = true
    } else if ((ch === 43 || ch === 45) && !inBlock && digits > 0) {
      begin = at + (negative ? -value : value)
      reversed = ch === 45
      inBlock = true
      value = 0
      digits = 0
    } else if (ch === 44 && inBlock && digits > 0) {
      before = begins.length === 0 ? run : before
      bp += value + (begins.length === 0 ? 0 : run)
      off.push(begins.length === 0 ? 0 : run)
      begins.push(begin)
      lengths.push(value)
      against.push(reversed)
      at = reversed ? begin - value : begin + value
      run = 0
      inBlock = false
      value = 0
      digits = 0
      negative = false
    } else if (ch === 44 && digits > 0 && !negative && run === 0 && value > 0) {
      run = value
      value = 0
      digits = 0
    } else {
      throw new WalkRowError(`not a blocks row's list: ${list.slice(0, 80)}`)
    }
  }
  if (begins.length === 0) {
    throw new WalkRowError(
      `a blocks row with no block: ${cols.join('\t').slice(0, 120)}`,
    )
  }
  return { before, after: run, begins, lengths, against, off, end: at, bp }
}

export type PieceBlocks = ReturnType<typeof parseBlocks>

const blockLow = (blocks: PieceBlocks, k: number) =>
  blocks.against[k]!
    ? blocks.begins[k]! - blocks.lengths[k]!
    : blocks.begins[k]!

const blockEnd = (blocks: PieceBlocks, k: number) =>
  blocks.against[k]!
    ? blocks.begins[k]! - blocks.lengths[k]!
    : blocks.begins[k]! + blocks.lengths[k]!

// Whether the departure before block k, from the block before it, can lie
// over `span`
const jumpOver = (
  blocks: PieceBlocks,
  k: number,
  span: { start: number; end: number },
) => {
  const leave = blockEnd(blocks, k - 1)
  const rejoin = blocks.begins[k]!
  return (
    Math.min(leave, rejoin) - 1 < span.end &&
    Math.max(leave, rejoin) + 1 > span.start
  )
}

// A record's own steps in an abridged piece: their stretch of the contig,
// whether they trail the piece, the record's end beside them, its leave where
// they trail and else its rejoin, and that end's edge key
export interface OwnRun {
  key: string
  start: number
  end: number
  trailing: boolean
  pos: number
  reversed: boolean
}

// Where a record's own steps lie in an abridged piece of its path, on their
// contig: nowhere, or in the piece's whole run off R at one end
export function ownStretch(
  piece: Pick<AbridgedPiece, 'name' | 'piece' | 'hapOffset' | 'hapBp'>,
  record: { hapStart: number; ownBp: number },
) {
  const pieceEnd = piece.hapOffset + piece.hapBp
  const start = Math.max(piece.hapOffset, record.hapStart)
  const end = Math.min(pieceEnd, record.hapStart + record.ownBp)
  const trailing = start > piece.hapOffset
  if (start < end && trailing === end < pieceEnd) {
    throw new WalkRowError(
      `the departure record of ${piece.name} with own steps over [${record.hapStart}, ${record.hapStart + record.ownBp}) of its contig lies ${trailing ? 'inside' : 'over all of'} piece ${piece.piece}, which covers [${piece.hapOffset}, ${pieceEnd}), where a record's own steps are a run at one end of a piece`,
    )
  }
  return start < end ? { start, end, trailing } : undefined
}

const pieceRow = (piece: AbridgedPiece, ids: number[], bp: number[]) => ({
  chunkStart: piece.chunkStart,
  name: piece.name,
  fragStart: piece.fragStart,
  hapOffset: piece.hapOffset,
  piece: piece.piece,
  n: ids.length,
  steps: {
    ids: Int32Array.from(ids),
    rev: new Uint8Array(ids.length),
    refStarts: new Float64Array(ids.length).fill(Number.NaN),
    against: new Uint8Array(ids.length),
    bp: Int32Array.from(bp),
  },
})

// The walk row of an abridged piece the cut builds no step of: its bp
export const bpRow = (piece: AbridgedPiece) =>
  pieceRow(piece, [SKIPPED], [piece.hapBp])

// The walk row of an abridged piece as its records' own runs, each a
// placeholder beside the step on R its record names, from its abridged row,
// the records and its tile's R alone
export function ownRunsRow(
  piece: AbridgedPiece,
  runs: OwnRun[],
  path: RefPath,
  placeholder: (start: number, end: number) => number,
) {
  const ids: number[] = []
  const bp: number[] = []
  const own: NonNullable<WalkRow['own']> = []
  const nodes = new Map<number, PlacedNode>()
  const add = (start: number, end: number, run?: OwnRun) => {
    const id = run ? placeholder(start, end) : SKIPPED
    if (run) {
      const atNodeEnd = run.trailing !== run.reversed
      const i = atNodeEnd ? path.endingAt(run.pos) : path.startingAt(run.pos)
      if (i === undefined) {
        throw new WalkRowError(
          `a departure record of ${piece.name} ${run.trailing ? 'leaves' : 'rejoins'} at ${run.pos}, where the reference path of the tile holding its own steps has no node ${atNodeEnd ? 'ending' : 'starting'}`,
        )
      }
      nodes.set(id, {
        refName: piece.name,
        start,
        end,
        rank: 1,
        passes: 0,
        placeholder: true,
      })
      own.push({
        key: run.key,
        leave: run.trailing,
        from: ids.length,
        to: ids.length + 1,
      })
    }
    ids.push(id)
    bp.push(end - start)
  }
  const leading = runs.find(run => !run.trailing)
  const trailing = runs.find(run => run.trailing)
  if (leading) {
    add(leading.start, leading.end, leading)
  }
  add(
    leading?.end ?? piece.hapOffset,
    trailing?.start ?? piece.hapOffset + piece.hapBp,
  )
  if (trailing) {
    add(trailing.start, trailing.end, trailing)
  }
  const row: WalkRow = {
    ...pieceRow(piece, ids, bp),
    abridged: piece.abridged,
    own,
  }
  return { row, nodes }
}

// The first and last blocks over `span`
export function blocksOver(
  blocks: PieceBlocks,
  span: { start: number; end: number },
) {
  let first = -1
  let last = -1
  for (let k = 0; k < blocks.begins.length; k++) {
    const lo = blockLow(blocks, k)
    if (lo < span.end && lo + blocks.lengths[k]! > span.start) {
      first = first < 0 ? k : first
      last = k
    }
  }
  return { first, last }
}

// The walk row of an abridged piece the cut does not load: blocks `core[0]`
// to `core[1]` whole, and of the rest its end steps and each departure that
// can lie over `over`, a SKIPPED step for what lies between. No departure
// across a SKIPPED step lies over `over`, so none a cut computes there is one
// it returns
export function unloadedRow(
  piece: AbridgedPiece,
  blocks: PieceBlocks,
  path: RefPath,
  over: { start: number; end: number },
  core: [number, number] | undefined,
  placeholder: (start: number, end: number) => number,
) {
  const { name, abridged } = piece
  const count = blocks.begins.length
  const from = new Int32Array(count)
  const to = new Int32Array(count)
  let low = Infinity
  let high = -Infinity
  for (let k = 0; k < count; k++) {
    const length = blocks.lengths[k]!
    const lo = blockLow(blocks, k)
    const i = path.startingAt(lo)
    const j = i === undefined ? undefined : path.past(i, length)
    if (
      i === undefined ||
      j === undefined ||
      j === i ||
      path.end(j - 1) !== lo + length
    ) {
      throw new WalkRowError(
        `block ${k} of ${name}'s blocks row, piece ${piece.piece}, covers no whole steps of its tile's reference path over [${lo}, ${lo + length})`,
      )
    }
    from[k] = i
    to[k] = j
    low = Math.min(low, lo)
    high = Math.max(high, lo + length)
  }
  if (low !== piece.span.start || high !== piece.span.end) {
    throw new WalkRowError(
      `${name}'s blocks row, piece ${piece.piece}, spans [${low}, ${high}) and its abridged row [${piece.span.start}, ${piece.span.end})`,
    )
  }
  const inCore = (k: number) =>
    core !== undefined && k >= core[0] && k <= core[1]
  const joined = (k: number) =>
    inCore(k - 1) || inCore(k) || jumpOver(blocks, k, over)
  const ids: number[] = []
  const rev: number[] = []
  const refStarts: number[] = []
  const against: number[] = []
  const bp: number[] = []
  const nodes = new Map<number, PlacedNode>()
  const { hapOffset } = piece
  let hap = hapOffset
  let skipped = 0
  let omitted = false
  const push = (
    id: number,
    reversed: number,
    refStart: number,
    opposed: number,
    length: number,
  ) => {
    ids.push(id)
    rev.push(reversed)
    refStarts.push(refStart)
    against.push(opposed)
    bp.push(length)
    hap += length
  }
  const flush = () => {
    if (skipped > 0 || omitted) {
      const length = skipped
      skipped = 0
      omitted = false
      push(SKIPPED, 0, Number.NaN, 0, length)
    }
  }
  const offR = (length: number) => {
    flush()
    const id = placeholder(hap, hap + length)
    nodes.set(id, {
      refName: name,
      start: hap,
      end: hap + length,
      rank: 1,
      passes: 0,
      placeholder: true,
    })
    push(id, 0, Number.NaN, 0, length)
  }
  const step = (i: number, reversed: boolean) => {
    flush()
    const h = path.handleAt(i, !reversed)
    push(
      idOf(h),
      h & 1,
      path.starts[i]!,
      reversed ? 1 : 0,
      path.end(i) - path.starts[i]!,
    )
  }
  const firstOf = (k: number) => (blocks.against[k]! ? to[k]! - 1 : from[k]!)
  const lastOf = (k: number) => (blocks.against[k]! ? from[k]! : to[k]! - 1)
  if (blocks.before > 0) {
    offR(blocks.before)
  }
  for (let k = 0; k < count; k++) {
    const reversed = blocks.against[k]!
    const run = blocks.off[k]!
    if (k > 0 && run > 0 && joined(k)) {
      offR(run)
    } else {
      skipped += run
    }
    const size = to[k]! - from[k]!
    const needFirst = k === 0 || joined(k)
    const needLast = k === count - 1 || joined(k + 1)
    if (inCore(k)) {
      for (
        let i = firstOf(k);
        reversed ? i >= from[k]! : i < to[k]!;
        i += reversed ? -1 : 1
      ) {
        step(i, reversed)
      }
    } else {
      const firstBp = path.end(firstOf(k)) - path.starts[firstOf(k)]!
      const lastBp = path.end(lastOf(k)) - path.starts[lastOf(k)]!
      if (needFirst || (size === 1 && needLast)) {
        step(firstOf(k), reversed)
      } else {
        skipped += firstBp
        omitted = true
      }
      if (size > 1) {
        skipped += blocks.lengths[k]! - firstBp - lastBp
        omitted ||= size > 2
        if (needLast) {
          step(lastOf(k), reversed)
        } else {
          skipped += lastBp
          omitted = true
        }
      }
    }
  }
  flush()
  if (blocks.after > 0) {
    offR(blocks.after)
  }
  const n = ids.length
  const row: WalkRow = {
    chunkStart: piece.chunkStart,
    name,
    fragStart: piece.fragStart,
    hapOffset,
    piece: piece.piece,
    n,
    steps: {
      ids: Int32Array.from(ids),
      rev: Uint8Array.from(rev),
      refStarts: Float64Array.from(refStarts),
      against: Uint8Array.from(against),
      bp: Int32Array.from(bp),
    },
    abridged,
  }
  return { row, nodes }
}

const refRows = (rows: string[][]) =>
  rows.filter(c => c[3] === 'ref').sort((a, b) => +a[6]! - +b[6]!)

// A chunk's R and its nodes' lengths, which its detail records share
export interface ChunkRef {
  path: RefPath
  lengths: Map<number, number>
}

export function chunkRef(rows: string[][]): ChunkRef {
  const path = new RefPath()
  const lengths = new Map<number, number>()
  for (const cols of refRows(rows)) {
    const handles = tokens(cols[9]!, { id: 0 }, [])
    const nodeLengths = numbers(cols[10]!)
    path.add(handles, nodeLengths, +cols[6]!)
    for (let k = 0; k < handles.length; k++) {
      lengths.set(idOf(handles[k]!), nodeLengths[k]!)
    }
  }
  return { path, lengths }
}

/**
 * Decodes one tile's rows, each split on tabs, the tile starting at
 * `tileStart` nominally. `keep` names the walk rows whose steps a cut wants,
 * undefined for every one. Given `ref`, the rows are a detail record decoded
 * against that tile's R, and the result leaves out R's rows, nodes, ends and
 * bases, which the tile's own decode gave.
 */
export function decodeChunk(
  rows: string[][],
  tileStart: number,
  keep: ((name: string) => boolean) | undefined,
  withSequences = false,
  ref?: ChunkRef,
): DecodedChunk {
  const path = ref?.path ?? new RefPath()
  const lengths = new Map<number, number>()
  const routes: Route[] = []
  const out: DecodedChunk = {
    sequences: new Map(),
    rows: [],
    abridged: [],
    nodes: new Map(),
    revisits: [],
    chunk: { path, routes, walks: [], lengths, refLengths: ref?.lengths },
  }
  const placed = (id: number) =>
    out.nodes.has(id) || ref?.lengths.has(id) === true
  const previous = { id: 0 }
  for (const cols of ref ? [] : refRows(rows)) {
    const handles = tokens(cols[9]!, { id: 0 }, [])
    const nodeLengths = numbers(cols[10]!)
    const passes = passCounts(cols, handles.length)
    const hapOffset = +cols[6]!
    const first = path.starts.length
    path.add(handles, nodeLengths, hapOffset)
    if (withSequences) {
      sliceBases(cols, 11, handles.map(idOf), nodeLengths, out.sequences)
    }
    let start = hapOffset
    for (let k = 0; k < handles.length; k++) {
      const id = idOf(handles[k]!)
      const length = nodeLengths[k]!
      lengths.set(id, length)
      const node = {
        refName: cols[4]!,
        start,
        end: start + length,
        rank: 0,
        passes: passes[k]!,
      }
      if (out.nodes.has(id)) {
        out.revisits.push({ id, ...node })
      } else {
        out.nodes.set(id, node)
      }
      start += length
    }
    if (keep === undefined || keep(cols[4]!)) {
      out.rows.push({
        chunkStart: tileStart,
        name: cols[4]!,
        fragStart: +cols[5]!,
        hapOffset,
        piece: +cols[7]!,
        n: handles.length,
        steps: {
          ids: Int32Array.from(handles, idOf),
          rev: Uint8Array.from(handles, h => h & 1),
          refStarts: Float64Array.from(path.starts.slice(first)),
          against: new Uint8Array(handles.length),
          bp: Int32Array.from(nodeLengths),
        },
      })
    }
  }
  for (const cols of rows) {
    if (cols[3] === 'nodes') {
      const ids = numbers(cols[4]!)
      const nodeLengths = numbers(cols[5]!)
      let id = 0
      for (let k = 0; k < ids.length; k++) {
        id = checkedId(id + ids[k]!)
        ids[k] = id
        lengths.set(id, nodeLengths[k]!)
      }
      if (withSequences) {
        sliceBases(cols, 6, ids, nodeLengths, out.sequences)
      }
    } else if (cols[3] === 'route') {
      const entry = tokens(cols[4]!, previous, [])
      const interior = tokens(cols[5]!, previous, [])
      const exit = tokens(cols[6]!, previous, [])
      routes.push({ entry: entry[0] ?? NONE, interior, exit: exit[0] ?? NONE })
    }
  }
  const routeOf = (index: number) => {
    const route = routes[index]
    if (!route) {
      throw new WalkRowError(`route ${index} is not in the chunk`)
    }
    return route
  }
  for (const cols of rows) {
    if (cols[3] === 'abridged') {
      const name = cols[4]!
      if (keep === undefined || keep(name)) {
        out.abridged.push(abridgedPiece(cols, tileStart))
      }
      continue
    }
    if (cols[3] !== 'walk') {
      continue
    }
    const choices = parseChoices(cols[10]!)
    out.chunk.walks.push({ cols, choices })
    const name = cols[4]!
    if (keep === undefined || keep(name)) {
      let interior = 0
      // the walk's nodes off R, by length; placeOffR gives their position
      for (let k = 0; k < choices.length; k += 4) {
        interior += routeOf(choices[k]!).interior.length
        for (const h of routeOf(choices[k]!).interior) {
          const id = idOf(h)
          if (!placed(id)) {
            const length = lengths.get(id) ?? 0
            out.nodes.set(id, {
              refName: name,
              start: 0,
              end: length,
              rank: 1,
              passes: 0,
            })
          }
        }
      }
      const steps = decodeWalk(
        path,
        routes,
        cols,
        choices,
        id => lengths.get(id) ?? ref?.lengths.get(id) ?? -1,
      )
      const n = steps.ids.length
      out.rows.push({
        chunkStart: tileStart,
        name,
        fragStart: +cols[5]!,
        hapOffset: +cols[6]!,
        piece: +cols[7]!,
        n,
        interior,
        steps,
        span: [+cols[1]!, +cols[2]!],
      })
    }
  }
  return out
}

/**
 * The positions of the cut's nodes off R, from the chunks a cut read: `kept`
 * holds the node ids the cut draws. A node off R takes its position from the
 * path whose name sorts first among those visiting it in the chunks read, at
 * its earliest offset, as format 1's node rows had it. A node on any read
 * chunk's R keeps its reference position, even where another chunk's route
 * steps on it.
 */
export function placeOffR(
  chunks: Chunk[],
  kept: Set<number>,
  place: (id: number, node: PlacedNode) => void,
) {
  const offR = new Set<number>()
  for (const id of kept) {
    if (!chunks.some(({ path }) => path.occurrences.has(id))) {
      offR.add(id)
    }
  }
  const best = new Map<number, PlacedNode>()
  for (const { path, routes, walks, lengths, refLengths } of chunks) {
    // routes holding a node of the cut off R
    const wanted = Uint8Array.from(routes, ({ interior }) =>
      interior.some(h => offR.has(idOf(h))) ? 1 : 0,
    )
    const lengthOf = (h: number) =>
      lengths.get(idOf(h)) ?? refLengths?.get(idOf(h)) ?? 0
    const bp = new Float64Array(routes.length).fill(-1)
    const routeBp = (index: number) => {
      if (bp[index]! < 0) {
        bp[index] = routes[index]!.interior.reduce((s, h) => s + lengthOf(h), 0)
      }
      return bp[index]!
    }
    // Walk rows run in path-name order, so a route's first taker is the
    // first-named path to reach its nodes here; each row is followed only as
    // far as the last route it is first to take
    const taken = new Uint8Array(routes.length)
    for (const { cols, choices } of walks) {
      let remaining = 0
      for (let k = 0; k < choices.length; k += 4) {
        const index = choices[k]!
        if (wanted[index]! !== 0 && taken[index] === 0) {
          remaining++
        }
      }
      if (remaining === 0) {
        continue
      }
      const name = cols[4]!
      let offset = +cols[6]!
      const [start, startCopy] = parseCopy(cols[9]!)
      traverse(path, routes, start, startCopy, choices, +cols[8]!, {
        step(h) {
          offset += lengthOf(h)
        },
        stretch(from, to, forward) {
          offset += forward
            ? path.bp(to + 1) - path.bp(from + 1)
            : path.bp(from) - path.bp(to)
        },
        route(index, flip) {
          if (wanted[index] === 0 || taken[index]! !== 0) {
            offset += routeBp(index)
            return true
          }
          taken[index] = 1
          const route = routes[index]!
          for (let k = 0; k < route.interior.length; k++) {
            const h = interiorAt(route, flip, k)
            const id = idOf(h)
            const length = lengthOf(h)
            if (offR.has(id)) {
              const held = best.get(id)
              if (
                !held ||
                name < held.refName ||
                (name === held.refName && offset < held.start)
              ) {
                best.set(id, {
                  refName: name,
                  start: offset,
                  end: offset + length,
                  rank: 1,
                  passes: 0,
                })
              }
            }
            offset += length
          }
          remaining--
          return remaining > 0
        },
      })
    }
  }
  for (const [id, node] of best) {
    place(id, node)
  }
}

// The reference's nodes in a chunk's ref rows, with their positions
export function referenceNodes(rows: string[][]) {
  const nodes: (PlacedNode & { id: number })[] = []
  for (const cols of rows) {
    if (cols[3] !== 'ref') {
      continue
    }
    const handles = tokens(cols[9]!, { id: 0 }, [])
    const nodeLengths = numbers(cols[10]!)
    const passes = passCounts(cols, handles.length)
    let start = +cols[6]!
    for (let k = 0; k < handles.length; k++) {
      const length = nodeLengths[k]!
      nodes.push({
        id: idOf(handles[k]!),
        refName: cols[4]!,
        start,
        end: start + length,
        rank: 0,
        passes: passes[k]!,
      })
      start += length
    }
  }
  return nodes
}
