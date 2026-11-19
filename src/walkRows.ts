import {
  RegionTooLargeError,
  bpString,
  panSNHaplotypeOrSample,
  panSNMatchesPrefix,
} from './util.ts'
import { SKIPPED } from './walkChunk.ts'
import {
  ATTACH_BP,
  adjacentSteps,
  attachesOwnSteps,
  departureStretch,
} from './walkDeparture.ts'
import { GAP_KINDS, walkGaps } from './walkGaps.ts'

import type { DepartureTable, GapTable, GraphTables } from './util.ts'
import type {
  DepartureEnds,
  DepartureRecord,
  Hole,
  WalkDeparture,
} from './walkDeparture.ts'
import type { WalkGap } from './walkGaps.ts'

// A walk row of gfa-to-tabix's format 6, its steps decoded with its chunk
// (walkChunk.ts). A haplotype path is cut into pieces, one per chunk it
// passes, and each piece's row takes the next piece index, so consecutive
// indices of one path join.
export interface WalkRow {
  name: string
  fragStart: number
  hapOffset: number
  piece: number
  n: number
  // of the n steps, those off R, in routes' interiors
  interior?: number
  // refStarts: where each step lies on the reference, NaN off R or unknown;
  // against: 1 where the step runs against R's orientation at its pass
  // bp: each step's bp, -1 where only its node states it
  steps: {
    ids: Int32Array
    rev: Uint8Array
    refStarts?: Float64Array
    against?: Uint8Array
    bp?: Int32Array
  }
  chunkStart: number
  // the piece's span on R, as its line states it
  span?: [number, number]
  // set on a piece built from its `blocks` row, which the cut did not load
  abridged?: Abridged
  own?: OwnSteps[]
}

// Steps `from` up to `to` as the own steps of the record whose leave, or
// rejoin, has the edge key `key`
export interface OwnSteps {
  key: string
  leave: boolean
  from: number
  to: number
}

// A piece the cut did not load: its steps in all and the detail unit that
// holds it
export interface Abridged {
  steps: number
  detail: { refName: string; start: number; name: string }
}

// Whether a step's id is a placeholder's: a run of an unloaded piece's
// steps off R, of known bp and no node identity
export const isPlaceholder = (id: number) => id < 0

export interface WalkNode {
  refName: string
  start: number
  end: number
  rank: number
  // how many times the reference's paths pass the node, from the ref rows the
  // cut read; 0 for a node it places off R, whether off the reference or on
  // it only where no ref row it read holds the node
  passes: number
  // the reference steps the node reversed at `start`, from the ref row that
  // places it; absent off the reference
  reversed?: boolean
  // on the reference inside the window plus the cut's context
  onWindow: boolean
  // no node of the graph: a run of an unloaded piece's steps off R, on
  // `refName`, its contig, from `start` to `end`
  placeholder?: true
}

export interface WalkWindow {
  refName: string
  start: number
  end: number
  context: number
  // keep each node's bases, from the SQ:Z: tags `--sequences` writes
  sequences?: boolean
}

export interface WalkHeader {
  // format:i:8; the reader refuses any other, and format 1 states none
  format?: number
  // the bp of reference a tile holds, tile:i:
  tile?: number
  references: string[]
  haplotypes: string[]
  holes: Hole[]
  // how many nodes the reference's paths pass more than once: 0 states
  // every position and departure in the file as a fact of the graph
  revisited?: number
  graph?: string
  leftout?: number
}

function integerTag(line: string, tag: string) {
  const match = new RegExp(`(?:^|\\t)${tag}:i:(\\d+)`).exec(line)
  return match ? +match[1]! : undefined
}

const INTEGER_TAGS: ('format' | 'tile' | 'revisited' | 'leftout')[] = [
  'format',
  'tile',
  'revisited',
  'leftout',
]

function stringTag(line: string, tag: string) {
  return new RegExp(`(?:^|\\t)${tag}:Z:([^\\t]+)`).exec(line)?.[1]
}

/**
 * What a walk file's header states: a `#walks` line with the format and the
 * tile size, `format:i:8\ttile:i:65536`, a `#reference\tGRCh38` line naming
 * the reference sample, a `#haplotype\tHG002#1` line per other haplotype with
 * rows and a `#hole` line per stretch of a reference sequence no line of the
 * reference covers
 */
export function walkHeader(lines: string[]): WalkHeader {
  const header: WalkHeader = {
    references: [],
    haplotypes: [],
    holes: [],
  }
  for (const line of lines) {
    const [key, value = '', start = '', end = ''] = line.split('\t')
    if (key === '#reference') {
      header.references.push(value)
    } else if (key === '#haplotype') {
      header.haplotypes.push(value)
    } else if (key === '#hole') {
      header.holes.push({ refName: value, start: +start, end: +end })
    } else if (key === '#walks') {
      for (const tag of INTEGER_TAGS) {
        const stated = integerTag(line, tag)
        if (stated !== undefined) {
          header[tag] ??= stated
        }
      }
      const graph = stringTag(line, 'graph')
      if (graph !== undefined) {
        header.graph ??= graph
      }
    }
  }
  return header
}

// Which walks a cut decodes: the reference's, and those `wanted` names by
// PanSN prefix at sample or haplotype depth; undefined decodes every one
export function walkNameFilter(
  wanted: string[] | undefined,
  referenceHaplotype: string | undefined,
) {
  return wanted === undefined || wanted.length === 0
    ? undefined
    : (name: string) =>
        panSNMatchesPrefix(name, referenceHaplotype) ||
        wanted.some(prefix => panSNMatchesPrefix(name, prefix))
}

/**
 * The walks a cut on `refName` decodes: its reference's and those `wanted`
 * names. A pggb P-line names a reference with no haplotype field
 * (`chm13#chr1`), whose walks are its sample's.
 */
export function walksKept(wanted: string[] | undefined, refName: string) {
  return walkNameFilter(wanted, panSNHaplotypeOrSample(refName))
}

// One stretch of a haplotype path: rows of consecutive piece index, decoded
export interface WalkRun {
  name: string
  fragStart: number
  hapOffset: number
  firstPiece: number
  // where each row after the first starts in `ids`, and the chunk each row
  // is filed under
  pieceStarts: number[]
  pieceChunks: number[]
  ids: Int32Array
  rev: Uint8Array
  // where each step lies on the reference, NaN off R or unknown, and whether
  // it runs against R's orientation at its pass
  refStarts: Float64Array
  against: Uint8Array
  // each step's bp, -1 where only its node states it
  bp?: Int32Array
  // the pieces the cut did not load, each with its steps' range
  unloaded?: { from: number; to: number; piece: Abridged }[]
  own?: OwnSteps[]
}

export function joinPieces(rows: WalkRow[]) {
  const byPath = new Map<string, WalkRow[]>()
  const runs: WalkRun[] = []
  for (const row of rows) {
    const key = `${row.name}\t${row.fragStart}`
    const list = byPath.get(key)
    if (list) {
      list.push(row)
    } else {
      byPath.set(key, [row])
    }
  }
  for (const pieces of byPath.values()) {
    pieces.sort((a, b) => a.piece - b.piece)
    let from = 0
    for (let i = 1; i <= pieces.length; i++) {
      if (
        i === pieces.length ||
        pieces[i]!.piece !== pieces[i - 1]!.piece + 1
      ) {
        runs.push(decodeRun(pieces.slice(from, i)))
        from = i
      }
    }
  }
  return runs
}

function decodeRun(pieces: WalkRow[]): WalkRun {
  const n = pieces.reduce((sum, p) => sum + p.n, 0)
  const ids = new Int32Array(n)
  const rev = new Uint8Array(n)
  const refStarts = new Float64Array(n).fill(Number.NaN)
  const against = new Uint8Array(n)
  const pieceStarts: number[] = []
  const pieceChunks: number[] = []
  const bp = pieces.some(p => p.steps.bp)
    ? new Int32Array(n).fill(-1)
    : undefined
  const unloaded: NonNullable<WalkRun['unloaded']> = []
  const own: OwnSteps[] = []
  let at = 0
  for (const p of pieces) {
    for (const steps of p.own ?? []) {
      own.push({ ...steps, from: at + steps.from, to: at + steps.to })
    }
    if (at > 0) {
      pieceStarts.push(at)
    }
    pieceChunks.push(p.chunkStart)
    if (p.steps.bp) {
      bp?.set(p.steps.bp, at)
    }
    if (p.abridged) {
      unloaded.push({ from: at, to: at + p.n, piece: p.abridged })
    }
    ids.set(p.steps.ids, at)
    rev.set(p.steps.rev, at)
    if (p.steps.refStarts) {
      refStarts.set(p.steps.refStarts, at)
    }
    if (p.steps.against) {
      against.set(p.steps.against, at)
    }
    at += p.n
  }
  const first = pieces[0]!
  return {
    name: first.name,
    fragStart: first.fragStart,
    hapOffset: first.hapOffset,
    firstPiece: first.piece,
    pieceStarts,
    pieceChunks,
    ids,
    rev,
    refStarts,
    against,
    ...(bp ? { bp } : {}),
    ...(unloaded.length > 0 ? { unloaded } : {}),
    ...(own.length > 0 ? { own } : {}),
  }
}

/**
 * One walk's stretch through the cut, from its first step on the drawn
 * stretch of the reference to its last, in consecutive pieces of the tiles
 * over the window and its context. `trimmedBefore` and `trimmedAfter` count
 * the bp of the walk those tiles hold but the cut kept out, before and after
 * it, so they depend on the tiles read.
 */
export interface WalkFragment {
  name: string
  fragStart: number
  hapStart: number
  hapEnd: number
  ids: Int32Array
  rev: Uint8Array
  // where each step lies on the reference at its pass, NaN off it, and
  // whether it runs against the reference's orientation there
  refStarts: Float64Array
  against: Uint8Array
  trimmedBefore: number
  trimmedAfter: number
  // the pieces with a step here that the cut did not load, in the walk's
  // order, and how many of the steps are theirs
  abridged?: Abridged[]
  unloadedSteps?: number
}

export interface WalkCutOptions {
  // keep the reference over the window and its context alone, and split each
  // walk where it steps onto the reference off them
  contained?: boolean
  // bp of each walk kept past either end of a contained span
  context?: number
  // the most bp of own sequence a departure gets its steps attached for
  attachBp?: number
  // a position with the reference's holes before it removed
  squeeze?: (bp: number) => number
  // the window and its context: every pass of a run's step over it widens
  // the drawn stretch to the step's node, and a departure that leaves or
  // rejoins on it gets its own steps. Not given, the reference nodes on the
  // window stand for it
  drawn?: { start: number; end: number }
  // runs of tiles read only for records' own steps, which give those steps
  // and nothing else
  attach?: WalkRun[]
  // the id of the placeholder for a walk's unloaded run over a stretch of
  // its contig; given it, a contained cut keeps `context` bp of a placeholder
  placeholder?: (name: string, start: number, end: number) => number
  // false lists no links
  links?: boolean
}

// A record whose own steps the cut attaches and may not hold: the cut lacks
// its leave piece or its rejoin piece, either of which can hold them
export interface Unattached {
  departure: Readonly<DepartureRecord>
}

// A placeholder a contained span's context ends inside: step `at`, of which
// the span keeps `keep` bp, its last where `tail`
interface Split {
  at: number
  keep: number
  tail: boolean
}

export const edgeKey = (
  name: string,
  pos: number,
  reversed: boolean,
  hap: number,
) => `${name}\t${pos}\t${reversed ? '-' : '+'}\t${hap}`

/**
 * The cut a window draws from its walks: the reference inside the window plus
 * its context, the drawn stretch, and each run from its first step on that
 * stretch to its last. The drawn stretch is the window and its context,
 * widened to the reference nodes over them, whole, and to every pass of a
 * run's step over them. Nothing is followed past the drawn stretch and no fragment runs on
 * through nodes other walks kept, so a haplotype's fragment does not depend on
 * which others are cut.
 *
 * A contained cut splits a run wherever it steps onto the reference off the
 * drawn stretch. Each piece between those steps that steps on the stretch is
 * kept from its first such step to its last, so a detour that comes back is
 * kept whole, and `context` bp of the walk past either end. At AMY1
 * haplotypes go from a 4 kb window to the other copies and back: the cut
 * above keeps 9,928 nodes there, the contained cut 182.
 *
 * Each departure comes back once: those inside one piece are computed from
 * its steps, and kept where their stretch overlaps the drawn stretch; those
 * between pieces, in different tiles or across W lines, are `records`, every
 * one of which comes back, as the reader read them over the window and its
 * context. A computed departure's leave and rejoin passes come from the ref
 * rows the cut read, 0 for a node it places off R. A departure within
 * `attachBp` of own sequence that leaves or rejoins on `drawn`, the base of
 * its leave step or of its rejoin step beside the breakpoint lying there,
 * gets its own steps, a 0 bp step included: a record from the leave and
 * rejoin pieces of `runs` or of `attach`, which give nothing else.
 * `unattached` lists those records lacking either piece, whose tiles the
 * reader fetches into `attach`. Any other departure states its ends and its
 * own bp alone.
 */
export function walkCut(
  readRuns: WalkRun[],
  nodes: Map<number, WalkNode>,
  refName: string,
  records: readonly Readonly<DepartureRecord>[],
  {
    contained = false,
    context = 0,
    attachBp = ATTACH_BP,
    squeeze = (bp: number) => bp,
    drawn,
    attach = [],
    placeholder,
    links: linked = true,
  }: WalkCutOptions = {},
) {
  // the drawn stretch: the window and its context, widened to the reference
  // nodes on it, so a window inside a hole in the reference has one
  let windowStart = drawn?.start ?? Infinity
  let windowEnd = drawn?.end ?? -Infinity
  for (const node of nodes.values()) {
    if (node.onWindow) {
      windowStart = Math.min(windowStart, node.start)
      windowEnd = Math.max(windowEnd, node.end)
    }
  }
  const length = (id: number) => {
    const node = nodes.get(id)
    return node === undefined ? 0 : node.end - node.start
  }
  const passesOf = (id: number) => nodes.get(id)?.passes ?? 0
  const onR = (run: WalkRun, i: number) => !Number.isNaN(run.refStarts[i])
  // bp of a run's steps before each, and after its last
  const prefixes = new Map<WalkRun, Float64Array>()
  const before = (run: WalkRun) => {
    let out = prefixes.get(run)
    if (!out) {
      out = new Float64Array(run.ids.length + 1)
      const stated = run.bp
      for (let i = 0; i < run.ids.length; i++) {
        const bp = stated === undefined ? -1 : stated[i]!
        out[i + 1] = out[i]! + (bp < 0 ? length(run.ids[i]!) : bp)
      }
      prefixes.set(run, out)
    }
    return out
  }
  // a node the reference passes twice is held at one pass, so the other
  // passes come from the steps placed on them
  if (drawn) {
    for (const run of readRuns) {
      const bp = before(run)
      for (let i = 0; i < run.ids.length; i++) {
        const start = run.refStarts[i]!
        if (start < drawn.end && start + bp[i + 1]! - bp[i]! > drawn.start) {
          windowStart = Math.min(windowStart, start)
          windowEnd = Math.max(windowEnd, start + bp[i + 1]! - bp[i]!)
        }
      }
    }
  }
  const attaches = (d: DepartureEnds) =>
    attachesOwnSteps(
      d,
      drawn ?? { start: windowStart, end: windowEnd },
      attachBp,
    )
  const stepOnWindow = (run: WalkRun, i: number) => {
    const start = run.refStarts[i]!
    const bp = before(run)
    return (
      !Number.isNaN(start) &&
      start < windowEnd &&
      start + bp[i + 1]! - bp[i]! > windowStart
    )
  }
  const spans: { run: WalkRun; first: number; last: number }[] = []
  // A contained run's spans: each stretch between steps onto the reference
  // off the window that has steps on it, from its first such step to its
  // last, and `context` bp of the walk past either: whole steps, but of a
  // placeholder only the bp left, which is where to split it
  const containedSpans = (run: WalkRun) => {
    const { ids } = run
    const bp = before(run)
    const found: typeof spans = []
    const splits: Split[] = []
    const take = (i: number, left: number, tail: boolean) => {
      const stepBp = bp[i + 1]! - bp[i]!
      if (stepBp > left && ids[i]! < 0 && ids[i] !== SKIPPED) {
        splits.push({ at: i, keep: left, tail })
      }
      return left - stepBp
    }
    let from = 0
    let first = -1
    let last = -1
    for (let i = 0; i <= ids.length; i++) {
      const on = i < ids.length && stepOnWindow(run, i)
      if (i === ids.length || ids[i] === SKIPPED || (onR(run, i) && !on)) {
        if (first >= 0) {
          for (let left = context; first > from && left > 0;) {
            left = take(--first, left, true)
          }
          for (let left = context; last < i - 1 && left > 0;) {
            left = take(++last, left, false)
          }
          found.push({ run, first, last })
        }
        from = i + 1
        first = -1
      } else if (on) {
        if (first < 0) {
          first = i
        }
        last = i
      }
    }
    return { found, splits }
  }
  // the links the steps returned walk, as pairs of handles, each once
  const links: number[] = []
  const linkedTo = new Map<number, Set<number>>()
  const link = (a: number, b: number) => {
    if (linked) {
      const held = linkedTo.get(a) ?? new Set<number>()
      linkedTo.set(a, held)
      if (!held.has(b)) {
        held.add(b)
        links.push(a, b)
      }
    }
  }
  // links each step of `run` from `from` to `to` to the one before it
  const walk = (run: WalkRun, from: number, to: number) => {
    const { ids, rev } = run
    for (let k = Math.max(from, 0) + 1; k <= to && k < ids.length; k++) {
      if (ids[k - 1] !== SKIPPED && ids[k] !== SKIPPED) {
        link(ids[k - 1]! * 2 + rev[k - 1]!, ids[k]! * 2 + rev[k]!)
      }
    }
  }
  // `run` with each split placeholder as two, the part within a span's
  // context and the rest
  const splitRun = (
    run: WalkRun,
    splits: Split[],
    mint: (name: string, start: number, end: number) => number,
  ) => {
    const bp = before(run)
    const n = run.ids.length + splits.length
    const ids = new Int32Array(n)
    const rev = new Uint8Array(n)
    const refStarts = new Float64Array(n).fill(Number.NaN)
    const against = new Uint8Array(n)
    const stepBp = new Int32Array(n)
    const shifted = new Int32Array(run.ids.length + 1)
    let s = 0
    let out = 0
    for (let i = 0; i < run.ids.length; i++) {
      shifted[i] = out
      const split = splits[s]
      const start = run.hapOffset + bp[i]!
      const end = run.hapOffset + bp[i + 1]!
      if (split?.at === i) {
        s++
        const edge = split.tail ? end - split.keep : start + split.keep
        const part = (a: number, b: number) => {
          const id = mint(run.name, a, b)
          nodes.set(id, {
            refName: run.name,
            start: a,
            end: b,
            rank: 1,
            passes: 0,
            onWindow: false,
            placeholder: true,
          })
          ids[out] = id
          stepBp[out++] = b - a
        }
        part(start, edge)
        part(edge, end)
      } else {
        ids[out] = run.ids[i]!
        rev[out] = run.rev[i]!
        refStarts[out] = run.refStarts[i]!
        against[out] = run.against[i]!
        stepBp[out++] = end - start
      }
    }
    shifted[run.ids.length] = out
    return {
      ...run,
      ids,
      rev,
      refStarts,
      against,
      bp: stepBp,
      pieceStarts: run.pieceStarts.map(at => shifted[at]!),
      unloaded: run.unloaded?.map(u => ({
        ...u,
        from: shifted[u.from]!,
        to: shifted[u.to]!,
      })),
      own: run.own?.map(steps => ({
        ...steps,
        from: shifted[steps.from]!,
        to: shifted[steps.to]!,
      })),
    }
  }
  const stretchSpan = (run: WalkRun) => {
    const bp = before(run)
    let first = -1
    let last = -1
    for (let i = 0; i < run.ids.length; i++) {
      const start = run.refStarts[i]!
      if (start < windowEnd && start + bp[i + 1]! - bp[i]! > windowStart) {
        if (first < 0) {
          first = i
        }
        last = i
      }
    }
    if (first >= 0) {
      spans.push({ run, first, last })
    }
  }
  const runs = readRuns.map(run => {
    const scan = contained ? containedSpans(run) : undefined
    const split =
      scan && placeholder && scan.splits.length > 0
        ? splitRun(run, scan.splits, placeholder)
        : run
    if (!scan) {
      stretchSpan(run)
    } else {
      for (const span of split === run
        ? scan.found
        : containedSpans(split).found) {
        spans.push(span)
      }
    }
    return split
  })
  const kept = new Set<number>()
  for (const [id, node] of nodes) {
    if (node.onWindow) {
      kept.add(id)
    }
  }
  for (const span of spans) {
    for (let i = span.first; i <= span.last; i++) {
      kept.add(span.run.ids[i]!)
    }
  }
  const runEnd = (run: WalkRun) => run.hapOffset + before(run)[run.ids.length]!
  // whether a fragment holds step `i` of `run`; a departure's own steps are
  // linked to its leave step and to its rejoin step where one does
  const spansOf = new Map<WalkRun, { first: number; last: number }[]>()
  for (const span of spans) {
    const held = spansOf.get(span.run) ?? []
    spansOf.set(span.run, held)
    held.push(span)
  }
  const returned = (run: WalkRun, i: number) =>
    spansOf.get(run)?.some(({ first, last }) => first <= i && i <= last) ===
    true
  const fragments: WalkFragment[] = []
  for (const { run, first, last } of spans) {
    const { ids, rev } = run
    const bp = before(run)
    walk(run, first, last)
    const hapStart = run.hapOffset + bp[first]!
    const hapEnd = run.hapOffset + bp[last + 1]!
    const held = (run.unloaded ?? []).filter(
      ({ from, to }) => from <= last && to > first,
    )
    const unloadedSteps = held.reduce(
      (sum, { from, to }) =>
        sum + Math.min(to, last + 1) - Math.max(from, first),
      0,
    )
    fragments.push({
      name: run.name,
      fragStart: run.fragStart,
      hapStart,
      hapEnd,
      ids: ids.subarray(first, last + 1),
      rev: rev.subarray(first, last + 1),
      refStarts: run.refStarts.subarray(first, last + 1),
      against: run.against.subarray(first, last + 1),
      trimmedBefore: hapStart - run.hapOffset,
      trimmedAfter: runEnd(run) - hapEnd,
      ...(held.length > 0
        ? { abridged: held.map(u => u.piece), unloadedSteps }
        : {}),
    })
  }
  // A step's breakpoint on R: where the path's sequence stops following R,
  // the node's end for a step along R and its start for one against it, and
  // the other way round for a rejoin
  const breakpoint = (run: WalkRun, i: number, leave: boolean) => {
    const reversed = run.against[i] === 1
    const start = run.refStarts[i]!
    return {
      pos:
        reversed === leave
          ? start
          : start + before(run)[i + 1]! - before(run)[i]!,
      reversed,
    }
  }
  const departures: WalkDeparture[] = []
  const attached = new Map<WalkRun, Set<number>>()
  const markAttached = (run: WalkRun, from: number, to: number) => {
    const marks = attached.get(run) ?? new Set<number>()
    attached.set(run, marks)
    for (let k = from; k < to; k++) {
      marks.add(k)
      kept.add(run.ids[k]!)
    }
  }
  const unloadedIn = (run: WalkRun, from: number, to: number) =>
    (run.unloaded ?? [])
      .filter(u => u.from < to && u.to > from)
      .map(u => u.piece)
  for (const run of runs) {
    if (run.name === refName) {
      continue
    }
    const { ids, rev, refStarts, against } = run
    const bp = before(run)
    const bounds = pieceBounds(run)
    for (let p = 0; p + 1 < bounds.length; p++) {
      let last = -1
      for (let i = bounds[p]!; i < bounds[p + 1]!; i++) {
        const start = refStarts[i]!
        if (Number.isNaN(start)) {
          continue
        }
        if (last >= 0) {
          const ownBp = bp[i]! - bp[last + 1]!
          const leaveReversed = against[last] === 1
          const rejoinReversed = against[i] === 1
          const leavePos = leaveReversed
            ? refStarts[last]!
            : refStarts[last]! + bp[last + 1]! - bp[last]!
          const rejoinPos = rejoinReversed ? start + bp[i + 1]! - bp[i]! : start
          const ends = {
            ownBp,
            leavePos,
            leaveReversed,
            rejoinPos,
            rejoinReversed,
          }
          if (!adjacentSteps(ends, squeeze)) {
            const [lo, hi] = departureStretch(ends)
            if (lo < windowEnd && hi > windowStart) {
              const d: WalkDeparture = {
                name: run.name,
                hapStart: run.hapOffset + bp[last + 1]!,
                ownBp,
                missingBp: 0,
                leavePos,
                leaveReversed,
                leavePasses: passesOf(ids[last]!),
                rejoinPos,
                rejoinReversed,
                rejoinPasses: passesOf(ids[i]!),
                fromSteps: true,
              }
              if (i > last + 1 && attaches(ends)) {
                d.steps = {
                  ids: ids.subarray(last + 1, i),
                  rev: rev.subarray(last + 1, i),
                }
                const held = unloadedIn(run, last + 1, i)
                if (held.length > 0) {
                  d.abridged = held
                }
                markAttached(run, last + 1, i)
                walk(
                  run,
                  returned(run, last) ? last : last + 1,
                  returned(run, i) ? i : i - 1,
                )
              }
              departures.push(d)
            }
          }
        }
        last = i
      }
    }
  }
  // A record's own steps: those after its leave step to the end of that
  // piece, then those before its rejoin step from the start of its piece. A
  // leave step is the last step on R of its piece and a rejoin step the
  // first of its, so the cut indexes those by name, breakpoint and contig
  // position
  // a piece's last step on R and where the piece ends, and its first and
  // where the piece starts
  const leaves = new Map<string, [WalkRun, number, number]>()
  const rejoins = new Map<string, [WalkRun, number, number]>()
  for (const run of [...runs, ...attach]) {
    const bp = before(run)
    const bounds = pieceBounds(run)
    for (let p = 0; p + 1 < bounds.length; p++) {
      let first = -1
      let last = -1
      for (let i = bounds[p]!; i < bounds[p + 1]!; i++) {
        if (onR(run, i)) {
          if (first < 0) {
            first = i
          }
          last = i
        }
      }
      if (first < 0) {
        continue
      }
      const leave = breakpoint(run, last, true)
      leaves.set(
        edgeKey(
          run.name,
          leave.pos,
          leave.reversed,
          run.hapOffset + bp[last + 1]!,
        ),
        [run, last, bounds[p + 1]!],
      )
      const rejoin = breakpoint(run, first, false)
      rejoins.set(
        edgeKey(
          run.name,
          rejoin.pos,
          rejoin.reversed,
          run.hapOffset + bp[first]!,
        ),
        [run, first, bounds[p]!],
      )
    }
  }
  for (const run of [...runs, ...attach]) {
    for (const { key, leave, from, to } of run.own ?? []) {
      if (leave) {
        leaves.set(key, [run, from - 1, to])
      } else {
        rejoins.set(key, [run, to, from])
      }
    }
  }
  const unattached: Unattached[] = []
  for (const record of records) {
    const d: WalkDeparture = { ...record }
    if (attaches(d)) {
      const ids: number[] = []
      const rev: number[] = []
      const held: Abridged[] = []
      const leave = leaves.get(
        edgeKey(d.name, d.leavePos, d.leaveReversed, d.hapStart),
      )
      if (leave) {
        const [run, i, end] = leave
        for (let k = i + 1; k < end; k++) {
          ids.push(run.ids[k]!)
          rev.push(run.rev[k]!)
        }
        markAttached(run, i + 1, end)
        walk(run, returned(run, i) ? i : i + 1, end - 1)
        for (const piece of unloadedIn(run, i + 1, end)) {
          held.push(piece)
        }
      }
      const rejoin = rejoins.get(
        edgeKey(d.name, d.rejoinPos, d.rejoinReversed, d.hapStart + d.ownBp),
      )
      if (rejoin) {
        const [run, i, start] = rejoin
        for (let k = start; k < i; k++) {
          ids.push(run.ids[k]!)
          rev.push(run.rev[k]!)
        }
        markAttached(run, start, i)
        walk(run, start, returned(run, i) ? i : i - 1)
        for (const piece of unloadedIn(run, start, i)) {
          held.push(piece)
        }
      }
      const oneRun =
        leave !== undefined &&
        rejoin !== undefined &&
        leave[0].fragStart === rejoin[0].fragStart &&
        pieceOf(leave[0], leave[1]) + 1 === pieceOf(rejoin[0], rejoin[1]) &&
        (leave[2] - 1 > leave[1] || returned(leave[0], leave[1])) &&
        (rejoin[2] < rejoin[1] || returned(rejoin[0], rejoin[1]))
      if (leave && rejoin && oneRun && ids.length > 0) {
        const [left, , end] = leave
        const [right, , start] = rejoin
        const last = left.ids[end - 1]
        const next = right.ids[start]
        if (
          last !== undefined &&
          next !== undefined &&
          last !== SKIPPED &&
          next !== SKIPPED
        ) {
          link(last * 2 + left.rev[end - 1]!, next * 2 + right.rev[start]!)
        }
      }
      if (ids.length > 0) {
        d.steps = { ids: Int32Array.from(ids), rev: Uint8Array.from(rev) }
        if (!oneRun) {
          d.breakAt = leave ? leave[2] - leave[1] - 1 : 0
        }
        if (held.length > 0) {
          d.abridged = held
        }
      }
      if (!leave || !rejoin) {
        unattached.push({ departure: record })
      }
    }
    departures.push(d)
  }
  departures.sort(
    (a, b) =>
      (a.name < b.name ? -1 : a.name > b.name ? 1 : 0) ||
      a.hapStart - b.hapStart ||
      a.leavePos - b.leavePos,
  )
  let gaps: WalkGap[] | undefined
  const gapsOf = () => {
    gaps ??= walkGaps({
      runs,
      attach,
      fragments,
      departures,
      refName,
      before,
      attached,
      breakpoint,
      spans,
    })
    return gaps
  }
  return { kept, fragments, departures, unattached, gapsOf, links }
}

const pieceOf = (run: WalkRun, step: number) =>
  run.firstPiece + run.pieceStarts.filter(start => start <= step).length

// Where each of a run's pieces starts, and where the last ends
function pieceBounds(run: WalkRun) {
  const bounds = [0]
  for (const start of run.pieceStarts) {
    bounds.push(start)
  }
  bounds.push(run.ids.length)
  return bounds
}

/**
 * Adds to `into` the nodes of decoded rows that walkCut keeps whatever else
 * the cut reads: those of `chunk.nodes` on the window, and those a row steps
 * on between two of its steps on the window, `start`-`end` with its context,
 * with no step onto the reference off the window between them in a contained
 * cut. A node placed later at a reference position off the window can split a
 * contained walk, so the count holds for the placements at the time. A
 * placeholder is no node, and is not counted.
 */
export function surelyKept(
  chunk: { nodes: Map<number, unknown>; rows: WalkRow[] },
  nodes: Map<number, WalkNode>,
  { contained, start, end }: { contained: boolean; start: number; end: number },
  into: Set<number>,
) {
  for (const id of chunk.nodes.keys()) {
    if (nodes.get(id)?.onWindow === true) {
      into.add(id)
    }
  }
  for (const { steps } of chunk.rows) {
    const { ids, refStarts } = steps
    let last = -1
    for (let i = 0; i < ids.length; i++) {
      const node = nodes.get(ids[i]!)
      const at = refStarts?.[i]
      const onR = at !== undefined && !Number.isNaN(at)
      if (
        node?.onWindow === true &&
        onR &&
        at < end &&
        at + node.end - node.start > start
      ) {
        if (last >= 0) {
          for (let j = last + 1; j < i; j++) {
            if (!isPlaceholder(ids[j]!)) {
              into.add(ids[j]!)
            }
          }
        }
        last = i
      } else if (contained && onR) {
        last = -1
      }
    }
  }
}

// A run of a walk's steps by its walk and its stretch of the contig, which no
// node's name is
export const runName = (name: string, start: number, end: number) =>
  `${name}:${start}-${end}`

// The places of a fragment's `n` steps it states: none from a host's fragment
// without `refStarts`, and none past its steps from one whose `refStarts` runs
// on
const statedPlaces = (placed: Float64Array | undefined, n: number) =>
  placed?.subarray(0, n) ?? new Float64Array()

// A node's name in the tables: its id, or for a placeholder its run's name
const nodeName = (id: number, node?: WalkNode) =>
  node?.placeholder ? runName(node.refName, node.start, node.end) : String(id)

// Whether placement `node` beats `held`: a pass of the reference over the
// window beats any other, then a reference position beats one off it, and
// among equals the path sorting first at its earliest start wins, so the order
// a cut reads chunks in does not matter
const outranks = (node: WalkNode, held: WalkNode) =>
  node.onWindow !== held.onWindow
    ? node.onWindow
    : node.rank === 0 &&
      (held.rank !== 0 ||
        node.refName < held.refName ||
        (node.refName === held.refName && node.start < held.start))

/**
 * The nodes and links of the chunks a cut reads, and the tables of the cut
 * they make with its walks. A node in two chunks is placed once.
 */
export class WalkGraph {
  readonly nodes = new Map<number, WalkNode>()
  readonly sequences = new Map<number, string>()
  private links: number[] = []

  private readonly window: WalkWindow
  private readonly revisited: number

  constructor(window: WalkWindow, revisited = Number.NaN) {
    this.window = window
    this.revisited = revisited
  }

  // A node with its position, unless the cut already holds it at one that
  // outranks it
  placeNode(id: number, node: Omit<WalkNode, 'onWindow'>) {
    const { window } = this
    const placed = {
      ...node,
      onWindow:
        node.rank === 0 &&
        node.refName === window.refName &&
        node.end > window.start - window.context &&
        node.start < window.end + window.context,
    }
    const held = this.nodes.get(id)
    if (!held || outranks(placed, held)) {
      this.nodes.set(id, placed)
    }
  }

  // A link between two endpoints, each id * 2 + 1 when reversed
  addLinkEnds(source: number, target: number) {
    this.links.push(source, target)
  }

  // Bases of nodes a chunk's rows give, for a window that asks for them
  addSequences(sequences: Map<number, string>) {
    if (this.window.sequences === true) {
      for (const [id, bases] of sequences) {
        this.sequences.set(id, bases)
      }
    }
  }

  /**
   * The cut as typed arrays, which the worker transfers rather than copies:
   * nodes by reference position, each link once in the order its L line
   * sorts, then a walk per fragment, the reference's first and the rest by
   * name, and the departures with their own steps. graphTablesGFA writes
   * them as the GFA they stand for, the same bytes for the same rows.
   */
  tables(
    kept: Set<number>,
    fragments: WalkFragment[],
    departures: WalkDeparture[],
    gaps: WalkGap[] = [],
  ): GraphTables {
    this.placeOffReference(kept, fragments, departures)
    const held: (WalkNode & { id: number })[] = []
    for (const id of kept) {
      const node = this.nodes.get(id)
      if (node) {
        held.push({ ...node, id })
      }
    }
    held.sort((x, y) =>
      x.refName < y.refName
        ? -1
        : x.refName > y.refName
          ? 1
          : x.start - y.start || x.id - y.id,
    )
    const index = new Map<number, number>()
    const names: string[] = []
    const refNames: string[] = []
    const refIndex = new Map<string, number>()
    const lengths = new Int32Array(held.length)
    const refs = new Int32Array(held.length)
    const starts = new Float64Array(held.length)
    const ranks = new Int32Array(held.length)
    const passes = new Int32Array(held.length)
    const placeholders = Uint8Array.from(held, node =>
      node.placeholder ? 1 : 0,
    )
    const reversedOnReference = Uint8Array.from(held, node =>
      node.reversed === true ? 1 : 0,
    )
    held.forEach((node, i) => {
      const { id, refName, start, end, rank, passes: n } = node
      index.set(id, i)
      names.push(nodeName(id, node))
      let ref = refIndex.get(refName)
      if (ref === undefined) {
        ref = refNames.length
        refIndex.set(refName, ref)
        refNames.push(refName)
      }
      lengths[i] = end - start
      refs[i] = ref
      starts[i] = start
      ranks[i] = rank
      passes[i] = n
    })
    // a kept node with no row of its own is named, and not drawn
    const nodeIndex = (id: number) => {
      let i = index.get(id)
      if (i === undefined) {
        i = names.length
        index.set(id, i)
        names.push(nodeName(id, this.nodes.get(id)))
      }
      return i
    }

    const byLine = new Map<string, number>()
    for (let i = 0; i < this.links.length; i += 2) {
      const a = this.links[i]!
      const b = this.links[i + 1]!
      const source = Math.floor(a / 2)
      const target = Math.floor(b / 2)
      if (kept.has(source) && kept.has(target)) {
        byLine.set(
          `${source}\t${a & 1 ? '-' : '+'}\t${target}\t${b & 1 ? '-' : '+'}`,
          i,
        )
      }
    }
    const lines = [...byLine.keys()].sort()
    const from = new Int32Array(lines.length)
    const to = new Int32Array(lines.length)
    const strands = new Uint8Array(lines.length)
    lines.forEach((line, k) => {
      const i = byLine.get(line)!
      const a = this.links[i]!
      const b = this.links[i + 1]!
      from[k] = nodeIndex(Math.floor(a / 2))
      to[k] = nodeIndex(Math.floor(b / 2))
      strands[k] = (a & 1) | ((b & 1) << 1)
    })

    const walks = this.referenceFirst(fragments)
    // a piece counts on the first walk holding it, else on the first
    // departure
    const counted = new Set<Abridged>()
    const stepsOf = (pieces: Abridged[] | undefined) => {
      let steps = pieces ? 0 : Number.NaN
      for (const piece of pieces ?? []) {
        if (!counted.has(piece)) {
          counted.add(piece)
          steps += piece.steps
        }
      }
      return steps
    }
    const abridgedSteps = Float64Array.from(walks, f => stepsOf(f.abridged))
    const offsets = new Int32Array(walks.length + 1)
    walks.forEach(({ ids }, p) => {
      offsets[p + 1] = offsets[p]! + ids.length
    })
    const steps = new Int32Array(offsets[walks.length]!)
    const reversed = new Uint8Array(steps.length)
    const refStarts = new Float64Array(steps.length).fill(Number.NaN)
    walks.forEach(({ ids, rev, refStarts: placed }, p) => {
      const at = offsets[p]!
      for (let j = 0; j < ids.length; j++) {
        steps[at + j] = nodeIndex(ids[j]!)
      }
      reversed.set(rev, at)
      refStarts.set(statedPlaces(placed, ids.length), at)
    })
    return {
      nodes: {
        names,
        lengths,
        refs,
        starts,
        ranks,
        refNames,
        passes,
        reversed: reversedOnReference,
        placeholders,
      },
      links: { from, to, strands },
      walks: {
        names: walks.map(f => f.name),
        starts: Float64Array.from(walks, f => f.hapStart),
        ends: Float64Array.from(walks, f => f.hapEnd),
        offsets,
        steps,
        reversed,
        refStarts,
        trimmedBefore: Float64Array.from(walks, f => f.trimmedBefore),
        trimmedAfter: Float64Array.from(walks, f => f.trimmedAfter),
        abridgedSteps,
      },
      departures: departureTable(departures, nodeIndex, stepsOf),
      gaps: gapTable(gaps, new Map(walks.map((f, p) => [f, p])), nodeIndex),
      revisited: this.revisited,
    }
  }

  private placeOffReference(
    kept: Set<number>,
    fragments: WalkFragment[],
    departures: WalkDeparture[],
  ) {
    const lengthOf = (id: number) => {
      const node = this.nodes.get(id)
      return node ? node.end - node.start : 0
    }
    const placed = new Map<number, WalkNode>()
    const visit = (name: string, ids: Int32Array, start: number) => {
      let at = start
      for (const id of ids) {
        const node = this.nodes.get(id)
        const length = lengthOf(id)
        const held = placed.get(id)
        if (
          node &&
          node.rank !== 0 &&
          !node.placeholder &&
          kept.has(id) &&
          (!held ||
            name < held.refName ||
            (name === held.refName && at < held.start))
        ) {
          placed.set(id, {
            ...node,
            refName: name,
            start: at,
            end: at + length,
          })
        }
        at += length
      }
    }
    for (const f of fragments) {
      visit(f.name, f.ids, f.hapStart)
    }
    for (const d of departures) {
      if (d.steps) {
        const { ids } = d.steps
        const leaving = ids.subarray(0, d.breakAt ?? ids.length)
        const rejoining = ids.subarray(leaving.length)
        const rejoiningBp = rejoining.reduce((sum, id) => sum + lengthOf(id), 0)
        visit(d.name, leaving, d.hapStart)
        visit(d.name, rejoining, d.hapStart + d.ownBp - rejoiningBp)
      }
    }
    for (const [id, node] of placed) {
      this.nodes.set(id, node)
    }
  }

  // The view anchors on the first walk when no walk names the assembly, as
  // gbz-base's cut has it; in file order that was whichever haplotype the
  // chunk listed first, CHM13 on a GRCh38 window
  private referenceFirst(fragments: WalkFragment[]) {
    const { refName } = this.window
    return [...fragments].sort(
      (a, b) =>
        Number(b.name === refName) - Number(a.name === refName) ||
        (a.name < b.name ? -1 : a.name > b.name ? 1 : a.hapStart - b.hapStart),
    )
  }
}

function gapTable(
  gaps: WalkGap[],
  walkIndex: Map<WalkFragment, number>,
  nodeIndex: (id: number) => number,
): GapTable {
  const segments = gaps.flatMap(g => g.segments)
  const segmentOffsets = new Int32Array(gaps.length + 1)
  gaps.forEach((g, k) => {
    segmentOffsets[k + 1] = segmentOffsets[k]! + g.segments.length
  })
  const stepOffsets = new Int32Array(segments.length + 1)
  segments.forEach((x, k) => {
    stepOffsets[k + 1] = stepOffsets[k]! + x.ids.length
  })
  const steps = new Int32Array(stepOffsets[segments.length]!)
  const reversed = new Uint8Array(steps.length)
  segments.forEach((x, k) => {
    const at = stepOffsets[k]!
    x.ids.forEach((id, j) => {
      steps[at + j] = nodeIndex(id)
    })
    reversed.set(x.rev, at)
  })
  return {
    names: gaps.map(g => g.name),
    walks: Int32Array.from(gaps, g =>
      g.fragment ? (walkIndex.get(g.fragment) ?? -1) : -1,
    ),
    edges: Uint8Array.from(gaps, g => g.edge),
    starts: Float64Array.from(gaps, g => g.start),
    ends: Float64Array.from(gaps, g => g.end),
    segmentOffsets,
    segments: {
      kinds: Uint8Array.from(segments, x => GAP_KINDS.indexOf(x.kind)),
      bp: Float64Array.from(segments, x => x.bp),
      missingBp: Float64Array.from(segments, x => x.missingBp),
      before: Float64Array.from(segments, x => x.before),
      after: Float64Array.from(segments, x => x.after),
      stepOffsets,
      steps,
      reversed,
    },
  }
}

function departureTable(
  departures: WalkDeparture[],
  nodeIndex: (id: number) => number,
  stepsOf: (pieces: Abridged[] | undefined) => number,
): DepartureTable {
  const n = departures.length
  const offsets = new Int32Array(n + 1)
  departures.forEach((d, k) => {
    offsets[k + 1] = offsets[k]! + (d.steps?.ids.length ?? 0)
  })
  const steps = new Int32Array(offsets[n]!)
  const reversed = new Uint8Array(steps.length)
  departures.forEach((d, k) => {
    if (d.steps) {
      const at = offsets[k]!
      for (let j = 0; j < d.steps.ids.length; j++) {
        steps[at + j] = nodeIndex(d.steps.ids[j]!)
      }
      reversed.set(d.steps.rev, at)
    }
  })
  return {
    names: departures.map(d => d.name),
    hapStarts: Float64Array.from(departures, d => d.hapStart),
    ownBp: Float64Array.from(departures, d => d.ownBp),
    missingBp: Float64Array.from(departures, d => d.missingBp),
    leavePos: Float64Array.from(departures, d => d.leavePos),
    leaveReversed: Uint8Array.from(departures, d => (d.leaveReversed ? 1 : 0)),
    leavePasses: Int32Array.from(departures, d => d.leavePasses),
    rejoinPos: Float64Array.from(departures, d => d.rejoinPos),
    rejoinReversed: Uint8Array.from(departures, d =>
      d.rejoinReversed ? 1 : 0,
    ),
    rejoinPasses: Int32Array.from(departures, d => d.rejoinPasses),
    fromSteps: Uint8Array.from(departures, d => (d.fromSteps ? 1 : 0)),
    breakAt: Int32Array.from(departures, d => d.breakAt ?? -1),
    offsets,
    steps,
    reversed,
    abridgedSteps: Float64Array.from(departures, d => stepsOf(d.abridged)),
  }
}

// How a notice names the budgets, such as the config slots that set them
export interface BudgetNames {
  steps: string
  bytes: string
}

export const BUDGET_NAMES: BudgetNames = {
  steps: 'the step budget',
  bytes: 'a byte budget of',
}

// The reader decodes a step along R several times as fast as one off it, so
// the read budget counts it as a quarter
const ALONG_R_PER_STEP = 4

export const readSteps = (row: Pick<WalkRow, 'n' | 'interior'>) => {
  const interior = row.interior ?? 0
  return interior + (row.n - interior) / ALONG_R_PER_STEP
}

/**
 * The zoom-in notice for a cut whose walk rows hold more steps than
 * `readPerKept` times `stepBudget`, or undefined when they fit. A window reads
 * whole tiles, so the span that fits is counted in tiles of `tile` bp from
 * the window's start; when the tiles a window there always reads are over on
 * their own, no zoom fits and the notice asks for fewer
 * haplotypes. A cut for every haplotype is also offered fewer where the
 * reference and one haplotype, at the walks' average steps, would fit the
 * whole window.
 */
export function stepBudgetError(
  rows: Pick<WalkRow, 'n' | 'interior' | 'chunkStart' | 'name'>[],
  stepBudget: number,
  readPerKept: number,
  window: { start: number; end: number },
  tile: number,
  filtered: boolean,
  names = BUDGET_NAMES,
) {
  const budget = stepBudget * readPerKept
  const over = readStepsOver(rows, budget, window, tile)
  const fewerFit =
    over !== undefined &&
    !filtered &&
    over.walks > 2 &&
    (2 * over.total) / over.walks <= budget
  return over
    ? new RegionTooLargeError(
        over.fitsBp > 0
          ? `Zoom in to about ${bpString(over.fitsBp)}${fewerFit ? ', or choose fewer haplotypes,' : ''} to see the graph`
          : `Too many haplotype steps here to draw ${filtered ? 'these haplotypes' : 'every haplotype'} (${Math.round(over.read).toLocaleString()} steps read, a step along the reference counting a quarter, over ${readPerKept} × ${names.steps} of ${stepBudget.toLocaleString()}); choose fewer haplotypes`,
        budget,
        window.end - window.start,
        Math.max(over.fitsBp, 1),
      )
    : undefined
}

export function countStepsError(
  rows: Pick<WalkRow, 'n' | 'interior' | 'chunkStart' | 'name'>[],
  stepBudget: number,
  readPerKept: number,
  window: { start: number; end: number },
  tile: number,
  names = BUDGET_NAMES,
) {
  const over = readStepsOver(rows, stepBudget * readPerKept, window, tile)
  return over
    ? new RegionTooLargeError(
        `Too many haplotype steps in this window to count: its walk rows hold ${Math.round(over.total).toLocaleString()} steps, a step along the reference counting a quarter, over ${readPerKept} × ${names.steps} of ${stepBudget.toLocaleString()}. ${over.fitsBp > 0 ? `Count a window of about ${bpString(over.fitsBp)} or less, or leave this one uncounted` : `Its first tile alone is over, so leave this window uncounted or raise ${names.steps}`}`,
        stepBudget * readPerKept,
        window.end - window.start,
        Math.max(over.fitsBp, 1),
      )
    : undefined
}

function readStepsOver(
  rows: Pick<WalkRow, 'n' | 'interior' | 'chunkStart' | 'name'>[],
  budget: number,
  window: { start: number; end: number },
  tile: number,
) {
  let total = 0
  const perChunk = new Map<number, number>()
  const walks = new Set<string>()
  for (const row of rows) {
    const steps = readSteps(row)
    total += steps
    perChunk.set(row.chunkStart, (perChunk.get(row.chunkStart) ?? 0) + steps)
    walks.add(row.name)
  }
  let read = 0
  let fitsBp = 0
  for (const cs of [...perChunk.keys()].sort((a, b) => a - b)) {
    read += perChunk.get(cs)!
    if (read > budget) {
      break
    }
    fitsBp = Math.max(fitsBp, cs + tile - window.start)
  }
  return total > budget ? { total, read, fitsBp, walks: walks.size } : undefined
}

/**
 * The zoom-in notice for a cut that keeps more steps than `budget`, the
 * fragments' and the departures' attached steps together, or undefined when
 * it fits. The span that fits is the window scaled by how far over the cut
 * is. A step built for a piece the cut did not load counts as any other.
 */
export function keptStepsError(
  fragments: WalkFragment[],
  departures: WalkDeparture[],
  budget: number,
  window: { start: number; end: number },
  filtered: boolean,
) {
  const total =
    fragments.reduce((sum, f) => sum + f.ids.length, 0) +
    departures.reduce((sum, d) => sum + (d.steps?.ids.length ?? 0), 0)
  return keptOverError(total, budget, window, filtered)
}

// keptStepsError's notice for a cut that keeps more nodes than `limit`, which
// it counted on the window before `readTo` when it stopped reading there
export function keptNodesError(
  nodes: number,
  limit: number,
  window: { start: number; end: number },
  filtered: boolean,
  readTo = window.end,
) {
  return keptOverError(nodes, limit, window, filtered, readTo)
}

function keptOverError(
  kept: number,
  limit: number,
  window: { start: number; end: number },
  filtered: boolean,
  readTo = window.end,
) {
  if (kept <= limit) {
    return undefined
  }
  const windowBp = window.end - window.start
  const readBp = Math.min(readTo, window.end) - window.start
  const fitsBp = Math.max(1, Math.floor((readBp * limit) / kept))
  return new RegionTooLargeError(
    `Zoom in to about ${bpString(fitsBp)}${filtered ? '' : ', or choose fewer haplotypes,'} to see the graph`,
    limit,
    windowBp,
    fitsBp,
  )
}

/**
 * The zoom-in notice for a cut whose walk read would fetch more than
 * `budget` compressed bytes, or undefined when it fits. `bytesTo(end)` is the
 * bytes of the tiles a read from the window's first base to `end` fetches;
 * the span that fits grows a tile at a time from the tile holding the
 * window's start, and when the tiles that reaches are over on their own, no
 * zoom fits. `attaching` marks bytes the tiles fetched for records' own
 * steps take past the budget.
 */
export function byteBudgetError(
  bytesTo: (end: number) => number,
  budget: number,
  window: { start: number; end: number },
  tile: number,
  names = BUDGET_NAMES,
  attaching = false,
) {
  const over = bytesOver(bytesTo, budget, window, tile)
  const windowBp = window.end - window.start
  return over === undefined
    ? undefined
    : over.fitsBp > 0
      ? new RegionTooLargeError(
          `Zoom in to about ${bpString(over.fitsBp)} to see the graph`,
          budget,
          windowBp,
          over.fitsBp,
          { bytes: over.bytes, zoomHelps: true },
        )
      : new RegionTooLargeError(
          attaching
            ? `Too many departures here whose own steps the cut must fetch, even zoomed in (${megabytes(over.least)} against ${names.bytes} ${megabytes(budget)}); choose fewer haplotypes or a smaller attachBp`
            : `Too much graph here to fetch, even zoomed in (${megabytes(over.least)} against ${names.bytes} ${megabytes(budget)}). A satellite or rDNA array usually causes this: the graph folds the repeat into a few nodes that each haplotype loops through many times`,
          budget,
          windowBp,
          1,
          { bytes: over.bytes },
        )
}

/**
 * byteBudgetError's refusal for counts, which draw nothing and choose no
 * haplotypes: what the read would fetch, and a narrower window to count
 * instead where one fits
 */
export function countBytesError(
  bytesTo: (end: number) => number,
  budget: number,
  window: { start: number; end: number },
  tile: number,
  names = BUDGET_NAMES,
  attaching = false,
) {
  const over = bytesOver(bytesTo, budget, window, tile)
  const windowBp = window.end - window.start
  const against = `against ${names.bytes} ${megabytes(budget)}`
  const own = attaching
    ? " with the tiles holding its departures' own steps"
    : ''
  const instead = attaching
    ? '. Count with a smaller attachBp, or leave this window uncounted'
    : ', so leave this window uncounted'
  return over === undefined
    ? undefined
    : new RegionTooLargeError(
        `Too many bytes to read to count this window: ${
          over.fitsBp === 0
            ? `${megabytes(over.least)} for its first tile alone${own}, ${against}${instead}`
            : `${megabytes(over.bytes)} ${against}. Count a window of about ${bpString(over.fitsBp)} or less, or leave this one uncounted`
        }`,
        budget,
        windowBp,
        Math.max(over.fitsBp, 1),
        over.fitsBp > 0
          ? { bytes: over.bytes, zoomHelps: true }
          : { bytes: over.bytes },
      )
}

// The widest span from the window's start, narrower than the window, whose
// read fits `budget`, 0 where none does; `bytesTo` never falls as its end
// grows
function narrowerFit(
  bytesTo: (end: number) => number,
  budget: number,
  window: { start: number; end: number },
) {
  let fits = window.start
  let over = window.end
  while (over - fits > 1) {
    const middle = Math.floor((fits + over) / 2)
    if (bytesTo(middle) <= budget) {
      fits = middle
    } else {
      over = middle
    }
  }
  return fits - window.start
}

// What a read past `budget` would fetch: `bytes` for the window, `least` for
// its first tile alone, and `fitsBp`, the span from the window's start whose
// read fits, 0 where none does; undefined for a read within the budget. The
// span grows a tile at a time, and where every tile fits and the window's
// read is still over, as the tiles fetched to attach can make it, it is the
// widest narrower one that fits, so it is never the window or more
function bytesOver(
  bytesTo: (end: number) => number,
  budget: number,
  window: { start: number; end: number },
  tile: number,
) {
  const bytes = bytesTo(window.end)
  let fitsBp = 0
  let least: number | undefined
  let more = bytes > budget
  for (
    let cs = Math.floor(window.start / tile) * tile;
    more && cs < window.end;
    cs += tile
  ) {
    const read = bytesTo(cs + 1)
    least ??= read
    more = read <= budget
    fitsBp = more ? cs + tile - window.start : fitsBp
  }
  const reaches = fitsBp >= window.end - window.start
  const narrower = reaches ? narrowerFit(bytesTo, budget, window) : fitsBp
  return bytes > budget
    ? {
        bytes,
        fitsBp: narrower,
        least:
          reaches && narrower === 0
            ? bytesTo(window.start + 1)
            : (least ?? bytes),
      }
    : undefined
}

export const megabytes = (bytes: number) => `${(bytes / 1e6).toFixed(1)} MB`
