import { panSNHaplotypeOrSample } from './util.ts'
import {
  NONE,
  RefPath,
  abridgedPiece,
  blockPositions,
  canonicalLink,
  checkedId,
  idOf,
  interiorAt,
  numbers,
  parseChoices,
  parseCopy,
  passCounts,
  refRows,
  tokens,
  traverse,
} from './walkChunk.ts'
import { WalkFormatError } from './walkFile.ts'

import type { AbridgedPiece, PieceBlocks, Route } from './walkChunk.ts'
import type { DepartureRecord } from './walkDeparture.ts'

export type CountReason = 'unloaded' | 'unattached' | 'revisited' | 'elsewhere'

/**
 * Why a count may fall short, one bit each; a count with none is exact.
 * `unloaded`: a haplotype over the window not counted there has a run off the
 * reference the call did not load. `unattached`: one has own sequence of a
 * departure over the window the call did not read. `revisited`: the reference
 * passes the node more than once, so a step on it can lie at another pass.
 * `elsewhere`: steps on it can lie in rows filed away from the window.
 */
export const COUNT_REASONS: Readonly<Record<CountReason, number>> =
  Object.freeze({
    unloaded: 1,
    unattached: 2,
    revisited: 4,
    elsewhere: 8,
  })

const REASONS: CountReason[] = [
  'unloaded',
  'unattached',
  'revisited',
  'elsewhere',
]

export const countReasons = (bits: number) =>
  REASONS.filter(reason => (bits & COUNT_REASONS[reason]) !== 0)

export interface CountColumns {
  haplotypes: Int32Array
  passes: Float64Array
  haplotypeReasons: Uint8Array
  passReasons: Uint8Array
}

/**
 * How many of a walk file's haplotypes step each node and link over a
 * window, and how often. `nodes.ids` ascend; links are handle pairs, each
 * `2 * id + reversed`, written the one way round a cut writes them, and
 * ascend by `from`, then `to`. Each count states in its reasons column why it
 * may fall short, 0 where it is exact
 */
export interface WalkCounts {
  haplotypes: {
    file: number
    window: number
    unloaded: string[]
    unattached: string[]
  }
  nodes: CountColumns & { ids: Int32Array }
  links: CountColumns & { from: Uint32Array; to: Uint32Array }
}

function popcount(word: number) {
  const pairs = word - ((word >>> 1) & 0x55555555)
  const nibbles = (pairs & 0x33333333) + ((pairs >>> 2) & 0x33333333)
  return Math.imul((nibbles + (nibbles >>> 4)) & 0x0f0f0f0f, 0x01010101) >>> 24
}

// A growing table of counted things, each with the haplotypes stepping it as
// a bitset and the steps on it
class Tally {
  readonly words: number
  size = 0
  bits: Uint32Array
  passes: Float64Array

  constructor(words: number) {
    this.words = words
    this.bits = new Uint32Array(256 * words)
    this.passes = new Float64Array(256)
  }

  add() {
    if (this.size === this.passes.length) {
      const bits = new Uint32Array(this.bits.length * 2)
      bits.set(this.bits)
      this.bits = bits
      const passes = new Float64Array(this.passes.length * 2)
      passes.set(this.passes)
      this.passes = passes
    }
    return this.size++
  }

  mark(k: number, haplotype: number) {
    const at = k * this.words + (haplotype >>> 5)
    this.bits[at] = this.bits[at]! | (1 << (haplotype & 31))
    this.passes[k] = this.passes[k]! + 1
  }

  count(k: number) {
    let n = 0
    for (let w = 0; w < this.words; w++) {
      n += popcount(this.bits[k * this.words + w]!)
    }
    return n
  }

  // whether a haplotype of `set` does not step k
  lacks(k: number, set: Uint32Array) {
    let lacking = false
    for (let w = 0; w < this.words; w++) {
      lacking ||= (set[w]! & ~this.bits[k * this.words + w]!) !== 0
    }
    return lacking
  }
}

const anyOf = (set: Uint32Array) => set.some(word => word !== 0)

// A route's nodes and links as counted, and its bp
interface RouteCounts {
  nodes: Int32Array
  links: Int32Array
  bp: number
}

export interface CountWalk {
  name: string
  fragStart: number
  piece: number
  hapOffset: number
  n: number
  interior: number
  span: { start: number; end: number }
  start: number
  startCopy: number
  choices: number[]
}

// A tile's R, its nodes and its links between neighbouring positions, which
// its detail units share
interface CountPath {
  path: RefPath
  nodes: Int32Array
  links: Int32Array
  lengths: Map<number, number>
}

// A tile's or detail unit's rows, ready to count
export interface CountUnit {
  start: number
  r: CountPath
  references: {
    name: string
    fragStart: number
    piece: number
    hapOffset: number
    from: number
    to: number
  }[]
  routes: Route[]
  counted: RouteCounts[]
  walks: CountWalk[]
  abridged: AbridgedPiece[]
}

// A piece's first and last handles, -1 where a run not loaded ends it
interface PieceEnds {
  haplotype: number
  first: number
  last: number
}

const lineKey = (name: string, fragStart: number) => `${name}\t${fragStart}`

/**
 * Counts, as it reads rows, the haplotypes stepping each node and link and
 * the steps on each, keeping no step: per node and link a bitset of
 * haplotypes and a pass count, and per tile its R's nodes and links once
 */
export class WalkTally {
  private readonly haplotypeIndex: Map<string, number>
  private readonly names: string[]
  private readonly nodes: Tally
  private readonly links: Tally
  private readonly nodeAt = new Map<number, number>()
  private readonly nodeIds: number[] = []
  // how many times the reference passes each node by a window tile's ref
  // row, 0 where none holds it, and whether such a pass lies on the drawn
  // stretch
  private readonly refPasses: number[] = []
  private readonly drawn: boolean[] = []
  // each link under the end it leaves, as an end's node index twice plus its
  // strand: the ends it reaches and its index, in turn
  private readonly linkAt: number[][] = []
  private readonly linkEnds: number[] = []
  private readonly window: Uint32Array
  private readonly unloaded: Uint32Array
  private readonly unattached: Uint32Array
  private readonly pieces = new Map<string, Map<number, PieceEnds>>()
  // the stretches of each contig the counted pieces cover, and whether the
  // call loaded every step of each
  private readonly stretches = new Map<
    string,
    { start: number; end: number; loaded: boolean }[]
  >()

  constructor(haplotypes: string[]) {
    this.names = haplotypes
    this.haplotypeIndex = new Map(haplotypes.map((name, k) => [name, k]))
    const words = Math.max(1, Math.ceil(haplotypes.length / 32))
    this.nodes = new Tally(words)
    this.links = new Tally(words)
    this.window = new Uint32Array(words)
    this.unloaded = new Uint32Array(words)
    this.unattached = new Uint32Array(words)
  }

  haplotype(name: string) {
    const listed = panSNHaplotypeOrSample(name)
    const k = this.haplotypeIndex.get(listed)
    if (k === undefined) {
      throw new WalkFormatError(
        `The walk file has rows of ${name}, and its header lists no haplotype ${listed}`,
      )
    }
    return k
  }

  private node(id: number) {
    let k = this.nodeAt.get(id)
    if (k === undefined) {
      k = this.nodes.add()
      this.nodeAt.set(id, k)
      this.nodeIds.push(id)
      this.refPasses.push(0)
      this.drawn.push(false)
    }
    return k
  }

  // The link from handle a to handle b, which is the link from b flipped to
  // a flipped, held under whichever way round leaves the smaller end
  private link(
    a: number,
    b: number,
    nodeA = this.node(idOf(a)),
    nodeB = this.node(idOf(b)),
  ) {
    const endA = 2 * nodeA + (a & 1)
    const endB = 2 * nodeB + (b & 1)
    const forward =
      endA < (endB ^ 1) || (endA === (endB ^ 1) && endB <= (endA ^ 1))
    const from = forward ? endA : endB ^ 1
    const to = forward ? endB : endA ^ 1
    const held = this.linkAt[from] ?? []
    this.linkAt[from] = held
    let k = -1
    for (let i = 0; i < held.length && k < 0; i += 2) {
      k = held[i] === to ? held[i + 1]! : -1
    }
    if (k < 0) {
      k = this.links.add()
      held.push(to, k)
      const [first, second] = canonicalLink(a, b)
      this.linkEnds.push(first, second)
    }
    return k
  }

  private pathLink(r: CountPath, i: number) {
    let k = r.links[i]!
    if (k < 0) {
      k = this.link(
        r.path.handles[i]!,
        r.path.handles[i + 1]!,
        r.nodes[i],
        r.nodes[i + 1],
      )
      r.links[i] = k
    }
    return k
  }

  private routesOf(rows: string[][], lengths: Map<number, number>) {
    const previous = { id: 0 }
    const routes: Route[] = []
    const counted: RouteCounts[] = []
    for (const cols of rows) {
      if (cols[3] === 'route') {
        const entry = tokens(cols[4]!, previous, [])[0] ?? NONE
        const interior = tokens(cols[5]!, previous, [])
        const exit = tokens(cols[6]!, previous, [])[0] ?? NONE
        const n = interior.length
        const nodes = new Int32Array(n)
        const links = new Int32Array(
          n - 1 + (entry === NONE ? 0 : 1) + (exit === NONE ? 0 : 1),
        )
        let at = 0
        let bp = 0
        let before = entry
        let beforeNode = entry === NONE ? -1 : this.node(idOf(entry))
        for (let k = 0; k < n; k++) {
          const h = interior[k]!
          const length = lengths.get(idOf(h))
          if (length === undefined) {
            throw new WalkFormatError(
              `a route steps on node ${idOf(h)}, which its unit's nodes row lacks`,
            )
          }
          bp += length
          const node = this.node(idOf(h))
          nodes[k] = node
          if (before !== NONE) {
            links[at++] = this.link(before, h, beforeNode, node)
          }
          before = h
          beforeNode = node
        }
        if (before !== NONE && exit !== NONE) {
          links[at++] = this.link(before, exit, beforeNode)
        }
        routes.push({ entry, interior, exit })
        counted.push({ nodes, links, bp })
      }
    }
    return { routes, counted }
  }

  private walksOf(rows: string[][], routes: Route[]) {
    const walks: CountWalk[] = []
    for (const cols of rows) {
      if (cols[3] === 'walk') {
        const choices = parseChoices(cols[10]!)
        let interior = 0
        for (let k = 0; k < choices.length; k += 4) {
          interior += routes[choices[k]!]?.interior.length ?? 0
        }
        const [start, startCopy] = parseCopy(cols[9]!)
        walks.push({
          name: cols[4]!,
          fragStart: +cols[5]!,
          hapOffset: +cols[6]!,
          piece: +cols[7]!,
          n: +cols[8]!,
          interior,
          span: { start: +cols[1]!, end: +cols[2]! },
          start,
          startCopy,
          choices,
        })
      }
    }
    return walks
  }

  private lengthsOf(rows: string[][], into: Map<number, number>) {
    for (const cols of rows) {
      if (cols[3] === 'nodes') {
        const deltas = numbers(cols[4]!)
        const lengths = numbers(cols[5]!)
        let id = 0
        for (let k = 0; k < deltas.length; k++) {
          id = checkedId(id + deltas[k]!)
          into.set(id, lengths[k]!)
        }
      }
    }
    return into
  }

  /**
   * A tile's rows, split on tabs, ready to count. A tile of the window
   * states how many times the reference passes each node of its R, which
   * makes a count on that node exact where its pass lies on the drawn
   * stretch; a tile fetched to attach states nothing of the kind
   */
  tile(rows: string[][], tileStart: number, ofWindow: boolean): CountUnit {
    const path = new RefPath()
    const lengths = new Map<number, number>()
    const references: CountUnit['references'] = []
    for (const cols of refRows(rows)) {
      const handles = tokens(cols[9]!, { id: 0 }, [])
      const nodeLengths = numbers(cols[10]!)
      const passes = passCounts(cols, handles.length)
      const from = path.handles.length
      path.add(handles, nodeLengths, +cols[6]!)
      for (let k = 0; k < handles.length; k++) {
        const id = idOf(handles[k]!)
        lengths.set(id, nodeLengths[k]!)
        const node = this.node(id)
        if (ofWindow) {
          this.refPasses[node] = Math.max(this.refPasses[node]!, passes[k]!)
        }
      }
      references.push({
        name: cols[4]!,
        fragStart: +cols[5]!,
        piece: +cols[7]!,
        hapOffset: +cols[6]!,
        from,
        to: path.handles.length,
      })
    }
    const r = {
      path,
      nodes: Int32Array.from(path.handles, h => this.node(idOf(h))),
      links: new Int32Array(path.handles.length).fill(-1),
      lengths: this.lengthsOf(rows, lengths),
    }
    const { routes, counted } = this.routesOf(rows, r.lengths)
    return {
      start: tileStart,
      r,
      references,
      routes,
      counted,
      walks: this.walksOf(rows, routes),
      abridged: rows
        .filter(cols => cols[3] === 'abridged')
        .map(cols => abridgedPiece(cols, tileStart)),
    }
  }

  // A detail unit's rows, decoded against its tile's R
  detail(rows: string[][], tile: CountUnit): CountUnit {
    const lengths = this.lengthsOf(rows, new Map(tile.r.lengths))
    const { routes, counted } = this.routesOf(rows, lengths)
    return {
      start: tile.start,
      r: tile.r,
      references: [],
      routes,
      counted,
      walks: this.walksOf(rows, routes),
      abridged: [],
    }
  }

  // Marks the nodes of a window tile's R whose pass lies on `drawn`
  placeReference(unit: CountUnit, drawn: { start: number; end: number }) {
    const { path, nodes } = unit.r
    for (let i = 0; i < path.handles.length; i++) {
      const node = nodes[i]!
      this.drawn[node] ||=
        path.starts[i]! < drawn.end && path.end(i) > drawn.start
    }
  }

  private piece(
    name: string,
    fragStart: number,
    piece: number,
    ends: PieceEnds,
  ) {
    const key = lineKey(name, fragStart)
    const line = this.pieces.get(key) ?? new Map<number, PieceEnds>()
    this.pieces.set(key, line)
    line.set(piece, ends)
    this.mark(this.window, ends.haplotype)
  }

  private stretch(name: string, start: number, end: number, loaded: boolean) {
    const held = this.stretches.get(name) ?? []
    this.stretches.set(name, held)
    held.push({ start, end, loaded })
  }

  private mark(set: Uint32Array, haplotype: number) {
    set[haplotype >>> 5] = set[haplotype >>> 5]! | (1 << (haplotype & 31))
  }

  // The reference's own piece in a window tile, counted on its haplotype
  countReference(
    unit: CountUnit,
    { name, fragStart, piece, from, to }: CountUnit['references'][number],
  ) {
    const haplotype = this.haplotype(name)
    const { r } = unit
    for (let i = from; i < to; i++) {
      this.nodes.mark(r.nodes[i]!, haplotype)
      if (i + 1 < to) {
        this.links.mark(this.pathLink(r, i), haplotype)
      }
    }
    this.piece(name, fragStart, piece, {
      haplotype,
      first: to > from ? r.path.handles[from]! : -1,
      last: to > from ? r.path.handles[to - 1]! : -1,
    })
  }

  /**
   * A walk row's steps as the row decodes them, counted on its haplotype
   * when `counting`, and its bp on its contig either way
   */
  walk(unit: CountUnit, walk: CountWalk, counting: boolean) {
    const haplotype = this.haplotype(walk.name)
    const { r, routes, counted } = unit
    const { path } = r
    let first = -1
    let last = -1
    let bp = 0
    let emitted = 0
    traverse(path, routes, walk.start, walk.startCopy, walk.choices, walk.n, {
      step: (handle, i) => {
        first = first < 0 ? handle : first
        last = handle
        bp += path.before[i + 1]! - path.before[i]!
        emitted++
        if (counting) {
          this.nodes.mark(r.nodes[i]!, haplotype)
        }
      },
      stretch: (from, to, forward) => {
        last = path.handleAt(to, forward)
        bp += forward
          ? path.before[to + 1]! - path.before[from + 1]!
          : path.before[from]! - path.before[to]!
        emitted += Math.abs(to - from)
        const low = Math.min(from, to)
        const high = counting ? Math.max(from, to) : low - 1
        for (let j = low; j <= high; j++) {
          if (j !== from) {
            this.nodes.mark(r.nodes[j]!, haplotype)
          }
          if (j < high) {
            this.links.mark(this.pathLink(r, j), haplotype)
          }
        }
      },
      route: (index, flip) => {
        const route = routes[index]!
        const { nodes, links, bp: routeBp } = counted[index]!
        const n = route.interior.length
        first = first < 0 && n > 0 ? interiorAt(route, flip, 0) : first
        last = n > 0 ? interiorAt(route, flip, n - 1) : last
        bp += routeBp
        emitted += n
        if (counting) {
          for (const node of nodes) {
            this.nodes.mark(node, haplotype)
          }
          for (const link of links) {
            this.links.mark(link, haplotype)
          }
        }
        return true
      },
    })
    if (emitted !== walk.n) {
      throw new WalkFormatError(
        `a walk row of ${walk.name} decodes to ${emitted} steps, not its ${walk.n}`,
      )
    }
    if (counting) {
      this.piece(walk.name, walk.fragStart, walk.piece, {
        haplotype,
        first,
        last,
      })
      this.stretch(walk.name, walk.hapOffset, walk.hapOffset + bp, true)
    }
    return bp
  }

  /**
   * An abridged piece the call does not load, from its blocks row: its steps
   * on R exactly, counted on its haplotype, and its runs off R as bp, which
   * leave the haplotype's counts off R lower bounds
   */
  blocks(unit: CountUnit, piece: AbridgedPiece, blocks: PieceBlocks) {
    const haplotype = this.haplotype(piece.name)
    const { r } = unit
    const { path } = r
    const { from, to } = blockPositions(piece, blocks, path)
    const count = from.length
    const firstOf = (k: number) =>
      blocks.against[k]!
        ? path.handleAt(to[k]! - 1, false)
        : path.handleAt(from[k]!, true)
    const lastOf = (k: number) =>
      blocks.against[k]!
        ? path.handleAt(from[k]!, false)
        : path.handleAt(to[k]! - 1, true)
    let offR = blocks.before + blocks.after
    for (let k = 0; k < count; k++) {
      for (let i = from[k]!; i < to[k]!; i++) {
        this.nodes.mark(r.nodes[i]!, haplotype)
        if (i + 1 < to[k]!) {
          this.links.mark(this.pathLink(r, i), haplotype)
        }
      }
      if (k > 0 && blocks.off[k] === 0) {
        this.links.mark(this.link(lastOf(k - 1), firstOf(k)), haplotype)
      }
      offR += blocks.off[k]!
    }
    if (offR > 0) {
      this.mark(this.unloaded, haplotype)
    }
    this.piece(piece.name, piece.fragStart, piece.piece, {
      haplotype,
      first: blocks.before > 0 ? -1 : firstOf(0),
      last: blocks.after > 0 ? -1 : lastOf(count - 1),
    })
    this.stretch(
      piece.name,
      piece.hapOffset,
      piece.hapOffset + piece.hapBp,
      offR === 0,
    )
  }

  // An abridged piece the call holds nothing of but its bp, which leaves its
  // haplotype's counts off R lower bounds
  unread(piece: AbridgedPiece) {
    const haplotype = this.haplotype(piece.name)
    this.mark(this.unloaded, haplotype)
    this.piece(piece.name, piece.fragStart, piece.piece, {
      haplotype,
      first: -1,
      last: -1,
    })
    this.stretch(
      piece.name,
      piece.hapOffset,
      piece.hapOffset + piece.hapBp,
      false,
    )
  }

  /**
   * Every departure record over the window: its haplotype is over the
   * window, and one whose own sequence, less what lies between W lines, the
   * counted pieces do not cover leaves that haplotype's counts off R lower
   * bounds
   */
  departures(records: readonly Readonly<DepartureRecord>[]) {
    for (const d of records) {
      const haplotype = this.haplotype(d.name)
      this.mark(this.window, haplotype)
      const need = d.ownBp - d.missingBp
      let loaded = 0
      let held = 0
      for (const { start, end, loaded: whole } of this.stretches.get(d.name) ??
        []) {
        const overlap = Math.max(
          0,
          Math.min(end, d.hapStart + d.ownBp) - Math.max(start, d.hapStart),
        )
        loaded += whole ? overlap : 0
        held += overlap
      }
      if (need > 0 && loaded < need && held < need) {
        this.mark(this.unattached, haplotype)
      }
    }
  }

  // The links between consecutive pieces of a line, where both are counted
  // and the call loaded the steps either side
  joinPieces() {
    for (const line of this.pieces.values()) {
      for (const [piece, ends] of line) {
        const next = line.get(piece + 1)
        if (next && ends.last >= 0 && next.first >= 0) {
          this.links.mark(this.link(ends.last, next.first), ends.haplotype)
        }
      }
    }
  }

  // Why the count of node k, or of a link with an end there, may fall short
  // where one of `missing` does not step it: 0 on a node of a window tile's
  // R passed once on the drawn stretch
  private nodeReasons(k: number, missing: (set: Uint32Array) => boolean) {
    const passes = this.refPasses[k]!
    return passes === 1
      ? this.drawn[k]!
        ? 0
        : COUNT_REASONS.elsewhere
      : COUNT_REASONS.elsewhere |
          (passes > 1 ? COUNT_REASONS.revisited : 0) |
          (missing(this.unloaded) ? COUNT_REASONS.unloaded : 0) |
          (missing(this.unattached) ? COUNT_REASONS.unattached : 0)
  }

  private columns(
    tally: Tally,
    order: number[],
    reasons: (k: number, missing: (set: Uint32Array) => boolean) => number,
  ): CountColumns {
    const file = this.names.length
    const held = new Map([
      [this.unloaded, anyOf(this.unloaded)],
      [this.unattached, anyOf(this.unattached)],
    ])
    const haplotypes = Int32Array.from(order, k => tally.count(k))
    return {
      haplotypes,
      passes: Float64Array.from(order, k => tally.passes[k]!),
      haplotypeReasons: Uint8Array.from(order, (k, i) =>
        haplotypes[i] === file
          ? 0
          : reasons(k, set => held.get(set) === true && tally.lacks(k, set)),
      ),
      passReasons: Uint8Array.from(order, k =>
        reasons(k, set => held.get(set) === true),
      ),
    }
  }

  private namesIn(set: Uint32Array) {
    return this.names.filter((_, k) => (set[k >>> 5]! & (1 << (k & 31))) !== 0)
  }

  result(): WalkCounts {
    this.joinPieces()
    const nodeOrder: number[] = []
    for (let k = 0; k < this.nodes.size; k++) {
      if (this.nodes.passes[k]! > 0) {
        nodeOrder.push(k)
      }
    }
    nodeOrder.sort((a, b) => this.nodeIds[a]! - this.nodeIds[b]!)
    const linkOrder: number[] = []
    for (let k = 0; k < this.links.size; k++) {
      if (this.links.passes[k]! > 0) {
        linkOrder.push(k)
      }
    }
    const fromOf = (k: number) => this.linkEnds[2 * k]!
    const toOf = (k: number) => this.linkEnds[2 * k + 1]!
    linkOrder.sort((a, b) => fromOf(a) - fromOf(b) || toOf(a) - toOf(b))
    const endNode = (handle: number) => this.nodeAt.get(idOf(handle))!
    const linkReasons = (k: number, missing: (set: Uint32Array) => boolean) =>
      this.nodeReasons(endNode(fromOf(k)), missing) |
      this.nodeReasons(endNode(toOf(k)), missing)
    return {
      haplotypes: {
        file: this.names.length,
        window: this.namesIn(this.window).length,
        unloaded: this.namesIn(this.unloaded),
        unattached: this.namesIn(this.unattached),
      },
      nodes: {
        ids: Int32Array.from(nodeOrder, k => this.nodeIds[k]!),
        ...this.columns(this.nodes, nodeOrder, (k, missing) =>
          this.nodeReasons(k, missing),
        ),
      },
      links: {
        from: Uint32Array.from(linkOrder, fromOf),
        to: Uint32Array.from(linkOrder, toOf),
        ...this.columns(this.links, linkOrder, linkReasons),
      },
    }
  }
}
