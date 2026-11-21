import { RefPath, idOf } from './walkChunk.ts'
import { WalkGraph, walkCut } from './walkRows.ts'

import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkNode, WalkRun } from './walkRows.ts'

const otherStrand = (handle: number) => (handle & 1 ? handle - 1 : handle + 1)

// A deterministic stream of numbers in [0, 1)
function random(seed: number) {
  let s = seed >>> 0
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0
    return s / 2 ** 32
  }
}

// RefPath.reach as it was written before its binary search: every occurrence
// of the target's node, the nearest on the way R is read
function reachByScan(path: RefPath, i: number, h: number, target: number) {
  const forward = path.handles[i] === h
  let found: number | undefined
  for (const j of path.occurrences.get(idOf(target)) ?? []) {
    if (forward && j >= i && path.handles[j] === target) {
      found = found === undefined || j < found ? j : found
    } else if (!forward && j <= i && otherStrand(path.handles[j]!) === target) {
      found = found === undefined || j > found ? j : found
    }
  }
  return found
}

const pathOf = (handles: number[]) => {
  const path = new RefPath()
  path.add(
    handles,
    handles.map(() => 1),
    0,
  )
  return path
}

describe('RefPath.reach', () => {
  // node 1 four times, twice each way, node 2 twice forward, node 3 once
  const path = pathOf([2, 4, 3, 6, 2, 4, 3, 2])

  test.each([
    { at: 0, target: 2, want: 0, why: 'the position itself' },
    { at: 1, target: 2, want: 4, why: 'the next forward copy past one' },
    { at: 5, target: 2, want: 7, why: 'the last copy' },
    { at: 0, target: 6, want: 3, why: 'a node R holds once' },
    { at: 4, target: 6, want: undefined, why: 'none left ahead' },
    { at: 0, target: 8, want: undefined, why: 'a node R lacks' },
    { at: 0, target: 5, want: undefined, why: 'a strand R never holds' },
  ])('read forward from $at, $target: $why', ({ at, target, want }) => {
    expect(path.reach(at, path.handles[at]!, target)).toBe(want)
    expect(reachByScan(path, at, path.handles[at]!, target)).toBe(want)
  })

  test.each([
    { at: 7, target: 3, want: 7, why: 'the position itself, flipped' },
    { at: 6, target: 3, want: 4, why: 'the copy before a reversed one' },
    { at: 3, target: 3, want: 0, why: 'the first copy' },
    { at: 7, target: 7, want: 3, why: 'a node R holds once' },
    { at: 2, target: 7, want: undefined, why: 'none left behind' },
    { at: 7, target: 9, want: undefined, why: 'a node R lacks' },
  ])('read backward from $at, $target: $why', ({ at, target, want }) => {
    const h = otherStrand(path.handles[at]!)
    expect(path.reach(at, h, target)).toBe(want)
    expect(reachByScan(path, at, h, target)).toBe(want)
  })

  test('finds what a scan of every occurrence finds, on paths that repeat a few nodes hundreds of times either way', () => {
    const next = random(7)
    const differ: string[] = []
    for (let p = 0; p < 20; p++) {
      const handles = Array.from(
        { length: 600 },
        () => 2 * Math.floor(next() * (2 + p)) + (next() < 0.3 ? 1 : 0),
      )
      const repeating = pathOf(handles)
      for (let q = 0; q < 300; q++) {
        const i = Math.floor(next() * handles.length)
        const h = next() < 0.5 ? handles[i]! : otherStrand(handles[i]!)
        const target = 2 * Math.floor(next() * (3 + p)) + (next() < 0.5 ? 1 : 0)
        const got = repeating.reach(i, h, target)
        const want = reachByScan(repeating, i, h, target)
        if (got !== want) {
          differ.push(`path ${p} from ${i} ${h} to ${target}: ${got} ${want}`)
        }
      }
    }
    expect(differ).toEqual([])
  })
})

const window = { refName: 'R', start: 0, end: 100, context: 0 }

// The L lines a cut's tables list, as `source\tstrand\ttarget\tstrand`, each
// node named by its id, as the tables order them
function tableLines(links: [number, number][], kept: number[]) {
  const graph = new WalkGraph(window)
  for (const [a, b] of links) {
    graph.addLinkEnds(a, b)
  }
  const { nodes, links: table } = graph.tables(new Set(kept), [], [])
  return Array.from(table.from, (from, k) => {
    const strands = table.strands[k]!
    return `${nodes.names[from]}\t${strands & 1 ? '-' : '+'}\t${nodes.names[table.to[k]!]}\t${strands & 2 ? '-' : '+'}`
  })
}

// The same links each once as their L lines, sorted as strings, as the tables
// sorted them before they keyed links by number
const lineOf = ([a, b]: [number, number]) =>
  `${Math.floor(a / 2)}\t${a & 1 ? '-' : '+'}\t${Math.floor(b / 2)}\t${b & 1 ? '-' : '+'}`
const sortedLines = (links: [number, number][], kept: number[]) =>
  [
    ...new Set(
      links
        .filter(
          ([a, b]) =>
            kept.includes(Math.floor(a / 2)) &&
            kept.includes(Math.floor(b / 2)),
        )
        .map(lineOf),
    ),
  ].sort()

describe('the tables order links as their L lines sort as strings', () => {
  const both = (ids: number[]) =>
    ids.flatMap(a =>
      ids.flatMap(b =>
        [0, 1, 2, 3].map((s): [number, number] => [
          2 * a + (s & 1),
          2 * b + (s >> 1),
        ]),
      ),
    )

  test.each([
    {
      ids: [9, 10, 1, 100, 99],
      why: 'a shorter id after a longer one it does not begin',
    },
    { ids: [12, 123, 1234, 124, 13], why: 'an id before every id it begins' },
    { ids: [0, 1, 10, -1, -10, -2], why: 'negatives first, then 0' },
    { ids: [-12, -123, -9, -90], why: 'negatives by their digits' },
    {
      ids: [2147483647, 214748364, 2147483648, 999999999, 1000000000],
      why: '10-digit ids',
    },
    {
      ids: [999999999999, 99999999999, 100000000000, -999999999999, -1],
      why: 'ids just below 1e12 either way',
    },
  ])('$why: $ids', ({ ids }) => {
    const links = both(ids)
    expect(tableLines(links, ids)).toEqual(sortedLines(links, ids))
  })

  test('the strands of a link order it after its ids, + before -, source before target', () => {
    const links: [number, number][] = [
      [3, 4],
      [2, 5],
      [3, 5],
      [2, 4],
      [2, 3],
    ]
    expect(tableLines(links, [1, 2])).toEqual([
      '1\t+\t1\t-',
      '1\t+\t2\t+',
      '1\t+\t2\t-',
      '1\t-\t2\t+',
      '1\t-\t2\t-',
    ])
  })

  test('a link between nodes the cut does not keep is left out', () => {
    expect(
      tableLines(
        [
          [2, 4],
          [2, 6],
        ],
        [1, 2],
      ),
    ).toEqual(['1\t+\t2\t+'])
  })

  test('each link once, of tens of thousands with repeats and handles alike in their low 32 bits, as the strings sort', () => {
    const next = random(11)
    const ids = [
      ...Array.from({ length: 300 }, (_, k) => k - 150),
      // handles of -5 and 2147483643 share their low 32 bits, as do those of
      // -6 and 2147483642, so they hash alike
      2147483643,
      2147483642,
      -2147483647,
      2147483647,
    ]
    const links: [number, number][] = Array.from({ length: 60000 }, () => {
      const a = ids[Math.floor(next() * ids.length)]!
      const b = ids[Math.floor(next() * ids.length)]!
      return [2 * a + (next() < 0.5 ? 1 : 0), 2 * b + (next() < 0.5 ? 1 : 0)]
    })
    const kept = ids.filter(id => id % 7 !== 0)
    const lines = tableLines([...links, ...links.slice(0, 5000)], kept)
    expect(lines.length).toBeGreaterThan(20000)
    expect(lines).toEqual(sortedLines(links, kept))
  })
})

// A run of `ids` all along R, one bp each from `start`, every step on the
// drawn stretch
function onR(name: string, ids: number[], start = 0): WalkRun {
  return {
    name,
    fragStart: 0,
    hapOffset: 0,
    firstPiece: 0,
    pieceStarts: [],
    pieceChunks: [0],
    ids: Int32Array.from(ids),
    rev: new Uint8Array(ids.length),
    refStarts: Float64Array.from(ids, (_, k) => start + k),
    against: new Uint8Array(ids.length),
    bp: Int32Array.from(ids, () => 1),
  }
}

describe('walkCut lists each link its steps walk once, in the order first walked', () => {
  test('of runs that walk thousands of links again and again, among handles alike in their low 32 bits', () => {
    const next = random(3)
    const pool = [
      1,
      2,
      3,
      -5,
      2147483643,
      -6,
      2147483642,
      ...Array.from({ length: 200 }, (_, k) => 10 + k),
    ]
    const runs = Array.from({ length: 30 }, (_, r) =>
      onR(
        `H${r}#1#c`,
        Array.from(
          { length: 2000 },
          () => pool[Math.floor(next() * pool.length)]!,
        ),
      ),
    )
    const nodes = new Map<number, WalkNode>()
    const { links } = walkCut(runs, nodes, 'R', [], {
      drawn: { start: 0, end: 5000 },
    })
    const seen = new Set<string>()
    const want: number[] = []
    for (const { ids } of runs) {
      for (let k = 1; k < ids.length; k++) {
        const key = `${2 * ids[k - 1]!} ${2 * ids[k]!}`
        if (!seen.has(key)) {
          seen.add(key)
          want.push(2 * ids[k - 1]!, 2 * ids[k]!)
        }
      }
    }
    expect(want.length / 2).toBeGreaterThan(30000)
    expect(links).toEqual(want)
  })
})

describe('walkCut sorts departures by name as strings, then contig position, then leave position', () => {
  const record = (
    name: string,
    hapStart: number,
    leavePos: number,
    rejoinPos: number,
  ): WalkDeparture => ({
    name,
    hapStart,
    ownBp: 5,
    missingBp: 0,
    leavePos,
    leaveReversed: false,
    leavePasses: 1,
    rejoinPos,
    rejoinReversed: false,
    rejoinPasses: 1,
    fromSteps: false,
  })
  // the order the comparator it replaced gave, a stable sort on name strings
  const byStrings = (records: WalkDeparture[]) =>
    records.toSorted(
      (a, b) =>
        (a.name < b.name ? -1 : a.name > b.name ? 1 : 0) ||
        a.hapStart - b.hapStart ||
        a.leavePos - b.leavePos,
    )
  const facts = (departures: WalkDeparture[]) =>
    departures.map(d => `${d.name} ${d.hapStart} ${d.leavePos} ${d.rejoinPos}`)
  const cut = (records: WalkDeparture[]) =>
    walkCut([], new Map(), 'R', records, {
      drawn: { start: 1000, end: 1001 },
      attachBp: 0,
    }).departures

  test('names whose numbers sort otherwise than their strings: H10 before H2, H1#10 before H1#9', () => {
    const records = [
      record('H2#1#c', 5, 50, 60),
      record('H10#1#c', 5, 50, 60),
      record('H1#9#c', 5, 50, 60),
      record('H1#10#c', 5, 50, 60),
      record('H9#1#c', 1, 50, 60),
      record('h1#1#c', 1, 50, 60),
      record('H#1#c', 1, 50, 60),
    ]
    expect(facts(cut(records))).toEqual([
      'H#1#c 1 50 60',
      'H1#10#c 5 50 60',
      'H1#9#c 5 50 60',
      'H10#1#c 5 50 60',
      'H2#1#c 5 50 60',
      'H9#1#c 1 50 60',
      'h1#1#c 1 50 60',
    ])
  })

  test('ties on all three keep the order given, among thousands in runs of a name', () => {
    const next = random(5)
    const names = ['H10#1#c', 'H2#1#c', 'H1#1#c', 'H1#2#c']
    const records = Array.from({ length: 4000 }, (_, k) =>
      record(
        names[Math.floor(k / 500) % names.length]!,
        Math.floor(next() * 20),
        Math.floor(next() * 3),
        k,
      ),
    )
    expect(facts(cut(records))).toEqual(facts(byStrings(records)))
  })
})
