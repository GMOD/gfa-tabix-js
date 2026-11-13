import { LocalFile } from 'generic-filehandle2'

import { walkFileOf } from './testing.ts'
import { RegionTooLargeError } from './util.ts'
import { RefPath, WalkRowError, decodeChunk, ownRunsRow } from './walkChunk.ts'
import { formsCut, walkDepartures } from './walkDepartures.ts'
import { WalkFormatError } from './walkFile.ts'
import { WalkReader } from './walkReader.ts'
import { isPlaceholder } from './walkRows.ts'

import type { WalkCutOptions } from './walkReader.ts'
import type { WalkFragment } from './walkRows.ts'

// What a cut states of an abridged piece it did not load, from the piece's
// `blocks` row (gfa-to-tabix's ADR 0012), against what the cut that loads it
// states. The cases are those of that repo's
// docs/handoff/2026-11-11-abridged-pieces.md, each on a fixture every walk of
// which the writer abridged.

const DATA = new URL('../test_data/', import.meta.url).pathname

type Cut = Awaited<ReturnType<WalkReader['cut']>>

// A reader of a fixture and the kinds of the units each of its reads fetched
function counted(name: string) {
  const path = `${DATA}${name}.walks.gz`
  const data = new LocalFile(path)
  const reader = WalkReader.fromFilehandles({
    data,
    directory: new LocalFile(`${path}.idx`),
  })
  const fetched: string[] = []
  const read = data.read.bind(data)
  const watch = async (refName: string) => {
    const units = await reader.units(refName)
    data.read = (length, position = 0, opts = {}) => {
      for (const u of units) {
        if (u.offset >= position && u.offset < position + length) {
          fetched.push(`${u.kind} ${u.name}`)
        }
      }
      return read(length, position, opts)
    }
  }
  return { reader, fetched, watch }
}

const reader = (name: string) => counted(name).reader

const all = () => true

// A fragment's steps: each on the reference as its strand, node and place,
// and each run off it as its bp
function steps(cut: Cut, f: WalkFragment) {
  const out: string[] = []
  let off = 0
  for (let i = 0; i <= f.ids.length; i++) {
    const offR = i < f.ids.length && Number.isNaN(f.refStarts[i])
    if (offR) {
      const node = cut.graph.nodes.get(f.ids[i]!)!
      off += node.end - node.start
    } else {
      if (off > 0) {
        out.push(`off:${off}`)
      }
      off = 0
      if (i < f.ids.length) {
        out.push(
          `${f.against[i] === 1 ? '<' : '>'}${f.ids[i]}@${f.refStarts[i]}`,
        )
      }
    }
  }
  return out.join(' ')
}

// Each fragment of a walk as its stretch of the contig and its steps
const rows = (cut: Cut, name: string) =>
  cut.fragments
    .filter(f => f.name === name)
    .map(f => [f.hapStart, f.hapEnd, steps(cut, f)])

const departures = (cut: Cut, name: string) =>
  cut.departures
    .filter(d => d.name === name)
    .map(d => ({
      hapStart: d.hapStart,
      ownBp: d.ownBp,
      missingBp: d.missingBp,
      leave: `${d.leavePos}${d.leaveReversed ? '-' : '+'}`,
      rejoin: `${d.rejoinPos}${d.rejoinReversed ? '-' : '+'}`,
      fromSteps: d.fromSteps,
    }))

const haplotypes = (cut: Cut, refName: string) =>
  walkDepartures(
    {
      fragments: cut.fragments,
      departures: cut.departures,
      nodes: cut.graph.nodes,
    },
    refName,
  ).sort((a, b) => (a.haplotype < b.haplotype ? -1 : 1))

const bpOf = (cut: Cut, ids: Int32Array) =>
  ids.reduce((sum, id) => {
    const node = cut.graph.nodes.get(id)!
    return sum + node.end - node.start
  }, 0)

// Each link by the node ends it joins, a node off the reference as `off`
function links(cut: Cut) {
  const { nodes, links: table } = cut.graph.tables(
    cut.kept,
    cut.fragments,
    cut.departures,
    cut.gaps,
  )
  const end = (index: number, right: boolean) =>
    nodes.ranks[index] === 0
      ? `${nodes.names[index]}${right ? ']' : '['}`
      : 'off'
  const out = new Set<string>()
  for (let k = 0; k < table.from.length; k++) {
    const strands = table.strands[k]!
    const a = end(table.from[k]!, (strands & 1) === 0)
    const b = end(table.to[k]!, (strands & 2) !== 0)
    if (a !== 'off' || b !== 'off') {
      out.add(a < b ? `${a} ${b}` : `${b} ${a}`)
    }
  }
  return [...out].sort()
}

// Everything a cut states that a cut loading the pieces must state too
const stated = (cut: Cut, refName: string) => ({
  fragments: cut.fragments
    .filter(f => f.name !== refName)
    .map(
      f =>
        `${f.name} ${f.hapStart}-${f.hapEnd} trimmed ${f.trimmedBefore} and ${f.trimmedAfter}: ${steps(cut, f)}`,
    )
    .sort(),
  departures: cut.departures.map(d => ({
    name: d.name,
    ...departures(cut, d.name).find(x => x.hapStart === d.hapStart),
    passes: [d.leavePasses, d.rejoinPasses],
    held: d.steps ? bpOf(cut, d.steps.ids) : 'none',
  })),
  gaps: cut.gaps.map(
    g =>
      `${g.name} edge ${g.edge} ${g.start}-${g.end} beside ${g.fragment?.hapStart}: ${g.segments.map(x => `${x.kind} ${x.bp} missing ${x.missingBp} ${x.before}..${x.after}`).join(', ')}`,
  ),
  kept: [...cut.kept]
    .filter(id => cut.graph.nodes.get(id)?.rank === 0)
    .sort((a, b) => a - b),
  links: links(cut),
  haplotypes: haplotypes(cut, refName),
})

// A step's name in the tables: its node's id, or a placeholder's walk and
// stretch of its contig
const nameOf = (cut: Cut, id: number) => {
  const node = cut.graph.nodes.get(id)
  return node?.placeholder
    ? `${node.refName}:${node.start}-${node.end}`
    : String(id)
}

const R = 'R#0#c'

// blocks.R.gfa: R over twelve 10 bp nodes from 0, in tiles of 60 bp, and
// every walk abridged
async function cuts(start: number, end: number, opts: WalkCutOptions = {}) {
  const walks = reader('blocks.R')
  const window = { start, end }
  return {
    unloaded: await walks.cut(R, window, { context: 0, ...opts }),
    loaded: await walks.cut(R, window, { context: 0, ...opts, expand: all }),
  }
}

// W HAIRPIN 1 c >7>8>9>10>21<10<9<8<7, 7 bp in 21: its first and last steps
// are both on 7, the first node of the tile at 60, and its abridged row
// spans 60-70
test.each([
  { start: 90, end: 100, row: [30, 57, '>10@90 off:7 <10@90'] },
  { start: 80, end: 90, row: [20, 67, '>9@80 >10@90 off:7 <10@90 <9@80'] },
  {
    start: 60,
    end: 70,
    row: [0, 87, '>7@60 >8@70 >9@80 >10@90 off:7 <10@90 <9@80 <8@70 <7@60'],
  },
])(
  'a hairpin that starts and ends on one node at its tile’s start is kept where it steps, both passes and the bp off R between them, cut $start-$end',
  async ({ start, end, row }) => {
    const { unloaded, loaded } = await cuts(start, end)
    expect(rows(unloaded, 'HAIRPIN#1#c')).toEqual([row])
    expect(stated(unloaded, R)).toEqual(stated(loaded, R))
  },
)

// W MIXED 1 c >8<11<10<9<8: its first step runs along R and the other four
// against it
test('a piece whose ends differ in strand is kept past its end steps, each step with its own strand', async () => {
  const { unloaded, loaded } = await cuts(100, 110)
  expect(rows(unloaded, 'MIXED#1#c')).toEqual([[10, 20, '<11@100']])
  const whole = await cuts(60, 120)
  expect(rows(whole.unloaded, 'MIXED#1#c')).toEqual([
    [0, 50, '>8@70 <11@100 <10@90 <9@80 <8@70'],
  ])
  const mixed = (cut: Cut) =>
    haplotypes(cut, R).find(h => h.haplotype === 'MIXED#1')
  expect(mixed(whole.unloaded)).toEqual(mixed(whole.loaded))
  expect(
    whole.unloaded.fragments
      .filter(f => f.name === 'MIXED#1#c')
      .map(f => [...f.against]),
  ).toEqual([[0, 1, 1, 1, 1]])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// W LOWFIRST 1 c >10>11>12>7>8>9: both end steps along R, at 90 and 80, and
// steps down to 60
test('a piece whose first step is not its lowest is kept below both end steps', async () => {
  const { unloaded, loaded } = await cuts(60, 70)
  expect(rows(unloaded, 'LOWFIRST#1#c')).toEqual([[30, 40, '>7@60']])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// W JUMP 1 c >7>8>22>12, 5 bp in 22: it has no step on 80-110
test('a stretch a piece jumps over holds its departure and no step of it', async () => {
  const { unloaded, loaded } = await cuts(90, 100)
  expect(rows(unloaded, 'JUMP#1#c')).toEqual([])
  expect(departures(unloaded, 'JUMP#1#c')).toEqual([
    {
      hapStart: 20,
      ownBp: 5,
      missingBp: 0,
      leave: '80+',
      rejoin: '110+',
      fromSteps: true,
    },
  ])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
  const whole = await cuts(60, 120)
  expect(rows(whole.unloaded, 'JUMP#1#c')).toEqual([
    [0, 35, '>7@60 >8@70 off:5 >12@110'],
  ])
})

// W PARTIAL 1 c >7>8>9>23>10>11>12, 50 bp in 23: 110 bp, 10 of them over
// 60-70
test('a window over part of a piece keeps that part, and counts the rest as trimmed', async () => {
  const { unloaded, loaded } = await cuts(60, 70)
  const [f] = unloaded.fragments.filter(x => x.name === 'PARTIAL#1#c')
  expect(rows(unloaded, 'PARTIAL#1#c')).toEqual([[0, 10, '>7@60']])
  expect([f!.trimmedBefore, f!.trimmedAfter]).toEqual([0, 100])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// W TWICE 1 c >8>9>10>24>8>9>10, 6 bp in 24: it passes 9 twice, with steps on
// 10 and 8 between
test('a contained cut splits an unloaded piece where it steps off the window, as it splits the loaded one', async () => {
  const { unloaded, loaded } = await cuts(80, 90, { contained: true })
  expect(rows(unloaded, 'TWICE#1#c')).toEqual([
    [10, 20, '>9@80'],
    [46, 56, '>9@80'],
  ])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// cover.R.gfa: R over twelve 1000 bp nodes in one tile. FAR walks 1 to 8,
// 300 bp of its own and back to 1, so its abridged row spans 0-1000; NEAR
// walks 5 to 9 and OFF 11 and 12
test('the forms cut covers a window with every haplotype stepping on it, whatever its end steps span', async () => {
  const region = { refName: 'c', start: 6000, end: 6500 }
  const forms = await formsCut(reader('cover.R'), region)
  const loaded = await reader('cover.R').cut(R, region, {
    links: false,
    expand: all,
  })
  expect(forms.filter(h => h.covered.length > 0)).toEqual(haplotypes(loaded, R))
  expect(forms.map(h => [h.haplotype, h.covered])).toEqual([
    ['FAR#1', [[5000, 8000]]],
    ['NEAR#1', [[5000, 8000]]],
    ['OFF#1', []],
  ])
})

// W MULTI 1 c >7>8>1>2>9>25>10>3>11>12, 3 bp in 25: pieces 0, 2 and 4 are in
// the tile at 60 and share its detail unit 60:1
test('each piece of a detail unit holding several is read from its own blocks row', async () => {
  const { unloaded, loaded } = await cuts(60, 120)
  expect(rows(unloaded, 'MULTI#1#c')).toEqual([
    [0, 20, '>7@60 >8@70'],
    [40, 63, '>9@80 off:3 >10@90'],
    [73, 93, '>11@100 >12@110'],
  ])
  expect(
    unloaded.fragments
      .filter(f => f.name === 'MULTI#1#c')
      .map(f => f.abridged!.map(a => [a.detail.name, a.steps])),
  ).toEqual([[['60:1', 2]], [['60:1', 3]], [['60:1', 2]]])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// W TAIL 1 c >5>6>26>7>8, 4 bp in 26, which trail its piece in the tile at
// 0; its record from 6 to 7 owns them
test('an unloaded piece joins the next piece read, the bp off R between them a placeholder its record holds', async () => {
  const { unloaded, loaded } = await cuts(50, 70)
  expect(rows(unloaded, 'TAIL#1#c')).toEqual([[10, 34, '>6@50 off:4 >7@60']])
  const [record] = unloaded.departures.filter(d => d.name === 'TAIL#1#c')
  expect(departures(unloaded, 'TAIL#1#c')).toEqual([
    {
      hapStart: 20,
      ownBp: 4,
      missingBp: 0,
      leave: '60+',
      rejoin: '60+',
      fromSteps: false,
    },
  ])
  expect([...record!.steps!.ids].map(isPlaceholder)).toEqual([true])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
  const gaps = (cut: Cut) =>
    cut.gaps
      .filter(g => g.name === 'TAIL#1#c')
      .map(g => [g.edge, g.start, g.end, g.segments.map(x => [x.kind, x.bp])])
  expect(gaps(unloaded)).toEqual([
    [1, 0, 10, [['read', 10]]],
    [2, 34, 44, [['read', 10]]],
  ])
  expect(gaps(unloaded)).toEqual(gaps(loaded))
})

// TAIL's piece in the tile at 0, >5>6>26, is not over 60-70: its record from
// 6 to 7 holds 26 as a placeholder, and no fragment holds any of it
test('a piece that reaches the cut only as placeholders among a record’s own steps is counted on the record, in the cut and in the tables', async () => {
  const { unloaded } = await cuts(60, 70)
  const [record] = unloaded.departures.filter(d => d.name === 'TAIL#1#c')
  expect([...record!.steps!.ids].map(isPlaceholder)).toEqual([true])
  expect(record!.abridged?.map(a => [a.detail.name, a.steps])).toEqual([
    ['0:0', 3],
  ])
  const counts = (cut: Cut) => {
    const tables = cut.graph.tables(
      cut.kept,
      cut.fragments,
      cut.departures,
      cut.gaps,
    )
    const of = (names: string[], count: Float64Array | undefined) =>
      names.flatMap((name, i) => (name === 'TAIL#1#c' ? [count?.[i]] : []))
    return {
      walks: of(tables.walks.names, tables.walks.abridgedSteps),
      departures: of(tables.departures.names, tables.departures.abridgedSteps),
    }
  }
  expect(counts(unloaded)).toEqual({ walks: [2], departures: [3] })
  expect(counts((await cuts(50, 70)).unloaded)).toEqual({
    walks: [5],
    departures: [0],
  })
})

// hairpin.R.gfa, the piece of gfa-to-tabix's walks-format.md: R over six
// 6 bp nodes, and W H 1 c >12>1>2<2<1>11>5>6>3>13>14 with 7 bp of its own
// before its first step on R, 10 bp between and 50 bp after its last
test('a piece that turns back, runs off R and jumps back is its list in path order, with its bp off R at either end', async () => {
  const walks = reader('hairpin.R')
  const window = { start: 0, end: 36 }
  const unloaded = await walks.cut(R, window, { context: 0 })
  const loaded = await walks.cut(R, window, { context: 0, expand: all })
  expect(rows(unloaded, 'H#1#c')).toEqual([
    [7, 59, '>1@0 >2@6 <2@6 <1@0 off:10 >5@24 >6@30 >3@12'],
  ])
  expect(
    unloaded.fragments
      .filter(f => f.name === 'H#1#c')
      .map(f => [f.trimmedBefore, f.trimmedAfter]),
  ).toEqual([[7, 50]])
  expect(links(unloaded)).toEqual(
    expect.arrayContaining(['2] 2]', '3[ 6]', '1[ off', '5[ off']),
  )
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// W TAIL 1 c >5>6>26>7>8: its piece in the tile at 0 ends with 4 bp of its
// own, and a window on 5 reads neither its next piece nor its record
test('a piece’s bp off R after its last step on R are read, the next piece unread', async () => {
  const { unloaded, loaded } = await cuts(40, 50)
  const [f] = unloaded.fragments.filter(x => x.name === 'TAIL#1#c')
  expect(rows(unloaded, 'TAIL#1#c')).toEqual([[0, 10, '>5@40']])
  expect([f!.trimmedBefore, f!.trimmedAfter]).toEqual([0, 14])
  expect(
    unloaded.gaps
      .filter(g => g.name === 'TAIL#1#c')
      .map(g => [g.edge, g.start, g.end, g.segments.map(x => [x.kind, x.bp])]),
  ).toEqual([[2, 10, 24, [['read', 14]]]])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// contigends.R.gfa: W L 1 c >11>12>13>3>4>5>6 starts with 30 bp of its own,
// and W M 1 c >1>2>3>11>12>13 ends with them, at their contigs' ends
test('a piece’s bp off R before its first step on R and after its last are read at a contig’s ends', async () => {
  const walks = reader('contigends.R')
  const window = { start: 12, end: 18 }
  const unloaded = await walks.cut(R, window, { context: 0 })
  const loaded = await walks.cut(R, window, { context: 0, expand: all })
  expect(
    unloaded.fragments
      .filter(f => f.name !== R)
      .map(f => [
        f.name,
        f.hapStart,
        f.hapEnd,
        f.trimmedBefore,
        f.trimmedAfter,
        steps(unloaded, f),
      ]),
  ).toEqual([
    ['L#1#c', 30, 36, 30, 6, '>3@12'],
    ['M#1#c', 12, 18, 0, 30, '>3@12'],
  ])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
  const whole = { start: 0, end: 36 }
  expect(stated(await walks.cut(R, whole), R)).toEqual(
    stated(await walks.cut(R, whole, { expand: all }), R),
  )
})

// ownends.R.gfa: W X 1 c >1>2>11>12>13>5>6 leaves R after 2 for 30 bp of its
// own, which trail its piece in the tile at 0; a window on 5 reads that tile
// for the record's own steps alone, and the blocks unit of its own tile at 20
test.each([
  { expand: undefined, held: ['placeholder'] },
  { expand: all, held: [11, 12, 13] },
])(
  'a record whose leave piece is abridged in a tile fetched for its own steps holds them, a placeholder from the abridged row and the record unless expanded, and that tile’s blocks unit is not read: $held',
  async ({ expand, held }) => {
    const { reader: walks, fetched, watch } = counted('ownends.R')
    await watch(R)
    const cut = await walks.cut(
      R,
      { start: 24, end: 30 },
      { context: 0, keep: name => name === 'X#1#c' || name === R, expand },
    )
    const [record] = cut.departures
    expect(departures(cut, 'X#1#c')).toEqual([
      {
        hapStart: 12,
        ownBp: 30,
        missingBp: 0,
        leave: '12+',
        rejoin: '24+',
        fromSteps: false,
      },
    ])
    expect(
      [...record!.steps!.ids].map(id =>
        isPlaceholder(id) ? 'placeholder' : id,
      ),
    ).toEqual(held)
    expect(bpOf(cut, record!.steps!.ids)).toBe(30)
    expect(fetched.filter(u => u.startsWith('blocks'))).toEqual(
      expand ? [] : ['blocks 20'],
    )
    expect(
      [...record!.steps!.ids].filter(isPlaceholder).map(id => nameOf(cut, id)),
    ).toEqual(expand ? [] : ['X#1#c:12-42'])
    expect(fetched.filter(u => u.startsWith('detail'))).toHaveLength(
      expand ? 2 : 0,
    )
  },
)

// The graphs of gfa-to-tabix's ADR 0013, each abridged whole, cut on a window
// where a record rejoins and whose tile lacks its leave piece, so its tile at
// 0 is fetched to attach its own steps there. In each pair of twins the tiles and the record differ in
// hapBp alone, and H or P has 10 or 5 bp, or 700 or 690, in its leave piece
// before a later W line or a step on another sequence. adv's CROSS, NOREF
// and PLAIN have rows alike but for path and hapBp. pickwrong's record owns
// the 7 bp that lead its contig's third W line, in a tile where its first W
// line's piece starts before the record does. twopieces' record owns 10 bp
// that trail one piece, in the tile fetched, and 7 that lead another, in the
// window's. The one blocks unit a cut reads is its window's tile's
test.each([
  {
    fixture: 'twinmidlineA.R',
    windowBlocks: ['blocks 24'],
    refName: 'R#0#chr1',
    window: { start: 24, end: 36 },
    held: { 'H#1#c': ['H#1#c:6-16'] },
    nodes: { 'H#1#c': ['11'] },
  },
  {
    fixture: 'twinmidlineB.R',
    windowBlocks: ['blocks 24'],
    refName: 'R#0#chr1',
    window: { start: 24, end: 36 },
    held: { 'H#1#c': ['H#1#c:6-11'] },
    nodes: { 'H#1#c': ['11'] },
  },
  {
    fixture: 'twinotherseqA.R',
    windowBlocks: ['blocks 24'],
    refName: 'R#0#chr1',
    window: { start: 24, end: 36 },
    held: { 'H#1#c': ['H#1#c:6-16'] },
    nodes: { 'H#1#c': ['11'] },
  },
  {
    fixture: 'twinotherseqB.R',
    windowBlocks: ['blocks 24'],
    refName: 'R#0#chr1',
    window: { start: 24, end: 36 },
    held: { 'H#1#c': ['H#1#c:6-11'] },
    nodes: { 'H#1#c': ['11'] },
  },
  {
    fixture: 'twinown700.R',
    windowBlocks: ['blocks 12'],
    refName: 'R#0#c',
    window: { start: 12, end: 24 },
    held: { 'P#1#c': ['P#1#c:12-712'] },
    nodes: { 'P#1#c': ['21'] },
  },
  {
    fixture: 'twinown690.R',
    windowBlocks: ['blocks 12'],
    refName: 'R#0#c',
    window: { start: 12, end: 24 },
    held: { 'P#1#c': ['P#1#c:12-702'] },
    nodes: { 'P#1#c': ['21'] },
  },
  {
    fixture: 'adv.R',
    windowBlocks: ['blocks 60'],
    refName: 'R#0#c1',
    window: { start: 90, end: 120 },
    held: {
      'CROSS#1#c': ['CROSS#1#c:30-37'],
      'NOREF#1#c': ['NOREF#1#c:30-37'],
      'PLAIN#1#c': ['PLAIN#1#c:30-36'],
    },
    nodes: { 'CROSS#1#c': ['21'], 'NOREF#1#c': ['21'], 'PLAIN#1#c': ['24'] },
  },
  {
    fixture: 'pickwrong.R',
    windowBlocks: ['blocks 12'],
    refName: 'R#0#chr1',
    window: { start: 12, end: 24 },
    held: { 'H#1#c': ['H#1#c:40-47'] },
    nodes: { 'H#1#c': ['13'] },
  },
  {
    fixture: 'twopieces.R',
    windowBlocks: ['blocks 24'],
    refName: 'R#0#chr1',
    window: { start: 24, end: 36 },
    held: { 'H#1#c': ['H#1#c:6-16', 'H#1#c:30-37'] },
    nodes: { 'H#1#c': ['11', '13'] },
  },
])(
  'a record’s own steps in a tile fetched to attach them come from the abridged row and the record, with no blocks unit read, as the expanded cut reads them: $fixture',
  async ({ fixture, refName, window, windowBlocks, held, nodes }) => {
    const ownSteps = (cut: Cut) =>
      Object.fromEntries(
        cut.departures
          .filter(d => d.steps)
          .map(d => [d.name, [...d.steps!.ids].map(id => nameOf(cut, id))]),
      )
    const { reader: walks, fetched, watch } = counted(fixture)
    await watch(refName)
    const tilesRead = () => fetched.filter(u => !u.startsWith('departures'))
    const unloaded = await walks.cut(refName, window, { context: 0 })
    expect(ownSteps(unloaded)).toEqual(held)
    expect(tilesRead()).toContain('tile 0')
    expect(tilesRead().filter(u => u.startsWith('blocks'))).toEqual(
      windowBlocks,
    )
    fetched.length = 0
    const loaded = await walks.cut(refName, window, { context: 0, expand: all })
    expect(ownSteps(loaded)).toEqual(nodes)
    expect(tilesRead().filter(u => u.startsWith('blocks'))).toEqual([])
    expect(loaded.placeholders).toBe(0)
    expect(stated(unloaded, refName)).toEqual(stated(loaded, refName))
  },
)

// twopiecesonetile.R.gfa: W H 1 c >1>11 and, from 30, >13>2, both pieces in
// the tile at 0, which the window reads with its blocks unit
test('a record’s own steps in two pieces of a window’s tile are each piece’s run at one end, by its blocks row and by its abridged row alike', async () => {
  const walks = reader('twopiecesonetile.R')
  const window = { start: 0, end: 12 }
  const unloaded = await walks.cut('R#0#chr1', window, { context: 0 })
  const loaded = await walks.cut('R#0#chr1', window, {
    context: 0,
    expand: all,
  })
  expect(
    unloaded.departures.map(d =>
      [...(d.steps?.ids ?? [])].map(id => nameOf(unloaded, id)),
    ),
  ).toEqual([['H#1#c:6-16', 'H#1#c:30-37']])
  expect(loaded.departures.map(d => [...(d.steps?.ids ?? [])])).toEqual([
    [11, 13],
  ])
  expect(stated(unloaded, 'R#0#chr1')).toEqual(stated(loaded, 'R#0#chr1'))
})

// edge.R.gfa: R as in blocks.R, W LEAD 1 c >1>2>27>7>8 with 5 bp in 27 and
// W TRAIL 1 c >7>8>28>1>2 with 3 bp in 28. Over 60-80 with 10 bp of context,
// contained, LEAD's piece in the tile at 0 lies off the drawn stretch and its
// last 5 bp lie within the context of its next piece's first step
test.each([{ attachBp: 1 }, { attachBp: undefined }])(
  'a contained cut builds the piece whose own steps its context reaches before a piece over the drawn stretch, whatever attachBp: $attachBp',
  async ({ attachBp }) => {
    const walks = reader('edge.R')
    const window = { start: 60, end: 80 }
    const opts = { context: 10, contained: true, attachBp }
    const unloaded = await walks.cut(R, window, opts)
    const loaded = await walks.cut(R, window, { ...opts, expand: all })
    expect(rows(loaded, 'LEAD#1#c')).toEqual([[20, 45, 'off:5 >7@60 >8@70']])
    expect(rows(loaded, 'TRAIL#1#c')).toEqual([[0, 23, '>7@60 >8@70 off:3']])
    expect(
      loaded.fragments.flatMap(f => [...f.ids]).filter(id => id < 0),
    ).toEqual([])
    expect(loaded.placeholders).toBe(0)
    const [lead] = unloaded.fragments.filter(f => f.name === 'LEAD#1#c')
    expect(nameOf(unloaded, lead!.ids[0]!)).toBe('LEAD#1#c:20-25')
    expect(stated(unloaded, R)).toEqual(stated(loaded, R))
  },
)

test('a cut that expands a walk attaches the steps its unexpanded cut holds as placeholders, tile by tile', async () => {
  const walks = reader('ownends.R')
  for (const start of [0, 6, 12, 18, 24, 30]) {
    const window = { start, end: start + 6 }
    expect(stated(await walks.cut(R, window, { context: 0 }), R)).toEqual(
      stated(await walks.cut(R, window, { context: 0, expand: all }), R),
    )
  }
})

// W PARTIAL 1 c >7>8>9>23>10>11>12, 50 bp in the one node 23 before 10, and
// W JUMP 1 c >7>8>22>12, 5 bp in 22 before 12
test('a contained cut keeps its context’s bp of a placeholder, where the loaded cut keeps the node that straddles it whole', async () => {
  const { unloaded, loaded } = await cuts(100, 110, {
    contained: true,
    context: 10,
  })
  expect(rows(unloaded, 'PARTIAL#1#c')).toEqual([
    [70, 110, 'off:10 >10@90 >11@100 >12@110'],
  ])
  expect(rows(loaded, 'PARTIAL#1#c')).toEqual([
    [30, 110, 'off:50 >10@90 >11@100 >12@110'],
  ])
  const [f] = unloaded.fragments.filter(x => x.name === 'PARTIAL#1#c')
  expect(nameOf(unloaded, f!.ids[0]!)).toBe('PARTIAL#1#c:70-80')
  expect([f!.trimmedBefore, f!.trimmedAfter]).toEqual([70, 0])
  const [d] = unloaded.departures.filter(x => x.name === 'PARTIAL#1#c')
  expect([...d!.steps!.ids].map(id => nameOf(unloaded, id))).toEqual([
    'PARTIAL#1#c:30-70',
    'PARTIAL#1#c:70-80',
  ])
  expect(
    unloaded.gaps
      .filter(g => g.name === 'PARTIAL#1#c')
      .map(g => [g.edge, g.start, g.end, g.segments.map(x => [x.kind, x.bp])]),
  ).toEqual([
    [
      1,
      0,
      70,
      [
        ['read', 30],
        ['steps', 40],
      ],
    ],
  ])
  const both = [stated(unloaded, R), stated(loaded, R)]
  const [a, b] = both.map(
    ({ departures: left, kept, links: joined, haplotypes: read }) => ({
      left,
      kept,
      joined,
      read,
    }),
  )
  expect(a).toEqual(b)
  const others = both.map(x =>
    x.fragments.filter(row => !row.startsWith('PARTIAL')),
  )
  expect(others[0]).toEqual(others[1])
})

test('a contained cut keeps a placeholder within its context whole, as the loaded cut keeps its nodes', async () => {
  const { unloaded, loaded } = await cuts(110, 120, {
    contained: true,
    context: 5,
  })
  expect(rows(unloaded, 'JUMP#1#c')).toEqual([[20, 35, 'off:5 >12@110']])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

test('a placeholder is no node: it is flagged, as long as its run, placed on its walk’s contig and named for both', async () => {
  const walks = reader('blocks.R')
  const window = { start: 60, end: 120 }
  const unloaded = await walks.cut(R, window, { context: 0 })
  const loaded = await walks.cut(R, window, { context: 0, expand: all })
  const placeholders = unloaded.fragments.flatMap(f =>
    [...f.ids].filter(isPlaceholder).map(id => ({ id, name: f.name })),
  )
  expect(placeholders.map(p => p.name).sort()).toEqual([
    'HAIRPIN#1#c',
    'JUMP#1#c',
    'MULTI#1#c',
    'PARTIAL#1#c',
    'TWICE#1#c',
  ])
  expect(new Set(placeholders.map(p => p.id)).size).toBe(placeholders.length)
  expect(placeholders.every(p => p.id < 0)).toBe(true)
  expect(
    placeholders.map(({ id, name }) => {
      const node = unloaded.graph.nodes.get(id)!
      return [
        name,
        node.placeholder,
        node.refName,
        node.start,
        node.end,
        node.rank,
        node.passes,
        node.onWindow,
      ]
    }),
  ).toEqual(
    expect.arrayContaining([
      ['HAIRPIN#1#c', true, 'HAIRPIN#1#c', 40, 47, 1, 0, false],
      ['PARTIAL#1#c', true, 'PARTIAL#1#c', 30, 80, 1, 0, false],
    ]),
  )
  const real = [...unloaded.graph.nodes].filter(([, n]) => !n.placeholder)
  expect(real.every(([id]) => id >= 0)).toBe(true)
  expect(loaded.fragments.some(f => f.ids.some(isPlaceholder))).toBe(false)
  expect(loaded.fragments.some(f => f.abridged)).toBe(false)
  const tables = unloaded.graph.tables(
    unloaded.kept,
    unloaded.fragments,
    unloaded.departures,
    unloaded.gaps,
  )
  const flagged = tables.nodes.names.filter(
    (_, i) => tables.nodes.placeholders![i] === 1,
  )
  expect(flagged).toEqual(
    expect.arrayContaining([
      'HAIRPIN#1#c:40-47',
      'JUMP#1#c:20-25',
      'MULTI#1#c:50-53',
      'PARTIAL#1#c:30-80',
      'TWICE#1#c:30-36',
    ]),
  )
  expect(flagged.some(name => /^-?\d+$/.test(name))).toBe(false)
})

// Every walk of blocks.R over 60-120 keeps placeholders for its runs off R
test('placeholders do not count towards the node limit, and a cut says how many it returns', async () => {
  const walks = reader('blocks.R')
  const window = { start: 60, end: 120 }
  const cut = await walks.cut(R, window, { context: 0 })
  const placeholders = [...cut.kept].filter(isPlaceholder).length
  const nodes = cut.kept.size - placeholders
  expect(placeholders).toBeGreaterThan(4)
  expect(cut.placeholders).toBe(placeholders)
  expect([...cut.kept].every(id => cut.graph.nodes.has(id))).toBe(true)
  const within = await walks.cut(R, window, { context: 0, nodeLimit: nodes })
  expect(within.kept.size).toBe(cut.kept.size)
  await expect(
    walks.cut(R, window, { context: 0, nodeLimit: nodes - 1 }),
  ).rejects.toThrow('Zoom in')
  const loaded = await walks.cut(R, window, { context: 0, expand: all })
  expect(loaded.placeholders).toBe(0)
})

test('a step built for an unloaded piece counts towards the step budget as any step does', async () => {
  const walks = reader('blocks.R')
  const window = { start: 60, end: 120 }
  const cut = await walks.cut(R, window, { context: 0 })
  const held =
    cut.fragments.reduce((sum, f) => sum + f.ids.length, 0) +
    cut.departures.reduce((sum, d) => sum + (d.steps?.ids.length ?? 0), 0)
  expect(
    cut.fragments.reduce((sum, f) => sum + (f.unloadedSteps ?? 0), 0),
  ).toBe(
    cut.fragments
      .filter(f => f.name !== R)
      .reduce((sum, f) => sum + f.ids.length, 0),
  )
  const within = await walks.cut(R, window, { context: 0, stepBudget: held })
  expect(within.fragments).toHaveLength(cut.fragments.length)
  await expect(
    walks.cut(R, window, { context: 0, stepBudget: held - 1 }),
  ).rejects.toThrow('Zoom in')
})

// ownends.R.gfa, X#1 and R alone over 24-30: the cut reads the tile at 20
// and its blocks unit, then the tile at 0 for X's record
test('the bytes of a tile fetched for own steps are checked with the blocks unit read before it, before it is fetched', async () => {
  const { reader: walks, fetched, watch } = counted('ownends.R')
  await watch(R)
  const units = await walks.units(R)
  const bytes = (kind: string, name: string) =>
    units.find(u => u.kind === kind && u.name === name)!.length
  const cut = (byteBudget: number) =>
    walks.cut(
      R,
      { start: 24, end: 30 },
      {
        context: 0,
        keep: name => name === 'X#1#c' || name === R,
        byteBudget,
      },
    )
  const read = () => fetched.filter(u => !u.startsWith('departures')).sort()
  const window = bytes('tile', '20') + bytes('blocks', '20')
  const whole = window + bytes('tile', '0')
  expect((await cut(whole)).departures).toHaveLength(1)
  expect(read()).toEqual(['blocks 20', 'tile 0', 'tile 20'])
  fetched.length = 0
  await expect(cut(whole - 1)).rejects.toBeInstanceOf(RegionTooLargeError)
  expect(read()).toEqual(['blocks 20', 'tile 20'])
  fetched.length = 0
  await expect(cut(window - 1)).rejects.toBeInstanceOf(RegionTooLargeError)
  expect(read()).toEqual(['tile 20'])
})

// One contig in two W lines of one tile, both abridged: >5>6 then 30 bp off
// R in node 20 to its line's end, 50 bp no line covers, and >1>2. Its record
// from 6 to 1 owns the 30 bp, which trail a piece with no step on a window
// over 1
test('a cut that expands a walk reads the piece a record attaches steps of, though it has no step on the window', async () => {
  const walks = WalkReader.fromFilehandles(
    walkFileOf([
      '#walks format:i:8 tile:i:60',
      '#reference R',
      '#haplotype W#1',
      'R#0#c 0 ref R#0#c 0 0 0 6 2,2,2,2,2,2 10,10,10,10,10,10 .',
      'R#0#c 0 abridged W#1#c 0 0 0 3 50 0',
      'R#0#c 0 abridged W#1#c 100 100 0 2 20 1',
      'R#0#c 0:0 40 60 nodes 20 30',
      'R#0#c 0:0 40 60 route 12 28 .',
      'R#0#c 0:0 40 60 walk W#1#c 0 0 0 3 10 0',
      'R#0#c 0:1 0 20 walk W#1#c 100 100 0 2 2 .',
      'R#0#c 0 40 60 blocks W#1#c 0 0 40+20,30',
      'R#0#c 0 0 20 blocks W#1#c 100 0 0+20',
      'R#0#c 0 60 departure W#1#c 20 80 50 60 + 1 0 + 1 0',
    ]),
  )
  const window = { start: 0, end: 10 }
  const loaded = await walks.cut(R, window, { context: 0, expand: all })
  const unloaded = await walks.cut(R, window, { context: 0 })
  const held = (cut: Cut) =>
    cut.departures.map(d => [
      d.hapStart,
      d.ownBp,
      d.missingBp,
      [...(d.steps?.ids ?? [])].map(id => nameOf(cut, id)),
    ])
  expect(held(loaded)).toEqual([[20, 80, 50, ['20']]])
  expect(held(unloaded)).toEqual([[20, 80, 50, ['W#1#c:20-50']]])
  expect(loaded.placeholders).toBe(0)
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// between.R.gfa: R as in blocks.R. W MID 1 c >7>8>1>27>2>9 has pieces 7,8 in
// the tile at 60, then 1,27,2 in the tile at 0, with 5 bp in 27, then 9 in the
// tile at 60; W AFTER 1 c >7>8>1>28>2 ends on its piece in the tile at 0, 23
// bp. A window on 60-90 with 10 bp of context reads the tile at 0 and draws
// none of 0-20
test('a cut that expands a walk loads its piece between two over the drawn stretch, which lies off it, and of one past its last states the bp, with no blocks unit read', async () => {
  const { reader: walks, fetched, watch } = counted('between.R')
  await watch(R)
  const window = { start: 60, end: 90 }
  const loaded = await walks.cut(R, window, { context: 10, expand: all })
  expect(fetched.sort()).toEqual([
    'departures .',
    'departures .',
    'detail 0:1',
    'detail 60:0',
    'detail 60:1',
    'tile 0',
    'tile 60',
  ])
  const [mid] = loaded.fragments.filter(f => f.name === 'MID#1#c')
  expect([...mid!.ids]).toEqual([7, 8, 1, 27, 2, 9])
  const [after] = loaded.fragments.filter(f => f.name === 'AFTER#1#c')
  expect([[...after!.ids], after!.hapEnd, after!.trimmedAfter]).toEqual([
    [7, 8],
    20,
    23,
  ])
  expect(loaded.fragments.some(f => f.abridged)).toBe(false)
  expect(loaded.placeholders).toBe(0)
  const unloaded = await walks.cut(R, window, { context: 10 })
  expect(rows(unloaded, 'MID#1#c')).toEqual([
    [0, 55, '>7@60 >8@70 >1@0 off:5 >2@10 >9@80'],
  ])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

// unbuilt.R.gfa: R as in blocks.R. W BACK 1 c >1>2>7>8>4>6 has pieces 1,2 and
// 4,6 in the tile at 0, both in one detail unit, and the second skips node 5.
// W LOOP 1 c >1>2>3>4>5>6>1 passes node 1 twice, so a cut of 0-10 keeps nodes
// 4 and 6 by LOOP's steps, and of BACK builds the first piece alone
test('a cut that expands a walk holds no link of a piece it does not build, where the piece’s detail unit holds one it builds', async () => {
  const walks = reader('unbuilt.R')
  const window = { start: 0, end: 10 }
  const loaded = await walks.cut(R, window, { context: 0, expand: all })
  const unloaded = await walks.cut(R, window, { context: 0 })
  expect(rows(loaded, 'BACK#1#c')).toEqual([[0, 10, '>1@0']])
  expect([...loaded.kept].sort((a, b) => a - b)).toEqual([1, 2, 3, 4, 5, 6])
  expect(links(loaded)).toEqual([
    '1[ 6]',
    '1] 2[',
    '2] 3[',
    '3] 4[',
    '4] 5[',
    '5] 6[',
  ])
  expect(links(unloaded)).toEqual(links(loaded))
})

test('a piece in two fragments of a contained cut counts its steps on the first', async () => {
  const { unloaded } = await cuts(80, 90, { contained: true })
  const { walks } = unloaded.graph.tables(
    unloaded.kept,
    unloaded.fragments,
    unloaded.departures,
  )
  expect(
    walks.names.flatMap((name, i) =>
      name === 'TWICE#1#c' ? [[walks.starts[i], walks.abridgedSteps![i]]] : [],
    ),
  ).toEqual([
    [10, 7],
    [46, 0],
  ])
})

describe('the units a cut reads', () => {
  const kinds = (fetched: string[]) =>
    fetched.map(u => u.split(' ')[0] ?? '').sort()

  test('a cut that leaves a kept piece unloaded reads its tile’s blocks unit, and no detail unit', async () => {
    const { reader: walks, fetched, watch } = counted('cover.R')
    await watch(R)
    await walks.cut(R, { start: 6000, end: 6500 })
    expect(fetched.sort()).toEqual(['blocks 0', 'tile 0'])
  })

  test('a cut that expands every walk reads the detail units of the pieces whose span lies over the drawn stretch, no other, and no blocks unit', async () => {
    const { reader: walks, fetched, watch } = counted('cover.R')
    await watch(R)
    const { fragments } = await walks.cut(
      R,
      { start: 6000, end: 6500 },
      { expand: all },
    )
    expect(fetched.sort()).toEqual(['detail 0:0', 'detail 0:1', 'tile 0'])
    expect(fragments.map(f => f.name).sort()).toEqual([
      'FAR#1#c',
      'NEAR#1#c',
      'R#0#c',
    ])
  })

  test('a cut that expands the walks it keeps reads their detail units by their abridged rows', async () => {
    const { reader: walks, fetched, watch } = counted('cover.R')
    await watch(R)
    const { fragments } = await walks.cut(
      R,
      { start: 6000, end: 6500 },
      { keep: name => name !== 'OFF#1#c', expand: all },
    )
    expect(fetched.sort()).toEqual(['detail 0:0', 'detail 0:1', 'tile 0'])
    expect(fragments.some(f => f.abridged)).toBe(false)
  })

  test('a cut that expands some walks reads the blocks unit for the others', async () => {
    const { reader: walks, fetched, watch } = counted('cover.R')
    await watch(R)
    const { fragments } = await walks.cut(
      R,
      { start: 6000, end: 6500 },
      { expand: name => name === 'FAR#1#c' },
    )
    expect(fetched.sort()).toEqual(['blocks 0', 'detail 0:0', 'tile 0'])
    expect(
      fragments
        .map(f => `${f.name} ${f.abridged ? 'unloaded' : 'loaded'}`)
        .sort(),
    ).toEqual(['FAR#1#c loaded', 'NEAR#1#c unloaded', 'R#0#c loaded'])
  })

  test.each([
    { expand: undefined, read: ['blocks', 'tile'] },
    { expand: all, read: ['tile'] },
  ])(
    'a cut whose tile abridges no piece over the window and its context reads no detail unit, and the tile’s blocks unit unless it expands every walk: $read',
    async ({ expand, read }) => {
      const { reader: walks, fetched, watch } = counted('cover.R')
      await watch(R)
      const { fragments } = await walks.cut(
        R,
        { start: 9000, end: 9100 },
        { context: 0, expand },
      )
      expect(kinds(fetched)).toEqual(read)
      expect(fragments.map(f => f.name)).toEqual(['R#0#c'])
    },
  )

  test('the blocks unit counts against the byte budget, and in the bytes a cut is said to read', async () => {
    const walks = reader('cover.R')
    const units = await walks.units(R)
    const bytes = (kind: string) =>
      units.filter(u => u.kind === kind).reduce((sum, u) => sum + u.length, 0)
    const window = { start: 6000, end: 6500 }
    expect(await walks.bytesFor(R, window)).toBe(
      bytes('tile') + bytes('blocks'),
    )
    await expect(
      walks.cut(R, window, { byteBudget: bytes('tile') }),
    ).rejects.toThrow('fetch')
    const { fragments } = await walks.cut(R, window, {
      byteBudget: bytes('tile') + bytes('blocks'),
    })
    expect(fragments).toHaveLength(3)
  })
})

describe('a file whose blocks unit is not what its tile needs', () => {
  const header = [
    '#walks format:i:8 tile:i:20',
    '#reference R',
    '#haplotype H#1',
  ]
  const abridged = [
    'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 5,5,5,5 .',
    'R#0#chr 0 abridged H#1#h 0 0 0 3 18 0',
    'R#0#chr 0:0 10 20 nodes 5 8',
    'R#0#chr 0:0 10 20 route 6 4 -2',
    'R#0#chr 0:0 10 20 walk H#1#h 0 0 0 3 6 0',
  ]
  const blocks = 'R#0#chr 0 10 20 blocks H#1#h 0 0 10+5,8,0+5'
  const open = (lines: string[]) =>
    WalkReader.fromFilehandles(walkFileOf([...header, ...lines]))
  const window = { start: 0, end: 20 }

  test('a sequence with its blocks unit is read', async () => {
    const cut = await open([...abridged, blocks]).cut('R#0#chr', window, {
      context: 0,
    })
    expect(rows(cut, 'H#1#h')).toEqual([[0, 18, '>3@10 off:8 >4@15']])
  })

  test('a tile with an abridged row and no blocks unit is refused by a cut that leaves the piece unloaded, and read by one that expands it', async () => {
    const walks = open(abridged)
    const error: unknown = await walks
      .cut('R#0#chr', window)
      .catch((e: unknown) => e)
    expect(error).toBeInstanceOf(WalkFormatError)
    expect(error).toHaveProperty(
      'message',
      "R#0#chr's unit table lists no blocks unit 0, which the abridged rows of that tile need",
    )
    const cut = await walks.cut('R#0#chr', window, { context: 0, expand: all })
    expect(rows(cut, 'H#1#h')).toEqual([[0, 18, '>3@10 off:8 >4@15']])
  })

  test.each([
    ['spanning other than its abridged row states', '5+5,8,5+5'],
    ['covering part of a node', '10+5,10,0+3'],
    ['that is no list', '10+5;0+5'],
    ['with a run off R after a colon', '10+5,0+5:8'],
    ['with two runs off R side by side', '10+5,3,5,0+5'],
    ['with no block', '18'],
  ])('a blocks row %s is refused', async (_, list) => {
    await expect(
      open([...abridged, blocks.replace('10+5,8,0+5', list)]).cut(
        'R#0#chr',
        window,
      ),
    ).rejects.toThrow(/blocks row/)
  })

  const nextPiece = (hapOffset: number, expand?: () => boolean) =>
    open([
      ...abridged,
      blocks,
      'R#0#chr 20 ref R#0#chr 0 20 1 4 10,2,2,2 5,5,5,5 .',
      `R#0#chr 20 walk H#1#h 0 ${hapOffset} 1 2 10 .`,
    ]).cut('R#0#chr', { start: 0, end: 40 }, { context: 0, expand })

  test('an abridged row that ends where the line’s next piece starts is read', async () => {
    expect((await nextPiece(18)).placeholders).toBe(1)
  })

  test.each([undefined, all])(
    'an abridged row that ends elsewhere than the line’s next piece starts is refused, by its hapBp: expand %o',
    async expand => {
      await expect(nextPiece(20, expand)).rejects.toThrow(
        "R#0#chr's abridged row of H#1#h, piece 0 in tile 0, ends at 18 on its contig, and the line's next piece starts at 20",
      )
    },
  )

  test('a blocks unit without the row of a piece its tile abridges is refused', async () => {
    await expect(
      open([...abridged, blocks.replace('H#1#h 0 0', 'H#1#h 0 1')]).cut(
        'R#0#chr',
        window,
      ),
    ).rejects.toThrow('holds no row for piece 0 of H#1#h')
  })
})

// One contig in two W lines of one tile, both abridged, as above: >5>6 and
// 30 bp off R to its line's end at 50, and >1>2 from 100. Its record leaves
// 6 at 20 on the contig and owns the 80 bp to 100
describe('a file whose abridged rows, blocks rows, detail rows and records disagree', () => {
  const lines = [
    '#walks format:i:8 tile:i:60',
    '#reference R',
    '#haplotype W#1',
    'R#0#c 0 ref R#0#c 0 0 0 6 2,2,2,2,2,2 10,10,10,10,10,10 .',
    'R#0#c 0 abridged W#1#c 0 0 0 3 50 0',
    'R#0#c 0 abridged W#1#c 100 100 0 2 20 1',
    'R#0#c 0:0 40 60 nodes 20 30',
    'R#0#c 0:0 40 60 route 12 28 .',
    'R#0#c 0:0 40 60 walk W#1#c 0 0 0 3 10 0',
    'R#0#c 0:1 0 20 walk W#1#c 100 100 0 2 2 .',
    'R#0#c 0 40 60 blocks W#1#c 0 0 40+20,30',
    'R#0#c 0 0 20 blocks W#1#c 100 0 0+20',
    'R#0#c 0 60 departure W#1#c 20 80 50 60 + 1 0 + 1 0',
  ]
  const cut = (edit: (line: string) => string, expand?: () => boolean) =>
    WalkReader.fromFilehandles(walkFileOf(lines.map(edit))).cut(
      R,
      { start: 0, end: 10 },
      { context: 0, expand },
    )
  const hapBp49 = (line: string) =>
    line.replace('abridged W#1#c 0 0 0 3 50 0', 'abridged W#1#c 0 0 0 3 49 0')
  const record = (hapStart: number, ownBp: number) => (line: string) =>
    line.replace(
      'departure W#1#c 20 80',
      `departure W#1#c ${hapStart} ${ownBp}`,
    )

  test('the file as written is read, unexpanded and expanded', async () => {
    expect((await cut(line => line)).placeholders).toBe(1)
    expect((await cut(line => line, all)).placeholders).toBe(0)
  })

  test('a blocks row that does not add up to its abridged row’s hapBp is refused', async () => {
    const error: unknown = await cut(hapBp49).catch((e: unknown) => e)
    expect(error).toBeInstanceOf(WalkFormatError)
    expect(error).toHaveProperty(
      'message',
      "R#0#c's blocks row of W#1#c, piece 0 in tile 0, adds up to 50 bp, and its abridged row states a hapBp of 49",
    )
  })

  test('a detail row whose steps do not add up to its abridged row’s hapBp is refused', async () => {
    await expect(cut(hapBp49, all)).rejects.toThrow(
      "R#0#c's detail unit 0:0 holds 50 bp of steps for piece 0 of W#1#c, and its abridged row states a hapBp of 49",
    )
  })

  test('a record whose own steps in a piece are not the run its blocks row ends with is refused', async () => {
    await expect(cut(record(25, 75))).rejects.toThrow(
      "R#0#c's blocks row of W#1#c, piece 0 in tile 0, ends with 30 bp off the reference, and a departure record's own steps there are 25 bp by its abridged row",
    )
  })

  test.each([
    { hapStart: 20, ownBp: 10, lies: 'inside piece 0, which covers [0, 50)' },
    {
      hapStart: 0,
      ownBp: 60,
      lies: 'over all of piece 0, which covers [0, 50)',
    },
  ])(
    'a record whose own steps lie $lies is refused, unexpanded and expanded',
    async ({ hapStart, ownBp, lies }) => {
      for (const expand of [undefined, all]) {
        const error: unknown = await cut(record(hapStart, ownBp), expand).catch(
          (e: unknown) => e,
        )
        expect(error).toBeInstanceOf(WalkRowError)
        expect(error).toHaveProperty(
          'message',
          expect.stringContaining(`of its contig lies ${lies}`),
        )
      }
    },
  )
})

// R over two 10 bp nodes: a record leaving along R at 15 names no node's end
test('own steps beside a leave or rejoin that is no edge of a node of the tile’s reference path are refused', () => {
  const path = new RefPath()
  path.add([2, 4], [10, 10], 0)
  const [piece] = decodeChunk(
    [
      'R#0#c\t0\t20\tref\tR#0#c\t0\t0\t0\t2\t2,2\t10,10\t.'.split('\t'),
      'R#0#c\t0\t20\tabridged\tW#1#c\t0\t0\t0\t3\t25\t0'.split('\t'),
    ],
    0,
    undefined,
  ).abridged
  const run = (pos: number, reversed: boolean) =>
    ownRunsRow(
      piece!,
      [{ key: 'k', start: 20, end: 25, trailing: true, pos, reversed }],
      path,
      () => -1,
    )
  const own = [{ key: 'k', leave: true, from: 1, to: 2 }]
  expect(run(20, false).row.own).toEqual(own)
  expect(run(10, true).row.own).toEqual(own)
  expect(() => run(15, false)).toThrow(
    'a departure record of W#1#c leaves at 15, where the reference path of the tile holding its own steps has no node ending',
  )
  expect(() => run(15, true)).toThrow('has no node starting')
})

// The walks holding a placeholder among their steps, a departure's or a
// gap's, whose walks' and departures' abridgedSteps state none of its steps
function uncounted(cut: Cut) {
  const tables = cut.graph.tables(
    cut.kept,
    cut.fragments,
    cut.departures,
    cut.gaps,
  )
  const holding = new Set<string>()
  const scan = (name: string, ids: Int32Array, from: number, to: number) => {
    for (let j = from; j < to; j++) {
      if (tables.nodes.placeholders![ids[j]!] === 1) {
        holding.add(name)
      }
    }
  }
  const { walks } = tables
  walks.names.forEach((name, i) => {
    scan(name, walks.steps, walks.offsets[i]!, walks.offsets[i + 1]!)
  })
  const records = tables.departures
  records.names.forEach((name, k) => {
    scan(name, records.steps, records.offsets[k]!, records.offsets[k + 1]!)
  })
  const { gaps } = tables
  gaps?.names.forEach((name, g) => {
    const { stepOffsets } = gaps.segments
    scan(
      name,
      gaps.segments.steps,
      stepOffsets[gaps.segmentOffsets[g]!]!,
      stepOffsets[gaps.segmentOffsets[g + 1]!]!,
    )
  })
  const said = new Map<string, number>()
  const sum = (names: string[], count: Float64Array | undefined) => {
    names.forEach((name, i) => {
      const n = count?.[i] ?? Number.NaN
      if (!Number.isNaN(n)) {
        said.set(name, (said.get(name) ?? 0) + n)
      }
    })
  }
  sum(walks.names, walks.abridgedSteps)
  sum(records.names, records.abridgedSteps)
  return [...holding].filter(name => !((said.get(name) ?? 0) > 0))
}

test.each([
  'abredge.R',
  'abridged.R',
  'blocks.R',
  'contigends.R',
  'hairpin.R',
  'interleaved.R',
  'ownends.R',
])(
  'every walk holding a placeholder in a cut of %s, in its steps, a departure’s or a gap’s, states its unloaded steps',
  async name => {
    const walks = reader(name)
    const missing: string[] = []
    for (const refName of await walks.referenceSequences()) {
      const length = await walks.sequenceLength(refName)
      for (let start = 0; start < length; start += 5) {
        for (const contained of [false, true]) {
          const cut = await walks.cut(
            refName,
            { start, end: start + 5 },
            { context: 0, contained },
          )
          for (const walk of uncounted(cut)) {
            missing.push(`${refName}:${start} ${walk}`)
          }
        }
      }
    }
    expect(missing).toEqual([])
  },
)
