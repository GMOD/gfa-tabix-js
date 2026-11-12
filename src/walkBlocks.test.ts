import { LocalFile } from 'generic-filehandle2'

import { walkFileOf } from './testing.ts'
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
      .map(f => f.abridged!.map(a => [a.detail.name, a.steps, a.bp])),
  ).toEqual([[['60:1', 2, 20]], [['60:1', 3, 23]], [['60:1', 2, 20]]])
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
// for the record's own steps alone
test.each([
  { expand: undefined, held: ['placeholder'] },
  { expand: all, held: [11, 12, 13] },
])(
  'a record whose leave piece is abridged in a tile fetched for its own steps holds them, a placeholder unless expanded: $held',
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
    expect(fetched.filter(u => u.startsWith('blocks')).sort()).toEqual([
      'blocks 0',
      'blocks 20',
    ])
    expect(fetched.filter(u => u.startsWith('detail'))).toHaveLength(
      expand ? 2 : 0,
    )
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

test('a placeholder is no node: it is flagged, as long as its run, placed on its walk’s contig and named for both, and one run keeps one id from cut to cut', async () => {
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
  const idsByRun = (cut: Cut) =>
    Object.fromEntries(
      [...cut.graph.nodes].flatMap(([id, node]) =>
        node.placeholder ? [[nameOf(cut, id), id]] : [],
      ),
    )
  const narrower = await walks.cut(R, { start: 90, end: 100 }, { context: 0 })
  expect(idsByRun(unloaded)).toMatchObject(idsByRun(narrower))
  expect(Object.keys(idsByRun(narrower)).length).toBeGreaterThan(3)
})

// Every walk of blocks.R over 60-120 keeps placeholders for its runs off R
test('placeholders do not count towards the node limit', async () => {
  const walks = reader('blocks.R')
  const window = { start: 60, end: 120 }
  const { kept, graph } = await walks.cut(R, window, { context: 0 })
  const placeholders = [...kept].filter(isPlaceholder).length
  const nodes = kept.size - placeholders
  expect(placeholders).toBeGreaterThan(4)
  expect([...kept].every(id => graph.nodes.has(id))).toBe(true)
  const within = await walks.cut(R, window, { context: 0, nodeLimit: nodes })
  expect(within.kept.size).toBe(kept.size)
  await expect(
    walks.cut(R, window, { context: 0, nodeLimit: nodes - 1 }),
  ).rejects.toThrow('Zoom in')
})

test('every fragment holding an unloaded piece says so, with the steps loading it reads', async () => {
  const { unloaded } = await cuts(90, 100)
  const tables = unloaded.graph.tables(
    unloaded.kept,
    unloaded.fragments,
    unloaded.departures,
  )
  const walks = Object.fromEntries(
    tables.walks.names.map((name, i) => [name, tables.walks.abridgedSteps![i]]),
  )
  expect(walks).toEqual({
    'R#0#c': NaN,
    'HAIRPIN#1#c': 9,
    'LOWFIRST#1#c': 6,
    'MIXED#1#c': 5,
    'MULTI#1#c': 3,
    'PARTIAL#1#c': 7,
    'TWICE#1#c': 7,
  })
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

  test('a cut that expands every walk reads the detail units of the pieces with a block over the drawn stretch, and no other', async () => {
    const { reader: walks, fetched, watch } = counted('cover.R')
    await watch(R)
    const { fragments } = await walks.cut(
      R,
      { start: 6000, end: 6500 },
      { expand: all },
    )
    expect(fetched.sort()).toEqual([
      'blocks 0',
      'detail 0:0',
      'detail 0:1',
      'tile 0',
    ])
    expect(fragments.map(f => f.name).sort()).toEqual([
      'FAR#1#c',
      'NEAR#1#c',
      'R#0#c',
    ])
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

  test.each([{ expand: undefined }, { expand: all }])(
    'a cut whose tile abridges no piece over the window and its context reads its blocks unit and no detail unit: %o',
    async ({ expand }) => {
      const { reader: walks, fetched, watch } = counted('cover.R')
      await watch(R)
      const { fragments } = await walks.cut(
        R,
        { start: 9000, end: 9100 },
        { context: 0, expand },
      )
      expect(kinds(fetched)).toEqual(['blocks', 'tile'])
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

describe('a file without blocks units', () => {
  const header = [
    '#walks format:i:7 tile:i:20',
    '#reference R',
    '#haplotype H#1',
  ]
  const plain = [
    'R#0#plain 0 ref R#0#plain 0 0 0 4 2,2,2,2 5,5,5,5 .',
    'R#0#plain 0 walk H#1#p 0 0 0 4 2 .',
  ]
  const abridged = [
    'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 5,5,5,5 .',
    'R#0#chr 0 abridged H#1#h 0 0 0 3 18 6 8 0',
    'R#0#chr 0:0 10 20 nodes 5 8',
    'R#0#chr 0:0 10 20 route 6 10 -4',
    'R#0#chr 0:0 10 20 walk H#1#h 0 0 0 3 6 0',
  ]
  const blocks = 'R#0#chr 0 10 20 blocks H#1#h 0 0 10+5,8,0+5'
  const open = (lines: string[]) =>
    WalkReader.fromFilehandles(walkFileOf([...header, ...lines]))
  const window = { start: 0, end: 20 }

  test('a sequence with a detail unit and no blocks unit for its tile is refused, saying which gfa-to-tabix to rebuild with', async () => {
    const walks = open([...plain, ...abridged])
    const error: unknown = await walks
      .cut('R#0#chr', window)
      .catch((e: unknown) => e)
    expect(error).toBeInstanceOf(WalkFormatError)
    expect(error).toHaveProperty(
      'message',
      "R#0#chr's unit table lists the detail unit 0:0 and no blocks unit for its tile, as a walk file written before gfa-to-tabix 4.1 does; @jbrowse/graph-walks reads an abridged piece's place on the reference from its blocks, so rebuild the walk file with gfa-to-tabix 4.1 or later",
    )
    await expect(
      walks.cut('R#0#chr', window, { expand: all }),
    ).rejects.toBeInstanceOf(WalkFormatError)
  })

  test('the same sequence with its blocks unit is read', async () => {
    const cut = await open([...plain, ...abridged, blocks]).cut(
      'R#0#chr',
      window,
      { context: 0 },
    )
    expect(rows(cut, 'H#1#h')).toEqual([[0, 18, '>3@10 off:8 >4@15']])
  })

  test.each([
    ['holding more bp than its abridged row states', '10+5,9,0+5'],
    ['starting on another step', '5+5,8,5+5'],
    ['ending on another strand', '10+5,8,5-5'],
    ['covering part of a node', '10+5,10,0+3'],
    ['that is no list', '10+5;0+5'],
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

  test.each([
    ['with its bp off R after a colon', '0 0 10+5,0+5:8'],
    ['with ownBefore and ownAfter columns', '0 0 0 0 10+5,0+5:8'],
    ['with ownBefore and ownAfter columns and no bp off R', '0 0 0 0 10+5,0+5'],
  ])(
    'a blocks row of an earlier shape, %s, is refused, saying to rebuild the file',
    async (_, columns) => {
      const error: unknown = await open([
        ...abridged,
        blocks.replace('0 0 10+5,8,0+5', columns),
      ])
        .cut('R#0#chr', window)
        .catch((e: unknown) => e)
      expect(error).toBeInstanceOf(WalkFormatError)
      expect(error).toHaveProperty(
        'message',
        expect.stringContaining(
          'rebuild the walk file with gfa-to-tabix 4.1 or later',
        ),
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
