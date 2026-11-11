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
  let own = 0
  for (let i = 0; i <= f.ids.length; i++) {
    const offR = i < f.ids.length && Number.isNaN(f.refStarts[i])
    if (offR) {
      const node = cut.graph.nodes.get(f.ids[i]!)!
      own += node.end - node.start
    } else {
      if (own > 0) {
        out.push(`own:${own}`)
      }
      own = 0
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

// Everything a cut states that a cut loading the pieces must state too
const stated = (cut: Cut, refName: string) => ({
  fragments: cut.fragments
    .filter(f => f.name !== refName)
    .map(f => `${f.name} ${f.hapStart}-${f.hapEnd} ${steps(cut, f)}`)
    .sort(),
  departures: cut.departures.flatMap(d => departures(cut, d.name)),
  haplotypes: haplotypes(cut, refName),
})

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
  { start: 90, end: 100, row: [30, 57, '>10@90 own:7 <10@90'] },
  { start: 80, end: 90, row: [20, 67, '>9@80 >10@90 own:7 <10@90 <9@80'] },
  {
    start: 60,
    end: 70,
    row: [0, 87, '>7@60 >8@70 >9@80 >10@90 own:7 <10@90 <9@80 <8@70 <7@60'],
  },
])(
  'a hairpin that starts and ends on one node at its tile’s start is kept where it steps, both passes and the own bp between them, cut $start-$end',
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
    [0, 35, '>7@60 >8@70 own:5 >12@110'],
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
    [40, 63, '>9@80 own:3 >10@90'],
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
test('an unloaded piece joins the next piece read, the own bp between them a placeholder its record holds', async () => {
  const { unloaded, loaded } = await cuts(50, 70)
  expect(rows(unloaded, 'TAIL#1#c')).toEqual([[10, 34, '>6@50 own:4 >7@60']])
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
// 6 bp nodes, and W H 1 c >1>2<2<1>11>5>6>3 with 10 bp in 11, whose blocks
// are 0+12,0-12,24+12:10,-24+6
test('a piece that turns back, runs own bp and jumps back is its blocks in path order', async () => {
  const walks = reader('hairpin.R')
  const window = { start: 0, end: 36 }
  const unloaded = await walks.cut(R, window, { context: 0 })
  const loaded = await walks.cut(R, window, { context: 0, expand: all })
  expect(rows(unloaded, 'H#1#c')).toEqual([
    [0, 52, '>1@0 >2@6 <2@6 <1@0 own:10 >5@24 >6@30 >3@12'],
  ])
  expect(stated(unloaded, R)).toEqual(stated(loaded, R))
})

test('a placeholder is no node: its id is below every node’s and its own, it is flagged, as long as its run and placed on its walk’s contig', async () => {
  const { unloaded, loaded } = await cuts(60, 120)
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
  const flagged = [...tables.nodes.placeholders!].flatMap((flag, i) =>
    flag === 1 ? [Number(tables.nodes.names[i])] : [],
  )
  expect(flagged.sort((a, b) => a - b)).toEqual(
    placeholders.map(p => p.id).sort((a, b) => a - b),
  )
  expect(tables.links.from).toHaveLength(tables.links.to.length)
  const linked = new Set([...tables.links.from, ...tables.links.to])
  expect([...linked].some(i => tables.nodes.placeholders![i] === 1)).toBe(false)
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

describe('the units a cut reads', () => {
  const kinds = (fetched: string[]) =>
    fetched.map(u => u.split(' ')[0] ?? '').sort()

  test('a cut that leaves a kept piece unloaded reads its tile’s blocks unit, and no detail unit', async () => {
    const { reader: walks, fetched, watch } = counted('cover.R')
    await watch(R)
    await walks.cut(R, { start: 6000, end: 6500 })
    expect(fetched.sort()).toEqual(['blocks 0', 'tile 0'])
  })

  test('a cut that expands every piece it keeps reads no blocks unit, and only the detail units whose extent overlaps the window and its context', async () => {
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
    'a cut whose tile abridges no piece over the window and its context reads neither a blocks unit nor a detail unit: %o',
    async ({ expand }) => {
      const { reader: walks, fetched, watch } = counted('cover.R')
      await watch(R)
      const { fragments } = await walks.cut(
        R,
        { start: 9000, end: 9100 },
        { context: 0, expand },
      )
      expect(kinds(fetched)).toEqual(['tile'])
      expect(fragments.map(f => f.name)).toEqual(['R#0#c'])
    },
  )

  test('the blocks unit counts against the byte budget', async () => {
    const walks = reader('cover.R')
    const units = await walks.units(R)
    const bytes = (kind: string) =>
      units.filter(u => u.kind === kind).reduce((sum, u) => sum + u.length, 0)
    const window = { start: 6000, end: 6500 }
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
  const blocks = 'R#0#chr 0 10 20 blocks H#1#h 0 0 10+5,0+5:8'
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

  test('its sequences with no detail unit read as before', async () => {
    const { fragments } = await open([...plain, ...abridged]).cut(
      'R#0#plain',
      window,
    )
    expect(fragments.map(f => f.name)).toEqual(['R#0#plain', 'H#1#p'])
  })

  test('the same sequence with its blocks unit is read', async () => {
    const cut = await open([...plain, ...abridged, blocks]).cut(
      'R#0#chr',
      window,
      { context: 0 },
    )
    expect(rows(cut, 'H#1#h')).toEqual([[0, 18, '>3@10 own:8 >4@15']])
  })

  test.each([
    ['holding more bp than its abridged row states', '10+5,0+5:9'],
    ['starting on another step', '5+5,5+5:8'],
    ['ending on another strand', '10+5,5-5:8'],
    ['covering part of a node', '10+5,0+3:10'],
    ['that is no list', '10+5;0+5'],
  ])('a blocks row %s is refused', async (_, list) => {
    await expect(
      open([...abridged, blocks.replace('10+5,0+5:8', list)]).cut(
        'R#0#chr',
        window,
      ),
    ).rejects.toThrow(/blocks row/)
  })

  test('a blocks unit without the row of a piece its tile abridges is refused', async () => {
    await expect(
      open([...abridged, blocks.replace('H#1#h 0 0', 'H#1#h 0 1')]).cut(
        'R#0#chr',
        window,
      ),
    ).rejects.toThrow('holds no row for piece 0 of H#1#h')
  })
})
