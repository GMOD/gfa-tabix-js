import { BlobFile, LocalFile } from 'generic-filehandle2'

import { bgzf, digestOf, graphTablesGFA, walkFileOf } from './testing.ts'
import { RegionTooLargeError, panSNHaplotype } from './util.ts'
import { WalkRowError, canonicalLink } from './walkChunk.ts'
import { departureKind, placementChoice } from './walkDeparture.ts'
import { WalkFile, WalkFormatError } from './walkFile.ts'
import {
  DepartureReadError,
  WALK_BUDGETS,
  WalkReader,
  readInOrder,
} from './walkReader.ts'
import { isPlaceholder, walkNameFilter, walksKept } from './walkRows.ts'

import type { AbridgedPiece } from './walkChunk.ts'
import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkUnit } from './walkFile.ts'
import type { WalkGap } from './walkGaps.ts'
import type { WalkFragment } from './walkRows.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

// the abridged pieces a cut makes a bp row of, those it did not build
const unbuilt = vi.hoisted(() => new Set<AbridgedPiece>())
vi.mock(import('./walkChunk.ts'), async importOriginal => {
  const chunk = await importOriginal()
  return {
    ...chunk,
    bpRow: (piece: AbridgedPiece) => {
      unbuilt.add(piece)
      return chunk.bpRow(piece)
    },
  }
})

function reader(name: string) {
  const path = `${DATA}${name}.walks.gz`
  return WalkReader.fromFilehandles({
    data: new LocalFile(path),
    directory: new LocalFile(`${path}.idx`),
  })
}

// The units of `name`'s sequence each read in `run` fetches, as
// `kind name`, or `departures start` for a departure block, named `.`, and
// how many reads it makes
async function unitsRead(
  name: string,
  refName: string,
  run: (walks: WalkReader) => Promise<unknown>,
) {
  const walks = reader(name)
  const units = await walks.units(refName)
  await walks.header()
  const read = vi.spyOn(LocalFile.prototype, 'read')
  try {
    await run(walks)
    return read.mock.calls.map(([length, position = 0]) =>
      units
        .filter(u => u.offset >= position && u.offset < position + length)
        .map(u => `${u.kind} ${u.kind === 'departures' ? u.start : u.name}`),
    )
  } finally {
    read.mockRestore()
  }
}

async function gfa(
  walks: WalkReader,
  start: number,
  end: number,
  haplotypes?: string[],
) {
  const refName = 'GRCh38#0#chr1'
  const { graph, kept, fragments, departures } = await walks.cut(
    refName,
    { start, end },
    { keep: walkNameFilter(haplotypes, panSNHaplotype(refName)) },
  )
  return graphTablesGFA(graph.tables(kept, fragments, departures))
}

const walkNames = (text: string) =>
  text
    .split('\n')
    .filter(l => l.startsWith('W\t'))
    .map(l => l.split('\t').slice(1, 4).join('#'))

const fragment = (f: WalkFragment) => [
  f.name,
  f.hapStart,
  f.hapEnd,
  [...f.ids].map((id, k) => (f.rev[k]! !== 0 ? '<' : '>') + id).join(''),
  f.trimmedBefore,
  f.trimmedAfter,
]

const departure = (d: WalkDeparture) => [
  d.name,
  d.hapStart,
  d.ownBp,
  d.missingBp,
  `${d.leavePos}${d.leaveReversed ? '-' : '+'}`,
  `${d.rejoinPos}${d.rejoinReversed ? '-' : '+'}`,
  departureKind(d),
  d.fromSteps ? 'steps' : 'record',
  d.steps
    ? d.steps.ids.some(isPlaceholder)
      ? d.steps.ids.length
      : [...d.steps.ids]
          .map((id, k) => (d.steps!.rev[k]! !== 0 ? '<' : '>') + id)
          .join('')
    : undefined,
]

// Four walks over a 3.6 kb reference in 1 kb chunks (header.GRCh38.gfa):
//   W GRCh38  0 chr1 >1>2>3>4>5>6
//   W HG002   1 chr1 >1>2>7>4>5>6   7 replaces 3
//   W HG002   2 chr1 >1>2>3>4>8>6   8 replaces 5
//   W HG00097 1 chr1 >1>2>3>4>5>6
test('the header names the format, the reference, the haplotypes, the tile size, the graph and the lines left out, read from byte 0 and no tile', async () => {
  const read = vi.spyOn(LocalFile.prototype, 'read')
  try {
    const header = await reader('header.GRCh38').header()
    expect(header).toMatchObject({
      format: 8,
      tile: 1000,
      references: ['GRCh38'],
      haplotypes: ['HG00097#1', 'HG002#1', 'HG002#2'],
      holes: [],
      graph: '84388ad9f4ab9098',
      leftout: 0,
    })
    expect(read).toHaveBeenCalled()
    expect(read.mock.calls.map(([, position = 0]) => position)).toEqual(
      read.mock.calls.map(() => 0),
    )
  } finally {
    read.mockRestore()
  }
})

// A walk file whose directory's first line is `walks`, with its header
// unit and data as walkFileOf writes them unless `data` is given
function directoryOf(walks: string, data = new Uint8Array()) {
  return WalkReader.fromFilehandles({
    data: new BlobFile(new Blob([data])),
    directory: new BlobFile(new Blob([`${walks}\n#header\t0\t0\n`])),
  })
}

test.each([
  { walks: '#walks', format: 1 },
  ...[2, 3, 4, 5, 6, 7].map(format => ({
    walks: `#walks\tformat:i:${format}\tchunk:i:1000`,
    format,
  })),
])(
  'a walk file whose directory states $walks is refused as format $format, with what to rebuild it with',
  async ({ walks, format }) => {
    const walkFile = directoryOf(walks)
    await expect(walkFile.header()).rejects.toThrow(WalkFormatError)
    await expect(
      walkFile.cut('GRCh38#0#chr1', { start: 0, end: 600 }),
    ).rejects.toThrow(
      `This walk file is format ${format}; @jbrowse/graph-walks reads format 8, so rebuild it with gfa-to-tabix 5.0 or later`,
    )
  },
)

test.each([9, 10])(
  'a format %i walk file is refused, asking for a newer reader',
  async format => {
    await expect(
      directoryOf(`#walks\tformat:i:${format}\ttile:i:1000`).header(),
    ).rejects.toThrow(
      `This walk file is format ${format}; this version of @jbrowse/graph-walks reads format 8, so update the viewer, or @jbrowse/graph-walks, to one that reads format ${format}`,
    )
  },
)

test('a file that is not a directory is refused', async () => {
  await expect(directoryOf('GRCh38#0#chr1\t0\t100').header()).rejects.toThrow(
    "This is not a walk file's directory",
  )
})

test('open finds the directory beside the data, before any query', async () => {
  const opened: string[] = []
  const walks = WalkReader.open(
    `${DATA}header.GRCh38.walks.gz?v=1#top`,
    url => {
      opened.push(url)
      return new LocalFile(url.replace(/[?#].*$/, ''))
    },
  )
  expect(opened).toEqual([
    `${DATA}header.GRCh38.walks.gz?v=1#top`,
    `${DATA}header.GRCh38.walks.gz.idx?v=1#top`,
  ])
  expect((await walks.header()).tile).toBe(1000)
})

// header.GRCh38's rows, with its directory changed by `edit`
async function edited(edit: (directory: string) => string) {
  const path = `${DATA}header.GRCh38.walks.gz`
  const directory = await new LocalFile(`${path}.idx`).readFile('utf8')
  return WalkReader.fromFilehandles({
    data: new LocalFile(path),
    directory: new BlobFile(new Blob([edit(directory)])),
  })
}

const region = { start: 0, end: 600 }
const entryPoints: [string, (walks: WalkReader) => Promise<unknown>][] = [
  ['header', walks => walks.header()],
  ['referenceSequences', walks => walks.referenceSequences()],
  ['sequenceFor', walks => walks.sequenceFor('GRCh38', 'chr1')],
  ['units', walks => walks.units('GRCh38#0#chr1')],
  ['sequenceLength', walks => walks.sequenceLength('GRCh38#0#chr1')],
  ['bytesFor', walks => walks.bytesFor('GRCh38#0#chr1', region)],
  ['referenceNodes', walks => walks.referenceNodes('GRCh38#0#chr1', region)],
  ['cut', walks => walks.cut('GRCh38#0#chr1', region)],
]

test.each(entryPoints)(
  'a directory and data that disagree are refused by %s',
  async (_, run) => {
    await expect(
      run(await edited(d => d.replace(/\tbytes:i:\d+/, '\tbytes:i:12'))),
    ).rejects.toThrow(
      /directory does not match its data: the directory states 12 bytes, the data holds \d+/,
    )
    await expect(
      run(await edited(d => d.replace(/\tbytes:i:\d+/, ''))),
    ).rejects.toThrow(
      /directory does not match its data: the directory states no bytes, the data holds \d+/,
    )
    await expect(
      run(await edited(d => d.replace('stay:i:311', 'stay:i:0'))),
    ).rejects.toThrow(/directory does not match its data/)
    await expect(run(await edited(d => d))).resolves.toBeDefined()
  },
)

test('a cut keeps the reference and the haplotypes asked for', async () => {
  const names = walkNames(
    await gfa(reader('header.GRCh38'), 1000, 2600, ['HG002#2']),
  )
  expect(names).toEqual(['GRCh38#0#chr1', 'HG002#2#chr1'])
})

test('a cut lists each departure over the window, with the own steps it holds', async () => {
  const { fragments, departures } = await reader('header.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 1000, end: 2600 },
    { context: 0 },
  )
  expect(fragments.map(fragment)).toEqual([
    ['GRCh38#0#chr1', 600, 3000, '>2>3>4>5', 600, 0],
    ['HG00097#1#chr1', 600, 3000, '>2>3>4>5', 600, 0],
    ['HG002#1#chr1', 600, 2950, '>2>7>4>5', 600, 0],
    ['HG002#2#chr1', 600, 2400, '>2>3>4', 600, 20],
  ])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 1200, 550, 0, '1200+', '1800+', 'forward', 'record', '>7'],
    ['HG002#2#chr1', 2400, 20, 0, '2400+', '3000+', 'forward', 'record', '>8'],
  ])
})

// An unchopped graph whose node 1 spans 0-2600, in tile 0 alone:
//   W GRCh38 0 chr1 >1>2>4
//   W HG002  1 chr1 >1>3>4   3 replaces 2
test('a node longer than a tile is read from the tile whose extent holds it', async () => {
  const reads = await unitsRead('maxnode.GRCh38', 'GRCh38#0#chr1', walks =>
    walks.cut('GRCh38#0#chr1', { start: 2100, end: 2500 }, { context: 0 }),
  )
  expect(reads.flat()).toContain('tile 0')
  const text = await gfa(reader('maxnode.GRCh38'), 2100, 2500)
  expect(text).toMatch(/^S\t1\t/m)
  expect(text).toMatch(/^W\tGRCh38\t0\tchr1\t0\t3800\t>1>2>4$/m)
  expect(walkNames(text)).toEqual(['GRCh38#0#chr1', 'HG002#1#chr1'])
})

// Tile 0's extent is 0-1200, as node 2 runs past it; tile 1000's 1200-2400.
// HG002#1 rejoins at 1800 and HG002#2 leaves at 2400, the two ends of node 4,
// which either window draws whole
test.each([
  {
    window: { start: 2100, end: 2600 },
    read: [['tile 1000'], ['tile 2000', 'departures 2000']],
    steps: [undefined, '>8'],
  },
  {
    window: { start: 1700, end: 2300 },
    read: [['tile 1000', 'departures 1000'], ['departures 2000'], ['tile 0']],
    steps: ['>7', undefined],
  },
])(
  'a cut reads the tiles and departure blocks its window overlaps, in one read where they touch, and a tile for the own steps of a record that leaves or rejoins on it alone: $window',
  async ({ window, read, steps }) => {
    let departures: WalkDeparture[] = []
    const reads = await unitsRead(
      'header.GRCh38',
      'GRCh38#0#chr1',
      async walks => {
        departures = (await walks.cut('GRCh38#0#chr1', window, { context: 0 }))
          .departures
      },
    )
    expect(reads).toEqual(read)
    expect(departures.map(departure)).toEqual([
      [
        'HG002#1#chr1',
        1200,
        550,
        0,
        '1200+',
        '1800+',
        'forward',
        'record',
        steps[0],
      ],
      [
        'HG002#2#chr1',
        2400,
        20,
        0,
        '2400+',
        '3000+',
        'forward',
        'record',
        steps[1],
      ],
    ])
  },
)

// A reference that passes nodes twice, in 1 kb chunks:
//   W GRCh38 0 chr1 >1>2>3>4>5>6>5>6>7>8>7>8>9>10
//   W HG002  1 chr1 >1>2>3>4>5>6>5>6>7>8>7>11>9>10   11 replaces the second 8
// 5 and 6 are 250 bp, both passes in chunk 2000; the rest are 500 bp, so 7
// and 8 pass in chunk 3000 and again in chunk 4000
test.each([
  [2500, 3000, '5@2500 6@2750', '>5>6'],
  [4000, 5000, '7@4000 8@4500', '>7>8'],
])(
  'a node the reference passes twice is placed at its pass over the window %i-%i',
  async (start, end, placed, steps) => {
    const { graph, kept, fragments, departures } = await reader(
      'revisit.GRCh38',
    ).cut('GRCh38#0#chr1', { start, end }, { context: 0 })
    const onWindow = [...graph.nodes]
      .filter(([, node]) => node.onWindow)
      .map(([id, node]) => `${id}@${node.start}`)
    expect(onWindow.join(' ')).toBe(placed)
    expect(graphTablesGFA(graph.tables(kept, fragments, departures))).toMatch(
      new RegExp(`^W\\tGRCh38\\t0\\tchr1\\t\\d+\\t\\d+\\t${steps}$`, 'm'),
    )
  },
)

// graph-walks 4's two-argument call must not compile, or a caller that
// missed the bump would draw no departure at all
test('tables takes the departures, so the call without them is a type error', async () => {
  const { graph, kept, fragments, departures } = await reader(
    'skip.GRCh38',
  ).cut('GRCh38#0#chr1', { start: 31, end: 60 }, { context: 0 })
  // @ts-expect-error departures are required
  expect(() => graph.tables(kept, fragments)).toThrow(TypeError)
  expect(graph.tables(kept, fragments, departures).departures.names).toEqual([
    'HG002#1#ctg',
  ])
})

test('a step on a second pass over the window departs from that pass', async () => {
  const { fragments, departures } = await reader('revisit.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 4000, end: 5000 },
    { context: 0 },
  )
  expect(fragments.map(fragment)).toEqual([
    ['GRCh38#0#chr1', 4000, 5000, '>7>8', 0, 0],
    ['HG002#1#chr1', 4000, 4500, '>7', 0, 500],
  ])
  expect(departures.map(departure)).toEqual([
    [
      'HG002#1#chr1',
      4500,
      500,
      0,
      '4500+',
      '5000+',
      'forward',
      'record',
      '>11',
    ],
  ])
})

// A reference that passes 7 and 8, 250 bp each, in the chunk before the
// window and again on it, with 9 and 10 between; the rest are 500 bp:
//   W GRCh38 0 chr1 >1>2>3>4>5>6>7>8>9>10>7>8>11>12
//   W HG002  1 chr1 >1>2>3>4>5>6>7>8>9>10>7>13>11>12   13 replaces the second 8
test("a cut leaves out the reference's other pass over the window's nodes", async () => {
  const { graph, kept, fragments, departures } = await reader(
    'farpass.GRCh38',
  ).cut('GRCh38#0#chr1', { start: 4500, end: 5000 }, { context: 0 })
  const text = graphTablesGFA(graph.tables(kept, fragments, departures))
  expect(text).toMatch(/^W\tGRCh38\t0\tchr1\t4500\t5000\t>7>8$/m)
  expect(text).toMatch(/^W\tHG002\t1\tchr1\t4500\t4750\t>7$/m)
  expect(text).not.toMatch(/^S\t(9|10)\t/m)
  const both = await reader('farpass.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 3900, end: 5000 },
    { context: 0 },
  )
  expect(
    graphTablesGFA(
      both.graph.tables(both.kept, both.fragments, both.departures),
    ),
  ).toMatch(/^S\t7\t\*\tLN:i:250\tSN:Z:GRCh38#0#chr1\tSO:i:4500\tSR:i:0$/m)
  expect(departures.map(departure)).toEqual([
    [
      'HG002#1#chr1',
      4750,
      250,
      0,
      '4750+',
      '5000+',
      'forward',
      'record',
      '>13',
    ],
  ])
})

// Twelve 500 bp reference nodes in 1 kb chunks (follow.GRCh38.gfa): HG002#1
// leaves chunk 1000 for 20 and 21, 2 kb each, and rejoins at 9, in chunk
// 4000; HG003#1's contig ends on 22
test('a walk leaving the window is cut there, and its record says where it rejoins', async () => {
  let cut: Awaited<ReturnType<WalkReader['cut']>> | undefined
  const reads = await unitsRead(
    'follow.GRCh38',
    'GRCh38#0#chr1',
    async walks => {
      cut = await walks.cut(
        'GRCh38#0#chr1',
        { start: 1000, end: 1500 },
        { context: 0 },
      )
    },
  )
  expect(reads).toEqual([['tile 1000', 'departures 1000']])
  const { fragments, departures } = cut!
  expect(fragments.map(fragment)).toEqual([
    ['GRCh38#0#chr1', 1000, 1500, '>3', 0, 500],
    ['HG002#1#chr1', 1000, 1500, '>3', 0, 4000],
    ['HG002#2#chr1', 1000, 1500, '>3', 0, 500],
    ['HG003#1#chr1', 1000, 1500, '>3', 0, 300],
  ])
  expect(departures.map(departure)).toEqual([
    [
      'HG002#1#chr1',
      1500,
      4000,
      0,
      '1500+',
      '4000+',
      'forward',
      'record',
      '>20>21',
    ],
  ])
})

test('the same graph in 10 kb chunks computes the departure from the steps of one piece', async () => {
  const { departures } = await reader('rejoin.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 4000, end: 4200 },
    { context: 0 },
  )
  expect(departures.map(departure)).toEqual([
    [
      'HG002#1#chr1',
      1500,
      4000,
      0,
      '1500+',
      '4000+',
      'forward',
      'steps',
      '>20>21',
    ],
  ])
})

// Twenty 500 bp reference nodes in 1 kb chunks (records.GRCh38.gfa):
//   W GRCh38 0 chr1 0    10000 >1>2>...>20
//   W HG002  1 chr1 0    2000  >1>2>3>4
//   W HG002  1 chr1 2500 3200  >30          visits no reference node
//   W HG002  1 chr1 4000 7000  >7>8>...>12
//   W HG003  1 chr1 0    10000 >1>2>3>19>20>4>5>...>18
const records = async (start: number, end: number, name: string) => {
  const { fragments, departures } = await reader('records.GRCh38').cut(
    'GRCh38#0#chr1',
    { start, end },
    { context: 0 },
  )
  return {
    fragments: fragments.filter(f => f.name === name).map(fragment),
    departures: departures.filter(d => d.name === name).map(departure),
  }
}

test("a departure across a contig's W lines counts the bp the graph lacks as missing", async () => {
  expect(await records(1000, 2000, 'HG002#1#chr1')).toEqual({
    fragments: [['HG002#1#chr1', 1000, 2000, '>3>4', 0, 0]],
    departures: [
      [
        'HG002#1#chr1',
        2000,
        2000,
        1300,
        '2000+',
        '3000+',
        'forward',
        'record',
        undefined,
      ],
    ],
  })
  expect(await records(3000, 3600, 'HG002#1#chr1')).toEqual({
    fragments: [['HG002#1#chr1', 4000, 5000, '>7>8', 0, 0]],
    departures: [
      [
        'HG002#1#chr1',
        2000,
        2000,
        1300,
        '2000+',
        '3000+',
        'forward',
        'record',
        undefined,
      ],
    ],
  })
})

test('a walk jumping away for eight chunks and back has a record each way', async () => {
  expect(await records(1000, 2000, 'HG003#1#chr1')).toEqual({
    fragments: [
      ['HG003#1#chr1', 1000, 1500, '>3', 0, 0],
      ['HG003#1#chr1', 2500, 3000, '>4', 0, 0],
    ],
    departures: [
      [
        'HG003#1#chr1',
        1500,
        0,
        0,
        '1500+',
        '9000+',
        'forward',
        'record',
        undefined,
      ],
      [
        'HG003#1#chr1',
        2500,
        0,
        0,
        '10000+',
        '1500+',
        'back',
        'record',
        undefined,
      ],
    ],
  })
})

test('a cut on a signal already aborted reads no row', async () => {
  const walks = reader('header.GRCh38')
  await walks.header()
  await walks.units('GRCh38#0#chr1')
  const read = vi.spyOn(LocalFile.prototype, 'read')
  const aborted = new AbortController()
  aborted.abort(new Error('panned away'))
  await expect(
    walks.cut(
      'GRCh38#0#chr1',
      { start: 1000, end: 2600 },
      { signal: aborted.signal },
    ),
  ).rejects.toThrow('panned away')
  expect(read).not.toHaveBeenCalled()
  read.mockRestore()
})

test.each([
  [
    'cut',
    (walks: WalkReader, r: { start: number; end: number }) =>
      walks.cut('GRCh38#0#chr1', r),
  ],
  [
    'referenceNodes',
    (walks: WalkReader, r: { start: number; end: number }) =>
      walks.referenceNodes('GRCh38#0#chr1', r),
  ],
  [
    'bytesFor',
    (walks: WalkReader, r: { start: number; end: number }) =>
      walks.bytesFor('GRCh38#0#chr1', r),
  ],
])('a range that is not finite queries nothing by %s', async (_, run) => {
  await expect(
    run(reader('header.GRCh38'), { start: NaN, end: 2600 }),
  ).rejects.toThrow("is not a finite range to query the graph's index for")
  await expect(
    run(reader('header.GRCh38'), { start: 0, end: Infinity }),
  ).rejects.toThrow("is not a finite range to query the graph's index for")
})

test('a cut over its step budget says how far to zoom in', async () => {
  const error = await reader('header.GRCh38')
    .cut('GRCh38#0#chr1', { start: 0, end: 3600 }, { stepBudget: 2 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(RegionTooLargeError)
  expect((error as RegionTooLargeError).regionTooLarge).toBe(true)
})

test('a cut over its departure budget fails without a zoom-in notice', async () => {
  const error = await reader('span.GRCh38')
    .cut('GRCh38#0#chr1', { start: 50, end: 55 }, { departureByteBudget: 1 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(DepartureReadError)
  expect(error).not.toBeInstanceOf(RegionTooLargeError)
  expect((error as DepartureReadError).message).toMatch(
    /^Too many departures from the reference over this window to read \([\d.]+ MB against a budget of 0\.0 MB\)/,
  )
  expect((error as DepartureReadError).limit).toBe(1)
})

// The chrM graph of HPRC's release 2, from `vg convert -f` of its GBZ and
// gfa-to-tabix --walks --refs GRCh38 --sequences
test("a cut keeps each node's bases when asked", async () => {
  const walks = reader('hprc-chrM.GRCh38')
  const { graph, kept, fragments, departures } = await walks.cut(
    'GRCh38#0#chrM',
    { start: 8000, end: 8300 },
    { sequences: true },
  )
  expect(
    walkNames(graphTablesGFA(graph.tables(kept, fragments, departures))),
  ).toHaveLength(44)
  for (const id of kept) {
    const node = graph.nodes.get(id)!
    expect(graph.sequences.get(id)).toMatch(/^[ACGTN]+$/)
    expect(graph.sequences.get(id)).toHaveLength(node.end - node.start)
  }
  expect(await walks.referenceSequences()).toEqual(['GRCh38#0#chrM'])
})

test('a contig names its sequence on the reference sample', async () => {
  const walks = reader('hprc-chrM.GRCh38')
  expect(await walks.sequenceFor('GRCh38', 'chrM')).toBe('GRCh38#0#chrM')
  expect(await walks.sequenceFor('hg38', 'GRCh38#0#chrM')).toBe('GRCh38#0#chrM')
  expect(await walks.sequenceFor('GRCh38', 'chr1')).toBeUndefined()
})

test("a region cut names the file's sequence for a contig, and keeps the haplotypes asked for", async () => {
  const cut = await reader('header.GRCh38').cutRegion(
    { refName: 'chr1', start: 1000, end: 2600 },
    { haplotypes: ['HG002#2'] },
  )
  expect(cut.refName).toBe('GRCh38#0#chr1')
  expect(cut.sample).toBe('GRCh38')
  expect(walkNames(graphTablesGFA(cut.tables))).toEqual([
    'GRCh38#0#chr1',
    'HG002#2#chr1',
  ])
  expect(cut.tables.departures.names).toEqual(['HG002#2#chr1'])
  const whole = await reader('header.GRCh38').cutRegion({
    refName: 'GRCh38#0#chr1',
    start: 1000,
    end: 2600,
  })
  expect(graphTablesGFA(whole.tables)).toBe(
    await gfa(reader('header.GRCh38'), 1000, 2600),
  )
})

test('a region cut on a sequence the file lacks names the ones it holds', async () => {
  await expect(
    reader('header.GRCh38').cutRegion({ refName: 'chr2', start: 0, end: 100 }),
  ).rejects.toThrow(
    'The walk file has no sequence chr2 on GRCh38; it holds GRCh38#0#chr1',
  )
})

test('a region cut over its node limit says how far to zoom in', async () => {
  const error = await reader('header.GRCh38')
    .cutRegion({ refName: 'chr1', start: 0, end: 3600 }, { nodeLimit: 4 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(RegionTooLargeError)
  expect(error).toMatchObject({ limit: 4, windowBp: 3600 })
  expect((error as Error).message).toMatch(
    /^Zoom in to about [\d.]+K?bp, or choose fewer haplotypes, to see the graph$/,
  )
})

// A walk file from rows as walkFileOf takes them, its units written in order
// or reversed, which notes the nominal start of each tile it serves
function inline(lines: string[], reversed = false, served: number[] = []) {
  const { data, directory, units } = walkFileOf(lines, { reversed })
  const read = data.read.bind(data)
  data.read = (length, position, opts) => {
    for (const unit of units) {
      if (
        unit.kind === 'tile' &&
        unit.offset >= position &&
        unit.offset < position + length
      ) {
        served.push(+unit.name)
      }
    }
    return read(length, position, opts)
  }
  return WalkReader.fromFilehandles({ data, directory })
}

// Eight 4 bp reference nodes in 16 bp chunks:
//   W R  0 chr >1>2>3>4>5>6>7>8
//   W H1 1 c   >1>2>6>3>4>5>6>7>8
//   W H2 1 c   >5>6>2>7>8
// `pieces` files each visit to the other chunk as a piece of its own. `stays`
// keeps the visit in its piece as a route interior, as gfa-to-tabix's --stay
// does on a reference that passes 2 and 6 twice, so chunk 0's route steps on
// 6, which chunk 16's ref row places, and chunk 16's on 2.
const tinyHeader = [
  '#walks format:i:8 tile:i:16',
  '#reference R',
  '#haplotype H1#1',
  '#haplotype H2#1',
]
const pieces = [
  ...tinyHeader,
  'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 .',
  'R#0#chr 0 walk H1#1#c 0 0 0 2 2 .',
  'R#0#chr 0 walk H1#1#c 0 12 2 2 6 .',
  'R#0#chr 0 walk H2#1#c 0 8 1 1 4 .',
  'R#0#chr 16 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 .',
  'R#0#chr 16 walk H1#1#c 0 8 1 1 12 .',
  'R#0#chr 16 walk H1#1#c 0 20 3 4 10 .',
  'R#0#chr 16 walk H2#1#c 0 0 0 2 10 .',
  'R#0#chr 16 walk H2#1#c 0 12 2 2 14 .',
]
const stays = [
  ...tinyHeader,
  'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 .',
  'R#0#chr 0 nodes 6 4',
  'R#0#chr 0 route 4 8 -6',
  'R#0#chr 0 walk H1#1#c 0 0 0 5 2 0',
  'R#0#chr 16 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 .',
  'R#0#chr 16 nodes 2 4',
  'R#0#chr 16 route 12 -8 10',
  'R#0#chr 16 walk H1#1#c 0 20 1 4 10 .',
  'R#0#chr 16 walk H2#1#c 0 0 0 5 10 0',
]

// Each window reads both chunks, the first as the chunk before it. In
// `stays` a route's step on a reference node is the own step of a departure
// inside one piece, so the cut keeps the node, at its reference position
describe.each([
  { start: 4, end: 20, haplotypes: ['H1'], node: 6, at: 20, step: '>6' },
  { start: 12, end: 32, haplotypes: ['H2'], node: 2, at: 4, step: '>2' },
  { start: 0, end: 32, haplotypes: undefined, node: 6, at: 20, step: '>6' },
])(
  'a reference node a route steps on, cut $start-$end for $haplotypes',
  ({ start, end, haplotypes, node, at, step }) => {
    test('keeps its reference position, in whichever order the chunks come', async () => {
      const cut = (reversed?: boolean) =>
        inline(stays, reversed)
          .cutRegion({ refName: 'chr', start, end }, { haplotypes, context: 0 })
          .then(({ tables }) => ({
            text: graphTablesGFA(tables),
            own: [...tables.departures.steps].map(i => tables.nodes.names[i]),
          }))
      const { text, own } = await cut()
      expect(text).toMatch(
        new RegExp(
          `^S\\t${node}\\t\\*\\tLN:i:4\\tSN:Z:R#0#chr\\tSO:i:${at}\\tSR:i:0$`,
          'm',
        ),
      )
      expect(own).toContain(String(node))
      expect(step).toBe(`>${node}`)
      expect(await cut(true)).toEqual({ text, own })
    })
  },
)

// The pieces join into one run each, so an uncontained cut keeps each walk
// from its first step on the window to its last, the visit to the other
// chunk included, and the jumps between the pieces come from the records
test('the same walks in pieces of their own cut whole, with the jumps as records', async () => {
  const departureRecords = [
    ...pieces,
    'R#0#chr 4 24 departure H1#1#c 8 0 0 8 + 1 20 + 1 .',
    'R#0#chr 8 24 departure H1#1#c 12 0 0 24 + 1 8 + 1 .',
    'R#0#chr 4 24 departure H2#1#c 8 0 0 24 + 1 4 + 1 .',
    'R#0#chr 4 28 departure H2#1#c 12 0 0 8 + 1 24 + 1 .',
  ]
  const { tables } = await inline(departureRecords).cutRegion(
    { refName: 'chr', start: 4, end: 20 },
    { context: 0 },
  )
  const text = graphTablesGFA(tables)
  expect(text.split('\n').filter(l => l.startsWith('W'))).toEqual([
    'W\tR\t0\tchr\t4\t20\t>2>3>4>5',
    'W\tH1\t1\tc\t4\t24\t>2>6>3>4>5',
    'W\tH2\t1\tc\t0\t12\t>5>6>2',
  ])
  expect(
    [...tables.departures.leavePos].map(
      (pos, k) =>
        `${tables.departures.names[k]} ${pos}->${tables.departures.rejoinPos[k]}`,
    ),
  ).toEqual(['H1#1#c 8->20', 'H1#1#c 24->8', 'H2#1#c 24->4', 'H2#1#c 8->24'])
})

test("a route's own steps inside one piece are a departure computed from the steps", async () => {
  const { departures } = await inline(stays).cut(
    'R#0#chr',
    { start: 4, end: 20 },
    { context: 0 },
  )
  expect(departures.map(departure)).toEqual([
    ['H1#1#c', 8, 4, 0, '8+', '8+', 'adjacent', 'steps', '>6'],
  ])
})

// Four 4 bp reference nodes, the second of which the reference passes twice,
// its other pass elsewhere; H#1 runs 9 between 2 and 3:
//   R  >1>2>3>4, passes 1,2,1,1
//   H  >1>2>9>3>4
test("a departure computed from steps takes its ends' pass counts from the ref rows", async () => {
  const { graph, kept, fragments, departures } = await inline([
    '#walks format:i:8 tile:i:16 revisited:i:1',
    '#reference R',
    '#haplotype H#1',
    'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 1,2,1,1',
    'R#0#chr 0 nodes 9 4',
    'R#0#chr 0 route 4 14 -12',
    'R#0#chr 0 walk H#1#h 0 0 0 5 2 0',
  ]).cut('R#0#chr', { start: 0, end: 16 }, { context: 0 })
  expect(
    departures.map(d => [
      d.leavePos,
      d.leavePasses,
      d.rejoinPos,
      d.rejoinPasses,
      d.fromSteps,
    ]),
  ).toEqual([[8, 2, 8, 1, true]])
  expect(departures.map(placementChoice)).toEqual([true])
  const tables = graph.tables(kept, fragments, departures)
  expect(
    tables.nodes.names.map((name, i) => `${name}:${tables.nodes.passes![i]}`),
  ).toEqual(['9:0', '1:1', '2:2', '3:1', '4:1'])
  expect([...tables.departures.leavePasses]).toEqual([2])
})

test('a cut keeps to WALK_BUDGETS unless its caller lifts them', async () => {
  const window = { start: 16, end: 32 }
  const defaults = { ...WALK_BUDGETS }
  try {
    WALK_BUDGETS.departureByteBudget = 1
    await expect(
      inline(copyJump).cut('R#0#chr', window),
    ).rejects.toBeInstanceOf(DepartureReadError)
    const lifted = await inline(copyJump).cut('R#0#chr', window, {
      departureByteBudget: Infinity,
    })
    expect(lifted.departures).toHaveLength(2)
    WALK_BUDGETS.departureByteBudget = defaults.departureByteBudget
    WALK_BUDGETS.byteBudget = 1
    await expect(
      inline(copyJump).cut('R#0#chr', window),
    ).rejects.toBeInstanceOf(RegionTooLargeError)
  } finally {
    Object.assign(WALK_BUDGETS, defaults)
  }
})

test('a ref row whose passes do not count its steps is refused', async () => {
  const error = await inline([
    '#walks format:i:8 tile:i:16',
    '#reference R',
    'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4',
  ])
    .cut('R#0#chr', { start: 0, end: 16 }, { context: 0 })
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(WalkRowError)
})

// A handle token is 2 × the id's difference from the one before, so the
// tokens 2 and 2 × (id − 1) are the ids 1 and `id`
describe('a node id at or past 2^31, the format’s limit, or below 0', () => {
  const tileOf = (...rows: string[]) =>
    inline([
      '#walks format:i:8 tile:i:16',
      '#reference R',
      'R#0#chr 0 ref R#0#chr 0 0 0 2 2,2 4,4 .',
      ...rows,
    ])
  const refOf = (id: number) =>
    inline([
      '#walks format:i:8 tile:i:16',
      '#reference R',
      `R#0#chr 0 ref R#0#chr 0 0 0 2 2,${2 * (id - 1)} 4,4 .`,
    ])
  const window = { start: 0, end: 16 }

  test('2^31 − 1 is read as itself, which is no placeholder', async () => {
    const cut = await refOf(2 ** 31 - 1).cut('R#0#chr', window, { context: 0 })
    expect([...cut.kept].sort((a, b) => a - b)).toEqual([1, 2 ** 31 - 1])
    expect([...cut.kept].some(isPlaceholder)).toBe(false)
  })

  test.each([
    { where: 'a ref row', file: refOf(2 ** 31), id: 2 ** 31 },
    { where: 'a ref row', file: refOf(2 ** 32), id: 2 ** 32 },
    { where: 'a ref row', file: refOf(-1), id: -1 },
    {
      where: 'a nodes row',
      file: tileOf(`R#0#chr 0 nodes 3,${2 ** 31 - 3} 4,4`),
      id: 2 ** 31,
    },
    {
      where: 'a route row',
      file: tileOf(`R#0#chr 0 route 4 ${2 ** 32 - 4} .`),
      id: 2 ** 31,
    },
  ])('is refused in $where: $id', async ({ file, id }) => {
    const error = await file
      .cut('R#0#chr', window, { context: 0 })
      .catch((e: unknown) => e)
    expect(error).toBeInstanceOf(WalkRowError)
    expect(String(error)).toContain(`node id ${id}`)
    expect(String(error)).toContain('2^31')
  })

  test('is a WalkRowError a consumer of the package can catch by type', async () => {
    const { WalkReader: Exported, WalkRowError: Refusal } =
      await import('./index.ts')
    const { data, directory } = walkFileOf([
      '#walks format:i:8 tile:i:16',
      '#reference R',
      `R#0#chr 0 ref R#0#chr 0 0 0 2 2,${2 ** 32} 4,4 .`,
    ])
    const error = await Exported.fromFilehandles({ data, directory })
      .cut('R#0#chr', window, { context: 0 })
      .catch((e: unknown) => e)
    expect(error).toBeInstanceOf(Refusal)
  })
})

// Twelve 4 bp reference nodes in 16 bp chunks; H#1 goes from 7 to 11, in
// chunk 32, and back to 6:
//   W R  0 chr >1>2>3>4>5>6>7>8>9>10>11>12
//   W H  1 h   >5>6>7>11>6>7>8
const copyJump = [
  '#walks format:i:8 tile:i:16',
  '#reference R',
  '#haplotype H#1',
  'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 .',
  'R#0#chr 16 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 .',
  'R#0#chr 16 walk H#1#h 0 0 0 3 10 .',
  'R#0#chr 16 walk H#1#h 0 16 2 3 12 .',
  'R#0#chr 32 ref R#0#chr 0 32 2 4 18,2,2,2 4,4,4,4 .',
  'R#0#chr 32 walk H#1#h 0 12 1 1 22 .',
  'R#0#chr 24 44 departure H#1#h 12 0 0 28 + 1 40 + 1 .',
  'R#0#chr 20 44 departure H#1#h 16 0 0 44 + 1 20 + 1 .',
]
const copyJumpCut = async (
  start: number,
  end: number,
  contained: boolean,
  walks = inline(copyJump),
) => {
  const { graph, kept, fragments, departures } = await walks.cut(
    'R#0#chr',
    { start, end },
    { context: 0, contained, nodeLimit: 5 },
  )
  return {
    text: graphTablesGFA(graph.tables(kept, fragments, departures)),
    walks: fragments.filter(f => f.name === 'H#1#h').map(fragment),
    departures: departures.map(departure),
  }
}

test('a contained cut splits a walk where it steps onto the reference off the window', async () => {
  const detour = await copyJumpCut(20, 36, false)
  expect(detour.walks).toEqual([['H#1#h', 4, 28, '>6>7>11>6>7>8', 4, 0]])
  const departureRecords = [
    ['H#1#h', 12, 0, 0, '28+', '40+', 'forward', 'record', undefined],
    ['H#1#h', 16, 0, 0, '44+', '20+', 'back', 'record', undefined],
  ]
  expect(detour.departures).toEqual(departureRecords)
  expect(detour.text).toMatch(/^S\t11\t/m)
  expect(detour.text).toMatch(/^W\tR\t0\tchr\t20\t36\t>6>7>8>9$/m)
  const contained = await copyJumpCut(20, 36, true)
  expect(contained.walks).toEqual([
    ['H#1#h', 4, 12, '>6>7', 4, 16],
    ['H#1#h', 16, 28, '>6>7>8', 16, 0],
  ])
  expect(contained.departures).toEqual(departureRecords)
  expect(contained.text).not.toMatch(/^S\t11\t/m)
})

test('a cut reads the tiles its window overlaps, and nothing else', async () => {
  const served: number[] = []
  const { walks, departures } = await copyJumpCut(
    20,
    28,
    true,
    inline(copyJump, false, served),
  )
  expect(walks).toEqual([
    ['H#1#h', 4, 12, '>6>7', 4, 0],
    ['H#1#h', 16, 24, '>6>7', 0, 4],
  ])
  expect(departures).toHaveLength(2)
  expect(new Set(served)).toEqual(new Set([16]))
})

// Eight chunks of four 4 bp reference nodes
const eightChunks = [
  '#walks format:i:8 tile:i:16',
  '#reference R',
  ...Array.from(
    { length: 8 },
    (_, c) =>
      `R#0#chr ${16 * c} ref R#0#chr 0 ${16 * c} ${c} 4 ${8 * c + 2},2,2,2 4,4,4,4 .`,
  ),
]

test('a cut over its node limit stops decoding at the tile that puts it over', async () => {
  const error = await inline(eightChunks)
    .cutRegion(
      { refName: 'chr', start: 0, end: 128 },
      { context: 0, nodeLimit: 6 },
    )
    .catch((e: unknown) => e)
  expect(error).toBeInstanceOf(RegionTooLargeError)
  expect(error).toMatchObject({ limit: 6, windowBp: 128, fitsBp: 24 })
  const whole = await inline(eightChunks).cutRegion(
    { refName: 'chr', start: 0, end: 128 },
    { context: 0, nodeLimit: 32 },
  )
  expect(whole.tables.nodes.names).toHaveLength(32)
})

// gfa-to-tabix's walks-two-refs.gfa on R in 10 bp chunks, once as written and
// once with --abridge-bytes 1, which abridges every walk in every chunk
async function twoRefs(name: string, expand?: (name: string) => boolean) {
  const { graph, kept, fragments, departures } = await reader(name).cut(
    'R#0#chr',
    { start: 0, end: 30 },
    { expand },
  )
  return {
    tables: graph.tables(kept, fragments, departures),
    fragments,
    departures,
  }
}

// H#1#h is >1>2 in tile 0, >3>7>4 in tile 10, with 3 bp of its own in 7,
// and >5 in tile 20
test('an unloaded walk cuts to the fragment its loaded one does, a placeholder for its run of own steps, and lists the pieces it did not load', async () => {
  const { tables, fragments } = await twoRefs('abridged.R')
  const loaded = await twoRefs('abridged.R', () => true)
  const steps = (f: WalkFragment) =>
    [...f.ids].map((id, k) =>
      isPlaceholder(id) ? 'placeholder' : `${id}@${f.refStarts[k]}`,
    )
  const h = fragments.filter(f => f.name === 'H#1#h')
  expect(
    h.map(f => [f.hapStart, f.hapEnd, f.trimmedBefore, f.trimmedAfter]),
  ).toEqual([[100, 133, 0, 0]])
  expect(h.map(steps)).toEqual([
    ['1@0', '2@6', '3@12', 'placeholder', '4@18', '5@24'],
  ])
  expect(
    h[0]!.abridged!.map(a => [a.steps, a.detail.refName, a.detail.name]),
  ).toEqual([
    [2, 'R#0#chr', '0:2'],
    [3, 'R#0#chr', '10:0'],
    [1, 'R#0#chr', '20:1'],
  ])
  const whole = loaded.fragments.filter(f => f.name === 'H#1#h')
  expect(whole.map(f => [f.hapStart, f.hapEnd, f.abridged])).toEqual([
    [100, 133, undefined],
  ])
  expect(whole.map(steps)).toEqual([
    ['1@0', '2@6', '3@12', '7@NaN', '4@18', '5@24'],
  ])
  const { nodes, walks } = tables
  const placeholders = [...nodes.placeholders!].flatMap((flag, i) =>
    flag === 1
      ? [
          {
            name: nodes.names[i],
            bp: nodes.lengths[i],
            on: nodes.refNames[nodes.refs[i]!],
            from: nodes.starts[i],
          },
        ]
      : [],
  )
  expect(placeholders).toContainEqual({
    name: 'H#1#h:118-121',
    bp: 3,
    on: 'H#1#h',
    from: 118,
  })
  expect(new Set(placeholders.map(p => p.name)).size).toBe(placeholders.length)
  expect(loaded.tables.nodes.placeholders!.some(flag => flag === 1)).toBe(false)
  walks.names.forEach((name, i) => {
    expect(walks.abridgedSteps![i]).toBe(
      name === 'R#0#chr'
        ? NaN
        : fragments
            .find(f => f.name === name && f.hapStart === walks.starts[i])!
            .abridged!.reduce((sum, a) => sum + a.steps, 0),
    )
  })
  expect(walks.abridgedSteps![walks.names.indexOf('H#1#h')]).toBe(6)
})

// R over four 5 bp nodes in 20 bp tiles, and H abridged from 3 at 10 to 4
// at 15 (abredge.R.gfa):
//   W R 0 chr >1>2>3>4
//   W H 1 h   >3>5>4   5 of 8 bp
test('an abridged piece starting where the drawn stretch ends is not cut', async () => {
  const { fragments } = await reader('abredge.R').cut(
    'R#0#chr',
    { start: 0, end: 10 },
    { context: 0 },
  )
  expect(fragments.map(f => f.name)).toEqual(['R#0#chr'])
})

test('an expanded walk is read from its detail records as the unabridged file cuts it', async () => {
  const all = () => true
  const expanded = await twoRefs('abridged.R', all)
  const whole = await twoRefs('unabridged.R')
  expect(graphTablesGFA(expanded.tables)).toBe(graphTablesGFA(whole.tables))
  expect(expanded.departures.map(departure)).toEqual(
    whole.departures.map(departure),
  )
  expect(whole.departures.map(departure)).toEqual([
    ['D#1#d', 24, 0, 0, '24+', '0+', 'back', 'record', undefined],
    ['H#1#h', 118, 3, 0, '18+', '18+', 'adjacent', 'steps', '>7'],
    ['I#1#i', 6, 0, 0, '6+', '18-', 'switch', 'record', undefined],
    ['I#1#i', 18, 0, 0, '6-', '18+', 'switch', 'record', undefined],
    ['Q#0#chr', 6, 6, 0, '6+', '12+', 'forward', 'record', '>6'],
  ])
})

function gridCut(
  name: string,
  refName: string,
  span: { start: number; end: number },
  mode: { contained: boolean; context: number },
) {
  return reader(name).cut(refName, span, { ...mode, expand: () => true })
}

function expandedFacts(
  { graph, kept, fragments, departures }: Awaited<ReturnType<typeof gridCut>>,
  gaps: WalkGap[],
) {
  const tables = graph.tables(kept, fragments, departures, gaps)
  return {
    gfa: graphTablesGFA(tables),
    tables,
    fragments: fragments.toSorted(
      (a, b) =>
        (a.name < b.name ? -1 : a.name > b.name ? 1 : 0) ||
        a.fragStart - b.fragStart ||
        a.hapStart - b.hapStart,
    ),
    departures,
    gaps,
  }
}

// The README's exception to an abridged file cutting as its twin: an unread
// segment's `before` or `after` beside an abridged piece the cut did not
// build, one ending where the segment starts or starting where it ends, is
// NaN where the twin states the breakpoint
function twinBreakpointsForNaN(
  gaps: WalkGap[],
  twin: WalkGap[],
  skipped: AbridgedPiece[],
) {
  const ending = new Set(
    skipped.map(p => `${p.name}\t${p.hapOffset + p.hapBp}`),
  )
  const starting = new Set(skipped.map(p => `${p.name}\t${p.hapOffset}`))
  return gaps.map((gap, g) => {
    const starts = gap.segments.map((_, s) =>
      gap.segments.slice(0, s).reduce((at, x) => at + x.bp, gap.start),
    )
    return {
      ...gap,
      segments: gap.segments.map((segment, s) => {
        const other = twin[g]?.segments[s]
        const start = starts[s]!
        const end = start + segment.bp
        return segment.kind === 'unread' && other?.kind === 'unread'
          ? {
              ...segment,
              before:
                Number.isNaN(segment.before) &&
                ending.has(`${gap.name}\t${start}`)
                  ? other.before
                  : segment.before,
              after:
                Number.isNaN(segment.after) &&
                starting.has(`${gap.name}\t${end}`)
                  ? other.after
                  : segment.after,
            }
          : segment
      }),
    }
  })
}

// Every window of a grid over each of `name`'s sequences, plain and
// contained, with and without context
async function* gridOf(name: string) {
  const walks = reader(name)
  for (const refName of await walks.referenceSequences()) {
    const length = await walks.sequenceLength(refName)
    const step = Math.ceil(length / 12)
    for (const mode of [
      { contained: false, context: 0 },
      { contained: true, context: 0 },
      { contained: false, context: step },
      { contained: true, context: step },
    ]) {
      for (let start = 0; start < length; start += step) {
        for (let end = start + step; end < length + step; end += step) {
          yield { refName, span: { start, end: Math.min(end, length) }, mode }
        }
      }
    }
  }
}

const where = ({
  refName,
  span,
  mode,
}: {
  refName: string
  span: { start: number; end: number }
  mode: { contained: boolean; context: number }
}) => `${refName}:${span.start}-${span.end} ${JSON.stringify(mode)}`

test.each([
  ['detailtwice.R', 'detailtwicewhole.R'],
  ['sortsfirst.R', 'sortsfirstwhole.R'],
  ['ownends.R', 'ownendswhole.R'],
  ['unbuilt.R', 'unbuiltwhole.R'],
  ['edge.R', 'edgewhole.R'],
  ['adv.R', 'advwhole.R'],
  ['gapnan.R', 'gapnanwhole.R'],
  ['secondline.R', 'secondlinewhole.R'],
])(
  'every window of %s, expanding every walk, cuts as %s, its graph built without abridging, does',
  async (abridged, whole) => {
    const differ: string[] = []
    for await (const at of gridOf(abridged)) {
      const { refName, span, mode } = at
      unbuilt.clear()
      const [a, b] = await Promise.all([
        gridCut(abridged, refName, span, mode),
        gridCut(whole, refName, span, mode),
      ])
      const skipped = [...unbuilt]
      if (
        digestOf(
          expandedFacts(a, twinBreakpointsForNaN(a.gaps, b.gaps, skipped)),
        ) !== digestOf(expandedFacts(b, b.gaps))
      ) {
        differ.push(where(at))
      }
    }
    expect(differ).toEqual([])
  },
)

// gapnan.R in 20 bp tiles: R walks 1 to 6 and H walks 1, 5, 2, 6 and 1 again,
// each 10 bp, abridged to pieces of one step. A cut of 0-10 builds H's pieces
// on 1 and leaves its piece on 2, at 10-20 in the same tile, unbuilt
test('an unread stretch beside an abridged piece the cut did not build has no breakpoint there, where its twin built without abridging states one', async () => {
  const gapsOf = async (name: string) => {
    const { gaps } = await reader(name).cut(
      'R#0#c',
      { start: 0, end: 10 },
      { context: 0, expand: () => true },
    )
    return gaps.map(gap => [
      gap.name,
      gap.start,
      gap.end,
      gap.segments.map(x => [x.kind, x.bp, x.before, x.after]),
    ])
  }
  expect(await gapsOf('gapnan.R')).toEqual([
    [
      'H#1#h',
      10,
      40,
      [
        ['unread', 10, 10, Number.NaN],
        ['read', 10, Number.NaN, Number.NaN],
        ['unread', 10, Number.NaN, 0],
      ],
    ],
  ])
  expect(await gapsOf('gapnanwhole.R')).toEqual([
    [
      'H#1#h',
      10,
      40,
      [
        ['unread', 10, 10, 10],
        ['read', 10, Number.NaN, Number.NaN],
        ['unread', 10, 20, 0],
      ],
    ],
  ])
})

// edge.R in 60 bp tiles: TRAIL walks 7 and 8, its own 28, then 1 and 2. A
// cut of 20-80 reads tile 0 and one of 60-80 does not; neither returns
// TRAIL's step on 1
test.each([
  ['edge.R', 'edgewhole.R'],
  ['adv.R', 'advwhole.R'],
  ['abut.R'],
  ['farrejoin.R'],
])(
  "every window of %s that returns a walk's same steps states the same break in each of its records",
  async (...names) => {
    const differ: string[] = []
    for (const name of names) {
      const breaks = new Map<string, number | undefined>()
      for await (const at of gridOf(name)) {
        const { fragments, departures } = await reader(name).cut(
          at.refName,
          at.span,
          { ...at.mode, expand: () => true },
        )
        for (const d of departures) {
          if (d.steps) {
            const key = JSON.stringify([
              d.name,
              d.hapStart,
              [...d.steps.ids],
              fragments
                .filter(f => f.name === d.name)
                .toSorted((a, b) => a.hapStart - b.hapStart)
                .map(f => [f.hapStart, [...f.ids], [...f.refStarts]]),
            ])
            if (breaks.has(key) && breaks.get(key) !== d.breakAt) {
              differ.push(
                `${name} ${d.name} at ${d.hapStart}: ${breaks.get(key)} and ${d.breakAt} at ${where(at)}`,
              )
            }
            breaks.set(key, d.breakAt)
          }
        }
      }
    }
    expect(differ).toEqual([])
  },
)

async function placedOf(
  name: string,
  node: string,
  start: number,
  end: number,
) {
  const { graph, kept, fragments, departures } = await reader(name).cut(
    'R#0#c',
    { start, end },
    { context: 0, expand: () => true },
  )
  const { nodes } = graph.tables(kept, fragments, departures)
  const k = nodes.names.indexOf(node)
  return [nodes.refNames[nodes.refs[k]!], nodes.starts[k]]
}

// detailtwice.R in 40 bp tiles: H walks 1, its own 20, 2, then 9 and 10 in
// tile 80, then 3, its own 20 again, and 4. Its two pieces in tile 0 are one
// detail unit
test('a node off the reference lies at the first step on it that the cut returns, as the file built without abridging has it', async () => {
  for (const name of ['detailtwice.R', 'detailtwicewhole.R']) {
    expect(await placedOf(name, '20', 20, 40)).toEqual(['H#1#h', 55])
    expect(await placedOf(name, '20', 0, 40)).toEqual(['H#1#h', 10])
  }
})

// sortsfirst.R in one 40 bp tile: A walks 1, its own 20 and 2, Z walks 3, its
// own 20 and 4, each abridged to a detail unit of its own
test('a node off the reference takes no place from a walk the cut returns no step of', async () => {
  for (const name of ['sortsfirst.R', 'sortsfirstwhole.R']) {
    expect(await placedOf(name, '20', 20, 40)).toEqual(['Z#1#z', 10])
    expect(await placedOf(name, '20', 0, 40)).toEqual(['A#1#a', 10])
  }
})

// secondline.R in 20 bp tiles: R walks 1 to 4, each 10 bp. H's first W line
// walks 1, 2 and its own 50 of 5 bp to 25, its second starts at 40 with its
// own 51 of 5 bp, then 3 and 4. Its record from 2 to 3 holds 50 and 51,
// breaking between them, with the 15 bp between the lines missing
test("a departure's steps that lead to its rejoin lie where its walk has them, past the sequence its record is missing", async () => {
  for (const name of ['secondline.R', 'secondlinewhole.R']) {
    expect(await placedOf(name, '50', 10, 30)).toEqual(['H#1#h', 20])
    expect(await placedOf(name, '51', 10, 30)).toEqual(['H#1#h', 40])
  }
})

test('the reference sequences leave out detail records and departures', async () => {
  expect(await reader('abridged.R').referenceSequences()).toEqual(['R#0#chr'])
})

// abridged.R, each read of a detail unit made through `around`
async function heldDetails(
  around: <T>(
    read: () => Promise<T>,
    detail: { start: number; signal?: AbortSignal },
  ) => Promise<T>,
) {
  const path = `${DATA}abridged.R.walks.gz`
  const units = await reader('abridged.R').units('R#0#chr')
  const details = new Set(
    units.filter(u => u.kind === 'detail').map(u => u.offset),
  )
  const data = new LocalFile(path)
  const read = data.read.bind(data)
  data.read = (length, position = 0, opts = {}) =>
    details.has(position)
      ? around(() => read(length, position, opts), {
          start: position,
          signal: opts.signal,
        })
      : read(length, position, opts)
  return WalkReader.fromFilehandles({
    data,
    directory: new LocalFile(`${path}.idx`),
  })
}

async function expandedCut(walks: WalkReader, signal?: AbortSignal) {
  const { graph, kept, fragments, departures } = await walks.cut(
    'R#0#chr',
    { start: 0, end: 30 },
    { expand: () => true, signal },
  )
  return {
    gfa: graphTablesGFA(graph.tables(kept, fragments, departures)),
    fragments,
  }
}

test('detail records read at once, finishing in reverse, cut as reads one at a time do', async () => {
  const abridgedPieces = (await twoRefs('abridged.R')).fragments.flatMap(f =>
    (f.abridged ?? []).map(a => `${a.detail.refName} ${a.detail.name}`),
  )
  const details = new Set(abridgedPieces).size
  expect(details).toBeLessThan(abridgedPieces.length)
  let queue = Promise.resolve()
  let reading = 0
  let most = 0
  const serial = await expandedCut(
    await heldDetails(read => {
      const turn = queue.then(async () => {
        most = Math.max(most, ++reading)
        const bytes = await read()
        reading--
        return bytes
      })
      queue = turn.then(() => undefined)
      return turn
    }),
  )
  expect(most).toBe(1)
  const waiting: (() => void)[] = []
  const finished: number[] = []
  const concurrent = await expandedCut(
    await heldDetails(async (read, { start }) => {
      const bytes = await read()
      await new Promise<void>(resolve => {
        waiting.push(() => {
          finished.push(start)
          resolve()
        })
        if (waiting.length === details) {
          for (const release of waiting.toReversed()) {
            release()
          }
        }
      })
      return bytes
    }),
  )
  expect(finished).toHaveLength(details)
  expect(finished).not.toEqual(finished.toSorted((a, b) => a - b))
  expect(concurrent).toEqual(serial)
})

test('aborting an expanding cut stops its detail reads', async () => {
  const controller = new AbortController()
  const signals: AbortSignal[] = []
  const walks = await heldDetails(
    (_read, { signal }) =>
      new Promise<never>((_resolve, reject) => {
        signals.push(signal!)
        signal!.addEventListener('abort', () => {
          reject(signal!.reason as Error)
        })
        if (signals.length === 3) {
          controller.abort(new Error('panned away'))
        }
      }),
  )
  await expect(expandedCut(walks, controller.signal)).rejects.toThrow(
    'panned away',
  )
  expect(signals).toHaveLength(3)
  expect(signals.every(signal => signal.aborted)).toBe(true)
})

test('readInOrder keeps at most its limit reading, and uses results in order', async () => {
  let reading = 0
  let most = 0
  const used: number[] = []
  await readInOrder(
    [5, 1, 4, 2, 3, 0],
    3,
    new AbortController().signal,
    async (ms: number) => {
      most = Math.max(most, ++reading)
      await new Promise(resolve => setTimeout(resolve, ms))
      reading--
      return ms * 10
    },
    (_ms, result) => used.push(result),
  )
  expect(most).toBe(3)
  expect(used).toEqual([50, 10, 40, 20, 30, 0])
})

// gfa-to-tabix --walks --chunk 10 of >1>2>9>10 against a reference >1 to >10
// of 10 bp nodes: HG002#1#ctg leaves the reference after 20 and rejoins at 80,
// and its record is filed over the stretch between
// The drawn stretch is the window and its context, so a cut keeps the same
// reference whether the context is asked for or the window widened by it, on
// both sides and past the chunks under the window itself
test('a cut reads the chunks under its context, not only under its window', async () => {
  const walks = reader('skip.GRCh38')
  const reference = async (start: number, end: number, context: number) =>
    (
      await walks.cut('GRCh38#0#chr1', { start, end }, { context })
    ).fragments.find(f => f.name === 'GRCh38#0#chr1')?.refStarts
  expect(await reference(32, 38, 5)).toEqual(await reference(27, 43, 0))
  expect(await reference(32, 38, 15)).toEqual(await reference(17, 53, 0))
})

test('a window inside a jump finds the jump in the departure read', async () => {
  const walks = reader('skip.GRCh38')
  const { fragments, departures } = await walks.cut(
    'GRCh38#0#chr1',
    { start: 31, end: 60 },
    { context: 0 },
  )
  expect(fragments.map(f => f.name)).toEqual(['GRCh38#0#chr1'])
  expect(departures.map(departure)).toEqual([
    ['HG002#1#ctg', 20, 0, 0, '20+', '80+', 'forward', 'record', undefined],
  ])
  const { tables } = await walks.cutRegion({
    refName: 'GRCh38#0#chr1',
    start: 31,
    end: 60,
  })
  expect(tables.departures.names).toEqual(['HG002#1#ctg'])
  expect([
    ...tables.departures.leavePos,
    ...tables.departures.rejoinPos,
  ]).toEqual([20, 80])
  const filtered = await walks.cut(
    'GRCh38#0#chr1',
    { start: 31, end: 60 },
    { keep: name => name === 'GRCh38#0#chr1' },
  )
  expect(filtered.departures).toEqual([])
})

// gfa-to-tabix --walks --chunk 20 of test_data/span.gfa, ten 10 bp reference
// nodes and three off it:
//   W GRCh38 0 chr1 >1>2>3>4>5>6>7>8>9>10
//   W HG002  1 chr1 >1>2>3>4>20>21>22>8>9>10   its own steps filed in chunk 20
//   W HG003  1 chr1 <10<9<8<22<21<20<4<3<2<1   its own steps filed in chunk 60
//   W HG004  1 chr1 <10<9<8<5<4<3<2<1          jumps from 70 to 50
// HG002#1 and HG003#1 leave or rejoin at 40 and at 70, and a cut of 50-55
// lies between
const spanCut = (
  haplotype: string,
  opts = {},
  walks = reader('span.GRCh38'),
  window = { start: 50, end: 55 },
) =>
  walks.cut('GRCh38#0#chr1', window, {
    context: 0,
    keep: walkNameFilter([haplotype], 'GRCh38'),
    ...opts,
  })
const atRejoin = { start: 70, end: 75 }

test('a record that jumps over the window keeps its ends and its own bp, and the cut fetches no tile for its own steps', async () => {
  let cut: Awaited<ReturnType<typeof spanCut>> | undefined
  const reads = await unitsRead('span.GRCh38', 'GRCh38#0#chr1', async walks => {
    cut = await spanCut('HG002#1', {}, walks)
  })
  expect(reads).toEqual([['tile 40', 'departures 40']])
  expect(cut?.fragments.map(f => f.name)).toEqual(['GRCh38#0#chr1'])
  expect(cut?.departures.map(departure)).toEqual([
    ['HG002#1#chr1', 40, 30, 0, '40+', '70+', 'forward', 'record', undefined],
  ])
  expect([...(cut?.kept ?? [])].sort((a, b) => a - b)).toEqual([6])
})

test("a record that rejoins on the window carries the walk's own steps, from the tile of its leave piece", async () => {
  let departures: WalkDeparture[] = []
  const reads = await unitsRead('span.GRCh38', 'GRCh38#0#chr1', async walks => {
    departures = (await spanCut('HG002#1', {}, walks, atRejoin)).departures
  })
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 40, 30, 0, '40+', '70+', 'forward', 'record', '>20>21>22'],
  ])
  expect(reads).toEqual([['tile 60', 'departures 60'], ['tile 20']])
})

test('a leave piece the read lacks is fetched for the own steps, once per tile, within the byte budget', async () => {
  const atLeave = { start: 35, end: 40 }
  const cut = (opts = {}, walks?: WalkReader) =>
    spanCut('HG003#1', opts, walks, atLeave)
  let departures: WalkDeparture[] = []
  const reads = await unitsRead('span.GRCh38', 'GRCh38#0#chr1', async walks => {
    departures = (await cut({}, walks)).departures
  })
  expect(departures.map(departure)).toEqual([
    ['HG003#1#chr1', 30, 30, 0, '70-', '40-', 'forward', 'record', '<22<21<20'],
  ])
  expect(reads).toEqual([['tile 20', 'departures 20'], ['tile 60']])
  const units = await reader('span.GRCh38').units('GRCh38#0#chr1')
  const bytes = units
    .filter(u => u.kind === 'tile' && (u.name === '20' || u.name === '60'))
    .reduce((sum, u) => sum + u.length, 0)
  expect((await cut({ byteBudget: bytes })).departures.map(departure)).toEqual(
    departures.map(departure),
  )
  await expect(cut({ byteBudget: bytes - 1 })).rejects.toBeInstanceOf(
    RegionTooLargeError,
  )
})

test('a record over the attach threshold carries no steps, and the cut fetches no tile for them', async () => {
  let cut: Awaited<ReturnType<typeof spanCut>> | undefined
  const reads = await unitsRead('span.GRCh38', 'GRCh38#0#chr1', async walks => {
    cut = await spanCut('HG002#1', { attachBp: 29 }, walks, atRejoin)
  })
  expect(reads).toEqual([['tile 60', 'departures 60']])
  expect(cut?.departures.map(departure)).toEqual([
    ['HG002#1#chr1', 40, 30, 0, '40+', '70+', 'forward', 'record', undefined],
  ])
  expect([...(cut?.kept ?? [])].sort((a, b) => a - b)).toEqual([8])
})

test('a direct jump across the window is a record with nothing of its own', async () => {
  const { fragments, departures } = await spanCut('HG004#1')
  expect(fragments.map(f => f.name)).toEqual(['GRCh38#0#chr1'])
  expect(departures.map(departure)).toEqual([
    ['HG004#1#chr1', 30, 0, 0, '70-', '50-', 'forward', 'record', undefined],
  ])
})

// Twelve 10 bp reference nodes in 20 bp chunks (attach.GRCh38.gfa); HG003#1
// leaves after 2, in chunk 0, for 21 and HG002#1 after 3, in chunk 20, for
// 20, and both rejoin at 11, in chunk 100:
//   W GRCh38 0 chr1 >1>2>3>4>5>6>7>8>9>10>11>12
//   W HG002  1 chr1 >1>2>3>20>11>12
//   W HG003  1 chr1 >1>2>21>11>12
test('two records whose leave pieces lie in adjacent chunks both get their own steps', async () => {
  const { departures } = await reader('attach.GRCh38').cut(
    'GRCh38#0#chr1',
    { start: 100, end: 120 },
    { context: 0 },
  )
  expect(departures.map(departure)).toEqual([
    ['HG002#1#chr1', 30, 10, 0, '30+', '100+', 'forward', 'record', '>20'],
    ['HG003#1#chr1', 20, 10, 0, '20+', '100+', 'forward', 'record', '>21'],
  ])
})

// An empty entry first in the table at the file's start, or after the
// first tile at that tile's offset, where it would shadow the tile in a read
// that took it
test.each([
  { at: 'first, at the file start', after: false },
  { at: "after the first tile, at the tile's offset", after: true },
])('a table entry of length 0 is never read: $at', async ({ after }) => {
  const { data, directory, units } = walkFileOf(copyJump)
  const bytes = new Uint8Array(await data.readFile())
  const end = bytes.length - 28
  const entry = (u: WalkUnit) => [
    '#entry',
    u.refName,
    u.start,
    u.end,
    u.kind,
    u.name,
    u.offset,
    u.length,
  ]
  const [first, ...rest] = units.map(entry)
  const empty = entry({
    ...units[0]!,
    offset: after ? units[0]!.offset : 0,
    length: 0,
  })
  const entries = [
    ['#unit', 'R#0#chr', 0, 48, 'table', '.'],
    ...(after ? [first!, empty] : [empty, first!]),
    ...rest,
  ]
  const table = bgzf(`${entries.map(e => e.join('\t')).join('\n')}\n`)
  const tableLength = table.reduce((sum, block) => sum + block.length, 0)
  const parts = [bytes.subarray(0, end), ...table, bytes.subarray(end)]
  const size = parts.reduce((sum, part) => sum + part.length, 0)
  const [walks = '', header = ''] = (await directory.readFile('utf8')).split(
    '\n',
  )
  const index = `${walks.replace(/\tbytes:i:\d+/, `\tbytes:i:${size}`)}\n${header}\nR#0#chr\t${end}\t${tableLength}\n`
  const headerReader = WalkReader.fromFilehandles({
    data: new BlobFile(new Blob(parts)),
    directory: new BlobFile(new Blob([index])),
  })
  expect(
    (await headerReader.units('R#0#chr')).filter(u => u.length === 0),
  ).toHaveLength(1)
  const cut = async (source: WalkReader) => {
    const { fragments, departures } = await source.cut(
      'R#0#chr',
      { start: 0, end: 48 },
      { context: 0 },
    )
    return [fragments.map(fragment), departures.map(departure)]
  }
  const whole = await cut(inline(copyJump))
  expect(whole[1]).toHaveLength(2)
  expect(await cut(headerReader)).toEqual(whole)
})

// header.GRCh38.gfa: HG002#1's record over 600-2400 has its own steps in the
// tile at 0, which a cut of 1700-2300, where it rejoins, fetches to attach
// them, and a cut of 2100-2600 does not
test("a sequence's length and the reference nodes over a window come from its tiles, and the bytes a cut reads from its table and its records' ownTiles", async () => {
  const walks = reader('header.GRCh38')
  const units = await walks.units('GRCh38#0#chr1')
  const tile = (name: string) =>
    units.find(u => u.kind === 'tile' && u.name === name)!
  expect(await walks.sequenceLength('GRCh38#0#chr1')).toBe(3600)
  expect(await walks.sequenceLength('GRCh38#0#chr2')).toBe(0)
  const departureBlock = (start: number) =>
    units.find(u => u.kind === 'departures' && u.start === start)!
  expect(
    await walks.bytesFor(
      'GRCh38#0#chr1',
      { start: 1700, end: 2300 },
      { context: 0 },
    ),
  ).toBe(
    tile('1000').length +
      departureBlock(1000).length +
      departureBlock(2000).length +
      tile('0').length,
  )
  expect(
    await walks.bytesFor(
      'GRCh38#0#chr1',
      { start: 2100, end: 2600 },
      { context: 0 },
    ),
  ).toBe(
    tile('1000').length + tile('2000').length + departureBlock(2000).length,
  )
  expect(
    (
      await walks.referenceNodes('GRCh38#0#chr1', { start: 1100, end: 2500 })
    ).map(n => `${n.id}@${n.start}-${n.end}`),
  ).toEqual(['2@600-1200', '3@1200-1800', '4@1800-2400', '5@2400-3000'])
})

// Ten 10 bp reference nodes in 20 bp tiles (jumpback.GRCh38.gfa). HG002#1's
// first W line ends on 9 and its second starts with 20 before rejoining at
// 3, so its record jumps back from 90 to 20 and its own step leads the
// rejoin piece, in tile 20:
//   W GRCh38 0 chr1 0  100 >1>2>3>4>5>6>7>8>9>10
//   W HG002  1 chr1 0  20  >8>9
//   W HG002  1 chr1 25 75  >20>3>4>5>6
// Its record's ownTiles names tile 20, which a cut over the leave at 90
// fetches and a cut over the rejoin at 20 holds; a cut of 50-60, which the
// jump back passes over, fetches nothing and states the record without them
test.each([
  { start: 80, end: 90, steps: '>20', fetched: ['tile 20'] },
  { start: 20, end: 30, steps: '>20', fetched: [] },
  { start: 50, end: 60, steps: undefined, fetched: [] },
])(
  'a jump back whose own steps lie in the rejoin piece gets them attached where it leaves or rejoins on the window, cut $start-$end',
  async ({ start, end, steps, fetched }) => {
    let departures: WalkDeparture[] = []
    const reads = await unitsRead(
      'jumpback.GRCh38',
      'GRCh38#0#chr1',
      async walks => {
        departures = (
          await walks.cut('GRCh38#0#chr1', { start, end }, { context: 0 })
        ).departures
      },
    )
    expect(departures.map(departure)).toEqual([
      ['HG002#1#chr1', 20, 15, 5, '90+', '20+', 'back', 'record', steps],
    ])
    expect(reads.slice(1).flat()).toEqual(fetched)
  },
)

test('a file past 4 GiB is read at its exact offsets', async () => {
  const gap = 2 ** 32 + 7
  const { data, directory, units } = walkFileOf(copyJump, { gap })
  expect(units[0]!.offset).toBeGreaterThan(2 ** 32)
  const walks = WalkReader.fromFilehandles({ data, directory })
  expect((await walks.units('R#0#chr')).map(u => u.offset)).toEqual(
    units.map(u => u.offset),
  )
  expect(await copyJumpCut(20, 36, false, walks)).toEqual(
    await copyJumpCut(20, 36, false),
  )
})

// ownends.R abridges every walk: X#1 leaves R after 2, at 12, takes 11, 12
// and 13 of its own and rejoins at 5, at 24, and both its pieces are
// abridged. Reading their detail units on every attach fetch cost pggb
// chrY's 1 kb windows 27 MB at p90, so only an expanded walk reads them
//   W R 0 c >1>2>3>4>5>6
//   W X 1 c >1>2>11>12>13>5>6
test.each([
  { expand: undefined, steps: 1 },
  { expand: (name: string) => name === 'X#1#c', steps: '>11>12>13' },
])(
  "an attach fetch reads an abridged piece's own steps where the cut expands the walk, and else holds a placeholder for them: $steps",
  async ({ expand, steps }) => {
    const { departures } = await reader('ownends.R').cut(
      'R#0#c',
      { start: 24, end: 30 },
      {
        context: 0,
        keep: name => name === 'X#1#c' || name === 'R#0#c',
        expand,
      },
    )
    expect(departures.map(departure)).toEqual([
      ['X#1#c', 12, 30, 0, '12+', '24+', 'forward', 'record', steps],
    ])
  },
)

// The reference passes 2 and 3 twice, 6 bp nodes in 10 bp tiles
// (loop.R.gfa): 2 at 6-12 and 18-24, 3 at 12-18 and 24-30
//   W R 0 chr >1>2>3>2>3>4
//   W A 1 a   >1>2>3>2>3>2>3>4
//   W F 1 f   >1>2>3>2>5>3>4   5 between the second passes of 2 and 3
test("a cut's drawn stretch takes every pass over the window, not one per node", async () => {
  const { fragments, departures } = await reader('loop.R').cut(
    'R#0#chr',
    { start: 8, end: 20 },
    { context: 0 },
  )
  expect(
    fragments
      .filter(f => f.name === 'R#0#chr' || f.name === 'F#1#f')
      .map(f => [f.name, f.hapStart, f.hapEnd, [...f.refStarts].join(',')]),
  ).toEqual([
    ['R#0#chr', 6, 24, '6,12,18'],
    ['F#1#f', 6, 24, '6,12,18'],
  ])
  expect(departures.filter(d => d.name === 'F#1#f').map(departure)).toEqual([
    ['F#1#f', 24, 2, 0, '24+', '24+', 'adjacent', 'record', undefined],
  ])
})

// F leaves at 24, the end of the second pass of 2, which a window to 20 draws
// whole and does not reach the end of
test.each([
  { end: 20, steps: undefined },
  { end: 24, steps: '>5' },
  { end: 25, steps: '>5' },
])(
  'a walk that leaves at the end of a node the window draws gets its own steps where the window holds the base it leaves from: to $end',
  async ({ end, steps }) => {
    const { fragments, departures } = await reader('loop.R').cut(
      'R#0#chr',
      { start: 8, end },
      { context: 0 },
    )
    expect(
      fragments
        .filter(f => f.name === 'F#1#f')
        .map(f => [...f.refStarts].slice(0, 3).join(',')),
    ).toEqual(['6,12,18'])
    expect(departures.filter(d => d.name === 'F#1#f').map(departure)).toEqual([
      ['F#1#f', 24, 2, 0, '24+', '24+', 'adjacent', 'record', steps],
    ])
  },
)

// drawnedge.R.gfa, one 10 kb tile over 10-20 kb that abridges G: G walks 11,
// to 16,000, then 300 bp of its own in place of 12, and rejoins at 13, at
// 17,000, all in one piece, so the cut computes the departure
//   W R 0 chr >1>...>25
//   W G 1 g   >5>6>7>8>9>10>11>200>13>14>15
test.each([
  { expand: undefined, own: 1 },
  { expand: () => true, own: '>200' },
])(
  'a departure inside one piece gets its own steps where it leaves on the window, and none where it jumps over it, expanded or not: $own',
  async ({ expand, own }) => {
    const g = async (start: number) => {
      const cut = await reader('drawnedge.R').cut(
        'R#0#chr',
        { start, end: start + 100 },
        { context: 100, expand },
      )
      return {
        departures: cut.departures
          .filter(d => d.name === 'G#1#g')
          .map(departure),
        gaps: cut.gaps
          .filter(gap => gap.name === 'G#1#g' && gap.start === 8000)
          .map(gap => gap.segments.map(x => `${x.kind} ${x.bp}`)),
        kept: cut.kept.size,
      }
    }
    const record = (steps: unknown) => [
      ['G#1#g', 8000, 300, 0, '16000+', '17000+', 'forward', 'steps', steps],
    ]
    expect(await g(16_400)).toEqual({
      departures: record(undefined),
      gaps: [['read 300']],
      kept: 1,
    })
    expect(await g(15_900)).toEqual({
      departures: record(own),
      gaps: [['steps 300', 'read 3000']],
      kept: 3,
    })
  },
)

// Twelve 4 bp reference nodes in 16 bp tiles. H's second piece goes to tile
// 32 and its record back to 20 names that tile; G's record names tile 0 and
// comes first in the block:
//   W R 0 chr >1>2>...>12
//   W H 1 h   >5>6>7>11>13>12>14>6>7>8   13, 14 own; pieces in tiles 16, 32, 16
//   W G 1 g   >1>15>6>7>8                15 own; pieces in tiles 0, 16
const farAttach = [
  '#walks format:i:8 tile:i:16',
  '#reference R',
  '#haplotype G#1',
  '#haplotype H#1',
  'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 .',
  'R#0#chr 0 nodes 15 4',
  'R#0#chr 0 route 2 28 .',
  'R#0#chr 0 walk G#1#g 0 0 0 2 2 0',
  'R#0#chr 16 ref R#0#chr 0 16 1 4 10,2,2,2 4,4,4,4 .',
  'R#0#chr 16 walk G#1#g 0 8 1 3 12 .',
  'R#0#chr 16 walk H#1#h 0 0 0 3 10 .',
  'R#0#chr 16 walk H#1#h 0 28 2 3 12 .',
  'R#0#chr 32 ref R#0#chr 0 32 2 4 18,2,2,2 4,4,4,4 .',
  'R#0#chr 32 nodes 13,1 4,4',
  'R#0#chr 32 route 22 4 -2',
  'R#0#chr 32 route 0 4 .',
  'R#0#chr 32 walk H#1#h 0 12 1 4 22 0,2',
  'R#0#chr 0 24 departure G#1#g 4 4 0 4 + 1 20 + 1 0',
  'R#0#chr 20 48 departure H#1#h 24 4 0 48 + 1 20 + 1 32',
  'R#0#chr 24 44 departure H#1#h 12 0 0 28 + 1 40 + 1 .',
]

test.each([
  { opts: {}, steps: '>14' },
  { opts: { keep: (name: string) => name !== 'G#1#g' }, steps: '>14' },
  { opts: { attachFetches: 0 }, steps: '>14' },
])(
  'a tile fetched to attach gives its record the own steps and nothing else: $opts',
  async ({ opts, steps }) => {
    const { fragments, departures } = await inline(farAttach).cut(
      'R#0#chr',
      { start: 20, end: 28 },
      { context: 0, ...opts },
    )
    expect(fragments.filter(f => f.name === 'H#1#h').map(fragment)).toEqual([
      ['H#1#h', 4, 12, '>6>7', 4, 0],
      ['H#1#h', 28, 36, '>6>7', 0, 4],
    ])
    expect(departures.filter(d => d.name === 'H#1#h').map(departure)).toEqual([
      ['H#1#h', 12, 0, 0, '28+', '40+', 'forward', 'record', undefined],
      ['H#1#h', 24, 4, 0, '48+', '20+', 'back', 'record', steps],
    ])
  },
)

// Forty 10 bp reference nodes in 10 bp tiles (attachmany.R.gfa). Each Hk,
// k from 1 to 10, leaves R after node k, at 10k, for its own 5 bp node 100+k
// and rejoins at 41-k, at 400-10k, so all ten jump over a cut in the middle,
// and their own steps lie in tiles 0 to 90:
//   W R  0 chr >1>2>...>40
//   W Hk 1 h   >1>...>k>(100+k)>(41-k)>...>40
const attachMany = { start: 195, end: 205 }
const attachManyRejoins = { start: 300, end: 400 }

const attachManyCut = async (
  window: { start: number; end: number },
  opts = {},
) => {
  let departures: WalkDeparture[] = []
  const reads = await unitsRead('attachmany.R', 'R#0#chr', async walks => {
    departures = (await walks.cut('R#0#chr', window, { context: 0, ...opts }))
      .departures
  })
  return {
    steps: departures.map(d => `${d.name} ${departure(d)[8]}`).sort(),
    ends: departures.map(d => `${d.name} ${d.leavePos} ${d.rejoinPos}`).sort(),
    ownBp: departures.map(d => d.ownBp),
    tiles: reads
      .flat()
      .filter(unit => unit.startsWith('tile'))
      .map(unit => Number(unit.slice(5)))
      .sort((a, b) => a - b),
  }
}
const tenth = <T>(of: (k: number) => T) =>
  Array.from({ length: 10 }, (_, k) => of(k + 1))

test('a cut attaches the own steps of every record that rejoins on the window, however many tiles hold them', async () => {
  const { steps, tiles } = await attachManyCut(attachManyRejoins)
  expect(steps).toEqual(tenth(k => `H${k}#1#h >${100 + k}`).sort())
  expect(tiles).toEqual([
    ...tenth(k => 10 * (k - 1)),
    ...tenth(k => 290 + 10 * k),
  ])
})

test('a cut fetches no tile for the own steps of records that jump over the window, and states each with its ends and its own bp', async () => {
  const { steps, ends, ownBp, tiles } = await attachManyCut(attachMany)
  expect(tiles).toEqual([190, 200])
  expect(steps).toEqual(tenth(k => `H${k}#1#h undefined`).sort())
  expect(ends).toEqual(tenth(k => `H${k}#1#h ${10 * k} ${400 - 10 * k}`).sort())
  expect(ownBp).toEqual(tenth(() => 5))
})

// Over 365-375, H1 and H2 jump over the window, H3 rejoins on it at 370, and
// H4 rejoins at 360, where the node the window draws from 360 to 370 starts
test('of the records over a window, a cut fetches a tile for those alone that leave or rejoin on it', async () => {
  const { steps, tiles } = await attachManyCut({ start: 365, end: 375 })
  expect(tiles).toEqual([20, 360, 370])
  expect(steps).toEqual([
    'H1#1#h undefined',
    'H2#1#h undefined',
    'H3#1#h >103',
    'H4#1#h undefined',
  ])
})

test('a cut whose attach fetches take it past the byte budget attaches none and says so', async () => {
  const walks = reader('attachmany.R')
  const tiles = (await walks.units('R#0#chr')).filter(u => u.kind === 'tile')
  const bytesOf = (starts: (start: number) => boolean) =>
    tiles.filter(u => starts(u.start)).reduce((sum, u) => sum + u.length, 0)
  const bytes = bytesOf(start => start < 100 || start >= 300)
  const cut = (byteBudget: number) =>
    walks.cut('R#0#chr', attachManyRejoins, { context: 0, byteBudget })
  expect((await cut(bytes)).departures.every(d => d.steps !== undefined)).toBe(
    true,
  )
  await expect(cut(bytes - 1)).rejects.toMatchObject({
    regionTooLarge: true,
    zoomHelps: true,
  })
  const firstTile = bytesOf(start => start < 100 || start === 300)
  const message: unknown = expect.stringMatching(
    /^Too many departures here whose own steps the cut must fetch.*choose fewer haplotypes or a smaller attachBp$/,
  )
  await expect(cut(firstTile - 1)).rejects.toMatchObject({
    regionTooLarge: true,
    zoomHelps: false,
    message,
  })
})

test('the byte budget counts the tiles of the context past the window', async () => {
  const walks = reader('attachmany.R')
  const opts = { context: 10 }
  const bytes = (await walks.units('R#0#chr'))
    .filter(u => u.kind === 'tile' && u.start >= 180 && u.start <= 210)
    .reduce((sum, u) => sum + u.length, 0)
  await walks.cut('R#0#chr', attachMany, { ...opts, byteBudget: bytes })
  await expect(
    walks.cut('R#0#chr', attachMany, { ...opts, byteBudget: bytes - 1 }),
  ).rejects.toBeInstanceOf(RegionTooLargeError)
})

// Four 6 bp reference nodes in 10 bp tiles (ownzero.R.gfa); G jumps back
// from 4 to 1 over node 9 of 0 bp, which the writer drops, so its record has
// no own step:
//   W R 0 chr 0 24 >1>2>3>4
//   W G 1 g 0 30 >1>2>3>4>9>1
test.each([
  { start: 0, end: 6 },
  { start: 0, end: 24 },
])(
  'a record over a 0 bp segment alone holds no step, cut $start-$end',
  async ({ start, end }) => {
    const { departures } = await reader('ownzero.R').cut(
      'R#0#chr',
      { start, end },
      { context: 0 },
    )
    expect(departures.map(departure)).toEqual([
      ['G#1#g', 24, 0, 0, '24+', '0+', 'back', 'record', undefined],
    ])
  },
)

// Two reference lines, 0-8 and 40-48, with a hole between, in 16 bp tiles
// (hole.R.gfa); H jumps from 1 to 4 across it:
//   W R 0 chr 0  8  >1>2
//   W R 0 chr 40 48 >3>4
//   W H 1 h   0  8  >1>4
test.each([
  { start: 16, end: 24 },
  { start: 0, end: 48 },
])(
  'a cut returns every record whose line overlaps its window, in a hole too: $start-$end',
  async ({ start, end }) => {
    const { departures } = await reader('hole.R').cut(
      'R#0#chr',
      { start, end },
      { context: 0 },
    )
    expect(departures.map(departure)).toEqual([
      ['H#1#h', 4, 0, 0, '4+', '44+', 'forward', 'record', undefined],
    ])
  },
)

// Twelve 4 bp reference nodes in one tile (excursion.R.gfa); H goes from 1
// to 9, jumps to 11 and comes back to 2, all in one piece:
//   W R 0 chr >1>2>...>12
//   W H 1 h   >1>9>11>2
test('a cut computes no departure outside its drawn stretch, even inside a fragment', async () => {
  const { fragments, departures } = await reader('excursion.R').cut(
    'R#0#chr',
    { start: 0, end: 8 },
    { context: 0 },
  )
  expect(fragments.filter(f => f.name === 'H#1#h').map(fragment)).toEqual([
    ['H#1#h', 0, 16, '>1>9>11>2', 0, 0],
  ])
  expect(departures.map(departure)).toEqual([
    ['H#1#h', 4, 0, 0, '4+', '32+', 'forward', 'steps', undefined],
    ['H#1#h', 12, 0, 0, '44+', '4+', 'back', 'steps', undefined],
  ])
})

test('referenceNodes gives every pass of the reference over the region', async () => {
  const nodes = await reader('loop.R').referenceNodes('R#0#chr', {
    start: 0,
    end: 36,
  })
  expect(nodes.map(n => `${n.id}@${n.start}-${n.end}`)).toEqual([
    '1@0-6',
    '2@6-12',
    '3@12-18',
    '2@18-24',
    '3@24-30',
    '4@30-45',
  ])
})

const linkLine = ([a, b]: [number, number]) =>
  `${a >> 1} ${a & 1 ? '-' : '+'} ${b >> 1} ${b & 1 ? '-' : '+'}`

test.each([
  { a: 14, b: 6, line: '7 + 3 +' },
  { a: 7, b: 15, line: '7 + 3 +' },
  { a: 16, b: 12, line: '8 + 6 +' },
  { a: 13, b: 17, line: '8 + 6 +' },
  { a: 14, b: 7, line: '3 + 7 -' },
  { a: 6, b: 15, line: '3 + 7 -' },
  { a: 15, b: 6, line: '3 - 7 +' },
])(
  'a link is written the way round with more + ends, the smaller first on a tie: $line',
  ({ a, b, line }) => {
    expect(linkLine(canonicalLink(a, b))).toBe(line)
  },
)

test("a cut's links run the way a walk along the reference takes them", async () => {
  const { graph, kept, fragments, departures } = await reader(
    'excursion.R',
  ).cut('R#0#chr', { start: 0, end: 8 }, { context: 0 })
  const links = graphTablesGFA(graph.tables(kept, fragments, departures))
    .split('\n')
    .filter(line => line.startsWith('L\t'))
  expect(links).toContain('L\t11\t+\t2\t+\t0M')
  expect(links).not.toContain('L\t2\t-\t11\t-\t0M')
})

// H's two W lines on one contig, each numbering its pieces from 0, over 10 bp
// nodes in 20 bp tiles (twolines.R.gfa):
//   W R 0 chr   0  40 >1>2>3>4
//   W H 1 h     0  40 >1>2<3>4
//   W H 1 h   100 140 >1>2>3>4
test("a link between two pieces of one W line holds though the contig's other W line numbers its pieces the same", async () => {
  const { graph, kept, fragments, departures } = await reader('twolines.R').cut(
    'R#0#chr',
    { start: 0, end: 40 },
    { context: 0 },
  )
  expect(
    graphTablesGFA(graph.tables(kept, fragments, departures)).split('\n'),
  ).toContain('L\t2\t+\t3\t-\t0M')
})

const placedSteps = (f: WalkFragment) =>
  [...f.ids]
    .map((id, k) => `${f.rev[k]! !== 0 ? '<' : '>'}${id}@${f.refStarts[k]}`)
    .join('')

// The reference passes 2 and 3 twice, 6 bp nodes in one tile (copies.R.gfa).
// F leaves from 2's second pass, `:1`, and lands on 3 at the pass as near
// the entry as the first, the later; U lands on 3's second pass, `@1`; S
// starts on 2's second pass, `start@1`:
//   W R 0 chr >1>2>3>2>3>4
//   W F 1 f   >1>2>3>2>8>3>4   8 of 2 bp
//   W U 1 u   >1>9>3>4         9 of 18 bp
//   W S 1 s   >2>3>4
test('a walk row takes the passes its :k, @k and start@k name, and the defaults where it names none', async () => {
  const { fragments } = await reader('copies.R').cut(
    'R#0#chr',
    { start: 0, end: 36 },
    { context: 0 },
  )
  expect(fragments.map(f => `${f.name} ${placedSteps(f)}`).sort()).toEqual([
    'F#1#f >1@0>2@6>3@12>2@18>8@NaN>3@24>4@30',
    'R#0#chr >1@0>2@6>3@12>2@18>3@24>4@30',
    'S#1#s >2@18>3@24>4@30',
    'U#1#u >1@0>9@NaN>3@24>4@30',
  ])
})

test("a cut's tables state the pass each step takes, NaN for a step off R", async () => {
  const { graph, kept, fragments, departures } = await reader('copies.R').cut(
    'R#0#chr',
    { start: 0, end: 36 },
    { context: 0 },
  )
  const { walks, nodes } = graph.tables(kept, fragments, departures)
  expect(
    walks.names.map((name, p) => {
      const from = walks.offsets[p]!
      const placed = Array.from(
        walks.steps.subarray(from, walks.offsets[p + 1]),
        (node, k) => `${nodes.names[node]}@${walks.refStarts![from + k]}`,
      )
      return `${name} ${placed.join(' ')}`
    }),
  ).toEqual([
    'R#0#chr 1@0 2@6 3@12 2@18 3@24 4@30',
    'F#1#f 1@0 2@6 3@12 2@18 8@NaN 3@24 4@30',
    'S#1#s 2@18 3@24 4@30',
    'U#1#u 1@0 9@NaN 3@24 4@30',
  ])
})

async function copiesPlaced(
  name: string,
  edit: (f: WalkFragment) => WalkFragment,
) {
  const { graph, kept, fragments, departures } = await reader('copies.R').cut(
    'R#0#chr',
    { start: 0, end: 36 },
    { context: 0 },
  )
  const { walks, nodes } = graph.tables(
    kept,
    fragments.map(f => (f.name === name ? edit({ ...f }) : f)),
    departures,
  )
  const p = walks.names.indexOf(name)
  const from = walks.offsets[p]!
  return Array.from(
    walks.steps.subarray(from, walks.offsets[p + 1]),
    (node, k) => `${nodes.names[node]}@${walks.refStarts![from + k]}`,
  ).join(' ')
}

test("a host's fragment with no refStarts gives its steps NaN places", async () => {
  expect(
    await copiesPlaced('F#1#f', f => {
      Reflect.deleteProperty(f, 'refStarts')
      return f
    }),
  ).toBe('1@NaN 2@NaN 3@NaN 2@NaN 8@NaN 3@NaN 4@NaN')
})

test("a host's fragment with ids trimmed and refStarts not gives the places of the steps it keeps", async () => {
  expect(
    await copiesPlaced('U#1#u', f => ({
      ...f,
      ids: f.ids.subarray(0, 3),
      rev: f.rev.subarray(0, 3),
    })),
  ).toBe('1@0 9@NaN 3@24')
})

test("a host's fragment with fewer refStarts than steps gives the rest NaN places, not 0", async () => {
  expect(
    await copiesPlaced('S#1#s', f => ({
      ...f,
      refStarts: f.refStarts.subarray(0, 2),
    })),
  ).toBe('2@18 3@24 4@NaN')
})

// copies.R's header states revisited:i:2, skip.GRCh38's revisited:i:0, and
// the inline file's none
test("a cut's tables carry the header's revisited count, NaN where it states none", async () => {
  const revisited = async (walks: WalkReader, refName: string, end: number) => {
    const { graph, kept, fragments, departures } = await walks.cut(
      refName,
      { start: 0, end },
      { context: 0 },
    )
    return graph.tables(kept, fragments, departures).revisited
  }
  expect(await revisited(reader('copies.R'), 'R#0#chr', 36)).toBe(2)
  expect(await revisited(reader('skip.GRCh38'), 'GRCh38#0#chr1', 100)).toBe(0)
  expect(
    await revisited(
      inline([
        '#walks format:i:8 tile:i:16',
        '#reference R',
        'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 4,4,4,4 .',
      ]),
      'R#0#chr',
      16,
    ),
  ).toBeNaN()
})

// The reference passes 2 and 3 twice, 6 bp nodes in one tile (early.R.gfa).
// E's route to 8 names no `:k`, so it leaves from the first pass of 2 its walk
// reaches, not the last:
//   W R 0 chr >1>2>3>2>3>4
//   W E 1 e   >1>2>8>3>2>3>4   8 of 2 bp
test('a route naming no :k leaves from the first pass of its entry the walk reaches', async () => {
  const { fragments } = await reader('early.R').cut(
    'R#0#chr',
    { start: 0, end: 36 },
    { context: 0 },
  )
  expect(fragments.map(f => `${f.name} ${placedSteps(f)}`)).toContain(
    'E#1#e >1@0>2@6>8@NaN>3@12>2@18>3@24>4@30',
  )
})

// Every walk abridged, in one tile (copyends.R.gfa), and each unloaded piece
// read from its blocks as its loaded one is. Y's end is 3 at its
// second pass, the one nearest its start, and V's is 3 at its second pass
// too, which `end@1` names as the first is nearer its start:
//   W R 0 chr >1>2>3>6>7>8>4>2>3>5
//   W Y 1 y   >4>2>3
//   W V 1 v   >1>2>3>6>7>8>4>2>3
test("an abridged row's end takes the pass end@k names, else the one nearest its start", async () => {
  const cut = (expand?: (name: string) => boolean) =>
    reader('copyends.R').cut(
      'R#0#chr',
      { start: 0, end: 60 },
      { context: 0, expand },
    )
  const { fragments } = await cut()
  expect(
    fragments
      .map(f => `${f.name} ${f.hapStart}-${f.hapEnd} ${placedSteps(f)}`)
      .sort(),
  ).toEqual([
    'R#0#chr 0-60 >1@0>2@6>3@12>6@18>7@24>8@30>4@36>2@42>3@48>5@54',
    'V#1#v 0-54 >1@0>2@6>3@12>6@18>7@24>8@30>4@36>2@42>3@48',
    'Y#1#y 0-18 >4@36>2@42>3@48',
  ])
  const expanded = await cut(() => true)
  expect(
    expanded.fragments.map(f => `${f.name} ${placedSteps(f)}`).sort(),
  ).toEqual([
    'R#0#chr >1@0>2@6>3@12>6@18>7@24>8@30>4@36>2@42>3@48>5@54',
    'V#1#v >1@0>2@6>3@12>6@18>7@24>8@30>4@36>2@42>3@48',
    'Y#1#y >4@36>2@42>3@48',
  ])
})

// 4 bp reference nodes in 8 bp tiles, with 9 of 0 bp on the reference at 12
// (zerostep.R.gfa), which the writer drops. G jumps from 2 to 9 and from 9
// to 6, one departure from 2 to 6 without it:
//   W R 0 chr >1>2>3>9>4>5>6>7
//   W G 1 g   >1>2>9>6>7
test('a walk through a 0 bp segment of the reference departs once over it', async () => {
  const { departures } = await reader('zerostep.R').cut(
    'R#0#chr',
    { start: 0, end: 28 },
    { context: 0 },
  )
  expect(departures.map(departure)).toEqual([
    ['G#1#g', 8, 0, 0, '8+', '20+', 'forward', 'record', undefined],
  ])
})

// hole.R.gfa in one 64 bp tile: H jumps from 1 to 4 across the hole, and J
// steps straight across it from 2 to 3, which the hole makes adjacent
test('steps either side of a hole in the reference depart nowhere', async () => {
  const { fragments, departures } = await reader('holetile.R').cut(
    'R#0#chr',
    { start: 0, end: 48 },
    { context: 0 },
  )
  expect(fragments.find(f => f.name === 'J#1#j')?.ids).toEqual(
    Int32Array.from([1, 2, 3, 4]),
  )
  expect(departures.map(departure)).toEqual([
    ['H#1#h', 4, 0, 0, '4+', '44+', 'forward', 'steps', undefined],
  ])
})

// Each link of a cut's tables, as its two nodes and their strands
function linkNames(cut: Awaited<ReturnType<WalkReader['cut']>>) {
  const { nodes, links } = cut.graph.tables(
    cut.kept,
    cut.fragments,
    cut.departures,
  )
  return Array.from(
    links.from,
    (from, k) =>
      `${nodes.names[from]}${links.strands[k]! & 1 ? '-' : '+'} ${nodes.names[links.to[k]!]}${links.strands[k]! & 2 ? '-' : '+'}`,
  ).sort()
}

// refhole.R.gfa: R is two W lines, >1>2 at 0-8 and >3>4 at 40-48, in one
// 64 bp tile, and H jumps from 1 to 4. No walk goes from 2 to 3, so the graph
// has no such link
test('a cut lists no link over a hole in the reference that no walk crosses, where both sides lie in one tile', async () => {
  const cut = await reader('refhole.R').cut(
    'R#0#chr',
    { start: 0, end: 64 },
    { context: 0 },
  )
  expect(linkNames(cut)).toEqual(['1+ 2+', '1+ 4+', '3+ 4+'])
})

// hprc-chrM.GRCh38.gfa: HPRC's chrM, 44 haplotypes, no walk abridged. The
// file's tile holds every haplotype's rows, whichever a cut keeps
test('a cut of four haplotypes lists the links its fragments walk, and between reference nodes no other', async () => {
  const walks = reader('hprc-chrM.GRCh38')
  const [refName = ''] = await walks.referenceSequences()
  const { haplotypes } = await walks.header()
  const cut = await walks.cut(
    refName,
    { start: 0, end: 16_000 },
    { keep: walksKept(haplotypes.slice(0, 4), refName) },
  )
  const walked = new Set<string>()
  for (const { ids, rev } of cut.fragments) {
    for (let i = 1; i < ids.length; i++) {
      const [a, b] = canonicalLink(
        ids[i - 1]! * 2 + rev[i - 1]!,
        ids[i]! * 2 + rev[i]!,
      )
      walked.add(`${a >> 1}${a & 1 ? '-' : '+'} ${b >> 1}${b & 1 ? '-' : '+'}`)
    }
  }
  const onReference = (name: string) =>
    cut.graph.nodes.get(Number(name.slice(0, -1)))?.rank === 0
  const listed = linkNames(cut)
  expect(cut.fragments.length).toBe(5)
  expect(listed).toEqual(expect.arrayContaining([...walked]))
  expect(
    listed.filter(
      link => link.split(' ').every(onReference) && !walked.has(link),
    ),
  ).toEqual([])
})

// adv.R.gfa: R over twelve 10 bp nodes. CROSS leaves after node 3 for 21,
// two nodes of another sequence and 22, and rejoins at node 10; NOREF's
// first line ends on 21 and a later one starts at node 10; PLAIN leaves
// after node 3 for 24 and rejoins at node 10. A cut holds 21 as an own step
// of CROSS and of NOREF, and neither walks from 21 to node 10
describe('the links at a departure record’s ends', () => {
  const all = () => true

  test('own steps that stop short of the rejoin node are not linked to it', async () => {
    const cut = await reader('adv.R').cut(
      'R#0#c1',
      { start: 84, end: 96 },
      { context: 0, expand: all },
    )
    expect(
      cut.departures.map(d => [d.name, [...(d.steps?.ids ?? [])]]),
    ).toEqual([
      ['CROSS#1#c', [21]],
      ['NOREF#1#c', [21]],
      ['PLAIN#1#c', [24]],
    ])
    expect(linkNames(cut)).toEqual(['24+ 10+', '9+ 10+'])
  })

  test('the node a record leaves is linked to its first own step, and its last to the node it rejoins, where the cut holds them', async () => {
    const cut = await reader('adv.R').cut(
      'R#0#c1',
      { start: 0, end: 15 },
      { context: 1000, expand: all },
    )
    const links = linkNames(cut)
    expect(links).toContain('3+ 21+')
    expect(links).toContain('3+ 24+')
    expect(links).toContain('24+ 10+')
    expect(links.filter(link => link.includes('21'))).toEqual(['3+ 21+'])
  })

  // jumpback.GRCh38.gfa: HG002's second W line starts on node 20, off the
  // reference, and rejoins at node 3
  test('own steps before the rejoin are linked to the node the record rejoins', async () => {
    const cut = await reader('jumpback.GRCh38').cut(
      'GRCh38#0#chr1',
      { start: 0, end: 16 },
      { context: 10 },
    )
    expect(linkNames(cut)).toContain('20+ 3+')
  })

  // abridged.R.gfa: H leaves R after node 3 for node 7, 3 bp, and rejoins at
  // node 4, inside one piece, so the cut computes the departure from H's steps
  test.each([
    { expand: all, links: ['3+ 7+', '7+ 4+'] },
    { expand: undefined, links: ['3+ H#1#h:118-121+', 'H#1#h:118-121+ 4+'] },
  ])(
    'a departure inside a piece links its own steps to the steps either side: $links',
    async ({ expand, links }) => {
      const cut = await reader('abridged.R').cut(
        'R#0#chr',
        { start: 10, end: 20 },
        { context: 0, ...(expand ? { expand } : {}) },
      )
      expect(linkNames(cut)).toEqual(expect.arrayContaining(links))
    },
  )

  // the same departure in a window that ends before node 4: H's fragment
  // stops at node 3, and D's, which goes round by node 4, keeps that node
  test.each([
    { expand: all, link: '3+ 7+', none: '7+ 4+' },
    {
      expand: undefined,
      link: '3+ H#1#h:118-121+',
      none: 'H#1#h:118-121+ 4+',
    },
  ])(
    'own steps are not linked to a rejoin step their walk’s fragment lacks, though another walk keeps its node: $none',
    async ({ expand, link, none }) => {
      const cut = await reader('abridged.R').cut(
        'R#0#chr',
        { start: 10, end: 18 },
        { context: 0, ...(expand ? { expand } : {}) },
      )
      const links = linkNames(cut)
      expect(cut.kept.has(4)).toBe(true)
      expect(links).toContain(link)
      expect(links).not.toContain(none)
    },
  )

  // abridged.R.gfa: Q leaves R after node 1 for 6 bp of its own and rejoins
  // at node 3; its own steps are one placeholder in a cut that expands none
  test('a placeholder for own steps is linked to the steps its record leaves and rejoins', async () => {
    const cut = await reader('abridged.R').cut(
      'R#0#chr',
      { start: 4, end: 16 },
      { context: 0 },
    )
    expect(linkNames(cut)).toEqual([
      '1+ 2+',
      '1+ 3-',
      '1+ Q#0#chr:6-12+',
      '2+ 3+',
      '3+ 4+',
      '4+ 1+',
      'Q#0#chr:6-12+ 3+',
    ])
  })

  // a window that starts after node 1: Q's fragment starts at node 3, and
  // D's keeps node 1
  test('a placeholder for own steps is not linked to a leave step its walk’s fragment lacks', async () => {
    const cut = await reader('abridged.R').cut(
      'R#0#chr',
      { start: 8, end: 16 },
      { context: 0 },
    )
    expect(cut.kept.has(1)).toBe(true)
    expect(linkNames(cut)).toEqual([
      '1+ 2+',
      '2+ 3+',
      '3+ 4+',
      '4+ 1+',
      'Q#0#chr:6-12+ 3+',
    ])
  })
})

// abut.R: ONE walks 1, 2, its own 20, then 4, 5, 6 as one W line; ABUT the
// same steps as two lines, the second starting on its contig where the first
// ends. Their departure records state the same
describe('a record between two W lines that abut', () => {
  const cutOf = async (haplotype: string, start: number, end: number) => {
    const cut = await reader('abut.R').cut(
      'R#0#c',
      { start, end },
      { context: 0, keep: walksKept([haplotype], 'R#0#c') },
    )
    return {
      links: linkNames(cut),
      departures: cut.departures.map(d => ({
        ...d,
        name: '',
        steps: [...(d.steps?.ids ?? [])],
      })),
    }
  }

  test('a cut holding both pieces links the own step to the rejoin step on one line, and states a break between two', async () => {
    const one = await cutOf('ONE#1', 10, 40)
    const abut = await cutOf('ABUT#1', 10, 40)
    expect(one.links).toContain('20+ 4+')
    expect(abut.links).not.toContain('20+ 4+')
    expect(one.departures.map(d => [d.steps, d.breakAt])).toEqual([
      [[20], undefined],
    ])
    expect(abut.departures.map(d => [d.steps, d.breakAt])).toEqual([[[20], 1]])
    expect(abut.links).toEqual(one.links.filter(link => link !== '20+ 4+'))
  })

  test('a cut holding the leave piece alone states the same of both', async () => {
    const one = await cutOf('ONE#1', 10, 20)
    const abut = await cutOf('ABUT#1', 10, 20)
    expect(one.departures.map(d => [d.steps, d.breakAt])).toEqual([[[20], 1]])
    expect(abut).toEqual(one)
  })
})

// farrejoin.R in 10 bp tiles: X walks 1, 2, its own 50, then 8, its own 51,
// then 3 and 4. A cut of [0, 30) reads the piece leaving 2 and fetches tile
// 70 for the own steps of the record rejoining at 3, and that tile holds the
// piece rejoining at 8 too. The cut returns no step on 8
test('a record whose rejoin piece comes in a fetch for another record breaks where the cut returns no step beside its own', async () => {
  const cut = await reader('farrejoin.R').cut(
    'R#0#c',
    { start: 0, end: 30 },
    { context: 0 },
  )
  expect(
    cut.departures.map(d => [d.hapStart, [...(d.steps?.ids ?? [])], d.breakAt]),
  ).toEqual([
    [20, [50], 1],
    [35, [51], undefined],
  ])
})

// A cut's links are those each of its walks' steps walk, so cutting several
// haplotypes lists the links of cutting each alone with the reference
describe('the links of a cut of several haplotypes are those of each cut alone', () => {
  test.each([
    'abridged.R',
    'unabridged.R',
    'adv.R',
    'ownends.R',
    'between.R',
    'blocks.R',
    'copies.R',
    'unbuilt.R',
    'attach.GRCh38',
    'jumpback.GRCh38',
    'abut.R',
  ])('%s', async name => {
    const walks = reader(name)
    const { haplotypes } = await walks.header()
    const differing: string[] = []
    for (const refName of await walks.referenceSequences()) {
      const end = (await walks.units(refName)).reduce(
        (most, unit) => Math.max(most, unit.end),
        0,
      )
      for (let start = 0; start < end; start += 6) {
        for (const context of [0, 3, 1000]) {
          for (const contained of [false, true]) {
            for (const expand of [undefined, () => true]) {
              const cutOf = async (kept: string[]) =>
                linkNames(
                  await walks.cut(
                    refName,
                    { start, end: start + 8 },
                    {
                      keep: walksKept(kept, refName),
                      context,
                      contained,
                      ...(expand ? { expand } : {}),
                    },
                  ),
                )
              const alone = new Set<string>()
              for (const haplotype of haplotypes) {
                for (const link of await cutOf([haplotype])) {
                  alone.add(link)
                }
              }
              const together = await cutOf(haplotypes)
              if (together.join(' | ') !== [...alone].sort().join(' | ')) {
                differing.push(
                  `${refName}:${start}-${start + 8} context ${context}${contained ? ' contained' : ''}${expand ? ' expanded' : ''}: ${together.join(' | ')} AGAINST ${[...alone].sort().join(' | ')}`,
                )
              }
            }
          }
        }
      }
    }
    expect(differing).toEqual([])
  })
})

// holetile.R's tile with its ref rows written last first: R is the rows in
// hapOffset order whatever order the tile lists them in
test("a tile's R joins its ref rows in hapOffset order", async () => {
  const { fragments } = await inline([
    '#walks format:i:8 tile:i:64',
    '#reference R',
    '#haplotype J#1',
    '#hole R#0#chr 8 40',
    'R#0#chr 0 ref R#0#chr 40 40 0 2 6,2 4,4 .',
    'R#0#chr 0 ref R#0#chr 0 0 0 2 2,2 4,4 .',
    'R#0#chr 0 walk J#1#j 0 0 0 4 2 .',
  ]).cut('R#0#chr', { start: 0, end: 48 }, { context: 0 })
  expect(fragments.map(f => `${f.name} ${placedSteps(f)}`).sort()).toEqual([
    'J#1#j >1@0>2@4>3@40>4@44',
    'R#0#chr >1@0>2@4',
    'R#0#chr >3@40>4@44',
  ])
})

// attach.GRCh38's HG003#1 leaves after 2, at 10-20, so its line starts at
// 10, in the departure block of 0-20
test.each([
  { end: 10, names: [] },
  { end: 11, names: ['HG003#1#chr1'] },
])(
  'a cut keeps the departure lines overlapping its window, of the blocks it reads: 0-$end',
  async ({ end, names }) => {
    const { departures } = await reader('attach.GRCh38').cut(
      'GRCh38#0#chr1',
      { start: 0, end },
      { context: 0 },
    )
    expect(departures.map(d => d.name)).toEqual(names)
  },
)

const lengthOf = (blocks: Uint8Array[]) =>
  blocks.reduce((sum, block) => sum + block.length, 0)
const REF_ROW = 'R#0#chr\t0\t4\tref\tR#0#chr\t0\t0\t0\t1\t2\t4\t.\n'
const TILE_HEAD = '#unit\tR#0#chr\t0\t4\ttile\t0\n'
const tableOf = (at: number, length: number, kind = 'tile', name = '0') =>
  `#unit\tR#0#chr\t0\t4\ttable\t.\n#entry\tR#0#chr\t0\t4\t${kind}\t${name}\t${at}\t${length}\n`

// The filehandles of a walk file whose one sequence, R#0#chr, has one unit
// over [0, 4), held in `unit`, and a table of `table`'s text listing it
function oneUnit(unit: Uint8Array<ArrayBuffer>[], table = tableOf) {
  const walks = '#walks\tformat:i:8\ttile:i:16'
  const header = bgzf(`${walks}\n#reference\tR\n`)
  const unitAt = lengthOf(header)
  const tableAt = unitAt + lengthOf(unit)
  const tableUnit = bgzf(table(unitAt, lengthOf(unit)))
  const parts = [...header, ...unit, ...tableUnit]
  const directory = `${walks}\tbytes:i:${lengthOf(parts)}\n#header\t0\t${unitAt}\nR#0#chr\t${tableAt}\t${lengthOf(tableUnit)}\n`
  return {
    data: new BlobFile(new Blob(parts)),
    directory: new BlobFile(new Blob([directory])),
  }
}

// A walk file whose one sequence, R#0#chr, has one tile, 0 over [0, 4)
const oneTile = (tile: Uint8Array<ArrayBuffer>[], table = tableOf) =>
  WalkReader.fromFilehandles(oneUnit(tile, table))

const cutAll = (walks: WalkReader) =>
  walks.cut('R#0#chr', { start: 0, end: 4 }, { context: 0 })

const TILE_AT =
  "The unit at byte 69 of the walk file, R#0#chr's tile 0 over [0, 4) as its table lists it,"

test('a one-tile file whose tile holds what its table says is read', async () => {
  const { fragments } = await cutAll(oneTile(bgzf(TILE_HEAD + REF_ROW)))
  expect(fragments.map(f => f.name)).toEqual(['R#0#chr'])
})

test.each([
  ['no #unit line', REF_ROW, `${TILE_AT} starts "R#0#chr\\t0\\t4\\tref`],
  [
    "another tile's #unit line, as a table whose entries swapped their byte ranges gives",
    TILE_HEAD.replace('tile\t0', 'tile\t16') + REF_ROW,
    `${TILE_AT} starts "#unit\\tR#0#chr\\t0\\t4\\ttile\\t16"`,
  ],
  [
    'a half row, as a table length ending at a member boundary inside the unit gives',
    `${TILE_HEAD}${REF_ROW}R#0#chr\t0\t4\twa`,
    `${TILE_AT} ends partway through a line`,
  ],
  [
    'a row kind the format does not give',
    `${TILE_HEAD}${REF_ROW}R#0#chr\t0\t4\tsketch\t1\n`,
    `${TILE_AT} holds a row that is not one of R#0#chr's ref, nodes, route, walk, abridged rows: "R#0#chr\\t0\\t4\\tsketch\\t1"`,
  ],
  [
    "another sequence's row",
    `${TILE_HEAD}${REF_ROW.replace(/^R#0#chr/, 'R#0#chr2')}`,
    `${TILE_AT} holds a row that is not one of R#0#chr's ref, nodes, route, walk, abridged rows: "R#0#chr2\\t0`,
  ],
])('a tile with %s is refused', async (_, text, message) => {
  await expect(cutAll(oneTile(bgzf(text)))).rejects.toThrow(message)
})

test.each([
  ['departures', '.', REF_ROW, "R#0#chr's departure rows"],
  ['detail', '0:0', REF_ROW, "R#0#chr's nodes, route, walk rows"],
  ['blocks', '0', REF_ROW, "R#0#chr's blocks rows"],
  [
    'tile',
    '0',
    'R#0#chr\t0\t4\tdeparture\n',
    "R#0#chr's ref, nodes, route, walk, abridged rows",
  ],
])(
  'a %s unit holding a row of a kind it does not hold is refused',
  async (kind, name, row, kinds) => {
    const unit = bgzf(`#unit\tR#0#chr\t0\t4\t${kind}\t${name}\n${row}`)
    const { data, directory } = oneUnit(
      unit,
      (at, length) =>
        `${tableOf(at, length, kind, name)}${kind === 'detail' ? `#entry\tR#0#chr\t0\t4\tblocks\t0\t${at}\t${length}\n` : ''}`,
    )
    const file = new WalkFile(data, directory)
    const units = await file.units('R#0#chr')
    await expect(file.read(units.filter(u => u.kind === kind))).rejects.toThrow(
      `holds a row that is not one of ${kinds}`,
    )
  },
)

test('a unit of a kind the format does not give is read whatever its rows', async () => {
  const unit = bgzf(`#unit\tR#0#chr\t0\t4\tsketch\t.\nanything\n`)
  const { data, directory } = oneUnit(unit, (at, length) =>
    tableOf(at, length, 'sketch', '.'),
  )
  const file = new WalkFile(data, directory)
  expect(await file.read(await file.units('R#0#chr'))).toEqual([['anything']])
})

test('a tile cut partway through a gzip member is refused', async () => {
  const [whole] = bgzf(TILE_HEAD + REF_ROW)
  await expect(
    cutAll(oneTile([whole!.subarray(0, whole!.length - 10)])),
  ).rejects.toThrow(`${TILE_AT} is not whole gzip members`)
})

test.each([
  [
    'a #unit line of another kind',
    (at: number, length: number) =>
      tableOf(at, length).replace('table', 'tile'),
    'of the walk file, starts "#unit\\tR#0#chr\\t0\\t4\\ttile\\t."',
  ],
  [
    'a line that is not an #entry',
    (at: number, length: number) =>
      tableOf(at, length).replace('#entry', 'entry'),
    "R#0#chr's unit table holds entry\tR#0#chr",
  ],
  [
    'a half line',
    (at: number, length: number) => tableOf(at, length).slice(0, -3),
    "R#0#chr's table ends partway through a line",
  ],
  [
    'an #entry line with a tag',
    (at: number, length: number) =>
      tableOf(at, length).replace(/\n$/, '\tzz:i:1\n'),
    "R#0#chr's unit table holds #entry\tR#0#chr",
  ],
  [
    'a #unit line with a tag',
    (at: number, length: number) =>
      tableOf(at, length).replace('table\t.\n', 'table\t.\tzz:i:1\n'),
    'of the walk file, starts "#unit\\tR#0#chr\\t0\\t4\\ttable\\t.\\tzz:i:1"',
  ],
])('a table with %s is refused', async (_, table, message) => {
  await expect(
    oneTile(bgzf(TILE_HEAD + REF_ROW), table).units('R#0#chr'),
  ).rejects.toThrow(message)
})

test("a reader that cannot learn the data's length skips the directory's bytes check", async () => {
  const path = `${DATA}header.GRCh38.walks.gz`
  const data = new LocalFile(path)
  data.stat = () => Promise.resolve({ size: 0 })
  const directory = (
    await new LocalFile(`${path}.idx`).readFile('utf8')
  ).replace(/\tbytes:i:\d+/, '\tbytes:i:12')
  const walks = WalkReader.fromFilehandles({
    data,
    directory: new BlobFile(new Blob([directory])),
  })
  expect((await walks.header()).tile).toBe(1000)
})

// A unit of a kind this reader does not know, listed in the table between
// the tile and the table, is a later addition: the cut reads around it, and
// its ref row of node 2 joins nothing
test('a unit of a kind the reader does not know is skipped', async () => {
  const unitBytes = (blocks: Uint8Array[]) =>
    blocks.reduce((sum, block) => sum + block.length, 0)
  const walks = '#walks\tformat:i:8\ttile:i:16'
  const header = bgzf(`${walks}\n#reference\tR\n`)
  const tile = bgzf(
    '#unit\tR#0#chr\t0\t4\ttile\t0\nR#0#chr\t0\t4\tref\tR#0#chr\t0\t0\t0\t1\t2\t4\t.\n',
  )
  const sketch = bgzf(
    '#unit\tR#0#chr\t0\t4\tsketch\t0\nR#0#chr\t0\t4\tref\tR#0#chr\t0\t0\t0\t1\t4\t4\t.\n',
  )
  const tileAt = unitBytes(header)
  const sketchAt = tileAt + unitBytes(tile)
  const tableAt = sketchAt + unitBytes(sketch)
  const table = bgzf(
    [
      '#unit\tR#0#chr\t0\t4\ttable\t.',
      `#entry\tR#0#chr\t0\t4\ttile\t0\t${tileAt}\t${unitBytes(tile)}`,
      `#entry\tR#0#chr\t0\t4\tsketch\t0\t${sketchAt}\t${unitBytes(sketch)}`,
      '',
    ].join('\n'),
  )
  const parts = [...header, ...tile, ...sketch, ...table]
  const directory = `${walks}\tbytes:i:${unitBytes(parts)}\n#header\t0\t${tileAt}\nR#0#chr\t${tableAt}\t${unitBytes(table)}\n`
  const walkFile = WalkReader.fromFilehandles({
    data: new BlobFile(new Blob(parts)),
    directory: new BlobFile(new Blob([directory])),
  })
  const read = vi.spyOn(BlobFile.prototype, 'read')
  try {
    const { fragments, kept } = await walkFile.cut(
      'R#0#chr',
      { start: 0, end: 4 },
      { context: 0 },
    )
    expect(fragments.map(fragment)).toEqual([['R#0#chr', 0, 4, '>1', 0, 0]])
    expect([...kept]).toEqual([1])
    expect(
      read.mock.calls.filter(
        ([length, position = 0]) =>
          position <= sketchAt && sketchAt < position + length,
      ),
    ).toEqual([])
  } finally {
    read.mockRestore()
  }
})

// `stays` with a tag the reader does not know on its header lines and every
// row, before or after an SQ:Z: tag, and `farAttach` likewise on its
// departure lines: each cuts as the rows without the tag do
const tagged = (lines: string[]) =>
  lines.map((line, k) =>
    line.startsWith('#walks')
      ? `${line} xx:i:1`
      : line.startsWith('#')
        ? `${line} xx:Z:a`
        : line.endsWith(' SQ:Z:AAAA')
          ? line.replace(' SQ:Z:AAAA', ' xx:Z:a SQ:Z:AAAA')
          : `${line} ${k % 2 ? 'xx:i:1' : 'yy:Z:b xx:i:1'}`,
  )

test('tags a reader does not know are ignored on any row and header line', async () => {
  const withBases = stays.map(line =>
    line.includes(' nodes ') ? `${line} SQ:Z:AAAA` : line,
  )
  const cut = async (lines: string[]) => {
    const { graph, kept, fragments, departures } = await inline(lines).cut(
      'R#0#chr',
      { start: 4, end: 20 },
      { context: 0, sequences: true },
    )
    return {
      gfa: graphTablesGFA(graph.tables(kept, fragments, departures)),
      bases: [...graph.sequences].sort(([a], [b]) => a - b),
      departures: departures.map(departure),
    }
  }
  expect(tagged(withBases)).not.toEqual(withBases)
  const plain = await cut(withBases)
  expect(plain.bases).toEqual([
    [2, 'AAAA'],
    [6, 'AAAA'],
  ])
  expect(await cut(tagged(withBases))).toEqual(plain)
  expect(await inline(tagged(stays)).header()).toEqual(
    await inline(stays).header(),
  )
  const attached = async (lines: string[]) => {
    const { fragments, departures } = await inline(lines).cut(
      'R#0#chr',
      { start: 20, end: 28 },
      { context: 0 },
    )
    return [fragments.map(fragment), departures.map(departure)]
  }
  expect(await attached(tagged(farAttach))).toEqual(await attached(farAttach))
})

// `stays` with the #walks line's graph:Z: and leftout:i: and without: the
// header returns them where the file states them, and both cut alike
test('a file that states its graph and the lines it leaves out cuts as one that does not', async () => {
  const stated = stays.map(line =>
    line.startsWith('#walks')
      ? `${line} graph:Z:6f8a0816aa66e2d6 leftout:i:1`
      : line,
  )
  const cut = async (lines: string[]) => {
    const { graph, kept, fragments, departures } = await inline(lines).cut(
      'R#0#chr',
      { start: 4, end: 20 },
      { context: 0 },
    )
    return graphTablesGFA(graph.tables(kept, fragments, departures))
  }
  const header = await inline(stated).header()
  expect(header).toMatchObject({ graph: '6f8a0816aa66e2d6', leftout: 1 })
  const plain = await inline(stays).header()
  expect(plain.graph).toBeUndefined()
  expect(plain.leftout).toBeUndefined()
  expect({ ...plain, graph: header.graph, leftout: 1 }).toEqual(header)
  expect(await cut(stated)).toEqual(await cut(stays))
})

// test_data/extension: one tile whose two walks the writer abridged, B#1#c
// with 50 and 10 bp off R, as written (plain) and with one later addition
// each: a tag on every row, a unit of a kind this reader does not know, and
// a header line whose first word it does not know
const extended = async (name: string, expand?: (name: string) => boolean) => {
  const walks = reader(`extension/${name}.R`)
  const cut = await walks.cut(
    'R#0#c',
    { start: 0, end: 36 },
    { context: 0, expand },
  )
  const { graph, kept, fragments, departures, gaps, placeholders } = cut
  return {
    header: await walks.header(),
    fragments,
    departures,
    gaps,
    placeholders,
    kept: [...kept],
    tables: graph.tables(kept, fragments, departures, gaps),
  }
}

test.each(
  ['tags', 'unit-kind', 'header-line'].flatMap(name => [
    { name, cut: 'unloaded', expand: undefined },
    { name, cut: 'expanded', expand: () => true },
  ]),
)(
  'a file with a later addition cuts as the file without it: $name, $cut',
  async ({ name, expand }) => {
    const plain = await extended('plain', expand)
    expect(plain.fragments.map(f => f.name).sort()).toEqual([
      'A#1#c',
      'B#1#c',
      'R#0#c',
    ])
    expect(plain.placeholders).toBe(expand ? 0 : 1)
    expect(await extended(name, expand)).toEqual(plain)
  },
)

// Every row kind, with a record: K takes a route of the tile, and H is
// abridged, with its detail and blocks units
const everyKind = [
  '#walks format:i:8 tile:i:20',
  '#reference R',
  '#haplotype H#1',
  '#haplotype K#1',
  'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 5,5,5,5 .',
  'R#0#chr 0 nodes 9 3',
  'R#0#chr 0 route 2 16 -10',
  'R#0#chr 0 walk K#1#k 0 0 0 3 2 0',
  'R#0#chr 0 abridged H#1#h 0 0 0 3 18 0',
  'R#0#chr 0:0 10 20 nodes 5 8',
  'R#0#chr 0:0 10 20 route 6 4 -2',
  'R#0#chr 0:0 10 20 walk H#1#h 0 0 0 3 6 0',
  'R#0#chr 0 10 20 blocks H#1#h 0 0 10+5,8,0+5',
  'R#0#chr 0 20 departure K#1#k 5 3 0 5 + 1 15 + 1 0',
]
const everyKindCut = (lines: string[], expand?: (name: string) => boolean) =>
  inline(lines).cut('R#0#chr', { start: 0, end: 20 }, { context: 0, expand })
const withRow = (kind: string, edit: (line: string) => string, from = 0) =>
  everyKind.map((line, k) =>
    k === everyKind.findIndex((l, i) => i >= from && l.includes(` ${kind} `))
      ? edit(line)
      : line,
  )
const withColumn = (kind: string, column: string, from = 0) =>
  withRow(kind, line => `${line} ${column}`, from)

const rowKinds = [
  { kind: 'ref', unit: 'tile' },
  { kind: 'nodes', unit: 'tile' },
  { kind: 'route', unit: 'tile' },
  { kind: 'walk', unit: 'tile' },
  { kind: 'abridged', unit: 'tile' },
  { kind: 'blocks', unit: 'blocks' },
  { kind: 'departure', unit: 'departures' },
  { kind: 'nodes', unit: 'detail', from: 9, expand: () => true },
  { kind: 'route', unit: 'detail', from: 9, expand: () => true },
  { kind: 'walk', unit: 'detail', from: 9, expand: () => true },
]

test.each(rowKinds)(
  'a column after the fixed columns of a $kind row of a $unit unit is ignored where it is a tag, and refused where it is not',
  async ({ kind, unit, from, expand }) => {
    const plain = await everyKindCut(everyKind, expand)
    const withTag = await everyKindCut(
      withColumn(kind, 'zz:i:1 Yy9:Z:a:b z:f:', from),
      expand,
    )
    expect(withTag.fragments).toEqual(plain.fragments)
    expect(withTag.departures).toEqual(plain.departures)
    const refused: unknown = await everyKindCut(
      withColumn(kind, 'junk', from),
      expand,
    ).catch((e: unknown) => e)
    expect(refused).toBeInstanceOf(WalkFormatError)
    expect(refused).toHaveProperty(
      'message',
      expect.stringMatching(
        new RegExp(
          `R#0#chr's ${unit} .* holds a row with a column after its fixed ones that is no tag, "junk"`,
        ),
      ),
    )
  },
)

test.each(rowKinds)(
  'a $kind row of a $unit unit without its last fixed column is a WalkRowError',
  async ({ kind, unit, from, expand }) => {
    const short = withRow(
      kind,
      line => line.slice(0, line.lastIndexOf(' ')),
      from,
    )
    const refused: unknown = await everyKindCut(short, expand).catch(
      (e: unknown) => e,
    )
    expect(refused).toBeInstanceOf(WalkRowError)
    expect(refused).toHaveProperty(
      'message',
      expect.stringMatching(
        new RegExp(
          `R#0#chr's ${unit} .* holds a row of fewer columns than the \\d+ its kind holds`,
        ),
      ),
    )
  },
)

test.each(['zz:1', '9z:i:1', ':i:1', 'zz:ii:1', 'z_z:i:1', 'zz'])(
  'a column after a row’s fixed columns that is no well-formed tag is refused: %s',
  async column => {
    await expect(everyKindCut(withColumn('ref', column))).rejects.toThrow(
      'that is no tag',
    )
  },
)

test('a directory line of a later kind is ignored, and a sequence line with a column more is refused', async () => {
  const { data, directory } = walkFileOf(everyKind)
  const text = await directory.readFile('utf8')
  const withLine = (edit: (text: string) => string) =>
    WalkReader.fromFilehandles({
      data,
      directory: new BlobFile(new Blob([edit(text)])),
    })
  expect(
    await withLine(t => `${t}#future\tx\ty\n`).referenceSequences(),
  ).toEqual(['R#0#chr'])
  await expect(
    withLine(t => t.replace(/(R#0#chr\t\d+\t\d+)\n/, '$1\tzz:i:1\n')).units(
      'R#0#chr',
    ),
  ).rejects.toThrow('A walk file directory line holds no byte range')
})

test("a tag after an abridged row's detail column names no other detail unit", async () => {
  const { fragments } = await inline([
    '#walks format:i:8 tile:i:20',
    '#reference R',
    '#haplotype H#1',
    'R#0#chr 0 ref R#0#chr 0 0 0 4 2,2,2,2 5,5,5,5 .',
    'R#0#chr 0 abridged H#1#h 0 0 0 3 18 0 xx:i:1',
    'R#0#chr 0:0 10 20 nodes 5 8',
    'R#0#chr 0:0 10 20 route 6 4 -2',
    'R#0#chr 0:0 10 20 walk H#1#h 0 0 0 3 6 0',
    'R#0#chr 0 10 20 blocks H#1#h 0 0 10+5,8,0+5 xx:i:1',
  ]).cut('R#0#chr', { start: 0, end: 20 }, { context: 0 })
  expect(
    fragments.flatMap(f => (f.abridged ?? []).map(a => a.detail.name)),
  ).toEqual(['0:0'])
})

// R's chr0 is two nodes of 0 bp: its tile states [0, 0), so a window from 0
// on reads it never, and the sequence's length is 0
test('a unit of a sequence of 0 bp, stating [0, 0), is read by no window', async () => {
  const served: number[] = []
  const walks = inline(
    [
      '#walks format:i:8 tile:i:16',
      '#reference R',
      '#haplotype H#1',
      'R#0#chr0 0 ref R#0#chr0 0 0 0 2 2,2 0,0 .',
      'R#0#chr0 0 walk H#1#h 0 0 0 2 2 .',
    ],
    false,
    served,
  )
  expect(
    (await walks.units('R#0#chr0')).map(u => [u.kind, u.start, u.end]),
  ).toEqual([['tile', 0, 0]])
  expect(await walks.sequenceLength('R#0#chr0')).toBe(0)
  const { fragments } = await walks.cut(
    'R#0#chr0',
    { start: 0, end: 16 },
    { context: 0 },
  )
  expect(fragments).toEqual([])
  expect(served).toEqual([])
})
