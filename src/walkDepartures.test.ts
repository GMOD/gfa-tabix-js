import { LocalFile } from 'generic-filehandle2'

import { AmbiguousSequenceError } from './refNames.ts'
import { describeForm, structuralForms } from './structuralForms.ts'
import { walkFileFrom, walkFileOf } from './testing.ts'
import { RegionTooLargeError, panSNHaplotypeOrSample } from './util.ts'
import { walkDepartures } from './walkDepartures.ts'
import { DepartureReadError, WalkFile } from './walkFile.ts'
import { walkCut } from './walkRows.ts'

import type { WalkDeparture } from './walkDeparture.ts'
import type { HaplotypeDepartures } from './walkDepartures.ts'
import type { WalkNode, WalkRun } from './walkRows.ts'

const R = 'GRCh38#0#chr1'

// R as nodes 1-6 of 100, 500, 100, 1000, 100 and 100 bp from 0, and nodes 50
// and up off R
const refLengths = [100, 500, 100, 1000, 100, 100]
const offR: Record<number, number> = { 50: 300, 51: 20, 52: 2000 }

function graphNodes() {
  const nodes = new Map<number, WalkNode>()
  let start = 0
  refLengths.forEach((length, i) => {
    nodes.set(i + 1, {
      refName: R,
      start,
      end: start + length,
      rank: 0,
      passes: 1,
      onWindow: true,
    })
    start += length
  })
  for (const [id, length] of Object.entries(offR)) {
    nodes.set(+id, {
      refName: 'off',
      start: 0,
      end: length,
      rank: 1,
      passes: 0,
      onWindow: false,
    })
  }
  return nodes
}

const refStartOf = (nodes: Map<number, WalkNode>, id: number) => {
  const node = nodes.get(id)
  return node?.rank === 0 ? node.start : NaN
}

// A walk's run through the cut: `+3` steps forward through node 3, `-3`
// against it, from `hapStart` on its contig. R holds every node forward and
// once unless `against` or `refStarts` say otherwise
function walk(
  name: string,
  steps: string,
  more: {
    hapStart?: number
    against?: number[]
    refStarts?: number[]
  } = {},
): WalkRun {
  const nodes = graphNodes()
  const parsed = steps
    .split(' ')
    .map(s => [+s.slice(1), s.startsWith('-')] as const)
  return {
    name,
    fragStart: 0,
    hapOffset: more.hapStart ?? 0,
    firstPiece: 0,
    pieceStarts: [],
    pieceChunks: [0],
    ids: Int32Array.from(parsed, ([id]) => id),
    rev: Uint8Array.from(parsed, ([, rev]) => (rev ? 1 : 0)),
    refStarts:
      more.refStarts === undefined
        ? Float64Array.from(parsed, ([id]) => refStartOf(nodes, id))
        : Float64Array.from(more.refStarts),
    against:
      more.against === undefined
        ? Uint8Array.from(parsed, ([, rev]) => (rev ? 1 : 0))
        : Uint8Array.from(more.against),
  }
}

// A record of `name`'s departure: `ownBp` of its own from `hapStart` on its
// contig, between the leave and rejoin breakpoints, `+700` along R at 700
// and `-600` against it
function record(
  name: string,
  hapStart: number,
  ownBp: number,
  leave: string,
  rejoin: string,
  missingBp = 0,
): WalkDeparture {
  return {
    name,
    hapStart,
    ownBp,
    missingBp,
    leavePos: +leave.slice(1),
    leaveReversed: leave.startsWith('-'),
    leavePasses: 1,
    rejoinPos: +rejoin.slice(1),
    rejoinReversed: rejoin.startsWith('-'),
    rejoinPasses: 1,
    fromSteps: false,
  }
}

// The cut of these runs over the whole of R, the departures inside each run
// computed from its steps and `records` those between runs, read for forms
function departuresOf(runs: WalkRun[], records: WalkDeparture[] = []) {
  const nodes = graphNodes()
  const { fragments, departures } = walkCut(runs, nodes, R, records)
  return walkDepartures({ fragments, departures, nodes }, R)
}

const plain = {
  invertedBp: 0,
  repeatedBp: 0,
}

test('a walk along R departs nowhere and covers what it reaches', () => {
  expect(departuresOf([walk('A#1#c', '+1 +2 +3')])).toEqual([
    { haplotype: 'A#1', covered: [[0, 700]], departures: [] },
  ])
})

test('a skipped stretch of R is a deletion with no ALT', () => {
  const [a] = departuresOf([walk('A#1#c', '+1 +3')])
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 600, altLength: 0, ...plain },
  ])
})

test('steps off R between abutting R steps insert at a point', () => {
  const [a] = departuresOf([walk('A#1#c', '+1 +50 +2')])
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 100, altLength: 300, ...plain },
  ])
})

test('steps off R in place of a stretch of R replace it', () => {
  const [a] = departuresOf([walk('A#1#c', '+1 +50 +3')])
  expect(a!.departures).toEqual([
    { refStart: 100, refEnd: 600, altLength: 300, ...plain },
  ])
})

test('a stretch of R walked against it is ALT and inverted', () => {
  const [a] = departuresOf([walk('A#1#c', '+2 -3 +4')])
  expect(a!.departures).toEqual([
    {
      refStart: 600,
      refEnd: 700,
      altLength: 100,
      invertedBp: 100,
      repeatedBp: 0,
    },
  ])
})

test('a step back along R is a tandem copy at the junction', () => {
  const [a] = departuresOf([walk('A#1#c', '+1 +2 +3 +2 +3 +4')])
  expect(a!.departures).toEqual([
    {
      refStart: 700,
      refEnd: 700,
      altLength: 600,
      invertedBp: 0,
      repeatedBp: 600,
      copied: [100, 700],
    },
  ])
  expect(a!.covered).toEqual([[0, 1700]])
})

test('a walk reading R mostly reversed is read the other way round', () => {
  const [a] = departuresOf([walk('A#1#c', '-4 -50 -2 -1')])
  expect(a!.departures).toEqual([
    { refStart: 600, refEnd: 700, altLength: 300, ...plain },
  ])
  expect(a!.covered).toEqual([[0, 1700]])
})

const inversion = {
  refStart: 100,
  refEnd: 1700,
  altLength: 1600,
  invertedBp: 1600,
  repeatedBp: 0,
}

test('a walk forward at its flanks reads an inversion over most of the window as inverted', () => {
  const haplotypes = departuresOf([
    walk('A#1#c', '+1 -4 -3 -2 +5'),
    walk('B#1#c', '+1 +2 +3 +4 +5'),
    walk('C#1#c', '-5 -4 -3 -2 -1'),
    walk('D#1#c', '-5 +2 +3 +4 -1'),
  ])
  expect(haplotypes).toEqual([
    { haplotype: 'A#1', covered: [[0, 1800]], departures: [inversion] },
    { haplotype: 'B#1', covered: [[0, 1800]], departures: [] },
    { haplotype: 'C#1', covered: [[0, 1800]], departures: [] },
    { haplotype: 'D#1', covered: [[0, 1800]], departures: [inversion] },
  ])
  const { forms } = structuralForms(haplotypes, {
    window: { start: 0, end: 1800 },
  })
  expect(forms.map(f => [describeForm(f), f.members])).toEqual([
    ['ref', ['B#1', 'C#1']],
    ['inv 1.6kb', ['A#1', 'D#1']],
  ])
})

test('one flank is enough to tell an inversion from a walk stored back to front', () => {
  const [a, b, c] = departuresOf([
    walk('A#1#c', '+1 -4 -3 -2'),
    walk('B#1#c', '+2 +3 +4 -1'),
    walk('C#1#c', '-4 -3 -2 +5'),
  ])
  expect(a).toEqual({ haplotype: 'A#1', covered: [[0, 100]], departures: [] })
  expect(b).toEqual({ haplotype: 'B#1', covered: [[0, 100]], departures: [] })
  expect(c).toEqual({
    haplotype: 'C#1',
    covered: [[1700, 1800]],
    departures: [],
  })
})

test('a walk wholly inside an inversion reads as one stored back to front', () => {
  expect(departuresOf([walk('A#1#c', '-4 -3 -2')])).toEqual([
    { haplotype: 'A#1', covered: [[100, 1700]], departures: [] },
  ])
})

test('a walk forward throughout that ends before where it began is not turned around', () => {
  const [a] = departuresOf([walk('A#1#c', '+3 +4 +1 +2')])
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[600, 1700]],
    departures: [
      {
        refStart: 1700,
        refEnd: 1700,
        altLength: 1700,
        invertedBp: 0,
        repeatedBp: 1700,
        copied: [0, 1700],
      },
    ],
  })
})

test('steps off R past either end of a walk are no departure', () => {
  expect(departuresOf([walk('A#1#c', '+51 +1 +2 +52')])).toEqual([
    { haplotype: 'A#1', covered: [[0, 600]], departures: [] },
  ])
  expect(departuresOf([walk('B#1#c', '+50 +52')])).toEqual([])
})

test("a haplotype's paths join into one list, the reference's are left out", () => {
  const result = departuresOf([
    walk(R, '+1 +2 +3 +4 +5 +6'),
    walk('A#1#c1', '+4 +6'),
    walk('A#1#c2', '+1 +50 +2'),
    walk('A#2#c1', '+1 +2'),
  ])
  expect(
    result.map(h => [h.haplotype, h.covered, h.departures.length]),
  ).toEqual([
    [
      'A#1',
      [
        [0, 600],
        [700, 1900],
      ],
      2,
    ],
    ['A#2', [[0, 600]], 0],
  ])
  expect(result[0]!.departures.map(d => [d.refStart, d.refEnd])).toEqual([
    [100, 100],
    [1700, 1800],
  ])
})

test('R that two paths of a haplotype cover, or a junction between them, has no call', () => {
  const [overlapping, abutting] = departuresOf([
    walk('A#1#c1', '+1 +2 +3 +4'),
    walk('A#1#c2', '+3 +4 +5'),
    walk('A#2#c1', '+1 +2'),
    walk('A#2#c2', '+3 +4'),
  ])
  expect(overlapping!.covered).toEqual([
    [0, 600],
    [1700, 1800],
  ])
  expect(abutting!.covered).toEqual([
    [0, 600],
    [600, 1700],
  ])
})

test('a departure with no step on the window is read from its record alone', () => {
  expect(departuresOf([], [record('A#1#c', 500, 0, '+100', '+1700')])).toEqual([
    {
      haplotype: 'A#1',
      covered: [[100, 1700]],
      departures: [{ refStart: 100, refEnd: 1700, altLength: 0, ...plain }],
    },
  ])
  expect(departuresOf([], [record('B#1#c', 500, 0, '-700', '-100')])).toEqual([
    {
      haplotype: 'B#1',
      covered: [[100, 700]],
      departures: [{ refStart: 100, refEnd: 700, altLength: 0, ...plain }],
    },
  ])
})

test("a departure across a contig's W lines counts the sequence the graph lacks as its own", () => {
  const [a] = departuresOf(
    [walk('A#1#c', '+1 +2'), walk('A#1#c', '+3 +4', { hapStart: 5600 })],
    [record('A#1#c', 600, 5000, '+600', '+600', 5000)],
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 1700]],
    departures: [{ refStart: 600, refEnd: 600, altLength: 5000, ...plain }],
  })
})

test("a contig's second visit to the window gives the R both reach no call", () => {
  const [a] = departuresOf([
    walk('A#1#c', '+1 +2 +3'),
    walk('A#1#c', '+2 +3 +4', { hapStart: 5000 }),
  ])
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [
      [0, 100],
      [700, 1700],
    ],
    departures: [],
  })
})

test('a walk along R through a node R holds reversed stays on R', () => {
  expect(
    departuresOf([walk('A#1#c', '-2 +3 +4', { against: [0, 0, 0] })]),
  ).toEqual([{ haplotype: 'A#1', covered: [[100, 1700]], departures: [] }])
  const [a] = departuresOf([walk('A#1#c', '+1 +50 -3', { against: [0, 0, 0] })])
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[0, 700]],
    departures: [{ refStart: 100, refEnd: 600, altLength: 300, ...plain }],
  })
})

test("a fragment's ends lie where their pass does, where R passes a node again", () => {
  const [a] = departuresOf(
    [walk('A#1#c', '+2 +4', { hapStart: 400, refStarts: [700, 1200] })],
    [record('A#1#c', 100, 300, '+100', '+700')],
  )
  expect(a).toEqual({
    haplotype: 'A#1',
    covered: [[100, 2200]],
    departures: [{ refStart: 100, refEnd: 700, altLength: 300, ...plain }],
  })
})

// A graph of P lines named `sample#contig`, with no haplotype field: the
// header lists H#c and H#d under H, and R#chr2 is the reference's
test('paths named with no haplotype field are their sample’s, so a forms cut lists the header’s haplotypes and no other', async () => {
  const walks = walkFileFrom(
    walkFileOf([
      '#walks format:i:8 tile:i:100',
      '#reference R',
      '#haplotype H',
      'R#chr1 0 ref R#chr1 0 0 0 3 2,2,2 6,6,6 .',
      'R#chr1 0 walk H#c 0 0 0 3 2 .',
      'R#chr1 0 walk H#d 0 0 0 1 6 .',
      'R#chr1 0 walk R#chr2 0 0 0 2 2 .',
    ]),
  )
  expect((await walks.getHeader()).haplotypes).toEqual(['H'])
  const forms = await walks.getFormsForRange('chr1', 0, 18)
  expect(forms.map(h => [h.haplotype, h.covered])).toEqual([['H', [[0, 12]]]])
  expect(panSNHaplotypeOrSample('H#c')).toBe('H')
  expect(panSNHaplotypeOrSample('H#1#c')).toBe('H#1')
  expect(panSNHaplotypeOrSample('H#1#c#0')).toBe('H#1')
  expect(panSNHaplotypeOrSample('chr1')).toBe('chr1')
})

describe('a forms cut', () => {
  const DATA = new URL('../test_data/', import.meta.url).pathname
  const reader = (name: string) => {
    const path = `${DATA}${name}.walks.gz`
    return new WalkFile({
      walksFilehandle: new LocalFile(path),
      idxFilehandle: new LocalFile(`${path}.idx`),
    })
  }
  const region = { refName: 'chr', start: 0, end: 30 }
  const byHaplotype = (hs: HaplotypeDepartures[]) =>
    Object.fromEntries(
      hs
        .filter(h => h.covered.length > 0 || h.departures.length > 0)
        .map(h => [h.haplotype, h]),
    )
  const net = (h: HaplotypeDepartures) =>
    h.departures.reduce(
      (sum, d) => sum + d.altLength - (d.refEnd - d.refStart),
      0,
    )
  // the departures of a cut that reads every abridged walk in full
  const loaded = async (
    name: string,
    window: { refName: string; start: number; end: number },
  ) => {
    const walks = reader(name)
    const { references } = await walks.getHeader()
    const refName = await walks.getSequenceName(references[0]!, window.refName)
    const cut = await walks.getCutForRange(refName!, window.start, window.end, {
      links: 'none',
      expand: () => true,
    })
    return byHaplotype(
      walkDepartures(
        {
          fragments: cut.fragments,
          departures: cut.departures,
          nodes: cut.graph.nodes,
        },
        refName!,
      ),
    )
  }
  const unloaded = async (
    name: string,
    window: { refName: string; start: number; end: number },
  ) =>
    byHaplotype(
      await reader(name).getFormsForRange(
        window.refName,
        window.start,
        window.end,
      ),
    )

  test('reads unloaded pieces as the unabridged file gives them', async () => {
    const whole = await unloaded('unabridged.R', region)
    expect(await unloaded('abridged.R', region)).toEqual(whole)
    expect(Object.keys(whole).sort()).toEqual(['D#1', 'H#1', 'I#1', 'Q#0'])
  })

  test('reads no detail unit', async () => {
    const path = `${DATA}abridged.R.walks.gz`
    const data = new LocalFile(path)
    const details = new Set(
      (await reader('abridged.R').getUnits('R#0#chr'))
        .filter(u => u.kind === 'detail')
        .map(u => u.offset),
    )
    const read = data.read.bind(data)
    let detailReads = 0
    data.read = (length, position = 0, opts = {}) => {
      detailReads += details.has(position) ? 1 : 0
      return read(length, position, opts)
    }
    await new WalkFile({
      walksFilehandle: data,
      idxFilehandle: new LocalFile(`${path}.idx`),
    }).getFormsForRange(region.refName, region.start, region.end)
    expect(details.size).toBeGreaterThan(0)
    expect(detailReads).toBe(0)
  })

  test('keeps to FORMS_STEP_BUDGET steps at least, whatever smaller budget its caller gives', async () => {
    expect(
      await reader('abridged.R').getFormsForRange(
        region.refName,
        region.start,
        region.end,
        { stepBudget: 3 },
      ),
    ).toEqual(
      await reader('abridged.R').getFormsForRange(
        region.refName,
        region.start,
        region.end,
      ),
    )
  })

  test('reads an unloaded walk as its loaded one when a piece of it lies in another tile', async () => {
    // X#1's pieces 0 and 2 are abridged in tile 0, and piece 1 between them
    // walks R in tile 10
    const window = { refName: 'c', start: 0, end: 36 }
    const whole = await loaded('interleaved.R', window)
    expect(net(whole['X#1']!)).toBe(84)
    expect(await unloaded('interleaved.R', window)).toEqual(whole)
  })

  // X#1 and Z#1 leave R from an abridged piece whose last steps are their
  // own, V#1 and T#1 rejoin across W lines into one whose first steps are
  // their own, Q#1 through a W line no piece holds, P#1 through a step on
  // another reference sequence, and S#1 rejoins reversed into a piece
  // against R
  test('reads each haplotype as its walks loaded when a departure covers own steps of an unloaded piece, or a piece runs against R', async () => {
    const window = { refName: 'c', start: 0, end: 36 }
    const whole = await loaded('ownends.R', window)
    expect(await unloaded('ownends.R', window)).toEqual(whole)
    expect(
      Object.keys(whole)
        .sort()
        .map(haplotype => [haplotype, net(whole[haplotype]!)]),
    ).toEqual([
      ['P#1', 70_110],
      ['Q#1', 70_028],
      ['S#1', 0],
      ['T#1', 70_008],
      ['V#1', 28],
      ['X#1', 18],
      ['Y#1', 40],
      ['Z#1', 10],
    ])
  })

  // L#1 starts with 30 bp of its own before R's 3 and M#1 ends with them
  // after it, at their contigs' ends, where no departure record covers them;
  // their pieces there are abridged
  test("reads each haplotype as its walks loaded when an unloaded piece holds its contig's first or last steps", async () => {
    const window = { refName: 'c', start: 0, end: 36 }
    const whole = await loaded('contigends.R', window)
    expect(whole).toEqual({
      'L#1': { haplotype: 'L#1', covered: [[12, 36]], departures: [] },
      'M#1': { haplotype: 'M#1', covered: [[0, 18]], departures: [] },
    })
    expect(await unloaded('contigends.R', window)).toEqual(whole)
    const fragments = async (expand?: () => boolean) => {
      const cut = await reader('contigends.R').getCutForRange(
        'R#0#c',
        window.start,
        window.end,
        { expand },
      )
      return cut.fragments.map(f => [
        f.name,
        f.hapStart,
        f.hapEnd,
        f.trimmedBefore,
        f.trimmedAfter,
        [...f.refStarts].join(','),
      ])
    }
    expect(await fragments()).toEqual(await fragments(() => true))
    expect(await fragments()).toEqual([
      ['R#0#c', 0, 36, 0, 0, '0,6,12,18,24,30'],
      ['M#1#c', 0, 18, 0, 30, '0,6,12'],
      ['L#1#c', 30, 54, 30, 0, '12,18,24,30'],
    ])
  })

  // R over five 2000 bp nodes in tile 0 and 1000 bp nodes from 10000
  // (drawnedge.R.gfa). H#1 walks R to 13000, runs 600 bp of its own, rejoins
  // R at 4000 and walks it again to 14000; its pieces in tile 10000 are
  // abridged. Over 10500-13500 the drawn stretch is 8000-15000, so the
  // piece after the rejoin starts outside it and the cut trims its first
  // steps. G#1's one abridged piece walks R from 10000 to 20000 and drops
  // 1000 bp of it for 300 of its own at 16000, past the drawn stretch
  test("chains a departure record's rejoin to an unloaded piece whose first steps lie outside the drawn stretch, and trims a piece reaching past it, as the loaded read does", async () => {
    const window = { refName: 'chr', start: 10_500, end: 13_500 }
    const whole = await loaded('drawnedge.R', window)
    expect(whole['H#1']).toEqual({
      haplotype: 'H#1',
      covered: [[8000, 14_000]],
      departures: [
        {
          refStart: 13_000,
          refEnd: 13_000,
          altLength: 9600,
          invertedBp: 0,
          repeatedBp: 9000,
          copied: [4000, 13_000],
        },
      ],
    })
    expect(whole['G#1']).toEqual({
      haplotype: 'G#1',
      covered: [[8000, 15_000]],
      departures: [],
    })
    expect(await unloaded('drawnedge.R', window)).toEqual(whole)
  })

  test('throws RegionTooLargeError past its byte budget', async () => {
    await expect(
      reader('abridged.R').getFormsForRange(
        region.refName,
        region.start,
        region.end,
        { byteBudget: 1 },
      ),
    ).rejects.toBeInstanceOf(RegionTooLargeError)
    expect(
      await reader('abridged.R').getFormsForRange(
        region.refName,
        region.start,
        region.end,
        { byteBudget: 1e9 },
      ),
    ).toEqual(
      await reader('abridged.R').getFormsForRange(
        region.refName,
        region.start,
        region.end,
      ),
    )
  })

  test('throws on a contig both haplotypes of a diploid reference sample hold', async () => {
    const window = { refName: 'chr1', start: 0, end: 40 }
    await expect(
      reader('dip.HG002').getFormsForRange(
        window.refName,
        window.start,
        window.end,
      ),
    ).rejects.toThrow(AmbiguousSequenceError)
    const named = await reader('dip.HG002').getFormsForRange(
      'HG002#1#chr1',
      window.start,
      window.end,
    )
    expect(named.map(h => h.haplotype)).toEqual(['HG00097#1'])
  })

  test('throws DepartureReadError past its departure byte budget', async () => {
    const window = { refName: 'chr1', start: 50, end: 55 }
    await expect(
      reader('span.GRCh38').getFormsForRange(
        window.refName,
        window.start,
        window.end,
        { departureByteBudget: 1 },
      ),
    ).rejects.toBeInstanceOf(DepartureReadError)
    expect(
      await reader('span.GRCh38').getFormsForRange(
        window.refName,
        window.start,
        window.end,
        {
          departureByteBudget: 1e9,
        },
      ),
    ).toEqual(
      await reader('span.GRCh38').getFormsForRange(
        window.refName,
        window.start,
        window.end,
      ),
    )
  })
})
