import { readdirSync } from 'node:fs'

import { LocalFile } from 'generic-filehandle2'

import { MissingSequenceError } from './refNames.ts'
import {
  bytesRead,
  endsAgainstGaps,
  firstDifference,
  readsOf,
  runsBreaks,
  unitsRead,
  walkFileFrom,
  walkFileOf,
} from './testing.ts'
import { panSNHaplotypeOrSample } from './util.ts'
import { WalkContainer } from './walkContainer.ts'
import { departureStretch, squeezer } from './walkDeparture.ts'
import { chainDepartures, walkDepartures } from './walkDepartures.ts'
import { WalkFile } from './walkFile.ts'
import {
  END_KINDS,
  RunsRefusedError,
  fragmentRuns,
  runDepartures,
  runEvents,
} from './walkRuns.ts'

import type { GraphTables } from './util.ts'
import type { Departure } from './walkDepartures.ts'
import type { WalkRuns } from './walkRuns.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

const FIXTURES = readdirSync(DATA)
  .filter(name => name.endsWith('.walks.gz'))
  .map(name => name.replace(/\.walks\.gz$/, ''))
  .sort()

type Cut = Awaited<ReturnType<WalkFile['getCutForRange']>>

const tablesOf = (cut: Cut) => cut.tables()

// A promise's value or its error, whichever it settles to
const settle = <T>(promise: Promise<T>) =>
  promise.then(
    value => ({ value, error: undefined }),
    (error: unknown) => ({ value: undefined, error }),
  )

// Where each piece of each line past its first starts on its contig, by walk
async function pieceStarts(name: string, refName: string) {
  const file = new WalkContainer(
    new LocalFile(`${DATA}${name}.walks.gz`),
    new LocalFile(`${DATA}${name}.walks.gz.idx`),
  )
  const units = (await file.units(refName)).filter(u => u.kind === 'tile')
  const starts = new Map<string, number[]>()
  for (const lines of await file.read(units)) {
    for (const line of lines) {
      const cols = line.split('\t')
      const kind = cols[3] ?? ''
      if (['walk', 'abridged', 'ref'].includes(kind) && Number(cols[7]) > 0) {
        const held = starts.get(cols[4]!) ?? []
        starts.set(cols[4]!, held)
        held.push(Number(cols[6]))
      }
    }
  }
  return starts
}

const CASES = [
  'single-base variant',
  'deletion',
  'inverted stretch',
  'jump back',
  'hole',
  'lead',
  'tail',
  'placeholder',
  'piece join',
  'record between pieces',
  'record at an end',
  'departing walk',
  'contained split',
  'beyond',
  'elsewhere',
  'edge',
  'ends',
] as const

type Case = (typeof CASES)[number]

const squeezeOf = (runs: WalkRuns) =>
  squeezer(
    runs.holes.map(h => ({ ...h, refName: runs.refName })),
    runs.refName,
  )

// The cases one window's runs hold
function casesOf(
  runs: WalkRuns,
  tables: GraphTables,
  starts: Map<string, number[]>,
  contained: boolean,
) {
  const held = new Set<Case>()
  const { walks, stretches } = runs
  if (runs.departing.names.length > 0) {
    held.add('departing walk')
  }
  const { departures } = tables
  const squeeze = squeezeOf(runs)
  walks.names.forEach((name, p) => {
    const from = walks.offsets[p]!
    const to = walks.offsets[p + 1]!
    const ways = new Set<number>()
    for (let k = from; k < to; k++) {
      ways.add(stretches.against[k]!)
      if (k > from && stretches.against[k] === stretches.against[k - 1]) {
        const against = stretches.against[k] === 1
        const prevEnd = against
          ? stretches.starts[k - 1]!
          : stretches.starts[k - 1]! + stretches.bp[k - 1]!
        const begin = against
          ? stretches.starts[k]! + stretches.bp[k]!
          : stretches.starts[k]!
        const sign = against ? -1 : 1
        const gap = sign * (begin - prevEnd)
        const squeezed = sign * (squeeze(begin) - squeeze(prevEnd))
        const own = stretches.offBefore[k]!
        if (own === 1 && gap === 1) {
          held.add('single-base variant')
        }
        if (own === 0 && squeezed > 0) {
          held.add('deletion')
        }
        if (gap !== squeezed) {
          held.add('hole')
        }
        if (gap < 0) {
          held.add('jump back')
        }
      }
      if (stretches.loaded[k] === 0) {
        held.add('placeholder')
      }
    }
    if (ways.size > 1) {
      held.add('inverted stretch')
    }
    if (walks.leadBp[p]! > 0) {
      held.add('lead')
    }
    if (walks.tailBp[p]! > 0) {
      held.add('tail')
    }
    if (walks.leadLoaded[p] === 0 || walks.tailLoaded[p] === 0) {
      held.add('placeholder')
    }
    const start = walks.starts[p]!
    const end = walks.ends[p]!
    if ((starts.get(name) ?? []).some(at => at > start && at < end)) {
      held.add('piece join')
    }
    const breaks = runsBreaks(runs, p)
    departures.names.forEach((d, k) => {
      const leave = departures.hapStarts[k]!
      const rejoin = leave + departures.ownBp[k]!
      const record = departures.fromSteps[k] === 0
      if (d === name && record && leave > start && rejoin < end) {
        held.add('record between pieces')
      }
      const atEnd =
        (departures.leavePos[k] === breaks.last &&
          leave === end - walks.tailBp[p]!) ||
        (departures.rejoinPos[k] === breaks.first &&
          rejoin === start + walks.leadBp[p]!)
      if (d === name && record && atEnd) {
        held.add('record at an end')
      }
    })
    const after = END_KINDS[walks.after.kinds[p]!]!
    held.add(after)
    held.add(END_KINDS[walks.before.kinds[p]!]!)
    if (
      contained &&
      (after === 'beyond' || after === 'elsewhere') &&
      walks.after.ownBp[p]! < walks.tailBp[p]! + walks.trimmedAfter[p]!
    ) {
      held.add('contained split')
    }
  })
  return held
}

const byPlace = (a: Departure, b: Departure) =>
  a.refStart - b.refStart || a.refEnd - b.refEnd

/**
 * walkDepartures of a cut against the departures its runs read, for each
 * haplotype with one walk and one fragment whose every departure in the cut
 * touches that fragment or its rejoining ends, and lies on the drawn
 * stretch where it is between two of its stretches, as a cut lists only
 * those
 */
function departuresAgainstRuns(cut: Cut, runs: WalkRuns, refName: string) {
  const tables = tablesOf(cut)
  const read = walkDepartures(
    {
      fragments: cut.fragments,
      departures: cut.departures,
      nodes: cut.graph.nodes,
    },
    refName,
  )
  const { walks, stretches, drawn } = runs
  let compared = 0
  const mismatches: string[] = []
  for (const h of read) {
    const fragments = walks.names.flatMap((name, p) =>
      panSNHaplotypeOrSample(name) === h.haplotype ? [p] : [],
    )
    const [p] = fragments
    if (fragments.length === 1 && p !== undefined) {
      const name = walks.names[p]!
      let onDrawn = true
      for (let k = walks.offsets[p]! + 1; k < walks.offsets[p + 1]!; k++) {
        const end = (j: number) => stretches.starts[j]! + stretches.bp[j]!
        const leaveReversed = stretches.against[k - 1] === 1
        const rejoinReversed = stretches.against[k] === 1
        const [lo, hi] = departureStretch({
          ownBp: stretches.offBefore[k]!,
          leavePos: leaveReversed ? stretches.starts[k - 1]! : end(k - 1),
          leaveReversed,
          rejoinPos: rejoinReversed ? end(k) : stretches.starts[k]!,
          rejoinReversed,
        })
        onDrawn &&= lo < drawn.end && hi > drawn.start
      }
      const rejoins = (kinds: Uint8Array) =>
        kinds[p]! <= END_KINDS.indexOf('elsewhere')
      const from =
        walks.starts[p]! -
        (rejoins(walks.before.kinds) ? walks.before.ownBp[p]! : 0)
      const to =
        walks.ends[p]! +
        (rejoins(walks.after.kinds) ? walks.after.ownBp[p]! : 0)
      const touching = tables.departures.names.every(
        (d, k) =>
          d !== name ||
          (tables.departures.hapStarts[k]! >= from &&
            tables.departures.hapStarts[k]! + tables.departures.ownBp[k]! <=
              to),
      )
      if (onDrawn && touching) {
        compared++
        const expected = chainDepartures(
          runEvents(runs, p, true),
        ).departures.toSorted(byPlace)
        const found = firstDifference(h.departures, expected)
        if (found !== undefined) {
          mismatches.push(`${h.haplotype}: ${found}`)
        }
      }
    }
  }
  return { compared, mismatches }
}

interface Sweep {
  fixtures: number
  cuts: number
  refused: string[]
  differences: string[]
  cases: Map<Case, number>
  byFixture: Map<string, Set<Case>>
  gaps: Map<string, number>
  gapMismatches: string[]
  departuresCompared: number
  departureMismatches: string[]
  overBytes: string[]
  details: string[]
}

async function sweep(): Promise<Sweep> {
  const out: Sweep = {
    fixtures: 0,
    cuts: 0,
    refused: [],
    differences: [],
    cases: new Map(CASES.map(c => [c, 0])),
    byFixture: new Map(),
    gaps: new Map(),
    gapMismatches: [],
    departuresCompared: 0,
    departureMismatches: [],
    overBytes: [],
    details: [],
  }
  for (const name of FIXTURES) {
    const data = new LocalFile(`${DATA}${name}.walks.gz`)
    const reads = readsOf(data)
    const reader = new WalkFile({
      walksFilehandle: data,
      idxFilehandle: new LocalFile(`${DATA}${name}.walks.gz.idx`),
    })
    const header = await reader.getHeader()
    const fixtureCases = new Set<Case>()
    out.byFixture.set(name, fixtureCases)
    out.fixtures++
    const refused = await reader
      .getRunsForRange((await reader.getReferenceSequences())[0]!, 0, 1)
      .then(
        () => false,
        (error: unknown) => error instanceof RunsRefusedError,
      )
    if (refused) {
      out.refused.push(name)
    }
    for (const refName of refused ? [] : await reader.getReferenceSequences()) {
      const starts = await pieceStarts(name, refName)
      const units = await reader.getUnits(refName)
      const length = await reader.getSequenceLength(refName)
      const regions = [length / 12, length / 3, length].flatMap(bp => {
        const width = Math.max(5, Math.round(bp))
        const step = Math.max(1, Math.round(width / 2))
        const held: { start: number; end: number }[] = []
        for (let start = 0; start < length; start += step) {
          held.push({ start, end: Math.min(length, start + width) })
        }
        return held
      })
      for (const region of regions) {
        for (const contained of [false, true]) {
          for (const context of [0, 10]) {
            const opts = { contained, context }
            const where = `${name} ${refName}:${region.start}-${region.end} ${contained ? 'contained' : 'plain'} context ${context}`
            const cut = await settle(
              reader.getCutForRange(refName, region.start, region.end, {
                ...opts,
                stepBudget: Infinity,
              }),
            )
            reads.length = 0
            const runs = await settle(
              reader.getRunsForRange(refName, region.start, region.end, opts),
            )
            const fetched = unitsRead(units, reads)
            const bytes = bytesRead(reads)
            out.cuts++
            if (fetched.some(u => u.kind === 'detail')) {
              out.details.push(where)
            }
            if (!cut.value || !runs.value) {
              if (cut.value || runs.value) {
                out.differences.push(
                  `${where}: ${String(cut.error ?? runs.error)}`,
                )
              }
            } else {
              const allowed = await reader.getBytesForRange(
                refName,
                region.start,
                region.end,
                { context },
              )
              if (bytes > allowed) {
                out.overBytes.push(`${where}: ${bytes} > ${allowed}`)
              }
              const c = cut.value
              const read = runs.value
              const tables = tablesOf(c)
              const expected = fragmentRuns(tables, refName, header.holes)
              const found = firstDifference(read, expected)
              if (found !== undefined) {
                out.differences.push(`${where}: ${found}`)
              }
              for (const k of casesOf(expected, tables, starts, contained)) {
                out.cases.set(k, out.cases.get(k)! + 1)
                fixtureCases.add(k)
              }
              const gapCheck = endsAgainstGaps(expected, tables)
              for (const [what, n] of gapCheck.checked) {
                out.gaps.set(what, (out.gaps.get(what) ?? 0) + n)
              }
              for (const m of gapCheck.mismatches) {
                out.gapMismatches.push(`${where}: ${m}`)
              }
              if (!contained) {
                const gate = departuresAgainstRuns(c, expected, refName)
                out.departuresCompared += gate.compared
                for (const m of gate.mismatches) {
                  out.departureMismatches.push(`${where}: ${m}`)
                }
              }
            }
          }
        }
      }
    }
  }
  return out
}

describe('the fixture sweep: every fixture and reference sequence, plain and contained, with context 0 and 10', () => {
  let result: Sweep
  beforeAll(async () => {
    result = await sweep()
    console.warn(
      `runs sweep: ${result.fixtures} fixtures, ${result.cuts} cuts, ${result.refused.length} fixtures refused (${result.refused.join(', ')})`,
    )
    console.warn(
      `cuts holding each case: ${[...result.cases].map(([c, n]) => `${c} ${n}`).join(', ')}`,
    )
    console.warn(
      `ends held to the gaps: ${[...result.gaps].map(([c, n]) => `${c} ${n}`).join(', ')}`,
    )
    console.warn(
      `haplotypes whose departures walkDepartures and the runs both read: ${result.departuresCompared}`,
    )
  }, 900_000)

  test('runs equal the cut read back as runs in every window', () => {
    expect(result.cuts).toBeGreaterThan(5000)
    expect(result.differences.slice(0, 20)).toEqual([])
  })

  test('runs refuse every fixture whose header does not state revisited:i:0, and only those', async () => {
    const stated = await Promise.all(
      FIXTURES.map(async name => {
        const { revisited } = await new WalkFile({
          walksPath: `${DATA}${name}.walks.gz`,
        }).getHeader()
        return revisited === 0 ? [] : [name]
      }),
    )
    expect(result.refused).toEqual(stated.flat())
    expect(result.refused).toContain('contained.R')
  })

  test('every case is held by some cut', () => {
    expect(
      [...result.cases].filter(([, n]) => n === 0).map(([c]) => c),
    ).toEqual([])
  })

  test.each<[string, Case[]]>([
    ['abut.R', ['record between pieces', 'piece join']],
    ['between.R', ['piece join', 'placeholder', 'record between pieces']],
    ['records.GRCh38', ['record at an end']],
    ['twolines.R', ['record between pieces', 'inverted stretch']],
    ['hole.R', ['hole']],
    ['holetile.R', ['hole']],
    ['refhole.R', ['hole']],
    ['excursion.R', ['contained split', 'jump back']],
    ['drawnedge.R', ['placeholder', 'edge']],
    ['abredge.R', ['placeholder', 'tail']],
    ['blocks.R', ['placeholder', 'piece join', 'inverted stretch']],
    ['abridged.R', ['placeholder', 'record between pieces']],
    ['hprc-chrM.GRCh38', ['single-base variant', 'deletion']],
  ])('%s holds its cases: %s', (name, cases) => {
    expect(
      cases.filter(c => result.byFixture.get(name)?.has(c) !== true),
    ).toEqual([])
  })

  test('each end equals the gap rows beside its fragment', () => {
    expect(result.gaps.size).toBeGreaterThanOrEqual(5)
    expect(result.gapMismatches.slice(0, 20)).toEqual([])
  })

  test('walkDepartures reads the departures the runs read, haplotype by haplotype', () => {
    expect(result.departuresCompared).toBeGreaterThan(100)
    expect(result.departureMismatches.slice(0, 20)).toEqual([])
  })

  test('no runs read a detail unit, or more bytes than bytesFor says of the window', () => {
    expect(result.details).toEqual([])
    expect(result.overBytes).toEqual([])
  })
})

// One tile of R, nodes 1 to 7 of 10, 1, 10, 10, 10, 10 and 10 bp from 0 to
// 61, and nodes 20, 21 and 22 off R of 1, 7 and 12 bp. SNP takes 20 for 2,
// DEL skips 4, INV takes 4 reversed, JMP goes back from 4 to 2, LT has 21
// before 1 and 22 after 7, EXC goes from 2 to 6 and back to 3, and PH is
// abridged: 1 to 3, a run of 15 bp, 5 to 7 and a run of 9 bp
const ONE_TILE = [
  '#walks format:i:8 tile:i:100 revisited:i:0',
  '#reference R',
  ...['SNP', 'DEL', 'INV', 'JMP', 'LT', 'EXC', 'PH'].map(
    h => `#haplotype ${h}#1`,
  ),
  'R#0#c 0 ref R#0#c 0 0 0 7 2,2,2,2,2,2,2 10,1,10,10,10,10,10 .',
  'R#0#c 0 nodes 20,1,1 1,7,12',
  'R#0#c 0 route 2 38 -34',
  'R#0#c 0 route 0 . 4',
  'R#0#c 0 route -4 . 3',
  'R#0#c 0 route 1 . 2',
  'R#0#c 0 route -2 . -4',
  'R#0#c 0 route . 38 -40',
  'R#0#c 0 route 12 30 .',
  'R#0#c 0 route -40 . 8',
  'R#0#c 0 route 0 . -6',
  'R#0#c 0 walk SNP#1#s 0 0 0 7 2 0',
  'R#0#c 0 walk DEL#1#d 0 0 0 6 2 2',
  'R#0#c 0 walk INV#1#i 0 0 0 7 2 4,2',
  'R#0#c 0 walk JMP#1#j 0 0 0 10 2 8',
  'R#0#c 0 walk LT#1#l 0 0 0 9 2 10,2',
  'R#0#c 0 walk EXC#1#e 0 0 0 5 2 14,2',
  'R#0#c 0 abridged PH#1#p 0 0 0 8 75 0',
  'R#0#c 0 0 61 blocks PH#1#p 0 0 0+21,15,10+30,9',
]

// R as two lines, 1 and 2 of 10 bp from 0 and 3 and 4 from 40, with a hole
// between; J walks 1 to 4 in one walk row
const HOLE = [
  '#walks format:i:8 tile:i:100 revisited:i:0',
  '#reference R',
  '#haplotype J#1',
  '#hole R#0#c 20 40',
  'R#0#c 0 ref R#0#c 0 0 0 2 2,2 10,10 .',
  'R#0#c 0 ref R#0#c 40 40 0 2 6,2 10,10 .',
  'R#0#c 0 walk J#1#j 0 0 0 4 2 .',
]

// R in two tiles of 30 bp, nodes 1 to 6 of 10 bp. JOIN walks 1 to 6 in a
// piece a tile; RECD walks 1 to 3 and node 30, 5 bp off R, in its first
// piece and 5 and 6 in its second, joined by a departure record
const TWO_TILES = [
  '#walks format:i:8 tile:i:30 revisited:i:0',
  '#reference R',
  '#haplotype JOIN#1',
  '#haplotype RECD#1',
  'R#0#c 0 ref R#0#c 0 0 0 3 2,2,2 10,10,10 .',
  'R#0#c 0 nodes 30 5',
  'R#0#c 0 route 6 54 .',
  'R#0#c 0 walk JOIN#1#a 0 0 0 3 2 .',
  'R#0#c 0 walk RECD#1#b 0 0 0 4 2 0',
  'R#0#c 30 ref R#0#c 0 30 1 3 8,2,2 10,10,10 .',
  'R#0#c 30 walk JOIN#1#a 0 30 1 3 8 .',
  'R#0#c 30 walk RECD#1#b 0 35 1 2 10 .',
  'R#0#c 20 50 departure RECD#1#b 30 5 0 30 + 1 40 + 1 0',
]

// A fragment's stretches as `[start, bp, way, offBefore, loaded]`, its lead
// and tail, and its ends as their kinds, breakpoints and own bp
function fragmentsOf(runs: WalkRuns) {
  const { walks, stretches } = runs
  return walks.names.map((name, p) => {
    const held: (string | number)[][] = []
    for (let k = walks.offsets[p]!; k < walks.offsets[p + 1]!; k++) {
      held.push([
        stretches.starts[k]!,
        stretches.bp[k]!,
        stretches.against[k] === 1 ? '-' : '+',
        stretches.offBefore[k]!,
        stretches.loaded[k]!,
      ])
    }
    const end = (ends: WalkRuns['walks']['before']) =>
      `${END_KINDS[ends.kinds[p]!]} ${ends.pos[p]} ${ends.ownBp[p]}`
    return {
      name,
      span: [walks.starts[p]!, walks.ends[p]!],
      stretches: held,
      lead: [walks.leadBp[p]!, walks.leadLoaded[p]!],
      tail: [walks.tailBp[p]!, walks.tailLoaded[p]!],
      before: end(walks.before),
      after: end(walks.after),
    }
  })
}

// The runs of a hand-built file over a window, held to its cut read back
async function handRuns(
  lines: string[],
  region: { start: number; end: number },
  opts: { contained?: boolean; context?: number } = {},
) {
  const reader = walkFileFrom(walkFileOf(lines))
  const cut = await reader.getCutForRange(
    'R#0#c',
    region.start,
    region.end,
    opts,
  )
  const runs = await reader.getRunsForRange(
    'R#0#c',
    region.start,
    region.end,
    opts,
  )
  const { holes } = await reader.getHeader()
  expect(runs).toEqual(fragmentRuns(tablesOf(cut), 'R#0#c', holes))
  return runs
}

const named = (runs: WalkRuns, name: string) =>
  fragmentsOf(runs).filter(f => f.name.startsWith(`${name}#`))

describe('hand-built fragments', () => {
  const whole = { start: 0, end: 61 }

  test('a single-base variant is a run of 1 bp between two stretches, and the reference one stretch', async () => {
    const runs = await handRuns(ONE_TILE, whole, { context: 0 })
    expect(named(runs, 'R')).toEqual([
      {
        name: 'R#0#c',
        span: [0, 61],
        stretches: [[0, 61, '+', 0, 1]],
        lead: [0, 1],
        tail: [0, 1],
        before: 'edge NaN NaN',
        after: 'edge NaN NaN',
      },
    ])
    expect(named(runs, 'SNP')[0]!.stretches).toEqual([
      [0, 10, '+', 0, 1],
      [11, 50, '+', 1, 1],
    ])
    expect(runs.reference.names).toEqual(['1', '2', '3', '4', '5', '6', '7'])
    expect(runs.drawn).toEqual({ refName: 'R#0#c', start: 0, end: 61 })
  })

  test('a deletion is two stretches with no run between, and an inverted stretch one taken against R', async () => {
    const runs = await handRuns(ONE_TILE, whole, { context: 0 })
    expect(named(runs, 'DEL')[0]!.stretches).toEqual([
      [0, 21, '+', 0, 1],
      [31, 30, '+', 0, 1],
    ])
    expect(named(runs, 'INV')[0]!.stretches).toEqual([
      [0, 21, '+', 0, 1],
      [21, 10, '-', 0, 1],
      [31, 30, '+', 0, 1],
    ])
  })

  test('a jump back is a stretch that begins behind where the one before ends', async () => {
    const runs = await handRuns(ONE_TILE, whole, { context: 0 })
    expect(named(runs, 'JMP')[0]!.stretches).toEqual([
      [0, 31, '+', 0, 1],
      [10, 51, '+', 0, 1],
    ])
  })

  test('a contained cut keeps a lead and a tail of whole loaded nodes past its context, and of a placeholder only the context', async () => {
    const contained = await handRuns(ONE_TILE, whole, {
      contained: true,
      context: 10,
    })
    expect(named(contained, 'LT')).toEqual([
      {
        name: 'LT#1#l',
        span: [0, 80],
        stretches: [[0, 61, '+', 0, 1]],
        lead: [7, 1],
        tail: [12, 1],
        before: 'edge NaN NaN',
        after: 'edge NaN NaN',
      },
    ])
    expect(named(contained, 'PH')[0]!.tail).toEqual([9, 0])
    const short = await handRuns(ONE_TILE, whole, {
      contained: true,
      context: 5,
    })
    expect(named(short, 'PH')[0]!.tail).toEqual([5, 0])
    const plain = await handRuns(ONE_TILE, whole, { context: 10 })
    expect(named(plain, 'LT')[0]!.span).toEqual([7, 68])
  })

  test("a placeholder is a run of an abridged piece's blocks row, not loaded", async () => {
    const runs = await handRuns(ONE_TILE, whole, { context: 0 })
    expect(named(runs, 'PH')[0]!.stretches).toEqual([
      [0, 21, '+', 0, 1],
      [31, 30, '+', 15, 0],
    ])
  })

  test('a contained walk is split at its step onto the reference off the window, and its ends say where it went', async () => {
    const region = { start: 0, end: 21 }
    const contained = await handRuns(ONE_TILE, region, {
      contained: true,
      context: 0,
    })
    expect(named(contained, 'EXC')).toEqual([
      {
        name: 'EXC#1#e',
        span: [0, 11],
        stretches: [[0, 11, '+', 0, 1]],
        lead: [0, 1],
        tail: [0, 1],
        before: 'edge NaN NaN',
        after: 'beyond 41 0',
      },
      {
        name: 'EXC#1#e',
        span: [21, 31],
        stretches: [[11, 10, '+', 0, 1]],
        lead: [0, 1],
        tail: [0, 1],
        before: 'elsewhere 51 0',
        after: 'edge NaN NaN',
      },
    ])
    const plain = await handRuns(ONE_TILE, region, { context: 0 })
    expect(named(plain, 'EXC')[0]!.stretches).toEqual([
      [0, 11, '+', 0, 1],
      [41, 10, '+', 0, 1],
      [11, 10, '+', 0, 1],
    ])
  })

  test('a stretch never spans a hole, though the walk follows R across it', async () => {
    const runs = await handRuns(HOLE, { start: 0, end: 60 }, { context: 0 })
    expect(named(runs, 'J')[0]!.stretches).toEqual([
      [0, 20, '+', 0, 1],
      [40, 20, '+', 0, 1],
    ])
    expect(runs.holes).toEqual([{ start: 20, end: 40 }])
    expect(runDepartures(runs, runs.walks.names.indexOf('J#1#j'))).toEqual([])
  })

  test('pieces join into one stretch, and a departure record between them is a run', async () => {
    const runs = await handRuns(
      TWO_TILES,
      { start: 0, end: 60 },
      {
        context: 0,
      },
    )
    expect(named(runs, 'JOIN')).toEqual([
      expect.objectContaining({
        span: [0, 60],
        stretches: [[0, 60, '+', 0, 1]],
      }),
    ])
    expect(named(runs, 'RECD')[0]!.stretches).toEqual([
      [0, 30, '+', 0, 1],
      [40, 20, '+', 5, 1],
    ])
  })

  test('a departure record between pieces is a fragment’s end where its rejoin lies past the window', async () => {
    const runs = await handRuns(
      TWO_TILES,
      { start: 0, end: 30 },
      {
        context: 0,
      },
    )
    expect(named(runs, 'RECD')[0]!.after).toBe('beyond 40 5')
    expect(named(runs, 'JOIN')[0]!.after).toBe('edge NaN NaN')
  })

  test('a walk that departs across the whole window holds no fragment, and its departures are listed as departing, computed in its piece or read from a record', async () => {
    const jumped = await handRuns(
      ONE_TILE,
      { start: 22, end: 30 },
      {
        context: 0,
      },
    )
    expect(named(jumped, 'DEL')).toEqual([])
    const { departing } = jumped
    expect(departing.names).toEqual(['DEL#1#d', 'PH#1#p'])
    expect([...departing.hapStarts]).toEqual([21, 21])
    expect([...departing.leavePos]).toEqual([21, 21])
    expect([...departing.rejoinPos]).toEqual([31, 31])
    expect([...departing.ownBp]).toEqual([0, 15])
    const recorded = await handRuns(
      TWO_TILES,
      { start: 31, end: 39 },
      {
        context: 0,
      },
    )
    expect(named(recorded, 'RECD')).toEqual([])
    expect(recorded.departing.names).toEqual(['RECD#1#b'])
    expect([
      recorded.departing.leavePos[0],
      recorded.departing.rejoinPos[0],
      recorded.departing.ownBp[0],
    ]).toEqual([30, 40, 5])
  })

  test('runDepartures reads where each fragment leaves a stretch along its way and rejoins the next', async () => {
    const runs = await handRuns(ONE_TILE, whole, { context: 0 })
    const departed = (name: string) =>
      runDepartures(
        runs,
        runs.walks.names.findIndex(n => n.startsWith(`${name}#`)),
      )
    const one = (
      leavePos: number,
      rejoinPos: number,
      ownBp: number,
      invertedBp: number,
      loaded = true,
    ) => [{ leavePos, rejoinPos, ownBp, invertedBp, loaded }]
    expect(departed('SNP')).toEqual(one(10, 11, 1, 0))
    expect(departed('DEL')).toEqual(one(21, 31, 0, 0))
    expect(departed('INV')).toEqual(one(21, 31, 0, 10))
    expect(departed('JMP')).toEqual(one(31, 10, 0, 0))
    expect(departed('PH')).toEqual(one(21, 31, 15, 0, false))
    expect(departed('R')).toEqual([])
  })
})

describe('getRunsForRange takes a range as the other range methods do', () => {
  const chrM = () =>
    new WalkFile({ walksPath: `${DATA}hprc-chrM.GRCh38.walks.gz` })

  test('a contig resolves on the reference sample to the file’s name for it, and the runs equal the cut’s read back', async () => {
    const file = chrM()
    const byContig = await file.getRunsForRange('chrM', 1000, 1500)
    expect(byContig.refName).toBe('GRCh38#0#chrM')
    expect(byContig).toEqual(
      await file.getRunsForRange('GRCh38#0#chrM', 1000, 1500),
    )
    const cut = await file.getCutForRange('chrM', 1000, 1500)
    expect(byContig).toEqual(fragmentRuns(cut.tables(), 'GRCh38#0#chrM'))
  })

  test('haplotypes keeps their walks beside the reference’s, as getGraphForRange does', async () => {
    const file = chrM()
    const runs = await file.getRunsForRange('chrM', 1000, 1500, {
      haplotypes: ['CHM13#0'],
    })
    const { tables } = await file.getGraphForRange('chrM', 1000, 1500, {
      haplotypes: ['CHM13#0'],
    })
    expect(runs).toEqual(fragmentRuns(tables, 'GRCh38#0#chrM'))
    expect(new Set(runs.walks.names.map(panSNHaplotypeOrSample))).toEqual(
      new Set(['GRCh38#0', 'CHM13#0']),
    )
  })

  test('a sequence the file lacks throws MissingSequenceError', async () => {
    const error: unknown = await chrM()
      .getRunsForRange('chr1', 0, 10)
      .catch((e: unknown) => e)
    expect(error).toBeInstanceOf(MissingSequenceError)
  })
})
