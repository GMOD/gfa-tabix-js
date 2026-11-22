import * as walks from './index.ts'
import {
  pairAlignments,
  pairCigar,
  reverseComplement,
} from './pairAlignment.ts'

import type { PairChain } from './pairAlignment.ts'

const SEQUENCES: Record<number, string> = {
  1: 'ACGTACGTAC',
  2: 'ACGATGCA',
  3: 'CCCCCCCCCC',
  4: 'ACGTTGCA',
  5: 'T'.repeat(24),
  6: 'CA'.repeat(15),
  7: 'GATTACAGAT'.repeat(4),
  8: 'TTGACCAGTA'.repeat(4),
  9: 'A',
  10: 'G',
  11: 'ACGATGCA',
}
const sequenceOf = (id: number) => SEQUENCES[id]!
const forward = (...ids: number[]) => ids.map(id => 2 * id)
const flipped = (...ids: number[]) => ids.map(id => 2 * id + 1)
const count = (edits: [string, number][], op: string) =>
  edits.reduce((sum, [o, len]) => sum + (o === op ? len : 0), 0)

describe('pairAlignments', () => {
  test('aligns the private stretch between two shared nodes base by base', () => {
    const [chain] = pairAlignments({
      query: forward(1, 4, 3),
      target: forward(1, 2, 3),
      sequenceOf,
      minMatch: 1,
    })
    expect(chain).toMatchObject({
      queryStart: 0,
      queryEnd: 28,
      targetStart: 0,
      targetEnd: 28,
      strand: '+',
    })
    expect(pairCigar(chain!.edits)).toBe('13=1X14=')
  })

  test('writes a flipped walk as a - record whose edits read along the target', () => {
    const [chain] = pairAlignments({
      query: flipped(3, 4, 1),
      target: forward(1, 2, 3),
      sequenceOf,
      minMatch: 1,
    })
    expect(chain).toMatchObject({
      queryStart: 0,
      queryEnd: 28,
      targetStart: 0,
      targetEnd: 28,
      strand: '-',
    })
    expect(pairCigar(chain!.edits)).toBe('13=1X14=')
  })

  test('writes two unrelated stretches as an insertion and a deletion', () => {
    const [chain] = pairAlignments({
      query: forward(7, 5, 8),
      target: forward(7, 6, 8),
      sequenceOf,
      minMatch: 1,
    })
    expect(pairCigar(chain!.edits)).toBe('40=24I30D40=')
    expect(chain!.sharedBases).toBe(80)
  })

  test('breaks a record at a private run longer than maxGap', () => {
    const chains = pairAlignments({
      query: forward(7, 5, 8),
      target: forward(7, 6, 8),
      sequenceOf,
      maxGap: 25,
      minMatch: 1,
    })
    expect(chains.map(chain => pairCigar(chain.edits))).toEqual(['40=', '40='])
    expect(
      chains.map(chain => chain.targetStart).sort((a, b) => a - b),
    ).toEqual([0, 70])
  })

  test('starts a run at the first target step when the query visits the last one just before', () => {
    const records = pairAlignments({
      query: forward(8, 5, 7),
      target: forward(7, 8),
      sequenceOf,
      minMatch: 1,
      bases: false,
    })
    expect(
      records.map(r => [r.queryStart, r.targetStart, pairCigar(r.edits)]),
    ).toEqual([
      [0, 40, '40='],
      [64, 0, '40='],
    ])
  })

  test('aligns a single-node target to the same node walked backwards', () => {
    const records = pairAlignments({
      query: flipped(7),
      target: forward(7),
      sequenceOf,
      minMatch: 1,
    })
    expect(records.map(r => [r.strand, pairCigar(r.edits)])).toEqual([
      ['-', '40='],
    ])
  })

  test.each([
    [
      'a SNP bubble',
      forward(1, 9, 3),
      forward(1, 10, 3),
      '10=1X10=',
      '10=1I1D10=',
    ],
    [
      'a SNP bubble on a flipped walk',
      flipped(3, 9, 1),
      forward(1, 10, 3),
      '10=1X10=',
      '10=1D1I10=',
    ],
    [
      'a stretch both walks hold on different nodes',
      forward(1, 4, 3),
      forward(1, 2, 3),
      '13=1X14=',
      '10=8I8D10=',
    ],
    [
      'the same bases on different nodes',
      forward(1, 11, 3),
      forward(1, 2, 3),
      '28=',
      '10=8I8D10=',
    ],
    [
      'two unrelated stretches',
      forward(7, 5, 8),
      forward(7, 6, 8),
      '40=24I30D40=',
      '40=24I30D40=',
    ],
    [
      'an insertion private to the query',
      forward(7, 5, 8),
      forward(7, 8),
      '40=24I40=',
      '40=24I40=',
    ],
    [
      'a deletion private to the target',
      forward(7, 8),
      forward(7, 5, 8),
      '40=24D40=',
      '40=24D40=',
    ],
  ])(
    'with bases: false, writes %s from the shared nodes alone',
    (_, query, target, compared, shared) => {
      const records = (bases: boolean) =>
        pairAlignments({ query, target, sequenceOf, minMatch: 1, bases }).map(
          ({ edits, sharedBases, ...span }) => ({
            span,
            cigar: pairCigar(edits),
            matches: count(edits, '='),
            sharedBases,
          }),
        )
      const withBases = records(true)
      const withoutBases = records(false)
      expect(withBases.map(r => r.cigar)).toEqual([compared])
      expect(withoutBases.map(r => r.cigar)).toEqual([shared])
      expect(withoutBases.map(r => r.span)).toEqual(withBases.map(r => r.span))
      expect(withoutBases.every(r => r.matches === r.sharedBases)).toBe(true)
    },
  )
})

// An LCG's low bits are its worst, and 'ACGT'[state % 4] off one writes poly-A
function randomSequence(length: number, seed: number) {
  let state = seed
  return Array.from({ length }, () => {
    state = (state + 0x6d2b79f5) | 0
    let t = Math.imul(state ^ (state >>> 15), 1 | state)
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t
    return 'ACGT'[((t ^ (t >>> 14)) >>> 0) % 4]!
  }).join('')
}

describe('pairAlignments past the exact alignment', () => {
  // long enough to be worth chaining across a 3 kb indel, as a real anchor is
  const FLANK = 200
  const unit = randomSequence(3000, 3)
  const sequences: Record<number, string> = {
    1: randomSequence(FLANK, 1),
    2: randomSequence(FLANK, 2),
    3: unit,
    4: reverseComplement(unit),
    5: randomSequence(60, 5),
    6: randomSequence(60, 6),
  }
  const lookup = (id: number) => sequences[id]!

  test('reports an inversion the graph holds as two unrelated nodes as a - record', () => {
    const records = pairAlignments({
      query: forward(1, 3, 2),
      target: forward(1, 4, 2),
      sequenceOf: lookup,
    })
    expect(
      records.map(r => [
        r.strand,
        r.queryStart,
        r.queryEnd,
        pairCigar(r.edits),
      ]),
    ).toEqual([
      ['+', 0, 3400, '200=3000I3000D200='],
      ['-', 200, 3200, '3000='],
    ])
    expect(records[1]).toMatchObject({ targetStart: 200, targetEnd: 3200 })
  })

  const throughSharedInversion = [...forward(1), ...flipped(3), ...forward(2)]
  const summary = (records: PairChain[]) =>
    records.map(r => [
      r.strand,
      r.queryStart,
      r.queryEnd,
      r.targetStart,
      r.targetEnd,
      pairCigar(r.edits),
      r.sharedBases,
    ])

  test('reports an inversion a walk takes through a shared node once, though its bases align again in the gap of the forward record over it', () => {
    expect(
      summary(
        pairAlignments({
          query: throughSharedInversion,
          target: forward(1, 3, 2),
          sequenceOf: lookup,
        }),
      ),
    ).toEqual([
      ['+', 0, 3400, 0, 3400, '200=3000I3000D200=', 400],
      ['-', 200, 3200, 200, 3200, '3000=', 3000],
    ])
  })

  test('with bases: false, reports the inversion a walk takes through a shared node and none that only the bases would find', () => {
    const noBases = (query: number[], target: number[]) =>
      summary(
        pairAlignments({ query, target, sequenceOf: lookup, bases: false }),
      )
    expect(noBases(forward(1, 3, 2), forward(1, 4, 2))).toEqual([
      ['+', 0, 3400, 0, 3400, '200=3000I3000D200=', 400],
    ])
    expect(noBases(throughSharedInversion, forward(1, 3, 2))).toEqual([
      ['+', 0, 3400, 0, 3400, '200=3000I3000D200=', 400],
      ['-', 200, 3200, 200, 3200, '3000=', 3000],
    ])
  })

  test('aligns a tandem copy by its bases where a walk revisits the node', () => {
    const [record, ...rest] = pairAlignments({
      query: forward(1, 3, 3, 2),
      target: forward(1, 3, 2),
      sequenceOf: lookup,
    })
    expect(rest).toEqual([])
    expect(count(record!.edits, '=')).toBe(3400)
    expect(count(record!.edits, 'I')).toBe(3000)
    expect(count(record!.edits, 'D') + count(record!.edits, 'X')).toBe(0)
    expect(record!.sharedBases).toBe(2 * FLANK)
  })

  test('keeps a chain that matches minMatch bases though a shorter one outscores it', () => {
    const lengths = [90, 200, 200, 60, 60, 50, 200]
    const nodes = Object.fromEntries(
      lengths.map((length, i) => [i + 1, randomSequence(length, 40 + i)]),
    )
    const [P, C, E, A, B, X, D] = [1, 2, 3, 4, 5, 6, 7]
    const records = pairAlignments({
      query: forward(P, C, E, A, B),
      target: forward(A, X, B, D, P),
      sequenceOf: id => nodes[id]!,
      bases: false,
    })
    expect(
      records.map(r => [r.queryStart, r.targetStart, pairCigar(r.edits)]),
    ).toEqual([[490, 0, '60=50D60=']])
  })

  test('leaves anchors worth fewer bases than the gap between them costs apart', () => {
    const records = pairAlignments({
      query: forward(5, 3, 3, 6),
      target: forward(5, 3, 6),
      sequenceOf: lookup,
      minMatch: 1,
    })
    expect(records.map(r => pairCigar(r.edits))).toEqual(['60=', '60='])
  })
})

describe('pairAlignments inside a large private stretch', () => {
  // anchors long enough to chain across the stretch between them
  const records = (query: string, target: string) => {
    const sequences: Record<number, string> = {
      1: randomSequence(300, 11),
      2: randomSequence(300, 12),
      3: query,
      4: target,
    }
    return pairAlignments({
      query: forward(1, 3, 2),
      target: forward(1, 4, 2),
      sequenceOf: id => sequences[id]!,
    })
  }

  test('aligns the gap between two chained k-mer runs by seeding it again', () => {
    // the target repeats the unit 66 times, past the 64 occurrences a k-mer may
    // seed from, so only the flanks chain; the gap between them holds one copy
    const unit = randomSequence(2100, 13)
    const left = randomSequence(1000, 14)
    const right = randomSequence(1000, 15)
    const [record, ...rest] = records(
      left + unit + right,
      unit.repeat(65) + left + unit + right,
    )
    expect(rest).toEqual([])
    expect(count(record!.edits, 'I')).toBe(0)
    expect(count(record!.edits, '=')).toBe(300 + 1000 + 2100 + 1000 + 300)
  })

  test('chains k-mer runs that match minMatch bases though a shorter one outscores them', () => {
    const shorter = randomSequence(99, 31)
    const a = randomSequence(60, 32)
    const b = randomSequence(60, 33)
    const [record] = records(
      randomSequence(1000, 35) +
        shorter +
        randomSequence(1000, 36) +
        a +
        b +
        randomSequence(1000, 37),
      a +
        randomSequence(300, 34) +
        b +
        randomSequence(1000, 38) +
        shorter +
        randomSequence(1000, 39),
    )
    expect(pairCigar(record!.edits)).toMatch(/^300=2099I60=300D6\d=/)
  })

  test('keeps an inversion nested in the gap of the forward chain over the same stretch', () => {
    const left = randomSequence(1000, 14)
    const inverted = randomSequence(1000, 15)
    const right = randomSequence(1000, 16)
    const rs = records(
      left + inverted + right,
      left + reverseComplement(inverted) + right,
    )
    // each flank may run a base into the inversion by chance, which costs the
    // inversion that base and no more
    expect(rs.map(r => r.strand)).toEqual(['+', '-'])
    const inversion = rs[1]!
    expect(inversion.queryStart).toBeGreaterThanOrEqual(1300)
    expect(inversion.queryEnd).toBeLessThanOrEqual(2300)
    expect(inversion.edits.every(([op]) => op === '=')).toBe(true)
    expect(count(inversion.edits, '=')).toBeGreaterThanOrEqual(990)
  })

  test('lets no later chain re-claim, on a paralog, what an earlier record of its orientation spans', () => {
    // flipped against the target, left-gap-right chains as one record, and the
    // gap's own copy sits past it in the target, where only a nested chain
    // could reach it
    const left = randomSequence(1000, 14)
    const gap = randomSequence(1000, 15)
    const other = randomSequence(1000, 16)
    const right = randomSequence(1000, 17)
    const rs = records(
      left + gap + right,
      reverseComplement(left + other + right + gap),
    )
    expect(rs.filter(r => r.strand === '-')).toHaveLength(1)
  })

  test('aligns a base in one record only, the one scoring higher over the stretch two records share', () => {
    // the target holds the query's copy inverted base for base, and forward
    // with every tenth base changed, which no k-mer seeds and only the exact
    // alignment of the forward record's gap reaches
    const unit = randomSequence(1000, 21)
    const diverged = unit.replace(/./g, (base, i: number) =>
      i % 10 === 5 ? (base === 'A' ? 'C' : 'A') : base,
    )
    const left = randomSequence(1000, 22)
    const right = randomSequence(1000, 23)
    const rs = records(
      left + unit + right,
      left + reverseComplement(unit) + diverged + right,
    )
    const alignedTimes = new Uint8Array(300 + 3000 + 300)
    for (const r of rs) {
      let qi = 0
      for (const [op, len] of r.edits) {
        if (op === '=' || op === 'X') {
          for (let k = 0; k < len; k++) {
            alignedTimes[
              r.strand === '-' ? r.queryEnd - 1 - qi - k : r.queryStart + qi + k
            ]! += 1
          }
        }
        qi += op === 'D' ? 0 : len
      }
    }
    expect(alignedTimes.reduce((most, n) => Math.max(most, n), 0)).toBe(1)
    const inversion = rs.find(r => r.strand === '-')!
    expect(count(inversion.edits, '=')).toBeGreaterThanOrEqual(990)
  })
})

describe('the package exports', () => {
  const lengths: Record<number, number> = { 1: 200, 2: 1, 3: 1, 4: 300, 5: 40 }
  const standIn = (id: number) => 'N'.repeat(lengths[id]!)
  const walk = (ids: number[]) => ids.map(id => 2 * id)

  test('aligns two walks on the nodes they share, from node lengths alone', () => {
    const chains = walks.pairAlignments({
      query: walk([1, 3, 5, 4]),
      target: walk([1, 2, 4]),
      sequenceOf: standIn,
      bases: false,
    })
    expect(chains).toHaveLength(1)
    const [chain] = chains
    expect(chain).toMatchObject({
      queryStart: 0,
      queryEnd: 541,
      targetStart: 0,
      targetEnd: 501,
      strand: '+',
      sharedBases: 500,
    })
    expect(walks.pairCigar(chain!.edits)).toBe('200=41I1D300=')
  })

  test('matches a walk that repeats a node against one that passes it once', () => {
    const [pairs, weight] = walks.weightedLcs(
      walk([1, 4, 5, 4, 2]),
      walk([1, 4, 2]),
      handle => lengths[handle / 2]!,
    )
    expect(pairs.map(([, b]) => b)).toEqual([0, 1, 2])
    expect(weight).toBe(501)
  })
})
