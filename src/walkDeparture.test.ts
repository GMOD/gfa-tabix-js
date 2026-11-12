import { decodeChunk, referenceNodes } from './walkChunk.ts'
import {
  adjacentSteps,
  departureKind,
  departureRowKey,
  departureStretch,
  parseDepartureRow,
  placementChoice,
  squeezer,
} from './walkDeparture.ts'

const ends = (
  leavePos: number,
  leave: '+' | '-',
  rejoinPos: number,
  rejoin: '+' | '-',
  ownBp = 0,
) => ({
  ownBp,
  leavePos,
  leaveReversed: leave === '-',
  rejoinPos,
  rejoinReversed: rejoin === '-',
})

test('a record row parses to its ends, own sequence and contig position', () => {
  const cols =
    'GRCh38#0#chr1~departure\t19\t81\tdeparture\tHG002#1#ctg\t20\t30\t12\t20\t+\t1\t80\t-\t3'.split(
      '\t',
    )
  expect(parseDepartureRow(cols)).toEqual({
    name: 'HG002#1#ctg',
    hapStart: 20,
    ownBp: 30,
    missingBp: 12,
    leavePos: 20,
    leaveReversed: false,
    leavePasses: 1,
    rejoinPos: 80,
    rejoinReversed: true,
    rejoinPasses: 3,
    fromSteps: false,
  })
  expect(placementChoice(parseDepartureRow(cols))).toBe(true)
  expect(departureRowKey(cols)).toBe(
    'departure\tHG002#1#ctg\t20\t30\t12\t20\t+\t1\t80\t-\t3',
  )
})

test('a record row with no pass counts is refused', () => {
  const cols =
    'GRCh38#0#chr1~departure\t19\t81\tdeparture\tHG002#1#ctg\t20\t30\t12\t20\t+\t80\t-'.split(
      '\t',
    )
  expect(() => parseDepartureRow(cols)).toThrow(/not a departure row/)
})

test('a record row passing an end zero times is refused, as a ref row is', () => {
  const row = (leavePasses: number, rejoinPasses: number) =>
    `GRCh38#0#chr1~departure\t19\t81\tdeparture\tHG002#1#ctg\t20\t30\t12\t20\t+\t${leavePasses}\t80\t-\t${rejoinPasses}`.split(
      '\t',
    )
  expect(() => parseDepartureRow(row(0, 1))).toThrow(/not a departure row/)
  expect(() => parseDepartureRow(row(1, 0))).toThrow(/not a departure row/)
  expect(parseDepartureRow(row(1, 1)).leavePasses).toBe(1)
})

test('a ref row with no pass counts is refused, and an abridged row that is none, that lacks its hapBp or states one below 1', () => {
  const ref = (more: string) =>
    `R#0#c\t0\t1\tref\tR#0#c\t0\t0\t0\t2\t2,2\t6,6${more}`.split('\t')
  expect(referenceNodes([ref('\t.')]).map(n => n.passes)).toEqual([1, 1])
  for (const lacking of ['', '\tSQ:Z:ACGTACGTACGT']) {
    expect(() => referenceNodes([ref(lacking)])).toThrow(
      "a ref row's passes do not count its 2 steps",
    )
  }
  const abridged = (nsteps: string, after = '\t12\t1') =>
    `R#0#c\t0\t6\tabridged\tX#1#c\t0\t0\t0\t${nsteps}${after}`.split('\t')
  for (const row of [
    abridged('x'),
    abridged('5', '\t1'),
    abridged('5', '\t0\t1'),
    abridged('5', '\t12.5\t1'),
  ]) {
    expect(() => decodeChunk([ref('\t.'), row], 0, undefined)).toThrow(
      /not an abridged row/,
    )
  }
  expect(
    decodeChunk([ref('\t.'), abridged('5')], 0, undefined).abridged.map(
      piece => piece.hapBp,
    ),
  ).toEqual([12])
})

test('a departure is a placement choice where either end lies on a node the reference passes more than once', () => {
  expect(placementChoice({ leavePasses: 1, rejoinPasses: 1 })).toBe(false)
  expect(placementChoice({ leavePasses: 2, rejoinPasses: 1 })).toBe(true)
  expect(placementChoice({ leavePasses: 1, rejoinPasses: 2 })).toBe(true)
})

test('a departure is a strand switch, a jump back, own sequence between adjacent steps or a jump forward', () => {
  expect(departureKind(ends(20, '+', 80, '-'))).toBe('switch')
  expect(departureKind(ends(80, '+', 20, '+'))).toBe('back')
  expect(departureKind(ends(20, '-', 80, '-'))).toBe('back')
  expect(departureKind(ends(20, '+', 20, '+', 30))).toBe('adjacent')
  expect(departureKind(ends(20, '+', 80, '+'))).toBe('forward')
  expect(departureKind(ends(80, '-', 20, '-'))).toBe('forward')
})

test('steps either side of a hole in the reference are adjacent', () => {
  const squeeze = squeezer(
    [
      { refName: 'R#0#chr', start: 12, end: 20 },
      { refName: 'R#0#chr2', start: 0, end: 50 },
    ],
    'R#0#chr',
  )
  expect([squeeze(5), squeeze(12), squeeze(20), squeeze(32)]).toEqual([
    5, 12, 12, 24,
  ])
  expect(adjacentSteps(ends(12, '+', 20, '+'), squeeze)).toBe(true)
  expect(adjacentSteps(ends(12, '+', 20, '+'), bp => bp)).toBe(false)
  expect(adjacentSteps(ends(12, '+', 20, '-'), squeeze)).toBe(false)
  expect(adjacentSteps(ends(12, '+', 20, '+', 8), squeeze)).toBe(false)
  expect(departureKind(ends(12, '+', 20, '+', 8), squeeze)).toBe('adjacent')
  expect(squeezer([], 'R#0#chr')(100)).toBe(100)
})

test('a departure covers the bases beside its breakpoints, as its record is filed', () => {
  expect(departureStretch(ends(20, '+', 80, '+'))).toEqual([19, 81])
  expect(departureStretch(ends(70, '-', 40, '-'))).toEqual([39, 71])
  expect(departureStretch(ends(70, '-', 50, '-'))).toEqual([49, 71])
  expect(departureStretch(ends(20, '+', 20, '+', 6))).toEqual([19, 21])
  expect(departureStretch(ends(0, '-', 0, '+'))).toEqual([0, 1])
})
