import { LocalFile } from 'generic-filehandle2'

import {
  cutFacts,
  digestOf,
  idsNotOneToOne,
  placeholderRuns,
} from './testing.ts'
import { WalkFile } from './walkFile.ts'

import type { WalkCutOptions } from './walkFile.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname
const R = 'R#0#c'

function open(fixture: string) {
  const path = `${DATA}${fixture}.walks.gz`
  return new WalkFile({
    walksFilehandle: new LocalFile(path),
    idxFilehandle: new LocalFile(`${path}.idx`),
  })
}

interface Window {
  start: number
  end: number
  opts: WalkCutOptions
}

const cutOf = (reader: WalkFile, { start, end, opts }: Window) =>
  reader.getCutForRange(R, start, end, opts)

// Windows of 5 and 20 bp along a fixture's R: plain, contained, which splits
// a placeholder at its context's end, and expanding some of the walks
async function sweep(reader: WalkFile) {
  const length = await reader.getSequenceLength(R)
  const windows: Window[] = []
  for (let start = 0; start < length; start += 5) {
    for (const width of [5, 20]) {
      for (const opts of [
        { context: 0 },
        { context: 10 },
        { context: 0, contained: true },
        { context: 10, contained: true },
        { context: 10, expand: (name: string) => name < 'N' },
      ]) {
        windows.push({ start, end: start + width, opts })
      }
    }
  }
  return windows
}

// Each cut keeps placeholders: every walk of blocks.R unloaded; PARTIAL's
// placeholder split at the context's end; MULTI alone loaded; and of
// ownends.R, a record's own steps from a tile fetched to attach them
const CUTS = [
  { fixture: 'blocks.R', start: 60, end: 120, opts: { context: 0 } },
  {
    fixture: 'blocks.R',
    start: 100,
    end: 110,
    opts: { context: 10, contained: true },
  },
  {
    fixture: 'blocks.R',
    start: 60,
    end: 80,
    opts: { context: 10, expand: (name: string) => name === 'MULTI#1#c' },
  },
  { fixture: 'ownends.R', start: 24, end: 30, opts: { context: 0 } },
  { fixture: 'interleaved.R', start: 0, end: 36, opts: { context: 0 } },
  { fixture: 'hairpin.R', start: 0, end: 36, opts: { context: 0 } },
]

describe.each(CUTS)(
  'the cut of $fixture over $start-$end, $opts',
  ({ fixture, ...window }) => {
    const facts = async (reader: WalkFile) =>
      cutFacts(await cutOf(reader, window))

    test('holds placeholders, and is the same on a second fresh reader', async () => {
      const fresh = await facts(open(fixture))
      expect(fresh.placeholders).toBeGreaterThan(0)
      expect(await facts(open(fixture))).toEqual(fresh)
    })

    test('is the same, ids and tables byte for byte, on a reader that has made other cuts', async () => {
      const fresh = await facts(open(fixture))
      const used = open(fixture)
      for (const other of (await sweep(used)).reverse()) {
        await cutOf(used, other)
      }
      const after = await facts(used)
      expect(after).toEqual(fresh)
      expect(digestOf(after)).toBe(digestOf(fresh))
      const again = await facts(used)
      expect(again).toEqual(fresh)
      expect(digestOf(again)).toBe(digestOf(fresh))
    })

    test('is the same, ids and tables byte for byte, made beside other cuts of its reader', async () => {
      const fresh = await facts(open(fixture))
      const shared = open(fixture)
      const others = (await sweep(shared))
        .reverse()
        .filter((_, k) => k % 7 === 0)
        .map(other => cutOf(shared, other))
      const [first, last] = await Promise.all([
        cutOf(shared, window),
        cutOf(shared, window),
        Promise.all(others),
      ])
      for (const beside of [cutFacts(first), cutFacts(last)]) {
        expect(beside).toEqual(fresh)
        expect(digestOf(beside)).toBe(digestOf(fresh))
      }
    })
  },
)

test.each([
  'blocks.R',
  'contigends.R',
  'hairpin.R',
  'interleaved.R',
  'ownends.R',
])(
  'in every cut of %s each placeholder id stands for one run, and each run has one id',
  async fixture => {
    const reader = open(fixture)
    const shared: string[] = []
    let placeholders = 0
    for (const window of await sweep(reader)) {
      const cut = await cutOf(reader, window)
      const runs = placeholderRuns(cut)
      placeholders += runs.length
      for (const line of idsNotOneToOne(runs)) {
        shared.push(`${window.start}-${window.end}: ${line}`)
      }
    }
    expect(shared).toEqual([])
    expect(placeholders).toBeGreaterThan(0)
  },
)

test('an id that stands for two runs and a run with two ids are each told, and a run met twice with its id is not', () => {
  const A = 'HG002#1#chr1:118-121'
  const B = 'HG002#1#chr1:300-304'
  expect(
    idsNotOneToOne([
      [-1, A],
      [-2, B],
      [-1, A],
    ]),
  ).toEqual([])
  expect(
    idsNotOneToOne([
      [-1, A],
      [-1, B],
    ]),
  ).toEqual([`-1 ${B}, beside -1 ${A}`])
  expect(
    idsNotOneToOne([
      [-1, A],
      [-2, A],
    ]),
  ).toEqual([`-2 ${A}, beside -1 ${A}`])
})
