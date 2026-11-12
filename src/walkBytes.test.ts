import { readdirSync } from 'node:fs'

import { LocalFile } from 'generic-filehandle2'

import { bytesRead, readsOf, unitsRead } from './testing.ts'
import { WalkReader } from './walkReader.ts'
import { walksKept } from './walkRows.ts'

import type { WalkCutOptions } from './walkReader.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

const walkFiles = (folder: string) =>
  readdirSync(`${DATA}${folder}`)
    .filter(file => file.endsWith('.walks.gz'))
    .map(file => `${folder}${file.slice(0, -'.walks.gz'.length)}`)

const FIXTURES = [...walkFiles(''), ...walkFiles('extension/')]

const CONTEXTS = [undefined, 0, 7]

async function opened(fixture: string) {
  const path = `${DATA}${fixture}.walks.gz`
  const data = new LocalFile(path)
  const reader = WalkReader.fromFilehandles({
    data,
    directory: new LocalFile(`${path}.idx`),
  })
  const header = await reader.header()
  const sequences = await Promise.all(
    (await reader.referenceSequences()).map(async refName => ({
      refName,
      units: await reader.units(refName),
      length: await reader.sequenceLength(refName),
    })),
  )
  return { reader, header, sequences, reads: readsOf(data) }
}

function windowsOf(length: number) {
  const windows = [{ start: 0, end: Math.max(length, 1) }]
  for (const width of [5, 20]) {
    const step = Math.max(width, Math.ceil(length / 8))
    for (let start = 0; start < length; start += step) {
      windows.push({ start, end: start + width })
    }
  }
  return windows
}

// Each window of each sequence of a fixture under each context, with what
// bytesFor says of it and the units a cut of it then reads
async function eachCut(
  fixture: string,
  options: (haplotypes: string[], refName: string) => WalkCutOptions,
) {
  const { reader, header, sequences, reads } = await opened(fixture)
  const cuts = []
  for (const { refName, units, length } of sequences) {
    for (const window of windowsOf(length)) {
      for (const context of CONTEXTS) {
        const said = await reader.bytesFor(refName, window, { context })
        reads.length = 0
        await reader.cut(refName, window, {
          context,
          ...options(header.haplotypes, refName),
        })
        const read = unitsRead(units, reads)
        const bytesOf = (kind: string) =>
          read
            .filter(u => u.kind === kind)
            .reduce((sum, u) => sum + u.length, 0)
        cuts.push({
          where: `${refName}:${window.start}-${window.end} context ${context}`,
          said,
          bytes: bytesRead(reads),
          tiles: bytesOf('tile'),
          blocks: bytesOf('blocks'),
          details: bytesOf('detail'),
          windowBlocks: units
            .filter(u => u.kind === 'blocks')
            .reduce((sum, u) => sum + u.length, 0),
        })
      }
    }
  }
  return cuts
}

test.each(FIXTURES)(
  'bytesFor is the bytes a cut that keeps every walk and expands none reads, on %s',
  async fixture => {
    const cuts = await eachCut(fixture, () => ({}))
    expect(cuts.length).toBeGreaterThan(0)
    expect(cuts.filter(cut => cut.bytes !== cut.said)).toEqual([])
  },
)

test.each(FIXTURES)(
  'a cut that keeps fewer walks reads no more than bytesFor, on %s',
  async fixture => {
    const cuts = await eachCut(fixture, (haplotypes, refName) => ({
      keep: walksKept(
        haplotypes.filter((_, i) => i % 2 === 0),
        refName,
      ),
    }))
    expect(cuts.filter(cut => cut.bytes > cut.said)).toEqual([])
  },
)

test.each(FIXTURES)(
  'a cut that expands every walk reads its detail units, no blocks unit, and else what bytesFor counts but the blocks units, on %s',
  async fixture => {
    const cuts = await eachCut(fixture, () => ({ expand: () => true }))
    expect(
      cuts.filter(
        cut =>
          cut.blocks > 0 ||
          cut.bytes - cut.details > cut.said ||
          cut.bytes - cut.details < cut.said - cut.windowBlocks,
      ),
    ).toEqual([])
  },
)

// ownends.R.gfa: over 24-30 lie the records of X#1, with 30 bp of its own in
// the tile at 0, and of Z#1, with 10 bp in the tile at 10
test('bytesFor counts the tiles a cut fetches for own steps within attachBp, and no blocks unit of theirs', async () => {
  const { reader, sequences, reads } = await opened('ownends.R')
  const [{ refName, units } = { refName: '', units: [] }] = sequences
  const window = { start: 24, end: 30 }
  const bytes = (kind: string, name: string) =>
    units
      .filter(u => u.kind === kind && u.name === name)
      .reduce((sum, u) => sum + u.length, 0)
  const departures = units
    .filter(u => u.kind === 'departures' && u.start < 30 && u.end > 24)
    .reduce((sum, u) => sum + u.length, 0)
  const windowBytes = bytes('tile', '20') + bytes('blocks', '20') + departures
  for (const { attachBp, tiles } of [
    { attachBp: undefined, tiles: ['0', '10'] },
    { attachBp: 29, tiles: ['10'] },
    { attachBp: 9, tiles: [] },
  ]) {
    const said = await reader.bytesFor(refName, window, {
      context: 0,
      attachBp,
    })
    expect(said).toBe(
      tiles.reduce((sum, name) => sum + bytes('tile', name), windowBytes),
    )
    reads.length = 0
    await reader.cut(refName, window, { context: 0, attachBp })
    expect(bytesRead(reads)).toBe(said)
  }
})

// attachmany.R.gfa: ten records, each with its own steps in one of the tiles
// at 0 to 90, leave at 10 to 100 and rejoin at 300 to 390. All ten lie over
// 195-205 and none leaves or rejoins there; over 365-375 lie four, and the
// one that rejoins at 370 names the tile at 20
test.each([
  { window: { start: 195, end: 205 }, own: [], tiles: ['190', '200'] },
  { window: { start: 365, end: 375 }, own: ['20'], tiles: ['360', '370'] },
  {
    window: { start: 300, end: 400 },
    own: Array.from({ length: 10 }, (_, k) => String(10 * k)),
    tiles: Array.from({ length: 10 }, (_, k) => String(300 + 10 * k)),
  },
])(
  'bytesFor counts a tile past the window’s for a record that leaves or rejoins on the window, and none for one that jumps over it: $window',
  async ({ window, own, tiles }) => {
    const { reader, sequences, reads } = await opened('attachmany.R')
    const [{ refName, units } = { refName: '', units: [] }] = sequences
    const bytes = (keep: (u: (typeof units)[number]) => boolean) =>
      units.filter(keep).reduce((sum, u) => sum + u.length, 0)
    const said = await reader.bytesFor(refName, window, { context: 0 })
    expect(said).toBe(
      bytes(u => u.kind === 'tile' && [...tiles, ...own].includes(u.name)) +
        bytes(
          u =>
            u.kind === 'departures' &&
            u.start < window.end &&
            u.end > window.start,
        ),
    )
    reads.length = 0
    await reader.cut(refName, window, { context: 0 })
    expect(bytesRead(reads)).toBe(said)
    expect(
      unitsRead(units, reads)
        .filter(u => u.kind === 'tile')
        .map(u => u.name)
        .sort((a, b) => Number(a) - Number(b)),
    ).toEqual([...own, ...tiles].sort((a, b) => Number(a) - Number(b)))
  },
)

test('bytesFor reads the departure blocks a cut would, and none past the departure budget', async () => {
  const { reader, sequences, reads } = await opened('ownends.R')
  const [{ refName, units } = { refName: '', units: [] }] = sequences
  const window = { start: 24, end: 30 }
  await reader.bytesFor(refName, window, { context: 0 })
  expect(unitsRead(units, reads).map(u => u.kind)).toEqual(['departures'])
  reads.length = 0
  await expect(
    reader.bytesFor(refName, window, { context: 0, departureByteBudget: 1 }),
  ).rejects.toThrow('Too many departures')
  expect(reads).toEqual([])
})
