import { LocalFile } from 'generic-filehandle2'
import { expect, test } from 'vitest'

import { rangesRead, readsOf } from './testing.ts'
import { RegionTooLargeError } from './util.ts'
import { countReasons } from './walkCounts.ts'
import { DepartureReadError, WalkReader } from './walkReader.ts'

import type { WalkCounts } from './walkCounts.ts'
import type { WalkCountOptions } from './walkReader.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

// counts.R.gfa: R walks nodes 1 to 8, 10 bp each, in tiles of 20 bp. A steps
// on 3 twice, around 9; B walks R backwards; C#1 and C#2 take 10 for 4, and
// F takes it backwards; D starts at 5; G jumps from 2 to 7; H steps on 12
// between 2 and 3, across a tile edge, and I between 3 and 4; J walks 3, 11
// and 4 twice; K steps on 13 for 2, its own run in tile 0's piece. Node 4 is
// stepped by R, A, B, H, I, J and K
const R = 'R#0#chr'

function countsOf(
  fixture: string,
  region: { start: number; end: number },
  opts: WalkCountOptions = {},
) {
  const path = `${DATA}${fixture}.R.walks.gz`
  return WalkReader.fromFilehandles({
    data: new LocalFile(path),
    directory: new LocalFile(`${path}.idx`),
  }).counts(R, region, { context: 0, ...opts })
}

const handle = (step: string) =>
  2 * Number(step.slice(1)) + (step.startsWith('<') ? 1 : 0)

function nodeCount({ nodes }: WalkCounts, id: number) {
  const k = nodes.ids.indexOf(id)
  return k < 0
    ? undefined
    : {
        haplotypes: nodes.haplotypes[k],
        passes: nodes.passes[k],
        haplotypeReasons: countReasons(nodes.haplotypeReasons[k]!),
        passReasons: countReasons(nodes.passReasons[k]!),
      }
}

// The counts of the link a walk takes from `a` to `b`, whichever way round
// the counts write it
function linkCount({ links }: WalkCounts, a: string, b: string) {
  const forward = [handle(a), handle(b)]
  const flipped = [handle(b) ^ 1, handle(a) ^ 1]
  const found: number[] = []
  links.from.forEach((from, k) => {
    const to = links.to[k]
    if (
      (from === forward[0] && to === forward[1]) ||
      (from === flipped[0] && to === flipped[1])
    ) {
      found.push(k)
    }
  })
  return found.map(k => ({
    haplotypes: links.haplotypes[k],
    passes: links.passes[k],
    haplotypeReasons: countReasons(links.haplotypeReasons[k]!),
    passReasons: countReasons(links.passReasons[k]!),
  }))
}

test('a haplotype counts once on a node however often it steps it, and every step is a pass', async () => {
  const counts = await countsOf('counts', { start: 0, end: 80 })
  expect(nodeCount(counts, 3)).toEqual({
    haplotypes: 10,
    passes: 12,
    haplotypeReasons: [],
    passReasons: [],
  })
  expect(nodeCount(counts, 11)).toMatchObject({ haplotypes: 1, passes: 2 })
  expect(nodeCount(counts, 9)).toMatchObject({ haplotypes: 1, passes: 1 })
})

test('a step against the reference counts on its node', async () => {
  const counts = await countsOf('counts', { start: 0, end: 80 })
  expect(nodeCount(counts, 5)).toMatchObject({ haplotypes: 11, passes: 11 })
  expect(nodeCount(counts, 10)).toMatchObject({ haplotypes: 3, passes: 3 })
})

test('a link counts each haplotype once whichever way it walks it, and is listed once', async () => {
  const counts = await countsOf('counts', { start: 0, end: 80 })
  expect(linkCount(counts, '>3', '>4')).toEqual([
    { haplotypes: 5, passes: 5, haplotypeReasons: [], passReasons: [] },
  ])
  expect(linkCount(counts, '>2', '>3')).toEqual([
    { haplotypes: 8, passes: 8, haplotypeReasons: [], passReasons: [] },
  ])
  expect(linkCount(counts, '>3', '>10')).toMatchObject([
    { haplotypes: 3, passes: 3 },
  ])
  expect(linkCount(counts, '>3', '>11')).toMatchObject([
    { haplotypes: 1, passes: 2 },
  ])
  expect(linkCount(counts, '>1', '>2')).toMatchObject([
    { haplotypes: 10, passes: 10 },
  ])
})

test('the haplotypes over the window are those with a piece there or a departure over it, out of every haplotype the file holds', async () => {
  const start = await countsOf('counts', { start: 0, end: 20 })
  expect(start.haplotypes).toEqual({
    file: 12,
    window: 11,
    unloaded: [],
    unattached: [],
  })
  const middle = await countsOf('counts', { start: 25, end: 55 })
  expect(middle.haplotypes.window).toBe(12)
  expect(nodeCount(middle, 4)).toEqual({
    haplotypes: 7,
    passes: 8,
    haplotypeReasons: [],
    passReasons: [],
  })
})

test('a count off the reference is never exact, and one the file holds every haplotype of is', async () => {
  const counts = await countsOf('counts', { start: 0, end: 80 })
  expect(nodeCount(counts, 10)).toMatchObject({
    haplotypeReasons: ['elsewhere'],
    passReasons: ['elsewhere'],
  })
})

test('a departure over the window whose own steps the counts do not read leaves counts off the reference lower bounds, and names its haplotype', async () => {
  const over = await countsOf('counts', { start: 25, end: 35 })
  expect(over.haplotypes.unattached).toEqual(['H#1', 'K#1'])
  expect(nodeCount(over, 12)).toEqual({
    haplotypes: 1,
    passes: 1,
    haplotypeReasons: ['unattached', 'elsewhere'],
    passReasons: ['unattached', 'elsewhere'],
  })
  const rejoining = await countsOf('counts', { start: 20, end: 30 })
  expect(rejoining.haplotypes.unattached).toEqual([])
  expect(nodeCount(rejoining, 12)).toMatchObject({ haplotypes: 2, passes: 2 })
  const past = await countsOf('counts', { start: 20, end: 30 }, { attachBp: 5 })
  expect(past.haplotypes.unattached).toEqual(['F#1', 'H#1', 'K#1'])
  expect(nodeCount(past, 12)).toMatchObject({
    haplotypes: 1,
    haplotypeReasons: ['unattached', 'elsewhere'],
  })
})

test('a piece of a window tile off the drawn stretch counts where it holds the own steps of a departure the counts attach', async () => {
  const counts = await countsOf('counts', { start: 15, end: 30 })
  expect(counts.haplotypes.unattached).toEqual([])
  expect(nodeCount(counts, 13)).toMatchObject({ haplotypes: 1, passes: 1 })
})

test('an abridged file counts the reference exactly from its blocks rows, off it lower bounds by unloaded, and expanded as its unabridged twin', async () => {
  const region = { start: 0, end: 80 }
  const twin = await countsOf('counts', region)
  const blocks = await countsOf('countsabridged', region)
  for (const id of [1, 2, 3, 4, 5, 6, 7, 8]) {
    expect(nodeCount(blocks, id)).toEqual(nodeCount(twin, id))
    expect(nodeCount(blocks, id)?.haplotypeReasons).toEqual([])
  }
  expect(linkCount(blocks, '>3', '>4')).toEqual(linkCount(twin, '>3', '>4'))
  expect(blocks.haplotypes.unloaded).toEqual([
    'A#1',
    'C#1',
    'F#1',
    'H#1',
    'I#1',
    'J#1',
    'K#1',
  ])
  expect(nodeCount(twin, 10)?.haplotypes).toBe(3)
  expect(nodeCount(blocks, 10)).toEqual({
    haplotypes: 1,
    passes: 1,
    haplotypeReasons: ['unloaded', 'elsewhere'],
    passReasons: ['unloaded', 'elsewhere'],
  })
  expect(nodeCount(blocks, 9)).toBeUndefined()
  const expanded = await countsOf('countsabridged', region, {
    expand: () => true,
  })
  expect(expanded).toEqual(twin)
})

test('counts refuse a window past the byte budgets, as a cut does', async () => {
  await expect(
    countsOf('counts', { start: 0, end: 80 }, { byteBudget: 1 }),
  ).rejects.toThrow(RegionTooLargeError)
  await expect(
    countsOf('counts', { start: 0, end: 80 }, { departureByteBudget: 1 }),
  ).rejects.toThrow(DepartureReadError)
  await expect(
    countsOf('counts', { start: 0, end: 80 }, { stepBudget: 1 }),
  ).rejects.toThrow(RegionTooLargeError)
})

test('a step budget that stops expanding leaves the rest unloaded, and says so', async () => {
  const region = { start: 0, end: 80 }
  const twin = await countsOf('counts', region)
  const stopped = await countsOf('countsabridged', region, {
    expand: () => true,
    stepBudget: 4,
  })
  expect(stopped.haplotypes.unloaded.length).toBeGreaterThan(0)
  const short = [...twin.nodes.ids].filter(
    id =>
      (nodeCount(stopped, id)?.haplotypes ?? 0) <
      (nodeCount(twin, id)?.haplotypes ?? 0),
  )
  expect(short.length).toBeGreaterThan(0)
  for (const id of short) {
    expect(nodeCount(stopped, id)?.haplotypeReasons ?? ['unloaded']).toContain(
      'unloaded',
    )
  }
  const fits = await countsOf('countsabridged', region, {
    expand: () => true,
    stepBudget: 1000,
  })
  expect(fits).toEqual(twin)
})

// between.R.gfa: MID#1 walks 7 and 8, jumps back to 1, takes 27 and jumps on
// to 9, each piece abridged; its two records own no steps, so the pieces
// either side of each only touch them
test('counts that expand every walk read the bytes an expanded cut of every walk reads, and no piece that only touches own steps', async () => {
  const path = `${DATA}between.R.walks.gz`
  const data = new LocalFile(path)
  const reader = WalkReader.fromFilehandles({
    data,
    directory: new LocalFile(`${path}.idx`),
  })
  await reader.units('R#0#c')
  const reads = readsOf(data)
  const region = { start: 20, end: 100 }
  const opts = { context: 0, expand: () => true }
  await reader.counts('R#0#c', region, opts)
  const counted = reads.splice(0)
  await reader.cut('R#0#c', region, opts)
  expect(rangesRead(counted)).toEqual(rangesRead(reads))
})
