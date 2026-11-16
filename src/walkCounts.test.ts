import { readFileSync } from 'node:fs'

import { LocalFile } from 'generic-filehandle2'
import { expect, test } from 'vitest'

import { rangesRead, readsOf, walkFileOf } from './testing.ts'
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

// twocontigs.R.gfa: R walks nodes 1 to 5, 10 bp each. H#1 walks 1 to 3 on
// contig a and 2 to 4 on contig b; H#2 takes 6 between 2 and 3 on both of
// its contigs; G walks R
test('a haplotype whose two contigs step a node counts once, and each step is a pass', async () => {
  const counts = await countsOf('twocontigs', { start: 0, end: 50 })
  expect(counts.haplotypes.file).toBe(4)
  expect(nodeCount(counts, 2)).toEqual({
    haplotypes: 4,
    passes: 6,
    haplotypeReasons: [],
    passReasons: [],
  })
  expect(linkCount(counts, '>2', '>3')).toEqual([
    { haplotypes: 3, passes: 4, haplotypeReasons: [], passReasons: [] },
  ])
  expect(nodeCount(counts, 6)).toEqual({
    haplotypes: 1,
    passes: 2,
    haplotypeReasons: ['elsewhere'],
    passReasons: ['elsewhere'],
  })
})

// contained.R.gfa, a random graph of the review of counts: R walks chrA
// twice, and H5#1 alone takes >10>13, written <13<10, off the reference
test('a contained count with context takes the piece before the drawn stretch that a contained cut keeps', async () => {
  const path = `${DATA}contained.R.walks.gz`
  const reader = WalkReader.fromFilehandles({
    data: new LocalFile(path),
    directory: new LocalFile(`${path}.idx`),
  })
  const opts = { context: 3, contained: true, expand: () => true }
  const region = { start: 44, end: 61 }
  const cut = await reader.cut('R#0#chrA', region, opts)
  const { nodes, links } = cut.graph.tables(
    cut.kept,
    cut.fragments,
    cut.departures,
  )
  const cutLinks = [...links.from].map((i, k) =>
    [nodes.names[i], nodes.names[links.to[k]!]].join(' and '),
  )
  expect(cutLinks).toContain('10 and 13')
  const counts = await reader.counts('R#0#chrA', region, opts)
  expect(linkCount(counts, '>10', '>13')).toMatchObject([
    { haplotypes: 1, passes: 1 },
  ])
})

test('counts past the step budget say what the window holds and what to do instead, offering no choice of haplotypes', async () => {
  await expect(
    countsOf('counts', { start: 0, end: 80 }, { stepBudget: 1 }),
  ).rejects.toMatchObject({
    message:
      'Too many haplotype steps in this window to count: its walk rows hold 29 steps, a step along the reference counting a quarter, over 8 × the step budget of 1. Count a window of about 20bp or less, or leave this one uncounted',
  })
  await expect(
    countsOf('counts', { start: 0, end: 80 }, { stepBudget: 0.5 }),
  ).rejects.toMatchObject({
    message:
      'Too many haplotype steps in this window to count: its walk rows hold 29 steps, a step along the reference counting a quarter, over 8 × the step budget of 0.5. Its first tile alone is over, so leave this window uncounted or raise the step budget',
  })
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

// R walks nodes 1 to 4, 5 bp each, in tiles of 10 bp; A walks them all and B
// the first two
const underHeader = (walks: string) => [
  walks,
  '#reference R',
  '#haplotype A#1',
  '#haplotype B#1',
  'R#0#chr 0 ref R#0#chr 0 0 0 2 2,2 5,5 .',
  'R#0#chr 0 walk A#1#a 0 0 0 2 2 .',
  'R#0#chr 0 walk B#1#b 0 0 0 2 2 .',
  'R#0#chr 10 ref R#0#chr 0 10 1 2 6,2 5,5 .',
  'R#0#chr 10 walk A#1#a 0 10 1 2 6 .',
]

test.each([
  ['#walks format:i:8 tile:i:10 revisited:i:0', []],
  ['#walks format:i:8 tile:i:10', ['elsewhere']],
  ['#walks format:i:8 tile:i:10 revisited:i:1', ['elsewhere']],
])(
  'under %s a count on a once-passed node of the drawn stretch has reasons %j, and a haplotype count of every haplotype the file holds is exact whatever the header states',
  async (walks, reasons) => {
    const counts = await WalkReader.fromFilehandles(
      walkFileOf(underHeader(walks)),
    ).counts(R, { start: 0, end: 20 }, { context: 0 })
    expect(counts.haplotypes.file).toBe(3)
    expect(nodeCount(counts, 3)).toEqual({
      haplotypes: 2,
      passes: 2,
      haplotypeReasons: reasons,
      passReasons: reasons,
    })
    expect(linkCount(counts, '>3', '>4')).toEqual([
      {
        haplotypes: 2,
        passes: 2,
        haplotypeReasons: reasons,
        passReasons: reasons,
      },
    ])
    expect(nodeCount(counts, 1)).toEqual({
      haplotypes: 3,
      passes: 3,
      haplotypeReasons: [],
      passReasons: reasons,
    })
  },
)

test('a ref row that states 2 passes leaves its node revisited, though the header states revisited:i:0', async () => {
  const lines = underHeader('#walks format:i:8 tile:i:10 revisited:i:0').map(
    line => (line.includes(' 0 ref ') ? line.replace(/\.$/, '2,1') : line),
  )
  const counts = await WalkReader.fromFilehandles(walkFileOf(lines)).counts(
    R,
    { start: 0, end: 20 },
    { context: 0 },
  )
  expect(nodeCount(counts, 1)?.passReasons).toEqual(['revisited', 'elsewhere'])
  expect(nodeCount(counts, 2)?.passReasons).toEqual([])
})

function plantedCounts(
  file: string,
  refName: string,
  region: { start: number; end: number },
) {
  const path = `${DATA}planted/${file}.walks.gz`
  return WalkReader.fromFilehandles({
    data: new LocalFile(path),
    directory: new LocalFile(`${path}.idx`),
  }).counts(refName, region, { context: 0 })
}

// The haplotypes of a planted GFA with a W line stepping node `id`
const steppingIn = (gfa: string, id: number) =>
  readFileSync(`${DATA}planted/${gfa}.gfa`, 'utf8')
    .split('\n')
    .filter(line => line.startsWith('W') && line.includes(`>${id}`)).length

// planted/shared.gfa: R walks 1 and 2 on c1 and 2 and 3 on c2, 5 bp each, and
// A walks 1, 2 and 3. gfa-to-tabix's broken.sh files A's step on 3 as the
// interior of a route in c1's tile and leaves c2's tile no row of A. The file
// states revisited:i:1, for node 2, and walk format 8 promises of such a file
// nothing about where a step on a once-passed node is filed
test('a file with a revisited node files a step on a once-passed node away from its pass, and the count there is a lower bound', async () => {
  const file = 'other-sequence-as-interior-revisited'
  expect(steppingIn('shared', 3)).toBe(2)
  const atPass = await plantedCounts(file, 'R#0#c2', { start: 0, end: 10 })
  expect(atPass.haplotypes.file).toBe(2)
  expect(nodeCount(atPass, 3)).toEqual({
    haplotypes: 1,
    passes: 1,
    haplotypeReasons: ['elsewhere'],
    passReasons: ['elsewhere'],
  })
  expect(nodeCount(atPass, 2)).toEqual({
    haplotypes: 1,
    passes: 1,
    haplotypeReasons: ['revisited', 'elsewhere'],
    passReasons: ['revisited', 'elsewhere'],
  })
  const away = await plantedCounts(file, 'R#0#c1', { start: 0, end: 10 })
  expect(nodeCount(away, 3)).toEqual({
    haplotypes: 1,
    passes: 1,
    haplotypeReasons: ['elsewhere'],
    passReasons: ['elsewhere'],
  })
  expect(nodeCount(away, 1)).toEqual({
    haplotypes: 2,
    passes: 2,
    haplotypeReasons: [],
    passReasons: ['elsewhere'],
  })
})

// planted/far.gfa: R walks 1 and 2 on c1 and 3 to 6 on c2, 5 bp each, in
// tiles of 10 bp; A and B walk all six. broken.sh files B's steps on 1 and 2
// as the interior of a route in c2's tile 10, and in the other file A's steps
// on 5 and 6 as one in c2's tile 0, each with no row at the nodes' own tile.
// Both state revisited:i:0, so neither is walk format 8: gfa-to-tabix's
// verify7.py and spec8.py refuse them, and the reader, which reads a window's
// tiles alone, cannot
test.each([
  ['other-sequence-as-interior', 'R#0#c1', { start: 0, end: 10 }, 1],
  ['other-tile-as-interior', 'R#0#c2', { start: 10, end: 20 }, 5],
])(
  '%s states revisited:i:0 and files a step on a reference node away from its pass, which the reader cannot tell: its count there is marked exact, one short',
  async (file, refName, region, id) => {
    expect(steppingIn('far', id)).toBe(3)
    const counts = await plantedCounts(file, refName, region)
    expect(counts.haplotypes.file).toBe(3)
    expect(nodeCount(counts, id)).toEqual({
      haplotypes: 2,
      passes: 2,
      haplotypeReasons: [],
      passReasons: [],
    })
  },
)

// counts.R's tiles are 191, 239, 196 and 177 bytes; a count of 20-30 fetches
// two more to attach own steps, 626 bytes in all
test('counts past the byte budget say what the read would fetch and what to count instead, with no graph to see and no haplotypes to choose', async () => {
  const refusal = (
    region: { start: number; end: number },
    byteBudget: number,
  ) =>
    countsOf('counts', region, { byteBudget }).catch((error: unknown) =>
      error instanceof RegionTooLargeError
        ? {
            message: error.message,
            fitsBp: error.fitsBp,
            zoomHelps: error.zoomHelps,
            bytes: error.bytes,
          }
        : error,
    )
  expect(await refusal({ start: 0, end: 80 }, 500)).toEqual({
    message:
      'Too many bytes to read to count this window: 0.0 MB against a byte budget of 0.0 MB. Count a window of about 40bp or less, or leave this one uncounted',
    fitsBp: 40,
    zoomHelps: true,
    bytes: 803,
  })
  expect(await refusal({ start: 0, end: 80 }, 150)).toEqual({
    message:
      'Too many bytes to read to count this window: 0.0 MB for its first tile alone, against a byte budget of 0.0 MB, so leave this window uncounted',
    fitsBp: 1,
    zoomHelps: false,
    bytes: 803,
  })
  expect(await refusal({ start: 20, end: 30 }, 300)).toEqual({
    message:
      "Too many bytes to read to count this window: 0.0 MB for its first tile alone with the tiles holding its departures' own steps, against a byte budget of 0.0 MB. Count with a smaller attachBp, or leave this window uncounted",
    fitsBp: 1,
    zoomHelps: false,
    bytes: 626,
  })
  expect(await refusal({ start: 20, end: 30 }, 500)).toEqual({
    message:
      "Too many bytes to read to count this window: 0.0 MB with the tiles holding its departures' own steps, against a byte budget of 0.0 MB. Count with a smaller attachBp, or leave this window uncounted",
    fitsBp: 20,
    zoomHelps: true,
    bytes: 626,
  })
})
