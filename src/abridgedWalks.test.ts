import { readFileSync } from 'node:fs'

import { LocalFile } from 'generic-filehandle2'
import { expect, test } from 'vitest'

import { abridgedText, abridgedWalks } from './abridgedWalks.ts'
import { WalkReader } from './walkReader.ts'

import type { WalkCutOptions } from './walkReader.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname
const README = new URL('../README.md', import.meta.url).pathname

function open(fixture: string) {
  const path = `${DATA}${fixture}.walks.gz`
  return WalkReader.fromFilehandles({
    data: new LocalFile(path),
    directory: new LocalFile(`${path}.idx`),
  })
}

async function cutOf(
  fixture: string,
  refName: string,
  region: { start: number; end: number },
  opts: WalkCutOptions = {},
) {
  const cut = await open(fixture).cut(refName, region, opts)
  const { graph, kept, fragments, departures } = cut
  return { cut, tables: graph.tables(kept, fragments, departures) }
}

// The bp of each walk's placeholders the cut keeps, by the walk each lies on
function placeholderBp(cut: Awaited<ReturnType<WalkReader['cut']>>) {
  const bp = new Map<string, number>()
  for (const id of cut.kept) {
    const node = cut.graph.nodes.get(id)
    if (node?.placeholder) {
      bp.set(node.refName, (bp.get(node.refName) ?? 0) + node.end - node.start)
    }
  }
  return Array.from(bp, ([name, total]) => ({ name, bp: total })).sort((a, b) =>
    a.name.localeCompare(b.name),
  )
}

test('each walk of abridged.R with a stretch off the reference the cut did not load, by its bp, and none it draws exactly', async () => {
  const { tables } = await cutOf('abridged.R', 'R#0#chr', {
    start: 0,
    end: 100,
  })
  expect(abridgedWalks(tables)).toEqual([
    { name: 'H#1#h', bp: 3 },
    { name: 'Q#0#chr', bp: 6 },
  ])
})

test('a walk counts the placeholders only its departure steps on', async () => {
  const { tables } = await cutOf(
    'ownends.R',
    'R#0#c',
    { start: 0, end: 30 },
    { context: 0 },
  )
  const v = abridgedWalks(tables).find(w => w.name === 'V#1#c')
  const onWalks = new Set<number>()
  tables.walks.names.forEach((name, p) => {
    if (name === 'V#1#c') {
      for (const i of tables.walks.steps.subarray(
        tables.walks.offsets[p],
        tables.walks.offsets[p + 1],
      )) {
        if (tables.nodes.placeholders?.[i] === 1) {
          onWalks.add(i)
        }
      }
    }
  })
  expect(onWalks.size).toBe(0)
  expect(v).toEqual({ name: 'V#1#c', bp: 20 })
})

const SWEPT = [
  { fixture: 'blocks.R', refName: 'R#0#c', lists: true },
  { fixture: 'between.R', refName: 'R#0#c', lists: true },
  { fixture: 'unbuilt.R', refName: 'R#0#c', lists: false },
  { fixture: 'ownends.R', refName: 'R#0#c', lists: true },
  { fixture: 'hairpin.R', refName: 'R#0#c', lists: true },
  { fixture: 'adv.R', refName: 'R#0#c1', lists: true },
  { fixture: 'abridged.R', refName: 'R#0#chr', lists: true },
]

const OPTIONS: WalkCutOptions[] = [
  { context: 0 },
  { context: 10 },
  { context: 0, contained: true },
  { context: 10, contained: true },
  { context: 10, expand: name => name < 'N' },
]

test.each(SWEPT)(
  'over windows of $fixture, a walk counts the bp of each placeholder the cut keeps on its contig once',
  async ({ fixture, refName, lists }) => {
    const reader = open(fixture)
    const length = await reader.sequenceLength(refName)
    let listed = 0
    for (let start = 0; start < length; start += 5) {
      for (const opts of OPTIONS) {
        const cut = await reader.cut(refName, { start, end: start + 20 }, opts)
        const tables = cut.graph.tables(cut.kept, cut.fragments, cut.departures)
        const found = abridgedWalks(tables)
        listed += found.length
        expect(found.toSorted((a, b) => a.name.localeCompare(b.name))).toEqual(
          placeholderBp(cut),
        )
      }
    }
    expect(listed > 0).toBe(lists)
  },
)

test('a cut of the unabridged twin, or one that expands every walk, lists none', async () => {
  const region = { start: 0, end: 100 }
  const twin = await cutOf('unabridged.R', 'R#0#chr', region)
  const expanded = await cutOf('abridged.R', 'R#0#chr', region, {
    expand: () => true,
  })
  expect(abridgedWalks(twin.tables)).toEqual([])
  expect(abridgedWalks(expanded.tables)).toEqual([])
})

test('ownends.R lists the walks with a stretch off the reference not loaded on the window, and not those it draws exactly or whose run off it lies past the window, as its twin draws none of them', async () => {
  const window = { start: 0, end: 30 }
  const { tables } = await cutOf('ownends.R', 'R#0#c', window, { context: 0 })
  const twin = await cutOf('ownendswhole.R', 'R#0#c', window, { context: 0 })
  expect(abridgedWalks(tables)).toEqual([
    { name: 'S#1#c', bp: 70_010 },
    { name: 'V#1#c', bp: 20 },
    { name: 'X#1#c', bp: 30 },
    { name: 'Y#1#c', bp: 40 },
    { name: 'Z#1#c', bp: 10 },
  ])
  expect(abridgedWalks(twin.tables)).toEqual([])
})

test('abridgedText says the bp off the reference not loaded in plain words', () => {
  expect(abridgedText({ name: 'HG00097#1#chr6', bp: 1 })).toBe(
    'HG00097#1#chr6: 1 bp off the reference not loaded',
  )
  expect(abridgedText({ name: 'HG00097#1#chr6', bp: 999 }, 'HG00097#1')).toBe(
    'HG00097#1: 999 bp off the reference not loaded',
  )
  expect(abridgedText({ name: 'a', bp: 1000 })).toBe(
    'a: 1.0 kb off the reference not loaded',
  )
  expect(abridgedText({ name: 'a', bp: 9_949 })).toBe(
    'a: 9.9 kb off the reference not loaded',
  )
  expect(abridgedText({ name: 'a', bp: 2_100_000 })).toBe(
    'a: 2100 kb off the reference not loaded',
  )
})

test('the README example is what ownends.R gives', async () => {
  const text = readFileSync(README, 'utf8')
  const example = text.slice(text.indexOf('abridgedWalks(tables).map('))
  const lines = example
    .slice(0, example.indexOf('```'))
    .split('\n')
    .flatMap(line => /^\/\/ '(.+)'$/.exec(line)?.[1] ?? [])
  const { tables } = await cutOf(
    'ownends.R',
    'R#0#c',
    { start: 0, end: 30 },
    { context: 0 },
  )
  const said = abridgedWalks(tables).map(walk => abridgedText(walk))
  expect(lines).toHaveLength(2)
  expect(said).toEqual(expect.arrayContaining(lines))
})
