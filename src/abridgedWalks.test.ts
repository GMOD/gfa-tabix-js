import { readFileSync } from 'node:fs'

import { LocalFile } from 'generic-filehandle2'
import { expect, test } from 'vitest'

import { abridgedText, abridgedWalks } from './abridgedWalks.ts'
import { WalkReader } from './walkReader.ts'

import type { WalkCutOptions } from './walkReader.ts'
import type { Abridged } from './walkRows.ts'

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

// Each walk's unloaded pieces as the cut returns them, each piece once, and
// the placeholders it keeps on the walk's contig
function piecesByWalk(cut: Awaited<ReturnType<WalkReader['cut']>>) {
  const pieces = new Map<string, Set<Abridged>>()
  for (const { name, abridged } of [...cut.fragments, ...cut.departures]) {
    for (const piece of abridged ?? []) {
      pieces.set(name, (pieces.get(name) ?? new Set()).add(piece))
    }
  }
  const bp = new Map<string, number>()
  for (const id of cut.kept) {
    const node = cut.graph.nodes.get(id)
    if (node?.placeholder) {
      bp.set(node.refName, (bp.get(node.refName) ?? 0) + node.end - node.start)
    }
  }
  return Array.from(pieces, ([name, held]) => {
    let steps = 0
    for (const piece of held) {
      steps += piece.steps
    }
    return { name, steps, bp: bp.get(name) ?? 0 }
  }).sort((a, b) => a.name.localeCompare(b.name))
}

test('each abridged walk of abridged.R, its steps on its walk and departures and its placeholders by their bp', async () => {
  const { tables } = await cutOf('abridged.R', 'R#0#chr', {
    start: 0,
    end: 100,
  })
  expect(abridgedWalks(tables)).toEqual([
    { name: 'D#1#d', steps: 9, bp: 0 },
    { name: 'H#1#h', steps: 6, bp: 3 },
    { name: 'I#1#i', steps: 5, bp: 0 },
    { name: 'Q#0#chr', steps: 5, bp: 6 },
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
  expect(v).toEqual({ name: 'V#1#c', steps: 7, bp: 20 })
})

const SWEPT = [
  { fixture: 'blocks.R', refName: 'R#0#c' },
  { fixture: 'between.R', refName: 'R#0#c' },
  { fixture: 'unbuilt.R', refName: 'R#0#c' },
  { fixture: 'ownends.R', refName: 'R#0#c' },
  { fixture: 'hairpin.R', refName: 'R#0#c' },
  { fixture: 'adv.R', refName: 'R#0#c1' },
  { fixture: 'abridged.R', refName: 'R#0#chr' },
]

const OPTIONS: WalkCutOptions[] = [
  { context: 0 },
  { context: 10 },
  { context: 0, contained: true },
  { context: 10, contained: true },
  { context: 10, expand: name => name < 'N' },
]

test.each(SWEPT)(
  'over windows of $fixture, a walk counts each unloaded piece the cut returns of it once, and the bp of its placeholders',
  async ({ fixture, refName }) => {
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
          piecesByWalk(cut),
        )
      }
    }
    expect(listed).toBeGreaterThan(0)
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

test('abridgedText says the steps and the bp off the reference in plain words', () => {
  expect(abridgedText({ name: 'HG00097#1#chr6', steps: 1, bp: 0 })).toBe(
    'HG00097#1#chr6: 1 step not loaded',
  )
  expect(
    abridgedText(
      { name: 'HG00097#1#chr6', steps: 140_212, bp: 999 },
      'HG00097#1',
    ),
  ).toBe('HG00097#1: 140,212 steps not loaded, 999 bp off the reference')
  expect(abridgedText({ name: 'a', steps: 2, bp: 1000 })).toBe(
    'a: 2 steps not loaded, 1.0 kb off the reference',
  )
  expect(abridgedText({ name: 'a', steps: 2, bp: 9_949 })).toBe(
    'a: 2 steps not loaded, 9.9 kb off the reference',
  )
  expect(abridgedText({ name: 'a', steps: 2, bp: 2_100_000 })).toBe(
    'a: 2 steps not loaded, 2100 kb off the reference',
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
