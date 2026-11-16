import { readFileSync, readdirSync } from 'node:fs'

import { LocalFile } from 'generic-filehandle2'
import { expect, test } from 'vitest'

import { rangesRead, readsOf } from './testing.ts'
import { RegionTooLargeError } from './util.ts'
import { COUNT_REASONS } from './walkCounts.ts'
import { WalkReader } from './walkReader.ts'

import type { CountColumns, WalkCounts } from './walkCounts.ts'
import type { WalkCountOptions } from './walkReader.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

const FIXTURES = ['', 'extension/'].flatMap(folder =>
  readdirSync(`${DATA}${folder}`)
    .filter(file => file.endsWith('.R.walks.gz'))
    .map(file => `${folder}${file.slice(0, -'.R.walks.gz'.length)}`),
)

const TWINS = new Map([
  ['abridged', 'unabridged'],
  ['advwhole', 'adv'],
  ['countsabridged', 'counts'],
  ['detailtwicewhole', 'detailtwice'],
  ['edgewhole', 'edge'],
  ['gapnanwhole', 'gapnan'],
  ['holetile', 'hole'],
  ['ownendswhole', 'ownends'],
  ['secondlinewhole', 'secondline'],
  ['sortsfirstwhole', 'sortsfirst'],
  ['unbuiltwhole', 'unbuilt'],
])

const gfaOf = (fixture: string) =>
  fixture.startsWith('extension/')
    ? 'extension/follow'
    : (TWINS.get(fixture) ?? fixture)

const SETTINGS: [string, WalkCountOptions][] = [
  ['unexpanded', {}],
  ['expanded', { expand: () => true }],
  ['expanded within 3 steps', { expand: () => true, stepBudget: 3 }],
]

const linkName = (a: number, b: number) => {
  const there = `${a},${b}`
  const back = `${b ^ 1},${a ^ 1}`
  return there < back ? there : back
}

const handleOf = (strand: string | undefined, id: number) =>
  2 * id + (strand === '<' ? 1 : 0)

const overlaps = (
  a: { start: number; end: number },
  b: { start: number; end: number },
) => a.start < b.end && a.end > b.start

// Each node's and each link's haplotypes, with how often each steps it, and
// each node's passes by the reference sequences, read from the GFA's W lines
// less their 0 bp steps, which gfa-to-tabix drops
function truthOf(gfa: string, references: Set<string>) {
  const lines = gfa.split('\n').map(line => line.split('\t'))
  const lengths = new Map<number, number>()
  for (const [kind, id, bases, ...tags] of lines) {
    if (kind === 'S') {
      const stated = tags.find(tag => tag.startsWith('LN:i:'))
      lengths.set(
        Number(id),
        bases === '*' ? Number(stated?.slice(5)) : (bases ?? '').length,
      )
    }
  }
  const nodes = new Map<number, Map<string, number>>()
  const links = new Map<string, Map<string, number>>()
  const passes = new Map<
    number,
    { refName: string; start: number; end: number }[]
  >()
  const step = <K>(
    into: Map<K, Map<string, number>>,
    key: K,
    haplotype: string,
  ) => {
    const held = into.get(key) ?? new Map<string, number>()
    into.set(key, held)
    held.set(haplotype, (held.get(haplotype) ?? 0) + 1)
  }
  for (const [kind, sample, hap, contig, start, , walk] of lines) {
    if (kind === 'W') {
      const haplotype = `${sample}#${hap}`
      const refName = `${haplotype}#${contig}`
      let at = Number(start)
      let before: number | undefined
      const steps = [...(walk ?? '').matchAll(/([<>])(\d+)/g)]
        .map(([, strand, digits]) => ({ strand, id: Number(digits) }))
        .filter(({ id }) => lengths.get(id)! > 0)
      for (const { strand, id } of steps) {
        const end = at + lengths.get(id)!
        const handle = handleOf(strand, id)
        step(nodes, id, haplotype)
        if (before !== undefined) {
          step(links, linkName(before, handle), haplotype)
        }
        if (references.has(refName)) {
          const held = passes.get(id) ?? []
          passes.set(id, held)
          held.push({ refName, start: at, end })
        }
        before = handle
        at = end
      }
    }
  }
  return { nodes, links, passes, lengths }
}

type Truth = ReturnType<typeof truthOf>

// Windows from every other base on, of a few sizes and the whole sequence,
// in units of a quarter of its shortest node on the reference
function windowsOf(refName: string, length: number, truth: Truth) {
  const shortest = [...truth.passes].reduce(
    (least, [id, held]) =>
      held.some(pass => pass.refName === refName)
        ? Math.min(least, truth.lengths.get(id)!)
        : least,
    length,
  )
  const unit = Math.max(1, Math.floor(shortest / 4))
  const windows = [{ start: 0, end: length }]
  for (let start = 0; start < length; start += 2 * unit) {
    for (const size of [1, 2, 5, 9, 17, 40]) {
      const end = start + size * unit
      if (end <= length + 3 * unit) {
        windows.push({ start, end })
      }
    }
  }
  return windows
}

const sumOf = (held: Map<string, number> | undefined) =>
  [...(held?.values() ?? [])].reduce((sum, n) => sum + n, 0)

// What is wrong with one count against its truth
function countProblems(
  columns: CountColumns,
  k: number,
  truth: Map<string, number> | undefined,
) {
  const haplotypes = truth?.size ?? 0
  const passes = sumOf(truth)
  return [
    columns.haplotypes[k]! > haplotypes ? 'haplotypes above the GFA' : '',
    columns.passes[k]! > passes ? 'passes above the GFA' : '',
    columns.haplotypeReasons[k] === 0 && columns.haplotypes[k] !== haplotypes
      ? `haplotypes ${columns.haplotypes[k]} exact, the GFA's ${haplotypes}`
      : '',
    columns.passReasons[k] === 0 && columns.passes[k] !== passes
      ? `passes ${columns.passes[k]} exact, the GFA's ${passes}`
      : '',
  ].filter(Boolean)
}

// What is wrong with counts over `drawn` against the GFA
function problemsOf(
  counts: WalkCounts,
  truth: Truth,
  refName: string,
  drawn: { start: number; end: number },
) {
  const { nodes, links } = counts
  const problems: string[] = []
  const passesOf = (id: number) => truth.passes.get(id) ?? []
  const onceOver = (id: number) => {
    const held = passesOf(id)
    return (
      held.length === 1 &&
      held[0]!.refName === refName &&
      overlaps(held[0]!, drawn)
    )
  }
  nodes.ids.forEach((id, k) => {
    const reasons = nodes.passReasons[k]!
    const held = passesOf(id)
    for (const problem of countProblems(nodes, k, truth.nodes.get(id))) {
      problems.push(`node ${id}: ${problem}`)
    }
    if (reasons === 0 && held.length !== 1) {
      problems.push(
        `node ${id}: passes exact, and the reference passes it ${held.length} times`,
      )
    }
    if (
      held.length > 1 &&
      held.some(pass => pass.refName === refName && overlaps(pass, drawn)) &&
      (reasons & COUNT_REASONS.revisited) === 0
    ) {
      problems.push(
        `node ${id}: passed again on the drawn stretch, not revisited`,
      )
    }
    if (held.length === 0 && (reasons & COUNT_REASONS.elsewhere) === 0) {
      problems.push(`node ${id}: off the reference, not elsewhere`)
    }
    if (k > 0 && nodes.ids[k - 1]! >= id) {
      problems.push(`node ${id}: after ${nodes.ids[k - 1]}`)
    }
  })
  const listed = new Map<string, number>()
  links.from.forEach((from, k) => {
    const to = links.to[k]!
    const name = linkName(from, to)
    listed.set(name, k)
    for (const problem of countProblems(links, k, truth.links.get(name))) {
      problems.push(`link ${from},${to}: ${problem}`)
    }
    const [lastFrom, lastTo] = [links.from[k - 1], links.to[k - 1]]
    if (
      lastFrom !== undefined &&
      lastTo !== undefined &&
      (lastFrom > from || (lastFrom === from && lastTo >= to))
    ) {
      problems.push(`link ${from},${to}: after ${lastFrom},${lastTo}`)
    }
  })
  const exact = (columns: CountColumns, k: number | undefined) =>
    k !== undefined &&
    columns.haplotypeReasons[k] === 0 &&
    columns.passReasons[k] === 0
  for (const id of truth.passes.keys()) {
    const k = nodes.ids.indexOf(id)
    if (onceOver(id) && !exact(nodes, k < 0 ? undefined : k)) {
      problems.push(`node ${id}: passed once on the drawn stretch, not exact`)
    }
  }
  for (const name of truth.links.keys()) {
    const [a, b] = name.split(',').map(Number)
    if (
      onceOver(a! >> 1) &&
      onceOver(b! >> 1) &&
      !exact(links, listed.get(name))
    ) {
      problems.push(
        `link ${name}: between two nodes passed once on the drawn stretch, not exact`,
      )
    }
  }
  return problems
}

// The nodes and links of a cut with no count
async function uncounted(
  reader: WalkReader,
  refName: string,
  region: { start: number; end: number },
  opts: WalkCountOptions,
  counts: WalkCounts,
) {
  const cut = await reader.cut(refName, region, opts)
  const { nodes, links } = cut.graph.tables(
    cut.kept,
    cut.fragments,
    cut.departures,
  )
  const ids = new Set(counts.nodes.ids)
  const counted = new Set(
    [...counts.links.from].map((from, k) =>
      linkName(from, counts.links.to[k]!),
    ),
  )
  const real = (i: number) => nodes.placeholders?.[i] !== 1
  const missing = nodes.names.flatMap((name, i) =>
    real(i) && !ids.has(Number(name)) ? [`cut node ${name}: no count`] : [],
  )
  links.from.forEach((i, k) => {
    const j = links.to[k]!
    const strands = links.strands[k]!
    const a = 2 * Number(nodes.names[i]) + (strands & 1)
    const b = 2 * Number(nodes.names[j]) + ((strands >> 1) & 1)
    if (real(i) && real(j) && !counted.has(linkName(a, b))) {
      missing.push(`cut link ${a},${b}: no count`)
    }
  })
  return missing
}

test.each(FIXTURES)(
  '%s: no count is above its GFA, an exact one equals it, the reasons say why, every node and link of a cut has one, and expanded they read what an expanded cut reads',
  async fixture => {
    const path = `${DATA}${fixture}.R.walks.gz`
    const data = new LocalFile(path)
    const reader = WalkReader.fromFilehandles({
      data,
      directory: new LocalFile(`${path}.idx`),
    })
    const sequences = await reader.referenceSequences()
    const reads = readsOf(data)
    const truth = truthOf(
      readFileSync(`${DATA}${gfaOf(fixture)}.R.gfa`, 'utf8'),
      new Set(sequences),
    )
    const problems: string[] = []
    for (const refName of sequences) {
      const length = await reader.sequenceLength(refName)
      for (const region of windowsOf(refName, length, truth)) {
        for (const context of [0, 3]) {
          for (const contained of [false, true]) {
            for (const [setting, more] of SETTINGS) {
              const opts = { context, contained, ...more }
              const where = `${refName}:${region.start}-${region.end}, context ${context}${contained ? ', contained' : ''}, ${setting}`
              const drawn = {
                start: Math.max(0, region.start - context),
                end: region.end + context,
              }
              try {
                reads.splice(0)
                const counts = await reader.counts(refName, region, opts)
                const counted = reads.splice(0)
                const found = problemsOf(counts, truth, refName, drawn)
                if (more.stepBudget === undefined) {
                  for (const problem of await uncounted(
                    reader,
                    refName,
                    region,
                    opts,
                    counts,
                  )) {
                    found.push(problem)
                  }
                  if (
                    more.expand &&
                    JSON.stringify(rangesRead(counted)) !==
                      JSON.stringify(rangesRead(reads))
                  ) {
                    found.push(
                      'reads other bytes than an expanded cut of every walk',
                    )
                  }
                }
                for (const problem of found) {
                  problems.push(`${where}: ${problem}`)
                }
              } catch (error) {
                if (
                  !(
                    more.stepBudget !== undefined &&
                    error instanceof RegionTooLargeError
                  )
                ) {
                  problems.push(`${where}: throws ${String(error)}`)
                }
              }
            }
          }
        }
      }
    }
    expect({ problems: problems.length, first: problems.slice(0, 8) }).toEqual({
      problems: 0,
      first: [],
    })
  },
  60_000,
)
