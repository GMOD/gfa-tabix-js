import { readFileSync } from 'node:fs'

import { LocalFile } from 'generic-filehandle2'
import { expect, test } from 'vitest'

import {
  CUT_PRESETS,
  WALK_BUDGETS,
  WALK_CONTEXT,
  WalkReader,
} from './walkReader.ts'

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

const cellValue = (cell: string) =>
  cell === 'Infinity'
    ? Infinity
    : cell === 'true' || cell === 'false'
      ? cell === 'true'
      : Number(cell.replaceAll(',', ''))

function readmePresets() {
  const text = readFileSync(README, 'utf8')
  const section = text.slice(
    text.indexOf('## Cut presets'),
    text.indexOf('## Fragments'),
  )
  const browser: Record<string, number | boolean> = {}
  const figure: Record<string, number | boolean> = {}
  for (const line of section.split('\n')) {
    const row = /^\| `(\w+)` +\| ([^|]+?) +\| ([^|]+?) +\|$/.exec(line)
    if (row) {
      browser[row[1]!] = cellValue(row[2]!)
      figure[row[1]!] = cellValue(row[3]!)
    }
  }
  return { browser, figure }
}

const BROWSER = {
  context: 1000,
  contained: false,
  attachBp: 65_536,
  stepBudget: 4_000_000,
  byteBudget: 8_000_000,
  departureByteBudget: 2_000_000,
}

const FIGURE = {
  ...BROWSER,
  stepBudget: Infinity,
  byteBudget: Infinity,
  departureByteBudget: Infinity,
}

const CUTS = [
  { fixture: 'ownends.R', refName: 'R#0#c', start: 0, end: 30 },
  { fixture: 'blocks.R', refName: 'R#0#c', start: 60, end: 120 },
  { fixture: 'abridged.R', refName: 'R#0#chr', start: 0, end: 12 },
  { fixture: 'attach.GRCh38', refName: 'GRCh38#0#chr1', start: 0, end: 40 },
  { fixture: 'records.GRCh38', refName: 'GRCh38#0#chr1', start: 0, end: 500 },
]

async function tablesOf(
  { fixture, refName, start, end }: (typeof CUTS)[number],
  opts: WalkCutOptions,
) {
  const { graph, kept, fragments, departures, gaps } = await open(fixture).cut(
    refName,
    { start, end },
    opts,
  )
  return graph.tables(kept, fragments, departures, gaps)
}

test('the presets are the ones the README tabulates', () => {
  expect(readmePresets()).toEqual(CUT_PRESETS)
})

test('the presets keep WALK_CONTEXT and differ only in their budgets', () => {
  expect(CUT_PRESETS).toEqual({ browser: BROWSER, figure: FIGURE })
  expect(CUT_PRESETS.figure.context).toBe(WALK_CONTEXT)
})

test('the presets, each preset and the budgets are frozen, so no caller changes them for another', () => {
  for (const shared of [
    CUT_PRESETS,
    CUT_PRESETS.browser,
    CUT_PRESETS.figure,
    WALK_BUDGETS,
  ]) {
    expect(Object.isFrozen(shared)).toBe(true)
  }
})

test.each(CUTS)(
  'a cut of $fixture with the browser preset is a cut with no options, and the figure preset draws the same graph',
  async cut => {
    const browser = await tablesOf(cut, CUT_PRESETS.browser)
    expect(await tablesOf(cut, {})).toEqual(browser)
    expect(await tablesOf(cut, CUT_PRESETS.figure)).toEqual(browser)
  },
)

test('the figure preset attaches no departure the browser preset leaves without steps', async () => {
  const cut = CUTS[0]!
  const figure = await tablesOf(cut, CUT_PRESETS.figure)
  const lifted = await tablesOf(cut, {
    ...CUT_PRESETS.figure,
    attachBp: Infinity,
  })
  expect(lifted.departures.steps.length).toBeGreaterThan(
    figure.departures.steps.length,
  )
})
