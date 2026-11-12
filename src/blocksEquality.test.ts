import { spawnSync } from 'node:child_process'

const ROOT = new URL('../', import.meta.url).pathname

const harness = (fixture: string, ...args: string[]) =>
  spawnSync(
    process.execPath,
    ['scripts/blocks-equality.ts', `test_data/${fixture}.walks.gz`, ...args],
    { cwd: ROOT, encoding: 'utf8' },
  )

// revref.gfa, reduced from pggb chrY: the reference walks 1-,2-,4-, and H#1
// leaves it after 1 for 3, so H's departure leaves node 1 by its left end
test('the equality harness walks a departure from the end of a node the reference passes reversed', () => {
  const run = harness(
    'revref.GRCh38',
    '--at',
    'GRCh38#0#chr:0-9,GRCh38#0#chr:3-5',
    '--context',
    '2',
  )
  expect(run.stdout).toContain('differences from the loaded cut: 0 in 0')
  expect(run.status).toBe(0)
})

// blocks.R and between.R hold a piece on 0-20 between two of its line over a
// window on 60-90, which a cut that expands its walk loads; every cut of
// either keeps placeholders, and the harness's reader makes each after the
// ones before it. Each has departures that leave or rejoin on a window and
// departures that jump over one
test.each(['blocks.R', 'between.R', 'ownends.R'])(
  'the equality harness finds no placeholder in a loaded cut of %s, each id on one run, each cut the same from a fresh reader, and own steps on the departures that leave or rejoin on the window alone',
  fixture => {
    const run = harness(
      fixture,
      '--windows',
      '30',
      '--width',
      '20',
      '--context',
      '10',
    )
    expect(run.stdout).toContain('differences from the loaded cut: 0 in 0')
    expect(run.stdout).not.toContain(' 0 placeholder steps')
    expect(run.stdout).toMatch(
      /[1-9]\d* departures that leave or rejoin there within the attach limit and [1-9]\d* that do not/,
    )
    expect(run.status).toBe(0)
  },
)

// The graphs of gfa-to-tabix's ADR 0013, where a record's own steps are told
// by hapBp alone: each window is one tile, so a tile is fetched to attach
// wherever a record leaves or rejoins on the window with own steps in
// another, and the harness holds the tiles fetched to those records, each
// placeholder to the run its blocks row gives, each expanding cut to no
// blocks unit, and each cut's bytes to bytesFor
test.each([
  'adv.R',
  'edge.R',
  'ownends.R',
  'pickwrong.R',
  'twinmidlineA.R',
  'twinmidlineB.R',
  'twinotherseqA.R',
  'twinown690.R',
  'twopieces.R',
])(
  'the equality harness finds each tile fetched for %s named by a record that leaves or rejoins on the window, each own run where its blocks row has it, no blocks unit read by an expanding cut, and bytesFor the bytes read',
  fixture => {
    const run = harness(
      fixture,
      '--windows',
      '12',
      '--width',
      '12',
      '--context',
      '0',
    )
    expect(run.stdout).toContain('differences from the loaded cut: 0 in 0')
    expect(run.stdout).toMatch(/[1-9]\d* tiles read past a window's/)
    expect(run.stdout).toMatch(
      /[1-9]\d* of them with own steps in an abridged piece; 0 pairs of cuts threw/,
    )
    expect(run.stdout).toContain('each held to bytesFor exactly')
    expect(run.status).toBe(0)
  },
)
