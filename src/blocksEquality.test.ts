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
// ones before it
test.each(['blocks.R', 'between.R', 'ownends.R'])(
  'the equality harness finds no placeholder in a loaded cut of %s, each id on one run, and each cut the same from a fresh reader',
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
    expect(run.status).toBe(0)
  },
)
