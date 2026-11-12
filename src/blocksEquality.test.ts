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
