// Breaks the runs reader on purpose and holds each break to failing: each
// mutant is a patch to a scratch copy of src/, which must fail the runs
// tests (the fixture sweep and the hand-built fragments) or, given a walk
// file, scripts/runs-equality.ts over the windows given. A patch whose
// anchor text is not in its file exactly once fails the run, as does a
// mutant that passes both. The control, every run's `loaded` inverted, must
// fail them too.
//
//   node scripts/runs-mutants.ts --scratch <dir> [--file <walks.gz>]
//     [--at <sequence>:<start>-<end>,...] [--always-script] [--only <name>]
//
// Without --always-script, the script runs only for a mutant the tests miss.
import { spawnSync } from 'node:child_process'
import {
  cpSync,
  existsSync,
  mkdirSync,
  readFileSync,
  rmSync,
  symlinkSync,
  writeFileSync,
} from 'node:fs'
import { join } from 'node:path'

interface Mutant {
  name: string
  file: string
  anchor: string
  patch: string
}

const RUNS = 'src/walkRuns.ts'

const MUTANTS: Mutant[] = [
  {
    name: 'a run dropped from a blocks row',
    file: RUNS,
    anchor: 'if (blocks.off[k]! > 0) {',
    patch: 'if (blocks.off[k]! > 0 && k !== 1) {',
  },
  {
    name: "a block's strand flipped",
    file: RUNS,
    anchor: 'const against = blocks.against[k]!',
    patch: 'const against = k === 1 ? !blocks.against[k]! : blocks.against[k]!',
  },
  {
    name: "a block's offset taken from the previous block's start",
    file: RUNS,
    anchor: 'const blocks = parseBlocks(cols)',
    patch: `const blocks = parseBlocks(cols)
  for (let k = blocks.begins.length - 1; k > 0; k--) {
    blocks.begins[k] = blocks.begins[k]! - blocks.lengths[k - 1]! * (blocks.against[k - 1] ? -1 : 1)
  }`,
  },
  {
    name: 'the bp at a piece join lost',
    file: RUNS,
    anchor: 'haps[i + 1] = haps[i]! + itemBp(item)',
    patch:
      'haps[i + 1] = haps[i]! + (i > 0 && items[i - 1]!.piece !== item.piece ? 0 : itemBp(item))',
  },
  {
    name: 'the bp at a piece join counted twice',
    file: RUNS,
    anchor: 'haps[i + 1] = haps[i]! + itemBp(item)',
    patch:
      'haps[i + 1] = haps[i]! + (i > 0 && items[i - 1]!.piece !== item.piece ? 2 : 1) * itemBp(item)',
  },
  {
    name: 'a record between pieces ignored',
    file: RUNS,
    anchor: `  for (const d of records) {
    leaves.set(`,
    patch: `  for (const d of records.filter(() => false)) {
    leaves.set(`,
  },
  {
    name: 'stretches merged across a hole',
    file: RUNS,
    anchor: 'if (path.starts[i] !== path.end(i - 1)) {',
    patch: 'if (path.starts[i] !== path.end(i - 1) && false) {',
  },
  {
    name: 'a contained walk not split at a reference step off the window',
    file: RUNS,
    anchor: `      if (on.first > on.last) {
        close({ i, k: 0 })
        from = { i: i + 1, k: 0 }
      } else {`,
    patch: `      if (on.first > on.last) {
      } else {`,
  },
  {
    name: 'a loaded tail clipped to exactly context',
    file: RUNS,
    anchor: 'bp += placeholder && size > left ? left : size',
    patch: 'bp += (placeholder || !backward) && size > left ? left : size',
  },
  {
    name: "two walks' names swapped",
    file: RUNS,
    anchor: `  const held = new Set(fragments.map(f => f.name))`,
    patch: `  const named = [...new Set(fragments.map(f => f.name).filter(n => n !== refName))].sort()
  for (const f of fragments) {
    f.name = f.name === named[0] ? named[1]! : f.name === named[1] ? named[0]! : f.name
  }
  const held = new Set(fragments.map(f => f.name))`,
  },
  {
    name: 'a departing walk dropped',
    file: RUNS,
    anchor: '[...jumping, ...recorded].filter(d => !held.has(d.name)),',
    patch:
      '[...jumping, ...recorded].filter(d => !held.has(d.name)).filter((d, _, all) => d.name !== all[0]?.name),',
  },
  {
    name: 'control: every run’s loaded inverted',
    file: RUNS,
    anchor:
      "this.items.push({ kind: 'off', piece: this.piece, bp: [bp], loaded })",
    patch:
      "this.items.push({ kind: 'off', piece: this.piece, bp: [bp], loaded: !loaded })",
  },
]

const flag = (name: string) => {
  const at = process.argv.indexOf(`--${name}`)
  return at < 0 ? undefined : process.argv[at + 1]
}
const root = new URL('../', import.meta.url).pathname
const scratch = flag('scratch') ?? ''
const file = flag('file')
const windows = flag('at') ?? ''
const always = process.argv.includes('--always-script')
const only = flag('only')

const print = (line: string) => {
  process.stdout.write(`${line}\n`)
}

const tally = (output: string) =>
  /Tests\s+(?:\S+\s+)?(\d+) failed/.exec(output)?.[1] ?? '0'

let survivors = 0
for (const mutant of MUTANTS.filter(
  m => only === undefined || m.name === only,
)) {
  const dir = join(scratch, mutant.name.replaceAll(/[^a-z0-9]+/gi, '-'))
  rmSync(dir, { recursive: true, force: true })
  mkdirSync(dir, { recursive: true })
  for (const part of [
    'src',
    'scripts',
    'package.json',
    'tsconfig.json',
    'vitest.config.ts',
  ]) {
    cpSync(join(root, part), join(dir, part), { recursive: true })
  }
  for (const part of ['node_modules', 'test_data']) {
    symlinkSync(join(root, part), join(dir, part))
  }
  const target = join(dir, mutant.file)
  const text = readFileSync(target, 'utf8')
  const count = text.split(mutant.anchor).length - 1
  if (count !== 1) {
    print(`${mutant.name}: ANCHOR found ${count} times in ${mutant.file}`)
    survivors++
  } else {
    writeFileSync(target, text.replace(mutant.anchor, mutant.patch))
    const tests = spawnSync(
      join(dir, 'node_modules/.bin/vitest'),
      ['run', 'src/walkRuns.test.ts'],
      { cwd: dir, encoding: 'utf8', env: { ...process.env, NO_COLOR: '1' } },
    )
    const testsFail = tests.status !== 0
    let line = `${mutant.name}: tests ${testsFail ? `fail (${tally(tests.stdout + tests.stderr)} failed)` : 'pass'}`
    let scriptFails = false
    if (file !== undefined && existsSync(file) && (always || !testsFail)) {
      const script = spawnSync(
        process.execPath,
        [
          '--max-old-space-size=36000',
          'scripts/runs-equality.ts',
          file,
          '--at',
          windows,
        ],
        { cwd: dir, encoding: 'utf8' },
      )
      scriptFails = script.status !== 0
      const summary =
        /compared: .*/.exec(script.stdout)?.[0] ??
        script.stderr.split('\n').find(l => l.includes('Error')) ??
        ''
      line += `; script ${scriptFails ? 'fails' : 'passes'} (${summary})`
    }
    print(line)
    if (!testsFail && !scriptFails) {
      survivors++
    }
  }
  rmSync(dir, { recursive: true, force: true })
}
print(
  `mutants that pass every check, or whose anchor is not found once: ${survivors}`,
)
process.exitCode = survivors > 0 ? 1 : 0
