// Holds a walk file's runs equal to its cuts read back as runs, over the
// windows it is given: WalkFile.getRunsForRange against fragmentRuns of
// getCutForRange's tables with
// the same options and no budget, every haplotype at once, plain and
// contained, with the browser's context. Each fragment end is held to the
// cut's gap rows beside it too, and no runs may read a detail unit. Prints
// each window's time and bytes for the runs and for the cut, every
// difference, and the ends of the walks `--probe` names at a position.
//
//   node scripts/runs-equality.ts <walks.gz or URL>
//     --at <sequence>:<start>-<end>,... [--context 1000]
//     [--probe <position>:<walk prefix>,<walk prefix>,...]
//
// Exits 1 on any difference.
import { LocalFile, RemoteFile } from 'generic-filehandle2'

import { END_KINDS, WalkFile, fragmentRuns } from '../src/index.ts'
import {
  bytesRead,
  endsAgainstGaps,
  firstDifference,
  readsOf,
  unitsRead,
} from '../src/testing.ts'

import { panSNHaplotypeOrSample } from '../src/util.ts'

import type { WalkRuns } from '../src/index.ts'

const flag = (name: string) => {
  const at = process.argv.indexOf(`--${name}`)
  return at < 0 ? undefined : process.argv[at + 1]
}
const location = process.argv[2] ?? ''
const context = Number(flag('context') ?? 1000)
const windows = (flag('at') ?? '')
  .split(',')
  .filter(Boolean)
  .map(at => {
    const [, sequence = '', start = '', end = ''] =
      /^(.*):(\d+)-(\d+)$/.exec(at.replaceAll('_', '')) ?? []
    return { sequence, start: Number(start), end: Number(end) }
  })
const [probeAt = '', probeNames = ''] = (flag('probe') ?? '').split(':')
const probe = {
  pos: Number(probeAt.replaceAll('_', '')),
  names: probeNames.split(',').filter(Boolean),
}

const print = (line: string) => {
  process.stdout.write(`${line}\n`)
}

const open = (url: string) =>
  /^https?:/.test(url) ? new RemoteFile(url) : new LocalFile(url)
const data = open(location)
const reads = readsOf(data)
const reader = new WalkFile({
  walksFilehandle: data,
  idxFilehandle: open(location.replace(/(?=$|[?#])/, '.idx')),
})
const header = await reader.getHeader()
const now = () => performance.now() / 1000

// What the walks `probe` names do at either end of their fragments around
// its position
function probeEnds(runs: WalkRuns) {
  const { walks } = runs
  const lines: string[] = []
  walks.names.forEach((name, p) => {
    const named = probe.names.some(prefix =>
      [name, panSNHaplotypeOrSample(name)].some(n => n.startsWith(prefix)),
    )
    if (named) {
      const end = (side: 'before' | 'after') => {
        const ends = walks[side]
        return `${side} ${END_KINDS[ends.kinds[p]!]} pos ${ends.pos[p]}${ends.reversed[p] === 1 ? '-' : '+'} own ${ends.ownBp[p]} missing ${ends.missingBp[p]}`
      }
      const first = walks.offsets[p]!
      const last = walks.offsets[p + 1]! - 1
      lines.push(
        `PROBE ${name} hap ${walks.starts[p]}-${walks.ends[p]} on R ${runs.stretches.starts[first]}-${runs.stretches.starts[last]! + runs.stretches.bp[last]!}: ${end('before')}; ${end('after')}`,
      )
    }
  })
  return lines
}

let differences = 0
let compared = 0
for (const { sequence, start, end } of windows) {
  const refName =
    (await reader.getSequenceName(header.references[0] ?? '', sequence)) ??
    sequence
  const units = await reader.getUnits(refName)
  for (const contained of [false, true]) {
    const where = `${refName}:${start}-${end} ${contained ? 'contained' : 'plain'} context ${context}`
    const opts = {
      contained,
      context,
      byteBudget: Infinity,
      departureByteBudget: Infinity,
    }
    reads.length = 0
    const t0 = now()
    const runs = await reader.getRunsForRange(sequence, start, end, opts)
    const runsS = now() - t0
    const runsBytes = bytesRead(reads)
    const details = unitsRead(units, reads).filter(u => u.kind === 'detail')
    reads.length = 0
    const t1 = now()
    const cut = await reader.getCutForRange(sequence, start, end, {
      ...opts,
      stepBudget: Infinity,
    })
    const tables = cut.tables()
    const cutS = now() - t1
    const cutBytes = bytesRead(reads)
    const allowed = await reader.getBytesForRange(sequence, start, end, {
      context,
      departureByteBudget: Infinity,
    })
    const expected = fragmentRuns(tables, refName, header.holes)
    const found = firstDifference(runs, expected)
    const gaps = endsAgainstGaps(expected, tables)
    compared++
    const kinds = (side: 'before' | 'after') =>
      END_KINDS.map(
        (kind, k) =>
          `${kind} ${runs.walks[side].kinds.filter(x => x === k).length}`,
      ).join(' ')
    print(
      `${where}: ${runs.walks.names.length} fragments, ${runs.stretches.starts.length} stretches, ${tables.walks.steps.length} steps in the cut; ends before: ${kinds('before')}; after: ${kinds('after')}`,
    )
    print(
      `  departing: ${new Set(runs.departing.names).size} walks with no fragment and ${runs.departing.names.length} departures over the window: ${[...new Set(runs.departing.names)].slice(0, 5).join(', ')}`,
    )
    print(
      `  runs ${runsS.toFixed(2)} s, ${runsBytes} bytes (bytesFor ${allowed}); cut ${cutS.toFixed(2)} s, ${cutBytes} bytes; peak RSS ${Math.round(process.resourceUsage().maxRSS / 1024)} MB`,
    )
    print(
      `  ends held to the gaps: ${[...gaps.checked].map(([what, n]) => `${what} ${n}`).join(', ')}`,
    )
    const problems = [
      ...(found === undefined ? [] : [`runs differ at ${found}`]),
      ...gaps.mismatches.slice(0, 5).map(m => `gap rows differ: ${m}`),
      ...(gaps.mismatches.length > 5
        ? [`and ${gaps.mismatches.length - 5} more gap rows`]
        : []),
      ...details.map(u => `a detail unit read: ${u.name}`),
      ...(runsBytes > allowed ? [`runs read ${runsBytes} > ${allowed}`] : []),
    ]
    for (const problem of problems) {
      print(`  DIFFERENCE ${problem}`)
    }
    differences += problems.length
    if (
      Number.isFinite(probe.pos) &&
      probe.pos >= start &&
      probe.pos < end &&
      !contained
    ) {
      for (const line of probeEnds(runs)) {
        print(`  ${line}`)
      }
    }
  }
}
print(
  `compared: ${compared} pairs of runs and cuts; differences: ${differences}`,
)
process.exitCode = differences > 0 ? 1 : 0
