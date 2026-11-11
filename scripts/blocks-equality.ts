// Holds a walk file's unloaded cuts equal to its loaded ones. Over random
// windows on its abridged tiles, it compares the cut that expands no walk,
// whose abridged pieces come from their `blocks` rows, with two cuts that
// expand every walk, each run of off-reference steps taken as one step of
// its bp:
//
// - `stated`: the loaded cut as far as a blocks row states a piece
//   (WalkCutOptions.expandStated). Any difference is a fault in the reader
//   or the file
// - `full`: the loaded cut whole. It also holds each loaded piece's own
//   steps after its last step on the reference, which no blocks row states,
//   so a difference here is either that or a fault
//
// Each is cut plain and contained. Compared: the fragments with their
// haplotype offsets, steps on the reference, strands and own-run bp, the bp
// trimmed before and after; the departures; the gaps; each haplotype's cover
// and departures; and the structural forms. Not compared: node ids off the
// reference, a departure's attached steps, and the flag on an own gap
// segment that its end piece is unloaded, which state what is loaded.
//
//   node scripts/blocks-equality.ts <walks.gz or URL> [--windows 300]
//     [--seed 1] [--width 1000] [--context 1000] [--contained-context 0]
//     [--sequence <name>] [--against stated|full] [--from 0]
//
// Prints every difference, then the counts. Exits 1 on a difference from
// the stated cut, 2 on differences from the full cut alone, 0 on none.
// `--against` makes one of the two comparisons, at two thirds of the reads,
// and `--from k` skips the first k windows of the seed's, to go on after a
// stop or to share the windows between runs. One cut is held at a time and
// each window's progress goes to stderr with the resident memory; a cut over
// the densest tiles holds about 1 GB, so give node a heap bound under the
// machine's, as `node --max-old-space-size=3072`.
// A contained cut keeps `--contained-context` bp of a walk past its steps on
// the window: above 0, an unloaded piece's own run there is one step, kept
// whole where the loaded cut keeps the nodes of its first bp, so the two
// differ.
import { LocalFile, RemoteFile } from 'generic-filehandle2'

import {
  WalkReader,
  isPlaceholder,
  structuralForms,
  walkDepartures,
} from '../src/index.ts'

import type {
  WalkDeparture,
  WalkFragment,
  WalkGap,
  WalkNode,
} from '../src/index.ts'

type Cut = Awaited<ReturnType<WalkReader['cut']>>
type Row = Record<string, string | number>
type Rows = Map<string, Row>

interface Window {
  refName: string
  start: number
  end: number
}

const flag = (name: string) => {
  const at = process.argv.indexOf(`--${name}`)
  return at < 0 ? undefined : process.argv[at + 1]
}
const path = process.argv[2] ?? ''
const windowCount = Number(flag('windows') ?? 300)
const seed = Number(flag('seed') ?? 1)
const width = Number(flag('width') ?? 1000)
const context = Number(flag('context') ?? 1000)
const containedContext = Number(flag('contained-context') ?? 0)
const onlySequence = flag('sequence')
const firstWindow = Number(flag('from') ?? 0)
const AGAINST = ['stated', 'full'] as const
const against = AGAINST.filter(
  name => flag('against') === undefined || flag('against') === name,
)

const print = (line: string) => {
  process.stdout.write(`${line}\n`)
}

// mulberry32
function random(from: number) {
  let state = from >>> 0
  return () => {
    state = (state + 0x6d2b79f5) >>> 0
    let t = state
    t = Math.imul(t ^ (t >>> 15), t | 1)
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61)
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}

const lengthIn = (nodes: Map<number, WalkNode>) => (id: number) => {
  const node = nodes.get(id)
  return node ? node.end - node.start : 0
}

// A fragment's steps, each on the reference as its place, strand and node,
// and each run off it as its bp
function fragmentSteps(f: WalkFragment, length: (id: number) => number) {
  const steps: string[] = []
  let own = 0
  let run = 0
  const flush = () => {
    if (run > 0) {
      steps.push(`own:${own}`)
    }
    own = 0
    run = 0
  }
  for (let i = 0; i < f.ids.length; i++) {
    const id = f.ids[i]!
    if (Number.isNaN(f.refStarts[i])) {
      own += length(id)
      run++
    } else {
      flush()
      steps.push(
        `${f.refStarts[i]}${f.against[i] === 1 ? '-' : '+'}${f.rev[i] === 1 ? '<' : '>'}${id}`,
      )
    }
  }
  flush()
  return steps
}

// A placeholder beside another step off the reference, which would leave
// one run of own steps as two
function strayPlaceholders(f: WalkFragment) {
  const offR = (k: number) =>
    k >= 0 && k < f.ids.length && Number.isNaN(f.refStarts[k])
  let stray = 0
  for (let i = 0; i < f.ids.length; i++) {
    stray +=
      isPlaceholder(f.ids[i]!) && (offR(i - 1) || offR(i + 1) || !offR(i))
        ? 1
        : 0
  }
  return stray
}

function fragmentRows(cut: Cut, refName: string) {
  const length = lengthIn(cut.graph.nodes)
  const rows: Rows = new Map()
  for (const f of cut.fragments) {
    if (f.name !== refName) {
      rows.set(`${f.name} line ${f.fragStart} from ${f.hapStart}`, {
        hapEnd: f.hapEnd,
        trimmedBefore: f.trimmedBefore,
        trimmedAfter: f.trimmedAfter,
        steps: fragmentSteps(f, length).join(' '),
      })
    }
  }
  return rows
}

const strand = (reversed: boolean) => (reversed ? '-' : '+')

function departureRows(departures: WalkDeparture[]) {
  const rows: Rows = new Map()
  for (const d of departures) {
    rows.set(
      `${d.name} at ${d.hapStart} leave ${d.leavePos}${strand(d.leaveReversed)} rejoin ${d.rejoinPos}${strand(d.rejoinReversed)}`,
      {
        ownBp: d.ownBp,
        missingBp: d.missingBp,
        leavePasses: d.leavePasses,
        rejoinPasses: d.rejoinPasses,
        from: d.fromSteps ? 'steps' : 'record',
      },
    )
  }
  return rows
}

function gapRows(gaps: WalkGap[]) {
  const rows: Rows = new Map()
  for (const g of gaps) {
    rows.set(`${g.name} edge ${g.edge} from ${g.start}`, {
      end: g.end,
      beside: g.fragment ? g.fragment.hapStart : 'none',
      segments: g.segments
        .map(
          s =>
            `${s.kind}:${s.bp}:missing${s.missingBp}:${s.before}..${s.after}`,
        )
        .join(' '),
    })
  }
  return rows
}

function formRows(cut: Cut, window: Window) {
  const haplotypes = walkDepartures(
    {
      fragments: cut.fragments,
      departures: cut.departures,
      nodes: cut.graph.nodes,
    },
    window.refName,
  )
  const rows: Rows = new Map()
  for (const h of haplotypes) {
    rows.set(h.haplotype, {
      covered: JSON.stringify(h.covered),
      departures: JSON.stringify(h.departures),
    })
  }
  const forms = structuralForms(haplotypes, {
    window: { start: window.start, end: window.end },
  })
  return {
    covers: rows,
    forms: new Map([
      [
        'forms',
        {
          sites: JSON.stringify(forms.sites),
          forms: JSON.stringify(
            forms.forms.map(f => ({ ...f, members: [...f.members].sort() })),
          ),
          unplaced: [...forms.unplaced].sort().join(' '),
        },
      ],
    ]),
  }
}

// Where two space-separated lists first differ, with the items around it
function firstDifference(a: string, b: string) {
  const x = a.split(' ')
  const y = b.split(' ')
  let at = 0
  while (at < x.length && at < y.length && x[at] === y[at]) {
    at++
  }
  const around = (items: string[]) =>
    items.slice(Math.max(0, at - 1), at + 3).join(' ') || 'nothing'
  return `${x.length} against ${y.length} items, at ${at}: ${around(x)} against ${around(y)}`
}

// Each key's row in `unloaded` and `loaded` where they differ, field by
// field, or where one lacks it
function differences(unloaded: Rows, loaded: Rows) {
  const keys = new Set(unloaded.keys())
  for (const key of loaded.keys()) {
    keys.add(key)
  }
  const out: string[] = []
  for (const key of [...keys].sort()) {
    const a = unloaded.get(key)
    const b = loaded.get(key)
    if (!a || !b) {
      out.push(
        `${key}: only ${a ? 'unloaded' : 'loaded'}: ${JSON.stringify(a ?? b).slice(0, 300)}`,
      )
    } else {
      const fields = Object.keys(a)
        .filter(field => a[field] !== b[field])
        .map(field => {
          const u = a[field]!
          const l = b[field]!
          return typeof u === 'string' &&
            typeof l === 'string' &&
            u.includes(' ')
            ? `${field} ${firstDifference(u, l)}`
            : `${field} ${u} against ${l}`
        })
      if (fields.length > 0) {
        out.push(`${key}: ${fields.join('; ')}`)
      }
    }
  }
  return out
}

const named = (fragments: WalkFragment[]) =>
  fragments.map(f => `${f.name} from ${f.hapStart}`)

// What a cut states, as rows to compare, with the fragments holding a piece
// it did not load, those with a placeholder out of place, and its placeholder
// steps; the cut itself is let go
function summary(cut: Cut, window: Window) {
  const { covers, forms } = formRows(cut, window)
  const walks = cut.fragments.filter(f => f.name !== window.refName)
  return {
    rows: {
      fragment: fragmentRows(cut, window.refName),
      departure: departureRows(cut.departures),
      gap: gapRows(cut.gaps),
      cover: covers,
      forms,
    },
    unloaded: named(
      walks.filter(f => f.abridged !== undefined || f.ids.some(isPlaceholder)),
    ),
    misplaced: named(walks.filter(f => strayPlaceholders(f) > 0)),
    placeholders: walks.reduce(
      (sum, f) => sum + f.ids.filter(isPlaceholder).length,
      0,
    ),
  }
}

const expand = () => true

const KINDS = ['fragment', 'departure', 'gap', 'cover', 'forms'] as const

const totals = {
  windows: 0,
  cuts: 0,
  fragments: 0,
  unloadedFragments: 0,
  placeholders: 0,
  departures: 0,
  gaps: 0,
  covers: 0,
}
const found = {
  stated: new Map<string, number>(),
  full: new Map<string, number>(),
}
const differing = { stated: new Set<string>(), full: new Set<string>() }

const count = (counts: Map<string, number>) =>
  [...counts.values()].reduce((sum, n) => sum + n, 0)

function report(
  name: 'stated' | 'full',
  kind: string,
  where: string,
  lines: string[],
) {
  for (const line of lines) {
    print(`DIFF ${name} ${kind} ${where} ${line}`)
  }
  if (lines.length > 0) {
    found[name].set(kind, (found[name].get(kind) ?? 0) + lines.length)
    differing[name].add(where)
  }
}

// A window's cut that expands nothing against those that expand every walk,
// one cut held at a time
async function compare(
  window: Window,
  contained: boolean,
  cut: (opts: {
    expand?: () => boolean
    expandStated?: boolean
  }) => Promise<Cut>,
) {
  const where = `${window.refName}:${window.start}-${window.end}${contained ? ' contained' : ''}`
  const unloaded = summary(await cut({}), window)
  for (const name of against) {
    const loaded = summary(
      await cut({ expand, expandStated: name === 'stated' }),
      window,
    )
    for (const kind of KINDS) {
      report(
        name,
        kind,
        where,
        differences(unloaded.rows[kind], loaded.rows[kind]),
      )
    }
    report(name, 'unloaded-piece-in-loaded-cut', where, loaded.unloaded)
  }
  report(
    against[0] ?? 'stated',
    'placeholder-misplaced',
    where,
    unloaded.misplaced,
  )
  totals.cuts++
  totals.fragments += unloaded.rows.fragment.size
  totals.unloadedFragments += unloaded.unloaded.length
  totals.placeholders += unloaded.placeholders
  totals.departures += unloaded.rows.departure.size
  totals.gaps += unloaded.rows.gap.size
  totals.covers += unloaded.rows.cover.size
}

async function abridgedTiles(reader: WalkReader) {
  const sequences = (await reader.referenceSequences()).filter(
    name => onlySequence === undefined || name === onlySequence,
  )
  const tiles: Window[] = []
  for (const refName of sequences) {
    const units = await reader.units(refName)
    const abridged = new Set(
      units.filter(u => u.kind === 'blocks').map(u => u.name),
    )
    for (const u of units) {
      if (u.kind === 'tile' && abridged.has(u.name)) {
        tiles.push({ refName, start: u.start, end: u.end })
      }
    }
  }
  return tiles
}

async function main() {
  const reader = WalkReader.open(path, location =>
    /^https?:/.test(location)
      ? new RemoteFile(location)
      : new LocalFile(location),
  )
  const tiles = await abridgedTiles(reader)
  print(
    `${path}: ${tiles.length} abridged tiles; windows ${firstWindow + 1} to ${windowCount} of ${width} bp, seed ${seed}, context ${context} plain and ${containedContext} contained`,
  )
  const next = random(seed)
  for (let k = 0; k < windowCount && tiles.length > 0; k++) {
    const tile = tiles[Math.floor(next() * tiles.length)]!
    const start =
      tile.start +
      Math.floor(next() * Math.max(1, tile.end - tile.start - width))
    const window = { refName: tile.refName, start, end: start + width }
    if (k >= firstWindow) {
      for (const contained of [false, true]) {
        await compare(window, contained, expansion =>
          reader.cut(window.refName, window, {
            context: contained ? containedContext : context,
            contained,
            links: false,
            stepBudget: Infinity,
            byteBudget: Infinity,
            departureByteBudget: Infinity,
            ...expansion,
          }),
        )
      }
      totals.windows++
      process.stderr.write(
        `window ${k + 1} of ${windowCount}, ${window.refName}:${window.start}-${window.end}: ${count(found.stated)} and ${count(found.full)} differences so far, ${Math.round(process.memoryUsage().rss / 1e6)} MB\n`,
      )
    }
  }
  print(
    `compared: ${totals.windows} windows, ${totals.cuts} cuts, ${totals.fragments} fragments (${totals.unloadedFragments} with unloaded pieces, ${totals.placeholders} placeholder steps), ${totals.departures} departures, ${totals.gaps} gaps, ${totals.covers} haplotype covers`,
  )
  for (const name of against) {
    print(
      `differences from the ${name} loaded cut: ${count(found[name])} in ${differing[name].size} cuts${[...found[name]].map(([kind, n]) => `, ${kind} ${n}`).join('')}`,
    )
  }
  process.exitCode = count(found.stated) > 0 ? 1 : count(found.full) > 0 ? 2 : 0
}

await main()
