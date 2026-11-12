// Holds a walk file's unloaded cuts equal to its loaded ones. Over random
// windows on its abridged tiles, it compares the cut that expands no walk,
// whose abridged pieces come from their `blocks` rows, with the cut that
// expands every walk: plain, contained, and contained with context.
//
// Compared, each row counted as often as it comes: the fragments with their
// haplotype offsets, bp trimmed and steps, a step on the reference as its
// place, strand and node and a run of steps off it as its bp; the
// departures with the bp of the steps they hold; the gaps; the reference
// nodes kept; the links between reference nodes, and each link with an end
// off the reference as often as the cut walks it; each haplotype's cover and
// departures; and the structural forms. A loaded cut that holds a
// placeholder anywhere, a pair of steps a cut walks with no link, and a walk
// with a placeholder anywhere in the unloaded cut's tables whose walks and
// departures state none of its unloaded steps are differences too. So are a
// placeholder id that stands for two runs of its cut, or a run of it with two
// ids, and an unloaded cut that a fresh reader of the file makes otherwise,
// in its tables or anything else it returns: the reader compared has by then
// made every cut before it. What is set aside is printed with the counts.
//
// A contained cut keeps `context` bp of a walk past its steps on the window.
// A loaded cut keeps whole nodes, so it runs on to the end of the node that
// straddles `context`; an unloaded one keeps `context` bp of a placeholder.
// There the two may differ by less than that node, and by nothing else: the
// fragments must agree once the loaded one's end is drawn in by that much,
// and the gaps must agree outside it.
//
// A run of steps off a tile's reference can hold a node the file's reference
// passes elsewhere, where the file's `revisited` is above 0. A loaded cut keeps
// such a node and its links, and an unloaded one cannot: its run is one
// placeholder. Where the unloaded cut keeps no such node and the loaded one
// keeps it by a step inside a run the unloaded one holds as a placeholder, the
// node is set aside from the reference nodes kept, and its links are compared
// as a node's off the reference; each one set aside is printed as ASIDE.
//
//   node scripts/blocks-equality.ts <walks.gz or URL> [--windows 300]
//     [--seed 1] [--width 1000] [--context 1000] [--sequence <name>]
//     [--from 0] [--haplotypes 64] [--at <sequence>:<start>-<end>,...]
//
// Prints every difference, then the counts, and exits 1 on any. `--from k`
// skips the first k windows of the seed's; `--at` cuts the windows it lists
// in place of random ones. Memory is bounded by one cut of
// `--haplotypes` haplotypes at a time: each window is cut once per that many
// of the file's haplotypes, as a walk's rows do not depend on the other
// walks cut; 0 cuts them all together. Progress goes to stderr.
import { LocalFile, RemoteFile } from 'generic-filehandle2'

import {
  WalkFile,
  WalkReader,
  isPlaceholder,
  structuralForms,
  walkDepartures,
  walksKept,
} from '../src/index.ts'

import {
  cutFacts,
  digestOf,
  idsNotOneToOne,
  placeholderRuns,
} from '../src/testing.ts'

import type { GraphTables, WalkFragment, WalkGap } from '../src/index.ts'

type Cut = Awaited<ReturnType<WalkReader['cut']>>
type Rows = Map<string, string[]>

interface Window {
  refName: string
  start: number
  end: number
}

interface Mode {
  name: string
  contained: boolean
  context: number
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
const onlySequence = flag('sequence')
const firstWindow = Number(flag('from') ?? 0)
const perCut = Number(flag('haplotypes') ?? 64)
const stated = (flag('at') ?? '')
  .split(',')
  .filter(Boolean)
  .map(at => {
    const [, refName = '', start = '', end = ''] =
      /^(.*):(\d+)-(\d+)$/.exec(at) ?? []
    return { refName, start: Number(start), end: Number(end) }
  })

const MODES: Mode[] = [
  { name: 'plain', contained: false, context },
  { name: 'contained', contained: true, context: 0 },
  { name: `contained+${context}`, contained: true, context },
]

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

const add = (rows: Rows, key: string, row: string) => {
  const held = rows.get(key)
  if (held) {
    held.push(row)
  } else {
    rows.set(key, [row])
  }
}

// A fragment as compared: its steps collapsed, and the bp of its first and
// last steps as cut
function fragmentOf(cut: Cut, f: WalkFragment) {
  const length = (i: number) => {
    const node = cut.graph.nodes.get(f.ids[i]!)
    return node ? node.end - node.start : 0
  }
  const steps: string[] = []
  let own = -1
  for (let i = 0; i <= f.ids.length; i++) {
    if (i < f.ids.length && Number.isNaN(f.refStarts[i])) {
      own = Math.max(own, 0) + length(i)
    } else {
      if (own >= 0) {
        steps.push(`off:${own}`)
      }
      own = -1
      if (i < f.ids.length) {
        steps.push(
          `${f.refStarts[i]}${f.against[i] === 1 ? '-' : '+'}${f.rev[i] === 1 ? '<' : '>'}${f.ids[i]}`,
        )
      }
    }
  }
  return {
    name: f.name,
    fragStart: f.fragStart,
    hapStart: f.hapStart,
    hapEnd: f.hapEnd,
    trimmedBefore: f.trimmedBefore,
    trimmedAfter: f.trimmedAfter,
    steps,
    firstBp: length(0),
    lastBp: length(f.ids.length - 1),
  }
}

type Fragment = ReturnType<typeof fragmentOf>

const fragmentRow = (f: Fragment) =>
  `${f.hapStart}-${f.hapEnd} trimmed ${f.trimmedBefore} and ${f.trimmedAfter}: ${f.steps.join(' ')}`

const strand = (reversed: boolean) => (reversed ? '-' : '+')

const gapRow = (g: WalkGap) =>
  `edge ${g.edge} ${g.start}-${g.end} beside ${g.fragment ? g.fragment.hapStart : 'none'}: ${g.segments
    .map(s => `${s.kind}:${s.bp}:missing${s.missingBp}:${s.before}..${s.after}`)
    .join(' ')}`

// A placeholder beside another step off the reference, or on it
function misplaced(f: WalkFragment) {
  const offR = (k: number) =>
    k >= 0 && k < f.ids.length && Number.isNaN(f.refStarts[k])
  return f.ids.some(
    (id, i) => isPlaceholder(id) && (offR(i - 1) || offR(i + 1) || !offR(i)),
  )
}

const END = 2 ** 25

// A link as a number: the two node ends it joins, either way round
const pair = (a: number, aRight: boolean, b: number, bRight: boolean) => {
  const x = a * 2 + (aRight ? 1 : 0)
  const y = b * 2 + (bRight ? 1 : 0)
  return x < y ? x * END + y : y * END + x
}

// A cut's links as rows: every link between two reference nodes, and each
// link with an end off the reference once for each time the cut walks it,
// as consecutive steps of a fragment or from a departure's leave node
// through its steps to its rejoin node, which ends it leaves and enters by
// the reference's strand at that pass. A link with an end off the reference
// that the cut walks nowhere joins nodes kept through other walks, and is
// set aside; a pair of steps a fragment or a computed departure walks with
// no link is `missing`. A step `inside` says lies in a run the unloaded cut
// holds as a placeholder is walked as a node off the reference, and its node
// is listed in `inside`
function linkRows(
  cut: Cut,
  window: Window,
  { nodes, links }: GraphTables,
  reversedOnR: (id: number, start: number) => boolean | undefined,
  inside: (name: string, id: number, start: number) => boolean,
) {
  const indexOfName = new Map(nodes.names.map((name, i) => [name, i]))
  const indexOf = (id: number | undefined) => {
    const node = id === undefined ? undefined : cut.graph.nodes.get(id)
    return id === undefined
      ? undefined
      : indexOfName.get(
          node?.placeholder
            ? `${node.refName}:${node.start}-${node.end}`
            : String(id),
        )
  }
  const length = (id: number) => {
    const node = cut.graph.nodes.get(id)
    return node ? node.end - node.start : 0
  }
  const onR = (index: number) =>
    nodes.ranks[index] === 0 && nodes.placeholders?.[index] !== 1
  const end = (index: number, right: boolean) =>
    `${nodes.names[index]}${right ? ']' : '['}`
  const held = new Set<number>()
  const rows: Rows = new Map()
  let offEnds = 0
  for (let k = 0; k < links.from.length; k++) {
    const strands = links.strands[k]!
    const from = links.from[k]!
    const to = links.to[k]!
    const fromRight = (strands & 1) === 0
    const toRight = (strands & 2) !== 0
    held.add(pair(from, fromRight, to, toRight))
    if (onR(from) && onR(to)) {
      const a = end(from, fromRight)
      const b = end(to, toRight)
      add(rows, 'between reference nodes', a < b ? `${a} ${b}` : `${b} ${a}`)
    } else {
      offEnds++
    }
  }
  const walkedOff = new Set<number>()
  const missing: string[] = []
  const insideNodes = new Set<number>()
  // a step as its node's index and whether it counts as on the reference
  const stepOf = (id: number | undefined, name = '', at = Number.NaN) => {
    const index = indexOf(id)
    const within = id !== undefined && !Number.isNaN(at) && inside(name, id, at)
    if (within) {
      insideNodes.add(id)
    }
    return {
      index,
      on: index !== undefined && onR(index) && !within,
    }
  }
  type Step = ReturnType<typeof stepOf>
  const walk = (
    a: Step,
    aRight: boolean,
    b: Step,
    bRight: boolean,
    where: string | undefined,
  ) => {
    if (a.index !== undefined && b.index !== undefined) {
      const key = pair(a.index, aRight, b.index, bRight)
      if (!held.has(key)) {
        if (where !== undefined) {
          missing.push(where)
        }
      } else if (a.on !== b.on) {
        walkedOff.add(key)
        add(
          rows,
          'walked, with an end off the reference',
          a.on ? end(a.index, aRight) : end(b.index, bRight),
        )
      } else if (!a.on) {
        walkedOff.add(key)
      }
    }
  }
  for (const f of cut.fragments) {
    let at = f.hapStart
    const stepAt = (i: number) =>
      stepOf(f.ids[i], f.name, Number.isNaN(f.refStarts[i]) ? at : Number.NaN)
    let last = stepAt(0)
    for (let i = 1; i < f.ids.length; i++) {
      at += length(f.ids[i - 1]!)
      const next = stepAt(i)
      walk(
        last,
        f.rev[i - 1] === 0,
        next,
        f.rev[i] !== 0,
        `${f.name} ${f.ids[i - 1]} to ${f.ids[i]}`,
      )
      last = next
    }
  }
  const startingAt = new Map<number, number>()
  const endingAt = new Map<number, number>()
  for (const id of cut.kept) {
    const node = cut.graph.nodes.get(id)
    if (
      node?.rank === 0 &&
      !node.placeholder &&
      node.refName === window.refName
    ) {
      startingAt.set(node.start, id)
      endingAt.set(node.end, id)
    }
  }
  const reversedAt = (id: number | undefined) =>
    id !== undefined &&
    reversedOnR(id, cut.graph.nodes.get(id)?.start ?? Number.NaN) === true
  for (const d of cut.departures) {
    if (d.steps && d.steps.ids.length > 0) {
      const { ids, rev } = d.steps
      const where = d.fromSteps
        ? `${d.name} departure at ${d.hapStart}`
        : undefined
      const leave = (d.leaveReversed ? startingAt : endingAt).get(d.leavePos)
      const rejoin = (d.rejoinReversed ? endingAt : startingAt).get(d.rejoinPos)
      const own: Step[] = []
      let at = d.hapStart
      for (const id of ids) {
        own.push(stepOf(id, d.name, at))
        at += length(id)
      }
      walk(
        stepOf(leave),
        reversedAt(leave) === d.leaveReversed,
        own[0]!,
        rev[0] !== 0,
        where,
      )
      for (let i = 1; i < ids.length; i++) {
        walk(own[i - 1]!, rev[i - 1] === 0, own[i]!, rev[i] !== 0, where)
      }
      walk(
        own.at(-1)!,
        rev.at(-1) === 0,
        stepOf(rejoin),
        reversedAt(rejoin) !== d.rejoinReversed,
        where,
      )
    }
  }
  return {
    rows,
    missing,
    setAside: offEnds - walkedOff.size,
    insideNodes,
    placeholderNames: nodes.names.filter(
      (_, i) => nodes.placeholders?.[i] === 1,
    ),
  }
}

// The reference's strand at each pass of its nodes, as `id@start`, from the
// ref rows of the tiles holding the positions asked for, each tile read once
function referenceStrands(file: WalkFile) {
  const tilesRead = new Set<number>()
  const reversed = new Map<string, boolean>()
  const cover = async (refName: string, positions: number[]) => {
    const units = (await file.units(refName)).filter(
      u =>
        u.kind === 'tile' &&
        !tilesRead.has(u.offset) &&
        positions.some(at => u.start <= at && at <= u.end),
    )
    if (tilesRead.size + units.length > TILES_HELD) {
      tilesRead.clear()
      reversed.clear()
    }
    const lines = await file.read(units)
    units.forEach((unit, k) => {
      tilesRead.add(unit.offset)
      for (const line of lines[k] ?? []) {
        const cols = line.split('\t')
        if (cols[3] === 'ref') {
          const lengths = (cols[10] ?? '').split(',').map(Number)
          let start = Number(cols[6])
          let id = 0
          for (const [j, token] of (cols[9] ?? '').split(',').entries()) {
            const t = Number(token)
            const r = t - 2 * Math.floor(t / 2)
            id += (t - r) / 2
            reversed.set(`${id}@${start}`, r === 1)
            start += lengths[j]!
          }
        }
      }
    })
  }
  return {
    cover,
    reversedOnR: (id: number, start: number) => reversed.get(`${id}@${start}`),
  }
}

type Strands = ReturnType<typeof referenceStrands>

const TILES_HELD = 4096

// Whether a step on `id` at `start` on its walk's contig lies inside a run
// the unloaded cut holds as a placeholder, on a node of the reference
const insideRun =
  (cut: Cut, unloaded: Unloaded | undefined) =>
  (name: string, id: number, start: number) => {
    const node = cut.graph.nodes.get(id)
    const end = start + (node ? node.end - node.start : 0)
    return (
      unloaded !== undefined &&
      node?.rank === 0 &&
      !node.placeholder &&
      (unloaded.runs.get(name) ?? []).some(([a, b]) => a <= start && end <= b)
    )
  }

// Each walk's stretches of its contig that the unloaded cut holds as
// placeholders, and the reference nodes it keeps
function unloadedRuns(cut: Cut) {
  const runs = new Map<string, [number, number][]>()
  for (const node of cut.graph.nodes.values()) {
    if (node.placeholder) {
      runs.set(node.refName, [
        ...(runs.get(node.refName) ?? []),
        [node.start, node.end],
      ])
    }
  }
  const keptOnReference = new Set<number>()
  for (const id of cut.kept) {
    const node = cut.graph.nodes.get(id)
    if (node?.rank === 0 && !node.placeholder) {
      keptOnReference.add(id)
    }
  }
  return { runs, keptOnReference }
}

type Unloaded = ReturnType<typeof unloadedRuns>

// The walks holding a placeholder among their steps, a departure's or a
// gap's, whose walks' and departures' abridgedSteps state no step
function uncounted({ nodes, walks, departures, gaps }: GraphTables) {
  const holding = new Set<string>()
  const scan = (name: string, steps: Int32Array, from: number, to: number) => {
    for (let j = from; j < to; j++) {
      if (nodes.placeholders?.[steps[j]!] === 1) {
        holding.add(name)
      }
    }
  }
  walks.names.forEach((name, i) => {
    scan(name, walks.steps, walks.offsets[i]!, walks.offsets[i + 1]!)
  })
  departures.names.forEach((name, k) => {
    scan(
      name,
      departures.steps,
      departures.offsets[k]!,
      departures.offsets[k + 1]!,
    )
  })
  gaps?.names.forEach((name, g) => {
    const { stepOffsets, steps } = gaps.segments
    scan(
      name,
      steps,
      stepOffsets[gaps.segmentOffsets[g]!]!,
      stepOffsets[gaps.segmentOffsets[g + 1]!]!,
    )
  })
  const counted = new Map<string, number>()
  const sum = (names: string[], steps: Float64Array | undefined) => {
    names.forEach((name, i) => {
      const n = steps?.[i] ?? Number.NaN
      if (!Number.isNaN(n)) {
        counted.set(name, (counted.get(name) ?? 0) + n)
      }
    })
  }
  sum(walks.names, walks.abridgedSteps)
  sum(departures.names, departures.abridgedSteps)
  return [...holding].filter(name => !((counted.get(name) ?? 0) > 0))
}

// Where a loaded cut holds a placeholder: in a fragment, a departure's
// steps, a gap's steps or the nodes it keeps
function placeholdersIn(cut: Cut) {
  const out: string[] = []
  for (const f of cut.fragments) {
    if (f.abridged !== undefined || f.ids.some(isPlaceholder)) {
      out.push(`fragment of ${f.name} from ${f.hapStart}`)
    }
  }
  for (const d of cut.departures) {
    if (d.abridged !== undefined || d.steps?.ids.some(isPlaceholder) === true) {
      out.push(`departure of ${d.name} at ${d.hapStart}`)
    }
  }
  for (const g of cut.gaps) {
    if (g.segments.some(x => x.ids.some(isPlaceholder))) {
      out.push(`gap of ${g.name} from ${g.start}`)
    }
  }
  let kept = 0
  for (const id of cut.kept) {
    kept += isPlaceholder(id) ? 1 : 0
  }
  if (kept > 0 || cut.placeholders > 0) {
    out.push(`${kept} kept, ${cut.placeholders} counted`)
  }
  return out
}

// A fragment whose unloadedSteps is not a count of its steps: one with no
// unloaded piece that states it, or one whose count is under its
// placeholders or over its steps
function unloadedStepsOff(cut: Cut) {
  return cut.fragments.flatMap(f => {
    const placeholders = f.ids.filter(isPlaceholder).length
    const fits =
      f.abridged === undefined
        ? f.unloadedSteps === undefined && placeholders === 0
        : f.unloadedSteps !== undefined &&
          f.unloadedSteps >= placeholders &&
          f.unloadedSteps <= f.ids.length
    return fits
      ? []
      : [`${f.name} from ${f.hapStart}: ${f.unloadedSteps} of ${f.ids.length}`]
  })
}

// What a cut states, as rows to compare; the cut itself is let go. Given
// what the unloaded cut of the same window holds, the reference nodes only a
// loaded cut can keep are set aside, with their links
function summary(
  cut: Cut,
  window: Window,
  strands: Strands,
  unloaded?: Unloaded,
) {
  const tables = cut.graph.tables(
    cut.kept,
    cut.fragments,
    cut.departures,
    cut.gaps,
  )
  const links = linkRows(
    cut,
    window,
    tables,
    strands.reversedOnR,
    insideRun(cut, unloaded),
  )
  const aside = new Set(
    [...links.insideNodes].filter(
      id => unloaded?.keptOnReference.has(id) === false,
    ),
  )
  const walks = cut.fragments.filter(f => f.name !== window.refName)
  const fragments = new Map<string, Fragment[]>()
  for (const f of walks) {
    const key = `${f.name} line ${f.fragStart}`
    const held = fragments.get(key) ?? []
    fragments.set(key, held)
    held.push(fragmentOf(cut, f))
  }
  const departures: Rows = new Map()
  for (const d of cut.departures) {
    const steps = d.steps
      ? d.steps.ids.reduce((sum, id) => {
          const node = cut.graph.nodes.get(id)
          return sum + (node ? node.end - node.start : 0)
        }, 0)
      : 'none'
    add(
      departures,
      `${d.name} at ${d.hapStart}`,
      `own ${d.ownBp} missing ${d.missingBp} leave ${d.leavePos}${strand(d.leaveReversed)}x${d.leavePasses} rejoin ${d.rejoinPos}${strand(d.rejoinReversed)}x${d.rejoinPasses} ${d.fromSteps ? 'steps' : 'record'} holding ${steps}`,
    )
  }
  const gaps = new Map<string, WalkGap[]>()
  for (const g of cut.gaps) {
    const held = gaps.get(g.name) ?? []
    gaps.set(g.name, held)
    held.push(g)
  }
  const haplotypes = walkDepartures(
    {
      fragments: cut.fragments,
      departures: cut.departures,
      nodes: cut.graph.nodes,
    },
    window.refName,
  )
  const covers: Rows = new Map()
  for (const h of haplotypes) {
    add(
      covers,
      h.haplotype,
      `${JSON.stringify(h.covered)} ${JSON.stringify(h.departures)}`,
    )
  }
  const forms = structuralForms(haplotypes, {
    window: { start: window.start, end: window.end },
  })
  const kept: Rows = new Map()
  for (const id of cut.kept) {
    const node = cut.graph.nodes.get(id)
    if (node?.rank === 0 && !node.placeholder && !aside.has(id)) {
      add(kept, 'on the reference', String(id))
    }
  }
  let offReference = 0
  for (const id of cut.kept) {
    offReference += cut.graph.nodes.get(id)?.rank === 0 ? 0 : 1
  }
  return {
    fragments,
    departures,
    gaps,
    covers,
    forms: new Map([
      [
        'forms',
        [
          JSON.stringify({
            sites: forms.sites,
            forms: forms.forms.map(f => ({
              ...f,
              members: [...f.members].sort(),
            })),
            unplaced: [...forms.unplaced].sort(),
          }),
        ],
      ],
    ]),
    kept,
    links: links.rows,
    linksMissing: links.missing,
    linksSetAside: links.setAside,
    offReference,
    placeholdersHeld: placeholdersIn(cut),
    unloadedStepsOff: unloadedStepsOff(cut),
    uncounted: uncounted(tables),
    runs: unloaded ? undefined : unloadedRuns(cut),
    aside: [...aside],
    insideNodes: links.insideNodes,
    placeholderNames: links.placeholderNames,
    idsNotOneToOne: idsNotOneToOne(placeholderRuns(cut)),
    unloaded: walks.filter(
      f => f.abridged !== undefined || f.ids.some(isPlaceholder),
    ).length,
    misplaced: walks.filter(misplaced).map(f => `${f.name} ${f.hapStart}`),
    placeholders: walks.reduce(
      (sum, f) => sum + f.ids.filter(isPlaceholder).length,
      0,
    ),
  }
}

// Where two lists first differ, with the items around
function firstDifference(a: string[], b: string[]) {
  let at = 0
  while (at < a.length && at < b.length && a[at] === b[at]) {
    at++
  }
  const around = (items: string[]) =>
    items
      .slice(Math.max(0, at - 1), at + 2)
      .join(' | ')
      .slice(0, 600) || 'nothing'
  return `${a.length} against ${b.length}, at ${at}: ${around(a)} AGAINST ${around(b)}`
}

// Each key whose rows differ, every row counted as often as it comes
function differences(unloaded: Rows, loaded: Rows) {
  const keys = new Set(unloaded.keys())
  for (const key of loaded.keys()) {
    keys.add(key)
  }
  const out: string[] = []
  for (const key of [...keys].sort()) {
    const a = [...(unloaded.get(key) ?? [])].sort()
    const b = [...(loaded.get(key) ?? [])].sort()
    if (a.length !== b.length || a.some((row, i) => row !== b[i])) {
      out.push(`${key}: ${firstDifference(a, b)}`)
    }
  }
  return out
}

// The bp an own run at a fragment's edge states, or -1 where the edge is on
// the reference
const edgeOwn = (step: string | undefined) =>
  step?.startsWith('off:') === true ? Number(step.slice(4)) : -1

// How far a loaded fragment's edge runs past the unloaded one's, where the
// loaded cut kept the node that straddles `context` whole: less than that
// node, from an unloaded edge at `context` or more. NaN is any other
// difference
function overrun(unloaded: number, loaded: number, nodeBp: number, bp: number) {
  const more = loaded - unloaded
  return more === 0
    ? 0
    : more > 0 && unloaded >= bp && more < nodeBp
      ? more
      : Number.NaN
}

// The loaded fragment with its ends drawn in by `head` and `tail` bp
const drawnIn = (f: Fragment, head: number, tail: number) => ({
  ...f,
  hapStart: f.hapStart + head,
  hapEnd: f.hapEnd - tail,
  trimmedBefore: f.trimmedBefore + head,
  trimmedAfter: f.trimmedAfter + tail,
  steps: f.steps.map((step, i) =>
    i === 0 && head > 0
      ? `off:${edgeOwn(step) - head}`
      : i === f.steps.length - 1 && tail > 0
        ? `off:${edgeOwn(step) - tail}`
        : step,
  ),
})

// A walk's contig as its gaps' segments and its fragments, in order, less
// the stretches in `left`, touching stretches of one label as one
function tape(
  gaps: WalkGap[],
  fragments: Fragment[],
  left: [number, number][],
) {
  const parts: { start: number; end: number; label: string }[] = []
  for (const f of fragments) {
    parts.push({ start: f.hapStart, end: f.hapEnd, label: 'fragment' })
  }
  for (const g of gaps) {
    let at = g.start
    for (const s of g.segments) {
      parts.push({
        start: at,
        end: at + s.bp,
        label: `${s.kind} missing ${s.missingBp} ${s.before}..${s.after}`,
      })
      at += s.bp
    }
  }
  const outside = parts.flatMap(part => {
    let pieces = [part]
    for (const [from, to] of left) {
      pieces = pieces.flatMap(x =>
        x.end <= from || x.start >= to
          ? [x]
          : [
              { ...x, end: Math.min(x.end, from) },
              { ...x, start: Math.max(x.start, to) },
            ].filter(y => y.end > y.start),
      )
    }
    return pieces
  })
  outside.sort((a, b) => a.start - b.start || a.end - b.end)
  const out: typeof parts = []
  for (const part of outside) {
    const last = out.at(-1)
    if (
      last?.label === part.label &&
      last.end >= part.start &&
      part.end > part.start
    ) {
      last.end = Math.max(last.end, part.end)
    } else {
      out.push({ ...part })
    }
  }
  return out.map(x => `${x.start}-${x.end} ${x.label}`)
}

type Summary = ReturnType<typeof summary>

// The fragments and gaps of the two cuts, walk by walk. Where a contained
// cut's context ends inside a node the loaded cut kept whole, the fragments
// are compared with the loaded one's end drawn in and the gaps outside that
// stretch
function fragmentsAndGaps(unloaded: Summary, loaded: Summary, mode: Mode) {
  const fragment: string[] = []
  const gap: string[] = []
  let overruns = 0
  const keys = new Set(unloaded.fragments.keys())
  for (const key of loaded.fragments.keys()) {
    keys.add(key)
  }
  const drawn = new Map<string, [number, number][]>()
  const drawnFragments = new Map<string, Fragment[]>()
  for (const key of [...keys].sort()) {
    const byStart = (a: Fragment, b: Fragment) => a.hapStart - b.hapStart
    const a = [...(unloaded.fragments.get(key) ?? [])].sort(byStart)
    const b = [...(loaded.fragments.get(key) ?? [])].sort(byStart)
    const adjusted = b.map((f, i) => {
      const u = a[i]
      const allowed = mode.contained && mode.context > 0 && u !== undefined
      const head = allowed
        ? overrun(
            edgeOwn(u.steps[0]),
            edgeOwn(f.steps[0]),
            f.firstBp,
            mode.context,
          )
        : 0
      const tail = allowed
        ? overrun(
            edgeOwn(u.steps.at(-1)),
            edgeOwn(f.steps.at(-1)),
            f.lastBp,
            mode.context,
          )
        : 0
      return {
        f,
        head: Number.isNaN(head) ? 0 : head,
        tail: Number.isNaN(tail) ? 0 : tail,
      }
    })
    const rowsA = a.map(fragmentRow)
    const rowsB = adjusted.map(x => fragmentRow(drawnIn(x.f, x.head, x.tail)))
    if (rowsA.length !== rowsB.length || rowsA.some((r, i) => r !== rowsB[i])) {
      fragment.push(`${key}: ${firstDifference(rowsA, rowsB)}`)
    }
    const name = a[0]?.name ?? b[0]?.name ?? ''
    for (const x of adjusted) {
      overruns += (x.head > 0 ? 1 : 0) + (x.tail > 0 ? 1 : 0)
      if (x.head > 0) {
        drawn.set(name, [
          ...(drawn.get(name) ?? []),
          [x.f.hapStart, x.f.hapStart + x.head],
        ])
      }
      if (x.tail > 0) {
        drawn.set(name, [
          ...(drawn.get(name) ?? []),
          [x.f.hapEnd - x.tail, x.f.hapEnd],
        ])
      }
    }
    drawnFragments.set(name, [
      ...(drawnFragments.get(name) ?? []),
      ...adjusted.map(x => drawnIn(x.f, x.head, x.tail)),
    ])
  }
  const names = new Set(unloaded.gaps.keys())
  for (const name of loaded.gaps.keys()) {
    names.add(name)
  }
  for (const name of [...names].sort()) {
    const a = unloaded.gaps.get(name) ?? []
    const b = loaded.gaps.get(name) ?? []
    const left = drawn.get(name)
    const walkFragments = (s: Summary) =>
      [...s.fragments].flatMap(([, fs]) => fs.filter(f => f.name === name))
    const rowsA = left
      ? tape(a, walkFragments(unloaded), left)
      : a.map(gapRow).sort()
    const rowsB = left
      ? tape(b, drawnFragments.get(name) ?? [], left)
      : b.map(gapRow).sort()
    if (rowsA.length !== rowsB.length || rowsA.some((r, i) => r !== rowsB[i])) {
      gap.push(`${name}: ${firstDifference(rowsA, rowsB)}`)
    }
  }
  return { fragment, gap, overruns }
}

const totals = {
  windows: 0,
  cuts: 0,
  fragments: 0,
  unloadedFragments: 0,
  placeholders: 0,
  departures: 0,
  gaps: 0,
  covers: 0,
  links: 0,
  kept: 0,
  overruns: 0,
}
// what the two cuts hold that is not compared row for row
const setAside = {
  unloadedLinks: 0,
  loadedLinks: 0,
  unloadedOffReference: 0,
  loadedOffReference: 0,
  insideNodes: 0,
  insideLinks: 0,
}

// The links between reference nodes the loaded cut's tables hold and the
// unloaded cut's do not, with an end on a node the loaded cut steps on inside
// a run the unloaded cut holds as a placeholder, taken out of `loaded`
function linksOnlyInside(unloaded: Rows, loaded: Rows, inside: Set<number>) {
  const key = 'between reference nodes'
  const theirs = new Set(unloaded.get(key) ?? [])
  const taken: string[] = []
  const kept: string[] = []
  for (const row of loaded.get(key) ?? []) {
    const ends = row.split(' ').map(end => Number(end.slice(0, -1)))
    if (!theirs.has(row) && ends.some(id => inside.has(id))) {
      taken.push(row)
    } else {
      kept.push(row)
    }
  }
  if (loaded.has(key)) {
    loaded.set(key, kept)
  }
  return taken
}
const found = new Map<string, number>()
const differing = new Set<string>()
// A cut's tables, and everything it returns with them, each as a digest
function digests(cut: Cut) {
  const facts = cutFacts(cut)
  return { tables: digestOf(facts.tables), cut: digestOf(facts) }
}

const count = (counts: Map<string, number>) =>
  [...counts.values()].reduce((sum, n) => sum + n, 0)

function report(kind: string, where: string, lines: string[]) {
  for (const line of lines) {
    print(`DIFF ${kind} ${where} ${line}`)
  }
  if (lines.length > 0) {
    found.set(kind, (found.get(kind) ?? 0) + lines.length)
    differing.add(where)
  }
}

const size = (rows: Rows) =>
  [...rows.values()].reduce((sum, list) => sum + list.length, 0)

const expand = () => true

// A window's cut that expands nothing against the one that expands every
// walk, and against itself made again by a fresh reader, one cut held at a
// time
async function compare(
  window: Window,
  where: string,
  mode: Mode,
  cut: (expanding: boolean, fresh?: boolean) => Promise<Cut>,
  revisited: boolean,
  strands: Strands,
) {
  const summarize = async (held: Cut, runs?: Unloaded) => {
    await strands.cover(
      window.refName,
      held.departures.flatMap(d => [d.leavePos, d.rejoinPos]),
    )
    return summary(held, window, strands, runs)
  }
  const digested = async (held: Cut) => ({
    digests: digests(held),
    ...(await summarize(held)),
  })
  const unloaded = await digested(await cut(false))
  const again = digests(await cut(false, true))
  const loaded = await summarize(
    await cut(true),
    revisited ? unloaded.runs : undefined,
  )
  const asideLinks = linksOnlyInside(
    unloaded.links,
    loaded.links,
    loaded.insideNodes,
  )
  if (loaded.aside.length > 0 || asideLinks.length > 0) {
    setAside.insideNodes += loaded.aside.length
    setAside.insideLinks += asideLinks.length
    print(
      `ASIDE reference-node-inside-unloaded-run ${where} nodes ${loaded.aside.join(' ')} links ${asideLinks.join(' | ')}`,
    )
  }
  const both = fragmentsAndGaps(unloaded, loaded, mode)
  report('fragment', where, both.fragment)
  report('gap', where, both.gap)
  report(
    'departures',
    where,
    differences(unloaded.departures, loaded.departures),
  )
  report('covers', where, differences(unloaded.covers, loaded.covers))
  report('forms', where, differences(unloaded.forms, loaded.forms))
  report('kept', where, differences(unloaded.kept, loaded.kept))
  report('links', where, differences(unloaded.links, loaded.links))
  report('link-missing-unloaded', where, unloaded.linksMissing)
  report('link-missing-loaded', where, loaded.linksMissing)
  report('placeholder-in-loaded-cut', where, loaded.placeholdersHeld)
  report('unloaded-steps', where, [
    ...unloaded.unloadedStepsOff,
    ...loaded.unloadedStepsOff,
  ])
  report('placeholder-misplaced', where, unloaded.misplaced)
  report('placeholder-uncounted', where, unloaded.uncounted)
  report('placeholder-ids-not-one-to-one', where, [
    ...unloaded.idsNotOneToOne,
    ...loaded.idsNotOneToOne,
  ])
  report(
    'cut-differs-on-a-fresh-reader',
    where,
    again.cut === unloaded.digests.cut
      ? []
      : [
          again.tables === unloaded.digests.tables
            ? 'the tables the same, the ids or rows the cut returns not'
            : 'the tables',
        ],
  )
  totals.cuts++
  totals.fragments += [...unloaded.fragments.values()].reduce(
    (sum, list) => sum + list.length,
    0,
  )
  totals.unloadedFragments += unloaded.unloaded
  totals.placeholders += unloaded.placeholders
  totals.departures += size(unloaded.departures)
  totals.gaps += [...unloaded.gaps.values()].reduce(
    (sum, list) => sum + list.length,
    0,
  )
  totals.covers += size(unloaded.covers)
  totals.links += size(unloaded.links)
  setAside.unloadedLinks += unloaded.linksSetAside
  setAside.loadedLinks += loaded.linksSetAside
  setAside.unloadedOffReference += unloaded.offReference
  setAside.loadedOffReference += loaded.offReference
  totals.kept += size(unloaded.kept)
  totals.overruns += both.overruns
  return unloaded.placeholderNames
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

const open = (location: string) =>
  /^https?:/.test(location) ? new RemoteFile(location) : new LocalFile(location)

async function main() {
  const file = new WalkFile(
    open(path),
    open(path.replace(/(?=$|[?#])/, '.idx')),
  )
  const reader = new WalkReader(file)
  const strands = referenceStrands(file)
  const tiles = await abridgedTiles(reader)
  const { haplotypes, revisited = 0 } = await reader.header()
  const batches =
    perCut > 0 && haplotypes.length > perCut
      ? Array.from({ length: Math.ceil(haplotypes.length / perCut) }, (_, b) =>
          haplotypes.slice(b * perCut, (b + 1) * perCut),
        )
      : [haplotypes]
  const windows: (Window | undefined)[] = []
  if (stated.length > 0) {
    for (const at of stated) {
      windows.push(at)
    }
  } else {
    const next = random(seed)
    for (let k = 0; k < windowCount && tiles.length > 0; k++) {
      const tile = tiles[Math.floor(next() * tiles.length)]!
      const start =
        tile.start +
        Math.floor(next() * Math.max(1, tile.end - tile.start - width))
      windows.push(
        k >= firstWindow
          ? { refName: tile.refName, start, end: start + width }
          : undefined,
      )
    }
  }
  print(
    `${path}: ${tiles.length} abridged tiles, ${haplotypes.length} haplotypes cut ${batches.length > 1 ? `${perCut} at a time` : 'together'}; ${stated.length > 0 ? `the ${stated.length} windows stated` : `windows ${firstWindow + 1} to ${windowCount} of ${width} bp, seed ${seed}`}, ${MODES.map(m => m.name).join(', ')}`,
  )
  for (const [k, window] of windows.entries()) {
    if (window) {
      for (const mode of MODES) {
        for (const [b, batch] of batches.entries()) {
          const where = `${window.refName}:${window.start}-${window.end} ${mode.name}${batches.length > 1 ? ` haplotypes ${b * perCut + 1}-${b * perCut + batch.length}` : ''}`
          const names = await compare(
            window,
            where,
            mode,
            (expanding, fresh = false) =>
              (fresh ? new WalkReader(file) : reader).cut(
                window.refName,
                window,
                {
                  keep:
                    batches.length > 1
                      ? walksKept(batch, window.refName)
                      : undefined,
                  context: mode.context,
                  contained: mode.contained,
                  stepBudget: Infinity,
                  byteBudget: Infinity,
                  departureByteBudget: Infinity,
                  ...(expanding ? { expand } : {}),
                },
              ),
            revisited > 0,
            strands,
          )
          report(
            'placeholder-name-as-a-node',
            where,
            names.filter(name => /^-?\d+$/.test(name)),
          )
        }
      }
      totals.windows++
      process.stderr.write(
        `window ${k + 1} of ${windows.length}, ${window.refName}:${window.start}-${window.end}: ${count(found)} differences so far, ${Math.round(process.memoryUsage().rss / 1e6)} MB\n`,
      )
    }
  }
  print(
    `compared: ${totals.windows} windows, ${totals.cuts} pairs of cuts, ${totals.fragments} fragments (${totals.unloadedFragments} with unloaded pieces, ${totals.placeholders} placeholder steps), ${totals.departures} departures, ${totals.gaps} gaps, ${totals.covers} haplotype covers, ${totals.kept} reference nodes kept, ${totals.links} links`,
  )
  print(
    `fragment ends where the loaded cut ran on to the end of the node straddling the context: ${totals.overruns}`,
  )
  print(
    `set aside, not compared: ${setAside.unloadedLinks} links of the unloaded cuts and ${setAside.loadedLinks} of the loaded with an end off the reference that the cut walks nowhere; the identity of the ${setAside.unloadedOffReference} nodes and placeholders off the reference the unloaded cuts keep and the ${setAside.loadedOffReference} nodes the loaded keep, whose bp are compared run by run; links between two nodes off the reference; each fragment's list of unloaded pieces and its unloadedSteps, which only an unloaded cut has`,
  )
  print(
    `set aside, by design: ${setAside.insideNodes} reference nodes and ${setAside.insideLinks} links between reference nodes that only a loaded cut keeps, by a step inside a run the unloaded cut holds as a placeholder, on a node the file's reference passes elsewhere (revisited:i:${revisited}); such a step's links are walked as a node's off the reference`,
  )
  print(
    `differences from the loaded cut: ${count(found)} in ${differing.size} pairs of cuts${[...found].map(([kind, n]) => `, ${kind} ${n}`).join('')}`,
  )
  process.exitCode = count(found) > 0 ? 1 : 0
}

await main()
