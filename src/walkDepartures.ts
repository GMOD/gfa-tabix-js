import { panSNHaplotypeOrSample } from './util.ts'

import type { WalkDeparture } from './walkDeparture.ts'
import type { WalkReader, WalkReaderHeader } from './walkReader.ts'
import type { WalkFragment, WalkNode } from './walkRows.ts'

/**
 * Where a haplotype leaves the reference path R and rejoins it: R's bases
 * `refStart`–`refEnd` against `altLength` of its own. `invertedBp` counts the
 * steps along R against its orientation. A step back along R walks a tandem
 * copy again: the departure is then a point at the junction, `repeatedBp`
 * counts the bp walked again, and `copied` spans them on R.
 */
export interface Departure {
  refStart: number
  refEnd: number
  altLength: number
  invertedBp: number
  repeatedBp: number
  copied?: [number, number]
}

/**
 * A `sample#haplotype`'s departures, or a sample's where its paths' names
 * have no haplotype field, by refStart. `covered` spans the R
 * exactly one of its paths reaches, each span within one path.
 */
export interface HaplotypeDepartures {
  haplotype: string
  covered: [number, number][]
  departures: Departure[]
}

// A cut's fragments, its departures and its nodes, as WalkReader.cut returns
// them
export interface DepartureSource {
  fragments: WalkFragment[]
  departures: WalkDeparture[]
  nodes: Map<number, WalkNode>
}

// What a path does along its contig: follows R over [lo, hi), along or
// against it, or runs `bp` of its own
type Event =
  | { kind: 'on'; lo: number; hi: number; reversed: boolean }
  | { kind: 'own'; bp: number }

// A point on a path's contig: where it is on R until (a fragment's last step
// or a departure's leave, `order` 0) or from (a fragment's first step or a
// departure's rejoin, 1). `own` is the bp of the path's own that follow a
// leave
interface Anchor {
  hap: number
  order: 0 | 1
  pos: number
  reversed: boolean
  exit?: true
  own?: number
}

// A fragment or departure as the bp of its contig it holds
type Item = { hap: [number, number] } & (
  | { fragment: WalkFragment }
  | { departure: WalkDeparture }
)

const onR = (start: number) => !Number.isNaN(start)

/**
 * Each haplotype's departures from R, `refName`, read off a cut: its
 * departures, every one over the drawn stretch, give where each path leaves
 * R and rejoins it and how much of its own it runs between, and its
 * fragments' first and last steps give where it is on R otherwise. A
 * fragment with the departures at its ends is a chain, as is a departure
 * touching no fragment, such as one jumping over the window; each reads in
 * its own orientation, mostly reversed reading reversed, and a stretch of R
 * walked against it reads as inverted, a jump back as a tandem copy, and
 * sequence between a contig's W lines as the path's own. R that two chains
 * of a haplotype reach, such as a contig's second visit to the window, has
 * no call. A piece the cut did not load reads as any other, from its steps
 * on R and the bp of its own between them. The walks on R's own haplotype
 * are left out.
 */
export function walkDepartures(
  { fragments, departures, nodes }: DepartureSource,
  refName: string,
) {
  const length = (id: number) => {
    const node = nodes.get(id)
    return node ? node.end - node.start : 0
  }
  const reference = panSNHaplotypeOrSample(refName)
  const byPath = new Map<string, Item[]>()
  const add = (name: string, item: Item) => {
    const items = byPath.get(name)
    if (items) {
      items.push(item)
    } else {
      byPath.set(name, [item])
    }
  }
  for (const f of fragments) {
    add(f.name, {
      hap: [f.hapStart - f.trimmedBefore, f.hapEnd + f.trimmedAfter],
      fragment: f,
    })
  }
  for (const d of departures) {
    add(d.name, { hap: [d.hapStart, d.hapStart + d.ownBp], departure: d })
  }
  const byHaplotype = new Map<string, HaplotypeDepartures>()
  for (const [name, items] of byPath) {
    const haplotype = panSNHaplotypeOrSample(name)
    if (haplotype === reference || name === refName) {
      continue
    }
    let h = byHaplotype.get(haplotype)
    if (!h) {
      h = { haplotype, covered: [], departures: [] }
      byHaplotype.set(haplotype, h)
    }
    for (const chain of pathChains(items, length)) {
      readPath(chain, h)
    }
  }
  const out = [...byHaplotype.values()]
  for (const h of out) {
    h.departures.sort((a, b) => a.refStart - b.refStart || a.refEnd - b.refEnd)
    h.covered = coveredOnce(h.covered)
  }
  return out
}

/**
 * A path's chains, each its events in the path's order: a fragment with the
 * departures leaving its last step, rejoining its first or lying inside it,
 * fragments that abut, and a departure touching no fragment on its own. What
 * the read did not see lies between chains, so each reads on its own, and R
 * that two chains reach has no call.
 */
function pathChains(items: Item[], length: (id: number) => number) {
  // items chain where one's contig stretch touches another's: a departure's
  // leave at a fragment's last step or inside it, its rejoin at the first,
  // or two fragments abutting. In start order, an item touches the chain so
  // far when it starts no later than the chain's furthest end. A fragment's
  // stretch holds the steps the cut trimmed, so a record's rejoin chains to
  // a piece whose first steps lie outside the drawn stretch
  const all = items
    .map(item => ({ item, edges: item.hap }))
    .sort((a, b) => a.edges[0] - b.edges[0])
  const chains: Item[][] = []
  let chainEnd = -Infinity
  for (const {
    item,
    edges: [start, end],
  } of all) {
    const chain = chains.at(-1)
    if (chain && start <= chainEnd) {
      chain.push(item)
    } else {
      chains.push([item])
    }
    chainEnd = Math.max(chainEnd, end)
  }
  return chains.map(chain => chainEvents(chain, length))
}

const on = (from: number, to: number, reversed: boolean): Event => ({
  kind: 'on',
  lo: Math.min(from, to),
  hi: Math.max(from, to),
  reversed,
})

function chainEvents(items: Item[], length: (id: number) => number) {
  const anchors: Anchor[] = []
  for (const item of items) {
    if ('fragment' in item) {
      const f = item.fragment
      const n = f.ids.length
      const sa = f.refStarts[0]!
      const sb = f.refStarts[n - 1]!
      const ra = f.against[0] === 1
      const rb = f.against[n - 1] === 1
      if (onR(sa)) {
        const pos = ra ? sa + length(f.ids[0]!) : sa
        anchors.push({ hap: f.hapStart, order: 1, pos, reversed: ra })
      }
      if (onR(sb)) {
        const pos = rb ? sb : sb + length(f.ids[n - 1]!)
        anchors.push({ hap: f.hapEnd, order: 0, pos, reversed: rb, exit: true })
      }
    } else {
      const d = item.departure
      anchors.push(
        {
          hap: d.hapStart,
          order: 0,
          pos: d.leavePos,
          reversed: d.leaveReversed,
          own: d.ownBp,
        },
        {
          hap: d.hapStart + d.ownBp,
          order: 1,
          pos: d.rejoinPos,
          reversed: d.rejoinReversed,
        },
      )
    }
  }
  anchors.sort((a, b) => a.hap - b.hap || a.order - b.order)
  const events: Event[] = []
  // the open stretch on R: where it began, and the furthest fragment end seen
  // along it
  let cur: { pos: number; reversed: boolean } | undefined
  let end: number | undefined
  const close = () => {
    if (cur) {
      events.push(on(cur.pos, end ?? cur.pos, cur.reversed))
    }
    cur = undefined
    end = undefined
  }
  for (const a of anchors) {
    if (cur && cur.reversed !== a.reversed) {
      close()
    }
    if (a.order === 1) {
      cur ??= { pos: a.pos, reversed: a.reversed }
    } else if (a.exit) {
      if (cur) {
        end = a.pos
      }
    } else {
      cur ??= { pos: a.pos, reversed: a.reversed }
      end = a.pos
      close()
      events.push({ kind: 'own', bp: a.own ?? 0 })
    }
  }
  close()
  return events
}

/**
 * A path's departures in the forward view: its stretches along R in order,
 * reversed where it walks R mostly against its orientation, each stretch
 * against it or of its own between two along it making one departure.
 */
function readPath(events: Event[], h: HaplotypeDepartures) {
  let bias = 0
  for (const e of events) {
    if (e.kind !== 'own') {
      bias += (e.reversed ? 1 : -1) * (e.hi - e.lo || 1)
    }
  }
  const flip = bias > 0
  let first = NaN
  let prevEnd = NaN
  let reach = -Infinity
  let alt = 0
  let inverted = 0
  for (const e of flip ? events.toReversed() : events) {
    if (e.kind === 'own') {
      alt += e.bp
    } else if (e.reversed !== flip) {
      alt += e.hi - e.lo
      inverted += e.hi - e.lo
    } else {
      const start = e.lo
      if (Number.isNaN(prevEnd)) {
        first = start
      } else if (start !== prevEnd || alt > 0) {
        const repeated = start < prevEnd ? prevEnd - start : 0
        h.departures.push({
          refStart: prevEnd,
          refEnd: repeated > 0 ? prevEnd : start,
          altLength: alt + repeated,
          invertedBp: inverted,
          repeatedBp: repeated,
          ...(repeated > 0 ? { copied: [start, prevEnd] } : {}),
        })
      }
      prevEnd = e.hi
      reach = Math.max(reach, e.hi)
      alt = 0
      inverted = 0
    }
  }
  if (!Number.isNaN(first)) {
    h.covered.push([first, reach])
  }
}

/**
 * The stretches of R that exactly one of a haplotype's paths covers, each
 * within one path. Where two paths overlap, such as a second contig over the
 * same locus, the haplotype holds two copies of that stretch and no call;
 * where one path ends and the next starts, no call crosses the junction.
 */
function coveredOnce(spans: [number, number][]) {
  const points = [...new Set(spans.flat())].sort((a, b) => a - b)
  const out: [number, number][] = []
  let lastPath = -1
  for (let k = 0; k + 1 < points.length; k++) {
    const a = points[k]!
    const b = points[k + 1]!
    const paths = spans.flatMap(([s, e], i) => (s <= a && e >= b ? [i] : []))
    const last = out.at(-1)
    if (paths.length !== 1) {
      lastPath = -1
    } else if (last && paths[0] === lastPath && last[1] === a) {
      last[1] = b
    } else {
      out.push([a, b])
      lastPath = paths[0]!
    }
  }
  return out
}

// The most steps a forms cut keeps
export const FORMS_STEP_BUDGET = 16_000_000

export interface FormsCutOptions {
  // the most steps the cut keeps, FORMS_STEP_BUDGET at least
  stepBudget?: number
  // the most compressed bytes the cut fetches from the region's tiles and
  // their blocks units; WALK_BUDGETS' by default
  byteBudget?: number
  // the most compressed bytes the cut's departure read fetches;
  // WALK_BUDGETS' by default
  departureByteBudget?: number
  signal?: AbortSignal
}

// Every haplotype's departures over a region, for structuralForms, from one
// cut without links that expands no walk. A haplotype the header names with
// no walk in the cut is listed with nothing covered.
export async function formsCut(
  reader: Pick<WalkReader, 'header' | 'sequenceFor' | 'cut'>,
  region: { refName: string; start: number; end: number },
  opts: FormsCutOptions = {},
) {
  const {
    stepBudget = FORMS_STEP_BUDGET,
    byteBudget,
    departureByteBudget,
    signal,
  } = opts
  const header: WalkReaderHeader = await reader.header({ signal })
  const sample = header.references[0] ?? ''
  const refName =
    (await reader.sequenceFor(sample, region.refName, { signal })) ??
    region.refName
  const cut = await reader.cut(refName, region, {
    links: 'none',
    stepBudget: Math.max(stepBudget, FORMS_STEP_BUDGET),
    byteBudget,
    departureByteBudget,
    signal,
  })
  const haplotypes = walkDepartures(
    {
      fragments: cut.fragments,
      departures: cut.departures,
      nodes: cut.graph.nodes,
    },
    refName,
  )
  const seen = new Set(haplotypes.map(h => h.haplotype))
  for (const haplotype of header.haplotypes) {
    if (!seen.has(haplotype)) {
      haplotypes.push({ haplotype, covered: [], departures: [] })
    }
  }
  return haplotypes
}
