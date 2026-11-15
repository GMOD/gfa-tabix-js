import type { GraphTables } from './util.ts'

export interface AbridgedWalk {
  name: string
  steps: number
  bp: number
}

export function abridgedWalks({ nodes, walks, departures }: GraphTables) {
  const byName = new Map<string, { steps: number; held: Set<number> }>()
  const add = (name: string, steps: Int32Array, abridged = Number.NaN) => {
    const walk = byName.get(name) ?? { steps: 0, held: new Set<number>() }
    walk.steps += Number.isNaN(abridged) ? 0 : abridged
    for (const i of steps) {
      if (nodes.placeholders?.[i] === 1) {
        walk.held.add(i)
      }
    }
    byName.set(name, walk)
  }
  walks.names.forEach((name, p) => {
    add(
      name,
      walks.steps.subarray(walks.offsets[p], walks.offsets[p + 1]),
      walks.abridgedSteps?.[p],
    )
  })
  departures.names.forEach((name, k) => {
    add(
      name,
      departures.steps.subarray(
        departures.offsets[k],
        departures.offsets[k + 1],
      ),
      departures.abridgedSteps?.[k],
    )
  })
  const found: AbridgedWalk[] = []
  for (const [name, { steps, held }] of byName) {
    let bp = 0
    for (const i of held) {
      bp += nodes.lengths[i]!
    }
    if (steps > 0 || bp > 0) {
      found.push({ name, steps, bp })
    }
  }
  return found
}

const groupDigits = (n: number) => String(n).replace(/\B(?=(\d{3})+$)/g, ',')

const readoutBp = (bp: number) =>
  bp < 1000 ? `${bp} bp` : `${(bp / 1000).toFixed(bp < 10_000 ? 1 : 0)} kb`

export function abridgedText({ name, steps, bp }: AbridgedWalk, label = name) {
  const notLoaded = `${label}: ${groupDigits(steps)} step${steps === 1 ? '' : 's'} not loaded`
  return bp > 0 ? `${notLoaded}, ${readoutBp(bp)} off the reference` : notLoaded
}
