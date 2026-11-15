import type { GraphTables } from './util.ts'

export interface AbridgedWalk {
  name: string
  bp: number
}

export function abridgedWalks({ nodes, walks, departures }: GraphTables) {
  const byName = new Map<string, Set<number>>()
  const add = (name: string, steps: Int32Array) => {
    const held = byName.get(name) ?? new Set<number>()
    for (const i of steps) {
      if (nodes.placeholders?.[i] === 1) {
        held.add(i)
      }
    }
    byName.set(name, held)
  }
  walks.names.forEach((name, p) => {
    add(name, walks.steps.subarray(walks.offsets[p], walks.offsets[p + 1]))
  })
  departures.names.forEach((name, k) => {
    add(
      name,
      departures.steps.subarray(
        departures.offsets[k],
        departures.offsets[k + 1],
      ),
    )
  })
  const found: AbridgedWalk[] = []
  for (const [name, held] of byName) {
    let bp = 0
    for (const i of held) {
      bp += nodes.lengths[i]!
    }
    if (bp > 0) {
      found.push({ name, bp })
    }
  }
  return found
}

const readoutBp = (bp: number) =>
  bp < 1000 ? `${bp} bp` : `${(bp / 1000).toFixed(bp < 10_000 ? 1 : 0)} kb`

export function abridgedText({ name, bp }: AbridgedWalk, label = name) {
  return `${label}: ${readoutBp(bp)} off the reference not loaded`
}
