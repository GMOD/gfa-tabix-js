import { createHash } from 'node:crypto'
import { crc32, deflateRawSync } from 'node:zlib'

import { BlobFile } from 'generic-filehandle2'

import { isPlaceholder, runName } from './walkRows.ts'

import type { GraphTables } from './util.ts'
import type { WalkFilehandle, WalkUnit } from './walkFile.ts'
import type { WalkReader } from './walkReader.ts'
import type { FilehandleOptions } from 'generic-filehandle2'
import type { Hash } from 'node:crypto'

type Cut = Awaited<ReturnType<WalkReader['cut']>>

// Every read `data` makes from here on, as its byte range
export function readsOf(data: WalkFilehandle) {
  const read = data.read.bind(data)
  const reads: { position: number; length: number }[] = []
  data.read = (length, position, opts) => {
    reads.push({ position, length })
    return read(length, position, opts)
  }
  return reads
}

// The units of a sequence's table that `reads` fetched, each as often as read
export const unitsRead = (
  units: WalkUnit[],
  reads: { position: number; length: number }[],
) =>
  reads.flatMap(({ position, length }) =>
    units.filter(u => u.offset >= position && u.offset < position + length),
  )

export const bytesRead = (reads: { length: number }[]) =>
  reads.reduce((sum, read) => sum + read.length, 0)

// The byte ranges `reads` cover, merged where they touch or overlap
export function rangesRead(reads: { position: number; length: number }[]) {
  const ranges: { start: number; end: number }[] = []
  for (const { position, length } of [...reads].sort(
    (a, b) => a.position - b.position,
  )) {
    const last = ranges.at(-1)
    if (last && position <= last.end) {
      last.end = Math.max(last.end, position + length)
    } else {
      ranges.push({ start: position, end: position + length })
    }
  }
  return ranges
}

// Each placeholder of a cut as its id and its run's name: a node flagged one
// by the stretch it states, and a placeholder step of a fragment by where the
// fragment walks it
export function placeholderRuns(cut: Cut) {
  const runs: [number, string][] = []
  for (const [id, node] of cut.graph.nodes) {
    if (node.placeholder) {
      runs.push([id, runName(node.refName, node.start, node.end)])
    }
  }
  for (const f of cut.fragments) {
    let at = f.hapStart
    for (const id of f.ids) {
      const node = cut.graph.nodes.get(id)
      const end = at + (node ? node.end - node.start : Number.NaN)
      if (isPlaceholder(id) && !Number.isNaN(end)) {
        runs.push([id, runName(f.name, at, end)])
      }
      at = end
    }
  }
  return runs
}

// The placeholders whose id stands for another run too, or whose run has
// another id too
export function idsNotOneToOne(runs: [number, string][]) {
  const runOf = new Map<number, string>()
  const idOf = new Map<string, number>()
  const out = new Set<string>()
  for (const [id, run] of runs) {
    const heldRun = runOf.get(id) ?? run
    const heldId = idOf.get(run) ?? id
    if (heldRun !== run || heldId !== id) {
      out.add(`${id} ${run}, beside ${heldId} ${heldRun}`)
    }
    runOf.set(id, heldRun)
    idOf.set(run, heldId)
  }
  return [...out]
}

// What a cut returns and the tables it makes
export function cutFacts(
  cut: Cut,
  tables = cut.graph.tables(cut.kept, cut.fragments, cut.departures, cut.gaps),
) {
  return {
    tables,
    kept: [...cut.kept],
    nodes: [...cut.graph.nodes],
    fragments: cut.fragments,
    departures: cut.departures,
    gaps: cut.gaps,
    placeholders: cut.placeholders,
  }
}

function hashInto(hash: Hash, value: unknown) {
  if (ArrayBuffer.isView(value)) {
    hash.update(`${value.constructor.name} ${value.byteLength}:`)
    hash.update(
      new Uint8Array(value.buffer, value.byteOffset, value.byteLength),
    )
  } else if (Array.isArray(value)) {
    hash.update(`[${value.length}:`)
    for (const item of value) {
      hashInto(hash, item)
    }
  } else if (typeof value === 'object' && value !== null) {
    for (const [key, item] of Object.entries(value)) {
      hash.update(`{${key}:`)
      hashInto(hash, item)
    }
    hash.update('}')
  } else {
    hash.update(`${typeof value} ${String(value)};`)
  }
}

// The same for two values only where they are the same, byte for byte
export function digestOf(value: unknown) {
  const hash = createHash('sha256')
  hashInto(hash, value)
  return hash.digest('hex')
}

// Each walk of a cut's tables, its name and then each step as its node's name
// and where the step lies on the reference
export const walkPlaces = ({ nodes, walks }: GraphTables) =>
  walks.names.map((name, p) => {
    const from = walks.offsets[p]!
    const steps = Array.from(
      walks.steps.subarray(from, walks.offsets[p + 1]),
      (node, k) => `${nodes.names[node]}@${walks.refStarts![from + k]}`,
    )
    return `${name} ${steps.join(' ')}`
  })

// The GFA a cut's tables stand for, as @jbrowse/bandage-core writes it
export function graphTablesGFA({ nodes, links, walks }: GraphTables) {
  const { names } = nodes
  const lines = ['H\tVN:Z:1.1']
  for (let i = 0; i < nodes.lengths.length; i++) {
    const ref = nodes.refs[i]!
    const stable =
      ref < 0
        ? ''
        : `\tSN:Z:${nodes.refNames[ref]}\tSO:i:${nodes.starts[i]}\tSR:i:${nodes.ranks[i]}`
    lines.push(`S\t${names[i]}\t*\tLN:i:${nodes.lengths[i]}${stable}`)
  }
  for (let i = 0; i < links.from.length; i++) {
    const strands = links.strands[i]!
    lines.push(
      `L\t${names[links.from[i]!]}\t${strands & 1 ? '-' : '+'}\t${names[links.to[i]!]}\t${strands & 2 ? '-' : '+'}\t0M`,
    )
  }
  for (let p = 0; p < walks.names.length; p++) {
    const [sample, haplotype, ...contig] = walks.names[p]!.split('#')
    let body = ''
    for (let s = walks.offsets[p]!; s < walks.offsets[p + 1]!; s++) {
      body += (walks.reversed[s]! !== 0 ? '<' : '>') + names[walks.steps[s]!]
    }
    lines.push(
      `W\t${sample}\t${haplotype}\t${contig.join('#')}\t${walks.starts[p]}\t${walks.ends[p]}\t${body}`,
    )
  }
  return lines.join('\n')
}

const BGZF_BLOCK = 0xff00
const hex = (text: string) =>
  Uint8Array.from(text.match(/../g) ?? [], byte => parseInt(byte, 16))
const BGZF_HEADER = hex('1f8b08040000000000ff060042430200')
const BGZF_EOF = hex('1f8b08040000000000ff0600424302001b0003000000000000000000')

// `text` as BGZF blocks, as gfa-to-tabix writes a unit
export function bgzf(text: string) {
  const data = new TextEncoder().encode(text)
  const blocks: Uint8Array<ArrayBuffer>[] = []
  for (let at = 0; at < data.length; at += BGZF_BLOCK) {
    const piece = data.subarray(at, at + BGZF_BLOCK)
    const deflated = deflateRawSync(piece)
    const block = new Uint8Array(26 + deflated.length)
    block.set(BGZF_HEADER)
    block.set(deflated, 18)
    const view = new DataView(block.buffer)
    view.setUint16(16, block.length - 1, true)
    view.setUint32(18 + deflated.length, crc32(piece), true)
    view.setUint32(22 + deflated.length, piece.length, true)
    blocks.push(block)
  }
  return blocks
}

const lengthOf = (blocks: Uint8Array[]) =>
  blocks.reduce((sum, block) => sum + block.length, 0)

// A file whose bytes past `from` lie `gap` further on than the blob holds them
class GappedFile extends BlobFile {
  private readonly from: number
  private readonly gap: number

  constructor(blob: Blob, from: number, gap: number) {
    super(blob)
    this.from = from
    this.gap = gap
  }

  override read(length: number, position: number, opts?: FilehandleOptions) {
    return super.read(
      length,
      position < this.from ? position : position - this.gap,
      opts,
    )
  }

  override async stat() {
    const { size } = await super.stat()
    return { size: size + this.gap }
  }
}

// A detail or blocks unit's key, `kind name`, as its kind, its name, its
// tile's start and its detail number, a blocks unit after its tile's details
function movedUnit(key: string) {
  const [kind = '', name = ''] = key.split(' ')
  const [tileStart, detail = Infinity] = name.split(':').map(Number)
  return { kind, name, tileStart: tileStart ?? 0, detail }
}

const spanOf = (cols: string[]) =>
  (cols[9] ?? '').split(',').reduce((end, n) => end + Number(n), +cols[5]!)

/**
 * A format 8 walk file in memory from `lines`, each space-separated: the
 * header's `#` lines; tile rows as `seq tile kind ...`, a tile's rows from
 * `kind` on under its nominal start; and departure lines as the file holds
 * them, `seq start end departure ...`, each filed in every block of the
 * tile size it overlaps; a detail unit's rows as `seq tile:detail start end
 * kind ...`; and a blocks unit's as `seq tile start end blocks ...`. A
 * tile's extent is its ref rows' span, and each of its lines states that
 * extent, but an abridged row its blocks row's span; every inline tile has
 * ref rows. A detail or blocks unit's extent
 * is its rows' span, and its units follow the sequence's tiles and
 * departure blocks. `reversed` writes a sequence's tiles and departure
 * blocks last first, and `gap` puts that many bytes the file never reads
 * between the header and the first unit.
 */
export function walkFileOf(
  lines: string[],
  { reversed = false, gap = 0 } = {},
) {
  const rows = lines.map(line => line.split(' '))
  const header = rows.filter(r => r[0]!.startsWith('#')).map(r => r.join('\t'))
  const tile = Number(/tile:i:(\d+)/.exec(header[0] ?? '')?.[1] ?? 65536)
  const body = rows.filter(r => !r[0]!.startsWith('#'))
  const sequences = [...new Set(body.map(r => r[0]!))]
  const parts = bgzf(`${header.join('\n')}\n`)
  const headerLength = lengthOf(parts)
  let offset = headerLength + gap
  const directory = [`#header\t0\t${headerLength}`]
  const units: WalkUnit[] = []
  const write = (text: string) => {
    const blocks = bgzf(text)
    const length = lengthOf(blocks)
    for (const block of blocks) {
      parts.push(block)
    }
    offset += length
    return { offset: offset - length, length }
  }
  for (const seq of sequences) {
    const items: { start: number; order: number; unit: () => WalkUnit }[] = []
    const tiles = new Map<number, string[][]>()
    const blocks = new Map<number, string[]>()
    // the rows of each detail and blocks unit, by kind and name
    const moved = new Map<string, string[][]>()
    for (const r of body.filter(row => row[0] === seq)) {
      if (r[1]!.includes(':') || r[4] === 'blocks') {
        const key = `${r[4] === 'blocks' ? 'blocks' : 'detail'} ${r[1]}`
        moved.set(key, [...(moved.get(key) ?? []), [seq, ...r.slice(2)]])
      } else if (Number.isNaN(Number(r[2]))) {
        tiles.set(+r[1]!, [...(tiles.get(+r[1]!) ?? []), r])
      } else {
        for (
          let at = Math.floor(+r[1]! / tile) * tile;
          at < +r[2]!;
          at += tile
        ) {
          blocks.set(at, [...(blocks.get(at) ?? []), r.join('\t')])
        }
      }
    }
    // an abridged row states its piece's span, which its blocks row gives
    const spans = new Map(
      body
        .filter(r => r[0] === seq && r[4] === 'blocks')
        .map(r => [[r[5], r[6], r[7]].join(' '), [r[2]!, r[3]!]]),
    )
    for (const [start, tileRows] of tiles) {
      const refs = tileRows.filter(r => r[2] === 'ref')
      const lo = refs.reduce((at, r) => Math.min(at, +r[5]!), Infinity)
      const hi = refs.reduce((at, r) => Math.max(at, spanOf(r)), -Infinity)
      items.push({
        start,
        order: 0,
        unit: () => {
          const text = [
            ['#unit', seq, lo, hi, 'tile', start],
            ...tileRows.map(r => [
              seq,
              ...(r[2] === 'abridged'
                ? (spans.get([r[3], r[4], r[6]].join(' ')) ?? [lo, hi])
                : [lo, hi]),
              ...r.slice(2),
            ]),
          ]
          const at = write(`${text.map(c => c.join('\t')).join('\n')}\n`)
          const unit = { refName: seq, start: lo, end: hi, kind: 'tile' }
          return { ...unit, name: String(start), ...at }
        },
      })
    }
    for (const [start, blockLines] of blocks) {
      items.push({
        start,
        order: 1,
        unit: () => {
          const head = ['#unit', seq, start, start + tile, 'departures', '.']
          const at = write(`${[head.join('\t'), ...blockLines].join('\n')}\n`)
          const unit = { refName: seq, start, end: start + tile }
          return { ...unit, kind: 'departures', name: '.', ...at }
        },
      })
    }
    items.sort((a, b) =>
      reversed
        ? b.start - a.start || b.order - a.order
        : a.start - b.start || a.order - b.order,
    )
    const held = items.map(item => item.unit())
    for (const [key, unitRows] of [...moved].sort(
      ([a], [b]) =>
        movedUnit(a).tileStart - movedUnit(b).tileStart ||
        movedUnit(a).detail - movedUnit(b).detail,
    )) {
      const { kind, name } = movedUnit(key)
      const lo = unitRows.reduce((at, r) => Math.min(at, +r[1]!), Infinity)
      const hi = unitRows.reduce((at, r) => Math.max(at, +r[2]!), -Infinity)
      const text = [['#unit', seq, lo, hi, kind, name], ...unitRows]
      const at = write(`${text.map(c => c.join('\t')).join('\n')}\n`)
      held.push({ refName: seq, start: lo, end: hi, kind, name, ...at })
    }
    for (const unit of held) {
      units.push(unit)
    }
    const end = held.reduce((at, u) => Math.max(at, u.end), 0)
    const table = write(
      `${[
        ['#unit', seq, 0, end, 'table', '.'],
        ...(gap > 0 && seq === sequences[0]
          ? [['#entry', seq, 0, 0, 'filler', '.', headerLength, gap]]
          : []),
        ...held.map(u => [
          '#entry',
          seq,
          u.start,
          u.end,
          u.kind,
          u.name,
          u.offset,
          u.length,
        ]),
      ]
        .map(c => c.join('\t'))
        .join('\n')}\n`,
    )
    directory.push(`${seq}\t${table.offset}\t${table.length}`)
  }
  parts.push(BGZF_EOF)
  const index = `${header[0]}\tbytes:i:${lengthOf(parts) + gap}\n${directory.join('\n')}\n`
  return {
    data: new GappedFile(new Blob(parts), headerLength, gap),
    directory: new BlobFile(new Blob([index])),
    units,
  }
}
