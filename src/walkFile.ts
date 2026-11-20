import { unzip } from '@gmod/bgzf-filehandle'

import { WalkRowError } from './walkChunk.ts'

import type { GenericFilehandle } from 'generic-filehandle2'

/**
 * Walk format 8 (gfa-to-tabix's ADRs 0010 and 0012): `<prefix>.<sample>.walks.gz`, a
 * series of gzip units, and its directory `<prefix>.<sample>.walks.gz.idx`.
 * The directory's first line is the data's `#walks` line plus
 * `bytes:i:<data length>`, then `#header 0 <length>` and one
 * `seq <offset> <length>` line per reference sequence, giving the unit that
 * lists the sequence's units: `#entry seq start end kind name offset length`.
 */
export const WALK_FORMAT = 8

// What the reader needs of a file: generic-filehandle2's handles are ones,
// and `stat` gives the data's length, from Content-Range on a RemoteFile's
// first read. Where it gives 0, as a RemoteFile does when CORS hides
// Content-Range, the reader skips checking the directory's `bytes`
export type WalkFilehandle = Pick<
  GenericFilehandle,
  'read' | 'readFile' | 'stat' | 'source'
>

// A walk file the reader refuses: one in a format other than 8, or one whose
// units do not hold what its directory and tables say
export class WalkFormatError extends Error {
  override name = 'WalkFormatError'
}

/**
 * One unit of a sequence's table: `start end` is a tile's true extent, a
 * departure block's span, a detail record's pieces' span on the reference,
 * or the span of a tile's abridged pieces' blocks; `name` a tile's nominal
 * start, a detail's `tile:detail`, a `blocks` unit's tile, or `.` for a
 * departure block; `offset length` its compressed bytes in the data file
 */
export interface WalkUnit {
  refName: string
  start: number
  end: number
  kind: string
  name: string
  offset: number
  length: number
}

interface Directory {
  walks: string[]
  bytes: number | undefined
  header: { offset: number; length: number }
  tables: Map<string, { offset: number; length: number }>
}

const decoder = new TextDecoder()

// The row kinds each unit kind holds; a unit kind not here is a later
// addition, whose rows a reader never reads
const ROW_KINDS = new Map([
  ['tile', ['ref', 'nodes', 'route', 'walk', 'abridged']],
  ['departures', ['departure']],
  ['detail', ['nodes', 'route', 'walk']],
  ['blocks', ['blocks']],
])

async function unitLines(bytes: Uint8Array, unit: string) {
  const text = await unzip(bytes).then(
    unzipped => decoder.decode(unzipped),
    (error: unknown) => {
      throw new WalkFormatError(`${unit} is not whole gzip members`, {
        cause: error,
      })
    },
  )
  if (!text.endsWith('\n')) {
    throw new WalkFormatError(
      `${unit} ends partway through a line: ${JSON.stringify(text.slice(-40))}`,
    )
  }
  return text.split('\n').filter(line => line !== '')
}

function unitName(unit: WalkUnit) {
  return `The unit at byte ${unit.offset} of the walk file, ${unit.refName}'s ${unit.kind} ${unit.name} over [${unit.start}, ${unit.end}) as its table lists it,`
}

function isTableHead(head: string, refName: string) {
  const [tag, seq, start, , kind, name, ...more] = head.split('\t')
  return (
    more.length === 0 &&
    [tag, seq, start, kind, name].join('\t') ===
      ['#unit', refName, '0', 'table', '.'].join('\t')
  )
}

// The columns each row kind holds; any after them are tags, `NAME:T:value`
const FIXED_COLUMNS = new Map([
  ['ref', 12],
  ['nodes', 6],
  ['route', 7],
  ['walk', 11],
  ['abridged', 11],
  ['blocks', 8],
  ['departure', 15],
])

const TAG = /^[A-Za-z][A-Za-z0-9]*:[A-Za-z]:/

// A column after a row's fixed ones that is no tag
function untagged(row: string) {
  const fixed = FIXED_COLUMNS.get(row.split('\t', 4)[3] ?? '') ?? Infinity
  let at = -1
  let found = 0
  while (found < fixed) {
    at = row.indexOf('\t', at + 1)
    found = at < 0 ? Infinity : found + 1
  }
  return found === fixed
    ? row
        .slice(at + 1)
        .split('\t')
        .find(col => !TAG.test(col))
    : undefined
}

// The fixed columns a row's kind holds, where the row holds fewer
function lacking(row: string) {
  const fixed = FIXED_COLUMNS.get(row.split('\t', 4)[3] ?? '')
  return fixed !== undefined && row.split('\t', fixed).length < fixed
    ? fixed
    : undefined
}

function checkUnit(unit: WalkUnit, head: string | undefined, rows: string[]) {
  const stated = [
    '#unit',
    unit.refName,
    unit.start,
    unit.end,
    unit.kind,
    unit.name,
  ].join('\t')
  const kinds = ROW_KINDS.get(unit.kind)
  let stray: string | undefined
  let extra: string | undefined
  let short: string | undefined
  if (kinds !== undefined) {
    for (const row of rows) {
      const cols = row.split('\t')
      const kind = cols[3] ?? ''
      const fixed = FIXED_COLUMNS.get(kind)
      if (
        stray === undefined &&
        (cols[0] !== unit.refName || !kinds.includes(kind))
      ) {
        stray = row
      }
      if (fixed !== undefined && cols.length !== fixed) {
        if (cols.length < fixed) {
          short ??= row
        } else if (extra === undefined) {
          for (let k = fixed; k < cols.length; k++) {
            if (!TAG.test(cols[k]!)) {
              extra = row
            }
          }
        }
      }
    }
  }
  if (head !== stated) {
    throw new WalkFormatError(
      `${unitName(unit)} starts ${JSON.stringify(head?.slice(0, 80) ?? '')}`,
    )
  } else if (stray !== undefined) {
    throw new WalkFormatError(
      `${unitName(unit)} holds a row that is not one of ${unit.refName}'s ${kinds?.join(', ')} rows: ${JSON.stringify(stray.slice(0, 80))}`,
    )
  } else if (extra !== undefined) {
    throw new WalkFormatError(
      `${unitName(unit)} holds a row with a column after its fixed ones that is no tag, ${JSON.stringify(untagged(extra)?.slice(0, 40))}: ${JSON.stringify(extra.slice(0, 80))}`,
    )
  } else if (short !== undefined) {
    throw new WalkRowError(
      `${unitName(unit)} holds a row of fewer columns than the ${lacking(short)} its kind holds: ${JSON.stringify(short.slice(0, 80))}`,
    )
  }
}

function integerTag(fields: string[], tag: string) {
  const field = fields.find(f => f.startsWith(`${tag}:i:`))
  return field === undefined ? undefined : +field.slice(tag.length + 3)
}

function formatError(format: number) {
  return new WalkFormatError(
    format > WALK_FORMAT
      ? `This walk file is format ${format}; this version of @jbrowse/graph-walks reads format ${WALK_FORMAT}, so update the viewer, or @jbrowse/graph-walks, to one that reads format ${format}`
      : `This walk file is format ${format}; @jbrowse/graph-walks reads format ${WALK_FORMAT}, so rebuild it with gfa-to-tabix 5.0 or later`,
  )
}

function parseDirectory(text: string): Directory {
  const [first = '', ...rest] = text.split('\n').filter(line => line !== '')
  const fields = first.split('\t')
  if (fields[0] !== '#walks') {
    throw new WalkFormatError(
      `This is not a walk file's directory: its first line starts ${JSON.stringify(first.slice(0, 40))}, not #walks`,
    )
  }
  const format = integerTag(fields, 'format') ?? 1
  if (format !== WALK_FORMAT) {
    throw formatError(format)
  }
  let header: Directory['header'] | undefined
  const tables = new Map<string, { offset: number; length: number }>()
  // a line whose first word starts with # and is not #header is a later
  // addition
  for (const line of rest) {
    const [name = '', offset, length, ...more] = line.split('\t')
    const range = { offset: Number(offset), length: Number(length) }
    const known = name === '#header' || !name.startsWith('#')
    if (
      known &&
      (more.length > 0 ||
        !Number.isSafeInteger(range.offset) ||
        !Number.isSafeInteger(range.length))
    ) {
      throw new WalkFormatError(
        `A walk file directory line holds no byte range, or more than one: ${line}`,
      )
    }
    if (name === '#header') {
      header = range
    } else if (known) {
      tables.set(name, range)
    }
  }
  if (!header) {
    throw new WalkFormatError('The walk file directory has no #header line')
  }
  return {
    walks: fields.filter(f => !f.startsWith('bytes:i:')),
    bytes: integerTag(fields, 'bytes'),
    header,
    tables,
  }
}

// A read shared by every caller, tried again by the next caller once it fails
export function once<T>(read: () => Promise<T>) {
  let held: Promise<T> | undefined
  return () => {
    if (!held) {
      const started = read()
      started.catch(() => {
        if (held === started) {
          held = undefined
        }
      })
      held = started
    }
    return held
  }
}

// Settles like `promise` unless `signal` aborts first
export function raceAbort<T>(promise: Promise<T>, signal?: AbortSignal) {
  return signal
    ? new Promise<T>((resolve, reject) => {
        signal.throwIfAborted()
        const onAbort = () => {
          reject(
            signal.reason instanceof Error
              ? signal.reason
              : new Error(String(signal.reason)),
          )
        }
        signal.addEventListener('abort', onAbort, { once: true })
        promise.then(resolve, reject).finally(() => {
          signal.removeEventListener('abort', onAbort)
        })
      })
    : promise
}

/**
 * A format 8 walk file read by byte ranges: the directory and the header
 * once, each sequence's table once, and units as asked, the ranges of units
 * that touch merged into one single-range read
 */
export class WalkFile {
  private readonly data: WalkFilehandle
  private readonly directory: WalkFilehandle
  private readonly tableReads = new Map<string, () => Promise<WalkUnit[]>>()

  constructor(data: WalkFilehandle, directory: WalkFilehandle) {
    this.data = data
    this.directory = directory
  }

  private readonly readDirectory = once(async () =>
    parseDirectory(await this.directory.readFile('utf8')),
  )

  // The header unit's lines, once the directory and the data agree
  readonly header = once(async () => {
    const directory = await this.readDirectory()
    const { offset, length } = directory.header
    const header = await unitLines(
      await this.data.read(length, offset),
      'The header unit',
    )
    const { size } = await this.data.stat()
    if (header[0] !== directory.walks.join('\t')) {
      throw new WalkFormatError(
        `The walk file's directory does not match its data: the directory states ${directory.walks.join(' ')}, the data ${header[0] ?? 'no #walks line'}`,
      )
    } else if (size > 0 && directory.bytes !== size) {
      throw new WalkFormatError(
        `The walk file's directory does not match its data: the directory states ${directory.bytes ?? 'no'} bytes, the data holds ${size}`,
      )
    }
    return header
  })

  // The reference sequences, in the directory's order, once the directory
  // and the data agree
  async sequences() {
    await this.header()
    return [...(await this.readDirectory()).tables.keys()]
  }

  // A sequence's units in file order, once the directory and the data
  // agree; none for a sequence the file lacks
  async units(refName: string) {
    let read = this.tableReads.get(refName)
    if (!read) {
      read = once(async () => {
        const range = (await this.readDirectory()).tables.get(refName)
        const [checked, table] = await Promise.allSettled([
          this.header(),
          range ? this.table(refName, range) : [],
        ])
        if (checked.status === 'rejected') {
          throw checked.reason
        } else if (table.status === 'rejected') {
          throw table.reason
        }
        return table.value
      })
      this.tableReads.set(refName, read)
    }
    return read()
  }

  private async table(
    refName: string,
    { offset, length }: { offset: number; length: number },
  ) {
    const [head = '', ...entries] = await unitLines(
      await this.data.read(length, offset),
      `${refName}'s table`,
    )
    if (!isTableHead(head, refName)) {
      throw new WalkFormatError(
        `${refName}'s table, at byte ${offset} of the walk file, starts ${JSON.stringify(head.slice(0, 80))}`,
      )
    }
    return entries.map(line => {
      const [
        tag,
        seq = '',
        start,
        end,
        kind = '',
        name = '',
        at,
        size,
        ...more
      ] = line.split('\t')
      const unit = {
        refName: seq,
        start: Number(start),
        end: Number(end),
        kind,
        name,
        offset: Number(at),
        length: Number(size),
      }
      if (
        tag !== '#entry' ||
        more.length > 0 ||
        seq !== refName ||
        ![unit.start, unit.end, unit.offset, unit.length].every(
          Number.isSafeInteger,
        )
      ) {
        throw new WalkFormatError(`${refName}'s unit table holds ${line}`)
      }
      return unit
    })
  }

  /**
   * The lines of each unit, without its `#unit` line, in the order given.
   * Units whose byte ranges touch are fetched by one read, and each is
   * gunzipped on its own.
   */
  async read(units: WalkUnit[], signal?: AbortSignal) {
    const byOffset = [...new Map(units.map(u => [u.offset, u])).values()]
      .filter(u => u.length > 0)
      .sort((a, b) => a.offset - b.offset)
    const runs: WalkUnit[][] = []
    for (const unit of byOffset) {
      const run = runs.at(-1)
      const last = run?.at(-1)
      if (run && last && last.offset + last.length === unit.offset) {
        run.push(unit)
      } else {
        runs.push([unit])
      }
    }
    const held = new Map<number, string[]>()
    await Promise.all(
      runs.map(async run => {
        const first = run[0]!
        const last = run.at(-1)!
        const length = last.offset + last.length - first.offset
        const bytes = await this.data.read(length, first.offset, { signal })
        for (const unit of run) {
          const at = unit.offset - first.offset
          const [head, ...rest] = await unitLines(
            bytes.subarray(at, at + unit.length),
            unitName(unit),
          )
          checkUnit(unit, head, rest)
          held.set(unit.offset, rest)
        }
      }),
    )
    return units.map(unit => held.get(unit.offset) ?? [])
  }
}
