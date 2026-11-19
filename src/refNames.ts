import { panSNContig, panSNHaplotype, panSNSample } from './util.ts'

// A tab can occur in neither half: both come out of BED columns
function qualifiedKey(assemblyName: string, refName: string) {
  return `${assemblyName}\t${refName}`
}

export class AmbiguousSequenceError extends Error {
  override name = 'AmbiguousSequenceError'
  readonly sequences: string[]

  constructor(assemblyName: string, refName: string, sequences: string[]) {
    super(
      `${refName} on ${assemblyName} is ${sequences.length} sequences of the file, ${sequences.join(', ')}; name one of them, or its haplotype: ${sequences.map(name => panSNHaplotype(name)).join(', ')}`,
    )
    this.sequences = sequences
  }
}

/**
 * A file's sequence names under each spelling a caller may ask by: the file's
 * own (`GRCh38#0#chr6`), and a PanSN name's contig under its sample (`chr6`
 * on `GRCh38`) and under its haplotype (`chr6` on `GRCh38#0`). A pggb name
 * with no haplotype field, `chm13#chr1`, is `chr1` on `chm13`. A spelling
 * two names share, `chr1` on `HG002` where the file holds `HG002#1#chr1` and
 * `HG002#2#chr1`, keeps both.
 */
export function refNameLookup(names: string[]) {
  const lookup = new Map<string, string[]>()
  const file = (key: string, name: string) => {
    const held = lookup.get(key)
    if (held === undefined) {
      lookup.set(key, [name])
    } else if (!held.includes(name)) {
      held.push(name)
    }
  }
  for (const name of names) {
    file(name, name)
    const sample = panSNSample(name)
    if (sample !== name) {
      const contig = panSNContig(name)
      file(qualifiedKey(sample, contig), name)
      const haplotype = panSNHaplotype(name)
      if (haplotype !== undefined) {
        file(qualifiedKey(haplotype, contig), name)
      }
    }
  }
  return lookup
}

/**
 * The file's name for `refName` on `assemblyName`, a sample or a haplotype,
 * else for `refName` as the file spells it; undefined when it holds neither.
 * Throws AmbiguousSequenceError where the spelling names more than one.
 */
export function resolveRefName(
  lookup: Map<string, string[]>,
  assemblyName: string,
  refName: string,
) {
  const held =
    lookup.get(qualifiedKey(assemblyName, refName)) ?? lookup.get(refName)
  if (held !== undefined && held.length > 1) {
    throw new AmbiguousSequenceError(assemblyName, refName, held)
  }
  return held?.[0]
}
