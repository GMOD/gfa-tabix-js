# API reference

Open a `WalkFile` and call `getGraphForRange`. This page covers that path: the
file, the method, its options, the tables it returns and the errors it throws.
The other methods are listed under [Going further](#going-further).

Coordinates are 0-based half-open throughout: `chr6:160,614,798-160,647,758` is
`160_614_797, 160_647_758`.

## Opening a walk file

```ts
import { WalkFile } from '@gmod/gfa-tabix'

const walks = new WalkFile({
  walksUrl:
    'https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38-f8.GRCh38.walks.gz',
})
```

| Option            | Effect                                                                           |
| ----------------- | -------------------------------------------------------------------------------- |
| `walksUrl`        | The walk file by URL, read by HTTP range requests                                |
| `walksPath`       | The walk file by local path                                                      |
| `walksFilehandle` | The walk file already open, as a generic-filehandle2 handle                      |
| `idxUrl`          | The directory by URL; by default `walksUrl` plus `.idx`, before any query string |
| `idxPath`         | The directory by local path; by default `walksPath` plus `.idx`                  |
| `idxFilehandle`   | The directory already open, which a `walksFilehandle` needs                      |

- A `WalkFile` takes exactly one of the three `walks` options and at most one of
  the three `idx` options, and throws a `TypeError` otherwise
- `WalkFile` reads walk format 8, which gfa-to-tabix 5.0 and later write
- jbrowse.org hosts the HPRC release 2.1 walk files, each with its `.idx`:
  `https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38-f8.GRCh38.walks.gz` and
  `https://jbrowse.org/demos/hprc/hprc-v2.1-mc-chm13-f8.CHM13.walks.gz`

## Reading a range

```ts
const { tables, refName, sample } = await walks.getGraphForRange(
  'chr6',
  160_614_797,
  160_647_758,
  { haplotypes: ['HG00097', 'HG00133'] },
)
```

- `getGraphForRange(refName, start, end, opts)` returns the graph over the range
  as [tables](#tables), with the sequence's full `refName` and the reference
  `sample`
- `haplotypes` is a list of PanSN prefixes such as `HG00097` or `HG002#1`. The
  cut keeps the reference and the haplotypes named, or every haplotype when the
  list is empty
- `sequences: true` adds `sequences`, each node's bases by id, from a file
  written with `--sequences`

### Naming the range

- `refName` is a contig of the reference sample, `chr6`, or the file's full name
  for it, `GRCh38#0#chr6`. The reference sample is the header's first reference
  unless `referenceSample` names another
- A sequence the file lacks throws a `MissingSequenceError`, whose `sequences`
  lists the ones it has
- A contig on more than one haplotype of the reference sample throws an
  `AmbiguousSequenceError`. In a file built with `--refs HG002`, `chr1` is both
  `HG002#1#chr1` and `HG002#2#chr1`, and `referenceSample: 'HG002#1'` or the
  full name picks one
- A file with a diploid reference sample has one coordinate system per reference
  haplotype. Where both pass a node, another haplotype's step there may be filed
  under either, so its walk can show in windows of one and in none of the other
- `refNameLookup` and `resolveRefName` resolve a contig the same way against any
  list of sequence names

## Options

| Option                | Effect                                                                                 | Default        |
| --------------------- | -------------------------------------------------------------------------------------- | -------------- |
| `haplotypes`          | The walks to keep beside the reference, by PanSN prefix                                | every walk     |
| `referenceSample`     | The sample whose contig `refName` names                                                | header's first |
| `context`             | bp of reference kept past each end of the range                                        | 1000           |
| `contained`           | `true` keeps each walk within the window; see [Contained cuts](cut.md#contained-cuts)  | false          |
| `attachBp`            | The most bp of own sequence a [departure](cut.md#departures) gets its steps for        | 65,536         |
| `stepBudget`          | The most steps a cut keeps                                                             | 4,000,000      |
| `byteBudget`          | The most compressed bytes the walk and attach reads fetch                              | 8 MB           |
| `departureByteBudget` | The most compressed bytes the departure read fetches                                   | 2 MB           |
| `nodeLimit`           | The most nodes a cut keeps, placeholders not among them                                | no limit       |
| `placeholderLimit`    | The most placeholders a cut keeps                                                      | no limit       |
| `links`               | `'walked'`, `'between'` or `'none'`; see [Links](#links)                               | `'walked'`     |
| `sequences`           | `true` returns each node's bases, from a file written with `--sequences`               | false          |
| `expand`              | Which abridged walks to read in full: PanSN prefixes, or a function of the walk's name | none           |
| `signal`              | An `AbortSignal` that stops the reads                                                  |                |
| `budgetNames`         | The words a `RegionTooLargeError` message uses for the step and byte budgets           |                |

A cut over a budget or limit throws a `RegionTooLargeError`, whose `fitsBp` is
about how wide a window would fit. [Large genomic regions](reads.md) covers the
error, placeholders and `expand`.

## Cut presets

```ts
import { CUT_PRESETS } from '@gmod/gfa-tabix'

const { tables } = await walks.getGraphForRange('chr6', start, end, {
  ...CUT_PRESETS.browser,
  haplotypes: ['HG00097'],
})
```

| Option                | `browser` | `figure` |
| --------------------- | --------- | -------- |
| `context`             | 1000      | 1000     |
| `contained`           | false     | false    |
| `attachBp`            | 65,536    | 65,536   |
| `stepBudget`          | 4,000,000 | Infinity |
| `byteBudget`          | 8,000,000 | Infinity |
| `departureByteBudget` | 2,000,000 | Infinity |

- `CUT_PRESETS.browser` is what a cut does with no options, for a track that
  follows a browser window
- `CUT_PRESETS.figure` lifts the budgets, for a static figure of one region
- `context`, `contained` and `attachBp` are the same in both, so a cut within
  the browser's budgets returns the same tables with either preset
- `WALK_BUDGETS` holds the default budgets and `attachBp`, and `WALK_CONTEXT`
  the default context. All of these objects are frozen, so spread a preset into
  the options as above

## Tables

`tables` holds the cut as typed arrays, one entry per node, link or walk.
`graphTablesGFA` in
[@jbrowse/bandage-core](https://github.com/GMOD/bandage-core) writes the nodes,
links and walks as GFA.

```ts
const { nodes, walks } = tables
walks.names.forEach((name, w) => {
  for (let k = walks.offsets[w]; k < walks.offsets[w + 1]; k++) {
    const node = walks.steps[k]
    console.log(name, nodes.names[node], walks.reversed[k], nodes.lengths[node])
  }
})
```

| Column                                      | Meaning                                                                               |
| ------------------------------------------- | ------------------------------------------------------------------------------------- |
| `nodes.names`, `nodes.lengths`              | Each node's name and bp, sorted by position                                           |
| `nodes.refs`, `nodes.refNames`              | The sequence the node is placed on, as an index into `refNames`                       |
| `nodes.starts`                              | The node's position on that sequence                                                  |
| `nodes.ranks`                               | 0 for a node on the reference                                                         |
| `nodes.passes`                              | How many times the reference passes the node                                          |
| `nodes.reversed`                            | 1 where the reference steps the node reversed                                         |
| `nodes.placeholders`                        | 1 for a [placeholder](reads.md#abridged-walks)                                        |
| `links.from`, `links.to`                    | The two nodes of each link, as node indices                                           |
| `links.strands`                             | Bit 0 is set where the `from` end is reversed, and bit 1 where the `to` end is        |
| `walks.names`                               | Each walk, `sample#haplotype#contig`, the reference first and the rest by name        |
| `walks.starts`, `walks.ends`                | The walk's stretch of its contig                                                      |
| `walks.offsets`                             | Walk `w` has the steps from `offsets[w]` up to `offsets[w + 1]`                       |
| `walks.steps`, `walks.reversed`             | Each step's node index, and 1 where the walk steps the node reversed                  |
| `walks.refStarts`                           | Each step's position on the reference, NaN for a step off it                          |
| `walks.trimmedBefore`, `walks.trimmedAfter` | bp of the walk the cut read past either end of what it kept                           |
| `departures`, `gaps`                        | One row per [departure](cut.md#departures) and per [gap](cut.md#gaps)                 |
| `revisited`                                 | How many nodes the reference passes more than once, NaN where the header has no count |

- A walk that leaves the part of the file the cut read and comes back has more
  than one entry in `walks`, one per [fragment](cut.md#fragments)
- A node the reference passes more than once is placed at its pass over the
  window, the earliest of several. A node off the reference is placed on the
  walk whose name sorts first among those that step it, so its place can change
  with the window and the walks kept
- Unless `revisited` is 0, a run of steps off the reference may be a copy of
  sequence the reference passes elsewhere

### Links

- `'walked'` lists a link for each two consecutive steps of a fragment or of a
  departure's own steps. A link that no kept walk takes is left out even when
  both its nodes are in the cut
- `'between'` adds the links other walks take between the cut's nodes, such as a
  deletion that no kept walk has. This list depends on which tiles the cut read
- Each link is written the way round with more `+` ends, and with the smaller
  first handle on a tie: a walk `>7>3` gives `L 7 + 3 +`

## The GFA tabix file

```ts
import { GfaTabixFile } from '@gmod/gfa-tabix'
import { TabixIndexedFile } from '@gmod/tabix'

const gfaTabix = new GfaTabixFile({
  graph: new TabixIndexedFile({ url: graphUrl, tbiUrl: `${graphUrl}.tbi` }),
})
const gfa = await gfaTabix.getGraphForRange(
  'GRCh38#0#chr6',
  160_614_797,
  160_647_758,
  { hops: 1 },
)
```

- `GfaTabixFile` takes `{ graph }`, format 2's one `<prefix>.graph.bed.gz`, or
  `{ segments, links }`, format 1's `<prefix>.segs.bed.gz` and
  `<prefix>.links.bed.gz`. Each is a `TabixIndexedFile` of
  [@gmod/tabix](https://github.com/GMOD/tabix-js), which the caller installs
- `getGraphForRange` returns the segments and links over the range as GFA text,
  sorted, with LN, SN, SO and SR tags on each segment. Both formats give the
  same text
- `hops` is for gfa-to-tabix's contig layout, where each round follows links out
  from the window. The default anchored layout needs none
- A segment tagged `mv:i:1` belongs to a component that gfa-to-tabix's
  `--max-component` left on separate coordinates, for a viewer to draw as a stub
- `parseSegmentRow`, `parseLinkLine` and `segmentSamples` read single rows of
  format 1, and `parseGraphHeader` and `parseGraphRow` read format 2's

## Errors

| Error                    | Thrown when                                                             |
| ------------------------ | ----------------------------------------------------------------------- |
| `WalkFormatError`        | The file is not walk format 8, or a unit does not match the directory   |
| `MissingSequenceError`   | The file has no sequence of that name                                   |
| `AmbiguousSequenceError` | The contig is on more than one haplotype of the reference sample        |
| `RegionTooLargeError`    | A cut is over a budget or limit ([Large genomic regions](reads.md))     |
| `DepartureReadError`     | The departure read is over `departureByteBudget`                        |
| `RunsRefusedError`       | `getRunsForRange` on a file whose header does not state `revisited:i:0` |

## File information

| Method                            | Returns                                                               |
| --------------------------------- | --------------------------------------------------------------------- |
| `getHeader()`                     | The reference sample, every haplotype with rows, and the tile size    |
| `getReferenceSequences()`         | The sequences a range can be on, in full PanSN: `GRCh38#0#chr6`       |
| `getSequenceName(sample, contig)` | The file's name for a contig on a sample, undefined where it has none |
| `getSequenceLength(refName)`      | Where the sequence's reference ends, 0 for a sequence the file lacks  |
| `getUnits(refName)`               | The sequence's unit table, empty for a sequence the file lacks        |

## Going further

| To get                                                    | Call                | Page                                            |
| --------------------------------------------------------- | ------------------- | ----------------------------------------------- |
| Fragments, departures and gaps as objects                 | `getCutForRange`    | [Cut objects](cut.md)                           |
| A region too large to load, or a file with abridged walks | `expand`, budgets   | [Large genomic regions](reads.md)               |
| Each walk as stretches of reference, with no step decoded | `getRunsForRange`   | [Runs](runs.md)                                 |
| How many haplotypes pass through each node and link       | `getCountsForRange` | [Counts](counts.md)                             |
| Haplotypes grouped by the structural changes they make    | `getFormsForRange`  | [Structural forms](structural-forms.md)         |
| The bytes a cut would read, before making it              | `getBytesForRange`  | [Units a cut reads](reads.md#units-a-cut-reads) |
| Every pass of the reference over a range                  | `getReferenceNodes` |                                                 |
