# API reference

Most callers open a `WalkFile` and call `getGraphForRange`. The rest of this
page is for looking things up.

- [Opening a walk file](#opening-a-walk-file)
- [Methods](#methods) and their [options](#options)
- What a cut returns: [tables](#tables), [fragments](#fragments),
  [departures](#departures) and [gaps](#gaps)
- [The GFA tabix file](#the-gfa-tabix-file)
- [Errors](#errors)

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

## Methods

Every range method takes `refName, start, end` and an options object.

| Method              | Returns                                                                                 | More                                            |
| ------------------- | --------------------------------------------------------------------------------------- | ----------------------------------------------- |
| `getGraphForRange`  | The graph over the range as [tables](#tables), with `refName`, `sample` and `sequences` | [Tables](#tables)                               |
| `getCutForRange`    | The cut as objects: `graph`, `kept`, `fragments`, `departures`, `gaps`, `placeholders`  | [Fragments](#fragments)                         |
| `getRunsForRange`   | Each walk as stretches of reference and the bp between them, with no step decoded       | [Runs](runs.md)                                 |
| `getCountsForRange` | How many of the file's haplotypes step each node and link                               | [Counts](counts.md)                             |
| `getFormsForRange`  | Each haplotype's departures, the input of `structuralForms`                             | [Structural forms](structural-forms.md)         |
| `getBytesForRange`  | The compressed bytes a cut of the range would read                                      | [Units a cut reads](reads.md#units-a-cut-reads) |
| `getReferenceNodes` | Every pass of the reference over the range                                              |                                                 |

- `getGraphForRange` takes `haplotypes`, a list of PanSN prefixes such as
  `HG00097` or `HG002#1`, and keeps the reference and every haplotype when the
  list is empty
- `getCutForRange` takes `keep` in place of `haplotypes`, a filter on
  `sample#haplotype#contig` that also decides whether the reference walk is
  kept. `walksKept(haplotypes, refName)` builds the filter `haplotypes` uses.
  The cut's `tables()` writes the same tables `getGraphForRange` returns

Methods that take no range:

| Method                            | Returns                                                               |
| --------------------------------- | --------------------------------------------------------------------- |
| `getHeader()`                     | The reference sample, every haplotype with rows, and the tile size    |
| `getReferenceSequences()`         | The sequences a range can be on, in full PanSN: `GRCh38#0#chr6`       |
| `getSequenceName(sample, contig)` | The file's name for a contig on a sample, undefined where it has none |
| `getSequenceLength(refName)`      | Where the sequence's reference ends, 0 for a sequence the file lacks  |
| `getUnits(refName)`               | The sequence's unit table, empty for a sequence the file lacks        |

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
| `keep`                | `getCutForRange`'s filter on `sample#haplotype#contig`, the reference's walk included  | every walk     |
| `referenceSample`     | The sample whose contig `refName` names                                                | header's first |
| `context`             | bp of reference kept past each end of the range                                        | 1000           |
| `contained`           | `true` keeps each walk within the window; see [Contained cuts](#contained-cuts)        | false          |
| `attachBp`            | The most bp of own sequence a [departure](#departures) gets its steps for              | 65,536         |
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

`getGraphForRange(...).tables` and `getCutForRange(...).tables()` hold the cut
as typed arrays. `graphTablesGFA` in
[@jbrowse/bandage-core](https://github.com/GMOD/bandage-core) writes the nodes,
links and walks as GFA.

| Table        | Columns                                                                                                             |
| ------------ | ------------------------------------------------------------------------------------------------------------------- |
| `nodes`      | `names`, `lengths`, `refs`, `starts`, `ranks`, `refNames`, `passes`, `reversed`, `placeholders`                     |
| `links`      | `from`, `to`, `strands`                                                                                             |
| `walks`      | `names`, `starts`, `ends`, `offsets` into `steps` and `reversed`, `refStarts`, `trimmedBefore`, `trimmedAfter`      |
| `departures` | One row per [departure](#departures), with own steps as node indices at `offsets[k]` to `offsets[k + 1]` of `steps` |
| `gaps`       | One row per [gap](#gaps), each walk by its index in `walks` and each step as a node index                           |
| `revisited`  | The header's count of nodes the reference passes more than once, NaN where the header states none                   |

- `walks` lists the reference first, then the other walks by name and `hapStart`
- `walks.refStarts` is each step's position on the reference, NaN for a step off
  it
- `nodes.passes` counts how many times the reference passes the node, and
  `nodes.reversed` is 1 where the reference steps it reversed
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

## Fragments

A fragment is one walk's stretch through the cut, from the walk's first step on
the window and its context to its last. A walk that leaves for a part of the
file the cut did not read and comes back has one fragment on each side.

| Field                           | Meaning                                                           |
| ------------------------------- | ----------------------------------------------------------------- |
| `name`                          | The walk, `sample#haplotype#contig`                               |
| `fragStart`                     | The start of the walk's W line                                    |
| `hapStart`, `hapEnd`            | The fragment's stretch of the contig                              |
| `ids`, `rev`                    | The steps: node ids and whether each is reversed                  |
| `refStarts`, `against`          | Each step's position on the reference, NaN off it, and its strand |
| `trimmedBefore`, `trimmedAfter` | bp of the walk the cut read past either end of the fragment       |

- A walk's fragments depend on the window, its context and the file, and are the
  same whichever other walks the cut keeps
- `fragments` come in no stated order

### Contained cuts

- By default a fragment keeps whatever the walk visits between its first and
  last step on the window, including steps on the reference far from it
- `contained: true` splits a walk wherever it steps onto the reference off the
  window and its context, and keeps `context` bp of the walk past either end of
  each piece. A detour that comes back to the window stays whole
- At AMY1, haplotypes go from a 4 kb window to the other copies of the gene and
  back. The default cut keeps 9,928 nodes there, and the contained cut 182

## Departures

A departure is where a walk leaves the reference at one breakpoint and rejoins
it at another. A cut returns each departure over the window and its context
once.

| Field                             | Meaning                                                                 |
| --------------------------------- | ----------------------------------------------------------------------- |
| `name`                            | The walk                                                                |
| `leavePos`, `rejoinPos`           | The breakpoints on the reference                                        |
| `leaveReversed`, `rejoinReversed` | Whether the step at each end runs against the reference                 |
| `leavePasses`, `rejoinPasses`     | How many times the reference passes each end's node                     |
| `hapStart`, `ownBp`               | Where the walk's own sequence starts on its contig and how far it runs  |
| `missingBp`                       | The part of `ownBp` between the contig's W lines, which the graph lacks |
| `steps`                           | The own steps, `ids` and `rev`, where the cut attaches them             |
| `breakAt`                         | The index where `steps` is not one run from the leave to the rejoin     |
| `fromSteps`                       | True where the cut computed the departure from a piece's steps          |

- `departureKind(d)` names what the departure does: `switch`, the strands
  differ; `back`, the walk jumps back over reference it walks again; `adjacent`,
  own sequence between adjacent steps; `forward`, a jump over reference the walk
  leaves out
- `placementChoice(d)` is true where either end's node is passed more than once.
  The writer chose which pass the step takes there, as on a collapsed repeat, so
  the departure may not be a jump the walk makes
- A departure gets `steps` where it leaves or rejoins on the window or its
  context and has at most `attachBp` of own sequence. Their nodes are in the cut
  and count against `stepBudget`
- A departure that jumps over the whole window, or has more than `attachBp` of
  own sequence, comes back with every other field and no `steps`
- Own steps in an abridged piece the cut did not load are a placeholder, which
  `expand` naming the walk reads in full

## Gaps

`gaps` lists what lies on each walk's contig between its fragments, past its
first and last fragment, and in its departures outside them.

- A gap has `name`, `start` and `end` on the contig, the `fragment` it follows,
  and `edge`: 0 between two fragments, 1 at the start of a walk row, 2 at its
  end, 3 for a departure outside the row
- Its `segments` are in contig order and their `bp` sum to `end - start`

| Segment kind | Meaning                                                                               |
| ------------ | ------------------------------------------------------------------------------------- |
| `steps`      | Own steps the cut attached to a departure, with `ids` and `rev`; their nodes are kept |
| `read`       | Steps the cut read and left out, such as reference steps off the window               |
| `own`        | A departure's own sequence the cut did not read, with its `missingBp`                 |
| `unread`     | Pieces the cut did not read, with the reference breakpoints `before` and `after` them |

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
