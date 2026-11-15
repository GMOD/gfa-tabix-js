# @jbrowse/graph-walks

@jbrowse/graph-walks reads a window of a pangenome graph from
[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix) files by range requests.

- **In:** a `gfa-to-tabix --walks` file, `<prefix>.<sample>.walks.gz` with its
  directory `<prefix>.<sample>.walks.gz.idx`, or a segment index,
  `<prefix>.segs.bed.gz` and `<prefix>.links.bed.gz` with theirs or format 2's
  `<prefix>.graph.bed.gz`, plus a reference region such as
  `chr6:160,614,798-160,647,758`
- **Out:** from a walk file, the reference over the region, every haplotype's
  walk through it, and the nodes and links they visit, as typed arrays; from a
  segment index, the segments and links under the region as GFA text
- The reader reads walk format 8, which gfa-to-tabix 5.0 and later write, and
  refuses an older format with a `WalkFormatError` asking for a rebuild, and a
  newer one asking for a newer reader
- It reads a walk file through
  [generic-filehandle2](https://github.com/GMOD/generic-filehandle2)'s handles
  and gunzips with
  [@gmod/bgzf-filehandle](https://github.com/GMOD/bgzf-filehandle); a segment
  index through any object with
  [@gmod/tabix](https://github.com/GMOD/tabix-js)'s `TabixIndexedFile` methods

The JBrowse graph plugin, BandageJS, ggbandage and the sequenceTubeMap fork read
through it.

## Install

```bash
npm install @jbrowse/graph-walks
```

A segment index also needs `@gmod/tabix`.

## Usage

```ts
import { WalkReader, walksKept } from '@jbrowse/graph-walks'

const reader = WalkReader.open(
  'https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38-f8.GRCh38.walks.gz',
)

const refName = 'GRCh38#0#chr6'
const { graph, kept, fragments, departures } = await reader.cut(
  refName,
  { start: 160_614_798, end: 160_647_758 },
  { keep: walksKept(['HG00097', 'HG00133'], refName) },
)
const tables = graph.tables(kept, fragments, departures)
```

- `WalkReader.open(url)` reads the walk file at `url` and its directory at `url`
  plus `.idx`, put before any query string or fragment, by HTTP range requests
  through generic-filehandle2's `RemoteFile`. Its second argument opens each
  location another way: `url => new LocalFile(url)` for a local file
- `WalkReader.fromFilehandles({ data, directory })` takes the two handles
  already open: `data` the `.walks.gz`, `directory` its `.idx`. Any object with
  a generic-filehandle2 handle's `read`, `readFile`, `stat` and `source` will do
- A file gfa-to-tabix 5.0 or later wrote, in walk format 8; the reader refuses
  any other format and says to rebuild the file. jbrowse.org hosts the HPRC
  release 2.1 walk files in format 8, each with its `.idx`:
  `https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38-f8.GRCh38.walks.gz` and
  `https://jbrowse.org/demos/hprc/hprc-v2.1-mc-chm13-f8.CHM13.walks.gz`
- `tables` lists the nodes with their lengths and reference positions, the links
  the cut's steps walk, each walk's steps, and every departure from the
  reference over the window ([Departures](#departures)). `graphTablesGFA` in
  [@jbrowse/bandage-core](https://github.com/GMOD/bandage-core) writes the
  nodes, links and walks as GFA
- **A cut lists the links its steps walk, and no other:** one for each two
  consecutive steps of a fragment, and one for each two consecutive own steps a
  departure holds. A departure's first own step is linked to its leave step
  where a fragment of its walk holds that step, the one on the reference that
  ends at the departure's `hapStart`, and its last own step to its rejoin step,
  the one that starts where its own bp end, likewise. A departure that leaves on
  the window and rejoins past its context so has no link to the node it rejoins,
  though another walk, or another pass of the reference, keeps that node in the
  cut. Where a departure holds the steps after its leave and the steps before
  its rejoin and not what lies between, the two stretches are not linked, and
  its `breakAt` says where they part. A link no walk of the cut takes is not
  listed though both its nodes are, so the links follow from the fragments and
  departures a cut returns, and are those of cutting each of its walks alone,
  however the file is tiled or abridged and whichever tiles the cut read. The
  reference is one of the cut's walks: a `keep` that leaves it out keeps the
  reference's nodes on the window and none of its links, and `walksKept` always
  keeps it
- A node the reference passes more than once, as in a collapsed repeat, takes
  the position of its pass over the window, and each walk's step on it the pass
  its row names: `tables.walks.refStarts` gives every step's position on the
  reference, NaN for a step off it, so a reader of the tables never picks a pass
  itself
- `tables.revisited` is the walk file header's `revisited`, how many nodes the
  reference passes more than once, NaN where the header states none. Unless it
  is 0, a run of steps off the reference may be a copy of reference sequence the
  reference passes elsewhere, not sequence the reference lacks, so a reader of
  the tables calls such a run the haplotype's own only where it is 0
- A node off the reference lies on the walk whose name sorts first among the
  fragments and departures with a step on it, at the first such step. Like the
  links, it follows from what the cut returns, however the file is tiled or
  abridged and whichever tiles the cut read. So its place moves with the window
  and the walks kept: a wider window can return an earlier step on it or a walk
  whose name sorts first, and a `keep` that leaves that walk out places it on
  another. A host should not hold a node's place fixed from one cut to the next
- `reader.header()` lists the reference sample and every haplotype with rows,
  and the tile size, `tile`, read from the walk file's header unit
  - `graph` is the hash of the GFA's segments, as its node and link index states
    it, and `leftout` how many of its lines have no rows in the file; each is
    absent where the file doesn't state it
- `reader.referenceSequences()` lists the sequences a region can be on, in full
  PanSN: `GRCh38#0#chr6`
- `reader.units(refName)` is a sequence's unit table,
  `reader.sequenceLength(refName)` where its reference ends, and
  `reader.referenceNodes(refName, region)` every pass of the reference over a
  region
- `reader.bytesFor(refName, region, { context, attachBp })` is the compressed
  bytes a cut of the region reads when it keeps every walk and expands none
  ([What a cut reads](#what-a-cut-reads))
- Links are written the way round with more `+` ends, as a walk along the
  reference takes them, and with the smaller first handle on a tie: a walk
  `>7>3` gives `L 7 + 3 +`

## Changes in 7

graph-walks 7 reads walk format 8 only (gfa-to-tabix's ADRs 0012 and 0013), and
reads an abridged piece it does not load from its tile's `blocks` unit: its
steps on the reference, exactly, and one placeholder step for each run of steps
off it. The two-step fragment that stood for such a piece, and everything that
netted it, is gone. What the `abridged` row alone settles, it reads no `blocks`
unit for: which pieces to build, where a piece ends on its contig, and a
record's own steps in a tile fetched to attach them.

- **Refused:** any format but 8, at open, with a `WalkFormatError` that says to
  rebuild the walk file with gfa-to-tabix 5.0 or later; `WALK_FORMAT` is 8. A
  cut that decodes a node id at or above 2^31, or below 0, which the format does
  not allow, throws a `WalkRowError` naming the limit rather than read it as
  another id. `WalkRowError`, which a cut throws for a row it cannot decode, is
  exported beside `WalkFormatError`. An `abridged` row without `hapBp`, as a
  format 8 file built before gfa-to-tabix 5.0.0 was released has it, is refused
  as no abridged row: rebuild the file
- **Refused where rows disagree:** a `blocks` row or a detail row that does not
  add up to its `abridged` row's `hapBp`; an `abridged` row that ends elsewhere
  than its line's next piece starts; a record whose own steps in a piece are not
  a run at one end of it, or not the run its `blocks` row has there. Each is
  checked on the rows a cut reads
- **Links:** a cut lists the links its steps walk and no other
  ([Usage](#usage)). graph-walks 6 also listed links of walks a cut left out,
  where the file kept their rows in a tile, and a link between the reference's
  neighbouring rows of a tile, which the graph may not have. A cut of a few
  haplotypes so lists fewer links than it did
- **Nodes off the reference:** each lies at a step the cut returns
  ([Usage](#usage)). graph-walks 6 placed one by every walk row of the tiles a
  cut read, so a walk the cut left out, or a piece of an abridged walk it did
  not build, could place it
- **`WalkFragment.abridged`** is a list, one `{ steps, detail }` per piece of
  the fragment the cut did not load, in the walk's order, and `unloadedSteps`
  counts the fragment's steps that are theirs. The fragment is an ordinary one:
  `ids`, `refStarts`, `against`, `hapStart`, `hapEnd`, `trimmedBefore` and
  `trimmedAfter` mean what they mean on any fragment, and an unloaded piece
  joins the pieces beside it
- **`WalkDeparture.abridged`** lists, the same way, the pieces the cut did not
  load with a step among the departure's `steps`. A piece can reach a cut there
  alone: a record's leave piece or rejoin piece with no step on the drawn
  stretch, in a tile of the window or in one fetched to attach its own steps
- **Placeholders:** a step with `isPlaceholder(id)`, an id below 0, is a run of
  steps off the tile's reference path that the cut did not load. Its `WalkNode`
  has `placeholder: true`, the run's bp as its length, and its place on its
  walk's contig. An id is its cut's alone: in a cut each run has one id and each
  id one run, and the same cut gives the same ids and the same tables whatever
  its reader cut before or cuts beside it, while another cut gives the run
  another id. A placeholder's name in the tables, `walk:start-end`, is its walk
  and its stretch of the contig whatever the cut, so a consumer that merges cuts
  matches placeholders by name. Placeholders are in `kept`, in a departure's
  `steps`, in a gap's `steps` segment and in the links
- **`GraphTables`:** `nodes.placeholders` flags placeholder nodes, each named
  for its walk and its stretch of the contig, `HG002#1#chr1:118-121`, which is
  no node's name. `walks.abridgedSteps` is the steps of a walk's unloaded
  pieces, each piece counted on the first walk holding it, NaN for a walk with
  none; `departures.abridgedSteps` counts, the same way, the pieces no walk
  holds. A gap's steps are its walk's departures' steps, so a walk's unloaded
  steps are its walks' `abridgedSteps` and its departures' together.
  `walks.abridgedBp` is gone: the placeholders' lengths are the bp not loaded
- **`GraphTables.revisited`** is required, where every other addition to the
  tables is optional, so a host that builds its own tables now has to pass it:
  how many nodes its reference passes more than once, or NaN where it states
  none ([Usage](#usage))
- **`walks.refStarts`** is each fragment's `refStarts`, step for step, as
  `graph.tables` reads them: NaN for a step of a fragment with no `refStarts`
  or with fewer than its `ids`, and none past its `ids` from one with more
- **Gaps:** `GapKind` and `GAP_KINDS` lose `abridged`, so `own` is kind 2 and
  `unread` kind 3 in a gap table, and a segment has no `abridged` flag. An
  unloaded piece's steps outside a fragment are `read`, as a loaded piece's are
- **Forms:** `Departure.abridged`, the `net` mark kind and
  `StructuralForms.abridged` are gone, since an unloaded piece gives the
  departures its loaded one gives. `formsCut` makes one cut and reads no detail
  unit: its `expand` and `partial` options and `ExpansionTooLargeError` are gone
- **`expand`** reads the detail units of the pieces the cut builds, and no
  other, told by the `abridged` rows and the departure records alone
  ([Abridged walks](#abridged-walks)). A cut reads no `blocks` unit for a walk
  it expands
- **When a cut reads a `blocks` unit:** for each tile under its window and
  context that abridges a walk it keeps and does not expand, and for nothing
  else ([What a cut reads](#what-a-cut-reads))
- **Attach from the row:** a tile fetched to attach a record's own steps gives
  them, for a piece the cut does not expand, as a placeholder from the piece's
  `abridged` row and the record, with no `blocks` unit read
- **Own steps where a walk leaves or rejoins:** a departure gets its own steps
  only where it leaves or rejoins the reference on the window or its context
  ([Departures](#departures)). One that jumps over them comes back with its ends
  and its own bp and no `steps`, keeps no node in the cut, and costs no read
  past its line in the departure blocks. graph-walks 6 gave every departure
  within `attachBp` over the window its steps, and fetched a tile for each
  record it lacked a piece of, however far off its ends lay.
  `leavesOrRejoinsOn(d, span)` tells the two apart
- **`bytesFor`** counts the tiles a cut fetches to attach, which it reads the
  window's departure blocks to learn, and takes `attachBp` and
  `departureByteBudget` as a cut does. It is the bytes a cut that keeps every
  walk and expands none reads
- **Later additions:** a tag, `NAME:T:value`, after a row's fixed columns is
  ignored on every row kind and on the `#walks` line; so are a header or
  directory line whose first word the reader does not know, and a unit of a kind
  it does not know. A column after a row's fixed columns that is no tag is
  refused, as is any column more on a `#unit`, `#entry` or directory sequence
  line
- **`WALK_BUDGETS`** is frozen, as `CUT_PRESETS` is
  ([Cut presets](#cut-presets))
- **Budgets:** a step built for an unloaded piece counts toward `stepBudget` as
  any step does. A placeholder does not count toward `nodeLimit`, which bounds
  the nodes a layout places; a cut's `placeholders` says how many it returned,
  for a consumer that lays them out. `byteBudget` counts the `blocks` units a
  cut reads, and `bytesFor` says them with everything else the cut reads
- **Contained cuts:** a fragment can end on a placeholder of exactly `context`
  bp, the part of a longer run within the cut's context, with the rest of the
  run as another placeholder outside it
- New: `isPlaceholder`, the `Abridged` type, `WalkNode.placeholder`,
  `WalkFragment.unloadedSteps`, `WalkDeparture.abridged`,
  `departures.abridgedSteps`, `placeholders` on a cut, and `abridgedWalks`,
  `abridgedText` and the `AbridgedWalk` type, which say a cut's runs off the
  reference not loaded by their bp ([Abridged walks](#abridged-walks))

## A region by contig

```ts
const { tables, refName, sample } = await reader.cutRegion(
  { refName: 'chr6', start: 160_614_798, end: 160_647_758 },
  { haplotypes: ['HG00097', 'HG00133'], nodeLimit: 100_000 },
)
```

- `cutRegion` resolves `chr6` to the reference sample's sequence,
  `GRCh38#0#chr6`, and keeps the haplotypes named, every one when none are
- `referenceSample` names another sample than the header's first reference. A
  sequence the file lacks throws, naming the ones it holds
- The cut options below apply too, with `haplotypes` in place of `keep`

## Cut options

| Option                | Effect                                                                                                |
| --------------------- | ----------------------------------------------------------------------------------------------------- |
| `keep`                | Which walks to decode, the reference's among them, by `sample#haplotype#contig`; every one by default |
| `context`             | bp of reference kept past each end of the region; default 1000, `WALK_CONTEXT`                        |
| `stepBudget`          | The most steps a cut keeps, departures' own steps included; 4,000,000 by default                      |
| `byteBudget`          | The most compressed bytes the walk and attach reads fetch; 8 MB by default                            |
| `departureByteBudget` | The most compressed bytes the departure read fetches; 2 MB by default                                 |
| `attachBp`            | The most bp of own sequence a departure gets its steps attached for; default 65,536                   |
| `attachFetches`       | Deprecated and ignored: the tiles of every record that gets its steps are fetched                     |
| `nodeLimit`           | The most nodes a cut keeps; no limit by default                                                       |
| `contained`           | `true` keeps each walk within the window; see below                                                   |
| `links`               | `false` lists no links, for a caller that draws no edges                                              |
| `sequences`           | Keep each node's bases in `graph.sequences`, from a `--sequences` file                                |
| `signal`              | An `AbortSignal` that stops the reads                                                                 |
| `expand`              | Which abridged walks to read in full, by `sample#haplotype#contig`; none by default                   |
| `budgetNames`         | The words a `RegionTooLargeError` message uses for the step and byte budgets                          |

- A cut over a budget or its node limit throws a `RegionTooLargeError`. The
  error's `fitsBp` holds about how wide a window would fit, and JBrowse shows it
  as a zoom-in notice
- Over the byte budget, the error's `bytes` is what the window would fetch, so a
  viewer with a fixed window can offer to load it anyway with `byteBudget`
  raised to `bytes`. `zoomHelps` is false where no narrower window fits
- The cut checks the node limit as it decodes each tile, against the nodes it
  keeps whatever else it reads, and stops decoding at the first tile over it
- A window past `departureByteBudget` throws a `DepartureReadError`, which is
  not a `RegionTooLargeError`: a departure's record is filed over its whole
  stretch, so zooming in does not lower that read. A collapsed repeat that
  haplotypes hop between the copies of is the usual cause
- `WALK_BUDGETS` holds the default budgets and `attachBp`, what a browser window
  holds; a caller that wants no limit passes `Infinity`
- The byte budgets bound what a cut fetches: the unit table gives the exact
  compressed bytes of the window's tiles and of its departure blocks before
  either is fetched, and the cut refuses either past its budget. The tiles
  fetched for records' own steps count against `byteBudget` with the window's,
  and a cut they would take past it throws the same `RegionTooLargeError` before
  fetching any of them, so a window attaches the own steps of every record that
  [gets them](#departures) or of none. The `blocks` units a cut reads count
  against `byteBudget` too. The detail records `expand` names count only against
  `stepBudget`, once read; the writer moves paths out of any tile whose rows
  would compress past its header's `abridgebytes`
- A narrow window reads whole tiles, so a cut reads several times the steps it
  keeps. The rows a cut decodes from the window's tiles, and from the detail
  records `expand` names once they are decoded, may hold up to 8 times
  `stepBudget` steps (`READ_STEPS_PER_KEPT`), a step along the reference
  counting a quarter; past that the cut throws a `RegionTooLargeError` before it
  cuts. The rows of tiles fetched for records' own steps are not counted there,
  and the own steps they give count against `stepBudget` with the rest

## Cut presets

```ts
import { CUT_PRESETS } from '@jbrowse/graph-walks'

const cut = await reader.cut(refName, region, {
  ...CUT_PRESETS.browser,
  keep: walksKept(['HG00097'], refName),
})
```

- `CUT_PRESETS.browser` is for a browser track at a window: the JBrowse graph
  plugin, BandageJS, jb2hubs and the sequenceTubeMap fork. It is what a cut does
  with no options
- `CUT_PRESETS.figure` is for a static figure of one region: ggbandage and
  ggtubemap. A figure waits for its one region and has no window to protect, so
  it lifts the step and byte budgets
- The cut rules, `context`, `contained` and `attachBp`, are the same in both, so
  one region draws one graph in every reader: a cut within the browser's budgets
  gives the same tables with either. A figure keeps `attachBp`, since a
  departure past it would bring nodes into the figure that the browser leaves
  out
- `WALK_CONTEXT` is the context both keep
- `CUT_PRESETS`, each preset and `WALK_BUDGETS` are frozen: a reader spreads
  one into its options, as above, so no reader changes them for another
- `keep`, `expand`, `nodeLimit`, `links` and `sequences` stay each reader's own:
  which walks it shows, which it reads in full, how many nodes its layout
  places, and what it draws

| Option                | `browser` | `figure` |
| --------------------- | --------- | -------- |
| `context`             | 1000      | 1000     |
| `contained`           | false     | false    |
| `attachBp`            | 65,536    | 65,536   |
| `stepBudget`          | 4,000,000 | Infinity |
| `byteBudget`          | 8,000,000 | Infinity |
| `departureByteBudget` | 2,000,000 | Infinity |

## Fragments

Each of `fragments` is one walk's stretch through the cut: `name`, `fragStart`
(the W line's start), `hapStart`/`hapEnd`, and its steps, from the walk's first
step on the drawn stretch of the reference, the window and its context, to its
last. `trimmedBefore` and `trimmedAfter` count the bp of the walk's run past
either end of the fragment, so they depend on the tiles the cut read. They
include own steps the cut attaches to a departure there, so they say how far the
run goes, not what lies there: [Gaps](#gaps) does. An abridged piece the cut
does not build counts there by its `abridged` row's `hapBp`.

- Fragments come from the tiles overlapping the drawn stretch alone. A tile the
  cut fetches for a record's own steps (see [Departures](#departures)) gives
  that record its steps and nothing else, so no fragment holds a step from it
- A walk's pieces in those tiles join where their piece indices are consecutive.
  Where a walk goes into a tile the cut did not read and comes back, it has a
  fragment either side, each ending at its last step on the stretch, and its
  departure records say where it leaves and rejoins
- A walk's fragments therefore depend on the window, its context and the file
  alone, not on which other walks the cut keeps
- `fragments` come in no stated order. `graph.tables` orders its `walks` itself:
  the reference first, then by name and `hapStart`
- `graph.tables(kept, fragments, departures)` keeps `trimmedBefore` and
  `trimmedAfter` as columns of `walks`

## Gaps

`gaps` states what lies on each walk's contig between its fragments, past its
first and last, and in its departures outside them, as
[ADR 0001](docs/adr/0001-a-cut-states-what-lies-between-fragments.md) decides.
Each has `name`, `edge` (0 between two fragments, 1 a row's start, 2 its end, 3
a departure outside the row), `start`/`end` on the contig, the `fragment` it
follows (or precedes, for a row's start), and `segments` in contig order whose
`bp` sum to `end - start`:

- `steps`: own steps the cut attached to a departure, with their `ids` and
  `rev`; these nodes are in the cut
- `read`: steps of the walk's runs over the window's tiles that the cut kept
  out, such as reference steps off the drawn stretch, or own steps of a
  departure over `attachBp` or of one that neither leaves nor rejoins on the
  window or its context
- `own`: a record's own sequence the cut did not read, with its `missingBp`, in
  no stated order
- `unread`: pieces the cut did not read, with the reference breakpoints the cut
  knows `before` and `after` them, NaN beside an abridged piece it read as its
  bp alone

A row end runs over the bp of the run past the fragment, and over the whole of
any record that reaches it. Tiles fetched for a record's own steps give `steps`
alone, and the cut fetches them for every record that
[gets its own steps](#departures) or refuses the window, so which records get
`steps` follows from each record and the window alone.
`graph.tables(kept, fragments, departures, gaps)` writes them as a `gaps` table,
each stretch's walk by its fragment's index in `walks`, and each step as a node
index.

## Departures

A departure is where a walk's consecutive steps on the reference are not
adjacent there: it leaves the reference at one breakpoint and rejoins it at
another, with `ownBp` of its own sequence between. `cut` returns each departure
once:

- Every departure record whose line overlaps the window and its context, from
  the departure blocks there. gfa-to-tabix writes one for every departure
  between two pieces of a walk. A record comes back even where no reference node
  lies under the window, as in a `#hole`
- Every departure inside one piece of a walk, which the cut computes from the
  piece's steps, whose stretch, from the base beside its leave to the base
  beside its rejoin, overlaps the drawn stretch: the reference's nodes over the
  window and its context, whole, widened by every pass of a walk's step over
  them. A walk's jump inside one piece far from the window, even inside its
  fragment, is not returned

```ts
import { departureKind } from '@jbrowse/graph-walks'

for (const d of departures) {
  // HG00097#1#chr6 leaves 160,620,000 and rejoins 160,674,000 (forward), 120,000 bp of its own
  console.log(
    d.name,
    d.leavePos,
    d.rejoinPos,
    departureKind(d),
    d.ownBp,
    d.steps?.ids.length,
  )
}
```

- `leavePos` and `rejoinPos` are breakpoints on the reference, where the walk's
  sequence stops following it and starts again; `leaveReversed` and
  `rejoinReversed` say the step runs against the reference's orientation at its
  pass
- `hapStart` is where the own sequence starts on the walk's contig and `ownBp`
  how far it runs, 0 for a direct jump; `missingBp` is the part of it between
  the contig's W lines, which the graph lacks
- `departureKind` derives what the departure does from its ends: `switch`, the
  strands differ; `back`, the walk jumps back over reference it walks again;
  `adjacent`, sequence of its own between adjacent steps; `forward`, a jump
  forward over reference it leaves out. The tools make no further call
- `leavePasses` and `rejoinPasses` say how many times the reference passes each
  end's node. `placementChoice(d)` is true where either is above 1: the writer
  chose which pass the step takes, so the departure is its choice, not a jump
  the walk makes, as where a walk follows one copy of a collapsed repeat the
  reference passes several times. A node's `passes` in `graph.nodes`, and the
  tables' `nodes.passes`, say the same of each step on the reference, from the
  ref rows the cut read; a node the cut places off the reference has 0, though
  it may lie on the reference outside those rows; the header's `revisited`
  counts the nodes passed more than once, 0 where every departure in the file is
  a fact of the graph
- **A departure carries its own steps in `steps` where it leaves or rejoins on
  the window or its context,** and has at most `attachBp` of own sequence. It
  leaves there when the base of its leave step beside the breakpoint lies in
  `[start − context, end + context)`, and rejoins there when that base of its
  rejoin step does: for a step along the reference, the base before `leavePos`
  and the base at `rejoinPos`; for one against it, the base at `leavePos` and
  the base before `rejoinPos`. `leavesOrRejoinsOn(d, span)` says so of any
  departure and stretch
- `steps` are the walk's steps between its leave and its rejoin, from the end of
  the leave piece and the start of the rejoin piece. A step of 0 bp is an own
  step, so a record with `ownBp` 0, or no more than its `missingBp`, can carry
  one; own sequence between the contig's W lines has no steps
- `breakAt` is set where the cut does not return `steps` as one run from the
  leave step to the rejoin step: the steps before that index follow the leave
  step and those from it lead to the rejoin step, and the cut's links do not
  take the walk from one to the other. A record's steps are one run where its
  leave and rejoin pieces are consecutive pieces of one W line and the cut
  returns, at each end, the end step in a fragment of its walk or an own step
  beside it; a computed departure's steps always are. So it is 0 where the cut
  returns neither the leave step nor an own step after it, and their count where
  it returns neither the rejoin step nor an own step before it. Like the links,
  it follows from the steps the cut returns and not from the tiles it read, so
  an abridged file cuts as the file built without abridging, and the same steps
  of a walk come back with the same break whatever the window. The `departures`
  table has it as `breakAt`, -1 where not set
- **A departure that jumps over the window and its context** comes back with
  everything its row states and no `steps`: its walk, `hapStart`, `ownBp`,
  `missingBp`, both breakpoints with their strands and passes, and so its
  `departureKind`. What a consumer does not get is the nodes of its own
  sequence, which join the graph at its two ends, off the window and its
  context: a cut that kept them drew them joined to nothing there
- The rule reads bases, not nodes. A cut draws whole the node that straddles an
  end of the window and its context, and a walk that leaves that node at its far
  end, past the stretch, gets no steps: a record states its breakpoints and not
  how far its end nodes run, and `bytesFor` has to tell the rule from the record
  alone. A wider `context` reaches it
- The rule is the same for a record and for a departure the cut computes inside
  one piece, so what a cut keeps does not turn on where the writer cut a walk
  into pieces
- **`attachBp` bounds what one departure brings into a cut.** A walk that leaves
  the reference on the window for megabases of its own would bring every node of
  them into the cut; past `attachBp` the departure states its bp and no steps,
  and the cut fetches no tile for them. It is not what keeps a window's reads
  down: the rule above is, with `byteBudget` past it
- The steps come from the pieces the cut read. For a record that gets its steps
  and whose leave or rejoin piece lies outside the read, the cut fetches the
  tiles its `ownTiles` column names that it has not read for that walk, and
  takes from them that record's own steps alone. It fetches them for every such
  record, or throws a `RegionTooLargeError` where they would take the cut past
  `byteBudget`, so whether a record gets its steps does not depend on the other
  walks the cut keeps, as gfa-to-tabix ADR 0009 asks
- A record that leaves on the window has its leave piece in the window's tiles,
  and one that rejoins there its rejoin piece, so a record costs a fetch only
  for the piece at its far end, and only where that piece holds own steps
- Own steps in a piece the cut did not load are a placeholder of their bp;
  `expand` naming the walk reads them from its detail unit. The piece holding
  them is the one whose stretch of the contig, `[hapOffset, hapOffset + hapBp)`
  by its `abridged` row, overlaps the record's `[hapStart, hapStart + ownBp)`,
  and the overlap is the placeholder's stretch
- `fromSteps` is true for a departure computed from a piece's steps
- Those steps' nodes are kept, so the cut draws a short insertion or replacement
  at the window's edge whole, and counts them against `stepBudget`
- `graph.tables(kept, fragments, departures)` keeps them as the `departures`
  table: `names`, `hapStarts`, `ownBp`, `missingBp`, `leavePos`,
  `leaveReversed`, `leavePasses`, `rejoinPos`, `rejoinReversed`, `rejoinPasses`,
  `fromSteps`, and each one's own steps as node indices at `offsets[k]` to
  `offsets[k + 1]` of `steps`, with `reversed`

## Abridged walks

- In a tile over its budget, gfa-to-tabix moves the costliest paths' rows to
  detail records and leaves one `abridged` row per moved piece, and a `blocks`
  unit for the tile with a row per piece: the piece in path order as its runs on
  the reference and the bp of its runs off it. A satellite array in a pggb graph
  is the usual case
- A cut builds the same pieces of an abridged walk whether or not it expands the
  walk, and tells which from the `abridged` rows and the departure records
  alone:
  - those whose span lies over the drawn stretch
  - those between two such pieces of a line, by piece number
  - those holding the own steps of a record that [gets them](#departures)
  - in a contained cut with context, the piece whose last steps the context
    reaches before a piece over the drawn stretch
- With `expand` naming the walk, the cut reads those pieces from their detail
  units and cuts them in full, as a file without abridging would. Otherwise it
  builds them from the tile's `blocks` unit
- **An abridged file, expanding every walk a cut keeps, cuts as its twin:** the
  file gfa-to-tabix writes from the same graph with the same `--refs`, `--tile`
  and `--stay`, and `--abridge-steps 0 --abridge-bytes 0`. Over any window both
  files cut within their budgets, the two give the same fragments, links, node
  places, departures with their steps and `breakAt`, and gaps, except:
  - an `unread` gap segment's `before` or `after` beside an abridged piece the
    cut did not build is NaN, where the twin states the breakpoint
    ([Gaps](#gaps))
  - `fragments` come in no stated order, so the two can list them differently;
    `graph.tables` orders the walks of both alike
- Every other piece of a tile under the window is its bp alone, from its
  `abridged` row's `hapBp`: it counts in the `trimmedBefore` or `trimmedAfter`
  of the fragment beside it and as `read` in a gap, and the cut states no step
  of it
- A piece built from its `blocks` row is cut as a loaded one is: kept where a
  step of it lies on the drawn stretch, trimmed to its first and last steps
  there, joined to the pieces beside it, split by a contained cut, linked, and
  read for departures, gaps and forms. Its steps on the reference are the
  reference's own, and each run off it is one placeholder step with that run's
  bp and no node
- A fragment's `abridged` lists the pieces in it the cut did not load, each with
  its `steps` and the `detail` unit that holds it: its sequence `refName`, the
  tile's nominal `start`, and its `name`, `tile:detail`. `steps` counts all the
  piece's steps, which reading the unit decodes
- `isPlaceholder(id)` tells a placeholder step, whose id is below every node's.
  Its node, in `graph.nodes` and in the tables' `nodes` with `placeholders` set,
  has the run's bp as its length, lies on the walk's contig and is named for
  both, `walk:start-end`. The id is its cut's alone, and the name the same in
  any cut. A run off the reference path can hold nodes of the reference, where
  the reference passes a node more than once, and two runs of equal bp may hold
  different sequence: only expanding tells
- A contained cut keeps `context` bp of a walk past its steps on the window. It
  keeps a loaded node whole, so it can run less than that node past `context`;
  of a placeholder it keeps exactly `context` bp, as a placeholder of its own,
  and the rest of the run as another
- A departure's `abridged` lists the pieces the cut did not load with a step
  among its `steps`. They can be the only steps of a piece a cut holds: a
  record's own steps in a piece with no step on the drawn stretch, in a tile of
  the window or in one fetched to attach them
- `graph.tables(kept, fragments, departures)` keeps the steps of a walk's
  unloaded pieces as `walks.abridgedSteps`, each piece on the first walk holding
  it, and of those no walk holds as `departures.abridgedSteps`, each on the
  first departure holding it, NaN for a walk or departure with none. A walk's
  unloaded steps are the two together
- `abridgedWalks(tables)` lists, as an `AbridgedWalk` `{ name, bp }` each, the
  walks with a run off the reference that the cut did not load, by full name, in
  the tables' order. `bp` is the lengths of the placeholders the walk's
  fragments and departures step on, each once: what the cut draws in place of
  sequence it did not load, which expanding the walk would draw as its steps. A
  walk the cut draws exactly is not listed, though pieces of it were not loaded,
  nor is one whose only run off the reference lies past the window, so a cut of
  the twin lists none. It counts no steps: a `blocks` row gives a run's bp and
  not its steps, and `walks.abridgedSteps` is what expanding would decode. A
  reader that lists haplotypes adds up their walks
- `abridgedText(walk, label)` says it in words, `label` the walk's name by
  default, with bp under a kilobase and kb above as bandage-core's walk row
  readouts write them:

```ts
abridgedWalks(tables).map(walk => abridgedText(walk))
// ownends.R over R#0#c:0-30 with no context, among others:
// 'S#1#c: 70 kb off the reference not loaded'
// 'V#1#c: 20 bp off the reference not loaded'
```

- The cut builds whole only what it draws of an unloaded piece: from the first
  block over the drawn stretch to the last, the piece's end steps, and each
  departure that can lie over the stretch. What lies between counts in bp
- `pnpm blocks-equality <walks.gz or URL>` holds the two readings equal on a
  walk file: over random windows on its abridged tiles, the cut that expands
  nothing against the cut that expands every walk, plain, contained and
  contained with context, and against the same cut from a fresh reader. It also
  holds own steps to the departures that leave or rejoin on the window and its
  context, which it tells from the records' columns itself, the tiles read past
  a window's to the ones those records name, each record's placeholders to the
  runs the tiles' `blocks` rows give, the expanding cut to no `blocks` unit, and
  the bytes read to `bytesFor` (`scripts/blocks-equality.ts`)

## What a cut reads

A cut reads these units, each once:

- **The window's tiles and departure blocks:** the `tile` and `departures` units
  whose `start end` overlap the window and its context
- **`blocks` units:** the unit of each of those tiles that abridges a walk the
  cut keeps and does not expand. It builds from it the pieces it draws
  ([Abridged walks](#abridged-walks)), and checks every row it reads against its
  `abridged` row's `hapBp`
- **Detail units:** with `expand`, those of the pieces it builds of the walks
  named, among the window's tiles and then among the tiles fetched to attach
- **Tiles fetched to attach:** for each record within `attachBp` that leaves or
  rejoins on the window or its context ([Departures](#departures)) and whose
  leave piece or rejoin piece the window's tiles lack, the tiles its `ownTiles`
  names. No `blocks` unit is read for them: an unexpanded piece's own run comes
  from its `abridged` row and the record

So a `blocks` unit is read for a tile under the window and its context that
abridges a walk kept and not expanded, and never otherwise: not to attach, not
to choose what to expand, and not for a walk the cut expands.

And a tile past the window's is read for a record that leaves or rejoins there,
and never for one that only lies over it:

- A record lies over every window between its ends, so on a collapsed repeat,
  where the placement hops between copies far apart, most of the records over a
  window end nowhere near it. Each costs its line in the departure blocks, which
  `departureByteBudget` bounds, and nothing more
- What a window reads past its own tiles, `blocks` units and departure blocks
  therefore follows what happens on the window, not how many walks pass over it

`reader.bytesFor(refName, region, { context, attachBp })` adds these up before a
cut is made:

- It counts the window's tiles, each one's `blocks` unit, the window's departure
  blocks, and the tiles named in `ownTiles`, beyond the window's, by the records
  within `attachBp` that leave or rejoin on the window or its context. It reads
  the departure blocks to learn those, as the cut does, and tells which records
  they are from each record's columns alone; it throws the cut's
  `DepartureReadError` past `departureByteBudget`
- **Exact** for a cut that keeps every walk and expands none, within its
  budgets, whatever its `context`, `attachBp` and `contained`
- **An upper bound with `keep` or `haplotypes`:** fewer walks leave out the
  `blocks` unit of a tile that abridges none of them and the tiles only other
  walks' records name, and add no read
- **With `expand`:** the cut reads no more than `bytesFor` but for its detail
  units, and no less than `bytesFor` without the `blocks` units. Which detail
  units it reads is in the tiles' `abridged` rows, which `bytesFor` does not
  read

## Moved components

- gfa-to-tabix 1.5's `--max-component` leaves a component over the limit on its
  own coordinates, so an anchored cut reaches it only through the links into it
- `rgfaCut` tags each segment of such a component that it holds `mv:i:1`
  (`MOVED_TAG`), for a viewer to draw as a stub

## Structural forms

```ts
import {
  describeForm,
  formTag,
  formsCut,
  structuralForms,
} from '@jbrowse/graph-walks'

const region = { refName: 'chr1', start: 109_680_000, end: 109_715_000 }
const haplotypes = await formsCut(reader, region)
const result = structuralForms(haplotypes, { window: region })
result.forms.map(f => `${f.members.length} ${describeForm(f)}`)
// ['306 del 18.4kb', '156 ref']
result.forms.map(f => formTag(f, result))
// ['del 18.4kb', 'ref']
```

- `formsCut` cuts the region uncontained and without links, and returns each
  haplotype's departures from the reference: exact REF span and ALT length,
  inverted bp, and bp walked again by a step back along the reference. It
  expands no walk: an abridged piece's `blocks` row gives the departures its
  detail record gives
- `formsCut` keeps to `FORMS_STEP_BUDGET` steps, or a larger `stepBudget`, and
  throws `RegionTooLargeError` past it
- `formsCut` takes `byteBudget`, and a region whose tiles the table puts past it
  throws `RegionTooLargeError` before any read
- `formsCut` reads each haplotype's departures off the cut's `departures`, every
  one over the region, and its fragments' ends, so a haplotype that leaves the
  region and rejoins the reference far away, or jumps over it whole, is a
  departure with no step on it; `walkDepartures` does the same for any cut
- `structuralForms` groups the haplotypes that make the same changes at every
  informative site. Each change is a `Mark`: an insertion, deletion,
  replacement, inversion, shift or duplication
- Each mark in a result names its `site` in `sites` and its `state` key there. A
  site is `edge`, `shared` (no form differs there) or informative, and
  `formStateAt(result, form, site)` gives a form's state key at any of them.
  `REFERENCE_STATE`, `NO_CALL_STATE` and `RARER_STATE` are the keys for no
  change, no call and a rare change among no-calls
- `placed` counts the haplotypes in a form, so a form's share is its members
  over `placed`
- `form.onlyShifts` is true where every change outside the shared sites is a
  shift and every site has a call: the reference's length with a stretch aligned
  elsewhere, as in AMY1's H3r
- `markTag`, `formTag` and `describeForm` label marks and forms, as
  [Labels](#labels) shows
- Every result holds plain data, so a worker can post it, and regrouping at
  another floor needs no new cut
- [docs/structural-forms.md](docs/structural-forms.md) gives the grouping rules
  and the measurements behind them

| Option          | Effect                                                                                                                    |
| --------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `window`        | The region; marks must touch it, and `auto` reads its width                                                               |
| `floorBp`       | The smallest change a form tells apart; `auto` gives 50 bp, 1 kb or 10 kb by width                                        |
| `minCarriers`   | The fewest haplotypes a state or form needs; 1% of them, at least one, by default                                         |
| `joinBp`        | Departures of one haplotype nearer than this on the reference join; default 50                                            |
| `exclude`       | Samples to leave out, such as a second reference carried as a haplotype; none by default                                  |
| `foldRareForms` | Fold states and forms under `minCarriers` into the nearest common one they agree with, else `~` or `other`; on by default |

### Labels

- `markTag(mark)` is the per-mark label and takes a `Mark` or a `SiteMark`. C4's
  common change gives `del 6.46kb->94bp`
- `formTag` names a form in about 20 characters for a row label: its largest
  change's `markTag`, then ` +N` counting the rest. Given the result, it leaves
  out the changes every form shares
- `describeForm` lists every change's `markTag`, largest first, comma-separated
- Sizes take 3 significant figures with the unit attached (`94bp`, `6.46kb`,
  `1.2Mb`), as `sizeTag` writes them
- A replacement whose smaller side is under 5% (`LOPSIDED_RATIO`) of its larger
  side stays `kind: 'replacement'` with its exact REF span and ALT length, and
  gains `lopsided: 'deletion' | 'insertion'`, so a caller can colour it as a
  loss or gain. Its label reads as one

| Change                         | Token                                   |
| ------------------------------ | --------------------------------------- |
| deletion, insertion            | `del 18.4kb`, `ins 22kb`                |
| lopsided replacement           | `del 6.46kb->94bp`, `ins 178bp->13.9kb` |
| other replacement              | `repl 32.7kb->58.7kb`                   |
| shift, inversion, duplication  | `shift 94.2kb`, `inv 12kb`, `dup 49kb`  |
| shift that changes length      | `shift 94.2kb net +2.81kb`              |
| no change, only shared changes | `ref`, `shared`                         |
| rare change, no call           | `rarer`, `no call`                      |
| more changes, in `formTag`     | `del 32.7kb +1`                         |

## Contained cuts

- By default a cut keeps each walk from its first step on the drawn stretch, the
  window and its context, to its last, whatever it visits between, and follows
  nothing past that stretch, so one haplotype's fragment never depends on which
  others are cut
- `contained: true` splits a walk wherever it steps onto the reference off the
  drawn stretch. Each piece between those steps that steps on the stretch is
  kept from its first such step to its last, so a detour that comes back stays
  whole, with `context` bp of the walk past either end
- A fragment counts the bp of its walk the cut read past either end as
  `trimmedBefore` and `trimmedAfter`
- At AMY1 haplotypes go from a 4 kb window to the other copies of the gene and
  back. The default cut keeps 9,928 nodes there, and the contained cut 182

## Naming the window

- `reader.sequenceFor(sample, refName)` gives the walk file's sequence for a
  contig on a sample: `chr6` on `GRCh38` is `GRCh38#0#chr6`, and a pggb file's
  `chm13#chr1` is `chr1` on `chm13`
- `panSNHaplotypeOrSample(name)` is what a walk file's header lists a path
  under: `HG002#1` for `HG002#1#chr6`, and the sample, `chm13`, for a name with
  no haplotype field, `chm13#chr6`. `walkDepartures` groups paths by it
- `walksKept(haplotypes, refName)` returns the `keep` filter for a cut on
  `refName`, which keeps the reference walk on `refName` and the walks of the
  haplotypes named
- `refNameLookup` and `resolveRefName` resolve a contig the same way against any
  file's sequence names

## Segment indexes

```ts
import { rgfaCut } from '@jbrowse/graph-walks'

const gfa = await rgfaCut(
  { segments: segsFile, links: linksFile },
  { refName: 'GRCh38#0#chr6', start: 160_614_798, end: 160_647_758 },
  { hops: 1 },
)
```

- `segsFile` and `linksFile` are `TabixIndexedFile`s, or any object with its
  `getLines`
- For gfa-to-tabix's format 2, pass `{ graph: graphFile }`, the one
  `<prefix>.graph.bed.gz`, with `getLines` and `getHeaderLines`. `rgfaCut` reads
  its layout from the header, refuses a format other than 2, and writes each
  node's `SB:H:` bits as the `SM:Z:` names format 1 states, in `#samples` order.
  A link end its row gives no carriers takes them from its node row in the same
  query, so a cut of either format gives the same segments and links
- gfa-to-tabix's default, anchored layout files every allele under the reference
  interval its bubble hangs from, so one read of each file holds the whole cut
- Its contig layout files each segment under the segment's own sequence. Each of
  `hops` rounds follows links out from the window, and a last read closes every
  allele reached
- `rgfaCut` writes each segment as `*` with LN/SN/SO/SR tags and the tag
  column's tags, sorted, so one region always gives the same text
- `parseSegmentRow`, `parseLinkLine` and `segmentSamples` read single rows;
  `parseGraphHeader`, `parseGraphRow` and `carrierNames` read format 2's

## How a cut works

- A tile holds the reference's path through it once, each way a haplotype leaves
  that path once, and every other path's piece as the list of ways it takes;
  gfa-to-tabix's
  [walks-format.md](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/walks-format.md)
  describes the rows, and its ADR 0010 the file
- At open the reader reads the directory, refusing any format but 8, then the
  header unit and the sequence's unit table at once, and refuses a directory
  whose first line differs from the header's `#walks` line or whose `bytes`
  differs from the data's length. Over HTTP the length comes from Content-Range,
  so a server must expose that header to cross-origin readers for the check to
  hold; where the reader cannot learn the length it skips the check. Every call
  that reads the file, `referenceSequences`, `units` and `referenceNodes` among
  them, waits on this check
- It refuses, with a `WalkFormatError`, a unit or table that does not hold what
  lists it: one that is not whole gzip members, one whose text ends partway
  through a line, a unit whose `#unit` line is not its table entry's, and a row
  of another sequence or of a kind its unit does not hold. A unit of a kind the
  format does not give it never reads
- A cut selects the `tile` and `departures` units whose `start end` overlap the
  window and its context, skipping any of length 0, merges the byte ranges that
  touch, and fetches each as one single-range request. It gunzips each unit,
  decodes each tile once, and keeps the departure lines overlapping the window
  and its context, each line once
- A departure inside one piece comes from the piece's steps, and one between
  pieces from its record, so each comes back once. A departure's own steps are
  attached where its own bp are within `attachBp` and it leaves or rejoins on
  the window or its context, from its leave piece's steps after the leave and
  its rejoin piece's before the rejoin. Where the read lacks either piece, the
  cut fetches the tiles the record's `ownTiles` names, every one wanted in one
  more read counted against `byteBudget`, and takes from them the records' own
  steps alone; a tile fetched again decodes only the walks it was not decoded
  for. An abridged piece the cut builds is read from its tile's `blocks` unit,
  from its `detail` unit where the cut expands its walk, or, in a tile fetched
  to attach, from its `abridged` row and the record

## License

MIT
