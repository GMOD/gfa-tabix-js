# Cutting a window

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
- A contig on more than one haplotype of the reference sample throws an
  `AmbiguousSequenceError` naming them: `chr1` in a file built with
  `--refs HG002`, which holds `HG002#1#chr1` and `HG002#2#chr1`.
  `referenceSample: 'HG002#1'`, or the full name `HG002#1#chr1`, picks one
- Such a file is two coordinate systems, one per reference haplotype. Where both
  pass a node, another haplotype's step there may be filed under either, so its
  walk can lie wholly under one and show in no window of the other
- The cut options below apply too, with `haplotypes` in place of `keep`

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
- Where a contig is on more than one haplotype of the sample, `sequenceFor` and
  `resolveRefName` throw an `AmbiguousSequenceError` with the names in
  `sequences`: `chr1` on `HG002` where the file holds `HG002#1#chr1` and
  `HG002#2#chr1`. `chr1` on `HG002#1`, or the full name, resolves

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
- `CUT_PRESETS`, each preset and `WALK_BUDGETS` are frozen: a reader spreads one
  into its options, as above, so no reader changes them for another
- `keep`, `expand`, `nodeLimit`, `placeholderLimit`, `links` and `sequences`
  stay each reader's own: which walks it shows, which it reads in full, how many
  nodes and placeholders its layout places, and what it draws

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
[ADR 0001](adr/0001-a-cut-states-what-lies-between-fragments.md) decides. Each
has `name`, `edge` (0 between two fragments, 1 a row's start, 2 its end, 3 a
departure outside the row), `start`/`end` on the contig, the `fragment` it
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
  beside its rejoin, overlaps the drawn stretch: the window and its context,
  widened to the reference's nodes over them, whole, and to every pass of a
  walk's step over them. A window inside a hole in the reference so returns the
  jumps across it, as a record would be. A walk's jump inside one piece far from
  the window, even inside its fragment, is not returned

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

## Links and node places

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
- **`links: 'between'` adds the links other walks take between the cut's
  nodes,** for a host that draws the graph around a few walks: a collapsed
  repeat's jumps between copies, or a deletion no kept walk has. It lists each
  link that a ref row, a route, a walk row stepping R across two ref rows, a
  blocks row, or two consecutive pieces of a line take in the tiles the cut
  reads, where both ends are among its nodes, and reads the window's `blocks`
  units for it, against `byteBudget`. None is outside the GFA. Unlike the
  default, the list depends on the tiles read: a link only a detail record
  holds, or only a tile the cut does not read, is missing, so counts, which
  fetch every record's own steps, can find one it lacks
- A node the reference passes more than once, as in a collapsed repeat, takes
  the position of its pass over the window and its context, the earliest where
  several are, else its earliest pass in the tiles the cut reads; and each
  walk's step on it the pass its row names: `tables.walks.refStarts` gives every
  step's position on the reference, NaN for a step off it, so a reader of the
  tables never picks a pass itself. A node with no pass in those tiles takes the
  first place a fragment or a departure's own steps hold it, by walk name and
  then position on the contig
- A node on the reference states the strand the reference steps it on:
  `tables.nodes.reversed` is 1 where the `ref` row that places the node steps it
  reversed, at the pass `nodes.starts` gives, and a node's `reversed` in
  `graph.nodes` and in `reader.referenceNodes` says the same. A cut's reference
  walk covers the drawn stretch, and a haplotype's walk can bring in a node the
  reference steps elsewhere in the tiles read, which no step of that walk gives
  the reference's strand of; the column gives it for every node with rank 0. It
  is 0 for a node placed off the reference, whose strand along its walk the
  walk's steps give
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
- Links are written the way round with more `+` ends, as a walk along the
  reference takes them, and with the smaller first handle on a tie: a walk
  `>7>3` gives `L 7 + 3 +`

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

## Moved components

- gfa-to-tabix 1.5's `--max-component` leaves a component over the limit on its
  own coordinates, so an anchored cut reaches it only through the links into it
- `rgfaCut` tags each segment of such a component that it holds `mv:i:1`
  (`MOVED_TAG`), for a viewer to draw as a stub
