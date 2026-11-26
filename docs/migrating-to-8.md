# Migrating to graph-walks 8

The names here are graph-walks 8's; [Migrating to 10](migrating-to-10.md) maps
them to @gmod/gfa-tabix's.

## Changes in 8

- **`links`** is `'walked'`, `'between'` or `'none'`: `'walked'` is graph-walks
  7's `true` and the default, `'none'` its `false`, and `'between'` adds the
  links other walks take between the cut's nodes
  ([Links and node places](cutting.md#links-and-node-places)). graph-walks 7
  dropped those links, so a cut of a few haplotypes drew a tandem array's copies
  as one straight run

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
- **Refused where a table leaves bytes out:** a sequence's table whose entries
  do not cover every byte from the end of the table before it, or of the header
  unit, to the table itself, in file order with no gap and no overlap. A table
  that left out a unit's entry read before as a sequence without that unit
- **Links:** a cut lists the links its steps walk and no other
  ([Links and node places](cutting.md#links-and-node-places)). graph-walks 6
  also listed links of walks a cut left out, where the file kept their rows in a
  tile, and a link between the reference's neighbouring rows of a tile, which
  the graph may not have. A cut of a few haplotypes so lists fewer links than it
  did
- **Nodes off the reference:** each lies at a step the cut returns
  ([Links and node places](cutting.md#links-and-node-places)). graph-walks 6
  placed one by every walk row of the tiles a cut read, so a walk the cut left
  out, or a piece of an abridged walk it did not build, could place it
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
- **A contig on two reference haplotypes:** `refNameLookup` maps each spelling
  to every name it fits, a `Map<string, string[]>`, and `resolveRefName`,
  `sequenceFor`, `cutRegion` and `formsCut` throw an `AmbiguousSequenceError`
  where it fits more than one. graph-walks 6 resolved `chr1` on `HG002` to the
  last of `HG002#1#chr1` and `HG002#2#chr1` without a word
  ([A region by contig](cutting.md#a-region-by-contig))
- **`GraphTables.revisited`** is required, where every other addition to the
  tables is optional, so a host that builds its own tables now has to pass it:
  how many nodes its reference passes more than once, or NaN where it states
  none ([Links and node places](cutting.md#links-and-node-places)). A host that
  builds its own `WalkGraph` passes it as the constructor's second argument,
  `new WalkGraph(window, revisited)`, and without it the tables' `revisited` is
  NaN
- **`walks.refStarts`** is each fragment's `refStarts`, step for step, as
  `graph.tables` reads them: NaN for a step of a fragment with no `refStarts` or
  with fewer than its `ids`, and none past its `ids` from one with more
- **Gaps:** `GapKind` and `GAP_KINDS` lose `abridged`, so `own` is kind 2 and
  `unread` kind 3 in a gap table, and a segment has no `abridged` flag. An
  unloaded piece's steps outside a fragment are `read`, as a loaded piece's are
- **Forms:** `Departure.abridged`, the `net` mark kind and
  `StructuralForms.abridged` are gone, since an unloaded piece gives the
  departures its loaded one gives. `formsCut` makes one cut and reads no detail
  unit: its `expand` and `partial` options and `ExpansionTooLargeError` are gone
- **`expand`** reads the detail units of the pieces the cut builds, and no
  other, told by the `abridged` rows and the departure records alone
  ([Abridged walks](reads.md#abridged-walks)). A cut reads no `blocks` unit for
  a walk it expands
- **When a cut reads a `blocks` unit:** for each tile under its window and
  context that abridges a walk it keeps and does not expand, and for nothing
  else ([What a cut reads](reads.md#what-a-cut-reads))
- **Attach from the row:** a tile fetched to attach a record's own steps gives
  them, for a piece the cut does not expand, as a placeholder from the piece's
  `abridged` row and the record, with no `blocks` unit read
- **Own steps where a walk leaves or rejoins:** a departure gets its own steps
  only where it leaves or rejoins the reference on the window or its context
  ([Departures](cutting.md#departures)). One that jumps over them comes back
  with its ends and its own bp and no `steps`, keeps no node in the cut, and
  costs no read past its line in the departure blocks. graph-walks 6 gave every
  departure within `attachBp` over the window its steps, and fetched a tile for
  each record it lacked a piece of, however far off its ends lay.
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
  ([Cut presets](cutting.md#cut-presets))
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
  reference not loaded by their bp ([Abridged walks](reads.md#abridged-walks))
