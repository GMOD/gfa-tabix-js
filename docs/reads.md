# Reads, budgets and abridged walks

## Budgets

- Over the byte budget, the error's `bytes` is what the window would fetch, so a
  viewer with a fixed window can offer to load it anyway with `byteBudget`
  raised to `bytes`. `zoomHelps` is false where no narrower window fits, and
  `fitsBp` is otherwise narrower than the window refused, also where the tiles
  fetched to attach are what takes the read over
- The cut checks the node limit as it decodes each tile, against the nodes it
  keeps whatever else it reads, and stops decoding at the first tile over it
- A placeholder is a node in the cut's tables that `nodeLimit` does not count
  ([Abridged walks](#abridged-walks)), so where the file abridges walks, a cut
  within its node limit can return far more table nodes than the limit: at a
  `nodeLimit` of 30,000, 20 kb of HLA class II in the HPRC file returns 214,636,
  of which 209,147 are placeholders. `placeholderLimit` bounds them. A cut that
  keeps more throws a `RegionTooLargeError` whose `limit` is the placeholder
  limit, when it builds the piece that takes it over, before it cuts. A reader
  that lays out every node of the tables passes both limits. Expanding a walk
  turns its placeholders into nodes, which `nodeLimit` counts
- A window past `departureByteBudget` throws a `DepartureReadError`, which is
  not a `RegionTooLargeError`: a departure's record is filed over its whole
  stretch, so zooming in does not lower that read. A collapsed repeat that
  haplotypes hop between the copies of is the usual cause
- The byte budgets bound what a cut fetches: the unit table gives the exact
  compressed bytes of the window's tiles and of its departure blocks before
  either is fetched, and the cut refuses either past its budget. The tiles
  fetched for records' own steps count against `byteBudget` with the window's,
  and a cut they would take past it throws the same `RegionTooLargeError` before
  fetching any of them, so a window attaches the own steps of every record that
  [gets them](cutting.md#departures) or of none. The `blocks` units a cut reads
  count against `byteBudget` too. The detail records `expand` names count only
  against `stepBudget`, once read; the writer moves paths out of any tile whose
  rows would compress past its header's `abridgebytes`
- A narrow window reads whole tiles, so a cut reads several times the steps it
  keeps. The rows a cut decodes from the window's tiles, and from the detail
  records `expand` names once they are decoded, may hold up to 8 times
  `stepBudget` steps (`READ_STEPS_PER_KEPT`), a step along the reference
  counting a quarter; past that the cut throws a `RegionTooLargeError` before it
  cuts. The rows of tiles fetched for records' own steps are not counted there,
  and the own steps they give count against `stepBudget` with the rest
- Before it decodes a row, the cut counts the steps it will surely keep of the
  window's tiles whose extent lies strictly inside the window and its context.
  When they are over `stepBudget` and the rows and the unit table show that no
  refusal the cut makes before it counts the steps it keeps can come, it throws
  there the `RegionTooLargeError` the finished cut would, the same in every
  field but `fitsBp` and the span its message names. So such a window costs its
  fetch, not its cut. The refusals it rules out are those of a malformed row for
  a route its tile lacks, a route without an entry or exit inside the row, fewer
  steps than its routes hold, or an `abridged` row's numbers; the node limit,
  which no more nodes than the tiles' `ref` rows and routes name can reach,
  where no detail record `expand` reads or tile fetched to attach can bring
  more; the placeholder limit, where no tile abridges a walk the cut keeps and
  does not expand and no tile is fetched to attach; the byte budget, over the
  `blocks` units the cut reads and every tile its records could fetch to attach;
  and the 8 times `stepBudget` the rows may hold, over the most the rows and the
  pieces built from them can. Where it can't rule one out, the finished cut
  refuses the window. A file malformed in any other way, such as a `ref` row
  whose passes or a `walk` row whose `nsteps` disagree with its steps, or a
  `blocks` or detail unit that disagrees with its rows, may be refused here as
  too large before the cut would find the fault, which a window narrow enough to
  cut then reports. The count takes each `ref` and `walk` row of a walk `keep`
  takes at its steps, less a walk row's steps off the reference before its first
  step on it and after its last. Every step a row places on the reference lies
  in its tile's extent, so a cut over that extent keeps the row from its first
  such step to its last and every step between, contained or not; the steps off
  the reference at a row's ends are kept only where a neighbouring piece reaches
  the window. An `abridged` row counts nothing there, since the cut builds its
  piece with a placeholder per run off the reference, in fewer steps than its
  `nsteps`, and neither do departures' own steps. So the count never passes the
  steps the cut's fragments keep. Its `fitsBp` is worked out tile by tile from
  the start of the window's context, each counted tile holding its count evenly
  over its extent, a tile not counted holding the counted tiles' steps per bp
  and a stretch no tile covers none: where those steps reach `stepBudget` is
  where the context of a window that fits ends, so it is narrower than the
  window refused, as the finished cut's is. Both are estimates, the finished
  cut's the window scaled by the steps it kept
- Where the rows' count fits `stepBudget` and the `nsteps` of the window's
  `abridged` rows could take it over, the cut reads the window's `blocks` units
  before it decodes, the same units it reads after decoding when it goes on, and
  adds what it surely keeps of each abridged piece it does not expand, in every
  tile of the window: the steps on the reference under each block of its
  `blocks` row that lies wholly inside the window and its context, and one
  placeholder for each run off the reference between two such blocks. The cut
  builds those blocks whole and keeps every step it places on that stretch with
  every step between two of them, so this count too never passes the steps the
  fragments keep. Over `stepBudget`, and with no earlier refusal to come, it
  refuses there: 90 kb of HLA class II in the HPRC file, where 332 of 464 walks
  are abridged, is refused in 0.7 s where the finished cut took 7 s to refuse
  it. A piece without a `blocks` row, or with one that does not add up to its
  `abridged` row or cover whole steps of the reference, leaves the window to the
  finished cut, which reports it. A piece the cut expands counts nothing, since
  the steps off the reference at its ends are in its detail record. For
  `fitsBp`, a tile at the edge of the window and its context holds its pieces'
  count and, as an estimate, its rows' count scaled to the part of its extent
  the cut draws

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
  - those holding the own steps of a record that
    [gets them](cutting.md#departures)
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
    ([Gaps](cutting.md#gaps))
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
- The tables keep the steps of a walk's unloaded pieces as
  `walks.abridgedSteps`, each piece on the first walk holding it, and of those
  no walk holds as `departures.abridgedSteps`, each on the first departure
  holding it, NaN for a walk or departure with none. A walk's unloaded steps are
  the two together
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
  the bytes read to `getBytesForRange` (`scripts/blocks-equality.ts`)

## What a cut reads

A cut reads these units, each once:

- **The window's tiles and departure blocks:** the `tile` and `departures` units
  whose `start end` overlap the window and its context
- **`blocks` units:** the unit of each of those tiles that abridges a walk the
  cut keeps and does not expand. It builds from it the pieces it draws
  ([Abridged walks](#abridged-walks)), and checks every row it reads against its
  `abridged` row's `hapBp`
- **Detail units:** with `expand`, those of the pieces it builds of the walks
  named, among the window's tiles and then among the tiles fetched to attach.
  Those whose bytes touch in the file are fetched by one read, up to 2 MB a
  read, and 16 reads at a time: expanding every walk over 90 kb of HLA class II
  in the HPRC file reads 519 detail units, back to back in two stretches of the
  file, by 4 reads, the same bytes
- **Tiles fetched to attach:** for each record within `attachBp` that leaves or
  rejoins on the window or its context ([Departures](cutting.md#departures)) and
  whose leave piece or rejoin piece the window's tiles lack, the tiles its
  `ownTiles` names. No `blocks` unit is read for them: an unexpanded piece's own
  run comes from its `abridged` row and the record

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

`walks.getBytesForRange(refName, start, end, { context, attachBp })` adds these
up before a cut is made:

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
- **With `expand`:** the cut reads no more than `getBytesForRange` but for its
  detail units, and no less than `getBytesForRange` without the `blocks` units.
  Which detail units it reads is in the tiles' `abridged` rows, which
  `getBytesForRange` does not read

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
  that reads the file, `getReferenceSequences`, `getUnits` and
  `getReferenceNodes` among them, waits on this check
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
