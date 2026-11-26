# Large genomic regions and abridged walks

A cut reads every selected haplotype's walk through the region, so the bytes to
download and the nodes to draw grow with the region's width and the number of
haplotypes. A cut over one of its limits throws. Where a stretch of the graph is
too large for any window, the file stores its costliest paths as abridged walks,
which a cut reads in full when `expand` names them.

## Large genomic regions

| Limit                 | Counts                                                           | Error                 |
| --------------------- | ---------------------------------------------------------------- | --------------------- |
| `byteBudget`          | Compressed bytes of the tiles and `blocks` units the cut fetches | `RegionTooLargeError` |
| `stepBudget`          | Steps the cut keeps, the own steps of departures included        | `RegionTooLargeError` |
| `nodeLimit`           | Nodes the cut keeps, placeholders not among them                 | `RegionTooLargeError` |
| `placeholderLimit`    | Placeholders the cut keeps                                       | `RegionTooLargeError` |
| `departureByteBudget` | Compressed bytes of the departure blocks over the window         | `DepartureReadError`  |

- A `RegionTooLargeError` has `limit`, the limit the cut passed, `fitsBp`, about
  how wide a window would fit, and `zoomHelps`, which is false where no narrower
  window fits
- Over `byteBudget` the error also has `bytes`, what the window would fetch. The
  cut knows the bytes from the file's unit table and throws before fetching, so
  a viewer with a fixed window can offer to load it anyway with `byteBudget`
  raised to `bytes`
- A window that is surely over `stepBudget` is refused before the cut decodes
  it, so refusing it takes the time of its fetch
- A narrow window reads whole tiles. A cut whose rows hold more than 8 times
  `stepBudget` steps throws a `RegionTooLargeError` too
- A `DepartureReadError` is a separate class because zooming in does not help: a
  departure's record is filed over its whole stretch. The usual cause is a
  collapsed repeat whose copies the haplotypes jump between
- `nodeLimit` does not count placeholders, so a cut of a file with abridged
  walks can return far more table nodes than the limit. At a `nodeLimit` of
  30,000, 20 kb of HLA class II in the HPRC file returns 214,636 nodes, of which
  209,147 are placeholders. A reader that lays out every node of the tables
  passes both limits

## Abridged walks

In a tile over its budget, gfa-to-tabix moves the rows of the costliest paths
out of the tile into detail records. A satellite array in a pggb graph is the
usual case.

- By default a cut returns an abridged piece with its steps on the reference and
  one placeholder step for each run off the reference. The placeholder has the
  run's bp and no node
- `expand` naming the walk reads the piece from its detail record. A cut that
  expands every walk it keeps returns what a file written without abridging
  returns, within its budgets
- `isPlaceholder(id)` identifies a placeholder step. Its node, in `graph.nodes`
  and in the tables' `nodes` with `placeholders` set, has the run's bp as its
  length and the name `walk:start-end`. The name is the same in every cut, and
  the id is valid in its cut alone
- Two placeholders of equal bp may stand for different sequence, and a
  placeholder can stand for nodes of the reference where the reference passes a
  node more than once. Expanding the walk is the way to tell
- A fragment's `abridged`, and a departure's, list the pieces the cut did not
  load, each with its `steps` and the `detail` unit that holds it. The tables
  have the step counts as `walks.abridgedSteps` and `departures.abridgedSteps`
- `abridgedWalks(tables)` lists each walk with a run off the reference that the
  cut did not load, as `{ name, bp }`, and `abridgedText(walk, label)` writes it
  as a sentence:

```ts
abridgedWalks(tables).map(walk => abridgedText(walk))
// ownends.R over R#0#c:0-30 with no context, among others:
// 'S#1#c: 70 kb off the reference not loaded'
// 'V#1#c: 20 bp off the reference not loaded'
```

`pnpm blocks-equality <walks.gz or URL>` checks a walk file: over random windows
on its abridged tiles, the cut that expands nothing must agree with the cut that
expands every walk.

## What a cut reads

A cut reads these units of the file, each once:

- **The window's tiles and departure blocks,** which are the units that overlap
  the window and its context
- **`blocks` units,** one for each of those tiles that abridges a walk the cut
  keeps and does not expand
- **Detail units,** with `expand`, for the pieces of the walks named
- **Tiles fetched to attach,** for each departure that leaves or rejoins on the
  window or its context ([Departures](api.md#departures)) and has its own steps
  in a tile outside the window's

A departure that only lies over the window, with both ends elsewhere, adds its
line in the departure blocks and no tile.

`getBytesForRange(refName, start, end, { context, attachBp })` adds these up
before a cut is made:

- The result is exact for a cut that keeps every walk and expands none
- With `keep` or `haplotypes` it is an upper bound, since fewer walks need fewer
  `blocks` units and attached tiles
- With `expand` the cut also reads its detail units, which `getBytesForRange`
  does not count
- It reads the departure blocks, so it throws a `DepartureReadError` past
  `departureByteBudget`

### Serving a walk file

- The reader fetches units by single-range HTTP requests, merging ranges that
  touch
- At open it checks the directory against the file's length, which over HTTP
  comes from `Content-Range`. A server must expose that header to cross-origin
  readers, or the reader skips the check
- A unit or table that does not match the directory throws a `WalkFormatError`
- gfa-to-tabix's
  [walks-format.md](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/walks-format.md)
  describes the file
