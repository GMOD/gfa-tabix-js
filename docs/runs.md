# Runs

`getRunsForRange` returns a range's walks as stretches of the reference and the
bp between them, without decoding a step or making a node.

```ts
import { WalkFile } from '@gmod/gfa-tabix'

const walks = new WalkFile({ walksUrl })
const runs = await walks.getRunsForRange('chr6', 32_510_000, 32_600_000, {
  contained: true,
})
```

- The runs are the walks of `getCutForRange` over the same range with the same
  `context`, `contained` and `attachBp`, for a cut that expands nothing: the
  reference's walk and those of the `haplotypes` asked for, every one by default
- The options are `haplotypes`, `context`, `contained`, `attachBp`,
  `byteBudget`, `departureByteBudget`, `referenceSample`, `signal` and
  `budgetNames`
- `getRunsForRange` has no step budget, node limit or placeholder limit. It
  throws past `byteBudget` or `departureByteBudget` as a cut does
- A file whose header does not state `revisited:i:0` throws a
  `RunsRefusedError`, since a run off the reference could then hold a node the
  reference passes. Cut the window for such a file

## Returned fields

A stretch is consecutive steps on the reference, each starting on the base where
the last ended, all in one direction. A run is what lies between two stretches.

| Field       | Contents                                                                                 |
| ----------- | ---------------------------------------------------------------------------------------- |
| `walks`     | One fragment per walk, in the order of the cut's tables, with `offsets` into `stretches` |
| `stretches` | `starts` and `bp` on the reference, `against`, `offBefore` and `loaded`                  |
| `reference` | The reference's nodes over its fragments                                                 |
| `drawn`     | The span of those nodes, NaN where the reference has no node over the window             |
| `holes`     | The file's holes on the sequence                                                         |
| `departing` | The departures over the window of walks that have no fragment there                      |
| `revisited` | Always 0                                                                                 |

`walks` columns:

- `names`, and `starts` and `ends` on the contig
- `trimmedBefore` and `trimmedAfter` count the walk row's bp past either end of
  the fragment. They include its steps on the reference, so they do not measure
  sequence off the reference
- `leadBp` and `tailBp` are the bp off the reference before the first stretch
  and after the last that a contained cut's context keeps, with `leadLoaded` and
  `tailLoaded`
- `before` and `after` are the fragment's [ends](#ends)

`stretches` columns:

- `offBefore` is the bp of the run before the stretch, 0 for a fragment's first
- `loaded` is 0 where that run is an abridged piece of the file
- Stretches either side of a hole follow on with `offBefore` 0, so take `holes`
  out before calling the gap between two stretches a deletion

`departing` has the columns of the cut's departure table: `names`, `hapStarts`,
`ownBp`, `missingBp`, `leavePos`, `leaveReversed`, `rejoinPos` and
`rejoinReversed`. It lets a host tell a walk that jumps over the whole window
from a haplotype with no walk there.

In a plain cut a walk that goes away and comes back has stretches off the
window. `contained: true` splits the walk there and keeps every stretch on the
window, which is what a consumer that knows the reference's nodes alone wants.

## Ends

Each fragment end has a kind, an index into `END_KINDS`:

| Kind        | Meaning                                                                                                          |
| ----------- | ---------------------------------------------------------------------------------------------------------------- |
| `beyond`    | The walk is next on the reference past the drawn stretch on that side, the same way round                        |
| `elsewhere` | Any other departure: reversed, back on the other side, or with its far end inside the drawn stretch              |
| `edge`      | No departure, and the fragment reaches the drawn stretch's edge, so whatever the walk does next lies off it      |
| `nowhere`   | No departure, and the end lies inside the drawn stretch: the walk has no step on this reference sequence past it |

- For `beyond` and `elsewhere`, `pos` is the departure's far breakpoint,
  `reversed` is 1 where that step runs against the reference, `ownBp` is the
  walk's bp between the fragment's stretch and that step, and `missingBp` is the
  part between the contig's W lines
- For `edge` and `nowhere`, `pos` is NaN and `reversed` and `missingBp` are 0.
  `ownBp` is NaN for `edge`, and for `nowhere` it is the rest of the run past
  the stretch, all of it off the reference

## Reading a cut back as runs

`fragmentRuns(tables, refName, holes)` converts a cut's tables to runs. `holes`
are the walk file header's, which the tables lack. It returns what
`getRunsForRange` returns for a cut that expands nothing, and it copies the
tables' `revisited`, which the caller checks is 0.

`scripts/runs-equality.ts` checks that `getRunsForRange` equals `fragmentRuns`
of the cut over every fixture window and over the hosted GRCh38 file.
