# Runs

`getRunsForRange(refName, start, end, opts)` returns a range's walks as runs: no
step is decoded and no node or placeholder made. It takes a range as the other
range methods do, resolving a contig on `referenceSample` and throwing
`MissingSequenceError` for a sequence the file lacks, and the options `signal`,
`haplotypes`, `context`, `contained`, `attachBp`, `byteBudget`,
`departureByteBudget` and `budgetNames`.

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
- It reads the window's tiles and departure blocks and those tiles' `blocks`
  units: no detail unit, and no tile fetched to attach
- It never refuses on a step budget, a node limit or a placeholder limit; it
  refuses a window past `byteBudget` or `departureByteBudget` as a cut does

## What runs hold

- A stretch is consecutive steps on the reference, each starting on the base the
  last ended on, all one way. It never spans a hole. A run is what lies between
  two stretches
- `walks` holds one fragment per walk of the cut's tables, in their order:
  `names`, `starts` and `ends` on the contig, `trimmedBefore` and
  `trimmedAfter`, `offsets` into `stretches`, `leadBp` and `tailBp`, the bp off
  the reference before the first stretch and after the last that a contained
  cut's context keeps, with `leadLoaded` and `tailLoaded`, and the `before` and
  `after` ends
- `stretches` holds `starts` and `bp` on the reference, `against`, `offBefore`,
  the bp of the run before, 0 for a fragment's first, and `loaded`, 0 where that
  run holds a placeholder. `loaded` reports the file's abridging, a piece the
  writer moved out of its tile, not which nodes a host loaded
- In a plain cut a fragment runs from its first step on the window to its last,
  so a walk that goes away and comes back has stretches off the drawn stretch
  and off the reference's nodes. A contained cut splits it there and keeps every
  stretch on the window: a fold of the runs, which knows the reference's nodes
  alone, uses `contained: true`
- Stretches either side of a hole follow on with `offBefore` 0, so a consumer
  takes `holes` out before it calls the gap between them a deletion
- `reference` lists the reference's nodes over its fragments and `drawn` spans
  them; `drawn` is NaN where the reference has no node over the window
- `holes` are the file's holes on the sequence, and `revisited` is always 0
- `departing` lists the departures over the drawn stretch of the walks with no
  fragment, such as one that departs across the whole window, as the cut's
  departure table lists them: `names`, `hapStarts`, `ownBp`, `missingBp`,
  `leavePos`, `leaveReversed`, `rejoinPos` and `rejoinReversed`, so a host can
  tell such a walk from one with no walk there. A walk with a fragment is left
  out, even where it jumps the window elsewhere

## Ends

An end (`WalkEnds`, its kind as an index into `END_KINDS`) comes from the
departure that leaves the fragment's last stretch or rejoins its first, computed
in its piece or read from a departure record between pieces:

| Kind        | Meaning                                                                                                          |
| ----------- | ---------------------------------------------------------------------------------------------------------------- |
| `beyond`    | The walk is next on the reference past the drawn stretch on that side, the same way round                        |
| `elsewhere` | Any other departure: reversed, back on the other side, or with its far end inside the drawn stretch              |
| `edge`      | No departure, and the fragment reaches the drawn stretch's edge, so whatever the walk does next lies off it      |
| `nowhere`   | No departure, and the end lies inside the drawn stretch: the walk has no step on this reference sequence past it |

- For `beyond` and `elsewhere`, `pos` is the departure's far breakpoint,
  `reversed` 1 where that step runs against the reference, `ownBp` the walk's bp
  between the fragment's stretch and that step, a lead or tail included, and
  `missingBp` the part between the contig's W lines
- For `edge` and `nowhere`, `pos` is NaN and `reversed` and `missingBp` are 0.
  `ownBp` is NaN for `edge`, and for `nowhere` the rest of the run the read
  holds past the stretch, all of it off the reference

## Reading a cut back as runs

`fragmentRuns(tables, refName, holes)` reads a cut's tables, from
`getCutForRange(...).tables()` or `getGraphForRange`, back as runs; `holes` are
the walk file header's, which the tables do not state. It is what
`getRunsForRange` returns for a cut that expands nothing, and the oracle the
checks below hold it to. It does not refuse: it copies the tables' `revisited`,
which a caller checks is 0.

## What runs refuse

A file whose header does not state `revisited:i:0` is refused with a
`RunsRefusedError`, since a run off the reference could then hold a node the
reference passes; a host cuts the window instead.

## How runs are checked

`getRunsForRange` equals `fragmentRuns` of the cut for every fixture window,
plain and contained, and over the hosted GRCh38 file
(`scripts/runs-equality.ts`), and `scripts/runs-mutants.ts` breaks the reader on
purpose and holds each break to failing.
