# Runs

`getRunsForRange(refName, start, end, opts)` returns the walks `getCutForRange`
returns for the same range and options, as runs: no step is decoded and no node
or placeholder made. It takes a range as the other range methods do, resolving a
contig on `referenceSample` and throwing `MissingSequenceError` for a sequence
the file lacks, and the options `signal`, `keep` or `haplotypes`, `context`,
`contained`, `byteBudget`, `departureByteBudget` and `budgetNames`. It reads the
window's tiles and departure blocks and those tiles' `blocks` units: no detail
unit, and no tile fetched to attach.

```ts
import { WalkFile } from '@gmod/gfa-tabix'

const walks = new WalkFile({ walksUrl })
const runs = await walks.getRunsForRange('chr6', 32_510_000, 32_600_000)
```

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
  run holds a piece the read did not load
- `reference` lists the reference's nodes over its fragments, `drawn` spans
  them, and `holes` are the file's holes on the sequence
- `departing` lists the departures over the drawn stretch of the walks with no
  fragment, such as one that departs across the whole window, as the cut's
  departure table lists them: `names`, `hapStarts`, `ownBp`, `missingBp`,
  `leavePos`, `leaveReversed`, `rejoinPos` and `rejoinReversed`, so a host can
  tell such a walk from one with no walk there

## Ends

An end (`WalkEnds`, its kind by `END_KINDS`) comes from the departure that
leaves the fragment's last stretch or rejoins its first, computed in its piece
or read from a departure record between pieces:

- `beyond`: the walk is next on the reference past the drawn stretch on that
  side, the same way round
- `elsewhere`: it is next on the reference reversed, or back on the other side
- `edge`: no departure, and the fragment reaches the drawn stretch's edge, so
  whatever the walk does next lies off it
- `ends`: no departure inside the drawn stretch, so the walk has no more steps
  on this reference sequence

`pos`, `reversed`, `ownBp` and `missingBp` are the departure's.

## What runs refuse

A file whose header does not state `revisited:i:0` is refused with a
`RunsRefusedError`, since a run off the reference could then hold a node the
reference passes; a host cuts the window instead. A window past the byte budgets
is refused as a cut refuses it.

## How runs are checked

`getRunsForRange` equals `getCutForRange`'s tables read back as runs for every
fixture window, plain and contained, and over the hosted GRCh38 file
(`scripts/runs-equality.ts`), and `scripts/runs-mutants.ts` breaks the reader on
purpose and holds each break to failing.
