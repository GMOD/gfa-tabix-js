# Cut objects

`getCutForRange` returns the cut as objects, for a caller that needs more than
the [tables](api.md#tables) of `getGraphForRange`: each walk's fragments, its
departures from the reference, and the gaps between its fragments.

```ts
const { fragments, departures, gaps, tables } = await walks.getCutForRange(
  'chr6',
  160_614_797,
  160_647_758,
)
```

- `getCutForRange` takes the [options](api.md#options) of `getGraphForRange`,
  with `keep` in place of `haplotypes`. `keep` is a filter on
  `sample#haplotype#contig` that also decides whether the reference walk is
  kept, and `walksKept(haplotypes, refName)` builds the filter `haplotypes` uses
- The cut also has `graph`, `kept`, `placeholders` and the `refName` it cut.
  `tables()` writes the tables `getGraphForRange` returns

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
