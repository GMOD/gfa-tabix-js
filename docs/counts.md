# Counts

`getCountsForRange` counts how many of the file's haplotypes step each node and
link over a region. It counts every walk of the file, where a cut returns the
walks it was asked to keep.

```ts
import { CUT_PRESETS, countReasons } from '@gmod/gfa-tabix'

const counts = await walks.getCountsForRange(
  'chr6',
  160_614_797,
  160_647_758,
  CUT_PRESETS.browser,
)
const k = counts.nodes.ids.indexOf(id)
counts.nodes.haplotypes[k] // haplotypes stepping the node
counts.nodes.passes[k] // their steps on it
countReasons(counts.nodes.haplotypeReasons[k]) // [] where exact
counts.haplotypes.window // haplotypes over the window
```

## What counts return

| Field                    | Meaning                                                                          |
| ------------------------ | -------------------------------------------------------------------------------- |
| `nodes.ids`              | Node ids, ascending                                                              |
| `links.from`, `links.to` | Link ends as handles, `2 × id + 1` where reversed, ascending by `from` then `to` |
| `haplotypes`             | How many haplotypes step the node or link                                        |
| `passes`                 | How many steps all of them take on it                                            |
| `haplotypeReasons`       | 0 where `haplotypes` is exact, else why it is a lower bound                      |
| `passReasons`            | The same for `passes`                                                            |
| `haplotypes.file`        | How many haplotypes the file has rows for                                        |
| `haplotypes.window`      | How many have a piece or a departure over the window and its context             |
| `haplotypes.unloaded`    | The haplotypes with a run off the reference the counts did not load              |
| `haplotypes.unattached`  | The haplotypes with a departure whose own sequence the counts did not read       |

- A haplotype is what the header lists a path under, such as `HG002#1`. One that
  steps a node three times adds one to `haplotypes` and three to `passes`, and
  its contigs count as one haplotype
- A link counts whichever way a walk takes it: `>3>4` and `<4<3` are one link
- `getCountsForRange` takes a cut's `context`, `contained`, `attachBp`, budgets
  and `expand`, and has no `keep`. Every node and link of a cut of the same
  region with the same options has a count
- `haplotypes.file` can be fewer than the graph's haplotypes, since a haplotype
  that steps no reference node has no row in the file

## Exact counts and lower bounds

A count is exact where its reason mask is 0. Otherwise it is a lower bound, and
`countReasons(mask)` names the `COUNT_REASONS` bits:

| Reason       | Meaning                                                                                                                |
| ------------ | ---------------------------------------------------------------------------------------------------------------------- |
| `elsewhere`  | The node is off the reference, or off the window, so a step on it can be in a part of the file the counts did not read |
| `unloaded`   | A haplotype over the window has a run off the reference the counts did not load: an abridged piece not expanded        |
| `unattached` | A haplotype over the window has a departure whose own sequence the counts did not read, such as one past `attachBp`    |
| `revisited`  | The reference passes the node more than once, so a step on it can be at a pass the counts did not read                 |

- A count is exact for a node the reference steps on the window and its context,
  and for a link between two such nodes, in a file whose header states
  `revisited:i:0`. `(await walks.getHeader()).revisited` is 0 for such a file
- In a file that states `revisited` above 0, or does not state it, every count
  has `elsewhere` and is a lower bound
- A haplotype count equal to `haplotypes.file` is exact in any file
- `expand: () => true` removes `unloaded` wherever the step budget holds
- An exact count is of the lines the file has. gfa-to-tabix's checks compare the
  file against its GFA

## What to say with a count

- **Mask 0:** "N of the M haplotypes in the file step this node", with M from
  `haplotypes.file`. `haplotypes.window` is not an M, since it depends on which
  tiles the window reads
- **Any other mask:** "at least N haplotypes step this node". A lower bound
  supports no ratio, percentage, complement, sum over nodes or comparison with
  another count, so draw it as a different kind of mark from an exact count
- **`passes`** counts steps: "stepped P times" with mask 0, and "stepped at
  least P times" otherwise
- A count is of the whole file, so the walks a cut shows may not add up to it. A
  haplotype can be in a node's count although its drawn fragment does not step
  the node, though never on a node marked exact

## Reads and budgets

- Counts read the units a cut that keeps every walk reads with the same options
  ([What a cut reads](reads.md#what-a-cut-reads))
- Counts decode each walk row once and build no fragment, table or graph. Memory
  is one bit per haplotype for each node and link, and a pass count
- `byteBudget` and `departureByteBudget` apply as in a cut, and the counts throw
  before reading. Walk rows past 8 times `stepBudget` steps throw a
  `RegionTooLargeError`, whose message is worded for counting
- With `expand`, counts stop expanding at the first detail unit that would pass
  the step budget, and the haplotypes of the pieces left are `unloaded`
- A file whose rows disagree with its header or its `abridged` rows throws a
  `WalkFormatError`
