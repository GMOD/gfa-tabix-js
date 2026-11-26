# Counts

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

- `getCountsForRange` counts every walk of the file over a region, where a cut
  draws the few it keeps: for each node and link the walks step, `haplotypes`,
  how many haplotypes step it, and `passes`, how many steps all of them take on
  it. A haplotype is what the header lists a path under, `HG002#1`, so a
  haplotype stepping a node three times counts once in the first and three times
  in the second, and its contigs count as one
- A link counts however a walk takes it: `>3>4` and `<4<3` are one link, with
  each haplotype once
- It takes a cut's `context`, `contained`, `attachBp`, budgets and `expand`, and
  no `keep`. Every node and link of a cut of the same region with the same
  options, whichever walks it keeps, has a count; counts that expand every walk
  cover a cut that expands any
- `nodes.ids` ascend. `links.from` and `links.to` are handles, `2 × id + 1`
  where reversed, written the way round a cut writes its links, and ascend by
  `from`, then `to`
- `haplotypes.file` is how many haplotypes the file holds: the reference's and
  each one a `#haplotype` line of the header lists, which are the haplotypes
  with rows in the file. A graph can hold more, since a haplotype whose every
  line steps no reference node has no row
- `haplotypes.window` is how many have a piece the counts read or a departure
  over the window and its context. It includes a haplotype that crosses the
  window with no step in it, and one whose record there holds sequence the graph
  lacks

What it counts:

- Each piece a cut builds, whole: those over the drawn stretch, those between
  two such pieces of a line, and, contained with context, the one before such a
  piece. A cut keeps a walk from its first step on the drawn stretch to its
  last; the counts also take the rest of the piece, such as a contig's run off
  the reference before its first step on the reference, which the cut reads as a
  `read` gap
- Each piece of a window's tile holding the own steps of a departure a cut
  attaches, and of a tile fetched to attach them, whole
- An abridged piece by its detail unit where `expand` names its walk, and else
  by its `blocks` row: its steps on the reference exactly, its runs off it as
  bp, which leave its haplotype's counts off the reference lower bounds

## Exact and lower bound

`haplotypeReasons[k]` and `passReasons[k]` are masks, 0 where the count is
exact, and else the `COUNT_REASONS` bits it may fall short by, which
`countReasons` names:

- **Exact, in a file whose header states `revisited:i:0`:** a node the reference
  steps on the drawn stretch, and a link between two such nodes. Every step on
  the node lies at its pass, in the tile holding it, and the counts read every
  piece with a step there, an abridged one's steps on the reference from its
  `blocks` row, so the count is of every line the file holds
- **What walk format 8 promises of a file that states `revisited:i:0`:** no
  route interior holds a reference node, so the tiles a window's read selects
  hold a row of every line of the file that steps a reference node over the
  window, and a `ref` row's `passes` counts the `ref` rows of the whole file
  ([What a reader relies on](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/walks-format.md#what-a-reader-relies-on),
  [ADR 0016](https://github.com/GMOD/gfa-to-tabix/blob/main/docs/adr/0016-format-8-is-frozen.md)).
  The exact counts rest on those two statements
- **A file that states `revisited` above 0, or none:** the format promises
  nothing about where such a file files a step on a node the reference passes
  once, so the counts mark no node or link of it exact by its `passes`: each
  takes `elsewhere`, and is a lower bound. gfa-to-tabix keeps the rule on every
  file it writes, but a file names no writer, so the counts go by the header:
  `(await walks.getHeader()).revisited` is 0 where a count can be exact by its
  `passes`
- **A haplotype count of `haplotypes.file`** is exact in any file, since the
  file holds no other haplotype. Its pass count is exact only by the rule above
- **A file that states `revisited:i:0` and holds a reference node in a route
  interior** is not walk format 8. The counts read a window's tiles and cannot
  tell, so a count they mark exact there falls short by each line filed away
  from the node's pass. gfa-to-tabix's `verify7.py` and `spec8.py` read the
  whole file and refuse it
- **An exact count is of the lines the file holds.** Whether the file holds
  every line of its graph that steps the reference is checked against the GFA,
  by gfa-to-tabix's checks, and cannot be read off the file
- **`elsewhere`:** any other node, and a link with an end there. A step off the
  reference is filed beside the step on it before it, which can lie anywhere: in
  a stretch of a haplotype that leaves and rejoins the reference away from the
  window, or in a W line with no step on the reference, which the file holds no
  row of (the header's `leftout` counts them). No row over the window points to
  either. A node on the reference off the drawn stretch is counted only from the
  pieces read. With no other reason, the count is exact over every walk the file
  places over the window: each piece the counts read and each departure over the
  window, every step loaded. Of the graph's haplotypes it is still a lower bound
- **`unloaded`:** a haplotype over the window holds a run off the reference the
  counts did not load, an abridged piece not expanded or one the step budget
  stopped. `haplotypes.unloaded` names them
- **`unattached`:** a haplotype over the window holds own sequence of a
  departure over it the counts did not read: one past `attachBp`, one that
  neither leaves nor rejoins on the window and its context, or own sequence in a
  W line with no step on the reference. `haplotypes.unattached` names them
- **`revisited`:** the reference passes the node more than once, so a step on it
  can lie at a pass, or off the reference in a tile, the counts did not read
- A haplotype count takes `unloaded` or `unattached` only where a haplotype they
  name does not step the node or link; a pass count takes them wherever any is
  named
- With mask 0, `haplotypes[k]` of the `haplotypes.file` haplotypes the file
  holds step it, each in a line the file holds, and no other of them does; with
  any other mask, at least that many step it and no share follows. Expanding
  every walk, `expand: () => true`, leaves no `unloaded` wherever the step
  budget holds
- A haplotype a cut shows can be in a node's count although its drawn fragment
  does not step the node, since a cut keeps a walk only from its first step on
  the drawn stretch and the counts take the whole piece. This never happens on a
  node marked exact: every step there lies on the drawn stretch

## What the numbers support

- **Mask 0:** "N of the M haplotypes the file holds step this node", with M
  `haplotypes.file`. M is the file's and can be fewer than the graph's.
  `haplotypes.window` is not an M, since which haplotypes it counts follows from
  the tiles the window reads
- **Any other mask:** "at least N haplotypes step this node" and nothing more:
  never a ratio, a percentage or a complement such as "M − N do not", never a
  sum over nodes or a comparison of two such counts, and drawn as a different
  kind of mark from an exact count
- **`passes`** is steps, not haplotypes: "stepped P times" with mask 0, and
  "stepped at least P times" with any other
- **A count is of the whole file.** It is not explained by the walks a cut
  shows, and should not be presented as if it were

## What counts cost

- **Reads:** the units a cut that keeps every walk reads with the same options.
  Without `expand`, `getBytesForRange`'s: the window's tiles, their `blocks`
  units and departure blocks in one read, then the tiles fetched to attach. With
  it, the window's tiles and departure blocks, then the detail units of the
  pieces it expands and the `blocks` units of the pieces it does not, all in one
  read of every run of touching units, then the tiles fetched to attach, and the
  detail units it expands of their pieces
- **Work:** each walk row it counts decoded once, with no step kept: a step
  along the reference is an array read and two bit sets, and a route's nodes and
  links are looked up once per tile. It grows with the steps of the pieces
  counted, as a cut of every walk does, and builds no fragment, table or graph
- **Memory:** per node and link a bitset of the file's haplotypes, a 32-bit word
  per 32 haplotypes, and a pass count
- **Budgets:** `byteBudget` over the tiles, `blocks` units and tiles fetched to
  attach, and `departureByteBudget`, as a cut: past either the counts throw
  before reading. The window's walk rows past `READ_STEPS_PER_KEPT` ×
  `stepBudget` steps, a step along the reference a quarter, throw a
  `RegionTooLargeError`. Expanding adds each detail unit's steps to them and
  stops at the first unit that would pass the budget, leaving that piece and
  every later one unloaded, so their haplotypes are `unloaded`
- **A refusal is worded for counts.** Past `byteBudget` or the step budget the
  `RegionTooLargeError` has a cut's fields and a message about counting, such as
  "Too many bytes to read to count this window: … Count a window of about 40bp
  or less, or leave this one uncounted"
- A row of a haplotype the header does not list throws a `WalkFormatError`, as
  do a missing `blocks` row, and a `blocks` row, a detail row or a line's next
  piece that does not add up to its `abridged` row
