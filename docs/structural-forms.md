# Structural forms from walks

- **In:** a `gfa-to-tabix --walks` file and a reference window
- **Out:** the window's forms: groups of haplotypes that make the same
  structural changes, each change with its kind, exact REF span, ALT length and
  carriers
- **Where:** `walkDepartures.ts` and `structuralForms.ts`. The graphgenomeviewer
  plugin (GPL) and jb2hubs call them
- **Measured on:** the published `hprc-v2.1-mc-grch38-f5.GRCh38.walks.bed.gz`,
  Node 24 on a laptop, 462 haplotypes once GRCh38 and CHM13 are left out

## Sites come from the walks

Each haplotype's departures from the reference path R give the sites. The
alternatives fall short:

| Source                               | Verdict                                                                                                                                                                       |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| bubbles file                         | Top-level bubbles only, with no per-haplotype allele; an rGFA's comes from `gfatools` on the SV graph, whose segment ids differ from the base-level graph the walk file holds |
| MinigraphBubbleAdapter / GetSubgraph | One feature per bubble with no genotypes                                                                                                                                      |
| **walks against R**                  | **Chosen**                                                                                                                                                                    |

`walkDepartures` reads each haplotype off the cut's departures and its
fragments' ends, in chains:

- A chain is a fragment with the departures leaving its last step, rejoining its
  first or lying inside it, fragments that abut on the contig, or a departure
  touching no fragment, such as a contig that jumps over the window or rejoins
  far outside it, which the window's chunks never held
- Each chain reads in its own orientation, since a graph may store a contig back
  to front. One that walks R only against it reads the other way round. One that
  walks R both ways reads the other way round only where it leaves R behind
  where it entered, measured in bp of R from its first step on R to its last. So
  an inversion that fills most of the window, between flanks walked forward,
  reads as an inversion (chr16:75,204,000-75,225,000 of the HPRC v2.1 GRCh38
  graph: 439 of 463 haplotypes, 16.5 kb). The bp walked each way decide only for
  a chain that leaves R where it entered
- A chain lying wholly inside an inverted block walks R only against it, as a
  contig stored back to front does, and the two cannot be told apart: it reads
  as reference there. A window that reaches either flank by one step tells them
  apart
- Between a departure's rejoin and the next leave the path follows R. A rejoin
  that is not where the last stretch along R ended closes a departure. Its REF
  span runs from that end to the rejoin, and the path's own bp between, sequence
  the graph lacks between a contig's W lines included, give its ALT length
- A stretch of R walked against it inside a forward chain counts as ALT and as
  inverted bp (AMY1: 317 haplotypes)
- A step back along R is a tandem copy. The departure is the zero-width
  junction, and the bp walked again count as ALT, as `repeatedBp`, with their
  span as `copied`
- R that two chains of one haplotype reach, such as a contig's second visit to
  the window or a second contig over the same locus, is no call there
- An abridged piece the cut did not load reads as its loaded one does. Its
  `blocks` row states its steps on R, each with its strand, and the bp of its
  own between them, so every departure inside it is exact in place and size
  (gfa-to-tabix's ADR 0012). Before that row, a reader had the piece's first and
  last steps on R and its bp alone, and netted the bp between them: a piece that
  turned back or jumped read as a block over R it never stepped on, and one
  stepping past both end steps was absent there

Walks give exact lengths. Each haplotype's summed departures across the truth
table's span match its assembly in 2,779 of 2,780 scored cases, median gap 0 bp;
the HPRC sidecar's equivalent is 99% with a median of 36 bp. The one miss,
NA19909#1 at AMY1 (walk −94,529 bp, assembly −101 bp), came from a second contig
over the AMY copies whose departures added to the first's; R that two walks of
one haplotype cover is now no call.

## Grouping rules

| Rule            | Decision                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `stateKey`      | Size change to 2 significant figures and the kind of change: one LPA copy number differs by a few bp (22,201 vs 22,204). A deletion keys as `-` and an insertion as `+` with its size, a replacement as `r-` or `r+`, a duplication as `d`, so a state is one kind of change and its name holds for every carrier: at chr22:37,509,713, 206 haplotypes leave out 119 bp and 4 walk 1 bp of their own in place of 120, which had shared `-120` and the name `del 119bp`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Detection floor | 50 bp. At a site no wider than 1 kb it applies to each haplotype's net change across the site, its plain pieces summed, since walks split a repeat's change into pieces: at rand29's VNTR the walk net matches the assembly for 408 of 408 haplotypes, the pieces that pass 50 bp alone for 235. A wider site keeps its changes apart: summed, the haplotypes without UGT2B17's 117 kb deletion added up small indels into 7 forms of −50 to −56 bp                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `minCarriers`   | max(1, round(n / 100)): 5 of 462                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `changeKind`    | On the two sides a mark's departures state: the bp of R they leave out and the bp of the haplotype's own they run. A deletion where its own side is 0 bp, an insertion where R's is, otherwise a replacement, at any ratio and with no size threshold. A mark is made only of departures that change R's length, walk it against it or walk it again, so a SNP beside a change neither widens it, names it, nor joins it to another: C4's common change is a 6,367 bp deletion beside six SNPs, drawn over the 6,367 bp it leaves out. R walked between two joined departures is on neither side                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Join            | Departures of one haplotype less than 50 bp of R apart join, keeping their net, once same-length departures are left out. One under 50 bp on both sides, which no mark states alone, joins only another such, so a small indel beside a change neither names nor widens it, while a repeat's change the walks split into small pieces still adds up. The test stays at 50 bp at every floor, so a name does not change with the zoom. At CYP2D6, 38 haplotypes run 12,156 and 1,564 bp of their own one base apart and leave out 1 bp of R 7 bp on: joined whole they read `repl 1bp->13.7kb`, and with the 1 bp apart `ins 13.7kb`. Over the 638 measured windows 493 of 21,606 replacement haplotype-marks are still deletions and insertions with no bp on one side joined, mostly small pieces that add up. The join replaces the sidecar's split-event merging: at GSTM1 walks give 306 carriers of an exact −18,445 bp deletion and 156 with no change, where the sidecar had 6 forms, two of them "shift" forms holding 82 haplotypes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Nesting         | Marks that `minCarriers` haplotypes share, two marks shared where each covers at least half of the other on R, are common, and overlapping common marks form one site, whose span is theirs; a haplotype's state there is only its own marks, so a parent deletion and a nested child never both count. UGT2B17's forms held 230/229/3 haplotypes from walks and 232/229/1 from the sidecar under the earlier union of every mark, not remeasured. A rarer mark states itself at each site it overlaps without widening it: at chr22:20.3-20.4 Mb one HG01175#1 departure over 108.7 kb of R had merged two sites of 461 and 460 haplotypes into one of 137 states. A summed change reaching past its site does not widen it either: at chr22:50.78 Mb a widened site held 1 kb insertions of 3 haplotypes whose state was read over the narrower span                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Place           | At a site wider than `SUMMED_SITE_BP` each change keys on its place as well as its size, `+14000@42125923`. A place is common spans every two of which cover half of each other, taken most carried span first, then by start and end: a common span joins the place with the most carriers all of whose spans it half-covers, then the earliest by start and end, or starts its own, and a rarer mark's place is its own span, so no one haplotype names a place or joins two, and the order the haplotypes come in changes no place. Three +5 kb replacements from one start over 10, 6 and 18 kb, two carriers each, had made `{a,b,c,d}` and `{e,f}` in one order and `{a,b,e,f}` and `{c,d}` in another while equal carriers and starts were left in input order. A place keys as its first span's start, or as that span where a place of the same size and kind starts there too: at LCR22A, 42 haplotypes insert 61 bp at the point 18,651,717, `+61@18651717`, and 21 replace 24 bp there with 85, `r+61@18651717`, where one key had held all 63. Insertions share a place only at one point. CYP2D6's 13.7 kb insertion has place 42,125,923 in 38 haplotypes and 42,134,095 in 4 (MC chr22 42.10-42.15 Mb, 1 kb floor); keyed on size alone they were one state, drawn at the first place for all 42. Places were once the first start among all the marks the half-cover relation chained together, so one haplotype's +1.6 kb replacement named the 38's place 42,125,912, and 10 haplotypes deleting 10-20 kb and 10 deleting 15.5-25.5 kb became one state of 21 when one more deleted 12.75-22.75 kb, covering half of each. A summed site keys on its net change alone, so a repeat's change keeps one state wherever the walks put it |
| Rare states     | A state under `minCarriers` folds into the common state nearest in net change among those it agrees with: change by change at a site keyed on place, each pair covering half of each other, net and largest change at a summed site, each pair under `agrees` (both under the floor, or apart by under the floor, a tenth of the larger, or the last figure the key gives the change's size). A change no common state agrees with is `~`, rarer. R's structure and a no call stay themselves. `foldRareForms: false` keeps every state, as a scorer wants every member's own sizes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Rare forms      | A haplotype whose form has under `minCarriers` joins the largest common form whose state at every site is its own or one its changes agree with, or `other`, the last form, key `*`, which states nothing of its members. Without it LPA at 50 bp has 48 single-haplotype rows                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Edge sites      | A site that most haplotypes in the window have no call at is `edge` and left out of the key. A haplotype's call reads its coverage over the site's part inside the window, so a departure followed past the edge does not make the others no-calls                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| No calls        | A haplotype whose walks leave a gap in a site's part inside the window has no call there, even with a mark at the site, as a walk running off R at AMY1's insertion junction has. R that two of its walks cover, such as a second contig, and a junction between two walks, such as a record break whose bases the graph lacks, are gaps too                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Shared sites    | A site where every form but `other` is in one state, a call, tells none apart: it is out of the key and `formTag` leaves its marks out                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |

### What a form states

A form states, at every site that is not an edge, one state of its members: its
marks there, a no call, or a rarer change. `formContradictions` checks that
statement against each member's own departures. It recomputes the member's
changes itself, joining its departures, giving each a kind, and adding up plain
pieces at a summed site, and shares with `structuralForms` only `agrees`, the
tolerance, `changeKind`, the naming rule, and the constants. A test holds it to
that import.

- `call`: the form states a call only where the member's walks span the site's
  part inside the window, and `hidden` a no call only where they do not
- `rarer`: a rarer change only where the member makes one
- At a summed site, `net` and `size`: the member's net and largest change there,
  its plain pieces added up into one change where they reach the floor, agree
  with the state's
- At a wider site the member's changes touching it pair one-to-one with the
  state's marks, each pair covering half of each other and agreeing in net
  change and size: `count` where the member makes more or fewer changes there,
  `place` where no pairing by half-cover exists, `net` or `size` where none of
  those also agrees
- `unsited`: every change of the member lies in a site

`other` and edge sites state nothing, so nothing there can contradict. It checks
what a form's keys state: its key at each informative site and, at a shared one,
the state of its marks there. A form's `marks` and `net`, which `describeForm`
and `formTag` draw, are built from those states and not checked again.
structuralForms meets it by construction. A state folds, and a haplotype joins a
form, only through the same `agrees`, change by change under half-cover at a
wider site. A place is spans every two of which cover half of each other, and
two places of one size never share a key, so every carrier of a state pairs with
the carrier it shows. What agrees with nothing is `~` or `other`.

Members a result contradicts, summed over windows, edf4acd's nearest fold
against the fold that replaced it, as the checker before this one read them; not
rerun on the gate:

| Windows                                 | Floor, fold     | Nearest fold | These rules |
| --------------------------------------- | --------------- | ------------ | ----------- |
| 28 MC chr22 windows, format 6           | auto, folded    | 315          | 0           |
| 28 MC chr22 windows, format 6           | 50 bp, folded   | 1,849        | 0           |
| 28 MC chr22 windows, format 6           | 50 bp, unfolded | 663          | 0           |
| jb2hubs gate, 19 loci and 40 random, f6 | auto, folded    | 417          | 0           |
| jb2hubs gate, 19 loci and 40 random, f6 | 50 bp, folded   | 1,073        | 0           |
| jb2hubs gate, 19 loci and 40 random, f6 | 50 bp, unfolded | 206          | 0           |

That checker missed what it claimed to check. It compared marks that merely
touched rather than half-covered, so a member of one deletion's state passed in
the state of a deletion 5.5 kb over. It looked only around the form's marks, so
a member with the form's deletion and two more changes elsewhere in the site
passed. It imported the marks, the fold and the state lookup it was checking. It
also reported a rare shift whose insertion lay at a summed site as a net change
it was not. The pairing checker on the 28 windows, members contradicted, all
`place`:

| Floor, fold         | auto | 50 bp | 1 kb | 10 kb |
| ------------------- | ---- | ----- | ---- | ----- |
| ec8ecd9, folded     | 5    | 24    | 6    | 0     |
| ec8ecd9, unfolded   | 5    | 30    | 7    | 0     |
| These rules, either | 0    | 0     | 0    | 0     |

Two causes underlay them. Places of one size that start together but do not
half-cover each other shared a key, such as a 61 bp insertion at LCR22A's
18,651,717 and a 24 bp replacement there gaining 61 bp. A place's spans half
covered its first span but not each other, so the carrier a state shows could
lie outside another carrier's half-cover. The fold also matched places by key,
not by half-cover. Over the 28 windows at the auto floor, folded, 5 more members
are now `~`, 237 against 232. Unfolded, a state split this way adds 2 forms, 441
against 439.

The cost is members that no common form fits: at the auto floor 70 in `other`
and 371 in a form with a rarer change, of 25,887 placed over the gate's 59
windows; at 50 bp folded, 561 and 502. They claim no size there. The gate's
claims, 50 bp unfolded, scored by `scoreForms.py --complete-claims`: 6,981 of
7,003 at the loci and 13,013 of 13,013 in the random windows, against 6,931 of
6,952 and 12,946 of 12,946 with the nearest fold, which scored 51 and 67 rare
carriers apart and missed 4 and 18 of them. The 28 windows are the `W` lines of
gfa-to-tabix's `mc22.forms-compare-review.txt`.

The nearest fold this replaced folded any rare state to the common state nearest
in net, a rare form to the form nearest in net, and a rare no call to the
majority. That kept the rows few and named, and at LPA, AMY1 and MHC raised the
prototype's claims over a majority fold (LPA 381 to 395 of 406, MHC 367 to 407);
the part of that gain that is honest stays, since a state within tolerance still
folds, and the rest becomes `~` or `other`, which claims no size.

### Mark kinds

| Kind                                   | Meaning                                                                |
| -------------------------------------- | ---------------------------------------------------------------------- |
| `insertion`, `deletion`, `replacement` | `changeKind` on the bp of R left out and the bp of the haplotype's own |
| `inversion`                            | Net change under the floor, inverted bp at or over it                  |
| `duplication`                          | A step back along R; `repeatedBp` and `copied` give the copy           |

A deletion and an insertion elsewhere stay two changes, whatever their sizes and
however near: the file states a loss here and a gain there. A `shift` kind that
paired them by a 10% size match within twice the size had said the lost sequence
was the gained one, which no row states; it is gone.

Every mark states its two sides: `refLength`, the bp of R its departures leave
out, and `altLength`, the bp they run, with `delta` their difference. A
deletion's `altLength` is 0 and an insertion's `refLength` is 0. No ratio or
size makes a replacement a deletion or an insertion: GSTT1's single departure of
37,159 bp of R against 6 bp of its own is `repl 37.2kb->6bp`. A caller colours a
replacement as a loss or gain by the sign of `delta`.

Duplication is a kind of its own because core JBrowse already has the term
(`DUP:TANDEM` folds to `DUP`, and VcfFeature maps it to `tandem_duplication`),
so the plugin can colour them as a gain. AMY1 has 231 back-steps of 50 bp or
more over 160 haplotypes, 19-65 kb each (median 48.6 kb); LPA and MHC have none,
since their extra copies show as insertions. Whether a copy reads as a back-step
depends on how the graph collapsed the repeat, so the key treats a duplication
like an insertion of its size. Moving its REF onto the copied span instead would
join neighbouring sites: AMY1 goes from 18 to 23 forms (ARI 0.915).

## Labels

`markTag` writes a mark as a short token, `formTag` a form as its largest token
with ` +N` for the rest, and `describeForm` every token, largest first; the
README's table lists them. `sizeTag` writes 3 significant figures. A mark gives
one size where a side is 0 bp and both otherwise, R's first, so no label leaves
a sum to its reader: `repl 6.35kb->26.4kb` at C4, and at AMY1 `dup 94.1kb` for a
step back and `dup 113kb->207kb` for one joined to a 94 kb deletion and 19 kb
walked reversed. Where two forms' row labels would read the same, `formTag`
names, after the largest change, the changes that tell them apart, each the
largest of those the most of the alike forms lack.

Row labels at the `auto` floor, 462 haplotypes once CHM13 is left out:

| Window | Haplotypes | `formTag`                            | Note                                           |
| ------ | ---------- | ------------------------------------ | ---------------------------------------------- |
| C4     | 217        | `del 6.37kb`                         | A 6,367 bp deletion beside six SNPs            |
| C4     | 95         | `ins 32.7kb +1`                      | And a 32.7 kb deletion                         |
| C4     | 13         | `del 32.7kb, repl 6.32kb->32.7kb`    | Both had read `del 32.7kb +1`                  |
| C4     | 9          | `del 32.7kb, del 6.37kb`             |                                                |
| CYP2D6 | 38         | `ins 13.7kb`                         | 12,156 and 1,564 bp one base apart             |
| GSTT1  | 66         | `repl 37.2kb->6bp`                   | One departure                                  |
| MHC    | 26         | `del 86kb +2`                        | One 86,048 bp departure                        |
| MHC    | 46         | `repl 3.77kb->48.1kb, ins 1.2kb +11` | Three forms had read `repl 3.77kb->48.1kb +12` |
| AMY1   | 164        | `del 94.2kb, repl 13bp->94.1kb`      | H3r                                            |
| AMY1   | 29         | `dup 113kb->207kb`                   |                                                |

## Display floor

`floorBp: 'auto'` takes the window's width / 100, about 10 px, snapped to 50 bp
under 25 kb, 1 kb from 25 to 300 kb and 10 kb above. Forms that become identical
at a higher floor merge.

Forms, every one with 5 or more haplotypes once rare forms fold, measured under
the nearest fold and not remeasured:

| Window      | 50 bp | 500 bp | 1 kb | 5 kb |
| ----------- | ----- | ------ | ---- | ---- |
| GSTM1 35 kb | 2     | 2      | 2    | 2    |
| C4 70 kb    | 11    | 10     | 10   | 10   |
| LPA 130 kb  | 26    | 30     | 30   | 30   |
| AMY1 150 kb | 13    | 14     | 14   | 14   |
| MHC 90 kb   | 29    | 21     | 14   | 9    |

At gene scale a fixed 1 kb floor leaves the track blank: in HP
72,040,000-72,045,000, a 300 bp deletion in 5 haplotypes falls to one "no
change" row, and a rand1 tile with a 51-55 bp insertion does the same. The auto
rule gives that 5 kb window 50 bp and the five windows above 1 kb.

## Abridged walks

`formsCut` expands no abridged walk. A piece the cut does not load comes from
its tile's `blocks` unit with its steps on R and the bp of its own between them,
and its departures are those its detail record gives: forms key on a departure's
place and size, never on its nodes.

- `scripts/blocks-equality.ts` holds the two equal on a walk file: over random
  windows on its abridged tiles, each haplotype's cover and departures and the
  structural forms from the cut that expands nothing against the cut that
  expands every walk
- Until graph-walks 7 the forms cut read abridged walks in full from their
  detail records, fewest steps first under `FORMS_STEP_BUDGET`, and left any
  past it as a `net` mark over its whole span. That cost a second cut over the
  same tiles, and a haplotype whose piece the first cut dropped was never
  expanded, so it came back with nothing covered

`structuralForms` takes under 0.1 s in every window, so the cut dominates.
`formsCut` spends about 1 s over a single expanded cut on its first, plain cut.

### What a `dp` tag must list

structuralForms reads more than the 50 bp departures, so a `dp` tag of those
alone would change the forms:

- Departures already joined at `joinBp` (50 bp), as `join` joins them, so each
  entry is one change with its exact net
- The nonzero pieces under 50 bp inside each change of 50 bp or more. At a site
  up to `SUMMED_SITE_BP` (1 kb) wide, structuralForms sums a haplotype's plain
  pieces, which lifts rand29's VNTR from 235 haplotypes matching the assembly to
  408 of 408
- Gate it on forms from `dp` equalling the expanded forms at MHC, KIR, LPA, AMY1
  and the jb2hubs gate windows

## Validation

Through `formsCut` and `structuralForms` against the published file, expanded,
scored against the assemblies (`truth.tsv`) and the sidecar (`forms-v3.jsonl`):

| Window | Walk net = assembly | Forms 50 bp / 1 kb | ARI vs sidecar, 50 bp | Sidecar forms |
| ------ | ------------------- | ------------------ | --------------------- | ------------- |
| GSTM1  | 405/405             | 2 / 2              | 0.86                  | 6             |
| C4     | 407/407             | 11 / 10            | 0.95                  | 15            |
| LPA    | 406/406             | 26 / 30            | 0.37                  | 179           |
| MHC    | 407/407             | 29 / 14            | 0.67                  | 120           |

MHC's net is scored on a cut widened to the truth table's span,
chr6:32,476,308-32,610,000. The ARI is low where the sidecar splits LPA into 179
forms and MHC into 120, most of them single haplotypes.

The jb2hubs gate (`printWalkForms.ts`, `scoreForms.py --complete-claims`) over
19 loci and 40 random 50 kb windows, at 50 bp with `foldRareForms: false`, on
the published format 5 file before the nearest fold was replaced (format 6
figures under [What a form states](#what-a-form-states)):

- Claims within 10% or 500 bp of the assembly: 6,910 of 6,931 at the loci and
  12,946 of 12,946 in the random windows. The sidecar sizes 5,528 and 12,860
- Per site, 50 bp to 1 kb, members both score: walks 78.8% against the sidecar's
  69.9% at the loci, 96.1% against 92.7% in the random windows
- Of the 21 misses, 20 are at KIR, where a walk and the assembly's flanks put a
  16 kb or 31 kb stretch on opposite sides of a span's edge among the KIR
  paralogs. The other, NA18508#1 at FCGR, walks one contig 82 kb longer than the
  flanks measure

## Consumers

- The plugin's `WalkFormsAdapter` feeds core `LinearMultiSampleVariantDisplay`,
  one row per form, pinned to a window, with the floor in its menu. It caches
  `HaplotypeDepartures` by walk file and pinned window, since a new floor
  regroups them without a new cut
- jb2hubs moves the HPRC and cattle pages to walk forms; cattle needs a walk
  file built from its P lines. Arabidopsis keeps the sidecar until its graph has
  paths
