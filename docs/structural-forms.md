# Structural forms

A form is a group of haplotypes that make the same structural changes over a
window. `getFormsForRange` reads each haplotype's departures from the reference,
and `structuralForms` groups them.

```ts
import { describeForm, formTag, structuralForms } from '@gmod/gfa-tabix'

const region = { refName: 'chr1', start: 109_680_000, end: 109_715_000 }
const haplotypes = await walks.getFormsForRange(
  region.refName,
  region.start,
  region.end,
)
const result = structuralForms(haplotypes, { window: region })
result.forms.map(f => `${f.members.length} ${describeForm(f)}`)
// ['306 del 18.4kb', '156 ref']
result.forms.map(f => formTag(f, result))
// ['del 18.4kb', 'ref']
```

- `getFormsForRange` returns, for each haplotype, its departures with exact REF
  span and ALT length, the bp it walks inverted, and the bp it walks again
- `getFormsForRange` keeps to `FORMS_STEP_BUDGET` steps, or a larger
  `stepBudget`, and takes `byteBudget`. It throws a `RegionTooLargeError` past
  either
- `getFormsForRange` expands no abridged walk, and returns the departures an
  expanded cut gives
- `walkDepartures` reads the same departures from a cut you already have
- A result is plain data, so a worker can post it. Grouping again at another
  floor needs no new cut, and `structuralForms` took under 0.1 s in every window
  we measured

## Options

| Option          | Effect                                                                             | Default          |
| --------------- | ---------------------------------------------------------------------------------- | ---------------- |
| `window`        | The region. Marks must touch it, and `auto` reads its width                        |                  |
| `floorBp`       | The smallest change a form tells apart                                             | `'auto'`         |
| `minCarriers`   | The fewest haplotypes a state or form needs                                        | 1%, at least one |
| `joinBp`        | Departures of one haplotype nearer than this on the reference join into one change | 50               |
| `exclude`       | Samples to leave out, such as a second reference stored as a haplotype             | none             |
| `foldRareForms` | Fold states and forms under `minCarriers` into a common one they agree with        | true             |

`floorBp: 'auto'` is 50 bp for a window under 25 kb, 1 kb from 25 to 300 kb and
10 kb above. Forms that become identical at a higher floor merge.

## What a result holds

| Field          | Contents                                                                                   |
| -------------- | ------------------------------------------------------------------------------------------ |
| `forms`        | Each form's `members`, `marks`, `net` change, and `states`, one key per informative site   |
| `sites`        | Each site's span, its `states` with their `carriers`, and whether it is `edge` or `shared` |
| `informative`  | The indices of the sites that tell forms apart                                             |
| `placed`       | How many haplotypes are in a form, so a form's share is its members over `placed`          |
| `unplaced`     | The haplotypes with no call in the window and no form                                      |
| `rareCarriers` | The haplotypes with a state or form under `minCarriers`                                    |

- A mark has `kind`, `refStart`, `refEnd`, `refLength`, `altLength` and `delta`,
  and in a result its `site` and its `state` key there
- `formStateAt(result, form, site)` is a form's state key at any site.
  `REFERENCE_STATE`, `NO_CALL_STATE` and `RARER_STATE` are the keys for no
  change, no call and a rare change that agrees with no common state
- The last form can be `other`, with key `OTHER_STATE`. It holds the haplotypes
  no common form fits, and has no marks

## Mark kinds

| Kind          | Meaning                                                                      | Label                 |
| ------------- | ---------------------------------------------------------------------------- | --------------------- |
| `deletion`    | Leaves out reference and has 0 bp of the haplotype's sequence in its place   | `del 18.4kb`          |
| `insertion`   | Has sequence of the haplotype and leaves out 0 bp of reference               | `ins 22kb`            |
| `replacement` | Has bp on both sides, at any ratio                                           | `repl 32.7kb->58.7kb` |
| `inversion`   | Net change under the floor, with inverted bp at or over it                   | `inv 12kb`            |
| `duplication` | A step back along the reference; `repeatedBp` and `copied` describe the copy | `dup 94.1kb`          |

- `refLength` is the bp of reference the mark leaves out, `altLength` the bp of
  the haplotype in their place, and `delta` their difference. A caller colours a
  replacement as a loss or a gain by the sign of `delta`
- No size threshold turns a replacement into a deletion or an insertion: GSTT1's
  departure of 37,159 bp of reference against 6 bp is `repl 37.2kb->6bp`
- A deletion and an insertion elsewhere are two changes, however near they are
- A mark is made of the departures that change the reference's length, walk it
  backwards or walk it again, so a SNP beside a deletion neither widens the
  deletion nor changes its label

## Labels

- `markTag(mark)` is the label of one mark, as in the table above
- `formTag(form, result)` names a form in about 20 characters for a row label:
  its largest change, then ` +N` for the rest, as in `del 32.7kb +1`. It leaves
  out the changes every form shares, and where two forms would get the same
  label it names the changes that tell them apart, as in
  `del 32.7kb, del 6.37kb`
- `describeForm(form)` lists every change, largest first
- `sizeTag` writes sizes to 3 significant figures: `94bp`, `6.46kb`, `1.2Mb`
- A form with no change is `ref`, one with only shared changes is `shared`, and
  the rare and no-call states are `rarer` and `no call`

## Grouping rules

| Rule         | Decision                                                                                                                          |
| ------------ | --------------------------------------------------------------------------------------------------------------------------------- |
| Floor        | Changes under 50 bp are left out. At a site up to 1 kb wide the floor applies to a haplotype's net change over the site           |
| Join         | Departures of one haplotype less than `joinBp` apart join. One under 50 bp on both sides joins only another such                  |
| Sites        | A mark that `minCarriers` haplotypes share is common, and overlapping common marks form one site                                  |
| States       | A state is one kind of change at one size, to 2 significant figures. At a site wider than 1 kb it also has a place                |
| Rare states  | A state under `minCarriers` folds into the nearest common state it agrees with, and is `rarer` where it agrees with none          |
| Rare forms   | A haplotype in a form under `minCarriers` joins the largest common form it agrees with at every site, and `other` where none fits |
| No calls     | A haplotype has no call at a site its walks do not span inside the window, or where two of its walks cover the same reference     |
| Edge sites   | A site where most haplotypes have no call is `edge` and left out of the form key                                                  |
| Shared sites | A site where every form has the same state is `shared`, left out of the key and out of `formTag`                                  |

Two marks are shared where each covers at least half of the other on the
reference. Two changes agree where both are under the floor, or they differ by
less than the floor or by less than a tenth of the larger.

`foldRareForms: false` keeps every haplotype's state at every site, which a
scorer comparing sizes wants.

## Limits

- A walk lying wholly inside an inverted block looks the same as a contig the
  graph stores back to front, so a window inside an inversion reports the
  inversion's haplotypes as reference. A window that reaches either flank by one
  step tells the two apart
- Whether an extra copy of a repeat is a `duplication` or an `insertion` depends
  on how the graph collapsed the repeat. AMY1's copies are duplications, and
  LPA's and MHC's are insertions

## Validation

Scored on the published HPRC v2.1 GRCh38 walk file against the assemblies, with
462 haplotypes once GRCh38 and CHM13 are left out:

| Window | Walk net = assembly | Forms at 50 bp / 1 kb |
| ------ | ------------------- | --------------------- |
| GSTM1  | 405/405             | 2 / 2                 |
| C4     | 407/407             | 11 / 10               |
| LPA    | 406/406             | 26 / 30               |
| MHC    | 407/407             | 29 / 14               |

[The design notes](../agent-docs/structural-forms-design.md) have the
measurement behind each grouping rule and the full validation.
