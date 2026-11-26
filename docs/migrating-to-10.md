# Migrating to @gmod/gfa-tabix 10

`@jbrowse/graph-walks` 9 is `@gmod/gfa-tabix` 10, the reader for gfa-to-tabix's
GFA tabix file and walk file, shaped like @gmod/bam: a class per file, opened
from one options object, with `get`-prefixed methods that take
`refName, start, end`.

```bash
npm uninstall @jbrowse/graph-walks
npm install @gmod/gfa-tabix
```

Coordinates stay 0-based half-open.

## Opening a walk file

| 9                                                    | 10                                                    |
| ---------------------------------------------------- | ----------------------------------------------------- |
| `WalkReader.open(url)`                               | `new WalkFile({ walksUrl: url })`                     |
| `WalkReader.open(path, p => new LocalFile(p))`       | `new WalkFile({ walksPath: path })`                   |
| `WalkReader.fromFilehandles({ data, directory })`    | `new WalkFile({ walksFilehandle, idxFilehandle })`    |
| `new WalkReader(new WalkFile(data, directory))`      | `new WalkFile({ walksFilehandle, idxFilehandle })`    |
| the directory at `url` plus `.idx`, before any query | the same by default; `idxUrl`, `idxPath` name another |
| `WalkFile`, the container class                      | internal, not exported                                |
| `WalkReaderHeader`                                   | `WalkFileHeader`                                      |

## Methods

| 9                                                       | 10                                                                        |
| ------------------------------------------------------- | ------------------------------------------------------------------------- |
| `reader.cutRegion({ refName, start, end }, opts)`       | `walks.getGraphForRange(refName, start, end, opts)`                       |
| `reader.cut(refName, { start, end }, opts)`             | `walks.getCutForRange(refName, start, end, opts)`                         |
| `graph.tables(kept, fragments, departures)` after `cut` | `cut.tables()`, which holds the gaps too                                  |
| `reader.counts(refName, { start, end }, opts)`          | `walks.getCountsForRange(refName, start, end, opts)`                      |
| `reader.bytesFor(refName, { start, end }, opts)`        | `walks.getBytesForRange(refName, start, end, opts)`                       |
| `reader.referenceNodes(refName, { start, end })`        | `walks.getReferenceNodes(refName, start, end)`                            |
| `formsCut(reader, { refName, start, end }, opts)`       | `walks.getFormsForRange(refName, start, end, opts)`                       |
| `reader.header()`                                       | `walks.getHeader()`                                                       |
| `reader.referenceSequences()`                           | `walks.getReferenceSequences()`                                           |
| `reader.sequenceFor(sample, contig)`                    | `walks.getSequenceName(sample, contig)`                                   |
| `reader.units(refName)`                                 | `walks.getUnits(refName)`                                                 |
| `reader.sequenceLength(refName)`                        | `walks.getSequenceLength(refName)`                                        |
| `rgfaCut(files, { refName, start, end }, { hops })`     | `new GfaTabixFile(files).getGraphForRange(refName, start, end, { hops })` |
| `RgfaFiles`                                             | `GfaTabixFileOptions`                                                     |

## Options and results

| 9                                                                    | 10                                                            |
| -------------------------------------------------------------------- | ------------------------------------------------------------- |
| `RegionCutOptions`                                                   | `GraphForRangeOptions`                                        |
| `WalkCutOptions` for `cut`                                           | `CutForRangeOptions`, `WalkCutOptions` with `referenceSample` |
| `WalkCountOptions` for `counts`                                      | `CountsForRangeOptions`                                       |
| `FormsCutOptions`                                                    | `FormsForRangeOptions`                                        |
| `expand: name => …`                                                  | the same, or a list of PanSN prefixes: `expand: ['HG00097']`  |
| `attachFetches`                                                      | gone; it was ignored                                          |
| `cutRegion`'s `{ tables, refName, sample }`                          | `getGraphForRange`'s, plus `sequences` with `sequences: true` |
| `cut`'s `{ graph, kept, fragments, departures, gaps, placeholders }` | the same, plus `refName` and `tables()`                       |

## Sequence names

- Every range method, `getCutForRange`, `getCountsForRange`, `getBytesForRange`,
  `getReferenceNodes` and `getFormsForRange` among them, takes a contig (`chr6`)
  as `getGraphForRange` does, resolved against `referenceSample`, as well as the
  file's own name
- A sequence the file lacks throws `MissingSequenceError`, with `sequences`
  listing those it holds, from every range method. `cutRegion` threw a plain
  `Error` with the same message, `cut` returned an empty cut, and `formsCut` cut
  the name as given

## Exports

The package exports its public API alone. Gone from the exports: `WalkGraph`,
`walkCut`, `walkHeader`, `joinPieces`, `squeezer`, `byteBudgetError`,
`stepBudgetError`, `keptNodesError`, `walkNameFilter`, `parseDepartureRow`,
`adjacentSteps`, `departureStretch`, `leavesOrRejoinsOn`,
`panSNHaplotypeOrSample`, `formsCut`, `rgfaCut`, `formContradictions`,
`carrierNames`, `MOVED_TAG`, `autoFloor`, `changeName`, `haplotypeMarks`,
`minCarriersFor`, and the types `WalkHeader`, `WalkRow`, `WalkRun`,
`WalkWindow`, `DepartureRecord`, `Hole`, `GraphHeader`, `RgfaLink` and
`StableSpan`.
