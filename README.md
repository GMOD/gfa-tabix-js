# @jbrowse/graph-walks

`@jbrowse/graph-walks` is the reader for
[gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix)'s pangenome files, as
[@gmod/bam](https://github.com/GMOD/bam-js) is for BAM.

It fetches only the part of the graph under a stretch of the reference, such as
one gene, by range requests from a web server, so a genome browser can draw it.

![A stretch of chromosome 6, and an arrow to each genome's route through it, as colored lines through shared pieces of DNA, one detouring through extra DNA and one skipping a part](docs/hero.svg)

A pangenome graph stores many people's genomes together: one reference genome,
with a branch wherever someone's DNA differs from it.

- **In:** a `gfa-to-tabix --walks` file, `<prefix>.<sample>.walks.gz` with its
  directory `<prefix>.<sample>.walks.gz.idx`, or a segment index,
  `<prefix>.segs.bed.gz` and `<prefix>.links.bed.gz` with theirs or format 2's
  `<prefix>.graph.bed.gz`, plus a reference region such as
  `chr6:160,614,798-160,647,758`
- **Out:** from a walk file, the reference over the region, every haplotype's
  walk through it, and the nodes and links they visit, as typed arrays; from a
  segment index, the segments and links under the region as GFA text
- The reader reads walk format 8, which gfa-to-tabix 5.0 and later write, and
  refuses an older format with a `WalkFormatError` asking for a rebuild, and a
  newer one asking for a newer reader
- It reads a walk file through
  [generic-filehandle2](https://github.com/GMOD/generic-filehandle2)'s handles
  and gunzips with
  [@gmod/bgzf-filehandle](https://github.com/GMOD/bgzf-filehandle); a segment
  index through any object with
  [@gmod/tabix](https://github.com/GMOD/tabix-js)'s `TabixIndexedFile` methods

The JBrowse graph plugin, BandageJS, jb2hubs and the sequenceTubeMap fork read
through it for a browser track, and ggbandage and ggtubemap for a static figure
([Cut presets](docs/cutting.md#cut-presets)).

## Install

```bash
npm install @jbrowse/graph-walks
```

A segment index also needs `@gmod/tabix`.

## Usage

```ts
import { WalkReader, walksKept } from '@jbrowse/graph-walks'

const reader = WalkReader.open(
  'https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38-f8.GRCh38.walks.gz',
)

const refName = 'GRCh38#0#chr6'
const { graph, kept, fragments, departures } = await reader.cut(
  refName,
  { start: 160_614_798, end: 160_647_758 },
  { keep: walksKept(['HG00097', 'HG00133'], refName) },
)
const tables = graph.tables(kept, fragments, departures)
```

- `WalkReader.open(url)` reads the walk file at `url` and its directory at `url`
  plus `.idx`, put before any query string or fragment, by HTTP range requests
  through generic-filehandle2's `RemoteFile`. Its second argument opens each
  location another way: `url => new LocalFile(url)` for a local file
- `WalkReader.fromFilehandles({ data, directory })` takes the two handles
  already open: `data` the `.walks.gz`, `directory` its `.idx`. Any object with
  a generic-filehandle2 handle's `read`, `readFile`, `stat` and `source` will do
- A file gfa-to-tabix 5.0 or later wrote, in walk format 8; the reader refuses
  any other format and says to rebuild the file. jbrowse.org hosts the HPRC
  release 2.1 walk files in format 8, each with its `.idx`:
  `https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38-f8.GRCh38.walks.gz` and
  `https://jbrowse.org/demos/hprc/hprc-v2.1-mc-chm13-f8.CHM13.walks.gz`
- `tables` lists the nodes with their lengths and reference positions, the links
  the cut's steps walk, each walk's steps, and every departure from the
  reference over the window ([Departures](docs/cutting.md#departures)).
  `graphTablesGFA` in
  [@jbrowse/bandage-core](https://github.com/GMOD/bandage-core) writes the
  nodes, links and walks as GFA
- `reader.header()` lists the reference sample and every haplotype with rows,
  and the tile size, `tile`, read from the walk file's header unit
  - `graph` is the hash of the GFA's segments, as its node and link index states
    it, and `leftout` how many of its lines have no rows in the file; each is
    absent where the file doesn't state it
- `reader.referenceSequences()` lists the sequences a region can be on, in full
  PanSN: `GRCh38#0#chr6`
- `reader.units(refName)` is a sequence's unit table,
  `reader.sequenceLength(refName)` where its reference ends, and
  `reader.referenceNodes(refName, region)` every pass of the reference over a
  region
- `reader.bytesFor(refName, region, { context, attachBp })` is the compressed
  bytes a cut of the region reads when it keeps every walk and expands none
  ([What a cut reads](docs/reads.md#what-a-cut-reads))
- `reader.counts(refName, region, opts)` is how many of the file's haplotypes
  step each node and link over a region, and how often, with each count exact or
  stated as a lower bound and why ([Counts](docs/counts.md))
- `reader.cutRegion` takes a contig such as `chr6` and haplotypes by name
  ([A region by contig](docs/cutting.md#a-region-by-contig)), `rgfaCut` cuts a
  segment index ([Segment indexes](docs/cutting.md#segment-indexes)), and
  `formsCut` with `structuralForms` groups the haplotypes by their structural
  changes ([Structural forms](docs/structural-forms.md#usage))

## Cut options

| Option                | Effect                                                                                                         |
| --------------------- | -------------------------------------------------------------------------------------------------------------- |
| `keep`                | Which walks to decode, the reference's among them, by `sample#haplotype#contig`; every one by default          |
| `context`             | bp of reference kept past each end of the region; default 1000, `WALK_CONTEXT`                                 |
| `stepBudget`          | The most steps a cut keeps, departures' own steps included; 4,000,000 by default                               |
| `byteBudget`          | The most compressed bytes the walk and attach reads fetch; 8 MB by default                                     |
| `departureByteBudget` | The most compressed bytes the departure read fetches; 2 MB by default                                          |
| `attachBp`            | The most bp of own sequence a departure gets its steps attached for; default 65,536                            |
| `attachFetches`       | Deprecated and ignored: the tiles of every record that gets its steps are fetched                              |
| `nodeLimit`           | The most nodes a cut keeps, placeholders not among them; no limit by default                                   |
| `placeholderLimit`    | The most placeholders a cut keeps; no limit by default                                                         |
| `contained`           | `true` keeps each walk within the window; see [Contained cuts](docs/cutting.md#contained-cuts)                 |
| `links`               | `'walked'` (default) the links the cut's walks take; `'between'` also other walks' between its nodes; `'none'` |
| `sequences`           | Keep each node's bases in `graph.sequences`, from a `--sequences` file                                         |
| `signal`              | An `AbortSignal` that stops the reads                                                                          |
| `expand`              | Which abridged walks to read in full, by `sample#haplotype#contig`; none by default                            |
| `budgetNames`         | The words a `RegionTooLargeError` message uses for the step and byte budgets                                   |

- A cut over a budget or its node limit throws a `RegionTooLargeError`. The
  error's `fitsBp` holds about how wide a window would fit, and JBrowse shows it
  as a zoom-in notice
- `WALK_BUDGETS` holds the default budgets and `attachBp`, what a browser window
  holds; a caller that wants no limit passes `Infinity`

## Docs

- [Cutting a window](docs/cutting.md): fragments, gaps, departures, links
- [Reads, budgets and abridged walks](docs/reads.md): what a cut fetches
- [Counts](docs/counts.md): haplotypes per node and link
- [Structural forms](docs/structural-forms.md): grouping haplotypes by change
- [Migrating to graph-walks 8](docs/migrating-to-8.md): changes in 7 and 8
- [ADR 0001](docs/adr/0001-a-cut-states-what-lies-between-fragments.md): gaps
  between fragments

## License

MIT
