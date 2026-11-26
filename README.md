# @gmod/gfa-tabix

The reader for [gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix)'s files, the
GFA tabix file and the walk file, as [@gmod/bam](https://github.com/GMOD/bam-js)
is for BAM.

![A stretch of chromosome 6 zooms out into each genome's route through it, as colored lines through shared pieces of DNA, one detouring through extra DNA and one skipping a part](docs/hero.svg)

## Install

```bash
npm install @gmod/gfa-tabix
```

The GFA tabix file also needs `@gmod/tabix`.

## Quick start

```ts
import { WalkFile } from '@gmod/gfa-tabix'

const walks = new WalkFile({
  walksUrl:
    'https://jbrowse.org/demos/hprc/hprc-v2.1-mc-grch38-f8.GRCh38.walks.gz',
})
const { tables } = await walks.getGraphForRange(
  'chr6',
  160_614_797,
  160_647_758,
  { haplotypes: ['HG00097', 'HG00133'] },
)
```

- `tables` holds the reference over the range, each haplotype's walk through it,
  and the nodes and links they visit, as typed arrays
- A URL is read by range requests in node and the browser, and `walksPath` reads
  a local file. The directory is the walk file plus `.idx`
- Coordinates are 0-based half-open, as in BED, so the samtools-style region
  `chr6:160,614,798-160,647,758` is `160_614_797, 160_647_758`
- `gfa-to-tabix --walks` writes the walk file, `<prefix>.<sample>.walks.gz`, and
  its directory

## Methods

Each `WalkFile` method below takes a contig and a range.

| Method              | Returns                                                      | Docs                                         |
| ------------------- | ------------------------------------------------------------ | -------------------------------------------- |
| `getGraphForRange`  | The graph and walks over the range, as tables                | [API reference](docs/api.md#tables)          |
| `getRunsForRange`   | Each walk as stretches of reference and the DNA between      | [Runs](docs/runs.md)                         |
| `getCountsForRange` | How many haplotypes pass through each node and link          | [Counts](docs/counts.md)                     |
| `getFormsForRange`  | The input for grouping haplotypes that make the same change  | [Structural forms](docs/structural-forms.md) |
| `getCutForRange`    | The cut as objects, for a caller that needs more than tables | [API reference](docs/api.md#fragments)       |

The [API reference](docs/api.md) lists every method, option and error.
[Large genomic regions and abridged walks](docs/reads.md) covers
`RegionTooLargeError` and what to do about a region too large to load.

## The GFA tabix file

`GfaTabixFile` reads gfa-to-tabix's segments and links, Tabix-indexed BED, as
GFA text ([API reference](docs/api.md#the-gfa-tabix-file)):

```ts
import { GfaTabixFile } from '@gmod/gfa-tabix'
import { TabixIndexedFile } from '@gmod/tabix'

const gfaTabix = new GfaTabixFile({
  graph: new TabixIndexedFile({ url: graphUrl, tbiUrl: `${graphUrl}.tbi` }),
})
const gfa = await gfaTabix.getGraphForRange('GRCh38#0#chr6', 0, 100_000)
```

## Words used here

- **Reference:** the one linear genome, such as GRCh38
- **Haplotype:** one copy of a chromosome from one person
- **Node and link:** a piece of DNA many haplotypes share, and the join to the
  next
- **Walk:** the route one haplotype takes through the nodes and links
- **Window:** a stretch of the reference, also called a region or range
- **Cut:** the nodes, links and walks under a window

## License

MIT
