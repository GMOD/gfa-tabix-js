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
const { tables, refName, sample } = await walks.getGraphForRange(
  'chr6',
  160_614_797,
  160_647_758,
  { haplotypes: ['HG00097', 'HG00133'] },
)
```

URLs read by range requests in node and the browser, `walksPath` reads a local
file, and the directory is the walk file plus `.idx` unless `idxUrl` or
`idxPath` says otherwise.

## Coordinates

Coordinates are 0-based half-open, as in BED, where a samtools-style region is
1-based closed: `chr6:160,614,798-160,647,758` is
`getGraphForRange('chr6', 160_614_797, 160_647_758)`.

## The walk file

`gfa-to-tabix --walks` writes `<prefix>.<sample>.walks.gz` and its directory
`<prefix>.<sample>.walks.gz.idx`. `WalkFile` reads the stretch of it under a
range: the reference, each haplotype's walk through it, and the nodes and links
they visit, as typed arrays. `refName` is a contig of the reference sample
(`chr6`) or the file's full name (`GRCh38#0#chr6`). `getCutForRange`,
`getCountsForRange` and `getFormsForRange` give the cut itself, per-node
haplotype counts and structural forms. `getRunsForRange` gives a range's walks
as runs of reference stretches, without decoding a step.

## The GFA tabix file

`GfaTabixFile` reads gfa-to-tabix's segments and links, Tabix-indexed BED, as
GFA text:

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
- **Window:** a stretch of the reference, also called a region
- **Cut:** the nodes, links and walks under a window

## Docs

Reference:

- [Cutting a window](docs/cutting.md): every method and option, and the GFA
  tabix file
- [Large genomic regions and abridged walks](docs/reads.md):
  `RegionTooLargeError`

Views of the haplotypes under a window:

- [Counts](docs/counts.md): haplotypes through each node and link
- [Runs](docs/runs.md): each walk as reference stretches and the DNA between
  them
- [Structural forms](docs/structural-forms.md): haplotypes that make the same
  change

## License

MIT
