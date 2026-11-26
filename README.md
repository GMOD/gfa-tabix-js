# @gmod/gfa-tabix

The reader for [gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix)'s files, the
GFA tabix file and the walk file, as [@gmod/bam](https://github.com/GMOD/bam-js)
is for BAM.

![A stretch of chromosome 6 zooms out into each genome's route through it, as colored lines through shared pieces of DNA, one detouring through extra DNA and one skipping a part](docs/hero.svg)

## Install

```bash
npm install @gmod/gfa-tabix
```

## Quick start

`WalkFile` reads the walk file, which holds each haplotype's route:

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

## The GFA tabix file

The GFA tabix file holds the graph alone, segments and links with no walks.
`GfaTabixFile` reads it as GFA text, for a viewer that draws the graph itself,
and needs `@gmod/tabix`. Its use is in the
[API reference](docs/api.md#the-gfa-tabix-file).

## Docs

- [API reference](docs/api.md)

## See also

- [gfa-to-tabix](https://github.com/GMOD/gfa-to-tabix), which writes the files
  this package reads

## License

MIT
