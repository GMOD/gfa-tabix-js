import { LocalFile } from 'generic-filehandle2'

import { WalkFile } from './walkFile.ts'

import type { WalkDeparture } from './walkDeparture.ts'

const DATA = new URL('../test_data/', import.meta.url).pathname

// Each call a cut makes to unattachedRecords, with what walkCut, as the cut
// ran it before the attach fetch until unattachedRecords replaced that cut,
// lists as unattached on the same runs, nodes, records and options
const calls = vi.hoisted(() => ({
  refName: '',
  seen: [] as { fresh: WalkDeparture[]; old: WalkDeparture[] }[],
}))

vi.mock(import('./walkRows.ts'), async importOriginal => {
  const rows = await importOriginal()
  return {
    ...rows,
    unattachedRecords: (
      ...[readRuns, nodes, records, options]: Parameters<
        typeof rows.unattachedRecords
      >
    ) => {
      const old = rows
        .walkCut(readRuns, nodes, calls.refName, records, options)
        .unattached.map(({ departure }) => departure)
      const fresh = rows.unattachedRecords(readRuns, nodes, records, options)
      calls.seen.push({ fresh, old })
      return fresh
    },
  }
})

function reader(name: string) {
  const path = `${DATA}${name}.walks.gz`
  return new WalkFile({
    walksFilehandle: new LocalFile(path),
    idxFilehandle: new LocalFile(`${path}.idx`),
  })
}

test.each(['attachmany.R', 'attach.GRCh38', 'ownends.R', 'farrejoin.R'])(
  'unattachedRecords finds, record for record, the ones the cut before the fetch found, over every window of %s',
  async name => {
    const walks = reader(name)
    const [refName = ''] = await walks.getReferenceSequences()
    calls.refName = refName
    calls.seen.length = 0
    const length = await walks.getSequenceLength(refName)
    for (let start = 0; start < length; start += 5) {
      for (const width of [1, 10, 50]) {
        for (const options of [
          { context: 0 },
          { context: 5 },
          { context: 5, contained: true },
          { context: 0, attachBp: 5 },
          { context: 20, contained: true, expand: () => true },
        ]) {
          await walks.getCutForRange(refName, start, start + width, options)
        }
      }
    }
    const differ = calls.seen.filter(
      ({ fresh, old }) =>
        fresh.length !== old.length || fresh.some((d, k) => d !== old[k]),
    )
    expect(differ).toEqual([])
    expect(
      calls.seen.filter(({ old }) => old.length > 0).length,
    ).toBeGreaterThan(20)
  },
)
