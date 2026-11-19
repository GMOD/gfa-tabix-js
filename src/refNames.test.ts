import {
  AmbiguousSequenceError,
  refNameLookup,
  resolveRefName,
} from './refNames.ts'

test('a PanSN name resolves under its sample and under its haplotype', () => {
  const lookup = refNameLookup([
    'GRCh38#0#chr6',
    'NA20809#2#CM094351.1',
    'chrM',
  ])
  expect(resolveRefName(lookup, 'GRCh38', 'chr6')).toBe('GRCh38#0#chr6')
  expect(resolveRefName(lookup, 'GRCh38#0', 'chr6')).toBe('GRCh38#0#chr6')
  expect(resolveRefName(lookup, 'NA20809', 'CM094351.1')).toBe(
    'NA20809#2#CM094351.1',
  )
  expect(resolveRefName(lookup, 'hg38', 'chrM')).toBe('chrM')
  expect(resolveRefName(lookup, 'hg38', 'GRCh38#0#chr6')).toBe('GRCh38#0#chr6')
  expect(resolveRefName(lookup, 'NA20809#1', 'CM094351.1')).toBeUndefined()
})

test('a pggb name with no haplotype field resolves under its sample', () => {
  const lookup = refNameLookup(['chm13#chr1'])
  expect(resolveRefName(lookup, 'chm13', 'chr1')).toBe('chm13#chr1')
  expect(resolveRefName(lookup, 'chm13#0', 'chr1')).toBeUndefined()
})

test('a contig on two haplotypes of a sample resolves on each haplotype and throws on the sample', () => {
  const lookup = refNameLookup(['HG002#1#chr1', 'HG002#2#chr1', 'HG002#1#chr2'])
  const onSample = () => resolveRefName(lookup, 'HG002', 'chr1')
  expect(onSample).toThrow(AmbiguousSequenceError)
  expect(onSample).toThrow(
    expect.objectContaining({
      message:
        'chr1 on HG002 is 2 sequences of the file, HG002#1#chr1, HG002#2#chr1; name one of them, or its haplotype: HG002#1, HG002#2',
      sequences: ['HG002#1#chr1', 'HG002#2#chr1'],
    }),
  )
  expect(resolveRefName(lookup, 'HG002#1', 'chr1')).toBe('HG002#1#chr1')
  expect(resolveRefName(lookup, 'HG002#2', 'chr1')).toBe('HG002#2#chr1')
  expect(resolveRefName(lookup, 'HG002', 'HG002#2#chr1')).toBe('HG002#2#chr1')
  expect(resolveRefName(lookup, 'hg38', 'HG002#1#chr1')).toBe('HG002#1#chr1')
  expect(resolveRefName(lookup, 'HG002', 'chr2')).toBe('HG002#1#chr2')
  expect(resolveRefName(lookup, 'HG002#2', 'chr2')).toBeUndefined()
})
