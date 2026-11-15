"""What a genotype says of the allele of its record, on an invented record of one alternate allele."""

from __future__ import annotations

import pathlib
from collections.abc import Sequence

import pysam
import pytest

from case_annotation import vcf_values

# Each genotype as a VCF writes it: its alleles as read, whether it carries the record's allele, and its zygosity.
GENOTYPES: dict[str, tuple[tuple[int | None, ...], bool | None, vcf_values.Zygosity | None]] = {
    '0/0': ((0, 0), False, 'hom_ref'),
    '0/1': ((0, 1), True, 'het'),
    '1/1': ((1, 1), True, 'hom_alt'),
    '0|1': ((0, 1), True, 'het'),
    # A haploid call, as on a male's X or Y.
    '1': ((1,), True, 'hom_alt'),
    '0': ((0,), False, 'hom_ref'),
    # No call.
    './.': ((None, None), None, None),
    '.': ((None,), None, None),
    # One allele missing: a copy of the allele is carried whatever the missing one is, and without a copy the missing
    # one may be the allele. Neither genotype says on how many copies the allele lies.
    './1': ((None, 1), True, None),
    '0/.': ((0, None), None, None),
    # A member with another allele of a site of several. The split of the site into a record per allele writes the
    # other allele as the reference, so such a member reads as the first two genotypes above: 0/2 becomes 0/0, and 1/2
    # becomes 1/0.
    '1/0': ((1, 0), True, 'het'),
    # Left as the index of an allele that the record lacks, the other allele is read as missing.
    '0/2': ((0, None), None, None),
    '1/2': ((1, None), True, None),
}


def _samples(tmp_path: pathlib.Path, keys: str, values: Sequence[str]) -> list[pysam.VariantRecordSample]:
    """A record of one alternate allele with a sample for each of the values, under the FORMAT keys given."""
    names = [f'G{number}' for number in range(len(values))]
    lines = [
        '##fileformat=VCFv4.2',
        '##contig=<ID=chr1,length=1000>',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="x">',
        '##FORMAT=<ID=FT,Number=1,Type=String,Description="x">',
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *names]),
        '\t'.join(['chr1', '5', '.', 'C', 'T', '50', 'PASS', '.', keys, *values]),
    ]
    path = tmp_path / 'genotypes.vcf'
    path.write_text('\n'.join(lines) + '\n')
    with pysam.VariantFile(str(path)) as vcf:
        (record,) = vcf
        return [record.samples[name] for name in names]


def test_a_genotype_says_whether_it_carries_the_allele_only_where_its_alleles_decide_it(tmp_path: pathlib.Path) -> None:
    samples = _samples(tmp_path, 'GT', list(GENOTYPES))
    for sample, (written, (alleles, carries, zygosity)) in zip(samples, GENOTYPES.items(), strict=True):
        read, gt = vcf_values.genotype(sample, 'here')
        assert gt == written.replace('2', '.')
        assert read == alleles, written
        assert vcf_values.carries(read) is carries, written
        assert vcf_values.zygosity(read) == zygosity, written


@pytest.mark.parametrize('written', GENOTYPES)
def test_the_zygosity_never_says_more_than_the_genotype_determines(written: str) -> None:
    alleles, _, _ = GENOTYPES[written]
    by_carries = {None: {None}, False: {'hom_ref'}, True: {'het', 'hom_alt', None}}
    assert vcf_values.zygosity(alleles) in by_carries[vcf_values.carries(alleles)]
    # A missing allele alone never makes a genotype one that lacks the allele, and it leaves the zygosity open.
    if None in alleles:
        assert vcf_values.carries(alleles) is not False
    assert (vcf_values.zygosity(alleles) is None) == (None in alleles)
    # A member who was called and lacks the allele is told by either field.
    assert (vcf_values.carries(alleles) is False) == (vcf_values.zygosity(alleles) == 'hom_ref')


def test_a_genotype_that_pysam_reads_without_alleles_is_refused(tmp_path: pathlib.Path) -> None:
    # GT stands after GQ, against the VCF specification, and pysam then returns no alleles for it.
    (sample,) = _samples(tmp_path, 'GQ:GT', ['40:0/1'])
    assert sample['GT'] == ()
    with pytest.raises(ValueError, match='here: G0 has a genotype without alleles'):
        vcf_values.genotype(sample, 'here')


def test_a_calls_filters_are_a_list_of_names_and_a_call_without_them_has_none(tmp_path: pathlib.Path) -> None:
    # The last call ends before FT, which htslib reads as the missing value that the second call writes.
    samples = _samples(tmp_path, 'GT:FT', ['0/1:PASS', '0/1:.', '0/1:LowDepth;LowGQ', '0/1:LowGQ', '0/1'])
    formats = frozenset({'GT', 'FT'})
    read = [vcf_values.names(sample, 'FT', formats, 'here') for sample in samples]
    assert read == [['PASS'], None, ['LowDepth', 'LowGQ'], ['LowGQ'], None]
    # A file whose header does not define the field has none either.
    assert vcf_values.names(samples[0], 'FT', frozenset({'GT'}), 'here') is None


@pytest.mark.parametrize('written', ['LowGQ;', ';LowGQ', 'LowDepth;;LowGQ'])
def test_a_list_of_names_with_an_empty_name_is_refused(tmp_path: pathlib.Path, written: str) -> None:
    (sample,) = _samples(tmp_path, 'GT:FT', [f'0/1:{written}'])
    with pytest.raises(ValueError, match=r'here: FT .* holds an empty name'):
        vcf_values.names(sample, 'FT', frozenset({'GT', 'FT'}), 'here')
