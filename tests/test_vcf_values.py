"""What a genotype says of the allele of its record, on an invented record of one alternate allele."""

from __future__ import annotations

import pathlib

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
    # one may be the allele.
    './1': ((None, 1), True, 'het'),
    '0/.': ((0, None), None, None),
    # A member with another allele of a site of several. The split of the site into a record per allele writes the
    # other allele as the reference, so such a member reads as the first two genotypes above: 0/2 becomes 0/0, and 1/2
    # becomes 1/0.
    '1/0': ((1, 0), True, 'het'),
    # Left as the index of an allele that the record lacks, the other allele is read as missing.
    '0/2': ((0, None), None, None),
    '1/2': ((1, None), True, 'het'),
}


def _samples(tmp_path: pathlib.Path) -> list[pysam.VariantRecordSample]:
    """A record of one alternate allele with a sample for each genotype, in the order of GENOTYPES."""
    names = [f'G{number}' for number in range(len(GENOTYPES))]
    lines = [
        '##fileformat=VCFv4.2',
        '##contig=<ID=chr1,length=1000>',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">',
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *names]),
        '\t'.join(['chr1', '5', '.', 'C', 'T', '50', 'PASS', '.', 'GT', *GENOTYPES]),
    ]
    path = tmp_path / 'genotypes.vcf'
    path.write_text('\n'.join(lines) + '\n')
    with pysam.VariantFile(str(path)) as vcf:
        (record,) = vcf
        return [record.samples[name] for name in names]


def test_a_genotype_says_whether_it_carries_the_allele_only_where_its_alleles_decide_it(tmp_path: pathlib.Path) -> None:
    for sample, (written, (alleles, carries, zygosity)) in zip(_samples(tmp_path), GENOTYPES.items(), strict=True):
        read, _ = vcf_values.genotype(sample)
        assert read == alleles, written
        assert vcf_values.carries(read) is carries, written
        assert vcf_values.zygosity(read) == zygosity, written


@pytest.mark.parametrize('written', GENOTYPES)
def test_the_zygosity_never_says_more_than_whether_the_allele_is_carried(written: str) -> None:
    alleles, _, _ = GENOTYPES[written]
    by_carries = {None: {None}, False: {'hom_ref'}, True: {'het', 'hom_alt'}}
    assert vcf_values.zygosity(alleles) in by_carries[vcf_values.carries(alleles)]
    # A missing allele alone never makes a genotype one that lacks the allele.
    if None in alleles:
        assert vcf_values.carries(alleles) is not False
