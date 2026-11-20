"""A union's values from public files that are joined to its alleles, from an invented union and invented files.

Every allele, frequency, feature and score here is invented. The files have the layout of the public ones: the INFO
fields of gnomAD's sites VCF of the mitochondrial genome, the nine fields of GENCODE's GTF of polyadenylation features,
and the header of APARENT2's table of effects. The alleles on chrM have a REF that is not the reference's sequence at
their position, and those on the other chromosomes lie in the first 10,000 bases, where the reference has no sequence.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import json
import pathlib
import struct
from collections.abc import Mapping, Sequence
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from case_annotation import cli, joined_values

type Allele = tuple[str, int, str, str]


def _float32(value: float) -> float:
    return struct.unpack('f', struct.pack('f', value))[0]


def _vid(allele: Allele) -> str:
    return '-'.join(str(part) for part in allele)


def _union(path: pathlib.Path, alleles: Sequence[Allele], samples: Sequence[str] = ()) -> pathlib.Path:
    """A sites VCF of the alleles, as the union is: one allele per record and no samples."""
    columns = [
        '#CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO',
        *(['FORMAT', *samples] if samples else []),
    ]
    lines = [
        '##fileformat=VCFv4.2',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        *(f'##contig=<ID={chrom}>' for chrom in dict.fromkeys(a[0] for a in alleles)),
        '\t'.join(columns),
        *(
            '\t'.join(
                [chrom, str(pos), '.', ref, alt, '.', '.', '.', *(['GT', *('0/1' for _ in samples)] if samples else [])]
            )
            for chrom, pos, ref, alt in alleles
        ),
    ]
    path.write_text('\n'.join(lines) + '\n')
    return path


def _convert(
    tmp_path: pathlib.Path,
    alleles: Sequence[Allele],
    files: Mapping[str, pathlib.Path],
    normalised: Mapping[str, pathlib.Path] | None = None,
) -> dict[str, dict[str, Any]]:
    """The table that the conversion writes for a union of the alleles, as rows by variant id.

    A file that is a VCF of alleles is given as its own normalised copy unless another is named.
    """
    sites = _union(tmp_path / 'union.vcf', alleles)
    if normalised is None:
        normalised = {name: file for name, file in files.items() if joined_values.SOURCES[name].normalised}
    joined_values.convert(sites, files, normalised, tmp_path / 'out')
    return {row['variant_id']: row for row in pq.read_table(tmp_path / 'out' / joined_values.TABLE_NAME).to_pylist()}


# ---- the bases that an allele changes ----

SPANS = {
    'a single-base change': ((100, 'A', 'G'), (99, 100)),
    'a change of two bases': ((100, 'AC', 'GT'), (99, 101)),
    'a deletion, without the base written before it': ((100, 'TCA', 'T'), (100, 102)),
    'an insertion, between its base and the next': ((100, 'T', 'TCA'), (100, 100)),
    'a deletion with an insertion': ((100, 'TCA', 'TG'), (100, 102)),
    'an indel that shares no first base': ((100, 'TC', 'G'), (99, 101)),
}


@pytest.mark.parametrize('case', SPANS)
def test_an_allele_changes_the_reference_bases_after_what_ref_and_alt_share(case: str) -> None:
    allele, span = SPANS[case]
    assert joined_values.changed_span(*allele) == span


# ---- gnomAD's mitochondrial sites ----

# The INFO fields that the columns come from, declared as gnomAD's file declares them.
_MITO_DEFINITIONS: dict[str, tuple[str, str]] = {
    'AN': ('1', 'Integer'),
    'AC_hom': ('1', 'Integer'),
    'AC_het': ('1', 'Integer'),
    'AF_hom': ('1', 'Float'),
    'AF_het': ('1', 'Float'),
    'max_hl': ('1', 'Float'),
    'hap_AF_hom': ('1', 'String'),
    'hap_AF_het': ('1', 'String'),
    'hapmax_AF_hom': ('1', 'String'),
    'hapmax_AF_het': ('1', 'String'),
    'faf_hapmax_hom': ('1', 'Float'),
    'hap_defining_variant': ('0', 'Flag'),
    'mitotip_score': ('1', 'Float'),
    'mitotip_trna_prediction': ('1', 'String'),
}
# What every record of gnomAD's file holds, with invented values. The lists by haplogroup are of three invented
# haplogroups, of which the third has no sample with a call.
_MITO_ALWAYS: dict[str, str] = {
    'AN': '50000',
    'AC_hom': '12',
    'AC_het': '3',
    'AF_hom': '0.00024',
    'AF_het': '6.00000e-05',
    'max_hl': '0.998000',
    'hap_AF_hom': '1.00000e-04|4.00000e-03|NaN',
    'hap_AF_het': '0.00000e+00|0.00000e+00|NaN',
    'faf_hapmax_hom': '0.00310000',
}


@dataclasses.dataclass(frozen=True)
class MitoRecord:
    """A record of the invented file: an allele, its FILTER, and the INFO fields it holds beside the usual ones.

    A field of the usual ones that is given as None is left out.
    """

    pos: int
    ref: str
    alt: str
    filters: str = 'PASS'
    info: Mapping[str, str | None] = dataclasses.field(default_factory=dict)
    chrom: str = 'chrM'

    def line(self) -> str:
        fields = {**_MITO_ALWAYS, **self.info}
        info = ';'.join(
            name if value == '' else f'{name}={value}' for name, value in fields.items() if value is not None
        )
        return '\t'.join([self.chrom, str(self.pos), '.', self.ref, self.alt, '.', self.filters, info])


def _mito(
    path: pathlib.Path,
    records: Sequence[MitoRecord],
    definitions: Mapping[str, tuple[str, str]] = _MITO_DEFINITIONS,
    samples: bool = False,
) -> pathlib.Path:
    lines = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=npg,Description="Invented filter">',
        '##FILTER=<ID=indel_stack,Description="Invented filter">',
        *(
            f'##INFO=<ID={name},Number={number},Type={kind},Description="Invented">'
            for name, (number, kind) in definitions.items()
        ),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##contig=<ID=chrM,length=16569>',
        '##contig=<ID=chr1,length=248956422>',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' + ('\tFORMAT\tS1' if samples else ''),
        *(record.line() + ('\tGT\t0/1' if samples else '') for record in records),
    ]
    path.write_text('\n'.join(lines) + '\n')
    return path


_MITO_RECORDS = (
    MitoRecord(3, 'G', 'C'),
    MitoRecord(
        7,
        'C',
        'G',
        'indel_stack;npg',
        {
            'AC_hom': '0',
            'AF_hom': '0.00000',
            'hapmax_AF_hom': 'L3',
            'hap_defining_variant': '',
            'mitotip_score': '12.8668',
            'mitotip_trna_prediction': 'possibly_benign',
        },
    ),
    # A deletion at its left-aligned position, and an insertion.
    MitoRecord(20, 'TCA', 'T', info={'hapmax_AF_het': 'H', 'hap_AF_hom': 'NaN|NaN|NaN'}),
    MitoRecord(30, 'G', 'GAA'),
)
_MITO_COLUMNS = joined_values.GNOMAD_MITO.columns


def test_an_allele_on_chrm_has_gnomads_values_as_the_file_gives_them(tmp_path: pathlib.Path) -> None:
    alleles = [
        ('chrM', 3, 'G', 'C'),
        ('chrM', 7, 'C', 'G'),
        ('chrM', 20, 'TCA', 'T'),
        ('chrM', 30, 'G', 'GAA'),
        # Another allele of a position that the file holds, the same deletion written two bases on, and the file's
        # first allele on another chromosome: the file lacks each.
        ('chrM', 3, 'G', 'T'),
        ('chrM', 22, 'ACA', 'A'),
        ('chr1', 3, 'G', 'C'),
    ]
    rows = _convert(tmp_path, alleles, {'gnomad_mito': _mito(tmp_path / 'mito.vcf', _MITO_RECORDS)})
    assert {vid: {c: row[c] for c in _MITO_COLUMNS} for vid, row in rows.items() if vid == 'chrM-3-G-C'} == {
        'chrM-3-G-C': {
            'gnomad_mito_ac_hom': 12,
            'gnomad_mito_ac_het': 3,
            'gnomad_mito_af_hom': _float32(0.00024),
            'gnomad_mito_af_het': _float32(0.00006),
            'gnomad_mito_an': 50000,
            'gnomad_mito_max_heteroplasmy': _float32(0.998),
            # The highest of the list by haplogroup, without the haplogroup that has no value.
            'gnomad_mito_hapmax_af_hom': _float32(0.004),
            'gnomad_mito_hapmax_af_het': 0.0,
            'gnomad_mito_hapmax_hom_haplogroup': None,
            'gnomad_mito_hapmax_het_haplogroup': None,
            'gnomad_mito_faf_hapmax_hom': _float32(0.0031),
            'gnomad_mito_filters': ['PASS'],
            'gnomad_mito_hap_defining': False,
            'gnomad_mito_mitotip_score': None,
            'gnomad_mito_mitotip_prediction': None,
        }
    }
    filtered = rows['chrM-7-C-G']
    assert filtered['gnomad_mito_filters'] == ['indel_stack', 'npg']
    assert (filtered['gnomad_mito_ac_hom'], filtered['gnomad_mito_af_hom']) == (0, 0.0)
    assert (filtered['gnomad_mito_hapmax_hom_haplogroup'], filtered['gnomad_mito_hap_defining']) == ('L3', True)
    assert filtered['gnomad_mito_mitotip_score'] == _float32(12.8668)
    assert filtered['gnomad_mito_mitotip_prediction'] == 'possibly_benign'
    # An indel is matched by its position and both alleles, as the file writes it after normalising.
    deletion = rows['chrM-20-TCA-T']
    assert (deletion['gnomad_mito_an'], deletion['gnomad_mito_hapmax_het_haplogroup']) == (50000, 'H')
    # No haplogroup has a value in its list, so there is no highest one.
    assert deletion['gnomad_mito_hapmax_af_hom'] is None
    assert rows['chrM-30-G-GAA']['gnomad_mito_an'] == 50000
    for lacking in ('chrM-3-G-T', 'chrM-22-ACA-A', 'chr1-3-G-C'):
        assert {rows[lacking][c] for c in _MITO_COLUMNS} == {None}


def test_gnomads_sites_are_read_from_the_normalised_copy_and_recorded_as_given(tmp_path: pathlib.Path) -> None:
    # The file as given writes the deletion at the right end of a repeat, and its normalised copy at the left end.
    given = _mito(tmp_path / 'given.vcf', [MitoRecord(22, 'ACA', 'A')])
    normalised = _mito(tmp_path / 'normalised.vcf', [MitoRecord(20, 'TCA', 'T')])
    alleles = [('chrM', 20, 'TCA', 'T'), ('chrM', 22, 'ACA', 'A')]
    rows = _convert(tmp_path, alleles, {'gnomad_mito': given}, {'gnomad_mito': normalised})
    assert rows['chrM-20-TCA-T']['gnomad_mito_an'] == 50000
    assert rows['chrM-22-ACA-A']['gnomad_mito_an'] is None
    (record,) = joined_values.read_info(tmp_path / 'out').files
    assert record == joined_values.JoinedFile(
        source='gnomad_mito',
        file='given.vcf',
        bytes=given.stat().st_size,
        sha256=hashlib.sha256(given.read_bytes()).hexdigest(),
    )


def _without(name: str) -> dict[str, tuple[str, str]]:
    return {field: declared for field, declared in _MITO_DEFINITIONS.items() if field != name}


MITO_REFUSED: dict[str, tuple[dict[str, Any], str]] = {
    'a header without one of the INFO fields': (
        {'records': _MITO_RECORDS[:1], 'definitions': _without('mitotip_score')},
        'no INFO field mitotip_score',
    ),
    'an INFO field declared as another type': (
        {'records': _MITO_RECORDS[:1], 'definitions': {**_MITO_DEFINITIONS, 'hapmax_AF_hom': ('1', 'Float')}},
        r"INFO hapmax_AF_hom is declared as \('Float', 1\), not as \('String', 1\)",
    ),
    'an INFO field declared with another number of values': (
        {'records': _MITO_RECORDS[:1], 'definitions': {**_MITO_DEFINITIONS, 'AN': ('.', 'Integer')}},
        'INFO AN is declared as',
    ),
    'a VCF with samples': ({'records': _MITO_RECORDS[:1], 'samples': True}, 'this VCF has some'),
    'a record on another chromosome': ({'records': [MitoRecord(3, 'G', 'C', chrom='chr1')]}, 'chr1:3: not on chrM'),
    'a record of two alleles': (
        {'records': [MitoRecord(3, 'G', 'C,T')]},
        'chrM:3 does not hold one reference and one alternate allele',
    ),
    'two records of one allele': (
        {'records': [MitoRecord(3, 'G', 'C'), MitoRecord(3, 'G', 'C')]},
        'chrM:3: a second record of the allele chrM-3-G-C',
    ),
    'a record without a field that every record has': (
        {'records': [MitoRecord(3, 'G', 'C', info={'AN': None})]},
        'chrM:3: no AN, which every record',
    ),
    'a count that is no number': (
        {'records': [MitoRecord(3, 'G', 'C', info={'AC_hom': 'many'})]},
        'chrM:3: no AC_hom, which every record',
    ),
    'a record without a list by haplogroup': (
        {'records': [MitoRecord(3, 'G', 'C', info={'hap_AF_het': None})]},
        'chrM:3: no hap_AF_het, which every record',
    ),
    'a list by haplogroup that holds no numbers': (
        {'records': [MitoRecord(3, 'G', 'C', info={'hap_AF_hom': '0.1|high'})]},
        r"chrM:3: hap_AF_hom is '0.1\|high', not numbers",
    ),
}


@pytest.mark.parametrize('case', MITO_REFUSED)
def test_a_file_that_is_not_gnomads_mitochondrial_sites_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    written, message = MITO_REFUSED[case]
    file = _mito(tmp_path / 'mito.vcf', **written)
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, [('chrM', 3, 'G', 'C')], {'gnomad_mito': file})
    assert not (tmp_path / 'out' / joined_values.INFO_NAME).exists()


def test_a_line_that_is_no_record_is_refused(tmp_path: pathlib.Path) -> None:
    file = _mito(tmp_path / 'mito.vcf', _MITO_RECORDS)
    file.write_text(file.read_text() + 'chrM\t40\t.\tT\n')
    with pytest.raises(ValueError, match=r'mito\.vcf: a record that does not parse'):
        _convert(tmp_path, [('chrM', 3, 'G', 'C')], {'gnomad_mito': file})


def test_a_file_that_is_no_vcf_is_refused_as_gnomads_sites(tmp_path: pathlib.Path) -> None:
    file = tmp_path / 'mito.vcf'
    file.write_text('chrom\tpos\nchrM\t3\n')
    with pytest.raises(ValueError, match=r'mito\.vcf: not a VCF'):
        _convert(tmp_path, [('chrM', 3, 'G', 'C')], {'gnomad_mito': file})


# ---- GENCODE's polyadenylation features ----

type Feature = tuple[str, str, int, int]

# An invented signal of six bases on each strand, a site, and a pseudogene's signal.
_FEATURES: tuple[Feature, ...] = (
    ('chr1', 'polyA_signal', 1001, 1006),
    ('chr1', 'polyA_site', 1020, 1021),
    ('chr1', 'pseudo_polyA', 1100, 1105),
    ('chr2', 'polyA_signal', 500, 505),
)


def _gtf_line(feature: Feature, strand: str = '+') -> str:
    chrom, kind, start, end = feature
    attributes = f'gene_id "{start}"; transcript_id "{start}"; gene_type "{kind.lower()}"; level 2;'
    return '\t'.join([chrom, 'HAVANA', kind, str(start), str(end), '.', strand, '.', attributes])


def _gtf(path: pathlib.Path, features: Sequence[Feature] = _FEATURES, extra: Sequence[str] = ()) -> pathlib.Path:
    """A GTF of the features with comment lines before them, gzipped where the path's name says so."""
    lines = ['##description: invented polyA features', '##format: gtf', *(_gtf_line(f) for f in features), *extra]
    text = ('\n'.join(lines) + '\n').encode()
    path.write_bytes(gzip.compress(text) if path.suffix == '.gz' else text)
    return path


# Alleles around the signal at chr1:1001-1006, each with whether it overlaps the signal.
OVERLAPS: dict[str, tuple[Allele, bool]] = {
    'a change of its first base': (('chr1', 1001, 'A', 'G'), True),
    'a change of its last base': (('chr1', 1006, 'A', 'G'), True),
    'a change of the base before it': (('chr1', 1000, 'A', 'G'), False),
    'a change of the base after it': (('chr1', 1007, 'A', 'G'), False),
    'a change of two bases that ends in it': (('chr1', 1000, 'AC', 'GT'), True),
    'a deletion of its first base, written from the base before': (('chr1', 1000, 'TA', 'T'), True),
    'a deletion after it, written from its last base': (('chr1', 1006, 'AC', 'A'), False),
    'a deletion that ends before it': (('chr1', 998, 'TCA', 'T'), False),
    'a deletion of all of it and more': (('chr1', 990, 'T' * 40, 'T'), True),
    'an insertion between two of its bases': (('chr1', 1003, 'A', 'ACC'), True),
    'an insertion after its first base': (('chr1', 1001, 'A', 'AC'), True),
    'an insertion before its last base': (('chr1', 1005, 'A', 'AC'), True),
    'an insertion before its first base': (('chr1', 1000, 'T', 'TC'), False),
    'an insertion after its last base': (('chr1', 1006, 'A', 'AC'), False),
    'a change in a site, which is no signal': (('chr1', 1020, 'A', 'G'), False),
    "a change in a pseudogene's signal": (('chr1', 1102, 'A', 'G'), False),
    'a change at its position on a chromosome without one': (('chr3', 1003, 'A', 'G'), False),
    'a change in the signal of another chromosome': (('chr2', 502, 'A', 'G'), True),
    'an allele on chrM': (('chrM', 1003, 'C', 'G'), False),
}


@pytest.mark.parametrize('name', ['polya.gtf', 'polya.gtf.gz'])
def test_an_allele_overlaps_a_signal_that_holds_a_base_it_changes(tmp_path: pathlib.Path, name: str) -> None:
    alleles = [allele for allele, _ in OVERLAPS.values()]
    rows = _convert(tmp_path, alleles, {'gencode_polya': _gtf(tmp_path / name)})
    found = {case: rows[_vid(allele)]['gencode_polya_signal'] for case, (allele, _) in OVERLAPS.items()}
    assert found == {case: overlaps for case, (_, overlaps) in OVERLAPS.items()}
    # The flag is false and never null for an allele that overlaps none.
    schema = pq.read_schema(tmp_path / 'out' / joined_values.TABLE_NAME)
    assert (schema.field('gencode_polya_signal').type, schema.field('gencode_polya_signal').nullable) == (
        pa.bool_(),
        False,
    )


def test_an_overlap_is_found_across_the_stretches_that_the_join_looks_in(tmp_path: pathlib.Path) -> None:
    # A signal over the end of a stretch of 1,024 bases, one long feature, and the same signal twice.
    features: list[Feature] = [
        ('chr1', 'polyA_signal', 2046, 2051),
        ('chr1', 'polyA_signal', 3000, 4000),
        ('chr1', 'polyA_signal', 9500, 9505),
        ('chr1', 'polyA_signal', 9500, 9505),
    ]
    alleles: dict[Allele, bool] = {
        ('chr1', 2046, 'A', 'G'): True,
        ('chr1', 2050, 'A', 'G'): True,
        ('chr1', 2048, 'A', 'AC'): True,
        ('chr1', 2052, 'A', 'G'): False,
        ('chr1', 3500, 'A', 'G'): True,
        # A deletion that starts three stretches before a signal and ends in it.
        ('chr1', 6400, 'T' * 3102, 'T'): True,
        ('chr1', 9505, 'A', 'G'): True,
        ('chr1', 9900, 'A', 'G'): False,
    }
    rows = _convert(tmp_path, list(alleles), {'gencode_polya': _gtf(tmp_path / 'polya.gtf', features)})
    assert {allele: rows[_vid(allele)]['gencode_polya_signal'] for allele in alleles} == alleles
    assert len(rows) == len(alleles)


GTF_REFUSED: dict[str, tuple[dict[str, Any], str]] = {
    'a line of eight fields': (
        {'extra': ['\t'.join(_gtf_line(_FEATURES[0]).split('\t')[:8])]},
        "not GENCODE's polyadenylation features, a GTF",
    ),
    'a feature of another kind': (
        {'extra': [_gtf_line(('chr1', 'exon', 3000, 3100))]},
        r"the features \['exon'\], and GENCODE's polyadenylation features are 'polyA_signal'",
    ),
    'a start that is no position': (
        {'extra': [_gtf_line(('chr1', 'polyA_signal', 3000, 3005)).replace('3000', '3e3')]},
        r"a start and an end after it, such as \[\('chr1', 'polyA_signal', '3e3', '3005'\)",
    ),
    'a feature without a chromosome': (
        {'extra': ['\t' + _gtf_line(('chr1', 'polyA_signal', 3000, 3005)).split('\t', 1)[1]]},
        r"a start and an end after it, such as \[\(None, 'polyA_signal', '3000', '3005'\)",
    ),
    'a feature without a kind': (
        {'extra': [_gtf_line(('chr1', 'polyA_signal', 3000, 3005)).replace('polyA_signal', '')]},
        r'the features \[None\], and GENCODE',
    ),
    'a start of 0': (
        {'features': [*_FEATURES, ('chr1', 'polyA_site', 0, 1)]},
        'features without a chromosome, a start and an end after it',
    ),
    'an end before its start': (
        {'features': [*_FEATURES, ('chr1', 'polyA_signal', 3005, 3000)]},
        'features without a chromosome, a start and an end after it',
    ),
    'no signal on a primary chromosome': (
        {'features': [('1', 'polyA_signal', 1001, 1006), ('chr1', 'polyA_site', 1020, 1021)]},
        'no polyA_signal on a primary chromosome; are its chromosomes named chr1',
    ),
}


@pytest.mark.parametrize('case', GTF_REFUSED)
def test_a_file_that_is_not_gencodes_polyadenylation_features_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    written, message = GTF_REFUSED[case]
    file = _gtf(tmp_path / 'polya.gtf.gz', **written)
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, [('chr1', 1001, 'A', 'G')], {'gencode_polya': file})


# ---- APARENT2's effects ----

_APARENT2_HEADER = (
    'chrom\tstrand\tgene\tvar_position_hg19\tvar_position_hg38\tref_nucleotide\tvar_nucleotide\tpas_type\t'
    'pas_position_hg19\tpas_position_hg38\tref_usage\tdelta_logodds_narrow\tdelta_usage_narrow\tdelta_logodds\t'
    'delta_usage'
)


def _effect(
    chrom: str,
    pos: object,
    ref: str,
    alt: str,
    delta: object,
    usage: object,
    *,
    site: int = 5000,
    strand: str = '+',
    gene: str = 'GENEA',
) -> str:
    """A row of the invented table: a change, the site it is scored for, the site's usage and the change's effect."""
    fields = [chrom, strand, gene, 90_000, pos, ref, alt, '3_most_exon', site + 80_000, site, usage, delta]
    return '\t'.join(str(field) for field in (*fields, 0.001, 0.25, 0.002))


# A change with a row for each of two sites and one for a gene on the other strand, a change with one row, and a change
# whose two lowest rows are equal.
_EFFECTS = (
    _effect('chr1', 4990, 'A', 'G', -3.19127, 0.47273),
    _effect('chr1', 4990, 'A', 'G', 0.25, 0.06413, site=5100),
    _effect('chr1', 4990, 'A', 'G', -0.00352, 0.9, site=4900, strand='-', gene='GENEB'),
    _effect('chr1', 4990, 'A', 'C', 0.01185, 0.47273),
    _effect('chr2', 700, 'T', 'C', -1.5, 0.2),
    _effect('chr2', 700, 'T', 'C', -1.5, 0.6, site=800),
    _effect('chr2', 700, 'T', 'C', 2e-05, 0.99, site=900),
)


def _effects(path: pathlib.Path, rows: Sequence[str] = _EFFECTS, header: str = _APARENT2_HEADER) -> pathlib.Path:
    text = ('\n'.join([header, *rows]) + '\n').encode()
    path.write_bytes(gzip.compress(text) if path.suffix == '.gz' else text)
    return path


@pytest.mark.parametrize('name', ['effects.csv.gz', 'effects.csv'])
def test_a_single_base_change_has_the_lowest_and_the_highest_effect_of_its_rows(
    tmp_path: pathlib.Path, name: str
) -> None:
    alleles = [
        ('chr1', 4990, 'A', 'G'),
        ('chr1', 4990, 'A', 'C'),
        ('chr2', 700, 'T', 'C'),
        # Without a row: another change at a scored position, an insertion there, a change far from a site, and the
        # scored change on another chromosome.
        ('chr1', 4990, 'A', 'T'),
        ('chr1', 4990, 'A', 'AG'),
        ('chr1', 9000, 'A', 'G'),
        ('chr3', 4990, 'A', 'G'),
    ]
    rows = _convert(tmp_path, alleles, {'aparent2': _effects(tmp_path / name)})
    columns = joined_values.APARENT2.columns
    assert columns == ('aparent2_min', 'aparent2_max', 'aparent2_ref_usage')
    found = {vid: tuple(row[c] for c in columns) for vid, row in rows.items()}
    assert found == {
        # The usage is that of the site of the lowest row.
        'chr1-4990-A-G': (_float32(-3.19127), 0.25, _float32(0.47273)),
        'chr1-4990-A-C': (_float32(0.01185), _float32(0.01185), _float32(0.47273)),
        # Of two sites with the lowest effect, the more used one.
        'chr2-700-T-C': (-1.5, _float32(2e-05), _float32(0.6)),
        'chr1-4990-A-T': (None, None, None),
        'chr1-4990-A-AG': (None, None, None),
        'chr1-9000-A-G': (None, None, None),
        'chr3-4990-A-G': (None, None, None),
    }
    schema = pq.read_schema(tmp_path / 'out' / joined_values.TABLE_NAME)
    assert {schema.field(c).type for c in columns} == {pa.float32()}


EFFECTS_REFUSED: dict[str, tuple[dict[str, Any], str]] = {
    'a header with another column name': (
        {'header': _APARENT2_HEADER.replace('var_position_hg38', 'var_position')},
        r"its header is \['chrom', 'strand', 'gene', 'var_position_hg19', 'var_position'",
    ),
    'a header without its last column': (
        {'header': _APARENT2_HEADER.rsplit('\t', 1)[0]},
        "its header is .*, and APARENT2's is",
    ),
    'a row of fourteen fields': (
        {'rows': [*_EFFECTS, _EFFECTS[0].rsplit('\t', 1)[0]]},
        "not APARENT2's table of effects",
    ),
    'a position that is no whole number': (
        {'rows': [*_EFFECTS, _effect('chr9', '4990.5', 'A', 'G', 0.1, 0.2)]},
        r"such as \[\['chr9', '4990.5', 'A', 'G', '0.2', '0.1'\]\]",
    ),
    'a position of 0': ({'rows': [_effect('chr9', 0, 'A', 'G', 0.1, 0.2)]}, 'rows without'),
    'a reference that is no base': ({'rows': [_effect('chr9', 4990, 'N', 'G', 0.1, 0.2)]}, 'rows without'),
    'a changed nucleotide of two bases': ({'rows': [_effect('chr9', 4990, 'A', 'GG', 0.1, 0.2)]}, 'rows without'),
    'an effect that is no number': ({'rows': [_effect('chr9', 4990, 'A', 'G', 'low', 0.2)]}, 'rows without'),
    'an effect that is not finite': ({'rows': [_effect('chr9', 4990, 'A', 'G', 'nan', 0.2)]}, 'rows without'),
    'an infinite effect': ({'rows': [_effect('chr9', 4990, 'A', 'G', '-inf', 0.2)]}, 'rows without'),
    'an empty usage': ({'rows': [_effect('chr9', 4990, 'A', 'G', 0.1, '')]}, 'rows without'),
    'an empty chromosome': ({'rows': [_effect('', 4990, 'A', 'G', 0.1, 0.2)]}, 'rows without'),
}


@pytest.mark.parametrize('case', EFFECTS_REFUSED)
def test_a_file_that_is_not_aparent2s_effects_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    # The malformed rows are of alleles that the union lacks, so a refusal shows that every row of the file is checked.
    written, message = EFFECTS_REFUSED[case]
    file = _effects(tmp_path / 'effects.csv.gz', **written)
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, [('chr1', 4990, 'A', 'G')], {'aparent2': file})


def test_a_table_of_effects_that_is_cut_off_is_refused(tmp_path: pathlib.Path) -> None:
    whole = _effects(tmp_path / 'whole.csv.gz', [_effect('chr1', 1000 + i, 'A', 'G', i / 7, 0.5) for i in range(5000)])
    cut = tmp_path / 'effects.csv.gz'
    cut.write_bytes(whole.read_bytes()[: whole.stat().st_size // 2])
    with pytest.raises(ValueError, match=r"effects\.csv\.gz: not APARENT2's table of effects"):
        _convert(tmp_path, [('chr1', 4990, 'A', 'G')], {'aparent2': cut})


def test_a_table_of_effects_that_is_cut_off_within_its_header_is_refused(tmp_path: pathlib.Path) -> None:
    whole = _effects(tmp_path / 'whole.csv.gz')
    cut = tmp_path / 'effects.csv.gz'
    cut.write_bytes(whole.read_bytes()[:20])
    with pytest.raises(ValueError, match=r"effects\.csv\.gz: not APARENT2's table of effects"):
        _convert(tmp_path, [('chr1', 4990, 'A', 'G')], {'aparent2': cut})


# ---- the conversion ----

_ALLELES: tuple[Allele, ...] = (
    ('chr1', 1003, 'A', 'G'),
    ('chr1', 4990, 'A', 'G'),
    ('chr2', 60, 'C', 'CT'),
    ('chrM', 3, 'G', 'C'),
)


def _files(tmp_path: pathlib.Path) -> dict[str, pathlib.Path]:
    return {
        'gnomad_mito': _mito(tmp_path / 'mito.vcf', _MITO_RECORDS),
        'gencode_polya': _gtf(tmp_path / 'polya.gtf.gz'),
        'aparent2': _effects(tmp_path / 'effects.csv.gz'),
    }


def test_the_table_holds_each_allele_of_the_union_once_with_the_columns_of_the_files_given(
    tmp_path: pathlib.Path,
) -> None:
    files = _files(tmp_path)
    rows = _convert(tmp_path, _ALLELES, files)
    # The rows are in the union's order.
    assert list(rows) == [_vid(allele) for allele in _ALLELES]
    table = pq.read_table(tmp_path / 'out' / joined_values.TABLE_NAME)
    assert table.schema.equals(joined_values.schema(tuple(joined_values.SOURCES.values())))
    held = {vid: (row['gnomad_mito_an'], row['gencode_polya_signal'], row['aparent2_min']) for vid, row in rows.items()}
    assert held == {
        'chr1-1003-A-G': (None, True, None),
        'chr1-4990-A-G': (None, False, _float32(-3.19127)),
        'chr2-60-C-CT': (None, False, None),
        'chrM-3-G-C': (50000, False, None),
    }
    # Each file is recorded as the run read it, in the order of the columns.
    written = json.loads((tmp_path / 'out' / joined_values.INFO_NAME).read_text())
    assert written == {
        'files': [
            {
                'source': name,
                'file': files[name].name,
                'bytes': files[name].stat().st_size,
                'sha256': hashlib.sha256(files[name].read_bytes()).hexdigest(),
            }
            for name in joined_values.SOURCES
        ]
    }
    assert joined_values.read_info(tmp_path / 'out').given() == tuple(joined_values.SOURCES.values())


SOME_FILES = {
    'no file': (),
    "gnomAD's sites alone": ('gnomad_mito',),
    'the two files on polyadenylation': ('gencode_polya', 'aparent2'),
}


@pytest.mark.parametrize('case', SOME_FILES)
def test_a_run_without_a_file_has_none_of_its_columns(tmp_path: pathlib.Path, case: str) -> None:
    given = {name: file for name, file in _files(tmp_path).items() if name in SOME_FILES[case]}
    rows = _convert(tmp_path, _ALLELES, given)
    sources = [joined_values.SOURCES[name] for name in SOME_FILES[case]]
    schema = pq.read_schema(tmp_path / 'out' / joined_values.TABLE_NAME)
    assert schema.names == ['variant_id', *(column for source in sources for column in source.columns)]
    assert sorted(rows) == sorted(_vid(allele) for allele in _ALLELES)
    info = joined_values.read_info(tmp_path / 'out')
    assert [file.source for file in info.files] == list(SOME_FILES[case])
    assert info.given() == tuple(sources)


def test_the_columns_are_32_bit_numbers_text_flags_and_one_list_of_names() -> None:
    kinds = {field.name: field.type for source in joined_values.SOURCES.values() for field in source.fields}
    assert kinds.pop('gnomad_mito_filters') == pa.list_(pa.string())
    assert set(kinds.values()) == {pa.int32(), pa.float32(), pa.string(), pa.bool_()}
    # A source's name starts each of its columns, so that a column says which file it is from.
    for source in joined_values.SOURCES.values():
        assert all(column.startswith(source.name) for column in source.columns), source.name


def test_a_file_of_no_source_is_refused(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match=r"\['clinvar'\] are no files that the tables read, which are gnomad_mito"):
        _convert(tmp_path, _ALLELES, {'clinvar': tmp_path / 'clinvar.vcf'}, {})


def test_a_vcf_of_alleles_without_its_normalised_copy_is_refused(tmp_path: pathlib.Path) -> None:
    files = _files(tmp_path)
    with pytest.raises(
        ValueError, match=r"rewrote them are \['gnomad_mito'\], and the normalised copies given are of \[\]"
    ):
        _convert(tmp_path, _ALLELES, files, {})
    assert not (tmp_path / 'out').exists()


def test_a_normalised_copy_of_a_file_that_is_no_vcf_of_alleles_is_refused(tmp_path: pathlib.Path) -> None:
    files = _files(tmp_path)
    with pytest.raises(ValueError, match=r"given are of \['gencode_polya', 'gnomad_mito'\]"):
        _convert(tmp_path, _ALLELES, files, {'gnomad_mito': files['gnomad_mito'], 'gencode_polya': files['aparent2']})
    # Also the copy of a file that the run was not given.
    with pytest.raises(
        ValueError, match=r"rewrote them are \[\], and the normalised copies given are of \['gnomad_mito'\]"
    ):
        _convert(tmp_path, _ALLELES, {'aparent2': files['aparent2']}, {'gnomad_mito': files['gnomad_mito']})


def test_a_union_with_samples_is_refused(tmp_path: pathlib.Path) -> None:
    sites = _union(tmp_path / 'union.vcf', _ALLELES, samples=('S1',))
    with pytest.raises(ValueError, match='the union is sites-only, and this VCF has samples'):
        joined_values.convert(sites, {}, {}, tmp_path / 'out')


def test_a_union_record_of_two_alleles_is_refused(tmp_path: pathlib.Path) -> None:
    sites = _union(tmp_path / 'union.vcf', [('chr1', 5, 'A', 'G,T')])
    with pytest.raises(ValueError, match='chr1:5 does not hold one reference and one alternate allele'):
        joined_values.convert(sites, {}, {}, tmp_path / 'out')


def test_a_source_that_gives_an_allele_two_rows_is_refused(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def doubled(_con: duckdb.DuckDBPyConnection, _path: pathlib.Path) -> str:
        return 'SELECT variant_id, 1 AS copies FROM alleles UNION ALL SELECT variant_id, 2 FROM alleles'

    source = joined_values.Source(name='doubled', fields=(pa.field('copies', pa.int32()),), read=doubled)
    monkeypatch.setattr(joined_values, 'SOURCES', {'doubled': source})
    with pytest.raises(ValueError, match='the joined values repeats alleles'):
        _convert(tmp_path, _ALLELES, {'doubled': tmp_path / 'union.vcf'})
    assert not (tmp_path / 'out' / joined_values.INFO_NAME).exists()


def test_a_union_that_repeats_an_allele_is_refused(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='the union repeats alleles'):
        _convert(tmp_path, [*_ALLELES, _ALLELES[0]], {})


# ---- the command ----


def _command(tmp_path: pathlib.Path, files: Mapping[str, pathlib.Path]) -> list[str]:
    sites = _union(tmp_path / 'union.vcf', _ALLELES)
    options = [option for name, file in files.items() for option in ('--file', name, str(file))]
    return ['joined-values', '--sites', str(sites), *options, '--output', str(tmp_path / 'out')]


def test_the_command_joins_the_files_it_is_given(tmp_path: pathlib.Path) -> None:
    files = _files(tmp_path)
    scratch = tmp_path / 'scratch'
    scratch.mkdir()
    normalised = ['--normalised', 'gnomad_mito', str(files['gnomad_mito']), '--scratch', str(scratch)]
    cli.main([*_command(tmp_path, files), *normalised])
    assert pq.read_table(tmp_path / 'out' / joined_values.TABLE_NAME).num_rows == len(_ALLELES)
    assert [file.source for file in joined_values.read_info(tmp_path / 'out').files] == list(joined_values.SOURCES)
    # Nothing is left in the folder for temporary files.
    assert not list(scratch.iterdir())


def test_the_command_writes_the_alleles_alone_for_a_run_without_a_file(tmp_path: pathlib.Path) -> None:
    cli.main(_command(tmp_path, {}))
    assert pq.read_schema(tmp_path / 'out' / joined_values.TABLE_NAME).names == ['variant_id']
    assert joined_values.read_info(tmp_path / 'out').files == ()


@pytest.mark.parametrize('option', ['--file', '--normalised'])
def test_the_command_refuses_an_option_given_twice_for_one_source(tmp_path: pathlib.Path, option: str) -> None:
    files = _files(tmp_path)
    twice = [option, 'gnomad_mito', str(files['gnomad_mito']), option, 'gnomad_mito', str(files['gnomad_mito'])]
    with pytest.raises(SystemExit, match=f'{option} is given twice for gnomad_mito'):
        cli.main([*_command(tmp_path, {}), *twice])
