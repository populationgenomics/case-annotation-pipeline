"""A union's values from public files that are joined to its alleles, from an invented union and invented files.

Every allele, frequency, feature, window, region, transcript and score here is invented. The files have the layout of
the public ones: the INFO fields of gnomAD's sites VCF of the mitochondrial genome, the nine fields of GENCODE's GTF of
polyadenylation features, and the headers of APARENT2's table of effects, of MitoTIP's table of scores and of gnomAD's
tables of genomic constraint by window, of regions of missense constraint and of MPC scores, and the columns of the
catalogue of branch points as this package prepares it. The alleles on chrM have a REF that is not the reference's
sequence at their position, and those on the other chromosomes lie in the first 10,000 bases, where the reference has
no sequence. The transcripts have ids that Ensembl has not given.
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
import pysam
import pytest

from case_annotation import branch_points, cli, joined_values

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


def _pairs(tmp_path: pathlib.Path) -> dict[tuple[str, str], dict[str, Any]]:
    """The table by allele and transcript that the conversion wrote, as rows by variant id and transcript id."""
    rows = pq.read_table(tmp_path / 'out' / joined_values.TRANSCRIPTS_NAME).to_pylist()
    return {(row.pop('variant_id'), row.pop('transcript_id')): row for row in rows}


def _table(path: pathlib.Path, header: str, rows: Sequence[str]) -> pathlib.Path:
    """A tab-separated table with a header, gzipped where the path's name says so."""
    text = ('\n'.join([header, *rows]) + '\n').encode()
    path.write_bytes(gzip.compress(text) if path.suffix == '.gz' else text)
    return path


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
        }
    }
    filtered = rows['chrM-7-C-G']
    assert filtered['gnomad_mito_filters'] == ['indel_stack', 'npg']
    assert (filtered['gnomad_mito_ac_hom'], filtered['gnomad_mito_af_hom']) == (0, 0.0)
    assert (filtered['gnomad_mito_hapmax_hom_haplogroup'], filtered['gnomad_mito_hap_defining']) == ('L3', True)
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
        {'records': _MITO_RECORDS[:1], 'definitions': _without('max_hl')},
        'no INFO field max_hl',
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
    return _table(path, header, rows)


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


# ---- MitoTIP's scores ----

_MITOTIP_HEADER = 'Position\trCRS\tAlt\tMitoTIP_Score\tQuartile\tCount\tPercentage\tMitomap_Status'


def _scored(pos: object, ref: str, alt: str, score: object, quartile: str = 'Q2') -> str:
    """A row of the invented table: a change of a position with its score, and nothing that a database says of it."""
    return '\t'.join(str(field) for field in (pos, ref, alt, score, quartile, 0, '0.000', 'Absent'))


# Two invented positions, each with the deletion of its base and its three changes, as the file has every position.
_SCORES = (
    _scored(3, 'G', ':', 20.3833, 'Q1'),
    _scored(3, 'G', 'A', 17.75, 'Q1'),
    _scored(3, 'G', 'C', 12.8668, 'Q2'),
    _scored(3, 'G', 'T', -5.89285, 'Q4'),
    _scored(7, 'C', ':', 9.5, 'Q3'),
    _scored(7, 'C', 'A', 8.44301, 'Q3'),
    _scored(7, 'C', 'G', 0, 'Q4'),
    _scored(7, 'C', 'T', 16.2594, 'Q1'),
)


def _mitotip(path: pathlib.Path, rows: Sequence[str] = _SCORES, header: str = _MITOTIP_HEADER) -> pathlib.Path:
    return _table(path, header, rows)


@pytest.mark.parametrize('name', ['mitotip_scores.txt', 'mitotip_scores.txt.gz'])
def test_a_single_base_change_on_chrm_has_mitotips_score_and_quartile(tmp_path: pathlib.Path, name: str) -> None:
    alleles = [
        ('chrM', 3, 'G', 'C'),
        ('chrM', 3, 'G', 'T'),
        ('chrM', 7, 'C', 'G'),
        # Without a row: a change at a position that the file lacks, a change with another REF than the file's, the
        # deletion of a scored base as a VCF writes it, an insertion at one, and a scored change on another
        # chromosome.
        ('chrM', 5, 'G', 'C'),
        ('chrM', 7, 'A', 'G'),
        ('chrM', 2, 'AG', 'A'),
        ('chrM', 3, 'G', 'GC'),
        ('chr1', 3, 'G', 'C'),
    ]
    rows = _convert(tmp_path, alleles, {'mitotip': _mitotip(tmp_path / name)})
    assert joined_values.MITOTIP.columns == ('mitotip_score', 'mitotip_quartile')
    assert {vid: (row['mitotip_score'], row['mitotip_quartile']) for vid, row in rows.items()} == {
        'chrM-3-G-C': (_float32(12.8668), 'Q2'),
        'chrM-3-G-T': (_float32(-5.89285), 'Q4'),
        # A score of 0 is a score.
        'chrM-7-C-G': (0.0, 'Q4'),
        'chrM-5-G-C': (None, None),
        'chrM-7-A-G': (None, None),
        'chrM-2-AG-A': (None, None),
        'chrM-3-G-GC': (None, None),
        'chr1-3-G-C': (None, None),
    }
    assert not _pairs(tmp_path)


MITOTIP_REFUSED: dict[str, tuple[dict[str, Any], str]] = {
    'a header with another column name': (
        {'header': _MITOTIP_HEADER.replace('MitoTIP_Score', 'Score')},
        r"its header is \['Position', 'rCRS', 'Alt', 'Score'.*, and MitoTIP's is \['Position'",
    ),
    'a row of seven fields': ({'rows': [*_SCORES, _SCORES[0].rsplit('\t', 1)[0]]}, "not MitoTIP's table of scores"),
    'a position that is no whole number': (
        {'rows': [*_SCORES, _scored('9.5', 'G', 'A', 1.5)]},
        r"rows without Position, rCRS, Alt, MitoTIP_Score, Quartile as .*such as \[\['9.5', 'G', 'A', '1.5', 'Q2'\]\]",
    ),
    'a position of 0': ({'rows': [_scored(0, 'G', 'A', 1.5)]}, 'rows without'),
    'a reference that is no base': ({'rows': [_scored(9, 'N', 'A', 1.5)]}, 'rows without'),
    'a change of two bases': ({'rows': [_scored(9, 'G', 'AC', 1.5)]}, 'rows without'),
    'a deletion written otherwise': ({'rows': [_scored(9, 'G', '-', 1.5)]}, 'rows without'),
    'a change to the reference base': ({'rows': [_scored(9, 'G', 'G', 1.5)]}, 'rows without'),
    'an empty change': ({'rows': [_scored(9, 'G', '', 1.5)]}, 'rows without'),
    'a score that is no number': ({'rows': [_scored(9, 'G', 'A', 'high')]}, 'rows without'),
    'a score that is not finite': ({'rows': [_scored(9, 'G', 'A', 'inf')]}, 'rows without'),
    'a quartile that is none of the four': ({'rows': [_scored(9, 'G', 'A', 1.5, 'Q5')]}, 'rows without'),
    'two rows of one change': (
        {'rows': [*_SCORES, _scored(7, 'C', 'T', 2.5)]},
        r"several rows for one change of a position, such as \[\(7, 'C', 'T'\)\]",
    ),
    'two rows of one deletion': ({'rows': [*_SCORES, _scored(7, 'C', ':', 2.5)]}, 'several rows for one change'),
}


@pytest.mark.parametrize('case', MITOTIP_REFUSED)
def test_a_file_that_is_not_mitotips_scores_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    written, message = MITOTIP_REFUSED[case]
    file = _mitotip(tmp_path / 'mitotip_scores.txt', **written)
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, [('chrM', 3, 'G', 'C')], {'mitotip': file})
    assert not (tmp_path / 'out' / joined_values.INFO_NAME).exists()


# ---- gnomAD's genomic constraint by window ----

_WINDOWS_HEADER = 'chrom\tstart\tend\telement_id\tpossible\texpected\tobserved\toe\tz'


def _window(chrom: str, start: object, end: object, z: object) -> str:
    """A row of the invented table: a window, 0-based and half-open, with its z score."""
    fields = (chrom, start, end, f'{chrom}-{start}-{end}', '2443.0', '159.09', '142.0', '0.89', z)
    return '\t'.join(str(field) for field in fields)


# Two neighbouring windows, one after a window that the file lacks, a window that ends before its 1,000 bases, the
# first window of a chromosome, and a window at the same place on another chromosome.
_WINDOWS = (
    _window('chr1', 2000, 3000, 1.3551011656250318),
    _window('chr1', 3000, 4000, -0.5),
    _window('chr1', 5000, 6000, 4.25),
    _window('chr1', 7000, 7400, 0),
    _window('chr2', 0, 1000, -9.95674792768871),
    _window('chr2', 2000, 3000, 2.18),
)


def _windows(path: pathlib.Path, rows: Sequence[str] = _WINDOWS, header: str = _WINDOWS_HEADER) -> pathlib.Path:
    return _table(path, header, rows)


# Alleles by position, each with the z score of its window: the window of the 1-based position p starts at the
# multiple of 1,000 at or below p - 1, so a position that is a multiple of 1,000 is the last base of a window.
WINDOWED: dict[str, tuple[Allele, float | None]] = {
    'the last base before a window': (('chr1', 2000, 'A', 'G'), None),
    'the first base of a window': (('chr1', 2001, 'A', 'G'), _float32(1.3551011656250318)),
    'a base within a window': (('chr1', 2500, 'A', 'G'), _float32(1.3551011656250318)),
    'the last base of a window, at a multiple of 1,000': (('chr1', 3000, 'A', 'G'), _float32(1.3551011656250318)),
    'the first base of the next window': (('chr1', 3001, 'A', 'G'), -0.5),
    'the last base of that window': (('chr1', 4000, 'A', 'G'), -0.5),
    'the first base after it, in a window that the file lacks': (('chr1', 4001, 'A', 'G'), None),
    'the last base of a window that the file lacks': (('chr1', 5000, 'A', 'G'), None),
    'the first base of a window after one that the file lacks': (('chr1', 5001, 'A', 'G'), 4.25),
    'a deletion, by the position that the VCF writes it at': (('chr1', 3000, 'ATT', 'A'), _float32(1.3551011656250318)),
    'an insertion, by its position': (('chr1', 3001, 'A', 'ACC'), -0.5),
    'the last base of a window that ends early, with a score of 0': (('chr1', 7400, 'A', 'G'), 0.0),
    'the base after a window that ends early': (('chr1', 7401, 'A', 'G'), None),
    'the first base of a chromosome': (('chr2', 1, 'A', 'G'), _float32(-9.95674792768871)),
    'the last base of the first window': (('chr2', 1000, 'A', 'G'), _float32(-9.95674792768871)),
    'a base of the same window on another chromosome': (('chr2', 2500, 'A', 'G'), _float32(2.18)),
    'a base on a chromosome without a window': (('chrX', 2500, 'A', 'G'), None),
    'a base on chrM': (('chrM', 2500, 'A', 'G'), None),
}


@pytest.mark.parametrize('name', ['windows.txt.gz', 'windows.txt'])
def test_an_allele_has_the_z_score_of_the_window_that_holds_its_position(tmp_path: pathlib.Path, name: str) -> None:
    alleles = [allele for allele, _ in WINDOWED.values()]
    rows = _convert(tmp_path, alleles, {'genomic_constraint': _windows(tmp_path / name)})
    found = {case: rows[_vid(allele)]['genomic_constraint_z'] for case, (allele, _) in WINDOWED.items()}
    assert found == {case: z for case, (_, z) in WINDOWED.items()}
    assert (
        pq.read_schema(tmp_path / 'out' / joined_values.TABLE_NAME).field('genomic_constraint_z').type == pa.float32()
    )


WINDOWS_REFUSED: dict[str, tuple[dict[str, Any], str]] = {
    'a header with another column name': (
        {'header': _WINDOWS_HEADER.replace('\tz', '\tz_score')},
        r"its header is \['chrom', 'start', 'end', .*'z_score'\], and gnomAD's is \['chrom'",
    ),
    'a row of eight fields': (
        {'rows': [*_WINDOWS, _WINDOWS[0].rsplit('\t', 1)[0]]},
        "not gnomAD's genomic constraint by window",
    ),
    'a chromosome without its prefix': (
        {'rows': [*_WINDOWS, _window('1', 8000, 9000, 1.5)]},
        r"rows without chrom, start, end, z as .*such as \[\['1', '8000', '9000', '1.5'\]\]",
    ),
    'a contig that is no primary chromosome': (
        {'rows': [_window('chr1_KI270706v1_random', 8000, 9000, 1.5)]},
        'rows without',
    ),
    'a start that is no multiple of 1,000': ({'rows': [_window('chr1', 8500, 9500, 1.5)]}, 'rows without'),
    'a start that is no whole number': ({'rows': [_window('chr1', '8e3', 9000, 1.5)]}, 'rows without'),
    'a start below 0': ({'rows': [_window('chr1', -1000, 0, 1.5)]}, 'rows without'),
    'an end at its start': ({'rows': [_window('chr1', 8000, 8000, 1.5)]}, 'rows without'),
    'an end more than 1,000 bases on': ({'rows': [_window('chr1', 8000, 9001, 1.5)]}, 'rows without'),
    'an empty end': ({'rows': [_window('chr1', 8000, '', 1.5)]}, 'rows without'),
    'a z score that is no number': ({'rows': [_window('chr1', 8000, 9000, 'NA')]}, 'rows without'),
    'a z score that is not finite': ({'rows': [_window('chr1', 8000, 9000, 'nan')]}, 'rows without'),
    'two rows of one window': (
        {'rows': [*_WINDOWS, _window('chr1', 5000, 6000, 1.5)]},
        r"several rows for one window, such as \[\('chr1', 5000\)\]",
    ),
    'two rows of one window with different ends': (
        {'rows': [*_WINDOWS, _window('chr1', 5000, 5500, 1.5)]},
        'several rows for one window',
    ),
}


@pytest.mark.parametrize('case', WINDOWS_REFUSED)
def test_a_file_that_is_not_gnomads_genomic_constraint_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    # The malformed rows are of windows that hold no allele of the union: every row of the file is checked.
    written, message = WINDOWS_REFUSED[case]
    file = _windows(tmp_path / 'windows.txt.gz', **written)
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, [('chr1', 2500, 'A', 'G')], {'genomic_constraint': file})
    assert not (tmp_path / 'out' / joined_values.INFO_NAME).exists()


def test_a_table_of_windows_that_is_cut_off_is_refused(tmp_path: pathlib.Path) -> None:
    whole = _windows(
        tmp_path / 'whole.txt.gz', [_window('chr1', 1000 * i, 1000 * i + 1000, i / 7) for i in range(5000)]
    )
    cut = tmp_path / 'windows.txt.gz'
    cut.write_bytes(whole.read_bytes()[: whole.stat().st_size // 2])
    with pytest.raises(ValueError, match=r"windows\.txt\.gz: not gnomAD's genomic constraint by window"):
        _convert(tmp_path, [('chr1', 2500, 'A', 'G')], {'genomic_constraint': cut})


# ---- gnomAD's regions of missense constraint ----

_MCR_HEADER = (
    'transcript\tn_total_transcript_regions\tchr\tregion_start\tregion_end\tregion_obs\tregion_exp\tregion_oe\t'
    'region_oe_chisq\tregion_oe_chisq_p\tis_high_coverage'
)
# Invented transcripts, whose ids Ensembl has not given.
_T1, _T2, _T3, _T4, _T5 = (f'ENST99{number:09d}' for number in range(1, 6))


def _region(
    transcript: str,
    regions: object,
    chrom: str,
    start: object,
    end: object,
    oe: object,
    p: object = '6.1148e-04',
    coverage: str = 'true',
) -> str:
    """A row of the invented table: a region of a transcript, 1-based and closed, with its ratio and its p-value."""
    fields = (transcript, regions, chrom, start, end, 40, '6.8322e+01', oe, '1.1741e+01', p, coverage)
    return '\t'.join(str(field) for field in fields)


# A transcript with three regions, of which the first two touch and the third follows a gap; a transcript with one
# region over parts of all three; a transcript on another chromosome; and one with a region of many stretches of the
# join. One p-value is below the smallest 32-bit float.
_REGIONS = (
    _region(_T1, 3, 'chr1', 1001, 1500, '0.0000e+00', '1.6359e-10', 'false'),
    _region(_T1, 3, 'chr1', 1501, 2000, '5.8546e-01'),
    _region(_T1, 3, 'chr1', 2301, 2600, '1.0553e-01', '7.2240e-256'),
    _region(_T2, 1, 'chr1', 1400, 2400, '3.6000e-01'),
    _region(_T3, 1, 'chr2', 500, 800, '4.8073e+00', '9.9987e-01'),
    _region(_T4, 2, 'chr2', 3000, 8999, '9.0000e-01'),
    _region(_T4, 2, 'chr2', 9000, 9000, '2.0000e-01'),
)


def _regions(path: pathlib.Path, rows: Sequence[str] = _REGIONS, header: str = _MCR_HEADER) -> pathlib.Path:
    return _table(path, header, rows)


# Alleles by position, each with the ratio of the region of each transcript that holds the position.
IN_REGIONS: dict[str, tuple[Allele, dict[str, float]]] = {
    'the base before a first region': (('chr1', 1000, 'A', 'G'), {}),
    'the first base of a region': (('chr1', 1001, 'A', 'G'), {_T1: 0.0}),
    'a base in the regions of two transcripts': (('chr1', 1450, 'A', 'C'), {_T1: 0.0, _T2: 0.36}),
    'the last base of a region': (('chr1', 1500, 'A', 'G'), {_T1: 0.0, _T2: 0.36}),
    'the first base of the region that follows it': (('chr1', 1501, 'A', 'G'), {_T1: 0.58546, _T2: 0.36}),
    'a deletion, by the position that the VCF writes it at': (('chr1', 1500, 'ATT', 'A'), {_T1: 0.0, _T2: 0.36}),
    'an insertion, by its position': (('chr1', 1501, 'A', 'ACC'), {_T1: 0.58546, _T2: 0.36}),
    'a base between two regions of a transcript': (('chr1', 2100, 'A', 'G'), {_T2: 0.36}),
    'the first base after the gap': (('chr1', 2301, 'A', 'G'), {_T1: 0.10553, _T2: 0.36}),
    'the last base of one transcript, after the end of the other': (('chr1', 2600, 'A', 'G'), {_T1: 0.10553}),
    'the base after a last region': (('chr1', 2601, 'A', 'G'), {}),
    'a base of a region on another chromosome': (('chr2', 650, 'T', 'C'), {_T3: 4.8073}),
    'that position on a chromosome without a region': (('chr3', 650, 'T', 'C'), {}),
    'the first base of a long region': (('chr2', 3000, 'A', 'G'), {_T4: 0.9}),
    'a base of it several stretches on': (('chr2', 7168, 'A', 'G'), {_T4: 0.9}),
    'the last base of a long region': (('chr2', 8999, 'A', 'G'), {_T4: 0.9}),
    'a region of one base': (('chr2', 9000, 'A', 'G'), {_T4: 0.2}),
    'the base after it': (('chr2', 9001, 'A', 'G'), {}),
}


@pytest.mark.parametrize('name', ['mcrs.tsv', 'mcrs.tsv.gz'])
def test_an_allele_has_the_region_of_each_transcript_that_holds_its_position(tmp_path: pathlib.Path, name: str) -> None:
    alleles = [allele for allele, _ in IN_REGIONS.values()]
    rows = _convert(tmp_path, alleles, {'mcr': _regions(tmp_path / name)})
    pairs = _pairs(tmp_path)
    by_allele: dict[str, dict[str, float]] = {_vid(allele): {} for allele in alleles}
    for (vid, transcript), values in pairs.items():
        by_allele[vid][transcript] = values['mcr_oe']
    expected = {
        _vid(allele): {transcript: _float32(oe) for transcript, oe in regions.items()}
        for allele, regions in IN_REGIONS.values()
    }
    assert by_allele == expected
    # The value of the allele is the lowest ratio among its regions, and none for an allele in no region.
    assert {vid: row['mcr_oe_min'] for vid, row in rows.items()} == {
        vid: min(regions.values(), default=None) for vid, regions in expected.items()
    }
    # A region gives its transcript's number of regions, its p-value in 64 bits and its coverage.
    assert pairs['chr1-1001-A-G', _T1] == {
        'mcr_oe': 0.0,
        'mcr_regions': 3,
        'mcr_oe_p': 1.6359e-10,
        'mcr_high_coverage': False,
    }
    assert pairs['chr1-2301-A-G', _T1]['mcr_oe_p'] == 7.224e-256
    assert pairs['chr1-2301-A-G', _T2] == {
        'mcr_oe': _float32(0.36),
        'mcr_regions': 1,
        'mcr_oe_p': 6.1148e-04,
        'mcr_high_coverage': True,
    }
    # The rows are in the union's order, and those of an allele in the order of the transcripts' ids.
    assert list(pairs) == [
        (_vid(allele), transcript) for allele, regions in IN_REGIONS.values() for transcript in sorted(regions)
    ]
    assert pq.read_schema(tmp_path / 'out' / joined_values.TRANSCRIPTS_NAME).equals(
        joined_values.transcripts_schema([joined_values.MCR])
    )


MCR_REFUSED: dict[str, tuple[dict[str, Any], str]] = {
    'a header with another column name': (
        {'header': _MCR_HEADER.replace('region_oe\t', 'oe\t')},
        r"its header is \['transcript', .*'oe', .*, and gnomAD's is \['transcript'",
    ),
    'a row of ten fields': (
        {'rows': [*_REGIONS, _REGIONS[0].rsplit('\t', 1)[0]]},
        "not gnomAD's regions of missense constraint",
    ),
    'a transcript with its version': (
        {'rows': [*_REGIONS, _region(f'{_T5}.2', 1, 'chr3', 100, 200, 0.5)]},
        r"rows without transcript, .* as an Ensembl transcript id without a version, .*such as \[\['ENST99000000005.2'",
    ),
    'a transcript that is not of Ensembl': (
        {'rows': [*_REGIONS, _region('NM_990000001', 1, 'chr3', 100, 200, 0.5)]},
        'rows without',
    ),
    'a number of regions of 0': ({'rows': [_region(_T5, 0, 'chr3', 100, 200, 0.5)]}, 'rows without'),
    'a chromosome without its prefix': ({'rows': [_region(_T5, 1, '3', 100, 200, 0.5)]}, 'rows without'),
    'a start that is no position': ({'rows': [_region(_T5, 1, 'chr3', '1e2', 200, 0.5)]}, 'rows without'),
    'an end before its start': ({'rows': [_region(_T5, 1, 'chr3', 200, 199, 0.5)]}, 'rows without'),
    'an end that is no position': ({'rows': [_region(_T5, 1, 'chr3', 100, 'end', 0.5)]}, 'rows without'),
    'a ratio that is no number': ({'rows': [_region(_T5, 1, 'chr3', 100, 200, 'NA')]}, 'rows without'),
    'a ratio below 0': ({'rows': [_region(_T5, 1, 'chr3', 100, 200, -0.5)]}, 'rows without'),
    'a ratio that is not finite': ({'rows': [_region(_T5, 1, 'chr3', 100, 200, 'inf')]}, 'rows without'),
    'a p-value above 1': ({'rows': [_region(_T5, 1, 'chr3', 100, 200, 0.5, p='1.5')]}, 'rows without'),
    'a p-value below 0': ({'rows': [_region(_T5, 1, 'chr3', 100, 200, 0.5, p='-0.1')]}, 'rows without'),
    'an empty p-value': ({'rows': [_region(_T5, 1, 'chr3', 100, 200, 0.5, p='')]}, 'rows without'),
    'a coverage that is neither true nor false': (
        {'rows': [_region(_T5, 1, 'chr3', 100, 200, 0.5, coverage='TRUE')]},
        'rows without',
    ),
    'two regions of a transcript that share a base': (
        {'rows': [*_REGIONS, _region(_T5, 2, 'chr3', 100, 200, 0.5), _region(_T5, 2, 'chr3', 200, 300, 0.7)]},
        rf"regions of one transcript that overlap or lie on two chromosomes, such as \[\('{_T5}', 'chr3', 200, 300\)\]",
    ),
    'a region of a transcript twice': (
        {'rows': [_region(_T5, 2, 'chr3', 100, 200, 0.5), _region(_T5, 2, 'chr3', 100, 200, 0.5)]},
        'regions of one transcript that overlap',
    ),
    'a region within another of its transcript': (
        {'rows': [_region(_T5, 2, 'chr3', 100, 900, 0.5), _region(_T5, 2, 'chr3', 300, 400, 0.7)]},
        'regions of one transcript that overlap',
    ),
    'a transcript on two chromosomes': (
        {'rows': [_region(_T5, 2, 'chr3', 100, 200, 0.5), _region(_T5, 2, 'chr4', 300, 400, 0.7)]},
        rf"or lie on two chromosomes, such as \[\('{_T5}', 'chr4', 300, 400\)\]",
    ),
    'a transcript with fewer regions than its rows say': (
        {'rows': [*_REGIONS, _region(_T5, 3, 'chr3', 100, 200, 0.5), _region(_T5, 3, 'chr3', 201, 300, 0.7)]},
        rf"another number of regions than n_total_transcript_regions says, .*such as \[\('{_T5}', 2, 3, 3\)\]",
    ),
    'a transcript whose rows give two numbers of regions': (
        {'rows': [_region(_T5, 2, 'chr3', 100, 200, 0.5), _region(_T5, 3, 'chr3', 201, 300, 0.7)]},
        rf"another number of regions than n_total_transcript_regions says, .*such as \[\('{_T5}', 2, 2, 3\)\]",
    ),
}


@pytest.mark.parametrize('case', MCR_REFUSED)
def test_a_file_that_is_not_gnomads_regions_of_missense_constraint_is_refused(
    tmp_path: pathlib.Path, case: str
) -> None:
    # The refused rows are of regions that hold no allele of the union: every row of the file is checked.
    written, message = MCR_REFUSED[case]
    file = _regions(tmp_path / 'mcrs.tsv', **written)
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, [('chr1', 1450, 'A', 'C')], {'mcr': file})
    assert not (tmp_path / 'out' / joined_values.INFO_NAME).exists()


# ---- gnomAD's MPC scores ----

_MPC_HEADER = 'locus\talleles\ttranscript\tmpc'


def _mpc_row(locus: str, alleles: str, transcript: str, mpc: object) -> str:
    """A row of the invented table: an allele as its locus and its list of two alleles, a transcript and a score."""
    return '\t'.join((locus, alleles, transcript, str(mpc)))


# An allele with a score for each of two transcripts, another change at its position, a change of two bases, an
# allele that lies in no region, and an allele on another chromosome. One score is the file's value of 6.
_MPC_ROWS = (
    _mpc_row('chr1:900', '["A","C"]', _T5, '1.2500e+00'),
    _mpc_row('chr1:1450', '["A","C"]', _T1, '2.5000e+00'),
    _mpc_row('chr1:1450', '["A","C"]', _T2, '8.8102e-01'),
    _mpc_row('chr1:1450', '["A","G"]', _T1, '6.0000e+00'),
    _mpc_row('chr1:1600', '["AC","GT"]', _T1, '5.5732e-05'),
    _mpc_row('chr2:650', '["T","C"]', _T3, '0.0000e+00'),
)


def _mpc(path: pathlib.Path, rows: Sequence[str] = _MPC_ROWS, header: str = _MPC_HEADER) -> pathlib.Path:
    return _table(path, header, rows)


@pytest.mark.parametrize('name', ['mpc.tsv.bgz', 'mpc.tsv'])
def test_an_allele_has_its_mpc_score_for_each_transcript_and_the_highest_of_them(
    tmp_path: pathlib.Path, name: str
) -> None:
    alleles = [
        ('chr1', 900, 'A', 'C'),
        ('chr1', 1450, 'A', 'C'),
        ('chr1', 1450, 'A', 'G'),
        ('chr1', 1600, 'AC', 'GT'),
        ('chr2', 650, 'T', 'C'),
        # Without a row: another change at a scored position, an insertion there, a scored change at another
        # position and on another chromosome.
        ('chr1', 1450, 'A', 'T'),
        ('chr1', 1450, 'A', 'AC'),
        ('chr1', 1451, 'A', 'C'),
        ('chr3', 1450, 'A', 'C'),
    ]
    file = _mpc(tmp_path / name)
    if name.endswith('.bgz'):
        # As gnomAD's file is: blocks of BGZF, which are the members of a gzip stream.
        pysam.tabix_compress(str(_mpc(tmp_path / 'plain.tsv')), str(file), force=True)
        assert file.read_bytes().count(b'\x1f\x8b\x08\x04') == 2
    rows = _convert(tmp_path, alleles, {'mpc': file})
    assert {vid: row['mpc_max'] for vid, row in rows.items()} == {
        'chr1-900-A-C': 1.25,
        'chr1-1450-A-C': 2.5,
        'chr1-1450-A-G': 6.0,
        'chr1-1600-AC-GT': _float32(5.5732e-05),
        # A score of 0 is a score.
        'chr2-650-T-C': 0.0,
        'chr1-1450-A-T': None,
        'chr1-1450-A-AC': None,
        'chr1-1451-A-C': None,
        'chr3-1450-A-C': None,
    }
    assert _pairs(tmp_path) == {
        ('chr1-900-A-C', _T5): {'mpc': 1.25},
        ('chr1-1450-A-C', _T1): {'mpc': 2.5},
        ('chr1-1450-A-C', _T2): {'mpc': _float32(0.88102)},
        ('chr1-1450-A-G', _T1): {'mpc': 6.0},
        ('chr1-1600-AC-GT', _T1): {'mpc': _float32(5.5732e-05)},
        ('chr2-650-T-C', _T3): {'mpc': 0.0},
    }


MPC_REFUSED: dict[str, tuple[dict[str, Any], str]] = {
    'a header with another column name': (
        {'header': _MPC_HEADER.replace('mpc', 'MPC')},
        r"its header is \['locus', 'alleles', 'transcript', 'MPC'\], and gnomAD's is \['locus'",
    ),
    'a row of three fields': (
        {'rows': [*_MPC_ROWS, _MPC_ROWS[0].rsplit('\t', 1)[0]]},
        "not gnomAD's table of MPC scores",
    ),
    'a locus without a position': (
        {'rows': [*_MPC_ROWS, _mpc_row('chr9', '["A","C"]', _T5, 1.5)]},
        r"rows without locus, alleles, transcript, mpc as .*such as \[\['chr9', '\[\"A\",\"C\"\]', 'ENST990",
    ),
    'a locus without the prefix of its chromosome': (
        {'rows': [_mpc_row('9:4990', '["A","C"]', _T5, 1.5)]},
        'rows without',
    ),
    'a locus at position 0': ({'rows': [_mpc_row('chr9:0', '["A","C"]', _T5, 1.5)]}, 'rows without'),
    'a locus with a range': ({'rows': [_mpc_row('chr9:4990-4991', '["A","C"]', _T5, 1.5)]}, 'rows without'),
    'alleles that are no list': ({'rows': [_mpc_row('chr9:4990', 'A,C', _T5, 1.5)]}, 'rows without'),
    'a list of three alleles': ({'rows': [_mpc_row('chr9:4990', '["A","C","G"]', _T5, 1.5)]}, 'rows without'),
    'an allele that is no base': ({'rows': [_mpc_row('chr9:4990', '["A","N"]', _T5, 1.5)]}, 'rows without'),
    'an empty allele': ({'rows': [_mpc_row('chr9:4990', '["A",""]', _T5, 1.5)]}, 'rows without'),
    'a transcript with its version': ({'rows': [_mpc_row('chr9:4990', '["A","C"]', f'{_T5}.3', 1.5)]}, 'rows without'),
    'an empty transcript': ({'rows': [_mpc_row('chr9:4990', '["A","C"]', '', 1.5)]}, 'rows without'),
    'a score that is no number': ({'rows': [_mpc_row('chr9:4990', '["A","C"]', _T5, 'NA')]}, 'rows without'),
    'a score that is not finite': ({'rows': [_mpc_row('chr9:4990', '["A","C"]', _T5, 'nan')]}, 'rows without'),
    'two rows of one allele of the union and one transcript': (
        {'rows': [*_MPC_ROWS, _mpc_row('chr1:1450', '["A","C"]', _T2, 1.5)]},
        rf"several rows for one allele and transcript, such as \[\('chr1', 1450, 'A', 'C', '{_T2}'\)\]",
    ),
}


@pytest.mark.parametrize('case', MPC_REFUSED)
def test_a_file_that_is_not_gnomads_mpc_scores_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    # The malformed rows are of alleles that the union lacks, so a refusal shows that every row of the file is checked.
    written, message = MPC_REFUSED[case]
    file = _mpc(tmp_path / 'mpc.tsv', **written)
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, [('chr1', 1450, 'A', 'C')], {'mpc': file})
    assert not (tmp_path / 'out' / joined_values.INFO_NAME).exists()


def test_a_table_of_mpc_scores_that_is_cut_off_is_refused(tmp_path: pathlib.Path) -> None:
    plain = _mpc(tmp_path / 'plain.tsv', [_mpc_row(f'chr1:{1000 + i}', '["A","C"]', _T1, i / 7) for i in range(5000)])
    whole = tmp_path / 'whole.tsv.bgz'
    pysam.tabix_compress(str(plain), str(whole), force=True)
    cut = tmp_path / 'mpc.tsv.bgz'
    cut.write_bytes(whole.read_bytes()[: whole.stat().st_size // 2])
    with pytest.raises(ValueError, match=r"mpc\.tsv\.bgz: not gnomAD's table of MPC scores"):
        _convert(tmp_path, [('chr1', 1450, 'A', 'C')], {'mpc': cut})


def test_an_allele_and_a_transcript_have_one_row_with_the_values_of_both_files(tmp_path: pathlib.Path) -> None:
    alleles = [('chr1', 900, 'A', 'C'), ('chr1', 1001, 'A', 'G'), ('chr1', 1450, 'A', 'C'), ('chr2', 650, 'T', 'C')]
    files = {'mcr': _regions(tmp_path / 'mcrs.tsv'), 'mpc': _mpc(tmp_path / 'mpc.tsv')}
    _convert(tmp_path, alleles, files)
    held = {pair: (row['mcr_oe'], row['mcr_regions'], row['mpc']) for pair, row in _pairs(tmp_path).items()}
    assert held == {
        # A score outside every region, a region without a score, and both.
        ('chr1-900-A-C', _T5): (None, None, 1.25),
        ('chr1-1001-A-G', _T1): (0.0, 3, None),
        ('chr1-1450-A-C', _T1): (0.0, 3, 2.5),
        ('chr1-1450-A-C', _T2): (_float32(0.36), 1, _float32(0.88102)),
        ('chr2-650-T-C', _T3): (_float32(4.8073), 1, 0.0),
    }
    assert pq.read_schema(tmp_path / 'out' / joined_values.TRANSCRIPTS_NAME).equals(
        joined_values.transcripts_schema([joined_values.MCR, joined_values.MPC])
    )


# ---- the catalogue of branch points ----

type Intron = tuple[str, str, str, int, int, int, int, str]

# The columns of a row of the paper's table by intron: gene, introns, kind and length of the intron, distance to the
# acceptor, rank, total and consensus.
_NEAR: Intron = ('INVENTED1', 'ENST01_IVS2', 'major', 900, -24, 1, 2, '1:YTNAY')
_FAR: Intron = ('INVENTED1', 'ENST02_IVS1,ENST03_IVS1', 'major', 4000, -30, 1, 1, '2:YTNA')
_MINOR: Intron = ('INVENTED2', 'ENST04_IVS7', 'minor', 85, -19, 2, 3, '0:none')


def _catalogued(
    chrom: str,
    pos: int,
    strand: str,
    kind: str = 'c',
    base: str = 'A',
    intron: Intron | None = None,
    was: int | None = None,
) -> str:
    """A row of the prepared catalogue: an invented branch point, which the lift moved from `was` if that is given."""
    sources = {'m': 'eBP_Invented,cBP_Invented', 'e': 'eBP_Invented', 'c': 'cBP_Invented'}[kind]
    name, grch37_name = (f'{kind}BP_{chrom}_{at}_{strand}_{base}' for at in (pos, pos if was is None else was))
    described = [''] * 8 if intron is None else [str(value) for value in intron]
    return '\t'.join([chrom, str(pos), strand, name, *described, sources, grch37_name])


# Invented branch points: one on each strand of chr1, far from the other; two that lie two bases apart, of which the
# first is in two introns and the second, which is no adenine, in none; and one on chr2.
_FORWARD = _catalogued('chr1', 1010, '+', 'm', intron=_NEAR, was=1510)
_REVERSE = _catalogued('chr1', 2000, '-', 'e', intron=_MINOR)
_BRANCH_POINTS = (
    _FORWARD,
    _REVERSE,
    _catalogued('chr1', 3010, '+', 'm', intron=_FAR),
    _catalogued('chr1', 3010, '+', 'm', intron=_NEAR),
    _catalogued('chr1', 3012, '+', base='C'),
    _catalogued('chr2', 500, '+'),
)
_CATALOGUE_HEADER = '\t'.join(branch_points.COLUMNS)
_PREAMBLE = (f'{branch_points.PREAMBLE}invented branch points', f'{branch_points.PREAMBLE}chain: invented.chain')


def _catalogue(
    path: pathlib.Path,
    rows: Sequence[str] = _BRANCH_POINTS,
    header: str = _CATALOGUE_HEADER,
    before: Sequence[str] = _PREAMBLE,
) -> pathlib.Path:
    """A catalogue with the lines that say what it was made from before its header, gzipped where the name says so."""
    text = ('\n'.join([*before, header, *rows]) + '\n').encode()
    path.write_bytes(gzip.compress(text) if path.suffix == '.gz' else text)
    return path


def _hits(rows: Mapping[str, Mapping[str, Any]], allele: Allele) -> list[tuple[str, list[int]]] | None:
    """The branch points of an allele's row with the positions that it changes, or None for an allele without one."""
    entries = rows[_vid(allele)]['branch_points']
    return None if entries is None else [(entry['name'], entry['hit']) for entry in entries]


# Alleles around the branch point at chr1:1010 of the forward strand, whose three bases are 1008 to 1010, and around
# the one at chr1:2000 of the reverse strand, whose three bases are 2002 down to 2000: each with the positions among
# -2, -1 and 0 that it changes, or None where it overlaps none of the three.
MOTIF_HITS: dict[str, tuple[Allele, list[int] | None]] = {
    'a change of the branch point': (('chr1', 1010, 'A', 'G'), [0]),
    'a change of the base before it': (('chr1', 1009, 'A', 'G'), [-1]),
    'a change of the base two before it': (('chr1', 1008, 'A', 'G'), [-2]),
    'a change of the base three before it': (('chr1', 1007, 'A', 'G'), None),
    'a change of the base after it': (('chr1', 1011, 'A', 'G'), None),
    'a change of two bases that ends at the base two before it': (('chr1', 1007, 'AC', 'GT'), [-2]),
    'a deletion of the three bases and more': (('chr1', 1000, 'T' * 20, 'T'), [-2, -1, 0]),
    'a deletion of the branch point and the base after it': (('chr1', 1009, 'TCA', 'T'), [0]),
    'a deletion that ends before the three bases': (('chr1', 1005, 'TCA', 'T'), None),
    'a deletion after the branch point, written from its base': (('chr1', 1010, 'AC', 'A'), None),
    'an insertion between the branch point and the base before it': (('chr1', 1009, 'A', 'ACC'), []),
    'an insertion between the two bases before it': (('chr1', 1008, 'A', 'AC'), []),
    'an insertion before the three bases': (('chr1', 1007, 'A', 'AC'), None),
    'an insertion after the branch point': (('chr1', 1010, 'A', 'AC'), None),
    'on the reverse strand, a change of the branch point': (('chr1', 2000, 'T', 'C'), [0]),
    'on the reverse strand, a change of the base before it': (('chr1', 2001, 'T', 'C'), [-1]),
    'on the reverse strand, a change of the base two before it': (('chr1', 2002, 'T', 'C'), [-2]),
    'on the reverse strand, a change of the base three before it': (('chr1', 2003, 'T', 'C'), None),
    'on the reverse strand, a change of the base after it': (('chr1', 1999, 'T', 'C'), None),
    'on the reverse strand, a deletion of the branch point and the base before it': (
        ('chr1', 1999, 'TCA', 'T'),
        [-1, 0],
    ),
    'on the reverse strand, an insertion between the branch point and the base after it': (
        ('chr1', 1999, 'T', 'TC'),
        None,
    ),
    'on the reverse strand, an insertion between the branch point and the base before it': (
        ('chr1', 2000, 'T', 'TC'),
        [],
    ),
    'a change at its position on a chromosome without one': (('chr3', 1010, 'A', 'G'), None),
    'an allele on chrM': (('chrM', 1010, 'A', 'G'), None),
}


@pytest.mark.parametrize('name', ['branch-points.tsv.gz', 'branch-points.tsv'])
def test_an_allele_has_the_branch_point_whose_three_bases_it_changes(tmp_path: pathlib.Path, name: str) -> None:
    alleles = [allele for allele, _ in MOTIF_HITS.values()]
    rows = _convert(tmp_path, alleles, {'branch_points': _catalogue(tmp_path / name)})
    found = {case: _hits(rows, allele) for case, (allele, _) in MOTIF_HITS.items()}
    assert found == {
        case: None if hit is None else [('mBP_chr1_1010_+_A' if allele[1] < 1500 else 'eBP_chr1_2000_-_A', hit)]
        for case, (allele, hit) in MOTIF_HITS.items()
    }
    # An allele without a branch point has no list, and an entry holds what the catalogue's row holds.
    schema = pq.read_schema(tmp_path / 'out' / joined_values.TABLE_NAME)
    assert schema.field('branch_points') == joined_values.BRANCH_POINTS.fields[0]
    assert schema.field('branch_points').nullable
    assert rows['chr1-1010-A-G']['branch_points'] == [
        {
            'name': 'mBP_chr1_1010_+_A',
            'pos': 1010,
            'strand': '+',
            'hit': [0],
            'acceptor_offset': -24,
            'rank': 1,
            'total': 2,
            'sources': ['eBP_Invented', 'cBP_Invented'],
            'consensus': '1:YTNAY',
            'gene': 'INVENTED1',
            'transcript_introns': ['ENST01_IVS2'],
            'intron_type': 'major',
            'intron_length': 900,
        }
    ]


def test_an_allele_has_an_entry_for_each_row_of_each_branch_point_that_it_changes(tmp_path: pathlib.Path) -> None:
    # Beside the invented ones, a branch point on each strand at one position, and one whose three bases lie over the
    # end of a stretch of 1,024 bases, which the join looks in.
    more = [_catalogued('chr2', 500, '-'), _catalogued('chr1', 2049, '+')]
    alleles: dict[Allele, list[tuple[str, list[int], int | None]] | None] = {
        # The base of a branch point in two introns, which is the base two before the next branch point: the entries
        # are in the order of position, and those of one branch point with the nearer acceptor first.
        ('chr1', 3010, 'A', 'G'): [
            ('mBP_chr1_3010_+_A', [0], -24),
            ('mBP_chr1_3010_+_A', [0], -30),
            ('cBP_chr1_3012_+_C', [-2], None),
        ],
        ('chr1', 3012, 'C', 'G'): [('cBP_chr1_3012_+_C', [0], None)],
        ('chr1', 3009, 'TAAC', 'T'): [
            ('mBP_chr1_3010_+_A', [0], -24),
            ('mBP_chr1_3010_+_A', [0], -30),
            ('cBP_chr1_3012_+_C', [-2, -1, 0], None),
        ],
        # The forward strand's entry stands before the reverse strand's.
        ('chr2', 500, 'A', 'G'): [('cBP_chr2_500_+_A', [0], None), ('cBP_chr2_500_-_A', [0], None)],
        ('chr2', 502, 'A', 'G'): [('cBP_chr2_500_-_A', [-2], None)],
        # An allele and three bases that share two stretches give one entry.
        ('chr1', 2040, 'T' * 20, 'T'): [('cBP_chr1_2049_+_A', [-2, -1, 0], None)],
        ('chr1', 2048, 'A', 'G'): [('cBP_chr1_2049_+_A', [-1], None)],
        ('chr1', 2050, 'A', 'G'): None,
    }
    file = _catalogue(tmp_path / 'branch-points.tsv', [*_BRANCH_POINTS, *more])
    rows = _convert(tmp_path, list(alleles), {'branch_points': file})
    found = {
        allele: None
        if (entries := rows[_vid(allele)]['branch_points']) is None
        else [(entry['name'], entry['hit'], entry['acceptor_offset']) for entry in entries]
        for allele in alleles
    }
    assert found == alleles
    # The entry of a branch point that the paper gives no intron has its name, its place and its sources alone.
    (alone,) = rows['chr1-3012-C-G']['branch_points']
    assert {name for name, value in alone.items() if value is None} == {
        'acceptor_offset',
        'rank',
        'total',
        'consensus',
        'gene',
        'transcript_introns',
        'intron_type',
        'intron_length',
    }
    assert alone['sources'] == ['cBP_Invented']
    # An entry gives the introns of one row as a list.
    assert rows['chr1-3010-A-G']['branch_points'][1]['transcript_introns'] == ['ENST02_IVS1', 'ENST03_IVS1']


def test_a_catalogue_is_read_with_and_without_the_lines_before_its_header(tmp_path: pathlib.Path) -> None:
    alleles = [('chr1', 1010, 'A', 'G')]
    for number, before in enumerate([(), _PREAMBLE[:1], _PREAMBLE]):
        folder = tmp_path / str(number)
        folder.mkdir()
        rows = _convert(folder, alleles, {'branch_points': _catalogue(folder / 'branch-points.tsv.gz', before=before)})
        assert _hits(rows, alleles[0]) == [('mBP_chr1_1010_+_A', [0])]


def _recatalogued(row: str, **columns: str) -> str:
    """A row of the catalogue with other text in the columns named."""
    fields = dict(zip(branch_points.COLUMNS, row.split('\t'), strict=True)) | columns
    return '\t'.join(fields.values())


_REFUSED_ROW = 'rows without chrom, pos, strand, name, gene, .*, sources, grch37_name as a primary chromosome'
CATALOGUE_REFUSED: dict[str, tuple[dict[str, Any], str]] = {
    'another header': (
        {'header': _CATALOGUE_HEADER.replace('acceptor_offset', 'dist_3ss')},
        r"its header is \['chrom', 'pos', .*'dist_3ss'.*, and the catalogue's is \['chrom', 'pos'",
    ),
    'a line before the header that does not start as the others': (
        {'before': ('# invented',)},
        r"its header is \['# invented'\], and the catalogue's is",
    ),
    'a row of thirteen fields': (
        {'rows': [_FORWARD.rsplit('\t', 1)[0]]},
        'not the prepared catalogue of branch points',
    ),
    'a chromosome that the tables lack': (
        {'rows': [_catalogued('chrUn_test', 10, '+')]},
        _REFUSED_ROW,
    ),
    'a position that is no number': (
        {'rows': [_recatalogued(_FORWARD, pos='1e3', name='mBP_chr1_1e3_+_A')]},
        _REFUSED_ROW,
    ),
    'a strand that is neither': (
        {'rows': [_recatalogued(_FORWARD, strand='.', name='mBP_chr1_1010_._A')]},
        _REFUSED_ROW,
    ),
    'a name of another position': ({'rows': [_recatalogued(_FORWARD, name='mBP_chr1_1011_+_A')]}, _REFUSED_ROW),
    'a name of the other strand': ({'rows': [_recatalogued(_FORWARD, name='mBP_chr1_1010_-_A')]}, _REFUSED_ROW),
    'a name without a kind': ({'rows': [_recatalogued(_FORWARD, name='BP_chr1_1010_+_A')]}, _REFUSED_ROW),
    'a name whose base is none': ({'rows': [_recatalogued(_FORWARD, name='mBP_chr1_1010_+_N')]}, _REFUSED_ROW),
    'no name as the paper wrote it': ({'rows': [_recatalogued(_FORWARD, grch37_name='')]}, _REFUSED_ROW),
    'a name as the paper wrote it of another kind': (
        {'rows': [_recatalogued(_FORWARD, grch37_name='cBP_chr1_1510_+_A')]},
        _REFUSED_ROW,
    ),
    'a name as the paper wrote it with another base': (
        {'rows': [_recatalogued(_FORWARD, grch37_name='mBP_chr1_1510_+_C')]},
        _REFUSED_ROW,
    ),
    'sources that are no names of sources': ({'rows': [_recatalogued(_FORWARD, sources='Invented')]}, _REFUSED_ROW),
    'no source': ({'rows': [_recatalogued(_FORWARD, sources='')]}, _REFUSED_ROW),
    'an intron without its gene': ({'rows': [_recatalogued(_FORWARD, gene='')]}, _REFUSED_ROW),
    'a gene alone': ({'rows': [_recatalogued(_catalogued('chr1', 1010, '+'), gene='INVENTED1')]}, _REFUSED_ROW),
    'an intron without its consensus': ({'rows': [_recatalogued(_FORWARD, consensus='')]}, _REFUSED_ROW),
    'introns that are no labels': ({'rows': [_recatalogued(_FORWARD, transcript_introns='ENST01')]}, _REFUSED_ROW),
    'another kind of intron': ({'rows': [_recatalogued(_FORWARD, intron_type='U12')]}, _REFUSED_ROW),
    'a length that is no number': ({'rows': [_recatalogued(_FORWARD, intron_length='long')]}, _REFUSED_ROW),
    'a distance that is not negative': ({'rows': [_recatalogued(_FORWARD, acceptor_offset='24')]}, _REFUSED_ROW),
    'a distance of 0': ({'rows': [_recatalogued(_FORWARD, acceptor_offset='-0')]}, _REFUSED_ROW),
    'a rank of 0': ({'rows': [_recatalogued(_FORWARD, rank='0')]}, _REFUSED_ROW),
    'a rank above its total': ({'rows': [_recatalogued(_FORWARD, rank='3')]}, _REFUSED_ROW),
    'a total that is no number': ({'rows': [_recatalogued(_FORWARD, total='two')]}, _REFUSED_ROW),
    'another consensus': ({'rows': [_recatalogued(_FORWARD, consensus='YTNAY')]}, _REFUSED_ROW),
    'a branch point twice in the same introns': (
        {'rows': [_FORWARD, _recatalogued(_FORWARD, acceptor_offset='-25')]},
        r"several rows for one branch point and its introns, such as \[\('mBP_chr1_1510_\+_A', \['ENST01_IVS2'\]\)\]",
    ),
    'a branch point without an intron twice': (
        {'rows': [_catalogued('chr2', 500, '+'), _catalogued('chr2', 500, '+')]},
        r"several rows for one branch point and its introns, such as \[\('cBP_chr2_500_\+_A', None\)\]",
    ),
    'no row': ({'rows': []}, 'no branch point, and the prepared catalogue holds some'),
}


@pytest.mark.parametrize('case', CATALOGUE_REFUSED)
def test_a_file_that_is_not_the_prepared_catalogue_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    written, message = CATALOGUE_REFUSED[case]
    file = _catalogue(tmp_path / 'branch-points.tsv.gz', **written)
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, [('chr1', 1010, 'A', 'G')], {'branch_points': file})


def test_a_refused_row_of_the_catalogue_is_shown_as_the_file_writes_it(tmp_path: pathlib.Path) -> None:
    row = _recatalogued(_FORWARD, rank='3')
    file = _catalogue(tmp_path / 'branch-points.tsv', [_REVERSE, row])
    with pytest.raises(ValueError, match=_REFUSED_ROW) as refused:
        _convert(tmp_path, [('chr1', 1010, 'A', 'G')], {'branch_points': file})
    assert str(row.split('\t')) in str(refused.value)
    assert 'eBP_chr1_2000' not in str(refused.value)


def test_a_catalogue_that_is_cut_off_is_refused(tmp_path: pathlib.Path) -> None:
    rows = [_catalogued('chr1', 3 * number, '+') for number in range(1, 2000)]
    whole = _catalogue(tmp_path / 'whole.tsv.gz', rows)
    cut = tmp_path / 'branch-points.tsv.gz'
    cut.write_bytes(whole.read_bytes()[: whole.stat().st_size // 2])
    with pytest.raises(ValueError, match=r'branch-points\.tsv\.gz: not the prepared catalogue of branch points'):
        _convert(tmp_path, [('chr1', 1010, 'A', 'G')], {'branch_points': cut})


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
        'mitotip': _mitotip(tmp_path / 'mitotip_scores.txt'),
        'genomic_constraint': _windows(tmp_path / 'windows.txt.gz'),
        'mcr': _regions(tmp_path / 'mcrs.tsv'),
        'mpc': _mpc(tmp_path / 'mpc.tsv'),
        'branch_points': _catalogue(tmp_path / 'branch-points.tsv.gz'),
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
    columns = ('gnomad_mito_an', 'mitotip_score', 'gencode_polya_signal', 'aparent2_min', 'mcr_oe_min')
    assert {vid: tuple(row[column] for column in columns) for vid, row in rows.items()} == {
        'chr1-1003-A-G': (None, None, True, None, 0.0),
        'chr1-4990-A-G': (None, None, False, _float32(-3.19127), None),
        'chr2-60-C-CT': (None, None, False, None, None),
        'chrM-3-G-C': (50000, _float32(12.8668), False, None, None),
    }
    by_transcript = pq.read_table(tmp_path / 'out' / joined_values.TRANSCRIPTS_NAME)
    assert by_transcript.schema.equals(joined_values.transcripts_schema(tuple(joined_values.SOURCES.values())))
    assert list(_pairs(tmp_path)) == [('chr1-1003-A-G', _T1)]
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
    'the MPC scores without the regions': ('mpc',),
    "gnomAD's windows and its regions": ('genomic_constraint', 'mcr'),
    'the branch points alone': ('branch_points',),
}


@pytest.mark.parametrize('case', SOME_FILES)
def test_a_run_without_a_file_has_none_of_its_columns(tmp_path: pathlib.Path, case: str) -> None:
    given = {name: file for name, file in _files(tmp_path).items() if name in SOME_FILES[case]}
    rows = _convert(tmp_path, _ALLELES, given)
    sources = [joined_values.SOURCES[name] for name in SOME_FILES[case]]
    schema = pq.read_schema(tmp_path / 'out' / joined_values.TABLE_NAME)
    assert schema.names == ['variant_id', *(column for source in sources for column in source.columns)]
    assert sorted(rows) == sorted(_vid(allele) for allele in _ALLELES)
    by_transcript = pq.read_schema(tmp_path / 'out' / joined_values.TRANSCRIPTS_NAME)
    assert by_transcript.names == [
        'variant_id',
        'transcript_id',
        *(column for source in sources for column in source.transcript_columns),
    ]
    info = joined_values.read_info(tmp_path / 'out')
    assert [file.source for file in info.files] == list(SOME_FILES[case])
    assert info.given() == tuple(sources)


def test_the_columns_are_32_bit_numbers_text_flags_two_lists_and_one_p_value() -> None:
    kinds = {
        field.name: field.type
        for source in joined_values.SOURCES.values()
        for field in (*source.fields, *source.transcript_fields)
    }
    assert kinds.pop('gnomad_mito_filters') == pa.list_(pa.string())
    assert kinds.pop('mcr_oe_p') == pa.float64()
    # The entries of a branch point are structs of the same kinds, with lists of text and one of 8-bit numbers.
    entry = kinds.pop('branch_points').value_type
    lists = {field.name: field.type.value_type for field in entry if isinstance(field.type, pa.ListType)}
    assert lists == {'hit': pa.int8(), 'sources': pa.string(), 'transcript_introns': pa.string()}
    assert {field.type for field in entry if field.name not in lists} == {pa.int32(), pa.string()}
    assert set(kinds.values()) == {pa.int32(), pa.float32(), pa.string(), pa.bool_()}
    # A source's name starts each of its columns, so that a column says which file it is from, and no column has the
    # name of another or of a key.
    for source in joined_values.SOURCES.values():
        assert all(column.startswith(source.name) for column in (*source.columns, *source.transcript_columns)), source
    names = [
        name for source in joined_values.SOURCES.values() for name in (*source.columns, *source.transcript_columns)
    ]
    assert len(set(names) | {'variant_id', 'transcript_id'}) == len(names) + 2


def test_a_reader_gives_values_by_transcript_where_its_source_has_such_columns(tmp_path: pathlib.Path) -> None:
    files = _files(tmp_path)
    sites = _union(tmp_path / 'union.vcf', _ALLELES)
    joined_values._write_alleles(sites, tmp_path / 'alleles.parquet')
    con = duckdb.connect()
    con.read_parquet(str(tmp_path / 'alleles.parquet')).create_view('alleles')
    for name, source in joined_values.SOURCES.items():
        values = source.read(con, files[name])
        assert (values.transcripts is not None) == bool(source.transcript_fields), name
        assert con.sql(values.alleles).columns == ['variant_id', *source.columns], name
        if values.transcripts is not None:
            assert con.sql(values.transcripts).columns == ['variant_id', 'transcript_id', *source.transcript_columns]


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
    def doubled(_con: duckdb.DuckDBPyConnection, _path: pathlib.Path) -> joined_values.Values:
        return joined_values.Values(
            'SELECT variant_id, 1 AS copies FROM alleles UNION ALL SELECT variant_id, 2 FROM alleles'
        )

    source = joined_values.Source(name='doubled', fields=(pa.field('copies', pa.int32()),), read=doubled)
    monkeypatch.setattr(joined_values, 'SOURCES', {'doubled': source})
    with pytest.raises(ValueError, match='the joined values repeats alleles'):
        _convert(tmp_path, _ALLELES, {'doubled': tmp_path / 'union.vcf'})
    assert not (tmp_path / 'out' / joined_values.INFO_NAME).exists()


def test_a_source_that_gives_a_transcript_of_an_allele_two_rows_is_refused(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def doubled(_con: duckdb.DuckDBPyConnection, _path: pathlib.Path) -> joined_values.Values:
        once = "SELECT variant_id, 'ENST99000000001' AS transcript_id, 1 AS copies FROM alleles"
        return joined_values.Values('SELECT variant_id, 1 AS copy FROM alleles', f'{once} UNION ALL {once}')

    source = joined_values.Source(
        name='doubled',
        fields=(pa.field('copy', pa.int32()),),
        read=doubled,
        transcript_fields=(pa.field('copies', pa.int32()),),
    )
    monkeypatch.setattr(joined_values, 'SOURCES', {'doubled': source})
    with pytest.raises(
        ValueError,
        match=rf"the values of doubled repeat a transcript of an allele, such as \[\('chr1-1003-A-G', '{_T1}'",
    ):
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
    by_transcript = pq.read_table(tmp_path / 'out' / joined_values.TRANSCRIPTS_NAME)
    assert (by_transcript.column_names, by_transcript.num_rows) == (['variant_id', 'transcript_id'], 0)
    assert joined_values.read_info(tmp_path / 'out').files == ()


@pytest.mark.parametrize('option', ['--file', '--normalised'])
def test_the_command_refuses_an_option_given_twice_for_one_source(tmp_path: pathlib.Path, option: str) -> None:
    files = _files(tmp_path)
    twice = [option, 'gnomad_mito', str(files['gnomad_mito']), option, 'gnomad_mito', str(files['gnomad_mito'])]
    with pytest.raises(SystemExit, match=f'{option} is given twice for gnomad_mito'):
        cli.main([*_command(tmp_path, {}), *twice])
