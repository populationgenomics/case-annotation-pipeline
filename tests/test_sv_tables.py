"""The rows of sv, from invented VCFs of a structural-variant caller and a copy-number caller."""

from __future__ import annotations

import pathlib
from collections.abc import Sequence
from typing import Any

import pytest
import structural_vcfs

from case_annotation import sv_tables

_PRIMARY = frozenset({'chr1', 'chr2', 'chr10', 'chrX', 'chrY'})
_SV_META = (*structural_vcfs.SV_INFO, *structural_vcfs.SV_FORMAT)
_CNV_META = (*structural_vcfs.CNV_INFO, *structural_vcfs.CNV_FORMAT)


def _rows(
    tmp_path: pathlib.Path,
    caller: sv_tables.Caller,
    records: Sequence[structural_vcfs.Record],
    meta: Sequence[str] | None = None,
    samples: Sequence[str] = structural_vcfs.SAMPLES,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    """The records' own columns of sv, and the members' calls of them, each by its record's id."""
    path = structural_vcfs.write(
        tmp_path / f'{caller}.vcf', meta or (_SV_META if caller == 'sv' else _CNV_META), records
    )
    rows = list(sv_tables.rows(path, caller, samples, _PRIMARY))
    return {site['sv_id']: site for site, _ in rows}, {site['sv_id']: calls for site, calls in rows}


GENOTYPE_CALLS = {
    '0/1': ((0, 1), (True, 'het')),
    '1/1': ((1, 1), (True, 'hom_alt')),
    '0/0': ((0, 0), (False, 'hom_ref')),
    './1': ((None, 1), (True, 'het')),
    # A genotype with an allele missing and no copy of the allele says neither that the member carries it nor that it
    # does not.
    '0/.': ((0, None), (None, None)),
    './.': ((None, None), (None, None)),
    # A haploid call, as on a male's X or Y.
    '1': ((1,), (True, 'hom_alt')),
    '0': ((0,), (False, 'hom_ref')),
    '.': ((None,), (None, None)),
}


@pytest.mark.parametrize('gt', GENOTYPE_CALLS)
def test_a_member_carries_a_structural_variant_by_its_genotype(gt: str) -> None:
    alleles, expected = GENOTYPE_CALLS[gt]
    assert sv_tables.genotype_call(alleles) == expected


_LOSS, _GAIN = ('<DEL>',), ('<DUP>',)
_BOTH = ('<DEL>', '<DUP>')
# The genotype, the copy number and the segment's alleles, then whether the member carries it and its zygosity.
COPY_NUMBER_CALLS = {
    'two copies where two are expected': ((None, None), 2, _LOSS, (False, 'hom_ref')),
    'one copy where two are expected': ((0, 1), 1, _LOSS, (True, 'het')),
    'no copy where two are expected': ((1, 1), 0, _LOSS, (True, 'hom_alt')),
    'three copies where two are expected': ((None, 1), 3, _GAIN, (True, None)),
    'four copies where two are expected': ((None, 1), 4, _GAIN, (True, None)),
    'a loss in a segment that another member gains': ((0, 1), 1, _BOTH, (True, 'het')),
    'a gain in a segment that another member loses': ((None, 2), 3, _BOTH, (True, None)),
    "one copy of a male's X": ((0,), 1, _LOSS, (False, 'hom_ref')),
    "no copy of a male's X": ((1,), 0, _LOSS, (True, 'hom_alt')),
    "two copies of a male's X": ((1,), 2, _GAIN, (True, None)),
    "one copy of a male's Y in a gained segment": ((0,), 1, _GAIN, (False, 'hom_ref')),
    # Without a copy number the record says neither, which a member at the expected copy number does not look like.
    "no copy number, as a female's on Y": ((None,), None, _LOSS, (None, None)),
    'no copy number of a diploid member': ((None, None), None, _GAIN, (None, None)),
}


@pytest.mark.parametrize('case', COPY_NUMBER_CALLS)
def test_a_member_carries_a_copy_number_segment_by_its_copy_number(case: str) -> None:
    alleles, cn, alts, expected = COPY_NUMBER_CALLS[case]
    assert sv_tables.copy_number_call(alleles, cn, alts, 'here') == expected


COPY_NUMBERS_REFUSED = {
    'a genotype of the reference without a copy number': ((0, 0), None, _LOSS, 'has no copy number'),
    'a negative copy number': ((0, 1), -1, _LOSS, 'a copy number of -1'),
    'an allele at the expected copy number': ((0, 1), 2, _LOSS, 'at the expected copy number, 2'),
    'no allele below the expected copy number': ((None, None), 1, _LOSS, 'names no allele'),
    'the loss named above the expected copy number': (
        (0, 1),
        3,
        _BOTH,
        r"is a <DUP>, and the genotype 0/1 names \['<DEL>'\]",
    ),
    'the gain named below the expected copy number': (
        (None, 2),
        1,
        _BOTH,
        r"is a <DEL>, and the genotype ./2 names \['<DUP>'\]",
    ),
    'one copy under a genotype of two lost': ((1, 1), 1, _LOSS, 'loses 1, and the genotype 1/1 names 2 lost'),
    'no copy under a genotype of one lost': ((0, 1), 0, _LOSS, 'loses 2, and the genotype 0/1 names 1 lost'),
    # A caller that writes a diploid member's one-copy loss as the haploid genotype 1.
    'one copy under a haploid genotype of the allele': ((1,), 1, _LOSS, 'at the expected copy number, 1'),
}


@pytest.mark.parametrize('case', COPY_NUMBERS_REFUSED)
def test_a_genotype_that_contradicts_its_copy_number_is_refused(case: str) -> None:
    alleles, cn, alts, message = COPY_NUMBERS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        sv_tables.copy_number_call(alleles, cn, alts, 'here')


BREAKENDS = {
    'G]chr10:900]': ('chr10', 900, 'left', 'left'),
    'G[chr10:900[': ('chr10', 900, 'left', 'right'),
    ']chr10:900]G': ('chr10', 900, 'right', 'left'),
    '[chr10:900[G': ('chr10', 900, 'right', 'right'),
    # Bases inserted at the junction, and a contig whose name holds colons.
    'GACT[HLA-A*01:01:01:01:55[': ('HLA-A*01:01:01:01', 55, 'left', 'right'),
}


@pytest.mark.parametrize('alt', BREAKENDS)
def test_a_breakends_allele_gives_its_mate_and_the_sides_that_are_joined(alt: str) -> None:
    chrom, pos, extends, mate_extends = BREAKENDS[alt]
    assert sv_tables.breakend(alt, 'here') == {
        'mate_chrom': chrom,
        'mate_pos': pos,
        'bnd_extends': extends,
        'mate_extends': mate_extends,
    }


@pytest.mark.parametrize('alt', ['G]chr10]', ']chr10:900]', 'G]chr10:900]T', 'G]chr10:900[', 'G]chr10:x]'])
def test_an_allele_that_is_no_breakend_with_a_mate_is_refused(alt: str) -> None:
    with pytest.raises(ValueError, match='is not a breakend with a mate'):
        sv_tables.breakend(alt, 'here')


def test_a_structural_variant_caller_gives_a_row_per_record_on_the_chromosomes_named(tmp_path: pathlib.Path) -> None:
    sites, calls = _rows(tmp_path, 'sv', structural_vcfs.SV_RECORDS)
    assert set(sites) == set(structural_vcfs.SV_IDS) == set(calls)
    deletion = sites['inv:DEL:1']
    assert deletion == {
        'chrom': 'chr1',
        'pos': 1000,
        'end_pos': 5000,
        'sv_id': 'inv:DEL:1',
        'caller': 'sv',
        'sv_type': 'DEL',
        'svlen': -4000,
        'ref': 'T',
        'alt': '<DEL>',
        'qual': 80.0,
        'filters': ['PASS'],
        'imprecise': True,
        'cipos': [-10, 10],
        'ciend': [-20, 20],
        'event_id': None,
    }
    # A record written as a sequence ends at the last base of its REF, and has both filters.
    sequence = sites['inv:DEL:2']
    assert (sequence['pos'], sequence['end_pos'], sequence['ref'], sequence['alt']) == (300, 303, 'ACGT', 'A')
    assert (sequence['filters'], sequence['imprecise'], sequence['cipos']) == (['MinQUAL', 'SampleFT'], False, None)
    # A record gives columns of sv, and never the link to a small variant, which the family's variants decide.
    assert all(set(site) <= set(sv_tables.RECORD_COLUMNS) for site in sites.values())
    assert [*sv_tables.RECORD_COLUMNS, 'small_variant_id'] == sv_tables.SV.names


def test_a_breakend_has_one_row_with_its_mate(tmp_path: pathlib.Path) -> None:
    sites, _ = _rows(tmp_path, 'sv', structural_vcfs.SV_RECORDS)
    mates = ('mate_id', 'mate_chrom', 'mate_pos', 'bnd_extends', 'mate_extends', 'event_id')
    assert [sites['inv:BND:1:0'][key] for key in mates] == ['inv:BND:1:1', 'chr10', 900, 'left', 'left', 'inv:event:1']
    assert [sites['inv:BND:1:1'][key] for key in mates] == ['inv:BND:1:0', 'chr2', 500, 'right', 'right', 'inv:event:1']
    assert (sites['inv:BND:1:0']['sv_type'], sites['inv:BND:1:0']['end_pos'], sites['inv:BND:1:0']['svlen']) == (
        'BND',
        500,
        None,
    )
    assert 'mate_chrom' not in sites['inv:DEL:1']


def test_each_members_call_of_a_structural_variant_holds_its_fields(tmp_path: pathlib.Path) -> None:
    _, calls = _rows(tmp_path, 'sv', structural_vcfs.SV_RECORDS)
    assert calls['inv:DEL:1']['S3'] == {
        'gt': '1/1',
        'carries': True,
        'zygosity': 'hom_alt',
        'ft': 'LowGQ',
        'gq': 30,
        'pl': [400, 40, 0],
        'pr_ref': 0,
        'pr_alt': 9,
        'sr_ref': 0,
        'sr_alt': 7,
        'vf_ref': 1,
        'vf_alt': 11,
        'cn': None,
        'sm': None,
        'bc': None,
        'pe_start': None,
        'pe_end': None,
        'qs': None,
        'dn': 'DeNovo',
        'dq': 3.5,
    }
    carried = {sv_id: [row[s]['carries'] for s in structural_vcfs.SAMPLES] for sv_id, row in calls.items()}
    assert carried['inv:DEL:1'] == [True, False, True]
    # A record without read counts or likelihoods, and a member without a genotype, who is not said to lack the allele.
    assert [calls['inv:DEL:2'][s]['zygosity'] for s in structural_vcfs.SAMPLES] == ['het', 'hom_ref', None]
    assert (calls['inv:DEL:2']['S1']['pl'], calls['inv:DEL:2']['S1']['pr_alt']) == (None, None)
    assert carried['inv:DEL:2'] == [True, False, None]
    assert calls['inv:DEL:2']['S3']['gt'] == './.'


def test_haploid_calls_of_a_structural_variant_on_a_males_x_and_y(tmp_path: pathlib.Path) -> None:
    _, calls = _rows(tmp_path, 'sv', structural_vcfs.SV_RECORDS)
    father = calls['inv:DEL:3']['S1']
    assert (father['gt'], father['carries'], father['zygosity'], father['pl']) == ('1', True, 'hom_alt', [80, 0])
    on_y = calls['inv:DUP:TANDEM:1']
    assert [(on_y[s]['gt'], on_y[s]['carries'], on_y[s]['zygosity']) for s in structural_vcfs.SAMPLES] == [
        ('0', False, 'hom_ref'),
        ('.', None, None),
        ('.', None, None),
    ]
    assert on_y['S1']['pl'] == [0, 70]


def test_a_copy_number_caller_gives_a_row_per_segment_away_from_the_reference(tmp_path: pathlib.Path) -> None:
    sites, calls = _rows(tmp_path, 'cnv', structural_vcfs.CNV_RECORDS)
    assert set(sites) == set(structural_vcfs.CNV_CALLS)
    assert 'inv:REF:1' not in sites
    assert {sv_id: site['sv_type'] for sv_id, site in sites.items()} == {
        'inv:LOSS:1': 'DEL',
        'inv:GAIN:1': 'DUP',
        'inv:CNV:1': 'CNV',
        'inv:LOSS:2': 'DEL',
        'inv:GAIN:2': 'DUP',
        'inv:GAIN:3': 'DUP',
        'inv:LOSS:3': 'DEL',
        'inv:GAIN:4': 'DUP',
    }
    loss, both = sites['inv:LOSS:1'], sites['inv:CNV:1']
    assert (loss['caller'], loss['pos'], loss['end_pos'], loss['svlen'], loss['alt']) == (
        'cnv',
        20_000,
        30_000,
        -10_000,
        '<DEL>',
    )
    # A segment of a loss and a gain has a length for each, and so none of its own.
    assert (both['alt'], both['svlen'], both['qual']) == ('<DEL>,<DUP>', None, None)
    for sv_id, expected in structural_vcfs.CNV_CALLS.items():
        found = tuple((calls[sv_id][s]['carries'], calls[sv_id][s]['zygosity']) for s in structural_vcfs.SAMPLES)
        assert found == expected, sv_id
    daughter = calls['inv:LOSS:1']['S3']
    values = ('gt', 'cn', 'sm', 'bc', 'pe_start', 'pe_end', 'qs', 'ft', 'pl', 'pr_ref')
    assert [daughter[key] for key in values] == ['1/1', 0, 0.75, 12, 3, 2, 40, 'PASS', None, None]
    # The genotype alone says nothing of a member at the expected copy number, or of how many copies one gained.
    assert [calls['inv:GAIN:1'][s]['gt'] for s in structural_vcfs.SAMPLES] == ['./1', './.', './1']
    assert [calls['inv:GAIN:1'][s]['cn'] for s in structural_vcfs.SAMPLES] == [3, 2, 4]
    on_y = calls['inv:LOSS:3']['S2']
    assert (on_y['gt'], on_y['cn'], on_y['sm'], on_y['carries']) == ('.', None, None, None)


def test_the_call_columns_follow_the_samples_given_whatever_the_vcfs_order(tmp_path: pathlib.Path) -> None:
    path = structural_vcfs.write(tmp_path / 'sv.vcf', _SV_META, structural_vcfs.SV_RECORDS)
    (_, calls), *_ = sv_tables.rows(path, 'sv', ('S3', 'S1', 'S2'), _PRIMARY)
    assert list(calls) == ['S3', 'S1', 'S2']
    # The members' columns follow the record's own, at the top level of the table.
    schema = sv_tables.schema(('S3', 'S1', 'S2'))
    assert schema.names == [*sv_tables.SV.names, 'S3', 'S1', 'S2']
    assert all(schema.field(s).type == sv_tables.CALL and not schema.field(s).nullable for s in ('S3', 'S1', 'S2'))


def _sv(alt: str, info: str, record_id: str = 'inv:1', ref: str = 'T') -> structural_vcfs.Record:
    return ('chr1', '1000', record_id, ref, alt, '80', 'PASS', info, 'GT', '0/1', '0/0', '0/0')


def _without(meta: Sequence[str], key: str) -> tuple[str, ...]:
    return tuple(line for line in meta if f'<ID={key},' not in line)


def _redefined(meta: Sequence[str], key: str, number: str, type_: str) -> tuple[str, ...]:
    kind = next(line for line in meta if f'<ID={key},' in line).split('=')[0].lstrip('#')
    return (*_without(meta, key), f'##{kind}=<ID={key},Number={number},Type={type_},Description="invented">')


# Per case: the caller, the records, the header lines, and the message.
RECORDS_REFUSED: dict[str, tuple[sv_tables.Caller, Sequence[structural_vcfs.Record], Sequence[str], str]] = {
    'a record without an ID': ('sv', [_sv('<DEL>', 'SVTYPE=DEL;END=2000', record_id='.')], _SV_META, 'without an ID'),
    'a structural variant of two alleles': (
        'sv',
        [_sv('<DEL>,<DUP>', 'SVTYPE=DEL;END=2000')],
        _SV_META,
        'holds one alternate allele',
    ),
    'a structural variant without an allele': (
        'sv',
        [_sv('.', 'SVTYPE=DEL;END=2000')],
        _SV_META,
        'one alternate allele',
    ),
    'a structural variant without its type': ('sv', [_sv('<DEL>', 'END=2000')], _SV_META, 'no SVTYPE'),
    'a breakend of another type': (
        'sv',
        [_sv('T]chr2:500]', 'SVTYPE=DEL')],
        _SV_META,
        'SVTYPE is DEL and the alternate',
    ),
    'a breakend type without a breakend': (
        'sv',
        [_sv('<DEL>', 'SVTYPE=BND')],
        _SV_META,
        'SVTYPE is BND and the alternate',
    ),
    'a breakend without a mate': ('sv', [_sv('T.', 'SVTYPE=INS')], _SV_META, 'not a sequence, a symbolic allele or a'),
    'a breakend with two mates': (
        'sv',
        [_sv('T]chr2:500]', 'SVTYPE=BND;MATEID=inv:2,inv:3')],
        _SV_META,
        'MATEID .* is not one name',
    ),
    'a confidence interval of three values': (
        'sv',
        [_sv('<DEL>', 'SVTYPE=DEL;END=2000;CIPOS=-1,0,1')],
        _redefined(_SV_META, 'CIPOS', '.', 'Integer'),
        'CIPOS .* is not an interval of two integers',
    ),
    'a length that is no integer': (
        'sv',
        [_sv('<DEL>', 'SVTYPE=DEL;END=2000;SVLEN=long')],
        _redefined(_SV_META, 'SVLEN', '.', 'String'),
        'SVLEN .* is not a length',
    ),
    'a copy-number segment of another allele': (
        'cnv',
        [structural_vcfs.cnv_record('chr1', 20_000, 'inv:1', '<CNV>', (('0/1', '1'),) * 3)],
        _CNV_META,
        "a copy-number segment's alleles are <DEL> and <DUP>",
    ),
    'a genotype that contradicts its copy number': (
        'cnv',
        [structural_vcfs.cnv_record('chr1', 20_000, 'inv:1', '<DEL>', (('0/1', '2'),) * 3)],
        _CNV_META,
        'chr1:20000: the genotype 0/1 names',
    ),
    'read counts for more than the two alleles': (
        'sv',
        [
            (
                'chr1',
                '1000',
                'inv:1',
                'T',
                '<DEL>',
                '80',
                'PASS',
                'SVTYPE=DEL;END=2000',
                'GT:PR',
                '0/1:5,5,1',
                '0/0:9,0',
                '0/0:9,0',
            )
        ],
        _SV_META,
        'PR .* is not a count for the reference and one for the allele',
    ),
    "likelihoods of a diploid call that are a haploid one's": (
        'sv',
        [
            (
                'chr1',
                '1000',
                'inv:1',
                'T',
                '<DEL>',
                '80',
                'PASS',
                'SVTYPE=DEL;END=2000',
                'GT:PL',
                '0/1:5,0',
                '0/0:0,9,90',
                '0/0:0,9,90',
            )
        ],
        _SV_META,
        'a call of 2 alleles has 3',
    ),
    'a copy number that is no integer': (
        'cnv',
        [structural_vcfs.cnv_record('chr1', 20_000, 'inv:1', '<DEL>', (('0/1', '1.5'),) * 3)],
        _redefined(_CNV_META, 'CN', '1', 'Float'),
        'CN is 1.5, not an integer',
    ),
}


@pytest.mark.parametrize('case', RECORDS_REFUSED)
def test_a_record_outside_the_contract_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    caller, records, meta, message = RECORDS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _rows(tmp_path, caller, records, meta)


def test_a_record_off_the_chromosomes_named_is_left_out_unread(tmp_path: pathlib.Path) -> None:
    off = ('chr1_KI270706v1_random', '50', '.', 'T', '<DEL>,<DUP>', '80', 'PASS', '.', 'GT', '0/1', '0/0', '0/0')
    assert _rows(tmp_path, 'sv', [off]) == ({}, {})


HEADERS_REFUSED: dict[str, tuple[sv_tables.Caller, Sequence[str], Sequence[str], str]] = {
    'a sample the family lacks': ('sv', _SV_META, ('S1', 'S2', 'S9'), r"its samples are \['S1', 'S2', 'S9'\]"),
    'a member without a column': ('sv', _SV_META, ('S1', 'S2'), r"the family's VCF has \['S1', 'S2', 'S3'\]"),
    'a member more than the family has': ('cnv', _CNV_META, ('S1', 'S2', 'S3', 'S4'), 'its samples are'),
    'no genotypes': ('sv', _without(_SV_META, 'GT'), structural_vcfs.SAMPLES, 'no GT among its FORMAT fields'),
    'no copy numbers': ('cnv', _without(_CNV_META, 'CN'), structural_vcfs.SAMPLES, 'no CN among its FORMAT fields'),
}


@pytest.mark.parametrize('case', HEADERS_REFUSED)
def test_a_vcf_that_does_not_fit_the_family_is_refused_by_its_header(tmp_path: pathlib.Path, case: str) -> None:
    caller, meta, samples, message = HEADERS_REFUSED[case]
    path = structural_vcfs.write(tmp_path / 'callers.vcf', meta, [], samples)
    with pytest.raises(ValueError, match=message):
        sv_tables.check(path, caller, structural_vcfs.SAMPLES)
    with pytest.raises(ValueError, match=message):
        list(sv_tables.rows(path, caller, structural_vcfs.SAMPLES, _PRIMARY))


def test_a_copy_number_vcf_is_no_structural_variant_callers_and_the_other_way_round(tmp_path: pathlib.Path) -> None:
    # Each caller's file given as the other's fails at its first record rather than filling the wrong columns.
    with pytest.raises(ValueError, match='holds one alternate allele'):
        _rows(tmp_path, 'sv', structural_vcfs.CNV_RECORDS, _CNV_META)
    path = structural_vcfs.write(tmp_path / 'swapped.vcf', _SV_META, structural_vcfs.SV_RECORDS)
    with pytest.raises(ValueError, match='no CN among its FORMAT fields'):
        sv_tables.check(path, 'cnv', structural_vcfs.SAMPLES)
