"""The rows of repeats, from invented repeat VCFs of the members of a trio."""

from __future__ import annotations

import pathlib
from collections.abc import Sequence
from typing import Any

import pytest
import structural_vcfs

from case_annotation import repeat_tables

_PRIMARY = frozenset({'chr1', 'chr2', 'chr10', 'chrX', 'chrY'})
_META = (*structural_vcfs.REPEAT_INFO, *structural_vcfs.REPEAT_FORMAT)


def _rows(tmp_path: pathlib.Path) -> dict[str, dict[str, Any]]:
    paths = structural_vcfs.write_repeats(tmp_path)
    return {row['repeat_id']: row for row in repeat_tables.rows(paths, structural_vcfs.SAMPLES, _PRIMARY)}


def test_a_repeat_has_one_row_with_a_column_per_member(tmp_path: pathlib.Path) -> None:
    rows = _rows(tmp_path)
    assert set(rows) == set(structural_vcfs.REPEAT_IDS)
    assert list(rows['RPTA']) == [*repeat_tables.LOCUS.names, *structural_vcfs.SAMPLES]
    assert {key: rows['RPTB_GCC'][key] for key in repeat_tables.LOCUS.names} == {
        'chrom': 'chr2',
        'pos': 700,
        'end_pos': 721,
        'repeat_id': 'RPTB_GCC',
        'repeat_unit': 'GCC',
        'ref_copies': 7,
    }


def test_a_members_column_holds_a_value_for_each_of_its_alleles(tmp_path: pathlib.Path) -> None:
    rows = _rows(tmp_path)
    assert rows['RPTA']['S3'] == {
        'gt': '1/2',
        'repcn': [12, 45],
        'repci_lower': [12, 40],
        'repci_upper': [12, 58],
        'so': ['SPANNING', 'INREPEAT'],
        'adsp': [5, 0],
        'adfl': [7, 9],
        'adir': [0, 4],
        'lc': 28.0,
        'filters': ['PASS'],
    }
    # A haploid genotype, as on a male's X, has one value of each.
    father = rows['RPTX']['S1']
    assert (father['gt'], father['repcn'], father['repci_upper'], father['so']) == ('1', [23], [23], ['SPANNING'])
    assert rows['RPTX']['S2']['repcn'] == [23, 30]


def test_a_member_without_a_call_of_a_repeat_has_no_values(tmp_path: pathlib.Path) -> None:
    mother = _rows(tmp_path)['RPTB_GCC']['S2']
    assert mother == {
        'gt': './.',
        'repcn': None,
        'repci_lower': None,
        'repci_upper': None,
        'so': None,
        'adsp': None,
        'adfl': None,
        'adir': None,
        'lc': 2.5,
        'filters': ['PASS'],
    }


def test_a_repeat_off_the_chromosomes_named_is_left_out(tmp_path: pathlib.Path) -> None:
    paths = structural_vcfs.write_repeats(tmp_path)
    rows = repeat_tables.rows(paths, structural_vcfs.SAMPLES, frozenset({'chr1', 'chrX'}))
    assert [row['repeat_id'] for row in rows] == ['RPTA', 'RPTX']


def _family(tmp_path: pathlib.Path, changed: Sequence[structural_vcfs.Record], meta: Sequence[str] = _META) -> None:
    """Reads the trio's repeat VCFs with the daughter's records replaced by those given."""
    paths = structural_vcfs.write_repeats(tmp_path, ('S1', 'S2'))
    daughter = structural_vcfs.write(tmp_path / 'S3.repeats.vcf', meta, changed, ('S3',))
    repeat_tables.rows((*paths, daughter), structural_vcfs.SAMPLES, _PRIMARY)


def _with(repeat_id: str, values: str, info: str | None = None) -> list[structural_vcfs.Record]:
    """The daughter's records, the one of a repeat with other values or another INFO."""
    records = []
    for record in structural_vcfs.repeat_records('S3'):
        if f'VARID={repeat_id}' in record[7]:
            records.append((*record[:7], info if info is not None else record[7], record[8], values))
        else:
            records.append(record)
    return records


def _redefined(key: str, number: str, type_: str) -> tuple[str, ...]:
    kept = tuple(line for line in _META if f'<ID={key},' not in line)
    kind = next(line for line in _META if f'<ID={key},' in line).split('=')[0].lstrip('#')
    return (*kept, f'##{kind}=<ID={key},Number={number},Type={type_},Description="invented">')


_RPTA = 'END=2030;REF=10;REPID=RPTA;RU=CAG;VARID=RPTA'
RECORDS_REFUSED: dict[str, tuple[list[structural_vcfs.Record], Sequence[str], str]] = {
    'one size for two alleles': (
        _with('RPTA', '1/2:SPANNING/INREPEAT:12:12-12/40-58:5/0:7/9:0/4:28.0'),
        _META,
        "REPCN '12' does not hold a value for each of the 2 alleles",
    ),
    'two sizes for a haploid genotype': (
        _with('RPTX', '1:SPANNING:23/23:23-23:7:9:0:32.0'),
        _META,
        "REPCN '23/23' does not hold a value for each of the 1 alleles",
    ),
    'a size missing beside another': (
        _with('RPTA', '1/2:SPANNING/INREPEAT:12/.:12-12/40-58:5/0:7/9:0/4:28.0'),
        _META,
        "REPCN '12/.' does not hold a value",
    ),
    'a size that is no integer': (
        _with('RPTA', '1/2:SPANNING/INREPEAT:12/many:12-12/40-58:5/0:7/9:0/4:28.0'),
        _META,
        "REPCN '12/many' does not hold integers",
    ),
    'a read count that is no integer': (
        _with('RPTA', '1/2:SPANNING/INREPEAT:12/45:12-12/40-58:5/0.5:7/9:0/4:28.0'),
        _META,
        "ADSP '5/0.5' does not hold integers",
    ),
    'a confidence interval without its upper bound': (
        _with('RPTA', '1/2:SPANNING/INREPEAT:12/45:12/40-58:5/0:7/9:0/4:28.0'),
        _META,
        "REPCI '12/40-58' does not hold an interval",
    ),
    'a repeat without its id': (
        _with('RPTA', structural_vcfs.REPEAT_CALLS['S3']['RPTA'][1], info='END=2030;REF=10;REPID=RPTA;RU=CAG'),
        _META,
        'no VARID',
    ),
    # With REF declared a string, every record's reference size is one.
    'a reference size that is no integer': (
        structural_vcfs.repeat_records('S3'),
        _redefined('REF', '1', 'String'),
        "REF '7' is not the reference's number of repeat units",
    ),
    'a repeat with two records': (
        [*structural_vcfs.repeat_records('S3'), structural_vcfs.repeat_records('S3')[0]],
        _META,
        "the repeat 'RPTB_GCC' has two records",
    ),
    'a repeat the other members lack': (
        _with('RPTA', structural_vcfs.REPEAT_CALLS['S3']['RPTA'][1], info=_RPTA.replace('VARID=RPTA', 'VARID=RPTC')),
        _META,
        r"not of one catalogue: 2 differ, such as \['RPTA', 'RPTC'\]",
    ),
    'a repeat fewer than the other members have': (
        structural_vcfs.repeat_records('S3')[:2],
        _META,
        r"not of one catalogue: 1 differ, such as \['RPTX'\]",
    ),
    'a repeat described otherwise than in the other members': (
        _with('RPTA', structural_vcfs.REPEAT_CALLS['S3']['RPTA'][1], info=_RPTA.replace('RU=CAG', 'RU=CTG')),
        _META,
        r"not of one catalogue: 1 differ, such as \['RPTA'\]",
    ),
}


@pytest.mark.parametrize('case', RECORDS_REFUSED)
def test_a_repeat_record_outside_the_contract_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    records, meta, message = RECORDS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _family(tmp_path, records, meta)


def test_the_first_members_empty_file_is_not_taken_for_any_catalogue(tmp_path: pathlib.Path) -> None:
    empty = structural_vcfs.write(tmp_path / 'S1.repeats.vcf', _META, [], ('S1',))
    others = structural_vcfs.write_repeats(tmp_path, ('S2', 'S3'))
    with pytest.raises(ValueError, match='not of one catalogue: 3 differ'):
        repeat_tables.rows((empty, *others), structural_vcfs.SAMPLES, _PRIMARY)


def _files(tmp_path: pathlib.Path, samples: Sequence[str]) -> tuple[pathlib.Path, ...]:
    """A repeat VCF per sample named, in the order given; a sample outside the trio has the daughter's records."""
    return tuple(
        structural_vcfs.write(
            tmp_path / f'file-{index}.vcf',
            _META,
            structural_vcfs.repeat_records(sample if sample in structural_vcfs.SAMPLES else 'S3'),
            (sample,),
        )
        for index, sample in enumerate(samples)
    )


FILES_REFUSED = {
    'some but not all members': (('S1', 'S3'), r"no repeat VCF of the samples \['S2'\]"),
    'one member alone': (('S2',), r"no repeat VCF of the samples \['S1', 'S3'\]"),
    'a member twice': (('S1', 'S2', 'S3', 'S1'), "holds the sample 'S1', as .*file-0.vcf does"),
    "another family's member": (('S1', 'S2', 'S3', 'S9'), r"hold the samples \['S9'\], which the family's VCF lacks"),
}


@pytest.mark.parametrize('case', FILES_REFUSED)
def test_repeat_vcfs_that_are_not_one_of_each_member_are_refused(tmp_path: pathlib.Path, case: str) -> None:
    samples, message = FILES_REFUSED[case]
    paths = _files(tmp_path, samples)
    with pytest.raises(ValueError, match=message):
        repeat_tables.by_sample(paths, structural_vcfs.SAMPLES)
    with pytest.raises(ValueError, match=message):
        repeat_tables.rows(paths, structural_vcfs.SAMPLES, _PRIMARY)


HEADERS_REFUSED = {
    'two samples in one file': (_META, ('S1', 'S2'), r"holds one member's genotypes, and this one has the samples"),
    'no genotypes': (tuple(line for line in _META if '<ID=GT,' not in line), ('S1',), 'no GT or REPCN'),
    'no repeat sizes': (tuple(line for line in _META if '<ID=REPCN,' not in line), ('S1',), 'no GT or REPCN'),
}


@pytest.mark.parametrize('case', HEADERS_REFUSED)
def test_a_file_that_is_no_members_repeat_vcf_is_refused_by_its_header(tmp_path: pathlib.Path, case: str) -> None:
    meta, samples, message = HEADERS_REFUSED[case]
    path = structural_vcfs.write(tmp_path / 'S1.repeats.vcf', meta, [], samples)
    with pytest.raises(ValueError, match=message):
        repeat_tables.by_sample((path,), ('S1',))
