"""STRchive's loci from an invented file, and the entry that each invented repeat gets."""

from __future__ import annotations

import datetime as dt
import pathlib
from collections.abc import Mapping, Sequence
from typing import Any

import invented_strchive
import pytest
import structural_vcfs

from case_annotation import repeat_tables, strchive

_RANGES = ('benign', 'intermediate', 'pathogenic')
_BOUNDS = tuple(f'{kind}_{bound}' for kind in _RANGES for bound in ('min', 'max'))
# The fields of an entry that the tables read.
_READ = (
    'id',
    'gene',
    'chrom',
    'start_hg38',
    'stop_hg38',
    'inheritance',
    'reference_motif_reference_orientation',
    'pathogenic_motif_reference_orientation',
    *_BOUNDS,
)


def _read(tmp_path: pathlib.Path, *loci: Mapping[str, Any]) -> tuple[strchive.Locus, ...]:
    return strchive.read(invented_strchive.write(tmp_path / 'loci.json', list(loci)))


def test_a_locus_is_read_with_its_place_its_motifs_and_its_values(tmp_path: pathlib.Path) -> None:
    first, second, third = strchive.read(invented_strchive.write(tmp_path / 'loci.json'))
    assert (first.id, first.chrom, first.start, first.stop, first.motifs) == (
        'INV1_GENEA',
        'chr1',
        1_990,
        2_040,
        ('AGC',),
    )
    assert first.columns == {
        'strchive_id': 'INV1_GENEA',
        'strchive_gene': 'GENEA',
        'strchive_inheritance': ['AD'],
        'strchive_start_hg38': 1_990,
        'strchive_stop_hg38': 2_040,
        'strchive_pathogenic_motif_reference_orientation': ['AGC'],
        'strchive_benign_min': 5,
        'strchive_benign_max': 20,
        'strchive_intermediate_min': 21,
        'strchive_intermediate_max': 35,
        'strchive_pathogenic_min': 36,
        'strchive_pathogenic_max': 200,
    }
    assert list(first.columns) == list(strchive.COLUMNS)
    # The reference motif and the pathogenic one, each once, and a range that the file leaves empty as nulls.
    assert second.motifs == ('GCN', 'NGC')
    assert (second.columns['strchive_intermediate_min'], second.columns['strchive_intermediate_max']) == (None, None)
    assert second.columns['strchive_inheritance'] == ['AD', 'AR']
    # A list without an entry is a null, as every list of the tables is.
    assert third.columns['strchive_inheritance'] is None


def test_a_bound_written_with_a_decimal_point_is_the_whole_number(tmp_path: pathlib.Path) -> None:
    # The releases before v2.12.0 write their bounds so.
    (locus,) = _read(tmp_path, invented_strchive.locus(benign_min=5.0, pathogenic_max=200.0))
    assert (locus.columns['strchive_benign_min'], locus.columns['strchive_pathogenic_max']) == (5, 200)
    assert all(type(locus.columns[f'strchive_{bound}']) is int for bound in _BOUNDS)


def test_the_motifs_compared_are_as_long_as_the_first_pathogenic_motif(tmp_path: pathlib.Path) -> None:
    # The ranges count units of the pathogenic motif. A reference motif of another length is of another unit, and the
    # entry's own motif_len is not read.
    (locus,) = _read(
        tmp_path,
        invented_strchive.locus(
            reference_motif_reference_orientation=['TTTC'],
            pathogenic_motif_reference_orientation=['CT', 'TC'],
            motif_len=4,
        ),
    )
    assert locus.motifs == ('CT', 'TC')
    assert locus.columns['strchive_pathogenic_motif_reference_orientation'] == ['CT', 'TC']


def test_a_range_of_one_size_is_read(tmp_path: pathlib.Path) -> None:
    (locus,) = _read(tmp_path, invented_strchive.locus(benign_min=4, benign_max=4, pathogenic_min=3, pathogenic_max=3))
    assert locus.columns['strchive_benign_min'] == locus.columns['strchive_benign_max'] == 4


def _without(key: str) -> list[dict[str, Any]]:
    return [{name: value for name, value in invented_strchive.locus().items() if name != key}]


FILES_REFUSED: Mapping[str, tuple[str | bytes, str]] = {
    'a table and no JSON': ('chrom\tstart\tstop\nchr1\t1990\t2040\n', 'not JSON'),
    'JSON that is cut off': ('[{"id": "INV1_GENEA", "gene": ', 'not JSON'),
    'bytes that are no text': (b'\xff\xfe\x00[', 'not JSON'),
    'an object and no list': ('{"loci": []}', "not STRchive's loci"),
    'a list without a locus': ('[]', "not STRchive's loci"),
    'a list that holds something other than objects': ('[{"id": "INV1_GENEA"}, "INV2_GENEB"]', "not STRchive's loci"),
}


@pytest.mark.parametrize('case', FILES_REFUSED)
def test_a_file_that_is_not_a_list_of_loci_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    content, message = FILES_REFUSED[case]
    path = tmp_path / 'loci.json'
    path.write_bytes(content.encode() if isinstance(content, str) else content)
    with pytest.raises(ValueError, match=f'loci.json: {message}'):
        strchive.read(path)


@pytest.mark.parametrize('key', _READ)
def test_a_locus_without_a_field_that_is_read_is_refused(tmp_path: pathlib.Path, key: str) -> None:
    with pytest.raises(ValueError, match=f'locus 1.*: no {key}, which every locus of STRchive has'):
        _read(tmp_path, *_without(key))


def _changed(**changed: Any) -> list[dict[str, Any]]:  # noqa: ANN401 (JSON values)
    return [invented_strchive.locus(**changed)]


_NOT_WHOLE = 'not a whole number from 0 to 2147483647'
LOCI_REFUSED: Mapping[str, tuple[Sequence[Mapping[str, Any]], str]] = {
    'an id that is no text': (_changed(id=5), 'locus 1: id is 5, not a text'),
    'an empty id': (_changed(id=''), "locus 1: id is '', not a text"),
    'a gene that is no text': (_changed(gene=None), r'locus 1 \(INV1_GENEA\): gene is None, not a text'),
    'a chromosome without the chr prefix': (_changed(chrom='1'), "chrom is '1', not a primary chromosome"),
    'a contig that is no primary chromosome': (_changed(chrom='chrUn_KI270302v1'), 'not a primary chromosome'),
    'a start with a fraction': (_changed(start_hg38=1_990.5), f'start_hg38 is 1990.5, {_NOT_WHOLE}'),
    'a start written as a text': (_changed(start_hg38='1990'), f"start_hg38 is '1990', {_NOT_WHOLE}"),
    'a stop that is true': (_changed(stop_hg38=True), f'stop_hg38 is True, {_NOT_WHOLE}'),
    'a stop that is missing a value': (_changed(stop_hg38=None), f'stop_hg38 is None, {_NOT_WHOLE}'),
    'a start below 0': (_changed(start_hg38=-1), f'start_hg38 is -1, {_NOT_WHOLE}'),
    'a stop at its start': (_changed(stop_hg38=1_990), 'start_hg38 1990 is not before stop_hg38 1990'),
    'a stop before its start': (_changed(stop_hg38=1_900), 'start_hg38 1990 is not before stop_hg38 1900'),
    'a bound with a fraction': (_changed(benign_max=20.5), f'benign_max is 20.5, {_NOT_WHOLE}'),
    'a bound written as a text': (_changed(pathogenic_min='36'), f"pathogenic_min is '36', {_NOT_WHOLE}"),
    'a bound below 0': (_changed(intermediate_min=-1), f'intermediate_min is -1, {_NOT_WHOLE}'),
    'a bound too large for its column': (_changed(pathogenic_max=2**31), f'pathogenic_max is 2147483648, {_NOT_WHOLE}'),
    **{
        f'a {kind} minimum above its maximum': (
            _changed(**{f'{kind}_min': 201, f'{kind}_max': 200}),
            f'{kind}_min 201 is above {kind}_max 200',
        )
        for kind in _RANGES
    },
    'motifs that are one text and no list': (
        _changed(pathogenic_motif_reference_orientation='AGC'),
        "pathogenic_motif_reference_orientation is 'AGC', not a list of texts",
    ),
    'a motif that is no text': (
        _changed(reference_motif_reference_orientation=['AGC', None]),
        'reference_motif_reference_orientation is .*, not a list of texts',
    ),
    'an empty motif': (_changed(pathogenic_motif_reference_orientation=['']), 'not a list of texts'),
    'a motif in lower case': (
        _changed(pathogenic_motif_reference_orientation=['agc']),
        r"pathogenic_motif_reference_orientation holds \['agc'\], which are no motifs",
    ),
    'a motif with a character that is no code': (
        _changed(reference_motif_reference_orientation=['AG-C']),
        r"reference_motif_reference_orientation holds \['AG-C'\], which are no motifs",
    ),
    'no pathogenic motif': (_changed(pathogenic_motif_reference_orientation=[]), 'no pathogenic motif'),
    'an inheritance that is one text and no list': (_changed(inheritance='AD'), "inheritance is 'AD', not a list"),
    'an inheritance that holds no text': (_changed(inheritance=[None]), 'inheritance is .*, not a list of texts'),
    'two loci of one id': (
        [invented_strchive.locus(), invented_strchive.locus(chrom='chr2')],
        r'locus 2 \(INV1_GENEA\): a second locus of this id',
    ),
}


@pytest.mark.parametrize('case', LOCI_REFUSED)
def test_a_locus_outside_what_strchive_publishes_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    loci, message = LOCI_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _read(tmp_path, *loci)


SAME = {
    'the same unit': ('CAG', 'CAG', True),
    'a rotation of the unit': ('CAG', 'AGC', True),
    'the other rotation': ('GCA', 'CAG', True),
    'a rotation of a longer unit': ('AAAAT', 'TAAAA', True),
    'another unit of the same length': ('GCC', 'CAG', False),
    'the unit of the other strand': ('CAG', 'CTG', False),
    'a unit that the code N of a motif can be': ('GCG', 'GCN', True),
    'a unit with a code, rotated': ('NGC', 'GCN', True),
    'a unit whose codes can be the motif': ('AARRG', 'AAAAG', True),
    'a unit whose codes can be another motif': ('AARRG', 'AAGGG', True),
    'a unit whose codes cannot be the motif': ('AARRG', 'ACAGG', False),
    'a unit that is the motif twice': ('CAGCAG', 'CAG', False),
    'a unit of another length': ('TTTC', 'CT', False),
    'a unit with a character that is no code': ('CAX', 'CAG', False),
    'a unit in lower case': ('cag', 'CAG', False),
}


@pytest.mark.parametrize('case', SAME)
def test_a_unit_is_of_a_motif_in_any_rotation_and_by_what_its_codes_stand_for(case: str) -> None:
    unit, motif, same = SAME[case]
    assert strchive.same_repeat(unit, motif) is same


def _repeat(**changed: Any) -> dict[str, Any]:  # noqa: ANN401 (a row's values)
    """A row of repeats as far as the match reads it: the repeat on chr1 of structural_vcfs."""
    return {'chrom': 'chr1', 'pos': 2_000, 'end_pos': 2_030, 'repeat_id': 'RPTA', 'repeat_unit': 'CAG'} | changed


def _matched(tmp_path: pathlib.Path, row: Mapping[str, Any], *loci: Mapping[str, Any]) -> str | None:
    (annotated,), matched = strchive.annotate([row], _read(tmp_path, *loci))
    assert matched == (annotated['strchive_id'] is not None)
    return annotated['strchive_id']


def test_each_repeat_gets_the_values_of_its_locus_and_one_without_a_locus_nulls(tmp_path: pathlib.Path) -> None:
    paths = structural_vcfs.write_repeats(tmp_path)
    rows = repeat_tables.rows(paths, structural_vcfs.SAMPLES, frozenset({'chr1', 'chr2', 'chrX'}))
    loci = strchive.read(invented_strchive.write(tmp_path / 'loci.json'))
    annotated, matched = strchive.annotate(rows, loci)
    assert [row['repeat_id'] for row in annotated] == [row['repeat_id'] for row in rows]
    by_id = {row['repeat_id']: row for row in annotated}
    assert {repeat: row['strchive_id'] for repeat, row in by_id.items() if row['strchive_id'] is not None} == (
        invented_strchive.MATCHED
    )
    assert matched == len(invented_strchive.MATCHED)
    by_locus = {locus.id: locus for locus in loci}
    for row, before in zip(annotated, rows, strict=True):
        # A row keeps what it held, and gains each column of STRchive.
        assert {name: row[name] for name in before} == before
        assert row.keys() - before.keys() == set(strchive.COLUMNS)
        expected = by_locus[row['strchive_id']].columns if row['strchive_id'] else dict.fromkeys(strchive.COLUMNS)
        assert {name: row[name] for name in strchive.COLUMNS} == expected
    assert {repeat for repeat, row in by_id.items() if row['strchive_id'] is None} == set(invented_strchive.UNMATCHED)


# A repeat's interval is the bases after `pos` up to `end_pos`, and a locus's those from its start to before its stop.
PLACES: Mapping[str, tuple[Mapping[str, Any], Mapping[str, Any], bool]] = {
    'the same interval': (_repeat(), invented_strchive.locus(start_hg38=2_000, stop_hg38=2_030), True),
    'a locus within the repeat': (_repeat(), invented_strchive.locus(start_hg38=2_010, stop_hg38=2_013), True),
    'a locus that holds the last base of the repeat': (
        _repeat(),
        invented_strchive.locus(start_hg38=2_029, stop_hg38=2_090),
        True,
    ),
    'a locus that begins after the last base of the repeat': (
        _repeat(),
        invented_strchive.locus(start_hg38=2_030, stop_hg38=2_090),
        False,
    ),
    'a locus that holds the first base of the repeat': (
        _repeat(),
        invented_strchive.locus(start_hg38=1_900, stop_hg38=2_001),
        True,
    ),
    'a locus that ends before the first base of the repeat': (
        _repeat(),
        invented_strchive.locus(start_hg38=1_900, stop_hg38=2_000),
        False,
    ),
    'a locus at the same place of another chromosome': (_repeat(chrom='chr2'), invented_strchive.locus(), False),
    'a repeat of another unit in the locus': (_repeat(repeat_unit='GCC'), invented_strchive.locus(), False),
    'a repeat whose file gives no unit': (_repeat(repeat_unit=None), invented_strchive.locus(), False),
    'a repeat of the pathogenic motif where the reference has another': (
        _repeat(repeat_unit='AAATG'),
        invented_strchive.locus(
            reference_motif_reference_orientation=['AAAAT'], pathogenic_motif_reference_orientation=['TGAAA']
        ),
        True,
    ),
    'a repeat of a reference motif of another length than the pathogenic one': (
        _repeat(repeat_unit='TTTC'),
        invented_strchive.locus(
            reference_motif_reference_orientation=['TTTC'], pathogenic_motif_reference_orientation=['CT']
        ),
        False,
    ),
}


@pytest.mark.parametrize('case', PLACES)
def test_a_repeat_is_of_the_locus_that_shares_a_base_and_a_unit_with_it(tmp_path: pathlib.Path, case: str) -> None:
    row, locus, same = PLACES[case]
    assert _matched(tmp_path, row, locus) == ('INV1_GENEA' if same else None)


def test_a_repeat_that_matches_two_loci_is_refused(tmp_path: pathlib.Path) -> None:
    loci = _read(tmp_path, invented_strchive.locus(), invented_strchive.locus(id='INV9_GENEA', start_hg38=2_020))
    with pytest.raises(ValueError, match=r"'RPTA' at chr1:2000-2030 matches the STRchive loci \['INV1_GENEA', 'INV9"):
        strchive.annotate([_repeat()], loci)
    # Each of two repeats can be of one locus, as a repeat and the repeat beside it in a catalogue are.
    rows = [_repeat(), _repeat(repeat_id='RPTA_2', pos=2_030, end_pos=2_036)]
    annotated, matched = strchive.annotate(rows, loci[:1])
    assert ([row['strchive_id'] for row in annotated], matched) == (['INV1_GENEA', 'INV1_GENEA'], 2)


def test_a_release_is_taken_with_its_name_and_its_date_as_stated(tmp_path: pathlib.Path) -> None:
    stated = strchive.stated(tmp_path / 'loci.json', 'v2.28.0', '2026-10-08')
    assert stated == strchive.Release(file=tmp_path / 'loci.json', name='v2.28.0', date=dt.date(2026, 10, 8))


STATED_REFUSED: Mapping[str, tuple[str, str, str]] = {
    'a name with a space': ('v 2.28.0', '2026-10-08', "'v 2.28.0' is no release's name"),
    'an empty name': ('', '2026-10-08', "'' is no release's name"),
    'a name that a shell would read': ('v2;rm', '2026-10-08', "is no release's name"),
    'a date in words': ('v2.28.0', '8 October 2026', "is '8 October 2026', not a date as YYYY-MM-DD"),
    'a date without its separators': ('v2.28.0', '20261008', "is '20261008', not a date as YYYY-MM-DD"),
    'a date with a day of one digit': ('v2.28.0', '2026-10-8', 'not a date as YYYY-MM-DD'),
    'a date with a time': ('v2.28.0', '2026-10-08T16:37:54Z', 'not a date as YYYY-MM-DD'),
    'a thirteenth month': ('v2.28.0', '2026-13-01', "is '2026-13-01', not a day of the calendar"),
    'a day that the month lacks': ('v2.28.0', '2026-02-30', 'not a day of the calendar'),
}


@pytest.mark.parametrize('case', STATED_REFUSED)
def test_a_release_that_is_not_stated_as_a_name_and_a_date_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    name, date, message = STATED_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        strchive.stated(tmp_path / 'loci.json', name, date)


def test_each_column_holds_the_field_of_its_name_under_the_prefix(tmp_path: pathlib.Path) -> None:
    entry = invented_strchive.locus()
    (locus,) = _read(tmp_path, entry)
    assert all(name.startswith('strchive_') for name in strchive.COLUMNS)
    assert locus.columns == {name: entry[name.removeprefix('strchive_')] for name in strchive.COLUMNS}
    # A repeat without a locus has a null in each.
    assert all(field.nullable for field in strchive.FIELDS)
