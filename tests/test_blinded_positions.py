from __future__ import annotations

import hashlib
import json
import pathlib
from collections.abc import Callable, Collection, Mapping, Sequence
from typing import Any

import duckdb
import invented_sv_annotation
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import test_family_tables

from case_annotation import annotation_tables, blinded_positions, cli, family_tables, sv_annotation, sv_sites

# Everything here is invented: the alleles, the records of ClinVar, the dbSNP ids and the positions.

_PRIMARY = family_tables.PRIMARY


def _file(
    path: pathlib.Path, small: Sequence[tuple[str, int, int]] = (), structural: Sequence[tuple[str, int, int]] = ()
) -> pathlib.Path:
    """A file of blinded positions with these intervals, each as its chromosome, its first and its last base."""

    def entries(intervals: Sequence[tuple[str, int, int]]) -> list[dict[str, Any]]:
        return [{'chrom': chrom, 'start': start, 'end': end} for chrom, start, end in intervals]

    held = {
        'assembly': 'GRCh38',
        'coordinates': '1-based, both ends inclusive',
        'small': entries(small),
        'structural': entries(structural),
    }
    path.write_text(json.dumps(held))
    return path


# ---- the file ----


def test_a_file_gives_the_intervals_of_its_two_lists(tmp_path: pathlib.Path) -> None:
    path = _file(tmp_path / 'positions.json', small=[('chr1', 100, 100), ('chrX', 7, 9)], structural=[('chrM', 1, 50)])
    positions = blinded_positions.read(path, _PRIMARY)
    assert [(i.chrom, i.start, i.end) for i in positions.small] == [('chr1', 100, 100), ('chrX', 7, 9)]
    assert [(i.chrom, i.start, i.end) for i in positions.structural] == [('chrM', 1, 50)]
    assert [(kind, place) for kind, place, _interval in positions.entries()] == [
        ('small', 0),
        ('small', 1),
        ('structural', 0),
    ]


def test_a_file_without_an_entry_is_read(tmp_path: pathlib.Path) -> None:
    positions = blinded_positions.read(_file(tmp_path / 'positions.json'), _PRIMARY)
    assert (positions.small, positions.structural) == ((), ())


_GOOD: Mapping[str, Any] = {
    'assembly': 'GRCh38',
    'coordinates': '1-based, both ends inclusive',
    'small': [{'chrom': 'chr1', 'start': 100, 'end': 100}],
    'structural': [{'chrom': 'chr2', 'start': 1001, 'end': 5000}],
}


def _with(key: str, value: object) -> dict[str, Any]:
    return {**_GOOD, key: value}


_REFUSED_FILES: Mapping[str, tuple[str, str]] = {
    'text that is no JSON': ('chr1\t100\t100\n', 'not a file of blinded positions'),
    'a list in place of the object': (json.dumps([_GOOD]), 'not a file of blinded positions'),
    'lists nested deeper than a reader of JSON follows': ('[' * 100_000, 'not a file of blinded positions'),
    'another assembly': (json.dumps(_with('assembly', 'GRCh37')), 'assembly'),
    'no assembly': (json.dumps({k: v for k, v in _GOOD.items() if k != 'assembly'}), 'assembly'),
    'other coordinates': (json.dumps(_with('coordinates', '0-based, half-open')), 'coordinates'),
    'no list of structural variants': (json.dumps({k: v for k, v in _GOOD.items() if k != 'structural'}), 'structural'),
    'a key that the file does not have': (json.dumps(_with('family', 'fam-1')), 'family'),
    'an entry with a key that it does not have': (
        json.dumps(_with('small', [{'chrom': 'chr1', 'start': 100, 'end': 100, 'alt': 'T'}])),
        'alt',
    ),
    'a list twice, of which a reader of JSON keeps the second': (
        json.dumps(_GOOD).replace('"structural"', '"small": [], "structural"'),
        r"the keys \['small'\] more than once",
    ),
    'an end twice in one entry': (
        json.dumps(_GOOD).replace('"end": 100', '"end": 100, "end": 200'),
        r"the keys \['end'\] more than once",
    ),
    'an entry without its end': (json.dumps(_with('small', [{'chrom': 'chr1', 'start': 100}])), 'end'),
    'an entry that ends before it starts': (
        json.dumps(_with('structural', [{'chrom': 'chr2', 'start': 5000, 'end': 1001}])),
        'ends before it starts',
    ),
    'an entry that starts before the first base': (
        json.dumps(_with('small', [{'chrom': 'chr1', 'start': 0, 'end': 1}])),
        'start',
    ),
    'a position written as text': (
        json.dumps(_with('small', [{'chrom': 'chr1', 'start': '100', 'end': 100}])),
        'start',
    ),
    'a position past the last that a table holds': (
        json.dumps(_with('small', [{'chrom': 'chr1', 'start': 100, 'end': 2**31}])),
        'end',
    ),
    'a position with a fraction': (json.dumps(_with('small', [{'chrom': 'chr1', 'start': 100, 'end': 100.5}])), 'end'),
    'a chromosome without its prefix': (
        json.dumps(_with('small', [{'chrom': '1', 'start': 100, 'end': 100}])),
        r"small\[0\] is on '1', which is no chromosome of the tables",
    ),
    'a contig that is no primary chromosome': (
        json.dumps(_with('structural', [{'chrom': 'chr1_KI270706v1_random', 'start': 1, 'end': 9}])),
        r"structural\[0\] is on 'chr1_KI270706v1_random'",
    ),
}


@pytest.mark.parametrize('case', _REFUSED_FILES)
def test_a_file_that_is_not_one_of_blinded_positions_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    text, message = _REFUSED_FILES[case]
    path = tmp_path / 'positions.json'
    path.write_text(text)
    with pytest.raises(ValueError, match=message) as refused:
        blinded_positions.read(path, _PRIMARY)
    assert str(path) in str(refused.value)
    # The good file, which each case changes in one place, is read.
    path.write_text(json.dumps(_GOOD))
    blinded_positions.read(path, _PRIMARY)


# ---- the rows of variants ----


def _known(tag: str) -> dict[str, Any]:
    """A value in every column that an entry empties, and a dbSNP id, each of which carries the tag."""
    record = {
        'accession': f'VCV-{tag}',
        'version': '1',
        'record_type': 'classified',
        'variant_type': 'SNV',
        'classification': f'Pathogenic-{tag}',
        'review_status': 'criteria provided, single submitter',
        'stars': 1,
        'last_evaluated': '2099-01-01',
        'conditions': [f'an invented condition of {tag}'],
    }
    return {
        'rsids': [f'rs-{tag}'],
        'clinvar_accession': record['accession'],
        'clinvar_classification': record['classification'],
        'clinvar_review_status': record['review_status'],
        'clinvar_stars': record['stars'],
        'clinvar_conditions': record['conditions'],
        'clinvar_last_evaluated': record['last_evaluated'],
        'clinvar_last_updated': '2099-02-02',
        'clinvar_pubmed_ids': [f'pmid-{tag}'],
        'clinvar_records': [record],
        'clinvar_other_alleles': [record | {'accession': f'VCV-other-{tag}', 'ref_allele': 'C', 'alt_allele': 'G'}],
    }


def test_the_invented_values_fill_every_column_that_an_entry_empties() -> None:
    known = _known('x')
    assert set(known) == {*blinded_positions.VARIANTS_COLUMNS, 'rsids'}
    assert all(value is not None for value in known.values())
    # The ClinVar columns alone: an allele keeps its dbSNP ids.
    assert set(blinded_positions.VARIANTS_COLUMNS) == {
        name for name in annotation_tables.VARIANTS.names if name.startswith('clinvar_')
    }
    assert len(blinded_positions.VARIANTS_COLUMNS) == 10
    assert set(blinded_positions.SV_COLUMNS) == {
        name for name in sv_annotation.COLUMNS.names if name.startswith('clinvar_')
    }


# A family's alleles around three invented variants: a single-base change at chr1:100, at a site of two alleles; an
# insertion written at chr2:200; and a deletion of the bases chr10:301 and chr10:302, written at chr10:300. Around
# each are substitutions, insertions and deletions at the bases that its entry has to reach and at the nearest
# bases that it must leave alone.
_AROUND_THE_BASE = (
    'chr1-90-ACCCCCCCCCCC-A',
    'chr1-98-G-GA',
    'chr1-98-G-T',
    'chr1-99-A-AT',
    'chr1-99-A-C',
    'chr1-99-AC-A',
    'chr1-100-C-CA',
    'chr1-100-C-G',
    'chr1-100-C-T',
    'chr1-100-CG-C',
    'chr1-101-G-A',
    'chr1-101-G-GT',
    'chr1-101-GT-G',
    'chr1-102-T-C',
)
_AROUND_THE_INSERTION = (
    'chr2-198-T-C',
    'chr2-198-T-TG',
    'chr2-199-G-C',
    'chr2-199-G-GT',
    'chr2-199-GA-G',
    'chr2-200-A-ACC',
    'chr2-200-A-G',
    'chr2-200-AT-A',
    'chr2-201-T-G',
    'chr2-201-T-TA',
    'chr2-201-TC-T',
    'chr2-202-C-CG',
    'chr2-202-C-T',
)
_AROUND_THE_DELETION = (
    'chr10-296-TTT-T',
    'chr10-297-TTA-T',
    'chr10-298-T-C',
    'chr10-298-T-TG',
    'chr10-299-A-AT',
    'chr10-299-A-G',
    'chr10-300-G-T',
    'chr10-300-GCA-G',
    'chr10-301-C-T',
    'chr10-302-A-AG',
    'chr10-302-A-T',
    'chr10-302-AC-A',
    'chr10-303-C-CT',
    'chr10-303-C-G',
    'chr10-303-CC-C',
    'chr10-304-C-G',
)
_ELSEWHERE = ('chrX-100-C-T',)
_FAMILY = (*_AROUND_THE_BASE, *_AROUND_THE_INSERTION, *_AROUND_THE_DELETION, *_ELSEWHERE)


def _family(
    tmp_path: pathlib.Path,
    alleles: Sequence[str] = _FAMILY,
    without_records: Collection[str] = (),
    sources: Sequence[str] = ('Ensembl', 'dbSNP', 'clinvar'),
) -> dict[str, Any]:
    """The arguments of assemble for a family with these alleles, in a union that has records of each of them.

    Every allele has a dbSNP id and a value in each column that an entry empties, which carry the allele's id. The
    alleles of `without_records` have the dbSNP id alone, as an allele that ClinVar does not hold.
    """
    inputs = test_family_tables._inputs(tmp_path)
    carried = (test_family_tables._REF, test_family_tables._REF, test_family_tables._HET)
    placed = [(allele, *allele.split('-')) for allele in alleles]
    test_family_tables._write_vcf(
        tmp_path / 'alleles.vcf.gz', [(chrom, int(pos), ref, alt, carried) for _allele, chrom, pos, ref, alt in placed]
    )
    union = tmp_path / 'alleles-union'
    (union / 'variants').mkdir(parents=True)
    (union / 'transcripts').mkdir()
    rows = [
        {'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt, 'variant_id': allele, 'worst_impact': 'modifier'}
        | ({'rsids': _known(allele)['rsids']} if allele in without_records else _known(allele))
        for allele, chrom, pos, ref, alt in placed
    ]
    pq.write_table(pa.Table.from_pylist(rows, schema=annotation_tables.VARIANTS), union / 'variants' / 'p0.parquet')
    pq.write_table(annotation_tables.TRANSCRIPTS.empty_table(), union / 'transcripts' / 'p0.parquet')
    (union / 'annotation.json').write_text(test_family_tables._info(*sources).model_dump_json())
    test_family_tables._write_joint(tmp_path / 'alleles-joint', dict.fromkeys(alleles))
    test_family_tables._write_joined(tmp_path / 'alleles-joined', alleles, ())
    test_family_tables._write_plugins(tmp_path / 'alleles-plugins', (), test_family_tables._plugins_info(None))
    return {
        **inputs,
        'vcf_path': tmp_path / 'alleles.vcf.gz',
        'annotation': union,
        'archive_values': tmp_path / 'alleles-joint',
        'joined': tmp_path / 'alleles-joined',
        'plugins': tmp_path / 'alleles-plugins',
    }


def _assemble(inputs: Mapping[str, Any], **changes: Any) -> family_tables.FamilyManifest:  # noqa: ANN401
    """Writes a family's folder from the inputs, with some of them changed."""
    return family_tables.assemble(**{**inputs, **changes})


def _emptied(path: pathlib.Path, key: str, columns: Sequence[str]) -> set[str]:
    """The rows of a table, by their key, that hold a null in every one of the columns; a row with some fails."""
    emptied = set()
    for row in pq.read_table(path).to_pylist():
        nulls = [row[column] is None for column in columns]
        assert all(nulls) or not any(nulls), row
        if all(nulls):
            emptied.add(row[key])
    return emptied


# What each entry reaches of the family's alleles: the alleles whose footprint shares a base with it. A footprint is
# the bases of the allele's REF, and for an insertion the base after them too.
_REACHED: Mapping[str, tuple[tuple[str, int, int], set[str]]] = {
    'a single base, at a site of two alleles': (
        ('chr1', 100, 100),
        # Both alleles of the site, the insertions written at the base and at the base before it, and the
        # deletions whose REF holds the base. Not the changes of the two neighbouring bases, the insertion written
        # at the next base or two bases before, or the deletion written at the next base.
        {
            'chr1-90-ACCCCCCCCCCC-A',
            'chr1-99-A-AT',
            'chr1-99-AC-A',
            'chr1-100-C-CA',
            'chr1-100-C-G',
            'chr1-100-C-T',
            'chr1-100-CG-C',
        },
    ),
    'an insertion, by its anchor base and the next base': (
        ('chr2', 200, 201),
        # The insertion, the changes of its anchor base and of the base after it, the insertions written one base
        # to either side, and the deletions whose REF holds one of the two bases. Not the change of the base before
        # the anchor, or anything written two bases from the anchor.
        {
            'chr2-199-G-GT',
            'chr2-199-GA-G',
            'chr2-200-A-ACC',
            'chr2-200-A-G',
            'chr2-200-AT-A',
            'chr2-201-T-G',
            'chr2-201-T-TA',
            'chr2-201-TC-T',
        },
    ),
    'an insertion, by its anchor base alone': (
        ('chr2', 200, 200),
        # An entry without the base after the anchor leaves the alleles written at that base alone, among them
        # the change of it.
        {'chr2-199-G-GT', 'chr2-199-GA-G', 'chr2-200-A-ACC', 'chr2-200-A-G', 'chr2-200-AT-A'},
    ),
    'a deletion, from its anchor base to its last base': (
        ('chr10', 300, 302),
        # The deletion, the changes of its anchor base and of the bases that it removes, the insertion written at
        # the base before it, and the insertion and the deletion written at its last base. Not the alleles written
        # at the base after its last base, the change of the base before it, or the deletion that ends there.
        {
            'chr10-299-A-AT',
            'chr10-300-G-T',
            'chr10-300-GCA-G',
            'chr10-301-C-T',
            'chr10-302-A-AG',
            'chr10-302-A-T',
            'chr10-302-AC-A',
        },
    ),
    'a base on another chromosome than any allele': (('chr3', 100, 100), set()),
    'a base far from any allele': (('chr1', 5000, 5000), set()),
}


@pytest.mark.parametrize('case', _REACHED)
def test_an_entry_of_a_small_variant_empties_the_alleles_whose_footprint_shares_a_base_with_it(
    tmp_path: pathlib.Path, case: str
) -> None:
    interval, reached = _REACHED[case]
    inputs = _family(tmp_path)
    _assemble(inputs, output=tmp_path / 'plain')
    manifest = _assemble(
        inputs, output=tmp_path / 'blinded', blinded=_file(tmp_path / 'positions.json', small=[interval])
    )
    assert _emptied(tmp_path / 'blinded' / 'variants.parquet', 'variant_id', blinded_positions.VARIANTS_COLUMNS) == (
        reached
    )
    assert reached <= set(_FAMILY)
    # Every other value of every row is the one of a run without the file, and so is every other table.
    plain, blinded = (pq.read_table(tmp_path / run / 'variants.parquet') for run in ('plain', 'blinded'))
    assert plain.drop_columns(list(blinded_positions.VARIANTS_COLUMNS)).equals(
        blinded.drop_columns(list(blinded_positions.VARIANTS_COLUMNS))
    )
    kept = [row for row in blinded.to_pylist() if row['variant_id'] not in reached]
    assert kept == [row for row in plain.to_pylist() if row['variant_id'] not in reached]
    assert all(row['clinvar_accession'] == f'VCV-{row["variant_id"]}' for row in kept)
    assert manifest.alleles == len(_FAMILY)


# The shapes of an allele, each as its REF and its ALT. The bases are invented and no reference is read.
_SHAPES: Mapping[str, tuple[str, str]] = {
    'a change of one base': ('A', 'G'),
    'a change of two bases in a row': ('AC', 'GT'),
    'a change of three bases in a row': ('ACG', 'TGA'),
    'an insertion of one base': ('A', 'AT'),
    'an insertion of two bases': ('A', 'ATT'),
    'an insertion whose first base is written in lower case': ('a', 'AT'),
    'a deletion of one base': ('AC', 'A'),
    'a deletion of three bases': ('ACGT', 'A'),
    'two bases replaced by one, after a shared base': ('ACG', 'AT'),
    'one base replaced by two, after a shared base': ('AC', 'ATT'),
    'two bases replaced by one': ('AC', 'G'),
    'one base replaced by two': ('A', 'GT'),
}


def _is_insertion(ref: str, alt: str) -> bool:
    return len(alt) > len(ref) and alt.upper().startswith(ref.upper())


def _begin(pos: int, ref: str, alt: str) -> int:
    """The base at which the Annotator places an allele: its `pos` after the leading bases that REF and ALT share."""
    shared = 0
    while shared < min(len(ref), len(alt)) and ref[shared].upper() == alt[shared].upper():
        shared += 1
    return pos + shared


def _footprint(pos: int, ref: str, alt: str) -> tuple[int, int]:
    """The first and the last base of an allele's footprint, which is the entry of a file for the allele."""
    return pos, pos + len(ref) - 1 + _is_insertion(ref, alt)


def test_the_rule_takes_every_pair_that_begins_at_one_base_or_shares_a_base_of_ref_and_little_beyond() -> None:
    # Measured, the Annotator lists a ClinVar record on a row where the two alleles begin at the same base. The
    # rule must take each such pair, whichever of the two the entry is of, and each pair whose REFs share a base.
    # It is written with footprints instead, which take one kind of pair more: an insertion with an allele whose
    # REF begins at the base after the insertion's REF and that begins a base later itself, as another insertion,
    # a deletion or a replacement written with the base before it does. Every pair of two shapes at up to six
    # bases apart is compared here.
    entry_pos = 100
    rows = [(shape, entry_pos + offset) for shape in _SHAPES for offset in range(-6, 7)]
    con = duckdb.connect()
    con.register(
        'a',
        pa.table(
            {
                'shape': pa.array([shape for shape, _pos in rows]),
                'chrom': pa.array(['chr1'] * len(rows)),
                'pos': pa.array([pos for _shape, pos in rows], pa.int32()),
                'ref': pa.array([_SHAPES[shape][0] for shape, _pos in rows]),
                'alt': pa.array([_SHAPES[shape][1] for shape, _pos in rows]),
            }
        ),
    )
    footprints = {shape: _footprint(entry_pos, *allele) for shape, allele in _SHAPES.items()}
    con.register(
        'e',
        pa.table(
            {
                'kind': pa.array(['small'] * len(footprints)),
                'shape': pa.array(list(footprints)),
                'chrom': pa.array(['chr1'] * len(footprints)),
                'first_base': pa.array([first for first, _last in footprints.values()], pa.int64()),
                'last_base': pa.array([last for _first, last in footprints.values()], pa.int64()),
            }
        ),
    )
    taken = set(
        con.sql(f'SELECT e.shape, a.shape, a.pos FROM a JOIN e ON {blinded_positions._SHARES_A_BASE}').fetchall()  # noqa: S608 (a condition of the package's own)
    )

    def spans(pos: int, ref: str) -> range:
        return range(pos, pos + len(ref))

    def overlap(one: range, other: range) -> bool:
        return one.start < other.stop and other.start < one.stop

    def adjacent(insertion: tuple[int, str, str], other: tuple[int, str, str]) -> bool:
        """Whether an insertion stands directly before the REF of an allele that begins a base later itself."""
        pos, ref, alt = insertion
        return _is_insertion(ref, alt) and other[0] == pos + len(ref) and _begin(*other) == other[0] + 1

    measured, beyond, apart = set(), set(), set()
    for entry_shape, entry_allele in _SHAPES.items():
        entry = (entry_pos, *entry_allele)
        for row_shape, pos in rows:
            row = (pos, *_SHAPES[row_shape])
            pair = (entry_shape, row_shape, pos)
            if _begin(*entry) == _begin(*row) or overlap(spans(*entry[:2]), spans(*row[:2])):
                measured.add(pair)
            elif adjacent(entry, row) or adjacent(row, entry):
                beyond.add(pair)
            if not overlap(range(entry_pos - 1, entry_pos + len(entry[1]) + 1), spans(*row[:2])):
                apart.add(pair)
    assert measured <= taken
    assert taken - measured == beyond
    # A row whose REF lies two bases or more from the REF of the entry's variant is never taken.
    assert apart
    assert not apart & taken
    # The pairs of the measurement that the REFs alone do not find, and one pair of each kind that the rule takes
    # beyond the measured ones.
    assert {
        ('an insertion of one base', 'a change of one base', entry_pos + 1),
        ('a change of one base', 'an insertion of one base', entry_pos - 1),
        ('a change of one base', 'an insertion whose first base is written in lower case', entry_pos - 1),
        ('an insertion of two bases', 'two bases replaced by one', entry_pos + 1),
    } <= measured
    assert {
        ('an insertion of one base', 'an insertion of two bases', entry_pos + 1),
        ('an insertion of one base', 'a deletion of one base', entry_pos + 1),
        ('a deletion of three bases', 'an insertion of one base', entry_pos - 1),
        ('an insertion of one base', 'two bases replaced by one, after a shared base', entry_pos + 1),
    } <= beyond


def test_an_entry_of_a_structural_variant_empties_no_row_of_variants(tmp_path: pathlib.Path) -> None:
    inputs = _family(tmp_path)
    over_everything = [(chrom, 1, 10_000) for chrom in ('chr1', 'chr2', 'chr10', 'chrX')]
    _assemble(inputs, blinded=_file(tmp_path / 'positions.json', structural=over_everything))
    assert _emptied(tmp_path / 'out' / 'variants.parquet', 'variant_id', blinded_positions.VARIANTS_COLUMNS) == set()


def _sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _everything(path: pathlib.Path) -> str:
    """What a reader gets from a Parquet file: its rows, its schema with its metadata, and its statistics."""
    file = pq.ParquetFile(path)
    return f'{file.read().to_pylist()} {file.schema_arrow.metadata} {file.metadata.to_dict()}'


def _folder(path: pathlib.Path) -> dict[str, bytes]:
    """The files of a family's folder, by name, each with its bytes."""
    return {file.name: file.read_bytes() for file in path.iterdir()}


def test_a_blinded_table_is_the_file_of_a_union_that_never_had_the_records(tmp_path: pathlib.Path) -> None:
    entries = [interval for interval, _reached in _REACHED.values()]
    reached = set().union(*(reached for _interval, reached in _REACHED.values()))
    (tmp_path / 'with').mkdir()
    (tmp_path / 'without').mkdir()
    blinded = _assemble(_family(tmp_path / 'with'), blinded=_file(tmp_path / 'positions.json', small=entries))
    never = _assemble(_family(tmp_path / 'without', without_records=reached))
    # The same bytes: data pages, dictionaries, statistics, Bloom filters and metadata.
    assert {t.table: t.sha256 for t in blinded.tables} == {t.table: t.sha256 for t in never.tables}
    path = tmp_path / 'with' / 'out' / 'variants.parquet'
    assert _sha256(path) == next(t.sha256 for t in blinded.tables if t.table == 'variants')
    # Nothing that a reader gets from the file names an emptied record, and the kept ones are there.
    read = _everything(path)
    kept = set(_FAMILY) - reached
    assert kept
    for allele in reached:
        assert not any(tag in read for tag in (f'VCV-{allele}', f'VCV-other-{allele}', f'pmid-{allele}'))
        assert f'rs-{allele}' in read
    for allele in kept:
        assert all(tag in read for tag in (f'VCV-{allele}', f'VCV-other-{allele}', f'rs-{allele}', f'pmid-{allele}'))
    # The two manifests are the same, key for key.
    assert blinded == never


def test_a_file_without_an_entry_leaves_the_folder_as_a_run_without_a_file_writes_it(tmp_path: pathlib.Path) -> None:
    inputs = _family(tmp_path)
    plain = _assemble(inputs, output=tmp_path / 'plain')
    empty = _assemble(inputs, output=tmp_path / 'empty', blinded=_file(tmp_path / 'positions.json'))
    assert plain == empty
    assert _folder(tmp_path / 'plain') == _folder(tmp_path / 'empty')


def test_the_manifest_holds_nothing_of_the_file(tmp_path: pathlib.Path) -> None:
    inputs = _family(tmp_path)
    one = _file(tmp_path / 'one-entry.json', small=[('chr1', 100, 100)])
    three = _file(
        tmp_path / 'three-entries.json',
        small=[('chr1', 100, 100), ('chr10', 77_777, 77_777)],
        structural=[('chr2', 88_888, 99_999)],
    )
    _assemble(inputs, output=tmp_path / 'one', blinded=one)
    _assemble(inputs, output=tmp_path / 'three', blinded=three)
    manifests = [(tmp_path / run / 'manifest.json').read_text() for run in ('one', 'three')]
    # Two files that empty the same rows give the same manifest: nothing in it depends on the file's name, on its
    # size, on its bytes or on how many entries it holds.
    assert manifests[0] == manifests[1]
    for path in (one, three):
        assert path.name not in manifests[1]
        assert path.stem not in manifests[1]
        assert hashlib.sha256(path.read_bytes()).hexdigest() not in manifests[1]
    assert not any(position in manifests[1] for position in ('77777', '88888', '99999'))
    # The folder holds the tables, the manifest and the notes, and no copy of the file.
    assert sorted(p.name for p in (tmp_path / 'three').iterdir()) == sorted(
        ['manifest.json', 'tables.md', 'variants.parquet', 'transcripts.parquet', 'roh.parquet', 'upd.parquet']
    )


def test_the_union_and_another_familys_tables_keep_what_one_familys_file_empties(tmp_path: pathlib.Path) -> None:
    inputs = _family(tmp_path)
    union = sorted((inputs['annotation'] / 'variants').iterdir())
    before = [_sha256(path) for path in union]
    positions = _file(tmp_path / 'positions.json', small=[('chr1', 100, 100)])
    # Two families of one union that hold the same alleles: one with a file, and one without.
    _assemble(inputs, output=tmp_path / 'with', blinded=positions)
    _assemble(inputs, output=tmp_path / 'without')
    assert [_sha256(path) for path in union] == before
    assert _emptied(tmp_path / 'with' / 'variants.parquet', 'variant_id', blinded_positions.VARIANTS_COLUMNS)
    assert not _emptied(tmp_path / 'without' / 'variants.parquet', 'variant_id', blinded_positions.VARIANTS_COLUMNS)


def test_the_step_logs_how_many_rows_each_entry_reached_and_an_entry_that_reached_none_is_no_error(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    positions = _file(
        tmp_path / 'positions.json', small=[('chr1', 100, 100), ('chr1', 5000, 5000)], structural=[('chr2', 1, 900)]
    )
    _assemble(_family(tmp_path), blinded=positions)
    # An entry by its list and its place in it, and no position.
    assert capsys.readouterr().out.splitlines() == [
        'blinded positions: small[0]: variants 7',
        'blinded positions: small[1]: variants 0',
        'blinded positions: structural[0]: no table that it could reach',
    ]
    (tmp_path / 'again').mkdir()
    _assemble(_family(tmp_path / 'again'))
    assert capsys.readouterr().out == ''


def test_a_blinded_row_keeps_its_dbsnp_ids(tmp_path: pathlib.Path) -> None:
    inputs = _family(tmp_path)
    _assemble(inputs, blinded=_file(tmp_path / 'positions.json', small=[interval for interval, _ in _REACHED.values()]))
    rows = pq.read_table(tmp_path / 'out' / 'variants.parquet').to_pylist()
    emptied = [row for row in rows if row['clinvar_records'] is None]
    assert len(emptied) > 10
    assert all(row['rsids'] == [f'rs-{row["variant_id"]}'] for row in rows)


def test_a_run_that_did_not_name_clinvar_has_no_column_to_empty(tmp_path: pathlib.Path) -> None:
    inputs = _family(tmp_path, without_records=_FAMILY, sources=('Ensembl', 'dbSNP'))
    positions = _file(tmp_path / 'positions.json', small=[('chr1', 100, 100)])
    manifest = _assemble(inputs, blinded=positions)
    assert set(blinded_positions.VARIANTS_COLUMNS) <= set(manifest.omitted_columns)
    rows = pq.read_table(tmp_path / 'out' / 'variants.parquet').to_pylist()
    assert all(row['rsids'] for row in rows)


def test_a_malformed_file_stops_the_step_before_it_writes_a_folder(tmp_path: pathlib.Path) -> None:
    inputs = _family(tmp_path)
    (tmp_path / 'positions.json').write_text('{"assembly": "GRCh38"}')
    with pytest.raises(ValueError, match='not a file of blinded positions'):
        _assemble(inputs, blinded=tmp_path / 'positions.json')
    assert not inputs['output'].exists()


# ---- the rows of sv ----

_SINGLE = 'criteria provided, single submitter'


def _record(accession: str, begin: int, end: int, kind: str, **extra: object) -> dict[str, Any]:
    return invented_sv_annotation.clinvar_entry(
        accession, begin, end, kind, 'Pathogenic', _SINGLE, '2099-01-01', **extra
    )


def _with_records(
    without: Collection[str] = (), begins: Mapping[tuple[str, str], int] | None = None
) -> Callable[[str, sv_sites.Site], dict[str, Any]]:
    """An invented annotation in which six records of the family's callers have a ClinVar record, each named after it.

    They are the structural-variant caller's symbolic deletion, the copy-number caller's loss, its gain and its
    segment of two alleles, which the Annotator matched by overlap, and the two short sequences, which it annotated
    as small variants and which hold a record of their own allele. The records of `without` have none. `begins`
    gives the first base that the Annotator annotated of an allele, by its record's ID and the allele, where that
    is another than the first base of the record.
    """

    def annotated(sv_id: str, site: sv_sites.Site) -> dict[str, Any]:
        whole = invented_sv_annotation.overlaps(1.0)
        kind = {'<DEL>': 'copy_number_loss', '<DUP>': 'copy_number_gain'}.get(site.alt, 'deletion')
        matched = {'inv:DEL:1': 'deletion', 'inv:LOSS:1': kind, 'inv:GAIN:1': kind, 'inv:CNV:1': kind}
        variant: dict[str, Any] = {'transcripts': [test_family_tables._GENE]}
        sources: dict[str, Any] = {}
        if sv_id in without:
            pass
        elif sv_id in matched:
            accession = f'VCV-{sv_id}-{site.alt}'
            assert site.end is not None
            sources['clinvar'] = [_record(accession, site.pos + 1, site.end, matched[sv_id], **whole)]
        elif sv_id in ('inv:INS:1', 'inv:DEL:2'):
            variant['clinvar'] = [
                _record(f'VCV-{sv_id}', site.pos + 1, site.pos + 1, 'insertion', isAlleleSpecific=True)
            ]
        if sv_id == 'inv:DEL:1':
            entry = invented_sv_annotation.gnomad_entry('inv-gnomad-1', 1001, 5000, 'deletion', 0.25, **whole)
            sources['gnomad'] = [entry]
        if (sv_id, site.alt) in (begins or {}):
            variant['begin'] = (begins or {})[sv_id, site.alt]
        return invented_sv_annotation.position(site, variant, **sources)

    return annotated


_WITH_RECORDS = {'inv:DEL:1', 'inv:LOSS:1', 'inv:GAIN:1', 'inv:CNV:1', 'inv:INS:1', 'inv:DEL:2'}

# What each file reaches of the family's records. inv:DEL:1 is the structural-variant caller's deletion of
# chr1:1001-5000, inv:LOSS:1 the copy-number caller's loss of chr1:20001-30000, inv:GAIN:1 its gain of
# chr2:40001-90000, and inv:CNV:1 its segment of two alleles at chr10:5001-15000. inv:INS:1 is the insertion
# chr10-7-G-GA, which is a row of variants, and inv:DEL:2 a deletion written at chr2:300 whose normalised allele is
# chr2-298-TGAC-T, which is no row of variants.
_REACHED_RECORDS: Mapping[str, tuple[dict[str, list[tuple[str, int, int]]], set[str]]] = {
    'the interval of a record of each caller': (
        {'structural': [('chr1', 1001, 5000), ('chr1', 20_001, 30_000)]},
        {'inv:DEL:1', 'inv:LOSS:1'},
    ),
    'an interval twice the length of a record': ({'structural': [('chr1', 1001, 9000)]}, {'inv:DEL:1'}),
    'an interval of half a record': ({'structural': [('chr1', 3001, 5000)]}, {'inv:DEL:1'}),
    'an interval of a base less than half a record': ({'structural': [('chr1', 3002, 5000)]}, set()),
    'an interval a base longer than twice a record': ({'structural': [('chr1', 1001, 9001)]}, set()),
    'an interval that covers two records and far more': ({'structural': [('chr1', 1, 1_000_000)]}, set()),
    'the interval of a record on another chromosome': ({'structural': [('chr2', 1001, 5000)]}, set()),
    'the interval of a segment of two alleles': ({'structural': [('chr10', 5001, 15_000)]}, {'inv:CNV:1'}),
    'the interval of a gain': ({'structural': [('chr2', 40_001, 90_000)]}, {'inv:GAIN:1'}),
    'a short sequence that is a row of variants': ({'small': [('chr10', 7, 7)]}, {'inv:INS:1'}),
    'a short sequence that is no row of variants, by its normalised allele': (
        {'small': [('chr2', 298, 301)]},
        {'inv:DEL:2'},
    ),
    'a short insertion, by the base after its anchor base': ({'small': [('chr10', 8, 8)]}, {'inv:INS:1'}),
    'a short insertion, by the base before its anchor base': ({'small': [('chr10', 6, 6)]}, set()),
    'a short insertion, two bases past its anchor base': ({'small': [('chr10', 9, 9)]}, set()),
    'a short deletion, by the last base of its normalised allele': ({'small': [('chr2', 301, 301)]}, {'inv:DEL:2'}),
    'a short deletion, by the base past its normalised allele': ({'small': [('chr2', 302, 302)]}, set()),
    'a short deletion, by the base before its normalised allele': ({'small': [('chr2', 297, 297)]}, set()),
    'a small variant inside the records matched by overlap': (
        {'small': [('chr1', 2000, 2000), ('chr1', 25_000, 25_000), ('chr10', 9000, 9000)]},
        set(),
    ),
    'a small variant as long as a record matched by overlap': ({'small': [('chr1', 1001, 5000)]}, set()),
    'the intervals of the two short sequences as structural variants': (
        {'structural': [('chr10', 7, 8), ('chr2', 299, 301)]},
        set(),
    ),
}


def _larger(tmp_path: pathlib.Path, without: Collection[str] = (), **changes: Any) -> dict[str, Any]:  # noqa: ANN401
    return test_family_tables._larger(tmp_path, annotated=_with_records(without), **changes)


@pytest.mark.parametrize('case', _REACHED_RECORDS)
def test_an_entry_empties_the_clinvar_columns_of_the_records_it_reaches(tmp_path: pathlib.Path, case: str) -> None:
    entries, reached = _REACHED_RECORDS[case]
    inputs = test_family_tables._inputs(tmp_path) | _larger(tmp_path)
    _assemble(inputs, output=tmp_path / 'plain')
    _assemble(inputs, output=tmp_path / 'blinded', blinded=_file(tmp_path / 'positions.json', **entries))
    plain, blinded = (pq.read_table(tmp_path / run / 'sv.parquet') for run in ('plain', 'blinded'))
    columns = list(blinded_positions.SV_COLUMNS)
    with_records = {row['sv_id'] for row in plain.to_pylist() if row['clinvar_records'] is not None}
    assert with_records == _WITH_RECORDS
    kept = {row['sv_id'] for row in blinded.to_pylist() if row['clinvar_records'] is not None}
    assert with_records - kept == reached
    # A reached record holds a null in each of the eight columns, and every other value is as without the file.
    for row in blinded.to_pylist():
        if row['sv_id'] in reached:
            assert all(row[column] is None for column in columns), row
    assert plain.drop_columns(columns).equals(blinded.drop_columns(columns))
    assert [r for r in blinded.to_pylist() if r['sv_id'] not in reached] == [
        r for r in plain.to_pylist() if r['sv_id'] not in reached
    ]


# What an entry reaches of two records whose annotated interval is not their own. The Annotator annotated the
# symbolic deletion inv:DEL:1, of chr1:1001-5000, from chr1:2001, and the gain of the segment inv:CNV:1, of
# chr10:5001-15000, from another base than its loss, so that the segment has no annotated interval.
_BEGINS = {('inv:DEL:1', '<DEL>'): 2001, ('inv:CNV:1', '<DUP>'): 5002}
_REACHED_BY_ANNOTATED: Mapping[str, tuple[tuple[str, int, int], set[str]]] = {
    'over half of the annotated interval and under half of the record': (('chr1', 2001, 3600), {'inv:DEL:1'}),
    'over half of the record and under half of the annotated interval': (('chr1', 1001, 3100), {'inv:DEL:1'}),
    'under half of both': (('chr1', 2001, 3400), set()),
    'the interval of a record without an annotated interval': (('chr10', 5001, 15_000), {'inv:CNV:1'}),
    'an interval before a record without an annotated interval': (('chr10', 101, 200), set()),
    'an interval after a record without an annotated interval': (('chr10', 20_001, 30_000), set()),
}


@pytest.mark.parametrize('case', _REACHED_BY_ANNOTATED)
def test_an_entry_is_compared_with_the_interval_that_the_annotator_annotated_where_a_record_has_one(
    tmp_path: pathlib.Path, case: str
) -> None:
    interval, reached = _REACHED_BY_ANNOTATED[case]
    inputs = test_family_tables._inputs(tmp_path) | test_family_tables._larger(
        tmp_path, annotated=_with_records(begins=_BEGINS)
    )
    _assemble(inputs, blinded=_file(tmp_path / 'positions.json', structural=[interval]))
    rows = {row['sv_id']: row for row in pq.read_table(tmp_path / 'out' / 'sv.parquet').to_pylist()}
    assert (rows['inv:DEL:1']['annotated_pos'], rows['inv:DEL:1']['annotated_end_pos']) == (2000, 5000)
    assert (rows['inv:CNV:1']['annotated_pos'], rows['inv:CNV:1']['annotated_end_pos']) == (None, 15_000)
    assert _WITH_RECORDS - {sv_id for sv_id, row in rows.items() if row['clinvar_records'] is not None} == reached


def test_a_blinded_sv_table_is_the_file_of_an_annotation_that_never_had_the_records(tmp_path: pathlib.Path) -> None:
    reached = {'inv:DEL:1', 'inv:LOSS:1', 'inv:CNV:1', 'inv:INS:1', 'inv:DEL:2'}
    positions = _file(
        tmp_path / 'positions.json',
        small=[('chr10', 7, 7), ('chr2', 298, 301)],
        structural=[('chr1', 1001, 5000), ('chr1', 20_001, 30_000), ('chr10', 5001, 15_000)],
    )
    (tmp_path / 'with').mkdir()
    (tmp_path / 'without').mkdir()
    blinded = _assemble(test_family_tables._inputs(tmp_path / 'with') | _larger(tmp_path / 'with'), blinded=positions)
    never = _assemble(test_family_tables._inputs(tmp_path / 'without') | _larger(tmp_path / 'without', without=reached))
    # The insertion is a row of variants too, which the file empties in the first run and the second does not.
    tables = ('sv', 'sv_transcripts', 'transcripts', 'repeats', 'roh', 'upd')
    assert {t.table: t.sha256 for t in blinded.tables if t.table in tables} == {
        t.table: t.sha256 for t in never.tables if t.table in tables
    }
    read = _everything(tmp_path / 'with' / 'out' / 'sv.parquet')
    assert 'VCV-inv:GAIN:1' in read
    assert 'inv-gnomad-1' in read
    assert not any(f'VCV-{sv_id}' in read for sv_id in reached)
    rows = {r['variant_id']: r for r in pq.read_table(tmp_path / 'with' / 'out' / 'variants.parquet').to_pylist()}
    assert rows['chr10-7-G-GA']['clinvar_accession'] is None
    assert rows['chr2-100-A-G']['clinvar_accession'] is not None


def test_the_step_logs_the_rows_of_both_tables_that_an_entry_reached(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    positions = _file(tmp_path / 'positions.json', small=[('chr10', 7, 7)], structural=[('chr1', 1001, 5000)])
    _assemble(test_family_tables._inputs(tmp_path) | _larger(tmp_path), blinded=positions)
    assert capsys.readouterr().out.splitlines() == [
        'blinded positions: small[0]: variants 1, sv 1',
        'blinded positions: structural[0]: sv 1',
    ]


def test_a_run_that_did_not_name_clinvar_has_no_column_of_sv_to_empty(tmp_path: pathlib.Path) -> None:
    sources = ('Ensembl', 'gnomad')
    inputs = test_family_tables._inputs(tmp_path, info=test_family_tables._info(*sources))
    larger = test_family_tables._larger(tmp_path, sources=sources)
    positions = _file(tmp_path / 'positions.json', structural=[('chr1', 1001, 5000)])
    manifest = _assemble(inputs | larger, blinded=positions)
    assert set(blinded_positions.SV_COLUMNS) <= set(manifest.sv_omitted_columns)


def test_a_file_without_an_entry_leaves_the_folder_of_a_family_with_larger_variants_as_it_is(
    tmp_path: pathlib.Path,
) -> None:
    inputs = test_family_tables._inputs(tmp_path) | _larger(tmp_path)
    _assemble(inputs, output=tmp_path / 'plain')
    _assemble(inputs, output=tmp_path / 'empty', blinded=_file(tmp_path / 'positions.json'))
    # The records' annotation goes through a view of its own in a run with a file, and sv has the same bytes.
    assert 'sv.parquet' in _folder(tmp_path / 'plain')
    assert _folder(tmp_path / 'plain') == _folder(tmp_path / 'empty')


# ---- a folder with a file beside one without ----


def _values(document: object, path: tuple[str | int, ...] = ()) -> dict[tuple[str | int, ...], object]:
    """Every value of a JSON document, by the keys and the places in lists that lead to it.

    An object and a list stand under their own path too, as the name of their type, so that an empty one differs
    from a null.
    """
    if isinstance(document, dict):
        children: list[tuple[str | int, object]] = list(document.items())
    elif isinstance(document, list):
        children = list(enumerate(document))
    else:
        return {path: document}
    found: dict[tuple[str | int, ...], object] = {path: type(document).__name__}
    for key, child in children:
        found |= _values(child, (*path, key))
    return found


def test_a_folder_differs_from_that_of_a_run_without_the_file_in_the_emptied_values_alone(
    tmp_path: pathlib.Path,
) -> None:
    inputs = test_family_tables._inputs(tmp_path) | _larger(tmp_path)
    positions = _file(
        tmp_path / 'given-positions.json',
        small=[('chr1', 5, 5), ('chr10', 7, 7)],
        structural=[('chr1', 1001, 5000)],
    )
    _assemble(inputs, output=tmp_path / 'without')
    _assemble(inputs, output=tmp_path / 'with', blinded=positions)
    without, given = _folder(tmp_path / 'without'), _folder(tmp_path / 'with')
    # The two tables that an entry empties values of, each with its key, the columns that an entry empties and the
    # rows that the entries reach.
    emptied = {
        'variants': ('variant_id', blinded_positions.VARIANTS_COLUMNS, {'chr1-5-C-T', 'chr10-7-G-GA'}),
        'sv': ('sv_id', blinded_positions.SV_COLUMNS, {'inv:DEL:1', 'inv:INS:1'}),
    }
    # The same files, and every file but the two tables and the manifest with the same bytes, the notes among them.
    assert without.keys() == given.keys()
    changed = {'manifest.json', *(f'{table}.parquet' for table in emptied)}
    assert {'tables.md', 'transcripts.parquet', 'sv_transcripts.parquet', 'repeats.parquet'} <= without.keys() - changed
    assert {name for name in without if without[name] != given[name]} == changed
    for table, (key, columns, reached) in emptied.items():
        before, after = (pq.ParquetFile(tmp_path / run / f'{table}.parquet') for run in ('without', 'with'))
        # The same columns, of the same types, and the same metadata, the file's own and the schema's.
        assert after.schema_arrow.equals(before.schema_arrow, check_metadata=True)
        assert after.metadata.metadata == before.metadata.metadata
        assert after.metadata.num_row_groups == before.metadata.num_row_groups
        rows, read = before.read().to_pylist(), after.read().to_pylist()
        assert [row[key] for row in read] == [row[key] for row in rows]
        differing = {
            (row[key], column)
            for row, other in zip(rows, read, strict=True)
            for column in row
            if row[column] != other[column]
        }
        # The values that differ are the ones that a reached row held in the columns that an entry empties, and
        # each is a null with the file.
        held = {
            (row[key], column) for row in rows if row[key] in reached for column in columns if row[column] is not None
        }
        assert differing == held
        assert {row for row, _column in differing} == reached
        assert all(row[column] is None for row in read for column in columns if row[key] in reached)
    # The manifests have the same keys at every level and lists of the same lengths. Their values differ in the size
    # and the checksum of the two tables alone, which follow from the emptied values.
    manifest = json.loads(without['manifest.json'])
    was, now = _values(manifest), _values(json.loads(given['manifest.json']))
    assert was.keys() == now.keys()
    places = {listed['table']: place for place, listed in enumerate(manifest['tables'])}
    follow = {('tables', places[table], key) for table in emptied for key in ('bytes', 'sha256')}
    unlike = {path for path in was if was[path] != now[path]}
    assert unlike <= follow
    assert {path for path in follow if path[-1] == 'sha256'} <= unlike
    # No file of the folder holds the word, or the name or the checksum of the positions' file, in its bytes as they
    # are stored. Of a table that covers its column names and its key-value metadata, and the comparison above its
    # values.
    marks = (b'blind', positions.stem.encode(), _sha256(positions).encode())
    for name, content in given.items():
        assert not any(mark in content.lower() for mark in marks), name


# ---- the commands ----


def test_the_commands_take_the_file_and_the_check_refuses_one_that_is_malformed(tmp_path: pathlib.Path) -> None:
    inputs = _family(tmp_path)
    positions = _file(tmp_path / 'positions.json', small=[('chr1', 100, 100)])
    check = [*test_family_tables._check_command(inputs), '--proband', 'S3']
    cli.main(check)
    cli.main([*check, '--blinded-positions', str(positions)])
    tables = [*test_family_tables._tables_command(inputs), '--proband', 'S3', '--blinded-positions', str(positions)]
    cli.main(tables)
    assert _emptied(inputs['output'] / 'variants.parquet', 'variant_id', blinded_positions.VARIANTS_COLUMNS)
    other = tmp_path / 'other.json'
    other.write_text(json.dumps(_with('assembly', 'GRCh37')))
    with pytest.raises(ValueError, match='not a file of blinded positions'):
        cli.main([*check, '--blinded-positions', str(other)])
