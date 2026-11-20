"""A union's values from echtvar archives, from an invented union as echtvar would annotate it with such archives."""

from __future__ import annotations

import copy
import hashlib
import json
import pathlib
import struct
import zipfile
from collections.abc import Mapping, Sequence
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from case_annotation import alphagenome_archives, cli, echtvar_values

# The fields of the archive, as echtvar embeds them: gnomAD's field, the archive's alias for it and how it is stored.
_CONFIG: list[dict[str, Any]] = [
    {'field': 'AC_joint', 'alias': 'gnomad_AC', 'missing_value': -2147483648, 'multiplier': 1, 'ftype': 'Integer'},
    {'field': 'AF_joint', 'alias': 'gnomad_AF', 'missing_value': 0, 'multiplier': 2000000, 'ftype': 'Float'},
    {
        'field': 'nhomalt_joint',
        'alias': 'gnomad_HomAlt',
        'missing_value': -2147483648,
        'multiplier': 1,
        'ftype': 'Integer',
    },
    {
        'field': 'fafmax_faf95_max_joint',
        'alias': 'gnomad_faf95max',
        'missing_value': 0,
        'multiplier': 2000000,
        'ftype': 'Float',
    },
]
for _field in _CONFIG:
    _field |= {'missing_string': 'MISSING', 'description': 'invented', 'zigzag': False, 'number': 'A'}
_GNOMAD = echtvar_values.GNOMAD_JOINT
_COLUMNS = tuple(_GNOMAD.columns)

_INFO = (
    '##INFO=<ID=gnomad_AC,Number=1,Type=Integer,Description="added by echtvar invented">',
    '##INFO=<ID=gnomad_AF,Number=1,Type=Float,Description="added by echtvar invented">',
    '##INFO=<ID=gnomad_HomAlt,Number=1,Type=Integer,Description="added by echtvar invented">',
    '##INFO=<ID=gnomad_faf95max,Number=1,Type=Float,Description="added by echtvar invented">',
)
# What echtvar writes for an allele its archive lacks: no count, and a frequency of 0.
_ABSENT = 'gnomad_AC=.;gnomad_AF=0;gnomad_HomAlt=.;gnomad_faf95max=0'
# Per record: chrom, pos, ref, alt and the INFO echtvar added.
_RECORDS = [
    ('chr1', 5, 'C', 'T', _ABSENT),
    # An allele gnomAD holds with a count of 0.
    ('chr1', 9, 'G', 'A', 'gnomad_AC=0;gnomad_AF=0;gnomad_HomAlt=0;gnomad_faf95max=0'),
    # Frequencies as 32-bit floats print: 0.001 and 0.0912345.
    ('chr2', 100, 'A', 'G', 'gnomad_AC=1600;gnomad_AF=0.0010000000475;gnomad_HomAlt=2;gnomad_faf95max=0.09123449773'),
    # An allele of more than four bases, which echtvar looks up in its table of long alleles.
    ('chr2', 200, 'CTGATGA', 'C', 'gnomad_AC=7;gnomad_AF=4e-06;gnomad_HomAlt=0;gnomad_faf95max=1e-06'),
    ('chrM', 7, 'G', 'A', _ABSENT),
]


def _write_vcf(
    path: pathlib.Path,
    records: Sequence[tuple[str, int, str, str, str]],
    info: Sequence[str] = (),
    samples: Sequence[str] = (),
) -> pathlib.Path:
    lines = ['##fileformat=VCFv4.2', *(f'##contig=<ID={c},length=1000>' for c in ('chr1', 'chr2', 'chrM')), *info]
    if samples:
        lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="x">')
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    lines.append('\t'.join([*columns, *(['FORMAT', *samples] if samples else [])]))
    for chrom, pos, ref, alt, fields in records:
        calls = ['GT', *('0/1' for _ in samples)] if samples else []
        lines.append('\t'.join([chrom, str(pos), '.', ref, alt, '.', '.', fields, *calls]))
    path.write_text('\n'.join(lines) + '\n')
    return path


def _write_archive(path: pathlib.Path, config: object = _CONFIG) -> pathlib.Path:
    """An archive with this configuration, as JSON unless it is given as bytes."""
    with zipfile.ZipFile(path, 'w') as z:
        z.writestr('echtvar/config.json', config if isinstance(config, bytes) else json.dumps(config))
        z.writestr('echtvar/chr1/0/var32.bin', b'invented')
    return path


def _convert(
    tmp_path: pathlib.Path,
    annotated: Sequence[tuple[str, int, str, str, str]] = _RECORDS,
    union: Sequence[tuple[str, int, str, str, str]] | None = None,
    config: object = _CONFIG,
    info: Sequence[str] = _INFO,
    samples: Sequence[str] = (),
) -> echtvar_values.ArchiveInfo:
    sites = [(c, p, r, a, '.') for c, p, r, a, _ in (annotated if union is None else union)]
    return echtvar_values.convert(
        _write_vcf(tmp_path / 'annotated.vcf', annotated, info, samples),
        _write_vcf(tmp_path / 'union.vcf', sites),
        {_GNOMAD.name: _write_archive(tmp_path / 'invented.zip', config)},
        'echtvar 0.0.0',
        tmp_path / 'out',
    ).gnomad_joint


def _rows(tmp_path: pathlib.Path) -> dict[str, dict[str, Any]]:
    table = pq.read_table(tmp_path / 'out' / echtvar_values.TABLE_NAME)
    assert table.schema.equals(echtvar_values.schema([_GNOMAD]))
    return {row.pop('variant_id'): row for row in table.to_pylist()}


def _float32(value: float) -> float:
    """The 32-bit float nearest to the value, as Python reads it from a 32-bit column."""
    return struct.unpack('f', struct.pack('f', value))[0]


def test_each_allele_has_the_values_the_archive_stores(tmp_path: pathlib.Path) -> None:
    _convert(tmp_path)
    rows = _rows(tmp_path)
    assert list(rows) == ['chr1-5-C-T', 'chr1-9-G-A', 'chr2-100-A-G', 'chr2-200-CTGATGA-C', 'chrM-7-G-A']
    none = dict.fromkeys(_COLUMNS)
    # An allele the archive lacks has no values; one it holds with a count of 0 keeps its 0s.
    assert rows['chr1-5-C-T'] == none
    assert rows['chrM-7-G-A'] == none
    assert rows['chr1-9-G-A'] == dict.fromkeys(_COLUMNS, 0)
    # A frequency is the 32-bit float echtvar returns, in a 32-bit column: 0.001 reads as 0.0010000000475.
    assert rows['chr2-100-A-G'] == {
        'gnomad_joint_af': _float32(0.001),
        'gnomad_joint_ac': 1600,
        'gnomad_joint_hom': 2,
        'gnomad_joint_faf95_max': _float32(0.0912345),
    }
    assert rows['chr2-200-CTGATGA-C'] == {
        'gnomad_joint_af': _float32(4e-06),
        'gnomad_joint_ac': 7,
        'gnomad_joint_hom': 0,
        'gnomad_joint_faf95_max': _float32(1e-06),
    }
    frequencies = {f.name: f.type for f in _GNOMAD.fields if pa.types.is_floating(f.type)}
    assert frequencies == {'gnomad_joint_af': pa.float32(), 'gnomad_joint_faf95_max': pa.float32()}


def test_the_archive_is_recorded_with_its_digest_and_fields(tmp_path: pathlib.Path) -> None:
    info = _convert(tmp_path)
    archive = tmp_path / 'invented.zip'
    assert info.archive == 'invented.zip'
    assert info.bytes == archive.stat().st_size
    assert info.sha256 == hashlib.sha256(archive.read_bytes()).hexdigest()
    assert info.echtvar_version == 'echtvar 0.0.0'
    assert [(f.field, f.alias, f.ftype, f.multiplier, f.missing_value) for f in info.fields] == [
        (c['field'], c['alias'], c['ftype'], c['multiplier'], c['missing_value']) for c in _CONFIG
    ]
    assert echtvar_values.read_info(tmp_path / 'out').gnomad_joint == info


def test_the_command_writes_the_table_and_the_record(tmp_path: pathlib.Path) -> None:
    annotated = _write_vcf(tmp_path / 'annotated.vcf', _RECORDS, _INFO)
    union = _write_vcf(tmp_path / 'union.vcf', [(c, p, r, a, '.') for c, p, r, a, _ in _RECORDS])
    archive = _write_archive(tmp_path / 'invented.zip')
    args = ['--annotated', str(annotated), '--sites', str(union), '--archive', _GNOMAD.name, str(archive)]
    cli.main(['echtvar-values', *args, '--echtvar-version', 'echtvar 0.0.0', '--output', str(tmp_path / 'out')])
    assert len(_rows(tmp_path)) == len(_RECORDS)
    assert echtvar_values.read_info(tmp_path / 'out').gnomad_joint.echtvar_version == 'echtvar 0.0.0'
    with pytest.raises(SystemExit, match='--archive is given twice for gnomad_joint'):
        cli.main(['echtvar-values', *args, *args[4:], '--echtvar-version', '0', '--output', str(tmp_path / 'again')])


def test_a_stand_in_no_frequency_can_be_leaves_no_value(tmp_path: pathlib.Path) -> None:
    config = copy.deepcopy(_CONFIG)
    config[1]['missing_value'] = -1
    record = ('chr1', 9, 'G', 'A', 'gnomad_AC=3;gnomad_AF=-1;gnomad_HomAlt=0;gnomad_faf95max=0')
    _convert(tmp_path, annotated=[record], config=config)
    assert _rows(tmp_path)['chr1-9-G-A'] == {
        'gnomad_joint_af': None,
        'gnomad_joint_ac': 3,
        'gnomad_joint_hom': 0,
        'gnomad_joint_faf95_max': 0,
    }


_OTHER = ('chr2', 300, 'T', 'A', _ABSENT)
ALLELES_REFUSED = {
    'an annotated allele the union lacks': (
        {'annotated': [*_RECORDS, _OTHER], 'union': _RECORDS},
        'echtvar output alleles the union lacks',
    ),
    'a union allele the annotated VCF lacks': (
        {'annotated': _RECORDS, 'union': [*_RECORDS, _OTHER]},
        'union alleles the echtvar output lacks',
    ),
    'an allele of other bases at the same position': (
        {'annotated': [('chr1', 5, 'C', 'G', _ABSENT)], 'union': [('chr1', 5, 'C', 'T', '.')]},
        'union alleles the echtvar output lacks',
    ),
    'an annotated allele repeated': (
        {'annotated': [*_RECORDS, _RECORDS[0]], 'union': _RECORDS},
        'the echtvar output repeats alleles',
    ),
    'a union allele repeated': ({'annotated': _RECORDS, 'union': [*_RECORDS, _RECORDS[0]]}, 'the union repeats'),
}


@pytest.mark.parametrize('case', ALLELES_REFUSED)
def test_annotated_alleles_other_than_the_unions_are_refused(tmp_path: pathlib.Path, case: str) -> None:
    changes, message = ALLELES_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, **changes)


def _config(index: int, **changes: object) -> list[dict[str, Any]]:
    config = copy.deepcopy(_CONFIG)
    config[index] |= changes
    return config


ARCHIVES_REFUSED = {
    'a field without its type': ([{k: v for k, v in f.items() if k != 'ftype'} for f in _CONFIG], 'not the fields'),
    'a field of a type echtvar has none of': (_config(0, ftype='Text'), 'not the fields'),
    'a configuration that is not a list of fields': ({'fields': _CONFIG}, 'not the fields'),
    'a configuration that is not JSON': (b'invented', 'not the fields'),
    "an archive without gnomAD's joint count": (_CONFIG[1:], "0 fields encode gnomAD's AC_joint"),
    'an archive of the genomes alone': (_config(1, field='AF_genomes'), "0 fields encode gnomAD's AF_joint"),
    'a field encoded twice': ([*_CONFIG, _config(2, alias='again')[2]], "2 fields encode gnomAD's nhomalt_joint"),
    'a frequency encoded as an integer': (_config(1, ftype='Integer'), 'AF_joint is encoded as Integer, not Float'),
    'a stand-in count that a count can be': (_config(0, missing_value=0), 'which a value of the field can be'),
}


@pytest.mark.parametrize('case', ARCHIVES_REFUSED)
def test_an_archive_without_the_joint_fields_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    config, message = ARCHIVES_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, config=config)


def test_a_file_that_is_no_echtvar_archive_is_refused(tmp_path: pathlib.Path) -> None:
    annotated = _write_vcf(tmp_path / 'annotated.vcf', _RECORDS, _INFO)
    not_a_zip = tmp_path / 'text.zip'
    not_a_zip.write_text('invented')
    with pytest.raises(ValueError, match=r'text\.zip: not a zip file'):
        echtvar_values.convert(annotated, annotated, {_GNOMAD.name: not_a_zip}, 'echtvar 0.0.0', tmp_path / 'out')
    without_config = tmp_path / 'other.zip'
    with zipfile.ZipFile(without_config, 'w') as z:
        z.writestr('readme.txt', 'invented')
    with pytest.raises(ValueError, match=r'other\.zip: no echtvar/config\.json'):
        echtvar_values.convert(annotated, annotated, {_GNOMAD.name: without_config}, 'echtvar 0.0.0', tmp_path / 'out')


def test_a_conversion_without_echtvars_version_is_refused(tmp_path: pathlib.Path) -> None:
    annotated = _write_vcf(tmp_path / 'annotated.vcf', _RECORDS, _INFO)
    union = _write_vcf(tmp_path / 'union.vcf', [(c, p, r, a, '.') for c, p, r, a, _ in _RECORDS])
    archive = _write_archive(tmp_path / 'invented.zip')
    with pytest.raises(ValueError, match='echtvar_version'):
        echtvar_values.convert(annotated, union, {_GNOMAD.name: archive}, '', tmp_path / 'out')


def _info_as(alias: str, number: str, kind: str) -> list[str]:
    declared = f'##INFO=<ID={alias},Number={number},Type={kind},Description="invented">'
    return [declared if f'ID={alias},' in line else line for line in _INFO]


ANNOTATED_REFUSED = {
    'a VCF with samples': ({'samples': ['S1']}, 'has samples'),
    'a VCF echtvar did not annotate': ({'info': _INFO[1:]}, 'no INFO field gnomad_AC'),
    'a record of two alleles': ({'annotated': [('chr1', 5, 'C', 'T,G', _ABSENT)]}, 'one alternate allele'),
    'a record without an alternate allele': ({'annotated': [('chr1', 5, 'C', '.', _ABSENT)]}, 'one alternate allele'),
    'a frequency on an allele the archive lacks': (
        {'annotated': [('chr1', 5, 'C', 'T', 'gnomad_AC=.;gnomad_AF=0.25;gnomad_HomAlt=.;gnomad_faf95max=0')]},
        r"lacks the allele, yet it has \['gnomad_AF'\]",
    ),
    'two values of a field': (
        {
            'annotated': [('chr1', 5, 'C', 'T', 'gnomad_AC=3,4;gnomad_AF=0;gnomad_HomAlt=0;gnomad_faf95max=0')],
            'info': _info_as('gnomad_AC', '.', 'Integer'),
        },
        'gnomad_AC has 2 values',
    ),
    'a count that is not an integer': (
        {
            'annotated': [('chr1', 5, 'C', 'T', 'gnomad_AC=3.5;gnomad_AF=0;gnomad_HomAlt=0;gnomad_faf95max=0')],
            'info': _info_as('gnomad_AC', '1', 'Float'),
        },
        'gnomad_AC is 3.5, not an integer',
    ),
    'a frequency that is not a number': (
        {
            'annotated': [('chr1', 5, 'C', 'T', 'gnomad_AC=3;gnomad_AF=rare;gnomad_HomAlt=0;gnomad_faf95max=0')],
            'info': _info_as('gnomad_AF', '1', 'String'),
        },
        "gnomad_AF is 'rare', not a number",
    ),
}


@pytest.mark.parametrize('case', ANNOTATED_REFUSED)
def test_an_annotated_vcf_echtvar_would_not_write_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    changes, message = ANNOTATED_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, **changes)


# ---- the archives of AlphaGenome's scores ----

_AVI = echtvar_values.ALPHAGENOME_AVI
_BREAKDOWN = echtvar_values.ALPHAGENOME_AVI_BREAKDOWN
_SPLICING = echtvar_values.ALPHAGENOME_SPLICING
_SCORES = (_AVI, _BREAKDOWN, _SPLICING)


def _embedded(source: echtvar_values.Source) -> list[dict[str, Any]]:
    """The configuration that echtvar embeds in an archive encoded with the configuration of the source's table."""
    table = alphagenome_archives.TABLES[source.name]
    return [
        {'missing_string': 'MISSING', 'ftype': 'Float', 'number': '1', 'description': 'invented', **field}
        for field in alphagenome_archives.config(table).values()
    ]


def _aliases(source: echtvar_values.Source) -> list[str]:
    return [field['alias'] for field in _embedded(source)]


def _added(source: echtvar_values.Source, *values: object) -> str:
    """What echtvar adds to a record from the source's archive: the values given, or none for an allele it lacks."""
    aliases = _aliases(source)
    return ';'.join(f'{alias}={value}' for alias, value in zip(aliases, values or ['.'] * len(aliases), strict=True))


# One contribution a step below a value of four digits, as echtvar cuts it, one below a millionth and one negative.
_SIXTEEN = ('0.2817', '0.0012929', '5.31e-05', '0', '-0.0004126', *['0.0125'] * 11)
# Per record: chrom, pos, ref, alt and what echtvar added from the archive of each source.
_SCORED = [
    # An allele that every archive holds; echtvar cut its raw score and its PHRED to the millionth below.
    (
        'chr1',
        5,
        'C',
        'T',
        {_AVI.name: ('-0.042449999', '97.313369'), _BREAKDOWN.name: _SIXTEEN, _SPLICING.name: ('1.52e-05',)},
    ),
    # A raw score and a PHRED of 0, which are values; no breakdown and no splicing score.
    ('chr1', 9, 'G', 'A', {_AVI.name: ('0', '0')}),
    # An insertion at a position whose single-base changes the archives hold, and an allele on chrM.
    ('chr1', 9, 'G', 'GA', {}),
    ('chrM', 7, 'G', 'A', {}),
]
_GNOMAD_HELD = 'gnomad_AC=7;gnomad_AF=4e-06;gnomad_HomAlt=0;gnomad_faf95max=1e-06'


def _scored_records(
    sources: Sequence[echtvar_values.Source], scored: Sequence[tuple[str, int, str, str, Any]] = _SCORED
) -> list[tuple[str, int, str, str, str]]:
    return [
        (c, p, r, a, ';'.join([_GNOMAD_HELD, *(_added(source, *values.get(source.name, ())) for source in sources)]))
        for c, p, r, a, values in scored
    ]


def _convert_scores(
    tmp_path: pathlib.Path,
    sources: Sequence[echtvar_values.Source] = _SCORES,
    scored: Sequence[tuple[str, int, str, str, Any]] = _SCORED,
    configs: Mapping[str, object] | None = None,
    declared: Sequence[echtvar_values.Source] | None = None,
) -> echtvar_values.Archives:
    """Converts a union that echtvar annotated with gnomAD's archive and those of the sources given."""
    records = _scored_records(sources, scored)
    info = [
        *_INFO,
        *(
            f'##INFO=<ID={alias},Number=1,Type=Float,Description="added by echtvar invented">'
            for source in (sources if declared is None else declared)
            for alias in _aliases(source)
        ),
    ]
    archives = {_GNOMAD.name: _write_archive(tmp_path / 'gnomad.zip')}
    for source in sources:
        config = (configs or {}).get(source.name, _embedded(source))
        archives[source.name] = _write_archive(tmp_path / f'{source.name}.zip', config)
    return echtvar_values.convert(
        _write_vcf(tmp_path / 'annotated.vcf', records, info),
        _write_vcf(tmp_path / 'union.vcf', [(c, p, r, a, '.') for c, p, r, a, _ in records]),
        archives,
        'echtvar 0.0.0',
        tmp_path / 'out',
    )


def _score_rows(tmp_path: pathlib.Path) -> dict[str, dict[str, Any]]:
    return {
        row.pop('variant_id'): row for row in pq.read_table(tmp_path / 'out' / echtvar_values.TABLE_NAME).to_pylist()
    }


def test_every_source_has_a_record_and_the_scores_come_from_the_fields_of_their_tables() -> None:
    assert list(echtvar_values.Archives.model_fields) == list(echtvar_values.SOURCES)
    columns = [column for source in echtvar_values.SOURCES.values() for column in source.columns]
    assert len(set(columns)) == len(columns)
    for source in _SCORES:
        table = alphagenome_archives.TABLES[source.name]
        assert {column.field for column in source.columns.values()} == set(table.kept)
        assert set(source.presence) == set(source.columns)
        # An archive that is encoded with the table's configuration is one that the source reads.
        fields = tuple(
            echtvar_values.ArchiveField(**{name: field[name] for name in echtvar_values.ArchiveField.model_fields})
            for field in _embedded(source)
        )
        assert (
            echtvar_values.column_fields(source, fields, pathlib.Path('invented.zip')).keys() == source.columns.keys()
        )
    # The two tables that are not for permissive use, and no other, give their values with a notice.
    assert {source.name for source in echtvar_values.SOURCES.values() if source.notice is not None} == {
        table.name
        for table in alphagenome_archives.TABLES.values()
        if any('description' in field for field in alphagenome_archives.config(table).values())
    }


def test_an_allele_has_the_scores_of_each_archive_that_holds_it(tmp_path: pathlib.Path) -> None:
    _convert_scores(tmp_path)
    rows = _score_rows(tmp_path)
    table = pq.read_schema(tmp_path / 'out' / echtvar_values.TABLE_NAME)
    assert table.equals(echtvar_values.schema([_GNOMAD, *_SCORES]))
    assert {table.field(column).type for source in _SCORES for column in source.columns} == {pa.float32()}
    held = rows['chr1-5-C-T']
    # A score is rounded to the decimals of its table, which gives the table's value back.
    assert (held['avi_raw'], held['avi_phred']) == (_float32(-0.04245), _float32(97.31337))
    assert [held[column] for column in list(_BREAKDOWN.columns)[:5]] == [
        _float32(v) for v in (0.2817, 0.001293, 0.000053, 0.0, -0.000413)
    ]
    assert held['alphagenome_splicing'] == _float32(0.000015)
    # A score of 0 is a score. The archives of the breakdown and of the splicing scores lack this allele.
    partly = rows['chr1-9-G-A']
    assert (partly['avi_raw'], partly['avi_phred']) == (0.0, 0.0)
    assert {partly[column] for source in (_BREAKDOWN, _SPLICING) for column in source.columns} == {None}
    # The archives hold single-base changes of chr1 to chrY: an insertion and an allele on chrM have no score.
    for lacking in ('chr1-9-G-GA', 'chrM-7-G-A'):
        assert {rows[lacking][column] for source in _SCORES for column in source.columns} == {None}
        assert rows[lacking]['gnomad_joint_ac'] == 7


RUNS = {
    'no archive of AlphaGenome': (),
    'the scores alone': (_AVI,),
    'the scores and their breakdown': (_AVI, _BREAKDOWN),
    'the splicing scores alone': (_SPLICING,),
    'all three': _SCORES,
}


@pytest.mark.parametrize('case', RUNS)
def test_a_run_has_the_columns_and_the_records_of_the_archives_it_was_given(tmp_path: pathlib.Path, case: str) -> None:
    sources = RUNS[case]
    info = _convert_scores(tmp_path, sources)
    assert pq.read_schema(tmp_path / 'out' / echtvar_values.TABLE_NAME).equals(
        echtvar_values.schema([_GNOMAD, *sources])
    )
    assert info.given() == (_GNOMAD, *sources)
    assert {name for name, record in info.records().items() if record is not None} == {
        _GNOMAD.name,
        *(source.name for source in sources),
    }
    assert echtvar_values.read_info(tmp_path / 'out') == info
    for source in sources:
        record = info.records()[source.name]
        assert record is not None
        archive = tmp_path / f'{source.name}.zip'
        assert (record.archive, record.bytes) == (archive.name, archive.stat().st_size)
        assert record.sha256 == hashlib.sha256(archive.read_bytes()).hexdigest()
        assert [field.alias for field in record.fields] == _aliases(source)
        # The notice is recorded for the archives whose terms ask for it, and the score's record has no such key.
        assert getattr(record, 'terms', None) == source.notice


def _changed(source: echtvar_values.Source, index: int, **changes: object) -> dict[str, object]:
    config = _embedded(source)
    config[index] |= changes
    return {source.name: config}


SCORE_ARCHIVES_REFUSED = {
    'a stand-in that a score can be': (
        {'configs': _changed(_AVI, 0, missing_value=-1)},
        'an allele it lacks gets the value -1 as its avi_raw, which a value of the field can be',
    ),
    'a score encoded as an integer': (
        {'configs': _changed(_AVI, 1, ftype='Integer')},
        'PHRED is encoded as Integer, not Float',
    ),
    'a breakdown without one of its contributions': (
        {'configs': {_BREAKDOWN.name: _embedded(_BREAKDOWN)[1:]}, 'declared': ()},
        "0 fields encode AlphaGenome's MERGED_SPLICING, and avi_shap_splicing needs one",
    ),
    'a breakdown whose field lacks the notice': (
        {'configs': _changed(_BREAKDOWN, 3, description='MAX_ABS_DNASE, invented')},
        'the description of MAX_ABS_DNASE does not end with the notice',
    ),
    'splicing scores without a description': (
        {'configs': {_SPLICING.name: [{k: v for k, v in _embedded(_SPLICING)[0].items() if k != 'description'}]}},
        'the description of alphagenome_splicing does not end with the notice',
    ),
    'two archives that add a field of the same alias': (
        {'configs': _changed(_SPLICING, 0, alias='avi_phred')},
        r'alphagenome_splicing\.zip and .*alphagenome_avi\.zip both add the INFO field avi_phred',
    ),
    "an archive that adds a field of gnomAD's alias beside its own": (
        {'configs': {_AVI.name: [*_embedded(_AVI), _embedded(_AVI)[0] | {'field': 'other', 'alias': 'gnomad_AF'}]}},
        r'alphagenome_avi\.zip and .*gnomad\.zip both add the INFO field gnomad_AF',
    ),
    'a raw score without its PHRED': (
        {'scored': [('chr1', 5, 'C', 'T', {_AVI.name: ('1.5', '.')})]},
        r"avi_phred says the archive lacks the allele, yet it has \['avi_raw'\]",
    ),
    'a breakdown without one contribution': (
        {'scored': [('chr1', 5, 'C', 'T', {_BREAKDOWN.name: ('.', *_SIXTEEN[1:])})]},
        'avi_shap_merged_splicing says the archive lacks the allele, yet it has',
    ),
    'a union that echtvar did not annotate with the archive': (
        {'declared': (_AVI, _SPLICING)},
        r'no INFO field avi_shap_merged_splicing, which alphagenome_avi_breakdown\.zip adds',
    ),
}


@pytest.mark.parametrize('case', SCORE_ARCHIVES_REFUSED)
def test_an_archive_or_a_union_that_is_not_what_a_source_describes_is_refused(
    tmp_path: pathlib.Path, case: str
) -> None:
    changes, message = SCORE_ARCHIVES_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _convert_scores(tmp_path, **changes)


def test_archives_other_than_a_runs_are_refused(tmp_path: pathlib.Path) -> None:
    annotated = _write_vcf(tmp_path / 'annotated.vcf', _RECORDS, _INFO)
    archive = _write_archive(tmp_path / 'invented.zip')
    with pytest.raises(ValueError, match=r"\['cadd'\] are no archives that the tables read"):
        echtvar_values.convert(annotated, annotated, {_GNOMAD.name: archive, 'cadd': archive}, '0', tmp_path / 'out')
    with pytest.raises(ValueError, match='no archive of gnomad_joint, which every run reads'):
        echtvar_values.convert(annotated, annotated, {_AVI.name: archive}, '0', tmp_path / 'out')
