"""gnomAD's joint values of a union, from an invented union as echtvar would annotate it with an invented archive."""

from __future__ import annotations

import copy
import hashlib
import json
import pathlib
import struct
import zipfile
from collections.abc import Sequence
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from case_annotation import cli, gnomad_joint

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
) -> gnomad_joint.GnomadJointInfo:
    sites = [(c, p, r, a, '.') for c, p, r, a, _ in (annotated if union is None else union)]
    return gnomad_joint.convert(
        _write_vcf(tmp_path / 'annotated.vcf', annotated, info, samples),
        _write_vcf(tmp_path / 'union.vcf', sites),
        _write_archive(tmp_path / 'invented.zip', config),
        'echtvar 0.0.0',
        tmp_path / 'out',
    )


def _rows(tmp_path: pathlib.Path) -> dict[str, dict[str, Any]]:
    table = pq.read_table(tmp_path / 'out' / gnomad_joint.TABLE_NAME)
    assert table.schema.equals(gnomad_joint.TABLE)
    return {row.pop('variant_id'): row for row in table.to_pylist()}


def _float32(value: float) -> float:
    """The 32-bit float nearest to the value, as Python reads it from a 32-bit column."""
    return struct.unpack('f', struct.pack('f', value))[0]


def test_each_allele_has_the_values_the_archive_stores(tmp_path: pathlib.Path) -> None:
    _convert(tmp_path)
    rows = _rows(tmp_path)
    assert list(rows) == ['chr1-5-C-T', 'chr1-9-G-A', 'chr2-100-A-G', 'chr2-200-CTGATGA-C', 'chrM-7-G-A']
    none = dict.fromkeys(gnomad_joint.COLUMNS)
    # An allele the archive lacks has no values; one it holds with a count of 0 keeps its 0s.
    assert rows['chr1-5-C-T'] == none
    assert rows['chrM-7-G-A'] == none
    assert rows['chr1-9-G-A'] == dict.fromkeys(gnomad_joint.COLUMNS, 0)
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
    frequencies = {f.name: f.type for f in gnomad_joint.TABLE if pa.types.is_floating(f.type)}
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
    assert gnomad_joint.read_info(tmp_path / 'out') == info


def test_the_command_writes_the_table_and_the_record(tmp_path: pathlib.Path) -> None:
    annotated = _write_vcf(tmp_path / 'annotated.vcf', _RECORDS, _INFO)
    union = _write_vcf(tmp_path / 'union.vcf', [(c, p, r, a, '.') for c, p, r, a, _ in _RECORDS])
    archive = _write_archive(tmp_path / 'invented.zip')
    args = ['--annotated', str(annotated), '--sites', str(union), '--archive', str(archive)]
    cli.main(['gnomad-joint', *args, '--echtvar-version', 'echtvar 0.0.0', '--output', str(tmp_path / 'out')])
    assert len(_rows(tmp_path)) == len(_RECORDS)
    assert gnomad_joint.read_info(tmp_path / 'out').echtvar_version == 'echtvar 0.0.0'


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
    'a stand-in count that a count can be': (_config(0, missing_value=0), 'which a count can be'),
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
        gnomad_joint.convert(annotated, annotated, not_a_zip, 'echtvar 0.0.0', tmp_path / 'out')
    without_config = tmp_path / 'other.zip'
    with zipfile.ZipFile(without_config, 'w') as z:
        z.writestr('readme.txt', 'invented')
    with pytest.raises(ValueError, match=r'other\.zip: no echtvar/config\.json'):
        gnomad_joint.convert(annotated, annotated, without_config, 'echtvar 0.0.0', tmp_path / 'out')


def test_a_conversion_without_echtvars_version_is_refused(tmp_path: pathlib.Path) -> None:
    annotated = _write_vcf(tmp_path / 'annotated.vcf', _RECORDS, _INFO)
    union = _write_vcf(tmp_path / 'union.vcf', [(c, p, r, a, '.') for c, p, r, a, _ in _RECORDS])
    archive = _write_archive(tmp_path / 'invented.zip')
    with pytest.raises(ValueError, match='echtvar_version'):
        gnomad_joint.convert(annotated, union, archive, '', tmp_path / 'out')


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
