"""The union of a group's sites of structural variants, from invented callers' records."""

from __future__ import annotations

import gzip
import pathlib
from collections.abc import Sequence
from typing import Any

import invented_sv_annotation
import pyarrow.parquet as pq
import pytest
import structural_vcfs

from case_annotation import cli, sv_sites

_DELETION = ('chr1', '1000', 'inv:DEL:1', 'T', '<DEL>', '.', '.', 'SVTYPE=DEL;END=5000;SVLEN=-4000')
_SEGMENT = ('chr10', '5000', 'inv:CNV:1', 'N', '<DEL>,<DUP>', '.', '.', 'SVTYPE=CNV;SVLEN=-10000,10000;END=15000')


def _union(tmp_path: pathlib.Path, families: dict[str, Sequence[Sequence[Sequence[str]]]]) -> sv_sites.Union:
    sites = {
        family: [
            invented_sv_annotation.write_sites(tmp_path / f'{family}.{index}.vcf.gz', records)
            for index, records in enumerate(callers)
        ]
        for family, callers in families.items()
    }
    return sv_sites.union(sites, tmp_path / 'union.vcf', tmp_path / 'map.parquet')


def _records(tmp_path: pathlib.Path) -> list[list[str]]:
    return [line.split('\t') for line in (tmp_path / 'union.vcf').read_text().splitlines() if not line.startswith('#')]


def _map(tmp_path: pathlib.Path) -> list[dict[str, Any]]:
    return pq.read_table(tmp_path / 'map.parquet').to_pylist()


def test_the_union_holds_each_allele_of_each_record_once(tmp_path: pathlib.Path) -> None:
    union = _union(tmp_path, {'fam-1': invented_sv_annotation.CALLERS})
    records = _records(tmp_path)
    kept = [
        r for r in (*structural_vcfs.SV_RECORDS, *structural_vcfs.CNV_RECORDS) if r[0].startswith('chr') and r[4] != '.'
    ]
    kept = [r for r in kept if r[0] != 'chr1_KI270706v1_random']
    alleles = [(r[2], alt) for r in kept for alt in r[4].split(',')]
    assert [(row['sv_id'], row['alt']) for row in _map(tmp_path)] == alleles
    assert (union.sites, union.alleles) == (len(alleles), len(alleles))
    assert len(records) == len(alleles)
    assert {r[4] for r in records} >= {'<DEL>', '<DUP>', '<DUP:TANDEM>', 'GA', '[chr2:500[T'}
    assert all(',' not in r[4] for r in records)
    # Sites-only: no samples, QUAL or FILTER, and INFO of the kept fields alone.
    assert {(len(r), r[5], r[6]) for r in records} == {(8, '.', '.')}
    assert {item.partition('=')[0] for r in records for item in r[7].split(';')} <= set(sv_sites.KEPT)
    # Every record but a breakend keeps the END its caller wrote.
    by_site = {r[2]: r[7] for r in records}
    info = {row['sv_id']: by_site[row['site_id']] for row in _map(tmp_path)}
    assert 'END=5000;' in info['inv:DEL:1']
    assert 'END=7;' in info['inv:INS:1']
    assert 'END=303;' in info['inv:DEL:2']
    assert 'END' not in info['inv:BND:1:0']
    assert union.with_end == sum('END=' in r[7] for r in records) == len(records) - 2


def test_a_record_of_two_alleles_becomes_a_site_of_each_with_its_own_length(tmp_path: pathlib.Path) -> None:
    _union(tmp_path, {'fam-1': [[_SEGMENT]]})
    rows = _map(tmp_path)
    assert [(row['sv_id'], row['alt']) for row in rows] == [('inv:CNV:1', '<DEL>'), ('inv:CNV:1', '<DUP>')]
    by_id = {r[2]: r for r in _records(tmp_path)}
    loss, gain = (by_id[row['site_id']] for row in rows)
    assert (loss[1], loss[4], loss[7]) == ('5000', '<DEL>', 'END=15000;SVTYPE=CNV;SVLEN=-10000')
    assert (gain[1], gain[4], gain[7]) == ('5000', '<DUP>', 'END=15000;SVTYPE=CNV;SVLEN=10000')


def test_a_site_has_the_id_of_its_content_whoever_holds_it(tmp_path: pathlib.Path) -> None:
    # Two families hold the deletion under different IDs and with the INFO fields in another order; the loss of one
    # family's two-allele segment is another family's segment of one allele; a longer deletion is another site.
    reordered = (*_DELETION[:2], 'inv:other:7', *_DELETION[3:7], 'SVLEN=-4000;END=5000;SVTYPE=DEL')
    longer = (*_DELETION[:2], 'inv:DEL:2', *_DELETION[3:7], 'SVTYPE=DEL;END=5001;SVLEN=-4001')
    loss = (*_SEGMENT[:2], 'inv:LOSS:9', 'N', '<DEL>', '.', '.', 'SVTYPE=CNV;SVLEN=-10000;END=15000')
    union = _union(tmp_path, {'fam-1': [[_DELETION, longer], [_SEGMENT]], 'fam-2': [[reordered, loss]]})
    sites = {(row['family'], row['sv_id'], row['alt']): row['site_id'] for row in _map(tmp_path)}
    assert sites[('fam-1', 'inv:DEL:1', '<DEL>')] == sites[('fam-2', 'inv:other:7', '<DEL>')]
    assert sites[('fam-1', 'inv:CNV:1', '<DEL>')] == sites[('fam-2', 'inv:LOSS:9', '<DEL>')]
    assert len(set(sites.values())) == 4
    assert (union.sites, union.alleles) == (4, 6)
    # The id says nothing of a family or of its record's ID, and the record under it has that content.
    found = sv_sites.read_union(_compressed(tmp_path))
    assert set(found) == set(sites.values())
    assert all(site_id == site.site_id and 'inv' not in site_id for site_id, site in found.items())


def _compressed(tmp_path: pathlib.Path) -> pathlib.Path:
    with gzip.open(tmp_path / 'union.vcf.gz', 'wt') as compressed:
        compressed.write((tmp_path / 'union.vcf').read_text())
    return tmp_path / 'union.vcf.gz'


def test_the_union_is_in_chromosome_order_and_declares_what_it_holds(tmp_path: pathlib.Path) -> None:
    _union(tmp_path, {'fam-1': invented_sv_annotation.CALLERS})
    records = _records(tmp_path)
    order = {contig: rank for rank, contig in enumerate(sv_sites.PRIMARY)}
    places = [(order[r[0]], int(r[1])) for r in records]
    assert places == sorted(places)
    assert {r[0] for r in records} == {'chr1', 'chr2', 'chr10', 'chrX', 'chrY'}
    header = [line for line in (tmp_path / 'union.vcf').read_text().splitlines() if line.startswith('#')]
    assert [line for line in header if line.startswith('##contig')] == [
        f'##contig=<ID={contig},length=200000000>' for contig in ('chr1', 'chr2', 'chr10', 'chrX', 'chrY', 'chrM')
    ]
    assert {line.split('ID=')[1].split(',')[0] for line in header if line.startswith('##INFO')} == set(sv_sites.KEPT)
    assert {line.split('ID=')[1].split(',')[0] for line in header if line.startswith('##ALT')} == {
        'DEL',
        'DUP',
        'DUP:TANDEM',
    }
    assert header[-1].split('\t') == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']


def _line(*columns: str) -> list[str]:
    return list(columns)


_OK = ('chr1', '1000', 'inv:DEL:1', 'T', '<DEL>', '.', '.')
RECORDS_REFUSED: dict[str, tuple[list[str], str]] = {
    'samples': ([*_OK, 'SVTYPE=DEL;END=5000', 'GT', '0/1'], 'the eight without samples'),
    'a contig off the primary chromosomes': (['chrUn', *_OK[1:], 'SVTYPE=DEL;END=5000'], 'not on a primary chromosome'),
    'a position that is no number': (['chr1', 'x', *_OK[2:], 'SVTYPE=DEL;END=5000'], 'is not a number'),
    'a QUAL': ([*_OK[:5], '50', '.', 'SVTYPE=DEL;END=5000'], 'have neither'),
    'a FILTER': ([*_OK[:5], '.', 'PASS', 'SVTYPE=DEL;END=5000'], 'have neither'),
    'an INFO field that the sites do not hold': ([*_OK, 'SVTYPE=DEL;END=5000;MATEID=inv:x'], 'INFO MATEID'),
    'an INFO field twice': ([*_OK, 'SVTYPE=DEL;END=5000;END=5001'], 'INFO END twice'),
    'no SVTYPE': ([*_OK, 'END=5000'], 'no SVTYPE'),
    'an SVTYPE without a value': ([*_OK, 'SVTYPE;END=5000'], 'no SVTYPE'),
    'an empty SVTYPE': ([*_OK, 'SVTYPE=;END=5000'], "SVTYPE '' is none of DEL"),
    'an SVTYPE that the Annotator is not known to annotate': ([*_OK, 'SVTYPE=CPX;END=5000'], "SVTYPE 'CPX' is none of"),
    'a copy-number segment of another allele': (
        [*_OK[:4], '<CNV>', '.', '.', 'SVTYPE=CNV;END=5000'],
        'the alleles of a copy-number segment are <DEL> and <DUP>, not <CNV>',
    ),
    'an END before the position': ([*_OK, 'SVTYPE=DEL;END=999'], "is not a position at or after the record's, 1000"),
    'IMPRECISE with a value': ([*_OK, 'SVTYPE=DEL;END=5000;IMPRECISE=1'], 'IMPRECISE is a flag'),
    'no alternate allele': ([*_OK[:4], '.', '.', '.', 'SVTYPE=DEL;END=5000'], 'are not one or several'),
    'an allele twice': ([*_OK[:4], '<DEL>,<DEL>', '.', '.', 'SVTYPE=CNV;END=5000'], 'are not one or several'),
    'one length for two alleles': (
        [*_OK[:4], '<DEL>,<DUP>', '.', '.', 'SVTYPE=CNV;END=5000;SVLEN=-4000'],
        'gives 1 lengths for 2 alleles',
    ),
    'a record without END that is no breakend': ([*_OK, 'SVTYPE=DEL;SVLEN=-4000'], 'no END, which every record but'),
    'a sequence without END': (
        [*_OK[:3], 'ACGT', 'A', '.', '.', 'SVTYPE=DEL;SVLEN=-3'],
        'no END, which every record but',
    ),
    'a breakend with END': ([*_OK[:4], 'T[chr2:500[', '.', '.', 'SVTYPE=BND;END=1000'], 'a breakend has no end'),
    'an END that is no position': ([*_OK, 'SVTYPE=DEL;END=later'], 'is not a position at or after'),
}


@pytest.mark.parametrize('case', RECORDS_REFUSED)
def test_a_record_outside_the_layout_of_a_familys_sites_is_refused(case: str) -> None:
    columns, message = RECORDS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        sv_sites.allele_sites(columns, 'sites.vcf: line 9')


def _write(
    path: pathlib.Path, lines: Sequence[str], contigs: str = '##contig=<ID=chr1,length=1000000>'
) -> pathlib.Path:
    with gzip.open(path, 'wt') as vcf:
        vcf.write(
            '\n'.join(['##fileformat=VCFv4.2', contigs, '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', *lines]) + '\n'
        )
    return path


_SITE = '\t'.join([*_OK, 'SVTYPE=DEL;END=5000'])


def test_a_record_without_an_id_is_refused(tmp_path: pathlib.Path) -> None:
    path = _write(tmp_path / 'sites.vcf.gz', [_SITE.replace('inv:DEL:1', '.')])
    with pytest.raises(ValueError, match='a record without an ID'):
        sv_sites.read_family(path)


FAMILIES_REFUSED: dict[str, tuple[list[list[str]], list[str], str]] = {
    "an ID twice in one caller's sites": ([[_SITE, _SITE.replace('1000', '2000')]], [], 'is that of several records'),
    "an ID in both callers' sites": ([[_SITE], [_SITE.replace('<DEL>', '<DUP>')]], [], 'is that of several records'),
    'a contig of two lengths': (
        [[_SITE], [_SITE.replace('inv:DEL:1', 'inv:DEL:2')]],
        ['##contig=<ID=chr1,length=1000000>', '##contig=<ID=chr1,length=999>'],
        'chr1 is 999 bases long, and 1000000 in another file',
    ),
    'a record on a contig that no header declares': (
        [[_SITE.replace('chr1\t', 'chr2\t')]],
        [],
        r"declare no contig \['chr2'\]",
    ),
}


@pytest.mark.parametrize('case', FAMILIES_REFUSED)
def test_sites_that_do_not_fit_together_are_refused(tmp_path: pathlib.Path, case: str) -> None:
    files, contigs, message = FAMILIES_REFUSED[case]
    paths = [
        _write(tmp_path / f'{index}.vcf.gz', lines, *([contigs[index]] if contigs else []))
        for index, lines in enumerate(files)
    ]
    with pytest.raises(ValueError, match=message):
        sv_sites.union({'fam-1': paths}, tmp_path / 'union.vcf', tmp_path / 'map.parquet')


def test_two_families_may_use_one_id_for_different_records(tmp_path: pathlib.Path) -> None:
    other = (*_DELETION[:1], '9000', *_DELETION[2:7], 'SVTYPE=DEL;END=9500')
    union = _union(tmp_path, {'fam-1': [[_DELETION]], 'fam-2': [[other]]})
    assert (union.sites, union.alleles) == (2, 2)
    assert len({row['site_id'] for row in _map(tmp_path)}) == 2


def _site_id(line: str) -> str:
    return invented_sv_annotation.site(line).site_id


UNIONS_REFUSED: dict[str, tuple[list[str], str]] = {
    "an ID that is not the id of the record's content": ([_SITE], "is not the id of the record's content"),
    'a site twice': ([_SITE.replace('inv:DEL:1', _site_id(_SITE))] * 2, 'or is there twice'),
    'a record of two alleles': (
        ['\t'.join([*_OK[:4], '<DEL>,<DUP>', '.', '.', 'SVTYPE=CNV;END=5000'])],
        'the union holds one for each record',
    ),
    'a record that lost its END': ([_SITE.replace(';END=5000', '')], 'no END, which every record but'),
}


@pytest.mark.parametrize('case', UNIONS_REFUSED)
def test_a_vcf_that_is_no_union_of_sites_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    lines, message = UNIONS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        sv_sites.read_union(_write(tmp_path / 'union.vcf.gz', lines))


def test_the_command_writes_the_union_and_its_map_and_refuses_an_empty_one(tmp_path: pathlib.Path) -> None:
    sites = invented_sv_annotation.write_sites(tmp_path / 'sv.vcf.gz', structural_vcfs.SV_RECORDS)
    segments = invented_sv_annotation.write_sites(tmp_path / 'cnv.vcf.gz', structural_vcfs.CNV_RECORDS)
    outputs = ['--output-vcf', str(tmp_path / 'union.vcf'), '--output-map', str(tmp_path / 'map.parquet')]
    cli.main(['sv-union', '--sites', 'fam-1', str(sites), '--sites', 'fam-1', str(segments), *outputs])
    assert {row['family'] for row in _map(tmp_path)} == {'fam-1'}
    assert len(_records(tmp_path)) == len(_map(tmp_path))
    empty = invented_sv_annotation.write_sites(tmp_path / 'empty.vcf.gz', [])
    with pytest.raises(SystemExit, match='give no record on the primary chromosomes'):
        cli.main(['sv-union', '--sites', 'fam-1', str(empty), *outputs])
