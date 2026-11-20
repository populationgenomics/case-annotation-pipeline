"""Runs of homozygosity and regions of uniparental disomy: what the tools are given, and their output as rows."""

from __future__ import annotations

import gzip
import importlib.metadata
import json
import pathlib
from collections.abc import Sequence

import invented_roh_upd
import pyarrow.parquet as pq
import pysam
import pytest

from case_annotation import cli, pedigree, roh_upd

_FAMILY = 'fam-1'
# A father, a mother and their two children, and an aunt, each by a PED line.
_LINES = {
    'S1': 'fam-1\tS1\t0\t0\t1\t1',
    'S2': 'fam-1\tS2\t0\t0\t2\t1',
    'S3': 'fam-1\tS3\tS1\tS2\t2\t2',
    'S4': 'fam-1\tS4\tS1\tS2\t1\t1',
    'S5': 'fam-1\tS5\t0\t0\t2\t1',
}


def _ped(tmp_path: pathlib.Path, members: Sequence[str] = tuple(_LINES)) -> pathlib.Path:
    path = tmp_path / 'family.ped'
    path.write_text(''.join(f'{_LINES[member]}\n' for member in members))
    return path


def _vcf(tmp_path: pathlib.Path, samples: Sequence[str], formats: Sequence[str] = ('GT', 'GQ')) -> pathlib.Path:
    """A VCF with these samples and FORMAT fields, and no record: the header is all that is read of it."""
    header = pysam.VariantHeader()
    header.contigs.add('chr1', length=1_000_000)
    for key in formats:
        kind = 'String' if key == 'GT' else 'Integer'
        header.add_meta('FORMAT', items=[('ID', key), ('Number', '1'), ('Type', kind), ('Description', 'invented')])
    for sample in samples:
        header.add_sample(sample)
    path = tmp_path / 'family.vcf.gz'
    with pysam.VariantFile(str(path), 'wz', header=header):
        pass
    return path


FAMILIES = {
    'a trio': (('S1', 'S2', 'S3'), [('S3', 'S2', 'S1')]),
    'two children, in the order of the samples': (('S4', 'S1', 'S2', 'S3'), [('S4', 'S2', 'S1'), ('S3', 'S2', 'S1')]),
    'a child and its mother': (('S2', 'S3'), []),
    'a child alone': (('S3',), []),
    'two parents and a relative without a child': (('S1', 'S2', 'S5'), []),
}


@pytest.mark.parametrize('case', FAMILIES)
def test_upd_assesses_each_child_whose_parents_are_both_sequenced(tmp_path: pathlib.Path, case: str) -> None:
    samples, expected = FAMILIES[case]
    found = roh_upd.trios(pedigree.read(_ped(tmp_path)), samples)
    assert [(trio.child, trio.mother, trio.father) for trio in found] == expected
    planned = roh_upd.plan(
        _FAMILY, _ped(tmp_path), _vcf(tmp_path, samples), invented_roh_upd.archive(tmp_path / 'a.zip')
    )
    assert planned.trios == found


def test_the_plan_gives_the_archives_name_of_the_frequency(tmp_path: pathlib.Path) -> None:
    archive = invented_roh_upd.archive(tmp_path / 'a.zip', alias='joint_af')
    planned = roh_upd.plan(_FAMILY, _ped(tmp_path), _vcf(tmp_path, ('S1', 'S2', 'S3')), archive)
    assert (planned.frequency.field, planned.frequency.alias) == ('AF_joint', 'joint_af')


ARCHIVES_REFUSED = {
    'a stand-in that both tools would read as a frequency': ({'missing_value': 1}, 'only where that value is 0'),
    'a stand-in that is written as a missing value': ({'missing_value': -1}, 'only where that value is 0'),
    'an alias that is no INFO key': ({'alias': 'AF joint'}, 'is no name of an INFO field'),
    'an alias that would end a command line': ({'alias': 'AF;x'}, 'is no name of an INFO field'),
    'an archive without the frequency': ({'field': 'AF_genomes'}, "0 fields encode gnomAD's AF_joint"),
}


@pytest.mark.parametrize('case', ARCHIVES_REFUSED)
def test_an_archive_that_the_tools_cannot_read_a_frequency_from_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    changed, message = ARCHIVES_REFUSED[case]
    archive = invented_roh_upd.archive(tmp_path / 'a.zip', **changed)
    with pytest.raises(ValueError, match=message):
        roh_upd.plan(_FAMILY, _ped(tmp_path), _vcf(tmp_path, ('S1', 'S2', 'S3')), archive)


def test_a_family_with_an_assessed_child_and_no_gq_is_refused(tmp_path: pathlib.Path) -> None:
    archive = invented_roh_upd.archive(tmp_path / 'a.zip')
    without = _vcf(tmp_path, ('S1', 'S2', 'S3'), formats=('GT',))
    with pytest.raises(ValueError, match=r'family\.vcf\.gz: no FORMAT field GQ, which upd needs of S3'):
        roh_upd.plan(_FAMILY, _ped(tmp_path), without, archive)
    # bcftools roh reads the genotypes alone, so a family that upd does not assess needs no GQ.
    duo = _vcf(tmp_path, ('S2', 'S3'), formats=('GT',))
    assert roh_upd.plan(_FAMILY, _ped(tmp_path), duo, archive).trios == ()


def test_the_plan_refuses_a_ped_file_that_lacks_a_sample(tmp_path: pathlib.Path) -> None:
    archive = invented_roh_upd.archive(tmp_path / 'a.zip')
    with pytest.raises(ValueError, match=r"lacks the VCF samples \['S5'\]"):
        roh_upd.plan(_FAMILY, _ped(tmp_path, ('S1', 'S2', 'S3')), _vcf(tmp_path, ('S1', 'S2', 'S3', 'S5')), archive)
    with pytest.raises(ValueError, match="is of family 'fam-1', not 'fam-2'"):
        roh_upd.plan('fam-2', _ped(tmp_path), _vcf(tmp_path, ('S1', 'S2', 'S3')), archive)


def test_the_plan_command_prints_the_frequencys_field_and_lists_the_children(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    arguments = ['roh-upd-plan', '--family', _FAMILY, '--ped', str(_ped(tmp_path))]
    arguments += ['--vcf', str(_vcf(tmp_path, ('S4', 'S1', 'S2', 'S3')))]
    arguments += ['--archive', str(invented_roh_upd.archive(tmp_path / 'a.zip')), '--trios', str(tmp_path / 'trios')]
    cli.main(arguments)
    assert capsys.readouterr().out == 'gnomad_AF\n'
    assert (tmp_path / 'trios').read_text() == 'S4\tS2\tS1\nS3\tS2\tS1\n'


def _roh(tmp_path: pathlib.Path, *lines: str) -> pathlib.Path:
    path = tmp_path / 'roh.txt'
    path.write_text(invented_roh_upd.HEADER + ''.join(f'{line}\n' for line in lines))
    return path


def test_a_run_is_a_row_with_the_first_and_the_last_position_of_the_run(tmp_path: pathlib.Path) -> None:
    lines = ('RG\tS3\tchr2\t1001\t3000\t2000\t41\t95.3', invented_roh_upd.run('S1', 'chr22', 7, 7, markers=1))
    rows = roh_upd.read_roh(_roh(tmp_path, *lines), ('S1', 'S2', 'S3'))
    assert rows[0] == {
        'member': 'S3',
        'chrom': 'chr2',
        'pos': 1001,
        'end_pos': 3000,
        'length': 2000,
        'markers': 41,
        'quality': 95.3,
    }
    # The length counts both ends, as the positions do.
    assert all(row['length'] == row['end_pos'] - row['pos'] + 1 for row in rows)
    assert (rows[1]['member'], rows[1]['length']) == ('S1', 1)
    assert roh_upd.read_roh(_roh(tmp_path), ('S1',)) == []


RUNS_REFUSED = {
    'a line of another kind': ('ST\tS3\tchr2\t1001\t0\t99.0', 'not a line of bcftools roh with the runs alone'),
    'a line without its quality': ('RG\tS3\tchr2\t1001\t3000\t2000\t41', 'not a line of bcftools roh'),
    'a sample that the VCF lacks': (invented_roh_upd.run('S9', 'chr2', 1, 9), "a run of 'S9'"),
    'a run on X': (invented_roh_upd.run('S3', 'chrX', 1, 9), "a run on 'chrX'"),
    'a run on a contig without its prefix': (invented_roh_upd.run('S3', '2', 1, 9), "a run on '2'"),
    'a start that is no integer': ('RG\tS3\tchr2\t1e3\t3000\t2000\t41\t95.3', "the start is '1e3'"),
    'a length that leaves an end out': ('RG\tS3\tchr2\t1001\t3000\t1999\t41\t95.3', 'is given the length 1999'),
    'a start before the first base': ('RG\tS3\tchr2\t0\t9\t10\t41\t95.3', 'is given the length 10'),
    'a quality that is no number': ('RG\tS3\tchr2\t1001\t3000\t2000\t41\thigh', "the quality is 'high'"),
}


@pytest.mark.parametrize('case', RUNS_REFUSED)
def test_output_that_is_not_bcftools_rohs_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    line, message = RUNS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        roh_upd.read_roh(_roh(tmp_path, invented_roh_upd.run('S3', 'chr1', 5, 90), line), ('S1', 'S2', 'S3'))


def _bed(tmp_path: pathlib.Path, *lines: str) -> pathlib.Path:
    path = tmp_path / 'regions.bed'
    path.write_text(''.join(f'{line}\n' for line in lines))
    return path


def test_a_region_is_a_row_with_both_pairs_of_bounds_and_without_upds_type(tmp_path: pathlib.Path) -> None:
    line = (
        'chr15\t1999\t9000\tORIGIN=PATERNAL;TYPE=ISODISOMY/DELETION;LOW_SIZE=7000;INF_SITES=31;SNPS=140;'
        'HET_HOM=0/96;OPP_SITES=2;START_LOW=1200;END_HIGH=9400;HIGH_SIZE=8200'
    )
    assert roh_upd.read_upd(_bed(tmp_path, line), 'S3') == [
        {
            'member': 'S3',
            'chrom': 'chr15',
            # A BED start is one below the first position.
            'pos': 2000,
            'end_pos': 9000,
            'outer_pos': 1200,
            'outer_end_pos': 9400,
            'origin': 'paternal',
            'informative_sites': 31,
            'snps': 140,
            'opposite_sites': 2,
        }
    ]
    assert list(roh_upd.read_upd(_bed(tmp_path, line), 'S3')[0]) == roh_upd.UPD.names
    assert roh_upd.read_upd(_bed(tmp_path, invented_roh_upd.region('chr2', 500, 9000)), 'S3')[0]['origin'] == 'maternal'
    assert roh_upd.read_upd(_bed(tmp_path), 'S3') == []


_REGION = invented_roh_upd.region('chr2', 500, 9000)
REGIONS_REFUSED = {
    'a line of three columns': ('chr2\t499\t9000', 'not a line of upd regions'),
    'a region without its origin': (_REGION.replace('ORIGIN=MATERNAL;', ''), 'are not those of upd regions'),
    'a value that upd does not write': (_REGION + ';MOSAIC=0.5', 'are not those of upd regions'),
    'a key without a value': (_REGION.replace('OPP_SITES=0', 'OPP_SITES'), 'are not those of upd regions'),
    'an origin of neither parent': (invented_roh_upd.region('chr2', 500, 9000, origin='BOTH'), "the origin is 'BOTH'"),
    'a region on X': (invented_roh_upd.region('chrX', 500, 9000), "a region on 'chrX'"),
    'a count that is no integer': (invented_roh_upd.region('chr2', 500, 9000, SNPS='many'), "SNPS is 'many'"),
    'an inner size that is not the distance of the inner bounds': (
        invented_roh_upd.region('chr2', 500, 9000, LOW_SIZE=8501),
        'are not the distances between the bounds',
    ),
    'an outer size that is not the distance of the outer bounds': (
        invented_roh_upd.region('chr2', 500, 9000, HIGH_SIZE=1),
        'are not the distances between the bounds',
    ),
}


@pytest.mark.parametrize('case', REGIONS_REFUSED)
def test_output_that_is_not_upds_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    line, message = REGIONS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        roh_upd.read_upd(_bed(tmp_path, _REGION, line), 'S3')


def test_the_folder_of_a_trio_holds_both_tables_and_the_record_of_what_ran(tmp_path: pathlib.Path) -> None:
    ped, vcf = _ped(tmp_path), _vcf(tmp_path, ('S1', 'S2', 'S3'))
    runs = (invented_roh_upd.run('S3', 'chr2', 400, 9500), invented_roh_upd.run('S1', 'chr1', 10, 80))
    folder = invented_roh_upd.folder(tmp_path / 'work', _FAMILY, ped, vcf, runs, {'S3': (_REGION,)}, sites=7)
    assert sorted(path.name for path in folder.iterdir()) == ['roh.parquet', 'roh_upd.json', 'upd.parquet']
    roh, upd = pq.read_table(folder / 'roh.parquet'), pq.read_table(folder / 'upd.parquet')
    assert (roh.schema, upd.schema) == (roh_upd.ROH, roh_upd.UPD)
    assert [(r['member'], r['chrom'], r['pos']) for r in roh.to_pylist()] == [('S3', 'chr2', 400), ('S1', 'chr1', 10)]
    assert [(r['member'], r['pos'], r['end_pos'], r['origin']) for r in upd.to_pylist()] == [
        ('S3', 500, 9000, 'maternal')
    ]
    info = roh_upd.read_info(folder)
    assert info == roh_upd.RohUpdInfo.model_validate_json((folder / 'roh_upd.json').read_text())
    assert (info.roh_members, info.sites) == (('S1', 'S2', 'S3'), 7)
    assert info.upd_children == (roh_upd.Trio(child='S3', mother='S2', father='S1'),)
    assert (info.upd_version, info.upd_commit) == invented_roh_upd.RAN['upd']
    assert (info.preparation, info.roh_command) == (
        invented_roh_upd.RAN['preparation'],
        invented_roh_upd.RAN['roh_command'],
    )
    archive = (tmp_path / 'work' / 'invented.zip').read_bytes()
    assert (info.frequency_archive.archive, info.frequency_archive.bytes) == ('invented.zip', len(archive))
    assert (info.frequency_archive.field.field, info.frequency_archive.field.alias) == ('AF_joint', 'gnomad_AF')


def test_an_assessed_child_without_a_region_has_a_table_without_rows(tmp_path: pathlib.Path) -> None:
    folder = invented_roh_upd.folder(tmp_path / 'work', _FAMILY, _ped(tmp_path), _vcf(tmp_path, ('S1', 'S2', 'S3')))
    assert pq.read_table(folder / 'upd.parquet').num_rows == 0
    assert pq.read_table(folder / 'roh.parquet').num_rows == 0
    assert [trio.child for trio in roh_upd.read_info(folder).upd_children] == ['S3']


def test_a_family_without_an_assessed_child_has_runs_and_no_table_of_regions(tmp_path: pathlib.Path) -> None:
    runs = (invented_roh_upd.run('S3', 'chr2', 400, 9500),)
    folder = invented_roh_upd.folder(tmp_path / 'work', _FAMILY, _ped(tmp_path), _vcf(tmp_path, ('S2', 'S3')), runs)
    assert sorted(path.name for path in folder.iterdir()) == ['roh.parquet', 'roh_upd.json']
    info = roh_upd.read_info(folder)
    assert (info.roh_members, info.upd_children) == (('S2', 'S3'), ())


CONVERSIONS_REFUSED = {
    'a number of sites below 0': ({'sites': -1}, '-1 is no number of records'),
    'regions of a member who is no child of two samples': (
        {'regions': {'S3': (), 'S1': ()}},
        r"regions of \['S1', 'S3'\] were given and \[\] as not assessed, and the children .* are \['S3'\]",
    ),
    'neither regions of the child nor the child as not assessed': ({'regions': {}}, r'regions of \[\] were given'),
    'a child with regions that is also given as not assessed': (
        {'regions': {'S3': ()}, 'not_assessed': ('S3',)},
        r"regions of \['S3'\] were given and \['S3'\] as not assessed",
    ),
    'a member who is no child of two samples as not assessed': (
        {'not_assessed': ('S3', 'S2')},
        r"\['S2', 'S3'\] as not assessed",
    ),
}


@pytest.mark.parametrize('case', CONVERSIONS_REFUSED)
def test_a_conversion_of_what_the_step_did_not_run_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    changed, message = CONVERSIONS_REFUSED[case]
    ped, vcf = _ped(tmp_path), _vcf(tmp_path, ('S1', 'S2', 'S3'))
    with pytest.raises(ValueError, match=message):
        invented_roh_upd.folder(tmp_path / 'work', _FAMILY, ped, vcf, **changed)
    assert not (tmp_path / 'work' / 'roh_upd').exists()


def _prepared(tmp_path: pathlib.Path, *records: tuple[str, tuple[str, ...]]) -> pathlib.Path:
    """A prepared VCF of three samples, as text: for each record its INFO and its three calls after `0/1`."""
    path = tmp_path / 'prepared.vcf.gz'
    with gzip.open(path, 'wt') as vcf:
        vcf.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n')
        for number, (info, calls) in enumerate(records, start=1):
            keys = 'GT:GQ' if all(calls) else 'GT'
            samples = '\t'.join(f'0/1:{call}' if call else '0/1' for call in calls)
            vcf.write(f'chr1\t{number * 100}\t.\tA\tC\t50\tPASS\t{info}\t{keys}\t{samples}\n')
    return path


_GOOD = ('gnomad_AF=0.1', ('60', '60', '60'))
_LAST = ('gnomad_AF=0', ('60', '60', '60'))
# Each case: the records of a prepared VCF, and whether upd keeps one of them at 0.05 and GQ 30.
SITES = {
    'a common site with high GQs, before another record': ((_GOOD, _LAST), True),
    'a site at both thresholds': ((('gnomad_AF=0.05', ('30', '30', '30')), _LAST), True),
    'such a site after sites that upd leaves out': ((_LAST, _LAST, _GOOD, _LAST), True),
    'a frequency in another notation': ((('gnomad_AF=5e-2;gnomad_AC=3', ('30', '31', '32')), _LAST), True),
    'no record': ((), False),
    'such a site as the only record, which upd never reads': ((_GOOD,), False),
    'such a site as the last record': ((_LAST, _GOOD), False),
    'a frequency below the threshold': ((('gnomad_AF=0.049', ('60', '60', '60')), _LAST), False),
    'an allele that the archive lacks': ((_LAST, _LAST), False),
    'a record without the frequency': ((('.', ('60', '60', '60')), _LAST), False),
    'a GQ below the threshold in the father': ((('gnomad_AF=0.1', ('29', '60', '60')), _LAST), False),
    'a GQ below the threshold in the child': ((('gnomad_AF=0.1', ('60', '60', '29')), _LAST), False),
    'a missing GQ': ((('gnomad_AF=0.1', ('60', '.', '60')), _LAST), False),
    'calls without GQ': ((('gnomad_AF=0.1', ('', '', '')), _LAST), False),
}


@pytest.mark.parametrize('case', SITES)
def test_upd_is_asked_for_only_where_it_would_keep_a_site(tmp_path: pathlib.Path, case: str) -> None:
    records, kept = SITES[case]
    assert roh_upd.keeps_a_site(_prepared(tmp_path, *records), 'gnomad_AF', 0.05, 30) is kept


def test_the_thresholds_are_the_ones_given(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path, ('gnomad_AF=0.01', ('20', '20', '20')), _LAST)
    assert not roh_upd.keeps_a_site(prepared, 'gnomad_AF', 0.05, 30)
    assert not roh_upd.keeps_a_site(prepared, 'gnomad_AF', 0.01, 30)
    assert roh_upd.keeps_a_site(prepared, 'gnomad_AF', 0.01, 20)
    # Another field's value is not the frequency.
    assert not roh_upd.keeps_a_site(prepared, 'gnomad_AC', 0.01, 20)


SITES_REFUSED = {
    'a frequency that is no number': (('gnomad_AF=common', ('60', '60', '60')), "gnomad_AF is 'common'"),
}


@pytest.mark.parametrize('case', SITES_REFUSED)
def test_a_prepared_vcf_that_upd_could_not_read_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    record, message = SITES_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        roh_upd.keeps_a_site(_prepared(tmp_path, record, _LAST), 'gnomad_AF', 0.05, 30)
    sites_only = tmp_path / 'sites.vcf.gz'
    with gzip.open(sites_only, 'wt') as vcf:
        vcf.write(
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC\t.\t.\t.\nchr1\t2\t.\tA\tC\t.\t.\t.\n'
        )
    with pytest.raises(ValueError, match='not a record of a VCF with samples'):
        roh_upd.keeps_a_site(sites_only, 'gnomad_AF', 0.05, 30)


def test_the_sites_command_prints_what_upd_would_keep(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    options = ['--af-tag', 'gnomad_AF', '--min-af', '0.05', '--min-gq', '30']
    cli.main(['roh-upd-sites', '--prepared', str(_prepared(tmp_path, _GOOD, _LAST)), *options])
    assert capsys.readouterr().out == 'kept\n'
    cli.main(['roh-upd-sites', '--prepared', str(_prepared(tmp_path, _GOOD)), *options])
    assert capsys.readouterr().out == 'none\n'


def test_a_child_that_upd_was_not_run_for_is_recorded_and_has_no_row(tmp_path: pathlib.Path) -> None:
    ped, vcf = _ped(tmp_path), _vcf(tmp_path, ('S4', 'S1', 'S2', 'S3'))
    # Without the list, the conversion refuses the family: no region was given for either child.
    with pytest.raises(ValueError, match=r'regions of \[\] were given and \[\] as not assessed'):
        invented_roh_upd.folder(tmp_path / 'refused', _FAMILY, ped, vcf, regions={})
    folder = invented_roh_upd.folder(tmp_path / 'work', _FAMILY, ped, vcf, sites=0, not_assessed=('S4', 'S3'))
    info = roh_upd.read_info(folder)
    assert (info.upd_children, info.sites) == ((), 0)
    assert [(entry.child, entry.mother, entry.father, entry.reason) for entry in info.upd_not_assessed] == [
        ('S4', 'S2', 'S1', 'no informative site'),
        ('S3', 'S2', 'S1', 'no informative site'),
    ]
    # The family has the table, which holds no row of a child that was not assessed.
    assert pq.read_table(folder / 'upd.parquet').num_rows == 0
    roh_upd.check(info, pedigree.read(ped), ('S4', 'S1', 'S2', 'S3'))
    # One child assessed and the other not, as the record can also say.
    mixed = invented_roh_upd.folder(
        tmp_path / 'mixed', _FAMILY, ped, vcf, regions={'S3': (_REGION,)}, not_assessed=('S4',)
    )
    info = roh_upd.read_info(mixed)
    assert [trio.child for trio in info.upd_children] == ['S3']
    assert [entry.child for entry in info.upd_not_assessed] == ['S4']
    assert {row['member'] for row in pq.read_table(mixed / 'upd.parquet').to_pylist()} == {'S3'}
    roh_upd.check(info, pedigree.read(ped), ('S4', 'S1', 'S2', 'S3'))
    with pytest.raises(ValueError, match='reason'):
        roh_upd.RohUpdInfo.model_validate(
            json.loads((mixed / 'roh_upd.json').read_text())
            | {'upd_not_assessed': [{'child': 'S4', 'mother': 'S2', 'father': 'S1', 'reason': 'too few reads'}]}
        )


def test_a_record_without_a_tools_version_is_refused(tmp_path: pathlib.Path) -> None:
    folder = invented_roh_upd.folder(tmp_path / 'work', _FAMILY, _ped(tmp_path), _vcf(tmp_path, ('S1', 'S2', 'S3')))
    record = json.loads((folder / 'roh_upd.json').read_text())
    changed = {
        'bcftools_version': '',
        'upd_version': '',
        'upd_commit': '0.1.1',
        'roh_command': '',
        'upd_command': '',
        'frequency_archive': record['frequency_archive'] | {'echtvar_version': ''},
    }
    for key, value in changed.items():
        with pytest.raises(ValueError, match=key):
            roh_upd.RohUpdInfo.model_validate(record | {key: value})
    # A commit is its 40 hexadecimal digits and nothing around them.
    with pytest.raises(ValueError, match='upd_commit'):
        roh_upd.RohUpdInfo.model_validate(record | {'upd_commit': 'commit ' + '0' * 40})


def test_a_folder_of_another_family_is_refused(tmp_path: pathlib.Path) -> None:
    ped = pedigree.read(_ped(tmp_path))
    folder = invented_roh_upd.folder(tmp_path / 'work', _FAMILY, _ped(tmp_path), _vcf(tmp_path, ('S1', 'S2', 'S3')))
    info = roh_upd.read_info(folder)
    roh_upd.check(info, ped, ('S1', 'S2', 'S3'))
    with pytest.raises(ValueError, match=r"runs of homozygosity are of \['S1', 'S2', 'S3'\]"):
        roh_upd.check(info, ped, ('S1', 'S2', 'S3', 'S4'))
    with pytest.raises(ValueError, match=r"runs of homozygosity are of \['S1', 'S2', 'S3'\]"):
        roh_upd.check(info, ped, ('S3', 'S2', 'S1'))
    # The same samples under a pedigree in which nobody is a child of two of them.
    unrelated = pedigree.read(_ped(tmp_path, ('S1', 'S2', 'S5')))
    renamed = info.model_copy(update={'roh_members': ('S1', 'S2', 'S5')})
    with pytest.raises(
        ValueError, match=r"upd assessed \['S3'\] and did not assess \[\], and the children .* are \[\]"
    ):
        roh_upd.check(renamed, unrelated, ('S1', 'S2', 'S5'))
    # A child that the record lists as assessed and as not assessed.
    twice = info.model_copy(
        update={
            'upd_not_assessed': (
                roh_upd.NotAssessed(**info.upd_children[0].model_dump(), reason='no informative site'),
            )
        }
    )
    with pytest.raises(ValueError, match=r"upd assessed \['S3'\] and did not assess \['S3'\]"):
        roh_upd.check(twice, ped, ('S1', 'S2', 'S3'))


class _Installed:
    """An installed distribution, as importlib.metadata gives one: its version and its record of where it came from."""

    version = '0.0'

    def __init__(self, origin: object) -> None:
        self._origin = origin

    def read_text(self, name: str) -> str | None:
        assert name == 'direct_url.json'
        return None if self._origin is None else json.dumps(self._origin)


def _install(monkeypatch: pytest.MonkeyPatch, origin: object) -> None:
    def distribution(name: str) -> _Installed:
        assert name == 'upd'
        return _Installed(origin)

    monkeypatch.setattr(importlib.metadata, 'distribution', distribution)


def test_the_installed_upd_is_recorded_by_its_version_and_its_commit(monkeypatch: pytest.MonkeyPatch) -> None:
    commit = 'ab' * 20
    _install(monkeypatch, {'url': roh_upd.UPD_REPOSITORY, 'vcs_info': {'vcs': 'git', 'commit_id': commit}})
    assert roh_upd.installed_upd() == ('0.0', commit)


UPD_REFUSED = {
    'a package from an index, as the other project of that name': None,
    'a commit of another repository': {
        'url': 'https://example.org/upd',
        'vcs_info': {'vcs': 'git', 'commit_id': 'ab' * 20},
    },
    'an archive of the repository, which names no commit': {'url': roh_upd.UPD_REPOSITORY, 'archive_info': {}},
}


@pytest.mark.parametrize('case', UPD_REFUSED)
def test_an_upd_that_is_not_from_a_commit_of_its_repository_is_refused(
    monkeypatch: pytest.MonkeyPatch, case: str
) -> None:
    _install(monkeypatch, UPD_REFUSED[case])
    with pytest.raises(ValueError, match='was not installed from a commit of'):
        roh_upd.installed_upd()


def test_tools_without_upd_beside_them_are_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    def distribution(name: str) -> None:
        raise importlib.metadata.PackageNotFoundError(name)

    monkeypatch.setattr(importlib.metadata, 'distribution', distribution)
    with pytest.raises(ValueError, match='no upd is installed beside these tools'):
        roh_upd.installed_upd()
