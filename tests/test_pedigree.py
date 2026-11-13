"""A family's PED file, on invented families."""

from __future__ import annotations

import pathlib

import pytest

from case_annotation import pedigree

_TRIO = 'fam-1\tS1\t0\t0\t1\t1\nfam-1\tS2\t0\t0\t2\t-9\nfam-1\tS3\tS1\tS2\t0\t2\n'


def _ped(tmp_path: pathlib.Path, text: str) -> pathlib.Path:
    path = tmp_path / 'family.ped'
    path.write_text(text)
    return path


def test_a_trio_reads_with_its_parents_sexes_and_phenotypes(tmp_path: pathlib.Path) -> None:
    ped = pedigree.read(_ped(tmp_path, '# a comment\n' + _TRIO))
    assert ped.family == 'fam-1'
    child = ped.member('S3')
    assert (child.father, child.mother, child.sex, child.affected) == ('S1', 'S2', 'unknown', True)
    assert ped.member('S2').affected is False  # -9, unknown, counts as unaffected


REJECTED = {
    'a line of five columns': ('fam-1\tS1\t0\t0\t1\n', 'six columns'),
    'two families': (_TRIO + 'fam-2\tS4\t0\t0\t1\t1\n', 'one family'),
    'a repeated individual': (_TRIO + 'fam-1\tS1\t0\t0\t1\t1\n', 'listed twice'),
    'a parent without a line': ('fam-1\tS3\tS1\t0\t1\t2\n', 'no line of its own'),
    'an individual who is their own parent': ('fam-1\tS1\tS1\t0\t1\t1\n', 'no line of its own'),
    'a female father': ('fam-1\tS1\t0\t0\t2\t1\nfam-1\tS3\tS1\t0\t1\t2\n', 'is female'),
    'a male mother': ('fam-1\tS2\t0\t0\t1\t1\nfam-1\tS3\t0\tS2\t1\t2\n', 'is male'),
    'an unknown sex code': ('fam-1\tS1\t0\t0\tM\t1\n', 'sex'),
    'an unknown phenotype code': ('fam-1\tS1\t0\t0\t1\t3\n', 'phenotype'),
    'a name with a space-free but unsafe character': ('fam/1\tS1\t0\t0\t1\t1\n', 'not a name'),
    'no individuals': ('# only a comment\n', 'no individuals'),
}


@pytest.mark.parametrize('case', REJECTED)
def test_a_malformed_ped_file_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    text, message = REJECTED[case]
    with pytest.raises(ValueError, match=message):
        pedigree.read(_ped(tmp_path, text))


def test_the_ped_file_must_be_the_familys(tmp_path: pathlib.Path) -> None:
    ped = pedigree.read(_ped(tmp_path, _TRIO))
    pedigree.check_family(ped, 'fam-1', ['S1', 'S2', 'S3'])
    with pytest.raises(ValueError, match='of family'):
        pedigree.check_family(ped, 'fam-2', ['S1', 'S2', 'S3'])


def test_every_vcf_sample_must_be_in_the_ped_file(tmp_path: pathlib.Path) -> None:
    ped = pedigree.read(_ped(tmp_path, _TRIO))
    with pytest.raises(ValueError, match='lacks the VCF samples'):
        pedigree.check_family(ped, 'fam-1', ['S1', 'S2', 'S3', 'S9'])


def test_a_vcf_without_samples_is_refused(tmp_path: pathlib.Path) -> None:
    ped = pedigree.read(_ped(tmp_path, _TRIO))
    with pytest.raises(ValueError, match='no samples'):
        pedigree.check_family(ped, 'fam-1', [])
