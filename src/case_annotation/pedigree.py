"""A family's PED file, as the family tables read it (docs/annotation.md, the input contract).

A PED file has six columns, separated by tabs or spaces: family, individual, father, mother, sex and phenotype. The
individual column holds the sample's name in the family's VCF, `0` stands for an absent parent, sex is `1` (male), `2`
(female) or `0` (unknown), and the phenotype is `2` (affected), `1` (unaffected), or `0` or `-9` (unknown). A line
that starts with `#` is a comment.
"""

from __future__ import annotations

import dataclasses
import pathlib
import re
from collections.abc import Collection
from typing import Literal

type Sex = Literal['male', 'female', 'unknown']

_SEXES: dict[str, Sex] = {'1': 'male', '2': 'female', '0': 'unknown'}
_AFFECTED = {'2': True, '1': False, '0': None, '-9': None}
# A name that each of the family's files, and a folder of its outputs, takes as it is.
NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9._-]*')


@dataclasses.dataclass(frozen=True)
class Member:
    sample: str
    father: str | None
    mother: str | None
    sex: Sex
    affected: bool | None  # None for an unknown phenotype


@dataclasses.dataclass(frozen=True)
class Pedigree:
    family: str
    members: tuple[Member, ...]

    def member(self, sample: str) -> Member:
        for m in self.members:
            if m.sample == sample:
                return m
        raise KeyError(f'{sample!r} is not in the family {self.family!r}')


def read(path: pathlib.Path) -> Pedigree:
    """The PED file's one family.

    Fails on a malformed line, a second family, a repeated or unknown individual, and a parent of the wrong sex.
    """
    rows: list[list[str]] = []
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split()
        if len(fields) != 6:
            raise ValueError(f'{path}:{number}: a PED line has six columns, and this one has {len(fields)}')
        rows.append(fields)
    if not rows:
        raise ValueError(f'{path}: no individuals')
    families = sorted({r[0] for r in rows})
    if len(families) != 1:
        raise ValueError(f'{path}: one family per PED file, and this one has {families}')
    members = []
    for family, sample, father, mother, sex, phenotype in rows:
        for name in (family, sample, father, mother):
            if not NAME.fullmatch(name):
                raise ValueError(f'{path}: {name!r} is not a name the outputs can carry')
        if sex not in _SEXES:
            raise ValueError(f'{path}: {sample}: sex {sex!r} is not 1, 2 or 0')
        if phenotype not in _AFFECTED:
            raise ValueError(f'{path}: {sample}: phenotype {phenotype!r} is not 2, 1, 0 or -9')
        members.append(
            Member(
                sample=sample,
                father=None if father == '0' else father,
                mother=None if mother == '0' else mother,
                sex=_SEXES[sex],
                affected=_AFFECTED[phenotype],
            )
        )
    pedigree = Pedigree(family=families[0], members=tuple(members))
    samples = [m.sample for m in members]
    if len(set(samples)) != len(samples):
        raise ValueError(f'{path}: an individual is listed twice')
    for m in members:
        for parent, wrong in ((m.father, 'female'), (m.mother, 'male')):
            if parent is None:
                continue
            if parent not in samples or parent == m.sample:
                raise ValueError(f'{path}: {m.sample}: the parent {parent!r} has no line of its own')
            if pedigree.member(parent).sex == wrong:
                raise ValueError(f'{path}: {m.sample}: the parent {parent!r} is {wrong}')
    return pedigree


def check_family(pedigree: Pedigree, family: str, vcf_samples: Collection[str]) -> None:
    """Fails unless the PED file is the family's and names every sample of its VCF.

    A PED individual the VCF lacks is a relative who was not sequenced.
    """
    if pedigree.family != family:
        raise ValueError(f'the PED file is of family {pedigree.family!r}, not {family!r}')
    if not vcf_samples:
        raise ValueError(f'the VCF of {family!r} has no samples')
    peds = {m.sample for m in pedigree.members}
    if missing := sorted(set(vcf_samples) - peds):
        raise ValueError(f'the PED file of {family!r} lacks the VCF samples {missing}')
