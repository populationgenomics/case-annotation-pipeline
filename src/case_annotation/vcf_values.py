"""Values of a VCF record as the family tables read them: its INFO, and one sample's FORMAT (docs/family-tables.md).

A value is None where the header does not define the field, the record lacks it or it is missing; a FORMAT value of
another type than the table's column fails.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import pysam

type Zygosity = Literal['hom_ref', 'het', 'hom_alt']


def info(record: pysam.VariantRecord, key: str) -> object:
    """An INFO value of a record, None where the header or the record lacks it."""
    # pysam fails on a key the header does not define.
    return record.info.get(key) if key in record.header.info else None


def value(sample: pysam.VariantRecordSample, key: str, formats: frozenset[str]) -> object:
    """A FORMAT value of one sample, None where the header or the record lacks it or it is missing."""
    if key not in formats:
        return None
    held = sample.get(key)
    if isinstance(held, tuple):
        return None if all(v is None for v in held) else held
    return None if held == '.' else held


def integer(sample: pysam.VariantRecordSample, key: str, formats: frozenset[str], where: str) -> int | None:
    held = value(sample, key, formats)
    if held is not None and not isinstance(held, int):
        raise ValueError(f'{where}: {key} is {held!r}, not an integer')
    return held


def number(sample: pysam.VariantRecordSample, key: str, formats: frozenset[str], where: str) -> float | None:
    held = value(sample, key, formats)
    if held is not None and not isinstance(held, int | float):
        raise ValueError(f'{where}: {key} is {held!r}, not a number')
    return None if held is None else float(held)


def text(sample: pysam.VariantRecordSample, key: str, formats: frozenset[str], where: str) -> str | None:
    held = value(sample, key, formats)
    if held is not None and not isinstance(held, str):
        raise ValueError(f'{where}: {key} is {held!r}, not a string')
    return held


def pair(
    sample: pysam.VariantRecordSample, key: str, formats: frozenset[str], where: str, of: str
) -> tuple[int | None, int | None]:
    """A field of two integers, such as the reads of the reference and of the allele; two Nones where it is missing.

    `of` says what the two are, for the message of a value that is not two integers.
    """
    held = value(sample, key, formats)
    if held is None:
        return None, None
    if not (isinstance(held, tuple) and len(held) == 2 and all(v is None or isinstance(v, int) for v in held)):
        raise ValueError(f'{where}: {key} {held} is not {of}')
    return held[0], held[1]


def likelihoods(
    sample: pysam.VariantRecordSample, formats: frozenset[str], ploidy: int, where: str
) -> list[int] | None:
    """The call's genotype likelihoods (PL); None where they are missing.

    The record holds one alternate allele, so a call of n alleles has n + 1 likelihoods: 0/0, 0/1 and 1/1 for a diploid
    call, 0 and 1 for a haploid one. In a record that also holds diploid calls, a haploid call's list ends early, as
    bcftools writes it, or is padded with missing values, which are dropped.
    """
    held = value(sample, 'PL', formats)
    if held is None:
        return None
    if not isinstance(held, tuple):
        raise ValueError(f'{where}: PL {held!r} is not a list of likelihoods')
    kept, padding = held[: ploidy + 1], held[ploidy + 1 :]
    if len(kept) != ploidy + 1 or any(v is not None for v in padding):
        raise ValueError(
            f'{where}: PL {held} has {len(held)} likelihoods, and a call of {ploidy} alleles has {ploidy + 1}'
        )
    if not all(isinstance(v, int) for v in kept):
        raise ValueError(f'{where}: PL {held!r} is not a list of integer likelihoods, each present')
    return list(kept)


def genotype(sample: pysam.VariantRecordSample) -> tuple[tuple[int | None, ...], str]:
    """The call's alleles, None for a missing one, and the genotype as a VCF writes it."""
    # htslib reads an allele index the record lacks as missing.
    alleles: tuple[int | None, ...] = sample['GT']
    separator = '|' if sample.phased and len(alleles) > 1 else '/'
    return alleles, separator.join('.' if a is None else str(a) for a in alleles)


def carries(alleles: Sequence[int | None]) -> bool | None:
    """Whether a genotype of a record with one alternate allele holds a copy of that allele; None where it does not say.

    True where an allele of the genotype is the alternate allele, also beside a missing one (./1). False where every
    allele is called and none is the alternate allele. None where an allele is missing and no called one is the
    alternate allele (./., a haploid ., 0/.): the missing allele may be the alternate one, so the genotype says neither.
    """
    if 1 in alleles:
        return True
    return None if None in alleles else False


def zygosity(alleles: Sequence[int | None]) -> Zygosity | None:
    """The zygosity of a genotype of a record with one alternate allele; None where carries is None.

    `hom_ref` for a genotype that does not hold the alternate allele, `hom_alt` where every allele is the alternate
    one, also in a haploid call, and `het` for any other genotype that holds it, also beside a missing allele (./1).
    """
    held = carries(alleles)
    if held is None:
        return None
    if not held:
        return 'hom_ref'
    return 'hom_alt' if all(a == 1 for a in alleles) else 'het'
