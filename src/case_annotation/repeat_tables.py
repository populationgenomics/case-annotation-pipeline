"""A family's repeat genotypes, from each member's repeat VCF (docs/family-tables.md).

A repeat genotyper such as DRAGEN's, which runs ExpansionHunter, sizes the short tandem repeats of a catalogue of loci
in one sample and writes one single-sample VCF, with a record per repeat of the catalogue whether or not the sample's
alleles differ from the reference. The members' files become one table, `repeats`. It has one row per repeat, keyed by
the catalogue's id of it (VARID, which DRAGEN's files give again as REPID), and, as in a VCF, one column per sequenced
member, which holds that member's two alleles, or its one where the caller calls it haploid.
"""

from __future__ import annotations

import pathlib
from collections.abc import Collection, Sequence
from typing import Any

import pyarrow as pa
import pysam

from case_annotation import vcf_values

_STRING = pa.string()
_INT = pa.int32()
_INTS = pa.list_(_INT)

# One member's genotype of one repeat: a column of the repeats table holds one per row. The lists hold a value per
# allele, in the genotype's order.
CALL = pa.struct(
    [
        pa.field('gt', _STRING, nullable=False),
        pa.field('repcn', _INTS),
        pa.field('repci_lower', _INTS),
        pa.field('repci_upper', _INTS),
        pa.field('so', pa.list_(_STRING)),
        pa.field('adsp', _INTS),
        pa.field('adfl', _INTS),
        pa.field('adir', _INTS),
        pa.field('lc', pa.float64()),
        pa.field('filters', pa.list_(_STRING)),
    ]
)
# The columns that describe the repeat, which every member's file must give alike; a member's column follows them.
LOCUS = pa.schema(
    [
        pa.field('chrom', _STRING, nullable=False),
        pa.field('pos', _INT, nullable=False),
        pa.field('end_pos', _INT, nullable=False),
        pa.field('repeat_id', _STRING, nullable=False),
        pa.field('repeat_unit', _STRING),
        pa.field('ref_copies', _INT),
    ]
)


def schema(samples: Sequence[str]) -> pa.Schema:
    """The repeats table of a family with these samples: the repeat, then a genotype column per sample."""
    return pa.schema([*LOCUS, *(pa.field(sample, CALL, nullable=False) for sample in samples)])


def _per_allele(
    sample: pysam.VariantRecordSample, key: str, formats: frozenset[str], ploidy: int, where: str
) -> list[str] | None:
    """A field that holds a value per allele, separated by `/`; None where every value is missing."""
    held = vcf_values.text(sample, key, formats, where)
    if held is None:
        return None
    parts = held.split('/')
    if all(part == '.' for part in parts):
        return None
    if len(parts) != ploidy or '.' in parts:
        raise ValueError(f'{where}: {key} {held!r} does not hold a value for each of the {ploidy} alleles')
    return parts


def _integers(parts: Sequence[str] | None, key: str, where: str) -> list[int] | None:
    if parts is None:
        return None
    if not all(part.isascii() and part.isdigit() for part in parts):
        raise ValueError(f'{where}: {key} {"/".join(parts)!r} does not hold integers')
    return [int(part) for part in parts]


def _intervals(parts: Sequence[str] | None, where: str) -> tuple[list[int] | None, list[int] | None]:
    """The lower and the upper bounds of the confidence intervals, each written as `lower-upper`."""
    if parts is None:
        return None, None
    bounds = [part.split('-') for part in parts]
    if not all(len(b) == 2 and all(v.isascii() and v.isdigit() for v in b) for b in bounds):
        raise ValueError(f'{where}: REPCI {"/".join(parts)!r} does not hold an interval, lower-upper, for each allele')
    return [int(lower) for lower, _ in bounds], [int(upper) for _, upper in bounds]


def _call(record: pysam.VariantRecord, sample: str, formats: frozenset[str], where: str) -> dict[str, Any]:
    """The member's genotype of a repeat, as a value of its column (CALL)."""
    held = record.samples[sample]
    alleles, gt = vcf_values.genotype(held, where)
    ploidy = len(alleles)
    lower, upper = _intervals(_per_allele(held, 'REPCI', formats, ploidy, where), where)
    return {
        'gt': gt,
        'repcn': _integers(_per_allele(held, 'REPCN', formats, ploidy, where), 'REPCN', where),
        'repci_lower': lower,
        'repci_upper': upper,
        'so': _per_allele(held, 'SO', formats, ploidy, where),
        'adsp': _integers(_per_allele(held, 'ADSP', formats, ploidy, where), 'ADSP', where),
        'adfl': _integers(_per_allele(held, 'ADFL', formats, ploidy, where), 'ADFL', where),
        'adir': _integers(_per_allele(held, 'ADIR', formats, ploidy, where), 'ADIR', where),
        'lc': vcf_values.number(held, 'LC', formats, where),
        'filters': list(record.filter.keys()),
    }


def _locus(record: pysam.VariantRecord, where: str) -> dict[str, Any]:
    repeat_id = vcf_values.info(record, 'VARID')
    if not isinstance(repeat_id, str):
        raise ValueError(f"{where}: no VARID, the catalogue's id of the repeat, which is its key in the table")
    ref_copies = vcf_values.info(record, 'REF')
    if ref_copies is not None and not isinstance(ref_copies, int):
        raise ValueError(f"{where}: REF {ref_copies!r} is not the reference's number of repeat units")
    return {
        'chrom': record.chrom,
        'pos': record.pos,
        'end_pos': record.stop,
        'repeat_id': repeat_id,
        'repeat_unit': vcf_values.info(record, 'RU'),
        'ref_copies': ref_copies,
    }


def _sample(path: pathlib.Path, header: pysam.VariantHeader) -> str:
    found = list(header.samples)
    if len(found) != 1:
        raise ValueError(f"{path}: a repeat VCF holds one member's genotypes, and this one has the samples {found}")
    if 'GT' not in header.formats or 'REPCN' not in header.formats:
        raise ValueError(f'{path}: no GT or REPCN among its FORMAT fields')
    return found[0]


def by_sample(paths: Sequence[pathlib.Path], samples: Collection[str]) -> dict[str, pathlib.Path]:
    """Each member's repeat VCF, by its sample, from the files' headers.

    Fails unless each file holds one sample and the files' samples are exactly the family's, one file each: every
    sequenced member has a column in the repeats table.
    """
    files: dict[str, pathlib.Path] = {}
    for path in paths:
        with pysam.VariantFile(str(path)) as vcf:
            sample = _sample(path, vcf.header)
        if sample in files:
            raise ValueError(f'{path}: holds the sample {sample!r}, as {files[sample]} does')
        files[sample] = path
    if extra := sorted(set(files) - set(samples)):
        raise ValueError(f"the repeat VCFs hold the samples {extra}, which the family's VCF lacks")
    if missing := sorted(set(samples) - set(files)):
        raise ValueError(f'no repeat VCF of the samples {missing}: a family brings one of every member, or none')
    return files


def rows(paths: Sequence[pathlib.Path], samples: Sequence[str], chromosomes: Collection[str]) -> list[dict[str, Any]]:
    """The rows of the repeats table, for the repeats on the chromosomes given, in the first sample's file's order.

    Fails unless every member's file holds the same repeats, each once and described alike: the files then come from
    one catalogue.
    """
    files = by_sample(paths, samples)
    loci: dict[str, dict[str, Any]] = {}
    calls: dict[str, dict[str, dict[str, Any]]] = {}
    for sample in samples:
        path = files[sample]
        found: dict[str, dict[str, Any]] = {}
        calls[sample] = {}
        with pysam.VariantFile(str(path)) as vcf:
            formats = frozenset(vcf.header.formats)
            for record in vcf:
                where = f'{path}: {record.chrom}:{record.pos}'
                locus = _locus(record, where)
                if locus['repeat_id'] in found:
                    raise ValueError(f'{where}: the repeat {locus["repeat_id"]!r} has two records')
                found[locus['repeat_id']] = locus
                calls[sample][locus['repeat_id']] = _call(record, sample, formats, where)
        if sample == samples[0]:
            loci = found
        elif found != loci:
            differing = sorted(k for k in found.keys() | loci.keys() if found.get(k) != loci.get(k))
            raise ValueError(
                f'{path}: its repeats are not those of {files[samples[0]]}, so the two are not of one catalogue: '
                f'{len(differing)} differ, such as {differing[:5]}'
            )
    return [
        locus | {sample: calls[sample][repeat_id] for sample in samples}
        for repeat_id, locus in loci.items()
        if locus['chrom'] in chromosomes
    ]
