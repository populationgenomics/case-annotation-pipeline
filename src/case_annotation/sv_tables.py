"""A family's structural variants and copy-number segments, from its callers' VCFs (docs/family-tables.md).

Two kinds of caller write such calls, and a family can bring a VCF of each. A structural-variant caller (`sv`) finds
junctions in the reads and writes one record per allele, as a sequence, a symbolic allele such as `<DEL>` or a breakend,
with a genotype per member. A copy-number caller (`cnv`) reads the depth and writes one record per segment, with a copy
number per member; DRAGEN's also writes the segments at the reference's copy number, which are left out here. The
records of both become rows of one table, `sv`, keyed by the caller's record id. A row says where the record lies, its
type and length, a breakend's mate, and its QUAL and FILTER, and then holds, as in a VCF and as the variants table does,
one column per sequenced member with that member's call of the record.

The family tables add to `sv` what Illumina's Annotator says of each record, and write a second table,
`sv_transcripts`, of the transcripts it lists for each (sv_annotation).

A member carries a structural variant by its genotype and a copy-number segment by its copy number: genotype_call and
copy_number_call define both. A member's call has the fields of its call of a small variant (family_tables.CALL)
wherever a field means the same thing, under the same name and type and in the same order.
"""

from __future__ import annotations

import pathlib
import re
from collections.abc import Collection, Iterator, Sequence
from typing import Any, Literal

import pyarrow as pa
import pysam

from case_annotation import vcf_values

type Caller = Literal['sv', 'cnv']

_STRING = pa.string()
_INT = pa.int32()

SV = pa.schema(
    [
        pa.field('chrom', _STRING, nullable=False),
        pa.field('pos', _INT, nullable=False),
        pa.field('end_pos', _INT, nullable=False),
        pa.field('sv_id', _STRING, nullable=False),
        pa.field('caller', _STRING, nullable=False),
        pa.field('sv_type', _STRING, nullable=False),
        pa.field('svlen', _INT),
        pa.field('ref', _STRING, nullable=False),
        pa.field('alt', _STRING, nullable=False),
        pa.field('qual', pa.float64()),
        pa.field('filters', pa.list_(_STRING)),
        pa.field('imprecise', pa.bool_(), nullable=False),
        pa.field('cipos', pa.list_(_INT)),
        pa.field('ciend', pa.list_(_INT)),
        pa.field('mate_id', _STRING),
        pa.field('mate_chrom', _STRING),
        pa.field('mate_pos', _INT),
        pa.field('bnd_extends', _STRING),
        pa.field('mate_extends', _STRING),
        pa.field('event_id', _STRING),
        pa.field('small_variant_id', _STRING),
    ]
)
# The columns that place and name a record, in every table with rows of a record.
RECORD_KEY = ('chrom', 'pos', 'end_pos', 'sv_id')
# The columns of sv that a record gives; small_variant_id comes from the family's variants table.
RECORD_COLUMNS = tuple(name for name in SV.names if name != 'small_variant_id')

# One member's call of one record: a member's column of the sv table holds one per row. gt, carries, zygosity, gq, pl,
# ft, dn and dq are the fields of a call of a small variant; the others are these callers' own.
CALL = pa.struct(
    [
        pa.field('gt', _STRING, nullable=False),
        pa.field('carries', pa.bool_()),
        pa.field('zygosity', _STRING),
        pa.field('gq', _INT),
        pa.field('pl', pa.list_(_INT)),
        pa.field('pr_ref', _INT),
        pa.field('pr_alt', _INT),
        pa.field('sr_ref', _INT),
        pa.field('sr_alt', _INT),
        pa.field('vf_ref', _INT),
        pa.field('vf_alt', _INT),
        pa.field('cn', _INT),
        pa.field('sm', pa.float64()),
        pa.field('bc', _INT),
        pa.field('pe_start', _INT),
        pa.field('pe_end', _INT),
        pa.field('qs', _INT),
        pa.field('ft', _STRING),
        pa.field('dn', _STRING),
        pa.field('dq', pa.float64()),
    ]
)

# The alleles of a copy-number segment, and the type of a segment by the alleles it holds.
_LOSS, _GAIN = '<DEL>', '<DUP>'
_SEGMENT_TYPES = {frozenset({_LOSS}): 'DEL', frozenset({_GAIN}): 'DUP', frozenset({_LOSS, _GAIN}): 'CNV'}

# A breakend's alternate allele (VCF 4.2, section 5.4): the bases at the position before or after the mate's position
# in brackets, which point to the side of it that is joined. A contig's name may hold colons, so the position is what
# follows the last one.
_BREAKEND = re.compile(
    r'(?P<before>[ACGTNacgtn]*)(?P<bracket>[\[\]])(?P<chrom>.+):(?P<pos>[0-9]+)(?P=bracket)(?P<after>[ACGTNacgtn]*)'
)
_BASES = re.compile(r'[ACGTNacgtn]+')


def schema(samples: Sequence[str]) -> pa.Schema:
    """The sv table of a family with these samples: the record's columns, then a call column per sample."""
    return pa.schema([*SV, *(pa.field(sample, CALL, nullable=False) for sample in samples)])


def transcripts_schema(columns: pa.Schema) -> pa.Schema:
    """The sv_transcripts table: the record's place and id, as in sv, then the columns of one of its transcripts."""
    return pa.schema([*(SV.field(name) for name in RECORD_KEY), *columns])


def genotype_call(alleles: Sequence[int | None]) -> tuple[bool | None, vcf_values.Zygosity | None]:
    """Whether a member carries a structural variant, and its zygosity, by its genotype.

    Both are those of a small variant's call (vcf_values.carries, vcf_values.zygosity). The member carries the variant
    if its genotype holds a copy of the allele, also beside a missing one, and does not if every allele is called and
    none is the allele; a genotype with an allele missing and no copy of the allele says neither, and has no zygosity.
    A haploid call of the allele, as on a male's X or Y, is `hom_alt`, and a haploid call of the reference `hom_ref`.
    """
    return vcf_values.carries(alleles), vcf_values.zygosity(alleles)


def copy_number_call(
    alleles: Sequence[int | None], cn: int | None, alts: Sequence[str], where: str
) -> tuple[bool | None, vcf_values.Zygosity | None]:
    """Whether a member carries a copy-number segment, and its zygosity, by its copy number.

    The copies expected are the alleles of the member's genotype: two, or one where the caller calls the member
    haploid, as on a male's X and Y outside the pseudo-autosomal regions. The member carries the segment if its copy
    number is another. The genotype alone does not say so: DRAGEN writes a member at the expected copy number as `./.`,
    and one above it as `./1`.

    The zygosity is `hom_ref` at the expected copy number, `hom_alt` at a copy number of 0 and `het` at any other below
    the expected one. Above it there is none, because the copy number does not say how the copies are split between the
    two chromosomes. A member without a copy number, such as a female on Y, has no call of the segment: the record
    says neither that the member carries it nor that it does not, and there is no zygosity.

    The genotype must say the same as the copy number, so that a caller with another convention fails here: an
    alternate allele exactly where the copy number is not the expected one, `<DEL>` below it and `<DUP>` above, and as
    many `<DEL>` alleles as copies are lost, so that `gt` and the zygosity agree.
    """
    expected = len(alleles)
    named = {alts[a - 1] for a in alleles if a}
    gt = '/'.join('.' if a is None else str(a) for a in alleles)
    if cn is None:
        if any(a is not None for a in alleles):
            raise ValueError(
                f'{where}: the genotype {gt} has no copy number, which says whether the segment is carried'
            )
        return None, None
    if cn < 0:
        raise ValueError(f'{where}: a copy number of {cn}')
    if cn == expected:
        if named:
            raise ValueError(f'{where}: the genotype {gt} names {sorted(named)} at the expected copy number, {cn}')
        return False, 'hom_ref'
    change = _LOSS if cn < expected else _GAIN
    if named != {change}:
        raise ValueError(
            f'{where}: a copy number of {cn} where {expected} are expected is a {change}, and the genotype {gt} names '
            f'{sorted(named) or "no allele"}'
        )
    if cn > expected:
        return True, None
    if (lost := sum(1 for a in alleles if a)) != expected - cn:
        raise ValueError(
            f'{where}: a copy number of {cn} where {expected} are expected loses {expected - cn}, and the genotype '
            f'{gt} names {lost} lost'
        )
    return True, 'hom_alt' if cn == 0 else 'het'


def breakend(alt: str, where: str) -> dict[str, Any]:
    """A breakend's mate columns, from its alternate allele.

    `bnd_extends` is the side of the record's position whose sequence the junction keeps: `left` where the allele
    starts with the position's bases (`t[p[`, `t]p]`), `right` where it ends with them (`]p]t`, `[p[t`). `mate_extends`
    is the side of the mate's position that is joined to it: `right` for `[` and `left` for `]`.
    """
    match = _BREAKEND.fullmatch(alt)
    if match is None or bool(match['before']) == bool(match['after']):
        raise ValueError(f'{where}: {alt!r} is not a breakend with a mate, such as G]chr2:300] or [chr2:300[G')
    return {
        'mate_chrom': match['chrom'],
        'mate_pos': int(match['pos']),
        'bnd_extends': 'left' if match['before'] else 'right',
        'mate_extends': 'right' if match['bracket'] == '[' else 'left',
    }


def is_sequence(alt: str) -> bool:
    """Whether an alternate allele is written as bases, and so can be the allele of a small variant."""
    return _BASES.fullmatch(alt) is not None


def _interval(record: pysam.VariantRecord, key: str, where: str) -> list[int] | None:
    held = vcf_values.info(record, key)
    if held is None:
        return None
    if not (isinstance(held, tuple) and len(held) == 2 and all(isinstance(v, int) for v in held)):
        raise ValueError(f'{where}: {key} {held!r} is not an interval of two integers')
    return list(held)


def _one(record: pysam.VariantRecord, key: str, where: str) -> str | None:
    """An INFO value that names one thing, None where the record has none; several fail."""
    held = vcf_values.info(record, key)
    if held is None or isinstance(held, str):
        return held
    if isinstance(held, tuple) and len(held) == 1 and isinstance(held[0], str):
        return held[0]
    raise ValueError(f'{where}: {key} {held!r} is not one name')


def _length(record: pysam.VariantRecord, where: str) -> int | None:
    """The record's SVLEN where it gives one length; a segment of a loss and a gain gives one for each."""
    held = vcf_values.info(record, 'SVLEN')
    if held is None:
        return None
    lengths = [v for v in (held if isinstance(held, tuple) else (held,)) if isinstance(v, int)]
    if not lengths or len(lengths) != (len(held) if isinstance(held, tuple) else 1):
        raise ValueError(f'{where}: SVLEN {held!r} is not a length')
    return lengths[0] if len(lengths) == 1 else None


def _site(record: pysam.VariantRecord, caller: Caller, where: str) -> dict[str, Any]:
    """A record's row of sv, without the link to a small variant."""
    alts = record.alts or ()
    if record.id is None:
        raise ValueError(f'{where}: a record without an ID, which is its key in the tables')
    mate: dict[str, Any] = {}
    if caller == 'cnv':
        if (sv_type := _SEGMENT_TYPES.get(frozenset(alts))) is None:
            raise ValueError(f"{where}: a copy-number segment's alleles are {_LOSS} and {_GAIN}, not {list(alts)}")
    else:
        if len(alts) != 1:
            raise ValueError(f'{where}: a structural variant record holds one alternate allele, not {list(alts)}')
        sv_type = vcf_values.info(record, 'SVTYPE')
        if not isinstance(sv_type, str):
            raise ValueError(f'{where}: no SVTYPE, the type of the structural variant')
        is_breakend = '[' in alts[0] or ']' in alts[0]
        if is_breakend != (sv_type == 'BND'):
            raise ValueError(f'{where}: SVTYPE is {sv_type} and the alternate allele is {alts[0]!r}')
        if is_breakend:
            mate = breakend(alts[0], where) | {'mate_id': _one(record, 'MATEID', where)}
        elif not (is_sequence(alts[0]) or (alts[0].startswith('<') and alts[0].endswith('>'))):
            raise ValueError(f'{where}: {alts[0]!r} is not a sequence, a symbolic allele or a breakend with a mate')
    return {
        'chrom': record.chrom,
        'pos': record.pos,
        # htslib's end of the record: END, or for a symbolic allele other than an insertion the position plus the
        # length in SVLEN where that lies further; without either, the last base of REF.
        'end_pos': record.stop,
        'sv_id': record.id,
        'caller': caller,
        'sv_type': sv_type,
        'svlen': _length(record, where),
        'ref': record.ref,
        'alt': ','.join(alts),
        'qual': record.qual,
        'filters': list(record.filter.keys()),
        'imprecise': vcf_values.info(record, 'IMPRECISE') is True,
        'cipos': _interval(record, 'CIPOS', where),
        'ciend': _interval(record, 'CIEND', where),
        'event_id': _one(record, 'EVENT', where),
        **mate,
    }


def _call(
    sample: pysam.VariantRecordSample, caller: Caller, alts: Sequence[str], formats: frozenset[str], where: str
) -> dict[str, Any]:
    """One member's call of a record, as a value of its column (CALL)."""
    alleles, gt = vcf_values.genotype(sample, where)
    cn = vcf_values.integer(sample, 'CN', formats, where)
    carries, zygosity = genotype_call(alleles) if caller == 'sv' else copy_number_call(alleles, cn, alts, where)
    counts = 'a count for the reference and one for the allele'
    pr_ref, pr_alt = vcf_values.pair(sample, 'PR', formats, where, counts)
    sr_ref, sr_alt = vcf_values.pair(sample, 'SR', formats, where, counts)
    vf_ref, vf_alt = vcf_values.pair(sample, 'VF', formats, where, counts)
    pe_start, pe_end = vcf_values.pair(sample, 'PE', formats, where, 'a count at the start and one at the end')
    return {
        'gt': gt,
        'carries': carries,
        'zygosity': zygosity,
        'gq': vcf_values.integer(sample, 'GQ', formats, where),
        'pl': vcf_values.likelihoods(sample, formats, len(alleles), where) if caller == 'sv' else None,
        'pr_ref': pr_ref,
        'pr_alt': pr_alt,
        'sr_ref': sr_ref,
        'sr_alt': sr_alt,
        'vf_ref': vf_ref,
        'vf_alt': vf_alt,
        'cn': cn,
        'sm': vcf_values.number(sample, 'SM', formats, where),
        'bc': vcf_values.integer(sample, 'BC', formats, where),
        'pe_start': pe_start,
        'pe_end': pe_end,
        'qs': vcf_values.integer(sample, 'QS', formats, where),
        'ft': vcf_values.text(sample, 'FT', formats, where),
        'dn': vcf_values.text(sample, 'DN', formats, where),
        'dq': vcf_values.number(sample, 'DQ', formats, where),
    }


def check_header(path: pathlib.Path, header: pysam.VariantHeader, caller: Caller, samples: Collection[str]) -> None:
    """Fails unless the VCF's samples are the family's and it has what the caller's calls are read from.

    Every sequenced member has a column in sv, so the samples are exactly the family's. A structural-variant
    caller's VCF needs genotypes, and a copy-number caller's needs copy numbers too.
    """
    found = list(header.samples)
    if set(found) != set(samples) or len(found) != len(samples):
        raise ValueError(f"{path}: its samples are {found}, and the family's VCF has {sorted(samples)}")
    needed = ('GT',) if caller == 'sv' else ('GT', 'CN')
    if missing := [key for key in needed if key not in header.formats]:
        raise ValueError(f'{path}: no {" or ".join(missing)} among its FORMAT fields')


def check(path: pathlib.Path, caller: Caller, samples: Collection[str]) -> None:
    """Fails unless a caller's VCF fits the family (check_header), from its header alone."""
    with pysam.VariantFile(str(path)) as vcf:
        check_header(path, vcf.header, caller, samples)


def rows(
    path: pathlib.Path, caller: Caller, samples: Sequence[str], chromosomes: Collection[str]
) -> Iterator[tuple[dict[str, Any], dict[str, Any]]]:
    """Each record's row of sv, for the records on the chromosomes given: its own columns, and each sample's call.

    A copy-number caller's segments at the reference's copy number, which have no alternate allele, are left out. The
    calls are in the order of the samples given.
    """
    with pysam.VariantFile(str(path)) as vcf:
        check_header(path, vcf.header, caller, samples)
        formats = frozenset(vcf.header.formats)
        for record in vcf:
            if record.chrom not in chromosomes or (caller == 'cnv' and not record.alts):
                continue
            where = f'{path}: {record.chrom}:{record.pos}'
            site = _site(record, caller, where)
            alts = record.alts or ()
            yield site, {s: _call(record.samples[s], caller, alts, formats, where) for s in samples}
