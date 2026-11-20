"""A family's runs of homozygosity and regions of uniparental disomy, from two published tools (docs/family-tables.md).

`bcftools roh` reports the runs of homozygosity of every sequenced member, and `upd` the regions in which a child with
both parents sequenced shows the alleles of one parent alone. Both run as they are published. This module says what
they are given and turns what they write into rows:

- `plan` reads a family's PED file, its VCF's header and the echtvar archive of gnomAD's joint release, and gives the
  INFO field that holds an allele's frequency and the children that `upd` assesses;
- the step (ROH_UPD in modules/family_tables.nf) prepares one VCF for both tools with bcftools and echtvar, and runs
  the tools on it;
- `keeps_a_site` says whether upd would keep any site of the prepared VCF, by upd's own rule. upd fails on a VCF of
  which it keeps none, as a call set of a few regions can be, so the step does not run it then, and the record lists
  the family's children as not assessed;
- `convert` reads the tools' output and writes a table of each, `roh` and `upd`, with a record of what ran.

The prepared VCF holds the family's PASS sites on the autosomes with one alternate allele of one base, genotypes
unphased, no INFO field of the family's VCF, and gnomAD's frequency of each allele. Both tools read a frequency of 0
as a site to leave out, which is what echtvar writes for an allele that the archive lacks: `plan` refuses an archive
that writes another stand-in. X, Y and the mitochondrial genome are left out, because a haploid call is no genotype
that the tools read, and written as a diploid one it makes a son's whole X a region of maternal disomy.

`upd` gives each region a type and counts of heterozygous and homozygous sites. It counts those sites up to the next
site that shows both parents' alleles, which lies beyond the region, so the type and the counts are wrong for a region
that covers part of a chromosome, and the tables leave them out. A run of homozygosity of the same child that overlaps
the region says what the type would have: the region is then an isodisomy or a deletion, and otherwise a heterodisomy.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import importlib.metadata
import json
import pathlib
import re
from collections.abc import Collection, Mapping, Sequence
from typing import Any, Literal

import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam

from case_annotation import echtvar_values, pedigree

ROH_NAME = 'roh.parquet'
UPD_NAME = 'upd.parquet'
INFO_NAME = 'roh_upd.json'
# Where upd is published. PyPI's package of the same name is another project, so the tools refuse an upd that was
# installed from anywhere else.
UPD_REPOSITORY = 'https://github.com/bjhall/upd'
# The chromosomes that both tools read.
AUTOSOMES = tuple(f'chr{number}' for number in range(1, 23))
# The column of the variants table whose field of the archive holds an allele's frequency.
_FREQUENCY = 'gnomad_joint_af'
# What a VCF's header takes as the ID of an INFO field, which goes into the tools' command lines.
_INFO_KEY = re.compile(r'[A-Za-z_][0-9A-Za-z_.]*')

_STRING = pa.string()
_INT = pa.int32()

# One run of homozygosity of one member, as `bcftools roh -O r` writes it: the first and the last position of the run,
# both part of it, its length in bases, the number of sites in it and the average of their phred scores.
ROH = pa.schema(
    [
        pa.field('member', _STRING, nullable=False),
        pa.field('chrom', _STRING, nullable=False),
        pa.field('pos', _INT, nullable=False),
        pa.field('end_pos', _INT, nullable=False),
        pa.field('length', _INT, nullable=False),
        pa.field('markers', _INT, nullable=False),
        pa.field('quality', pa.float32(), nullable=False),
    ]
)
# One region of one child, as `upd regions` writes it. `pos` and `end_pos` are the first and the last site that show
# the alleles of one parent alone, and `outer_pos` and `outer_end_pos` upd's START_LOW and END_HIGH.
UPD = pa.schema(
    [
        pa.field('member', _STRING, nullable=False),
        pa.field('chrom', _STRING, nullable=False),
        pa.field('pos', _INT, nullable=False),
        pa.field('end_pos', _INT, nullable=False),
        pa.field('outer_pos', _INT, nullable=False),
        pa.field('outer_end_pos', _INT, nullable=False),
        pa.field('origin', _STRING, nullable=False),
        pa.field('informative_sites', _INT, nullable=False),
        pa.field('snps', _INT, nullable=False),
        pa.field('opposite_sites', _INT, nullable=False),
    ]
)

type Origin = Literal['maternal', 'paternal']
_ORIGINS: Mapping[str, Origin] = {'MATERNAL': 'maternal', 'PATERNAL': 'paternal'}
# The keys of the fourth column of a line of `upd regions`, in the order upd writes them.
_REGION_KEYS = (
    'ORIGIN',
    'TYPE',
    'LOW_SIZE',
    'INF_SITES',
    'SNPS',
    'HET_HOM',
    'OPP_SITES',
    'START_LOW',
    'END_HIGH',
    'HIGH_SIZE',
)


class Trio(pydantic.BaseModel, frozen=True, extra='forbid'):
    """A child whose parents are both sequenced, by the three samples' names in the family's VCF."""

    child: str
    mother: str
    father: str


# Why upd did not assess a child: it would have kept no site of the prepared VCF (keeps_a_site).
type Reason = Literal['no informative site']


class NotAssessed(Trio, frozen=True, extra='forbid'):
    """A child whose parents are both sequenced and that upd was not run for, with the reason."""

    reason: Reason


class FrequencyArchive(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The echtvar archive that gave each site its frequency, and the echtvar that read it."""

    archive: str  # the file's name
    bytes: int
    sha256: str
    echtvar_version: str = pydantic.Field(min_length=1)
    field: echtvar_values.ArchiveField  # the field whose alias the tools read as the frequency


class RohUpdInfo(pydantic.BaseModel, frozen=True, extra='forbid'):
    """What made a family's roh and upd tables."""

    bcftools_version: str = pydantic.Field(min_length=1)
    upd_version: str = pydantic.Field(min_length=1)  # as the installed upd states it
    upd_commit: str = pydantic.Field(pattern=r'^[0-9a-f]{40}$')  # the commit of upd's repository it was installed from
    # The commands that wrote the VCF both tools read, in order, each without its files.
    preparation: tuple[str, ...]
    roh_command: str = pydantic.Field(min_length=1)  # without its files
    upd_command: str = pydantic.Field(min_length=1)  # without its files and its three samples
    frequency_archive: FrequencyArchive
    sites: int  # the records of the prepared VCF, which is 0 where the family has no such site
    roh_members: tuple[str, ...]  # every sequenced member, in the VCF's order
    # The children that upd assessed, in the VCF's order.
    upd_children: tuple[Trio, ...]
    # The children with both parents sequenced that upd was not run for, each with the reason, in the VCF's order. The
    # family has the upd table where it has a child of either list, and the table has no row of a child of this one.
    upd_not_assessed: tuple[NotAssessed, ...]


@dataclasses.dataclass(frozen=True)
class Plan:
    """What the step gives the tools: the INFO field of an allele's frequency, and the children upd assesses."""

    frequency: echtvar_values.ArchiveField
    trios: tuple[Trio, ...]


def trios(ped: pedigree.Pedigree, samples: Sequence[str]) -> tuple[Trio, ...]:
    """The family's children whose parents are both sequenced, in the order of the VCF's samples."""
    sequenced = set(samples)
    found = []
    for sample in samples:
        member = ped.member(sample)
        if member.mother in sequenced and member.father in sequenced:
            found.append(Trio(child=sample, mother=str(member.mother), father=str(member.father)))
    return tuple(found)


def frequency_field(archive: pathlib.Path) -> echtvar_values.ArchiveField:
    """The archive's field of gnomAD's joint allele frequency, which the tools read from each record.

    Raises:
        ValueError: if the file is no archive of gnomAD's joint release as the tables read it, if it writes another
            stand-in than 0 for the frequency of an allele it lacks, or if the field's alias is no INFO key.
    """
    fields = echtvar_values.column_fields(echtvar_values.GNOMAD_JOINT, echtvar_values.read_archive(archive), archive)
    field = fields[_FREQUENCY]
    if field.missing_value != 0:
        raise ValueError(
            f'{archive}: an allele it lacks gets the value {field.missing_value} as its {field.alias}, and the '
            'tools leave such a site out only where that value is 0'
        )
    if not _INFO_KEY.fullmatch(field.alias):
        raise ValueError(f'{archive}: {field.alias!r} is no name of an INFO field')
    return field


def _samples(vcf_path: pathlib.Path) -> tuple[list[str], frozenset[str]]:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        return list(vcf.header.samples), frozenset(vcf.header.formats)


def plan(family: str, ped_path: pathlib.Path, vcf_path: pathlib.Path, archive: pathlib.Path) -> Plan:
    """What the tools are given for a family, from its PED file, its VCF's header and the archive.

    Raises:
        ValueError: if the PED file is not the family's or lacks a sample of the VCF, if the archive is refused
            (frequency_field), or if the family has a child with both parents sequenced and its VCF no GQ: upd reads
            a call without GQ as one of quality 0 and would then report no region whatever the genotypes.
    """
    samples, formats = _samples(vcf_path)
    ped = pedigree.read(ped_path)
    pedigree.check_family(ped, family, samples)
    found = trios(ped, samples)
    if found and 'GQ' not in formats:
        raise ValueError(
            f'{vcf_path}: no FORMAT field GQ, which upd needs of {found[0].child} and its parents: it leaves out a '
            'site where a call has a GQ below 30, and reads a call without GQ as one of GQ 0'
        )
    return Plan(frequency=frequency_field(archive), trios=found)


def installed_upd() -> tuple[str, str]:
    """The version that the installed upd states, and the commit of its repository that it was installed from.

    Raises:
        ValueError: if no upd is installed beside this package, or one that did not come from a commit of upd's
            repository, as PyPI's package of the same name does not.
    """
    try:
        distribution = importlib.metadata.distribution('upd')
    except importlib.metadata.PackageNotFoundError as missing:
        raise ValueError('no upd is installed beside these tools, as the tools image has it') from missing
    origin = json.loads(distribution.read_text('direct_url.json') or '{}')
    commit = origin.get('vcs_info', {}).get('commit_id')
    if origin.get('url') != UPD_REPOSITORY or commit is None:
        raise ValueError(
            f'the installed upd {distribution.version} was not installed from a commit of {UPD_REPOSITORY}: '
            "PyPI's package of that name is another project"
        )
    return distribution.version, commit


def _integer(text: str, what: str, where: str) -> int:
    if not re.fullmatch(r'-?[0-9]+', text):
        raise ValueError(f'{where}: {what} is {text!r}, not an integer')
    return int(text)


def read_roh(path: pathlib.Path, samples: Sequence[str]) -> list[dict[str, Any]]:
    """The runs of homozygosity in the output of `bcftools roh -O r`, one row for each of its RG lines.

    Raises:
        ValueError: on a line that is no comment and no RG line of eight fields, on a run of a sample that the VCF
            lacks or off the autosomes, and on a length that is not the run's positions counted with both ends, which
            is how the table states them.
    """
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        where = f'{path}:{number}'
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) != 8 or fields[0] != 'RG':
            raise ValueError(f'{where}: not a line of bcftools roh with the runs alone: {line!r}')
        _, member, chrom, start, end, length, markers, quality = fields
        if member not in samples:
            raise ValueError(f'{where}: a run of {member!r}, which is no sample of the prepared VCF')
        if chrom not in AUTOSOMES:
            raise ValueError(f'{where}: a run on {chrom!r}, and the tools are given the autosomes alone')
        row = {
            'member': member,
            'chrom': chrom,
            'pos': _integer(start, 'the start', where),
            'end_pos': _integer(end, 'the end', where),
            'length': _integer(length, 'the length', where),
            'markers': _integer(markers, 'the number of markers', where),
        }
        if row['length'] != row['end_pos'] - row['pos'] + 1 or row['pos'] < 1:
            raise ValueError(f'{where}: the run from {start} to {end} is given the length {length}')
        try:
            row['quality'] = float(quality)
        except ValueError as malformed:
            raise ValueError(f'{where}: the quality is {quality!r}, not a number') from malformed
        rows.append(row)
    return rows


def read_upd(path: pathlib.Path, child: str) -> list[dict[str, Any]]:
    """The regions of one child in the output of `upd regions`, one row for each of its lines.

    A line is a BED interval, whose start is one below the first position of the region, with upd's values of the
    region as KEY=VALUE pairs in a fourth column. upd's TYPE and HET_HOM are read and not kept (see the module's
    opening).

    Raises:
        ValueError: on a line that is not of that form, with other keys than upd's, off the autosomes or with a size
            that is not the distance between its bounds, which is why the table can leave the sizes out.
    """
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        where = f'{path}:{number}'
        fields = line.split('\t')
        if len(fields) != 4:
            raise ValueError(f'{where}: not a line of upd regions, which has four columns: {line!r}')
        chrom, start, end, pairs = fields
        parts = [pair.partition('=') for pair in pairs.split(';')]
        if tuple(key for key, _, _ in parts) != _REGION_KEYS or not all(separator for _, separator, _ in parts):
            raise ValueError(f'{where}: the values of the region are not those of upd regions: {pairs!r}')
        values = {key: value for key, _, value in parts}
        if chrom not in AUTOSOMES:
            raise ValueError(f'{where}: a region on {chrom!r}, and the tools are given the autosomes alone')
        if values['ORIGIN'] not in _ORIGINS:
            raise ValueError(f'{where}: the origin is {values["ORIGIN"]!r}, not MATERNAL or PATERNAL')
        row = {
            'member': child,
            'chrom': chrom,
            'pos': _integer(start, 'the start', where) + 1,
            'end_pos': _integer(end, 'the end', where),
            'outer_pos': _integer(values['START_LOW'], 'START_LOW', where),
            'outer_end_pos': _integer(values['END_HIGH'], 'END_HIGH', where),
            'origin': _ORIGINS[values['ORIGIN']],
            'informative_sites': _integer(values['INF_SITES'], 'INF_SITES', where),
            'snps': _integer(values['SNPS'], 'SNPS', where),
            'opposite_sites': _integer(values['OPP_SITES'], 'OPP_SITES', where),
        }
        sizes = (_integer(values['LOW_SIZE'], 'LOW_SIZE', where), _integer(values['HIGH_SIZE'], 'HIGH_SIZE', where))
        if sizes != (row['end_pos'] - row['pos'], row['outer_end_pos'] - row['outer_pos']):
            raise ValueError(f'{where}: the sizes {sizes} are not the distances between the bounds: {line!r}')
        rows.append(row)
    return rows


def _kept(line: str, tag: str, min_af: float, min_gq: int, where: str) -> bool:
    """Whether upd keeps a record of the prepared VCF, by the rule of its code at the commit that the image installs.

    upd keeps a record whose frequency is at least `min_af` and whose GQ is at least `min_gq` in every sample of the
    VCF, the ones outside the trio too. It reads a record without the frequency as one of frequency 0, and a call
    without GQ, or with one that is no integer, as one of GQ 0.
    """
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 10:
        raise ValueError(f'{where}: not a record of a VCF with samples: {line!r}')
    info = {key: value for key, _, value in (pair.partition('=') for pair in fields[7].split(';'))}
    try:
        frequency = float(info.get(tag) or 0)
    except ValueError as malformed:
        raise ValueError(f'{where}: the frequency {tag} is {info[tag]!r}, not a number') from malformed
    if frequency < min_af:
        return False
    keys = fields[8].split(':')
    for call in fields[9:]:
        try:
            quality = int(dict(zip(keys, call.split(':'), strict=False)).get('GQ', '0'))
        except ValueError:
            quality = 0
        if quality < min_gq:
            return False
    return True


def keeps_a_site(prepared: pathlib.Path, tag: str, min_af: float, min_gq: int) -> bool:
    """Whether upd, with these thresholds, would keep at least one site of the prepared VCF.

    upd fails with a TypeError on a VCF of which it keeps no site, so the step asks before it runs upd. The answer
    follows upd's own reading of the file: its text, line by line, with the rule of _kept. upd never reads the last
    record of a file, so a site counts only where a record follows it. The reading stops at the first site that
    counts, which a whole genome has among its first records.

    Raises:
        ValueError: if a record is not one of a VCF with samples, or its frequency is no number.
    """
    with gzip.open(prepared, 'rt') as vcf:
        records = ((number, line) for number, line in enumerate(vcf, start=1) if not line.startswith('#'))
        previous = None
        for record in records:
            # `previous` has a record after it.
            if previous is not None and _kept(previous[1], tag, min_af, min_gq, f'{prepared}:{previous[0]}'):
                return True
            previous = record
    return False


def _digest(path: pathlib.Path) -> str:
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def convert(
    family: str,
    ped_path: pathlib.Path,
    prepared: pathlib.Path,
    archive: pathlib.Path,
    roh: pathlib.Path,
    regions: Mapping[str, pathlib.Path],
    output: pathlib.Path,
    *,
    not_assessed: Collection[str] = (),
    sites: int,
    bcftools_version: str,
    echtvar_version: str,
    upd: tuple[str, str],
    preparation: Sequence[str],
    roh_command: str,
    upd_command: str,
) -> RohUpdInfo:
    """Writes the family's runs and regions, and the record of what made them, into the output folder.

    The inputs are the family's PED file, the prepared VCF that both tools read, the archive that gave its sites their
    frequency, the output of `bcftools roh` and that of `upd regions` for each assessed child, by the child's sample.
    `not_assessed` are the children that the step did not run upd for, because upd would have kept no site
    (keeps_a_site). The folder gets the runs as ROH_NAME, the regions as UPD_NAME where the family has a child with
    both parents sequenced, assessed or not, and the record as INFO_NAME. `upd` is the installed upd's version and
    commit (installed_upd), and the commands are what the step ran.

    Raises:
        ValueError: if the number of sites is below 0, if the children with regions and those not assessed are not
            together exactly the children that the family's pedigree and the VCF's samples give, each once, or if a
            tool's output is not what the tool writes.
    """
    samples, _formats = _samples(prepared)
    ped = pedigree.read(ped_path)
    pedigree.check_family(ped, family, samples)
    if sites < 0:
        raise ValueError(f'{prepared}: {sites} is no number of records')
    with_parents = trios(ped, samples)
    children = [trio.child for trio in with_parents]
    if sorted([*regions, *not_assessed]) != sorted(children):
        raise ValueError(
            f'regions of {sorted(regions)} were given and {sorted(not_assessed)} as not assessed, and the children '
            f'of {family!r} with both parents sequenced are {sorted(children)}'
        )
    assessed = tuple(trio for trio in with_parents if trio.child in regions)
    field = frequency_field(archive)
    runs = read_roh(roh, samples)
    found = [row for trio in assessed for row in read_upd(regions[trio.child], trio.child)]
    output.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(runs, schema=ROH), output / ROH_NAME, compression='zstd')
    if with_parents:
        pq.write_table(pa.Table.from_pylist(found, schema=UPD), output / UPD_NAME, compression='zstd')
    info = RohUpdInfo(
        bcftools_version=bcftools_version,
        upd_version=upd[0],
        upd_commit=upd[1],
        preparation=tuple(preparation),
        roh_command=roh_command,
        upd_command=upd_command,
        frequency_archive=FrequencyArchive(
            archive=archive.name,
            bytes=archive.stat().st_size,
            sha256=_digest(archive),
            echtvar_version=echtvar_version,
            field=field,
        ),
        sites=sites,
        roh_members=tuple(samples),
        upd_children=assessed,
        upd_not_assessed=tuple(
            NotAssessed(**trio.model_dump(), reason='no informative site')
            for trio in with_parents
            if trio.child not in regions
        ),
    )
    (output / INFO_NAME).write_text(info.model_dump_json(indent=2) + '\n')
    return info


def read_info(folder: pathlib.Path) -> RohUpdInfo:
    """The record in a folder `convert` wrote."""
    return RohUpdInfo.model_validate_json((folder / INFO_NAME).read_text())


def check(info: RohUpdInfo, ped: pedigree.Pedigree, samples: Sequence[str]) -> None:
    """Fails unless a folder's record is of a family with this pedigree and these samples.

    Raises:
        ValueError: if the runs are of other members than the family's samples, or the children that upd assessed
            and those it did not are not together the children that the pedigree and the samples give, each once.
    """
    if list(info.roh_members) != list(samples):
        raise ValueError(
            f'the runs of homozygosity are of {list(info.roh_members)}, and the samples of {ped.family!r} are '
            f'{list(samples)}'
        )
    recorded = [*info.upd_children, *(Trio(**entry.model_dump(exclude={'reason'})) for entry in info.upd_not_assessed)]
    expected = trios(ped, samples)
    if len(recorded) != len(expected) or set(recorded) != set(expected):
        raise ValueError(
            f'upd assessed {[trio.child for trio in info.upd_children]} and did not assess '
            f'{[entry.child for entry in info.upd_not_assessed]}, and the children of {ped.family!r} with both '
            f'parents sequenced are {[trio.child for trio in expected]}'
        )
