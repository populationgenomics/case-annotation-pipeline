"""STRchive's size ranges for the repeats of a family's repeats table (docs/annotation.md).

STRchive curates the tandem repeats that cause disease. A release publishes its loci as one JSON file, a list with an
entry for each locus: where the repeat lies on GRCh38, its motifs, and the numbers of repeat units between which an
allele counts as benign, as intermediate and as pathogenic. This module reads that file, checks it, and gives each row
of `repeats` the entry of its repeat as the columns FIELDS, which hold the entry's values as STRchive gives them. No
column says whether a member's allele is pathogenic: a reader compares the sizes with the ranges.

A row is matched to an entry by what both publish about the repeat, and by no name of it, so that the rule holds for
any genotyper's catalogue. The two lie on one chromosome, their intervals share a base, and the row's repeat unit is
one of the entry's motifs (same_repeat). A row that matches no entry has nulls, and a row that matches two is refused.
One more column says whether the row's interval is the entry's own: only then do a member's sizes and the entry's
ranges count the same stretch.

The file states neither its release nor a date. STRchive is a dated source, so a run that gives the file states both
(Release), and a family's manifest records them as stated.
"""

from __future__ import annotations

import dataclasses
import datetime as dt
import json
import pathlib
import re
from collections.abc import Mapping, Sequence
from typing import Any

import pyarrow as pa
import pydantic

from case_annotation import pedigree, sv_sites

# The source's name wherever a run meets it: the workflow's param, the manifest's key and the source of its columns in
# the contract.
NAME = 'strchive'

_STRING = pa.string()
_INT = pa.int32()
_INT_MAX = 2**31 - 1
_RANGES = ('benign', 'intermediate', 'pathogenic')
_BOUNDS = tuple(f'{kind}_{bound}' for kind in _RANGES for bound in ('min', 'max'))
# Whether a row's interval is its entry's, start and stop alike: the one column that no field of the file gives.
SAME_INTERVAL = 'strchive_same_interval'
# The columns of repeats that hold an entry's values, each named as STRchive's field with the prefix. The six bounds
# count repeat units of the row's repeat unit's length (Locus.motifs). A catalogue can draw the repeat's ends
# elsewhere than STRchive does. A member's size and a bound then count units of one length over different stretches,
# and nothing here converts one to the other: the entry's interval is a column, and SAME_INTERVAL says where the two
# intervals are one.
FIELDS = (
    pa.field('strchive_id', _STRING),
    pa.field('strchive_gene', _STRING),
    pa.field('strchive_evidence', _STRING),
    pa.field('strchive_inheritance', pa.list_(_STRING)),
    pa.field('strchive_start_hg38', _INT),
    pa.field('strchive_stop_hg38', _INT),
    pa.field(SAME_INTERVAL, pa.bool_()),
    pa.field('strchive_pathogenic_motif_reference_orientation', pa.list_(_STRING)),
    *(pa.field(f'strchive_{bound}', _INT) for bound in _BOUNDS),
)
COLUMNS = tuple(field.name for field in FIELDS)
# STRchive's levels of evidence for the link between a locus and its disease, as its site lists its classifications.
# The releases up to v2.28.0 use all of them but `No Known Relationship`.
LEVELS = frozenset(
    {'Definitive', 'Strong', 'Moderate', 'Limited', 'Disputed', 'Refuted', 'No Known Relationship', 'Provisional'}
)

# The bases that a code of a motif or of a repeat unit stands for. STRchive writes N in a motif, and a genotyper's
# catalogue any IUPAC code in a unit.
_BASES: Mapping[str, frozenset[str]] = {
    code: frozenset(bases)
    for code, bases in {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': 'AG',
        'Y': 'CT',
        'S': 'CG',
        'W': 'AT',
        'K': 'GT',
        'M': 'AC',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT',
    }.items()
}
# How a run writes the date of a release, which goes into a command line as the release's name does.
DATE = re.compile(r'[0-9]{4}-[0-9]{2}-[0-9]{2}')


@dataclasses.dataclass(frozen=True)
class Release:
    """A release of STRchive as a run gives it: the file of its loci, and its name and date as the run states them."""

    file: pathlib.Path
    name: str
    date: dt.date


def stated(file: pathlib.Path, name: str, date: str) -> Release:
    """The release that a run gives, from the name and the date as the run writes them.

    Raises:
        ValueError: if the name is not one that a command line takes as it is, or the date is not a day of the
            calendar written as YYYY-MM-DD.
    """
    if not pedigree.NAME.fullmatch(name):
        raise ValueError(f"the STRchive release {name!r} is no release's name, such as v2.28.0")
    if not DATE.fullmatch(date):
        raise ValueError(f'the date of the STRchive release is {date!r}, not a date as YYYY-MM-DD')
    try:
        day = dt.date.fromisoformat(date)
    except ValueError as no_day:
        raise ValueError(f'the date of the STRchive release is {date!r}, not a day of the calendar') from no_day
    return Release(file=file, name=name, date=day)


class StrchiveFile(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The STRchive file that a family's repeats were matched with, as the run read it."""

    file: str  # the file's name
    bytes: int
    sha256: str
    # The release and the day it was published, as the run stated them: the file holds neither.
    release: str
    release_date: dt.date
    # The rows of repeats with an entry of the file, and those without one.
    matched_repeats: int
    unmatched_repeats: int


@dataclasses.dataclass(frozen=True)
class Locus:
    """One entry of STRchive's file: where its repeat lies, what a repeat unit is compared with, and its values.

    `start` and `stop` are the entry's interval on GRCh38, 0-based and half-open as STRchive gives it. `motifs` are
    its reference motif and its pathogenic motifs in the reference's orientation, those as long as its first
    pathogenic motif: STRchive counts a size in units of that length, so a row whose unit is one of them has sizes in
    the unit of the entry's ranges. The length is not read from the entry's `motif_len`, which in every release up to
    v2.28.0 is another number for some locus. `columns` are the entry's values of FIELDS, all but SAME_INTERVAL, which
    is of a row and its entry together.
    """

    chrom: str
    start: int
    stop: int
    motifs: tuple[str, ...]
    columns: Mapping[str, Any]

    @property
    def id(self) -> str:
        return self.columns['strchive_id']


def _field(entry: Mapping[str, Any], key: str, where: str) -> Any:  # noqa: ANN401 (a JSON value)
    if key not in entry:
        raise ValueError(f'{where}: no {key}, which every locus of STRchive has')
    return entry[key]


def _text(entry: Mapping[str, Any], key: str, where: str) -> str:
    value = _field(entry, key, where)
    if not isinstance(value, str) or not value:
        raise ValueError(f'{where}: {key} is {value!r}, not a text')
    return value


def _whole(entry: Mapping[str, Any], key: str, where: str) -> int:
    """A whole number of the entry. The releases before v2.12.0 write a bound as 40.0, which JSON reads as 40 too."""
    value = _field(entry, key, where)
    # A JSON true is an int to Python.
    numeric = isinstance(value, int | float) and not isinstance(value, bool)
    if not numeric or not 0 <= value <= _INT_MAX or not float(value).is_integer():
        raise ValueError(f'{where}: {key} is {value!r}, not a whole number from 0 to {_INT_MAX}')
    return int(value)


def _texts(entry: Mapping[str, Any], key: str, where: str) -> tuple[str, ...]:
    value = _field(entry, key, where)
    if not isinstance(value, list) or not all(isinstance(item, str) and item for item in value):
        raise ValueError(f'{where}: {key} is {value!r}, not a list of texts')
    return tuple(value)


def _motifs(entry: Mapping[str, Any], key: str, where: str) -> tuple[str, ...]:
    motifs = _texts(entry, key, where)
    if unknown := sorted(motif for motif in motifs if not set(motif) <= _BASES.keys()):
        raise ValueError(f'{where}: {key} holds {unknown}, which are no motifs of bases and IUPAC codes')
    return motifs


def _evidence(entry: Mapping[str, Any], where: str) -> str | None:
    """The entry's level of evidence, which the file gives as a list of one level, or of none.

    The releases before v2.18.0 have no such field, and `read` refuses a file in which only some entries have it.
    """
    if 'evidence' not in entry:
        return None
    levels = _texts(entry, 'evidence', where)
    if len(levels) > 1 or not set(levels) <= LEVELS:
        raise ValueError(
            f'{where}: evidence is {list(levels)}, not one of the levels {", ".join(sorted(LEVELS))}, or no level'
        )
    return levels[0] if levels else None


def _locus(entry: Mapping[str, Any], where: str) -> Locus:
    locus_id = _text(entry, 'id', where)
    where = f'{where} ({locus_id})'
    chrom = _text(entry, 'chrom', where)
    if chrom not in sv_sites.PRIMARY:
        raise ValueError(f'{where}: chrom is {chrom!r}, not a primary chromosome named as chr1, chr2, ...')
    start, stop = _whole(entry, 'start_hg38', where), _whole(entry, 'stop_hg38', where)
    if start >= stop:
        raise ValueError(f'{where}: start_hg38 {start} is not before stop_hg38 {stop}')
    bounds: dict[str, int | None] = {
        bound: None if _field(entry, bound, where) is None else _whole(entry, bound, where) for bound in _BOUNDS
    }
    for kind in _RANGES:
        low, high = bounds[f'{kind}_min'], bounds[f'{kind}_max']
        if low is not None and high is not None and low > high:
            raise ValueError(f'{where}: {kind}_min {low} is above {kind}_max {high}')
    pathogenic = _motifs(entry, 'pathogenic_motif_reference_orientation', where)
    if not pathogenic:
        raise ValueError(f'{where}: no pathogenic motif, whose length is the unit of its ranges')
    reference = _motifs(entry, 'reference_motif_reference_orientation', where)
    inheritance = _texts(entry, 'inheritance', where)
    return Locus(
        chrom=chrom,
        start=start,
        stop=stop,
        motifs=tuple(m for m in dict.fromkeys((*reference, *pathogenic)) if len(m) == len(pathogenic[0])),
        columns={
            'strchive_id': locus_id,
            'strchive_gene': _text(entry, 'gene', where),
            'strchive_evidence': _evidence(entry, where),
            # A list without an entry is a null, as every list of the tables is.
            'strchive_inheritance': list(inheritance) or None,
            'strchive_start_hg38': start,
            'strchive_stop_hg38': stop,
            'strchive_pathogenic_motif_reference_orientation': list(pathogenic) or None,
            **{f'strchive_{bound}': value for bound, value in bounds.items()},
        },
    )


def read(path: pathlib.Path) -> tuple[Locus, ...]:
    """The loci of a STRchive file, in the file's order.

    Raises:
        ValueError: if the file is not JSON or not a list of objects with at least one, a locus lacks a field that is
            read, its id, gene or chromosome is no text or the chromosome no primary one, a coordinate or a bound of a
            range is no whole number, its start is not before its stop, the minimum of a range is above its maximum,
            its motifs or its inheritance are no list of texts, a motif holds something other than bases and IUPAC
            codes, it has no pathogenic motif, its evidence is not one of STRchive's levels or none, some loci have an
            evidence and others lack the field, or two loci share an id.
    """
    try:
        entries = json.loads(path.read_bytes())
    # A JSONDecodeError and a UnicodeDecodeError are both ValueErrors.
    except ValueError as unreadable:
        raise ValueError(f'{path}: not JSON: {unreadable}') from unreadable
    if not isinstance(entries, list) or not entries or not all(isinstance(entry, dict) for entry in entries):
        raise ValueError(f"{path}: not STRchive's loci, which are a list with an object for each locus")
    # A release has the field in every entry or, before v2.18.0, in none, and its column is then null in every row.
    without = [number for number, entry in enumerate(entries, start=1) if 'evidence' not in entry]
    if without and len(without) < len(entries):
        raise ValueError(f'{path}: locus {without[0]} has no evidence, which other loci of the file have')
    loci: dict[str, Locus] = {}
    for number, entry in enumerate(entries, start=1):
        locus = _locus(entry, f'{path}: locus {number}')
        if locus.id in loci:
            raise ValueError(f'{path}: locus {number} ({locus.id}): a second locus of this id')
        loci[locus.id] = locus
    return tuple(loci.values())


def same_repeat(unit: str, motif: str) -> bool:
    """Whether a repeat unit and a motif describe one repeat.

    They do where they have one length and, in one of the unit's rotations, each base of the unit can be the motif's
    base at its place: a repeat of CAG is one of AGC and of GCA too, and a catalogue may start it at any of the three.
    A code stands for its bases, so GCG can be GCN. Both are read on the reference's strand, as STRchive gives the
    motif and a repeat VCF the unit. A unit with a character that is no base and no code is of no motif.
    """
    if len(unit) != len(motif):
        return False
    return any(
        all(_BASES.get(base, frozenset()) & _BASES[other] for base, other in zip(rotated, motif, strict=True))
        for rotated in (unit[shift:] + unit[:shift] for shift in range(len(unit)))
    )


def _matches(row: Mapping[str, Any], locus: Locus) -> bool:
    """Whether a row of repeats is of a locus.

    A repeat VCF gives the base before the repeat as POS and the repeat's last base as END, so `pos` and `end_pos`
    are the repeat's interval as 0-based and half-open, which is how STRchive gives its own.
    """
    return (
        row['chrom'] == locus.chrom
        and row['pos'] < locus.stop
        and row['end_pos'] > locus.start
        and row['repeat_unit'] is not None
        and any(same_repeat(row['repeat_unit'], motif) for motif in locus.motifs)
    )


def annotate(rows: Sequence[Mapping[str, Any]], loci: Sequence[Locus]) -> tuple[list[dict[str, Any]], int]:
    """The rows of repeats, each with the values of its entry of STRchive, and the number of rows that have one.

    A row with an entry also says whether its interval is the entry's. A row without one has a null in each of the
    columns.

    Raises:
        ValueError: if a row matches two entries, which no release of STRchive that was read gives a repeat.
    """
    annotated = []
    matched = 0
    for row in rows:
        found = [locus for locus in loci if _matches(row, locus)]
        if len(found) > 1:
            raise ValueError(
                f'the repeat {row["repeat_id"]!r} at {row["chrom"]}:{row["pos"]}-{row["end_pos"]} matches the '
                f'STRchive loci {[locus.id for locus in found]}, and a repeat is of one locus'
            )
        matched += len(found)
        values: Mapping[str, Any] = dict.fromkeys(COLUMNS)
        if found:
            same = (row['pos'], row['end_pos']) == (found[0].start, found[0].stop)
            values = {**found[0].columns, SAME_INTERVAL: same}
        annotated.append({**row, **values})
    return annotated, matched
