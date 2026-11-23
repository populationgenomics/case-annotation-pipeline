"""The positions at which a family's tables say nothing of ClinVar (docs/annotation.md, Blinded positions).

A run can name a file of such positions for a family. The family's tables then hold a null in the ClinVar columns of
the rows in reach of an entry, as a row without a record holds one, so that a reader of the tables cannot tell a
blinded row from a row that never had a record. The union's tables are not changed.

The file is a JSON object that states its assembly and its coordinates and holds two lists of intervals:

- `small`, the reference bases of a small variant, from the first to the last base of the allele's REF as the tables
  hold the allele: an insertion's entry is its anchor base, and a deletion's begins with its anchor base. An entry
  empties the row of `variants` of every allele whose REF lies within REACH bases of it, and the row of `sv` of every
  record written as a sequence whose allele does, after the normalisation of the family's VCF.
- `structural`, the bases of a structural variant. An entry empties the rows of `sv` whose record it overlaps by the
  rule that the Annotator matches ClinVar's records by (sv_annotation.RECIPROCAL_OVERLAP), and no row of `variants`.

The functions of the two rules read views of the connection that family_tables.assemble made, and name them.
"""

from __future__ import annotations

import json
import pathlib
from collections.abc import Collection, Sequence
from typing import Any, Literal, Self

import duckdb
import pyarrow as pa
import pydantic

from case_annotation import annotation_tables, sv_annotation

# How many bases beside the REF of an allele an entry of a small variant still reaches. The tables write an insertion
# or a deletion with the base before it, and the Annotator lists ClinVar's records of other alleles "at the same
# position" without saying which base that is for such an allele. If it is the first base that the allele changes,
# the record of an insertion written at a position stands on a single-base change one base on, and that of a
# single-base change on an insertion written one base before it. One base on either side covers both readings.
REACH = 1

# The columns that an entry empties, by table: the ClinVar columns of both. An allele keeps its dbSNP ids. An allele
# that a population database holds nearly always has one, so a row without one would stand out among its
# neighbours, and an id gives a reader nothing that the allele's position and bases do not give.
VARIANTS_COLUMNS = annotation_tables.SOURCE_COLUMNS['clinvar']
SV_COLUMNS = sv_annotation.SOURCE_COLUMNS['clinvar']

type Kind = Literal['small', 'structural']

# The last position that a table can hold, whose `pos` and `end_pos` are 32-bit integers.
_LAST = 2**31 - 1


class Interval(pydantic.BaseModel, frozen=True, extra='forbid', strict=True):
    """The reference bases from `start` to `end`, counted from 1, on a chromosome named as the tables name it."""

    chrom: str
    start: int = pydantic.Field(ge=1, le=_LAST)
    end: int = pydantic.Field(le=_LAST)

    @pydantic.model_validator(mode='after')
    def _in_order(self) -> Self:
        if self.start > self.end:
            raise ValueError(f'an interval from {self.start} to {self.end}, which ends before it starts')
        return self


class BlindedPositions(pydantic.BaseModel, frozen=True, extra='forbid', strict=True):
    """A family's file of blinded positions. The two stated values are the only ones that the tools read."""

    assembly: Literal['GRCh38']
    coordinates: Literal['1-based, both ends inclusive']
    small: tuple[Interval, ...]
    structural: tuple[Interval, ...]

    def entries(self) -> list[tuple[Kind, int, Interval]]:
        """Each entry with its list's name and its place in the list, counted from 0."""
        kinds: tuple[tuple[Kind, tuple[Interval, ...]], ...] = (('small', self.small), ('structural', self.structural))
        return [(kind, place, interval) for kind, intervals in kinds for place, interval in enumerate(intervals)]


def _once(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """The object of a JSON document's keys and values; fails for a key that the object holds twice.

    A JSON reader keeps the last value of a repeated key, so a second `small` would replace the first list
    without a word, and the positions of the first would stay in the tables.
    """
    held = dict(pairs)
    if len(held) != len(pairs):
        keys = [key for key, _value in pairs]
        raise ValueError(f'the keys {sorted({key for key in keys if keys.count(key) > 1})} more than once')
    return held


def read(path: pathlib.Path, chromosomes: Collection[str]) -> BlindedPositions:
    """A family's blinded positions from their file.

    Raises:
        ValueError: if the file is not the object described above, with a key that it does not know or a key
            twice among the refused ones, or an entry lies on a chromosome that is none of those given.
    """
    text = path.read_bytes()
    try:
        json.loads(text, object_pairs_hook=_once)
        positions = BlindedPositions.model_validate_json(text)
    except (ValueError, RecursionError) as error:
        raise ValueError(f'{path}: not a file of blinded positions: {error}') from error
    for kind, place, interval in positions.entries():
        if interval.chrom not in chromosomes:
            raise ValueError(
                f'{path}: {kind}[{place}] is on {interval.chrom!r}, which is no chromosome of the tables, such as chr1'
            )
    return positions


def _register(con: duckdb.DuckDBPyConnection, positions: BlindedPositions) -> None:
    """Gives the connection the entries as the table `blinded_entries`, one row for each."""
    entries = positions.entries()
    con.register(
        'blinded_entries',
        pa.table(
            {
                'kind': pa.array([kind for kind, _place, _interval in entries], pa.string()),
                'place': pa.array([place for _kind, place, _interval in entries], pa.int32()),
                'chrom': pa.array([interval.chrom for _kind, _place, interval in entries], pa.string()),
                'first_base': pa.array([interval.start for _kind, _place, interval in entries], pa.int64()),
                'last_base': pa.array([interval.end for _kind, _place, interval in entries], pa.int64()),
            }
        ),
    )


# Whether an allele `a`, with its `chrom`, its `pos` and its `ref`, lies in reach of an entry `e` of a small variant.
_IN_REACH = (
    "e.kind = 'small' AND e.chrom = a.chrom "
    f'AND a.pos <= e.last_base + {REACH} AND a.pos + length(a.ref) - 1 >= e.first_base - {REACH}'
)


def _overlapping(pos: str, end_pos: str) -> str:
    """Whether a record whose bases run from pos + 1 to end_pos and an entry `e` overlap as the Annotator's matches do.

    The bases that the two share are at least RECIPROCAL_OVERLAP of the longer of the two. A record without a base,
    as a breakend or an insertion is, shares none.
    """
    shared = f'least({end_pos}, e.last_base) - greatest({pos} + 1, e.first_base) + 1'
    longer = f'greatest({end_pos} - {pos}, e.last_base - e.first_base + 1)'
    return f'({shared}) >= {sv_annotation.RECIPROCAL_OVERLAP} * ({longer})'


def _emptied(source: str, key: str, columns: Sequence[str], matches: str) -> str:
    """A query of the source's rows, with the columns null in the rows whose key is among the matches."""
    replaced = ', '.join(f'CASE WHEN b.{key} IS NULL THEN s.{column} END AS {column}' for column in columns)
    return (
        f'SELECT s.* REPLACE ({replaced}) FROM {source} s '  # noqa: S608 (names of this package's own)
        f'LEFT JOIN (SELECT DISTINCT {key} FROM {matches}) b USING ({key})'
    )


def blind_variants(con: duckdb.DuckDBPyConnection, positions: BlindedPositions | None) -> str:
    """The name of the view that a family reads the union's variants from: the union's own for one without positions.

    The connection has the family's `records` and the union's variants as `union_variants`. For a family with
    positions it gets `blinded_entries`, then `blinded_variants`, each allele of the family with each entry that
    reaches it, and a view that holds the union's rows with VARIANTS_COLUMNS null for those alleles.
    """
    if positions is None:
        return 'union_variants'
    _register(con, positions)
    con.execute(
        'CREATE TEMP TABLE blinded_variants AS SELECT a.variant_id, e.kind, e.place '  # noqa: S608 (a condition of this module's own)
        'FROM (SELECT variant_id, chrom, pos, ref FROM union_variants SEMI JOIN records USING (variant_id)) a '
        f'JOIN blinded_entries e ON {_IN_REACH}'
    )
    con.sql(_emptied('union_variants', 'variant_id', VARIANTS_COLUMNS, 'blinded_variants')).create_view(
        'blinded_union_variants'
    )
    return 'blinded_union_variants'


def blind_records(con: duckdb.DuckDBPyConnection, columns: Collection[str]) -> str:
    """Makes the view of the annotation of a family's records that one with blinded positions reads; returns its name.

    The connection has `blinded_entries`, the family's `sv_records`, their annotation as `sv_annotated`, and
    `sv_alleles`, the normalised allele of each record that is written as a sequence. It gets `blinded_sv`, each record
    with each entry that reaches it, and the view, which holds the annotation with those of SV_COLUMNS that are
    among `columns`, the ones the run did not leave out, null for those records.

    An entry of a structural variant is compared with the record's own interval and with the one that the Annotator
    annotated, which can begin elsewhere for a record written as a sequence. A record of two alleles that were
    annotated over two intervals has none, and is compared by its own alone: DuckDB's least and greatest skip a
    null, so a comparison with a missing interval would hold for an entry anywhere. The joins read a record's own
    columns only, since a member's column may bear any other name.
    """
    annotated = (
        'n.annotated_pos IS NOT NULL AND n.annotated_end_pos IS NOT NULL '
        f'AND {_overlapping("n.annotated_pos", "n.annotated_end_pos")}'
    )
    con.execute(
        'CREATE TEMP TABLE blinded_sv AS '  # noqa: S608 (conditions of this module's own)
        'SELECT r.sv_id, e.kind, e.place FROM (SELECT sv_id, chrom, pos, end_pos FROM sv_records) r '
        'JOIN (SELECT sv_id, annotated_pos, annotated_end_pos FROM sv_annotated) n USING (sv_id) '
        "JOIN blinded_entries e ON e.kind = 'structural' AND e.chrom = r.chrom "
        f'AND ({_overlapping("r.pos", "r.end_pos")} OR ({annotated})) '
        f'UNION ALL SELECT a.sv_id, e.kind, e.place FROM sv_alleles a JOIN blinded_entries e ON {_IN_REACH}'
    )
    emptied = [column for column in SV_COLUMNS if column in columns]
    if not emptied:
        return 'sv_annotated'
    con.sql(_emptied('sv_annotated', 'sv_id', emptied, 'blinded_sv')).create_view('blinded_sv_annotated')
    return 'blinded_sv_annotated'


def log(con: duckdb.DuckDBPyConnection, positions: BlindedPositions | None, *, with_records: bool) -> None:
    """Prints one line for each entry of a family with positions, with how many rows of each table the entry reached.

    The lines are for the log of the step that writes the family's tables, which is no file of the family's folder.
    They name an entry by its list and its place in it and not by its position. `with_records` says whether the
    family has an `sv` table, and so whether the connection has `blinded_sv`.
    """
    if positions is None:
        return
    counted = {'variants': 'blinded_variants', **({'sv': 'blinded_sv'} if with_records else {})}
    rows: dict[tuple[str, int], dict[str, int]] = {}
    for table, matches in counted.items():
        query = (
            f'SELECT e.kind, e.place, count(m.place) FROM blinded_entries e '  # noqa: S608 (a table of this module's own)
            f'LEFT JOIN {matches} m USING (kind, place) GROUP BY ALL'
        )
        for kind, place, reached in con.sql(query).fetchall():
            rows.setdefault((kind, place), {})[table] = reached
    for (kind, place), reached in sorted(rows.items()):
        # An entry of a structural variant reaches no row of variants by its rule.
        tables = [table for table in reached if kind == 'small' or table != 'variants']
        counts = ', '.join(f'{table} {reached[table]}' for table in tables)
        print(f'blinded positions: {kind}[{place}]: {counts or "no table that it could reach"}')
