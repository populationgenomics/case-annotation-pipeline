"""Values of a sites union that are joined to its alleles from public files, as Parquet tables (docs/annotation.md).

Some values of an allele come from a public file that is neither a source of Illumina's bundle nor an echtvar archive:
a VCF of sites, an annotation of intervals, a table of scores, or the catalogue of branch points that this package
prepares from a paper's datasets (branch_points). This module joins each such file to the union's alleles with DuckDB,
by allele or by position, and writes one row per allele, keyed by variant id as `variants` is, with a record of each
file beside it. Some files give a value for an allele and one transcript. Those values go into a second table, keyed by
variant id and by Ensembl's transcript id without its version, with a row only for an allele and a transcript that has
a value.

A Source describes one file: its columns and the function that reads it. A run gives any of the files and none is
needed, so the tables of a run hold the columns of the files that it was given. A reader checks the file as it reads
it and fails for one that is not what its source publishes: a header or an INFO definition that differs, a row that
does not parse, or two rows that would give an allele two values. A further file is one more reader and one more
Source.

A source that is a VCF of alleles is read as `bcftools norm` rewrote it with the options that normalise a family's VCF,
so that its indels are split and left-aligned as the union's are. The record of the file is of the file as the run was
given it.

An allele is matched to an interval by the reference bases that it changes (changed_span). An insertion changes none
and lies between two bases, so it overlaps an interval that holds both of them, and none that it only touches.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import itertools
import math
import pathlib
import tempfile
from collections.abc import Callable, Iterator, Mapping, Sequence
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam

from case_annotation import annotation_tables, branch_points, sv_sites

TABLE_NAME = 'joined_values.parquet'
TRANSCRIPTS_NAME = 'joined_transcript_values.parquet'
INFO_NAME = 'joined_values.json'
_BATCH_ROWS = 100_000

_STRING = pa.string()
_INT = pa.int32()
_FLOAT = pa.float32()
# What a row of each table is of. A transcript's id is Ensembl's without its version, as the files give it.
_ALLELE_KEY = (pa.field('variant_id', _STRING, nullable=False),)
_TRANSCRIPT_KEY = (*_ALLELE_KEY, pa.field('transcript_id', _STRING, nullable=False))

# The union's alleles as the readers meet them, as the view `alleles`. `record` numbers them in the union's order,
# which is the order of the table's rows. changed_start and changed_end are the reference bases that the allele
# changes, as a half-open interval of 0-based positions (changed_span).
_ALLELES = pa.schema(
    [
        pa.field('record', pa.int64(), nullable=False),
        pa.field('variant_id', _STRING, nullable=False),
        pa.field('chrom', _STRING, nullable=False),
        pa.field('pos', _INT, nullable=False),
        pa.field('ref', _STRING, nullable=False),
        pa.field('alt', _STRING, nullable=False),
        pa.field('changed_start', _INT, nullable=False),
        pa.field('changed_end', _INT, nullable=False),
    ]
)
# An overlap is looked for among the intervals that share a stretch of this many bases with the allele, so that the
# join is on a chromosome and a stretch and compares an allele with the few intervals near it.
_STRETCH = 1024
_BASES = "('A', 'C', 'G', 'T')"
_PRIMARY = '(' + ', '.join(f"'{chrom}'" for chrom in sv_sites.PRIMARY) + ')'
_ENSEMBL_TRANSCRIPT = 'ENST[0-9]+'


class JoinedFile(pydantic.BaseModel, frozen=True, extra='forbid'):
    """A file that a run joined to the union's alleles, as the run read it."""

    source: str  # the name of its Source, which is the name of the param that gave it
    file: str  # the file's name
    bytes: int
    sha256: str


class JoinedFiles(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The files that a run joined, in the order of their columns in the tables; none in a run that was given none."""

    files: tuple[JoinedFile, ...]

    def given(self) -> tuple[Source, ...]:
        """The sources that the run read a file of, in the tables' order."""
        return tuple(SOURCES[file.source] for file in self.files)


@dataclasses.dataclass(frozen=True)
class Values:
    """A file's values of the union's alleles, as queries on the connection that its reader read it into.

    `alleles` gives `variant_id` and the source's columns of `variants` for each allele that has a value, each allele
    at most once. `transcripts` gives `variant_id`, `transcript_id` and the source's columns of `transcripts` for each
    allele and transcript that has a value, each pair at most once; a source without such columns has none.
    """

    alleles: str
    transcripts: str | None = None


@dataclasses.dataclass(frozen=True)
class Source:
    """What the tables take from one public file.

    `name` is the file's name wherever a run meets it: the workflow's param, the `source` of its record in a family's
    manifest, and the source of its columns in the contract.

    `fields` are its columns of `variants`, and `transcript_fields` its columns of `transcripts`, which a file has
    that gives a value for an allele and one transcript.

    `read` takes a connection that has the union's alleles as the view `alleles` (_ALLELES) and the file. It checks the
    file, and returns the queries of its values.
    """

    name: str
    fields: tuple[pa.Field, ...]
    read: Callable[[duckdb.DuckDBPyConnection, pathlib.Path], Values]
    transcript_fields: tuple[pa.Field, ...] = ()
    # Whether the file is a VCF of alleles, which is read as bcftools norm rewrote it.
    normalised: bool = False
    # Whether the columns are flags, which are false and not null for an allele without a row.
    flags: bool = False

    @property
    def columns(self) -> tuple[str, ...]:
        return tuple(field.name for field in self.fields)

    @property
    def transcript_columns(self) -> tuple[str, ...]:
        return tuple(field.name for field in self.transcript_fields)


def changed_span(pos: int, ref: str, alt: str) -> tuple[int, int]:
    """The reference bases that an allele changes, as a half-open interval of 0-based positions.

    The bases that REF and ALT start with alike are not changed, which leaves out the base that a VCF writes before an
    insertion or a deletion. An insertion changes no base: its interval is empty and lies between the two bases that
    the inserted sequence goes between.
    """
    shared = sum(1 for _ in itertools.takewhile(lambda bases: bases[0] == bases[1], zip(ref, alt, strict=False)))
    return pos - 1 + shared, pos - 1 + len(ref)


def _overlaps(intervals: str, selected: str) -> str:
    """A query of what is selected of each allele of the union, as `a`, and each interval of a table that it overlaps.

    The table holds `chrom` and the half-open interval `start0` to `end0` of 0-based positions, and the interval is
    `i` in what is selected. An allele overlaps an interval that holds a base it changes. An insertion overlaps an
    interval that holds the base on each side of it, which the same comparison gives for its empty interval, and none
    that begins or ends where it lies. An allele and an interval that share several stretches meet once in each, so
    the query gives what is selected once, however often it is found.
    """
    return (
        f'SELECT DISTINCT {selected} FROM ('  # noqa: S608 (a table name and column names of this module's own)
        'SELECT variant_id, chrom, changed_start, changed_end, '
        f'unnest(range(changed_start // {_STRETCH}, greatest(changed_end - 1, changed_start) // {_STRETCH} + 1)) '
        'AS stretch FROM alleles) a '
        f'JOIN (SELECT *, unnest(range(start0 // {_STRETCH}, (end0 - 1) // {_STRETCH} + 1)) '
        f'AS stretch FROM {intervals}) i '
        'ON a.chrom = i.chrom AND a.stretch = i.stretch AND a.changed_start < i.end0 AND a.changed_end > i.start0'
    )


def _overlapping(intervals: str) -> str:
    """A query of the variant ids of the union's alleles that overlap an interval of a table (_overlaps)."""
    return _overlaps(intervals, 'a.variant_id')


def _gzipped(path: pathlib.Path) -> bool:
    with path.open('rb') as f:
        return f.read(2) == b'\x1f\x8b'


def _text_rows(
    path: pathlib.Path, columns: Sequence[str], *, header: bool, comment: str | None = None, skip: int = 0
) -> str:
    """A table function that reads a tab-separated file, gzipped or not, with every field as text.

    No field is quoted, and a row with another number of fields than the columns fails when the rows are read. The
    file's first `skip` lines are not read.
    """
    declared = ', '.join(f"'{name}': 'VARCHAR'" for name in columns)
    options = [
        "delim = '\\t'",
        f'header = {str(header).lower()}',
        "quote = ''",
        "escape = ''",
        'auto_detect = false',
        'strict_mode = true',
        f"compression = '{'gzip' if _gzipped(path) else 'none'}'",
        f'columns = {{{declared}}}',
    ]
    if comment is not None:
        options.append(f"comment = '{comment}'")
    if skip:
        options.append(f'skip = {skip}')
    quoted = str(path).replace("'", "''")
    return f"read_csv('{quoted}', {', '.join(options)})"


def _create(con: duckdb.DuckDBPyConnection, table: str, query: str, path: pathlib.Path, what: str) -> None:
    """Reads a file's rows into a temporary table.

    Raises:
        ValueError: if DuckDB cannot read the file as the query describes it, as for a row with too few fields or a
            gzip stream that is cut off.
    """
    try:
        con.execute(f'CREATE TEMP TABLE {table} AS {query}')
    except duckdb.Error as unreadable:
        raise ValueError(f'{path}: not {what}: {unreadable}') from unreadable


def _whole(text: str, least: int = 1) -> str:
    """A whole number from a text column that is at least the least one, or null for any other text."""
    number = f'TRY_CAST({text} AS INTEGER)'
    return f"CASE WHEN regexp_full_match({text}, '[0-9]+') AND {number} >= {least} THEN {number} END"


def _finite(text: str) -> str:
    """A finite number from a text column, or null for any other text."""
    number = f'TRY_CAST({text} AS DOUBLE)'
    return f'CASE WHEN isfinite({number}) THEN {number} END'


@dataclasses.dataclass(frozen=True)
class _Layout:
    """A tab-separated table with a header line, as its source publishes it."""

    what: str  # what a message calls the file, such as "APARENT2's table of effects"
    whose: str  # whose header it is, such as "APARENT2's"
    header: tuple[str, ...]
    # What each line before the header starts with, in a file that says of itself what it was made from.
    preamble: str | None = None


def _read_rows(
    con: duckdb.DuckDBPyConnection,
    path: pathlib.Path,
    layout: _Layout,
    table: str,
    *,
    parsed: Mapping[str, str],
    malformed: str,
    shown: Sequence[str],
    required: str,
    keep: str | None = None,
) -> None:
    """Reads the rows of a file of a layout into a temporary table, and fails for a file that is not of the layout.

    The file is read once. `parsed` gives each column of the table as an expression on the file's columns, which are
    text. `malformed` is the condition on the table's columns under which a row is refused, and a row for which it is
    null is refused too. The table holds every row, or with `keep` the rows that meet that condition, so that the
    whole file is checked while only those rows are held. A refusal shows the `shown` columns of some refused rows as
    the file writes them, and says that they are not what is `required`. The lines before the header that start with
    the layout's preamble, where it has one, are passed over.

    Raises:
        ValueError: if the header is not the layout's, a row has another number of fields, a gzip stream is cut off,
            or a row is malformed.
    """
    before = 0
    try:
        with gzip.open(path, 'rb') if _gzipped(path) else path.open('rb') as f:
            header = f.readline().decode('utf-8', 'replace').rstrip('\r\n')
            while layout.preamble is not None and header.startswith(layout.preamble):
                before += 1
                header = f.readline().decode('utf-8', 'replace').rstrip('\r\n')
    # A gzip stream that ends within the header raises an EOFError, which is no OSError.
    except (OSError, EOFError) as unreadable:
        raise ValueError(f'{path}: not {layout.what}: {unreadable}') from unreadable
    if header != '\t'.join(layout.header):
        raise ValueError(f'{path}: its header is {header.split(chr(9))}, and {layout.whose} is {list(layout.header)}')
    columns = ', '.join(parsed)
    selected = ', '.join(f'{expression} AS {column}' for column, expression in parsed.items())
    # `text` holds a malformed row's fields as the file writes them, for the message.
    text = ', '.join(f'"{name}"' for name in shown)
    _create(
        con,
        table,
        'WITH parsed AS ('  # noqa: S608 (a file's path, quoted, and expressions of this module)
        f'SELECT {selected}, [{text}] AS text FROM {_text_rows(path, layout.header, header=True, skip=before)}'
        f'), checked AS (SELECT {columns}, text, coalesce({malformed}, true) AS malformed FROM parsed) '
        f'SELECT {columns}, malformed, CASE WHEN malformed THEN text END AS text FROM checked'
        + ('' if keep is None else f' WHERE malformed OR ({keep})'),
        path,
        layout.what,
    )
    if refused := con.sql(f'SELECT text FROM {table} WHERE malformed LIMIT 5').fetchall():  # noqa: S608 (a table of this module)
        raise ValueError(f'{path}: rows without {", ".join(shown)} as {required}, such as {[r[0] for r in refused]}')


def _unique(con: duckdb.DuckDBPyConnection, path: pathlib.Path, table: str, key: Sequence[str], what: str) -> None:
    """Fails if two rows of a table that a file was read into share a key, which would give an allele two values."""
    if repeated := con.sql(
        f'SELECT {", ".join(key)} FROM {table} GROUP BY ALL HAVING count(*) > 1 ORDER BY ALL LIMIT 5'  # noqa: S608 (names of this module)
    ).fetchall():
        raise ValueError(f'{path}: several rows for one {what}, such as {repeated}')


# gnomAD's sites VCF of the mitochondrial genome (release 3.1). A homoplasmic call is one at a heteroplasmy level of
# 0.95 or more, and a heteroplasmic one at 0.10 up to 0.95, as the file's header defines its counts.
_MITO_CHROM = 'chrM'
# The INFO fields that the columns come from, each with the type and the number of values that the header must declare,
# which are then the type and the number of what htslib returns for a record. htslib returns no value for one that is
# not of the declared type. The lists by haplogroup are one string, with `|` between the values.
_MITO_INFO: Mapping[str, tuple[str, int]] = {
    'AN': ('Integer', 1),
    'AC_hom': ('Integer', 1),
    'AC_het': ('Integer', 1),
    'AF_hom': ('Float', 1),
    'AF_het': ('Float', 1),
    'max_hl': ('Float', 1),
    'hap_AF_hom': ('String', 1),
    'hap_AF_het': ('String', 1),
    'hapmax_AF_hom': ('String', 1),
    'hapmax_AF_het': ('String', 1),
    'faf_hapmax_hom': ('Float', 1),
    'hap_defining_variant': ('Flag', 0),
}
# The columns that hold an INFO field as the file gives it, and whether every record has the field.
_MITO_COPIED: Mapping[str, tuple[str, bool]] = {
    'gnomad_mito_ac_hom': ('AC_hom', True),
    'gnomad_mito_ac_het': ('AC_het', True),
    'gnomad_mito_af_hom': ('AF_hom', True),
    'gnomad_mito_af_het': ('AF_het', True),
    'gnomad_mito_an': ('AN', True),
    'gnomad_mito_max_heteroplasmy': ('max_hl', True),
    'gnomad_mito_hapmax_hom_haplogroup': ('hapmax_AF_hom', False),
    'gnomad_mito_hapmax_het_haplogroup': ('hapmax_AF_het', False),
    'gnomad_mito_faf_hapmax_hom': ('faf_hapmax_hom', True),
}
# The columns that hold the highest of a list by haplogroup.
_MITO_HIGHEST: Mapping[str, str] = {
    'gnomad_mito_hapmax_af_hom': 'hap_AF_hom',
    'gnomad_mito_hapmax_af_het': 'hap_AF_het',
}
_MITO_FIELDS = (
    pa.field('gnomad_mito_ac_hom', _INT),
    pa.field('gnomad_mito_ac_het', _INT),
    pa.field('gnomad_mito_af_hom', _FLOAT),
    pa.field('gnomad_mito_af_het', _FLOAT),
    pa.field('gnomad_mito_an', _INT),
    pa.field('gnomad_mito_max_heteroplasmy', _FLOAT),
    pa.field('gnomad_mito_hapmax_af_hom', _FLOAT),
    pa.field('gnomad_mito_hapmax_af_het', _FLOAT),
    pa.field('gnomad_mito_hapmax_hom_haplogroup', _STRING),
    pa.field('gnomad_mito_hapmax_het_haplogroup', _STRING),
    pa.field('gnomad_mito_faf_hapmax_hom', _FLOAT),
    pa.field('gnomad_mito_filters', pa.list_(_STRING)),
    pa.field('gnomad_mito_hap_defining', pa.bool_()),
)


def _records(vcf: pysam.VariantFile, path: pathlib.Path) -> Iterator[pysam.VariantRecord]:
    """The records of a VCF; a line that htslib cannot read as a record fails under the file's name."""
    try:
        yield from vcf
    except OSError as unreadable:
        raise ValueError(f'{path}: a record that does not parse: {unreadable}') from unreadable


def _highest(listed: object, field: str, at: str) -> float | None:
    """The highest of the frequencies that gnomAD lists by haplogroup; none where no haplogroup has one.

    gnomAD writes NaN for a haplogroup without a sample that has a call at the site.
    """
    try:
        values = [float(value) for value in str(listed).split('|')]
    except ValueError as malformed:
        raise ValueError(f'{at}: {field} is {listed!r}, not numbers with | between them') from malformed
    return max((value for value in values if not math.isnan(value)), default=None)


def _gnomad_mito(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> Values:
    """Reads gnomAD's mitochondrial sites, normalised as the union is, and matches them by allele.

    Raises:
        ValueError: if the header lacks an INFO field that the columns come from or declares it otherwise, the VCF has
            samples, a line is no record, a record is not on chrM, does not hold one alternate allele, or lacks a
            field that every record of gnomAD's file has or holds no value of the field's type in it, a list by
            haplogroup holds something other than numbers, or two records hold one allele.
    """
    declared = pa.schema([*_ALLELE_KEY, *_MITO_FIELDS])
    rows: dict[str, dict[str, Any]] = {}
    try:
        vcf = pysam.VariantFile(str(path))
    except (OSError, ValueError) as unreadable:
        raise ValueError(f'{path}: not a VCF: {unreadable}') from unreadable
    with vcf:
        if len(vcf.header.samples):
            raise ValueError(f"{path}: gnomAD's sites have no samples, and this VCF has some")
        for name, (kind, number) in _MITO_INFO.items():
            if name not in vcf.header.info:
                raise ValueError(f"{path}: no INFO field {name}, which gnomAD's mitochondrial sites VCF defines")
            found = (vcf.header.info[name].type, vcf.header.info[name].number)
            if found != (kind, number):
                raise ValueError(f'{path}: INFO {name} is declared as {found}, not as {(kind, number)}')
        for record in _records(vcf, path):
            at = f'{path}: {record.chrom}:{record.pos}'
            if record.chrom != _MITO_CHROM:
                raise ValueError(f'{at}: not on {_MITO_CHROM}, which is all that the file holds')
            vid = annotation_tables.record_id(record, str(path))
            if vid in rows:
                raise ValueError(f'{at}: a second record of the allele {vid}')
            row: dict[str, Any] = {'variant_id': vid}
            for column, (name, required) in _MITO_COPIED.items():
                value: object = record.info.get(name)
                if value is None and required:
                    raise ValueError(f"{at}: no {name}, which every record of gnomAD's file has")
                row[column] = value
            for column, name in _MITO_HIGHEST.items():
                if record.info.get(name) is None:
                    raise ValueError(f"{at}: no {name}, which every record of gnomAD's file has")
                row[column] = _highest(record.info.get(name), name, at)
            row['gnomad_mito_filters'] = list(record.filter.keys())
            row['gnomad_mito_hap_defining'] = 'hap_defining_variant' in record.info
            rows[vid] = row
    con.register('gnomad_mito_rows', pa.Table.from_pylist(list(rows.values()), schema=declared))
    return Values('SELECT * FROM gnomad_mito_rows')


# GENCODE's polyadenylation features, a GTF. A polyA_signal is the signal before a site where a transcript is cut and
# polyadenylated, nearly always six bases. The file also holds the sites themselves and signals of pseudogenes.
_GTF_COLUMNS = ('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes')
_POLYA_FEATURES = ('polyA_signal', 'polyA_site', 'pseudo_polyA')
_POLYA_SIGNAL = 'polyA_signal'


def _gencode_polya(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> Values:
    """Reads GENCODE's polyadenylation features and finds the alleles that overlap a polyA_signal.

    Two features can hold the same bases, and the file has such pairs. They give no allele a second row, because the
    result is a flag.

    Raises:
        ValueError: if a line is not the nine fields of a GTF, a feature is none of the three kinds that the file
            holds, a feature's start or end is no position or its end lies before its start, or the file has no
            polyA_signal on a primary chromosome.
    """
    rows = _text_rows(path, _GTF_COLUMNS, header=False, comment='#')
    _create(
        con,
        'gencode_polya_rows',
        f'SELECT seqname AS chrom, feature, start AS start_text, "end" AS end_text, '  # noqa: S608 (a file's path, quoted)
        f'{_whole("start")} AS first, {_whole('"end"')} AS last FROM {rows}',
        path,
        "GENCODE's polyadenylation features, a GTF",
    )
    kinds = ', '.join(f"'{kind}'" for kind in _POLYA_FEATURES)
    if other := con.sql(
        f'SELECT DISTINCT feature FROM gencode_polya_rows WHERE feature IS NULL OR feature NOT IN ({kinds}) '  # noqa: S608 (constants of this module)
        'ORDER BY 1 LIMIT 5'
    ).fetchall():
        raise ValueError(
            f"{path}: the features {[r[0] for r in other]}, and GENCODE's polyadenylation features are {kinds}"
        )
    if malformed := con.sql(
        'SELECT chrom, feature, start_text, end_text FROM gencode_polya_rows '
        'WHERE chrom IS NULL OR first IS NULL OR last IS NULL OR last < first LIMIT 5'
    ).fetchall():
        raise ValueError(f'{path}: features without a chromosome, a start and an end after it, such as {malformed}')
    primary = ', '.join(f"'{chrom}'" for chrom in sv_sites.PRIMARY)
    con.execute(
        'CREATE TEMP TABLE gencode_polya_signals AS SELECT chrom, first - 1 AS start0, last AS end0 '  # noqa: S608 (constants of this module)
        f"FROM gencode_polya_rows WHERE feature = '{_POLYA_SIGNAL}' AND chrom IN ({primary})"
    )
    if not con.sql('SELECT 1 FROM gencode_polya_signals LIMIT 1').fetchall():
        raise ValueError(
            f'{path}: no {_POLYA_SIGNAL} on a primary chromosome; are its chromosomes named chr1, chr2, ...?'
        )
    return Values(f'SELECT variant_id, true AS gencode_polya_signal FROM ({_overlapping("gencode_polya_signals")})')  # noqa: S608 (a query of this module)


# APARENT2's precomputed effects of the single-base changes near the polyadenylation sites of PolyA_DB: one row for a
# change and a site. A change near two sites, or near sites of genes on both strands, has several rows. The positions
# are 1-based and the bases are of the plus strand, also for a gene on the minus strand.
_APARENT2 = _Layout(
    what="APARENT2's table of effects",
    whose="APARENT2's",
    header=(
        'chrom',
        'strand',
        'gene',
        'var_position_hg19',
        'var_position_hg38',
        'ref_nucleotide',
        'var_nucleotide',
        'pas_type',
        'pas_position_hg19',
        'pas_position_hg38',
        'ref_usage',
        'delta_logodds_narrow',
        'delta_usage_narrow',
        'delta_logodds',
        'delta_usage',
    ),
)


def _aparent2(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> Values:
    """Reads APARENT2's effects and gives each single-base change of the union the summary of its rows.

    A row is kept where it is of an allele of the union. The lowest and the highest change of the log odds are over an
    allele's rows, and the reference usage is of the site of the lowest row: where several rows are the lowest, of the
    most used of their sites.

    Raises:
        ValueError: if the header is not APARENT2's, a row has another number of fields, or a row lacks a chromosome,
            a GRCh38 position, a base as its reference or its changed nucleotide, or a finite number as its reference
            usage or its change of the log odds.
    """
    _read_rows(
        con,
        path,
        _APARENT2,
        'aparent2_rows',
        parsed={
            'chrom': 'chrom',
            'pos': _whole('var_position_hg38'),
            'ref': 'ref_nucleotide',
            'alt': 'var_nucleotide',
            'delta': _finite('delta_logodds_narrow'),
            'usage': _finite('ref_usage'),
        },
        malformed=(
            f'chrom IS NULL OR pos IS NULL OR ref NOT IN {_BASES} OR alt NOT IN {_BASES} '
            'OR delta IS NULL OR usage IS NULL'
        ),
        shown=('chrom', 'var_position_hg38', 'ref_nucleotide', 'var_nucleotide', 'ref_usage', 'delta_logodds_narrow'),
        required='a chromosome, a position, two bases and two finite numbers',
        keep=(
            '(chrom, pos, ref, alt) IN (SELECT (chrom, pos, ref, alt) FROM alleles WHERE len(ref) = 1 AND len(alt) = 1)'
        ),
    )
    return Values(
        'SELECT a.variant_id, s.aparent2_min, s.aparent2_max, s.aparent2_ref_usage FROM ('
        'SELECT chrom, pos, ref, alt, min(delta) AS aparent2_min, max(delta) AS aparent2_max, '
        'arg_min(usage, (delta, -usage)) AS aparent2_ref_usage FROM aparent2_rows GROUP BY ALL'
        ') s JOIN alleles a USING (chrom, pos, ref, alt)'
    )


# MitoTIP's scores of every change of a base of a mitochondrial tRNA gene: one row for a position and a change, which
# is one of the three other bases or, written as `:`, the deletion of the base. The positions and the reference bases
# are of the revised Cambridge Reference Sequence, which is chrM of GRCh38. A higher score is a change that MitoTIP
# takes as more likely to be pathogenic, and Q1 is the quartile of the highest scores.
_MITOTIP = _Layout(
    what="MitoTIP's table of scores",
    whose="MitoTIP's",
    header=('Position', 'rCRS', 'Alt', 'MitoTIP_Score', 'Quartile', 'Count', 'Percentage', 'Mitomap_Status'),
)
_MITOTIP_DELETION = ':'
_MITOTIP_QUARTILES = "('Q1', 'Q2', 'Q3', 'Q4')"


def _mitotip(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> Values:
    """Reads MitoTIP's scores and matches the single-base changes by allele.

    The row of a deletion is checked and not joined. MitoTIP scores the deletion of each base of a run of equal bases
    by itself, and a VCF writes all of them as one allele, so such an allele has no row of its own.

    Raises:
        ValueError: if the header is not MitoTIP's, a row has another number of fields or lacks a position, a base as
            its reference, another base or `:` as its change, a finite number as its score or one of Q1 to Q4 as its
            quartile, or two rows hold one change of one position.
    """
    _read_rows(
        con,
        path,
        _MITOTIP,
        'mitotip_rows',
        parsed={
            'pos': _whole('"Position"'),
            'ref': '"rCRS"',
            'alt': '"Alt"',
            'score': _finite('"MitoTIP_Score"'),
            'quartile': '"Quartile"',
        },
        malformed=(
            f"pos IS NULL OR ref NOT IN {_BASES} OR (alt NOT IN {_BASES} AND alt != '{_MITOTIP_DELETION}') "
            f'OR alt = ref OR score IS NULL OR quartile NOT IN {_MITOTIP_QUARTILES}'
        ),
        shown=('Position', 'rCRS', 'Alt', 'MitoTIP_Score', 'Quartile'),
        required=f'a position, a base, another base or {_MITOTIP_DELETION}, a finite number and one of Q1 to Q4',
    )
    _unique(con, path, 'mitotip_rows', ('pos', 'ref', 'alt'), 'change of a position')
    return Values(
        'SELECT a.variant_id, m.score AS mitotip_score, m.quartile AS mitotip_quartile FROM mitotip_rows m '  # noqa: S608 (a constant of this module)
        f"JOIN alleles a ON a.chrom = '{_MITO_CHROM}' AND a.pos = m.pos AND a.ref = m.ref AND a.alt = m.alt"
    )


# gnomAD's constraint of the genome in windows of 1,000 bases, from gnomAD v3.1.2: one row for a window of an autosome
# that passed gnomAD's quality control, with the z score of its depletion of variants. A window's start and end are
# 0-based and half-open and its start is a multiple of 1,000, so the window of the 1-based position p starts at the
# highest multiple of 1,000 that is at most p - 1.
_WINDOW = 1000
_WINDOWS = _Layout(
    what="gnomAD's genomic constraint by window",
    whose="gnomAD's",
    header=('chrom', 'start', 'end', 'element_id', 'possible', 'expected', 'observed', 'oe', 'z'),
)


def _genomic_constraint(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> Values:
    """Reads gnomAD's windows and gives each allele the z score of the window that holds its position.

    Raises:
        ValueError: if the header is not gnomAD's, a row has another number of fields or lacks a primary chromosome,
            a start that is a multiple of 1,000, an end after it and at most 1,000 bases on, or a finite number as its
            z score, or two rows hold one window.
    """
    _read_rows(
        con,
        path,
        _WINDOWS,
        'genomic_constraint_rows',
        parsed={
            'chrom': 'chrom',
            'start0': _whole('"start"', least=0),
            'end0': _whole('"end"'),
            'z': _finite('z'),
        },
        malformed=(
            f'chrom NOT IN {_PRIMARY} OR start0 IS NULL OR start0 % {_WINDOW} != 0 OR end0 IS NULL '
            f'OR end0 <= start0 OR end0 > start0 + {_WINDOW} OR z IS NULL'
        ),
        shown=('chrom', 'start', 'end', 'z'),
        required=(
            f'a primary chromosome, a start that is a multiple of {_WINDOW}, an end after it and at most {_WINDOW} '
            'bases on, and a finite number'
        ),
    )
    _unique(con, path, 'genomic_constraint_rows', ('chrom', 'start0'), 'window')
    return Values(
        'SELECT a.variant_id, w.z AS genomic_constraint_z FROM alleles a JOIN genomic_constraint_rows w '  # noqa: S608 (a constant of this module)
        f'ON a.chrom = w.chrom AND ((a.pos - 1) // {_WINDOW}) * {_WINDOW} = w.start0 AND a.pos <= w.end0'
    )


# gnomAD's regions of missense constraint (v4.1.1): one row for a region of a transcript's coding sequence, with the
# ratio of the missense variants that gnomAD observed in it to those that it expected. A region's start and end are
# 1-based and closed, and the region spans the introns between its coding bases. The regions of one transcript do not
# overlap and those of two transcripts do, so a region is matched by its transcript as well as by position. A
# transcript with a single region has its whole coding sequence as that region.
_MCR = _Layout(
    what="gnomAD's regions of missense constraint",
    whose="gnomAD's",
    header=(
        'transcript',
        'n_total_transcript_regions',
        'chr',
        'region_start',
        'region_end',
        'region_obs',
        'region_exp',
        'region_oe',
        'region_oe_chisq',
        'region_oe_chisq_p',
        'is_high_coverage',
    ),
)


def _mcr(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> Values:
    """Reads gnomAD's regions and gives each allele the region of each transcript that holds its position.

    An allele is matched by its position alone, an insertion or a deletion by the position that the VCF writes it at.
    The value of an allele is the lowest ratio among its regions.

    Raises:
        ValueError: if the header is not gnomAD's, a row has another number of fields or lacks an Ensembl transcript
            id without a version, a number of regions, a primary chromosome, a start and an end as positions in that
            order, a finite ratio that is not below 0, a p-value from 0 to 1 or `true` or `false` as its coverage,
            two regions of one transcript overlap or lie on two chromosomes, or a transcript has another number of
            regions than its rows say.
    """
    _read_rows(
        con,
        path,
        _MCR,
        'mcr_rows',
        parsed={
            'transcript_id': 'transcript',
            'regions': _whole('n_total_transcript_regions'),
            'chrom': '"chr"',
            'first': _whole('region_start'),
            'last': _whole('region_end'),
            'oe': _finite('region_oe'),
            'p': _finite('region_oe_chisq_p'),
            'high_coverage': "CASE is_high_coverage WHEN 'true' THEN true WHEN 'false' THEN false END",
        },
        malformed=(
            f"NOT regexp_full_match(transcript_id, '{_ENSEMBL_TRANSCRIPT}') OR regions IS NULL "
            f'OR chrom NOT IN {_PRIMARY} OR first IS NULL OR last IS NULL OR last < first '
            'OR oe IS NULL OR oe < 0 OR p IS NULL OR p < 0 OR p > 1 OR high_coverage IS NULL'
        ),
        shown=(
            'transcript',
            'n_total_transcript_regions',
            'chr',
            'region_start',
            'region_end',
            'region_oe',
            'region_oe_chisq_p',
            'is_high_coverage',
        ),
        required=(
            'an Ensembl transcript id without a version, a number of regions, a primary chromosome, a start and an '
            'end after it, a ratio of 0 or more, a p-value from 0 to 1, and true or false'
        ),
    )
    if apart := con.sql(
        'SELECT transcript_id, chrom, first, last FROM ('
        'SELECT *, lag(chrom) OVER regions AS chrom_before, lag(last) OVER regions AS last_before FROM mcr_rows '
        'WINDOW regions AS (PARTITION BY transcript_id ORDER BY chrom, first, last)'
        ') WHERE chrom_before != chrom OR last_before >= first ORDER BY ALL LIMIT 5'
    ).fetchall():
        raise ValueError(f'{path}: regions of one transcript that overlap or lie on two chromosomes, such as {apart}')
    if miscounted := con.sql(
        'SELECT transcript_id, count(*), min(regions), max(regions) FROM mcr_rows GROUP BY 1 '
        'HAVING count(*) != min(regions) OR count(*) != max(regions) ORDER BY 1 LIMIT 5'
    ).fetchall():
        raise ValueError(
            f'{path}: transcripts with another number of regions than n_total_transcript_regions says, as their id, '
            f'their rows and the lowest and the highest number that the rows give, such as {miscounted}'
        )
    con.execute(
        'CREATE TEMP TABLE mcr_values AS SELECT a.variant_id, r.transcript_id, r.oe AS mcr_oe, '  # noqa: S608 (a constant of this module)
        'r.regions AS mcr_regions, r.p AS mcr_oe_p, r.high_coverage AS mcr_high_coverage '
        f'FROM (SELECT variant_id, chrom, pos, (pos - 1) // {_STRETCH} AS stretch FROM alleles) a '
        f'JOIN (SELECT *, unnest(range((first - 1) // {_STRETCH}, (last - 1) // {_STRETCH} + 1)) AS stretch '
        'FROM mcr_rows) r ON a.chrom = r.chrom AND a.stretch = r.stretch AND a.pos BETWEEN r.first AND r.last'
    )
    return Values(
        alleles='SELECT variant_id, min(mcr_oe) AS mcr_oe_min FROM mcr_values GROUP BY 1',
        transcripts='SELECT * FROM mcr_values',
    )


# gnomAD's MPC score of missense changes (v4.1.1): one row for an allele and a transcript, with the allele as a locus
# such as `chr1:65568`, whose position is 1-based, and its reference and changed bases as a list such as `["A","C"]`.
_MPC = _Layout(
    what="gnomAD's table of MPC scores",
    whose="gnomAD's",
    header=('locus', 'alleles', 'transcript', 'mpc'),
)
_MPC_LOCUS = '^(.+):([0-9]+)$'
_MPC_ALLELES = r'^\["([ACGT]+)","([ACGT]+)"\]$'


def _mpc(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> Values:
    """Reads gnomAD's MPC scores and gives each allele of the union its score for each transcript.

    A row is kept where it is of an allele of the union, and two rows for one allele and transcript are looked for
    among the kept rows. The value of an allele is the highest score among its transcripts.

    Raises:
        ValueError: if the header is not gnomAD's, a row has another number of fields or lacks a locus on a primary
            chromosome, two alleles of bases, an Ensembl transcript id without a version or a finite number as its
            score, or two rows hold one allele of the union and one transcript.
    """
    _read_rows(
        con,
        path,
        _MPC,
        'mpc_rows',
        parsed={
            'chrom': f"regexp_extract(locus, '{_MPC_LOCUS}', 1)",
            'pos': _whole(f"regexp_extract(locus, '{_MPC_LOCUS}', 2)"),
            'ref': f"regexp_extract(alleles, '{_MPC_ALLELES}', 1)",
            'alt': f"regexp_extract(alleles, '{_MPC_ALLELES}', 2)",
            'transcript_id': 'transcript',
            'score': _finite('mpc'),
        },
        # regexp_extract gives an empty string for a text that is not of the pattern.
        malformed=(
            f"chrom NOT IN {_PRIMARY} OR pos IS NULL OR ref = '' OR alt = '' "
            f"OR NOT regexp_full_match(transcript_id, '{_ENSEMBL_TRANSCRIPT}') OR score IS NULL"
        ),
        shown=('locus', 'alleles', 'transcript', 'mpc'),
        required=(
            'a locus on a primary chromosome, a list of two alleles, an Ensembl transcript id without a version and '
            'a finite number'
        ),
        keep='(chrom, pos, ref, alt) IN (SELECT (chrom, pos, ref, alt) FROM alleles)',
    )
    _unique(con, path, 'mpc_rows', ('chrom', 'pos', 'ref', 'alt', 'transcript_id'), 'allele and transcript')
    con.execute(
        'CREATE TEMP TABLE mpc_values AS SELECT a.variant_id, m.transcript_id, m.score AS mpc '
        'FROM mpc_rows m JOIN alleles a USING (chrom, pos, ref, alt)'
    )
    return Values(
        alleles='SELECT variant_id, max(mpc) AS mpc_max FROM mpc_values GROUP BY 1',
        transcripts='SELECT * FROM mpc_values',
    )


# The prepared catalogue of branch points (branch_points): one row for a branch point and an intron of the paper's
# gene models, or one row for a branch point that the paper gives no intron. An allele is matched to the three bases at
# positions -2 to 0 of a branch point, along the branch point's strand.
_CATALOGUE = _Layout(
    what='the prepared catalogue of branch points',
    whose="the catalogue's",
    header=branch_points.COLUMNS,
    preamble=branch_points.PREAMBLE,
)
_BRANCH_POINT = pa.struct(
    [
        pa.field('name', _STRING, nullable=False),
        pa.field('pos', _INT, nullable=False),
        pa.field('strand', _STRING, nullable=False),
        pa.field('hit', pa.list_(pa.int8()), nullable=False),
        pa.field('acceptor_offset', _INT),
        pa.field('rank', _INT),
        pa.field('total', _INT),
        pa.field('sources', pa.list_(_STRING), nullable=False),
        pa.field('consensus', _STRING),
        pa.field('gene', _STRING),
        pa.field('transcript_introns', pa.list_(_STRING)),
        pa.field('intron_type', _STRING),
        pa.field('intron_length', _INT),
    ]
)


def _listed(text: str, entry: str) -> str:
    """A list from a text column that holds entries of a pattern with commas between, or null for any other text."""
    return f"CASE WHEN regexp_full_match({text}, '{entry}(,{entry})*') THEN string_split({text}, ',') END"


def _one_of(values: Sequence[str]) -> str:
    return '(' + ', '.join(f"'{value}'" for value in values) + ')'


def _branch_points(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> Values:
    """Reads the prepared catalogue and gives each allele the branch points whose motif it changes.

    An allele has an entry for each row of the catalogue whose three bases it overlaps (_overlaps), with the positions
    among -2, -1 and 0 whose base it changes. An insertion between two of the three bases changes none of them and has
    an entry without such a position. The entries of an allele are in the order of position, strand and name, and
    those of one branch point with the intron whose acceptor is nearest first.

    Raises:
        ValueError: if the header is not the catalogue's, a row has another number of fields or lacks a primary
            chromosome, a position, a strand, a name of that chromosome, position and strand, a name as the paper
            wrote it with the same kind and base, or a list of sources; fills some of the columns of an intron and
            not all, or fills them with something other than a list of a transcript's introns, a kind of intron, a
            length, a negative distance to the acceptor, a rank that is at most its total and a class of consensus;
            two rows hold one branch point in the same introns; or the file holds no row.
    """
    quoted = {column: f'"{column}"' for column in branch_points.COLUMNS}
    # A row of the paper's table by intron fills each of these columns, and a row without an intron none.
    intron_columns = branch_points.INTRON_COLUMNS
    filled = ' + '.join(f'({quoted[column]} IS NOT NULL)::INTEGER' for column in intron_columns)
    _read_rows(
        con,
        path,
        _CATALOGUE,
        'branch_point_rows',
        parsed={
            'chrom': 'chrom',
            'pos': _whole('pos'),
            'strand': 'strand',
            'name': '"name"',
            'grch37_name': 'grch37_name',
            'sources': _listed('sources', branch_points.SOURCE),
            'gene': 'gene',
            'transcript_introns': _listed('transcript_introns', branch_points.TRANSCRIPT_INTRON),
            'intron_type': 'intron_type',
            'intron_length': _whole('intron_length'),
            'acceptor_offset': (
                "CASE WHEN regexp_full_match(acceptor_offset, '-[1-9][0-9]*') "
                'THEN TRY_CAST(acceptor_offset AS INTEGER) END'
            ),
            'rank': _whole('"rank"'),
            'total': _whole('total'),
            'consensus': 'consensus',
            'filled': filled,
        },
        malformed=(
            f"chrom NOT IN {_PRIMARY} OR pos IS NULL OR strand NOT IN ('+', '-') "
            f"""OR NOT regexp_full_match("name", '{branch_points.NAME}') """
            f"""OR NOT regexp_full_match(grch37_name, '{branch_points.NAME}') """
            """OR "name" != concat("name"[1], 'BP_', chrom, '_', pos, '_', strand, '_', "name"[-1]) """
            'OR "name"[1] != grch37_name[1] OR "name"[-1] != grch37_name[-1] OR sources IS NULL '
            f'OR filled NOT IN (0, {len(intron_columns)}) '
            f'OR (filled = {len(intron_columns)} AND (transcript_introns IS NULL '
            f'OR intron_type NOT IN {_one_of(branch_points.INTRON_TYPES)} OR intron_length IS NULL '
            'OR acceptor_offset IS NULL OR "rank" IS NULL OR total IS NULL OR "rank" > total '
            f'OR consensus NOT IN {_one_of(branch_points.CONSENSUS)}))'
        ),
        shown=branch_points.COLUMNS,
        required=(
            'a primary chromosome, a position, a strand, a name of the three with a kind and a base, the name as the '
            'paper wrote it with that kind and base, and sources, and then either nothing or a gene, the introns of '
            'transcripts, a kind of intron, a length, a negative distance to the acceptor, a rank that is at most '
            'its total, and a class of consensus'
        ),
    )
    _unique(con, path, 'branch_point_rows', ('grch37_name', 'transcript_introns'), 'branch point and its introns')
    # The three bases of a row as an interval: the branch point's base and the two before it on its strand.
    con.execute(
        'CREATE TEMP TABLE branch_point_motifs AS SELECT *, row_number() OVER () AS entry, '
        "CASE strand WHEN '+' THEN 1 ELSE -1 END AS step, "
        "CASE strand WHEN '+' THEN pos - 3 ELSE pos - 1 END AS start0, "
        "CASE strand WHEN '+' THEN pos ELSE pos + 2 END AS end0 FROM branch_point_rows"
    )
    if not con.sql('SELECT 1 FROM branch_point_motifs LIMIT 1').fetchall():
        raise ValueError(f'{path}: no branch point, and the prepared catalogue holds some')
    # The 1-based position p is of a base that the allele changes where changed_start < p <= changed_end.
    hit = (
        'CAST(list_filter([-2, -1, 0], lambda k: o.changed_start < m.pos + k * m.step '
        'AND m.pos + k * m.step <= o.changed_end) AS TINYINT[])'
    )
    fields = {
        'name': 'm."name"',
        'pos': 'm.pos',
        'strand': 'm.strand',
        'hit': hit,
        **{name: f'm."{name}"' for name in _BRANCH_POINT.names[4:]},
    }
    entry = ', '.join(f"'{name}': {value}" for name, value in fields.items())
    found = _overlaps('branch_point_motifs', 'a.variant_id, a.changed_start, a.changed_end, i.entry')
    return Values(
        f'SELECT o.variant_id, list({{{entry}}} ORDER BY m.pos, m.strand, m."name", '  # noqa: S608 (names of this module's own)
        'm.acceptor_offset DESC NULLS LAST, m.transcript_introns, m.grch37_name) AS branch_points '
        f'FROM ({found}) o JOIN branch_point_motifs m USING (entry) GROUP BY o.variant_id'
    )


GNOMAD_MITO = Source(name='gnomad_mito', fields=_MITO_FIELDS, read=_gnomad_mito, normalised=True)
GENCODE_POLYA = Source(
    name='gencode_polya',
    fields=(pa.field('gencode_polya_signal', pa.bool_(), nullable=False),),
    read=_gencode_polya,
    flags=True,
)
APARENT2 = Source(
    name='aparent2',
    fields=(
        pa.field('aparent2_min', _FLOAT),
        pa.field('aparent2_max', _FLOAT),
        pa.field('aparent2_ref_usage', _FLOAT),
    ),
    read=_aparent2,
)

MITOTIP = Source(
    name='mitotip',
    fields=(pa.field('mitotip_score', _FLOAT), pa.field('mitotip_quartile', _STRING)),
    read=_mitotip,
)
GENOMIC_CONSTRAINT = Source(
    name='genomic_constraint',
    fields=(pa.field('genomic_constraint_z', _FLOAT),),
    read=_genomic_constraint,
)
# A p-value of the regions can be far below the smallest 32-bit float, so it is a 64-bit one.
MCR = Source(
    name='mcr',
    fields=(pa.field('mcr_oe_min', _FLOAT),),
    read=_mcr,
    transcript_fields=(
        pa.field('mcr_oe', _FLOAT),
        pa.field('mcr_regions', _INT),
        pa.field('mcr_oe_p', pa.float64()),
        pa.field('mcr_high_coverage', pa.bool_()),
    ),
)
MPC = Source(
    name='mpc',
    fields=(pa.field('mpc_max', _FLOAT),),
    read=_mpc,
    transcript_fields=(pa.field('mpc', _FLOAT),),
)
BRANCH_POINTS = Source(
    name='branch_points',
    fields=(pa.field('branch_points', pa.list_(_BRANCH_POINT)),),
    read=_branch_points,
)

# Every file that a run can join, in the order of their columns in the tables. A run reads any of them, or none.
SOURCES: Mapping[str, Source] = {
    source.name: source
    for source in (GNOMAD_MITO, MITOTIP, GENCODE_POLYA, APARENT2, GENOMIC_CONSTRAINT, MCR, MPC, BRANCH_POINTS)
}


def schema(sources: Sequence[Source]) -> pa.Schema:
    """The table by allele of a run that read these sources: the variant id, then each source's columns."""
    return pa.schema([*_ALLELE_KEY, *(field for source in sources for field in source.fields)])


def transcripts_schema(sources: Sequence[Source]) -> pa.Schema:
    """The table by allele and transcript of such a run: the two ids, then each source's columns of `transcripts`."""
    return pa.schema([*_TRANSCRIPT_KEY, *(field for source in sources for field in source.transcript_fields)])


def _write_alleles(sites: pathlib.Path, path: pathlib.Path) -> None:
    """Writes the union's alleles as the readers meet them (_ALLELES); a record of several alleles fails."""
    with pysam.VariantFile(str(sites)) as vcf, pq.ParquetWriter(path, _ALLELES, compression='zstd') as writer:
        if len(vcf.header.samples):
            raise ValueError(f'{sites}: the union is sites-only, and this VCF has samples')
        rows: list[dict[str, Any]] = []
        for number, record in enumerate(vcf):
            vid = annotation_tables.record_id(record, str(sites))
            chrom, pos, ref, alt = record.chrom, record.pos, str(record.ref), str(record.alts[0] if record.alts else '')
            start, end = changed_span(pos, ref, alt)
            rows.append(
                {
                    'record': number,
                    'variant_id': vid,
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'changed_start': start,
                    'changed_end': end,
                }
            )
            if len(rows) == _BATCH_ROWS:
                writer.write_table(pa.Table.from_pylist(rows, schema=_ALLELES))
                rows.clear()
        writer.write_table(pa.Table.from_pylist(rows, schema=_ALLELES))


def _write_rows(con: duckdb.DuckDBPyConnection, query: str, path: pathlib.Path, declared: pa.Schema) -> None:
    """Writes a query's rows as a Parquet file of the declared schema, whose columns the query gives in order."""
    reader = con.sql(query).to_arrow_reader(_BATCH_ROWS)
    with pq.ParquetWriter(path, declared, compression='zstd') as writer:
        for batch in reader:
            writer.write_table(pa.Table.from_batches([batch]).cast(declared))


def _pairs(con: duckdb.DuckDBPyConnection, queries: Mapping[str, str], declared: pa.Schema) -> str:
    """A query of the rows of the table by transcript, from each source's query of its values by transcript.

    An allele and a transcript have one row, with the values of every file that has some for the two.

    Raises:
        ValueError: if a source's query gives a transcript of an allele twice.
    """
    for name, query in queries.items():
        if repeated := con.sql(
            f'SELECT variant_id, transcript_id FROM ({query}) GROUP BY ALL HAVING count(*) > 1 ORDER BY ALL LIMIT 5'  # noqa: S608 (a query of this module)
        ).fetchall():
            raise ValueError(f'the values of {name} repeat a transcript of an allele, such as {repeated}')
    if not queries:
        return 'SELECT variant_id, NULL::VARCHAR AS transcript_id FROM alleles WHERE false'
    # Each source has a null in the columns of the others, and any_value takes the value that is not null.
    values = ', '.join(f'any_value(t.{name}) AS {name}' for name in declared.names[len(_TRANSCRIPT_KEY) :])
    return (
        f'SELECT t.variant_id, t.transcript_id, {values} '  # noqa: S608 (queries of this module)
        f'FROM ({" UNION ALL BY NAME ".join(f"({query})" for query in queries.values())}) t '
        'JOIN alleles a USING (variant_id) GROUP BY a.record, t.variant_id, t.transcript_id '
        'ORDER BY a.record, t.transcript_id'
    )


def _digest(path: pathlib.Path) -> str:
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def convert(
    sites: pathlib.Path,
    files: Mapping[str, pathlib.Path],
    normalised: Mapping[str, pathlib.Path],
    output: pathlib.Path,
    *,
    scratch: pathlib.Path | None = None,
) -> JoinedFiles:
    """Writes the union's values from the files, and the files' records, into the output folder.

    The inputs are the union, each file by the name of its source (SOURCES), and for each source that is a VCF of
    alleles the same file as `bcftools norm` rewrote it against the union's reference, one allele per record and
    left-aligned. The rows of both tables are in the union's order, and those of an allele in the table by transcript
    in the order of the transcripts' ids. Fails unless the alleles of the table by allele are exactly the union's,
    each once, and unless the table by transcript holds each transcript of an allele once.

    Raises:
        ValueError: if a file is of no source, a source that is a VCF of alleles comes without its normalised copy or
            another source with one, or a file is not what its source describes.
    """
    if unknown := sorted(files.keys() - SOURCES.keys()):
        raise ValueError(f'{unknown} are no files that the tables read, which are {", ".join(SOURCES)}')
    given = [source for name, source in SOURCES.items() if name in files]
    needed = {source.name for source in given if source.normalised}
    if needed != normalised.keys():
        raise ValueError(
            f'the files read as bcftools norm rewrote them are {sorted(needed)}, and the normalised copies given are '
            f'of {sorted(normalised)}'
        )
    output.mkdir(parents=True)
    table = output / TABLE_NAME
    with tempfile.TemporaryDirectory(dir=scratch) as tmp:
        work = pathlib.Path(tmp)
        _write_alleles(sites, work / 'alleles.parquet')
        con = duckdb.connect(config={'temp_directory': str(work / 'duckdb')})
        con.read_parquet(str(work / 'alleles.parquet')).create_view('alleles')
        selected = ['a.variant_id']
        joins = []
        by_transcript: dict[str, str] = {}
        for source in given:
            values = source.read(con, normalised[source.name] if source.normalised else files[source.name])
            joins.append(f'LEFT JOIN ({values.alleles}) {source.name} USING (variant_id)')
            selected += [
                f'coalesce({column}, false) AS {column}' if source.flags else column for column in source.columns
            ]
            if values.transcripts is not None:
                by_transcript[source.name] = values.transcripts
        _write_rows(
            con,
            f'SELECT {", ".join(selected)} FROM alleles a {" ".join(joins)} ORDER BY a.record',  # noqa: S608 (queries of this module)
            table,
            schema(given),
        )
        _write_rows(
            con,
            _pairs(con, by_transcript, transcripts_schema(given)),
            output / TRANSCRIPTS_NAME,
            transcripts_schema(given),
        )
        con.close()
        union = pq.read_table(work / 'alleles.parquet', columns=['variant_id'])
    annotation_tables.check_alleles(union, table, 'joined values')
    info = JoinedFiles(
        files=tuple(
            JoinedFile(
                source=source.name,
                file=files[source.name].name,
                bytes=files[source.name].stat().st_size,
                sha256=_digest(files[source.name]),
            )
            for source in given
        )
    )
    (output / INFO_NAME).write_text(info.model_dump_json(indent=2) + '\n')
    return info


def read_info(folder: pathlib.Path) -> JoinedFiles:
    """The files' records in a folder `convert` wrote."""
    return JoinedFiles.model_validate_json((folder / INFO_NAME).read_text())
