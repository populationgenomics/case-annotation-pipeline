"""Values of a sites union that are joined to its alleles from public files, as a Parquet table (docs/annotation.md).

Some values of an allele come from a small public file that is neither a source of Illumina's bundle nor an echtvar
archive: a VCF of sites, an annotation of intervals or a table of scores. This module joins each such file to the
union's alleles with DuckDB, by allele or by position, and writes one row per allele, keyed by variant id as `variants`
is, with a record of each file beside it.

A Source describes one file: its columns and the function that reads it. A run gives any of the files and none is
needed, so the table of a run holds the columns of the files that it was given. A reader checks the file as it reads
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

from case_annotation import annotation_tables, sv_sites

TABLE_NAME = 'joined_values.parquet'
INFO_NAME = 'joined_values.json'
_BATCH_ROWS = 100_000

_STRING = pa.string()
_INT = pa.int32()
_FLOAT = pa.float32()

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
class Source:
    """What the tables take from one public file.

    `name` is the file's name wherever a run meets it: the workflow's param, the `source` of its record in a family's
    manifest, and the source of its columns in the contract.

    `read` takes a connection that has the union's alleles as the view `alleles` (_ALLELES) and the file. It checks the
    file, and returns a query that gives `variant_id` and the source's columns for each allele that has a value, each
    allele at most once.
    """

    name: str
    fields: tuple[pa.Field, ...]
    read: Callable[[duckdb.DuckDBPyConnection, pathlib.Path], str]
    # Whether the file is a VCF of alleles, which is read as bcftools norm rewrote it.
    normalised: bool = False
    # Whether the columns are flags, which are false and not null for an allele without a row.
    flags: bool = False

    @property
    def columns(self) -> tuple[str, ...]:
        return tuple(field.name for field in self.fields)


def changed_span(pos: int, ref: str, alt: str) -> tuple[int, int]:
    """The reference bases that an allele changes, as a half-open interval of 0-based positions.

    The bases that REF and ALT start with alike are not changed, which leaves out the base that a VCF writes before an
    insertion or a deletion. An insertion changes no base: its interval is empty and lies between the two bases that
    the inserted sequence goes between.
    """
    shared = sum(1 for _ in itertools.takewhile(lambda bases: bases[0] == bases[1], zip(ref, alt, strict=False)))
    return pos - 1 + shared, pos - 1 + len(ref)


def _overlapping(intervals: str) -> str:
    """A query of the variant ids of the union's alleles that overlap an interval of a table.

    The table holds `chrom` and the half-open interval `start0` to `end0` of 0-based positions. An allele overlaps an
    interval that holds a base it changes. An insertion overlaps an interval that holds the base on each side of it,
    which the same comparison gives for its empty interval, and none that begins or ends where it lies.
    """
    return (
        'SELECT DISTINCT a.variant_id FROM ('  # noqa: S608 (a table name of this module's own)
        'SELECT variant_id, chrom, changed_start, changed_end, '
        f'unnest(range(changed_start // {_STRETCH}, greatest(changed_end - 1, changed_start) // {_STRETCH} + 1)) '
        'AS stretch FROM alleles) a '
        f'JOIN (SELECT chrom, start0, end0, unnest(range(start0 // {_STRETCH}, (end0 - 1) // {_STRETCH} + 1)) '
        f'AS stretch FROM {intervals}) i '
        'ON a.chrom = i.chrom AND a.stretch = i.stretch AND a.changed_start < i.end0 AND a.changed_end > i.start0'
    )


def _gzipped(path: pathlib.Path) -> bool:
    with path.open('rb') as f:
        return f.read(2) == b'\x1f\x8b'


def _text_rows(path: pathlib.Path, columns: Sequence[str], *, header: bool, comment: str | None = None) -> str:
    """A table function that reads a tab-separated file, gzipped or not, with every field as text.

    No field is quoted, and a row with another number of fields than the columns fails when the rows are read.
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


def _whole(text: str) -> str:
    """A whole number above 0 from a text column, or null for any other text."""
    number = f'TRY_CAST({text} AS INTEGER)'
    return f"CASE WHEN regexp_full_match({text}, '[0-9]+') AND {number} > 0 THEN {number} END"


def _finite(text: str) -> str:
    """A finite number from a text column, or null for any other text."""
    number = f'TRY_CAST({text} AS DOUBLE)'
    return f'CASE WHEN isfinite({number}) THEN {number} END'


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
    'mitotip_score': ('Float', 1),
    'mitotip_trna_prediction': ('String', 1),
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
    'gnomad_mito_mitotip_score': ('mitotip_score', False),
    'gnomad_mito_mitotip_prediction': ('mitotip_trna_prediction', False),
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
    pa.field('gnomad_mito_mitotip_score', _FLOAT),
    pa.field('gnomad_mito_mitotip_prediction', _STRING),
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


def _gnomad_mito(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> str:
    """Reads gnomAD's mitochondrial sites, normalised as the union is, and matches them by allele.

    Raises:
        ValueError: if the header lacks an INFO field that the columns come from or declares it otherwise, the VCF has
            samples, a line is no record, a record is not on chrM, does not hold one alternate allele, or lacks a
            field that every record of gnomAD's file has or holds no value of the field's type in it, a list by
            haplogroup holds something other than numbers, or two records hold one allele.
    """
    declared = pa.schema([pa.field('variant_id', _STRING, nullable=False), *_MITO_FIELDS])
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
    return 'SELECT * FROM gnomad_mito_rows'


# GENCODE's polyadenylation features, a GTF. A polyA_signal is the signal before a site where a transcript is cut and
# polyadenylated, nearly always six bases. The file also holds the sites themselves and signals of pseudogenes.
_GTF_COLUMNS = ('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes')
_POLYA_FEATURES = ('polyA_signal', 'polyA_site', 'pseudo_polyA')
_POLYA_SIGNAL = 'polyA_signal'


def _gencode_polya(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> str:
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
    return f'SELECT variant_id, true AS gencode_polya_signal FROM ({_overlapping("gencode_polya_signals")})'  # noqa: S608 (a query of this module)


# APARENT2's precomputed effects of the single-base changes near the polyadenylation sites of PolyA_DB: one row for a
# change and a site. A change near two sites, or near sites of genes on both strands, has several rows. The positions
# are 1-based and the bases are of the plus strand, also for a gene on the minus strand.
_APARENT2_HEADER = (
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
)


def _aparent2(con: duckdb.DuckDBPyConnection, path: pathlib.Path) -> str:
    """Reads APARENT2's effects and gives each single-base change of the union the summary of its rows.

    The file is read once. A row is kept where it is of an allele of the union, or malformed, so that the whole file
    is checked while only the union's rows are held. The lowest and the highest change of the log odds are over an
    allele's rows, and the reference usage is of the site of the lowest row: where several rows are the lowest, of the
    most used of their sites.

    Raises:
        ValueError: if the header is not APARENT2's, a row has another number of fields, or a row lacks a chromosome,
            a GRCh38 position, a base as its reference or its changed nucleotide, or a finite number as its reference
            usage or its change of the log odds.
    """
    try:
        with gzip.open(path, 'rb') if _gzipped(path) else path.open('rb') as f:
            header = f.readline().decode('utf-8', 'replace').rstrip('\r\n')
    # A gzip stream that ends within the header raises an EOFError, which is no OSError.
    except (OSError, EOFError) as unreadable:
        raise ValueError(f"{path}: not APARENT2's table of effects: {unreadable}") from unreadable
    if header != '\t'.join(_APARENT2_HEADER):
        raise ValueError(f"{path}: its header is {header.split(chr(9))}, and APARENT2's is {list(_APARENT2_HEADER)}")
    shown = ('chrom', 'var_position_hg38', 'ref_nucleotide', 'var_nucleotide', 'ref_usage', 'delta_logodds_narrow')
    # `text` holds a malformed row's fields as the file writes them, for the message.
    _create(
        con,
        'aparent2_rows',
        'WITH parsed AS ('  # noqa: S608 (a file's path, quoted, and constants of this module)
        f'SELECT chrom, {_whole("var_position_hg38")} AS pos, ref_nucleotide AS ref, var_nucleotide AS alt, '
        f'{_finite("delta_logodds_narrow")} AS delta, {_finite("ref_usage")} AS usage, '
        f'[{", ".join(shown)}] AS text FROM {_text_rows(path, _APARENT2_HEADER, header=True)}'
        '), checked AS ('
        'SELECT chrom, pos, ref, alt, delta, usage, text, coalesce('
        f'chrom IS NULL OR pos IS NULL OR ref NOT IN {_BASES} OR alt NOT IN {_BASES} '
        'OR delta IS NULL OR usage IS NULL, true) AS malformed FROM parsed'
        ') SELECT chrom, pos, ref, alt, delta, usage, malformed, CASE WHEN malformed THEN text END AS text '
        'FROM checked WHERE malformed OR (chrom, pos, ref, alt) IN ('
        'SELECT (chrom, pos, ref, alt) FROM alleles WHERE len(ref) = 1 AND len(alt) = 1)',
        path,
        "APARENT2's table of effects",
    )
    if malformed := con.sql('SELECT text FROM aparent2_rows WHERE malformed LIMIT 5').fetchall():
        raise ValueError(
            f'{path}: rows without {", ".join(shown)} as a chromosome, a position, two bases and two finite numbers, '
            f'such as {[r[0] for r in malformed]}'
        )
    return (
        'SELECT a.variant_id, s.aparent2_min, s.aparent2_max, s.aparent2_ref_usage FROM ('
        'SELECT chrom, pos, ref, alt, min(delta) AS aparent2_min, max(delta) AS aparent2_max, '
        'arg_min(usage, (delta, -usage)) AS aparent2_ref_usage FROM aparent2_rows GROUP BY ALL'
        ') s JOIN alleles a USING (chrom, pos, ref, alt)'
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

# Every file that a run can join, in the order of their columns in the tables. A run reads any of them, or none.
SOURCES: Mapping[str, Source] = {source.name: source for source in (GNOMAD_MITO, GENCODE_POLYA, APARENT2)}


def schema(sources: Sequence[Source]) -> pa.Schema:
    """The table of a run that read these sources: the variant id, then each source's columns."""
    return pa.schema(
        [pa.field('variant_id', _STRING, nullable=False), *(field for source in sources for field in source.fields)]
    )


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
    left-aligned. The table's rows are in the union's order. Fails unless its alleles are exactly the union's, each
    once.

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
    declared = schema(given)
    with tempfile.TemporaryDirectory(dir=scratch) as tmp:
        work = pathlib.Path(tmp)
        _write_alleles(sites, work / 'alleles.parquet')
        con = duckdb.connect(config={'temp_directory': str(work / 'duckdb')})
        con.read_parquet(str(work / 'alleles.parquet')).create_view('alleles')
        selected = ['a.variant_id']
        joins = []
        for source in given:
            query = source.read(con, normalised[source.name] if source.normalised else files[source.name])
            joins.append(f'LEFT JOIN ({query}) {source.name} USING (variant_id)')
            selected += [
                f'coalesce({column}, false) AS {column}' if source.flags else column for column in source.columns
            ]
        rows = f'SELECT {", ".join(selected)} FROM alleles a {" ".join(joins)} ORDER BY a.record'  # noqa: S608 (queries of this module)
        reader = con.sql(rows).to_arrow_reader(_BATCH_ROWS)
        with pq.ParquetWriter(table, declared, compression='zstd') as writer:
            for batch in reader:
                writer.write_table(pa.Table.from_batches([batch]).cast(declared))
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
