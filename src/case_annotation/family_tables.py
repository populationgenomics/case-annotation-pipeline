"""A family's annotation tables: its alleles from the union's tables, with its members' calls (docs/family-tables.md).

The input is the family's normalised VCF (one allele per record, primary chromosomes only, each allele carried by some
member), its PED file, the proband that the run's sheet names for it, if any, and the tables and the values from the
echtvar archives and from the joined files of the union that its alleles were annotated in, and the runs of
homozygosity and regions of uniparental disomy that two tools found in its own VCF (roh_upd). The output folder holds
three Parquet files, a manifest and the notes on the tables:

- `variants.parquet`: one row per allele. It holds the union's `variants` row of the allele, less the columns of the
  dated sources the run did not name, its values from the archives that the run read (gnomAD's joint values, and
  AlphaGenome's scores where the run names their archives), its values from the public files that the run joined
  (joined_values), the family's own QUAL and FILTER and, as in a VCF, one column per sequenced member, named as the
  member's sample in the VCF and in the VCF's order, with that member's own call of the allele;
- `transcripts.parquet`: the union's `transcripts` rows of the family's alleles, each with the values that the joined
  files give for its allele and its transcript;
- `roh.parquet`: one row per run of homozygosity of a sequenced member, on the autosomes;
- `manifest.json`: the tables, the pedigree with the roles of the named proband and its parents, what annotated the
  alleles (the Annotator's run, each echtvar archive that the run read and each file that it joined), and what found
  the runs and the regions;
- `tables.md`: a copy of the notes on the tables in this package (NOTES), which the manifest names.

A family whose sheet row also names VCFs of its structural variants, its copy-number segments or its members' repeats
(Structural) gets up to three more tables, which the manifest lists with the others:

- `sv.parquet`: one row per record of the family's structural-variant and copy-number callers, keyed by the caller's
  record id, with the annotation of its site in the group's union of such records (SvAnnotation), less the columns of
  the dated sources the run did not name, and, as in `variants`, one column per sequenced member with that member's
  call of the record (sv_tables). A structural variant written as a sequence that is the same allele as a row of
  `variants`, after the normalisation of the family's VCF, names that row;
- `sv_transcripts.parquet`: the transcripts the Annotator lists for each allele of each record;
- `repeats.parquet`: one row per repeat of the repeat genotyper's catalogue, with STRchive's size ranges of the repeat
  where the run gives a release of STRchive (strchive), and one column per sequenced member that holds the member's
  genotype of it (repeat_tables).

A family with a child whose parents are both sequenced gets one more:

- `upd.parquet`: one row per region in which such a child shows the alleles of one parent alone, on the autosomes.

Each table is sorted by chromosome (1 to 22, X, Y, M), then position, and written with zstd in row groups of 100,000
rows, with Bloom filters on the columns that look up a variant, a gene or a transcript. The columns whose values rarely
repeat within a row group are written without a dictionary, each with the encoding that _ENCODINGS gives it.
"""

from __future__ import annotations

import dataclasses
import hashlib
import importlib.metadata
import importlib.resources
import pathlib
import tempfile
from collections.abc import Collection, Mapping, Sequence
from typing import Any, Literal

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam

from case_annotation import (
    annotation_tables,
    echtvar_values,
    joined_values,
    pedigree,
    repeat_tables,
    roh_upd,
    strchive,
    sv_annotation,
    sv_tables,
    vcf_values,
)

# The version of these tables' layout: it names the output folder, and changes whenever a column is renamed, retyped
# or removed, or its meaning changes. A table, a column or a manifest key that is added leaves it, because a reader of
# the Parquet tables by column name, or of the manifest as plain JSON, reads what it knew as before, and the manifest
# lists the tables that a folder holds. A reader that validates the manifest with an earlier build's FamilyManifest,
# which forbids unknown keys, refuses the new key; nothing in this repository reads a manifest back.
# It is the one line of tables-version.txt in this package, so that a caller outside Python can read the same value.
SCHEMA_VERSION = importlib.resources.files('case_annotation').joinpath('tables-version.txt').read_text().strip()
# The notes on the tables (tables.md in this package): what a reader of the tables cannot see in them and gets wrong
# without an error. They list no columns: a table's schema and the contract's tables.json do. A family's folder gets a
# copy under the same name, so that a reader of the folder has the notes on the tables it reads.
NOTES = 'tables.md'
# The release of the pipeline that writes a folder: the version of this package as installed, which the tests hold to
# pyproject.toml, to pipeline-version.txt, to nextflow.config and to the contract's workflow.json. A run takes its tools
# from an image, and the family-tables command refuses to write for a workflow of another version, which the workflow
# reads from pipeline-version.txt, so the value is the workflow's as well.
PIPELINE_VERSION = importlib.metadata.version('case-annotation-tools')
ROW_GROUP_ROWS = 100_000
PRIMARY = tuple(f'chr{c}' for c in (*range(1, 23), 'X', 'Y', 'M'))
_CHROM_RANK = {c: i for i, c in enumerate(PRIMARY)}

_STRING = pa.string()
_INT = pa.int32()

# One member's call of one allele: a member's column of the variants table holds one per row. Whether the member
# carries the allele is read from its genotype alone, and is null where the genotype does not say (vcf_values.carries).
CALL = pa.struct(
    [
        pa.field('gt', _STRING, nullable=False),
        pa.field('carries', pa.bool_()),
        pa.field('zygosity', _STRING),
        pa.field('gq', _INT),
        pa.field('dp', _INT),
        pa.field('pl', pa.list_(_INT)),
        pa.field('ad_ref', _INT),
        pa.field('ad_alt', _INT),
        pa.field('allele_balance', pa.float64()),
        pa.field('ps', _INT),
        pa.field('ft', pa.list_(_STRING)),
        pa.field('dn', _STRING),
        pa.field('dq', pa.float64()),
        pa.field('msc', _INT),
    ]
)
# The family's own record of an allele: the QUAL and FILTER of its record of the VCF.
_RECORD = (pa.field('qual', pa.float64()), pa.field('filters', pa.list_(_STRING)))
# The columns of the variants table that describe the allele, in the table's order: the union's annotation of it,
# its values from each echtvar archive and from each joined file, and the family's record. The members' columns follow
# them.
ALLELE = pa.schema(
    [
        *annotation_tables.VARIANTS,
        *(field for source in echtvar_values.SOURCES.values() for field in source.fields),
        *(field for source in joined_values.SOURCES.values() for field in source.fields),
        *_RECORD,
    ]
)
# The columns of the variants table that a run can leave out, by the source that it then did not name: a dated source
# of the annotation, by its name in a run config, or an echtvar archive or a joined file that a run need not read, by
# its param's name.
VARIANT_SOURCES: Mapping[str, tuple[str, ...]] = {
    **annotation_tables.SOURCE_COLUMNS,
    **{
        source.name: tuple(source.columns)
        for source in echtvar_values.SOURCES.values()
        if source is not echtvar_values.REQUIRED
    },
    **{source.name: source.columns for source in joined_values.SOURCES.values()},
}
# The columns of the transcripts table, in the table's order: the union's annotation of an allele's transcript, then
# the values that each joined file gives for the two. The file declares every column nullable (TABLES).
TRANSCRIPT = pa.schema(
    [
        *(field.with_nullable(True) for field in annotation_tables.TRANSCRIPTS),
        *(field for source in joined_values.SOURCES.values() for field in source.transcript_fields),
    ]
)
# The columns of the transcripts table that a run can leave out, by the joined file that it then was not given.
TRANSCRIPT_SOURCES: Mapping[str, tuple[str, ...]] = {
    source.name: source.transcript_columns for source in joined_values.SOURCES.values() if source.transcript_fields
}


def variants_schema(samples: Sequence[str], omitted: Collection[str] = ()) -> pa.Schema:
    """The variants table of a VCF with these samples, in a run that leaves the omitted columns out.

    The allele's columns come first, then a call column per sample, in the VCF's order.
    """
    return pa.schema(
        [
            *(field for field in ALLELE if field.name not in omitted),
            *(pa.field(sample, CALL, nullable=False) for sample in samples),
        ]
    )


def _records_schema(samples: Sequence[str]) -> pa.Schema:
    """What read_calls takes from each record of the VCF: the allele's id, the record's QUAL and FILTER, the calls."""
    return pa.schema(
        [
            pa.field('variant_id', _STRING, nullable=False),
            *_RECORD,
            *(pa.field(sample, CALL, nullable=False) for sample in samples),
        ]
    )


def sv_schema(samples: Sequence[str], omitted: Collection[str] = ()) -> pa.Schema:
    """The sv table of a family with these samples, in a run that leaves the omitted columns out.

    The record's columns come first, then the annotation's of its site, then a call column per sample.
    """
    return pa.schema(
        [
            *sv_tables.SV,
            *(field for field in sv_annotation.COLUMNS if field.name not in omitted),
            *(pa.field(sample, sv_tables.CALL, nullable=False) for sample in samples),
        ]
    )


# The columns of the repeats table that describe the repeat, in the table's order: what the genotyper's catalogue
# says of it, then STRchive's values of it. The members' columns follow them.
REPEAT = pa.schema([*repeat_tables.LOCUS, *strchive.FIELDS])
# The columns of the repeats table that a run can leave out, by the file that it then was not given, under its param's
# name. STRchive is a dated source, and a run includes it only by giving a release of it.
REPEAT_SOURCES: Mapping[str, tuple[str, ...]] = {strchive.NAME: strchive.COLUMNS}


def repeats_schema(samples: Sequence[str], omitted: Collection[str] = ()) -> pa.Schema:
    """The repeats table of a family with these samples, in a run that leaves the omitted columns out.

    The repeat's columns come first, then a genotype column per sample.
    """
    return pa.schema(
        [
            *(field for field in REPEAT if field.name not in omitted),
            *(pa.field(sample, repeat_tables.CALL, nullable=False) for sample in samples),
        ]
    )


# The tables that hold one column per sequenced member, named as the member's sample, each with the columns that say
# what a row is of, which come before the members'.
ROW_COLUMNS: Mapping[str, pa.Schema] = {
    'variants': ALLELE,
    'sv': sv_schema(()),
    'repeats': REPEAT,
}


@dataclasses.dataclass(frozen=True)
class Table:
    """One of a family's tables as its file declares it, whatever the family's samples and the run's sources.

    `columns` are the columns that say what a row is of, in the file's order. In a table that holds the members' calls
    a column of each sequenced member follows them, named as the member's sample, in the order of the family's VCF and
    never null, and `member` is the type of each. `sources` gives the columns that hold the data of a source that a
    run can leave out, by the source's name: a dated source of the annotation by its name in a run config, and an
    echtvar archive, a joined file or a release of STRchive that a run need not read by its param's name. The table of
    a run that did not name the source lacks them.
    """

    columns: pa.Schema
    member: pa.StructType | None
    sources: Mapping[str, tuple[str, ...]]


# Every table that a family's folder can hold, in the manifest's order. transcripts is written as DuckDB returns the
# union's rows of the family's alleles, with no declared schema to cast them to (_write), so its file declares every
# column nullable, also the ones that hold a value in every row.
TABLES: Mapping[str, Table] = {
    'variants': Table(ALLELE, CALL, VARIANT_SOURCES),
    'transcripts': Table(TRANSCRIPT, None, TRANSCRIPT_SOURCES),
    'sv': Table(sv_schema(()), sv_tables.CALL, sv_annotation.SOURCE_COLUMNS),
    'sv_transcripts': Table(sv_tables.transcripts_schema(sv_annotation.TRANSCRIPT_COLUMNS), None, {}),
    'repeats': Table(REPEAT, repeat_tables.CALL, REPEAT_SOURCES),
    'roh': Table(roh_upd.ROH, None, {}),
    'upd': Table(roh_upd.UPD, None, {}),
}
# The table that only a family with a child whose parents are both sequenced has (roh_upd.trios).
TRIO_TABLE = 'upd'


def check_column_names(samples: Sequence[str]) -> None:
    """Fails unless each sample can name its column in every table that holds one column per member (ROW_COLUMNS).

    DuckDB matches column names without regard to case, so a sample may not share a name, in any case, with another
    sample or with a column of one of those tables that says what a row is of. Every such column counts: also one that
    a run leaves out because it did not name the column's source, and those of a table that the family does not get
    because its sheet row names no VCF for it. A sample name is then refused or accepted whatever the run names, and a
    query that joins two of a family's tables finds no member's column under the name of a column of the other.
    """
    seen: dict[str, tuple[str, str]] = {}
    for table, columns in ROW_COLUMNS.items():
        for name in columns.names:
            seen.setdefault(name.lower(), (name, table))
    for sample in samples:
        if (clash := seen.get(sample.lower())) is not None:
            name, table = clash
            raise ValueError(f'the VCF sample {sample!r} would share the column name {name!r} in the {table} table')
        seen[sample.lower()] = (sample, 'variants')


# The columns each table looks values up by, each by its path in the file, and about how many distinct values a row
# group holds in each. hgnc_ids is a list, and its filter is on the list's values.
_BLOOM: Mapping[str, Mapping[str, int]] = {
    'variants': {
        'variant_id': ROW_GROUP_ROWS,
        'hgnc_ids.list.element': 5_000,
    },
    'transcripts': {
        'variant_id': 20_000,
        'transcript_id': 20_000,
        'gene_symbol': 5_000,
        'gene_id': 5_000,
        'hgnc_id': 5_000,
    },
    'sv': {'sv_id': ROW_GROUP_ROWS},
    'sv_transcripts': {
        'sv_id': 20_000,
        'transcript_id': 20_000,
        'gene_symbol': 5_000,
        'gene_id': 5_000,
        'hgnc_id': 5_000,
    },
    'repeats': {'repeat_id': 1_000},
    'roh': {},
    'upd': {},
}
# The columns each table writes without a dictionary, each by its path in the file, and the encoding of its values. A
# row group holds few repeated values of each, so a dictionary would hold nearly every value once and add an index per
# row. pos rises within a chromosome, and DELTA_BINARY_PACKED stores the small differences between neighbours. The
# other columns hold plain values: the HGVS strings, the dbSNP ids, gnomAD's frequencies and its allele and homozygote
# counts, and the variant ids of variants. DELTA_BYTE_ARRAY would make those ids an eighth smaller again, and DuckDB
# took a quarter longer over a join of the table's rows on ids stored that way.
# Every other column keeps its dictionary, the fields of the members' calls and the HGNC ids among them. gnomAD's
# allele numbers repeat within a row group. In transcripts a variant id repeats over an allele's rows: plain ids would
# be a quarter smaller, and a join on them takes DuckDB twice as long.
# In the tables of the larger variants end_pos rises nearly as pos does and is stored the same way, which takes under
# half of what a dictionary of the ends does. A record id and an allele's sequence hardly ever repeat, and are plain.
# The annotation's columns of sv follow the same two rules: the interval the Annotator annotated and the intervals of
# the gnomAD matches rise as pos does, and the ids, frequencies, counts and overlaps of the matches hardly repeat.
# roh and upd hold a few thousand rows at most, and every column of theirs keeps its dictionary.
_ENCODINGS: Mapping[str, Mapping[str, str]] = {
    'variants': {
        'pos': 'DELTA_BINARY_PACKED',
        'variant_id': 'PLAIN',
        'hgvsg': 'PLAIN',
        'rsids.list.element': 'PLAIN',
        'genes.list.element.mane_select_hgvsc': 'PLAIN',
        **{
            f'gnomad_{source}_{value}': 'PLAIN'
            for source in ('genomes', 'exomes', 'joint')
            for value in ('af', 'ac', 'hom')
        },
        'gnomad_joint_faf95_max': 'PLAIN',
    },
    'transcripts': {'pos': 'DELTA_BINARY_PACKED', 'hgvsc': 'PLAIN'},
    'sv': {
        'pos': 'DELTA_BINARY_PACKED',
        'end_pos': 'DELTA_BINARY_PACKED',
        'sv_id': 'PLAIN',
        'ref': 'PLAIN',
        'alt': 'PLAIN',
        'small_variant_id': 'PLAIN',
        'annotated_pos': 'DELTA_BINARY_PACKED',
        'annotated_end_pos': 'DELTA_BINARY_PACKED',
        **{f'gnomad_sv_{value}': 'PLAIN' for value in ('id', 'af', 'ac', 'overlap', 'max_af')},
        **{f'gnomad_sv_matches.list.element.{value}': 'DELTA_BINARY_PACKED' for value in ('pos', 'end_pos')},
        **{
            f'gnomad_sv_matches.list.element.{value}': 'PLAIN'
            for value in ('id', 'af', 'ac', 'reciprocal_overlap', 'annotation_overlap')
        },
        **{f'gnomad_genomes_{value}': 'PLAIN' for value in ('af', 'ac', 'hom')},
    },
    'sv_transcripts': {'pos': 'DELTA_BINARY_PACKED', 'end_pos': 'DELTA_BINARY_PACKED'},
    'repeats': {'pos': 'DELTA_BINARY_PACKED', 'end_pos': 'DELTA_BINARY_PACKED'},
    'roh': {},
    'upd': {},
}


class TableFile(pydantic.BaseModel, frozen=True, extra='forbid'):
    table: str
    path: str  # relative to the manifest
    rows: int
    bytes: int
    sha256: str


class NotesFile(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The notes on the tables in a family's folder, named as a table's file is (TableFile)."""

    path: str  # relative to the manifest
    bytes: int
    sha256: str


# The proband that the run's sheet names (member_roles) and that member's parents by the PED file.
type Role = Literal['proband', 'father', 'mother']


class FamilyMember(pydantic.BaseModel, frozen=True, extra='forbid'):
    sample_id: str
    father: str | None
    mother: str | None
    sex: pedigree.Sex
    affected: bool | None  # None where the PED file gives the phenotype as unknown
    sequenced: bool  # a sample of the family's VCF, with a column in each table that holds the members' calls
    role: Role | None


class FamilyManifest(pydantic.BaseModel, frozen=True, extra='forbid'):
    schema_version: str
    # The release of the pipeline that wrote the folder (PIPELINE_VERSION), the same in every run of one release.
    pipeline_version: str
    family: str
    members: tuple[FamilyMember, ...]
    # variants and transcripts, then sv, sv_transcripts and repeats where the family has them, then roh, and upd
    # where the family has it
    tables: tuple[TableFile, ...]
    # The notes on the tables (NOTES), which every folder holds beside its manifest.
    notes: NotesFile
    # The callers whose records sv holds; none in a family without the table. The tables and their rows do not say
    # it: a caller that reports nothing leaves no row.
    sv_callers: tuple[sv_tables.Caller, ...]
    alleles: int
    annotation: annotation_tables.AnnotationInfo
    gnomad_joint: echtvar_values.ArchiveInfo  # the archive the gnomad_joint_ columns come from
    # The archives of AlphaGenome's scores that the avi_ columns and alphagenome_splicing come from, each None in a
    # run that was not given it. The two that are not for permissive use carry the notice their terms ask for.
    alphagenome_avi: echtvar_values.ArchiveInfo | None
    alphagenome_avi_breakdown: echtvar_values.NoticedArchiveInfo | None
    alphagenome_splicing: echtvar_values.NoticedArchiveInfo | None
    # The public files that the run joined to the alleles, each with the name of its source, in the order of their
    # columns in variants; none in a run that was given none.
    joined_files: tuple[joined_values.JoinedFile, ...]
    # The variants columns of the sources the run did not name: its dated sources, then the archives it need not read,
    # then the files it need not join.
    omitted_columns: tuple[str, ...]
    # The transcripts columns of the joined files that the run was not given.
    transcripts_omitted_columns: tuple[str, ...]
    # What annotated the family's structural variants and copy-number segments, and the sv columns of the dated
    # sources that run did not name; None and no columns in a family without the sv table.
    sv_annotation: annotation_tables.AnnotationInfo | None
    sv_omitted_columns: tuple[str, ...]
    # The release of STRchive that the strchive_ columns of repeats come from, and the repeats columns that a run
    # without one left out; None and no columns in a family without the repeats table.
    strchive: strchive.StrchiveFile | None
    repeats_omitted_columns: tuple[str, ...]
    # What found the runs of roh and the regions of upd, with the members and the children it read, and the children
    # that upd could not assess. A family whose record names no child in either list has no upd table. An upd table
    # has no row of an assessed child without a region, and none of a child that was not assessed.
    roh_upd: roh_upd.RohUpdInfo


@dataclasses.dataclass(frozen=True)
class _Call:
    """One member's call of one record's allele, as the member's column of the variants table holds it."""

    gt: str
    carries: bool | None  # by the genotype alone; None where it does not say (vcf_values.carries)
    zygosity: vcf_values.Zygosity | None
    dp: int | None
    gq: int | None
    pl: list[int] | None
    ad_ref: int | None
    ad_alt: int | None
    ps: int | None
    ft: list[str] | None
    dn: str | None
    dq: float | None
    msc: int | None

    @property
    def allele_balance(self) -> float | None:
        """The fraction of the reference's and the allele's reads that show the allele; None without such reads."""
        if self.ad_ref is None or self.ad_alt is None or self.ad_ref + self.ad_alt == 0:
            return None
        return self.ad_alt / (self.ad_ref + self.ad_alt)

    def column(self) -> dict[str, Any]:
        """The call as a value of its member's column (CALL)."""
        return {
            'gt': self.gt,
            'carries': self.carries,
            'zygosity': self.zygosity,
            'gq': self.gq,
            'dp': self.dp,
            'pl': self.pl,
            'ad_ref': self.ad_ref,
            'ad_alt': self.ad_alt,
            'allele_balance': self.allele_balance,
            'ps': self.ps,
            'ft': self.ft,
            'dn': self.dn,
            'dq': self.dq,
            'msc': self.msc,
        }


def _call(sample: pysam.VariantRecordSample, formats: frozenset[str], where: str) -> _Call:
    alleles, gt = vcf_values.genotype(sample, where)
    ad_ref, ad_alt = vcf_values.pair(sample, 'AD', formats, where, 'a depth of the reference and one of the allele')
    return _Call(
        gt=gt,
        carries=vcf_values.carries(alleles),
        zygosity=vcf_values.zygosity(alleles),
        dp=vcf_values.integer(sample, 'DP', formats, where),
        gq=vcf_values.integer(sample, 'GQ', formats, where),
        pl=vcf_values.likelihoods(sample, formats, len(alleles), where),
        ad_ref=ad_ref,
        ad_alt=ad_alt,
        ps=vcf_values.integer(sample, 'PS', formats, where),
        ft=vcf_values.names(sample, 'FT', formats, where),
        dn=vcf_values.text(sample, 'DN', formats, where),
        dq=vcf_values.number(sample, 'DQ', formats, where),
        msc=vcf_values.integer(sample, 'MSC', formats, where),
    )


class _Writer:
    """Rows into a Parquet file, in batches."""

    def __init__(self, path: pathlib.Path, schema: pa.Schema, batch: int = ROW_GROUP_ROWS) -> None:
        self._schema = schema
        self._columns: dict[str, list[Any]] = {f.name: [] for f in schema}
        self._batch = batch
        self._writer = pq.ParquetWriter(path, schema, compression='zstd')
        self.rows = 0

    def add(self, row: Mapping[str, Any]) -> None:
        for name, column in self._columns.items():
            column.append(row.get(name))
        self.rows += 1
        if self.rows % self._batch == 0:
            self._flush()

    def _flush(self) -> None:
        self._writer.write_table(pa.table(self._columns, schema=self._schema))
        for column in self._columns.values():
            column.clear()

    def close(self) -> None:
        self._flush()
        self._writer.close()


def read_calls(vcf_path: pathlib.Path, ped: pedigree.Pedigree, records_path: pathlib.Path) -> list[str]:
    """Writes each record's allele with its QUAL, its FILTER and the members' calls, unsorted.

    Returns the VCF's samples in order.
    """
    with pysam.VariantFile(str(vcf_path)) as vcf:
        samples = list(vcf.header.samples)
        pedigree.check_family(ped, ped.family, samples)
        check_column_names(samples)
        formats = frozenset(vcf.header.formats)
        if 'GT' not in formats:
            raise ValueError(f'{vcf_path}: no GT, so no genotypes')
        records = _Writer(records_path, _records_schema(samples))
        for record in vcf:
            where = f'{vcf_path}: {record.chrom}:{record.pos}'
            if record.chrom not in _CHROM_RANK:
                raise ValueError(f'{where}: not on a primary chromosome; the normalised VCF holds only those')
            vid = annotation_tables.record_id(record, str(vcf_path))
            calls = {s: _call(record.samples[s], formats, where) for s in samples}
            if not any(c.carries for c in calls.values()):
                raise ValueError(f'{where}: no member carries the allele; the normalised VCF holds only carried ones')
            records.add(
                {
                    'variant_id': vid,
                    'qual': record.qual,
                    'filters': list(record.filter.keys()),
                    **{sample: call.column() for sample, call in calls.items()},
                }
            )
        records.close()
    return samples


def _encodings(table: str, omitted: Collection[str]) -> dict[str, str]:
    """The table's columns without a dictionary and their encodings (_ENCODINGS), less those of the omitted columns."""
    return {path: encoding for path, encoding in _ENCODINGS[table].items() if path.split('.')[0] not in omitted}


def _write(
    con: duckdb.DuckDBPyConnection,
    query: str,
    path: pathlib.Path,
    table: str,
    metadata: Mapping[str, str],
    encodings: Mapping[str, str],
    declared: pa.Schema | None = None,
) -> int:
    """Writes a query's rows in row groups of ROW_GROUP_ROWS, with the table's Bloom filters; returns the rows.

    The columns named in encodings are written without a dictionary, each with its encoding, and every other column
    with one.

    The rows are cast to the declared schema where one is given, since DuckDB's Arrow output drops nullability and
    renames a list's values. The file declares no sorting columns: Parquet's would say that chrom is in ascending string
    order, and the rows are in chromosome order, which puts chr2 before chr10.
    """
    reader = con.sql(query).to_arrow_reader(ROW_GROUP_ROWS)
    target = declared if declared is not None else reader.schema
    schema = target.with_metadata({**metadata, 'table': table})
    # The filter of a column that the run leaves out is left out with it.
    bloom = {c: {'ndv': n, 'fpp': 0.01} for c, n in _BLOOM[table].items() if c.split('.')[0] in schema.names}
    rows = 0
    buffered: list[pa.RecordBatch] = []
    with pq.ParquetWriter(
        path,
        schema,
        compression='zstd',
        write_page_index=True,
        bloom_filter_options=bloom,
        **annotation_tables.column_options(schema, encodings),
    ) as writer:

        def flush(final: bool) -> None:
            nonlocal buffered
            whole = pa.Table.from_batches(buffered, schema=reader.schema)
            cut = whole.num_rows if final else whole.num_rows - whole.num_rows % ROW_GROUP_ROWS
            if cut:
                writer.write_table(whole.slice(0, cut).cast(schema), ROW_GROUP_ROWS)
            buffered = whole.slice(cut).to_batches()

        for batch in reader:
            buffered.append(batch)
            rows += batch.num_rows
            if sum(b.num_rows for b in buffered) >= ROW_GROUP_ROWS:
                flush(final=False)
        flush(final=True)
    return rows


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as f:
        while chunk := f.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def write_notes(output: pathlib.Path) -> NotesFile:
    """Copies the notes on the tables into a family's folder, byte for byte, and returns what names them there."""
    path = output / NOTES
    path.write_bytes(importlib.resources.files('case_annotation').joinpath(NOTES).read_bytes())
    return NotesFile(path=path.name, bytes=path.stat().st_size, sha256=_sha256(path))


@dataclasses.dataclass(frozen=True)
class Structural:
    """The VCFs of a family's larger variants that the run's sheet names beside the family's VCF, each of them optional.

    The structural-variant caller's and the copy-number caller's VCFs hold every sequenced member, and the repeat VCFs
    one member each; a family brings a repeat VCF of every sequenced member or of none.
    """

    sv_vcf: pathlib.Path | None
    cnv_vcf: pathlib.Path | None
    repeat_vcfs: tuple[pathlib.Path, ...]

    @property
    def callers(self) -> tuple[tuple[sv_tables.Caller, pathlib.Path], ...]:
        """The callers whose records the sv table holds, each with its VCF."""
        given: tuple[tuple[sv_tables.Caller, pathlib.Path | None], ...] = (('sv', self.sv_vcf), ('cnv', self.cnv_vcf))
        return tuple((caller, vcf) for caller, vcf in given if vcf is not None)

    def check(self, samples: Sequence[str]) -> None:
        """Fails unless each VCF holds the family's samples, from the VCFs' headers alone."""
        for caller, vcf in self.callers:
            sv_tables.check(vcf, caller, samples)
        if self.repeat_vcfs:
            repeat_tables.by_sample(self.repeat_vcfs, samples)


@dataclasses.dataclass(frozen=True)
class SvAnnotation:
    """The annotation of a group's structural variants and copy-number segments, as a family's sv table reads it.

    `tables` is the folder that sv_annotation.convert wrote for the union of the group's sites, and `site_map` the map
    from each family's records to the union's sites (sv_sites.union).
    """

    tables: pathlib.Path
    site_map: pathlib.Path

    def info(self) -> annotation_tables.AnnotationInfo:
        return annotation_tables.AnnotationInfo.model_validate_json((self.tables / 'annotation.json').read_text())


def _check_contigs(path: pathlib.Path) -> None:
    """Fails if a VCF has records and none on a primary chromosome, as one with contigs named 1, 2, ... has."""
    with pysam.VariantFile(str(path)) as vcf:
        contigs = {record.chrom for record in vcf}
    if contigs and not contigs & _CHROM_RANK.keys():
        raise ValueError(f'{path}: no record on the primary chromosomes; are its contigs named chr1, chr2, ...?')


def _read_structural(structural: Structural, samples: Sequence[str], records_path: pathlib.Path) -> set[str]:
    """Writes the callers' records as rows of sv, unsorted and without the link to a small variant.

    Returns the ids of the structural-variant caller's records that are written as a sequence, which are the ones a
    small variant can be the same allele as.
    """
    declared = sv_tables.schema(samples)
    records = _Writer(records_path, declared.remove(declared.get_field_index('small_variant_id')))
    sequences: set[str] = set()
    for caller, vcf in structural.callers:
        before = records.rows
        for site, calls in sv_tables.rows(vcf, caller, samples, _CHROM_RANK):
            records.add(site | calls)
            if caller == 'sv' and sv_tables.is_sequence(site['alt']):
                sequences.add(site['sv_id'])
        if records.rows == before:
            _check_contigs(vcf)
    records.close()
    return sequences


def read_sv_alleles(path: pathlib.Path, sequences: Collection[str]) -> pa.Table:
    """The variant id of each structural variant written as a sequence, from its normalised record.

    The file holds those records of the family's structural-variant VCF, on the primary chromosomes, after the
    normalisation that the family's own VCF gets, each under its ID. It fails unless its records are exactly the
    sequences given, each once, so that the alleles are this VCF's.
    """
    ids: dict[str, str] = {}
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            if record.id is None or record.id in ids:
                raise ValueError(f'{path}: {record.chrom}:{record.pos} has no ID, or the ID of another record')
            ids[record.id] = annotation_tables.record_id(record, str(path))
    if extra := sorted(ids.keys() - set(sequences))[:5]:
        raise ValueError(f"{path}: holds records that are no sequence of the family's structural variants: {extra}")
    if missing := sorted(set(sequences) - ids.keys())[:5]:
        raise ValueError(f"{path}: lacks sequences of the family's structural variants, such as {missing}")
    return pa.table({'sv_id': pa.array(list(ids.keys()), _STRING), 'variant_id': pa.array(list(ids.values()), _STRING)})


def _annotate_records(
    con: duckdb.DuckDBPyConnection, family: str, annotated: SvAnnotation, omitted: Collection[str]
) -> tuple[pa.Schema, str]:
    """Registers the annotation of the family's records; returns its columns and the query of sv_transcripts.

    The connection has the family's `sv_records`, and gets the family's rows of the map as `sv_site_map` and its
    records' annotation as `sv_annotated`. A record of one allele has the annotation of its site. A record of two
    alleles has that of both (sv_annotation.record_of_alleles). The joins with the map and the union read a record's
    own columns only: a member's column may bear a name of theirs, such as `site_id`. Fails unless the map's rows of
    the family are exactly the alleles of its records, and the union's annotation holds each of their sites.
    """
    site_map = pq.read_table(annotated.site_map, filters=[('family', '=', family)]).drop_columns(['family'])
    con.register('sv_site_map', site_map)
    con.read_parquet(str(annotated.tables / 'sv' / '*.parquet')).create_view('union_sv')
    con.read_parquet(str(annotated.tables / 'sv_transcripts' / '*.parquet')).create_view('union_sv_transcripts')
    con.sql("SELECT sv_id, unnest(string_split(alt, ',')) AS alt FROM sv_records").create_view('sv_record_alleles')
    for problem, query in {
        "records that the map of the group's sites lacks": (
            'SELECT DISTINCT sv_id FROM sv_record_alleles ANTI JOIN sv_site_map USING (sv_id, alt)'
        ),
        "records of the map of the group's sites that the family's callers lack": (
            'SELECT DISTINCT sv_id FROM sv_site_map ANTI JOIN sv_record_alleles USING (sv_id, alt)'
        ),
        "records whose site the group's annotation lacks": (
            'SELECT DISTINCT sv_id FROM sv_site_map ANTI JOIN union_sv USING (site_id)'
        ),
    }.items():
        if found := [r[0] for r in con.sql(f'{query} ORDER BY 1 LIMIT 5').fetchall()]:
            raise ValueError(f'{annotated.site_map}: {problem}, such as {found}')
    columns = pa.schema([field for field in sv_annotation.COLUMNS if field.name not in omitted])
    names = ', '.join(columns.names)
    several = 'SELECT sv_id FROM sv_site_map GROUP BY 1 HAVING count(*) > 1'
    # The alleles of a record of several, in the record's order.
    alleles = con.sql(
        f'SELECT m.sv_id, {names} FROM sv_site_map m JOIN union_sv USING (site_id) '  # noqa: S608 (column names of this package's own)
        'JOIN (SELECT sv_id, alt FROM sv_records) r USING (sv_id) '
        f"WHERE m.sv_id IN ({several}) ORDER BY m.sv_id, list_position(string_split(r.alt, ','), m.alt)"
    ).to_arrow_table()
    by_record: dict[str, list[dict[str, Any]]] = {}
    for row in alleles.to_pylist():
        by_record.setdefault(row['sv_id'], []).append(row)
    merged = [{'sv_id': sv_id} | sv_annotation.record_of_alleles(rows) for sv_id, rows in by_record.items()]
    con.register(
        'sv_annotated_alleles',
        pa.Table.from_pylist(merged, schema=pa.schema([pa.field('sv_id', _STRING, nullable=False), *columns])),
    )
    con.sql(
        f'SELECT m.sv_id, {names} FROM sv_site_map m JOIN union_sv USING (site_id) '  # noqa: S608 (column names of this package's own)
        f'WHERE m.sv_id NOT IN ({several}) UNION ALL BY NAME SELECT * FROM sv_annotated_alleles'
    ).create_view('sv_annotated')
    transcripts = (
        'SELECT r.chrom, r.pos, r.end_pos, r.sv_id, t.* EXCLUDE (site_id) '
        'FROM (SELECT chrom, pos, end_pos, sv_id FROM sv_records) r '
        'JOIN sv_site_map m USING (sv_id) JOIN union_sv_transcripts t USING (site_id) '
        'ORDER BY chrom_rank(r.chrom), r.pos, r.end_pos, r.sv_id, t.variant_type, t.transcript_id'
    )
    return columns, transcripts


def _structural_queries(
    con: duckdb.DuckDBPyConnection,
    family: str,
    structural: Structural,
    sv_alleles: pathlib.Path | None,
    annotated: SvAnnotation | None,
    sv_omitted: Collection[str],
    loci: Sequence[strchive.Locus] | None,
    samples: Sequence[str],
    work: pathlib.Path,
) -> tuple[dict[str, tuple[str, pa.Schema]], tuple[int, int] | None]:
    """The query and the schema of each table of the family's larger variants, by the table's name.

    Reads the VCFs into views of the connection, which also has the family's `records`. `annotated` is the annotation
    of the group's structural variants, which a family with a caller's VCF has, and `sv_omitted` the columns of sv
    that the family's table leaves out. `loci` are those of the release of STRchive that the run gives, or None in a
    run without one, whose repeats table lacks STRchive's columns. Also returns how many of the family's repeats have
    an entry of STRchive and how many have none, or None where nothing was matched: in a run without a release and in
    a family without repeat VCFs.
    """
    queries: dict[str, tuple[str, pa.Schema]] = {}
    if structural.callers:
        if annotated is None:
            raise ValueError("a family's records of structural variants need the annotation of the group's sites")
        sequences = _read_structural(structural, samples, work / 'sv.parquet')
        con.read_parquet(str(work / 'sv.parquet')).create_view('sv_records')
        if repeated := con.sql(
            'SELECT sv_id FROM sv_records GROUP BY 1 HAVING count(*) > 1 ORDER BY 1 LIMIT 5'
        ).fetchall():
            raise ValueError(f"the family's callers give several records one ID, such as {[r[0] for r in repeated]}")
        empty = pa.array([], _STRING)
        alleles = pa.table({'sv_id': empty, 'variant_id': empty})
        if sv_alleles is not None:
            alleles = read_sv_alleles(sv_alleles, sequences)
        con.register('sv_alleles', alleles)
        annotation, transcripts = _annotate_records(con, family, annotated, sv_omitted)
        columns = ', '.join(f'r.{column}' for column in sv_tables.RECORD_COLUMNS)
        # A record names a row of variants only where the family's variants table holds its allele. What a record
        # holds after its own columns are the members' calls, as sv ends.
        records = (
            f'SELECT {columns}, a.variant_id AS small_variant_id, {", ".join(f"n.{c}" for c in annotation.names)}, '  # noqa: S608 (column names of this package's own)
            f'r.* EXCLUDE ({", ".join(sv_tables.RECORD_COLUMNS)}) '
            'FROM sv_records r JOIN sv_annotated n USING (sv_id) '
            'LEFT JOIN (SELECT * FROM sv_alleles SEMI JOIN records USING (variant_id)) a USING (sv_id) '
            'ORDER BY chrom_rank(r.chrom), r.pos, r.end_pos, r.sv_id'
        )
        queries |= {
            'sv': (records, sv_schema(samples, sv_omitted)),
            'sv_transcripts': (transcripts, sv_tables.transcripts_schema(sv_annotation.TRANSCRIPT_COLUMNS)),
        }
    elif annotated is not None:
        raise ValueError("the annotation of the group's structural variants, for a family without a caller's VCF")
    counts = None
    if structural.repeat_vcfs:
        rows = repeat_tables.rows(structural.repeat_vcfs, samples, _CHROM_RANK)
        if not rows:
            for vcf in structural.repeat_vcfs:
                _check_contigs(vcf)
        if loci is not None:
            rows, matched = strchive.annotate(rows, loci)
            counts = (matched, len(rows) - matched)
        schema = repeats_schema(samples, () if loci is not None else strchive.COLUMNS)
        con.register('repeat_rows', pa.Table.from_pylist(rows, schema=schema))
        queries['repeats'] = ('SELECT * FROM repeat_rows ORDER BY chrom_rank(chrom), pos, end_pos, repeat_id', schema)
    return queries, counts


def assemble(
    family: str,
    vcf_path: pathlib.Path,
    ped_path: pathlib.Path,
    annotation: pathlib.Path,
    archive_values: pathlib.Path,
    joined: pathlib.Path,
    output: pathlib.Path,
    proband: str | None,
    structural: Structural,
    segments: pathlib.Path,
    *,
    sv_alleles: pathlib.Path | None = None,
    sv_annotated: SvAnnotation | None = None,
    ranges: strchive.Release | None = None,
    scratch: pathlib.Path | None = None,
) -> FamilyManifest:
    """Writes the family's tables and manifest into the output folder.

    The inputs are its normalised VCF, its PED file, the union's tables (annotation_tables.convert's output folder),
    the union's values from the echtvar archives (echtvar_values.convert's output folder), its values from the joined
    files (joined_values.convert's output folder), the sample id of the proband that the run's sheet names, or None for
    a family in which it designates none (member_roles), and the VCFs of its larger variants that the sheet names, with
    the normalised sequences of its structural variants where it names their VCF (read_sv_alleles) and the annotation
    of the group's structural variants where it names a caller's VCF (SvAnnotation). That annotation must be of the
    same catalog and sources as the union's, since one run makes both. `ranges` is the release of STRchive that the
    run gives, whose file is read and checked for every family, and matched with the repeats of a family that has
    them. `segments` is the folder of the family's runs of homozygosity and regions of uniparental disomy
    (roh_upd.convert's output folder), which must be of this family's samples and pedigree.
    """
    ped = pedigree.read(ped_path)
    if ped.family != family:
        raise ValueError(f'{ped_path}: the PED file is of family {ped.family!r}, not {family!r}')
    roles = member_roles(ped, _vcf_samples(vcf_path), proband)
    structural.check(_vcf_samples(vcf_path))
    if (structural.sv_vcf is None) != (sv_alleles is None):
        raise ValueError("a family's structural-variant VCF and its normalised sequences are given together")
    info = annotation_tables.AnnotationInfo.model_validate_json((annotation / 'annotation.json').read_text())
    sv_info = None if sv_annotated is None else sv_annotated.info()
    if sv_info is not None and (sv_info.catalog, sv_info.sources) != (info.catalog, info.sources):
        raise ValueError(
            f"the group's structural variants were annotated with catalog {sv_info.catalog} and the sources "
            f'{list(sv_info.sources)}, and its alleles with catalog {info.catalog} and {list(info.sources)}'
        )
    archives = echtvar_values.read_info(archive_values)
    joined_files = joined_values.read_info(joined)
    segment_info = roh_upd.read_info(segments)
    loci = None if ranges is None else strchive.read(ranges.file)
    output.mkdir(parents=True)
    with tempfile.TemporaryDirectory(dir=scratch) as tmp:
        work = pathlib.Path(tmp)
        samples = read_calls(vcf_path, ped, work / 'records.parquet')
        roh_upd.check(segment_info, ped, samples)
        con = duckdb.connect(config={'temp_directory': str(work / 'duckdb')})
        rank = ' '.join(f"WHEN '{c}' THEN {i}" for c, i in _CHROM_RANK.items())
        con.execute(f'CREATE MACRO chrom_rank(c) AS CASE c {rank} END')
        con.read_parquet(str(work / 'records.parquet')).create_view('records')
        con.read_parquet(str(annotation / 'variants' / '*.parquet')).create_view('union_variants')
        con.read_parquet(str(annotation / 'transcripts' / '*.parquet')).create_view('union_transcripts')
        con.read_parquet(str(archive_values / echtvar_values.TABLE_NAME)).create_view('union_archives')
        con.read_parquet(str(joined / joined_values.TABLE_NAME)).create_view('union_joined')
        con.read_parquet(str(joined / joined_values.TRANSCRIPTS_NAME)).create_view('union_joined_transcripts')
        for problem, query in {
            'the VCF repeats alleles': 'SELECT variant_id FROM records GROUP BY 1 HAVING count(*) > 1',
            "alleles the union's annotation lacks": (
                'SELECT variant_id FROM records ANTI JOIN union_variants USING (variant_id)'
            ),
            "alleles the union's values from the echtvar archives lack": (
                'SELECT variant_id FROM records ANTI JOIN union_archives USING (variant_id)'
            ),
            "alleles the union's values from the joined files lack": (
                'SELECT variant_id FROM records ANTI JOIN union_joined USING (variant_id)'
            ),
        }.items():
            if found := [r[0] for r in con.sql(f'{query} ORDER BY 1 LIMIT 5').fetchall()]:
                raise ValueError(f'{vcf_path}: {problem}, such as {found}')
        named = {
            *info.sources,
            *(source.name for source in archives.given()),
            *(source.name for source in joined_files.given()),
        }
        omitted = tuple(
            column for source, columns in VARIANT_SOURCES.items() if source not in named for column in columns
        )
        from_archives = [f'a.{column}' for source in archives.given() for column in source.columns]
        from_files = [f'j.{column}' for source in joined_files.given() for column in source.columns]
        transcripts_omitted = tuple(
            column for source, columns in TRANSCRIPT_SOURCES.items() if source not in named for column in columns
        )
        by_transcript = [f'x.{column}' for source in joined_files.given() for column in source.transcript_columns]
        # A joined file has a transcript under Ensembl's id without its version, and a row of the union under the id
        # with its version. A row of RefSeq's has another id and so no value.
        joined_transcripts = (
            'LEFT JOIN union_joined_transcripts x '
            "ON x.variant_id = t.variant_id AND x.transcript_id = split_part(t.transcript_id, '.', 1) "
            if by_transcript
            else ''
        )
        sv_omitted = tuple(
            column
            for source, columns in sv_annotation.SOURCE_COLUMNS.items()
            if sv_info is not None and source not in sv_info.sources
            for column in columns
        )
        # The union's variants hold the columns of every dated source, and its values from the archives and from the
        # joined files only the columns of the archives and files that the run read.
        dated = [column for column in omitted if column in annotation_tables.VARIANTS.names]
        exclude = f' EXCLUDE ({", ".join(dated)})' if dated else ''
        own = {'schema_version': SCHEMA_VERSION, 'family': family}
        # Each archive's record under its source's name, as in the manifest, and a JSON null for one the run lacks.
        records = {
            name: 'null' if record is None else record.model_dump_json() for name, record in archives.records().items()
        }
        # The joined files' records as in the manifest: a JSON list, empty for a run that joined none.
        joined_record = '[' + ','.join(file.model_dump_json() for file in joined_files.files) + ']'
        annotated = {**own, 'annotation': info.model_dump_json(), **records, 'joined_files': joined_record}
        # A record's columns after its id are the family's QUAL and FILTER and the members' calls, as variants ends.
        queries = {
            'variants': (
                f'SELECT v.*{exclude}, {", ".join([*from_archives, *from_files])}, r.* EXCLUDE (variant_id) '  # noqa: S608 (column names of this package's own)
                'FROM union_variants v JOIN records r USING (variant_id) JOIN union_archives a USING (variant_id) '
                'JOIN union_joined j USING (variant_id) ORDER BY chrom_rank(v.chrom), v.pos, v.ref, v.alt'
            ),
            'transcripts': (
                f'SELECT {", ".join(["t.*", *by_transcript])} '  # noqa: S608 (column names of this package's own)
                f'FROM union_transcripts t SEMI JOIN records USING (variant_id) {joined_transcripts}'
                'ORDER BY chrom_rank(t.chrom), t.pos, t.variant_id, t.transcript_id'
            ),
        }
        declared = {'variants': variants_schema(samples, omitted)}
        larger, counts = _structural_queries(
            con, family, structural, sv_alleles, sv_annotated, sv_omitted, loci, samples, work
        )
        ranges_record = None
        if ranges is not None and counts is not None:
            ranges_record = strchive.StrchiveFile(
                file=ranges.file.name,
                bytes=ranges.file.stat().st_size,
                sha256=_sha256(ranges.file),
                release=ranges.name,
                release_date=ranges.date,
                matched_repeats=counts[0],
                unmatched_repeats=counts[1],
            )
        repeats_omitted = REPEAT_SOURCES[strchive.NAME] if 'repeats' in larger and ranges is None else ()
        # A table of the larger variants carries the provenance of what annotated it: the annotation of the group's
        # structural variants for the two tables that hold it, and for the repeats the release of STRchive, as in the
        # manifest, or a JSON null in a run without one.
        metadata = dict.fromkeys(queries, annotated)
        sv_metadata = own if sv_info is None else {**own, 'annotation': sv_info.model_dump_json()}
        repeats_metadata = {
            **own,
            strchive.NAME: 'null' if ranges_record is None else ranges_record.model_dump_json(),
        }
        for table, (query, schema) in larger.items():
            queries[table], declared[table] = query, schema
            metadata[table] = sv_metadata if table in ('sv', 'sv_transcripts') else repeats_metadata
        # The runs and the regions as the step wrote them, sorted as the other tables are. The file of the regions is
        # there where the family has a child with both parents sequenced, whether or not upd could assess it.
        segment_tables = {'roh': (roh_upd.ROH_NAME, roh_upd.ROH)}
        if segment_info.upd_children or segment_info.upd_not_assessed:
            segment_tables[TRIO_TABLE] = (roh_upd.UPD_NAME, roh_upd.UPD)
        for table, (name, schema) in segment_tables.items():
            con.read_parquet(str(segments / name)).create_view(f'{table}_rows')
            queries[table] = f'SELECT * FROM {table}_rows ORDER BY chrom_rank(chrom), pos, end_pos, member'  # noqa: S608 (a table name of this package's own)
            declared[table] = schema
            metadata[table] = {**own, 'roh_upd': segment_info.model_dump_json()}
        files = []
        for table, query in queries.items():
            path = output / f'{table}.parquet'
            left_out = {
                'sv': sv_omitted,
                'transcripts': transcripts_omitted,
                'repeats': repeats_omitted,
                'roh': (),
                'upd': (),
            }.get(table, omitted)
            rows = _write(con, query, path, table, metadata[table], _encodings(table, left_out), declared.get(table))
            files.append(
                TableFile(table=table, path=path.name, rows=rows, bytes=path.stat().st_size, sha256=_sha256(path))
            )
        con.close()
    sequenced = set(samples)
    manifest = FamilyManifest(
        schema_version=SCHEMA_VERSION,
        pipeline_version=PIPELINE_VERSION,
        family=family,
        members=tuple(
            FamilyMember(
                sample_id=m.sample,
                father=m.father,
                mother=m.mother,
                sex=m.sex,
                affected=m.affected,
                sequenced=m.sample in sequenced,
                role=roles.get(m.sample),
            )
            for m in ped.members
        ),
        tables=tuple(files),
        notes=write_notes(output),
        sv_callers=tuple(caller for caller, _vcf in structural.callers),
        alleles=next(f.rows for f in files if f.table == 'variants'),
        annotation=info,
        gnomad_joint=archives.gnomad_joint,
        alphagenome_avi=archives.alphagenome_avi,
        alphagenome_avi_breakdown=archives.alphagenome_avi_breakdown,
        alphagenome_splicing=archives.alphagenome_splicing,
        joined_files=joined_files.files,
        omitted_columns=omitted,
        transcripts_omitted_columns=transcripts_omitted,
        sv_annotation=sv_info,
        sv_omitted_columns=sv_omitted,
        strchive=ranges_record,
        repeats_omitted_columns=repeats_omitted,
        roh_upd=segment_info,
    )
    (output / 'manifest.json').write_text(manifest.model_dump_json(indent=2) + '\n')
    return manifest


def member_roles(ped: pedigree.Pedigree, sequenced: Collection[str], proband: str | None) -> dict[str, Role]:
    """The roles in a family whose proband the run's sheet names, by sample: the proband and its parents.

    Nothing is derived. With no proband (None: the sheet designates none in this family) nobody has a role. A named
    proband must be an individual of the PED file and a sample of the VCF, and a parent of the proband by the PED file
    has a role only where sequenced, so that each role names a column of the variants table. The proband's phenotype is
    not checked: the sheet says who the proband is, and whether that member is affected is the PED file's to say.
    """
    if proband is None:
        return {}
    if proband not in {m.sample for m in ped.members}:
        raise ValueError(f'the proband {proband!r} of {ped.family!r} is not an individual of its PED file')
    if proband not in sequenced:
        raise ValueError(
            f'the proband {proband!r} of {ped.family!r} is not a sample of its VCF, so the variants table would have '
            'no column for it'
        )
    member = ped.member(proband)
    roles: dict[str, Role] = {proband: 'proband'}
    if member.father is not None and member.father in sequenced:
        roles[member.father] = 'father'
    if member.mother is not None and member.mother in sequenced:
        roles[member.mother] = 'mother'
    return roles


def _vcf_samples(vcf_path: pathlib.Path) -> list[str]:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        return list(vcf.header.samples)


def check_inputs(
    family: str, ped_path: pathlib.Path, vcf_path: pathlib.Path, proband: str | None, structural: Structural
) -> None:
    """Fails unless the family's inputs fit together, before anything is annotated.

    The PED file is the family's and names each sample of the VCF, each sample can name its column in the tables that
    hold one per member (check_column_names), and the proband that the run's sheet names, if any, is an individual of
    the PED file and a sample of the VCF. The VCFs of the family's larger variants hold exactly the VCF's samples, the
    repeat VCFs one sample each.
    """
    ped = pedigree.read(ped_path)
    samples = _vcf_samples(vcf_path)
    pedigree.check_family(ped, family, samples)
    check_column_names(samples)
    member_roles(ped, samples, proband)
    structural.check(samples)
