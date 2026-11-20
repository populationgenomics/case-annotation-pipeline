"""Values of a sites union from echtvar archives, as a Parquet table (docs/annotation.md).

Some values of an allele come from no source of Illumina's bundle: gnomAD's joint values, which the bundle lacks, and
AlphaGenome's scores. Each comes from an echtvar archive. `echtvar anno` adds the fields of every archive that it is
given to the union's records as INFO fields, matching each record's position, REF and ALT exactly. This module reads
those records and writes one row per allele, keyed by variant id as `variants` is, with what it read about each archive
beside it.

A Source describes what the tables take from one archive: its columns, the field of the encoded VCF that each column
comes from, and the columns that say whether the archive holds an allele. The archive of gnomAD's joint release is read
in every run, and an archive of AlphaGenome's scores where the run names one.

echtvar writes a stand-in value for a field of an allele its archive lacks, and for a field the encoded VCF did not
have on an allele. gnomAD's allele count has a stand-in that is not a count, so it tells an allele the archive lacks,
whose row holds no values, from one gnomAD counted 0 times. A frequency's stand-in is 0, so a frequency of 0 on an
allele the archive holds is gnomAD's 0 or a value gnomAD did not give. The VCFs that the AlphaGenome archives encode
hold every field on every record, and each field's stand-in is a missing value, so an allele has all of an archive's
values or none.

echtvar returns a frequency as a 32-bit float, and the table's frequency columns are 32-bit floats that hold it as it
is. It stores a score as a whole number of steps one decimal finer than the columns keep, cut down and not rounded, so
a score is rounded to the decimals of its source before it goes into its 32-bit column.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import pathlib
import zipfile
from collections.abc import Mapping, Sequence
from typing import Any, Literal

import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam

from case_annotation import alphagenome_archives, annotation_tables

TABLE_NAME = 'echtvar_values.parquet'
INFO_NAME = 'echtvar_values.json'
_BATCH_ROWS = 100_000


class ArchiveField(pydantic.BaseModel, frozen=True, extra='forbid'):
    """One field of an echtvar archive, as the configuration embedded in the archive records it."""

    field: str  # the INFO field of the VCF the archive encodes
    alias: str  # the INFO field echtvar adds to an annotated record
    ftype: Literal['Integer', 'Float', 'Categorical', 'Flag']
    # echtvar stores a float as a whole number of 1/multiplier. Its encoder multiplies the float by the multiplier in
    # 32-bit arithmetic and cuts the product to a whole number, which is the exact product's whole part or one below
    # it. It returns that number over the multiplier, again as a 32-bit float.
    multiplier: int
    missing_value: int
    missing_string: str

    def is_stand_in(self, value: object) -> bool:
        """Whether an annotated record's value of the field is what echtvar writes where the archive has none.

        A stand-in that is htslib's missing integer or its missing float reads as no value.
        """
        if value is None:
            return True
        if isinstance(value, str):
            return value == self.missing_string
        return isinstance(value, int | float) and float(value) == float(self.missing_value)


class ArchiveInfo(pydantic.BaseModel, frozen=True, extra='forbid'):
    """An echtvar archive that a run read values from, and the echtvar that read it."""

    archive: str  # the file's name
    bytes: int
    sha256: str
    echtvar_version: str = pydantic.Field(min_length=1)
    fields: tuple[ArchiveField, ...]


class NoticedArchiveInfo(ArchiveInfo, frozen=True, extra='forbid'):
    """An archive whose values are under terms that ask for a notice with every copy of them."""

    terms: str  # the notice, which the archive's embedded configuration holds in the description of each field read


@dataclasses.dataclass(frozen=True)
class Column:
    """A column of the table: the INFO field of the encoded VCF that it comes from, and the type echtvar stores."""

    field: str
    ftype: Literal['Integer', 'Float']

    @property
    def arrow(self) -> pa.DataType:
        return pa.float32() if self.ftype == 'Float' else pa.int32()


@dataclasses.dataclass(frozen=True)
class Source:
    """What the tables take from one echtvar archive.

    `name` is the archive's name wherever a run meets it: the workflow's param, the key of its record in a family's
    manifest, and the source of its columns in the contract. The archive's own name for each field, its alias, comes
    from the archive.
    """

    name: str
    origin: str  # whose values the archive holds, for a message
    columns: Mapping[str, Column]
    # The columns whose stand-in says that the archive lacks the allele. Such an allele has no value in any column,
    # and an allele with a value in one of these has all of them.
    presence: tuple[str, ...]
    # Whether a value can be negative. A stand-in below 0 then cannot tell an allele that the archive lacks, and only
    # a missing value can.
    negative: bool = False
    # The decimals a float is rounded to, where the archive stores it in finer steps than its source writes.
    decimals: int | None = None
    # The notice that the terms of the archive's values ask for, where they ask for one.
    notice: str | None = None

    @property
    def fields(self) -> tuple[pa.Field, ...]:
        """The source's columns as the tables declare them."""
        return tuple(pa.field(name, column.arrow) for name, column in self.columns.items())


GNOMAD_JOINT = Source(
    name='gnomad_joint',
    origin='gnomAD',
    columns={
        'gnomad_joint_af': Column('AF_joint', 'Float'),
        'gnomad_joint_ac': Column('AC_joint', 'Integer'),
        'gnomad_joint_hom': Column('nhomalt_joint', 'Integer'),
        'gnomad_joint_faf95_max': Column('fafmax_faf95_max_joint', 'Float'),
    },
    presence=('gnomad_joint_ac',),
)

# The breakdown's columns, each after the input of the score whose contribution it holds, in the order of DeepMind's
# table.
_CONTRIBUTIONS: Mapping[str, str] = {
    'avi_shap_splicing': 'MERGED_SPLICING',
    'avi_shap_atac': 'MAX_ABS_ATAC',
    'avi_shap_contact_maps': 'MAX_ABS_CONTACT_MAPS',
    'avi_shap_dnase': 'MAX_ABS_DNASE',
    'avi_shap_chip_tf': 'MAX_ABS_CHIP_TF',
    'avi_shap_chip_histone': 'MAX_ABS_CHIP_HISTONE',
    'avi_shap_cage': 'MAX_ABS_CAGE',
    'avi_shap_procap': 'MAX_ABS_PROCAP',
    'avi_shap_rna_seq': 'MAX_ABS_RNA_SEQ',
    'avi_shap_polyadenylation': 'MAX_ABS_POLYADENYLATION',
    'avi_shap_alphamissense': 'ALPHAMISSENSE',
    'avi_shap_cactus_241way': 'CACTUS_241_WAY',
    'avi_shap_protein_termination': 'PROTEIN_TERMINATION',
    'avi_shap_start_lost': 'START_LOST',
    'avi_shap_stop_lost': 'STOP_LOST',
    'avi_shap_phastcons_470way': 'PHASTCONS_470_WAY',
}


def _scores(name: str, columns: Mapping[str, str], decimals: int, notice: str | None = None) -> Source:
    """An archive of one of AlphaGenome's tables (alphagenome_archives), every field of which every row holds."""
    return Source(
        name=name,
        origin='AlphaGenome',
        columns={column: Column(field, 'Float') for column, field in columns.items()},
        presence=tuple(columns),
        negative=True,
        decimals=decimals,
        notice=notice,
    )


ALPHAGENOME_AVI = _scores('alphagenome_avi', {'avi_phred': 'PHRED', 'avi_raw': 'raw_score'}, decimals=5)
ALPHAGENOME_AVI_BREAKDOWN = _scores(
    'alphagenome_avi_breakdown', _CONTRIBUTIONS, decimals=6, notice=alphagenome_archives.NOTICE
)
ALPHAGENOME_SPLICING = _scores(
    'alphagenome_splicing',
    {'alphagenome_splicing': 'alphagenome_splicing'},
    decimals=6,
    notice=alphagenome_archives.NOTICE,
)

# Every archive that a run can read, in the order of their columns in the tables. A run reads the first and any of
# the others.
SOURCES: Mapping[str, Source] = {
    source.name: source for source in (GNOMAD_JOINT, ALPHAGENOME_AVI, ALPHAGENOME_AVI_BREAKDOWN, ALPHAGENOME_SPLICING)
}
REQUIRED = GNOMAD_JOINT


class Archives(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The archive of each source (SOURCES) that a run read, and None for one that the run was not given."""

    gnomad_joint: ArchiveInfo
    alphagenome_avi: ArchiveInfo | None
    alphagenome_avi_breakdown: NoticedArchiveInfo | None
    alphagenome_splicing: NoticedArchiveInfo | None

    def given(self) -> tuple[Source, ...]:
        """The sources that the run read an archive of, in the tables' order."""
        return tuple(source for name, source in SOURCES.items() if getattr(self, name) is not None)

    def records(self) -> dict[str, ArchiveInfo | None]:
        """Each source's record by the source's name, which is its key in a manifest."""
        return {name: getattr(self, name) for name in SOURCES}


def schema(sources: Sequence[Source]) -> pa.Schema:
    """The table of a run that read these sources: the variant id, then each source's columns."""
    return pa.schema(
        [pa.field('variant_id', pa.string(), nullable=False), *(field for source in sources for field in source.fields)]
    )


def _embedded(archive: pathlib.Path) -> list[dict[str, Any]]:
    """The configuration embedded in an echtvar archive: an object for each of its fields.

    Raises:
        ValueError: if the file is not a zip archive, or holds no echtvar configuration or a malformed one.
    """
    try:
        with zipfile.ZipFile(archive) as z:
            embedded = z.read('echtvar/config.json')
    except zipfile.BadZipFile as bad:
        raise ValueError(f'{archive}: not a zip file, so not an echtvar archive') from bad
    except KeyError as missing:
        raise ValueError(f'{archive}: no echtvar/config.json, so not an echtvar archive') from missing
    try:
        config = json.loads(embedded)
        for field in config:
            ArchiveField(**{name: field[name] for name in ArchiveField.model_fields})
    # A pydantic.ValidationError and a json.JSONDecodeError are ValueErrors.
    except (KeyError, TypeError, ValueError) as malformed:
        raise ValueError(f'{archive}: its echtvar/config.json is not the fields of an echtvar archive') from malformed
    return config


def read_archive(archive: pathlib.Path) -> tuple[ArchiveField, ...]:
    """The fields of an echtvar archive, from the configuration embedded in it.

    Raises:
        ValueError: if the file is not a zip archive, or holds no echtvar configuration or a malformed one.
    """
    return tuple(
        ArchiveField(**{name: field[name] for name in ArchiveField.model_fields}) for field in _embedded(archive)
    )


def _never_a_value(source: Source, field: ArchiveField) -> bool:
    """Whether the stand-in of a field is no value that the field can hold, so that it marks an absent allele.

    echtvar writes a float field whose stand-in has the bits of htslib's missing float as a missing value.
    """
    if field.ftype == 'Float' and field.missing_value == alphagenome_archives.FLOAT_MISSING:
        return True
    return not source.negative and field.missing_value < 0


def column_fields(source: Source, fields: tuple[ArchiveField, ...], archive: pathlib.Path) -> dict[str, ArchiveField]:
    """The archive's field for each column of a source.

    Raises:
        ValueError: if the archive lacks one of the source's fields, holds one twice or as another type, or writes
            a stand-in that a value could be for a field that says whether the archive holds an allele.
    """
    chosen = {}
    for name, column in source.columns.items():
        found = [f for f in fields if f.field == column.field]
        if len(found) != 1:
            raise ValueError(
                f"{archive}: {len(found)} fields encode {source.origin}'s {column.field}, and {name} needs one"
            )
        if found[0].ftype != column.ftype:
            raise ValueError(f'{archive}: {column.field} is encoded as {found[0].ftype}, not {column.ftype}')
        chosen[name] = found[0]
    for name in source.presence:
        if not _never_a_value(source, chosen[name]):
            raise ValueError(
                f'{archive}: an allele it lacks gets the value {chosen[name].missing_value} as its '
                f'{chosen[name].alias}, which a value of the field can be'
            )
    return chosen


def read_notice(source: Source, archive: pathlib.Path) -> str | None:
    """The notice of a source whose terms ask for one, as the archive's embedded configuration holds it.

    Raises:
        ValueError: if the description of a field that the source reads does not end with the notice: the archive
            would then annotate a VCF without it, and it was not built with the configuration of the source's table.
    """
    if source.notice is None:
        return None
    read = {column.field for column in source.columns.values()}
    for field in _embedded(archive):
        if field['field'] in read and not str(field.get('description', '')).endswith(source.notice):
            raise ValueError(
                f'{archive}: the description of {field["field"]} does not end with the notice that the terms of '
                f"{source.origin}'s values ask for: {source.notice!r}"
            )
    return source.notice


def _scalar(record: pysam.VariantRecord, alias: str, at: str) -> object:
    """A record's value of an INFO field with one value per record, whether or not it is declared per allele."""
    value: object = record.info.get(alias)
    if isinstance(value, tuple):
        if len(value) != 1:
            raise ValueError(f'{at}: {alias} has {len(value)} values, and the record one allele')
        return value[0]
    return value


def _stored(source: Source, field: ArchiveField, value: object, at: str) -> int | float | None:
    """The value the archive stores for an allele it holds, from an annotated record's value.

    In a source without negative values, a stand-in that no value can be (no value, or a negative one) gives None, and
    a stand-in of 0 stays 0. A float of a source without `decimals` is the 32-bit float echtvar returns, unchanged:
    htslib reads it exactly, and the table's column holds it exactly. A float of a source with them is rounded.
    """
    if value is None or (field.is_stand_in(value) and field.missing_value < 0 and not source.negative):
        return None
    if field.ftype == 'Float':
        if not isinstance(value, int | float):
            raise ValueError(f'{at}: {field.alias} is {value!r}, not a number')
        return float(value) if source.decimals is None else round(float(value), source.decimals)
    if not isinstance(value, int):
        raise ValueError(f'{at}: {field.alias} is {value!r}, not an integer')
    return value


def _values(
    record: pysam.VariantRecord, source: Source, columns: Mapping[str, ArchiveField], at: str
) -> dict[str, Any]:
    """A record's values of one source's columns, and none where the archive lacks the record's allele."""
    values = {column: _scalar(record, field.alias, at) for column, field in columns.items()}
    missing = [columns[c].alias for c in source.presence if columns[c].is_stand_in(values[c])]
    if missing:
        held = sorted(columns[c].alias for c in columns if not columns[c].is_stand_in(values[c]))
        if held:
            raise ValueError(f'{at}: {missing[0]} says the archive lacks the allele, yet it has {held}')
        return {}
    return {column: _stored(source, field, values[column], at) for column, field in columns.items()}


def _digest(archive: pathlib.Path) -> str:
    with archive.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def convert(
    annotated: pathlib.Path,
    sites: pathlib.Path,
    archives: Mapping[str, pathlib.Path],
    echtvar_version: str,
    output: pathlib.Path,
) -> Archives:
    """Writes the union's values from the archives, and the archives' records, into the output folder.

    The inputs are the union as one `echtvar anno` call annotated it with every archive, the union itself and each
    archive by the name of its source (SOURCES). Fails unless the annotated alleles are exactly the union's, each once.

    Raises:
        ValueError: if an archive is of no source, the archive of gnomAD's joint release is not among them, two fields
            of the archives share an alias, or an archive or the annotated VCF is not what its source describes.
    """
    if unknown := sorted(archives.keys() - SOURCES.keys()):
        raise ValueError(f'{unknown} are no archives that the tables read, which are {", ".join(SOURCES)}')
    if REQUIRED.name not in archives:
        raise ValueError(f'no archive of {REQUIRED.name}, which every run reads')
    given = [source for name, source in SOURCES.items() if name in archives]
    fields = {source.name: read_archive(archives[source.name]) for source in given}
    # echtvar adds every field of every archive, and a later archive's field replaces one of the same alias.
    seen: dict[str, str] = {}
    for source in given:
        for field in fields[source.name]:
            if (other := seen.setdefault(field.alias, source.name)) != source.name:
                raise ValueError(
                    f'{archives[source.name]} and {archives[other]} both add the INFO field {field.alias}, and echtvar '
                    'would keep one of the two'
                )
    columns = {source.name: column_fields(source, fields[source.name], archives[source.name]) for source in given}
    notices = {source.name: read_notice(source, archives[source.name]) for source in given}
    output.mkdir(parents=True)
    table = output / TABLE_NAME
    declared = schema(given)
    with pysam.VariantFile(str(annotated)) as vcf, pq.ParquetWriter(table, declared, compression='zstd') as writer:
        if len(vcf.header.samples):
            raise ValueError(f'{annotated}: the union is sites-only, and this VCF has samples')
        for source in given:
            for field in columns[source.name].values():
                if field.alias not in vcf.header.info:
                    raise ValueError(
                        f'{annotated}: no INFO field {field.alias}, which {archives[source.name].name} adds'
                    )
        rows: list[dict[str, Any]] = []
        for record in vcf:
            at = f'{annotated}: {record.chrom}:{record.pos}'
            row: dict[str, Any] = {'variant_id': annotation_tables.record_id(record, str(annotated))}
            for source in given:
                row |= _values(record, source, columns[source.name], at)
            rows.append(row)
            if len(rows) == _BATCH_ROWS:
                writer.write_table(pa.Table.from_pylist(rows, schema=declared))
                rows.clear()
        writer.write_table(pa.Table.from_pylist(rows, schema=declared))
    annotation_tables.check_alleles(annotation_tables.read_union_sites(sites), table, 'echtvar output')
    records: dict[str, ArchiveInfo | None] = dict.fromkeys(SOURCES)
    for source in given:
        archive = archives[source.name]
        entry: dict[str, Any] = {
            'archive': archive.name,
            'bytes': archive.stat().st_size,
            'sha256': _digest(archive),
            'echtvar_version': echtvar_version,
            'fields': fields[source.name],
        }
        notice = notices[source.name]
        records[source.name] = (
            ArchiveInfo.model_validate(entry)
            if notice is None
            else NoticedArchiveInfo.model_validate(entry | {'terms': notice})
        )
    info = Archives.model_validate(records)
    (output / INFO_NAME).write_text(info.model_dump_json(indent=2) + '\n')
    return info


def read_info(folder: pathlib.Path) -> Archives:
    """The archives' records in a folder `convert` wrote."""
    return Archives.model_validate_json((folder / INFO_NAME).read_text())
