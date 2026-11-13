"""gnomAD's joint values of a sites union, from an echtvar archive, as a Parquet table (docs/annotation.md).

Illumina's bundle carries gnomAD's genomes and exomes separately and no joint values. gnomAD's own come from an echtvar
archive of its joint release: `echtvar anno` adds the archive's fields to the union's records as INFO fields, matching
each record's position, REF and ALT exactly. This module reads those records and writes one row per allele, keyed by
variant id as `variants` is, with what it read about the archive beside it.

echtvar writes a stand-in value for a field of an allele its archive lacks, and for a field the encoded VCF did not
have on an allele. The allele count's stand-in is not a count, so it tells an allele the archive lacks, whose row holds
no values, from one gnomAD counted 0 times. A frequency's stand-in is 0, so a frequency of 0 on an allele the archive
holds is gnomAD's 0 or a value gnomAD did not give.

echtvar returns a frequency as a 32-bit float, and the table's frequency columns are 32-bit floats that hold it as it
is.
"""

from __future__ import annotations

import hashlib
import json
import pathlib
import zipfile
from collections.abc import Mapping
from typing import Any, Literal

import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam

from case_annotation import annotation_tables

TABLE_NAME = 'gnomad_joint.parquet'
INFO_NAME = 'gnomad_joint.json'
_BATCH_ROWS = 100_000

# The table's columns, and for each the INFO field of gnomAD's joint sites VCF that the archive must have encoded and
# the type echtvar stores it as. The archive's own name for each, its alias, comes from the archive.
_SOURCES: Mapping[str, tuple[str, Literal['Integer', 'Float']]] = {
    'gnomad_joint_af': ('AF_joint', 'Float'),
    'gnomad_joint_ac': ('AC_joint', 'Integer'),
    'gnomad_joint_hom': ('nhomalt_joint', 'Integer'),
    'gnomad_joint_faf95_max': ('fafmax_faf95_max_joint', 'Float'),
}
# The column whose stand-in says that the archive lacks the allele.
_PRESENCE = 'gnomad_joint_ac'

TABLE = pa.schema(
    [
        pa.field('variant_id', pa.string(), nullable=False),
        pa.field('gnomad_joint_af', pa.float32()),
        pa.field('gnomad_joint_ac', pa.int32()),
        pa.field('gnomad_joint_hom', pa.int32()),
        pa.field('gnomad_joint_faf95_max', pa.float32()),
    ]
)
COLUMNS = tuple(_SOURCES)


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

        A stand-in that is htslib's missing integer reads as no value.
        """
        if value is None:
            return True
        if isinstance(value, str):
            return value == self.missing_string
        return isinstance(value, int | float) and float(value) == float(self.missing_value)


class GnomadJointInfo(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The echtvar archive the joint values come from, and the echtvar that read it."""

    archive: str  # the file's name
    bytes: int
    sha256: str
    echtvar_version: str = pydantic.Field(min_length=1)
    fields: tuple[ArchiveField, ...]


def read_archive(archive: pathlib.Path) -> tuple[ArchiveField, ...]:
    """The fields of an echtvar archive, from the configuration embedded in it.

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
        return tuple(ArchiveField(**{name: field[name] for name in ArchiveField.model_fields}) for field in config)
    # A pydantic.ValidationError and a json.JSONDecodeError are ValueErrors.
    except (KeyError, TypeError, ValueError) as malformed:
        raise ValueError(f'{archive}: its echtvar/config.json is not the fields of an echtvar archive') from malformed


def column_fields(fields: tuple[ArchiveField, ...], archive: pathlib.Path) -> dict[str, ArchiveField]:
    """The archive's field for each column of the table.

    Raises:
        ValueError: if the archive lacks one of gnomAD's joint fields, holds one twice or as another type, or writes a
            stand-in for the allele count that could be a count.
    """
    chosen = {}
    for column, (source, ftype) in _SOURCES.items():
        found = [f for f in fields if f.field == source]
        if len(found) != 1:
            raise ValueError(f"{archive}: {len(found)} fields encode gnomAD's {source}, and {column} needs one")
        if found[0].ftype != ftype:
            raise ValueError(f'{archive}: {source} is encoded as {found[0].ftype}, not {ftype}')
        chosen[column] = found[0]
    if chosen[_PRESENCE].missing_value >= 0:
        raise ValueError(
            f'{archive}: an allele it lacks gets the count {chosen[_PRESENCE].missing_value}, which a count can be'
        )
    return chosen


def _scalar(record: pysam.VariantRecord, alias: str, at: str) -> object:
    """A record's value of an INFO field with one value per record, whether or not it is declared per allele."""
    value: object = record.info.get(alias)
    if isinstance(value, tuple):
        if len(value) != 1:
            raise ValueError(f'{at}: {alias} has {len(value)} values, and the record one allele')
        return value[0]
    return value


def _stored(field: ArchiveField, value: object, at: str) -> int | float | None:
    """The count or frequency the archive stores for an allele it holds, from an annotated record's value.

    A stand-in that no count or frequency can be (no value, or a negative one) gives None; a stand-in of 0 stays 0.
    A frequency is the 32-bit float echtvar returns, unchanged: htslib reads it exactly, and the table's column holds
    it exactly.
    """
    if value is None or (field.is_stand_in(value) and field.missing_value < 0):
        return None
    if field.ftype == 'Float':
        if not isinstance(value, int | float):
            raise ValueError(f'{at}: {field.alias} is {value!r}, not a number')
        return float(value)
    if not isinstance(value, int):
        raise ValueError(f'{at}: {field.alias} is {value!r}, not an integer')
    return value


def _row(record: pysam.VariantRecord, columns: Mapping[str, ArchiveField], where: str) -> dict[str, Any]:
    at = f'{where}: {record.chrom}:{record.pos}'
    values = {column: _scalar(record, field.alias, at) for column, field in columns.items()}
    row: dict[str, Any] = {'variant_id': annotation_tables.record_id(record, where)}
    if columns[_PRESENCE].is_stand_in(values[_PRESENCE]):
        held = sorted(columns[c].alias for c in columns if not columns[c].is_stand_in(values[c]))
        if held:
            raise ValueError(f'{at}: {columns[_PRESENCE].alias} says the archive lacks the allele, yet it has {held}')
        return row
    return row | {column: _stored(field, values[column], at) for column, field in columns.items()}


def convert(
    annotated: pathlib.Path, sites: pathlib.Path, archive: pathlib.Path, echtvar_version: str, output: pathlib.Path
) -> GnomadJointInfo:
    """Writes the union's joint values and the archive's record into the output folder.

    The inputs are the union as `echtvar anno` annotated it with the archive, the union itself and the archive. Fails
    unless the annotated alleles are exactly the union's, each once.
    """
    fields = read_archive(archive)
    columns = column_fields(fields, archive)
    output.mkdir(parents=True)
    table = output / TABLE_NAME
    with pysam.VariantFile(str(annotated)) as vcf, pq.ParquetWriter(table, TABLE, compression='zstd') as writer:
        if len(vcf.header.samples):
            raise ValueError(f'{annotated}: the union is sites-only, and this VCF has samples')
        for field in columns.values():
            if field.alias not in vcf.header.info:
                raise ValueError(f'{annotated}: no INFO field {field.alias}, which {archive.name} adds')
        rows: list[dict[str, Any]] = []
        for record in vcf:
            rows.append(_row(record, columns, str(annotated)))
            if len(rows) == _BATCH_ROWS:
                writer.write_table(pa.Table.from_pylist(rows, schema=TABLE))
                rows.clear()
        writer.write_table(pa.Table.from_pylist(rows, schema=TABLE))
    annotation_tables.check_alleles(annotation_tables.read_union_sites(sites), table, 'echtvar output')
    with archive.open('rb') as f:
        sha256 = hashlib.file_digest(f, 'sha256').hexdigest()
    info = GnomadJointInfo(
        archive=archive.name,
        bytes=archive.stat().st_size,
        sha256=sha256,
        echtvar_version=echtvar_version,
        fields=fields,
    )
    (output / INFO_NAME).write_text(info.model_dump_json(indent=2) + '\n')
    return info


def read_info(folder: pathlib.Path) -> GnomadJointInfo:
    """The archive's record in a folder `convert` wrote."""
    return GnomadJointInfo.model_validate_json((folder / INFO_NAME).read_text())
