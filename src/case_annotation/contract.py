"""The pipeline's published contract (contract/): what another program can hold its own code to without this source.

`write` writes every file of contract/ but its README. Three of them come from the code that writes a family's folder:
the JSON Schema of the manifest from its model, the columns of the tables from their schemas, and the tables' version.
The families sheet and the two workflows are read by Nextflow, so what `families_sheet` and `workflows` say of them is
declared here, and tests/test_contract.py holds each statement to the .nf source. The two examples are invented.
"""

from __future__ import annotations

import json
import pathlib
from collections.abc import Mapping, Sequence
from typing import Any, Literal

import pyarrow as pa
import pydantic.json_schema

from case_annotation import (
    alphagenome_archives,
    annotation_tables,
    echtvar_values,
    family_tables,
    joined_values,
    pedigree,
)

# A family's id, a sample's, a group's and a catalog's: each names files and folders and goes into command lines.
NAME = pedigree.NAME.pattern
# A source, a category or a catalog as Illumina names it, which may hold spaces and brackets.
ILLUMINA_NAME = r'[A-Za-z0-9][A-Za-z0-9 ._()-]*'
# What separates the VCFs of a row's repeat_vcfs, and what the index beside a VCF adds to the VCF's name.
SEPARATOR = ';'
INDEX_SUFFIX = '.tbi'

# The sheet's columns that give a family each table of its larger variants: a row that names a VCF in one of them gets
# the table.
TABLE_COLUMNS: Mapping[str, tuple[str, ...]] = {
    'sv': ('sv_vcf', 'cnv_vcf'),
    'sv_transcripts': ('sv_vcf', 'cnv_vcf'),
    'repeats': ('repeat_vcfs',),
}

# The folder of the Annotator's output with one catalog. The output does not change with the tables' version, so the
# folder does not name it: a run of any version of the tables can reuse what the folder holds.
_ANNOTATION_FOLDER = 'illumina-{catalog}'
# The folder of one version of the tables made with one catalog, which a family's folder and its group's both end in.
_TABLES_FOLDER = f'tables-{{tables_version}}-{_ANNOTATION_FOLDER}'
# Where everything that holds the alleles of all of a group's families is published.
_GROUP_FOLDER = 'annotation-unions/{group}'
_DATA_FOLDER = 'references/dragen-annotation/{annotator_version}/{catalog}'

_TYPE_NAMES: Mapping[pa.DataType, str] = {
    pa.string(): 'string',
    pa.bool_(): 'bool',
    pa.int8(): 'int8',
    pa.int32(): 'int32',
    pa.float32(): 'float32',
    pa.float64(): 'float64',
}


def _without_refusal(node: Any) -> Any:  # noqa: ANN401 (a JSON value)
    """A JSON Schema without `additionalProperties: false`, at every level."""
    if isinstance(node, dict):
        return {k: _without_refusal(v) for k, v in node.items() if not (k == 'additionalProperties' and v is False)}
    if isinstance(node, list):
        return [_without_refusal(v) for v in node]
    return node


def manifest_schema() -> dict[str, Any]:
    """The JSON Schema of a family's manifest.json, from FamilyManifest, for a reader of manifests.

    The model refuses a key that it does not know, so that a misspelt key fails where the manifest is written. A key
    that is added keeps the tables' version, so a reader that holds this schema must accept a manifest with one: the
    schema allows unknown keys in every object. Every key that the schema describes is required.
    """
    schema = _without_refusal(family_tables.FamilyManifest.model_json_schema())
    return {'$schema': pydantic.json_schema.GenerateJsonSchema.schema_dialect, **schema}


def families_sheet() -> dict[str, Any]:
    """The families sheet as family_tables.nf reads it: a CSV with a header line and one row for each family.

    `column_required` says whether the header must hold the column, and `value_required` whether a row must hold a
    value in it. A `pattern` is matched against the whole value. A column with an `index_suffix` names VCFs that each
    have an index beside them, named as the VCF with the suffix, and one with a `separator` names several.
    """

    def column(
        name: str, *, column_required: bool = True, value_required: bool = True, **rules: object
    ) -> dict[str, Any]:
        return {'name': name, 'column_required': column_required, 'value_required': value_required, **rules}

    larger = {'column_required': False, 'value_required': False, 'index_suffix': INDEX_SUFFIX}
    return {
        'format': 'csv',
        'header': True,
        'columns': [
            column('family', pattern=NAME, unique=True),
            column('vcf'),
            column('ped'),
            column(
                'proband',
                value_required=False,
                pattern=NAME,
                empty_means='no proband is designated in the family, so no member has a role in the manifest',
            ),
            column('sv_vcf', **larger),
            column('cnv_vcf', **larger),
            column('repeat_vcfs', **larger, separator=SEPARATOR),
        ],
        'other_columns': 'refused',
    }


def _param(name: str, kind: str, *default: object, required: bool = False, **rules: str) -> dict[str, Any]:
    """One param as a workflow declares it: `name: Type` (required), `name: Type?` or `name: Type = default`.

    Its rules are a `pattern`, which is matched against the whole value and against each value of a list, and an
    `index_suffix`, for a file with an index beside it that is named as the file with the suffix.
    """
    declared: dict[str, Any] = {'name': name, 'type': kind, 'required': required}
    if default:
        (declared['default'],) = default
    return declared | rules


def _family_tables_workflow() -> dict[str, Any]:
    annotator = ('annotator_container', 'annotator', 'annotator_version', 'license_secret')
    return {
        'params': [
            _param('families', 'Path', required=True),
            _param('group', 'String', required=True, pattern=NAME),
            _param('reference', 'Path', required=True, index_suffix='.fai'),
            _param('catalog', 'String', required=True, pattern=NAME),
            _param('annotation_data', 'Path', required=True),
            _param('sources', 'List<String>', [], pattern=ILLUMINA_NAME),
            _param('workers', 'Integer', 8),
            *(_param(name, 'String') for name in annotator),
            _param('annotated', 'Path'),
            _param('gnomad_joint', 'Path', required=True),
            *(
                _param(source.name, 'Path')
                for source in echtvar_values.SOURCES.values()
                if source is not echtvar_values.REQUIRED
            ),
            *(_param(name, 'Path') for name in joined_values.SOURCES),
            _param('tools_container', 'String', required=True),
            _param('output_dir', 'String'),
            _param('arGuid', 'String'),
        ],
        # A run annotates its group's unions or reuses an earlier run's annotation of them, and says which by the
        # params it gives. `annotated` takes the folder that the annotating run published as `annotator`.
        'annotation': {
            'annotate': {'needs': ['sources', *annotator], 'refuses': ['annotated']},
            'reuse': {
                'needs': ['annotated'],
                'refuses': ['sources'],
                'reads': [
                    {'path': 'out/{group}.json.gz'},
                    {'path': 'run.json'},
                    {'path': 'out/{group}.sv.json.gz', 'some_row_names_any_of': list(TABLE_COLUMNS['sv'])},
                ],
            },
        },
        'publishes': {
            'family': {
                'folder': f'annotations/{{family}}/{_TABLES_FOLDER}',
                'files': [
                    *({'name': f'{table}.parquet'} for table in family_tables.TABLES if table not in TABLE_COLUMNS),
                    {'name': 'manifest.json'},
                    {'name': family_tables.NOTES},
                    *(
                        {'name': f'{table}.parquet', 'row_names_any_of': list(columns)}
                        for table, columns in TABLE_COLUMNS.items()
                    ),
                ],
            },
            'union': {'folder': f'{_GROUP_FOLDER}/{_TABLES_FOLDER}'},
            'annotator': {'folder': f'{_GROUP_FOLDER}/{_ANNOTATION_FOLDER}'},
        },
    }


def _annotator_workflow() -> dict[str, Any]:
    return {
        'params': [
            _param('mode', 'String', required=True),
            _param('annotator_container', 'String', required=True),
            _param('annotator', 'String', required=True),
            _param('annotator_version', 'String', required=True, pattern=NAME),
            _param('license_secret', 'String', required=True),
            _param('assembly', 'String', 'GRCh38', pattern=NAME),
            _param('sources', 'List<String>', [], pattern=ILLUMINA_NAME),
            _param('catalogs', 'List<String>', [], pattern=ILLUMINA_NAME),
            _param('catalog', 'String', pattern=NAME),
            _param('categories', 'List<String>', ['basic', 'premium'], pattern=ILLUMINA_NAME),
            _param('exclude_sources', 'List<String>', [], pattern=ILLUMINA_NAME),
            _param('annotation_data', 'Path'),
            _param('input', 'Path'),
            _param('name', 'String', pattern=NAME),
            _param('regions', 'Path'),
            _param('label', 'String', pattern=NAME),
            _param('workers', 'Integer', 8),
            _param('tools_container', 'String'),
            _param('output_dir', 'String'),
            _param('arGuid', 'String'),
        ],
        # The params that a run of each value of `mode` needs a value of beside the required ones, a list at least
        # one, and where the run publishes. The data directory that a download publishes is what `annotation_data`
        # takes.
        'modes': {
            'list': {'needs': [], 'publishes': {'folder': 'annotation/list/{date}'}},
            'download': {
                'needs': ['catalog', 'categories'],
                'publishes': {'folder': _DATA_FOLDER, 'data_directory': f'{_DATA_FOLDER}/{{assembly}}'},
            },
            'annotate': {
                'needs': ['catalog', 'annotation_data', 'input', 'name', 'sources'],
                'needs_with': {'regions': ['tools_container']},
                'publishes': {
                    'folder': f'annotations/{{name}}/{_ANNOTATION_FOLDER}',
                    'folder_with_label': f'annotations/{{name}}/{_ANNOTATION_FOLDER}-{{label}}',
                },
            },
        },
    }


def workflows() -> dict[str, Any]:
    """The two workflows as a launcher meets them: their params, the names they accept and the folders they publish.

    A param is `required` where the workflow declares it without a default, and has a `default` where it declares one.
    A folder is relative to `output_dir`, with the run's values in braces: its params, a family's id from the sheet,
    the tables' version, and the date on which a listing run started.
    """
    return {
        'pipeline_version': family_tables.PIPELINE_VERSION,
        'tables_version': family_tables.SCHEMA_VERSION,
        'workflows': {'family_tables.nf': _family_tables_workflow(), 'annotator.nf': _annotator_workflow()},
    }


def _type(arrow: pa.DataType, nullable: bool) -> dict[str, Any]:
    if isinstance(arrow, pa.ListType):
        element = arrow.value_field
        return {'type': 'list', 'nullable': nullable, 'element': _type(element.type, element.nullable)}
    if isinstance(arrow, pa.StructType):
        return {'type': 'struct', 'nullable': nullable, 'fields': [_column(field) for field in arrow.fields]}
    if arrow not in _TYPE_NAMES:
        raise ValueError(f'no name in the contract for the Arrow type {arrow}')
    return {'type': _TYPE_NAMES[arrow], 'nullable': nullable}


def _column(field: pa.Field) -> dict[str, Any]:
    return {'name': field.name, **_type(field.type, field.nullable)}


def tables() -> dict[str, Any]:
    """The columns of each table as its Parquet file declares them (family_tables.TABLES).

    A type is the name of an Arrow type, or `list` with its `element`, or `struct` with its `fields`. `nullable` is
    what the file declares. A column with a `source` holds the data of a source that a run can leave out, and the table
    of a run that did not name the source lacks the column: a dated source, by its name in a run config, or an echtvar
    archive or a joined file that a run need not read, by its param's name. `member` is the type of the column that each
    sequenced member has after the others, named as the member's sample.
    """
    described = []
    for name, table in family_tables.TABLES.items():
        source = {column: source for source, columns in table.sources.items() for column in columns}
        if unknown := sorted(source.keys() - set(table.columns.names)):
            raise ValueError(f'{name}: a source is given the columns {unknown}, which the table lacks')
        entry: dict[str, Any] = {
            'name': name,
            'file': f'{name}.parquet',
            'columns': [
                _column(field) | ({'source': source[field.name]} if field.name in source else {})
                for field in table.columns
            ],
        }
        if table.member is not None:
            entry['member'] = _type(table.member, nullable=False)
        described.append(entry)
    return {'tables_version': family_tables.SCHEMA_VERSION, 'tables': described}


_EXAMPLE_SAMPLES = ('S0101', 'S0102', 'S0103')
_EXAMPLE_FOLDER = 'families/family-a'
# Two invented families: a trio with a proband and the VCFs of its larger variants, and a family that names neither.
_EXAMPLE_ROWS: Sequence[Mapping[str, str]] = (
    {
        'family': 'family-a',
        'vcf': f'{_EXAMPLE_FOLDER}/family-a.vcf.gz',
        'ped': f'{_EXAMPLE_FOLDER}/family-a.ped',
        'proband': 'S0103',
        'sv_vcf': f'{_EXAMPLE_FOLDER}/family-a.sv.vcf.gz',
        'cnv_vcf': f'{_EXAMPLE_FOLDER}/family-a.cnv.vcf.gz',
        'repeat_vcfs': SEPARATOR.join(f'{_EXAMPLE_FOLDER}/{sample}.repeats.vcf.gz' for sample in _EXAMPLE_SAMPLES),
    },
    {
        'family': 'family-b',
        'vcf': 'families/family-b/family-b.vcf.gz',
        'ped': 'families/family-b/family-b.ped',
        'proband': '',
        'sv_vcf': '',
        'cnv_vcf': '',
        'repeat_vcfs': '',
    },
)


# The four fields of gnomAD's joint release that the tables read, as an archive's configuration could encode them.
_EXAMPLE_ARCHIVE: Sequence[tuple[str, str, Literal['Integer', 'Float'], int, int]] = (
    ('AC_joint', 'gnomad_AC', 'Integer', 1, -2147483648),
    ('AF_joint', 'gnomad_AF', 'Float', 2_000_000, 0),
    ('nhomalt_joint', 'gnomad_HomAlt', 'Integer', 1, -2147483648),
    ('fafmax_faf95_max_joint', 'gnomad_faf95max', 'Float', 2_000_000, 0),
)


def _example_scores(name: str) -> dict[str, Any]:
    """An invented record of an archive of one of AlphaGenome's tables, with the fields of the table's configuration."""
    table = alphagenome_archives.TABLES[name]
    return {
        'archive': f'{name}.invented.zip',
        'bytes': 1_000_000 * len(table.kept),
        'sha256': f'{len(table.kept) % 16:x}' * 64,
        'echtvar_version': '0.2.4',
        'fields': tuple(
            echtvar_values.ArchiveField(
                field=field['field'],
                alias=field['alias'],
                ftype='Float',
                multiplier=field['multiplier'],
                missing_value=field['missing_value'],
                missing_string='MISSING',
            )
            for field in alphagenome_archives.config(table).values()
        ),
    }


def example_sheet() -> str:
    """An invented families sheet, whose first row is the family of the example manifest."""
    names = [column['name'] for column in families_sheet()['columns']]
    return ''.join(','.join(row) + '\n' for row in (names, *([row[name] for name in names] for row in _EXAMPLE_ROWS)))


def _example_annotation() -> annotation_tables.AnnotationInfo:
    return annotation_tables.AnnotationInfo(
        annotator='Invented Annotator 1.0',
        json_schema_version=6,
        assembly='GRCh38',
        catalog='209901',
        sources=('Ensembl', 'clinvar', 'gnomad'),
        source_versions=(
            annotation_tables.SourceVersion(source='Ensembl', annotation_type='GeneModels', version='1'),
            annotation_tables.SourceVersion(source='clinvar', annotation_type='SmallVariant', version='20990101'),
            annotation_tables.SourceVersion(source='clinvar', annotation_type='StructuralVariant', version='20990101'),
            annotation_tables.SourceVersion(source='gnomad', annotation_type='SmallVariant', version='1'),
            annotation_tables.SourceVersion(source='gnomad', annotation_type='StructuralVariant', version='1'),
        ),
        data_sources=(
            annotation_tables.DataSource(name='Ensembl', version='1', release_date=None),
            annotation_tables.DataSource(name='ClinVar', version='20990101', release_date='2099-01-01'),
            annotation_tables.DataSource(name='gnomAD', version='1', release_date='2099-01-01'),
        ),
        omitted=(
            annotation_tables.OmittedSource(source='REVEL', reason='current'),
            annotation_tables.OmittedSource(source='dbSNP', reason='dated'),
        ),
    )


def example_manifest() -> family_tables.FamilyManifest:
    """An invented manifest: a trio with the tables of its larger variants.

    Its run did not name dbSNP, read AlphaGenome's scores and splicing scores and no breakdown of the scores, and
    joined gnomAD's mitochondrial sites and GENCODE's polyadenylation features and no effects of APARENT2.
    """

    def member(
        sample: str, father: str | None, mother: str | None, sex: pedigree.Sex, role: family_tables.Role
    ) -> family_tables.FamilyMember:
        return family_tables.FamilyMember(
            sample_id=sample,
            father=father,
            mother=mother,
            sex=sex,
            affected=role == 'proband',
            sequenced=True,
            role=role,
        )

    father, mother, child = _EXAMPLE_SAMPLES
    sizes = {'variants': 5_000_000, 'transcripts': 55_000_000, 'sv': 12_000, 'sv_transcripts': 30_000, 'repeats': 60}
    annotation = _example_annotation()
    return family_tables.FamilyManifest(
        schema_version=family_tables.SCHEMA_VERSION,
        pipeline_version=family_tables.PIPELINE_VERSION,
        family=_EXAMPLE_ROWS[0]['family'],
        members=(
            member(father, None, None, 'male', 'father'),
            member(mother, None, None, 'female', 'mother'),
            member(child, father, mother, 'female', 'proband'),
        ),
        tables=tuple(
            family_tables.TableFile(
                table=table,
                path=f'{table}.parquet',
                rows=rows,
                bytes=rows * 40,
                sha256=f'{number:x}' * 64,
            )
            for number, (table, rows) in enumerate(sizes.items(), start=1)
        ),
        notes=family_tables.NotesFile(path=family_tables.NOTES, bytes=6_000, sha256='f' * 64),
        sv_callers=('sv', 'cnv'),
        alleles=sizes['variants'],
        annotation=annotation,
        gnomad_joint=echtvar_values.ArchiveInfo(
            archive='gnomad-joint.invented.zip',
            bytes=1_000_000,
            sha256='a' * 64,
            echtvar_version='0.2.4',
            fields=tuple(
                echtvar_values.ArchiveField(
                    field=field,
                    alias=alias,
                    ftype=ftype,
                    multiplier=multiplier,
                    missing_value=missing_value,
                    missing_string='MISSING',
                )
                for field, alias, ftype, multiplier, missing_value in _EXAMPLE_ARCHIVE
            ),
        ),
        alphagenome_avi=echtvar_values.ArchiveInfo(**_example_scores('alphagenome_avi')),
        alphagenome_avi_breakdown=None,
        alphagenome_splicing=echtvar_values.NoticedArchiveInfo(
            **_example_scores('alphagenome_splicing'), terms=alphagenome_archives.NOTICE
        ),
        joined_files=tuple(
            joined_values.JoinedFile(source=source.name, file=file, bytes=size, sha256=digit * 64)
            for source, file, size, digit in (
                (joined_values.GNOMAD_MITO, 'gnomad-mito.invented.vcf.bgz', 5_000_000, 'b'),
                (joined_values.GENCODE_POLYA, 'gencode-polya.invented.gtf.gz', 2_000_000, 'c'),
            )
        ),
        omitted_columns=(
            *family_tables.VARIANT_SOURCES['dbSNP'],
            *family_tables.VARIANT_SOURCES['alphagenome_avi_breakdown'],
            *family_tables.VARIANT_SOURCES['aparent2'],
        ),
        sv_annotation=annotation,
        sv_omitted_columns=(),
    )


def files() -> dict[str, str]:
    """Each file of the contract by its path in contract/, with its text."""

    def document(value: object) -> str:
        return json.dumps(value, indent=2) + '\n'

    return {
        'manifest.schema.json': document(manifest_schema()),
        'families-sheet.json': document(families_sheet()),
        'workflow.json': document(workflows()),
        'tables.json': document(tables()),
        'tables-version': family_tables.SCHEMA_VERSION + '\n',
        'examples/manifest.json': example_manifest().model_dump_json(indent=2) + '\n',
        'examples/families.csv': example_sheet(),
    }


def write(out: pathlib.Path) -> None:
    """Writes the contract's files into the folder, replacing those of the same names and leaving any other file."""
    for name, text in files().items():
        path = out / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
