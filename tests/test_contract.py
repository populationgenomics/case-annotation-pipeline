"""The published contract (contract/), held to the code that writes a family's folder and to the workflows' source.

The workflows are Nextflow's, so the tests read their source by regular expression: a statement of the contract that a
workflow no longer makes fails here.
"""

from __future__ import annotations

import csv
import importlib.resources
import io
import json
import pathlib
import re
import string
import tomllib
from collections.abc import Iterator, Mapping, Sequence
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import structural_vcfs
import test_family_tables

import case_annotation
from case_annotation import cli, contract, echtvar_values, family_tables, joined_values, strchive

REPOSITORY = pathlib.Path(__file__).parent.parent
CONTRACT = REPOSITORY / 'contract'
FAMILY_TABLES = (REPOSITORY / 'family_tables.nf').read_text('utf-8')
ANNOTATOR = (REPOSITORY / 'annotator.nf').read_text('utf-8')
FAMILY_MODULES = (REPOSITORY / 'modules' / 'family_tables.nf').read_text('utf-8')
ANNOTATOR_MODULES = (REPOSITORY / 'modules' / 'annotator.nf').read_text('utf-8')
INCLUDING = (REPOSITORY / 'tests' / 'workflows' / 'include' / 'main.nf').read_text('utf-8')
STUB_RUNS = (REPOSITORY / '.github' / 'workflows' / 'nextflow.yml').read_text('utf-8')

SHEET = contract.families_sheet()
COLUMNS: Mapping[str, Mapping[str, Any]] = {column['name']: column for column in SHEET['columns']}
WORKFLOWS: Mapping[str, Mapping[str, Any]] = contract.workflows()['workflows']
FAMILY_FILES: Sequence[Mapping[str, Any]] = WORKFLOWS['family_tables.nf']['publishes']['family']['files']
TABLES: Mapping[str, Mapping[str, Any]] = {table['name']: table for table in contract.tables()['tables']}


def _files(folder: pathlib.Path) -> dict[str, str]:
    return {str(path.relative_to(folder)): path.read_text('utf-8') for path in folder.rglob('*') if path.is_file()}


def _names(listed: str) -> list[str]:
    """The names of a Groovy list of quoted strings, from the text between its brackets."""
    return [name.strip().strip("'") for name in listed.split(',')]


def _function(text: str, name: str) -> str:
    """The body of a function at the top level of a workflow's source."""
    opening = f'\ndef {name}('
    assert text.count(opening) == 1
    start = text.index(opening)
    return text[start : text.index('\n}\n', start)]


def _block(text: str, opening: str) -> str:
    """A block at the top level of a workflow's source, from the line that opens it to the one that closes it."""
    assert text.count(f'\n{opening} {{\n') == 1
    start = text.index(f'\n{opening} {{\n')
    return text[start : text.index('\n}\n', start)]


def _code(text: str) -> str:
    """A workflow's source without its comment lines."""
    return '\n'.join(line for line in text.splitlines() if not line.lstrip().startswith('//'))


# ---- the folder ----


def test_the_committed_contract_is_the_one_the_command_writes(tmp_path: pathlib.Path) -> None:
    cli.main(['contract', '--out', str(tmp_path / 'contract')])
    committed = _files(CONTRACT)
    assert committed.pop('README.md')
    assert _files(tmp_path / 'contract') == committed, 'write it again: uv run case-annotation contract --out contract'


def test_writing_the_contract_again_replaces_its_files_and_leaves_the_others(tmp_path: pathlib.Path) -> None:
    contract.write(tmp_path)
    (tmp_path / 'README.md').write_text('kept\n')
    (tmp_path / 'tables-version').write_text('an earlier version\n')
    contract.write(tmp_path)
    assert _files(tmp_path) == {**contract.files(), 'README.md': 'kept\n'}


def test_the_contracts_versions_are_the_tables_and_the_packages() -> None:
    version = importlib.resources.files(case_annotation).joinpath('tables-version.txt').read_text()
    assert contract.files()['tables-version'] == version
    assert contract.tables()['tables_version'] == contract.workflows()['tables_version'] == version.strip()
    project = tomllib.loads((REPOSITORY / 'pyproject.toml').read_text('utf-8'))['project']
    assert contract.workflows()['pipeline_version'] == project['version']


# ---- the manifest ----


def _objects(node: object) -> Iterator[dict[str, Any]]:
    """Every JSON object of a document, the document itself among them."""
    if isinstance(node, dict):
        yield node
    if isinstance(node, dict | list):
        for value in node.values() if isinstance(node, dict) else node:
            yield from _objects(value)


def test_the_published_schema_is_the_models_without_its_refusal_of_unknown_keys() -> None:
    published = contract.manifest_schema()
    assert not any('additionalProperties' in node for node in _objects(published))
    # The model refuses an unknown key in every object that it describes, and the published schema differs in nothing
    # else.
    described = [node for node in _objects(published) if 'properties' in node]
    assert len(described) > 1
    for node in described:
        node['additionalProperties'] = False
    assert published.pop('$schema') == 'https://json-schema.org/draft/2020-12/schema'
    assert published == family_tables.FamilyManifest.model_json_schema()


def test_every_key_that_the_schema_describes_is_required() -> None:
    schema = contract.manifest_schema()
    assert set(schema['properties']) == set(family_tables.FamilyManifest.model_fields)
    # In the manifest and in every object within it: no key is left out for a value that stands for nothing.
    described = [node for node in _objects(schema) if 'properties' in node]
    assert len(described) > 1
    for node in described:
        assert set(node['required']) == set(node['properties']), node.get('title')
        assert not any('default' in key for key in node['properties'].values()), node.get('title')


def test_a_members_phenotype_is_required_and_may_be_null() -> None:
    member = contract.manifest_schema()['$defs']['FamilyMember']
    assert 'affected' in member['required']
    assert {option['type'] for option in member['properties']['affected']['anyOf']} == {'boolean', 'null'}


def test_the_example_manifest_is_one_that_the_model_reads_and_writes() -> None:
    text = (CONTRACT / 'examples' / 'manifest.json').read_text('utf-8')
    manifest = family_tables.FamilyManifest.model_validate_json(text)
    assert manifest.model_dump_json(indent=2) + '\n' == text
    assert manifest.schema_version == family_tables.SCHEMA_VERSION
    # What a run's manifest always holds: the tables in their order, the alleles as the rows of variants, each role on
    # a sequenced member, and as omitted columns those of a dated source that the run did not name.
    assert [table.table for table in manifest.tables] == list(family_tables.TABLES)
    assert manifest.alleles == manifest.tables[0].rows
    assert all(member.sequenced for member in manifest.members if member.role is not None)
    named = (
        set(manifest.annotation.sources)
        | {name for name in echtvar_values.SOURCES if getattr(manifest, name) is not None}
        | {file.source for file in manifest.joined_files}
        | ({strchive.NAME} if manifest.strchive is not None else set())
    )
    for omitted, table in (
        (manifest.omitted_columns, 'variants'),
        (manifest.transcripts_omitted_columns, 'transcripts'),
        (manifest.sv_omitted_columns, 'sv'),
        (manifest.repeats_omitted_columns, 'repeats'),
    ):
        sources = family_tables.TABLES[table].sources
        assert omitted == tuple(c for source, columns in sources.items() if source not in named for c in columns)
    # Each of its archives holds the fields that a run reads the source's values from, and it shows an archive that
    # the run read, one with a notice and one that the run was not given.
    for name, source in echtvar_values.SOURCES.items():
        record = getattr(manifest, name)
        assert record is None or echtvar_values.column_fields(source, record.fields, pathlib.Path(record.archive))
    assert manifest.alphagenome_avi is not None
    assert manifest.alphagenome_avi_breakdown is None
    assert manifest.alphagenome_splicing is not None
    assert manifest.alphagenome_splicing.terms == echtvar_values.ALPHAGENOME_SPLICING.notice
    # It shows files that the run joined and one that it was not given, in the order of their columns, and of the
    # files with columns in transcripts one of each.
    joined = [file.source for file in manifest.joined_files]
    assert joined == [name for name in joined_values.SOURCES if name in joined]
    assert 0 < len(joined) < len(joined_values.SOURCES)
    assert 0 < len(manifest.transcripts_omitted_columns) < len(test_family_tables._JOINED_TRANSCRIPT_COLUMNS)
    assert manifest.sv_annotation is not None
    assert (manifest.sv_annotation.catalog, manifest.sv_annotation.sources) == (
        manifest.annotation.catalog,
        manifest.annotation.sources,
    )
    # It shows a release of STRchive, whose counts are of the rows of repeats.
    assert manifest.strchive is not None
    repeats = next(table.rows for table in manifest.tables if table.table == 'repeats')
    assert manifest.strchive.matched_repeats + manifest.strchive.unmatched_repeats == repeats


# ---- the families sheet ----


def test_the_sheets_columns_are_the_ones_the_workflow_reads() -> None:
    (accepted,) = re.findall(r'\ndef sheetColumns\(\) \{\n\s*return \[([^\]]+)\]\n\}', FAMILY_TABLES)
    assert list(COLUMNS) == _names(accepted)
    # A sheet with any other column is refused, and the workflow reads no other column of a row.
    assert SHEET['other_columns'] == 'refused'
    assert re.search(
        r'def unknown = rows\[0\]\.keySet\(\) - sheetColumns\(\)\n\s+if \(unknown\) \{\n\s+error\(', FAMILY_TABLES
    )
    read = set(re.findall(r'\brow\.(\w+)', FAMILY_TABLES)) - {'containsKey'}
    assert read == set(COLUMNS)
    assert (SHEET['format'], SHEET['header']) == ('csv', True)
    assert 'params.families.splitCsv(header: true)' in FAMILY_TABLES


def test_the_sheets_required_columns_and_values_are_the_ones_the_workflow_asks_for() -> None:
    check = _function(FAMILY_TABLES, 'checkFamily')
    # The columns that need a value in every row, and the one that needs only to be there.
    (valued,) = re.findall(r'\[([^\]]+)\]\.each \{ column ->\s+if \(!row\[column\]\)', check)
    present = re.findall(r"!row\.containsKey\('(\w+)'\)", check)
    assert [name for name, column in COLUMNS.items() if column['value_required']] == _names(valued)
    assert [name for name, column in COLUMNS.items() if column['column_required']] == [*_names(valued), *present]
    # A family is listed once.
    assert 'if (ids.unique(false).size() != ids.size())' in FAMILY_TABLES
    assert 'def ids = rows.collect { row -> checkFamily(row) }' in FAMILY_TABLES
    assert check.rstrip().endswith('return row.family')
    assert [name for name, column in COLUMNS.items() if column.get('unique')] == ['family']


def test_the_sheets_patterns_are_the_ones_the_workflow_matches() -> None:
    assert 'checkIds(row.family, row.proband)' in _function(FAMILY_TABLES, 'checkFamily')
    check = _function(FAMILY_TABLES, 'checkIds')
    (family,) = re.findall(r'if \(!\(family ==~ /(.+?)/\)\)', check)
    # A proband is matched only where the row names one.
    (proband,) = re.findall(r'if \(proband && !\(proband ==~ /(.+?)/\)\)', check)
    patterns = {name: column['pattern'] for name, column in COLUMNS.items() if 'pattern' in column}
    assert patterns == {'family': family, 'proband': proband}


def test_an_empty_proband_reaches_the_tools_as_no_proband() -> None:
    assert 'empty_means' in COLUMNS['proband']
    assert not COLUMNS['proband']['value_required']
    assert 'return proband ? "--proband ${proband}" : \'--no-proband\'' in _function(FAMILY_MODULES, 'probandOption')
    # The option under which the tools write no role, in the two commands that take a row's proband.
    for command in ('family-check', 'family-tables'):
        stated = f'{command} --family ${{family}} --ped ${{ped}} --vcf ${{vcf}} ${{probandOption(proband)}}'
        assert f'case-annotation {stated}' in FAMILY_MODULES
        options = ['--family', 'f', '--ped', 'p', '--vcf', 'v', '--no-proband', '--pipeline-version', 'p']
        if command == 'family-tables':
            options += ['--annotation', 'a', '--echtvar-values', 'e', '--joined-values', 'j', '--transcript-rules', 'v']
            options += ['--roh-upd', 'r']
            options += ['--schema-version', 's', '--output', 'o']
        assert cli.parser().parse_args([command, *options]).proband is None


def test_the_workflow_reads_the_vcfs_of_the_larger_variants_as_the_sheet_names_them() -> None:
    # One VCF in each of sv_vcf and cnv_vcf, and several in repeat_vcfs, between separators.
    larger = 'indexed([row.sv_vcf]), indexed([row.cnv_vcf]), indexed(repeatVcfs(row))'
    assert f'tuple(row.family, file(row.vcf), file(row.ped), row.proband, {larger})' in FAMILY_TABLES
    separated = {name: column['separator'] for name, column in COLUMNS.items() if 'separator' in column}
    assert separated == {'repeat_vcfs': contract.SEPARATOR}
    assert f"row.repeat_vcfs.split('{contract.SEPARATOR}', -1)" in _function(FAMILY_TABLES, 'repeatVcfs')
    # Each of them with its index beside it, which the sheet's check looks for and a task is given. The family's own
    # VCF is read without one.
    indexed = {name: column['index_suffix'] for name, column in COLUMNS.items() if 'index_suffix' in column}
    assert indexed == dict.fromkeys(('sv_vcf', 'cnv_vcf', 'repeat_vcfs'), contract.INDEX_SUFFIX)
    check = _function(FAMILY_TABLES, 'checkFamily')
    assert '([row.sv_vcf, row.cnv_vcf] + repeatVcfs(row))' in check
    assert f'if (!file("${{vcf}}{contract.INDEX_SUFFIX}").exists())' in check
    beside = f'[file(vcf), file("${{vcf}}{contract.INDEX_SUFFIX}")]'
    assert f'return vcfs.findAll {{ vcf -> vcf }}.collectMany {{ vcf -> {beside} }}' in _function(
        FAMILY_TABLES, 'indexed'
    )


def test_the_example_sheet_meets_the_sheets_rules() -> None:
    text = (CONTRACT / 'examples' / 'families.csv').read_text('utf-8')
    reader = csv.DictReader(io.StringIO(text))
    rows = list(reader)
    assert reader.fieldnames == list(COLUMNS)
    for name, column in COLUMNS.items():
        values = [row[name] for row in rows]
        assert all(values) or not column['value_required']
        assert all(re.fullmatch(column['pattern'], value) for value in values if value and 'pattern' in column)
        assert len(set(values)) == len(values) or not column.get('unique')
    # It shows a row that names everything and one that leaves every optional value empty.
    optional = [name for name, column in COLUMNS.items() if not column['value_required']]
    assert {tuple(bool(row[name]) for name in optional) for row in rows} == {
        (True,) * len(optional),
        (False,) * len(optional),
    }
    # Its first row is the family of the example manifest, whose folder holds the files that the row's columns give.
    manifest = family_tables.FamilyManifest.model_validate_json(
        (CONTRACT / 'examples' / 'manifest.json').read_text('utf-8')
    )
    assert manifest.family == rows[0]['family']
    assert [member.sample_id for member in manifest.members if member.role == 'proband'] == [rows[0]['proband']]
    assert {table.path for table in manifest.tables} | {'manifest.json', manifest.notes.path} == _folder(
        [name for name in COLUMNS if rows[0][name]]
    )
    assert len(rows[0]['repeat_vcfs'].split(COLUMNS['repeat_vcfs']['separator'])) == len(manifest.members)


# ---- the workflows' params ----


def _declared(text: str) -> list[dict[str, Any]]:
    """Each param that a workflow declares, as the contract gives it.

    The `params` block runs from the line `params {` to the next line that is `}` alone, and each declaration has the
    form `name: Type`, `name: Type?` or `name: Type = default`.
    """
    lines = text.splitlines()
    assert lines.count('params {') == 1
    start = lines.index('params {')
    declared = []
    for line in lines[start + 1 : lines.index('}', start)]:
        stripped = line.strip()
        if not stripped or stripped.startswith('//'):
            continue
        match = re.fullmatch(r'(\w+): ([A-Za-z]+(?:<[A-Za-z, ]+>)?)(\?)?(?: = (.+))?', stripped)
        assert match is not None, f'not a param declaration: {line!r}'
        name, kind, optional, default = match.groups()
        assert optional is None or default is None
        param = {'name': name, 'type': kind, 'required': optional is None and default is None}
        if default is not None:
            param['default'] = json.loads(default.replace("'", '"'))
        declared.append(param)
    return declared


def _params(workflow: str) -> list[dict[str, Any]]:
    return WORKFLOWS[workflow]['params']


@pytest.mark.parametrize(('workflow', 'text'), [('family_tables.nf', FAMILY_TABLES), ('annotator.nf', ANNOTATOR)])
def test_the_params_are_the_ones_each_workflow_declares(workflow: str, text: str) -> None:
    declaration = ('name', 'type', 'required', 'default')
    assert [{k: v for k, v in param.items() if k in declaration} for param in _params(workflow)] == _declared(text)


def test_the_reference_has_the_index_the_workflow_looks_for_beside_it() -> None:
    indexed = {param['name']: param['index_suffix'] for param in _params('family_tables.nf') if 'index_suffix' in param}
    (name, suffix), *others = indexed.items()
    assert not others
    beside = f'settings.{name}.resolveSibling("${{settings.{name}.name}}{suffix}")'
    assert re.search(rf'def (\w+) = {re.escape(beside)}\n\s+if \(!\1\.exists\(\)\) \{{\n\s+error\(', FAMILY_TABLES)
    assert FAMILY_TABLES.count('.resolveSibling(') == 1


def test_the_patterns_of_the_family_tables_params_are_the_ones_the_workflow_matches() -> None:
    ((names, name),) = re.findall(r'\[([^\]]+)\]\.each \{ v ->\s+if \(!\(v ==~ /(.+?)/\)\)', FAMILY_TABLES)
    (source,) = re.findall(r'settings\.sources\.each \{ v ->\s+if \(!\(v ==~ /(.+?)/\)\)', FAMILY_TABLES)
    # A setting that a run can leave out is matched where the run gives it.
    stated = dict(re.findall(r'if \(!\(settings\.(\w+) ==~ /(.+?)/\)\) \{\n\s+error\(', FAMILY_TABLES))
    matched = dict.fromkeys(re.findall(r'settings\.(\w+)', names), name) | {'sources': source} | stated
    assert {param['name']: param['pattern'] for param in _params('family_tables.nf') if 'pattern' in param} == matched
    assert set(matched) == {'group', 'catalog', 'sources', 'strchive_release', 'strchive_release_date'}
    assert (stated['strchive_release'], stated['strchive_release_date']) == (contract.NAME, strchive.DATE.pattern)


def test_a_release_of_strchive_comes_with_its_name_and_its_date() -> None:
    ((file, release, date),) = WORKFLOWS['family_tables.nf']['given_together']
    assert file == strchive.NAME
    declared = {param['name']: param for param in _params('family_tables.nf')}
    assert [declared[name]['type'] for name in (file, release, date)] == ['Path', 'String', 'String']
    assert not any(declared[name]['required'] for name in (file, release, date))
    # The workflow refuses the file without the two, and either of the two without the file.
    included = _code(_block(FAMILY_TABLES, 'workflow ANNOTATE_FAMILIES'))
    assert re.search(
        rf'if \(settings\.{file} != null\) \{{\s+'
        rf'if \(settings\.{release} == null \|\| settings\.{date} == null\) \{{\s+error\([^\n]+\)\s+\}}',
        included,
    )
    assert re.search(rf'else if \(settings\.{release} != null \|\| settings\.{date} != null\) \{{\s+error\(', included)
    # The step that writes a family's tables takes the three, and gives the tools an option for each.
    given = f"settings.{file} ?: [],\n        settings.{release} ?: '',\n        settings.{date} ?: '',\n"
    assert given in included
    assert f"path {file}, stageAs: '{file}/*'\n    val {release}\n    val {date}\n" in FAMILY_MODULES
    options = ' '.join(f'--{name.replace("_", "-")} ${{{name}}}' for name in (file, release, date))
    assert f'def ranges = {file} ? "{options}" : \'\'' in FAMILY_MODULES
    # The stub runs give the three together, and each refusal its own message.
    assert f'--{file} /src/stub/strchive.json --{release} v0.0.0 --{date} 2099-01-01' in STUB_RUNS
    assert f"refused 'give {release} and {date}'" in STUB_RUNS
    assert f"refused 'the run gives no {file}'" in STUB_RUNS


def test_the_patterns_of_the_annotator_params_are_the_ones_the_workflow_matches() -> None:
    checks = re.findall(
        r'[\[(](params\.[^\])]+)[\])]\s+\.findAll \{ v -> v != null \}\s+\.each \{ v ->\s+if \(!\(v ==~ /(.+?)/\)\)',
        ANNOTATOR,
    )
    assert len(checks) == 2
    matched = {name: pattern for names, pattern in checks for name in re.findall(r'params\.(\w+)', names)}
    assert {param['name']: param['pattern'] for param in _params('annotator.nf') if 'pattern' in param} == matched
    assert set(matched.values()) == {contract.NAME, contract.ILLUMINA_NAME}


def test_a_run_of_the_family_tables_annotates_or_reuses_as_the_workflow_asks() -> None:
    ways = WORKFLOWS['family_tables.nf']['annotation']
    found = re.search(
        r'if \(settings\.annotated != null\) \{\s+if \(settings\.sources\) \{\s+error\([^\n]+\)\s+\}\s+\}\s+'
        r'else \{\s+if \(!settings\.sources\) \{\s+error\([^\n]+\)\s+\}\s+'
        r'\[([^\]]+)\]\.each \{ name ->\s+if \(settings\[name\] == null\) \{\s+error\(',
        FAMILY_TABLES,
    )
    assert found is not None
    # The entry workflow asks for the Annotator's image, which the config reads and ANNOTATE_FAMILIES does not take.
    (image,) = re.findall(
        r"if \(params\.annotated == null && params\.sources && params\.(\w+) == null\) \{\s+error\('give \1, which",
        FAMILY_TABLES,
    )
    assert ways['annotate'] == {'needs': ['sources', image, *_names(found[1])], 'refuses': ['annotated']}
    assert {key: ways['reuse'][key] for key in ('needs', 'refuses')} == {'needs': ['annotated'], 'refuses': ['sources']}
    assert list(ways) == ['annotate', 'reuse']
    # Each param that a way names is declared, and none of them is required of every run.
    declared = {param['name']: param['required'] for param in _params('family_tables.nf')}
    named = {name for way in ways.values() for name in (*way['needs'], *way['refuses'])}
    assert named <= declared.keys()
    assert not any(declared[name] for name in named)


def test_a_reuse_reads_the_files_the_workflow_resolves() -> None:
    reads = WORKFLOWS['family_tables.nf']['annotation']['reuse']['reads']
    group = {'group': '${settings.group}'}
    union, config, structural = (read['path'].format(**group) for read in reads)
    assert f'settings.annotated.resolve("{union}")' in FAMILY_TABLES
    assert f"settings.annotated.resolve('{config}')" in FAMILY_TABLES
    # The annotation of the structural variants is named after the group too, and read where a row names a VCF of them.
    assert 'def svName(group) {\n    return "${group}.sv"\n}' in FAMILY_MODULES
    named = 'out/${svName(settings.group)}.json.gz'
    assert named.replace('${svName(settings.group)}', '${settings.group}.sv') == structural
    assert f'settings.annotated.resolve("{named}")' in FAMILY_TABLES
    assert FAMILY_TABLES.count('settings.annotated.resolve(') == len(reads)
    conditions = [read.get('some_row_names_any_of') for read in reads]
    assert conditions[:2] == [None, None]
    # The entry workflow gives the condition as annotates_sv, the second of the three values that the workflow takes.
    some_row = f'rows.any {{ row -> {" || ".join(f"row.{name}" for name in conditions[2])} }}'
    assert f'ANNOTATE_FAMILIES(\n        families,\n        {some_row},\n' in FAMILY_TABLES
    assert re.search(r'take:\n    families\n    annotates_sv\n    settings\n', FAMILY_TABLES)
    assert re.search(r'if \(annotates_sv\) \{\s+def reused_sv = settings\.annotated\.resolve\(', FAMILY_TABLES)


# ---- the workflow in another pipeline ----


def test_the_included_workflow_reads_nothing_of_the_launched_pipeline() -> None:
    # Params, the project's folder and the manifest are those of the pipeline that was launched, which is another one
    # where a pipeline includes the workflow. Only the entry workflow's own parts read params: its declaration, the
    # entry workflow, the output block and the check of a sheet's row.
    included = FAMILY_TABLES
    for part in (
        _block(FAMILY_TABLES, 'params'),
        _block(FAMILY_TABLES, 'workflow'),
        _block(FAMILY_TABLES, 'output'),
        _function(FAMILY_TABLES, 'checkFamily'),
    ):
        assert part in included
        included = included.replace(part, '')
    assert '\nworkflow ANNOTATE_FAMILIES {\n' in included
    assert 'checkFamily(' not in included
    for text in (included, FAMILY_MODULES, ANNOTATOR_MODULES):
        assert not re.search(r'\bparams\b|\bprojectDir\b|\bworkflow\.', _code(text))
    # Nothing in the script reads the project's folder or the manifest, the entry workflow included.
    assert not re.search(r'\bprojectDir\b|\bworkflow\.', _code(FAMILY_TABLES))


def test_the_settings_of_the_included_workflow_are_the_entrys_params() -> None:
    (listed,) = re.findall(r'\ndef settingNames\(\) \{\n\s*return \[([^\]]+)\]\n\}', FAMILY_TABLES)
    settings = _names(listed)
    assert len(set(settings)) == len(settings)
    assert set(settings) < {param['name'] for param in _params('family_tables.nf')}
    assert 'settingNames().collectEntries { name -> [name, params[name]] }' in FAMILY_TABLES
    # The workflow refuses any other set of names, and reads each setting.
    included = _code(_block(FAMILY_TABLES, 'workflow ANNOTATE_FAMILIES'))
    assert 'if (settings.keySet() != settingNames().toSet()) {\n        error(' in included
    every_run, annotating = re.findall(r'\[([^\]]+)\]\.each \{ name ->\s+if \(settings\[name\] == null\)', included)
    read = set(re.findall(r'\bsettings\.(\w+)', included)) - {'keySet'}
    assert read | set(_names(every_run)) | set(_names(annotating)) == set(settings)
    # A params block refuses a launch without a required param and fills in a default. A map does neither, so the
    # workflow refuses a null for each setting whose param is required or has a default.
    declared = {param['name']: param for param in _params('family_tables.nf')}
    needed = [name for name in settings if declared[name]['required'] or 'default' in declared[name]]
    assert _names(every_run) == needed
    # The pipeline of the stub runs gives exactly these.
    (given,) = re.findall(r'\n    def settings = \[\n(.*?)\n    \]\n', INCLUDING, flags=re.DOTALL)
    assert re.findall(r'^        (\w+): ', given, flags=re.MULTILINE) == settings


def test_a_family_of_an_including_pipeline_is_checked_as_it_arrives() -> None:
    # The entry workflow checks its sheet before the run starts. The workflow holds each family that it is given to
    # the sheet's patterns, and to the run's steps for the structural variants.
    check = _function(FAMILY_TABLES, 'checkRecord')
    assert 'if (proband == null) {\n        error(' in check
    assert 'checkIds(family, proband)' in check
    assert 'if (!annotates_sv && (sv || cnv)) {\n        error(' in check
    assert check.rstrip().endswith('return record')
    assert 'FAMILY_CHECK(families.map { record -> checkRecord(record, annotates_sv) }, ' in FAMILY_TABLES
    # Nothing is checked of the families together. A collected channel also ends when another error stops the run,
    # and a check of what it then holds would report an error of its own in place of that one.
    assert 'sv_group = sv_sites.toList().map { cut -> svGroup(cut) }' in FAMILY_TABLES
    assert FAMILY_TABLES.count('.toList()') == 1


def test_the_including_pipeline_of_the_stub_runs_publishes_as_the_entry_does() -> None:
    assert _code(_block(INCLUDING, 'output')) == _code(_block(FAMILY_TABLES, 'output'))
    assert 'nextflow -C tests/workflows/include/nextflow.config run tests/workflows/include/main.nf -stub-run' in (
        STUB_RUNS
    )


# ---- the published folders ----


def test_the_folders_of_the_family_tables_are_the_ones_the_workflow_publishes_to() -> None:
    publishes = WORKFLOWS['family_tables.nf']['publishes']
    values = {
        'family': '${family}',
        'group': '${group}',
        'tables_version': '${tablesVersion()}',
        'catalog': '${catalog}',
    }
    # The workflow's own templates, with its expressions in place of the values. Each folder is a function, which the
    # output block of a pipeline that includes the workflow calls too.
    tables_folder = contract._TABLES_FOLDER.format(**values)
    assert f'def tablesFolder(catalog) {{\n    return "{tables_folder}"\n}}' in FAMILY_TABLES
    assert 'return file("${moduleDir}/src/case_annotation/tables-version.txt").text.trim()' in FAMILY_TABLES
    family, union = (
        publishes[name]['folder'].format(**values).replace(tables_folder, '${tablesFolder(catalog)}')
        for name in ('family', 'union')
    )
    assert all('${tablesFolder(catalog)}' in folder for folder in (family, union))
    assert f'def familyFolder(family, catalog) {{\n    return "{family}"\n}}' in FAMILY_TABLES
    assert f'def unionFolder(group, catalog) {{\n    return "{union}"\n}}' in FAMILY_TABLES
    annotator = publishes['annotator']['folder'].format(**values)
    assert f'def annotatorFolder(group, catalog) {{\n    return "{annotator}"\n}}' in FAMILY_TABLES
    # The output block gives each function the run's values.
    assert 'family_tables {\n        path { r -> familyFolder(r.family, params.catalog) }' in FAMILY_TABLES
    assert 'union {\n        path unionFolder(params.group, params.catalog)\n        overwrite true' in FAMILY_TABLES
    assert 'annotator {\n        path annotatorFolder(params.group, params.catalog)' in FAMILY_TABLES


def test_the_stub_runs_whose_logs_are_read_for_tasks_write_one_line_for_each_task() -> None:
    # Nextflow's other log is a progress table, which cuts the name of a process short where it is long, as it is
    # under the workflow's name. A check that looks for a task's whole name then fails, and one that looks for the
    # absence of a name passes whatever ran.
    for run in ('family_run', 'included_run'):
        (body,) = re.findall(rf'\n          {run}\(\) \{{\n(.*?)\n          \}}\n', STUB_RUNS, flags=re.DOTALL)
        assert 'docker run --rm -e NXF_ANSI_LOG=false ' in body


def _named(folder: str) -> set[str]:
    """The run's values that a folder of the contract names in braces."""
    return {name for _text, name, _spec, _conversion in string.Formatter().parse(folder) if name is not None}


def test_the_folder_of_the_annotators_output_does_not_name_the_tables_version() -> None:
    publishes = WORKFLOWS['family_tables.nf']['publishes']
    family, union, annotator = (publishes[name]['folder'] for name in ('family', 'union', 'annotator'))
    # The tables change with their version, and the Annotator's output of a group's unions with a catalog does not. A
    # folder of it that named the version would hide it from a run of the next version, which looks for an annotation
    # to reuse where the contract says one is published, and would then annotate again.
    assert _named(family) == {'family', 'tables_version', 'catalog'}
    assert _named(union) == {'group', 'tables_version', 'catalog'}
    assert _named(annotator) == {'group', 'catalog'}
    # It lies beside the folder of the union's tables, and has the name that annotator.nf gives the folder of an
    # annotation.
    group, _, name = annotator.rpartition('/')
    assert union.rpartition('/')[0] == group
    assert WORKFLOWS['annotator.nf']['modes']['annotate']['publishes']['folder'].rpartition('/')[2] == name


def test_the_stub_runs_look_for_the_folders_that_the_contract_gives() -> None:
    publishes = WORKFLOWS['family_tables.nf']['publishes']
    # What every stub run of family_tables.nf is given, and the tables' version.
    (given,) = re.findall(
        r'nextflow run family_tables\.nf -stub-run -profile stub \\\n((?: +--\w+ \S+ \\\n)+)', STUB_RUNS
    )
    values = dict(re.findall(r'--(\w+) (\S+)', given)) | {'tables_version': family_tables.SCHEMA_VERSION}
    union, annotator = (publishes[name]['folder'].format(**values) for name in ('union', 'annotator'))
    assert f'\n          union={union}\n' in STUB_RUNS
    assert f'\n          annotation={annotator}\n' in STUB_RUNS
    # The first run's files are looked for in both, and a later run reuses what the first one published.
    assert '"${union}"/{' in STUB_RUNS
    assert '"${annotation}"/{run.json,' in STUB_RUNS
    assert '--annotated "/src/${results}/${annotation}"' in STUB_RUNS
    # A family's folder is written out wherever a check looks into one.
    family = publishes['family']['folder'].format(family='fam-a', **values)
    assert family.startswith('annotations/fam-a/')
    assert set(re.findall(r'annotations/fam-[^/\s"]+/([^/\s"]+)', STUB_RUNS)) == {family.rpartition('/')[2]}


# The key of a file that a family has by its pedigree: a child whose parents are both samples of its VCF.
BY_PEDIGREE = 'some_child_has_both_parents_sequenced'


def _folder(named: Sequence[str], child: bool = True) -> set[str]:
    """The files of a family's folder by the contract, for a row that names a VCF in the columns given.

    `child` says whether a sample of the family's VCF is the child of two others.
    """
    return {
        file['name']
        for file in FAMILY_FILES
        if (
            child
            if BY_PEDIGREE in file
            else 'row_names_any_of' not in file or set(file['row_names_any_of']) & set(named)
        )
    }


def test_the_files_of_a_familys_folder_are_the_ones_the_workflow_publishes() -> None:
    always = [file['name'] for file in FAMILY_FILES if file.keys() == {'name'}]
    larger = [file for file in FAMILY_FILES if 'row_names_any_of' in file]
    trio = [file for file in FAMILY_FILES if BY_PEDIGREE in file]
    assert [*({'name': name} for name in always), *larger, *trio] == FAMILY_FILES
    assert trio == [{'name': f'{contract.TRIO_TABLE}.parquet', BY_PEDIGREE: True}]
    stems = [file['name'].removesuffix('.parquet') for file in (*larger, *trio)]
    assert (
        f'{", ".join(f"path({name!r})" for name in always)}, '
        f"path('{{{','.join(stems)}}}.parquet', arity: '0..{len(stems)}')"
    ) in FAMILY_MODULES
    names = {'manifest.json': 'manifest', family_tables.NOTES: 'notes'}
    published = ', '.join(names.get(name) or name.removesuffix('.parquet') for name in always)
    assert f'[family: family, files: [{published}] + larger]' in FAMILY_TABLES
    # The stub knows no pedigree, so it gives the table of none.
    (stub,) = re.findall(r'\n    touch variants\.parquet [^\n]*\n', FAMILY_MODULES)
    assert stub == f"\n    touch {' '.join(always)} ${{larger.join(' ')}}\n"
    # The task's inputs by the sheet's columns, and the tables that the stub gives for them.
    inputs = {'sv_vcf': 'sv', 'cnv_vcf': 'cnv', 'repeat_vcfs': 'repeats'}
    by_columns: dict[tuple[str, ...], list[str]] = {}
    for file in larger:
        by_columns.setdefault(tuple(file['row_names_any_of']), []).append(file['name'])
    assert set().union(*by_columns) == set(inputs)
    assert (
        ' + '.join(
            f'({" + ".join(f"vcfs({inputs[column]})" for column in columns)} ? {names} : [])'
            for columns, names in by_columns.items()
        )
        in FAMILY_MODULES
    )


# A row by the columns in which it names a VCF, and whether the family's child has both parents sequenced.
FOLDERS = {
    'a row that names no VCF of the larger variants': ((), True),
    'a structural-variant VCF': (('sv_vcf',), True),
    'a copy-number VCF': (('cnv_vcf',), True),
    'repeat VCFs': (('repeat_vcfs',), True),
    'every VCF of the larger variants': (('sv_vcf', 'cnv_vcf', 'repeat_vcfs'), True),
    'a child whose father is not sequenced': ((), False),
    'such a child with every VCF of the larger variants': (('sv_vcf', 'cnv_vcf', 'repeat_vcfs'), False),
}
# The invented trio's PED file with the father as an individual who is no parent of the child.
_NO_FATHER = 'fam-1\tS1\t0\t0\t1\t1\nfam-1\tS2\t0\t0\t2\t1\nfam-1\tS3\t0\tS2\t2\t2\n'


def _assemble(
    tmp_path: pathlib.Path,
    named: Sequence[str],
    sources: Sequence[str] = test_family_tables._SOURCES,
    scores: Sequence[echtvar_values.Source] = test_family_tables._SCORE_SOURCES,
    joined: Sequence[joined_values.Source] = test_family_tables._JOINED_SOURCES,
    ranges: bool = True,
    child: bool = True,
) -> family_tables.FamilyManifest:
    """Writes the folder of the invented family of test_family_tables, for a row that names the columns given.

    Its run named the sources given, read the archives of the scores given, joined the files given and, unless `ranges`
    is false, gave a release of STRchive. Unless `child` is false, its child has both parents sequenced.
    """
    larger = test_family_tables._larger(
        tmp_path,
        sv=structural_vcfs.SV_RECORDS if 'sv_vcf' in named else None,
        alleles=structural_vcfs.SV_ALLELES if 'sv_vcf' in named else None,
        cnv=structural_vcfs.CNV_RECORDS if 'cnv_vcf' in named else None,
        repeats='repeat_vcfs' in named,
        sources=sources,
        ranges=ranges,
    )
    inputs = test_family_tables._inputs(
        tmp_path,
        info=test_family_tables._info(*sources),
        scores=scores,
        joined=joined,
        ped=test_family_tables._PED if child else _NO_FATHER,
    )
    return family_tables.assemble(**{**inputs, **larger})


@pytest.mark.parametrize('case', FOLDERS)
def test_a_familys_folder_holds_the_files_the_contract_gives_for_its_row(tmp_path: pathlib.Path, case: str) -> None:
    named, child = FOLDERS[case]
    manifest = _assemble(tmp_path, named, child=child)
    assert {path.name for path in (tmp_path / 'out').iterdir()} == _folder(named, child)
    assert (f'{contract.TRIO_TABLE}.parquet' in _folder(named, child)) == child == bool(manifest.roh_upd.upd_children)
    assert [f'{table.table}.parquet' for table in manifest.tables] == [table.path for table in manifest.tables]


def _mode_branches() -> dict[str, str]:
    """The source of each mode of annotator.nf, by the mode: from its test of `params.mode` to the next one."""
    starts = [(found[1], found.start()) for found in re.finditer(r"if \(params\.mode == '(\w+)'\) \{", ANNOTATOR)]
    ends = [*(start for _mode, start in starts[1:]), ANNOTATOR.index('error("mode ${params.mode}:')]
    return {mode: ANNOTATOR[start:end] for (mode, start), end in zip(starts, ends, strict=True)}


def test_each_mode_of_the_annotator_needs_the_params_the_workflow_requires() -> None:
    modes = WORKFLOWS['annotator.nf']['modes']
    branches = _mode_branches()
    assert list(modes) == list(branches)
    assert f'give {", ".join(list(modes)[:-1])} or {list(modes)[-1]}' in ANNOTATOR
    assert re.findall(r'require\(params\.(\w+)', _function(ANNOTATOR, 'requireCatalog')) == ['catalog']
    declared = {param['name']: param for param in _params('annotator.nf')}
    for mode, branch in branches.items():
        # What the branch requires outright, in its order: the catalog, each param that it must be given, and each
        # that an empty value of stops.
        outright = (
            r"\n        (?:(requireCatalog)\(\)|require\(params\.(\w+), '\2'\)|if \(!params\.(\w+)\) \{\n +error\()"
        )
        needs = ['catalog' if catalog else given or valued for catalog, given, valued in re.findall(outright, branch)]
        assert modes[mode]['needs'] == needs
        # And what it requires with another param.
        nested = re.findall(r"\n        if \(params\.(\w+)\) \{\n            require\(params\.(\w+), '\2'\)", branch)
        assert modes[mode].get('needs_with', {}) == {given: [needed] for given, needed in nested}
        assert all(not declared[name]['required'] for name in needs)


def test_the_folders_of_the_annotator_are_the_ones_the_workflow_publishes_to() -> None:
    modes = WORKFLOWS['annotator.nf']['modes']
    values = {
        'date': '${workflow.start.toLocalDate()}',
        'annotator_version': '${params.annotator_version}',
        'catalog': '${params.catalog}',
        'assembly': '${params.assembly}',
        'name': '${params.name}',
        'label': '${params.label}',
    }
    listed = modes['list']['publishes']['folder'].format(**values)
    assert f'lists {{\n        path "{listed}"' in ANNOTATOR
    # A download publishes its data directory, a folder named after the assembly, and its log into one folder.
    download = modes['download']['publishes']
    for output in ('bundle', 'download_log'):
        assert f'{output} {{\n        path "{download["folder"].format(**values)}"' in ANNOTATOR
    assert download['data_directory'] == f'{download["folder"]}/{{assembly}}'
    assert 'path assembly, emit: bundle' in ANNOTATOR_MODULES
    annotate = modes['annotate']['publishes']
    labelled = "${params.label ? '-' + params.label : ''}"
    assert f'annotation {{\n        path "{annotate["folder"].format(**values)}{labelled}"' in ANNOTATOR
    assert annotate['folder_with_label'] == f'{annotate["folder"]}-{{label}}'


# ---- the tables ----


def _arrow(described: Mapping[str, Any]) -> pa.DataType:
    """The Arrow type that the contract describes, as a reader of tables.json builds it."""
    if described['type'] == 'list':
        element = described['element']
        return pa.list_(pa.field('element', _arrow(element), nullable=element['nullable']))
    if described['type'] == 'struct':
        return pa.struct([pa.field(f['name'], _arrow(f), nullable=f['nullable']) for f in described['fields']])
    return pa.type_for_alias(described['type'])


def _schema(table: Mapping[str, Any], samples: Sequence[str], sources: Sequence[str]) -> pa.Schema:
    """The schema of a table's file by the contract, for a family's samples and the sources that its run named."""
    columns = [column for column in table['columns'] if 'source' not in column or column['source'] in sources]
    members = [{'name': sample, **table['member']} for sample in samples if 'member' in table]
    return pa.schema([pa.field(c['name'], _arrow(c), nullable=c['nullable']) for c in (*columns, *members)])


_EVERY_FILE = test_family_tables._JOINED_SOURCES
# The sources that a run named, the archives of scores that it read, the files that it joined, and whether it gave a
# release of STRchive.
RUNS = {
    'every source, every archive, every file and STRchive': (
        test_family_tables._SOURCES,
        test_family_tables._SCORE_SOURCES,
        _EVERY_FILE,
        True,
    ),
    'no dated source of the annotation': (('Ensembl',), test_family_tables._SCORE_SOURCES, _EVERY_FILE, True),
    'no archive of scores': (test_family_tables._SOURCES, (), _EVERY_FILE, True),
    'the scores alone': (test_family_tables._SOURCES, test_family_tables._SCORE_SOURCES[:1], _EVERY_FILE, True),
    'no joined file': (test_family_tables._SOURCES, test_family_tables._SCORE_SOURCES, (), True),
    'one joined file': (test_family_tables._SOURCES, test_family_tables._SCORE_SOURCES, (joined_values.MCR,), True),
    'no release of STRchive': (test_family_tables._SOURCES, test_family_tables._SCORE_SOURCES, _EVERY_FILE, False),
    'nothing that a run can leave out': (('Ensembl',), (), (), False),
}


@pytest.mark.parametrize('case', RUNS)
def test_the_columns_of_each_table_are_the_ones_its_file_declares(tmp_path: pathlib.Path, case: str) -> None:
    sources, scores, joined, ranges = RUNS[case]
    manifest = _assemble(tmp_path, FOLDERS['every VCF of the larger variants'][0], sources, scores, joined, ranges)
    assert [table.table for table in manifest.tables] == list(TABLES)
    samples = [member.sample_id for member in manifest.members if member.sequenced]
    # A run names a dated source of the annotation in its sources, and an archive, a joined file or a release of
    # STRchive by giving it.
    named = (
        *sources,
        *(source.name for source in scores),
        *(source.name for source in joined),
        *((strchive.NAME,) if ranges else ()),
    )
    for name, table in TABLES.items():
        written = pq.read_schema(tmp_path / 'out' / table['file'])
        assert written.equals(_schema(table, samples, named), check_metadata=False), name
    # The columns that the contract gives a source are the ones the manifest lists as left out.
    for name, omitted in (
        ('variants', manifest.omitted_columns),
        ('transcripts', manifest.transcripts_omitted_columns),
        ('sv', manifest.sv_omitted_columns),
        ('repeats', manifest.repeats_omitted_columns),
    ):
        left_out = [c['name'] for c in TABLES[name]['columns'] if 'source' in c and c['source'] not in named]
        assert tuple(sorted(left_out)) == tuple(sorted(omitted))
    assert bool(manifest.omitted_columns) == (
        (sources, scores, joined) != (test_family_tables._SOURCES, test_family_tables._SCORE_SOURCES, _EVERY_FILE)
    )
    assert bool(manifest.transcripts_omitted_columns) == (joined != _EVERY_FILE)
    assert bool(manifest.sv_omitted_columns) == (sources == ('Ensembl',))
    assert bool(manifest.repeats_omitted_columns) == (manifest.strchive is None) == (not ranges)
    with_sources = {'variants', 'transcripts', 'sv', 'repeats'}
    assert all('source' not in c for name in set(TABLES) - with_sources for c in TABLES[name]['columns'])
    # The columns of STRchive are those that its release gives, and no other column of repeats has a source.
    assert [c['name'] for c in TABLES['repeats']['columns'] if 'source' in c] == list(strchive.COLUMNS)
    assert {c['source'] for c in TABLES['repeats']['columns'] if 'source' in c} == {strchive.NAME}
    # An archive's columns are those of its source, and the archive that every run reads is no column's source.
    by_source: dict[str, list[str]] = {}
    for column in TABLES['variants']['columns']:
        by_source.setdefault(column.get('source', ''), []).append(column['name'])
    for source in echtvar_values.SOURCES.values():
        expected = [] if source is echtvar_values.REQUIRED else list(source.columns)
        assert by_source.get(source.name, []) == expected
    # A joined file's columns are those of its source, in variants and in transcripts, and a run can leave each file
    # out. The transcript rules read the catalogue of branch points as well, and their columns of the gained `AG`s are
    # left out with it.
    by_transcript: dict[str, list[str]] = {}
    for column in TABLES['transcripts']['columns']:
        by_transcript.setdefault(column.get('source', ''), []).append(column['name'])
    for file in joined_values.SOURCES.values():
        gained = file is joined_values.BRANCH_POINTS
        assert by_source[file.name] == [*file.columns, *(['ag_gain_high_risk'] if gained else [])]
        assert by_transcript.get(file.name, []) == [*file.transcript_columns, *(['ag_gains'] if gained else [])]


def test_a_table_with_the_members_columns_reserves_the_names_of_its_other_columns() -> None:
    with_members = {name: table.columns for name, table in family_tables.TABLES.items() if table.member is not None}
    assert with_members == dict(family_tables.ROW_COLUMNS)
    assert {name for name, table in TABLES.items() if 'member' in table} == set(with_members)
    assert all(not table['member']['nullable'] for table in TABLES.values() if 'member' in table)


def test_each_type_of_the_contract_is_named_as_arrow_names_it() -> None:
    assert all(pa.type_for_alias(name) == arrow for arrow, name in contract._TYPE_NAMES.items())


def test_a_type_that_the_contract_has_no_name_for_is_refused() -> None:
    with pytest.raises(ValueError, match='no name in the contract for the Arrow type int64'):
        contract._column(pa.field('count', pa.int64()))


def test_a_source_of_a_column_that_its_table_lacks_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    table = family_tables.Table(family_tables.TABLES['repeats'].columns, None, {'clinvar': ('clinvar_accession',)})
    monkeypatch.setattr(family_tables, 'TABLES', {'repeats': table})
    with pytest.raises(ValueError, match=r"repeats: a source is given the columns \['clinvar_accession'\]"):
        contract.tables()
