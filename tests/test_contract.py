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
from case_annotation import cli, contract, family_tables, gnomad_joint

REPOSITORY = pathlib.Path(__file__).parent.parent
CONTRACT = REPOSITORY / 'contract'
FAMILY_TABLES = (REPOSITORY / 'family_tables.nf').read_text('utf-8')
ANNOTATOR = (REPOSITORY / 'annotator.nf').read_text('utf-8')
FAMILY_MODULES = (REPOSITORY / 'modules' / 'family_tables.nf').read_text('utf-8')
ANNOTATOR_MODULES = (REPOSITORY / 'modules' / 'annotator.nf').read_text('utf-8')
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
    named = set(manifest.annotation.sources)
    for omitted, table in ((manifest.omitted_columns, 'variants'), (manifest.sv_omitted_columns, 'sv')):
        sources = family_tables.TABLES[table].sources
        assert omitted == tuple(c for source, columns in sources.items() if source not in named for c in columns)
    # Its archive holds the fields that a run reads gnomAD's joint values from.
    assert gnomad_joint.column_fields(manifest.gnomad_joint.fields, pathlib.Path(manifest.gnomad_joint.archive))
    assert manifest.sv_annotation is not None
    assert (manifest.sv_annotation.catalog, manifest.sv_annotation.sources) == (
        manifest.annotation.catalog,
        manifest.annotation.sources,
    )


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
    check = _function(FAMILY_TABLES, 'checkFamily')
    (family,) = re.findall(r'if \(!\(row\.family ==~ /(.+?)/\)\)', check)
    # A proband is matched only where the row names one.
    (proband,) = re.findall(r'if \(row\.proband && !\(row\.proband ==~ /(.+?)/\)\)', check)
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
            options += ['--annotation', 'a', '--gnomad-joint', 'g', '--schema-version', 's', '--output', 'o']
        assert cli.parser().parse_args([command, *options]).proband is None


def test_the_workflow_reads_the_vcfs_of_the_larger_variants_as_the_sheet_names_them() -> None:
    # One VCF in each of sv_vcf and cnv_vcf, and several in repeat_vcfs, between separators.
    assert 'tuple(row.family, indexed([row.sv_vcf]), indexed([row.cnv_vcf]), indexed(repeatVcfs(row)))' in FAMILY_TABLES
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
    assert 'tuple(row.family, file(row.vcf), file(row.ped), row.proband)' in FAMILY_TABLES


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
    beside = f'params.{name}.resolveSibling("${{params.{name}.name}}{suffix}")'
    assert re.search(rf'def (\w+) = {re.escape(beside)}\n\s+if \(!\1\.exists\(\)\) \{{\n\s+error\(', FAMILY_TABLES)
    assert FAMILY_TABLES.count('.resolveSibling(') == 1


def test_the_patterns_of_the_family_tables_params_are_the_ones_the_workflow_matches() -> None:
    ((names, name),) = re.findall(r'\[([^\]]+)\]\.each \{ v ->\s+if \(!\(v ==~ /(.+?)/\)\)', FAMILY_TABLES)
    (source,) = re.findall(r'params\.sources\.each \{ v ->\s+if \(!\(v ==~ /(.+?)/\)\)', FAMILY_TABLES)
    matched = dict.fromkeys(re.findall(r'params\.(\w+)', names), name) | {'sources': source}
    assert {param['name']: param['pattern'] for param in _params('family_tables.nf') if 'pattern' in param} == matched
    assert set(matched) == {'group', 'catalog', 'sources'}


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
        r'if \(params\.annotated != null\) \{\s+if \(params\.sources\) \{\s+error\([^\n]+\)\s+\}\s+\}\s+'
        r'else \{\s+if \(!params\.sources\) \{\s+error\([^\n]+\)\s+\}\s+'
        r'\[([^\]]+)\]\.each \{ name ->\s+if \(params\[name\] == null\) \{\s+error\(',
        FAMILY_TABLES,
    )
    assert found is not None
    assert ways['annotate'] == {'needs': ['sources', *_names(found[1])], 'refuses': ['annotated']}
    assert {key: ways['reuse'][key] for key in ('needs', 'refuses')} == {'needs': ['annotated'], 'refuses': ['sources']}
    assert list(ways) == ['annotate', 'reuse']
    # Each param that a way names is declared, and none of them is required of every run.
    declared = {param['name']: param['required'] for param in _params('family_tables.nf')}
    named = {name for way in ways.values() for name in (*way['needs'], *way['refuses'])}
    assert named <= declared.keys()
    assert not any(declared[name] for name in named)


def test_a_reuse_reads_the_files_the_workflow_resolves() -> None:
    reads = WORKFLOWS['family_tables.nf']['annotation']['reuse']['reads']
    group = {'group': '${params.group}'}
    union, config, structural = (read['path'].format(**group) for read in reads)
    assert f'params.annotated.resolve("{union}")' in FAMILY_TABLES
    assert f"params.annotated.resolve('{config}')" in FAMILY_TABLES
    # The annotation of the structural variants is named after the group too, and read where a row names a VCF of them.
    assert 'def svName(group) {\n    return "${group}.sv"\n}' in FAMILY_MODULES
    named = 'out/${svName(params.group)}.json.gz'
    assert named.replace('${svName(params.group)}', '${params.group}.sv') == structural
    assert f'params.annotated.resolve("{named}")' in FAMILY_TABLES
    assert FAMILY_TABLES.count('params.annotated.resolve(') == len(reads)
    conditions = [read.get('some_row_names_any_of') for read in reads]
    assert conditions[:2] == [None, None]
    assert f'def annotates_sv = rows.any {{ row -> {" || ".join(f"row.{name}" for name in conditions[2])} }}' in (
        FAMILY_TABLES
    )
    assert re.search(r'if \(annotates_sv\) \{\s+def reused_sv = params\.annotated\.resolve\(', FAMILY_TABLES)


# ---- the published folders ----


def test_the_folders_of_the_family_tables_are_the_ones_the_workflow_publishes_to() -> None:
    publishes = WORKFLOWS['family_tables.nf']['publishes']
    values = {
        'family': '${r.family}',
        'group': '${params.group}',
        'tables_version': '${tablesVersion()}',
        'catalog': '${params.catalog}',
    }
    # The workflow's own templates, with its expressions in place of the values.
    tables_folder = contract._TABLES_FOLDER.format(**values)
    assert f'def tablesFolder() {{\n    return "{tables_folder}"\n}}' in FAMILY_TABLES
    assert 'return file("${projectDir}/src/case_annotation/tables-version.txt").text.trim()' in FAMILY_TABLES
    family, union = (
        publishes[name]['folder'].format(**values).replace(tables_folder, '${tablesFolder()}')
        for name in ('family', 'union')
    )
    assert all('${tablesFolder()}' in folder for folder in (family, union))
    assert f'family_tables {{\n        path {{ r -> "{family}" }}' in FAMILY_TABLES
    assert f'union {{\n        path "{union}"' in FAMILY_TABLES
    annotator = publishes['annotator']['folder'].format(**values)
    assert f'annotator {{\n        path "{annotator}"' in FAMILY_TABLES


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


def _folder(named: Sequence[str]) -> set[str]:
    """The files of a family's folder by the contract, for a row that names a VCF in the columns given."""
    return {
        file['name']
        for file in FAMILY_FILES
        if 'row_names_any_of' not in file or set(file['row_names_any_of']) & set(named)
    }


def test_the_files_of_a_familys_folder_are_the_ones_the_workflow_publishes() -> None:
    always = [file['name'] for file in FAMILY_FILES if 'row_names_any_of' not in file]
    larger = [file for file in FAMILY_FILES if 'row_names_any_of' in file]
    assert [*({'name': name} for name in always), *larger] == FAMILY_FILES
    stems = [file['name'].removesuffix('.parquet') for file in larger]
    assert (
        f'{", ".join(f"path({name!r})" for name in always)}, '
        f"path('{{{','.join(stems)}}}.parquet', arity: '0..{len(larger)}')"
    ) in FAMILY_MODULES
    assert '[family: family, files: [variants, transcripts, manifest, notes] + larger]' in FAMILY_TABLES
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


FOLDERS = {
    'a row that names no VCF of the larger variants': (),
    'a structural-variant VCF': ('sv_vcf',),
    'a copy-number VCF': ('cnv_vcf',),
    'repeat VCFs': ('repeat_vcfs',),
    'every VCF of the larger variants': ('sv_vcf', 'cnv_vcf', 'repeat_vcfs'),
}


def _assemble(
    tmp_path: pathlib.Path, named: Sequence[str], sources: Sequence[str] = test_family_tables._SOURCES
) -> family_tables.FamilyManifest:
    """Writes the folder of the invented family of test_family_tables, for a row that names the columns given."""
    larger = test_family_tables._larger(
        tmp_path,
        sv=structural_vcfs.SV_RECORDS if 'sv_vcf' in named else None,
        alleles=structural_vcfs.SV_ALLELES if 'sv_vcf' in named else None,
        cnv=structural_vcfs.CNV_RECORDS if 'cnv_vcf' in named else None,
        repeats='repeat_vcfs' in named,
        sources=sources,
    )
    inputs = test_family_tables._inputs(tmp_path, info=test_family_tables._info(*sources))
    return family_tables.assemble(**{**inputs, **larger})


@pytest.mark.parametrize('case', FOLDERS)
def test_a_familys_folder_holds_the_files_the_contract_gives_for_its_row(tmp_path: pathlib.Path, case: str) -> None:
    manifest = _assemble(tmp_path, FOLDERS[case])
    assert {path.name for path in (tmp_path / 'out').iterdir()} == _folder(FOLDERS[case])
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


@pytest.mark.parametrize('sources', [test_family_tables._SOURCES, ('Ensembl',)], ids=['every source', 'no dated one'])
def test_the_columns_of_each_table_are_the_ones_its_file_declares(
    tmp_path: pathlib.Path, sources: Sequence[str]
) -> None:
    manifest = _assemble(tmp_path, FOLDERS['every VCF of the larger variants'], sources)
    assert [table.table for table in manifest.tables] == list(TABLES)
    samples = [member.sample_id for member in manifest.members if member.sequenced]
    for name, table in TABLES.items():
        written = pq.read_schema(tmp_path / 'out' / table['file'])
        assert written.equals(_schema(table, samples, sources), check_metadata=False), name
    # The columns that the contract gives a source are the ones the manifest lists as left out.
    for name, omitted in (('variants', manifest.omitted_columns), ('sv', manifest.sv_omitted_columns)):
        left_out = [c['name'] for c in TABLES[name]['columns'] if 'source' in c and c['source'] not in sources]
        assert tuple(sorted(left_out)) == tuple(sorted(omitted))
        assert bool(left_out) == (sources == ('Ensembl',))
    assert all('source' not in c for name in set(TABLES) - {'variants', 'sv'} for c in TABLES[name]['columns'])


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
