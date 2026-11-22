from __future__ import annotations

import importlib.resources
import pathlib
import re
import tomllib

import pytest

import case_annotation
from case_annotation import cli, family_tables


def test_the_package_declares_its_inline_types() -> None:
    assert (pathlib.Path(case_annotation.__file__).parent / 'py.typed').is_file()


def test_the_command_line_requires_a_subcommand() -> None:
    with pytest.raises(SystemExit):
        cli.main([])


def test_the_tables_version_is_the_one_line_of_its_file() -> None:
    text = importlib.resources.files(case_annotation).joinpath('tables-version.txt').read_text()
    assert re.fullmatch(r'\d+\.\d+\.\d+\n', text)
    assert text.removesuffix('\n') == family_tables.SCHEMA_VERSION


def test_the_pipelines_version_is_the_one_line_of_its_file_and_the_packages() -> None:
    repository = pathlib.Path(__file__).parent.parent
    text = (repository / 'pipeline-version.txt').read_text('utf-8')
    assert re.fullmatch(r'\d+\.\d+\.\d+\n', text)
    version = text.removesuffix('\n')
    assert version == tomllib.loads((repository / 'pyproject.toml').read_text('utf-8'))['project']['version']
    # A manifest gives the installed package's version, which is the same one where the environment is in step.
    assert version == family_tables.PIPELINE_VERSION
    # Nextflow reads the manifest of a project that it pulls from the config's text, with no folder to find the file
    # in, so the config holds the version written out.
    config = (repository / 'nextflow.config').read_text('utf-8')
    (manifest,) = re.findall(r'^manifest \{\n(.*?)^\}', config, flags=re.MULTILINE | re.DOTALL)
    assert re.findall(r"^    version = '([^']+)'$", manifest, flags=re.MULTILINE) == [version]
    # The workflow reads the file from its own folder, which it also finds where another pipeline includes it, and
    # gives the version to the command that writes a family's folder, which refuses another. The workflow's test
    # gives the versions that this commit has.
    workflow = (repository / 'family_tables.nf').read_text('utf-8')
    assert 'def pipelineVersion() {\n    return file("${moduleDir}/pipeline-version.txt").text.trim()\n}' in workflow
    assert '        tablesVersion(),\n        pipelineVersion(),\n    )' in workflow
    modules = (repository / 'modules' / 'family_tables.nf').read_text('utf-8')
    assert modules.count('--pipeline-version ${pipeline_version}') == 2
    test = (repository / 'tests' / 'workflows' / 'family_tables.nf').read_text('utf-8')
    assert f"        '{family_tables.SCHEMA_VERSION}',\n        '{version}',\n    )" in test
    # Every call of the family check in the workflow's test gives the version too.
    given = re.escape(version)
    one_line = re.findall(rf"FAMILY_CHECK\(families\.join\(structural\), '{given}'\)", test)
    refused = re.findall(rf"FAMILY_CHECK_\w+_REFUSED\(\n        channel\.of\(.*\),\n        '{given}',\n    \)", test)
    assert len(one_line) + len(refused) == len(re.findall(r'= FAMILY_CHECK\w*\(', test)) == 5


def test_the_first_step_of_a_run_refuses_tools_of_another_release_before_the_annotator_runs() -> None:
    repository = pathlib.Path(__file__).parent.parent
    workflow = (repository / 'family_tables.nf').read_text('utf-8')
    # FAMILY_CHECK takes the workflow's version, and what the Annotator reads comes from families that it has checked:
    # the union of the sites of the normalised VCFs, and the union of the sites of the larger variants. A run that
    # reuses an annotation starts with the same step.
    order = (
        'checked = FAMILY_CHECK(families.map { record -> checkRecord(record, annotates_sv) }, pipelineVersion())',
        'passed = families.join(checked.map { f -> [f] })',
        'normalised = NORMALISE_FAMILY(\n        passed.map',
        'union = UNION_SITES(settings.group, normalised.sites.flatten().collect())',
        'sv_sites = SV_SITES(\n            passed.filter',
        'sv_union = SV_UNION(settings.group, sv_group.map',
        'annotated = ILLUMINA_ANNOTATE(\n            union.map',
        'sv_annotated = ILLUMINA_ANNOTATE_SV(\n                sv_union.union.map',
    )
    places = [workflow.find(step) for step in order]
    assert -1 not in places, [step for step, place in zip(order, places, strict=True) if place < 0]
    assert places == sorted(places)
    assert 'sv_group = sv_sites.toList()' in workflow
    assert workflow.count('FAMILY_CHECK(') == 1
