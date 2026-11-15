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


def test_the_pipelines_version_is_the_packages() -> None:
    repository = pathlib.Path(__file__).parent.parent
    config = (repository / 'nextflow.config').read_text('utf-8')
    (manifest,) = re.findall(r'^manifest \{\n(.*?)^\}', config, flags=re.MULTILINE | re.DOTALL)
    (version,) = re.findall(r"^    version = '([^']+)'$", manifest, flags=re.MULTILINE)
    assert version == tomllib.loads((repository / 'pyproject.toml').read_text('utf-8'))['project']['version']
    # A manifest gives the installed package's version, which is the same one where the environment is in step.
    assert version == family_tables.PIPELINE_VERSION
