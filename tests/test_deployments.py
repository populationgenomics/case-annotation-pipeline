"""Each deployment's params files (deployments/<organisation>/), held to the params that their workflows declare."""

from __future__ import annotations

import pathlib
import re

import pytest

from case_annotation import contract

REPOSITORY = pathlib.Path(__file__).parent.parent
# A params file is named after its access level and its workflow, as test.annotator.params.yaml is for annotator.nf.
PARAMS_FILES = sorted((REPOSITORY / 'deployments').glob('*/*.params.yaml'))


def _params(path: pathlib.Path) -> list[tuple[str, str]]:
    """The params of a file that holds one on each line that is no comment, each with its value as written."""
    lines = [line for line in path.read_text('utf-8').splitlines() if line and not line.startswith('#')]
    found = [re.fullmatch(r'(\w+): (\S+)', line) for line in lines]
    assert all(found), f'{path}: a line that is no param with one value'
    return [(match[1], match[2]) for match in found if match is not None]


def test_a_deployment_has_a_params_file() -> None:
    assert PARAMS_FILES


@pytest.mark.parametrize('path', PARAMS_FILES, ids=lambda path: str(path.relative_to(REPOSITORY)))
def test_a_deployments_params_file_gives_params_that_its_workflow_declares(path: pathlib.Path) -> None:
    _level, workflow, _params_suffix, _yaml = path.name.split('.')
    workflows = contract.workflows()['workflows']
    declared = {param['name']: param for param in workflows[f'{workflow}.nf']['params']}
    given = _params(path)
    names = [name for name, _value in given]
    assert len(set(names)) == len(names)
    assert set(names) <= declared.keys()
    # A run of family_tables.nf starts from a copy of the Annotator's file, so that file gives it no param that it
    # does not declare.
    assert workflow != 'annotator' or set(names) <= {param['name'] for param in workflows['family_tables.nf']['params']}
    # A value that names files and folders has the form its workflow accepts.
    assert all(re.fullmatch(declared[name]['pattern'], value) for name, value in given if 'pattern' in declared[name])
