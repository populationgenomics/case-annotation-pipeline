from __future__ import annotations

import pathlib

import case_annotation


def test_the_package_declares_its_inline_types() -> None:
    assert (pathlib.Path(case_annotation.__file__).parent / 'py.typed').is_file()
