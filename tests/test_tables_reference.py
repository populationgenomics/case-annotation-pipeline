"""The notes on the tables (tables.md in the package), held to the tables' version, to their size and to the schemas."""

from __future__ import annotations

import importlib.resources
import re
from collections.abc import Callable, Mapping
from typing import Any

import pytest

import case_annotation
from case_annotation import contract, family_tables

NOTES = importlib.resources.files(case_annotation).joinpath(family_tables.REFERENCE).read_text('utf-8')

# The most that the file holds, in characters: a reader takes it whole at the start of its work on a family, so it
# holds what a query gets wrong without an error and no list of columns.
LIMIT = 8_000
# How the notes write a member's column, which each family names as its sample.
MEMBER = '<member>'

# A name written with its table, or with MEMBER alone, wherever it stands between backticks: `table.column`,
# `table.column.field` for a field of a struct or of a list's entries, `table.<member>.field` and `<member>.field` for
# a field of a member's column. A `*` ends a name that stands for every name it starts.
_TABLES = '|'.join(sorted(family_tables.TABLES, key=len, reverse=True))
_NAME = re.compile(rf'(?<![\w.>])(?:(?P<table>{_TABLES})\.|(?={re.escape(MEMBER)}\.))(?P<path>[\w<>*]+(?:\.[\w*]+)*)')


def _fields(described: Mapping[str, Any]) -> dict[str, Any] | None:
    """The fields of a struct, or of the entries of a list of structs, by name; None for any other type."""
    held = described['element'] if described['type'] == 'list' else described
    return {field['name']: field for field in held['fields']} if held['type'] == 'struct' else None


def by_the_schemas() -> dict[str, dict[str, Any]]:
    """Each table's columns by name, as contract.tables describes them, with its members' column under MEMBER."""
    tables = {}
    for table in contract.tables()['tables']:
        columns: dict[str, Any] = {column['name']: column for column in table['columns']}
        if 'member' in table:
            columns[MEMBER] = table['member']
        tables[table['name']] = columns
    return tables


def _held(name: str, among: Mapping[str, Any]) -> bool:
    return any(known.startswith(name[:-1]) for known in among) if name.endswith('*') else name in among


def name_problems(text: str) -> list[str]:
    """The names between backticks that no table has as the notes write them.

    A name written with its table is looked up in that table, and one written with MEMBER alone in the members' columns
    of every table that has them. A bare word that is a column or a field of some table is a problem too: without its
    table nothing can hold it to a schema.
    """
    tables = by_the_schemas()
    every = {name for columns in tables.values() for name in columns} | {
        field for columns in tables.values() for column in columns.values() for field in _fields(column) or ()
    }
    problems = []
    for span in re.findall(r'`([^`]+)`', text):
        if re.fullmatch(r'[a-z_][a-z0-9_]*', span) and span in every and span not in tables:
            problems.append(f'{span} is written without its table')
        for match in _NAME.finditer(span):
            parts = match['path'].split('.')
            if len(parts) > 2:
                problems.append(f'{match[0]} is deeper than a field')
                continue
            column, field = parts[0], (parts[1] if len(parts) == 2 else None)
            candidates = [match['table']] if match['table'] else [t for t in tables if MEMBER in tables[t]]
            if match['table'] is None and field is None:
                continue
            found = False
            for table in candidates:
                columns = tables[table]
                if not _held(column, columns):
                    continue
                if field is None:
                    found = True
                elif not column.endswith('*') and (fields := _fields(columns[column])) is not None:
                    found = found or _held(field, fields)
            if not found:
                problems.append(f'{match[0]} is no column or field of {" or ".join(candidates)}')
    return problems


def test_every_name_the_notes_write_as_a_column_or_a_field_is_one_of_its_table() -> None:
    assert name_problems(NOTES) == []
    # The notes do write such names, of every table.
    named = {match['table'] for span in re.findall(r'`([^`]+)`', NOTES) for match in _NAME.finditer(span)}
    assert named >= set(family_tables.TABLES)


NAMES_WRONG = {
    'a column that the table lacks': (
        '`variants.gene_symbol`',
        'variants.gene_symbol is no column or field of variants',
    ),
    'a column of another table': ('`transcripts.sv_type`', 'transcripts.sv_type is no column or field of transcripts'),
    'a field that the entries lack': ('`variants.genes.biotype`', 'variants.genes.biotype is no column or field'),
    'a field of a column that has none': ('`variants.pos.start`', 'variants.pos.start is no column or field'),
    "a field that the table's members lack": (
        '`repeats.<member>.carries`',
        'repeats.<member>.carries is no column or field of repeats',
    ),
    "a field that no table's members have": ('`<member>.genotype`', '<member>.genotype is no column or field'),
    'a group of columns that the table lacks': ('`sv.gnomad_joint_*`', 'sv.gnomad_joint_* is no column or field of sv'),
    'a name inside a condition': ("`sv.kind = 'DEL'`", 'sv.kind is no column or field of sv'),
    'a column written without its table': ('`zygosity`', 'zygosity is written without its table'),
    'a field of an entry written without its table': (
        '`ensembl_gene_id`',
        'ensembl_gene_id is written without its table',
    ),
    'a name deeper than a field': ('`variants.genes.gene_symbol.first`', 'is deeper than a field'),
}


@pytest.mark.parametrize('case', NAMES_WRONG)
def test_a_name_that_no_table_has_is_caught(case: str) -> None:
    written, expected = NAMES_WRONG[case]
    found = name_problems(f'{NOTES}\n- Also {written}.\n')
    assert len(found) == 1, found
    assert expected in found[0]


def version_problems(text: str, version: str) -> list[str]:
    """What is wrong with the first line of the notes, which is their title and gives the tables' version."""
    first = text.partition('\n')[0]
    return [] if first.startswith('# ') and f'`{version}`' in first else [f'the first line does not give {version}']


def test_the_first_line_gives_the_tables_version() -> None:
    version = importlib.resources.files(case_annotation).joinpath('tables-version.txt').read_text().strip()
    assert version_problems(NOTES, version) == []


WRONG_VERSIONS: Mapping[str, Callable[[str], tuple[str, str]]] = {
    'the notes of other tables': lambda text: (text, '9.9.9'),
    'a version further down': lambda text: (f'# Notes\n\n{text}', family_tables.SCHEMA_VERSION),
}


@pytest.mark.parametrize('case', WRONG_VERSIONS)
def test_notes_that_do_not_open_with_the_tables_version_are_caught(case: str) -> None:
    text, version = WRONG_VERSIONS[case](NOTES)
    assert version_problems(text, version) == [f'the first line does not give {version}']


def too_long(text: str) -> bool:
    return len(text) > LIMIT


def test_the_notes_are_short_enough_to_be_read_whole() -> None:
    assert not too_long(NOTES)
    assert too_long(NOTES + 'x' * (LIMIT - len(NOTES) + 1))
    assert not too_long('x' * LIMIT)


def links_out(text: str) -> list[str]:
    """The links of a text that do not lead to one of its own headings."""
    headings = re.findall(r'^#+ (.+)$', text, flags=re.MULTILINE)
    anchors = {re.sub(r'[^\w\- ]', '', heading).lower().replace(' ', '-') for heading in headings}
    targets = re.findall(r'\]\(([^)]*)\)', text)
    return [target for target in targets if not (target.startswith('#') and target[1:] in anchors)]


def test_the_notes_stand_alone_in_a_familys_folder() -> None:
    # The file is copied into a folder that holds no other document, so a link could only lead to one of its headings.
    assert links_out(NOTES) == []


@pytest.mark.parametrize(
    'target', ['../../docs/family-tables.md', 'manifest.json', '#a-heading-the-file-lacks', 'https://example.org']
)
def test_a_link_out_of_the_notes_is_caught(target: str) -> None:
    assert links_out(f'{NOTES}\nSee [this]({target}).\n') == [target]
