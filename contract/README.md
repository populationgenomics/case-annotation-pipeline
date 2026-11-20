# The pipeline's contract

This folder states the pipeline's interface in files that a program can read, so that another repository can test its
own code against the interface without reading this repository's source. Such a repository writes families sheets,
launches the workflows or reads a family's folder. It copies this folder at the tag of the version it runs.

| File                                             | What it states                                                                                                                                                                                 |
| ------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [manifest.schema.json](manifest.schema.json)     | the JSON Schema of a family's `manifest.json`                                                                                                                                                  |
| [families-sheet.json](families-sheet.json)       | the columns of a families sheet in their order, which of them a sheet must have and a row must fill, the patterns of the names, and the rules for the VCFs of the larger variants              |
| [workflow.json](workflow.json)                   | for `family_tables.nf` and `annotator.nf`: each param, the names a run accepts, the folders a run publishes to, the files of a family's folder, and the files a reused annotation is read from |
| [tables.json](tables.json)                       | the columns of the five tables, each with its Arrow type and whether the file declares it nullable                                                                                             |
| [tables-version](tables-version)                 | the tables' version, in one line                                                                                                                                                               |
| [examples/manifest.json](examples/manifest.json) | a manifest as a run writes it                                                                                                                                                                  |
| [examples/families.csv](examples/families.csv)   | a families sheet with two rows                                                                                                                                                                 |

Both examples are invented: the families, the samples, the paths, the sizes and the checksums belong to no real case.

The notes on the tables, [`tables.md`](../src/case_annotation/tables.md), say what the tables do not show and a query
gets wrong without an error, and every family's folder has a copy of them under that name.
[Reading the family tables](../docs/family-tables.md) says what each manifest key means, and
[Annotation](../docs/annotation.md) what each param does.

## Reading the files

The manifest's schema accepts a key that it does not describe. The code that writes a manifest refuses unknown keys, but
a key that is added later keeps the tables' version, so a reader must take a manifest with such a key. Every key that
the schema describes is required: a manifest holds each of them, with an empty list or a null where a family has nothing
for it.

In `families-sheet.json`, `column_required` says whether a sheet's header must hold the column and `value_required`
whether a row must hold a value in it. A `pattern` is matched against the whole value. A column with an `index_suffix`
names VCFs that each have an index beside them, named as the VCF with the suffix. The `separator` stands between the
VCFs of one field.

In `workflow.json`, a param is `required` where its workflow declares it without a default, and a `pattern` on a list
param holds for each of its values. A param with an `index_suffix` names a file that has an index beside it, as in the
sheet. A mode's `needs` are the params that a run of the mode must have a value of beside the required ones, and for a
list at least one value. The params of one list of `given_together` are given together or not at all, as the file of a
release of STRchive is with the release's name and date. The workflows check more than the file states, such as the
catalog alias `LATEST` and the index beside `input` in a run with `regions`;
[Checks that stop a run](../docs/annotation.md#checks-that-stop-a-run) lists them. A folder is relative to `output_dir`,
and the names in braces are the run's values: its params, a family's id from the sheet, the tables' version, and the
date on which a listing run started. A file with `row_names_any_of` is in a family's folder when the family's row names
a VCF in one of those columns. A file with `some_row_names_any_of` is read when any row of the sheet does. The folder of
`annotator` is where a run of `family_tables.nf` that annotates publishes the Annotator's output, and what a later run
of the group with that catalog takes as `annotated`. It names no version of the tables, because the output does not
change with it.

In `tables.json`, a `type` is the name of an Arrow type, `list` with its `element`, or `struct` with its `fields`. The
tables `variants`, `sv` and `repeats` end with one column for each sequenced member, named as the member's sample, and
`member` gives the type of that column once. A column with a `source` holds the data of a source that a run can leave
out, and the table of a run that did not name the source lacks the column. The source is a dated source of the
annotation, under its name in a run config, or an echtvar archive of AlphaGenome's scores, a joined public file or a
release of STRchive, under the name of its param. The manifest lists the columns that a run left out. The file of
`transcripts` declares every column nullable, also the columns that hold a value in every row.

## When the contract changes

`manifest.schema.json` and `tables.json` change with the tables' version: without a new version, no column and no
manifest key is renamed, retyped or removed. A column, a table or a manifest key can be added under the same version,
and then these two files gain it. `families-sheet.json` and `workflow.json` change with the pipeline's version, which
`workflow.json` records as `pipeline_version` and a release tag names.

## Writing the files

`case-annotation contract` writes every file of this folder except this one:

```sh
uv run case-annotation contract --out contract
```

The schema, the tables' columns and the version come from the code that writes a family's folder. What the files say of
the sheet and the workflows is declared in `src/case_annotation/contract.py`, and `tests/test_contract.py` holds each
statement to the workflows' source. The same tests write the folder again and compare it with the committed files, so a
change to the interface that leaves this folder as it was fails the `tests` check.
