# case-annotation-pipeline

Annotates the variant calls of one case into Parquet tables that an analyst or an agent can query directly. A case is
one family: its members and everything sequenced for them.

The pipeline takes a sheet with one row for each family. A row names the family's joint VCF and its pedigree file and,
where the family has them, its VCFs of structural variants, copy number segments and tandem repeats. For each family the
pipeline writes one folder of tables, with a manifest that describes them.

## The two workflows

Both are Nextflow workflows at the root of the repository.

| Workflow                             | What a run does                                                                                                                                                                            |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [annotator.nf](annotator.nf)         | runs Illumina's Annotator in one of three modes: `list` asks Illumina which catalogs and source versions it hosts, `download` fetches one catalog's data, and `annotate` annotates one VCF |
| [family_tables.nf](family_tables.nf) | annotates the alleles of a group of families once, and their structural variants in a second run, and writes each family's folder of tables                                                |

A run of `family_tables.nf` reads the data that a `download` published, so an operator who wants family tables runs
`list`, `download` and then `family_tables.nf`.

The pipeline's own steps are the Python package `case-annotation-tools` in `src/case_annotation/`, which installs the
`case-annotation` command with one subcommand for each step. They run in the image `case-annotation-tools`, which
[containers/tools/Dockerfile](containers/tools/Dockerfile) builds. The Annotator is Illumina's licensed software and
runs in an image that the organisation running the pipeline builds itself.

## The documents

| Document                                             | Its reader                                                                    | What it holds                                                                                                                                  |
| ---------------------------------------------------- | ----------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| [Reading the family tables](docs/family-tables.md)   | an analyst or an agent who queries the tables                                 | one VCF record followed into its rows, the usual rare-disease queries with DuckDB from Python, and the folder's files and manifest             |
| [Notes on the tables](src/case_annotation/tables.md) | an agent that queries a family's tables; a copy lies in every family's folder | what the tables do not show and a query gets wrong without an error                                                                            |
| [Annotation](docs/annotation.md)                     | the operator who runs the pipeline                                            | the Annotator's modes and its licence, which sources are dated, what a families sheet holds, the steps, what a run publishes and what stops it |
| [The Annotator's image](docs/annotator-image.md)     | whoever builds the image of the Annotator                                     | what the image has to hold, and how a task fetches the licence key                                                                             |
| [Glossary](docs/glossary.md)                         | everyone                                                                      | each word that the documents use                                                                                                               |
| [The contract](contract/README.md)                   | a program that launches or reads a run                                        | the manifest's JSON Schema, and the sheet's columns, the workflows' params and the tables' columns as JSON                                     |

`deployments/<organisation>/` holds an organisation's own params files and its runbook.
[deployments/cpg/](deployments/cpg/README.md) is that of the Centre for Population Genomics.

## Development

The Python package and its checks need [uv](https://docs.astral.sh/uv/). Set up the environment and the commit hooks
once:

```sh
uv sync --locked --group lint
uv run pre-commit install
```

A pull request has to pass two checks, which run the same commands as these:

```sh
uv run pre-commit run --all-files
uv run pytest
```

A change to the interface changes [the contract](contract/README.md) too.
`uv run case-annotation contract --out contract` writes it again, and `pytest` fails while the committed folder differs
from what the command writes.

The `nextflow` workflow runs on a pull request too. It runs both workflows in stub mode, and it runs their processes in
Docker on invented inputs, which are the tests in `tests/workflows/`. The comment that opens each test gives its
command.

[CLAUDE.md](CLAUDE.md) holds the rules of the repository, and [REVIEW.md](REVIEW.md) says what a review of a change
reports.
