# case-annotation-pipeline

Annotates the variant calls of one case into Parquet tables that an analyst or an agent can query directly. A case is
one family: its members and everything sequenced for them.

The pipeline takes a sheet with one row for each family. A row names the family's joint VCF and its pedigree file and,
where the family has them, its VCFs of structural variants, copy number segments and tandem repeats. For each family the
pipeline writes one folder of tables, with a manifest that describes them.

The pipeline's own steps are the Python package `case-annotation-tools` in `src/case_annotation/`, which installs the
`case-annotation` command with one subcommand for each step. The steps run in the tools image, which
[containers/tools](containers/tools/README.md) says how to build. The workflows and the documentation follow in the next
changes.

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

[CLAUDE.md](CLAUDE.md) holds the rules of the repository, and [REVIEW.md](REVIEW.md) says what a review of a change
reports.
