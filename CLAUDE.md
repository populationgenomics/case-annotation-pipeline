# Rules for this repository

## Commands

```sh
uv sync --locked --group lint
uv run pre-commit install
uv run pytest
uv run pytest tests/test_package.py
uv run pytest tests/test_package.py::test_the_package_declares_its_inline_types
uv run pre-commit run --all-files
```

## This repository will be public

The repository is private for now and will be made public later. Publishing exposes its whole history, so every file,
commit message, branch name, pull request, issue and review comment written now has to be fit to publish. Rewriting
history before publication would not remove anything: pull requests keep their own commits on GitHub, along with their
review comments.

- No patient data, and nothing derived from a real case: no real variant tied to a phenotype, no case or sample
  accession, no text from a real referral or report. Fixtures are invented and say so.
- No names of case sources, cohorts or programmes, and nothing internal to any organisation that runs the pipeline:
  projects, datasets, systems, infrastructure (buckets, cloud projects, hosts), private repositories, documents or
  discussions. The tools take such values as data; they never appear as literals, examples, comments or test values.
- The exception is a deployment's own configuration: `deployments/<organisation>/` may name that organisation's own
  infrastructure (cloud projects, buckets, registries, regions, groups) and the public data it draws on. Nothing in it
  names an individual or holds patient data, details of a real case or a secret.
- Give the reasoning in the text itself instead of pointing at a private document, thread or ticket.
- When unsure whether something may be public, leave it out and ask.

## Design rule

- A table's bytes change only with the tables' version, the one line of `src/case_annotation/tables-version.txt`. The
  version changes when a column of a table is renamed, retyped or removed, or when its meaning changes. A table, a
  column or a manifest key that is added leaves the version as it is, because a reader of the Parquet tables by column
  name, or of the manifest as plain JSON, finds what it knew as before.

## Working norms

- Push back when a design, name or approach seems worse, with reasoning; offer a better alternative with its trade-offs.
  The author makes the final call.
- Investigate before producing: read the code and verify constraints rather than assuming a convention applies.
- Fail loudly: raise on missing or malformed input; never fall back to a default to carry on.
- Explain a non-obvious change before showing the diff. Ask when unsure rather than assume intent.

## Code style

- The package supports Python 3.12 and later; `from __future__ import annotations` in every module.
- Import modules, not symbols (`import pathlib`, then `pathlib.Path`), except names from `typing` and `collections.abc`.
  No `TYPE_CHECKING` blocks.
- Models are frozen (`frozen=True` on every class, which pyright checks) and use tuples for sequences.
- Comments only for a non-obvious mechanism or constraint; no history narration.
- Tests assert invariants and test each validator with an input that violates it, not today's values.
- Never bypass the dependency age gate (`exclude-newer`, Dependabot cooldown), and give no package an exception to it.
- A workflow pins each action to a commit SHA, with the version in a comment (`@<sha> # v6.1.0`), because a tag can be
  moved to other code. Dependabot updates both.

## Pull requests

- A change made by an agent always goes out as a pull request, never as a push to `main`.
- The `lint` and `tests` workflows pass before it merges. A pull request lists their checks as `lint / pre-commit` and
  `tests / pytest`, and the API names the check runs `pre-commit` and `pytest`. They run
  `uv run pre-commit run --all-files` and `uv run pytest`, so run both before pushing.
- An agent reviews its own pull request by `REVIEW.md`, fixes what the review finds and reviews again. It stops after
  two reviews in a row with no Important finding, and then says that the pull request is ready for a maintainer.
- A maintainer decides each merge, and every merge is a squash, so that `main` holds one commit for each pull request.
  An agent merges a pull request, or arms its auto-merge, only when a maintainer asks it directly to merge that pull
  request as it stands. Passing checks and an approval are not that request, and the request does not replace the checks
  above.
- Before reviewing a change, read `REVIEW.md`: what a review reports as Important, what as a Nit, and what not at all.
