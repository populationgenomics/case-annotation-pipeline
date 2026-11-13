# Review instructions

A Claude review of a change to this repository, most often a pull request, sorts what it finds as below. The rules that
it names are in `CLAUDE.md`.

## Important

The change should not merge until this is fixed. Only these three classes are Important:

- **A public-repository rule is broken**: data from a real case, a name of a case source or cohort, something internal
  to an organisation outside its deployment's own configuration, or a secret.
- **A correctness bug**: code that gives a wrong result or fails, or a validator that accepts invalid input.
- **A hard rule of the repository is broken**: the dependency age gate is bypassed, or a new validator has no test with
  an input that violates it.

## Nit

Everything else worth saying. Report at most three Nits for each review; when there are more, keep the three that matter
most.

## No Important finding

A review that found nothing Important says so. This is a comment and not a GitHub approval: a maintainer still decides
the merge.

## A later review

A pull request can be reviewed more than once, for example after each push to it. A review that finds an earlier review
on the same pull request is a later review, and it reports less:

- It does not repeat a point that is already there.
- It checks each Important finding of the earlier reviews against the change as it now stands. A finding that the change
  still has is still open. The review does not post it again: its summary names it, with its file and line. When the
  author has answered the finding in its thread, the summary says that it was answered, and a maintainer decides it.
- It says that there is no Important finding only when it found none and none is still open.
- It reports no new Nit. The author has answered the first review's Nits, and new ones after each push would never end.

## Not reported

Formatting, prose style, and anything else that `ruff`, `pyright` and `mdformat` already enforce.

## Where the words come from

Important and Nit are the words of the managed Code Review service, so that every kind of Claude review uses the same
two. The service keeps its own format in what it posts: this file changes what a review reports, not how the report
looks.
