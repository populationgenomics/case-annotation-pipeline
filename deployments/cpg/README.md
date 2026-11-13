# CPG's deployment

The Centre for Population Genomics (CPG) runs the pipeline through Seqera, in its `rd-benchmark` dataset, launched with
its [analysis-runner](https://github.com/populationgenomics/analysis-runner). This document is the runbook for an
operator at CPG: what is in place, how a workflow is launched, and the params of each workflow with where its output
lands. [Annotation](../../docs/annotation.md) explains what the workflows do.

## What is in place

| Input                                | Where                                                                             |
| ------------------------------------ | --------------------------------------------------------------------------------- |
| the Annotator's image                | `australia-southeast1-docker.pkg.dev/cpg-common/images-dev/dragen`                |
| the tools image                      | `australia-southeast1-docker.pkg.dev/cpg-common/images-dev/case-annotation-tools` |
| the licence key                      | the Secret Manager secret `dragen-api-key` in the dataset's project               |
| the annotation data                  | `gs://cpg-rd-benchmark-test/references/dragen-annotation/4.0.0/202610/GRCh38/`    |
| GRCh38                               | `gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta`           |
| the archive of gnomAD's joint values | `gs://cpg-common-test/gnomad/echtvar_simple/gnomad_4.1_whole_genome.zip`          |
| outputs                              | `gs://cpg-rd-benchmark-test`                                                      |

- The Annotator's image is CPG's DRAGEN image, which holds the Annotator at `/opt/dragen/share/nirvana/Annotator` with
  `gcp-secret` and `jq` ([The Annotator's image](../../docs/annotator-image.md)). It is built by hand, because DRAGEN's
  installer is behind an Illumina login, and runs at the test access level pull it from `images-dev`.

- The tools image is built from [containers/tools/Dockerfile](../../containers/tools/Dockerfile) at a release tag, and
  the Dockerfile shows the command. The params file pins it by digest, so a new tag means a new build and a new digest
  there. Until the image of the first tag is pushed, the line holds a placeholder in place of the digest, and a task
  that needs the tools image fails to pull one: every run of `family_tables.nf`, and an `annotate` run with `regions`.

- The licence is CPG's DRAGEN licence, whose API key the Annotator takes as an API key file. The dataset's accounts can
  read the key, in the project `rd-benchmark-220215`. The dataset's analysis group can add and disable versions, and
  cannot read them:

  ```sh
  gcloud secrets versions add dragen-api-key --project=rd-benchmark-220215 --data-file=-
  ```

- The annotation data is the October 2026 catalog, `202610`, as a `download` run published it
  ([The Annotator](#the-annotator)).

- The archive is CPG's copy of the echtvar archive of gnomAD v4.1's joint release, which has the size of the public file
  ([gnomAD's joint values](../../docs/annotation.md#gnomads-joint-values)). Each family's manifest records the SHA-256
  that its run read.

[test.annotator.params.yaml](test.annotator.params.yaml) holds the images, the Annotator's path and version, the
licence's secret and the output bucket for runs at the test access level, with the images pinned by digest. A run copies
it and adds its own params. No other access level has a params file yet.

## Launching a workflow

A Seqera launch takes one params file and one entry script, at one commit:

```sh
analysis-runner seqera --dataset rd-benchmark --access-level test \
    --repository https://github.com/populationgenomics/case-annotation-pipeline \
    --main-script <annotator.nf or family_tables.nf> --commit-id <commit> \
    --params-file <run>.params.yaml
```

The analysis-runner launches with Nextflow's v2 syntax parser, which the pipeline needs. Access levels above `test` need
the commit to be on `main`. The dataset has to allow this repository for the analysis-runner, and Seqera has to be able
to clone the repository while it is private.

Where the analysis-runner is not installed, `uvx --exclude-newer '7 days' analysis-runner@3.5.0` runs it in place of
`analysis-runner` in these commands. The option is the repository's dependency age gate: `uvx` does not read `[tool.uv]`
in `pyproject.toml`, so without it the tool's own dependencies would resolve to releases younger than seven days.

A run's input files are given as `gs://` paths that the dataset's accounts can read.

## The Annotator

`annotator.nf` runs the Annotator of the DRAGEN image in one of three modes ([Modes](../../docs/annotation.md#modes)).
Only the dataset's accounts can read the licence key, so the workflow runs from Seqera. A run copies the params file,
adds its mode and its own params, and launches:

```sh
cp deployments/cpg/test.annotator.params.yaml <run>.params.yaml
printf 'mode: list\n' >> <run>.params.yaml
analysis-runner seqera --dataset rd-benchmark --access-level test \
    --repository https://github.com/populationgenomics/case-annotation-pipeline \
    --main-script annotator.nf --commit-id <commit> \
    --params-file <run>.params.yaml
```

The table gives, for each mode, the params a run adds and where its output lands.

| Mode       | Params the run adds                                      | Output                                                                            |
| ---------- | -------------------------------------------------------- | --------------------------------------------------------------------------------- |
| `list`     | none                                                     | `gs://cpg-rd-benchmark-test/annotation/list/<date>/`                              |
| `download` | `catalog`, `exclude_sources`                             | `gs://cpg-rd-benchmark-test/references/dragen-annotation/4.0.0/<catalog>/GRCh38/` |
| `annotate` | `catalog`, `annotation_data`, `input`, `name`, `sources` | `gs://cpg-rd-benchmark-test/annotations/<name>/illumina-<catalog>/`               |

A `list` run writes Illumina's listing of the catalogs it hosts and of each source's versions into one folder for each
day, which is where a new catalog first shows.

The catalog `202610` was downloaded with the germline sources of the basic and premium categories, leaving out COSMIC's
somatic data:

```yaml
mode: download
catalog: '202610'
exclude_sources: [cosmic]
```

A second download of a catalog replaces the data directory of the first, and an annotation that reads the folder
meanwhile finds it incomplete. The `runConfigs/bundle.json` of the downloaded data lists the sources it holds:

```sh
gcloud storage cat gs://cpg-rd-benchmark-test/references/dragen-annotation/4.0.0/202610/GRCh38/runConfigs/bundle.json \
    | jq -c '.data.versionSelector | keys'
```

One VCF is annotated with the sources that the run names:

```yaml
mode: annotate
catalog: '202610'
annotation_data: gs://cpg-rd-benchmark-test/references/dragen-annotation/4.0.0/202610/GRCh38
input: gs://<bucket>/<path>/<name>.vcf.gz
name: <name>
sources: [<source>, <source>, ...]
```

A run over a few regions takes `regions`, a BED file in a bucket, and a `label` for its folder. It cuts the VCF in the
tools image, so it needs the digest in the params file.

## Family tables

`family_tables.nf` annotates the union of a group of families' sites once and writes a folder of Parquet tables for each
family ([Family tables](../../docs/annotation.md#family-tables)). This repository has no launcher of its own for it. A
run needs two things from whoever launches it: a families sheet and a params file.

The sheet is a CSV in a bucket with one row for each family, which names the family's VCF, its PED file and its proband,
and the VCFs of its larger variants where it has them
([the input contract](../../docs/annotation.md#the-input-contract)).
[contract/families-sheet.json](../../contract/families-sheet.json) states its columns and rules for a program that
writes one, and [contract/examples/families.csv](../../contract/examples/families.csv) is an invented one.

The params are the ones that [contract/workflow.json](../../contract/workflow.json) lists for `family_tables.nf`. Every
line of `test.annotator.params.yaml` is one of them, and its tools image needs its digest filled in. A run adds CPG's
reference and archive, and its own sheet, group and catalog:

```sh
cp deployments/cpg/test.annotator.params.yaml <run>.params.yaml
cat >> <run>.params.yaml <<'EOF'
reference: gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta
gnomad_joint: gs://cpg-common-test/gnomad/echtvar_simple/gnomad_4.1_whole_genome.zip
families: gs://<bucket>/<path>/families.csv
group: <group>
catalog: '202610'
annotation_data: gs://cpg-rd-benchmark-test/references/dragen-annotation/4.0.0/202610/GRCh38
sources: [<source>, <source>, ...]
EOF
analysis-runner seqera --dataset rd-benchmark --access-level test \
    --repository https://github.com/populationgenomics/case-annotation-pipeline \
    --main-script family_tables.nf --commit-id <commit> \
    --params-file <run>.params.yaml
```

> [!WARNING]
>
> A run with `sources` starts the Annotator, which is billed: once over the group's alleles, and a second time over its
> structural variants where a row names an `sv_vcf` or a `cnv_vcf`.

With `sv_vcf` or `cnv_vcf` in a sheet, `sources` should name the sources that hold structural variants (`gnomad`,
`oneKg` and `clinvar`).

A run that failed after the Annotator finished is launched again with `annotated` in place of `sources`, which reads the
annotation that the first run published and starts no Annotator:

```yaml
annotated: gs://cpg-rd-benchmark-test/annotation-unions/<group>/tables-<tables version>-illumina-<catalog>/annotator
```

That folder can be reused when it holds `out/<group>.json.gz` and `run.json`, and `out/<group>.sv.json.gz` as well when
a family of the group names an `sv_vcf` or a `cnv_vcf`.

| Output                                          | Where                                                                                              |
| ----------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| each family's folder                            | `gs://cpg-rd-benchmark-test/annotations/<family>/tables-<tables version>-illumina-<catalog>/`      |
| the union, its annotation and its gnomAD values | `gs://cpg-rd-benchmark-test/annotation-unions/<group>/tables-<tables version>-illumina-<catalog>/` |

The union's folder holds every family's alleles and records, so only a family's own folder is given to whoever analyses
that family ([What a run publishes](../../docs/annotation.md#what-a-run-publishes)).
