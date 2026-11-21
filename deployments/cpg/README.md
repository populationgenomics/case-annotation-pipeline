# CPG's deployment

The Centre for Population Genomics (CPG) runs the pipeline through Seqera, in its `rd-benchmark` dataset, launched with
its [analysis-runner](https://github.com/populationgenomics/analysis-runner). This document is the runbook for an
operator at CPG: what is in place, how a workflow is launched, and the params of each workflow with where its output
lands. [Annotation](../../docs/annotation.md) explains what the workflows do.

## What is in place

| Input                                | Where                                                                                         |
| ------------------------------------ | --------------------------------------------------------------------------------------------- |
| the Annotator's image                | `australia-southeast1-docker.pkg.dev/cpg-common/images-dev/dragen`                            |
| the tools image                      | `australia-southeast1-docker.pkg.dev/cpg-common/images-dev/case-annotation-tools`             |
| the licence key                      | the Secret Manager secret `dragen-api-key` in the dataset's project                           |
| the annotation data                  | `gs://cpg-rd-benchmark-test/references/dragen-annotation/4.0.0/202610/GRCh38/`                |
| GRCh38                               | `gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta`                       |
| the archive of gnomAD's joint values | `gs://cpg-common-test/gnomad/echtvar_simple/gnomad_4.1_whole_genome.zip`                      |
| the gene models                      | `gs://cpg-rd-benchmark-test/references/annotation/gene_models/ensembl-115_refseq-RS_2025_08/` |
| the archives of AlphaGenome's scores | `gs://cpg-common-test/references/alphagenome/echtvar/`                                        |
| the joined files and STRchive's loci | `gs://cpg-rd-benchmark-test/references/annotation/`                                           |
| outputs                              | `gs://cpg-rd-benchmark-test`                                                                  |

- The Annotator's image is CPG's DRAGEN image, which holds the Annotator at `/opt/dragen/share/nirvana/Annotator` with
  `gcp-secret` and `jq` ([The Annotator's image](../../docs/annotator-image.md)). It is built by hand, because DRAGEN's
  installer is behind an Illumina login, and runs at the test access level pull it from `images-dev`.

- The tools image is built from [containers/tools/Dockerfile](../../containers/tools/Dockerfile) at a release tag, and
  the Dockerfile shows the command. The params file pins it by digest, so a new tag means a new build and a new digest
  there. The pinned image was built at the tag `v0.5.2`, which is the commit `3aa097702d5e1657486cea8cbea22a2671d1188f`,
  and the registry tags it with both. The image has to be of the release of the workflow that a run launches:
  `family_tables.nf` gives its own version to its first step and to the step that writes a family's folder, and each
  refuses to run with tools of another version, the first before anything is annotated. A run of a later release
  therefore needs the image built at that release's tag, and the pin moved to it.

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

- The gene models are the folder that `case-annotation gene-models` wrote from Ensembl's release 115 and RefSeq's
  release `GCF_000001405.40-RS_2025_08`, which are the releases of the catalog `202610`
  ([The gene models](../../docs/annotation.md#the-gene-models)). Its `gene-models.json` records the public files with
  their SHA-256.

- The three archives of AlphaGenome's scores are echtvar archives that were built from DeepMind's downloads
  ([AlphaGenome's archives](../../docs/annotation.md#alphagenomes-archives)). Each has a folder of its own, which also
  holds its SHA-256 and a `receipt.json` that records the download it was built from. The breakdown archive
  (`alphagenome_avi_breakdown`) and the splicing archive (`alphagenome_splicing`) are for non-commercial use only, and
  the folder of each holds a `TERMS.txt` with the notice that DeepMind's terms ask for
  ([The terms of AlphaGenome's files](../../docs/annotation.md#the-terms-of-alphagenomes-files)).

- The seven public files of [Joined files](../../docs/annotation.md#joined-files) are copies as their sources publish
  them, each of the size that the section gives. The catalogue of branch points was prepared from the two datasets of
  Zhang et al. with UCSC's chain file `hg19ToHg38.over.chain.gz`, and the `branch-points.json` beside it records the
  three files with their SHA-256.
  [The catalogue of branch points](../../docs/annotation.md#the-catalogue-of-branch-points) quotes the terms of that
  chain file.

- STRchive's loci are the file of the release `v2.28.0`, which STRchive published on 8 October 2026
  ([STRchive's ranges of the repeats](../../docs/annotation.md#strchives-ranges-of-the-repeats)). No other release is in
  the bucket.

[test.annotator.params.yaml](test.annotator.params.yaml) holds the images, the Annotator's path and version, the
licence's secret and the output bucket for runs at the test access level, with the images pinned by digest. A run copies
it and adds its own params. The reference files that only `family_tables.nf` takes are not in it, because a run of
`annotator.nf` stops at a param that the workflow does not declare. Their lines are in [Family tables](#family-tables).
No other access level has a params file yet.

## Launching a workflow

A Seqera launch takes one params file and one entry script, at one commit:

```sh
analysis-runner seqera --dataset rd-benchmark --access-level test \
    --repository https://github.com/populationgenomics/case-annotation-pipeline \
    --main-script <annotator.nf or family_tables.nf> --commit-id <commit> \
    --params-file <run>.params.yaml
```

The analysis-runner launches with Nextflow's v2 syntax parser, which the pipeline needs. Access levels above `test` need
the commit to be on `main`. The dataset has to allow this repository for the analysis-runner. The repository is public,
so Seqera needs no credential to clone it.

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
tools image, which the params file pins as built at the tag `v0.5.2`.

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
line of `test.annotator.params.yaml` is one of them, and its tools image is the one built at the tag `v0.5.2`. A run
adds CPG's reference, its archives, its gene models and its other reference files, and its own sheet, group and catalog:

```sh
cp deployments/cpg/test.annotator.params.yaml <run>.params.yaml
cat >> <run>.params.yaml <<'EOF'
reference: gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta
gnomad_joint: gs://cpg-common-test/gnomad/echtvar_simple/gnomad_4.1_whole_genome.zip
gene_models: gs://cpg-rd-benchmark-test/references/annotation/gene_models/ensembl-115_refseq-RS_2025_08
alphagenome_avi: gs://cpg-common-test/references/alphagenome/echtvar/alphagenome_avi/alphagenome_avi.zip
alphagenome_avi_breakdown: gs://cpg-common-test/references/alphagenome/echtvar/alphagenome_avi_breakdown/alphagenome_avi_breakdown.zip
alphagenome_splicing: gs://cpg-common-test/references/alphagenome/echtvar/alphagenome_splicing/alphagenome_splicing.zip
gnomad_mito: gs://cpg-rd-benchmark-test/references/annotation/gnomad_mito/v3.1/gnomad.genomes.v3.1.sites.chrM.vcf.bgz
mitotip: gs://cpg-rd-benchmark-test/references/annotation/mitotip/2020-04-27/mitotip_scores.txt
gencode_polya: gs://cpg-rd-benchmark-test/references/annotation/gencode_polya/v50/gencode.v50.polyAs.gtf.gz
aparent2: gs://cpg-rd-benchmark-test/references/annotation/aparent2/2022-11-12/aparent2_ism_scores_polyadb_v3.csv.gz
genomic_constraint: gs://cpg-rd-benchmark-test/references/annotation/genomic_constraint/v3.1.2/constraint_z_genome_1kb.qc.download.txt.gz
mcr: gs://cpg-rd-benchmark-test/references/annotation/mcr/v4.1.1/gnomad_v4.1.1_non_outlier_mcrs.tsv
mpc: gs://cpg-rd-benchmark-test/references/annotation/mpc/v4.1.1/gnomad_v4.1.1_mpc.tsv.bgz
branch_points: gs://cpg-rd-benchmark-test/references/annotation/branch_points/pnas.2211194119-ucsc-hg19ToHg38/branch-points.tsv.gz
strchive: gs://cpg-rd-benchmark-test/references/annotation/strchive/v2.28.0/STRchive-loci-v2.28.0.json
strchive_release: v2.28.0
strchive_release_date: '2026-10-08'
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

Of these files, a run needs `reference`, `gnomad_joint` and `gene_models`. It can leave out any of the three AlphaGenome
archives and any of the eight files from `gnomad_mito` to `branch_points`, and its tables then lack the columns of what
it left out. The three `strchive` lines are given together or left out together, and the date is a quoted string.

> [!IMPORTANT]
>
> A run that gives `alphagenome_avi_breakdown` or `alphagenome_splicing` writes tables under DeepMind's terms for
> non-commercial use: `variants` then has the `avi_shap_` columns or `alphagenome_splicing`. Tables with those columns
> must not be handed to a commercial organisation. A run that leaves both lines out writes tables without those columns,
> and `alphagenome_avi` alone is offered for commercial use as well.

> [!WARNING]
>
> A run with `sources` starts the Annotator, which is billed: once over the group's alleles, and a second time over its
> structural variants where a row names an `sv_vcf` or a `cnv_vcf`.

With `sv_vcf` or `cnv_vcf` in a sheet, `sources` should name the sources that hold structural variants (`gnomad`,
`oneKg` and `clinvar`).

A run that failed after the Annotator finished is launched again with `annotated` in place of `sources`, which reads the
annotation that the first run published and starts no Annotator. So is a run of the same group and catalog at a commit
with another version of the tables, because the folder does not name that version
([Reusing an annotation](../../docs/annotation.md#reusing-an-annotation)):

```yaml
annotated: gs://cpg-rd-benchmark-test/annotation-unions/<group>/illumina-<catalog>
```

That folder can be reused when it holds `out/<group>.json.gz` and `run.json`, and `out/<group>.sv.json.gz` as well when
a family of the group names an `sv_vcf` or a `cnv_vcf`.

A release before `0.3.0` published the Annotator's output to
`gs://cpg-rd-benchmark-test/annotation-unions/<group>/tables-<tables version>-illumina-<catalog>/annotator/`. Such a
folder is reused by giving it as `annotated`, or by copying its contents to the folder above, where a launcher that
fills in the contract's folder looks for it.

| Output                                      | Where                                                                                              |
| ------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| each family's folder                        | `gs://cpg-rd-benchmark-test/annotations/<family>/tables-<tables version>-illumina-<catalog>/`      |
| the union, its tables and its gnomAD values | `gs://cpg-rd-benchmark-test/annotation-unions/<group>/tables-<tables version>-illumina-<catalog>/` |
| the Annotator's output for the group        | `gs://cpg-rd-benchmark-test/annotation-unions/<group>/illumina-<catalog>/`                         |

Both folders of the group hold every family's alleles and records, so only a family's own folder is given to whoever
analyses that family ([What a run publishes](../../docs/annotation.md#what-a-run-publishes)).
