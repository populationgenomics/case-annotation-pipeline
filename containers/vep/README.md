# The VEP image

The step that runs Ensembl VEP on a union's alleles runs in this image
([Upstream reading frames and nonsense-mediated decay](../../docs/annotation.md#upstream-reading-frames-and-nonsense-mediated-decay)).
It is the official image of Ensembl VEP 116.3, pinned by its digest, with one plugin added. The [Dockerfile](Dockerfile)
pins both.

| What                             | From                                                                                                                                                   | Licence    |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------- |
| Ensembl VEP and the plugin `NMD` | `ensemblorg/ensembl-vep:release_116.3`, which holds the plugins of [Ensembl/VEP_plugins](https://github.com/Ensembl/VEP_plugins) at release 116        | Apache 2.0 |
| the plugin `UTRannotator`        | [ImperialCardioGenetics/UTRannotator](https://github.com/ImperialCardioGenetics/UTRannotator) at the commit `3d9d466e356752aa69d61b62da493edddf374363` | MIT        |
| `uORF_5UTR_GRCh38_PUBLIC.txt`    | the same commit: the plugin's list of upstream reading frames that ribosome profiling has seen translated, which its authors curated from sorfs.org    | MIT        |

UTRannotator is the work of Xiaolei Zhang, Matthew Wakeling, James Ware and Nicola Whiffin ("Annotating high-impact 5′
untranslated region variants with the UTRannotator", Bioinformatics,
[doi:10.1093/bioinformatics/btaa783](https://doi.org/10.1093/bioinformatics/btaa783)). The build fetches its seven
source files, its data file and its licence from the commit above, checks each against its SHA-256, and puts them under
`/opt/utrannotator` unchanged, with the licence as `/opt/utrannotator/LICENSE`. This repository holds no copy of them.

The collection that VEP ships has a plugin of nearly the same name, `UTRAnnotator`, which Ensembl maintains. The image
holds it, and the pipeline does not run it: for an allele and transcript with several kinds of effect it reports one,
and which one can differ between two runs of the same input. The original plugin reports every kind.

Build the image from this folder, for linux/amd64, and tag it with the commit that it was built from:

```sh
commit=$(git rev-parse HEAD)
docker build --platform linux/amd64 \
  --label org.opencontainers.image.revision="${commit}" \
  -t <registry>/case-annotation-vep:"${commit}" containers/vep
```

The build reads no file of the repository, so the image changes only when the Dockerfile does. It ends by starting VEP
and compiling both plugins, so an image that cannot run them is never tagged. Pin the image by its digest where a run
names it, as `vep_container`.

The image runs as root, as the pipeline's other images do. The official image's user is `vep`, which may not be able to
write into a task's work directory.
