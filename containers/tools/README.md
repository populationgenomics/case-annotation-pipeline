# The tools image

The pipeline's own steps run in this image. It holds the `case-annotation` command, installed with this repository's
package from `uv.lock`, the three programs that the steps call beside it, bcftools, echtvar and
[upd](https://github.com/bjhall/upd), and one data file that a step reads. The [Dockerfile](Dockerfile) pins the
versions of bcftools and echtvar and the commit of the data file.

upd is installed from `uv.lock` too, into the environment of the package. PyPI's package named `upd` is another project,
so `pyproject.toml` gives upd's Git repository as its source, with the commit of its tag 0.1.1, and the lock holds that
commit and the versions of upd's two dependencies. upd is in the dependency group `image` and is no dependency of the
package, whose metadata would otherwise name a requirement that an index resolves to the other project. A newer upd is a
new commit in `pyproject.toml` and a new lock (`uv lock`), under the same age gate as every other package.

The data file belongs to UTRannotator, the plugin of Ensembl VEP whose rules `case-annotation vep-plugins` applies
([Upstream reading frames and nonsense-mediated decay](../../docs/annotation.md#upstream-reading-frames-and-nonsense-mediated-decay)).
The step reads it for the `evidence` of an upstream reading frame.

| What                          | From                                                                                                                                                                                                                                                                                       | Licence |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------- |
| `uORF_5UTR_GRCh38_PUBLIC.txt` | [ImperialCardioGenetics/UTRannotator](https://github.com/ImperialCardioGenetics/UTRannotator) at the commit `3d9d466e356752aa69d61b62da493edddf374363`: the plugin's list of upstream reading frames that ribosome profiling has seen translated, which its authors curated from sorfs.org | MIT     |

UTRannotator is the work of Xiaolei Zhang, Matthew Wakeling, James Ware and Nicola Whiffin ("Annotating high-impact 5′
untranslated region variants with the UTRannotator", Bioinformatics,
[doi:10.1093/bioinformatics/btaa783](https://doi.org/10.1093/bioinformatics/btaa783)). The build fetches the file and
the plugin's licence from the commit above, checks each against its SHA-256, and puts them under `/opt/utrannotator`
unchanged, with the licence as `/opt/utrannotator/LICENSE`. The image gives the file's path as `UORF_EVIDENCE`. This
repository holds no copy of the file.

Build the image from the repository's root, for linux/amd64, and tag it with the commit that it was built from:

```sh
commit=$(git rev-parse HEAD)
docker build --platform linux/amd64 -f containers/tools/Dockerfile \
  --label org.opencontainers.image.revision="${commit}" \
  -t <registry>/case-annotation-tools:"${commit}" .
```

Build from a clean checkout of that commit, because the image takes `src/`, `pyproject.toml`, `uv.lock` and `README.md`
from the working tree. Behind a proxy that intercepts TLS, add `--secret id=ca-bundle,src=<the proxy's CA bundle>`: the
downloads trust it, and the image does not keep it.

The build ends by starting each tool and looking for the data file, so an image whose tools do not start, or that lacks
the file, is never tagged. On Apple silicon the build runs under emulation, and the Dockerfile's opening comment says
which emulators can run DuckDB.

Build the image again whenever `src/`, `pyproject.toml` or `uv.lock` changes, and pin the new image by its digest where
a run names it.
