# The tools image

The pipeline's own steps run in this image. It holds the `case-annotation` command, installed with this repository's
package from `uv.lock`, and the two programs that the steps call beside it, bcftools and echtvar. The
[Dockerfile](Dockerfile) pins the version of each.

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

The build ends by starting each tool, so an image whose tools do not start is never tagged. On Apple silicon the build
runs under emulation, and the Dockerfile's opening comment says which emulators can run DuckDB.

Build the image again whenever `src/`, `pyproject.toml` or `uv.lock` changes, and pin the new image by its digest where
a run names it.
