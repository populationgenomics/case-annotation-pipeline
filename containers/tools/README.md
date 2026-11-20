# The tools image

The pipeline's own steps run in this image. It holds the `case-annotation` command, installed with this repository's
package from `uv.lock`, and the three programs that the steps call beside it: bcftools, echtvar and
[upd](https://github.com/bjhall/upd). The [Dockerfile](Dockerfile) pins the versions of bcftools and echtvar.

upd is installed from `uv.lock` too, into the environment of the package. PyPI's package named `upd` is another project,
so `pyproject.toml` gives upd's Git repository as its source, with the commit of its tag 0.1.1, and the lock holds that
commit and the versions of upd's two dependencies. upd is in the dependency group `image` and is no dependency of the
package, whose metadata would otherwise name a requirement that an index resolves to the other project. A newer upd is a
new commit in `pyproject.toml` and a new lock (`uv lock`), under the same age gate as every other package.

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
