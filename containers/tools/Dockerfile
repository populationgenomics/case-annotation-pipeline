# The tools image for the pipeline's own steps: the case-annotation command (this repository's package, as locked),
# bcftools and echtvar. It builds from the repository root, which holds the lock:
#
#   docker build --platform linux/amd64 -f containers/tools/Dockerfile -t <registry>/case-annotation-tools:<commit> .
#
# amd64 only, because the pipeline's VMs are. Behind a proxy that intercepts TLS, add
# `--secret id=ca-bundle,src=<the proxy's CA bundle>`: the downloads trust it, and the image does not keep it.
#
# DuckDB 2's linux/amd64 wheel is compiled for x86-64-v3, so `import duckdb` stops with an illegal instruction on a
# processor that lacks AVX2, BMI2 or FMA. The pipeline's VMs have them. On Apple silicon this image builds and runs under
# emulation, and the emulator has to provide them too, or the check that ends this build fails. Under Colima, Rosetta
# provides them (`colima start --vz-rosetta`, for a VM of type vz): with DuckDB 2.0.0 on macOS 27, this build passes
# there. The QEMU 7.0.0 that Colima registers otherwise, up to at least Colima 0.10.3, has no AVX. Of QEMU's later
# user-mode emulators, 8.0.4 and 8.1.5 ran the wheel, while 9.2.2, 10.0.4 and 10.2.3 crashed as soon as pyarrow or
# DuckDB was imported.

FROM ghcr.io/astral-sh/uv:0.13.8@sha256:5cc1805e2426484abca4ad46ea783bab281fdb225fd82d8ce6d66a81254bac85 AS uv

# The package without its sources or its dev groups. Its metadata names the README, so the build reads that too.
FROM python:3.12.16-slim-trixie@sha256:b18b261e1b4974b4f321451ad5251ded85abb22c2b8b262da527dda1d031ab30 AS python-packages
COPY --from=uv /uv /usr/local/bin/uv
ENV UV_PYTHON_DOWNLOADS=never \
    UV_PYTHON=/usr/local/bin/python3.12 \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv
WORKDIR /src
COPY pyproject.toml uv.lock README.md ./
COPY src src
RUN --mount=type=secret,id=ca-bundle,required=false \
    if [ -f /run/secrets/ca-bundle ]; then export SSL_CERT_FILE=/run/secrets/ca-bundle; fi \
    && uv sync --frozen --no-default-groups --no-editable

FROM python:3.12.16-slim-trixie@sha256:b18b261e1b4974b4f321451ad5251ded85abb22c2b8b262da527dda1d031ab30
# bcftools as Debian packages it, with Debian's htslib, which reads gs:// and https:// inputs. Compiling a newer release
# instead is not reliable on Apple silicon, where QEMU's x86 emulation crashes gcc now and then. The exact version fails
# the build once Debian replaces it, so that an update is a deliberate edit here. procps: Nextflow reads a task's
# resource use with ps.
RUN apt-get update -q \
    && apt-get install -q -y --no-install-recommends bcftools=1.21-1 procps \
    && rm -rf /var/lib/apt/lists/*
# echtvar's static release binary, checked against the SHA-256 GitHub records for it.
ADD --checksum=sha256:2c425bad9a4e6aeeade5eeec229e7f887d7385642a0aada3a03909b6303f1e09 \
    https://github.com/brentp/echtvar/releases/download/v0.2.4/echtvar-v0.2.4-x86_64-unknown-linux-musl \
    /usr/local/bin/echtvar
RUN chmod 0755 /usr/local/bin/echtvar
COPY --from=python-packages /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

# An image whose tools do not start is never pushed.
RUN case-annotation --help >/dev/null && bcftools version && echtvar --version
