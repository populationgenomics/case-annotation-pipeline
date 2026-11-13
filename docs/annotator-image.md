# The Annotator's image

The pipeline runs Illumina's Annotator in an image that this repository does not build. The Annotator is Illumina's
licensed software, so the organisation that runs the pipeline builds an image that holds it, keeps the image in a
private registry and names it in a run's params. This page says what the image has to hold, for the person who builds
it. [Annotation](annotation.md#annotating-with-illuminas-annotator) explains what the Annotator's tasks do.

## The params that name the image

Both workflows take the image and what is in it as params.

| Param                 | Holds                                                                                                             |
| --------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `annotator_container` | the image                                                                                                         |
| `annotator`           | the Annotator's path in the image, such as `/opt/dragen/share/nirvana/Annotator`                                  |
| `annotator_version`   | the Annotator's version, such as `4.0.0`                                                                          |
| `license_secret`      | the secret version that holds Illumina's licence key, as `projects/<project>/secrets/<secret>/versions/<version>` |

`annotator.nf` needs the four in every mode. `family_tables.nf` needs them in a run that annotates, and not in one that
reuses an earlier run's annotation (`annotated`), which starts no Annotator.

The image runs the three processes that start the Annotator: `ILLUMINA_LIST`, `ILLUMINA_DOWNLOAD` and
`ILLUMINA_ANNOTATE`. Every other step runs in the image `case-annotation-tools`, which
[containers/tools/Dockerfile](../containers/tools/Dockerfile) builds.

## What the image holds

| In the image                              | What the tasks need of it                                                                                                                                                                                           |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| the Annotator                             | an executable at the path that `annotator` names. A download and an annotation first run it with `--version`, and stop unless what it prints holds `DRAGEN Annotation <annotator_version>-`                         |
| `gcp-secret`                              | a script on the `PATH` that takes the value of `license_secret` as its one argument and prints the licence key ([The helper that fetches the key](#the-helper-that-fetches-the-key))                                |
| `jq`                                      | on the `PATH`. A download joins the categories' run configs into one with it, and an annotation cuts the bundle's run config down to the sources named                                                              |
| bash                                      | the tasks' scripts use arrays, `[[ ]]` and `pipefail`                                                                                                                                                               |
| GNU coreutils, findutils and grep, an awk | `mktemp`, `tee`, `date`, `sort`, `cp -t` and `du --exclude`; `find -printf` and `xargs -0 -r`; `grep -F -f`; and `awk`, with which a task redacts the key from the Annotator's output before it is written or shown |

The Annotator caches its runtime and its licence under `HOME`. A task points `HOME` at a scratch folder under `/tmp` and
deletes the folder when it ends, so the image needs no home directory prepared and nothing in `/tmp` kept.

An image that meets the table is built in three steps.

1. Start from an operating system that DRAGEN supports, such as one compatible with RHEL 9. Its base image has bash and
   the GNU tools of the last two rows.
1. Run Illumina's installer of DRAGEN Software Mode, which holds the Annotator. DRAGEN 4.6 installs DRAGEN Annotation
   4.0 at `share/nirvana/Annotator` under its own folder.
1. Install `jq` and copy the helper onto the `PATH`.

> [!IMPORTANT]
>
> The installer and the image hold Illumina's licensed software. Keep the installer out of every repository, and push
> the image to a private registry only.

The check that the tasks make can be made by hand before a run. This command has to print a line that holds
`DRAGEN Annotation 4.0.0-` for `annotator_version: 4.0.0`:

```sh
docker run --rm <image> /opt/dragen/share/nirvana/Annotator --version
```

## The helper that fetches the key

A task calls the helper once, when it starts:

```sh
gcp-secret <license_secret> > /tmp/dragen-api-key
```

The helper prints the payload of that secret version and nothing else, and it exits with an error where it cannot read
the secret, which stops the task. For Google Secret Manager it is a few lines, which read the secret with the VM's
service account through the metadata server and need `curl`, `jq` and `base64` in the image:

```sh
#!/usr/bin/env bash
# Prints the payload of a Secret Manager secret version, named as projects/<project>/secrets/<secret>/versions/<version>,
# using the VM's service account through the metadata server.
set -euo pipefail

version=${1:?usage: gcp-secret projects/<project>/secrets/<secret>/versions/<version>}

token=$(curl -sSf -H 'Metadata-Flavor: Google' \
    http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token | jq -er .access_token)
curl -sSf -H "Authorization: Bearer ${token}" "https://secretmanager.googleapis.com/v1/${version}:access" \
    | jq -er .payload.data | base64 -d
```

The service account that the tasks run as needs permission to read that secret version, and the VMs need outbound access
to the secret manager. A listing and a download also fetch from Illumina's servers.

The workflows call the helper by this one name and pass `license_secret` on as it is given. With another secret manager,
a script of the same name that takes the value of `license_secret` and prints the key serves in its place.

## The licence key

The key is fetched at run time and stays on the VM of the task that fetched it.

- It is in no image and no params file. A params file holds the name of the secret, not its value.
- A task writes it to `/tmp/dragen-api-key`, readable by its owner only, hands it to the Annotator as an API key file
  (`--license.apiKeyFile`), and deletes it when the task ends.
- It is never written to a published file. A task replaces each line of the key in the Annotator's output with
  `[redacted]`, and fails if a file under 10 MB among those it is about to publish holds a line of the key
  ([The licence](annotation.md#the-licence)).
- The key needs a line of 8 characters or more, since those lines are what the redaction and the check look for. A task
  given a key without one stops before it runs the Annotator.
