# A stand-in for the Annotator's image without the Annotator: a base with jq, for tests/workflows/annotator.nf. CI
# builds it as annotator-test before the test runs:
#
#   docker build -t annotator-test tests/workflows/annotator
FROM almalinux:9.8-20261002@sha256:9819dc675b67b595c2b59e42be7763fca1a8bb217fa5944e04daa22e9a64db16

RUN dnf install -q -y jq && dnf clean all
