# The image of the step that runs Ensembl VEP (docs/annotation.md): the official VEP image with the UTRannotator plugin
# beside it. It builds from this folder, which holds nothing that the build reads:
#
#   docker build --platform linux/amd64 -t <registry>/case-annotation-vep:<commit> containers/vep
#
# The pipeline's VMs are amd64. The base is pinned by the digest of its index, which also holds an arm64 image, so the
# build runs without emulation on Apple silicon too.
#
# VEP's own collection has a plugin of the same name in other letters, UTRAnnotator, which this image does not run: it
# reports one kind of effect for an allele and transcript where several apply, and which one differs between runs. The
# plugin here is the original, of Zhang, Wakeling, Ware and Whiffin, at one commit of its repository. Its files are
# fetched from that commit, each checked against its SHA-256, with its MIT licence beside them (README.md).

FROM ensemblorg/ensembl-vep:release_116.3@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3

# Where the plugin lies, which the step gives VEP as --dir_plugins, and the commit that it is of, which the step
# writes into its record.
ENV UTRANNOTATOR_DIR=/opt/utrannotator \
    UTRANNOTATOR_COMMIT=3d9d466e356752aa69d61b62da493edddf374363

# The base image's user is vep. The pipeline's other images run as root, which can write into a task's work directory
# however the executor mounts it, so this one does too.
USER root

# The plugin's folders, made before the files: a folder that ADD makes for a file gets the file's mode, which no user
# but root can enter.
RUN mkdir -p "${UTRANNOTATOR_DIR}/src"

# The plugin's module, the six files that it loads from src/, the file of upstream reading frames that were seen
# translated, and the licence. A file from a URL is readable by root alone unless its mode is given.
ADD --chmod=0644 --checksum=sha256:fd51f3f336e330d655be045cd6a79a58372ca2ae984c129e00a3e8b38e3dd658 \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/UTRannotator.pm \
    ${UTRANNOTATOR_DIR}/UTRannotator.pm
ADD --chmod=0644 --checksum=sha256:d6ebc67da174787e19ad18726f83b56629860a195e7b4352cb38a3b1d7391df5 \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/src/five_prime_UTR_utils.pl \
    ${UTRANNOTATOR_DIR}/src/five_prime_UTR_utils.pl
ADD --chmod=0644 --checksum=sha256:3e1576a35be1b029772ef940d471df53e700b320da7e5f4391f0625d56debd9c \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/src/uAUG_gained.pl \
    ${UTRANNOTATOR_DIR}/src/uAUG_gained.pl
ADD --chmod=0644 --checksum=sha256:84ed7f088f939ff750782bc2c7d3d417d2a1278d078349abcfe673f5fb3366c5 \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/src/uAUG_lost.pl \
    ${UTRANNOTATOR_DIR}/src/uAUG_lost.pl
ADD --chmod=0644 --checksum=sha256:68ca27d46dac0715f8c4b3098277d760eabd8718a6fda2fc90c332b0b70847ee \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/src/uSTOP_gained.pl \
    ${UTRANNOTATOR_DIR}/src/uSTOP_gained.pl
ADD --chmod=0644 --checksum=sha256:f9e57fd15aa45b5a38917424efdbc835eb97a0da383e581c93d94d16421ccf18 \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/src/uSTOP_lost.pl \
    ${UTRANNOTATOR_DIR}/src/uSTOP_lost.pl
ADD --chmod=0644 --checksum=sha256:766a4f2970ca3eb6a68bc8cd97e5ad2cb1beafc81e65d92dd307c809ea2b140c \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/src/uframeshift.pl \
    ${UTRANNOTATOR_DIR}/src/uframeshift.pl
ADD --chmod=0644 --checksum=sha256:e220023a747465dcdf6725ce809511bccc20b7a9990288da91032a1a7f51ebce \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/uORF_5UTR_GRCh38_PUBLIC.txt \
    ${UTRANNOTATOR_DIR}/uORF_5UTR_GRCh38_PUBLIC.txt
ADD --chmod=0644 --checksum=sha256:dfcec4f825cb1811176fb8f50eb67e84dfc58c6a2fdd8849b214ca74036d5721 \
    https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/${UTRANNOTATOR_COMMIT}/LICENSE \
    ${UTRANNOTATOR_DIR}/LICENSE

# An image whose VEP does not start, or cannot compile the plugin and the NMD plugin of its own collection, is never
# pushed. Perl finds both through the paths that the base image sets.
RUN vep --help >/dev/null \
    && perl -I"${UTRANNOTATOR_DIR}" -e 'use UTRannotator; use NMD; print "the plugins compile\n"'
