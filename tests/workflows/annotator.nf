// Runs the Annotator's processes on an invented stand-in for the Annotator (annotator/Annotator) and an invented key
// (annotator/gcp-secret), in a stand-in for the Annotator's image (annotator/Dockerfile), and checks what they write:
// the download's run config and data directory, the annotation's run config, and that the key appears in no output.
// Refused inputs must fail, and the key check must name a file it cannot read and find the key in a file of nearly
// 10 MB. Run with Docker from the repository's root, with this folder's config alone (-C), because a run from the root
// would otherwise read the pipeline's nextflow.config over it:
//   docker build -t annotator-test tests/workflows/annotator
//   nextflow -C tests/workflows/nextflow.config run tests/workflows/annotator.nf

include { ILLUMINA_LIST ; ILLUMINA_DOWNLOAD ; ILLUMINA_ANNOTATE ; prelude } from '../../modules/annotator'
include { ILLUMINA_LIST as LIST_SHORT_KEY_REFUSED } from '../../modules/annotator'
include { ILLUMINA_DOWNLOAD as DOWNLOAD_UNKNOWN_EXCLUDED_REFUSED } from '../../modules/annotator'
include { ILLUMINA_DOWNLOAD as DOWNLOAD_VERSION_CONFLICT_REFUSED } from '../../modules/annotator'
include { ILLUMINA_DOWNLOAD as DOWNLOAD_KEY_IN_DATA_REFUSED } from '../../modules/annotator'
include { ILLUMINA_DOWNLOAD as DOWNLOAD_CREDENTIALS_FILE_REFUSED } from '../../modules/annotator'
include { ILLUMINA_DOWNLOAD as DOWNLOAD_OTHER_ANNOTATOR_REFUSED } from '../../modules/annotator'
include { ILLUMINA_ANNOTATE as ANNOTATE_MISSING_SOURCE_REFUSED } from '../../modules/annotator'
include { ILLUMINA_ANNOTATE as ANNOTATE_OTHER_CATALOG_REFUSED } from '../../modules/annotator'
include { ILLUMINA_ANNOTATE as ANNOTATE_KEY_IN_OUTPUT_REFUSED } from '../../modules/annotator'

workflow {
    // Where nextflow.config mounts annotator/ in the container.
    def annotator = '/invented/Annotator'
    def secret = 'projects/invented/secrets/invented/versions/1'
    def categories = ['basic', 'premium']
    vcf = channel.of('##fileformat=VCFv4.2').collectFile(name: 'invented.vcf', newLine: true).map { f -> tuple('invented', f) }

    listed = ILLUMINA_LIST(channel.of(['catalogs', ['--catalogs']]), annotator, secret)
    downloaded = ILLUMINA_DOWNLOAD('202610', 'GRCh38', categories, ['cosmic'], annotator, '4.0.0', secret)
    annotated = ILLUMINA_ANNOTATE(vcf, downloaded.bundle, '202610', ['clinvar', 'spliceAI'], 2, annotator, '4.0.0', secret)
    CHECK(listed, downloaded.bundle, downloaded.log, annotated.annotation, annotated.config, annotated.log)
    KEY_CHECK_UNREADABLE(annotator, secret)
    KEY_CHECK_LARGE(annotator, secret)

    // Refused: a key too short to redact; an excluded source no category has; two categories pinning one source at
    // different versions; setup writing the key outside license/, or a credentials file; an Annotator of another
    // version; a source the bundle lacks; a bundle of another catalog; the Annotator writing the key beside its output.
    short_key = LIST_SHORT_KEY_REFUSED(channel.of(['catalogs', ['--catalogs']]), annotator, 'projects/invented/secrets/short/versions/1')
    credentials_file = DOWNLOAD_CREDENTIALS_FILE_REFUSED('stray', 'GRCh38', categories, [], annotator, '4.0.0', secret)
    unknown_excluded = DOWNLOAD_UNKNOWN_EXCLUDED_REFUSED('202610', 'GRCh38', categories, ['cosmik'], annotator, '4.0.0', secret)
    version_conflict = DOWNLOAD_VERSION_CONFLICT_REFUSED('conflict', 'GRCh38', categories, [], annotator, '4.0.0', secret)
    key_in_data = DOWNLOAD_KEY_IN_DATA_REFUSED('leak', 'GRCh38', categories, [], annotator, '4.0.0', secret)
    other_annotator = DOWNLOAD_OTHER_ANNOTATOR_REFUSED('202610', 'GRCh38', categories, [], annotator, '4.0.1', secret)
    missing_source = ANNOTATE_MISSING_SOURCE_REFUSED(vcf, downloaded.bundle, '202610', ['clinvar', 'cosmic'], 2, annotator, '4.0.0', secret)
    other_catalog = ANNOTATE_OTHER_CATALOG_REFUSED(vcf, downloaded.bundle, '202608', ['clinvar'], 2, annotator, '4.0.0', secret)
    key_in_output = ANNOTATE_KEY_IN_OUTPUT_REFUSED(vcf.map { _name, f -> tuple('leak', f) }, downloaded.bundle, '202610', ['clinvar'], 2, annotator, '4.0.0', secret)
    unknown_excluded.bundle
        .map { _b -> 'unknown excluded source' }
        .mix(
            short_key.map { _l -> 'short key' },
            credentials_file.bundle.map { _b -> 'credentials file' },
            version_conflict.bundle.map { _b -> 'version conflict' },
            key_in_data.bundle.map { _b -> 'key in the data' },
            other_annotator.bundle.map { _b -> 'other Annotator version' },
            missing_source.config.map { _c -> 'missing source' },
            other_catalog.config.map { _c -> 'other catalog' },
            key_in_output.config.map { _c -> 'key in the output' },
        )
        .toList()
        .map { accepted ->
            if (accepted) {
                error("should have been refused, but were accepted: ${accepted}")
            }
        }
}

// The bundle's run config holds the basic and premium sources less COSMIC, without the config's own license and data
// directory; setup's five configs moved aside; no license/ folder. The annotation's run config holds the sources named.
// The key, which the stand-in prints in its listing and annotation, is redacted, and appears in no output (only in
// this task's own script).
process CHECK {
    label 'illumina'

    input:
    path listing
    path bundle, stageAs: 'bundle'
    path download_log
    path annotation, stageAs: 'out/*'
    path run_config
    path annotate_log

    script:
    """
    set -o pipefail
    test "\$(jq -c '.data.versionSelector | keys' bundle/runConfigs/bundle.json)" = '["clinvar","gnomad","omim","spliceAI"]'
    test "\$(jq -c '[has("license"), (.data | has("directory"))]' bundle/runConfigs/bundle.json)" = '[false,false]'
    test "\$(ls bundle/runConfigs)" = bundle.json
    test "\$(ls bundle/runConfigs.setup | wc -l)" = 5
    test ! -e bundle/license
    test "\$(jq -c '.data.versionSelector | keys' ${run_config})" = '["clinvar","spliceAI"]'
    test "\$(zcat out/invented.json.gz | jq -c '.header.dataSources | keys')" = '["clinvar","spliceAI"]'
    grep -q 'listing with key \\[redacted\\]' ${listing}
    grep -q 'annotating with key \\[redacted\\]' ${annotate_log}
    if grep -rlF --exclude='.command.*' invented-key-0123456789abcdef -R .; then
        echo 'the key is in an output' >&2
        exit 1
    fi
    """
}

// The key check on a folder holding a file it cannot read: one at a path longer than the 4,096 bytes a program can
// open a file by, which stops root as well. The check is called in a condition, where the shell's -e does not apply,
// so it has to fail on its own, and its message has to name the file.
process KEY_CHECK_UNREADABLE {
    label 'illumina'

    input:
    val annotator
    val license_secret

    script:
    """
    ${prelude(annotator, license_secret, 'check.log')}
    mkdir "\${scratch}/folder"
    echo 'invented, without the key' > "\${scratch}/folder/readable.log"
    refuse_key_in "\${scratch}/folder"
    (
        cd "\${scratch}/folder"
        level=\$(printf 'd%.0s' {1..200})
        for _ in {1..25}; do
            mkdir "\${level}"
            cd "\${level}"
        done
        echo 'invented, without the key' > unreadable.log
    )
    if refuse_key_in "\${scratch}/folder" 2> refusal.txt; then
        echo 'the key check passed a folder holding a file it cannot read' >&2
        exit 1
    fi
    grep -q 'could not read .*/unreadable.log to check it for the key' refusal.txt
    """
}

// The key check on a file of just over 10,000,000 bytes that ends with the key: under the check's limit, and over the
// 9 MiB up to which a limit given to find in whole MiB would reach.
process KEY_CHECK_LARGE {
    label 'illumina'

    input:
    val annotator
    val license_secret

    script:
    """
    ${prelude(annotator, license_secret, 'check.log')}
    mkdir "\${scratch}/folder"
    { head -c 10000000 /dev/zero | tr '\\0' '\\n'; cat /tmp/dragen-api-key; } > "\${scratch}/folder/large.log"
    if refuse_key_in "\${scratch}/folder" 2> refusal.txt; then
        echo 'the key check passed a file of nearly 10 MB that holds the key' >&2
        exit 1
    fi
    grep -q 'files holding the key: .*/large.log' refusal.txt
    """
}
