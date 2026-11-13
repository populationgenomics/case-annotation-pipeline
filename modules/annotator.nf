// Illumina's standalone DRAGEN Annotator (DRAGEN Annotation 4.0, Illumina Connected Annotations), in an image of the
// deployment's that also holds gcp-secret and jq. Each task fetches the licence key into the VM's /tmp (fetchKey),
// hands it to the Annotator as an API key file, and deletes it when the task ends. The Annotator takes other kinds of
// credential too, but the redaction and the checks below look for the key file's lines, which would miss a secret value
// inside a credentials JSON.

// One listing of the data Illumina hosts, `Annotator list` with the listing's arguments, written to <name>.txt.
process ILLUMINA_LIST {
    tag "${name}"
    label 'illumina'

    input:
    tuple val(name), val(args)
    val annotator
    val license_secret

    output:
    path "${name}.txt"

    script:
    """
    ${prelude(annotator, license_secret, "${name}.txt")}
    run_annotator list ${quoted(args)} "\${license[@]}"
    """

    stub:
    """
    echo 'stub listing: ${args.join(' ')}' > ${name}.txt
    """
}

// Downloads one catalog's data for one assembly: `setup` writes a run config per category, the categories' sources
// less the excluded ones go into one run config, runConfigs/bundle.json, and `download` fetches them. The data
// directory, without the license/ folder setup writes the credentials into, is the task's output, named after the
// assembly; setup's own run configs move to runConfigs.setup/, so that runConfigs/ holds only the bundle's.
// The stub gives its folder a fixed time, as an object store gives every folder: a second stub run's folder then
// looks to Nextflow like the one the first run published, which is what a run into a bucket meets.
process ILLUMINA_DOWNLOAD {
    tag "${catalog} ${assembly}"
    label 'illumina'
    // Bound by Illumina's servers, not the CPUs. The data is downloaded to the boot disk, then written to the work
    // directory, which stages each file on the same disk.
    machineType 'n2d-standard-8'
    cpus 8
    memory 30.GB
    disk 500.GB

    input:
    val catalog
    val assembly
    val categories
    val exclude_sources
    val annotator
    val annotator_version
    val license_secret

    output:
    path assembly, emit: bundle
    path "${assembly}.download.log", emit: log

    script:
    def configs = categories.collect { c -> "\"\${scratch}/data/runConfigs/${assembly}.${catalog}.${c}.json\"" }.join(' ')
    """
    ${prelude(annotator, license_secret, "${assembly}.download.log")}
    check_version ${annotator_version}
    run_annotator setup --non-interactive --data.directory "\${scratch}/data" --data.catalog '${catalog}' \\
        --data.assemblies '${assembly}' --data.autoDownload false --license.autoDetect false "\${license[@]}"

    # One run config from the categories' configs, with the excluded sources left out. A source and type that two
    # categories pin at different versions, or an excluded source no category has, fails. The config's own license and
    # data directory go, since they name files on this VM.
    jq -s --argjson exclude '${jsonList(exclude_sources)}' '
        (map(.data.versionSelector | to_entries[] | .key as \$s | .value | to_entries[]
            | {source: \$s, type: .key, version: .value})) as \$all
        | (\$exclude - (\$all | map(.source))) as \$unknown
        | if \$unknown != [] then error("no category has the excluded sources \\(\$unknown)") else . end
        | (\$all | group_by([.source, .type]) | map(
            if (map(.version) | unique | length) > 1
            then error("\\(.[0].source) \\(.[0].type): the categories pin different versions") else .[0] end
          ) | map(select(.source as \$s | any(\$exclude[]; . == \$s) | not))) as \$kept
        | .[0] | del(.license) | del(.data.directory)
        | .data.versionSelector = (reduce \$kept[] as \$e ({}; .[\$e.source][\$e.type] = \$e.version))
    ' ${configs} > "\${scratch}/bundle.json"
    mkdir "\${scratch}/data/runConfigs.setup"
    mv "\${scratch}"/data/runConfigs/*.json "\${scratch}/data/runConfigs.setup/"
    mv "\${scratch}/bundle.json" "\${scratch}/data/runConfigs/bundle.json"
    jq . "\${scratch}/data/runConfigs/bundle.json" | tee -a ${assembly}.download.log

    run_annotator download --config "\${scratch}/data/runConfigs/bundle.json" --data.directory "\${scratch}/data" \\
        --parallel.workers ${task.cpus} "\${license[@]}"
    run_annotator version-validate --config "\${scratch}/data/runConfigs/bundle.json" --data.directory "\${scratch}/data"

    # Nothing that leaves the VM may hold the credentials: no file by its name, and no small file (configs, logs)
    # holding a line of the key.
    find "\${scratch}/data" -mindepth 1 -maxdepth 1 ! -name license -print0 > "\${scratch}/entries"
    found=\$(find "\${scratch}/data" -path "\${scratch}/data/license" -prune -o \\( -iname '*credential*' -o -iname '*.lic' \\) -print)
    if [[ -n \${found} ]]; then
        echo "credentials files outside license/: \${found}" >&2
        exit 1
    fi
    refuse_key_in "\${scratch}/data" license
    (cd "\${scratch}/data" && find . -path ./license -prune -o -type f -printf '%s\\t%p\\n' | sort -k2) >> ${assembly}.download.log
    du -sh --exclude=license "\${scratch}/data" | tee -a ${assembly}.download.log

    mkdir ${assembly}
    xargs -0 -r cp -r -t ${assembly} < "\${scratch}/entries"
    """

    stub:
    """
    mkdir -p ${assembly}/runConfigs
    echo '{"data": {"catalog": "${catalog}", "versionSelector": {"stub": {"SmallVariant": "1"}}}}' > ${assembly}/runConfigs/bundle.json
    echo 'stub download' > ${assembly}.download.log
    touch -t 200001010000 ${assembly}
    """
}

// Annotates one VCF with a downloaded bundle, its run config cut down to the sources named, into Illumina's JSON (and
// the .jsi index beside it, where the Annotator writes one). An invalid licence for any source fails the task rather
// than leaving the source out, as does a small output file holding the key.
process ILLUMINA_ANNOTATE {
    tag "${name}"
    label 'illumina'
    // The bundle is copied to the local SSDs, mounted at /tmp, rather than read through the bucket's mount.
    machineType 'n2d-highmem-16'
    cpus 16
    memory 120.GB
    disk [request: 375.GB, type: 'local-ssd']

    input:
    tuple val(name), path(vcf, stageAs: 'input/*')
    path data, stageAs: 'bundle'
    val catalog
    val sources
    val workers
    val annotator
    val annotator_version
    val license_secret

    output:
    path 'out/*', emit: annotation
    path 'run.json', emit: config
    path "${name}.annotate.log", emit: log

    script:
    """
    ${prelude(annotator, license_secret, "${name}.annotate.log")}
    check_version ${annotator_version}
    # The bundle's run config with only the sources named; a source the bundle lacks, or another catalog, fails.
    jq --argjson sources '${jsonList(sources)}' --arg catalog '${catalog}' '
        if .data.catalog != \$catalog then error("the bundle is of catalog \\(.data.catalog), not \\(\$catalog)") else . end
        | (\$sources - (.data.versionSelector | keys)) as \$missing
        | if \$missing != [] then error("the bundle has no \\(\$missing)") else . end
        | .data.versionSelector |= with_entries(select(.key as \$k | any(\$sources[]; . == \$k)))
    ' bundle/runConfigs/bundle.json > run.json
    jq . run.json | tee -a ${name}.annotate.log

    started=\$(date +%s)
    cp -rL bundle "\${scratch}/data"
    echo "copied the bundle in \$((\$(date +%s) - started)) s" | tee -a ${name}.annotate.log
    started=\$(date +%s)
    # The Annotator writes into the output directory but does not create it.
    mkdir "\${scratch}/out"
    run_annotator annotate --config "\${PWD}/run.json" --input.file "\${PWD}/${vcf}" --output.directory "\${scratch}/out" \\
        --output.prefix ${name} --output.format json --data.directory "\${scratch}/data" --parallel.workers ${workers} \\
        --annotationOptions.enableDq true --license.throwInvalidLicenseError true "\${license[@]}"
    echo "annotated in \$((\$(date +%s) - started)) s" | tee -a ${name}.annotate.log
    ls -l "\${scratch}/out" | tee -a ${name}.annotate.log
    refuse_key_in "\${scratch}/out"
    mkdir out
    cp -r "\${scratch}"/out/. out/
    """

    stub:
    """
    cp bundle/runConfigs/bundle.json run.json
    mkdir out
    touch out/${name}.json.gz out/${name}.json.gz.jsi
    echo 'stub annotation' > ${name}.annotate.log
    """
}

// The records of a VCF in the regions of a BED file, for a small annotation run.
process SUBSET_VCF {
    tag "${name}"
    label 'tools'

    input:
    tuple val(name), path(vcf, stageAs: 'input/*'), path(index, stageAs: 'input/*')
    path regions

    output:
    tuple val(name), path("${name}.vcf.gz")

    script:
    """
    bcftools view --no-version -R ${regions} -Oz -o ${name}.vcf.gz ${vcf}
    echo "\$(bcftools view -H ${name}.vcf.gz | wc -l) records in the regions"
    """

    stub:
    """
    touch ${name}.vcf.gz
    """
}

// The command that fetches the licence key into /tmp/dragen-api-key, which only the task's user can read.
def fetchKey(api_key_secret) {
    return "(umask 077 && gcp-secret ${api_key_secret} > /tmp/dragen-api-key)"
}

// The shell every Annotator task starts with: the key fetched into /tmp, a scratch folder for the Annotator removed
// with the key when the task ends, and run_annotator, which runs the Annotator in the scratch folder with HOME there,
// so that whatever it caches (its .NET runtime, its premium licence) stays on the VM rather than landing in the work
// directory. The key's lines of 8 characters or more, which the redaction and the checks look for, are written to
// the scratch folder first; a key with none fails before the Annotator runs. The Annotator's output goes to the log and
// the task's output with any of those lines replaced by [redacted], in case an error message quotes it.
def prelude(annotator, license_secret, log) {
    return """
set -o pipefail
${fetchKey(license_secret)}
scratch=\$(mktemp -d /tmp/annotator.XXXXXX)
trap 'rm -rf /tmp/dragen-api-key "\${scratch}"' EXIT
# Where the Annotator caches its premium licence under HOME; it does not create the folder itself.
mkdir -m 700 "\${scratch}/.ilmnAnnotations"
annotator_bin=${annotator}
log=${log}
touch "\${log}"
license=(--license.apiKeyFile /tmp/dragen-api-key)
""" + preludeFunctions() + 'key_lines\n'
}

def preludeFunctions() {
    return '''
key_lines() {
    awk '{ sub(/\\r$/, "") } length($0) >= 8' /tmp/dragen-api-key > "${scratch}/key-lines"
    if [[ ! -s ${scratch}/key-lines ]]; then
        echo 'the key has no line of 8 or more characters to redact' >&2
        return 1
    fi
}
redact() {
    awk -v key_file="${scratch}/key-lines" '
        BEGIN {
            while ((getline line < key_file) > 0) keys[++n] = line
            if (!n) { print "no key lines to redact" > "/dev/stderr"; exit 1 }
        }
        {
            for (i = 1; i <= n; i++)
                while ((p = index($0, keys[i])) > 0)
                    $0 = substr($0, 1, p - 1) "[redacted]" substr($0, p + length(keys[i]))
            print
        }'
}
run_annotator() {
    local status=0
    (cd "${scratch}" && HOME="${scratch}" "${annotator_bin}" "$@") > "${scratch}/output" 2>&1 || status=$?
    redact < "${scratch}/output" | tee -a "${log}"
    return "${status}"
}
# Fails if a file under 10 MB in the folder holds a line of the key: configs and logs that could quote it, rather than
# the data. A second argument names a subfolder left out, the data directory's license/. A folder that cannot be listed
# and a file that cannot be read have not been checked, so each fails here with a message naming it, rather than through
# the shell's -e, which does not apply where the check is called in a condition. The size limit is in bytes (10 MiB):
# find rounds a size up to whole units before comparing, so -10M would leave out every file over 9 MiB.
refuse_key_in() {
    local file status found=() pruned=()
    if [[ -n ${2:-} ]]; then
        pruned=(-path "$1/$2" -prune -o)
    fi
    if ! find "$1" "${pruned[@]}" -type f -size -10485760c -print0 > "${scratch}/checked"; then
        echo "could not list every file in $1 to check it for the key" >&2
        return 1
    fi
    while IFS= read -r -d '' file; do
        # grep exits 1 for a file without a line of the key, and 2 for one it cannot read.
        status=0
        grep -qF -f "${scratch}/key-lines" -- "${file}" || status=$?
        case ${status} in
            0) found+=("${file}") ;;
            1) ;;
            *)
                echo "could not read ${file} to check it for the key" >&2
                return 1
                ;;
        esac
    done < "${scratch}/checked"
    if ((${#found[@]})); then
        echo "files holding the key: ${found[*]}" >&2
        return 1
    fi
}
check_version() {
    local banner
    banner=$(cd "${scratch}" && HOME="${scratch}" "${annotator_bin}" --version 2>&1 || true)
    if ! grep -qF "DRAGEN Annotation $1-" <<< "${banner}"; then
        echo "the image's Annotator is not version $1:" >&2
        echo "${banner}" >&2
        return 1
    fi
}
'''
}

// Each argument in single quotes; the workflow refuses values that hold one.
def quoted(args) {
    return args.collect { a -> "'${a}'" }.join(' ')
}

// A list of strings as a JSON array, for jq; the workflow refuses values that hold a quote.
def jsonList(values) {
    return '[' + values.collect { v -> "\"${v}\"" }.join(',') + ']'
}
