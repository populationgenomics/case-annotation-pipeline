// Runs Illumina's standalone DRAGEN Annotator (docs/annotation.md) in one of three modes. list asks Illumina
// which catalogs it hosts and which versions of each data source, and publishes the listings under
// annotation/list/<date>/. download fetches one catalog's data for one assembly, the sources of the categories named
// less those excluded, and publishes the data directory under references/dragen-annotation/<version>/<catalog>/
// <assembly>/. annotate annotates one VCF with such a directory and the sources named, and publishes Illumina's JSON
// under annotations/<name>/illumina-<catalog>/.

include { ILLUMINA_LIST ; ILLUMINA_DOWNLOAD ; ILLUMINA_ANNOTATE ; SUBSET_VCF } from './modules/annotator'

params {
    // list, download or annotate.
    mode: String
    // The image holding the Annotator, gcp-secret and jq, the Annotator's path in it, and its version, which the
    // download and annotate tasks check and which names the downloaded data's folder.
    annotator_container: String
    annotator: String
    annotator_version: String
    // Secret Manager secret version holding the licence key, which the Annotator takes as an API key file:
    // projects/<project>/secrets/<secret>/versions/<version>.
    license_secret: String
    // The assembly listed, downloaded or annotated against.
    assembly: String = 'GRCh38'
    // The data sources: for list, a filter on those listed, all without it; for annotate, the sources annotated with,
    // each of which the downloaded data must hold.
    sources: List<String> = []
    // list: catalogs, by id or name, whose data sources are also listed one catalog at a time.
    catalogs: List<String> = []
    // download and annotate: the catalog's id, such as 202610, rather than an alias such as LATEST, which moves.
    catalog: String?
    // download: the categories whose run configs setup writes and whose sources are downloaded, and sources left out.
    categories: List<String> = ['basic', 'premium']
    exclude_sources: List<String> = []
    // annotate: the data directory a download published; the VCF, bgzipped, with its .tbi beside it if regions is
    // given; the name of the outputs and of their folder; a BED file the VCF is cut down to first, for a small run; a
    // suffix to the output folder, to keep such a run apart; and the Annotator's worker processes.
    annotation_data: Path?
    input: Path?
    name: String?
    regions: Path?
    label: String?
    workers: Integer = 8
    // annotate with regions: the tools image (containers/tools), which cuts the VCF.
    tools_container: String?
    // Read by nextflow.config.
    output_dir: String?
    // Some launchers set it on every launch; undeclared, it would fail the run.
    arGuid: String?
}

workflow {

    main:
    // The values go into the Annotator's command line in single quotes and into jq as JSON strings.
    (params.sources + params.catalogs + params.categories + params.exclude_sources)
        .findAll { v -> v != null }
        .each { v ->
            if (!(v ==~ /[A-Za-z0-9][A-Za-z0-9 ._()-]*/)) {
                error("${v}: a name may hold letters, digits, spaces and . _ ( ) - only")
            }
        }
    // The Annotator's version, the assembly, the catalog id, the output's name and its label become file and folder
    // names.
    [params.annotator_version, params.assembly, params.catalog, params.name, params.label]
        .findAll { v -> v != null }
        .each { v ->
            if (!(v ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                error("${v}: a version, assembly, catalog id, name or label may hold letters, digits and . _ - only")
            }
        }

    if (params.mode == 'list') {
        def filters = params.sources ? ['--data.sources', params.sources.join(',')] : []
        def listings = [
            ['catalogs', ['--catalogs']],
            [params.assembly, ['--data.assembly', params.assembly] + filters],
        ] + params.catalogs.collect { catalog ->
            ["${params.assembly}.catalog-${fileName(catalog)}", ['--data.assembly', params.assembly, '--data.catalog', catalog] + filters]
        }
        lists = ILLUMINA_LIST(
            channel.fromList(listings),
            params.annotator,
            params.license_secret,
        )
        bundle = channel.empty()
        download_log = channel.empty()
        annotation = channel.empty()
    }
    else if (params.mode == 'download') {
        requireCatalog()
        if (!params.categories) {
            error('download: name at least one category')
        }
        downloaded = ILLUMINA_DOWNLOAD(
            params.catalog,
            params.assembly,
            params.categories,
            params.exclude_sources,
            params.annotator,
            params.annotator_version,
            params.license_secret,
        )
        lists = channel.empty()
        bundle = downloaded.bundle
        download_log = downloaded.log
        annotation = channel.empty()
    }
    else if (params.mode == 'annotate') {
        requireCatalog()
        require(params.annotation_data, 'annotation_data')
        require(params.input, 'input')
        require(params.name, 'name')
        if (!params.sources) {
            error('annotate: name the sources to annotate with')
        }
        if (params.regions) {
            require(params.tools_container, 'tools_container')
            def index = params.input.resolveSibling("${params.input.name}.tbi")
            if (!index.exists()) {
                error("${params.input}: no .tbi beside it, which cutting it to regions needs")
            }
            vcf = SUBSET_VCF(channel.value(tuple(params.name, params.input, index)), params.regions)
        }
        else {
            vcf = channel.value(tuple(params.name, params.input))
        }
        annotated = ILLUMINA_ANNOTATE(
            vcf,
            params.annotation_data,
            params.catalog,
            params.sources,
            params.workers,
            params.annotator,
            params.annotator_version,
            params.license_secret,
        )
        lists = channel.empty()
        bundle = channel.empty()
        download_log = channel.empty()
        annotation = annotated.annotation.flatten().mix(annotated.config, annotated.log)
    }
    else {
        error("mode ${params.mode}: give list, download or annotate")
    }

    publish:
    lists = lists
    bundle = bundle
    download_log = download_log
    annotation = annotation
}

// Illumina's listing changes as it adds data, so each day's listing has a folder of its own. The downloaded data's
// folder is the data directory itself, named after the assembly.
// Nextflow's default replaces a published file whose size or time differs, and compares a folder by its own, which an
// object store gives every folder alike: a second download of a catalog would keep the first one's data directory
// beside its own log. The data directory is therefore always replaced, which removes an earlier one whole and leaves
// none of its files among the new ones.
output {
    lists {
        path "annotation/list/${workflow.start.toLocalDate()}"
    }
    bundle {
        path "references/dragen-annotation/${params.annotator_version}/${params.catalog}"
        overwrite true
    }
    download_log {
        path "references/dragen-annotation/${params.annotator_version}/${params.catalog}"
    }
    annotation {
        path "annotations/${params.name}/illumina-${params.catalog}${params.label ? '-' + params.label : ''}"
    }
}

def require(value, name) {
    if (value == null) {
        error("${params.mode}: give ${name}")
    }
}

// The catalog's id, which names the downloaded data's folder; LATEST is an alias that moves to each new catalog.
def requireCatalog() {
    require(params.catalog, 'catalog')
    if (params.catalog.equalsIgnoreCase('LATEST')) {
        error("${params.mode}: give the catalog's id rather than LATEST, which moves")
    }
}

// A catalog's id or name as part of a file name.
def fileName(value) {
    return value.replaceAll(/[^A-Za-z0-9._-]/, '_')
}
