// Writes the family tables (docs/annotation.md) of the families that share one source selection: each family's VCF is
// normalised, the union of their sites is annotated once with Illumina's Annotator and converted into tables, the
// union's gnomAD joint values are read from an echtvar archive, and each family's tables are cut from both. The union
// and what annotates it reveal every family's alleles, so they are published under annotation-unions/<group>/ and never
// beside a family's tables, which go to annotations/<family>/. The Annotator's output has a folder of its own there,
// which a later run can reuse in place of annotating again. A family whose row also names VCFs of its larger variants
// (its structural variants, its copy-number segments or its members' repeats) gets tables of those too. The structural
// variants and copy-number segments of the group's families go the same way as their alleles: the union of their sites
// is annotated once, in a second Annotator run with the same catalog and sources, and each family's tables read its
// records' annotation from it. Nothing annotates the repeats.

include { FAMILY_CHECK ; NORMALISE_FAMILY ; NORMALISE_SV ; UNION_SITES ; SV_SITES ; SV_UNION ; GNOMAD_JOINT ; ANNOTATION_TABLES ; SV_ANNOTATION_TABLES ; FAMILY_TABLES ; groupAnnotation ; svName ; svGroup } from './modules/family_tables'
include { ILLUMINA_ANNOTATE } from './modules/annotator'
include { ILLUMINA_ANNOTATE as ILLUMINA_ANNOTATE_SV } from './modules/annotator'

params {
    // A CSV with the columns family, vcf, ped and proband: each family's id, its VCF with genotypes, its PED file,
    // whose individuals are the VCF's samples, and the sample id of its proband, left empty for a family in which none
    // is designated (docs/annotation.md, the input contract). Nothing derives a proband, so a sheet without the column
    // stops the run.
    // Three more columns are optional, and so is each of their fields: sv_vcf and cnv_vcf, the VCFs of the family's
    // structural-variant caller and of its copy-number caller, whose samples are the family's, and repeat_vcfs, its
    // members' repeat VCFs, one of each sequenced member, separated by semicolons. Each is bgzipped with its .tbi
    // beside it. A family gets the tables of what its row names. A sheet with any other column stops the run.
    families: Path
    // The name of the families' source selection, which names the union and its folder. It is internal: where a
    // reader of the tables is not to learn the selection, pick a name that does not describe it.
    group: String
    // The reference the VCFs are called against, with its .fai beside it; the normalisation left-aligns against it.
    reference: Path
    // As for annotator.nf's annotate mode: the catalog's id, the data directory its download published, the sources to
    // annotate with, the Annotator's worker processes, its image, path and version, and the licence's secret.
    catalog: String
    annotation_data: Path
    sources: List<String> = []
    workers: Integer = 8
    annotator_container: String?
    annotator: String?
    annotator_version: String?
    license_secret: String?
    // An earlier run's annotation of the same union (the folder that run published the Annotator's output in, holding
    // out/<group>.json.gz and run.json), used instead of annotating again; the conversion fails unless its alleles are
    // this union's. For a group with a family that names sv_vcf or cnv_vcf the folder also holds
    // out/<group>.sv.json.gz, the annotation of the union of the group's structural variants, which is reused in the
    // same way.
    annotated: Path?
    // An echtvar archive of gnomAD's joint release, which holds gnomAD's own joint allele count, frequency, homozygote
    // count and filtering allele frequency (docs/annotation.md); echtvar reads it afresh in every run, also one that
    // reuses an annotation.
    gnomad_joint: Path
    tools_container: String
    output_dir: String?
    arGuid: String?
}

workflow {

    main:
    // The Annotator's image is the config's to read (nextflow.config), so the named workflow does not take it.
    if (params.annotated == null && params.sources && params.annotator_container == null) {
        error('give annotator_container, which annotating needs')
    }
    def rows = params.families.splitCsv(header: true)
    if (!rows) {
        error("${params.families}: no families")
    }
    // An optional column under a misspelt name would be read as absent, and its families would get no table of it
    // while the run ended well.
    def unknown = rows[0].keySet() - sheetColumns()
    if (unknown) {
        error("${params.families}: the columns ${unknown} are none of a families sheet, which has ${sheetColumns().join(', ')}")
    }
    def ids = rows.collect { row -> checkFamily(row) }
    if (ids.unique(false).size() != ids.size()) {
        error("${params.families}: a family is listed twice")
    }
    families = channel.fromList(rows).map { row -> tuple(row.family, file(row.vcf), file(row.ped), row.proband) }
    // The VCFs of each family's larger variants with their indexes: an empty list for what its row does not name.
    structural = channel.fromList(rows)
        .map { row ->
            tuple(row.family, indexed([row.sv_vcf]), indexed([row.cnv_vcf]), indexed(repeatVcfs(row)))
        }
    annotated = ANNOTATE_FAMILIES(
        families,
        structural,
        rows.any { row -> row.sv_vcf || row.cnv_vcf },
        [
            group: params.group,
            reference: params.reference,
            catalog: params.catalog,
            annotation_data: params.annotation_data,
            sources: params.sources,
            workers: params.workers,
            annotator: params.annotator,
            annotator_version: params.annotator_version,
            license_secret: params.license_secret,
            annotated: params.annotated,
            gnomad_joint: params.gnomad_joint,
            pipeline_version: pipelineVersion(),
        ],
    )

    publish:
    union = annotated.union
    annotator = annotated.annotator
    family_tables = annotated.family_tables
}

// The steps of the family tables, for a caller that has the families as channels: this script's entry workflow, which
// reads them from a sheet, or another workflow that makes their files in the same run. `families` holds each family as
// tuple(family, vcf, ped, proband), and `structural` each family as tuple(family, sv, cnv, repeats), where each of the
// three is the VCFs with their indexes, or an empty list. `annotates_sv` says whether any family has a structural-variant
// or a copy-number VCF: it decides which processes the run has, so the caller knows it before the run starts. `options`
// holds what the entry workflow takes as params (group, reference, catalog, annotation_data, sources, workers,
// annotator, annotator_version, license_secret, annotated, gnomad_joint), and pipeline_version, the version of this
// pipeline, which is not the caller's `workflow.manifest.version` when another pipeline calls it.
workflow ANNOTATE_FAMILIES {

    take:
    families
    structural
    annotates_sv
    options

    main:
    // The names become file and folder names and go into command lines.
    [options.group, options.catalog].each { v ->
        if (!(v ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
            error("${v}: a group or catalog id may hold letters, digits and . _ - only")
        }
    }
    options.sources.each { v ->
        if (!(v ==~ /[A-Za-z0-9][A-Za-z0-9 ._()-]*/)) {
            error("${v}: a source may hold letters, digits, spaces and . _ ( ) - only")
        }
    }
    if (options.annotated != null) {
        if (options.sources) {
            error('give sources to annotate, or annotated to reuse an annotation, not both')
        }
    }
    else {
        if (!options.sources) {
            error('name the sources to annotate with, or give annotated')
        }
        ['annotator', 'annotator_version', 'license_secret'].each { name ->
            if (options[name] == null) {
                error("give ${name}, which annotating needs")
            }
        }
    }
    def reference_index = options.reference.resolveSibling("${options.reference.name}.fai")
    if (!reference_index.exists()) {
        error("${options.reference}: no .fai beside it")
    }
    def bundle_config = options.annotation_data.resolve('runConfigs/bundle.json')

    // The first step in the tools image, which refuses an image of another release than this workflow's. Every step
    // before the Annotator's reads a family that it has checked.
    checked = FAMILY_CHECK(families.join(structural), options.pipeline_version)
    reference = channel.value(tuple(options.reference, reference_index))
    normalised = NORMALISE_FAMILY(
        families.join(checked.map { f -> [f] }).map { family, vcf, _ped, _proband -> tuple(family, vcf) },
        reference,
    )
    sv_alleles = NORMALISE_SV(
        structural.join(checked.map { f -> [f] }).filter { _family, sv, _cnv, _repeats -> sv }.map { family, sv, _cnv, _repeats -> tuple(family, sv) },
        reference,
    )
    union = UNION_SITES(options.group, normalised.sites.flatten().collect())
    // The union of the group's structural variants and copy-number segments, where a family names a VCF of either.
    if (annotates_sv) {
        sv_sites = SV_SITES(
            structural.join(checked.map { f -> [f] }).filter { _family, sv, cnv, _repeats -> sv || cnv }.map { family, sv, cnv, _repeats -> tuple(family, sv, cnv) }
        )
        sv_group = sv_sites.toList().map { cut -> svGroup(cut) }
        sv_union = SV_UNION(options.group, sv_group.map { g -> g.families }, sv_group.map { g -> g.sites })
    }

    joint = GNOMAD_JOINT(union, options.gnomad_joint)

    if (options.annotated != null) {
        def reused = options.annotated.resolve("out/${options.group}.json.gz")
        if (!reused.exists()) {
            error("${options.group}: ${options.annotated} holds no out/${options.group}.json.gz, the group's annotation")
        }
        json = channel.value(reused)
        run_config = channel.value(file(options.annotated.resolve('run.json'), checkIfExists: true))
        annotator = channel.empty()
        if (annotates_sv) {
            def reused_sv = options.annotated.resolve("out/${svName(options.group)}.json.gz")
            if (!reused_sv.exists()) {
                error("${options.group}: ${options.annotated} holds no out/${svName(options.group)}.json.gz, the annotation of the group's structural variants")
            }
            sv_json = channel.value(reused_sv)
        }
    }
    else {
        annotated = ILLUMINA_ANNOTATE(
            union.map { group, vcf, _index -> tuple(group, vcf) },
            options.annotation_data,
            options.catalog,
            options.sources,
            options.workers,
            options.annotator,
            options.annotator_version,
            options.license_secret,
        )
        json = annotated.annotation.map { outputs -> groupAnnotation(outputs, options.group) }.first()
        run_config = annotated.config
        annotator = annotated.annotation.flatten().mix(annotated.config, annotated.log)
        if (annotates_sv) {
            // The second run has the first one's catalog and sources, so its run config is the same file, and the
            // first run's is the one published and read.
            sv_annotated = ILLUMINA_ANNOTATE_SV(
                sv_union.union.map { _group, vcf, _index -> tuple(svName(options.group), vcf) },
                options.annotation_data,
                options.catalog,
                options.sources,
                options.workers,
                options.annotator,
                options.annotator_version,
                options.license_secret,
            )
            sv_json = sv_annotated.annotation.map { outputs -> groupAnnotation(outputs, svName(options.group)) }
            annotator = annotator.mix(sv_annotated.annotation.flatten(), sv_annotated.log)
        }
    }
    tables = ANNOTATION_TABLES(union.combine(json).combine(run_config), bundle_config, options.catalog).first()
    // A group without structural variants has no annotation of them, and its families read none.
    sv_tables = channel.value([])
    sv_map = channel.value([])
    if (annotates_sv) {
        sv_tables = SV_ANNOTATION_TABLES(sv_union.union.combine(sv_json).combine(run_config), bundle_config, options.catalog).first()
        sv_map = sv_union.map
    }
    // The PED file comes from the sheet rather than through NORMALISE_FAMILY: a task's input passed on as its output
    // is a link into that task's work directory, which another VM cannot follow to the file. The VCFs of the larger
    // variants come from the sheet too. A family without a structural-variant VCF has no normalised sequences.
    written = FAMILY_TABLES(
        normalised.family.join(families.map { family, _vcf, ped, proband -> [family, ped, proband] }).join(structural).join(sv_alleles, remainder: true).map { family, vcf, index, ped, proband, sv, cnv, repeats, alleles ->
            tuple(family, vcf, index, ped, proband, sv, cnv, repeats, alleles ?: [])
        },
        tables,
        joint,
        sv_tables,
        sv_map,
        tablesVersion(),
        options.pipeline_version,
    )
    sv_published = annotates_sv ? sv_union.union.map { _group, vcf, index -> [vcf, index] }.flatten().mix(sv_tables, sv_map) : channel.empty()


    emit:
    union = union.map { _group, vcf, index -> [vcf, index] }.flatten().mix(tables, joint, sv_published)
    annotator = annotator
    family_tables = written.map { family, variants, transcripts, manifest, notes, larger ->
        [family: family, files: [variants, transcripts, manifest, notes] + larger]
    }
}

// A family's tables go under the family's id; the folder names the tables' version and the catalog.
// The union's tables and its gnomAD joint values are folders. Nextflow's default replaces a published file whose
// size or time differs, and compares a folder by its own, which an object store gives every folder alike: a run
// into a folder that an earlier run wrote would keep that run's tables beside the families' new ones. The union is
// therefore always replaced, which removes an earlier folder whole and leaves none of its parts among the new ones.
// The Annotator's output does not change with the tables' version, so its folder names the catalog alone, as the folder
// of an annotation of annotator.nf does, and a run of any version of the tables can take it as annotated.
output {
    union {
        path unionFolder(params.group, params.catalog)
        overwrite true
    }
    annotator {
        path annotatorFolder(params.group, params.catalog)
    }
    family_tables {
        path { r -> familyFolder(r.family, params.catalog) }
    }
}

// The three folders under the output directory, for this script's output block and for that of a pipeline that calls
// ANNOTATE_FAMILIES and publishes what it emits.
def unionFolder(group, catalog) {
    return "annotation-unions/${group}/${tablesFolder(catalog)}"
}

def annotatorFolder(group, catalog) {
    return "annotation-unions/${group}/illumina-${catalog}"
}

def familyFolder(family, catalog) {
    return "annotations/${family}/${tablesFolder(catalog)}"
}

// The version of this pipeline, also in a run that another pipeline launched, whose `workflow.manifest.version` is
// that pipeline's: the manifest's version in the nextflow.config beside this script, read as text.
def pipelineVersion() {
    def found = file("${moduleDir}/nextflow.config").readLines().findAll { line -> line ==~ /    version = '[^']+'/ }
    if (found.size() != 1) {
        error("${moduleDir}/nextflow.config: expected one manifest version, found ${found.size()}")
    }
    return found[0].split("'")[1]
}

// The version of the tables' layout: the one line of src/case_annotation/tables-version.txt, which the tools read too
// (family_tables.SCHEMA_VERSION) and check against the version a task gives them. It changes when a column of a table
// is renamed, retyped or removed or its meaning changes, and not with a table or a column that is added.
def tablesVersion() {
    return file("${moduleDir}/src/case_annotation/tables-version.txt").text.trim()
}

def tablesFolder(catalog) {
    return "tables-${tablesVersion()}-illumina-${catalog}"
}

// The columns of the families sheet: four that every sheet has, and the optional ones of the larger variants.
def sheetColumns() {
    return ['family', 'vcf', 'ped', 'proband', 'sv_vcf', 'cnv_vcf', 'repeat_vcfs']
}

// A row of the families sheet: its family id, which names the family's outputs, its two files, and its proband's
// sample id, which goes into command lines. The id may be empty, for a family in which no proband is designated, but
// the sheet needs the column and the row the field: read as empty, their absence would pass for such a family.
def checkFamily(row) {
    ['family', 'vcf', 'ped'].each { column ->
        if (!row[column]) {
            error("${params.families}: a row without ${column}: ${row}")
        }
    }
    if (!row.containsKey('proband')) {
        error("${params.families}: no column proband, which holds each family's proband, or nothing where none is designated")
    }
    // splitCsv gives null for a field the row lacks and an empty string for an empty one.
    if (row.proband == null) {
        error("${params.families}: a row that ends before proband, which is left empty where none is designated: ${row}")
    }
    if (!(row.family ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error("${row.family}: a family id may hold letters, digits and . _ - only")
    }
    if (row.proband && !(row.proband ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error("${row.family}: the proband ${row.proband} is not a sample id, which may hold letters, digits and . _ - only")
    }
    // The VCFs of the larger variants are optional, so nothing is refused for their absence. A named one is read with
    // its index, and a missing index would otherwise fail a task after the families' checks.
    ([row.sv_vcf, row.cnv_vcf] + repeatVcfs(row))
        .findAll { vcf -> vcf }
        .each { vcf ->
            if (!file("${vcf}.tbi").exists()) {
                error("${row.family}: ${vcf}: no .tbi beside it")
            }
        }
    return row.family
}

// The repeat VCFs a row names, one of each sequenced member of the family, separated by semicolons; none where the
// sheet has no such column or the row leaves it empty.
def repeatVcfs(row) {
    if (!row.repeat_vcfs) {
        return []
    }
    def vcfs = row.repeat_vcfs.split(';', -1).collect { vcf -> vcf.trim() }
    if (vcfs.any { vcf -> !vcf }) {
        error("${row.family}: an empty entry among its repeat_vcfs, which are separated by single semicolons: ${row.repeat_vcfs}")
    }
    return vcfs
}

// Each named VCF with the index beside it, as a process stages them; nothing for a VCF the row does not name.
def indexed(vcfs) {
    return vcfs.findAll { vcf -> vcf }.collectMany { vcf -> [file(vcf), file("${vcf}.tbi")] }
}
