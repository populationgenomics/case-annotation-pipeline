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
//
// The steps are the workflow ANNOTATE_FAMILIES, which takes the families as a channel and its settings as values, so
// that another pipeline can include it and give it families that it makes in the same run (docs/annotation.md,
// Including the workflow in another pipeline). The entry workflow reads the families from a sheet and the settings
// from params, and publishes what ANNOTATE_FAMILIES emits.

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
    // The config reads the Annotator's image (nextflow.config), so ANNOTATE_FAMILIES does not take it.
    if (params.annotated == null && params.annotator_container == null) {
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
    // The VCFs of each family's larger variants with their indexes: an empty list for what its row does not name.
    families = channel.fromList(rows)
        .map { row ->
            tuple(row.family, file(row.vcf), file(row.ped), row.proband, indexed([row.sv_vcf]), indexed([row.cnv_vcf]), indexed(repeatVcfs(row)))
        }
    annotated = ANNOTATE_FAMILIES(
        families,
        rows.any { row -> row.sv_vcf || row.cnv_vcf },
        settingNames().collectEntries { name -> [name, params[name]] },
    )

    publish:
    union = annotated.union
    annotator = annotated.annotator
    family_tables = annotated.family_tables
}

// The steps of the family tables, for the entry workflow above and for a pipeline that includes this workflow.
//
// `families` holds each family once, as tuple(family, vcf, ped, proband, sv, cnv, repeats): its id, its VCF, its PED
// file, its proband's sample id (empty where none is designated), and the VCFs of its structural-variant caller, of
// its copy-number caller and of its members' repeats. Each of the last three is a list of the VCFs with their .tbi
// indexes, and an empty list for what the family does not have.
//
// `annotates_sv` is a Boolean, true where any family has a structural-variant or a copy-number VCF. The steps of their
// union are part of the run or not when the workflow is built, before a family arrives, so a caller knows the value
// when the run starts.
//
// `settings` is a map of the values that the entry workflow takes as the params of the same names (settingNames),
// each with the type that its param declares. They are known when the run starts too: the workflow checks them and
// looks for the reference's index and a reused annotation's files before any task runs.
//
// It emits what the output block below publishes: `union`, each file and folder of the group's two unions;
// `annotator`, each file of the Annotator's output, and nothing in a run that reuses an annotation; and
// `family_tables`, a record of each family with its id as `family` and the files of its folder as `files`.
workflow ANNOTATE_FAMILIES {
    take:
    families
    annotates_sv
    settings

    main:
    // A setting under a misspelt name would be read as absent, and a run that was to reuse an annotation would
    // annotate again.
    if (settings.keySet() != settingNames().toSet()) {
        error("ANNOTATE_FAMILIES takes the settings ${settingNames().join(', ')}, and was given ${settings.keySet().join(', ')}")
    }
    // The names become file and folder names and go into command lines.
    [settings.group, settings.catalog].each { v ->
        if (!(v ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
            error("${v}: a group or catalog id may hold letters, digits and . _ - only")
        }
    }
    settings.sources.each { v ->
        if (!(v ==~ /[A-Za-z0-9][A-Za-z0-9 ._()-]*/)) {
            error("${v}: a source may hold letters, digits, spaces and . _ ( ) - only")
        }
    }
    if (settings.annotated != null) {
        if (settings.sources) {
            error('give sources to annotate, or annotated to reuse an annotation, not both')
        }
    }
    else {
        if (!settings.sources) {
            error('name the sources to annotate with, or give annotated')
        }
        ['annotator', 'annotator_version', 'license_secret'].each { name ->
            if (settings[name] == null) {
                error("give ${name}, which annotating needs")
            }
        }
    }
    def reference_index = settings.reference.resolveSibling("${settings.reference.name}.fai")
    if (!reference_index.exists()) {
        error("${settings.reference}: no .fai beside it")
    }
    def bundle_config = settings.annotation_data.resolve('runConfigs/bundle.json')

    // The first step in the tools image, which refuses an image of another release than this workflow's. A family is
    // checked as it arrives too (checkRecord), and every step before the Annotator's reads a family that has passed
    // both.
    checked = FAMILY_CHECK(families.map { record -> checkRecord(record, annotates_sv) }, pipelineVersion())
    passed = families.join(checked.map { f -> [f] })
    reference = channel.value(tuple(settings.reference, reference_index))
    normalised = NORMALISE_FAMILY(
        passed.map { family, vcf, _ped, _proband, _sv, _cnv, _repeats -> tuple(family, vcf) },
        reference,
    )
    sv_alleles = NORMALISE_SV(
        passed.filter { _family, _vcf, _ped, _proband, sv, _cnv, _repeats -> sv }.map { family, _vcf, _ped, _proband, sv, _cnv, _repeats -> tuple(family, sv) },
        reference,
    )
    union = UNION_SITES(settings.group, normalised.sites.flatten().collect())
    // The union of the group's structural variants and copy-number segments, where a family has a VCF of either.
    if (annotates_sv) {
        sv_sites = SV_SITES(
            passed.filter { _family, _vcf, _ped, _proband, sv, cnv, _repeats -> sv || cnv }.map { family, _vcf, _ped, _proband, sv, cnv, _repeats -> tuple(family, sv, cnv) }
        )
        sv_group = sv_sites.toList().map { cut -> svGroup(checkSvSites(cut)) }
        sv_union = SV_UNION(settings.group, sv_group.map { g -> g.families }, sv_group.map { g -> g.sites })
    }

    joint = GNOMAD_JOINT(union, settings.gnomad_joint)

    if (settings.annotated != null) {
        def reused = settings.annotated.resolve("out/${settings.group}.json.gz")
        if (!reused.exists()) {
            error("${settings.group}: ${settings.annotated} holds no out/${settings.group}.json.gz, the group's annotation")
        }
        json = channel.value(reused)
        run_config = channel.value(file(settings.annotated.resolve('run.json'), checkIfExists: true))
        annotator = channel.empty()
        if (annotates_sv) {
            def reused_sv = settings.annotated.resolve("out/${svName(settings.group)}.json.gz")
            if (!reused_sv.exists()) {
                error("${settings.group}: ${settings.annotated} holds no out/${svName(settings.group)}.json.gz, the annotation of the group's structural variants")
            }
            sv_json = channel.value(reused_sv)
        }
    }
    else {
        annotated = ILLUMINA_ANNOTATE(
            union.map { group, vcf, _index -> tuple(group, vcf) },
            settings.annotation_data,
            settings.catalog,
            settings.sources,
            settings.workers,
            settings.annotator,
            settings.annotator_version,
            settings.license_secret,
        )
        json = annotated.annotation.map { outputs -> groupAnnotation(outputs, settings.group) }.first()
        run_config = annotated.config
        annotator = annotated.annotation.flatten().mix(annotated.config, annotated.log)
        if (annotates_sv) {
            // The second run has the first one's catalog and sources, so its run config is the same file, and the
            // first run's is the one published and read.
            sv_annotated = ILLUMINA_ANNOTATE_SV(
                sv_union.union.map { _group, vcf, _index -> tuple(svName(settings.group), vcf) },
                settings.annotation_data,
                settings.catalog,
                settings.sources,
                settings.workers,
                settings.annotator,
                settings.annotator_version,
                settings.license_secret,
            )
            sv_json = sv_annotated.annotation.map { outputs -> groupAnnotation(outputs, svName(settings.group)) }
            annotator = annotator.mix(sv_annotated.annotation.flatten(), sv_annotated.log)
        }
    }
    tables = ANNOTATION_TABLES(union.combine(json).combine(run_config), bundle_config, settings.catalog).first()
    // A group without structural variants has no annotation of them, and its families read none.
    sv_tables = channel.value([])
    sv_map = channel.value([])
    if (annotates_sv) {
        sv_tables = SV_ANNOTATION_TABLES(sv_union.union.combine(sv_json).combine(run_config), bundle_config, settings.catalog).first()
        sv_map = sv_union.map
    }
    // The PED file comes from the family as it was given rather than through NORMALISE_FAMILY: a task's input passed on
    // as its output is a link into that task's work directory, which another VM cannot follow to the file. The VCFs of
    // the larger variants come from there too. A family without a structural-variant VCF has no normalised sequences.
    written = FAMILY_TABLES(
        normalised.family.join(families.map { family, _vcf, ped, proband, sv, cnv, repeats -> [family, ped, proband, sv, cnv, repeats] }).join(sv_alleles, remainder: true).map { family, vcf, index, ped, proband, sv, cnv, repeats, alleles ->
            tuple(family, vcf, index, ped, proband, sv, cnv, repeats, alleles ?: [])
        },
        tables,
        joint,
        sv_tables,
        sv_map,
        tablesVersion(),
        pipelineVersion(),
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

// The three folders under the output directory, for the output block above and for that of a pipeline that includes
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

def tablesFolder(catalog) {
    return "tables-${tablesVersion()}-illumina-${catalog}"
}

// The version of the tables' layout: the one line of src/case_annotation/tables-version.txt, which the tools read too
// (family_tables.SCHEMA_VERSION) and check against the version a task gives them. It changes when a column of a table
// is renamed, retyped or removed or its meaning changes, and not with a table or a column that is added. The file is
// found from this script's folder, since projectDir is the including pipeline's where one includes this script.
def tablesVersion() {
    return file("${moduleDir}/src/case_annotation/tables-version.txt").text.trim()
}

// The version of this pipeline: the one line of pipeline-version.txt, which the tools check against their own
// (family_tables.PIPELINE_VERSION). `workflow.manifest.version` is the launched pipeline's, which is another one where
// a pipeline includes this script, so no step is given it.
def pipelineVersion() {
    return file("${moduleDir}/pipeline-version.txt").text.trim()
}

// The settings of ANNOTATE_FAMILIES, each of which the entry workflow takes as the param of the same name.
def settingNames() {
    return ['group', 'reference', 'catalog', 'annotation_data', 'sources', 'workers', 'annotator', 'annotator_version', 'license_secret', 'annotated', 'gnomad_joint']
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
    checkIds(row.family, row.proband)
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

// A family's id and its proband's sample id, which go into command lines; the id names the family's outputs too.
def checkIds(family, proband) {
    if (!(family ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error("${family}: a family id may hold letters, digits and . _ - only")
    }
    if (proband && !(proband ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error("${family}: the proband ${proband} is not a sample id, which may hold letters, digits and . _ - only")
    }
}

// A family as ANNOTATE_FAMILIES takes it, checked when it arrives: the entry workflow has checked its sheet before the
// run starts, and an including pipeline's families come during the run. A structural-variant or copy-number VCF in a
// run that was built without the steps of their union would reach FAMILY_TABLES with no annotation to read
// (checkSvSites has the opposite case).
def checkRecord(record, annotates_sv) {
    def (family, _vcf, _ped, proband, sv, cnv, _repeats) = record
    checkIds(family, proband)
    if (!annotates_sv && (sv || cnv)) {
        error("${family}: a structural-variant or copy-number VCF in a run with annotates_sv false, which has no step that annotates it")
    }
    return record
}

// The sites of the families with a structural-variant or a copy-number VCF, in a run that was built with the steps of
// their union. Without such a family the union would be empty, and the Annotator would still run over it.
def checkSvSites(sites) {
    if (!sites) {
        error('annotates_sv is true, and no family has a structural-variant or a copy-number VCF')
    }
    return sites
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
