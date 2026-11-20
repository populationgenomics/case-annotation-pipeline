// A pipeline that includes ANNOTATE_FAMILIES, as another project does (docs/annotation.md, Including the workflow in
// another pipeline). It lies in a folder of its own, so that projectDir is not the repository's root, and the run has
// the config beside it and nothing of the repository's nextflow.config. It makes one family's files in a process and
// gives the family to the workflow as a channel value, with no families sheet, and reuses an annotation. The stub-runs
// job runs it in stub mode (.github/workflows/nextflow.yml):
//
//   nextflow -C tests/workflows/include/nextflow.config run tests/workflows/include/main.nf -stub-run ...

include { ANNOTATE_FAMILIES ; unionFolder ; annotatorFolder ; familyFolder } from '../../../family_tables.nf'

params {
    // The family's id and its proband's sample id, and whether its files hold a structural-variant VCF. A run without
    // proband gives the workflow a family without one, which it refuses.
    family: String
    proband: String?
    with_sv: Boolean = false
    // What ANNOTATE_FAMILIES takes beside the family: whether the run has the steps of the structural variants, and
    // its settings. The run reuses an annotation, so it names no source and nothing of the Annotator's. It can name
    // any of the archives of AlphaGenome's scores and any of the public files that are joined to the alleles.
    annotates_sv: Boolean = false
    group: String
    reference: Path
    catalog: String
    annotation_data: Path
    annotated: Path
    gnomad_joint: Path
    alphagenome_avi: Path?
    alphagenome_avi_breakdown: Path?
    alphagenome_splicing: Path?
    gnomad_mito: Path?
    gencode_polya: Path?
    aparent2: Path?
    // A name to add to the settings, and the name of a setting to give as null, for the runs that ANNOTATE_FAMILIES
    // refuses for them.
    unknown_setting: String?
    null_setting: String?
    output_dir: String
}

workflow {

    main:
    def settings = [
        group: params.group,
        reference: params.reference,
        catalog: params.catalog,
        annotation_data: params.annotation_data,
        sources: [],
        workers: 8,
        annotator: null,
        annotator_version: null,
        license_secret: null,
        annotated: params.annotated,
        gnomad_joint: params.gnomad_joint,
        alphagenome_avi: params.alphagenome_avi,
        alphagenome_avi_breakdown: params.alphagenome_avi_breakdown,
        alphagenome_splicing: params.alphagenome_splicing,
        gnomad_mito: params.gnomad_mito,
        gencode_polya: params.gencode_polya,
        aparent2: params.aparent2,
    ]
    families = MADE_FAMILY(params.with_sv).map { vcf, ped, sv -> tuple(params.family, vcf, ped, params.proband, sv, [], []) }
    annotated = ANNOTATE_FAMILIES(
        families,
        params.annotates_sv,
        settings + (params.unknown_setting ? [(params.unknown_setting): true] : [:]) + (params.null_setting ? [(params.null_setting): null] : [:]),
    )

    publish:
    union = annotated.union
    annotator = annotated.annotator
    family_tables = annotated.family_tables
}

// The folders of family_tables.nf's own output block, and its replacing of the union.
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

// Stands for the steps of an including pipeline that write a family's files during the run. The files are empty, as
// the stub runs' are.
process MADE_FAMILY {
    input:
    val with_sv

    output:
    tuple path('made.vcf.gz'), path('made.ped'), path('made.sv.vcf.gz{,.tbi}', arity: '0..2')

    script:
    """
    touch made.vcf.gz made.ped ${with_sv ? 'made.sv.vcf.gz made.sv.vcf.gz.tbi' : ''}
    """
}
