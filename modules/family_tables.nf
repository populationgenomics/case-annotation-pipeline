// The steps of the family tables (family_tables.nf, docs/annotation.md) other than the Annotator's own: checking a
// family's inputs, normalising its VCF and the sequences of its structural variants, taking the union of the families'
// sites, reading the union's values from echtvar archives, joining public files to the union's alleles, applying the
// rules of two plugins of Ensembl VEP to the union's alleles in exons, converting the union's annotation, finding each
// family's runs of homozygosity and regions of uniparental disomy, and writing each family's tables. The structural
// variants and copy-number segments of a group have a union of their own, which SV_SITES cuts from each family's VCFs,
// SV_UNION joins and SV_ANNOTATION_TABLES converts the annotation of.
// NORMALISE_FAMILY, NORMALISE_SV, UNION_SITES and SV_SITES run bcftools alone, ECHTVAR_VALUES echtvar and
// case-annotation, SV_UNION, JOINED_VALUES and VEP_PLUGIN_TABLE case-annotation and bcftools, ROH_UPD bcftools,
// echtvar, upd and case-annotation, and the others case-annotation, all in the tools image (label 'tools').
//
// A family's structural-variant VCF, its copy-number VCF and its members' repeat VCFs are optional. A process takes
// each as the files the sheet names with their indexes, staged in a folder of their own, or as an empty list.

// Fails unless the PED file is the family's and names every sample of its VCF, each sample can name a column of the
// variants table (no two names, and no name and a column's of variants, sv or repeats, that differ only in case), and
// the proband the sheet names, if any, is an individual of the PED file and a sample of the VCF, before anything is
// annotated. The structural-variant and copy-number VCFs must hold exactly the VCF's samples, and the repeat VCFs one
// of them each and all of them together. It is the first step of a run in the tools image, so it also refuses an
// image of another release of the pipeline than the workflow's (pipeline_version), before the Annotator runs.
process FAMILY_CHECK {
    tag "${family}"
    label 'tools'

    input:
    tuple val(family), path(vcf, stageAs: 'input/*'), path(ped, stageAs: 'input/*'), val(proband), path(sv, stageAs: 'sv/*'), path(cnv, stageAs: 'cnv/*'), path(repeats, stageAs: 'repeats/*')
    val pipeline_version

    output:
    val family

    script:
    """
    case-annotation family-check --family ${family} --ped ${ped} --vcf ${vcf} ${probandOption(proband)} \\
        ${structuralOptions(sv, cnv, repeats)} --pipeline-version ${pipeline_version}
    """

    stub:
    """
    true
    """
}

// The family's VCF on the primary chromosomes, one allele per record, left-aligned against the reference, each allele
// carried by some member (a copy of it in the genotype, also beside a missing allele, as in ./1), with no spanning
// deletion (*) or symbolic allele; and its sites without samples, IDs, QUAL, FILTER or INFO. The sites and the calls
// come from the same normalised records, so their alleles match. An allele that two records hold after normalising
// fails: each would give its members two calls, and keeping one record would drop the other's calls.
process NORMALISE_FAMILY {
    tag "${family}"
    label 'tools'

    input:
    tuple val(family), path(vcf, stageAs: 'input/*')
    tuple path(reference), path(reference_index)

    output:
    tuple val(family), path("${family}.vcf.gz"), path("${family}.vcf.gz.tbi"), emit: family
    tuple path("${family}.sites.vcf.gz"), path("${family}.sites.vcf.gz.tbi"), emit: sites

    script:
    """
    set -o pipefail
    chromosomes=\$(printf 'chr%s,' {1..22} X Y M)
    bcftools view --no-version -t "\${chromosomes%,}" -Ou ${vcf} \\
        | bcftools norm --no-version -m -any -f ${reference} --check-ref e -Ou \\
        | bcftools view --no-version -i 'GT~"1" && ALT!="*" && ALT!~"^<"' -Oz -o ${family}.vcf.gz
    bcftools index -t ${family}.vcf.gz
    repeated=\$(bcftools query -f '%CHROM:%POS %REF>%ALT\\n' ${family}.vcf.gz | sort | uniq -d | awk 'NR <= 5')
    if [[ -n \${repeated} ]]; then
        echo "${family}: alleles that more than one record holds after normalising: \${repeated}" >&2
        exit 1
    fi
    records=\$(bcftools index -n ${family}.vcf.gz)
    echo "\${records} alleles"
    if [[ \${records} -eq 0 ]]; then
        echo "${family}: no allele on the primary chromosomes; are its contigs named chr1, chr2, ...?" >&2
        exit 1
    fi
    bcftools view --no-version -G -Ou ${family}.vcf.gz \\
        | bcftools annotate --no-version -x ID,QUAL,FILTER,INFO -Oz -o ${family}.sites.vcf.gz
    bcftools index -t ${family}.sites.vcf.gz
    """

    stub:
    """
    touch ${family}.vcf.gz ${family}.vcf.gz.tbi ${family}.sites.vcf.gz ${family}.sites.vcf.gz.tbi
    """
}

// The family's structural variants that are written as a sequence, on the primary chromosomes, normalised as
// NORMALISE_FAMILY normalises the family's VCF: each under its ID, without samples, QUAL, FILTER or INFO. A record whose
// normalised allele is a row of the family's variants table is that allele, and FAMILY_TABLES links the two. Symbolic
// alleles and breakends, which have no sequence to compare, are left out.
process NORMALISE_SV {
    tag "${family}"
    label 'tools'

    input:
    tuple val(family), path(sv, stageAs: 'sv/*')
    tuple path(reference), path(reference_index)

    output:
    tuple val(family), path("${family}.sv-alleles.vcf.gz")

    script:
    """
    set -o pipefail
    chromosomes=\$(printf 'chr%s,' {1..22} X Y M)
    bcftools view --no-version -G -t "\${chromosomes%,}" -e 'ALT~"^<" || ALT~"[][]"' -Ou ${vcfs(sv)[0]} \\
        | bcftools annotate --no-version -x QUAL,FILTER,INFO -Ou \\
        | bcftools norm --no-version -m -any -f ${reference} --check-ref e -Oz -o ${family}.sv-alleles.vcf.gz
    echo "\$(bcftools view --no-version -H ${family}.sv-alleles.vcf.gz | wc -l) structural variants written as a sequence"
    """

    stub:
    """
    touch ${family}.sv-alleles.vcf.gz
    """
}

// The sites of a family's structural variants and copy-number segments: the records of each of its callers on the
// primary chromosomes, without the copy-number segments at the reference's copy number, under the callers' IDs, without
// samples, QUAL or FILTER, and with INFO cut down to the six fields the union keeps. bcftools cuts them and leaves
// INFO END as the caller wrote it: the Annotator takes a record's extent from it, and annotates a deletion without it
// as a single base. tests/workflows/family_tables.nf compares the END of each record of the union with its caller's.
process SV_SITES {
    tag "${family}"
    label 'tools'

    input:
    tuple val(family), path(sv, stageAs: 'sv/*'), path(cnv, stageAs: 'cnv/*')

    output:
    tuple val(family), path("${family}.*.sites.vcf.gz")

    script:
    def callers = vcfs(sv).collect { f -> "sv ${f}" } + vcfs(cnv).collect { f -> "cnv ${f}" }
    """
    set -o pipefail
    chromosomes=\$(printf 'chr%s,' {1..22} X Y M)
    cut_sites() {
        # A caller's records with an alternate allele: its name, its VCF.
        local kept=(-t "\${chromosomes%,}" -e 'ALT="."')
        bcftools view --no-version -G "\${kept[@]}" -Ou "\$2" \\
            | bcftools annotate --no-version -x 'QUAL,FILTER,^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/CIPOS,INFO/CIEND,INFO/IMPRECISE' \\
                -Oz -o ${family}.\$1.sites.vcf.gz
        echo "\$(bcftools view --no-version -H ${family}.\$1.sites.vcf.gz | wc -l) records of \$2"
    }
    ${callers.collect { c -> "cut_sites ${c}" }.join('\n    ')}
    """

    stub:
    def names = (vcfs(sv) ? ['sv'] : []) + (vcfs(cnv) ? ['cnv'] : [])
    """
    touch ${names.collect { n -> "${family}.${n}.sites.vcf.gz" }.join(' ')}
    """
}

// The union of a group's sites of structural variants and copy-number segments: one record per allele, a record that
// several families hold written once, each under an id of its content, and the map from the families' records to the
// sites. case-annotation writes the records as text and bcftools compresses them. A site's id is a digest of its
// content, END among it, and SV_ANNOTATION_TABLES refuses a record of the compressed union whose ID is not that of its
// content. Like the union of the alleles, it reveals every family's records and stays internal.
process SV_UNION {
    tag "${group}"
    label 'tools'

    input:
    val group
    val families
    path sites, stageAs: 'sites/*'

    output:
    tuple val(group), path("${group}.sv.vcf.gz"), path("${group}.sv.vcf.gz.tbi"), emit: union
    path "${group}.sv.map.parquet", emit: map

    script:
    """
    set -o pipefail
    case-annotation sv-union ${families.collect { family, name -> "--sites ${family} sites/${name}" }.join(' ')} \\
        --output-vcf union.vcf --output-map ${group}.sv.map.parquet
    bcftools view --no-version -Oz -o ${group}.sv.vcf.gz union.vcf
    bcftools index -t ${group}.sv.vcf.gz
    """

    stub:
    """
    touch ${group}.sv.vcf.gz ${group}.sv.vcf.gz.tbi ${group}.sv.map.parquet
    """
}

// The union of the families' sites, each allele once. It reveals every family's alleles, so it stays internal.
process UNION_SITES {
    tag "${group}"
    label 'tools'

    input:
    val group
    path sites, stageAs: 'sites/*'

    output:
    tuple val(group), path("${group}.vcf.gz"), path("${group}.vcf.gz.tbi")

    script:
    """
    bcftools concat --no-version -a -d exact -Oz -o ${group}.vcf.gz \$(ls sites/*.vcf.gz | sort)
    bcftools index -t ${group}.vcf.gz
    echo "\$(bcftools index -n ${group}.vcf.gz) alleles in the union"
    """

    stub:
    """
    touch ${group}.vcf.gz ${group}.vcf.gz.tbi
    """
}

// The values of the union's alleles that come from echtvar archives: gnomAD's joint values, from an archive of
// gnomAD's joint release, and AlphaGenome's scores, from up to three archives that a run need not have
// (docs/annotation.md). An archive that the run does not have is an empty list. One echtvar call adds the fields of
// every archive to each record, matching position, REF and ALT exactly; it takes one allele per record, as the union
// has them. case-annotation writes the values as a table keyed by variant id, and records each archive's size, SHA-256
// and fields with echtvar's version. It fails unless the table's alleles are exactly the union's, and where two
// archives add a field of the same name, which echtvar does not check. The archives are copied to the task's disk,
// which archivesDisk sizes from them. The annotated BCF gets an index only because htslib reports a file it opens
// without one.
// The stub gives its folder a fixed time, as an object store gives every folder: a second stub run's folder then
// looks to Nextflow like the one the first run published, which is what a run into a bucket meets.
process ECHTVAR_VALUES {
    tag "${group}"
    label 'tools'
    machineType 'n2d-standard-4'
    cpus 4
    memory 16.GB
    disk { archivesDisk([gnomad_joint, alphagenome_avi, alphagenome_avi_breakdown, alphagenome_splicing]) }

    input:
    tuple val(group), path(sites), path(sites_index)
    path gnomad_joint, stageAs: 'gnomad_joint/*'
    path alphagenome_avi, stageAs: 'alphagenome_avi/*'
    path alphagenome_avi_breakdown, stageAs: 'alphagenome_avi_breakdown/*'
    path alphagenome_splicing, stageAs: 'alphagenome_splicing/*'

    output:
    path 'echtvar_values'

    script:
    // Each archive of the run under the name that case-annotation knows it by, which is its param's.
    def archives = [gnomad_joint: gnomad_joint, alphagenome_avi: alphagenome_avi, alphagenome_avi_breakdown: alphagenome_avi_breakdown, alphagenome_splicing: alphagenome_splicing].findAll { _source, archive -> archive }
    """
    echo "${archives.size()} archives, on a disk of ${task.disk}"
    version=\$(echtvar --version)
    echtvar anno ${archives.collect { _source, archive -> "-e ${archive}" }.join(' ')} ${sites} annotated.bcf
    bcftools index annotated.bcf
    case-annotation echtvar-values --annotated annotated.bcf --sites ${sites} \\
        ${archives.collect { source, archive -> "--archive ${source} ${archive}" }.join(' ')} \\
        --echtvar-version "\${version#echtvar }" --output echtvar_values
    """

    stub:
    """
    mkdir echtvar_values
    touch echtvar_values/echtvar_values.parquet echtvar_values/echtvar_values.json
    touch -t 200001010000 echtvar_values
    """
}

// The values of the union's alleles that come from public files which a run joins by allele or by position, each of
// which a run need not have (docs/annotation.md): gnomAD's sites VCF of the mitochondrial genome, MitoTIP's table of
// scores, GENCODE's polyadenylation features, APARENT2's table of effects on polyadenylation, gnomAD's tables of
// genomic constraint by window, of regions of missense constraint and of MPC scores, and the prepared catalogue of
// branch points. A file that the run does not have is an empty list. gnomAD's VCF is matched by position, REF and
// ALT, so bcftools first normalises it as NORMALISE_FAMILY normalises a family's VCF, one allele per record and
// left-aligned against the reference, and stops at a REF that is not the reference's. case-annotation writes the
// values as a table keyed by variant id, and those that a file gives for an allele and one transcript as a second
// table keyed by variant id and transcript id, and records each file's name, size and SHA-256 as the run was given
// it. It fails for a file that is not what its source publishes, and unless the first table's alleles are exactly the
// union's. The step runs in every run, also in one without a file, whose first table then holds the alleles alone and
// whose second has no row. The normalised VCF gets an index only because htslib reports a file it opens without one.
// The stub gives its folder a fixed time, as the stub of ECHTVAR_VALUES does.
process JOINED_VALUES {
    tag "${group}"
    label 'tools'
    machineType 'n2d-standard-4'
    cpus 4
    memory 16.GB
    disk 50.GB

    input:
    tuple val(group), path(sites), path(sites_index)
    tuple path(reference), path(reference_index)
    path gnomad_mito, stageAs: 'gnomad_mito/*'
    path mitotip, stageAs: 'mitotip/*'
    path gencode_polya, stageAs: 'gencode_polya/*'
    path aparent2, stageAs: 'aparent2/*'
    path genomic_constraint, stageAs: 'genomic_constraint/*'
    path mcr, stageAs: 'mcr/*'
    path mpc, stageAs: 'mpc/*'
    path branch_points, stageAs: 'branch_points/*'

    output:
    path 'joined_values'

    script:
    // Each file of the run under the name that case-annotation knows it by, which is its param's.
    def given = [gnomad_mito: gnomad_mito, mitotip: mitotip, gencode_polya: gencode_polya, aparent2: aparent2, genomic_constraint: genomic_constraint, mcr: mcr, mpc: mpc, branch_points: branch_points].findAll { _source, joined -> joined }
    // The files that are VCFs of alleles, each of which is read as bcftools norm rewrote it.
    def allele_vcfs = given.findAll { source, _joined -> source == 'gnomad_mito' }
    """
    set -o pipefail
    echo "${given.size()} files"
    normalised() {
        # A VCF of alleles as the union's are: its source's name, the VCF.
        bcftools norm --no-version -m -any -f ${reference} --check-ref e -Oz -o "\$1.normalised.vcf.gz" "\$2"
        bcftools index -t "\$1.normalised.vcf.gz"
    }
    ${allele_vcfs.collect { source, joined -> "normalised ${source} ${joined}" }.join('\n    ')}
    case-annotation joined-values --sites ${sites} \\
        ${given.collect { source, joined -> "--file ${source} ${joined}" }.join(' ')} \\
        ${allele_vcfs.collect { source, _joined -> "--normalised ${source} ${source}.normalised.vcf.gz" }.join(' ')} \\
        --output joined_values --scratch /tmp
    """

    stub:
    """
    mkdir joined_values
    touch joined_values/joined_values.parquet joined_values/joined_transcript_values.parquet joined_values/joined_values.json
    touch -t 200001010000 joined_values
    """
}

// What the rules of VEP's UTRannotator and NMD plugins say of the union's alleles, as a table by allele and
// transcript, with a record of what the values come from. case-annotation applies the rules to the transcripts of the
// prepared gene models (case-annotation gene-models) and to the reference, and the step runs no VEP.
// bcftools first cuts the union to the alleles that a rule can evaluate: those that overlap an exon of a gene model
// file, each under its variant id, with the positive controls of the prepared gene models among them. A variant in a
// 5' UTR and one that ends a reading frame lie in an exon, so the cut leaves out no allele that a rule would evaluate,
// and the command reads a few alleles in a hundred. The regions do not overlap one another, and an allele that the cut
// still writes twice stops the step, since an allele's id is the key of its rows in the table.
// The command fails unless each positive control holds what its rule gives it: gene models and a reference that do
// not fit one another would otherwise leave every row empty, like a union without a variant in a 5' UTR.
// A staged input can be a link into a bucket mount, as on Google Batch, and reading the reference and the gene models
// through such a mount, in many small reads, is slow. The step copies both into a folder of the task's own on the
// VM's disk, once, in one sequential read each, and the commands are given the copies and read nothing else of the two.
// The stub gives its folder a fixed time, as the stub of ECHTVAR_VALUES does.
process VEP_PLUGIN_TABLE {
    tag "${group}"
    label 'tools'
    machineType 'n2d-standard-8'
    cpus 8
    memory 16.GB
    disk 50.GB

    input:
    tuple val(group), path(sites), path(sites_index)
    tuple path(reference), path(reference_index)
    path gene_models, stageAs: 'gene-models'

    output:
    path 'vep_plugins'

    script:
    """
    set -o pipefail
    copies=\$(mktemp -d /tmp/plugin-rules.XXXXXX)
    trap 'rm -rf "\${copies}"' EXIT
    mkdir "\${copies}/gene-models"
    cp -L ${reference} ${reference_index} "\${copies}/"
    cp -L gene-models/* "\${copies}/gene-models/"
    bcftools view --no-version -R "\${copies}/gene-models/exons.bed" -Ou ${sites} \\
        | bcftools annotate --no-version --set-id '%CHROM-%POS-%REF-%ALT' -Ov -o exons.vcf
    echo "\$(grep -vc '^#' exons.vcf) alleles of the union overlap an exon"
    grep -v '^#' "\${copies}/gene-models/controls.vcf" >> exons.vcf
    bcftools sort -Oz -o ${group}.exons.vcf.gz exons.vcf
    bcftools index -t ${group}.exons.vcf.gz
    repeated=\$(bcftools query -f '%ID\\n' ${group}.exons.vcf.gz | sort | uniq -d | awk 'NR <= 5')
    if [[ -n \${repeated} ]]; then
        echo "${group}: alleles that the cut to the exons wrote more than once: \${repeated}" >&2
        exit 1
    fi
    case-annotation vep-plugins --given ${group}.exons.vcf.gz --gene-models "\${copies}/gene-models" \\
        --reference "\${copies}/${reference}" --uorf-evidence "\${UORF_EVIDENCE}" --workers ${task.cpus} \\
        --output vep_plugins
    """

    stub:
    """
    mkdir vep_plugins
    touch vep_plugins/vep_plugins.parquet vep_plugins/vep_plugins.json
    touch -t 200001010000 vep_plugins
    """
}

// The union's variants and transcripts tables from its Annotator JSON, checked against its sites. The JSON is read
// twice: its genes section, which ends it, gives the HGNC ids that the rows of its positions hold. CPU-bound: one
// process inflates the JSON and the others parse its lines.
// The stub gives its folder a fixed time, as the stub of ECHTVAR_VALUES does.
process ANNOTATION_TABLES {
    tag "${group}"
    label 'tools'
    machineType 'n2d-standard-16'
    cpus 16
    memory 64.GB
    disk 200.GB

    input:
    tuple val(group), path(sites), path(sites_index), path(annotation, stageAs: 'annotation/*'), path(run_config, stageAs: 'run.json')
    path bundle_config, stageAs: 'bundle.json'
    val catalog

    output:
    path 'tables'

    script:
    """
    case-annotation annotation-tables --json ${annotation} --sites ${sites} --run-config ${run_config} \\
        --bundle-config ${bundle_config} --catalog ${catalog} --workers ${task.cpus - 1} --output tables
    """

    stub:
    """
    mkdir -p tables/variants tables/transcripts
    echo '{}' > tables/annotation.json
    touch -t 200001010000 tables
    """
}

// The sv and sv_transcripts tables of the union of a group's structural variants, from its Annotator JSON, checked
// against its sites. One process reads the JSON, which holds far fewer positions than the union of the alleles.
// The stub gives its folder a fixed time, as the stub of ECHTVAR_VALUES does.
process SV_ANNOTATION_TABLES {
    tag "${group}"
    label 'tools'
    machineType 'n2d-standard-4'
    cpus 4
    memory 16.GB
    disk 50.GB

    input:
    tuple val(group), path(sites), path(sites_index), path(annotation, stageAs: 'annotation/*'), path(run_config, stageAs: 'run.json')
    path bundle_config, stageAs: 'bundle.json'
    val catalog

    output:
    path 'sv_annotation'

    script:
    """
    case-annotation sv-annotation-tables --json ${annotation} --sites ${sites} --run-config ${run_config} \\
        --bundle-config ${bundle_config} --catalog ${catalog} --output sv_annotation
    """

    stub:
    """
    mkdir -p sv_annotation/sv sv_annotation/sv_transcripts
    echo '{}' > sv_annotation/annotation.json
    touch -t 200001010000 sv_annotation
    """
}

// A family's runs of homozygosity and regions of uniparental disomy, from two published tools that run as they are
// (docs/annotation.md): `bcftools roh` for every sequenced member, and `upd regions` for each child whose parents are
// both sequenced. case-annotation first reads the pedigree for those children and the archive of gnomAD's joint release
// for the INFO field of an allele's frequency, and refuses a family with such a child whose VCF has no GQ. Both tools
// read one prepared VCF, cut from the family's VCF as the sheet names it: the PASS sites on chr1 to chr22, each
// alternate allele that no member carries removed, then the sites with one alternate allele of one base, without INFO
// and with every genotype unphased, and with gnomAD's frequency of each allele from the archive. The alleles are
// removed in a call of their own, since one call of bcftools view counts a site's alleles before it removes any. upd
// stops at a site of several alleles, at an INFO header line with a key beside ID, Number, Type and Description and at
// a frequency of `.`, and reads a phased genotype as no call. X, Y and M are left out: a haploid call is no genotype
// that upd reads.
// upd keeps the sites with a frequency of at least 0.05 and a GQ of at least 30 in every sample, which are its own
// thresholds and are given to it from the array that the question below takes them from. It fails with a TypeError on
// a VCF of which it keeps no site, as a call set of a few regions can be, and as one without a site on chr1 to chr22
// is. case-annotation therefore reads the prepared VCF by upd's rule first, and upd is not run where it would keep
// none: the family's children are then recorded as not assessed, and the family gets its tables. A tool that fails for
// any other reason fails the task, and its output is then never read.
// case-annotation writes the tools' output as tables, with a record of the commands that ran, which the step passes
// from the arrays it ran them with. The archive is copied to the task's disk, which archivesDisk sizes. The prepared
// VCF gets an index only because htslib reports a file it opens without one.
process ROH_UPD {
    tag "${family}"
    label 'tools'
    disk { archivesDisk([gnomad_joint]) }

    input:
    tuple val(family), path(vcf, stageAs: 'input/*'), path(ped, stageAs: 'input/*')
    path gnomad_joint, stageAs: 'gnomad_joint/*'

    output:
    tuple val(family), path('roh_upd')

    script:
    """
    set -o pipefail
    tag=\$(case-annotation roh-upd-plan --family ${family} --ped ${ped} --vcf ${vcf} --archive ${gnomad_joint} --trios trios.tsv)
    autosomes=\$(printf 'chr%s,' {1..22})
    kept=(view -f PASS,. --trim-alt-alleles -t "\${autosomes%,}")
    single=(view -m2 -M2 -v snps)
    bare=(annotate -x INFO)
    unphased=(-t a -n u)
    bcftools "\${kept[@]}" --no-version -Ou ${vcf} \\
        | bcftools "\${single[@]}" --no-version -Ou \\
        | bcftools "\${bare[@]}" --no-version -Ou \\
        | bcftools +setGT --no-version -Ob -o prepared.bcf -- "\${unphased[@]}"
    sites=\$(bcftools view --no-version -H prepared.bcf | wc -l)
    echo "\${sites} sites for bcftools roh and upd"
    echtvar_version=\$(echtvar --version)
    echtvar anno -e ${gnomad_joint} prepared.bcf prepared.vcf.gz
    bcftools index -t prepared.vcf.gz
    roh=(roh --AF-tag "\${tag}" -G30 -O r)
    bcftools "\${roh[@]}" -o roh.txt prepared.vcf.gz
    thresholds=(--af-tag "\${tag}" --min-af 0.05 --min-gq 30)
    called=("\${thresholds[@]}" regions)
    sites_of_upd=\$(case-annotation roh-upd-sites --prepared prepared.vcf.gz "\${thresholds[@]}")
    mkdir upd
    regions=()
    while IFS=\$'\\t' read -r child mother father; do
        if [[ \${sites_of_upd} == none ]]; then
            echo "${family}: upd would keep no site, so it is not run for \${child}"
            regions+=(--not-assessed "\${child}")
            continue
        fi
        upd --vcf prepared.vcf.gz --proband "\${child}" --mother "\${mother}" --father "\${father}" "\${called[@]}" --out "upd/\${child}.bed"
        regions+=(--regions "\${child}" "upd/\${child}.bed")
    done < trios.tsv
    case-annotation roh-upd --family ${family} --ped ${ped} --archive ${gnomad_joint} --prepared prepared.vcf.gz \\
        --sites "\${sites}" --roh roh.txt "\${regions[@]}" \\
        --bcftools-version "\$(bcftools --version | awk 'NR == 1 { print \$2 }')" \\
        --echtvar-version "\${echtvar_version#echtvar }" \\
        --preparation "bcftools \${kept[*]}" --preparation "bcftools \${single[*]}" \\
        --preparation "bcftools \${bare[*]}" --preparation "bcftools +setGT -- \${unphased[*]}" \\
        --preparation "echtvar anno" \\
        --roh-command "bcftools \${roh[*]}" --upd-command "upd \${called[*]}" --output roh_upd
    """

    stub:
    """
    mkdir roh_upd
    touch roh_upd/roh.parquet roh_upd/roh_upd.json
    """
}

// One family's variants table, which holds a column of each member's calls, its transcripts table, its manifest and
// the notes on the tables (tables.md, a copy of the file in the tools' package), from the union's tables, its values
// from the echtvar archives and from the joined files and its table of VEP's plugins. The manifest's roles are those of
// the proband the sheet names
// and of its parents; nobody has one where the sheet names none. DuckDB sorts the union's tables down to the family's
// alleles, spilling to the boot disk.
// A family with a structural-variant or a copy-number VCF also gets sv, which holds a column of each member's calls
// too, and sv_transcripts, and one with repeat VCFs repeats; sv_alleles is NORMALISE_SV's output for a family with a
// structural-variant VCF. sv_annotation and sv_map are the tables and the map of the union of the group's structural
// variants, which a family with either VCF reads its records' annotation from, and empty lists in a group without
// such a family. strchive is the file of the loci of one release of STRchive, with the release's name and date as the
// run states them, and an empty list and empty strings in a run without one: case-annotation checks the file for every
// family and adds its size ranges to the repeats table of a family that has one.
// Every family also gets roh, its members' runs of homozygosity, and one with a child whose parents are both sequenced
// upd, such a child's regions of uniparental disomy, both from ROH_UPD's folder of the family (segments). The stub
// cannot read a pedigree, so it writes no upd.
process FAMILY_TABLES {
    tag "${family}"
    label 'tools'
    machineType 'n2d-highmem-8'
    cpus 8
    memory 64.GB
    disk 200.GB

    input:
    tuple val(family), path(vcf), path(vcf_index), path(ped, stageAs: 'input/*'), val(proband), path(sv, stageAs: 'sv/*'), path(cnv, stageAs: 'cnv/*'), path(repeats, stageAs: 'repeats/*'), path(sv_alleles, stageAs: 'sv-alleles/*'), path(segments)
    path tables
    path echtvar_values
    path joined_values
    path vep_plugins
    path sv_annotation, stageAs: 'sv-union/*'
    path sv_map, stageAs: 'sv-union/*'
    path strchive, stageAs: 'strchive/*'
    val strchive_release
    val strchive_release_date
    val schema_version
    val pipeline_version

    output:
    tuple val(family), path('variants.parquet'), path('transcripts.parquet'), path('roh.parquet'), path('manifest.json'), path('tables.md'), path('{sv,sv_transcripts,repeats,upd}.parquet', arity: '0..4')

    script:
    // Only a family with a caller's VCF reads the annotation of the group's structural variants.
    def annotated = vcfs(sv) + vcfs(cnv) ? "--sv-annotation ${sv_annotation} --sv-map ${sv_map}" : ''
    def ranges = strchive ? "--strchive ${strchive} --strchive-release ${strchive_release} --strchive-release-date ${strchive_release_date}" : ''
    """
    case-annotation family-tables --family ${family} --ped ${ped} --vcf ${vcf} ${probandOption(proband)} \\
        ${structuralOptions(sv, cnv, repeats)} ${vcfs(sv_alleles).collect { f -> "--sv-alleles ${f}" }.join(' ')} \\
        ${annotated} ${ranges} \\
        --annotation ${tables} --echtvar-values ${echtvar_values} --joined-values ${joined_values} \\
        --vep-plugins ${vep_plugins} --roh-upd ${segments} \\
        --schema-version ${schema_version} --pipeline-version ${pipeline_version} --output out --scratch /tmp
    mv out/* .
    """

    stub:
    def larger = (vcfs(sv) + vcfs(cnv) ? ['sv.parquet', 'sv_transcripts.parquet'] : []) + (vcfs(repeats) ? ['repeats.parquet'] : [])
    """
    touch variants.parquet transcripts.parquet roh.parquet manifest.json tables.md ${larger.join(' ')}
    """
}

// How case-annotation takes a family's proband as the sheet states it: the sample id, or that none is designated where
// the sheet's value is empty.
def probandOption(proband) {
    return proband ? "--proband ${proband}" : '--no-proband'
}

// The VCFs among the files of an optional input, which holds each VCF with its index or nothing, in the order of
// their names.
def vcfs(files) {
    return [files].flatten().findAll { f -> !f.name.endsWith('.tbi') }.sort { f -> f.name }
}

// How case-annotation takes the VCFs of a family's larger variants: an option for each one the sheet names.
def structuralOptions(sv, cnv, repeats) {
    def options = vcfs(sv).collect { f -> "--sv-vcf ${f}" }
    options += vcfs(cnv).collect { f -> "--cnv-vcf ${f}" }
    options += vcfs(repeats).collect { f -> "--repeat-vcf ${f}" }
    return options.join(' ')
}

// The disk of ECHTVAR_VALUES: 50 GB for the image, the union and its annotated copy, and the size of every archive
// of the run, each of which is copied to the disk, rounded up to whole GB. An archive that the run does not have is an
// empty list among them.
def archivesDisk(archives) {
    def bytes = [archives].flatten().sum { archive -> archive.size() }
    return 1.GB * (50 + (bytes + 1024 * 1024 * 1024 - 1).intdiv(1024 * 1024 * 1024))
}

// The group's annotation among the files an Annotator run wrote, which is the one named after the group. A run that
// wrote none stops here with the group's name: left out silently, the group's tables would never be written and the
// workflow would still end well. The run over the group's structural variants is named svName(group).
def groupAnnotation(outputs, group) {
    def files = [outputs].flatten()
    def matching = files.findAll { f -> f.name == "${group}.json.gz" }
    if (matching.size() != 1) {
        error("${group}: the Annotator wrote no ${group}.json.gz for the group, but ${files.collect { f -> f.name }}")
    }
    return matching[0]
}

// The name of the Annotator run over the union of a group's structural variants, which names its JSON and its log
// beside those of the run over the group's alleles.
def svName(group) {
    return "${group}.sv"
}

// What SV_UNION takes of the families' sites: each family's id with the name of each file of its sites, and the files,
// both in the order of the families and the names, so that the task's command is the same in every run.
def svGroup(families) {
    def pairs = families
        .collectMany { family, files -> [files].flatten().collect { f -> [family, f] } }
        .sort { a, b -> a[0] <=> b[0] ?: a[1].name <=> b[1].name }
    return [families: pairs.collect { family, f -> [family, f.name] }, sites: pairs.collect { _family, f -> f }]
}
