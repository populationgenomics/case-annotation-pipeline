// Runs the family tables' processes on two invented families, wired as family_tables.nf wires them, and checks what
// they write: FAMILY_CHECK, NORMALISE_FAMILY and UNION_SITES; ANNOTATION_TABLES on an invented annotation of the union,
// which INVENTED_ANNOTATION writes in the Annotator's layout, with bundle.json and run.json standing in for the
// configs; ECHTVAR_VALUES on invented echtvar archives, which INVENTED_ARCHIVE and INVENTED_ALPHAGENOME_ARCHIVE
// encode; and FAMILY_TABLES, which reads
// each family's PED file and proband as family_tables.nf passes them. The sheet's proband is the daughter in family A
// and empty in family B, whose affected son therefore has no role. Everything under family_tables/ is invented, the PED
// files too: family-a.vcf holds a record of two alleles that both members carry, one of two alleles of which one nobody
// carries, a spanning deletion beside an allele, a deletion written at the right end of a repeat, a record nobody
// carries, an allele that the mother has no call of, a record that only a call of ./1 carries, a deletion of six bases,
// a record off the primary chromosomes and a record of two alleles on X with genotype likelihoods, called in a haploid
// father and a diploid mother and daughter; family-b.vcf holds an allele family A shares, the same deletion
// left-aligned, and one of its own. gnomad-joint.vcf stands in for gnomAD's joint sites VCF, with invented counts and
// frequencies under gnomAD's field names, and echtvar-config.json encodes it as the public archive of
// docs/annotation.md is encoded: it holds one of the two alleles at a site, an allele counted 0 times, the deletion at
// its left-aligned position, the deletion of six bases (echtvar keeps alleles of more than four bases in a table of
// their own), an allele whose frequency is a usual threshold, a filtering allele frequency that is a multiple of
// 1/2000000 and that echtvar's 32-bit arithmetic stores one step below, and another allele at a position of the
// union's. reference.fa holds an N at chr2:300 and the IUPAC code R at chr2:320, which bcftools reads as N: ref-n.vcf,
// whose REF is N at both, is normalised, and a REF of N at a base (ref-n-at-base.vcf) and a base as REF at the N and at
// the R (ref-base-at-n.vcf, ref-base-at-iupac.vcf) are refused. Family A also brings the VCFs of its larger variants,
// which COMPRESSED gives the indexes a run's sheet names them with, and family B brings none. family-a.sv.vcf holds a
// deletion written as a sequence at another place of a repeat than family-a.vcf writes it, the deletion of six bases as
// family-a.vcf has it, a deletion that no small variant is, a symbolic deletion, a breakend and its mate, a record off
// the primary chromosomes and a symbolic deletion on X, where the father is called haploid. family-a.cnv.vcf holds a
// segment at the reference's copy number, a loss, a segment that the father loses and the mother gains, and a gain on
// X. A1, A2 and A3.repeats.vcf hold each member's genotypes of two invented repeats. strchive.json stands in for the
// loci of a release of STRchive, with three invented loci under STRchive's field names: one that the repeat on chr1 is
// of, with a wider interval and another rotation of its unit, one that overlaps the repeat on X with another motif, and
// one that no repeat lies in. The checks refuse repeat VCFs of
// two of the three members, a structural-variant VCF of another family's samples, and a structural variant whose REF is
// not the reference's (sv-ref-mismatch.vcf).
// The structural variants and copy-number segments go the way family_tables.nf sends them: SV_SITES cuts family A's
// sites, SV_UNION writes their union and its map, INVENTED_SV_ANNOTATION writes an invented annotation of the union in
// the Annotator's layout, SV_ANNOTATION_TABLES converts it, and FAMILY_TABLES reads family A's records' annotation from
// it. A deletion written as a sequence without INFO END (sv-no-end.vcf) is refused when the union is written.
// The ##source line of the VCFs of the larger variants names build/tests/family.nf, which means this file.
// alphagenome_avi.tsv, alphagenome_avi_breakdown.tsv and alphagenome_splicing.tsv stand in for DeepMind's three tables
// of AlphaGenome's scores (docs/annotation.md, AlphaGenome's archives), with invented values under DeepMind's column
// names. INVENTED_ALPHAGENOME_ARCHIVE builds an archive of each as the document says, with the conversion and the
// encode configuration of the tools, and checks that the archive returns the table's rows. The tables hold the three
// single-base changes at each of a few positions of the union: every change at one position has a score, a breakdown
// and a splicing score; other positions lack the breakdown, the splicing score or both; the position of a deletion
// has scores of its single-base changes; and other positions of the union have no row. One raw score is 0, another is
// negative, a PHRED is 0, one is above 67, where a 32-bit float is coarser than the millionths that the archive
// stores, and one contribution is in scientific notation. ECHTVAR_VALUES reads the three archives with gnomAD's in one
// echtvar call, and a second time gnomAD's alone, as in a run that names no archive of the scores.
// gencode_polya.gtf and aparent2.tsv stand in for GENCODE's polyadenylation features and for APARENT2's table of
// effects, with invented features and values in the layout of the public files. JOINED_VALUES joins them to the union:
// three alleles change a base of a signal, one of them a deletion whose last base lies in it, and a change in a site
// and one in a pseudogene's signal overlap none; one allele has a row for each of two sites, and others one row or
// none. genomic_constraint.tsv, mcr.tsv and mpc.tsv stand in for gnomAD's tables of genomic constraint by window, of
// regions of missense constraint and of MPC scores, with invented windows, regions and scores under gnomAD's column
// names and the invented transcripts of the annotation. The union's alleles on chr1 and chr2 lie in a window each and
// those on chrX in none; ENST01, the transcript of every row of the annotation, has two regions on chr1, and ENST02,
// which the annotation lacks, one on chr2 and a score, so its values are in variants and on no row of transcripts.
// branch_points.tsv stands in for the prepared catalogue of branch points, with invented branch points under the
// catalogue's column names and two lines before its header: a single-base change at one branch point, a deletion that
// ends at one of the reverse strand, and a deletion over two branch points that lie two bases apart, of which the first
// is in two introns and the second in none.
// The run names no file of gnomAD's mitochondrial sites and none of MitoTIP's scores, so the families' tables lack
// those columns. gnomad_mito.vcf and mitotip.tsv stand in for those files, with invented values under the field names
// of gnomAD's file and the column names of MitoTIP's, and are joined to mito-union.vcf, an invented union on an
// invented chrM of 120 bases (mito-reference.fa): gnomad_mito.vcf writes a deletion at the right end of a repeat,
// which the union holds at the left end, as a family's normalised VCF would, and mitotip.tsv scores the changes of one
// base of the union and the deletion of a base that the union's deletion removes.
// The step that applies the rules of VEP's two plugins goes the way family_tables.nf sends it too. gene-models/ holds
// one invented gene of four exons on chr1, on the forward strand, in the layout of Ensembl's GFF3 and in that of
// NCBI's RefSeq GFF, and an invented assembly report. Its 5' UTR is the first exon and the start of the second, and
// holds one start codon, at chr1:82, whose reading frame runs into the coding sequence in frame. INVENTED_GENE_MODELS
// prepares the folder that a run is given, and VEP_PLUGIN_TABLE cuts the union to the gene's exons, adds the positive
// controls and writes the table of the plugins' rules for both files. Of the union's alleles, the two at chr1:100 lie
// in the 5' UTR without an effect, chr1-200-G-A makes a stop codon in the upstream reading frame, chr1-300-T-A is a
// missense variant, which neither plugin's rules evaluate, and chr1-400-TCA-T is a frameshift in the third exon, which
// the NMD plugin's rule evaluates and does not flag. A folder of gene models whose record names the controls on
// another transcript is refused, which is what a run meets whose gene models are not as their record describes them.
// Sites that hold an allele of the exons twice are refused before the rules are applied to them.
// Refused inputs must fail, and so must Annotator outputs without the group's annotation. The case-annotation and
// echtvar steps run in the tools image, built first. Run with Docker from the repository's root, with this folder's
// config alone (-C), because a run from the root would otherwise read the pipeline's nextflow.config over it:
//   docker build -f containers/tools/Dockerfile -t case-annotation-tools-test .
//   nextflow -C tests/workflows/nextflow.config run tests/workflows/family_tables.nf

include { FAMILY_CHECK ; NORMALISE_FAMILY ; NORMALISE_SV ; UNION_SITES ; SV_SITES ; SV_UNION ; ECHTVAR_VALUES ; JOINED_VALUES ; VEP_PLUGIN_TABLE ; ANNOTATION_TABLES ; SV_ANNOTATION_TABLES ; ROH_UPD ; FAMILY_TABLES ; groupAnnotation ; svName ; svGroup ; archivesDisk } from '../../modules/family_tables'
include { ROH_UPD as ROH_UPD_TRIO } from '../../modules/family_tables'
include { ROH_UPD as ROH_UPD_DUO } from '../../modules/family_tables'
include { ROH_UPD as ROH_UPD_NO_GQ_REFUSED } from '../../modules/family_tables'
include { ROH_UPD as ROH_UPD_NO_SITE } from '../../modules/family_tables'
include { ROH_UPD as ROH_UPD_NO_KEPT_SITE } from '../../modules/family_tables'
include { ROH_UPD as ROH_UPD_FAILING_UPD_REFUSED } from '../../modules/family_tables'
include { ECHTVAR_VALUES as ECHTVAR_VALUES_GNOMAD_ALONE } from '../../modules/family_tables'
include { JOINED_VALUES as JOINED_VALUES_MITO } from '../../modules/family_tables'
include { VEP_PLUGIN_TABLE as VEP_PLUGIN_TABLE_CONTROL_REFUSED } from '../../modules/family_tables'
include { VEP_PLUGIN_TABLE as VEP_PLUGIN_TABLE_REPEATED_REFUSED } from '../../modules/family_tables'
include { SV_SITES as SV_SITES_NO_END } from '../../modules/family_tables'
include { SV_UNION as SV_UNION_NO_END_REFUSED } from '../../modules/family_tables'
include { NORMALISE_SV as NORMALISE_SV_REF_MISMATCH_REFUSED } from '../../modules/family_tables'
include { FAMILY_CHECK as FAMILY_CHECK_REPEATS_REFUSED } from '../../modules/family_tables'
include { FAMILY_CHECK as FAMILY_CHECK_SV_SAMPLES_REFUSED } from '../../modules/family_tables'
include { NORMALISE_FAMILY as NORMALISE_UNPREFIXED_REFUSED } from '../../modules/family_tables'
include { NORMALISE_FAMILY as NORMALISE_REF_MISMATCH_REFUSED } from '../../modules/family_tables'
include { NORMALISE_FAMILY as NORMALISE_REPEATED_REFUSED } from '../../modules/family_tables'
include { NORMALISE_FAMILY as NORMALISE_REF_N } from '../../modules/family_tables'
include { NORMALISE_FAMILY as NORMALISE_REF_N_AT_BASE_REFUSED } from '../../modules/family_tables'
include { NORMALISE_FAMILY as NORMALISE_REF_BASE_AT_N_REFUSED } from '../../modules/family_tables'
include { NORMALISE_FAMILY as NORMALISE_REF_BASE_AT_IUPAC_REFUSED } from '../../modules/family_tables'
include { FAMILY_CHECK as FAMILY_CHECK_PROBAND_REFUSED } from '../../modules/family_tables'

workflow {
    def fixture = "${projectDir}/family_tables"
    reference = channel.value(tuple(file("${fixture}/reference.fa"), file("${fixture}/reference.fa.fai")))
    // The sheet's rows: family A names its proband, and family B designates none.
    families = channel.of(['family-a', 'A3'], ['family-b', ''])
        .map { f, proband ->
            tuple(f, file("${fixture}/${f}.vcf"), file("${fixture}/${f}.ped"), proband)
        }

    // The VCFs of family A's larger variants, each with its index, as family_tables.nf stages what a sheet names: the
    // structural-variant caller's, the copy-number caller's, and the three members' repeat VCFs. Family B has none.
    // The VCF whose deletion has no END is compressed with them, for the union that is refused below.
    compressed = COMPRESSED(
        channel.of(['sv', 'family-a.sv.vcf'], ['cnv', 'family-a.cnv.vcf'], ['repeats', 'A1.repeats.vcf'], ['repeats', 'A2.repeats.vcf'], ['repeats', 'A3.repeats.vcf'], ['no-end', 'sv-no-end.vcf']).map { kind, name -> tuple(kind, file("${fixture}/${name}")) }
    )
    larger = compressed
        .filter { kind, _vcf, _index -> kind != 'no-end' }
        .map { kind, vcf, index -> tuple(kind, [vcf, index]) }
        .groupTuple()
        .toList()
        .map { kinds ->
            def files = kinds.collectEntries { kind, pairs -> [kind, pairs.flatten()] }
            tuple('family-a', files.sv, files.cnv, files.repeats)
        }
    structural = larger.mix(channel.of(tuple('family-b', [], [], [])))

    checked = FAMILY_CHECK(families.join(structural), '0.5.0')
    normalised = NORMALISE_FAMILY(
        families.join(checked.map { f -> [f] }).map { f, vcf, _ped, _proband -> tuple(f, vcf) },
        reference,
    )
    sv_alleles = NORMALISE_SV(
        structural.join(checked.map { f -> [f] }).filter { _f, sv, _cnv, _repeats -> sv }.map { f, sv, _cnv, _repeats -> tuple(f, sv) },
        reference,
    )
    union = UNION_SITES('union', normalised.sites.flatten().collect())
    CHECK(
        normalised.family.filter { f -> f[0] == 'family-a' }.map { _f, vcf, index -> tuple(vcf, index) },
        normalised.sites.filter { s -> s[0].name == 'family-a.sites.vcf.gz' },
        union.map { _group, vcf, index -> tuple(vcf, index) },
        file("${fixture}/expected.family-a.txt"),
        file("${fixture}/expected.union.txt"),
    )
    annotation = INVENTED_ANNOTATION(union)
    tables = ANNOTATION_TABLES(
        union.combine(annotation).combine(channel.of(file("${fixture}/run.json"))),
        file("${fixture}/bundle.json"),
        '209901',
    ).first()
    // An archive of each of the two invented VCFs of gnomAD's joint values: the one the families' alleles are read
    // from, and the frequencies of the sites of the trio that ROH_UPD is tested on.
    archives = INVENTED_ARCHIVE(
        channel.of(['invented', 'gnomad-joint.vcf'], ['frequencies', 'roh-upd-gnomad.vcf']).map { name, vcf -> tuple(name, file("${fixture}/${vcf}")) },
        file("${fixture}/echtvar-config.json"),
    )
    archive = archives.filter { name, _zip -> name == 'invented' }.map { _name, zip -> zip }.first()
    frequencies = archives.filter { name, _zip -> name == 'frequencies' }.map { _name, zip -> zip }.first()
    // An archive of each of the three invented tables of AlphaGenome's scores, built and checked.
    scores = INVENTED_ALPHAGENOME_ARCHIVE(
        channel.of('alphagenome_avi', 'alphagenome_avi_breakdown', 'alphagenome_splicing').map { table -> tuple(table, file("${fixture}/${table}.tsv")) },
        file("${fixture}/reference.fa.fai"),
    )
    avi = scores.filter { table, _zip -> table == 'alphagenome_avi' }.map { _table, zip -> zip }.first()
    breakdown = scores.filter { table, _zip -> table == 'alphagenome_avi_breakdown' }.map { _table, zip -> zip }.first()
    splicing = scores.filter { table, _zip -> table == 'alphagenome_splicing' }.map { _table, zip -> zip }.first()
    joint = ECHTVAR_VALUES(union, archive, avi, breakdown, splicing)
    // A run that names no archive of the scores gives the step an empty list for each.
    CHECK_GNOMAD_ALONE(ECHTVAR_VALUES_GNOMAD_ALONE(union, archive, [], [], []))
    // The step's disk is 50 GB and the size of its archives, rounded up to whole GB, and an archive that the run does
    // not have adds nothing.
    def small = file("${fixture}/reference.fa")
    if (archivesDisk([small, [], [], []]) != 51.GB || archivesDisk([small, small, small, small]) != 51.GB) {
        error("the disk of four small archives is ${archivesDisk([small, small, small, small])}, and of one ${archivesDisk([small, [], [], []])}")
    }
    // The two files on polyadenylation, gnomAD's three tables of constraint and the catalogue of branch points,
    // joined to the union, and no file of gnomAD's mitochondrial sites or of MitoTIP's scores.
    polya = file("${fixture}/gencode_polya.gtf")
    effects = file("${fixture}/aparent2.tsv")
    windows = file("${fixture}/genomic_constraint.tsv")
    regions = file("${fixture}/mcr.tsv")
    scores = file("${fixture}/mpc.tsv")
    catalogue = file("${fixture}/branch_points.tsv")
    joined = JOINED_VALUES(union, reference, [], [], polya, effects, windows, regions, scores, catalogue)
    // gnomAD's mitochondrial sites and MitoTIP's scores, joined to an invented union on an invented chrM.
    mito = file("${fixture}/gnomad_mito.vcf")
    mitotip = file("${fixture}/mitotip.tsv")
    CHECK_JOINED_MITO(
        JOINED_VALUES_MITO(
            INVENTED_MITO_UNION(file("${fixture}/mito-union.vcf")),
            channel.value(tuple(file("${fixture}/mito-reference.fa"), file("${fixture}/mito-reference.fa.fai"))),
            mito,
            mitotip,
            [],
            [],
            [],
            [],
            [],
            [],
        ),
        mito,
        mitotip,
    )
    // The rules of VEP's two plugins, applied to the union's alleles in the exons of the invented gene.
    gene_models = INVENTED_GENE_MODELS(file("${fixture}/gene-models"), reference)
    vep_plugins = VEP_PLUGIN_TABLE(union, reference, gene_models)
    CHECK_VEP_PLUGINS(vep_plugins)
    // The union of the group's structural variants and copy-number segments, which are family A's, and its tables.
    sv_sites = SV_SITES(
        structural.join(checked.map { f -> [f] }).filter { _f, sv, cnv, _repeats -> sv || cnv }.map { f, sv, cnv, _repeats -> tuple(f, sv, cnv) }
    )
    sv_group = sv_sites.toList().map { cut -> svGroup(cut) }
    sv_union = SV_UNION('union', sv_group.map { g -> g.families }, sv_group.map { g -> g.sites })
    CHECK_SV_UNION(sv_union.union.map { _group, vcf, index -> tuple(vcf, index) }, sv_union.map)
    sv_annotation = INVENTED_SV_ANNOTATION(sv_union.union.map { _group, vcf, _index -> tuple(svName('union'), vcf) })
    sv_tables = SV_ANNOTATION_TABLES(
        sv_union.union.combine(sv_annotation).combine(channel.of(file("${fixture}/run.json"))),
        file("${fixture}/bundle.json"),
        '209901',
    ).first()
    // Each family's runs of homozygosity and regions of uniparental disomy, from its VCF as the sheet names it.
    segments = ROH_UPD(
        families.join(checked.map { f -> [f] }).map { f, vcf, ped, _proband -> tuple(f, vcf, ped) },
        archive,
    )
    written = FAMILY_TABLES(
        normalised.family.join(families.map { f, _vcf, ped, proband -> [f, ped, proband] }).join(structural).join(sv_alleles, remainder: true).map { f, vcf, index, ped, proband, sv, cnv, repeats, alleles ->
            tuple(f, vcf, index, ped, proband, sv, cnv, repeats, alleles ?: [])
        }.join(segments),
        tables,
        joint,
        joined,
        vep_plugins,
        sv_tables,
        sv_union.map,
        file("${fixture}/strchive.json"),
        'v0.0.0',
        '2099-01-01',
        '0.5.0',
        '0.5.0',
    )
    // What FAMILY_TABLES gives of a family: its id, variants, transcripts, roh, the manifest, the notes, and the
    // tables that only some families have.
    CHECK_TABLES(written.filter { w -> w[0] == 'family-a' }.map { w -> tuple(w[0], w[1], w[2], w[4], w[5]) }, archive)
    CHECK_NO_ROLES(written.filter { w -> w[0] == 'family-b' }.map { w -> w[4] })
    CHECK_JOINED(written.filter { w -> w[0] == 'family-a' }.map { w -> tuple(w[1], w[2], w[4]) }, polya, effects, windows, regions, scores, catalogue)
    CHECK_LARGER(written.filter { w -> w[0] == 'family-a' }.map { w -> tuple(w[1], w[4], w[6]) })
    CHECK_SEGMENTS(written.filter { w -> w[0] == 'family-a' }.map { w -> tuple(w[3], w[4], w[6]) }, archive)
    // A family whose row names no such VCF gets no table of them. Family B's son has both parents sequenced, so it
    // gets the table of his regions of uniparental disomy.
    written
        .filter { w -> w[0] == 'family-b' }
        .map { w ->
            def others = [w[6]].flatten().collect { f -> f.name }
            if (others != ['upd.parquet']) {
                error("family B brought no VCF of its larger variants and has a child with both parents, and has the tables ${others}")
            }
        }

    // ROH_UPD on an invented trio whose child has a region of each kind, and on the mother and the child alone.
    trio = file("${fixture}/roh-upd-trio.vcf")
    trio_ped = file("${fixture}/roh-upd-trio.ped")
    derived = DERIVED(
        channel.of(['duo', 'view -s T2,T3'], ['no-gq', 'annotate -x FORMAT/GQ'], ['only-x', 'view -t chrX']).map { name, command -> tuple(name, command, trio) }
    )
    CHECK_ROH_UPD(
        ROH_UPD_TRIO(channel.of(tuple('trio-t', trio, trio_ped)), frequencies).map { _f, folder -> folder },
        ROH_UPD_DUO(derived.filter { name, _vcf -> name == 'duo' }.map { _name, vcf -> tuple('trio-t', vcf, trio_ped) }, frequencies).map { _f, folder -> folder },
    )

    // A REF of N is normalised where the reference has an N or an IUPAC code.
    ref_n = NORMALISE_REF_N(channel.of(tuple('ref-n', file("${fixture}/ref-n.vcf"))), reference)
    CHECK_REF_N(ref_n.family.map { _f, vcf, index -> tuple(vcf, index) })

    // Where upd would keep no site it is not run, and the family gets its tables: the trio with the archive of the
    // families' alleles, which gives none of its sites a frequency that upd keeps, and the trio's record on X alone,
    // which leaves no site on chr1 to chr22.
    CHECK_NOT_ASSESSED(
        ROH_UPD_NO_KEPT_SITE(channel.of(tuple('trio-t', trio, trio_ped)), archive).map { _f, folder -> folder },
        ROH_UPD_NO_SITE(derived.filter { name, _vcf -> name == 'only-x' }.map { _name, vcf -> tuple('trio-t', vcf, trio_ped) }, frequencies).map { _f, folder -> folder },
    )

    // The group's annotation is the Annotator's output named after the group, alone or among others; outputs without
    // it are refused with the group's name.
    def outputs = [file('union.json.gz'), file('union.json.gz.jsi')]
    if (groupAnnotation(outputs, 'union') != outputs[0] || groupAnnotation(outputs[0], 'union') != outputs[0]) {
        error("the group's annotation is not the output named after the group")
    }
    def refusal = null
    try {
        groupAnnotation([file('other.json.gz'), file('union.json')], 'union')
    }
    catch (e: Exception) {
        // Nextflow wraps the error of an included function in the exception of the call.
        refusal = e.cause != null ? e.cause.message : e.message
    }
    if (refusal == null || !refusal.startsWith('union: the Annotator wrote no union.json.gz')) {
        error("outputs without the group's annotation were not refused by the group's name: ${refusal}")
    }

    // Refused: contigs not named chr1, chr2, ..., which leave no allele; a REF that is not the reference's; an allele
    // that two records hold; a REF of N where the reference has a base; a base as REF where the reference has an N, and
    // where it has an IUPAC code that the base fits; a proband that is not an individual of the family's PED file;
    // repeat VCFs of two of the family's three members; a structural-variant VCF whose samples are another family's; a
    // structural variant written as a sequence whose REF is not the reference's; a union of sites with a deletion that
    // has no INFO END, which the Annotator would annotate as a single base; for ROH_UPD, a trio whose VCF has no GQ,
    // and an upd that fails after it has written a region, which failing-upd stands in for (nextflow.config); a folder
    // of gene models whose record names the positive controls on a transcript that the gene models do not hold; sites
    // with an allele of the exons twice.
    unprefixed = NORMALISE_UNPREFIXED_REFUSED(channel.of(tuple('unprefixed', file("${fixture}/unprefixed.vcf"))), reference)
    mismatch = NORMALISE_REF_MISMATCH_REFUSED(channel.of(tuple('mismatch', file("${fixture}/ref-mismatch.vcf"))), reference)
    repeated = NORMALISE_REPEATED_REFUSED(channel.of(tuple('repeated', file("${fixture}/repeated.vcf"))), reference)
    n_at_base = NORMALISE_REF_N_AT_BASE_REFUSED(
        channel.of(tuple('n-at-base', file("${fixture}/ref-n-at-base.vcf"))),
        reference,
    )
    base_at_n = NORMALISE_REF_BASE_AT_N_REFUSED(
        channel.of(tuple('base-at-n', file("${fixture}/ref-base-at-n.vcf"))),
        reference,
    )
    base_at_iupac = NORMALISE_REF_BASE_AT_IUPAC_REFUSED(
        channel.of(tuple('base-at-iupac', file("${fixture}/ref-base-at-iupac.vcf"))),
        reference,
    )
    other_proband = FAMILY_CHECK_PROBAND_REFUSED(
        channel.of(tuple('family-a', file("${fixture}/family-a.vcf"), file("${fixture}/family-a.ped"), 'B3', [], [], [])),
        '0.5.0',
    )
    two_repeats = FAMILY_CHECK_REPEATS_REFUSED(
        channel.of(tuple('family-a', file("${fixture}/family-a.vcf"), file("${fixture}/family-a.ped"), 'A3', [], [], [file("${fixture}/A1.repeats.vcf"), file("${fixture}/A2.repeats.vcf")])),
        '0.5.0',
    )
    other_samples = FAMILY_CHECK_SV_SAMPLES_REFUSED(
        channel.of(tuple('family-b', file("${fixture}/family-b.vcf"), file("${fixture}/family-b.ped"), '', [file("${fixture}/family-a.sv.vcf")], [], [])),
        '0.5.0',
    )
    sv_mismatch = NORMALISE_SV_REF_MISMATCH_REFUSED(
        channel.of(tuple('sv-mismatch', [file("${fixture}/sv-ref-mismatch.vcf")])),
        reference,
    )
    no_end = SV_SITES_NO_END(
        compressed.filter { kind, _vcf, _index -> kind == 'no-end' }.map { _kind, vcf, index -> tuple('no-end', [vcf, index], []) }
    )
    no_end_union = SV_UNION_NO_END_REFUSED(
        'no-end',
        no_end.map { f, sites -> svGroup([[f, sites]]).families },
        no_end.map { f, sites -> svGroup([[f, sites]]).sites },
    )
    no_gq = ROH_UPD_NO_GQ_REFUSED(
        derived.filter { name, _vcf -> name == 'no-gq' }.map { _name, vcf -> tuple('trio-t', vcf, trio_ped) },
        frequencies,
    )
    failing_upd = ROH_UPD_FAILING_UPD_REFUSED(channel.of(tuple('trio-t', trio, trio_ped)), frequencies)
    other_control = VEP_PLUGIN_TABLE_CONTROL_REFUSED(union, reference, INVENTED_OTHER_CONTROLS(gene_models))
    repeated_sites = VEP_PLUGIN_TABLE_REPEATED_REFUSED(INVENTED_REPEATED_SITES(), reference, gene_models)
    unprefixed.sites
        .map { _s -> 'unprefixed contigs' }
        .mix(
            no_gq.map { _f, _folder -> 'a trio without GQ' },
            failing_upd.map { _f, _folder -> 'an upd that fails' },
            no_end_union.map.map { _m -> 'a deletion without END' },
            other_control.map { _t -> "controls on a transcript that the gene models do not hold" },
            repeated_sites.map { _t -> 'an allele twice among the sites in the exons' },
            two_repeats.map { _f -> "repeat VCFs of two of three members" },
            other_samples.map { _f -> "a structural-variant VCF of another family's samples" },
            sv_mismatch.map { _a -> "a structural variant whose REF is not the reference's" },
            mismatch.sites.map { _s -> 'a REF mismatch' },
            repeated.sites.map { _s -> 'a repeated allele' },
            n_at_base.sites.map { _s -> 'a REF of N at a base' },
            base_at_n.sites.map { _s -> 'a base as REF at an N' },
            base_at_iupac.sites.map { _s -> 'a base as REF at an IUPAC code' },
            other_proband.map { _f -> "another family's proband" },
        )
        .toList()
        .map { accepted ->
            if (accepted) {
                error("should have been refused, but were accepted: ${accepted}")
            }
        }
}

// A VCF made from another with one bcftools command.
process DERIVED {
    tag "${name}"

    input:
    tuple val(name), val(command), path(vcf)

    output:
    tuple val(name), path("${name}.vcf")

    script:
    """
    bcftools ${command} --no-version -o ${name}.vcf ${vcf}
    """
}

// The folder of gene models that a run is given, prepared from the invented gene models as from the public files,
// which are gzipped.
process INVENTED_GENE_MODELS {
    input:
    path invented, stageAs: 'invented'
    tuple path(reference), path(reference_index)

    output:
    path 'gene-models'

    script:
    """
    gzip -c invented/ensembl.gff3 > invented.gff3.gz
    gzip -c invented/refseq.gff > invented.gff.gz
    case-annotation gene-models --ensembl invented.gff3.gz --refseq invented.gff.gz \\
        --assembly-report invented/assembly-report.txt --reference ${reference} --output gene-models
    """
}

// A copy of the folder of gene models whose record names each control on a transcript that the invented gene models do
// not hold. Its other files are those of the prepared folder.
process INVENTED_OTHER_CONTROLS {
    input:
    path prepared, stageAs: 'prepared'

    output:
    path 'gene-models'

    script:
    """
    mkdir gene-models
    cp -L prepared/* gene-models/
    sed -e 's/"ENST01.1"/"ENST02.1"/' -e 's/"NM_01.1"/"NM_02.1"/' prepared/gene-models.json \\
        > gene-models/gene-models.json
    if cmp -s prepared/gene-models.json gene-models/gene-models.json; then
        echo "the record has neither transcript of the invented gene models" >&2
        exit 1
    fi
    """
}

// Sites that hold one allele in the first exon of the invented gene twice, as no union does.
process INVENTED_REPEATED_SITES {
    output:
    tuple val('repeated'), path('repeated.vcf.gz'), path('repeated.vcf.gz.tbi')

    script:
    """
    printf '##fileformat=VCFv4.2\\n##contig=<ID=chr1,length=600>\\n' > repeated.vcf
    printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' >> repeated.vcf
    printf 'chr1\\t100\\t.\\tG\\tA\\t.\\t.\\t.\\n' >> repeated.vcf
    printf 'chr1\\t100\\t.\\tG\\tA\\t.\\t.\\t.\\n' >> repeated.vcf
    bcftools view --no-version -Oz -o repeated.vcf.gz repeated.vcf
    bcftools index -t repeated.vcf.gz
    """
}

// The table of the plugins' rules for the union, and its record. The table has a row for an allele and transcript that
// a rule evaluated, for the transcript of each gene model file, and none for a control: the alleles in the 5' UTR have
// the counts of the reading frames that the UTR holds, which is one that runs into the coding sequence in frame; the
// allele that makes a stop codon in that frame has one record of a gained stop, with the plugin's values; the
// frameshift is evaluated by the NMD plugin's rule and not flagged; the missense variant has no row. The record holds
// the release of VEP that the step follows, the commits that the rules are restated from, the plugin's file of
// translated reading frames as the image holds it, and each gene model file with the invented file that it was
// prepared from, and no count of the union's.
process CHECK_VEP_PLUGINS {
    input:
    path vep_plugins

    script:
    """
    python3 - <<'EOF'
import json
import re

import pyarrow.parquet as pq

from case_annotation import vep_plugins

table = pq.read_table('${vep_plugins}/vep_plugins.parquet')
assert table.schema.equals(vep_plugins.TABLE), table.schema
rows = table.to_pylist()
keys = [(r['variant_id'], r['transcript_id'], r['source']) for r in rows]
evaluated = ['chr1-100-G-A', 'chr1-100-G-C', 'chr1-200-G-A', 'chr1-400-TCA-T']
assert keys == [(a, t, s) for a in evaluated for t, s in (('ENST01.1', 'Ensembl'), ('NM_01.1', 'RefSeq'))], keys
counts = ('utr5_existing_uorfs', 'utr5_existing_oorfs_out_of_frame', 'utr5_existing_oorfs_in_frame')
for row in rows:
    in_utr = row['variant_id'] != 'chr1-400-TCA-T'
    assert tuple(row[c] for c in counts) == ((0, 0, 1) if in_utr else (None, None, None)), row
    assert row['nmd_escaping'] is (None if in_utr else False), row
    if row['variant_id'] != 'chr1-200-G-A':
        assert (row['utr5_consequences'], row['utr5_uorfs']) == (None, None), row
        continue
    assert row['utr5_consequences'] == ['uSTOP_gained'], row
    (record,) = row['utr5_uorfs']
    assert {field: value for field, value in record.items() if value is not None} == {
        'consequence': 'uSTOP_gained', 'ref_type': 'InFrame_oORF', 'kozak_context': 'GTGATGC',
        'kozak_strength': 'Moderate', 'ref_start_distance_to_cds': 117, 'new_stop_distance_to_cds': 12,
        'evidence': False}, record
info = json.load(open('${vep_plugins}/vep_plugins.json'))
assert info['vep_release'] == '116.3', info
commits = (info['utrannotator_commit'], info['nmd_plugin_commit'])
assert all(re.fullmatch('[0-9a-f]{40}', commit) for commit in commits), info
assert (info['uorf_evidence']['name'], info['uorf_evidence']['bytes']) == ('uORF_5UTR_GRCh38_PUBLIC.txt', 2376724), info
assert info['uorf_evidence']['sha256'] == 'e220023a747465dcdf6725ce809511bccc20b7a9990288da91032a1a7f51ebce', info
models = [(m['source'], m['file']['name'], m['prepared_from']['name']) for m in info['gene_models']]
assert models == [('Ensembl', 'ensembl.gff3.gz', 'invented.gff3.gz'),
                  ('RefSeq', 'refseq.gff.gz', 'invented.gff.gz')], models
assert sorted(info) == ['gene_models', 'nmd_plugin_commit', 'uorf_evidence', 'utrannotator_commit', 'vep_release'], info
stated = [m['prepared_from']['stated'] for m in info['gene_models']]
assert [len(lines) for lines in stated] == [2, 3] and all(line.startswith('#!') for lines in stated for line in lines), info
EOF
    """
}

// A VCF bgzipped, with its index beside it.
process COMPRESSED {
    input:
    tuple val(kind), path(vcf)

    output:
    tuple val(kind), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
    """
    bcftools view --no-version -Oz -o ${vcf}.gz ${vcf}
    bcftools index -t ${vcf}.gz
    """
}

// The union of family A's structural variants and copy-number segments, and its map: one record per allele, so the
// segment of two alleles has two; each under an id of 32 hexadecimal digits that the map gives the record's allele;
// no samples, QUAL or FILTER; every record but the two breakends with the END its caller wrote, also the deletions
// written as a sequence; and none of a breakend's MATEID, which names a record of the family's VCF.
process CHECK_SV_UNION {
    input:
    tuple path(union), path(union_index)
    path map

    script:
    """
    set -o pipefail
    printf 'chr1\\t200\\t<DEL>\\t400\\nchr1\\t300\\tT[chr2:200[\\t.\\nchr1\\t406\\tA\\t408\\nchr1\\t500\\t<DEL>\\t560\\nchr2\\t90\\tG\\t96\\nchr2\\t97\\tG\\t100\\nchr2\\t200\\t]chr1:300]C\\t.\\nchr2\\t300\\t<DEL>\\t450\\nchr2\\t300\\t<DUP>\\t450\\nchrX\\t200\\t<DEL>\\t260\\nchrX\\t300\\t<DUP>\\t500\\n' > expected.txt
    bcftools query -f '%CHROM\\t%POS\\t%ALT\\t%INFO/END\\n' ${union} | sort -k1,1 -k2,2n -k3,3 | diff expected.txt -
    test -z "\$(bcftools query -l ${union})"
    test "\$(bcftools query -f '%QUAL%FILTER\\n' ${union} | sort -u)" = '..'
    test -z "\$(bcftools view -H ${union} | grep -E 'MATEID|EVENT|REFLEN' || true)"
    test "\$(bcftools query -f '%INFO/SVLEN\\n' -i 'POS=300 && CHROM="chr2"' ${union} | sort | tr '\\n' ' ')" = '-150 150 '
    python3 - <<'EOF'
import re

import pyarrow.parquet as pq
import pysam

ids = [record.id for record in pysam.VariantFile('${union}')]
assert len(ids) == len(set(ids)) == 11, ids
assert all(re.fullmatch('[0-9a-f]{32}', i) for i in ids), ids
rows = pq.read_table('${map}').to_pylist()
assert {r['family'] for r in rows} == {'family-a'}, rows
assert sorted(r['site_id'] for r in rows) == sorted(ids), rows
alleles = sorted((r['sv_id'], r['alt']) for r in rows if r['sv_id'] in ('inv:CNV:1', 'inv:DEL:1'))
assert alleles == [('inv:CNV:1', '<DEL>'), ('inv:CNV:1', '<DUP>'), ('inv:DEL:1', 'A')], alleles
EOF
    """
}

// Family A's alleles split, left-aligned, carried and without the spanning deletion, with their genotypes; its sites
// with no samples, ID, QUAL, FILTER or INFO; the union holding each allele of the two families once.
process CHECK {
    container 'quay.io/biocontainers/bcftools:1.24--h118bc1c_4@sha256:186c45de3059dbb1837006c27788284c89d17f8e693e7804624b6cf6ca5caca0'

    input:
    tuple path(vcf), path(vcf_index)
    tuple path(sites), path(sites_index)
    tuple path(union), path(union_index)
    path expected_family
    path expected_union

    script:
    """
    set -o pipefail
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ${vcf} | diff ${expected_family} -
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${union} | diff ${expected_union} -
    test -z "\$(bcftools query -l ${sites})"
    test "\$(bcftools query -f '%ID%QUAL%FILTER%INFO\\n' ${sites} | sort -u)" = '....'
    test -z "\$(bcftools view -h ${sites} | grep -E '^##(INFO|FORMAT)=')"
    """
}

// The alleles normalised from ref-n.vcf: its two records, each with N as its REF.
process CHECK_REF_N {
    input:
    tuple path(vcf), path(vcf_index)

    script:
    """
    set -o pipefail
    printf 'chr2\\t300\\tN\\tA\\nchr2\\t320\\tN\\tA\\n' > expected.txt
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${vcf} | diff expected.txt -
    """
}

// The union's annotation in the Annotator's layout: a header without samples, one position line per allele with one
// invented transcript, and the genes section, which gives the transcript's invented gene an invented HGNC id.
process INVENTED_ANNOTATION {
    input:
    tuple val(group), path(union), path(union_index)

    output:
    path "${group}.json.gz"

    script:
    """
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${union} > alleles.tsv
    python3 - <<'EOF'
import gzip
import json

transcript = {'transcript': 'ENST01.1', 'source': 'Ensembl', 'hgnc': 'GENEA', 'geneId': 'ENSG01',
              'consequence': ['intron_variant']}
header = {'annotator': 'Invented Annotator', 'genomeAssembly': 'GRCh38', 'schemaVersion': 6, 'dataSources': []}
positions = []
for line in open('alleles.tsv'):
    chrom, pos, ref, alt = line.rstrip('\\n').split('\\t')
    positions.append(json.dumps({'chromosome': chrom, 'position': int(pos), 'refAllele': ref, 'altAlleles': [alt],
                                 'variants': [{'transcripts': [transcript]}]}))
with gzip.open('${group}.json.gz', 'wt') as f:
    f.write('{"header":' + json.dumps(header) + ',"positions":[\\n')
    f.write(',\\n'.join(positions) + '\\n')
    f.write('],"genes":[\\n{"name":"GENEA","ensemblGeneId":"ENSG01","hgncId":90001}\\n]}\\n')
EOF
    """
}

// The annotation of the union of the structural variants in the Annotator's layout: a header without samples, one
// position line per site under the site's id, and the genes section. Each site has an invented transcript. A sequence
// of at most 50 bases is annotated as a small variant, and the deletion of six bases has gnomAD's values of the allele
// itself. The symbolic deletion on chromosome 1 matched two entries of gnomAD's structural variants, the more frequent
// one at the lower overlap, and each allele of the segment of two alleles matched an entry of its own direction.
process INVENTED_SV_ANNOTATION {
    input:
    tuple val(name), path(union)

    output:
    path "${name}.json.gz"

    script:
    """
    python3 - <<'EOF'
import gzip
import json

from case_annotation import sv_sites

TYPES = {'DEL': 'deletion', 'INS': 'insertion', 'DUP': 'tandem_duplication', 'BND': 'translocation'}
COPIES = {'<DEL>': 'copy_number_loss', '<DUP>': 'copy_number_gain'}


def gnomad(variant_id, begin, end, kind, af, overlap):
    return {'chromosome': '1', 'begin': begin, 'end': end, 'variantId': variant_id, 'variantType': kind, 'allAf': af,
            'allAc': 5, 'allAn': 100, 'allHc': 0, 'reciprocalOverlap': overlap, 'annotationOverlap': overlap}


header = {'annotator': 'Invented Annotator', 'genomeAssembly': 'GRCh38', 'schemaVersion': 6, 'dataSources': []}
positions = []
for site_id, site in sv_sites.read_union('${union}').items():
    breakend = sv_sites.is_breakend(site.alt)
    structural = breakend or site.alt.startswith('<') or max(len(site.ref), len(site.alt)) > 50
    kind = COPIES[site.alt] if site.sv_type == 'CNV' else TYPES[site.sv_type]
    transcript = {'transcript': f'ENST-{kind}.1', 'source': 'Ensembl', 'hgnc': 'GENEA', 'geneId': 'ENSG01',
                  'consequence': ['intron_variant', 'feature_truncation'], 'impact': 'high'}
    variant = {'vid': 'x', 'chromosome': site.chrom, 'begin': site.pos if breakend else site.pos + 1,
               'end': site.pos if breakend else site.end, 'refAllele': site.ref, 'altAllele': site.alt,
               'variantType': kind, 'transcripts': [transcript]}
    position = {'chromosome': site.chrom, 'position': site.pos, 'id': site_id, 'refAllele': site.ref,
                'altAlleles': [site.alt], 'vcfInfo': {'SVTYPE': site.sv_type}, 'variants': [variant]}
    if structural:
        variant['isStructuralVariant'] = True
        if not breakend:
            position['svEnd'] = site.end
    if (site.chrom, site.pos, site.alt) == ('chr1', 500, '<DEL>'):
        position['gnomad'] = [gnomad('inv-common', 480, 570, 'deletion', 0.25, 0.65934),
                              gnomad('inv-best', 501, 560, 'deletion', 0.001, 1.0)]
    if (site.chrom, site.pos) == ('chr2', 300):
        position['gnomad'] = [gnomad(f'inv-{kind}', 301, 450, kind, 0.5, 1.0)]
    if (site.chrom, site.pos) == ('chr2', 90):
        variant['gnomad'] = {'allAf': 0.125, 'allAc': 1, 'allAn': 8, 'allHc': 0}
    positions.append(json.dumps(position))
with gzip.open('${name}.json.gz', 'wt') as f:
    f.write('{"header":' + json.dumps(header) + ',"positions":[\\n')
    f.write(',\\n'.join(positions) + '\\n')
    f.write('],"genes":[\\n{"name":"GENEA","ensemblGeneId":"ENSG01","hgncId":90001}\\n]}\\n')
EOF
    """
}

// An echtvar archive of the invented gnomAD VCF, with the fields of the public archive.
process INVENTED_ARCHIVE {
    tag "${name}"

    input:
    tuple val(name), path(vcf)
    path config

    output:
    tuple val(name), path("${name}.zip")

    script:
    """
    echtvar encode ${name}.zip ${config} ${vcf}
    """
}

// An echtvar archive of one of the invented tables of AlphaGenome's scores, built as docs/annotation.md says: the
// table's rows as a VCF, encoded with the table's configuration. The archive is then checked as the document says:
// echtvar annotates the rows' alleles with it, and the values that it returns are compared with the rows.
process INVENTED_ALPHAGENOME_ARCHIVE {
    tag "${table}"

    input:
    tuple val(table), path(tsv)
    path reference_index

    output:
    tuple val(table), path("${table}.zip")

    script:
    """
    set -o pipefail
    case-annotation alphagenome-config --table ${table} > config.json
    case-annotation alphagenome-vcf --table ${table} --reference-index ${reference_index} < ${tsv} > rows.vcf
    echtvar encode ${table}.zip config.json rows.vcf
    bcftools annotate --no-version -x INFO -Ob -o alleles.bcf rows.vcf
    echtvar anno -e ${table}.zip alleles.bcf annotated.bcf
    bcftools index annotated.bcf
    case-annotation alphagenome-check --table ${table} --annotated annotated.bcf < ${tsv}
    """
}

// Family A's tables: its eleven alleles, each with its invented transcript, its gene's HGNC id, the entry of that gene
// and, after the allele's columns of variants, a call column per member in the VCF's order, which says that the member
// carries the allele where its genotype holds a copy of it, that it does not where every allele is called and none is
// the allele, and neither for the mother's missing call; the manifest gives the daughter, whom the sheet names, the
// proband's role and her parents theirs.
// Each split allele of the record on X has each member's likelihoods of its own genotypes, two for the haploid father
// and three for the others, and an allele of a record without PL has none. An allele the archive holds has the counts
// it stores, and each frequency as the 32-bit float echtvar returns: the multiple of 1/2000000 at or below the invented
// value, or the one below that for 0.062502, and 0 for a filtering allele frequency the invented VCF lacks. An allele
// the archive lacks, such as the other allele of a site or of a position it holds, has none; the manifest records the
// archive by its size and SHA-256, with its fields and echtvar's version. The folder's notes on the tables are the
// file of the tools' package, byte for byte, and the manifest names it with its size and SHA-256. The manifest gives
// the pipeline's version as the installed package has it.
// An allele has AlphaGenome's scores from each archive that holds it, as 32-bit floats rounded to the decimals of the
// invented tables, which gives the tables' values back: the score and its PHRED from the first archive, the sixteen
// contributions from the second and the splicing score from the third. The contributions of an allele add up to its
// raw score and 0.048, as the invented tables have them. The two deletions and the alleles at positions without a row
// have none, and the alleles at chr1:200 and on chrX those of some archives. The manifest and the table's metadata
// record each archive, the two that are not for permissive use with their notice.
// Each row of transcripts has what the rules of VEP's plugins say of its allele on the invented transcript, which both
// the invented annotation and Ensembl's layout of the invented gene models call ENST01.1: the counts of the 5' UTR's
// reading frames for the three alleles in the UTR, the record of a gained stop for chr1-200-G-A, the false of the NMD
// plugin's rule for the frameshift, and nothing for the alleles that no rule evaluated. The plugins' six columns end
// with nmd_escaping, and the five columns that no plugin gives follow them and end the table. The allele's row of
// variants has the kind of the effect. The manifest counts the four rows of the rules' table on the RefSeq transcript,
// which the invented annotation lacks, and holds the record of the run.
process CHECK_TABLES {
    input:
    tuple val(family), path(variants), path(transcripts), path(manifest), path(notes)
    path archive

    script:
    """
    python3 - <<'EOF'
import hashlib
import importlib.metadata
import importlib.resources
import json
import struct

import pyarrow as pa
import pyarrow.parquet as pq


def float32(value):
    return struct.unpack('f', struct.pack('f', value))[0]


manifest = json.load(open('${manifest}'))
notes = open('${notes}', 'rb').read()
assert notes == importlib.resources.files('case_annotation').joinpath('tables.md').read_bytes()
assert manifest['notes'] == {
    'path': 'tables.md', 'bytes': len(notes), 'sha256': hashlib.sha256(notes).hexdigest()}, manifest
assert manifest['pipeline_version'] == importlib.metadata.version('case-annotation-tools'), manifest
rows = {t['table']: t['rows'] for t in manifest['tables']}
assert {t: rows[t] for t in ('variants', 'transcripts')} == {'variants': 11, 'transcripts': 11}, rows
assert manifest['family'] == '${family}', manifest
roles = {m['sample_id']: m['role'] for m in manifest['members']}
assert roles == {'A1': 'father', 'A2': 'mother', 'A3': 'proband'}, roles
variants = pq.read_table('${variants}')
assert variants.column_names[:5] == ['chrom', 'pos', 'ref', 'alt', 'variant_id'], variants.column_names
assert variants.column_names[-5:] == ['qual', 'filters', 'A1', 'A2', 'A3'], variants.column_names
assert variants.num_rows == 11, variants.num_rows
pl = {(r['variant_id'], s): (r[s]['gt'], r[s]['pl']) for r in variants.to_pylist() for s in ('A1', 'A2', 'A3')}
carries = {key: r[key[1]]['carries'] for r in variants.to_pylist() for key in pl if key[0] == r['variant_id']}
zygosity = {key: r[key[1]]['zygosity'] for r in variants.to_pylist() for key in pl if key[0] == r['variant_id']}


def carried(gt):
    alleles = gt.split('/')
    return True if '1' in alleles else None if '.' in alleles else False


assert all(carries[key] is carried(gt) for key, (gt, _) in pl.items()), carries
assert set(carries.values()) == {True, False, None}, carries
# The mother has no call of chr2-60-C-A, which is not a call of the reference.
assert (pl[('chr2-60-C-A', 'A2')][0], carries[('chr2-60-C-A', 'A2')], zygosity[('chr2-60-C-A', 'A2')]) == (
    './.', None, None), pl
# She has the other allele of the site at chr1:100, which the split writes as the reference in this allele's record.
assert (pl[('chr1-100-G-C', 'A2')][0], carries[('chr1-100-G-C', 'A2')], zygosity[('chr1-100-G-C', 'A2')]) == (
    '0/0', False, 'hom_ref'), pl
# A genotype with an allele missing has no zygosity, also the father's ./1 of chr2-80-C-T, which carries the allele.
assert all((zygosity[key] is None) == ('.' in gt) for key, (gt, _) in pl.items()), zygosity
assert (pl[('chr2-80-C-T', 'A1')][0], carries[('chr2-80-C-T', 'A1')], zygosity[('chr2-80-C-T', 'A1')]) == (
    './1', True, None), pl
assert all((zygosity[key] == 'hom_ref') == (carries[key] is False) for key in carries), zygosity
assert pl[('chrX-100-T-A', 'A1')] == ('0', [90, 80]), pl
assert pl[('chrX-100-T-A', 'A2')] == ('0/1', [70, 0, 60]), pl
assert pl[('chrX-100-T-A', 'A3')] == ('1/0', [73, 62, 63]), pl
assert pl[('chrX-100-T-C', 'A1')] == ('1', [90, 0]), pl
assert pl[('chrX-100-T-C', 'A2')] == ('0/0', [70, 71, 72]), pl
assert pl[('chrX-100-T-C', 'A3')] == ('0/1', [73, 51, 52]), pl
assert pl[('chr2-60-C-A', 'A1')] == ('0/1', None), pl
assert variants.schema.field('hgnc_ids').type == pa.list_(pa.int32()), variants.schema
assert variants.column('hgnc_ids').to_pylist() == [[90001]] * 11, variants.column('hgnc_ids')
genes = [[(g['gene_symbol'], g['ensembl_gene_id'], g['hgnc_id'], g['worst_consequence'], g['worst_transcript_id'],
           g['mane_select_transcript_id']) for g in entries]
         for entries in variants.column('genes').to_pylist()]
assert genes == [[('GENEA', 'ENSG01', 90001, 'intron_variant', 'ENST01.1', None)]] * 11, genes
assert not {'gene_symbol', 'ensembl_gene_id', 'worst_transcript_id'} & set(variants.column_names), variants.column_names
transcripts = pq.read_table('${transcripts}')
assert transcripts.schema.field('hgnc_id').type == pa.int32(), transcripts.schema
assert transcripts.column('hgnc_id').to_pylist() == [90001] * 11, transcripts.column('hgnc_id')
columns = ('gnomad_joint_af', 'gnomad_joint_ac', 'gnomad_joint_hom', 'gnomad_joint_faf95_max')
kinds = [variants.schema.field(c).type for c in columns]
assert kinds == [pa.float32(), pa.int32(), pa.int32(), pa.float32()], kinds
joint = {r['variant_id']: tuple(r[c] for c in columns) for r in variants.to_pylist()}
absent = (None, None, None, None)
assert joint == {
    'chr1-100-G-A': (float32(1.85e-05), 3, 1, float32(4e-06)),
    'chr1-100-G-C': absent,
    'chr1-200-G-A': (0.0, 0, 0, 0.0),
    'chr1-300-T-A': absent,
    'chr1-400-TCA-T': (float32(0.074321), 120000, 4500, float32(0.0625015)),
    'chr2-40-C-T': (0.25, 250000, 30000, 0.25),
    'chr2-60-C-A': absent,
    'chr2-80-C-T': absent,
    'chr2-90-GCGTAAA-G': (float32(4e-06), 7, 0, float32(1e-06)),
    'chrX-100-T-A': absent,
    'chrX-100-T-C': (float32(0.001), 25, 1, float32(0.0008)),
}, joint
recorded = manifest['gnomad_joint']
archive = open('${archive}', 'rb').read()
assert (recorded['archive'], recorded['bytes']) == ('invented.zip', len(archive)), recorded
assert recorded['sha256'] == hashlib.sha256(archive).hexdigest(), recorded
assert recorded['echtvar_version'] == '0.2.4', recorded
fields = [(f['field'], f['alias'], f['ftype'], f['multiplier'], f['missing_value']) for f in recorded['fields']]
assert fields == [
    ('AC_joint', 'gnomad_AC', 'Integer', 1, -2147483648),
    ('AF_joint', 'gnomad_AF', 'Float', 2000000, 0),
    ('nhomalt_joint', 'gnomad_HomAlt', 'Integer', 1, -2147483648),
    ('fafmax_faf95_max_joint', 'gnomad_faf95max', 'Float', 2000000, 0),
], fields
shown = ('avi_phred', 'avi_raw', 'avi_shap_atac', 'avi_shap_procap', 'alphagenome_splicing')
scores = {r['variant_id']: tuple(r[c] for c in shown) for r in variants.to_pylist()}
none = (None,) * len(shown)
expected = {
    'chr1-100-G-A': (1.20713, -0.04245, 0.001293, 0.000053, 0.002873),
    'chr1-100-G-C': (0.0, -1.25109, 0.00062, 0.000441, 0.000015),
    'chr1-200-G-A': (97.31337, 2.44971, None, None, None),
    'chr1-300-T-A': none,
    'chr1-400-TCA-T': none,
    'chr2-40-C-T': none,
    'chr2-60-C-A': (48.90215, 2.75691, 0.003406, 0.007602, 1.207),
    'chr2-80-C-T': none,
    'chr2-90-GCGTAAA-G': none,
    'chrX-100-T-A': (2.50794, -0.00597, 0.002019, 0.002793, None),
    'chrX-100-T-C': (20.00001, 0.55782, 0.01204, 0.09413, None),
}
assert scores == {
    v: tuple(None if x is None else float32(x) for x in values) for v, values in expected.items()}, scores
contributions = [c for c in variants.column_names if c.startswith('avi_shap_')]
assert len(contributions) == 16, contributions
added = ('avi_phred', 'avi_raw', *contributions, 'alphagenome_splicing')
at = variants.column_names.index('gnomad_joint_faf95_max')
assert tuple(variants.column_names[at + 1:at + 1 + len(added)]) == added, variants.column_names
assert {variants.schema.field(c).type for c in added} == {pa.float32()}, variants.schema
broken_down = [r for r in variants.to_pylist() if r['avi_shap_atac'] is not None]
assert len(broken_down) == 5, broken_down
assert all(abs(sum(r[c] for c in contributions) - r['avi_raw'] - 0.048) < 0.00001 for r in broken_down), broken_down
for name in ('alphagenome_avi', 'alphagenome_avi_breakdown', 'alphagenome_splicing'):
    recorded = manifest[name]
    assert (recorded['archive'], recorded['echtvar_version']) == (name + '.zip', '0.2.4'), recorded
    assert ('terms' in recorded) == (name != 'alphagenome_avi'), recorded
    assert json.loads(variants.schema.metadata[name.encode()]) == recorded, variants.schema.metadata
notice = 'By using this information, you agree to AlphaGenome Output Terms of Use found at '
assert manifest['alphagenome_splicing']['terms'].startswith(notice), manifest
assert not set(added) & set(manifest['omitted_columns']), manifest
assert variants.column_names[-6] == 'utr5_consequences', variants.column_names
kinds = {r['variant_id']: r['utr5_consequences'] for r in variants.to_pylist() if r['utr5_consequences']}
assert kinds == {'chr1-200-G-A': ['uSTOP_gained']}, kinds
plugins = ('utr5_consequences', 'utr5_uorfs', 'utr5_existing_uorfs', 'utr5_existing_oorfs_out_of_frame',
           'utr5_existing_oorfs_in_frame', 'nmd_escaping')
beside = ('utr5_cds_kozak_context', 'utr5_cds_kozak_strength', 'utr5_cds_kozak_context_alt',
          'utr5_cds_kozak_strength_alt', 'nmd_escape_reason')
assert transcripts.column_names[-11:] == [*plugins, *beside], transcripts.column_names
said = {r['variant_id']: tuple(r[c] for c in plugins) for r in transcripts.to_pylist()
        if any(r[c] is not None for c in plugins)}
assert set(said) == {'chr1-100-G-A', 'chr1-100-G-C', 'chr1-200-G-A', 'chr1-400-TCA-T'}, said
assert said['chr1-100-G-A'] == said['chr1-100-G-C'] == (None, None, 0, 0, 1, None), said
assert said['chr1-400-TCA-T'] == (None, None, None, None, None, False), said
gained = said['chr1-200-G-A']
assert (gained[0], gained[2:]) == (['uSTOP_gained'], (0, 0, 1, None)), gained
(record,) = gained[1]
assert (record['consequence'], record['ref_type'], record['new_stop_distance_to_cds'], record['evidence']) == (
    'uSTOP_gained', 'InFrame_oORF', 12, False), record
assert record['type'] is None and record['distance_to_cds'] is None, record
ran = manifest['vep_plugins']
assert ran['run']['vep_release'] == '116.3', ran
assert sorted(ran['run']) == ['gene_models', 'nmd_plugin_commit', 'uorf_evidence', 'utrannotator_commit',
                              'vep_release'], ran
assert ran['matching'] == [
    {'source': 'Ensembl', 'rows_without_transcript': 0, 'rows_of_other_ids': 0, 'utr5_transcripts_without_values': 0},
    {'source': 'RefSeq', 'rows_without_transcript': 4, 'rows_of_other_ids': 0, 'utr5_transcripts_without_values': 0},
], ran
assert json.loads(variants.schema.metadata[b'vep_plugins']) == ran['run'], variants.schema.metadata
EOF
    """
}

// An invented union on chrM, bgzipped with its index, as UNION_SITES writes a union.
process INVENTED_MITO_UNION {
    input:
    path vcf

    output:
    tuple val('mito'), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
    """
    bcftools view --no-version -Oz -o ${vcf}.gz ${vcf}
    bcftools index -t ${vcf}.gz
    """
}

// Family A's values from the two files on polyadenylation and from gnomAD's three tables of constraint, after the
// columns of the echtvar archives: the flag of a signal, which is true for the three alleles that change a base of one
// and false for every other allele; APARENT2's lowest and highest effect with the usage of the lowest row's site, for
// the four single-base changes that the invented table scores; the z score of the window of each allele on chr1 and
// chr2; the lowest ratio among the regions that hold the allele's position; and the highest MPC score of the allele.
// The last two also count a transcript that the annotation lacks. transcripts holds the ratio, the number of regions,
// the p-value, the coverage and the MPC score on the row of the allele's transcript, which the files write without its
// version, and nothing of the transcript that the annotation lacks. After them variants holds the branch points whose
// three bases an allele changes, as a list: one entry for the single-base change at a branch point, with the intron
// that the invented catalogue gives it; one for the deletion whose last base is a branch point of the reverse strand,
// which the catalogue gives no intron; and three for the deletion over two branch points, two of them for the two
// introns of the first. Every other allele has no list. The manifest and the tables' metadata record the six files by
// name, size and SHA-256, and the manifest lists the columns of gnomAD's mitochondrial sites and of MitoTIP's scores,
// which the run did not name, as left out, after those of the dated sources.
process CHECK_JOINED {
    input:
    tuple path(variants), path(transcripts), path(manifest)
    path polya
    path effects
    path windows
    path regions
    path scores
    path catalogue

    script:
    """
    python3 - <<'EOF'
import hashlib
import json

import pyarrow as pa
import pyarrow.parquet as pq

variants = pq.read_table('${variants}')
added = ('gencode_polya_signal', 'aparent2_min', 'aparent2_max', 'aparent2_ref_usage', 'genomic_constraint_z',
         'mcr_oe_min', 'mpc_max')
at = variants.column_names.index('alphagenome_splicing')
assert tuple(variants.column_names[at + 1:at + 1 + len(added)]) == added, variants.column_names
after = variants.column_names[at + 1 + len(added):at + 3 + len(added)]
assert after == ['branch_points', 'utr5_consequences'], variants.column_names
kinds = [variants.schema.field(c).type for c in added]
assert kinds == [pa.bool_()] + [pa.float32()] * 6, kinds
assert not variants.schema.field('gencode_polya_signal').nullable, variants.schema
rows = {r['variant_id']: tuple(r[c] for c in added) for r in variants.to_pylist()}
assert rows == {
    'chr1-100-G-A': (False, -2.25, 0.75, 0.5, 2.5, 0.25, 3.25),
    'chr1-100-G-C': (False, 0.03125, 0.03125, 0.5, 2.5, 0.25, 0.5),
    'chr1-200-G-A': (False, None, None, None, 2.5, 0.75, None),
    'chr1-300-T-A': (True, None, None, None, 2.5, 0.75, None),
    'chr1-400-TCA-T': (True, None, None, None, 2.5, None, None),
    'chr2-40-C-T': (False, None, None, None, -1.25, None, None),
    'chr2-60-C-A': (False, -0.5, -0.5, 0.25, -1.25, 0.125, None),
    'chr2-80-C-T': (True, None, None, None, -1.25, 0.125, None),
    'chr2-90-GCGTAAA-G': (False, None, None, None, -1.25, None, None),
    'chrX-100-T-A': (False, None, None, None, None, None, None),
    'chrX-100-T-C': (False, 1.5, 1.5, 0.0625, None, None, 6.0),
}, rows
entry = variants.schema.field('branch_points').type.value_type
assert [(f.name, str(f.type), f.nullable) for f in entry] == [
    ('name', 'string', False), ('pos', 'int32', False), ('strand', 'string', False),
    ('hit', 'list<element: int8>', False), ('acceptor_offset', 'int32', True), ('rank', 'int32', True),
    ('total', 'int32', True), ('sources', 'list<element: string>', False), ('consensus', 'string', True),
    ('gene', 'string', True), ('transcript_introns', 'list<element: string>', True), ('intron_type', 'string', True),
    ('intron_length', 'int32', True)], entry
found = {r['variant_id']: r['branch_points'] for r in variants.to_pylist() if r['branch_points'] is not None}
alone = {'acceptor_offset': None, 'rank': None, 'total': None, 'consensus': None, 'gene': None,
         'transcript_introns': None, 'intron_type': None, 'intron_length': None}
both = ['eBP_Invented', 'cBP_Invented']
assert found == {
    'chr1-200-G-A': [
        {'name': 'mBP_chr1_200_+_G', 'pos': 200, 'strand': '+', 'hit': [0], 'acceptor_offset': -24, 'rank': 1,
         'total': 1, 'sources': both, 'consensus': '1:YTNAY', 'gene': 'INVENTED1',
         'transcript_introns': ['ENST01_IVS2'], 'intron_type': 'major', 'intron_length': 900}],
    'chr1-400-TCA-T': [
        {'name': 'cBP_chr1_402_-_A', 'pos': 402, 'strand': '-', 'hit': [0], 'sources': ['cBP_Invented'], **alone}],
    'chr2-90-GCGTAAA-G': [
        {'name': 'eBP_chr2_92_+_A', 'pos': 92, 'strand': '+', 'hit': [-1, 0], 'acceptor_offset': -20, 'rank': 1,
         'total': 2, 'sources': ['eBP_Invented'], 'consensus': '0:none', 'gene': 'INVENTED2',
         'transcript_introns': ['ENST01_IVS4'], 'intron_type': 'minor', 'intron_length': 85},
        {'name': 'eBP_chr2_92_+_A', 'pos': 92, 'strand': '+', 'hit': [-1, 0], 'acceptor_offset': -31, 'rank': 2,
         'total': 2, 'sources': ['eBP_Invented'], 'consensus': '3:TNA', 'gene': 'INVENTED2',
         'transcript_introns': ['ENST02_IVS1', 'ENST03_IVS1'], 'intron_type': 'major', 'intron_length': 4000},
        {'name': 'cBP_chr2_94_+_C', 'pos': 94, 'strand': '+', 'hit': [-2, -1, 0], 'sources': ['cBP_Invented'],
         **alone}],
}, found
transcripts = pq.read_table('${transcripts}')
by_transcript = ('mcr_oe', 'mcr_regions', 'mcr_oe_p', 'mcr_high_coverage', 'mpc')
plugins = transcripts.column_names.index('utr5_consequences')
assert tuple(transcripts.column_names[plugins - len(by_transcript):plugins]) == by_transcript, transcripts.column_names
kinds = [transcripts.schema.field(c).type for c in by_transcript]
assert kinds == [pa.float32(), pa.int32(), pa.float64(), pa.bool_(), pa.float32()], kinds
held = {(r['variant_id'], r['transcript_id']): tuple(r[c] for c in by_transcript) for r in transcripts.to_pylist()}
nothing = (None, None, None, None, None)
assert held == {
    ('chr1-100-G-A', 'ENST01.1'): (0.25, 2, 2.1e-60, True, 2.5),
    ('chr1-100-G-C', 'ENST01.1'): (0.25, 2, 2.1e-60, True, 0.5),
    ('chr1-200-G-A', 'ENST01.1'): (0.75, 2, 0.125, False, None),
    ('chr1-300-T-A', 'ENST01.1'): (0.75, 2, 0.125, False, None),
    ('chr1-400-TCA-T', 'ENST01.1'): nothing,
    ('chr2-40-C-T', 'ENST01.1'): nothing,
    ('chr2-60-C-A', 'ENST01.1'): nothing,
    ('chr2-80-C-T', 'ENST01.1'): nothing,
    ('chr2-90-GCGTAAA-G', 'ENST01.1'): nothing,
    ('chrX-100-T-A', 'ENST01.1'): nothing,
    ('chrX-100-T-C', 'ENST01.1'): (None, None, None, None, 6.0),
}, held
manifest = json.load(open('${manifest}'))
files = {'gencode_polya': '${polya}', 'aparent2': '${effects}', 'genomic_constraint': '${windows}',
         'mcr': '${regions}', 'mpc': '${scores}', 'branch_points': '${catalogue}'}
assert manifest['joined_files'] == [
    {'source': source, 'file': name, 'bytes': len(open(name, 'rb').read()),
     'sha256': hashlib.sha256(open(name, 'rb').read()).hexdigest()}
    for source, name in files.items()], manifest
for table in (variants, transcripts):
    assert json.loads(table.schema.metadata[b'joined_files']) == manifest['joined_files'], table.schema.metadata
left_out = manifest['omitted_columns'][-15:]
assert [c for c in left_out if c.startswith('gnomad_mito_')] == left_out[:13], manifest
assert left_out[13:] == ['mitotip_score', 'mitotip_quartile'], manifest
assert manifest['transcripts_omitted_columns'] == [], manifest
assert not [c for c in variants.column_names if c.startswith(('gnomad_mito_', 'mitotip_'))], variants.column_names
EOF
    """
}

// The values of the invented union on chrM from the invented files of gnomAD's mitochondrial sites and of MitoTIP's
// scores: the single-base change has both files' values, the deletion gnomAD's values of the record that the file
// writes six bases to the right and no score, although MitoTIP's file scores the deletion of a base that it removes,
// and the allele that both files lack none. No value is by transcript. The records are of the files as they were given.
process CHECK_JOINED_MITO {
    input:
    path values
    path mito
    path mitotip

    script:
    """
    python3 - <<'EOF'
import hashlib
import json
import struct

import pyarrow.parquet as pq


def float32(value):
    return struct.unpack('f', struct.pack('f', value))[0]


table = pq.read_table('${values}/joined_values.parquet')
assert table.column_names[0] == 'variant_id' and len(table.column_names) == 16, table.column_names
assert all(c.startswith('gnomad_mito_') for c in table.column_names[1:14]), table.column_names
assert table.column_names[14:] == ['mitotip_score', 'mitotip_quartile'], table.column_names
rows = {r.pop('variant_id'): {c.removeprefix('gnomad_mito_'): v for c, v in r.items()} for r in table.to_pylist()}
assert list(rows) == ['chrM-10-A-C', 'chrM-12-C-T', 'chrM-30-GCA-G'], rows
assert rows['chrM-10-A-C'] == {
    'ac_hom': 8, 'ac_het': 2, 'af_hom': float32(0.002), 'af_het': float32(0.0005), 'an': 4000,
    'max_heteroplasmy': 0.5, 'hapmax_af_hom': 0.125, 'hapmax_af_het': 0.03125, 'hapmax_hom_haplogroup': 'B',
    'hapmax_het_haplogroup': 'B', 'faf_hapmax_hom': 0.0625, 'filters': ['PASS'], 'hap_defining': True,
    'mitotip_score': 12.5, 'mitotip_quartile': 'Q3'}, rows
assert set(rows['chrM-12-C-T'].values()) == {None}, rows
deletion = rows['chrM-30-GCA-G']
assert (deletion['an'], deletion['ac_hom'], deletion['af_hom'], deletion['filters']) == (4000, 0, 0.0, ['npg']), rows
assert (deletion['hap_defining'], deletion['hapmax_hom_haplogroup']) == (False, None), rows
assert (deletion['mitotip_score'], deletion['mitotip_quartile']) == (None, None), rows
by_transcript = pq.read_table('${values}/joined_transcript_values.parquet')
assert (by_transcript.column_names, by_transcript.num_rows) == (['variant_id', 'transcript_id'], 0), by_transcript
files = {'gnomad_mito': '${mito}', 'mitotip': '${mitotip}'}
assert json.load(open('${values}/joined_values.json')) == {'files': [
    {'source': source, 'file': name, 'bytes': len(open(name, 'rb').read()),
     'sha256': hashlib.sha256(open(name, 'rb').read()).hexdigest()}
    for source, name in files.items()]}
EOF
    """
}

// The values of a run that names no archive of AlphaGenome's scores: gnomAD's joint columns alone, with the record of
// gnomAD's archive and of no other.
process CHECK_GNOMAD_ALONE {
    input:
    path values

    script:
    """
    python3 - <<'EOF'
import json

import pyarrow.parquet as pq

columns = pq.read_schema('${values}/echtvar_values.parquet').names
assert columns == ['variant_id', 'gnomad_joint_af', 'gnomad_joint_ac', 'gnomad_joint_hom',
                   'gnomad_joint_faf95_max'], columns
records = json.load(open('${values}/echtvar_values.json'))
assert len(records) == 4, records
assert [name for name, record in records.items() if record is not None] == ['gnomad_joint'], records
EOF
    """
}

// Family B's manifest: the sheet names no proband, so nobody has a role, although the PED file has an affected son with
// both parents. The family has no repeats table, so its manifest records no release of STRchive and no column left out
// of that table.
process CHECK_NO_ROLES {
    input:
    path manifest

    script:
    """
    python3 - <<'EOF'
import json

manifest = json.load(open('${manifest}'))
members = manifest['members']
assert {m['sample_id']: m['role'] for m in members} == {'B1': None, 'B2': None, 'B3': None}, members
assert [m['sample_id'] for m in members if m['affected']] == ['B3'], members
assert (manifest['strchive'], manifest['repeats_omitted_columns']) == (None, []), manifest
EOF
    """
}

// Family A's tables of its larger variants. sv holds both callers' records in one order, without the record off the
// primary chromosomes and the segment at the reference's copy number, and after a record's own columns and the
// annotation's a call column per member in the order of the family's VCF, as variants does, with the fields that the
// two tables share under the same names. Each record of sv has the annotation of its site: the HGNC id of its
// transcript's gene; the symbolic deletion its best gnomAD match, which overlaps it most, the highest frequency among
// its matches and both of them in full; the deletion of six bases, annotated as a small variant, gnomAD's values of
// the allele itself; the segment of two alleles the matches of both alleles and no best match. sv has no ClinVar
// column, since the run named no ClinVar, and the manifest lists the columns left out. sv_transcripts holds each
// allele's transcript with its gene's HGNC id, two for the segment of two alleles. The deletion written at another
// place of the repeat and the deletion of six bases each name their allele in variants, and the third deletion, the
// symbolic alleles, the breakends and the segments name none. Each breakend has its mate's place. The father carries
// the deletion on his one X by its haploid genotype, and the gain there by a copy number of 2 where 1 is expected,
// with no zygosity; the mother, at 2 copies, carries neither segment. repeats holds each member's sizes of both
// repeats, one for the father on X, and STRchive's values of the repeat on chr1, which is of an invented locus whose
// interval is two units longer than the repeat's. The repeat on X is of none, and the manifest records the release as
// the run stated it, with both counts.
process CHECK_LARGER {
    input:
    tuple path(variants), path(manifest), path(larger)

    script:
    """
    python3 - <<'EOF'
import json
import struct

import pyarrow as pa
import pyarrow.parquet as pq


def float32(value):
    return struct.unpack('f', struct.pack('f', value))[0]


manifest = json.load(open('${manifest}'))
rows = {t['table']: t['rows'] for t in manifest['tables']}
assert rows == {'variants': 11, 'transcripts': 11, 'sv': 10, 'sv_transcripts': 11, 'repeats': 2, 'roh': 0, 'upd': 0}, rows
assert manifest['sv_callers'] == ['sv', 'cnv'], manifest
assert manifest['schema_version'] == '0.5.0', manifest
assert manifest['sv_annotation']['sources'] == manifest['annotation']['sources'] == ['Ensembl', 'gnomad'], manifest
table = pq.read_table('sv.parquet')
assert table.column_names[:4] == ['chrom', 'pos', 'end_pos', 'sv_id'], table.column_names
assert table.column_names[-4:] == ['gnomad_exomes_failed_filter', 'A1', 'A2', 'A3'], table.column_names
at = table.column_names.index('small_variant_id')
assert table.column_names[at + 1:at + 5] == ['annotated_as', 'annotated_pos', 'annotated_end_pos', 'hgnc_ids']
shared = ('gt', 'carries', 'zygosity', 'gq', 'pl', 'ft', 'dn', 'dq')
small = pq.read_schema('${variants}').field('A1').type
assert all(table.schema.field('A1').type.field(f).equals(small.field(f)) for f in shared), table.schema
assert manifest['sv_omitted_columns'] == [f'clinvar_{c}' for c in (
    'accession', 'classification', 'review_status', 'stars', 'conditions', 'last_evaluated', 'overlap', 'records')]
assert not [c for c in table.column_names if c.startswith('clinvar')], table.column_names
assert table.schema.field('gnomad_sv_af').type == pa.float32(), table.schema
assert table.schema.field('hgnc_ids').equals(pq.read_schema('${variants}').field('hgnc_ids')), table.schema
sv = table.to_pylist()
order = ['inv:LOSS:1', 'inv:BND:1:0', 'inv:DEL:1', 'inv:DEL:4', 'inv:DEL:2', 'inv:DEL:3', 'inv:BND:1:1', 'inv:CNV:1',
         'inv:DEL:5', 'inv:GAIN:1']
assert [r['sv_id'] for r in sv] == order, sv
by_id = {r['sv_id']: r for r in sv}
linked = {r['sv_id']: r['small_variant_id'] for r in sv if r['small_variant_id'] is not None}
assert linked == {'inv:DEL:1': 'chr1-400-TCA-T', 'inv:DEL:2': 'chr2-90-GCGTAAA-G'}, linked
assert set(linked.values()) <= {v['variant_id'] for v in pq.read_table('${variants}').to_pylist()}
# The record keeps the place its caller wrote it at.
assert (by_id['inv:DEL:1']['pos'], by_id['inv:DEL:1']['end_pos'], by_id['inv:DEL:1']['ref']) == (406, 408, 'ACA')
assert {r['sv_id']: (r['caller'], r['sv_type']) for r in sv if r['caller'] == 'cnv'} == {
    'inv:LOSS:1': ('cnv', 'DEL'), 'inv:CNV:1': ('cnv', 'CNV'), 'inv:GAIN:1': ('cnv', 'DUP')}, sv
assert {tuple(r['hgnc_ids']) for r in sv} == {(90001,)}, sv
treated = {r['sv_id']: r['annotated_as'] for r in sv}
assert {i for i, how in treated.items() if how == 'small_variant'} == {'inv:DEL:1', 'inv:DEL:2', 'inv:DEL:3'}, treated
assert set(treated.values()) == {'structural', 'small_variant'}, treated
matched = by_id['inv:DEL:4']
assert (matched['gnomad_sv_id'], matched['gnomad_sv_overlap'], matched['gnomad_sv_distance']) == ('inv-best', 1.0, None)
assert (matched['gnomad_sv_af'], matched['gnomad_sv_max_af']) == (float32(0.001), 0.25), matched
assert [(m['id'], m['pos'], m['end_pos']) for m in matched['gnomad_sv_matches']] == [
    ('inv-best', 500, 560), ('inv-common', 479, 570)], matched
assert (matched['annotated_pos'], matched['annotated_end_pos']) == (500, 560), matched
small = by_id['inv:DEL:2']
assert (small['gnomad_genomes_af'], small['gnomad_genomes_an'], small['gnomad_sv_id']) == (0.125, 8, None), small
segment = by_id['inv:CNV:1']
assert [(m['id'], m['variant_type']) for m in segment['gnomad_sv_matches']] == [
    ('inv-copy_number_loss', 'copy_number_loss'), ('inv-copy_number_gain', 'copy_number_gain')], segment
assert (segment['gnomad_sv_id'], segment['gnomad_sv_af'], segment['gnomad_sv_max_af']) == (None, None, None), segment
assert (segment['annotated_as'], segment['annotated_pos'], segment['annotated_end_pos']) == ('structural', 300, 450)
assert by_id['inv:BND:1:0']['gnomad_sv_matches'] is None, by_id['inv:BND:1:0']
transcripts = pq.read_table('sv_transcripts.parquet')
assert transcripts.column_names[:6] == ['chrom', 'pos', 'end_pos', 'sv_id', 'variant_type', 'transcript_id']
assert transcripts.schema.field('hgnc_id').type == pa.int32(), transcripts.schema
transcripts = transcripts.to_pylist()
assert list(dict.fromkeys(t['sv_id'] for t in transcripts)) == order, transcripts
assert [(t['variant_type'], t['transcript_id']) for t in transcripts if t['sv_id'] == 'inv:CNV:1'] == [
    ('copy_number_gain', 'ENST-copy_number_gain.1'), ('copy_number_loss', 'ENST-copy_number_loss.1')], transcripts
assert {tuple(t['consequences']) for t in transcripts} == {('feature_truncation', 'intron_variant')}, transcripts
assert {(t['gene_symbol'], t['hgnc_id']) for t in transcripts} == {('GENEA', 90001)}, transcripts
assert all((t['chrom'], t['pos'], t['end_pos']) == (by_id[t['sv_id']]['chrom'], by_id[t['sv_id']]['pos'],
                                                    by_id[t['sv_id']]['end_pos']) for t in transcripts)
mate = by_id['inv:BND:1:0']
assert (mate['mate_id'], mate['mate_chrom'], mate['mate_pos'], mate['bnd_extends'], mate['mate_extends']) == (
    'inv:BND:1:1', 'chr2', 200, 'left', 'right'), mate
symbolic = by_id['inv:DEL:4']
assert (symbolic['end_pos'], symbolic['svlen'], symbolic['filters'], symbolic['imprecise'], symbolic['cipos']) == (
    560, -60, ['MinQUAL'], True, [-5, 5]), symbolic
seen = {s: (by_id['inv:DEL:5'][s]['gt'], by_id['inv:DEL:5'][s]['carries'], by_id['inv:DEL:5'][s]['zygosity'])
        for s in ('A1', 'A2', 'A3')}
assert seen == {'A1': ('1', True, 'hom_alt'), 'A2': ('0/1', True, 'het'), 'A3': ('0/0', False, 'hom_ref')}, seen
assert by_id['inv:DEL:5']['A1']['pl'] == [70, 0], by_id['inv:DEL:5']
copies = {sv_id: {s: (by_id[sv_id][s]['gt'], by_id[sv_id][s]['cn'], by_id[sv_id][s]['carries'],
                      by_id[sv_id][s]['zygosity']) for s in ('A1', 'A2', 'A3')}
          for sv_id in ('inv:LOSS:1', 'inv:GAIN:1')}
assert copies == {
    'inv:LOSS:1': {'A1': ('0/1', 1, True, 'het'), 'A2': ('./.', 2, False, 'hom_ref'), 'A3': ('1/1', 0, True, 'hom_alt')},
    'inv:GAIN:1': {'A1': ('1', 2, True, None), 'A2': ('./.', 2, False, 'hom_ref'), 'A3': ('./1', 3, True, None)},
}, copies
repeats = pq.read_table('repeats.parquet')
ranges = ['strchive_id', 'strchive_gene', 'strchive_evidence', 'strchive_inheritance', 'strchive_start_hg38',
          'strchive_stop_hg38', 'strchive_units_offset', 'strchive_pathogenic_motif_reference_orientation',
          'strchive_benign_min', 'strchive_benign_max', 'strchive_intermediate_min', 'strchive_intermediate_max',
          'strchive_pathogenic_min', 'strchive_pathogenic_max']
assert repeats.column_names == ['chrom', 'pos', 'end_pos', 'repeat_id', 'repeat_unit', 'ref_copies', *ranges, 'A1',
                                'A2', 'A3'], repeats.column_names
assert repeats.schema.field('strchive_pathogenic_min').type == pa.int32(), repeats.schema
assert repeats.schema.field('strchive_units_offset').type == pa.float32(), repeats.schema
assert json.loads(repeats.schema.metadata[b'strchive']) == manifest['strchive'], repeats.schema.metadata
repeats = {r['repeat_id']: r for r in repeats.to_pylist()}
assert [repeats['RPTA'][c] for c in ranges] == ['INVA_GENEA', 'GENEA', 'Definitive', ['AD'], 47, 83, 2.0, ['AGC'],
                                                5, 20, 21, 35, 36, 200], repeats
assert [repeats['RPTX_CGG'][c] for c in ranges] == [None] * len(ranges), repeats
assert {k: v for k, v in manifest['strchive'].items() if k not in ('bytes', 'sha256')} == {
    'file': 'strchive.json', 'release': 'v0.0.0', 'release_date': '2099-01-01', 'matched_repeats': 1,
    'unmatched_repeats': 1}, manifest
assert manifest['repeats_omitted_columns'] == [], manifest
# The daughter's longer allele of the repeat lies in the locus's pathogenic range.
assert repeats['RPTA']['strchive_pathogenic_min'] <= max(repeats['RPTA']['A3']['repcn']) <= 200, repeats
assert list(repeats) == ['RPTA', 'RPTX_CGG'], repeats
assert (repeats['RPTX_CGG']['repeat_unit'], repeats['RPTX_CGG']['ref_copies']) == ('CGG', 20), repeats
sizes = {s: repeats['RPTX_CGG'][s]['repcn'] for s in ('A1', 'A2', 'A3')}
assert sizes == {'A1': [23], 'A2': [23, 30], 'A3': [23, 23]}, sizes
assert repeats['RPTA']['A3']['repci_upper'] == [12, 58], repeats
assert repeats['RPTA']['A3']['so'] == ['SPANNING', 'INREPEAT'], repeats
EOF
    """
}

// Family A's runs of homozygosity and regions of uniparental disomy, as its folder holds them. The family's VCF has
// four sites that both tools read, one of them with a frequency that they keep, so both tables are there without a
// row: roh as a table of every family, and upd because the daughter has both parents sequenced. The manifest lists
// both after the other tables and records what ran: the versions of bcftools, echtvar and upd, the commit that upd was
// installed from, the commands, the archive with the field of the frequency, the members and the assessed child. Each
// table carries the same record.
process CHECK_SEGMENTS {
    input:
    tuple path(roh), path(manifest), path(larger)
    path archive

    script:
    """
    python3 - <<'EOF'
import hashlib
import json
import re

import pyarrow.parquet as pq

manifest = json.load(open('${manifest}'))
assert [t['table'] for t in manifest['tables']][-2:] == ['roh', 'upd'], manifest
roh, upd = pq.read_table('${roh}'), pq.read_table('upd.parquet')
assert roh.column_names == ['member', 'chrom', 'pos', 'end_pos', 'length', 'markers', 'quality'], roh.column_names
assert upd.column_names == ['member', 'chrom', 'pos', 'end_pos', 'outer_pos', 'outer_end_pos', 'origin',
                            'informative_sites', 'snps', 'opposite_sites'], upd.column_names
assert (roh.num_rows, upd.num_rows) == (0, 0), (roh, upd)
record = manifest['roh_upd']
assert (record['bcftools_version'], record['upd_version']) == ('1.21', '0.1'), record
assert re.fullmatch('[0-9a-f]{40}', record['upd_commit']), record
assert record['preparation'] == [
    'bcftools view -f PASS,. --trim-alt-alleles -t ' + ','.join(f'chr{n}' for n in range(1, 23)),
    'bcftools view -m2 -M2 -v snps', 'bcftools annotate -x INFO', 'bcftools +setGT -- -t a -n u', 'echtvar anno'], record
assert record['roh_command'] == 'bcftools roh --AF-tag gnomad_AF -G30 -O r', record
assert record['upd_command'] == 'upd --af-tag gnomad_AF --min-af 0.05 --min-gq 30 regions', record
archive = open('${archive}', 'rb').read()
used = record['frequency_archive']
assert (used['archive'], used['bytes'], used['sha256']) == (
    'invented.zip', len(archive), hashlib.sha256(archive).hexdigest()), used
assert (used['echtvar_version'], used['field']['field'], used['field']['alias']) == (
    '0.2.4', 'AF_joint', 'gnomad_AF'), used
# The site with an allele that nobody carries, two sites of one allele and the site that the father has a call of
# ./1 at; the indels, the sites of two carried alleles and the site on X are left out.
assert record['sites'] == 4, record
assert record['roh_members'] == ['A1', 'A2', 'A3'], record
assert record['upd_children'] == [{'child': 'A3', 'mother': 'A2', 'father': 'A1'}], record
assert record['upd_not_assessed'] == [], record
for table in (roh, upd):
    assert json.loads(table.schema.metadata[b'roh_upd']) == record, table.schema.metadata
EOF
    """
}

// What ROH_UPD writes for the invented trio (roh-upd-trio.vcf, with the frequencies of roh-upd-gnomad.vcf) and for its
// mother and child alone. On chr1 the child has the mother's allele twice at 22 sites from 1000 to 3000, between two
// sites that show both parents' alleles: upd reports the region with its outer bounds, and bcftools roh a run of the
// child over it. On chr2 the child has both of the mother's alleles at 16 sites from 1000 to 2500: upd reports the
// region up to its last site that shows the mother's alleles alone, and no run of the child lies over it. Read
// together, the first region is an isodisomy or a deletion and the second a heterodisomy. The regions are whole
// although sites inside them show both parents' alleles: a site that the archive lacks, one with a frequency below
// 0.05, one that failed a filter, an insertion, a site of three carried alleles and one with a GQ below 30 are not
// read. The sites of the trio's VCF with phased genotypes, with an allele that nobody carries and with an INFO value
// count among the 22, and the VCF's INFO header line with keys of its own does not stop upd. The mother and the child
// alone have the child's run and no table of regions.
process CHECK_ROH_UPD {
    input:
    path trio, stageAs: 'trio/*'
    path duo, stageAs: 'duo/*'

    script:
    """
    python3 - <<'EOF'
import json
import os

import duckdb
import pyarrow.parquet as pq

assert sorted(os.listdir('trio/roh_upd')) == ['roh.parquet', 'roh_upd.json', 'upd.parquet']
roh = pq.read_table('trio/roh_upd/roh.parquet').to_pylist()
assert all(r['length'] == r['end_pos'] - r['pos'] + 1 for r in roh), roh
(run,) = [r for r in roh if r['member'] == 'T3']
assert run['chrom'] == 'chr1' and run['pos'] <= 1000 and run['end_pos'] >= 3000 and run['markers'] >= 22, run
regions = pq.read_table('trio/roh_upd/upd.parquet').to_pylist()
assert regions == [
    {'member': 'T3', 'chrom': 'chr1', 'pos': 1000, 'end_pos': 3000, 'outer_pos': 900, 'outer_end_pos': 3099,
     'origin': 'maternal', 'informative_sites': 22, 'snps': 29, 'opposite_sites': 0},
    {'member': 'T3', 'chrom': 'chr2', 'pos': 1000, 'end_pos': 2400, 'outer_pos': 900, 'outer_end_pos': 2599,
     'origin': 'maternal', 'informative_sites': 8, 'snps': 21, 'opposite_sites': 0},
], regions
record = json.load(open('trio/roh_upd/roh_upd.json'))
assert record['roh_members'] == ['T1', 'T2', 'T3'], record
assert record['upd_children'] == [{'child': 'T3', 'mother': 'T2', 'father': 'T1'}], record
assert record['upd_not_assessed'] == [], record
# The share of each region that a run of the same member covers.
query = (
    'SELECT u.chrom, coalesce(sum(least(u.end_pos, r.end_pos) - greatest(u.pos, r.pos) + 1) '
    'FILTER (WHERE r.pos IS NOT NULL), 0) / (u.end_pos - u.pos + 1) '
    "FROM 'trio/roh_upd/upd.parquet' u LEFT JOIN 'trio/roh_upd/roh.parquet' r "
    'ON r.member = u.member AND r.chrom = u.chrom AND r.pos <= u.end_pos AND r.end_pos >= u.pos '
    'GROUP BY u.member, u.chrom, u.pos, u.end_pos ORDER BY u.chrom')
covered = duckdb.sql(query).fetchall()
assert covered == [('chr1', 1.0), ('chr2', 0.0)], covered
assert sorted(os.listdir('duo/roh_upd')) == ['roh.parquet', 'roh_upd.json']
alone = json.load(open('duo/roh_upd/roh_upd.json'))
assert (alone['roh_members'], alone['upd_children']) == (['T2', 'T3'], []), alone
runs = [r for r in pq.read_table('duo/roh_upd/roh.parquet').to_pylist() if r['member'] == 'T3']
assert [(r['chrom'], r['pos'] <= 1000, r['end_pos'] >= 3000) for r in runs] == [('chr1', True, True)], runs
EOF
    """
}

// What ROH_UPD writes where upd would keep no site, in which case upd is not run. The invented trio with the archive of
// the families' alleles has 74 sites, none of them with a frequency that upd keeps: with its own frequencies the same
// trio is assessed (CHECK_ROH_UPD). The trio's record on X alone leaves no site at all. In both the child is recorded
// as not assessed, with the reason, and not as assessed, and the family has the table of the regions without a row.
// Neither has a run of homozygosity, which bcftools roh needs a frequency for.
process CHECK_NOT_ASSESSED {
    input:
    path rare, stageAs: 'rare/*'
    path empty, stageAs: 'empty/*'

    script:
    """
    python3 - <<'EOF'
import json
import os

import pyarrow.parquet as pq

sites = {}
for name in ('rare', 'empty'):
    folder = f'{name}/roh_upd'
    assert sorted(os.listdir(folder)) == ['roh.parquet', 'roh_upd.json', 'upd.parquet'], os.listdir(folder)
    record = json.load(open(f'{folder}/roh_upd.json'))
    assert record['upd_children'] == [], record
    assert record['upd_not_assessed'] == [
        {'child': 'T3', 'mother': 'T2', 'father': 'T1', 'reason': 'no informative site'}], record
    assert record['roh_members'] == ['T1', 'T2', 'T3'], record
    assert pq.read_table(f'{folder}/upd.parquet').num_rows == 0
    assert pq.read_table(f'{folder}/roh.parquet').num_rows == 0
    sites[name] = record['sites']
assert sites == {'rare': 74, 'empty': 0}, sites
EOF
    """
}
