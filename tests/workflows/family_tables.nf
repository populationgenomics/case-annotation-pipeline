// Runs the family tables' processes on two invented families, wired as family_tables.nf wires them, and checks what
// they write: FAMILY_CHECK, NORMALISE_FAMILY and UNION_SITES; ANNOTATION_TABLES on an invented annotation of the union,
// which INVENTED_ANNOTATION writes in the Annotator's layout, with bundle.json and run.json standing in for the
// configs; GNOMAD_JOINT on an invented echtvar archive, which INVENTED_ARCHIVE encodes; and FAMILY_TABLES, which reads
// each family's PED file and proband as family_tables.nf passes them. The sheet's proband is the daughter in family A
// and empty in family B, whose affected son therefore has no role. Everything under family_tables/ is invented, the PED
// files too: family-a.vcf holds a record of two alleles that both members carry, one of two alleles of which one nobody
// carries, a spanning deletion beside an allele, a deletion written at the right end of a repeat, a record nobody
// carries, a record that only a call of ./1 carries, a deletion of six bases, a record off the primary chromosomes and
// a record of two alleles on X with genotype likelihoods, called in a haploid father and a diploid mother and daughter;
// family-b.vcf holds an allele family A shares, the same deletion left-aligned, and one of its own. gnomad-joint.vcf
// stands in for gnomAD's joint sites VCF, with invented counts and frequencies under gnomAD's field names, and
// echtvar-config.json encodes it as the public archive of docs/annotation.md is encoded: it holds one of the two
// alleles at a site, an allele counted 0 times, the deletion at its left-aligned position, the deletion of six bases
// (echtvar keeps alleles of more than four bases in a table of their own), an allele whose frequency is a usual
// threshold, a filtering allele frequency that is a multiple of 1/2000000 and that echtvar's 32-bit arithmetic stores
// one step below, and another allele at a position of the union's. reference.fa holds an N at chr2:300 and the IUPAC
// code R at chr2:320, which bcftools reads as N: ref-n.vcf, whose REF is N at both, is normalised, and a REF of N at a
// base (ref-n-at-base.vcf) and a base as REF at the N and at the R (ref-base-at-n.vcf, ref-base-at-iupac.vcf) are
// refused. Family A also brings the VCFs of its larger variants, which COMPRESSED gives the indexes a run's sheet names
// them with, and family B brings none. family-a.sv.vcf holds a deletion written as a sequence at another place of a
// repeat than family-a.vcf writes it, the deletion of six bases as family-a.vcf has it, a deletion that no small
// variant is, a symbolic deletion, a breakend and its mate, a record off the primary chromosomes and a symbolic
// deletion on X, where the father is called haploid. family-a.cnv.vcf holds a segment at the reference's copy number, a
// loss, a segment that the father loses and the mother gains, and a gain on X. A1, A2 and A3.repeats.vcf hold each
// member's genotypes of two invented repeats. The checks refuse repeat VCFs of two of the three members, a
// structural-variant VCF of another family's samples, and a structural variant whose REF is not the reference's
// (sv-ref-mismatch.vcf).
// The structural variants and copy-number segments go the way family_tables.nf sends them: SV_SITES cuts family A's
// sites, SV_UNION writes their union and its map, INVENTED_SV_ANNOTATION writes an invented annotation of the union in
// the Annotator's layout, SV_ANNOTATION_TABLES converts it, and FAMILY_TABLES reads family A's records' annotation from
// it. A deletion written as a sequence without INFO END (sv-no-end.vcf) is refused when the union is written.
// The ##source line of the VCFs of the larger variants names build/tests/family.nf, which means this file.
// Refused inputs must fail, and so must Annotator outputs without the group's annotation. The case-annotation and
// echtvar steps run in the tools image, built first. Run with Docker from the repository's root, with this folder's
// config alone (-C), because a run from the root would otherwise read the pipeline's nextflow.config over it:
//   docker build -f containers/tools/Dockerfile -t case-annotation-tools-test .
//   nextflow -C tests/workflows/nextflow.config run tests/workflows/family_tables.nf

include { FAMILY_CHECK ; NORMALISE_FAMILY ; NORMALISE_SV ; UNION_SITES ; SV_SITES ; SV_UNION ; GNOMAD_JOINT ; ANNOTATION_TABLES ; SV_ANNOTATION_TABLES ; FAMILY_TABLES ; groupAnnotation ; svName ; svGroup } from '../../modules/family_tables'
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

    checked = FAMILY_CHECK(families.join(structural))
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
    archive = INVENTED_ARCHIVE(file("${fixture}/gnomad-joint.vcf"), file("${fixture}/echtvar-config.json"))
    joint = GNOMAD_JOINT(union, archive)
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
    written = FAMILY_TABLES(
        normalised.family.join(families.map { f, _vcf, ped, proband -> [f, ped, proband] }).join(structural).join(sv_alleles, remainder: true).map { f, vcf, index, ped, proband, sv, cnv, repeats, alleles ->
            tuple(f, vcf, index, ped, proband, sv, cnv, repeats, alleles ?: [])
        },
        tables,
        joint,
        sv_tables,
        sv_union.map,
        '0.4.0',
    )
    CHECK_TABLES(written.filter { w -> w[0] == 'family-a' }.map { w -> w[0..3] }, archive)
    CHECK_NO_ROLES(written.filter { w -> w[0] == 'family-b' }.map { w -> w[3] })
    CHECK_LARGER(written.filter { w -> w[0] == 'family-a' }.map { w -> tuple(w[1], w[3], w[4]) })
    // A family whose row names no such VCF gets no table of them.
    written
        .filter { w -> w[0] == 'family-b' }
        .map { w ->
            if (w[4]) {
                error("family B brought no VCF of its larger variants, and has the tables ${w[4]}")
            }
        }

    // A REF of N is normalised where the reference has an N or an IUPAC code.
    ref_n = NORMALISE_REF_N(channel.of(tuple('ref-n', file("${fixture}/ref-n.vcf"))), reference)
    CHECK_REF_N(ref_n.family.map { _f, vcf, index -> tuple(vcf, index) })

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
    // has no INFO END, which the Annotator would annotate as a single base.
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
        channel.of(tuple('family-a', file("${fixture}/family-a.vcf"), file("${fixture}/family-a.ped"), 'B3', [], [], []))
    )
    two_repeats = FAMILY_CHECK_REPEATS_REFUSED(
        channel.of(tuple('family-a', file("${fixture}/family-a.vcf"), file("${fixture}/family-a.ped"), 'A3', [], [], [file("${fixture}/A1.repeats.vcf"), file("${fixture}/A2.repeats.vcf")]))
    )
    other_samples = FAMILY_CHECK_SV_SAMPLES_REFUSED(
        channel.of(tuple('family-b', file("${fixture}/family-b.vcf"), file("${fixture}/family-b.ped"), '', [file("${fixture}/family-a.sv.vcf")], [], []))
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
    unprefixed.sites
        .map { _s -> 'unprefixed contigs' }
        .mix(
            no_end_union.map.map { _m -> 'a deletion without END' },
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
    input:
    path vcf
    path config

    output:
    path 'invented.zip'

    script:
    """
    echtvar encode invented.zip ${config} ${vcf}
    """
}

// Family A's tables: its ten alleles, each with its invented transcript, its gene's HGNC id and, after the allele's
// columns of variants, a call column per member in the VCF's order, which says that the member carries the allele
// where its genotype holds a copy of it; the manifest gives the daughter, whom the sheet names, the proband's role and
// her parents theirs.
// Each split allele of the record on X has each member's likelihoods of its own genotypes, two for the haploid father
// and three for the others, and an allele of a record without PL has none. An allele the archive holds has the counts
// it stores, and each frequency as the 32-bit float echtvar returns: the multiple of 1/2000000 at or below the invented
// value, or the one below that for 0.062502, and 0 for a filtering allele frequency the invented VCF lacks. An allele
// the archive lacks, such as the other allele of a site or of a position it holds, has none; the manifest records the
// archive by its size and SHA-256, with its fields and echtvar's version.
process CHECK_TABLES {
    input:
    tuple val(family), path(variants), path(transcripts), path(manifest)
    path archive

    script:
    """
    python3 - <<'EOF'
import hashlib
import json
import struct

import pyarrow as pa
import pyarrow.parquet as pq


def float32(value):
    return struct.unpack('f', struct.pack('f', value))[0]


manifest = json.load(open('${manifest}'))
rows = {t['table']: t['rows'] for t in manifest['tables']}
assert {t: rows[t] for t in ('variants', 'transcripts')} == {'variants': 10, 'transcripts': 10}, rows
assert manifest['family'] == '${family}', manifest
roles = {m['sample_id']: m['role'] for m in manifest['members']}
assert roles == {'A1': 'father', 'A2': 'mother', 'A3': 'proband'}, roles
variants = pq.read_table('${variants}')
assert variants.column_names[:5] == ['chrom', 'pos', 'ref', 'alt', 'variant_id'], variants.column_names
assert variants.column_names[-5:] == ['qual', 'filters', 'A1', 'A2', 'A3'], variants.column_names
assert variants.num_rows == 10, variants.num_rows
pl = {(r['variant_id'], s): (r[s]['gt'], r[s]['pl']) for r in variants.to_pylist() for s in ('A1', 'A2', 'A3')}
carries = {key: r[key[1]]['carries'] for r in variants.to_pylist() for key in pl if key[0] == r['variant_id']}
assert all(carries[key] == ('1' in gt) for key, (gt, _) in pl.items()), carries
assert any(carries.values()) and not all(carries.values()), carries
assert pl[('chrX-100-T-A', 'A1')] == ('0', [90, 80]), pl
assert pl[('chrX-100-T-A', 'A2')] == ('0/1', [70, 0, 60]), pl
assert pl[('chrX-100-T-A', 'A3')] == ('1/0', [73, 62, 63]), pl
assert pl[('chrX-100-T-C', 'A1')] == ('1', [90, 0]), pl
assert pl[('chrX-100-T-C', 'A2')] == ('0/0', [70, 71, 72]), pl
assert pl[('chrX-100-T-C', 'A3')] == ('0/1', [73, 51, 52]), pl
assert pl[('chr2-60-C-A', 'A1')] == ('0/1', None), pl
assert variants.schema.field('hgnc_ids').type == pa.list_(pa.int32()), variants.schema
assert variants.column('hgnc_ids').to_pylist() == [[90001]] * 10, variants.column('hgnc_ids')
transcripts = pq.read_table('${transcripts}')
assert transcripts.schema.field('hgnc_id').type == pa.int32(), transcripts.schema
assert transcripts.column('hgnc_id').to_pylist() == [90001] * 10, transcripts.column('hgnc_id')
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
EOF
    """
}

// Family B's manifest: the sheet names no proband, so nobody has a role, although the PED file has an affected son with
// both parents.
process CHECK_NO_ROLES {
    input:
    path manifest

    script:
    """
    python3 - <<'EOF'
import json

members = json.load(open('${manifest}'))['members']
assert {m['sample_id']: m['role'] for m in members} == {'B1': None, 'B2': None, 'B3': None}, members
assert [m['sample_id'] for m in members if m['affected']] == ['B3'], members
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
// repeats, one for the father on X.
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
assert rows == {'variants': 10, 'transcripts': 10, 'sv': 10, 'sv_transcripts': 11, 'repeats': 2}, rows
assert manifest['sv_callers'] == ['sv', 'cnv'], manifest
assert manifest['schema_version'] == '0.4.0', manifest
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
assert repeats.column_names == ['chrom', 'pos', 'end_pos', 'repeat_id', 'repeat_unit', 'ref_copies', 'A1', 'A2',
                                'A3'], repeats.column_names
repeats = {r['repeat_id']: r for r in repeats.to_pylist()}
assert list(repeats) == ['RPTA', 'RPTX_CGG'], repeats
assert (repeats['RPTX_CGG']['repeat_unit'], repeats['RPTX_CGG']['ref_copies']) == ('CGG', 20), repeats
sizes = {s: repeats['RPTX_CGG'][s]['repcn'] for s in ('A1', 'A2', 'A3')}
assert sizes == {'A1': [23], 'A2': [23, 30], 'A3': [23, 23]}, sizes
assert repeats['RPTA']['A3']['repci_upper'] == [12, 58], repeats
assert repeats['RPTA']['A3']['so'] == ['SPANNING', 'INREPEAT'], repeats
EOF
    """
}
