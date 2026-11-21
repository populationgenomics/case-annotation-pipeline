"""The table of the plugins' rules: one allele on one transcript, and a whole run on the invented gene models.

The transcripts of invented_transcripts and the gene models of tests/workflows/family_tables/gene-models/ are
invented, and so is every allele here.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import json
import pathlib
import re
from collections.abc import Sequence
from typing import Any

import invented_transcripts
import pyarrow.parquet as pq
import pydantic
import pysam
import pytest

from case_annotation import cli, gene_models, nmd, plugin_gates, transcript_rules, utr5

FIXTURE = pathlib.Path(__file__).parent / 'workflows' / 'family_tables'
REFERENCE = FIXTURE / 'reference.fa'
STRANDS = pytest.mark.parametrize('strand', [1, -1])
UTR = len(invented_transcripts.UTR5)
COUNTS = {'utr5_existing_uorfs': 1, 'utr5_existing_oorfs_out_of_frame': 0, 'utr5_existing_oorfs_in_frame': 0}
OWN_START = {'utr5_cds_kozak_context': 'GCCATGG', 'utr5_cds_kozak_strength': 'Strong'}


def _annotated(tmp_path: pathlib.Path, invented: invented_transcripts.Invented) -> transcript_rules.Annotated:
    with pysam.FastaFile(str(invented_transcripts.reference(tmp_path, invented.genome))) as reference:
        return transcript_rules.Annotated(invented.transcript, reference)


def _columns(
    tmp_path: pathlib.Path, invented: invented_transcripts.Invented, allele: plugin_gates.Allele
) -> dict[str, Any] | None:
    return transcript_rules.transcript_columns(_annotated(tmp_path, invented), allele, frozenset())


@STRANDS
def test_an_allele_in_the_5_prime_utr_has_the_counts_and_the_own_start_codon(
    tmp_path: pathlib.Path, strand: int
) -> None:
    gene = invented_transcripts.gene(strand)
    assert _columns(tmp_path, gene, gene.allele(3, 1, 'T')) == COUNTS | OWN_START


@STRANDS
def test_an_allele_with_an_effect_has_its_records_and_their_kinds(tmp_path: pathlib.Path, strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    columns = _columns(tmp_path, gene, gene.allele(8, 1, 'C'))
    assert columns is not None
    assert columns['utr5_consequences'] == ['uAUG_lost']
    (record,) = columns['utr5_uorfs']
    assert (record['consequence'], record['cap_distance_to_start'], record['distance_to_cds']) == ('uAUG_lost', 6, 19)
    assert set(record) == {'consequence'} | transcript_rules.KINDS['uAUG_lost']
    assert 'nmd_escaping' not in columns


@STRANDS
def test_a_record_has_the_evidence_of_its_start_codon_s_place_on_the_genome(
    tmp_path: pathlib.Path, strand: int
) -> None:
    gene = invented_transcripts.gene(strand)
    seen = frozenset({gene.transcript.genomic(7)})
    columns = transcript_rules.transcript_columns(_annotated(tmp_path, gene), gene.allele(8, 1, 'C'), seen)
    assert columns is not None
    assert columns['utr5_uorfs'][0]['evidence'] is True


@STRANDS
def test_the_records_of_a_row_are_sorted_by_kind_and_from_the_cap(tmp_path: pathlib.Path, strand: int) -> None:
    two_frames = invented_transcripts.laid_out(('GGCACCATGGCCATGGCTCCGTAAGC' + 'ATGGCTTAA',), (27, 35), strand=strand)
    columns = _columns(tmp_path, two_frames, two_frames.allele(18, 0, 'ATGC'))
    assert columns is not None
    assert columns['utr5_consequences'] == ['uAUG_gained', 'uFrameShift']
    kinds = [(record['consequence'], record.get('ref_start_distance_to_cds')) for record in columns['utr5_uorfs']]
    assert kinds == [('uAUG_gained', None), ('uFrameShift', 20), ('uFrameShift', 14)]


@STRANDS
def test_an_allele_that_no_rule_evaluates_has_no_row(tmp_path: pathlib.Path, strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    intron_start, _end = gene.transcript.introns[1]
    in_intron = plugin_gates.Allele(intron_start + 12, intron_start + 12, 'C', 'A')
    assert _columns(tmp_path, gene, in_intron) is None
    assert _columns(tmp_path, gene, gene.coding_allele(40, 1, 'A')) is None  # a missense change
    beyond = plugin_gates.Allele(gene.transcript.end + 20, gene.transcript.end + 20, 'N', 'A')
    assert _columns(tmp_path, gene, beyond) is None


@STRANDS
def test_an_allele_that_ends_the_reading_frame_has_the_prediction_and_its_reason(
    tmp_path: pathlib.Path, strand: int
) -> None:
    gene = invented_transcripts.gene(strand)
    expected = {'first': 'cds_start+100', 'middle': None, 'penultimate': 'E-2_end-50', 'last': 'E-1'}
    for where, codon in invented_transcripts.CAG_CODONS.items():
        columns = _columns(tmp_path, gene, gene.coding_allele(codon * 3 + 1, 1, 'T'))
        assert columns == {'nmd_escaping': expected[where] is not None, 'nmd_escape_reason': expected[where]}
    shifted = _columns(tmp_path, gene, gene.coding_allele(150, 0, 'G'))
    assert shifted == {'nmd_escaping': False, 'nmd_escape_reason': None}
    assert expected.keys() == invented_transcripts.CAG_CODONS.keys()
    assert {reason for reason in expected.values() if reason} <= set(nmd.REASONS)


@STRANDS
def test_a_splice_site_allele_without_both_coding_ends_is_not_evaluated(tmp_path: pathlib.Path, strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    start, end = gene.transcript.introns[1]
    at_site = plugin_gates.Allele(start, start, gene.genome[start - 1], 'A' if gene.genome[start - 1] != 'A' else 'C')
    assert _columns(tmp_path, gene, at_site) is None
    whole = plugin_gates.Allele(start - 2, end + 2, gene.genome[start - 3 : end + 2], '')
    assert _columns(tmp_path, gene, whole) == {'nmd_escaping': True, 'nmd_escape_reason': 'cds_start+100'}


@STRANDS
def test_an_allele_across_the_end_of_the_utr_has_the_counts_and_no_record(tmp_path: pathlib.Path, strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    assert _columns(tmp_path, gene, gene.allele(UTR - 1, 4, '')) == COUNTS | OWN_START
    # An insertion directly before the first coding base is in the UTR for VEP, and the plugin's rules do not read
    # it. It changes the context of the own start codon.
    assert _columns(tmp_path, gene, gene.allele(UTR + 1, 0, 'T')) == COUNTS | OWN_START | {
        'utr5_cds_kozak_context_alt': 'CCTATGG',
        'utr5_cds_kozak_strength_alt': 'Moderate',
    }


@STRANDS
def test_a_deletion_from_one_exon_of_the_utr_into_the_next_has_no_record(tmp_path: pathlib.Path, strand: int) -> None:
    """The plugin takes as many bases from the spliced UTR as the deletion has on the genome, the intron's included."""
    gene = invented_transcripts.gene(strand)
    start, end = gene.transcript.introns[0]
    across = plugin_gates.Allele(start - 3, end + 2, gene.genome[start - 4 : end + 2], '')
    assert transcript_rules.utr_change(gene.transcript, across, gene.sequences.utr5) is None
    assert _columns(tmp_path, gene, across) == COUNTS | OWN_START


@STRANDS
def test_an_allele_as_a_change_of_the_utr_is_on_the_transcript_s_strand(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    assert transcript_rules.utr_change(gene.transcript, gene.allele(8, 1, 'C'), gene.sequences.utr5) == utr5.Change(
        7, 1, 'C'
    )
    assert transcript_rules.utr_change(gene.transcript, gene.allele(8, 0, 'CA'), gene.sequences.utr5) == utr5.Change(
        7, 0, 'CA'
    )
    assert transcript_rules.utr_change(gene.transcript, gene.allele(8, 3, ''), gene.sequences.utr5) == utr5.Change(
        7, 3, ''
    )
    # The UTR's last five bases are in the second exon.
    assert transcript_rules.utr_change(gene.transcript, gene.allele(22, 2, 'TT'), gene.sequences.utr5) == utr5.Change(
        21, 2, 'TT'
    )
    assert transcript_rules.utr_change(gene.transcript, gene.coding_allele(5, 1, 'A'), gene.sequences.utr5) is None


@STRANDS
def test_an_insertion_before_the_first_coding_base_is_read_where_it_can_be_written_in_the_utr(
    tmp_path: pathlib.Path, strand: int
) -> None:
    """The plugin reads an insertion by its following base, which is a coding one here, and gives no record.

    `C` after the UTR's last base, a `C`, gives the sequence that `C` before that base gives, which the plugin reads.
    """
    frame = invented_transcripts.laid_out(('GGCACCATGGCC' + 'ATGGCTGAGCTGTAAGTAGC',), (13, 32), strand=strand)
    utr = frame.sequences.utr5
    at_the_end, inside = frame.allele(13, 0, 'C'), frame.allele(12, 0, 'C')
    assert transcript_rules.utr_change(frame.transcript, at_the_end, utr) == utr5.Change(12, 0, 'C')
    assert transcript_rules.utr_change(frame.transcript, inside, utr) == utr5.Change(11, 0, 'C')
    written = [_columns(tmp_path, frame, allele) for allele in (at_the_end, inside)]
    assert written[0] is not None
    assert written[0] == written[1]
    assert written[0]['utr5_consequences'] == ['uFrameShift']
    # An insertion there that ends with another base can only be written before the coding base.
    assert transcript_rules.utr_change(frame.transcript, frame.allele(13, 0, 'CT'), utr) is None
    assert transcript_rules.utr_change(frame.transcript, frame.allele(14, 0, 'C'), utr) is None


@STRANDS
def test_a_change_of_the_own_start_codon_s_context_gives_the_context_with_the_allele(
    tmp_path: pathlib.Path, strand: int
) -> None:
    gene = invented_transcripts.gene(strand)
    before = _columns(tmp_path, gene, gene.allele(UTR - 2, 1, 'T'))
    assert before == COUNTS | OWN_START | {
        'utr5_cds_kozak_context_alt': 'TCCATGG',
        'utr5_cds_kozak_strength_alt': 'Moderate',
    }
    # The base behind the start codon is a coding one: the row has the own start codon's columns and no count.
    behind = _columns(tmp_path, gene, gene.coding_allele(4, 1, 'C'))
    assert behind == OWN_START | {'utr5_cds_kozak_context_alt': 'GCCATGC', 'utr5_cds_kozak_strength_alt': 'Moderate'}
    assert _columns(tmp_path, gene, gene.coding_allele(2, 1, 'C')) is None  # in the start codon


def test_a_transcript_without_a_whole_start_codon_has_no_context_of_one(tmp_path: pathlib.Path) -> None:
    other_start = invented_transcripts.laid_out(('GGCACCGCC' + 'CTGGCTTAA' + 'GG',), (10, 18))
    columns = _columns(tmp_path, other_start, other_start.allele(3, 1, 'T'))
    assert columns is not None
    assert (columns['utr5_cds_kozak_context'], columns['utr5_cds_kozak_strength']) == (None, None)
    assert _annotated(tmp_path, other_start).own_start is None
    padded = invented_transcripts.laid_out(('GGCACCGCC' + 'ATGGCTTAA' + 'GG',), (10, 18), phase=1)
    assert _annotated(tmp_path, padded).own_start is None


# A whole run


def _gzipped(path: pathlib.Path, text: str) -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='utf-8') as f:
        f.write(text)
    return path


def _prepared(tmp_path: pathlib.Path) -> pathlib.Path:
    models = FIXTURE / 'gene-models'
    public: dict[gene_models.Source, pathlib.Path] = {
        'Ensembl': _gzipped(tmp_path / 'invented.gff3.gz', (models / 'ensembl.gff3').read_text('utf-8')),
        'RefSeq': _gzipped(tmp_path / 'invented.gff.gz', (models / 'refseq.gff').read_text('utf-8')),
    }
    gene_models.prepare(public, models / 'assembly-report.txt', REFERENCE, tmp_path / 'gene-models')
    return tmp_path / 'gene-models'


ALLELES = (
    ('chr1', 100, 'G', 'A'),
    ('chr1', 100, 'G', 'C'),
    ('chr1', 200, 'G', 'A'),
    ('chr1', 300, 'T', 'A'),
    ('chr1', 400, 'TCA', 'T'),
)


def _given(
    tmp_path: pathlib.Path,
    prepared: pathlib.Path,
    alleles: Sequence[tuple[str, int, str, str]] = ALLELES,
    controls: bool = True,
) -> pathlib.Path:
    """The alleles, each under its variant id, and the controls of the prepared gene models, sorted and indexed."""
    records = [(chrom, pos, f'{chrom}-{pos}-{ref}-{alt}', ref, alt) for chrom, pos, ref, alt in alleles]
    if controls:
        for line in (prepared / gene_models.CONTROLS_NAME).read_text().splitlines():
            if not line.startswith('#'):
                chrom, pos, name, ref, alt = line.split('\t')[:5]
                records.append((chrom, int(pos), name, ref, alt))
    plain = tmp_path / 'given.vcf'
    head = '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=600>\n##contig=<ID=chr2,length=600>\n'
    lines = [f'{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t.\t.\t.\n' for chrom, pos, name, ref, alt in sorted(records)]
    plain.write_text(head + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' + ''.join(lines))
    pysam.tabix_index(str(plain), preset='vcf', force=True)
    return pathlib.Path(f'{plain}.gz')


def _evidence(
    tmp_path: pathlib.Path, text: str = 'chr\tpos\tgene\tstrand\ttype\tstop_pos\nchr1\t82\tGENEA\t+\tuORF\t0\n'
) -> pathlib.Path:
    path = tmp_path / 'invented-uorfs.txt'
    path.write_text(text)
    return path


@pytest.mark.parametrize('workers', [1, 2])
def test_a_run_writes_a_row_for_each_allele_and_transcript_that_a_rule_evaluated(
    tmp_path: pathlib.Path, workers: int
) -> None:
    prepared = _prepared(tmp_path)
    output = tmp_path / 'out'
    info = transcript_rules.annotate(
        _given(tmp_path, prepared), prepared, REFERENCE, _evidence(tmp_path), output, workers
    )
    table = pq.read_table(output / transcript_rules.TABLE_NAME)
    assert table.schema.equals(transcript_rules.TABLE)
    rows = table.to_pylist()
    evaluated = ['chr1-100-G-A', 'chr1-100-G-C', 'chr1-200-G-A', 'chr1-400-TCA-T']
    keys = [(row['variant_id'], row['transcript_id'], row['source']) for row in rows]
    assert keys == [(a, t, s) for a in evaluated for t, s in (('ENST01.1', 'Ensembl'), ('NM_01.1', 'RefSeq'))]
    for row in rows:
        in_utr = row['variant_id'] != 'chr1-400-TCA-T'
        counts = (
            row['utr5_existing_uorfs'],
            row['utr5_existing_oorfs_out_of_frame'],
            row['utr5_existing_oorfs_in_frame'],
        )
        assert counts == ((0, 0, 1) if in_utr else (None, None, None))
        assert (row['nmd_escaping'], row['nmd_escape_reason']) == ((None, None) if in_utr else (False, None))
        assert (row['utr5_cds_kozak_context'] is not None) is in_utr
        if row['variant_id'] != 'chr1-200-G-A':
            assert (row['utr5_consequences'], row['utr5_uorfs']) == (None, None)
            continue
        assert row['utr5_consequences'] == ['uSTOP_gained']
        (record,) = row['utr5_uorfs']
        assert {field: value for field, value in record.items() if value is not None} == {
            'consequence': 'uSTOP_gained',
            'ref_type': 'InFrame_oORF',
            'kozak_context': 'GTGATGC',
            'kozak_strength': 'Moderate',
            'ref_start_distance_to_cds': 117,
            'new_stop_distance_to_cds': 12,
            'evidence': True,
        }
    assert info == transcript_rules.read_info(output)
    assert (info.utrannotator_commit, info.nmd_plugin_commit) == (utr5.SOURCE_COMMIT, nmd.SOURCE_COMMIT)
    assert info.vep_release == transcript_rules.VEP_RELEASE
    assert info.uorf_evidence.name == 'invented-uorfs.txt'
    assert [(m.source, m.file.name, m.prepared_from.name) for m in info.gene_models] == [
        ('Ensembl', 'ensembl.gff3.gz', 'invented.gff3.gz'),
        ('RefSeq', 'refseq.gff.gz', 'invented.gff.gz'),
    ]
    written = json.loads((output / transcript_rules.INFO_NAME).read_text())
    assert (info.branch_points, written['branch_points']) == (None, None)
    assert sorted(written) == [
        'branch_points',
        'gene_models',
        'nmd_plugin_commit',
        'uorf_evidence',
        'utrannotator_commit',
        'vep_release',
    ]


# The invented catalogue of the workflow's test gives the third intron of the invented gene a branch point at chr1:486,
# 24 bases before the fourth exon, and the intron holds an `AG` behind it. A `G` at chr1:502 makes a second one.
CATALOGUE = FIXTURE / 'branch_points.tsv'
GAINED = {
    'branch_point_pos': 486,
    'branch_point_distance': 16,
    'acceptor_distance': -7,
    'yag': True,
    'high_risk': True,
    'existing_ag': True,
    'score': 4,
}
BESIDE_EXONS = (*ALLELES, ('chr1', 496, 'C', 'G'), ('chr1', 502, 'T', 'G'), ('chr1', 506, 'A', 'ACAGAG'))


@pytest.mark.parametrize('workers', [1, 2])
def test_a_run_with_the_catalogue_of_branch_points_writes_the_ags_that_an_allele_gains(
    tmp_path: pathlib.Path, workers: int
) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared, alleles=BESIDE_EXONS)
    output = tmp_path / 'out'
    info = transcript_rules.annotate(
        given, prepared, REFERENCE, _evidence(tmp_path), output, workers, branch_points=CATALOGUE
    )
    table = pq.read_table(output / transcript_rules.TABLE_NAME)
    assert table.schema.equals(transcript_rules.TABLE)
    gained = {
        (row['variant_id'], row['transcript_id'], row['source']): row[transcript_rules.AG_GAINS.name]
        for row in table.to_pylist()
        if row[transcript_rules.AG_GAINS.name] is not None
    }
    # The substitution gains one `AG` on the gene's transcript of each file, and the insertion two. The third allele
    # of the intron makes no `AG` and has no row.
    transcripts = (('ENST01.1', 'Ensembl'), ('NM_01.1', 'RefSeq'))
    assert set(gained) == {
        (allele, *transcript) for allele in ('chr1-502-T-G', 'chr1-506-A-ACAGAG') for transcript in transcripts
    }
    for transcript in transcripts:
        assert gained[('chr1-502-T-G', *transcript)] == [GAINED]
        two = gained[('chr1-506-A-ACAGAG', *transcript)]
        assert [(gain['branch_point_distance'], gain['acceptor_distance']) for gain in two] == [(23, -5), (25, -3)]
    assert all(row['variant_id'] != 'chr1-496-C-G' for row in table.to_pylist())
    # The rows of the two plugins' rules are those of a run without the catalogue.
    transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'without', workers)
    without = pq.read_table(tmp_path / 'without' / transcript_rules.TABLE_NAME)
    assert without.column(transcript_rules.AG_GAINS.name).null_count == without.num_rows
    assert table.filter(table.column(transcript_rules.AG_GAINS.name).is_null()).equals(without)
    # The record names the catalogue as the run read it.
    assert info.branch_points == transcript_rules.FileRecord(
        name=CATALOGUE.name, bytes=CATALOGUE.stat().st_size, sha256=hashlib.sha256(CATALOGUE.read_bytes()).hexdigest()
    )
    assert transcript_rules.read_info(tmp_path / 'without').branch_points is None


def test_the_plugins_rules_read_the_alleles_of_the_exon_regions_with_or_without_the_catalogue(
    tmp_path: pathlib.Path,
) -> None:
    """The gene's first intron has 12 bases and lies in its 5' UTR, where VEP's rules count the UTR's reading frames.

    The cut to the exons leaves an allele inside that intron out, so a run without the catalogue never reads it. A
    run with the catalogue is given it, and reads it for the gained `AG` alone.
    """
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared, alleles=(('chr1', 100, 'G', 'A'), ('chr1', 185, 'CC', 'AG')))
    lines = CATALOGUE.read_text().splitlines()
    row = 'chr1\t{pos}\t+\tcBP_chr1_{pos}_+_G\t\t\t\t\t\t\t\t\tcBP_Invented\tcBP_chr1_{pos}_+_G'
    rows = {}
    # A branch point ten bases before the second exon, in the intron, and one 13 bases before it, in the first exon,
    # which is no branch point of the intron.
    for where, pos in (('inside', 183), ('outside', 180)):
        catalogue = tmp_path / f'{where}.tsv'
        catalogue.write_text('\n'.join([*lines[:3], row.format(pos=pos), *lines[3:]]) + '\n')
        output = tmp_path / where
        transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), output, branch_points=catalogue)
        rows[where] = pq.read_table(output / transcript_rules.TABLE_NAME).to_pylist()
    transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'without')
    without = pq.read_table(tmp_path / 'without' / transcript_rules.TABLE_NAME).to_pylist()
    assert [row['variant_id'] for row in without] == ['chr1-100-G-A'] * 2
    assert rows['outside'] == without
    in_intron = [row for row in rows['inside'] if row['variant_id'] == 'chr1-185-CC-AG']
    assert [row for row in rows['inside'] if row not in in_intron] == without
    assert [(row['transcript_id'], row['source']) for row in in_intron] == [
        ('ENST01.1', 'Ensembl'),
        ('NM_01.1', 'RefSeq'),
    ]
    values = {'branch_point_pos': 183, 'branch_point_distance': 3, 'acceptor_distance': -6}
    for row in in_intron:
        assert row[transcript_rules.AG_GAINS.name] == [GAINED | values | {'high_risk': False, 'existing_ag': False}]
        assert all(row[column] is None for column in transcript_rules.COLUMNS.names[:-1])


# An invented gene of the reverse strand on chr2, whose last exon lies inside a run of thirteen `CT` from chr2:91 to
# chr2:116, with a `G` on either side. Its intron ends at chr2:107, and a branch point lies at chr2:130.
REVERSE_GENE = (
    '2\tinvented\tgene\t101\t300\t.\t-\t.\tID=gene:ENSG02;Name=GENEB;biotype=protein_coding;gene_id=ENSG02;version=1\n'
    '2\tinvented\tmRNA\t101\t300\t.\t-\t.\tID=transcript:ENST02;Parent=gene:ENSG02;Name=GENEB-201;'
    'biotype=protein_coding;transcript_id=ENST02;version=1\n'
    '2\tinvented\texon\t200\t300\t.\t-\t.\tParent=transcript:ENST02;Name=ENSE05;exon_id=ENSE05;rank=1;version=1\n'
    '2\tinvented\texon\t101\t106\t.\t-\t.\tParent=transcript:ENST02;Name=ENSE06;exon_id=ENSE06;rank=2;version=1\n'
    '2\tinvented\tCDS\t210\t290\t.\t-\t0\tID=CDS:ENSP02;Parent=transcript:ENST02;protein_id=ENSP02\n'
)


def test_an_insertion_that_a_record_writes_before_a_transcript_gains_the_ag_of_its_writing_in_the_intron(
    tmp_path: pathlib.Path,
) -> None:
    """One more `CT` of the run is one more `AG` on the gene's strand, behind the branch point of its last intron.

    The record writes it at the run's first base on the genome, which lies before the gene, and no base of the record
    is the gene's.
    """
    with pysam.FastaFile(str(REFERENCE)) as fixture:
        sequences = {name: fixture.fetch(name) for name in fixture.references}
    sequences['chr2'] = sequences['chr2'][:89] + 'G' + 'CT' * 13 + 'G' + sequences['chr2'][117:]
    reference = tmp_path / 'with-a-run.fa'
    reference.write_text(''.join(f'>{name}\n{bases}\n' for name, bases in sequences.items()))
    pysam.faidx(str(reference))
    models = FIXTURE / 'gene-models'
    public: dict[gene_models.Source, pathlib.Path] = {
        'Ensembl': _gzipped(tmp_path / 'invented.gff3.gz', (models / 'ensembl.gff3').read_text('utf-8') + REVERSE_GENE),
        'RefSeq': _gzipped(tmp_path / 'invented.gff.gz', (models / 'refseq.gff').read_text('utf-8')),
    }
    prepared = tmp_path / 'gene-models'
    gene_models.prepare(public, models / 'assembly-report.txt', reference, prepared)
    catalogue = tmp_path / 'catalogue.tsv'
    row = 'chr2\t130\t-\tcBP_chr2_130_-_A\t\t\t\t\t\t\t\t\tcBP_Invented\tcBP_chr2_130_-_A'
    catalogue.write_text('\n'.join([*CATALOGUE.read_text().splitlines()[:3], row]) + '\n')
    gained = GAINED | {'branch_point_pos': 130, 'branch_point_distance': 15, 'acceptor_distance': -10}
    # As the record writes it, before the gene, and as the same change is written inside the intron.
    for record in (('chr2', 90, 'G', 'GCT'), ('chr2', 114, 'T', 'TCT')):
        given = _given(tmp_path, prepared, alleles=(record,))
        output = tmp_path / f'at-{record[1]}'
        transcript_rules.annotate(given, prepared, reference, _evidence(tmp_path), output, branch_points=catalogue)
        rows = pq.read_table(output / transcript_rules.TABLE_NAME).to_pylist()
        assert [(row['transcript_id'], row[transcript_rules.AG_GAINS.name]) for row in rows] == [('ENST02.1', [gained])]
        assert all(rows[0][column] is None for column in transcript_rules.COLUMNS.names[:-1])


def test_a_record_is_in_the_exon_regions_by_the_bases_of_its_ref() -> None:
    exons = [(68, 181), (191, 261)]
    starts = [start for start, _end in exons]
    in_exons = {
        (pos, ref) for pos, ref in ((69, 'A'), (181, 'A'), (67, 'ACG'), (182, 'ACGTACGTACG'), (192, 'AC'), (261, 'A'))
    }
    outside = {(68, 'A'), (67, 'AC'), (182, 'A'), (182, 'ACGTACGTAC'), (262, 'A'), (5, 'A')}
    for pos, ref in in_exons | outside:
        assert transcript_rules._in_exons(starts, exons, pos, ref) is ((pos, ref) in in_exons), (pos, ref)


EXONS_REFUSED = {
    'a line without an end': 'chr1\t68\n',
    'a start that is no number': 'chr1\tx\t181\n',
    'an end at its start': 'chr1\t181\t181\n',
    'a region before the one above it': 'chr1\t191\t261\nchr1\t68\t181\n',
    'a region that meets the one above it': 'chr1\t68\t181\nchr1\t181\t261\n',
}


@pytest.mark.parametrize('case', EXONS_REFUSED)
def test_regions_of_the_exons_that_are_not_as_prepared_are_refused(tmp_path: pathlib.Path, case: str) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared)
    (prepared / gene_models.EXONS_NAME).write_text(EXONS_REFUSED[case])
    with pytest.raises(ValueError, match=r'exons\.bed: line [12] is'):
        transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out')


def test_a_catalogue_that_is_not_the_prepared_one_is_refused(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    other = tmp_path / 'other.tsv'
    other.write_text('chrom\tpos\nchr1\t486\n')
    with pytest.raises(ValueError, match='its header is not that of the prepared catalogue of branch points'):
        transcript_rules.annotate(
            _given(tmp_path, prepared), prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out', branch_points=other
        )
    assert not (tmp_path / 'out').exists()


def test_the_command_takes_the_catalogue_of_branch_points(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared, alleles=BESIDE_EXONS)
    arguments = ['--given', str(given), '--gene-models', str(prepared), '--reference', str(REFERENCE)]
    arguments += ['--uorf-evidence', str(_evidence(tmp_path)), '--output', str(tmp_path / 'out')]
    cli.main(['transcript-rules', *arguments, '--branch-points', str(CATALOGUE)])
    table = pq.read_table(tmp_path / 'out' / transcript_rules.TABLE_NAME)
    assert table.num_rows - table.column(transcript_rules.AG_GAINS.name).null_count == 4


def test_two_runs_of_the_same_inputs_write_the_same_table(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    given, evidence = _given(tmp_path, prepared), _evidence(tmp_path)
    transcript_rules.annotate(given, prepared, REFERENCE, evidence, tmp_path / 'one')
    transcript_rules.annotate(given, prepared, REFERENCE, evidence, tmp_path / 'two', workers=2)
    one, two = (tmp_path / name / transcript_rules.TABLE_NAME for name in ('one', 'two'))
    assert pq.read_table(one).equals(pq.read_table(two))


def test_a_union_without_an_allele_in_an_exon_writes_an_empty_table(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    transcript_rules.annotate(
        _given(tmp_path, prepared, alleles=()), prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out'
    )
    table = pq.read_table(tmp_path / 'out' / transcript_rules.TABLE_NAME)
    assert (table.num_rows, table.schema.equals(transcript_rules.TABLE)) == (0, True)


def test_sites_without_the_controls_are_refused(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared, controls=False)
    with pytest.raises(ValueError, match=r'positive controls .* are not among its records'):
        transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out')


def test_a_record_under_another_id_than_its_allele_s_is_refused(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared)
    with gzip.open(given, 'rt') as lines:
        text = lines.read().replace('chr1-200-G-A', 'chr1-200-G-T')
    (tmp_path / 'given.vcf').write_text(text)
    pysam.tabix_index(str(tmp_path / 'given.vcf'), preset='vcf', force=True)
    with pytest.raises(ValueError, match='not the variant id of its allele'):
        transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out')


def test_a_repeated_id_is_refused(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared, alleles=(*ALLELES, ALLELES[0]))
    with pytest.raises(ValueError, match='the ID of another record'):
        transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out')


def test_controls_on_a_transcript_that_the_gene_models_lack_are_refused(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    record = prepared / gene_models.RECORD_NAME
    record.write_text(record.read_text().replace('"ENST01.1"', '"ENST02.1"'))
    with pytest.raises(ValueError, match=r'control-UTRannotator-Ensembl has 0 rows on the transcript ENST02\.1'):
        transcript_rules.annotate(
            _given(tmp_path, prepared), prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out'
        )


def test_a_control_without_its_value_is_refused() -> None:
    control = gene_models.Control(
        id='control-NMD-Ensembl', plugin='NMD', transcript_id='ENST01.1', chrom='chr1', pos=5, ref='A', alt='AC'
    )
    row = {'transcript_id': 'ENST01.1', 'nmd_escaping': False}
    with pytest.raises(ValueError, match='lacks what the rules of NMD give it'):
        transcript_rules.check_control(control, [row], 'invented')
    transcript_rules.check_control(control, [row | {'nmd_escaping': True}], 'invented')
    gained = control.model_copy(update={'plugin': 'UTRannotator'})
    with pytest.raises(ValueError, match='lacks what the rules of UTRannotator give it'):
        transcript_rules.check_control(gained, [row | {'utr5_consequences': ['uAUG_lost']}], 'invented')
    transcript_rules.check_control(gained, [row | {'utr5_consequences': ['uAUG_gained']}], 'invented')


def test_a_gene_model_file_that_is_not_the_recorded_one_is_refused(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared)
    with (prepared / gene_models.FILES['RefSeq']).open('ab') as f:
        f.write(b'\0')
    with pytest.raises(ValueError, match=r'not the\s+file that gene-models\.json describes'):
        transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out')


def test_a_file_that_is_not_the_plugin_s_list_of_upstream_orfs_is_refused(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    evidence = _evidence(tmp_path, 'chr1\t82\n')
    with pytest.raises(ValueError, match='not with its header'):
        transcript_rules.annotate(_given(tmp_path, prepared), prepared, REFERENCE, evidence, tmp_path / 'out')


def test_a_record_of_two_alleles_is_refused(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared)
    with gzip.open(given, 'rt') as lines:
        text = re.sub(r'(control-NMD-Ensembl\t\w+\t\w+)', r'\1,T', lines.read())
    assert ',T\t' in text
    (tmp_path / 'given.vcf').write_text(text)
    pysam.tabix_index(str(tmp_path / 'given.vcf'), preset='vcf', force=True)
    with pytest.raises(ValueError, match='is no record of one allele under its id'):
        transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out')


def test_a_transcript_that_both_gene_model_files_give_a_row_is_refused(tmp_path: pathlib.Path) -> None:
    models = FIXTURE / 'gene-models'
    refseq = (models / 'refseq.gff').read_text('utf-8').replace('transcript_id=NM_01.1', 'transcript_id=ENST01.1')
    public: dict[gene_models.Source, pathlib.Path] = {
        'Ensembl': _gzipped(tmp_path / 'invented.gff3.gz', (models / 'ensembl.gff3').read_text('utf-8')),
        'RefSeq': _gzipped(tmp_path / 'invented.gff.gz', refseq),
    }
    prepared = tmp_path / 'gene-models'
    gene_models.prepare(public, models / 'assembly-report.txt', REFERENCE, prepared)
    with pytest.raises(ValueError, match=r'two rows for chr1-100-G-A and ENST01\.1'):
        transcript_rules.annotate(
            _given(tmp_path, prepared), prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out'
        )


def test_an_allele_with_a_base_that_is_not_one_of_four_has_no_row(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    given = _given(tmp_path, prepared, alleles=(('chr1', 100, 'G', 'N'), ('chr1', 200, 'G', '*')))
    transcript_rules.annotate(given, prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out')
    assert pq.read_table(tmp_path / 'out' / transcript_rules.TABLE_NAME).num_rows == 0


def test_a_transcript_without_a_coding_sequence_is_not_annotated(tmp_path: pathlib.Path) -> None:
    gene = invented_transcripts.gene()
    non_coding = dataclasses.replace(gene.transcript, coding=None)
    with (
        pysam.FastaFile(str(invented_transcripts.reference(tmp_path, gene.genome))) as reference,
        pytest.raises(ValueError, match='has no coding sequence'),
    ):
        transcript_rules.Annotated(non_coding, reference)
    with pytest.raises(ValueError, match='has no coding sequence'):
        plugin_gates.Placed(non_coding, gene.allele(3, 1, 'T'), gene.sequences)
    assert transcript_rules.utr_change(non_coding, gene.allele(3, 1, 'T'), gene.sequences.utr5) is None


def test_the_command_needs_at_least_one_process(tmp_path: pathlib.Path) -> None:
    prepared = _prepared(tmp_path)
    arguments = [
        '--given',
        str(_given(tmp_path, prepared)),
        '--gene-models',
        str(prepared),
        '--reference',
        str(REFERENCE),
    ]
    arguments += ['--uorf-evidence', str(_evidence(tmp_path)), '--output', str(tmp_path / 'out')]
    with pytest.raises(ValueError, match='at least one process'):
        cli.main(['transcript-rules', *arguments, '--workers', '0'])
    cli.main(['transcript-rules', *arguments])
    assert pq.read_table(tmp_path / 'out' / transcript_rules.TABLE_NAME).num_rows == 8


@pytest.mark.parametrize(
    'field',
    [
        {'utrannotator_commit': '3d9d466'},
        {'nmd_plugin_commit': 'release/116'},
        {'vep_release': ''},
        {'uorf_evidence': {'name': 'invented-uorfs.txt', 'bytes': 1, 'sha256': 'e3b0'}},
        {'vep_version': '116.3'},
    ],
)
def test_a_record_with_a_malformed_value_or_another_key_is_refused(field: dict[str, Any]) -> None:
    record = {
        'utrannotator_commit': '0' * 40,
        'nmd_plugin_commit': '1' * 40,
        'vep_release': '116.3',
        'uorf_evidence': {'name': 'invented-uorfs.txt', 'bytes': 1, 'sha256': '2' * 64},
        'gene_models': [],
        'branch_points': None,
    }
    transcript_rules.TranscriptRulesInfo.model_validate(record)
    with pytest.raises(pydantic.ValidationError):
        transcript_rules.TranscriptRulesInfo.model_validate(record | field)


def test_a_run_writes_its_progress_to_the_log(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(transcript_rules, 'PROGRESS_EVERY', 2)
    prepared = _prepared(tmp_path)
    transcript_rules.annotate(_given(tmp_path, prepared), prepared, REFERENCE, _evidence(tmp_path), tmp_path / 'out')
    lines = capsys.readouterr().out.splitlines()
    # Nine alleles on chr1, the five of the union and the four controls, read with each gene model file.
    of_chr1 = [line for line in lines if line.startswith('chr1, ')]
    assert [line.split(' after ')[0] for line in of_chr1] == [
        f'chr1, {source}: {done} of 9 alleles' for source in ('Ensembl', 'RefSeq') for done in (2, 4, 6, 8, 9)
    ]
    done = [line for line in lines if ' chromosomes done, ' in line]
    assert len(done) == len(gene_models.CHROMOSOMES)
    assert done[-1].startswith('25 of 25 chromosomes done, 9 alleles, after ')
