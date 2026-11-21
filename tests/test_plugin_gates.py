"""Which alleles VEP gives to its two plugins, on the invented gene of invented_transcripts, on either strand."""

from __future__ import annotations

import invented_transcripts
import pytest

from case_annotation import plugin_gates, transcript_models

STRANDS = pytest.mark.parametrize('strand', [1, -1])
UTR = len(invented_transcripts.UTR5)
CODING = len(invented_transcripts.CODING)
# Positions along the spliced transcript, counted from 1: the first exon is 1 to 20, the second 21 to 51 with the
# coding sequence from 26, the third 52 to 301 and the fourth from 302.


def placed(invented: invented_transcripts.Invented, allele: plugin_gates.Allele) -> plugin_gates.Placed:
    return plugin_gates.Placed(invented.transcript, allele, invented.sequences)


def intron_allele(
    invented: invented_transcripts.Invented, intron: int, offset: int, removed: int, inserted: str
) -> plugin_gates.Allele:
    """An allele at a base of an intron, counted from 1 along the transcript; an insertion stands before the base."""
    transcript = invented.transcript
    start, end = transcript.introns[intron]
    # On the reverse strand the allele's lowest base is its last on the transcript.
    backwards = end - offset + 1 - (removed - 1 if removed else -1)
    first = start + offset - 1 if transcript.strand == 1 else backwards
    ref = invented.genome[first - 1 : first - 1 + removed]
    alt = inserted if transcript.strand == 1 else transcript_models.reverse_complement(inserted)
    return plugin_gates.Allele(start=first, end=first + removed - 1, ref=ref, alt=alt)


def test_an_insertion_or_deletion_is_held_without_the_base_before_it() -> None:
    of_record = plugin_gates.Allele.of_record
    assert of_record(100, 'G', 'A') == plugin_gates.Allele(100, 100, 'G', 'A')
    assert of_record(100, 'G', 'GAT') == plugin_gates.Allele(101, 100, '', 'AT')
    assert of_record(100, 'GAT', 'G') == plugin_gates.Allele(101, 102, 'AT', '')
    assert of_record(100, 'GAT', 'GC') == plugin_gates.Allele(101, 102, 'AT', 'C')
    # A substitution of several bases keeps a shared first base, and so does a change whose first bases differ.
    assert of_record(100, 'GA', 'GC') == plugin_gates.Allele(100, 101, 'GA', 'GC')
    assert of_record(100, 'CC', 'A') == plugin_gates.Allele(100, 101, 'CC', 'A')


def test_the_stretches_where_two_alleles_differ_are_told_apart() -> None:
    assert plugin_gates.Allele(100, 103, 'ACGT', 'ATGA').differing == ((101, 101), (103, 103))
    assert plugin_gates.Allele(100, 103, 'ACGT', 'AG').differing == ((101, 103),)
    assert plugin_gates.Allele(100, 100, 'A', 'TG').differing == ((100, 100),)
    assert plugin_gates.Allele(101, 100, '', 'TG').differing == ((101, 100),)


def test_codons_are_translated_by_the_standard_code_or_the_mitochondrial_one() -> None:
    assert plugin_gates.translate('ATGGCTTAA', 1) == 'MA*'
    assert plugin_gates.translate('TGAAGAATA', 1) == '*RI'
    assert plugin_gates.translate('TGAAGAATA', 2) == 'W*M'
    # Every reading of CTN is leucine. TAN is tyrosine or a stop.
    assert plugin_gates.translate('CTNTANNNA', 1) == 'LXX'
    assert plugin_gates.translate('ATGGC', 1) == 'M'


def test_a_transcripts_protein_ends_before_its_stop_codon_and_begins_with_methionine() -> None:
    sequences = transcript_models.Sequences(spliced='', utr5='', coding='CTGGCTTAA', utr3='')
    assert plugin_gates.protein(sequences) == 'MA'
    without_stop = transcript_models.Sequences(spliced='', utr5='', coding='GCTGCTGC', utr3='')
    assert plugin_gates.protein(without_stop) == 'AA'
    # By the mitochondrial code ATA begins a protein, AGA ends it and TGA is tryptophan.
    mitochondrial = transcript_models.Sequences(spliced='', utr5='', coding='ATATGAAGA', utr3='')
    assert plugin_gates.protein(mitochondrial, 2) == 'MW'
    assert plugin_gates.protein(mitochondrial) == 'I*R'


@STRANDS
def test_an_allele_in_an_exon_before_the_coding_sequence_is_in_the_5_prime_utr(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    for cdna in (1, 10, 20, 21, UTR):
        assert placed(
            gene, gene.allele(cdna, 1, 'A' if gene.sequences.spliced[cdna - 1] != 'A' else 'C')
        ).five_prime_utr_variant
    assert placed(gene, gene.allele(10, 0, 'TT')).five_prime_utr_variant
    assert placed(gene, gene.allele(10, 3, '')).five_prime_utr_variant


@STRANDS
def test_an_allele_in_the_coding_sequence_or_an_intron_is_not_in_the_5_prime_utr(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    assert not placed(gene, gene.coding_allele(1, 1, 'C')).five_prime_utr_variant
    assert not placed(gene, gene.coding_allele(50, 1, 'A')).five_prime_utr_variant
    assert not placed(gene, intron_allele(gene, 0, 10, 1, 'A')).five_prime_utr_variant
    assert not placed(gene, gene.allele(UTR + CODING + 3, 1, 'A')).five_prime_utr_variant


@STRANDS
def test_the_ends_of_the_5_prime_utr_are_vep_s(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    # An insertion directly before the first coding base is in the UTR, and one behind that base is not.
    assert placed(gene, gene.allele(UTR + 1, 0, 'T')).five_prime_utr_variant
    assert not placed(gene, gene.allele(UTR + 2, 0, 'T')).five_prime_utr_variant
    # A deletion across the end of the UTR is in it. An insertion before the transcript's first base is not.
    assert placed(gene, gene.allele(UTR - 1, 4, '')).five_prime_utr_variant
    assert not placed(gene, gene.allele(1, 0, 'T')).five_prime_utr_variant
    # An insertion between an exon of the UTR and its intron is in the UTR.
    assert placed(gene, intron_allele(gene, 0, 1, 0, 'T')).five_prime_utr_variant


@STRANDS
def test_a_deletion_of_the_whole_transcript_has_none_of_the_terms(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    transcript = gene.transcript
    whole = plugin_gates.Allele(
        transcript.start - 2, transcript.end + 2, gene.genome[transcript.start - 3 : transcript.end + 2], ''
    )
    found = placed(gene, whole)
    terms = (found.five_prime_utr_variant, found.stop_gained, found.frameshift_variant, found.splice_site_variant)
    assert terms == (False, False, False, False)


@STRANDS
def test_a_substitution_that_makes_a_stop_codon_is_a_gained_stop(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    for codon in invented_transcripts.CAG_CODONS.values():
        found = placed(gene, gene.coding_allele(codon * 3 + 1, 1, 'T'))
        assert (found.stop_gained, found.frameshift_variant) == (True, False)
        assert (found.cds_start, found.cds_end) == (codon * 3 + 1, codon * 3 + 1)
    assert not placed(gene, gene.coding_allele(6, 1, 'C')).stop_gained  # GCT to GCC
    assert not placed(gene, gene.coding_allele(4, 1, 'A')).stop_gained  # GCT to ACT


@STRANDS
def test_a_change_of_the_stop_codon_is_no_gained_stop(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    assert not placed(gene, gene.coding_allele(CODING, 1, 'G')).stop_gained  # TAA to TAG
    assert not placed(gene, gene.coding_allele(CODING - 2, 1, 'C')).stop_gained  # TAA to CAA


@STRANDS
def test_an_insertion_or_deletion_in_the_coding_sequence_is_a_frameshift_unless_of_whole_codons(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    for removed, inserted, shifted in (
        (0, 'G', True),
        (0, 'GG', True),
        (0, 'GGC', False),
        (1, '', True),
        (3, '', False),
    ):
        found = placed(gene, gene.coding_allele(40, removed, inserted))
        assert found.frameshift_variant is shifted
    insertion = placed(gene, gene.coding_allele(40, 0, 'G'))
    assert (insertion.cds_start, insertion.cds_end) == (40, 39)
    # A substitution of two bases is none, and neither is an insertion in the UTR.
    assert not placed(gene, gene.coding_allele(40, 2, 'AA')).frameshift_variant
    assert not placed(gene, gene.allele(10, 0, 'G')).frameshift_variant


@STRANDS
def test_an_insertion_that_makes_a_stop_codon_in_frame_is_a_gained_stop(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    found = placed(gene, gene.coding_allele(40, 0, 'TAA'))
    assert (found.stop_gained, found.frameshift_variant) == (True, False)


@STRANDS
def test_an_allele_that_runs_out_of_the_coding_sequence_has_one_coding_end_and_shifts_no_frame(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    # The second exon ends at the 26th coding base. The deletion takes that base and the intron's first.
    transcript = gene.transcript
    low = min(transcript.genomic(51), transcript.introns[1][0 if strand == 1 else 1])
    across = plugin_gates.Allele(low, low + 1, gene.genome[low - 1 : low + 1], '')
    found = placed(gene, across)
    assert (found.cds_start, found.cds_end) == (26, None)
    assert (found.frameshift_variant, found.splice_site_variant) == (False, True)
    out_of_utr = placed(gene, gene.allele(UTR - 1, 4, ''))
    assert (out_of_utr.cds_start, out_of_utr.cds_end) == (None, 2)
    assert not out_of_utr.frameshift_variant


@STRANDS
def test_the_first_and_last_two_bases_of_an_intron_are_its_splice_sites(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    length = gene.transcript.introns[1][1] - gene.transcript.introns[1][0] + 1
    for offset, at_site in ((1, True), (2, True), (3, False), (length - 2, False), (length - 1, True), (length, True)):
        found = placed(gene, intron_allele(gene, 1, offset, 1, 'A'))
        assert found.splice_site_variant is at_site
        assert (found.cds_start, found.cds_end) == (None, None)
    assert not placed(gene, gene.allele(51, 1, 'T')).splice_site_variant  # the exon's last base
    # An insertion between the first two bases of the intron is at the site, and one before the first is not.
    assert placed(gene, intron_allele(gene, 1, 2, 0, 'T')).splice_site_variant
    assert not placed(gene, intron_allele(gene, 1, 1, 0, 'T')).splice_site_variant


@STRANDS
def test_a_deletion_of_a_whole_intron_is_at_its_splice_sites_and_has_both_coding_ends(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    start, end = gene.transcript.introns[1]
    whole = plugin_gates.Allele(start - 2, end + 2, gene.genome[start - 3 : end + 2], '')
    found = placed(gene, whole)
    assert found.splice_site_variant
    assert (found.cds_start, found.cds_end) == (25, 28)


def test_an_intron_of_a_few_bases_has_no_splice_sites() -> None:
    short = invented_transcripts.laid_out(('GGCACC' + 'ATGGCTGCT', 'GCTGCTGCTTAAGG'), (7, 27), introns=('GTCCCCCCAG',))
    start, _end = short.transcript.introns[0]
    assert short.transcript.has_frameshift_intron
    assert not placed(short, plugin_gates.Allele(start, start, 'G', 'A')).splice_site_variant


def test_an_allele_in_a_last_codon_that_the_coding_sequence_does_not_finish_shifts_no_frame() -> None:
    unfinished = invented_transcripts.laid_out(('GGCACC' + 'ATGGCTGCTGC',), (7, 17))
    assert placed(unfinished, unfinished.coding_allele(5, 1, '')).frameshift_variant
    assert not placed(unfinished, unfinished.coding_allele(10, 1, '')).frameshift_variant


def test_a_coding_phase_moves_the_coding_positions() -> None:
    padded = invented_transcripts.laid_out(('GCTCAGGCTGCTTAAGG',), (1, 15), phase=1)
    found = placed(padded, padded.allele(4, 1, 'G'))
    # With one `N` before it the sequence reads NGC TCA GGC, so the fourth base is the second of the second codon,
    # and a `G` there makes TGA.
    assert (found.cds_start, found.cds_end, found.stop_gained) == (5, 5, True)
    whole = invented_transcripts.laid_out(('GCTCAGGCTGCTTAAGG',), (1, 15))
    found = placed(whole, whole.allele(4, 1, 'G'))
    assert (found.cds_start, found.cds_end, found.stop_gained) == (4, 4, False)
