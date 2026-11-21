"""The NMD plugin's rule, on invented transcripts: a gene of four exons on either strand, and one of one exon."""

from __future__ import annotations

import invented_transcripts
import pytest

from case_annotation import nmd

STRANDS = pytest.mark.parametrize('strand', [1, -1])


def _base(invented: invented_transcripts.Invented, cdna: int) -> int:
    return invented.transcript.genomic(cdna)


@STRANDS
def test_the_last_exon_is_an_escape_with_its_untranslated_part(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    first_of_last = gene.transcript.exon_cdna_starts[-1]
    for cdna in (first_of_last, first_of_last + 50, gene.transcript.length):
        assert nmd.escape_reason(gene.transcript, _base(gene, cdna), 500) == 'E-1'
    assert nmd.escape_reason(gene.transcript, _base(gene, first_of_last - 60), 500) is None


@STRANDS
def test_the_last_52_bases_of_the_exon_before_the_last_are_an_escape(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    last_of_penultimate = gene.transcript.exon_cdna_starts[-1] - 1
    assert nmd.escape_reason(gene.transcript, _base(gene, last_of_penultimate), 500) == 'E-2_end-50'
    assert nmd.escape_reason(gene.transcript, _base(gene, last_of_penultimate - 51), 500) == 'E-2_end-50'
    assert nmd.escape_reason(gene.transcript, _base(gene, last_of_penultimate - 52), 500) is None


@STRANDS
def test_the_first_101_coding_bases_are_an_escape(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    middle = _base(gene, 150)
    assert nmd.escape_reason(gene.transcript, middle, 101) == 'cds_start+100'
    assert nmd.escape_reason(gene.transcript, middle, 102) is None


def test_a_transcript_of_one_exon_is_an_escape_before_any_other_reason() -> None:
    single = invented_transcripts.laid_out(('GGCACC' + 'ATGGCTGCTTAA' + 'GG',), (7, 18))
    assert nmd.escape_reason(single.transcript, _base(single, 9), 3) == 'intronless'


@STRANDS
def test_the_reasons_come_in_the_plugin_s_order(strand: int) -> None:
    gene = invented_transcripts.gene(strand)
    in_last_exon = _base(gene, gene.transcript.length - 5)
    in_reach = _base(gene, gene.transcript.exon_cdna_starts[-1] - 5)
    assert nmd.escape_reason(gene.transcript, in_last_exon, 50) == 'E-1'
    assert nmd.escape_reason(gene.transcript, in_reach, 50) == 'E-2_end-50'
    assert set(nmd.REASONS) == {'intronless', 'E-1', 'E-2_end-50', 'cds_start+100'}
