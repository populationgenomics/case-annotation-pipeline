"""The reading of a gene model file into transcripts, on invented gene models.

Every feature line here is invented. The lines are in the layout of Ensembl's GFF3 (`ENSEMBL`) and of NCBI's RefSeq
GFF (`REFSEQ`): one gene of the forward strand and one of the reverse strand.
"""

from __future__ import annotations

import dataclasses
import pathlib

import invented_transcripts
import pysam
import pytest

from case_annotation import gene_models, transcript_models

NAMES = {'1': '1', 'chr1': 'chr1', 'INV_000001.1': 'INV_000001.1'}
ENSEMBL = """\
1\tinvented\tgene\t70\t580\t.\t+\t.\tID=gene:ENSG01;biotype=protein_coding;gene_id=ENSG01;version=1
1\tinvented\tmRNA\t70\t580\t.\t+\t.\tID=transcript:ENST01;Parent=gene:ENSG01;biotype=protein_coding;transcript_id=ENST01;version=1
1\tinvented\tfive_prime_UTR\t70\t180\t.\t+\t.\tParent=transcript:ENST01
1\tinvented\texon\t70\t180\t.\t+\t.\tParent=transcript:ENST01;Name=ENSE01;rank=1
1\tinvented\texon\t193\t260\t.\t+\t.\tParent=transcript:ENST01;Name=ENSE02;rank=2
1\tinvented\texon\t290\t480\t.\t+\t.\tParent=transcript:ENST01;Name=ENSE03;rank=3
1\tinvented\tCDS\t211\t260\t.\t+\t0\tID=CDS:ENSP01;Parent=transcript:ENST01
1\tinvented\tCDS\t290\t450\t.\t+\t1\tID=CDS:ENSP01;Parent=transcript:ENST01
1\tinvented\tgene\t700\t990\t.\t-\t.\tID=gene:ENSG02;biotype=protein_coding;gene_id=ENSG02;version=3
1\tinvented\tmRNA\t700\t990\t.\t-\t.\tID=transcript:ENST02;Parent=gene:ENSG02;biotype=protein_coding;transcript_id=ENST02;version=4
1\tinvented\texon\t900\t990\t.\t-\t.\tParent=transcript:ENST02;Name=ENSE04;rank=1
1\tinvented\texon\t700\t800\t.\t-\t.\tParent=transcript:ENST02;Name=ENSE05;rank=2
1\tinvented\tCDS\t900\t960\t.\t-\t0\tID=CDS:ENSP02;Parent=transcript:ENST02
1\tinvented\tCDS\t750\t800\t.\t-\t2\tID=CDS:ENSP02;Parent=transcript:ENST02
1\tinvented\tncRNA_gene\t1100\t1200\t.\t+\t.\tID=gene:ENSG03;biotype=lncRNA;gene_id=ENSG03;version=1
1\tinvented\tlnc_RNA\t1100\t1200\t.\t+\t.\tID=transcript:ENST03;Parent=gene:ENSG03;biotype=lncRNA;transcript_id=ENST03;version=1
1\tinvented\texon\t1100\t1200\t.\t+\t.\tParent=transcript:ENST03;Name=ENSE06;rank=1
"""
REFSEQ = """\
INV_000001.1\tinvented\tgene\t70\t480\t.\t+\t.\tID=gene-GENEA;Name=GENEA;gbkey=Gene
INV_000001.1\tinvented\tmRNA\t70\t480\t.\t+\t.\tID=rna-NM_01.1;Parent=gene-GENEA;gbkey=mRNA;transcript_id=NM_01.1
INV_000001.1\tinvented\texon\t70\t180\t.\t+\t.\tID=exon-NM_01.1-1;Parent=rna-NM_01.1
INV_000001.1\tinvented\texon\t193\t260\t.\t+\t.\tID=exon-NM_01.1-2;Parent=rna-NM_01.1
INV_000001.1\tinvented\texon\t290\t480\t.\t+\t.\tID=exon-NM_01.1-3;Parent=rna-NM_01.1
INV_000001.1\tinvented\tCDS\t211\t260\t.\t+\t0\tID=cds-NP_01.1;Parent=rna-NM_01.1
INV_000001.1\tinvented\tCDS\t290\t450\t.\t+\t1\tID=cds-NP_01.1;Parent=rna-NM_01.1
INV_000001.1\tinvented\tgene\t600\t690\t.\t+\t.\tID=gene-SEGMENT;Name=SEGMENT;gbkey=Gene
INV_000001.1\tinvented\tV_gene_segment\t600\t690\t.\t+\t.\tID=id-SEGMENT;Parent=gene-SEGMENT;gbkey=V_segment
INV_000001.1\tinvented\texon\t600\t690\t.\t+\t.\tID=exon-id-SEGMENT-1;Parent=id-SEGMENT
INV_000001.1\tinvented\tCDS\t610\t690\t.\t+\t0\tID=cds-SEGMENT;Parent=id-SEGMENT
"""


def _prepared(folder: pathlib.Path, text: str, name: str = 'models.gff') -> pathlib.Path:
    """The feature lines as a prepared gene model file: sorted, compressed and indexed."""
    lines = sorted(text.splitlines(), key=lambda line: (line.split('\t')[0], *map(int, line.split('\t')[3:5])))
    plain = folder / name
    plain.write_text('\n'.join(lines) + '\n')
    pysam.tabix_compress(str(plain), f'{plain}.gz', force=True)
    pysam.tabix_index(f'{plain}.gz', preset='gff', force=True)
    return pathlib.Path(f'{plain}.gz')


def _read(folder: pathlib.Path, text: str) -> dict[str, transcript_models.Transcript]:
    transcripts = transcript_models.read(_prepared(folder, text), 'chr1', NAMES)
    return {transcript.transcript_id: transcript for transcript in transcripts}


def _edited(text: str, old: str, new: str) -> str:
    assert old in text
    return text.replace(old, new)


def test_a_transcript_has_its_exons_in_its_own_order_and_its_coding_ends(tmp_path: pathlib.Path) -> None:
    read = _read(tmp_path, ENSEMBL)
    assert sorted(read) == ['ENST01.1', 'ENST02.4']
    forward, reverse = read['ENST01.1'], read['ENST02.4']
    assert (forward.chrom, forward.strand, forward.start, forward.end) == ('chr1', 1, 70, 580)
    assert forward.exons == ((70, 180), (193, 260), (290, 480))
    # 111 bases of the first exon and 18 of the second precede the coding sequence, which ends 30 bases early.
    assert forward.coding == transcript_models.Coding(cdna_start=130, cdna_end=340, phase=0, codon_table=1)
    assert forward.coding_region == (211, 450)
    assert reverse.exons == ((900, 990), (700, 800))
    assert reverse.coding == transcript_models.Coding(cdna_start=31, cdna_end=142, phase=0, codon_table=1)
    assert reverse.coding_region == (750, 960)
    assert (reverse.introns, forward.introns) == (((801, 899),), ((181, 192), (261, 289)))
    assert forward.has_frameshift_intron
    assert not reverse.has_frameshift_intron


def test_a_transcript_without_a_coding_sequence_is_not_read(tmp_path: pathlib.Path) -> None:
    assert 'ENST03.1' not in _read(tmp_path, ENSEMBL)


def test_ncbis_transcripts_take_their_ends_from_their_exons_and_any_type_can_code(tmp_path: pathlib.Path) -> None:
    read = _read(tmp_path, REFSEQ)
    assert sorted(read) == ['NM_01.1', 'id-SEGMENT']
    assert (read['NM_01.1'].start, read['NM_01.1'].end) == (70, 480)
    assert read['id-SEGMENT'].coding == transcript_models.Coding(cdna_start=11, cdna_end=91, phase=0, codon_table=1)


def test_a_transcript_of_a_file_that_ranks_its_exons_keeps_the_ends_of_its_own_line(tmp_path: pathlib.Path) -> None:
    wider = _edited(ENSEMBL, 'mRNA\t70\t580', 'mRNA\t60\t580')
    assert _read(tmp_path, wider)['ENST01.1'].start == 60
    narrower = _edited(REFSEQ, 'mRNA\t70\t480', 'mRNA\t60\t480')
    assert _read(tmp_path, narrower)['NM_01.1'].start == 70


def test_the_phase_of_the_first_coding_line_is_the_number_of_bases_that_the_first_codon_lacks(
    tmp_path: pathlib.Path,
) -> None:
    for gff_phase, lacking in (('1', 2), ('2', 1)):
        text = _edited(ENSEMBL, 'CDS\t211\t260\t.\t+\t0', f'CDS\t211\t260\t.\t+\t{gff_phase}')
        coding = _read(tmp_path, text)['ENST01.1'].coding
        assert coding is not None
        assert coding.phase == lacking


def test_a_coding_line_that_begins_in_no_exon_is_left_out(tmp_path: pathlib.Path) -> None:
    text = _edited(ENSEMBL, 'CDS\t290\t450', 'CDS\t285\t450')
    coding = _read(tmp_path, text)['ENST01.1'].coding
    assert coding is not None
    assert (coding.cdna_start, coding.cdna_end) == (130, 179)


def test_a_transcript_whose_unranked_exons_overlap_is_left_out(tmp_path: pathlib.Path) -> None:
    text = _edited(REFSEQ, 'exon\t193\t260', 'exon\t180\t260')
    assert sorted(_read(tmp_path, text)) == ['id-SEGMENT']


def test_a_transcript_whose_gene_is_of_a_type_that_is_not_read_is_left_out(tmp_path: pathlib.Path) -> None:
    text = _edited(ENSEMBL, '1\tinvented\tgene\t70\t580', '1\tinvented\tbiological_region\t70\t580')
    assert sorted(_read(tmp_path, text)) == ['ENST02.4']


def test_a_transcript_of_the_mitochondrial_genome_has_its_own_genetic_code(tmp_path: pathlib.Path) -> None:
    for name, text in (
        ('MT', ENSEMBL.replace('1\tinvented', 'MT\tinvented')),
        ('INV_M.1', REFSEQ.replace('INV_000001.1', 'INV_M.1')),
    ):
        names = {name: name, 'chrM': 'chrM'}
        first, *_others = transcript_models.read(_prepared(tmp_path, text), 'chrM', names)
        assert first.coding is not None
        assert first.coding.codon_table == 2


@pytest.mark.parametrize(
    ('old', 'new', 'message'),
    [
        (
            '1\tinvented\tCDS\t290\t450\t.\t+\t1\tID=CDS:ENSP01',
            '1\tinvented\tmRNA\t290\t450\t.\t+\t1\tID=other',
            'child',
        ),
        ('Name=ENSE02;rank=2', 'Name=ENSE02;rank=5', 'ranks'),
        ('Name=ENSE02;rank=2', 'Name=ENSE02', 'some of its exons'),
        ('1\tinvented\texon\t70\t180\t.\t+', '1\tinvented\texon\t70\t180\t.\t.', 'strand'),
        ('ID=gene:ENSG03;biotype=lncRNA;gene_id=ENSG03;version=1', 'ID=gene:ENSG03\tmore', 'columns'),
    ],
)
def test_a_file_that_vep_would_not_read_fails(tmp_path: pathlib.Path, old: str, new: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _read(tmp_path, _edited(ENSEMBL, old, new))


def test_a_file_that_names_a_chromosome_in_two_ways_fails(tmp_path: pathlib.Path) -> None:
    both = ENSEMBL + REFSEQ.replace('gene-', 'gene-other-').replace('rna-', 'rna-other-')
    with pytest.raises(ValueError, match='several ways'):
        transcript_models.read(_prepared(tmp_path, both), 'chr1', NAMES)


def test_a_chromosome_that_the_file_does_not_hold_has_no_transcripts(tmp_path: pathlib.Path) -> None:
    assert transcript_models.read(_prepared(tmp_path, ENSEMBL), 'chr2', {'2': '2', 'chr2': 'chr2'}) == []


def test_the_synonyms_give_every_chromosome_its_names(tmp_path: pathlib.Path) -> None:
    pairs = [
        (chrom.removeprefix('chr'), chrom, f'INV_{number}.1') for number, chrom in enumerate(gene_models.CHROMOSOMES)
    ]
    path = tmp_path / 'synonyms.txt'
    path.write_text(''.join(f'{a}\t{b}\n{a}\t{c}\n{b}\t{c}\n' for a, b, c in pairs))
    names = transcript_models.sequence_names(path)
    assert list(names) == list(gene_models.CHROMOSOMES)
    assert sorted(names['chrX']) == sorted(['X', 'chrX', 'INV_22.1'])
    path.write_text(''.join(f'{a}\t{b}\n' for a, b, _c in pairs[1:]))
    with pytest.raises(ValueError, match="no other name of \\['chr1'\\]"):
        transcript_models.sequence_names(path)
    path.write_text('1 chr1\n')
    with pytest.raises(ValueError, match='no pair of names'):
        transcript_models.sequence_names(path)


@pytest.mark.parametrize('strand', [1, -1])
def test_a_transcripts_sequences_are_read_on_its_strand(tmp_path: pathlib.Path, strand: int) -> None:
    invented = invented_transcripts.gene(strand)
    with pysam.FastaFile(str(invented_transcripts.reference(tmp_path, invented.genome.lower()))) as reference:
        assert transcript_models.sequences(invented.transcript, reference) == invented.sequences
    assert invented.sequences.utr5 == invented_transcripts.UTR5
    assert invented.sequences.coding == invented_transcripts.CODING
    assert invented.sequences.utr3 == invented_transcripts.UTR3


def test_a_coding_phase_pads_the_coding_sequence(tmp_path: pathlib.Path) -> None:
    invented = invented_transcripts.laid_out(('GCTGCTGCTTAAGG',), (1, 12), phase=2)
    with pysam.FastaFile(str(invented_transcripts.reference(tmp_path, invented.genome))) as reference:
        assert transcript_models.sequences(invented.transcript, reference).coding == 'NNGCTGCTGCTTAA'


def test_a_reference_that_does_not_hold_an_exon_fails(tmp_path: pathlib.Path) -> None:
    invented = invented_transcripts.gene()
    with (
        pysam.FastaFile(str(invented_transcripts.reference(tmp_path, invented.genome[:300]))) as reference,
        pytest.raises(ValueError, match='the reference holds'),
    ):
        transcript_models.sequences(invented.transcript, reference)


@pytest.mark.parametrize('strand', [1, -1])
def test_a_position_of_the_genome_and_one_of_the_transcript_give_each_other(strand: int) -> None:
    transcript = invented_transcripts.gene(strand).transcript
    for cdna in range(1, transcript.length + 1):
        assert transcript.cdna(transcript.genomic(cdna)) == cdna
    in_intron = transcript.introns[0][0]
    assert transcript.cdna(in_intron) is None
    with pytest.raises(ValueError, match='has no base'):
        transcript.genomic(transcript.length + 1)


def test_a_span_is_cut_into_its_pieces_inside_and_outside_the_exons() -> None:
    Segment = transcript_models.Segment  # noqa: N806
    forward = invented_transcripts.gene(1).transcript
    reverse = invented_transcripts.gene(-1).transcript
    # The first exon is 101 to 120 on the forward strand, and the second begins at 153.
    assert forward.segments(105, 110) == (Segment(5, 10, mapped=True),)
    assert forward.segments(118, 155) == (Segment(18, 20, True), Segment(121, 152, False), Segment(21, 23, True))
    assert forward.segments(99, 102) == (Segment(99, 100, False), Segment(1, 2, True))
    # On the reverse strand the first exon is 574 to 593, and the pieces come in the transcript's order.
    assert reverse.segments(572, 580) == (Segment(14, 20, True), Segment(572, 573, False))


def test_an_insertion_lies_between_its_two_neighbours_or_beside_the_one_in_an_exon() -> None:
    Segment = transcript_models.Segment  # noqa: N806
    forward = invented_transcripts.gene(1).transcript
    reverse = invented_transcripts.gene(-1).transcript
    assert forward.segments(106, 105) == (Segment(6, 5, True),)
    assert forward.segments(121, 120) == (Segment(21, 20, True),)  # behind the first exon's last base
    assert forward.segments(153, 152) == (Segment(21, 20, True),)  # before the second exon's first base
    assert forward.segments(101, 100) == (Segment(1, 0, True),)  # before the transcript's first base
    assert forward.segments(130, 129) == (Segment(130, 129, False),)  # in the intron
    assert reverse.segments(590, 589) == (Segment(5, 4, True),)
    assert reverse.segments(574, 573) == (Segment(21, 20, True),)  # behind the first exon's last base


def test_the_transcripts_that_hold_a_base_are_found_once_each() -> None:
    first = invented_transcripts.gene(1).transcript
    wide = invented_transcripts.laid_out(('A' * 40_000,), (1, 30), transcript_id='ENST91.1').transcript
    overlapping = transcript_models.Overlapping([first, wide])
    assert overlapping.of(150, 150) == [first, wide]
    assert overlapping.of(30_000, 39_000) == [wide]
    assert overlapping.of(50_000, 50_001) == []
    assert overlapping.of(90, 100) == []


@pytest.mark.parametrize('strand', [1, -1])
def test_the_exons_and_introns_near_a_span_are_found_with_a_reach(strand: int) -> None:
    transcript = invented_transcripts.gene(strand).transcript
    for start, end in transcript.exons:
        assert transcript.exons_near(end, end) == [(start, end)]
        assert transcript.exons_near(end + 1, end + 1) == []
        assert transcript.exons_near(end + 3, end + 5, reach=3) == [(start, end)]
        assert transcript.exons_near(start - 5, start - 4, reach=3) == []
    first, second = sorted(transcript.introns)[:2]
    assert transcript.introns_near(first[0] - 3, first[0] - 3, reach=3) == [first]
    assert transcript.introns_near(first[0] - 4, first[0] - 4, reach=3) == []
    assert sorted(transcript.introns_near(first[1], second[0])) == [first, second]
    assert sorted(transcript.exons_near(transcript.start, transcript.end)) == sorted(transcript.exons)


def test_exons_that_overlap_or_lie_in_one_another_are_all_found() -> None:
    """A file that ranks its exons can give exons that overlap, and a long exon can hide a later short one."""
    exons = ((100, 400), (150, 160), (170, 180), (390, 500), (600, 610))
    nested = invented_transcripts.gene().transcript
    nested = dataclasses.replace(nested, exons=exons, start=100, end=610, coding=None)
    for low in range(90, 620, 7):
        for high in (low, low + 3, low + 40):
            for reach in (0, 3, 12):
                met = [(start, end) for start, end in exons if start - reach <= high and low <= end + reach]
                assert sorted(nested.exons_near(low, high, reach)) == met, (low, high, reach)
