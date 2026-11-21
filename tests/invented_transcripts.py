"""Invented transcripts for the tests of the plugins' rules: a small genome laid out around designed exons.

Nothing here is a real gene. `GENE` has four exons. Its 5' UTR is the first exon and the first five bases of the
second, with one upstream ORF, and its coding sequence of 360 bases has a `CAG` codon, which one substitution turns
into a stop codon, in each stretch that the NMD rule tells apart.
"""

from __future__ import annotations

import dataclasses
import pathlib
from collections.abc import Sequence

import pysam

from case_annotation import plugin_gates, transcript_models

CHROM = 'chr1'
_FIRST = 101

# The codons of the coding sequence, counted from 0, that read CAG: in the first 101 coding bases, in the middle of
# the third exon, in that exon's last 52 bases and in the last exon.
CAG_CODONS = {'first': 10, 'middle': 50, 'penultimate': 80, 'last': 100}
_CODONS = ['ATG'] + ['GCT'] * 118 + ['TAA']
for _codon in CAG_CODONS.values():
    _CODONS[_codon] = 'CAG'
CODING = ''.join(_CODONS)
UTR5 = 'GGCACCATGGCTCCGTAAGC' + 'GCGCC'
UTR3 = 'GGCCTTGGAACC'
_INTRON = 'GTAAGT' + 'C' * 24 + 'AG'
EXONS = (UTR5[:20], UTR5[20:] + CODING[:26], CODING[26:276], CODING[276:] + UTR3)


@dataclasses.dataclass(frozen=True)
class Invented:
    """A transcript with its sequences and the genome that holds it."""

    transcript: transcript_models.Transcript
    sequences: transcript_models.Sequences
    genome: str

    def allele(self, cdna: int, removed: int, inserted: str) -> plugin_gates.Allele:
        """The allele that replaces `removed` bases from a position of the spliced transcript with `inserted`.

        The position counts from 1, and an insertion stands before its base. The allele is given as VEP holds it.
        """
        transcript = self.transcript
        if transcript.strand == 1:
            start = transcript.genomic(cdna)
            alt = inserted
        else:
            start = transcript.genomic(cdna + removed - 1) if removed else transcript.genomic(cdna) + 1
            alt = transcript_models.reverse_complement(inserted)
        return plugin_gates.Allele(
            start=start, end=start + removed - 1, ref=self.genome[start - 1 : start - 1 + removed], alt=alt
        )

    def coding_allele(self, position: int, removed: int, inserted: str) -> plugin_gates.Allele:
        """The same for a position of the coding sequence, counted from 1."""
        coding = self.transcript.coding
        assert coding is not None
        return self.allele(coding.cdna_start + position - 1, removed, inserted)


def laid_out(
    exons: Sequence[str],
    coding: tuple[int, int],
    strand: int = 1,
    phase: int = 0,
    introns: Sequence[str] | None = None,
    transcript_id: str = 'ENST90.1',
) -> Invented:
    """A transcript of the given exons, on either strand of a genome that is made for it.

    `coding` gives the first and the last coding base along the spliced transcript, counted from 1.
    """
    between = list(introns) if introns is not None else [_INTRON] * (len(exons) - 1)
    layout, places, offset = '', [], 0
    for index, exon in enumerate(exons):
        places.append((offset, offset + len(exon) - 1))
        layout += exon + (between[index] if index < len(between) else '')
        offset = len(layout)
    last = _FIRST + len(layout) - 1
    coordinates = tuple(
        (_FIRST + start, _FIRST + end) if strand == 1 else (last - end, last - start) for start, end in places
    )
    genome = 'N' * (_FIRST - 1) + (layout if strand == 1 else transcript_models.reverse_complement(layout)) + 'N' * 60
    spliced = ''.join(exons)
    return Invented(
        transcript=transcript_models.Transcript(
            transcript_id=transcript_id,
            chrom=CHROM,
            strand=strand,
            start=_FIRST,
            end=last,
            exons=coordinates,
            coding=transcript_models.Coding(cdna_start=coding[0], cdna_end=coding[1], phase=phase, codon_table=1),
        ),
        sequences=transcript_models.Sequences(
            spliced=spliced,
            utr5=spliced[: coding[0] - 1],
            coding='N' * phase + spliced[coding[0] - 1 : coding[1]],
            utr3=spliced[coding[1] :],
        ),
        genome=genome,
    )


def gene(strand: int = 1) -> Invented:
    return laid_out(EXONS, (len(UTR5) + 1, len(UTR5) + len(CODING)), strand=strand)


def reference(folder: pathlib.Path, genome: str) -> pathlib.Path:
    """A FASTA file of the genome as its one sequence, with its index."""
    path = folder / 'invented.fa'
    path.write_text(f'>{CHROM}\n{genome}\n')
    pysam.faidx(str(path))
    return path
