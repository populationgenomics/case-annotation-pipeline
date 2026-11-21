"""Which alleles Ensembl VEP gives to its UTRannotator and NMD plugins, decided as VEP decides it (docs/annotation.md).

A plugin of VEP runs on an allele and a transcript only when VEP has given the pair one of the plugin's consequence
terms. UTRannotator runs for `5_prime_UTR_variant`. The NMD plugin runs for `stop_gained`, `frameshift_variant`,
`splice_donor_variant` and `splice_acceptor_variant`, and then only when both ends of the allele have a position in
the coding sequence. The tables tell an allele that a plugin evaluated and did not flag from one that it never saw,
so this module decides the same five terms, and no other, in the way VEP 116.3 does. It is a restatement in Python of
the predicates in `Bio::EnsEMBL::Variation::Utils::VariationEffect`, `BaseVariationFeatureOverlapAllele`,
`BaseTranscriptVariation`, `BaseTranscriptVariationAllele` and `TranscriptVariationAllele` of ensembl-variation, as
the VEP image of release 116.3 holds them (Apache 2.0, licenses/ENSEMBL-NOTICE.txt).

Two things are not VEP's.

- VEP reads the 3' UTR of a reverse-strand transcript of a GFF file from the wrong bases, the span mirrored about the
  transcript, and with them decides whether an insertion or deletion at a stop codon keeps the stop. This module reads
  the transcript's own 3' UTR on both strands.
- VEP translates a mitochondrial transcript by the vertebrate mitochondrial code when the genome's name is `MT`. With
  a reference that names it `chrM` it uses the standard code, in which `TGA` ends a protein and `AGA` does not. Here
  a transcript of `chrM` is translated by the mitochondrial code (transcript_models).
"""

from __future__ import annotations

import dataclasses
import functools
import itertools
import re

from case_annotation import transcript_models

_UNAMBIGUOUS = re.compile(r'[ACGT]*')
_BASES = 'TCAG'
# The amino acids of the 64 codons in the order TTT, TTC, TTA, TTG, TCT, ..., GGG: NCBI's translation tables 1 (the
# standard code) and 2 (the vertebrate mitochondrial code).
_TABLES = {
    1: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    2: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
}
_CODONS = [''.join(codon) for codon in itertools.product(_BASES, repeat=3)]
_AMINO_ACIDS = {table: dict(zip(_CODONS, acids, strict=True)) for table, acids in _TABLES.items()}
# The start codons of each code. A protein's first amino acid is methionine where the coding sequence begins with one.
_STARTS = {1: frozenset({'TTG', 'CTG', 'ATG'}), 2: frozenset({'ATT', 'ATC', 'ATA', 'ATG', 'GTG'})}


@dataclasses.dataclass(frozen=True)
class Allele:
    """An allele as VEP holds it: on the forward strand, an insertion or deletion without the base before it.

    An insertion has an empty `ref` and an end one below its start, and lies between the bases `end` and `start`.
    """

    start: int
    end: int
    ref: str
    alt: str

    @classmethod
    def of_record(cls, pos: int, ref: str, alt: str) -> Allele:
        """The allele of a VCF record with one alternate allele.

        VEP drops the first base of both alleles when their lengths differ and they share it, and no other base.
        """
        if len(ref) != len(alt) and ref[:1] == alt[:1]:
            return cls(start=pos + 1, end=pos + len(ref) - 1, ref=ref[1:], alt=alt[1:])
        return cls(start=pos, end=pos + len(ref) - 1, ref=ref, alt=alt)

    @property
    def low(self) -> int:
        return min(self.start, self.end)

    @property
    def high(self) -> int:
        return max(self.start, self.end)

    @functools.cached_property
    def differing(self) -> tuple[tuple[int, int], ...]:
        """The stretches of the genome where the two alleles differ, each as a start and an end.

        VEP compares the alleles base by base when both have more than one, so a substitution of several bases with
        an unchanged base in it gives several stretches. Any other allele gives its whole span.
        """
        if len(self.ref) <= 1 or len(self.alt) <= 1:
            return ((self.start, self.end),)
        differing = [index for index, (ref, alt) in enumerate(itertools.zip_longest(self.ref, self.alt)) if ref != alt]
        runs: list[tuple[int, int]] = []
        for index in differing:
            if runs and runs[-1][1] == self.start + index - 1:
                runs[-1] = (runs[-1][0], self.start + index)
            else:
                runs.append((self.start + index, self.start + index))
        return tuple(runs)


def _overlap(start: int, end: int, other_start: int, other_end: int) -> bool:
    return end >= other_start and start <= other_end


def translate(sequence: str, table: int) -> str:
    """The amino acids of a sequence of whole codons, `*` for a stop codon.

    A codon with an `N` gives the amino acid that every reading of it gives, or `X` when the readings differ.
    """
    acids = _AMINO_ACIDS[table]
    peptide = []
    for index in range(0, len(sequence) - 2, 3):
        codon = sequence[index : index + 3]
        acid = acids.get(codon)
        if acid is None:
            readings = {
                acids.get(''.join(reading), 'X')
                for reading in itertools.product(*(_BASES if base == 'N' else base for base in codon))
            }
            acid = readings.pop() if len(readings) == 1 else 'X'
        peptide.append(acid)
    return ''.join(peptide)


def protein(sequences: transcript_models.Sequences, table: int = 1) -> str:
    """A transcript's protein as VEP holds it beside the transcript: without its last stop codon."""
    return _protein(sequences.coding, table)


@functools.lru_cache(maxsize=2048)
def _protein(coding: str, table: int) -> str:
    """The protein of a coding sequence, kept for the alleles that follow on the same transcript."""
    coding = coding[: len(coding) - len(coding) % 3]
    peptide = translate(coding, table)
    if coding and _AMINO_ACIDS[table].get(coding[-3:]) == '*':
        peptide = peptide[:-1]
    if peptide and peptide[0] != 'M' and coding[:3] in _STARTS[table]:
        peptide = 'M' + peptide[1:]
    return peptide


class Placed:
    """An allele on a transcript with a coding sequence: where it lies, and the consequence terms that VEP gives it.

    Each term is decided as VEP's predicate of the same name decides it, behind the conditions that VEP checks before
    it calls a predicate.
    """

    def __init__(
        self, transcript: transcript_models.Transcript, allele: Allele, sequences: transcript_models.Sequences
    ) -> None:
        if transcript.coding is None or transcript.coding_region is None:
            raise ValueError(f'{transcript.transcript_id} has no coding sequence')
        self.transcript = transcript
        self.allele = allele
        self.sequences = sequences
        self._coding = transcript.coding
        self._coding_region = transcript.coding_region
        self._alt = allele.alt if transcript.strand == 1 else transcript_models.reverse_complement(allele.alt)
        self._insertion = allele.start == allele.end + 1

    # Where the allele lies

    @functools.cached_property
    def within_transcript(self) -> bool:
        return _overlap(self.allele.start, self.allele.end, self.transcript.start, self.transcript.end)

    @functools.cached_property
    def _ablates(self) -> bool:
        """Whether the allele deletes the whole transcript, for which VEP gives `transcript_ablation` and no more."""
        allele, transcript = self.allele, self.transcript
        covers = allele.start <= transcript.start and allele.end >= transcript.end
        return covers and len(allele.ref) > len(allele.alt)

    @functools.cached_property
    def _in_exon(self) -> bool:
        stretch = 12 if self.transcript.has_frameshift_intron else 0
        low, high = self.allele.low, self.allele.high
        return bool(self.transcript.exons_near(low, high, stretch))

    @functools.cached_property
    def _at_intron_boundary(self) -> bool:
        low, high = self.allele.low, self.allele.high
        return any(
            (start - 3 <= high and low <= start + 7) or (end - 7 <= high and low <= end + 3)
            for start, end in self.transcript.introns_near(low, high, 3)
        )

    @functools.cached_property
    def _overlaps_coding_region(self) -> bool:
        return _overlap(self.allele.low, self.allele.high, *self._coding_region)

    @functools.cached_property
    def _outside_coding_region(self) -> bool:
        """VEP's `utr`: the allele has a base, or a neighbour, outside the coding region."""
        first, last = self._coding_region
        return not self._overlaps_coding_region or self.allele.low < first or self.allele.high > last

    @functools.cached_property
    def segments(self) -> tuple[transcript_models.Segment, ...]:
        return self.transcript.segments(self.allele.start, self.allele.end)

    @functools.cached_property
    def _coding_segments(self) -> tuple[transcript_models.Segment, ...]:
        """The allele's pieces by coding position, counted from 1 at the first coding base, without the phase."""
        first, last = self._coding.cdna_start, self._coding.cdna_end
        pieces = []
        for piece in self.segments:
            if not piece.mapped or piece.end < first or piece.start > last:
                pieces.append(dataclasses.replace(piece, mapped=False))
                continue
            if piece.start < first:
                pieces.append(transcript_models.Segment(piece.start, first - 1, mapped=False))
            pieces.append(
                transcript_models.Segment(
                    max(piece.start, first) - first + 1, min(piece.end, last) - first + 1, mapped=True
                )
            )
            if piece.end > last:
                pieces.append(transcript_models.Segment(last + 1, piece.end, mapped=False))
        return tuple(pieces)

    @functools.cached_property
    def cds_start(self) -> int | None:
        """The coding position of the allele's first base on the transcript, or None when it has none."""
        first = self._coding_segments[0]
        return first.start + self._coding.phase if first.mapped else None

    @functools.cached_property
    def cds_end(self) -> int | None:
        last = self._coding_segments[-1]
        return last.end + self._coding.phase if last.mapped else None

    @functools.cached_property
    def _cdna_ends(self) -> tuple[int | None, int | None]:
        first, last = self.segments[0], self.segments[-1]
        return (first.start if first.mapped else None, last.end if last.mapped else None)

    @functools.cached_property
    def _codons(self) -> tuple[int, int] | None:
        """The first and last codon that the allele lies in, counted from 1, or None without both coding ends."""
        if self.cds_start is None or self.cds_end is None:
            return None
        return ((self.cds_start + 2) // 3, (self.cds_end + 2) // 3)

    @functools.cached_property
    def _intron_effects(self) -> tuple[bool, bool, bool]:
        """Whether the allele lies in a frameshift intron, and at the lower and the upper splice site of an intron."""
        low, high = self.allele.low, self.allele.high
        in_frameshift_intron = lower_site = upper_site = False
        for region_start, region_end in self.allele.differing:
            for intron in self.transcript.introns_near(low, high, 3):
                start, end = intron
                if transcript_models.is_frameshift_intron(intron) and _overlap(region_start, region_end, start, end):
                    in_frameshift_intron = True
                    continue
                if (start - 3 <= high and low <= start + 7) or (end - 7 <= high and low <= end + 3):
                    lower_site = lower_site or _overlap(region_start, region_end, start, start + 1)
                    upper_site = upper_site or _overlap(region_start, region_end, end - 1, end)
        return in_frameshift_intron, lower_site, upper_site

    @functools.cached_property
    def _within_cdna(self) -> bool:
        length = self.transcript.length
        if any(piece.mapped and piece.end > 0 and piece.start <= length for piece in self.segments):
            return True
        return self._intron_effects[0] and self.within_transcript

    # The consequence terms

    @functools.cached_property
    def five_prime_utr_variant(self) -> bool:
        if not self.within_transcript or self._ablates or not self._in_exon or not self._outside_coding_region:
            return False
        allele, transcript = self.allele, self.transcript
        first, last = self._coding_region
        if transcript.strand == 1:
            before = (self._insertion and allele.start == first) or _overlap(
                allele.start, allele.end, transcript.start, first - 1
            )
        else:
            before = (self._insertion and allele.end == last) or _overlap(
                allele.start, allele.end, last + 1, transcript.end
            )
        return before and self._within_cdna

    @functools.cached_property
    def splice_site_variant(self) -> bool:
        """Whether VEP gives `splice_donor_variant` or `splice_acceptor_variant`."""
        if not self.within_transcript or self._ablates or not self._at_intron_boundary:
            return False
        _in_frameshift_intron, lower_site, upper_site = self._intron_effects
        return lower_site or upper_site

    @functools.cached_property
    def _is_coding(self) -> bool:
        """VEP's `coding`: the allele overlaps the coding region and an exon."""
        return self.within_transcript and not self._ablates and self._overlaps_coding_region and self._in_exon

    def _codon(self, *, reference: bool) -> str | None:
        """The codons that the allele lies in, with the reference allele or with the alternate one; `-` for none."""
        codons = self._codons
        if codons is None or not codons[0] or not codons[1] or self.cds_start is None or self.cds_end is None:
            return None
        first_base, last_base = codons[0] * 3 - 2, codons[1] * 3
        length = last_base - first_base + 1
        if reference:
            codon = self.sequences.coding[first_base - 1 : first_base - 1 + length]
        else:
            # The coding sequence with the allele, and the 3' UTR behind it, from the first of these codons on: the
            # bases before the allele, the allele, and as many bases behind it as the codons still need.
            length = max(length + len(self._alt) - (self.cds_end - self.cds_start + 1), 0)
            coding = self.sequences.coding
            head = coding[first_base - 1 : self.cds_start - 1] + self._alt
            behind = coding[self.cds_end : self.cds_end + max(length - len(head), 0)]
            behind += self.sequences.utr3[: max(length - len(head) - len(behind), 0)]
            codon = (head + behind)[:length]
        return codon or '-'

    def _peptide(self, *, reference: bool) -> str | None:
        if not reference and not _UNAMBIGUOUS.fullmatch(self.allele.alt):
            return None
        if reference and not _UNAMBIGUOUS.fullmatch(self.allele.ref):
            return None
        codon = self._codon(reference=reference)
        if codon is None:
            return None
        if codon == '-':
            return '-'
        whole = len(codon) - len(codon) % 3
        peptide = translate(codon[:whole], self._coding.codon_table)
        if len(codon) > whole and peptide != '*':
            peptide += 'X'
        return peptide or '-'

    @functools.cached_property
    def _peptides(self) -> tuple[str, str] | None:
        """The amino acids of the allele's codons, reference and alternate, or None when VEP has not both."""
        alternate = self._peptide(reference=False)
        reference = self._peptide(reference=True) if alternate else None
        if not alternate or not reference:
            return None
        return ('' if reference == '-' else reference, '' if alternate == '-' else alternate)

    @functools.cached_property
    def _partial_codon(self) -> bool:
        """Whether the allele begins in a last codon that the coding sequence does not complete."""
        first = self._coding_segments[0]
        if not first.mapped:
            return False
        codon = (first.start + self._coding.phase + 2) // 3
        remaining = len(self.sequences.coding) - (codon * 3 - 3)
        return 0 < remaining < 3

    def _stop_altered(self, *, overlaps_stop: bool) -> bool:
        """Whether an insertion or deletion at the stop codon leaves no stop codon where the coding sequence ends."""
        allele = self.allele
        if not _UNAMBIGUOUS.fullmatch(allele.alt) or not overlaps_stop or len(allele.ref) == len(allele.alt):
            return False
        cdna_start, cdna_end = self._cdna_ends
        if not cdna_start or not cdna_end or not self.cds_start:
            return False
        coding = self.sequences.coding
        changed = coding + self.sequences.utr3
        changed = changed[: self.cds_start - 1] + self._alt + changed[self.cds_start - 1 + cdna_end - cdna_start + 1 :]
        if len(changed) < len(coding):
            return True
        return translate(changed[len(coding) - 3 : len(coding)], self._coding.codon_table) != '*'

    @functools.cached_property
    def _overlaps_stop_by_cdna(self) -> bool:
        cdna_start, cdna_end = self._cdna_ends
        if not cdna_start or not cdna_end:
            return False
        return _overlap(cdna_start, cdna_end, self._coding.cdna_end - 2, self._coding.cdna_end)

    @functools.cached_property
    def _overlaps_stop_by_genome(self) -> bool:
        """VEP's test for an insertion or deletion, which moves an insertion's start by the inserted length."""
        start, end = self.allele.start, self.allele.end
        if self._insertion and re.fullmatch(r'[ACGTN]+', self._alt):
            start += len(self._alt) * self.transcript.strand
        first, last = self._coding_region
        stop = (last - 2, last) if self.transcript.strand == 1 else (first, first + 2)
        return _overlap(start, end, *stop)

    @functools.cached_property
    def _stop_lost(self) -> bool:
        if self._partial_codon:
            return False
        peptides = self._peptides
        if peptides is not None and 'X' not in peptides[1]:
            return '*' in peptides[0] and '*' not in peptides[1]
        return self._stop_altered(overlaps_stop=self._overlaps_stop_by_cdna)

    @functools.cached_property
    def _stop_retained(self) -> bool:
        if self._partial_codon or self._stop_lost:
            return False
        peptides = self._peptides
        if peptides is not None and peptides[1] and 'X' not in peptides[1]:
            return self._same_protein(*peptides)
        allele = self.allele
        return (
            len(allele.ref) != len(allele.alt)
            and self._overlaps_stop_by_genome
            and not self._stop_altered(overlaps_stop=self._overlaps_stop_by_genome)
        )

    def _same_protein(self, reference: str, alternate: str) -> bool:
        """VEP's `ref_eq_alt_sequence`: whether the protein with the alternate amino acids still ends as it did."""
        if reference == 'X' and alternate == 'X':
            return False
        codons = self._codons
        if codons is None:
            return False
        whole = protein(self.sequences, self._coding.codon_table)
        if codons[0] > len(whole) and alternate.startswith('*'):
            return True
        if codons[0] - 1 > len(whole):
            return False
        changed = whole[: codons[0] - 1] + alternate + whole[max(codons[1], codons[0] - 1) :]
        beyond = changed[len(whole) :]
        if whole == changed[: len(whole)] and beyond.startswith('*'):
            return True
        return '*' in reference and reference.index('*') == alternate.find('*')

    @functools.cached_property
    def stop_gained(self) -> bool:
        if not self._is_coding or self._stop_retained or self._stop_lost:
            return False
        peptides = self._peptides
        return peptides is not None and '*' in peptides[1] and '*' not in peptides[0]

    @functools.cached_property
    def frameshift_variant(self) -> bool:
        allele = self.allele
        if not self._is_coding or len(allele.ref) == len(allele.alt):
            return False
        if self._partial_codon or self._stop_retained or self.cds_start is None or self.cds_end is None:
            return False
        if self._peptides is not None and self._peptides[0].startswith('*'):
            return False
        return abs(len(allele.alt) - (self.cds_end - self.cds_start + 1)) % 3 != 0
