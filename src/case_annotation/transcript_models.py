"""The transcripts of a prepared gene model file, built as Ensembl VEP builds them (docs/annotation.md).

The rules of the two plugins (utr5, nmd) and the conditions under which VEP gives an allele to a plugin
(plugin_gates) are stated over a transcript: its exons, its coding start and end and its spliced sequence. The tables
held the plugins' values as VEP 116.3 computed them from a GFF file and a FASTA file, without its cache, so this
module reads a gene model file the way VEP's GFF reader does. What follows from that reading, and is easy to take
for an oversight:

- A transcript is a feature whose parent is a gene, and it has a coding sequence when CDS features name it as their
  parent. Its type does not matter, so a `V_gene_segment` with CDS lines is read like an mRNA.
- The coding sequence runs from the start of the first CDS line to the end of the last, in the transcript's order. A
  gap or an overlap between two CDS lines of one exon is not looked at.
- A CDS line belongs to the exon that holds its first base on the genome. A CDS line that begins in no exon is left
  out.
- The phase of the first CDS line pads the coding sequence with one or two `N`, and moves a coding position by as
  much.
- A transcript of the file that gives its exons a `rank` (Ensembl's) keeps the start and end of its own line. Any
  other transcript (NCBI's) takes them from its exons, and is left out when two of its exons overlap.

The reading follows `Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF` and `GFF` of ensembl-vep and
`Bio::EnsEMBL::Transcript` and `TranscriptMapper` of ensembl, as release 116.3 of the VEP image holds them (Apache
2.0, licenses/ENSEMBL-NOTICE.txt). No code of theirs is copied.
"""

from __future__ import annotations

import bisect
import dataclasses
import functools
import itertools
import pathlib
import re
from collections.abc import Iterator, Mapping, Sequence

import pysam

from case_annotation import gene_models

# The feature types that VEP's GFF reader keeps. A line of another type is dropped before the features are linked, so
# a feature whose parent is of another type has no parent.
_TYPES = frozenset(
    {
        'aberrant_processed_transcript',
        'CDS',
        'C_gene_segment',
        'D_gene_segment',
        'exon',
        'gene',
        'J_gene_segment',
        'lnc_RNA',
        'lincRNA',
        'lincRNA_gene',
        'miRNA',
        'miRNA_gene',
        'mRNA',
        'mt_gene',
        'ncRNA',
        'ncRNA_gene',
        'nc_primary_transcript',
        'NMD_transcript_variant',
        'primary_transcript',
        'processed_pseudogene',
        'processed_transcript',
        'protein_coding_gene',
        'pseudogene',
        'pseudogenic_transcript',
        'RNA',
        'rRNA',
        'rRNA_gene',
        'snoRNA',
        'snoRNA_gene',
        'snRNA',
        'snRNA_gene',
        'supercontig',
        'transcript',
        'tRNA',
        'VD_gene_segment',
        'V_gene_segment',
    }
)
_VERSIONED = re.compile(r'\.\d+$')
_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')
# An intron of at most this many bases beyond its first is one that a gene build put in to mend a reading frame. VEP
# does not call an allele in it intronic or at a splice site.
_FRAMESHIFT_INTRON = 12
# The mitochondrial genome under the pipeline's name, whose transcripts are translated by the vertebrate mitochondrial
# code. VEP means to do the same and looks for the name `MT`, which it does not find when the reference names the
# genome `chrM`, so it translates them by the standard code (plugin_gates).
_MITOCHONDRIAL = 'chrM'


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclasses.dataclass(frozen=True)
class Coding:
    """Where a transcript's coding sequence lies."""

    cdna_start: int  # the first coding base, counted from 1 along the spliced transcript
    cdna_end: int  # the last one
    # The bases of the first codon that the transcript lacks: 0, 1 or 2. The coding sequence is padded with as many
    # `N`, and a coding position counts them.
    phase: int
    codon_table: int  # 1 for the standard genetic code, 2 for the vertebrate mitochondrial one


@dataclasses.dataclass(frozen=True)
class Segment:
    """A piece of an allele's span on a transcript, in the transcript's order.

    A piece inside an exon has `mapped` set and gives its positions along the spliced transcript, counted from 1. A
    piece outside the exons gives its positions on the genome. An insertion lies between two bases and has a start
    one above its end.
    """

    start: int
    end: int
    mapped: bool


@dataclasses.dataclass(frozen=True)
class Transcript:
    """One transcript of a gene model file."""

    transcript_id: str  # with its version, as `transcripts` holds it
    chrom: str  # as the pipeline's VCFs and reference name it
    strand: int  # 1 or -1
    start: int
    end: int
    exons: tuple[tuple[int, int], ...]  # start and end on the genome, in the transcript's order
    coding: Coding | None

    @functools.cached_property
    def exon_cdna_starts(self) -> tuple[int, ...]:
        """The position along the spliced transcript of each exon's first base, counted from 1."""
        starts, position = [], 1
        for start, end in self.exons:
            starts.append(position)
            position += end - start + 1
        return tuple(starts)

    @functools.cached_property
    def length(self) -> int:
        return sum(end - start + 1 for start, end in self.exons)

    @functools.cached_property
    def introns(self) -> tuple[tuple[int, int], ...]:
        """Start and end on the genome of each intron that has a base, in the transcript's order."""
        found = []
        for (start, end), (next_start, next_end) in itertools.pairwise(self.exons):
            intron = (end + 1, next_start - 1) if self.strand == 1 else (next_end + 1, start - 1)
            if intron[0] <= intron[1]:
                found.append(intron)
        return tuple(found)

    @functools.cached_property
    def has_frameshift_intron(self) -> bool:
        return any(is_frameshift_intron(intron) for intron in self.introns)

    @functools.cached_property
    def coding_region(self) -> tuple[int, int] | None:
        """The first and last base of the coding sequence on the genome, the lower one first."""
        if self.coding is None:
            return None
        ends = (self.genomic(self.coding.cdna_start), self.genomic(self.coding.cdna_end))
        return (min(ends), max(ends))

    @functools.cached_property
    def _ascending(self) -> tuple[tuple[int, int, int], ...]:
        """The exons by their place on the genome: start, end and index in the transcript's order."""
        return tuple(sorted((start, end, index) for index, (start, end) in enumerate(self.exons)))

    def genomic(self, cdna: int) -> int:
        """The base of the genome at a position along the spliced transcript."""
        index = bisect.bisect_right(self.exon_cdna_starts, cdna) - 1
        if index < 0 or cdna > self.length:
            raise ValueError(f'{self.transcript_id} has no base {cdna}')
        start, end = self.exons[index]
        offset = cdna - self.exon_cdna_starts[index]
        return start + offset if self.strand == 1 else end - offset

    def cdna(self, position: int) -> int | None:
        """The position along the spliced transcript of a base of the genome, or None for a base outside the exons."""
        ascending = self._ascending
        at = bisect.bisect_right(ascending, (position, float('inf'), 0)) - 1
        if at < 0 or position > ascending[at][1]:
            return None
        start, end, index = ascending[at]
        offset = position - start if self.strand == 1 else end - position
        return self.exon_cdna_starts[index] + offset

    def segments(self, start: int, end: int) -> tuple[Segment, ...]:
        """The pieces of a span of the genome on this transcript, as VEP's mapper gives them.

        An insertion, whose start is one above its end, is mapped by its two neighbouring bases. With both in one
        exon it lies between their two positions. With one of them in an exon it lies directly before or after that
        position, and the other neighbour is dropped. With none it is one piece outside the exons.
        """
        if start == end + 1:
            return self._insertion(end, start)
        pieces: list[Segment] = []
        position = start
        for exon_start, exon_end, _index in self._ascending:
            if exon_end < position:
                continue
            if exon_start > end:
                break
            if exon_start > position:
                pieces.append(Segment(position, exon_start - 1, mapped=False))
                position = exon_start
            last = min(end, exon_end)
            first_cdna, last_cdna = self.cdna(position), self.cdna(last)
            if first_cdna is None or last_cdna is None:
                raise AssertionError(f'{self.transcript_id}: {position} to {last} is not in one exon')
            pieces.append(Segment(min(first_cdna, last_cdna), max(first_cdna, last_cdna), mapped=True))
            position = last + 1
        if position <= end:
            pieces.append(Segment(position, end, mapped=False))
        return tuple(pieces if self.strand == 1 else reversed(pieces))

    def _insertion(self, before: int, after: int) -> tuple[Segment, ...]:
        lower, upper = self.cdna(before), self.cdna(after)
        if lower is None and upper is None:
            return (Segment(after, before, mapped=False),)
        if lower is not None and upper is not None and abs(lower - upper) == 1:
            return (Segment(max(lower, upper), min(lower, upper), mapped=True),)
        pieces = []
        # On the forward strand the insertion follows the lower neighbour and precedes the upper one. On the reverse
        # strand it is the other way round.
        if lower is not None:
            pieces.append(
                Segment(lower + 1, lower, mapped=True) if self.strand == 1 else Segment(lower, lower - 1, True)
            )
        if upper is not None:
            pieces.append(
                Segment(upper, upper - 1, mapped=True) if self.strand == 1 else Segment(upper + 1, upper, True)
            )
        return tuple(pieces if self.strand == 1 else reversed(pieces))


def is_frameshift_intron(intron: tuple[int, int]) -> bool:
    return intron[1] - intron[0] <= _FRAMESHIFT_INTRON


@dataclasses.dataclass(frozen=True)
class Sequences:
    """A transcript's sequence, read from the reference in capitals and on the transcript's strand."""

    spliced: str
    utr5: str  # the bases before the coding sequence; empty for a transcript without one
    coding: str  # from the first coding base to the last, behind the `N` of the coding phase
    utr3: str  # the bases after the coding sequence


def sequences(transcript: Transcript, reference: pysam.FastaFile) -> Sequences:
    """The sequences of a transcript.

    Raises:
        ValueError: if the reference does not hold an exon of the transcript in full.
    """
    parts = []
    for start, end in transcript.exons:
        bases = reference.fetch(transcript.chrom, start - 1, end).upper()
        if len(bases) != end - start + 1:
            raise ValueError(
                f'the reference holds {len(bases)} bases of {transcript.chrom}:{start}-{end}, an exon of '
                f'{transcript.transcript_id}'
            )
        parts.append(bases if transcript.strand == 1 else reverse_complement(bases))
    spliced = ''.join(parts)
    coding = transcript.coding
    if coding is None:
        return Sequences(spliced=spliced, utr5='', coding='', utr3='')
    return Sequences(
        spliced=spliced,
        utr5=spliced[: coding.cdna_start - 1],
        coding='N' * coding.phase + spliced[coding.cdna_start - 1 : coding.cdna_end],
        utr3=spliced[coding.cdna_end :],
    )


@dataclasses.dataclass
class _Feature:
    type: str
    start: int
    end: int
    strand: int
    phase: str
    attributes: dict[str, str]
    line: str
    children: list[_Feature] = dataclasses.field(default_factory=list)

    @property
    def id(self) -> str | None:
        attributes = self.attributes
        return attributes.get('ID') or attributes.get('Name') or attributes.get('id') or attributes.get('name')

    @property
    def parents(self) -> list[str]:
        named = self.attributes.get('Parent') or self.attributes.get('parent')
        return named.split(',') if named else []


def _features(lines: Iterator[str], where: str) -> list[_Feature]:
    """The features of one sequence that VEP keeps, in the file's order and each line once."""
    features, seen = [], set()
    for line in lines:
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9:
            raise ValueError(f'{where}: a line of {len(fields)} columns, and a GFF feature line has nine')
        if fields[2] not in _TYPES or line in seen:
            continue
        seen.add(line)
        if fields[6] not in ('+', '-'):
            raise ValueError(f'{where}: a {fields[2]} line with the strand {fields[6]!r}')
        attributes = {}
        for pair in fields[8].split(';'):
            key, _, value = pair.partition('=')
            attributes[key] = value.split('=', 1)[0]
        features.append(
            _Feature(
                type=fields[2],
                start=int(fields[3]),
                end=int(fields[4]),
                strand=1 if fields[6] == '+' else -1,
                phase=fields[7],
                attributes=attributes,
                line=line,
            )
        )
    return features


def _linked(features: Sequence[_Feature]) -> list[tuple[_Feature, _Feature | None]]:
    """The transcript features of one sequence, each with its gene, after VEP's linking of children to parents.

    A feature with a parent is a child of the last feature before it that has the parent's id, or else of the last
    feature of the sequence that has it. A feature without a parent is a gene when its type says so, and its children
    are its transcripts. Any other feature without a parent is a transcript of its own.
    """
    by_id: dict[str, _Feature] = {}
    orphans = []
    for feature in features:
        for parent_id in feature.parents:
            if parent_id in by_id:
                by_id[parent_id].children.append(feature)
            elif parent_id:
                orphans.append(feature)
        if feature.id:
            by_id[feature.id] = feature
    for feature in orphans:
        for parent_id in feature.parents:
            parent = by_id.get(parent_id)
            if parent is not None and all(child is not feature for child in parent.children):
                parent.children.append(feature)
    found: list[tuple[_Feature, _Feature | None]] = []
    for feature in by_id.values():
        if any(feature.parents):
            continue
        if 'gene' in feature.type or feature.type in ('processed_transcript', 'RNA'):
            found.extend((child, feature) for child in feature.children)
        else:
            found.append((feature, None))
    return found


def _biotype(feature: _Feature, gene: _Feature | None) -> str | None:
    attributes = feature.attributes
    stated = attributes.get('biotype') or attributes.get('transcript_type') or attributes.get('transcript_biotype')
    if stated:
        return stated
    kind = feature.type.lower()
    if kind == 'mrna':
        return 'protein_coding'
    if kind == 'ncrna':
        return attributes.get('ncrna_class') or None
    if segment := re.fullmatch(r'([a-z]+)_gene_segment', kind):
        return f'IG_{segment[1].upper()}_gene'
    if gene is not None and gene.attributes.get('description', '').startswith('microRNA'):
        return 'miRNA'
    return attributes.get('gbkey') or None


def _ensembl_phase(gff_phase: str) -> int:
    """The bases of the first codon that a CDS line lacks, from the bases that the line gives before its first codon."""
    return {'0': 0, '1': 2, '2': 1}.get(gff_phase, 0)


def _transcript(feature: _Feature, gene: _Feature | None, chrom: str) -> Transcript | None:
    """The transcript of a feature, or None for a feature that VEP builds no transcript from.

    Raises:
        ValueError: if the feature has a child that is no exon, CDS, stop codon or microRNA, for which VEP stops.
    """
    if not feature.children or not _biotype(feature, gene):
        return None
    name = re.sub(r'^(gene|transcript):', '', feature.attributes.get('transcript_id') or feature.id or '', flags=re.I)
    version = feature.attributes.get('version') or feature.attributes.get('transcript_version')
    if version and not _VERSIONED.search(name):
        name = f'{name}.{version}'
    exons: list[_Feature] = []
    cds_by_start: dict[int, _Feature] = {}
    for child in feature.children:
        kind = child.type.lower()
        if kind == 'exon':
            exons.append(child)
        elif kind in ('cds', 'stop_codon'):
            cds_by_start[child.start] = child
        elif kind != 'mirna':
            raise ValueError(f'{name} has a child of the type {child.type}, which VEP does not read')
    if not exons:
        return None
    exons.sort(key=lambda exon: exon.start, reverse=feature.strand == -1)
    ranked = all('rank' in exon.attributes for exon in exons)
    if ranked:
        ranks = [int(exon.attributes['rank']) for exon in exons]
        if ranks != list(range(1, len(exons) + 1)):
            raise ValueError(f'{name}: the ranks of its exons are {ranks}, not their order on its strand')
    elif any('rank' in exon.attributes for exon in exons):
        raise ValueError(f'{name}: some of its exons have a rank and some have none')
    else:
        for exon, following in itertools.pairwise(exons):
            apart = following.start > exon.end if feature.strand == 1 else following.end < exon.start
            if not apart:
                return None
    ordered_cds: list[tuple[_Feature, int]] = []
    for index, exon in enumerate(exons):
        starts = sorted((s for s in cds_by_start if exon.start <= s <= exon.end), reverse=feature.strand == -1)
        ordered_cds.extend((cds_by_start.pop(s), index) for s in starts)
    coding = None
    if ordered_cds:
        starts, position = [], 1
        for exon in exons:
            starts.append(position)
            position += exon.end - exon.start + 1
        (first, first_exon), (last, last_exon) = ordered_cds[0], ordered_cds[-1]
        if feature.strand == 1:
            cdna_start = starts[first_exon] + first.start - exons[first_exon].start
            cdna_end = starts[last_exon] + last.end - exons[last_exon].start
        else:
            cdna_start = starts[first_exon] + exons[first_exon].end - first.end
            cdna_end = starts[last_exon] + exons[last_exon].end - last.start
        coding = Coding(
            cdna_start=cdna_start,
            cdna_end=cdna_end,
            phase=_ensembl_phase(first.phase),
            codon_table=2 if chrom == _MITOCHONDRIAL else 1,
        )
    return Transcript(
        transcript_id=name,
        chrom=chrom,
        strand=feature.strand,
        start=feature.start if ranked else min(exon.start for exon in exons),
        end=feature.end if ranked else max(exon.end for exon in exons),
        exons=tuple((exon.start, exon.end) for exon in exons),
        coding=coding,
    )


def sequence_names(synonyms: pathlib.Path) -> dict[str, dict[str, str]]:
    """For each chromosome of the tables, the names that the chromosome synonyms give it, each with itself.

    Raises:
        ValueError: if a line is no pair of names, or a chromosome of the tables is in no pair.
    """
    names: dict[str, set[str]] = {chrom: {chrom} for chrom in gene_models.CHROMOSOMES}
    for number, line in enumerate(synonyms.read_text('utf-8').splitlines(), start=1):
        pair = line.split('\t')
        if len(pair) != 2:
            raise ValueError(f'{synonyms}: line {number} is no pair of names')
        for name, other in (pair, reversed(pair)):
            if name in names:
                names[name].add(other)
    if alone := [chrom for chrom, known in names.items() if len(known) == 1]:
        raise ValueError(f'{synonyms}: no other name of {alone}')
    return {chrom: {name: name for name in sorted(known)} for chrom, known in names.items()}


def read(path: pathlib.Path, chrom: str, names: Mapping[str, str]) -> list[Transcript]:
    """The transcripts with a coding sequence that a prepared gene model file gives one chromosome.

    `names` holds the chromosome's names (sequence_names). The file names it by one of them, or holds nothing of it.

    Raises:
        ValueError: if the file names the chromosome in two ways, or a line of it is not what VEP reads (_features,
            _transcript).
    """
    with pysam.TabixFile(str(path)) as gff:
        held = [name for name in names if name in gff.contigs]
        if len(held) > 1:
            raise ValueError(f'{path} names {chrom} in several ways: {held}')
        if not held:
            return []
        features = _features(gff.fetch(held[0]), f'{path}: {held[0]}')
    transcripts = []
    for feature, gene in _linked(features):
        transcript = _transcript(feature, gene, chrom)
        if transcript is not None and transcript.coding is not None:
            transcripts.append(transcript)
    return transcripts


class Overlapping:
    """The transcripts of one chromosome by the bases that they span."""

    _BIN = 16_384

    def __init__(self, transcripts: Sequence[Transcript]) -> None:
        self._bins: dict[int, list[Transcript]] = {}
        for transcript in transcripts:
            for index in range(transcript.start // self._BIN, transcript.end // self._BIN + 1):
                self._bins.setdefault(index, []).append(transcript)

    def of(self, start: int, end: int) -> list[Transcript]:
        """The transcripts that hold a base from `start` to `end`, in the order they were given, each once."""
        found: dict[int, Transcript] = {}
        for index in range(start // self._BIN, end // self._BIN + 1):
            for transcript in self._bins.get(index, ()):
                if transcript.start <= end and start <= transcript.end:
                    found.setdefault(id(transcript), transcript)
        return list(found.values())
