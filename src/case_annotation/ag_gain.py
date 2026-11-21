"""The AG dinucleotides that an allele creates between a branch point and its acceptor (docs/annotation.md).

The spliceosome takes the first `AG` behind a branch point as the acceptor, so the stretch between a branch point and
the acceptor of its intron normally holds none. Zhang et al. describe a rule for the variants that create one there
and call them AG-gain variants (Zhang P, et al. Genome-wide detection of human intronic AG-gain variants located
between splicing branchpoints and canonical splice acceptor sites. Proc Natl Acad Sci U S A. 2023;120(46):e2314225120,
doi:10.1073/pnas.2314225120, CC BY 4.0). This module states that rule in Python from the article's text, its Box 1
and its Figure 2. The authors' tool AGAIN is under other terms (CC BY-NC-ND 4.0) and was not read.

The rule is stated over the end of one intron on its strand, with the intron's first branch point: the catalogued
branch point nearest the acceptor. The article's zone 1 runs from the base behind the branch point to the fourth base
before the exon, which leaves the acceptor's own three bases out. An allele gains an `AG` when the changed sequence
has one, both bases in the zone, that the reference does not have. Each gained `AG` has:

- its distance from the branch point, counted to its `G`, so the first base behind the branch point is at 1;
- its distance to the acceptor, a negative number: the bases of the intron behind its `G`, so that the last `AG` of
  the zone is at -3;
- whether a `C` or a `T` stands before it, which makes it the `YAG` of an acceptor;
- whether it lies in the article's high-risk region, at least 8 bases from the branch point;
- whether the reference already has an `AG` in the zone;
- the article's score: 2 for a gained `AG` of zone 1, and one more each for the high-risk region, for a `YAG` and
  for a zone without an `AG` in the reference, so 2 to 5.

The article's tool also reports an `AG` gained between the second branch point and the first (its zone 2, with the
score 1). The published results that the rule was held against have none, so that part is left out.

The article gives one example of each kind of insertion and deletion and no rule for one that can be written in
several places. Here the reference and the changed sequence are compared as utr5 compares them, so that one changed
sequence has one answer: an `AG` is gained when no way of writing the allele makes it one of the reference's, and one
that the allele repeats is counted at the first copy. The distances are counted in the changed sequence.
"""

from __future__ import annotations

import bisect
import dataclasses
import gzip
import pathlib
import re
from collections.abc import Iterable
from typing import Any

import pysam

from case_annotation import branch_points, gene_models, plugin_gates, transcript_models

_GAINED = 'AG'
_WHOLE = re.compile(r'[1-9][0-9]*')
_PYRIMIDINES = 'CT'
# The acceptor's own bases, the last three of the intron, which the zone leaves out.
ACCEPTOR_BASES = 3
# The first distance from the branch point, to the `G` of a gained `AG`, that lies in the high-risk region.
HIGH_RISK_FROM = 8
# A branch point counts for an intron from the 4th base before the exon behind it to the 99th: the distances of the
# catalogue's own rows by intron (branch_points).
NEAREST = 4
REACH = 99
# The score of a gained `AG` of zone 1 before the three things that each add one to it.
_BASE_SCORE = 2
# The bases of the reference that are read at a time beside an insertion or deletion (written_span).
_CHUNK = 32


@dataclasses.dataclass(frozen=True)
class Gain:
    """One `AG` that an allele creates between a branch point and its acceptor."""

    branch_point_distance: int
    acceptor_distance: int
    yag: bool
    high_risk: bool
    existing_ag: bool
    score: int


def _shared(reference: str, changed: str, first: int, removed: int) -> tuple[int, int]:
    """How many bases a sequence shares with its changed form at their beginning, and how many at their end.

    The change replaces `removed` bases from `first`. Among repeated bases an insertion or deletion can be written at
    several places, and the two stretches then overlap.
    """
    shortest = min(len(reference), len(changed))
    prefix = first
    while prefix < shortest and reference[prefix] == changed[prefix]:
        prefix += 1
    suffix = len(reference) - first - removed
    while suffix < shortest and reference[-1 - suffix] == changed[-1 - suffix]:
        suffix += 1
    return prefix, suffix


def _places(sequence: str, first: int, last: int) -> list[int]:
    """Where an `AG` begins in a sequence, with its `A` at `first` or later and its `G` at `last` or earlier."""
    found, at = [], sequence.find(_GAINED, first, last + 1)
    while at >= 0:
        found.append(at)
        at = sequence.find(_GAINED, at + 1, last + 1)
    return found


def gains(sequence: str, branch_point: int, behind: int, position: int, removed: int, inserted: str) -> list[Gain]:
    """The `AG`s that an allele gains in zone 1 of an intron, from the branch point on.

    `sequence` is a stretch of the reference on the intron's strand that holds the branch point at the index
    `branch_point`, every base from there to the intron's last, and `behind` bases of the exon after it. The allele
    replaces `removed` bases from the index `position` with `inserted`, and an insertion removes none and stands
    before the base at `position`. Bases before the branch point and bases of the exon are there for an allele that a
    record writes outside the zone and that can be written inside it.

    The stretch is read back from the exon, which keeps its place in the changed sequence. An allele is read when it
    can be written with the exon after it, and when no way of writing it changes the branch point or lies before it:
    an allele that can be written before the branch point leaves everything from the branch point to the exon as it
    was. Any other allele gains nothing here.
    """
    changed = sequence[:position] + inserted + sequence[position + removed :]
    if changed == sequence:
        return []
    prefix, suffix = _shared(sequence, changed, position, removed)
    # Where the shared end begins in each sequence. The branch point is before every way of writing the allele when
    # it lies in the shared beginning and not in the shared end.
    reference_end, changed_end = len(sequence) - suffix, len(changed) - suffix
    if prefix <= branch_point or reference_end <= branch_point or suffix < behind:
        return []
    # An `AG` of the changed sequence is the reference's when it lies in the shared end, or in the shared beginning
    # without also being a copy of one in the shared end, or where a replacement of bases leaves it in place.
    pure = prefix + suffix >= min(len(sequence), len(changed)) and len(sequence) != len(changed)

    def inherited(at: int) -> bool:
        in_place = not pure and sequence[at : at + 2] == _GAINED and at < reference_end
        return at >= changed_end or (at + 2 <= prefix and at < reference_end) or in_place

    zone_end = len(sequence) - behind - ACCEPTOR_BASES - 1
    existing_ag = bool(_places(sequence, branch_point + 1, zone_end))
    intron_end = len(changed) - behind
    found = []
    for at in _places(changed, branch_point + 1, intron_end - ACCEPTOR_BASES - 1):
        if inherited(at):
            continue
        distance = at + 1 - branch_point
        yag = changed[at - 1] in _PYRIMIDINES
        high_risk = distance >= HIGH_RISK_FROM
        found.append(
            Gain(
                branch_point_distance=distance,
                acceptor_distance=at + 2 - intron_end,
                yag=yag,
                high_risk=high_risk,
                existing_ag=existing_ag,
                score=_BASE_SCORE + high_risk + yag + (not existing_ag),
            )
        )
    return found


class BranchPoints:
    """The catalogued branch points of one chromosome: their positions on each strand."""

    def __init__(self, forward: Iterable[int], reverse: Iterable[int]) -> None:
        self._positions = {1: sorted(set(forward)), -1: sorted(set(reverse))}

    def first(self, strand: int, exon_start: int, intron_length: int) -> int | None:
        """The first branch point of an intron: the catalogued one nearest its acceptor, or None for an intron without.

        `exon_start` is the first base of the exon behind the intron. A branch point counts from NEAREST to REACH
        bases before that exon, and inside the intron.
        """
        reach = min(REACH, intron_length)
        positions = self._positions[strand]
        if strand == 1:
            at = bisect.bisect_right(positions, exon_start - NEAREST) - 1
            return positions[at] if at >= 0 and positions[at] >= exon_start - reach else None
        at = bisect.bisect_left(positions, exon_start + NEAREST)
        return positions[at] if at < len(positions) and positions[at] <= exon_start + reach else None


def read_catalogue(path: pathlib.Path) -> dict[str, BranchPoints]:
    """The branch points of the prepared catalogue (branch_points), gzipped or not, for each chromosome of the tables.

    The step that joins the catalogue to the alleles checks every column of it (joined_values). This reader takes a
    row's chromosome, position and strand, and checks those.

    Raises:
        ValueError: if the file cannot be read, its header is not the catalogue's, a row has another number of
            fields, a chromosome that is none of the tables', no whole position or no strand, or the file holds no
            row.
    """
    found: dict[str, dict[str, list[int]]] = {chrom: {'+': [], '-': []} for chrom in gene_models.CHROMOSOMES}
    columns = branch_points.COLUMNS
    at_chrom, at_pos, at_strand = (columns.index(name) for name in ('chrom', 'pos', 'strand'))
    rows = 0
    try:
        with path.open('rb') as f:
            gzipped = f.read(2) == b'\x1f\x8b'
        with gzip.open(path, 'rt', encoding='utf-8') if gzipped else path.open(encoding='utf-8') as lines:
            header = lines.readline()
            while header.startswith(branch_points.PREAMBLE):
                header = lines.readline()
            if header.rstrip('\n').split('\t') != list(columns):
                raise ValueError(f'{path}: its header is not that of the prepared catalogue of branch points')
            for line in lines:
                rows += 1
                fields = line.rstrip('\n').split('\t')
                if len(fields) != len(columns):
                    raise ValueError(f'{path}: a row of {len(fields)} fields, and the catalogue has {len(columns)}')
                chrom, pos, strand = fields[at_chrom], fields[at_pos], fields[at_strand]
                if chrom not in found or _WHOLE.fullmatch(pos) is None or strand not in found[chrom]:
                    raise ValueError(
                        f'{path}: a row without one of chr1 to chrM, a whole position and a strand: {fields[:4]}'
                    )
                found[chrom][strand].append(int(pos))
    except (OSError, EOFError, UnicodeDecodeError) as unreadable:
        raise ValueError(f'{path}: not the prepared catalogue of branch points: {unreadable}') from unreadable
    if not rows:
        raise ValueError(f'{path}: no branch point, and the prepared catalogue holds some')
    return {chrom: BranchPoints(strands['+'], strands['-']) for chrom, strands in found.items()}


def written_span(allele: plugin_gates.Allele, reference: pysam.FastaFile, chrom: str) -> tuple[int, int]:
    """The lowest and the highest base of the genome that an allele reaches over every way of writing it.

    An insertion or a deletion among repeated bases can be written at several places, and a record writes it at the
    leftmost on the genome, which is the far end of the repeat for a transcript of the reverse strand. Any other
    allele has the span of its record.
    """
    unit = allele.alt or allele.ref
    if allele.ref and allele.alt:
        return allele.low, allele.high
    length = reference.get_reference_length(chrom)

    def repeats(first: int, step: int) -> int:
        """How many bases from `first` on, in the direction `step`, go on with the unit's bases in a circle."""
        count = 0
        while 1 <= first + step * count <= length:
            low = max(first + step * count, 1) if step == 1 else max(first - count - _CHUNK + 1, 1)
            high = min(low + _CHUNK - 1, length) if step == 1 else first - count
            bases = reference.fetch(chrom, low - 1, high).upper()
            for base in bases if step == 1 else reversed(bases):
                if base != unit[count % len(unit) if step == 1 else -1 - count % len(unit)]:
                    return count
                count += 1
        return count

    return allele.low - repeats(allele.start - 1, -1), allele.high + repeats(allele.end + 1, 1)


def of_transcript(
    transcript: transcript_models.Transcript,
    allele: plugin_gates.Allele,
    span: tuple[int, int],
    reference: pysam.FastaFile,
    catalogued: BranchPoints,
) -> list[dict[str, Any]]:
    """The `AG`s that an allele gains on a transcript, each with the position of its branch point on the genome.

    `span` is the allele's reach over every way of writing it (written_span), by which the introns are found that it
    can be written in. Each of them is read with its first branch point (BranchPoints.first), and an intron without
    one gives nothing. The entries are in the transcript's order.
    """
    forward = transcript.strand == 1
    found: list[dict[str, Any]] = []
    for start, end in transcript.introns_near(*span, reach=1):
        exon_start = end + 1 if forward else start - 1
        branch_point = catalogued.first(transcript.strand, exon_start, end - start + 1)
        if branch_point is None or span[1] < min(branch_point, exon_start) or span[0] > max(branch_point, exon_start):
            continue
        # The stretch from the branch point to the intron's last base, on the genome, with the allele as written.
        low = min(branch_point if forward else start, allele.start)
        high = max(end if forward else branch_point, allele.end)
        bases = reference.fetch(transcript.chrom, low - 1, high).upper()
        if len(bases) != high - low + 1:
            raise ValueError(f'the reference holds {len(bases)} bases of {transcript.chrom}:{low}-{high}')
        if forward:
            entries = gains(bases, branch_point - low, high - end, allele.start - low, len(allele.ref), allele.alt)
        else:
            entries = gains(
                transcript_models.reverse_complement(bases),
                high - branch_point,
                start - low,
                high - allele.end,
                len(allele.ref),
                transcript_models.reverse_complement(allele.alt),
            )
        found.extend({'branch_point_pos': branch_point} | dataclasses.asdict(entry) for entry in entries)
    return sorted(
        found, key=lambda entry: (entry['branch_point_pos'] * transcript.strand, entry['branch_point_distance'])
    )
