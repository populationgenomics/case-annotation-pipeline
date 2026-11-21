"""What an allele in a 5' UTR does to the upstream reading frames of a transcript (docs/annotation.md).

These are the rules of the UTRannotator plugin for Ensembl VEP, by Xiaolei Zhang, Matthew Wakeling, James Ware and
Nicola Whiffin (https://github.com/ImperialCardioGenetics/UTRannotator, MIT, licenses/UTRannotator-LICENSE.txt),
restated in Python from the plugin's seven Perl files at the commit SOURCE_COMMIT: `UTRannotator.pm` and, under
`src/`, `five_prime_UTR_utils.pl`, `uAUG_gained.pl`, `uAUG_lost.pl`, `uSTOP_gained.pl`, `uSTOP_lost.pl` and
`uframeshift.pl`. The Rust port of Ensembl's edition of the plugin in vep-rs
(https://github.com/natera-open-source/vep-rs, Apache 2.0, tag v0.6.0, `crates/vep-builtins/src/plugins/
utr_annotator.rs`) and its pages on that edition's defects were read for the logic. No code of either is copied.

An upstream start codon is an `ATG` in the 5' UTR. Its reading frame is an upstream ORF (`uORF`) when a stop codon
ends it inside the UTR, and otherwise an overlapping ORF that runs into the coding sequence, in its frame
(`InFrame_oORF`) or out of it (`OutOfFrame_oORF`). The plugin has five analyses, and `effects` gives the records of
every one that applies:

- `uAUG_gained`: the allele makes a start codon;
- `uAUG_lost`: it removes one;
- `uSTOP_gained`: a substitution puts a stop codon into a reading frame, before the frame's own end;
- `uSTOP_lost`: it removes the stop codon of an upstream ORF;
- `uFrameShift`: an insertion or deletion inside a reading frame changes the frame's length by no multiple of three.

The plugin decides an insertion or deletion by the place where the record writes it, so one changed sequence can get
different records. Here the reference UTR and the changed UTR are compared, so that one changed sequence gets one set
of records. A start codon is lost only when no way of writing the allele leaves its three bases alone, and gained
only when no way of writing the allele makes it one of the reference's. A stop codon is lost when a way of writing
the allele changes it and the frame does not end at its place afterwards. A start codon that the allele repeats, or
removes with its repeat, is counted at the first copy. The records also differ from the plugin's in a few places
where the plugin's own rule fails. For a substitution of one base the only such place is an allele at the first two
bases of the UTR. docs/annotation.md lists every difference with an example, and tests/test_utr5.py pins each.

One thing here is not the plugin's: the Kozak context and strength of the transcript's own start codon, with the
reference allele and with an allele that changes the context (kozak, own_start_change), by the rule that the plugin
applies to an upstream start codon.
"""

from __future__ import annotations

import bisect
import dataclasses
import re
from collections.abc import Callable, Iterable
from typing import Any

# The commit of https://github.com/ImperialCardioGenetics/UTRannotator that the rules were read at.
SOURCE_COMMIT = '3d9d466e356752aa69d61b62da493edddf374363'

_START = 'ATG'
_STOPS = frozenset({'TAA', 'TAG', 'TGA'})
_STOP_CODON = re.compile(r'(?=TAA|TAG|TGA)')


@dataclasses.dataclass(frozen=True)
class Change:
    """An allele inside a 5' UTR, on the transcript's strand.

    `position` counts the UTR's bases from 0 at its first. The allele replaces `removed` bases from there with
    `inserted`. An insertion removes none and stands before the base at `position`.
    """

    position: int
    removed: int
    inserted: str


def kozak(sequence: str, start: int) -> tuple[str | None, str | None]:
    """The Kozak context of the start codon at `start` in a sequence, and its strength.

    The context is the three bases before the codon, the codon and the base after it. It is strong with a purine
    three bases before the codon and a `G` after it, moderate with one of the two and weak with neither. A codon
    without three bases before it or one after it in the sequence has neither.
    """
    if start < 3 or start + 3 >= len(sequence):
        return None, None
    context = sequence[start - 3 : start + 4]
    purine, guanine = context[0] in 'AG', context[6] == 'G'
    return context, 'Strong' if purine and guanine else 'Moderate' if purine or guanine else 'Weak'


def _starts(sequence: str) -> list[int]:
    found, at = [], sequence.find(_START)
    while at >= 0:
        found.append(at)
        at = sequence.find(_START, at + 1)
    return found


class ReadingFrames:
    """A 5' UTR with the coding sequence behind it: its start codons, and where each one's reading frame ends."""

    def __init__(
        self, utr: str, coding: str, coding_stops: tuple[int | None, int | None, int | None] | None = None
    ) -> None:
        self.utr = utr
        self.coding = coding
        self.starts = _starts(utr)
        self._stops: tuple[list[int], list[int], list[int]] = ([], [], [])
        for match in _STOP_CODON.finditer(utr):
            self._stops[match.start() % 3].append(match.start())
        self._coding_stops = coding_stops

    @property
    def coding_stops(self) -> tuple[int | None, int | None, int | None]:
        """The first stop codon of the coding sequence at a place with each remainder on division by three."""
        if self._coding_stops is None:
            first: list[int | None] = [None, None, None]
            for match in _STOP_CODON.finditer(self.coding):
                if first[match.start() % 3] is None:
                    first[match.start() % 3] = match.start()
                    if None not in first:
                        break
            self._coding_stops = (first[0], first[1], first[2])
        return self._coding_stops

    def stop_in_utr(self, start: int) -> int | None:
        """The first stop codon in the frame of `start` that lies behind it and wholly in the UTR."""
        stops = self._stops[start % 3]
        at = bisect.bisect_right(stops, start)
        return stops[at] if at < len(stops) else None

    def stop(self, start: int) -> int | None:
        """The first stop codon in the frame of `start` behind it, in the UTR or in the coding sequence."""
        in_utr = self.stop_in_utr(start)
        if in_utr is not None:
            return in_utr
        length = len(self.utr)
        # A codon that begins in the UTR's last two bases ends in the coding sequence.
        for at in (length - 2, length - 1):
            if at > start and (at - start) % 3 == 0 and self.utr[at:] + self.coding[: 3 - (length - at)] in _STOPS:
                return at
        in_coding = self.coding_stops[(start - length) % 3]
        return None if in_coding is None else length + in_coding

    def type(self, start: int) -> str:
        """The type of the reading frame of the start codon at `start`."""
        if self.stop_in_utr(start) is not None:
            return 'uORF'
        return 'OutOfFrame_oORF' if (len(self.utr) - start) % 3 else 'InFrame_oORF'

    def existing(self) -> dict[str, int]:
        """The reading frames that the UTR has, as the plugin counts them: one for each start codon in it."""
        types = [self.type(start) for start in self.starts]
        return {
            'utr5_existing_uorfs': types.count('uORF'),
            'utr5_existing_oorfs_out_of_frame': types.count('OutOfFrame_oORF'),
            'utr5_existing_oorfs_in_frame': types.count('InFrame_oORF'),
        }


def _shared(reference: str, changed: str, first: int, removed: int) -> tuple[int, int]:
    """How many bases a sequence shares with its changed form at their beginning, and how many at their end.

    The change replaces `removed` bases from `first`. For an insertion or deletion among repeated bases the two
    stretches can overlap. A codon that lies wholly in either is one that the allele can be written without
    touching.
    """
    shortest = min(len(reference), len(changed))
    prefix = first
    while prefix < shortest and reference[prefix] == changed[prefix]:
        prefix += 1
    suffix = len(reference) - first - removed
    while suffix < shortest and reference[-1 - suffix] == changed[-1 - suffix]:
        suffix += 1
    return prefix, suffix


@dataclasses.dataclass(frozen=True)
class _Compared:
    """A reference UTR beside its changed form: what the two share, and what that says of a codon.

    Among repeated bases an insertion or deletion can be written at several places, and the shared beginning and
    the shared end then overlap.
    """

    utr: str
    changed: str
    prefix: int  # the bases that the two share at their beginning: the last place where the change can begin
    suffix: int  # the bases that they share at their end

    @property
    def same_length(self) -> bool:
        return len(self.utr) == len(self.changed)

    @property
    def pure(self) -> bool:
        """Whether the changed UTR is the reference with bases put in at one place, or taken out at one place.

        The alleles of a record do not say: a record can write a deletion as a longer allele replaced by a shorter.
        """
        return self.prefix + self.suffix >= min(len(self.utr), len(self.changed)) and not self.same_length

    @property
    def reference_end(self) -> int:
        """Where the shared end begins in the reference.

        A codon of the reference from there on stands whole in the changed UTR with everything behind it, so its
        reading frame is as it was.
        """
        return len(self.utr) - self.suffix

    @property
    def changed_end(self) -> int:
        return len(self.changed) - self.suffix

    def stays(self, start: int) -> bool:
        """Whether a start codon of the reference stands at its own place, where no codon of the shared end does."""
        return start + 3 <= self.prefix and start < self.changed_end

    def in_place(self, start: int) -> bool:
        """Whether a replacement of bases by others leaves a start codon at a place where both sequences have one.

        The codon at that place must be of neither shared end: a codon of the shared end is the reference's codon
        from further on or further back, and is not counted a second time.
        """
        both = self.utr[start : start + 3] == _START and self.changed[start : start + 3] == _START
        return not self.pure and both and start < self.reference_end and start < self.changed_end

    def inherited(self, start: int) -> bool:
        """Whether a start codon of the changed UTR is one of the reference's.

        It is one of the shared end, or one of the shared beginning whose codon in the reference is not also found
        in the shared end, or one that a replacement leaves in place.
        """
        of_beginning = start + 3 <= self.prefix and start < self.reference_end
        return start >= self.changed_end or of_beginning or self.in_place(start)

    def stop_touched(self, codon: int) -> bool:
        """Whether a way of writing the allele changes a base of a stop codon of the reference or parts it."""
        if not self.pure:
            return codon + 3 > self.prefix and codon < self.reference_end
        deleted = len(self.utr) - len(self.changed)
        if deleted < 0:
            return self.reference_end <= codon + 2 and self.prefix >= codon + 1
        # The first place where the deletion can begin is as many bases before the shared end as it removes.
        return self.reference_end - deleted <= codon + 2 and self.prefix + deleted > codon


def effects(frames: ReadingFrames, change: Change, evidence: Callable[[int], bool]) -> list[dict[str, Any]]:
    """The records of an allele in a 5' UTR: one for each start codon and analysis that applies.

    `evidence` says of a place in the UTR whether a start codon there was seen translated (uorf_evidence). The
    records are in no particular order.
    """
    utr = frames.utr
    changed = utr[: change.position] + change.inserted + utr[change.position + change.removed :]
    if changed == utr:
        return []
    compared = _Compared(utr, changed, *_shared(utr, changed, change.position, change.removed))
    prefix, same_length, reference_end = compared.prefix, compared.same_length, compared.reference_end
    after = ReadingFrames(changed, frames.coding, frames.coding_stops)
    shift = (len(changed) - len(utr)) % 3 != 0

    def of_start(start: int) -> dict[str, Any]:
        context, strength = kozak(utr, start)
        return {'kozak_context': context, 'kozak_strength': strength, 'evidence': evidence(start)}

    def distance(stop: int | None, start: int) -> int | None:
        return None if stop is None else stop - start

    records: list[dict[str, Any]] = []
    for start in after.starts:
        if compared.inherited(start):
            continue
        context, strength = kozak(changed, start)
        records.append(
            {
                'consequence': 'uAUG_gained',
                'type': after.type(start),
                'kozak_context': context,
                'kozak_strength': strength,
                'cap_distance_to_start': start,
                'distance_to_cds': len(changed) - start,
                'distance_to_stop': distance(after.stop(start), start),
            }
        )
    for start in frames.starts:
        if start >= reference_end:
            break
        if not compared.stays(start) and not compared.in_place(start):
            records.append(
                {
                    'consequence': 'uAUG_lost',
                    'type': frames.type(start),
                    'cap_distance_to_start': start,
                    'distance_to_cds': len(utr) - start,
                    'distance_to_stop': distance(frames.stop(start), start),
                }
                | of_start(start)
            )
        stop = frames.stop_in_utr(start)
        if stop is not None and compared.stop_touched(stop):
            # The stop codon stays where a stop codon still stands at its place and the frame up to it is the same.
            if changed[stop : stop + 3] in _STOPS and (same_length or prefix >= stop):
                continue
            # It also stays behind an insertion or deletion of whole codons that can be written before it.
            if same_length or shift or stop < reference_end:
                other = after.stop_in_utr(start) if changed[start : start + 3] == _START else None
                records.append(
                    {
                        'consequence': 'uSTOP_lost',
                        'alt_stop': other is not None,
                        'alt_stop_distance_to_cds': None if other is None else len(changed) - other,
                        'frame_with_cds': 'outOfFrame' if (len(changed) - start) % 3 else 'inFrame',
                    }
                    | of_start(start)
                )
            continue
        if not compared.stays(start):
            continue
        frame_end = len(utr) - 1 if stop is None else stop - 1
        if same_length:
            new_stop = after.stop_in_utr(start)
            if reference_end - 1 <= frame_end and new_stop is not None and new_stop < frame_end:
                records.append(
                    {
                        'consequence': 'uSTOP_gained',
                        'ref_type': frames.type(start),
                        'ref_start_distance_to_cds': len(utr) - start,
                        'new_stop_distance_to_cds': len(changed) - new_stop,
                    }
                    | of_start(start)
                )
            continue
        # The frame is shifted when the allele can only be written between its start codon and its stop codon.
        if shift and (stop is None or stop >= reference_end):
            reference_stop, changed_stop = frames.stop(start), after.stop(start)
            records.append(
                {
                    'consequence': 'uFrameShift',
                    'ref_type': frames.type(start),
                    'ref_type_length': None if reference_stop is None else reference_stop - start + 3,
                    'ref_start_distance_to_cds': len(utr) - start,
                    'alt_type': after.type(start),
                    'alt_type_length': None if changed_stop is None else changed_stop - start + 3,
                }
                | of_start(start)
            )
    return records


def own_start_change(spliced: str, start: int, first: int, removed: int, inserted: str) -> tuple[str, str] | None:
    """The Kozak context and strength of a transcript's own start codon with an allele, where the allele changes it.

    `spliced` is the transcript's sequence and `start` the place of its start codon in it, counted from 0. The
    allele replaces `removed` bases from `first` with `inserted`. None for an allele that leaves the context as it
    is, for one that cannot be written without touching the start codon, and where the changed transcript has no
    three bases before the codon.
    """
    if first >= start + 4 or first + removed <= start - 3:
        return None
    changed = spliced[:first] + inserted + spliced[first + removed :]
    prefix, suffix = _shared(spliced, changed, first, removed)
    # The own start codon is the one with the coding sequence behind it, so the shared end is looked at first: among
    # repeated bases the codon at its old place can be one that `effects` counts as gained.
    if start >= len(spliced) - suffix:
        moved = start + len(changed) - len(spliced)
    elif prefix >= start + 3:
        moved = start
    else:
        return None
    context, strength = kozak(changed, moved)
    if context is None or strength is None or context == kozak(spliced, start)[0]:
        return None
    return context, strength


_EVIDENCE_HEADER = 'chr\tpos\t'


def uorf_evidence(lines: Iterable[str]) -> frozenset[tuple[str, int]]:
    """The start codons that the plugin's file lists as seen translated: chromosome and place of the codon's `A`.

    The file is tab-separated. Its first line names its columns, and every other line begins with a chromosome
    under its `chr` name and a position.

    Raises:
        ValueError: if the first line is not the header, a line has no chromosome and whole position, or the file
            lists nothing.
    """
    listed = set()
    for number, line in enumerate(lines, start=1):
        fields = line.rstrip('\n').split('\t')
        if number == 1:
            if not line.startswith(_EVIDENCE_HEADER):
                raise ValueError(f'the file of translated upstream ORFs begins {line[:40]!r}, not with its header')
            continue
        if len(fields) < 2 or not fields[0].startswith('chr') or not fields[1].isdigit():
            raise ValueError(f'line {number} of the file of translated upstream ORFs has no chromosome and position')
        listed.add((fields[0], int(fields[1])))
    if not listed:
        raise ValueError('the file of translated upstream ORFs lists none')
    return frozenset(listed)
