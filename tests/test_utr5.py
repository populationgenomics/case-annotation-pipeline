"""The rules of UTRannotator as utr5 restates them, on invented sequences.

Every sequence here is invented. `UTR` holds two upstream ORFs, and `CODING` is a coding sequence that begins with
its start codon and has a stop codon in each reading frame.
"""

from __future__ import annotations

import random
from collections.abc import Iterator
from typing import Any

import pytest

from case_annotation import utr5

#      0         1         2         3         4
#      01234567890123456789012345678901234567890
UTR = 'GGCACCATGGCTCCGTAAGCGCCATGCCTTGAGCGGTCACC'
# An upstream ORF from the start codon at 6 to the stop codon at 15, and one from 23 to 29.
CODING = 'ATGGCTGAGCTGTAAGTAGC'


def records(utr: str, position: int, removed: int, inserted: str, coding: str = CODING) -> list[dict[str, Any]]:
    """The records of a change, each without its empty fields, sorted by kind and text."""
    found = utr5.effects(
        utr5.ReadingFrames(utr, coding), utr5.Change(position, removed, inserted), lambda _start: False
    )
    brief = [{key: value for key, value in record.items() if value is not None} for record in found]
    return sorted(brief, key=lambda record: (record['consequence'], str(record)))


def kinds(utr: str, position: int, removed: int, inserted: str) -> list[str]:
    return [record['consequence'] for record in records(utr, position, removed, inserted)]


LOST_FIRST_START = {
    'consequence': 'uAUG_lost',
    'type': 'uORF',
    'cap_distance_to_start': 6,
    'distance_to_cds': 35,
    'distance_to_stop': 9,
    'kozak_context': 'ACCATGG',
    'kozak_strength': 'Strong',
    'evidence': False,
}
FIRST_FRAME_SHIFTED = {
    'consequence': 'uFrameShift',
    'ref_type': 'uORF',
    'ref_type_length': 12,
    'ref_start_distance_to_cds': 35,
    'alt_type': 'uORF',
    'alt_type_length': 27,
    'kozak_context': 'ACCATGG',
    'kozak_strength': 'Strong',
    'evidence': False,
}


@pytest.mark.parametrize(
    ('sequence', 'start', 'expected'),
    [
        ('GCCATGG', 3, ('GCCATGG', 'Strong')),
        ('ACCATGG', 3, ('ACCATGG', 'Strong')),
        ('GCCATGC', 3, ('GCCATGC', 'Moderate')),
        ('TCCATGG', 3, ('TCCATGG', 'Moderate')),
        ('TCCATGC', 3, ('TCCATGC', 'Weak')),
        ('CCATGG', 2, (None, None)),
        ('GCCATG', 3, (None, None)),
    ],
)
def test_a_kozak_context_is_strong_with_a_purine_before_and_a_g_after(
    sequence: str, start: int, expected: tuple[str | None, str | None]
) -> None:
    assert utr5.kozak(sequence, start) == expected


def test_a_reading_frame_ends_at_the_first_stop_codon_in_its_frame() -> None:
    frames = utr5.ReadingFrames(UTR, CODING)
    assert frames.starts == [6, 23]
    assert (frames.stop_in_utr(6), frames.stop_in_utr(23)) == (15, 29)
    assert (frames.type(6), frames.type(23)) == ('uORF', 'uORF')
    assert frames.existing() == {
        'utr5_existing_uorfs': 2,
        'utr5_existing_oorfs_out_of_frame': 0,
        'utr5_existing_oorfs_in_frame': 0,
    }


def test_a_frame_without_a_stop_in_the_utr_ends_in_the_coding_sequence() -> None:
    in_frame = utr5.ReadingFrames('CCATGGCC', CODING)
    assert (in_frame.type(2), in_frame.stop_in_utr(2), in_frame.stop(2)) == ('InFrame_oORF', None, 8 + 12)
    out_of_frame = utr5.ReadingFrames('CCATGGC', CODING)
    assert (out_of_frame.type(2), out_of_frame.stop(2)) == ('OutOfFrame_oORF', 7 + 16)
    # A stop codon whose first two bases end the UTR and whose last begins the coding sequence.
    across = utr5.ReadingFrames('CCATGGCCTA', CODING)
    assert (across.type(2), across.stop(2)) == ('OutOfFrame_oORF', 8)
    assert in_frame.existing()['utr5_existing_oorfs_in_frame'] == 1
    assert out_of_frame.existing()['utr5_existing_oorfs_out_of_frame'] == 1


def test_a_frame_without_any_stop_has_no_end() -> None:
    assert utr5.ReadingFrames('CCATGGCC', 'ATGGCC').stop(2) is None


def test_a_substitution_that_makes_a_start_codon_is_a_gain() -> None:
    assert records(UTR, 28, 1, 'A') == [
        {
            'consequence': 'uAUG_gained',
            'type': 'OutOfFrame_oORF',
            'kozak_context': 'GCCATGA',
            'kozak_strength': 'Moderate',
            'cap_distance_to_start': 28,
            'distance_to_cds': 13,
            'distance_to_stop': 18,
        }
    ]


def test_a_substitution_in_a_start_codon_is_its_loss() -> None:
    assert records(UTR, 7, 1, 'C') == [LOST_FIRST_START]


def test_a_substitution_in_a_stop_codon_is_its_loss_and_names_the_next_stop() -> None:
    assert records(UTR, 16, 1, 'C') == [
        {
            'consequence': 'uSTOP_lost',
            'alt_stop': False,
            'frame_with_cds': 'outOfFrame',
            'kozak_context': 'ACCATGG',
            'kozak_strength': 'Strong',
            'evidence': False,
        }
    ]
    # With a later stop codon in the frame, the record says where it is.
    longer = 'GGCACCATGGCTCCGTAAGCGTGACC'
    (lost,) = records(longer, 16, 1, 'C')
    assert (lost['alt_stop'], lost['alt_stop_distance_to_cds'], lost['frame_with_cds']) == (True, 5, 'outOfFrame')


def test_a_substitution_that_leaves_another_stop_codon_at_the_place_is_no_loss() -> None:
    assert records(UTR, 17, 1, 'G') == []


def test_a_substitution_that_makes_a_stop_codon_in_a_frame_is_a_gained_stop() -> None:
    longer = 'GGCACCATGGCTCAGCCGTAAGCG'
    assert records(longer, 12, 1, 'T') == [
        {
            'consequence': 'uSTOP_gained',
            'ref_type': 'uORF',
            'ref_start_distance_to_cds': 18,
            'new_stop_distance_to_cds': 12,
            'kozak_context': 'ACCATGG',
            'kozak_strength': 'Strong',
            'evidence': False,
        }
    ]


def test_an_insertion_or_deletion_in_a_frame_shifts_it_unless_it_is_of_whole_codons() -> None:
    assert records(UTR, 10, 0, 'A') == [FIRST_FRAME_SHIFTED]
    assert records(UTR, 10, 2, '') == [FIRST_FRAME_SHIFTED | {'alt_type_length': 24}]
    assert records(UTR, 10, 0, 'AAA') == []
    assert records(UTR, 9, 3, '') == []


def test_a_change_far_behind_a_start_codon_meets_its_frame_and_no_frame_that_ended_before() -> None:
    """A frame that is open over thousands of bases is shifted, and one that a stop codon ended is left alone."""
    filler = 'GCC' * 2000
    utr = 'GGCATGGCTTAAGC' + 'ATGG' + filler + 'CC'
    (shifted,) = records(utr, len(utr) - 2, 0, 'C')
    assert (shifted['consequence'], shifted['ref_start_distance_to_cds']) == ('uFrameShift', len(utr) - 14)
    assert (shifted['ref_type'], shifted['alt_type']) == ('InFrame_oORF', 'OutOfFrame_oORF')
    assert records(utr, 20, 1, 'T') == []
    made = records(utr, len(utr) - 4, 1, 'A')  # GCC CC to GCA CC makes no codon of note
    assert made == []
    (stopped,) = records(utr, len(utr) - 6, 3, 'TGA')
    assert (stopped['consequence'], stopped['new_stop_distance_to_cds']) == ('uSTOP_gained', 6)


def test_a_deletion_shifts_the_frame_up_to_the_base_before_its_stop_codon() -> None:
    (shifted,) = records(UTR, 13, 2, '')
    assert shifted == FIRST_FRAME_SHIFTED | {'alt_type_length': 24}
    (lost,) = records(UTR, 14, 2, '')
    assert (lost['consequence'], lost['alt_stop'], lost['alt_stop_distance_to_cds']) == ('uSTOP_lost', True, 12)


def test_whole_codons_put_into_a_stop_codon_or_taken_out_of_it_are_its_loss() -> None:
    """Three bases that can only be written inside the stop codon part it, although they shift no frame."""
    (inserted,) = records(UTR, 16, 0, 'CCC')
    assert (inserted['consequence'], inserted['alt_stop'], inserted['frame_with_cds']) == (
        'uSTOP_lost',
        False,
        'outOfFrame',
    )
    (deleted,) = records('GGCACCATGGCTCCGTAAGCGCC', 16, 3, '')
    assert (deleted['consequence'], deleted['alt_stop']) == ('uSTOP_lost', False)


def test_an_insertion_or_a_deletion_that_makes_a_start_codon_is_a_gain() -> None:
    """The distances of a gained start codon are those of the changed UTR."""
    expected = {'consequence': 'uAUG_gained', 'cap_distance_to_start': 6, 'distance_to_cds': 8, 'distance_to_stop': 24}
    for utr, removed, inserted in (('GGCACCAGGCTCC', 0, 'T'), ('GGCACCACTGGCTCC', 1, '')):
        (gained,) = records(utr, 7, removed, inserted)
        assert {key: gained[key] for key in expected} == expected
        assert (gained['type'], gained['kozak_context']) == ('OutOfFrame_oORF', 'ACCATGG')


def test_a_lost_start_codon_whose_frame_runs_into_the_coding_sequence_has_its_stop_there() -> None:
    (lost,) = records('CCATGGCC', 3, 1, 'C')
    assert (lost['consequence'], lost['type'], lost['distance_to_cds'], lost['distance_to_stop']) == (
        'uAUG_lost',
        'InFrame_oORF',
        6,
        18,
    )


def test_a_shifted_frame_that_loses_its_stop_runs_into_the_coding_sequence() -> None:
    (shifted,) = records('CCATGGCTTAACC', 6, 0, 'G')
    assert (shifted['ref_type'], shifted['ref_type_length']) == ('uORF', 9)
    assert (shifted['alt_type'], shifted['alt_type_length']) == ('InFrame_oORF', 27)


def test_a_record_says_whether_its_start_codon_was_seen_translated() -> None:
    frames = utr5.ReadingFrames(UTR, CODING)
    (lost,) = utr5.effects(frames, utr5.Change(7, 1, 'C'), lambda start: start == 6)
    assert lost['evidence'] is True


def test_an_allele_that_changes_nothing_has_no_record() -> None:
    assert records(UTR, 7, 1, 'T') == []


# Where the records are not the plugin's. Each test says what the plugin gives for the same change.


def test_an_insertion_directly_before_a_start_codon_does_not_remove_it() -> None:
    """The plugin reads the codon at its old place, finds the inserted bases there and reports the start lost."""
    assert records(UTR, 6, 0, 'CC') == []


def test_an_insertion_that_can_be_written_before_a_start_codon_does_not_remove_it() -> None:
    """`CA` after the codon's `A` gives the sequence that `AC` before the codon gives.

    The plugin reports the start lost for both ways of writing it, and for the first one a start gained as well.
    """
    assert records(UTR, 7, 0, 'CA') == records(UTR, 6, 0, 'AC') == []


def test_an_insertion_that_leaves_a_start_codon_where_it_stood_does_not_make_it() -> None:
    """`G` between the codon's `T` and `G` gives the sequence that `G` after the codon gives.

    The plugin reports the codon gained. It is the reference's own, and what the insertion does is shift its frame.
    """
    assert records(UTR, 8, 0, 'G') == records(UTR, 9, 0, 'G') == [FIRST_FRAME_SHIFTED]


def test_a_substitution_that_leaves_a_start_codon_where_it_stood_does_not_make_it() -> None:
    """The plugin reports a gain for any start codon that the changed bases reach, the reference's own too."""
    assert records(UTR, 5, 5, 'TATGA') == []


def test_a_replacement_of_unequal_length_that_leaves_a_start_codon_where_it_stood_does_not_make_it() -> None:
    assert records('CCATGCC', 1, 5, 'TATGCA') == []


def test_a_replacement_does_not_keep_a_start_codon_by_one_that_comes_from_further_on() -> None:
    """`CATGC` replaced by `G` takes out the first of two start codons, and the second moves to its place."""
    two = 'GGCATGCATGCC'
    lost = records(two, 2, 5, 'G')
    assert [(record['consequence'], record['cap_distance_to_start']) for record in lost] == [('uAUG_lost', 3)]
    gained = records('GGGATGCC', 2, 1, 'CATGC')
    assert [(record['consequence'], record['cap_distance_to_start']) for record in gained] == [('uAUG_gained', 3)]


def test_gained_and_lost_start_codons_account_for_the_change_in_their_number() -> None:
    """Each start codon of the changed UTR is the reference's or gained, and each of the reference is kept or lost."""
    letters = random.Random(11)  # noqa: S311 (invented sequences, no secret)
    for _ in range(4000):
        utr = ''.join(letters.choice('ATGC'[: letters.randint(3, 4)]) for _ in range(letters.randint(3, 16)))
        position = letters.randrange(len(utr))
        removed = letters.randint(0, min(6, len(utr) - position))
        inserted = ''.join(letters.choice('ATG') for _ in range(letters.randint(0, 6)))
        changed = utr[:position] + inserted + utr[position + removed :]
        found = kinds(utr, position, removed, inserted)
        difference = found.count('uAUG_gained') - found.count('uAUG_lost')
        assert difference == changed.count('ATG') - utr.count('ATG'), (utr, position, removed, inserted)


def test_a_deletion_that_a_record_writes_as_a_replacement_has_the_deletion_s_records() -> None:
    """`ATGG` replaced by `G` takes out `ATG`, as the deletion of the first `ATG` does."""
    repeated = 'GATGATGGG'
    assert repeated[:4] + 'G' + repeated[8:] == repeated[:0] + repeated[3:]
    lost = records(repeated, 0, 3, '')
    assert [(record['consequence'], record['cap_distance_to_start']) for record in lost] == [('uAUG_lost', 1)]
    assert records(repeated, 4, 4, 'G') == lost
    assert records(repeated, 1, 4, 'A') == lost


def test_every_start_codon_that_an_allele_makes_is_a_gain() -> None:
    """The plugin reports the first one alone."""
    gained = records(UTR, 3, 0, 'ATGCATGC')
    assert [(g['consequence'], g['cap_distance_to_start'], g['type']) for g in gained] == [
        ('uAUG_gained', 3, 'OutOfFrame_oORF'),
        ('uAUG_gained', 7, 'uORF'),
    ]


def test_an_allele_at_the_first_bases_of_the_utr_is_not_read_with_its_last_bases() -> None:
    """The plugin looks for a new start codon from two bases before the allele.

    At the UTR's first two bases that reaches its last two, by Perl's negative index. Here the UTR ends with `A` and
    the changed UTR begins `TG`, and the plugin reports a start codon one base before the UTR.
    """
    assert records('CGCCGGCACCCGGA', 0, 1, 'T') == []


def test_a_repeated_start_codon_is_gained_or_lost_once_and_leaves_the_other_frame_alone() -> None:
    """One more or one fewer copy of `ATGC` can be written at either copy.

    The reading frame of a start codon is what follows it, so the copy with the unchanged bases behind it stays and
    the first copy is the one gained or lost. The plugin, given the insertion behind the first copy, reports the
    second copy gained and the first copy's frame shifted.
    """
    repeated = 'GGCACCATGCATGCCGTAAGCG'
    for position in (6, 10, 14):
        assert [(r['consequence'], r['cap_distance_to_start']) for r in records(repeated, position, 0, 'ATGC')] == [
            ('uAUG_gained', 6)
        ]
    for position in (6, 10):
        assert [(r['consequence'], r['cap_distance_to_start']) for r in records(repeated, position, 4, '')] == [
            ('uAUG_lost', 6)
        ]


def test_an_insertion_directly_before_a_stop_codon_shifts_the_frame_or_does_nothing() -> None:
    """The plugin reports the stop lost for both, and names the same stop codon as the next one."""
    assert records(UTR, 15, 0, 'C') == [FIRST_FRAME_SHIFTED]
    assert records(UTR, 15, 0, 'CCC') == []


def test_an_insertion_in_a_stop_codon_is_its_loss() -> None:
    (lost,) = records(UTR, 16, 0, 'C')
    assert (lost['consequence'], lost['alt_stop'], lost['alt_stop_distance_to_cds']) == ('uSTOP_lost', True, 12)


def test_a_change_that_leaves_the_frame_to_end_at_the_stop_codon_s_place_is_no_loss() -> None:
    """The plugin reports every deletion at a stop codon as its loss, and names the stop codon at the same place."""
    # ATG CCG TAA TAA without `GTA` reads ATG CCA TAA, and ATG CTG TAA GTA G without `CTGT` reads ATG AAG TAG.
    assert records('GGCATGCCGTAATAAGC', 8, 3, '') == []
    assert records('GGCATGCTGTAAGTAGC', 6, 4, '') == []
    assert records('GGCATGCAGTGATAGTC', 7, 6, 'GATA') == []


def test_a_deletion_that_leaves_a_stop_codon_at_the_place_is_no_loss() -> None:
    """`TAGAA` without `AG` reads `TAA` where `TAG` stood. The plugin reports every deletion at a stop codon lost."""
    assert records('GGCACCATGGCTCCGTAGAAGCGCC', 16, 2, '') == []


#           0         1         2
#           012345678901234567890123
STOP_RUN = 'CATGCCCCTTTGACCGGCCTAGCC'
# The frame of the start codon at 1 ends at the stop codon `TGA` at 10, behind three `T`.


def test_a_deletion_that_can_be_written_into_a_stop_codon_is_its_loss_however_it_is_written() -> None:
    """Two `T` fewer can be written before the stop codon or with its first base.

    The plugin reports a shifted frame for the first way and a lost stop for the second.
    """
    written = [records(STOP_RUN, position, 2, '') for position in (8, 9)]
    assert written[0] == written[1]
    assert [record['consequence'] for record in written[0]] == ['uSTOP_lost']


def test_whole_codons_that_can_be_written_before_a_stop_codon_keep_the_stop() -> None:
    assert kinds(STOP_RUN, 7, 3, '') == []
    # Three `T` put in behind the stop codon's own `T` can be written before the stop codon as well.
    assert kinds('GGCATGCTTTAAGC', 10, 0, 'TTT') == kinds('GGCATGCTTTAAGC', 7, 0, 'TTT') == []


def test_the_strength_of_a_start_codon_without_a_context_is_empty() -> None:
    """A start codon at the UTR's last three bases has no base behind it in the UTR.

    The plugin then leaves the strength of the start codon that it looked at before, here `Strong`, in the record.
    """
    lost = records('GGCACCATGGCCGCCATG', 5, 13, '')
    assert [(record['cap_distance_to_start'], record.get('kozak_strength')) for record in lost] == [
        (15, None),
        (6, 'Strong'),
    ]


@pytest.mark.parametrize(
    ('utr', 'unit'),
    [(UTR, 'C'), (UTR, 'CA'), ('GGCACCATGCATGCCGTAAGCG', 'ATGC'), (STOP_RUN, 'T'), (STOP_RUN, 'TT'), (UTR, 'GT')],
)
def test_one_changed_sequence_has_one_set_of_records(utr: str, unit: str) -> None:
    """Every way of writing an insertion or a deletion that gives the same sequence gives the same records."""
    by_changed: dict[str, list[list[dict[str, Any]]]] = {}
    for position in range(len(utr) + 1):
        by_changed.setdefault(utr[:position] + unit + utr[position:], []).append(records(utr, position, 0, unit))
        if utr[position : position + len(unit)] == unit:
            deleted = utr[:position] + utr[position + len(unit) :]
            by_changed.setdefault(deleted, []).append(records(utr, position, len(unit), ''))
    assert any(len(written) > 1 for written in by_changed.values())
    for written in by_changed.values():
        assert all(found == written[0] for found in written)


# The transcript's own start codon

#          0         1
#          0123456789012345678
SPLICED = 'GGCGCCGCCATGGCTGAGC'
OWN_START = 9


def test_the_own_start_codon_has_the_context_that_an_upstream_one_would_have() -> None:
    assert utr5.kozak(SPLICED, OWN_START) == ('GCCATGG', 'Strong')


@pytest.mark.parametrize(
    ('first', 'removed', 'inserted', 'expected'),
    [
        (6, 1, 'T', ('TCCATGG', 'Moderate')),  # the base three before the codon
        (12, 1, 'C', ('GCCATGC', 'Moderate')),  # the base after it
        (7, 1, 'T', ('GTCATGG', 'Strong')),  # a base of the context that the strength does not read
        (9, 0, 'T', ('CCTATGG', 'Moderate')),  # an insertion directly before the codon
        (12, 0, 'A', ('GCCATGA', 'Moderate')),  # an insertion directly after it
        (7, 2, '', ('CCGATGG', 'Moderate')),  # a deletion of the two bases before it
    ],
)
def test_an_allele_that_changes_the_context_of_the_own_start_codon_gives_both(
    first: int, removed: int, inserted: str, expected: tuple[str, str]
) -> None:
    assert utr5.own_start_change(SPLICED, OWN_START, first, removed, inserted) == expected


@pytest.mark.parametrize(
    ('first', 'removed', 'inserted'),
    [
        (5, 1, 'T'),  # before the context
        (13, 1, 'A'),  # behind it
        (10, 1, 'C'),  # in the start codon
        (8, 2, 'TT'),  # across the context and the start codon
        (6, 0, 'GCC'),  # a copy of the three bases before the codon, which leaves the same three before it
        (6, 3, ''),  # one of the two `GCC` before the codon
    ],
)
def test_an_allele_that_leaves_the_context_or_touches_the_start_codon_gives_none(
    first: int, removed: int, inserted: str
) -> None:
    assert utr5.own_start_change(SPLICED, OWN_START, first, removed, inserted) is None


def test_the_own_start_codon_is_the_copy_with_the_coding_sequence_behind_it() -> None:
    """Among repeated bases a start codon can stand at the own start codon's old place and be another one.

    A copy of `CATG` before the UTR's last base puts a start codon where the own one stood. The own one is the later
    copy, and the earlier one is an upstream start codon that the allele makes.
    """
    spliced = 'GGCGCGCC' + 'ATGGCTGCT'
    assert utr5.own_start_change(spliced, 8, 7, 0, 'CATG') == ('TGCATGG', 'Moderate')
    (gained,) = records(spliced[:8], 7, 0, 'CATG', coding=spliced[8:])
    assert (gained['consequence'], gained['cap_distance_to_start'], gained['kozak_context']) == (
        'uAUG_gained',
        8,
        'GCCATGC',
    )
    # Without the UTR's own `ATG`, the two that follow it are still the first codons of the coding sequence.
    spliced = 'GCCATG' + 'ATGATGGGG'
    assert utr5.own_start_change(spliced, 6, 3, 3, '') == ('GCCATGA', 'Moderate')
    (lost,) = records(spliced[:6], 3, 3, '', coding=spliced[6:])
    assert (lost['consequence'], lost['cap_distance_to_start']) == ('uAUG_lost', 3)


def test_the_context_with_an_allele_is_the_same_however_the_allele_is_written() -> None:
    """One more `ATG` among three can be written before, between or behind them."""
    spliced = 'GGCGCC' + 'ATGATGATGGCTTAA'
    written = {utr5.own_start_change(spliced, 6, first, 0, 'ATG') for first in (6, 9, 12, 15)}
    assert written == {('ATGATGA', 'Moderate')}
    # One `ATG` fewer leaves the same three bases before the first of them, and the same base behind it.
    assert {utr5.own_start_change(spliced, 6, first, 3, '') for first in (6, 9, 12)} == {None}


def test_a_change_that_leaves_no_three_bases_before_the_own_start_codon_gives_none() -> None:
    assert utr5.own_start_change('GCCATGGCT', 3, 0, 2, '') is None


# The file of translated upstream ORFs

EVIDENCE = [
    'chr\tpos\tgene\tstrand\ttype\tstop_pos\n',
    'chr1\t100\tINVENTED\t+\tuORF\t130\n',
    'chr2\t7\tOTHER\t-\tuORF\t1\n',
]


def test_the_file_of_translated_upstream_orfs_gives_their_start_codons() -> None:
    assert utr5.uorf_evidence(EVIDENCE) == {('chr1', 100), ('chr2', 7)}


@pytest.mark.parametrize(
    ('lines', 'message'),
    [
        (EVIDENCE[1:], 'not with its header'),
        ([EVIDENCE[0], '1\t100\tINVENTED\n'], 'line 2'),
        ([EVIDENCE[0], 'chr1\tabout 100\n'], 'line 2'),
        ([EVIDENCE[0], 'chr1\n'], 'line 2'),
        (EVIDENCE[:1], 'lists none'),
    ],
)
def test_a_file_that_is_not_the_list_of_translated_upstream_orfs_fails(lines: list[str], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        utr5.uorf_evidence(lines)


# Reading the changed UTR's frames from the reference's


def _random_changes(seed: int, count: int) -> Iterator[tuple[str, utr5.Change]]:
    letters = random.Random(seed)  # noqa: S311 (invented sequences, no secret)
    for _ in range(count):
        alphabet = letters.choice(('ATG', 'ATGA', 'ATGC', 'TAG'))
        utr = ''.join(letters.choice(alphabet) for _ in range(letters.randint(0, 24)))
        position = letters.randint(0, len(utr))
        removed = letters.randint(0, min(5, len(utr) - position))
        copied = utr[position : position + letters.randint(1, 5)]
        inserted = letters.choice((copied, ''.join(letters.choice(alphabet) for _ in range(letters.randint(0, 5)))))
        yield utr, utr5.Change(position, removed, inserted)


def test_the_frames_of_a_changed_utr_are_those_of_reading_it_through() -> None:
    """The changed UTR's stop codons and gained start codons are read from the reference's, without reading it."""
    for utr, change in _random_changes(5, 6000):
        changed = utr[: change.position] + change.inserted + utr[change.position + change.removed :]
        compared = utr5._Compared(utr, changed, *utr5._shared(utr, changed, change.position, change.removed))
        frames, through = utr5.ReadingFrames(utr, CODING), utr5.ReadingFrames(changed, CODING)
        after = utr5._ChangedFrames(frames, compared)
        for start in through.starts:
            assert after.stop_in_utr(start) == through.stop_in_utr(start), (utr, change, start)
            assert (after.stop(start), after.type(start)) == (through.stop(start), through.type(start))
        assert after.starts_made() == [start for start in through.starts if not compared.inherited(start)], (
            utr,
            change,
        )


def test_the_start_codons_that_are_passed_over_have_no_record(monkeypatch: pytest.MonkeyPatch) -> None:
    """Only the start codons whose frame the change can meet are looked at, and looking at every one finds no more."""
    found = [
        (utr, change, records(utr, change.position, change.removed, change.inserted))
        for utr, change in _random_changes(6, 6000)
    ]
    monkeypatch.setattr(
        utr5, '_starts_met', lambda frames, compared: [s for s in frames.starts if s < compared.reference_end]
    )
    for utr, change, met in found:
        assert records(utr, change.position, change.removed, change.inserted) == met, (utr, change)
    assert sum(bool(met) for _utr, _change, met in found) > 500
