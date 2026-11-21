"""The rule of the gained `AG`: one stretch of an intron, the catalogue's branch points, and an allele on a transcript.

Every sequence, branch point and allele here is invented.
"""

from __future__ import annotations

import gzip
import pathlib
from collections.abc import Sequence

import invented_transcripts
import pysam
import pytest

from case_annotation import ag_gain, branch_points, gene_models, plugin_gates, transcript_models

STRANDS = pytest.mark.parametrize('strand', [1, -1])

# The end of an invented intron: four bases, the branch point, a zone of 21 bases and the acceptor's three. The zone
# holds an `A` at 5 and a `G` at 15, counted from the branch point, and no `AG`.
BEFORE = 'TGCT'
ZONE = 'CTTCATTCCTTTCTGTTTCTT'
STRETCH = BEFORE + 'A' + ZONE + 'CAG'
BRANCH_POINT = len(BEFORE)
# The exon behind the intron.
EXON = 'GTCCATGGCTGCTTAA'


def _gains(
    position: int,
    removed: int,
    inserted: str,
    stretch: str = STRETCH,
    behind: int = 0,
    branch_point: int = BRANCH_POINT,
) -> list[ag_gain.Gain]:
    """The gains of the allele that replaces `removed` bases from a distance behind the branch point."""
    return ag_gain.gains(stretch, branch_point, behind, branch_point + position, removed, inserted)


def _gain(distance: int, to_acceptor: int, *, yag: bool, existing_ag: bool = False) -> ag_gain.Gain:
    high_risk = distance >= 8
    return ag_gain.Gain(
        branch_point_distance=distance,
        acceptor_distance=to_acceptor,
        yag=yag,
        high_risk=high_risk,
        existing_ag=existing_ag,
        score=2 + high_risk + yag + (not existing_ag),
    )


def test_a_substitution_that_makes_an_ag_is_a_gain_with_its_distances_to_the_g() -> None:
    # The `G` put behind the zone's `A`, six bases behind the branch point, with 18 bases of the intron behind it.
    assert _gains(6, 1, 'G') == [_gain(6, -18, yag=True)]
    # The `A` put before the zone's `G`, which stands 15 bases behind the branch point.
    assert _gains(14, 1, 'A') == [_gain(15, -9, yag=True)]


def test_a_gain_is_high_risk_from_eight_bases_behind_the_branch_point() -> None:
    (near,) = _gains(6, 0, 'AG')
    (far,) = _gains(7, 0, 'AG')
    assert (near.branch_point_distance, near.high_risk) == (7, False)
    assert (far.branch_point_distance, far.high_risk) == (8, True)
    # The first stands behind an `A` and the second behind a `T`, which the score counts apart.
    assert far.score - far.yag == near.score - near.yag + 1


def test_a_gain_is_a_yag_behind_a_pyrimidine_alone() -> None:
    (behind_c,) = _gains(14, 1, 'A')
    (behind_g,) = _gains(13, 2, 'GA')
    assert (behind_c.yag, behind_g.yag) == (True, False)
    assert behind_c.score == behind_g.score + 1


def test_a_zone_that_holds_an_ag_in_the_reference_lowers_the_score_by_one() -> None:
    with_ag = STRETCH[: BRANCH_POINT + 2] + 'AG' + STRETCH[BRANCH_POINT + 4 :]
    assert ag_gain.gains(with_ag, BRANCH_POINT, 0, BRANCH_POINT + 14, 1, 'A') == [
        _gain(15, -9, yag=True, existing_ag=True)
    ]
    (without,) = _gains(14, 1, 'A')
    assert without.score == 5


def test_the_zone_begins_behind_the_branch_point() -> None:
    # An `AG` whose `A` is the first base behind the branch point is the first of the zone.
    assert _gains(1, 2, 'AG') == [_gain(2, -22, yag=False)]
    # A `G` behind the branch point makes no `AG` of the zone: its `A` is the branch point itself.
    assert _gains(1, 1, 'G') == []


def test_the_zone_ends_before_the_acceptor_s_three_bases() -> None:
    assert _gains(20, 2, 'AG') == [_gain(21, -3, yag=True)]
    assert _gains(21, 2, 'AG') == []


def test_a_substitution_gains_exactly_the_ags_of_the_zone_that_the_reference_lacks() -> None:
    """Every substitution of one base in the stretch, against a count of the `AG`s in both sequences."""

    def in_zone(sequence: str) -> set[int]:
        return {at for at in range(BRANCH_POINT + 1, len(sequence) - 4) if sequence[at : at + 2] == 'AG'}

    for index, base in enumerate(STRETCH):
        for other in set('ACGT') - {base}:
            changed = STRETCH[:index] + other + STRETCH[index + 1 :]
            expected = in_zone(changed) - in_zone(STRETCH) if index > BRANCH_POINT else set()
            found = ag_gain.gains(STRETCH, BRANCH_POINT, 0, index, 1, other)
            assert {BRANCH_POINT + gain.branch_point_distance - 1 for gain in found} == expected, (index, other)


def test_every_gain_holds_the_rule_s_own_arithmetic() -> None:
    alleles = [(position, removed, inserted) for position in range(1, 24) for removed, inserted in ((1, 'A'), (1, 'G'))]
    alleles += [(position, 0, inserted) for position in range(1, 25) for inserted in ('AG', 'CAGT', 'G', 'A')]
    alleles += [(position, removed, '') for position in range(1, 22) for removed in (1, 2, 3)]
    seen = 0
    for position, removed, inserted in alleles:
        intron = len(STRETCH) - removed + len(inserted)
        for gain in _gains(position, removed, inserted):
            seen += 1
            assert gain.branch_point_distance >= 2
            assert gain.acceptor_distance <= -3
            assert gain.branch_point_distance - gain.acceptor_distance == intron - BRANCH_POINT - 1
            assert gain.high_risk is (gain.branch_point_distance >= ag_gain.HIGH_RISK_FROM)
            assert gain.score == 2 + gain.high_risk + gain.yag + (not gain.existing_ag)
    assert seen >= 40


def test_an_insertion_gains_the_ag_that_it_holds_or_completes() -> None:
    # An `AG` inside the inserted bases, a `G` put behind an `A`, and an `A` put before a `G`.
    assert _gains(10, 0, 'CAGT') == [_gain(12, -16, yag=True)]
    assert _gains(6, 0, 'G') == [_gain(6, -19, yag=True)]
    assert _gains(15, 0, 'A') == [_gain(16, -9, yag=True)]


def test_an_insertion_can_gain_several_ags() -> None:
    assert _gains(10, 0, 'AGCAG') == [_gain(11, -18, yag=True), _gain(14, -15, yag=True)]


def test_a_deletion_gains_the_ag_that_it_joins() -> None:
    # The bases between the zone's `A` and its `G` are removed.
    assert _gains(6, 9, '') == [_gain(6, -9, yag=True)]
    assert _gains(6, 8, '') == []


def _writings(stretch: str, position: int, removed: int, inserted: str) -> list[tuple[int, int, str]]:
    """Every way of writing an insertion or a deletion of the stretch: the place, the bases removed and those put in."""
    changed = stretch[:position] + inserted + stretch[position + removed :]
    return [
        (at, removed, changed[at : at + len(inserted)])
        for at in range(len(stretch) + 1)
        if stretch[:at] + changed[at : at + len(inserted)] + stretch[at + removed :] == changed
    ]


@pytest.mark.parametrize('unit', ['AG', 'A', 'G', 'T', 'TAG', 'CAGT', 'CTT', 'TTCA'])
def test_an_insertion_has_one_answer_however_it_is_written(unit: str) -> None:
    whole = STRETCH + EXON
    several = 0
    for position in range(len(whole) + 1):
        writings = _writings(whole, position, 0, unit)
        answers = {
            tuple(ag_gain.gains(whole, BRANCH_POINT, len(EXON), at, removed, inserted))
            for at, removed, inserted in writings
        }
        assert len(answers) == 1, (position, unit, writings)
        several += len(writings) > 1
    assert several


@pytest.mark.parametrize('removed', [1, 2, 3])
def test_a_deletion_has_one_answer_however_it_is_written(removed: int) -> None:
    whole = 'TGCTAGCTTCATTATTGTTCAGCAG' + 'CAGGTCC'
    several = 0
    for position in range(len(whole) - removed + 1):
        writings = _writings(whole, position, removed, '')
        answers = {tuple(ag_gain.gains(whole, BRANCH_POINT, 7, at, removed, '')) for at, _removed, _ in writings}
        assert len(answers) == 1, (position, writings)
        several += len(writings) > 1
    assert several


def test_an_ag_that_an_insertion_repeats_is_counted_once_at_the_first_copy() -> None:
    repeated = BEFORE + 'A' + 'CTTCTTTCAGTTTCTCTT' + 'CAG'
    # One more `AG` behind the one at 9 and 10: written before it or behind it, it is the copy at 10 that is gained.
    for position in (9, 11):
        assert _gains(position, 0, 'AG', stretch=repeated) == [_gain(10, -13, yag=True, existing_ag=True)]


def test_an_allele_that_can_be_written_before_the_branch_point_gains_nothing() -> None:
    # The branch point and the `G` behind it are repeated: written behind them, the allele would put an `AG` into the
    # zone, and written before the branch point it leaves everything from the branch point on as it was.
    stretch = BEFORE + 'A' + 'GCTTCTTTCCTTTCTCTTTCTT' + 'CAG'
    assert _writings(stretch, BRANCH_POINT + 2, 0, 'AG') == [
        (BRANCH_POINT, 0, 'AG'),
        (BRANCH_POINT + 1, 0, 'GA'),
        (BRANCH_POINT + 2, 0, 'AG'),
    ]
    assert _gains(2, 0, 'AG', stretch=stretch) == []
    assert _gains(0, 0, 'AG', stretch=stretch) == []
    # So does an insertion behind the branch point that ends with the branch point's base: it is `ACAG` before it.
    assert _writings(STRETCH, BRANCH_POINT + 1, 0, 'CAGA')[0] == (BRANCH_POINT, 0, 'ACAG')
    assert _gains(1, 0, 'CAGA') == _gains(0, 0, 'ACAG') == []
    # The same bases without that end can only be written behind the branch point.
    assert _gains(1, 0, 'CAGT') == [_gain(3, -25, yag=True)]


def test_an_allele_that_changes_the_branch_point_or_lies_before_it_gains_nothing() -> None:
    assert _gains(0, 1, 'G') == []
    assert _gains(-2, 0, 'CAGT') == []
    assert _gains(-1, 3, 'AG') == []


def test_an_allele_that_changes_the_exon_gains_nothing() -> None:
    whole = STRETCH + EXON
    # A deletion from the zone into the exon, and a substitution of the exon's first base.
    assert _gains(20, 6, 'AG', stretch=whole, behind=len(EXON)) == []
    assert ag_gain.gains(whole, BRANCH_POINT, len(EXON), len(STRETCH), 1, 'A') == []


def test_an_insertion_at_the_acceptor_that_leaves_one_ag_gains_nothing() -> None:
    """`CAG` with `GGG` put behind its `A` still ends in the one `AG`, with the inserted bases in the exon."""
    whole = STRETCH + EXON
    assert whole[len(STRETCH) - 3 : len(STRETCH) + 1] == 'CAGG'
    assert ag_gain.gains(whole, BRANCH_POINT, len(EXON), len(STRETCH) - 1, 0, 'GGG') == []


def test_an_insertion_at_the_acceptor_that_repeats_its_ag_gains_the_first_copy() -> None:
    whole = STRETCH + EXON
    # `AGTC` behind the acceptor's `C` makes `CAGTCAG`: the first `AG` is new, four bases before the exon's own.
    assert ag_gain.gains(whole, BRANCH_POINT, len(EXON), len(STRETCH) - 2, 0, 'AGTC') == [_gain(24, -4, yag=True)]


def test_the_first_branch_point_is_the_catalogued_one_nearest_the_acceptor() -> None:
    catalogued = ag_gain.BranchPoints(forward=[880, 960, 975], reverse=[1030, 1045, 1120])
    assert catalogued.first(1, 1000, 500) == 975
    assert catalogued.first(-1, 1000, 500) == 1030


def test_a_branch_point_counts_from_4_to_99_bases_before_the_exon() -> None:
    for distance, counts in ((3, False), (4, True), (99, True), (100, False)):
        forward = ag_gain.BranchPoints(forward=[1000 - distance], reverse=[])
        reverse = ag_gain.BranchPoints(forward=[], reverse=[1000 + distance])
        assert (forward.first(1, 1000, 500) is not None) is counts
        assert (reverse.first(-1, 1000, 500) is not None) is counts


def test_a_branch_point_outside_the_intron_or_on_the_other_strand_does_not_count() -> None:
    catalogued = ag_gain.BranchPoints(forward=[970], reverse=[1030])
    assert catalogued.first(1, 1000, 30) == 970
    assert catalogued.first(1, 1000, 29) is None
    assert catalogued.first(-1, 1000, 29) is None
    assert ag_gain.BranchPoints(forward=[], reverse=[970]).first(1, 1000, 500) is None


# The prepared catalogue


def _catalogue(
    tmp_path: pathlib.Path, rows: Sequence[Sequence[str]], *, header: Sequence[str] | None = None
) -> pathlib.Path:
    lines = [f'{branch_points.PREAMBLE}an invented catalogue', '\t'.join(header or branch_points.COLUMNS)]
    lines += ['\t'.join(row) for row in rows]
    path = tmp_path / branch_points.FILE_NAME
    with gzip.open(path, 'wt', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    return path


def _row(chrom: str = 'chr1', pos: str = '500', strand: str = '+') -> list[str]:
    name = f'cBP_{chrom}_{pos}_{strand}_A'
    return [chrom, pos, strand, name, *[''] * len(branch_points.INTRON_COLUMNS), 'cBP_Invented', name]


def test_the_catalogue_gives_each_chromosome_its_branch_points_by_strand(tmp_path: pathlib.Path) -> None:
    rows = [_row(pos='470'), _row(pos='470'), _row(pos='530', strand='-'), _row(chrom='chrX', pos='480')]
    read = ag_gain.read_catalogue(_catalogue(tmp_path, rows))
    assert read.keys() == set(gene_models.CHROMOSOMES)
    assert (read['chr1'].first(1, 500, 100), read['chr1'].first(-1, 500, 100)) == (470, 530)
    assert (read['chrX'].first(1, 500, 100), read['chr2'].first(1, 500, 100)) == (480, None)


def test_a_catalogue_that_is_not_compressed_is_read_as_well(tmp_path: pathlib.Path) -> None:
    plain = tmp_path / 'plain.tsv'
    plain.write_text(gzip.decompress(_catalogue(tmp_path, [_row(pos='470')]).read_bytes()).decode())
    assert ag_gain.read_catalogue(plain)['chr1'].first(1, 500, 100) == 470


CATALOGUES_REFUSED = {
    'another header': ({'header': ('chrom', 'pos', 'strand')}, [['chr1', '500', '+']], 'its header is not that'),
    'a row with another number of fields': ({}, [_row()[:-1]], 'a row of 13 fields'),
    'a chromosome that is none of the tables': ({}, [_row(chrom='chr1_alt')], 'a row without one of chr1 to chrM'),
    'a position that is no whole number': ({}, [_row(pos='5e2')], 'a row without one of chr1 to chrM'),
    'a position of zero': ({}, [_row(pos='0')], 'a row without one of chr1 to chrM'),
    'a strand that is none': ({}, [_row(strand='.')], 'a row without one of chr1 to chrM'),
    'no row': ({}, [], 'no branch point'),
}


@pytest.mark.parametrize('case', CATALOGUES_REFUSED)
def test_a_file_that_is_not_the_prepared_catalogue_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    options, rows, message = CATALOGUES_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        ag_gain.read_catalogue(_catalogue(tmp_path, rows, **options))


def test_a_catalogue_whose_compressed_stream_is_cut_off_is_refused(tmp_path: pathlib.Path) -> None:
    path = _catalogue(tmp_path, [_row(pos=str(pos)) for pos in range(100, 4000)])
    path.write_bytes(path.read_bytes()[:-200])
    with pytest.raises(ValueError, match='not the prepared catalogue of branch points'):
        ag_gain.read_catalogue(path)


# An allele on a transcript


def _invented(strand: int) -> invented_transcripts.Invented:
    """A transcript of two exons around an intron that ends in STRETCH."""
    intron = 'GTAAGT' + 'C' * 30 + STRETCH
    return invented_transcripts.laid_out(('GGCACC' + 'ATGGCT', EXON), (7, 22), strand=strand, introns=(intron,))


def _place(invented: invented_transcripts.Invented, offset: int) -> int:
    """The base of the genome at an offset from the second exon's first base, along the transcript."""
    transcript = invented.transcript
    return transcript.exons[1][0] + offset if transcript.strand == 1 else transcript.exons[1][1] - offset


def _catalogued(invented: invented_transcripts.Invented) -> ag_gain.BranchPoints:
    place = _place(invented, BRANCH_POINT - len(STRETCH))
    return ag_gain.BranchPoints([place], []) if invented.transcript.strand == 1 else ag_gain.BranchPoints([], [place])


def _allele(invented: invented_transcripts.Invented, offset: int, removed: int, inserted: str) -> plugin_gates.Allele:
    """The allele that replaces `removed` bases from an offset from the exon with `inserted`, on the transcript.

    An insertion stands before the base at the offset.
    """
    if invented.transcript.strand == 1:
        start, alt = _place(invented, offset), inserted
    else:
        start, alt = _place(invented, offset + removed - 1), transcript_models.reverse_complement(inserted)
    return plugin_gates.Allele(
        start=start, end=start + removed - 1, ref=invented.genome[start - 1 : start - 1 + removed], alt=alt
    )


def _found(
    tmp_path: pathlib.Path,
    invented: invented_transcripts.Invented,
    allele: plugin_gates.Allele,
    catalogued: ag_gain.BranchPoints | None = None,
) -> list[dict[str, object]]:
    path = invented_transcripts.reference(tmp_path, invented.genome)
    with pysam.FastaFile(str(path)) as reference:
        span = ag_gain.written_span(allele, reference, invented_transcripts.CHROM)
        return ag_gain.of_transcript(invented.transcript, allele, span, reference, catalogued or _catalogued(invented))


@STRANDS
def test_an_allele_on_either_strand_has_the_gain_of_its_intron_with_the_branch_point_s_place(
    tmp_path: pathlib.Path, strand: int
) -> None:
    invented = _invented(strand)
    # The zone's `A` stands five bases behind the branch point, and the allele puts a `G` behind it.
    allele = _allele(invented, BRANCH_POINT + 6 - len(STRETCH), 1, 'G')
    (entry,) = _found(tmp_path, invented, allele)
    assert entry == {'branch_point_pos': _place(invented, BRANCH_POINT - len(STRETCH))} | {
        'branch_point_distance': 6,
        'acceptor_distance': -18,
        'yag': True,
        'high_risk': False,
        'existing_ag': False,
        'score': 4,
    }


@STRANDS
def test_the_gains_of_an_allele_are_in_the_transcript_s_order(tmp_path: pathlib.Path, strand: int) -> None:
    invented = _invented(strand)
    found = _found(tmp_path, invented, _allele(invented, BRANCH_POINT + 10 - len(STRETCH), 0, 'AGCAG'))
    assert [entry['branch_point_distance'] for entry in found] == [11, 14]


@STRANDS
def test_an_intron_without_a_branch_point_on_its_strand_gives_nothing(tmp_path: pathlib.Path, strand: int) -> None:
    invented = _invented(strand)
    allele = _allele(invented, BRANCH_POINT + 6 - len(STRETCH), 1, 'G')
    place = _place(invented, BRANCH_POINT - len(STRETCH))
    other_strand = ag_gain.BranchPoints([], [place]) if strand == 1 else ag_gain.BranchPoints([place], [])
    assert _found(tmp_path, invented, allele, other_strand) == []
    assert _found(tmp_path, invented, allele, ag_gain.BranchPoints([], [])) == []


@STRANDS
def test_an_allele_outside_the_stretch_gives_nothing(tmp_path: pathlib.Path, strand: int) -> None:
    invented = _invented(strand)
    # Before the branch point, in the exon behind the intron, and at the donor's end of the intron.
    for offset in (-len(STRETCH), 3, -len(STRETCH) - 30):
        base = 'A' if invented.genome[_place(invented, offset) - 1] != 'A' else 'C'
        assert _found(tmp_path, invented, _allele(invented, offset, 1, base)) == []


@STRANDS
def test_an_insertion_written_in_the_exon_has_the_answer_of_its_writing_in_the_intron(
    tmp_path: pathlib.Path, strand: int
) -> None:
    """A record writes an insertion at the leftmost place on the genome, which is in the exon on one of the strands."""
    intron = 'GTAAGT' + 'C' * 30 + STRETCH
    exon = 'CAGGTCCATGGCTTAA'
    invented = invented_transcripts.laid_out(('GGCACC' + 'ATGGCT', exon), (7, 22), strand=strand, introns=(intron,))
    # One more `CAG` at the acceptor, written before the acceptor's own and behind the exon's first three bases.
    answers = [_found(tmp_path, invented, _allele(invented, offset, 0, 'CAG')) for offset in (-3, 0, 3)]
    assert answers[0] == answers[1] == answers[2]
    (entry,) = answers[0]
    assert (entry['branch_point_distance'], entry['acceptor_distance'], entry['yag']) == (24, -3, True)


# The ends of two invented introns, each with the exon behind it: one ends in an acceptor's `CAG`, and the other in
# an `A` before an exon that begins with `G`, so that an `AG` lies across the two.
BOUNDARIES = {'an acceptor': (STRETCH, EXON), 'an AG across the boundary': (STRETCH[:-3] + 'CCA', EXON)}


@STRANDS
@pytest.mark.parametrize('boundary', BOUNDARIES)
@pytest.mark.parametrize('unit', ['GTTT', 'G', 'GT', 'AG', 'CAG', 'AGTC', 'TTCA', 'CA', 'A'])
def test_an_insertion_beside_the_exon_has_one_answer_on_a_transcript_however_it_is_written(
    tmp_path: pathlib.Path, strand: int, boundary: str, unit: str
) -> None:
    """Every insertion of a unit from 12 bases before the exon to 8 bases into it, on both strands.

    The stretch that is read for a record ends where the allele's reach ends, so a record written in the intron is
    read with the bases of the exon that it can also be written behind.
    """
    end, exon = BOUNDARIES[boundary]
    invented = invented_transcripts.laid_out(
        ('GGCACC' + 'ATGGCT', exon), (7, 22), strand=strand, introns=('GTAAGT' + 'C' * 30 + end,)
    )
    whole = end + exon
    several = 0
    for offset in range(-12, 9):
        at = len(end) + offset
        changed = whole[:at] + unit + whole[at:]
        writings = [
            (other - len(end), changed[other : other + len(unit)])
            for other in range(len(whole) + 1)
            if whole[:other] + changed[other : other + len(unit)] + whole[other:] == changed
        ]
        answers = [_found(tmp_path, invented, _allele(invented, place, 0, bases)) for place, bases in writings]
        assert all(answer == answers[0] for answer in answers), (offset, writings, answers)
        several += len(writings) > 1
    assert several or unit in ('CA', 'AGTC')


@STRANDS
def test_an_insertion_at_an_ag_across_the_boundary_gains_nothing(tmp_path: pathlib.Path, strand: int) -> None:
    """`GTTT` behind the intron's last `A` is `TTTG` behind the exon's first `G`: the `AG` is the reference's."""
    end, exon = BOUNDARIES['an AG across the boundary']
    invented = invented_transcripts.laid_out(
        ('GGCACC' + 'ATGGCT', exon), (7, 22), strand=strand, introns=('GTAAGT' + 'C' * 30 + end,)
    )
    for offset, bases in ((0, 'GTTT'), (1, 'TTTG'), (2, 'TTGT')):
        assert _found(tmp_path, invented, _allele(invented, offset, 0, bases)) == []


@STRANDS
def test_a_reference_without_the_whole_end_of_an_intron_is_refused(tmp_path: pathlib.Path, strand: int) -> None:
    invented = _invented(strand)
    allele = _allele(invented, BRANCH_POINT + 6 - len(STRETCH), 1, 'G')
    # A reference that ends a few bases behind the allele, inside the stretch that is read for it on either strand.
    cut = invented.genome[: allele.start + (10 if strand == 1 else 3)]
    path = tmp_path / 'short.fa'
    path.write_text(f'>{invented_transcripts.CHROM}\n{cut}\n')
    pysam.faidx(str(path))
    with pysam.FastaFile(str(path)) as reference:
        span = ag_gain.written_span(allele, reference, invented_transcripts.CHROM)
        with pytest.raises(ValueError, match='the reference holds'):
            ag_gain.of_transcript(invented.transcript, allele, span, reference, _catalogued(invented))


@STRANDS
def test_the_span_of_an_insertion_or_deletion_reaches_over_the_bases_that_it_repeats(
    tmp_path: pathlib.Path, strand: int
) -> None:
    invented = invented_transcripts.laid_out(('GGCACC' + 'ATGTTTTTTGCACACACTAA',), (7, 26), strand=strand)
    path = invented_transcripts.reference(tmp_path, invented.genome)
    run = sorted(invented.transcript.genomic(cdna) for cdna in (10, 15))
    repeat = sorted(invented.transcript.genomic(cdna) for cdna in (17, 23))
    with pysam.FastaFile(str(path)) as reference:

        def span(allele: plugin_gates.Allele) -> tuple[int, int]:
            return ag_gain.written_span(allele, reference, invented_transcripts.CHROM)

        # One more base of a run, and one base fewer, wherever in the run the record writes it.
        for cdna in (10, 12, 15):
            assert span(invented.allele(cdna, 0, 'T')) == (run[0] - 1, run[1] + 1)
            assert span(invented.allele(cdna, 1, '')) == (run[0], run[1])
        # One more unit of a repeat of two bases, which the base behind its last unit goes on with.
        assert span(invented.allele(19, 0, 'CA')) == (repeat[0] - 1, repeat[1] + 1)
        # An insertion that repeats nothing, and a substitution.
        alone = invented.allele(12, 0, 'G')
        assert span(alone) == (alone.low, alone.high)
        assert span(invented.allele(12, 1, 'G')) == (invented.transcript.genomic(12),) * 2


def test_the_span_of_a_run_at_the_end_of_a_sequence_stops_there(tmp_path: pathlib.Path) -> None:
    path = tmp_path / 'ends.fa'
    path.write_text(f'>{invented_transcripts.CHROM}\n{"A" * 70}CG{"T" * 70}\n')
    pysam.faidx(str(path))
    with pysam.FastaFile(str(path)) as reference:
        first = plugin_gates.Allele(start=3, end=2, ref='', alt='A')
        last = plugin_gates.Allele(start=140, end=140, ref='T', alt='')
        assert ag_gain.written_span(first, reference, invented_transcripts.CHROM) == (0, 71)
        assert ag_gain.written_span(last, reference, invented_transcripts.CHROM) == (73, 142)
