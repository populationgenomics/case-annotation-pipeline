"""The conversion of an Annotator JSON of structural variants, on invented positions in the Annotator's layout."""

from __future__ import annotations

import copy
import gzip
import json
import pathlib
from collections.abc import Callable, Sequence
from typing import Any

import duckdb
import invented_sv_annotation as invented
import pyarrow as pa
import pytest

from case_annotation import annotation_tables, cli, sv_annotation, sv_sites

# One site of each kind: a deletion and a tandem duplication written as symbolic alleles, a copy-number loss and a
# gain, an insertion, a breakend, an inversion, a sequence of more than 50 bases and two short sequences.
_DELETION = invented.site('chr1\t1000\t.\tT\t<DEL>\t.\t.\tSVTYPE=DEL;END=5000;SVLEN=-4000')
_DUPLICATION = invented.site('chr1\t7000\t.\tT\t<DUP:TANDEM>\t.\t.\tSVTYPE=DUP;END=9000;SVLEN=2000')
_LOSS = invented.site('chr2\t1000\t.\tN\t<DEL>\t.\t.\tSVTYPE=CNV;END=5000;SVLEN=-4000')
_GAIN = invented.site('chr2\t1000\t.\tN\t<DUP>\t.\t.\tSVTYPE=CNV;END=5000;SVLEN=4000')
_INSERTION = invented.site('chr1\t6000\t.\tA\t<INS>\t.\t.\tSVTYPE=INS;END=6000;SVLEN=300')
_BREAKEND = invented.site('chr1\t8000\t.\tG\tG]chr2:900]\t.\t.\tSVTYPE=BND')
_INVERSION = invented.site('chr1\t20000\t.\tC\t<INV>\t.\t.\tSVTYPE=INV;END=30000;SVLEN=10000')
_LONG = invented.site('chr1\t40000\t.\tA' + 'CGT' * 20 + '\tA\t.\t.\tSVTYPE=DEL;END=40060;SVLEN=-60')
_SHORT = invented.site('chr1\t50000\t.\tACGT\tA\t.\t.\tSVTYPE=DEL;END=50003;SVLEN=-3')
_SHORT_INSERTION = invented.site('chr1\t51000\t.\tA\tACGT\t.\t.\tSVTYPE=INS;END=51000;SVLEN=3')
_SITES = (
    _DELETION,
    _DUPLICATION,
    _LOSS,
    _GAIN,
    _INSERTION,
    _BREAKEND,
    _INVERSION,
    _LONG,
    _SHORT,
    _SHORT_INSERTION,
)
_BY_ID = {site.site_id: site for site in _SITES}


def _rows(position: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    row, transcripts, _, _ = sv_annotation.position_rows(position, _BY_ID, invented.HGNC)
    return row, transcripts


def _row(position: dict[str, Any]) -> dict[str, Any]:
    return _rows(position)[0]


def test_a_position_gives_its_site_and_how_the_annotator_treated_it() -> None:
    treated = {site.site_id: _row(invented.position(site)) for site in _SITES}
    assert {(row['variant_type'], row['annotated_as']) for row in treated.values()} == {
        ('deletion', 'structural'),
        ('tandem_duplication', 'structural'),
        ('copy_number_loss', 'structural'),
        ('copy_number_gain', 'structural'),
        ('insertion', 'structural'),
        ('translocation', 'structural'),
        ('inversion', 'structural'),
        ('deletion', 'small_variant'),
        ('insertion', 'small_variant'),
    }
    assert treated[_LONG.site_id]['annotated_as'] == 'structural'
    assert treated[_SHORT.site_id]['annotated_as'] == 'small_variant'
    assert all(
        (row['chrom'], row['pos']) == (_BY_ID[site_id].chrom, _BY_ID[site_id].pos) for site_id, row in treated.items()
    )
    # A position without a match of any source has no value of its columns.
    assert all(row.get('gnomad_sv_id') is None and row.get('clinvar_records') is None for row in treated.values())
    pa.Table.from_pylist(list(treated.values()), schema=sv_annotation.SV)


def test_the_annotated_interval_is_the_annotators_in_the_convention_of_pos_and_end_pos() -> None:
    def interval(site: sv_sites.Site, **variant: int) -> tuple[int, int]:
        row = _row(invented.position(site, variant))
        return row['annotated_pos'], row['annotated_end_pos']

    # The bases concerned run from annotated_pos + 1 to annotated_end_pos, as from pos + 1 to end_pos.
    assert interval(_DELETION) == (1000, 5000)
    assert interval(_INSERTION, end=6000) == (6000, 6000)
    assert interval(_BREAKEND) == (8000, 8000)
    assert interval(_SHORT) == (50000, 50003)
    # A sequence that the Annotator left-aligned starts earlier and still ends at the site's END.
    assert interval(_LONG, begin=39991) == (39990, 40060)
    # A small variant has the interval of its trimmed allele.
    assert interval(_SHORT, begin=50002, end=50003) == (50001, 50003)


def test_the_best_gnomad_match_of_a_record_with_an_extent_overlaps_it_most() -> None:
    entry, overlaps = invented.gnomad_entry, invented.overlaps
    entries = [
        entry('g-common', 900, 5400, 'deletion', 0.3, **overlaps(0.8, 0.9)),
        entry('g-best', 1001, 5000, 'copy_number_loss', 0.001, failedFilter=True, **overlaps(1.0)),
        entry('g-rare', 1001, 4000, 'deletion', None, **overlaps(0.75, 1.0)),
    ]
    row = _row(invented.position(_DELETION, gnomad=entries))
    assert (row['gnomad_sv_id'], row['gnomad_sv_af'], row['gnomad_sv_failed_filter']) == ('g-best', 0.001, True)
    assert (row['gnomad_sv_ac'], row['gnomad_sv_an'], row['gnomad_sv_hom']) == (5, 100, 1)
    assert (row['gnomad_sv_overlap'], row['gnomad_sv_distance']) == (1.0, None)
    assert row['gnomad_sv_max_af'] == 0.3
    assert [m['id'] for m in row['gnomad_sv_matches']] == ['g-best', 'g-common', 'g-rare']
    common = row['gnomad_sv_matches'][1]
    assert (common['pos'], common['end_pos'], common['variant_type']) == (899, 5400, 'deletion')
    assert (common['reciprocal_overlap'], common['annotation_overlap'], common['failed_filter']) == (0.8, 0.9, False)
    assert row['gnomad_sv_matches'][2]['af'] is None


def test_matches_that_overlap_alike_are_ordered_by_frequency_then_id() -> None:
    entry, overlaps = invented.gnomad_entry, invented.overlaps
    entries = [
        entry('g-b', 1001, 5000, 'deletion', 0.01, **overlaps(1.0)),
        entry('g-unknown', 1001, 5000, 'deletion', None, **overlaps(1.0)),
        entry('g-a', 1001, 5000, 'deletion', 0.01, **overlaps(1.0)),
        entry('g-frequent', 1001, 5000, 'deletion', 0.2, **overlaps(1.0)),
    ]
    expected = ['g-frequent', 'g-a', 'g-b', 'g-unknown']
    for order in (entries, entries[::-1]):
        row = _row(invented.position(_DELETION, gnomad=order))
        assert [m['id'] for m in row['gnomad_sv_matches']] == expected


@pytest.mark.parametrize('site', [_INSERTION, _BREAKEND], ids=['an insertion', 'a breakend'])
def test_the_best_gnomad_match_of_a_record_at_one_position_is_the_nearest(site: sv_sites.Site) -> None:
    entry = invented.gnomad_entry
    entries = [
        entry('g-far', site.pos + 301, site.pos + 300, 'insertion', 0.4),
        entry('g-near', site.pos - 9, site.pos - 10, 'insertion', 0.002),
        entry('g-mobile', site.pos - 41, site.pos - 40, 'mobile_element_insertion', 0.1),
    ]
    row = _row(invented.position(site, gnomad=entries))
    assert (row['gnomad_sv_id'], row['gnomad_sv_af']) == ('g-near', 0.002)
    assert (row['gnomad_sv_overlap'], row['gnomad_sv_distance']) == (None, 10)
    assert row['gnomad_sv_max_af'] == 0.4
    assert [m['id'] for m in row['gnomad_sv_matches']] == ['g-near', 'g-mobile', 'g-far']
    assert [m['reciprocal_overlap'] for m in row['gnomad_sv_matches']] == [None, None, None]
    # An entry at the position itself is no bases away.
    at = _row(invented.position(site, gnomad=[entry('g-at', site.pos + 1, site.pos, 'insertion', 0.1)]))
    assert at['gnomad_sv_distance'] == 0


def test_the_other_frequency_sources_give_their_best_match() -> None:
    overlaps = invented.overlaps
    onekg = [
        {'begin': 1001, 'end': 4500, 'id': 'esv1', 'variantType': 'copy_number_loss', 'allAf': 0.5, **overlaps(0.875)},
        {'begin': 1001, 'end': 5000, 'id': 'esv2', 'variantType': 'copy_number_loss', 'allAf': 0.02, 'allAc': 100,
         'allAn': 5008, **overlaps(1.0)},
    ]  # fmt: skip
    cnv = [
        {'begin': 900, 'end': 5000, 'variantId': 'c1', 'variantType': 'copy_number_loss', 'allSf': 0.0001,
         **overlaps(0.9)},
        {'begin': 1001, 'end': 5000, 'variantId': 'c2', 'variantType': 'copy_number_loss', 'allSf': 0.00002,
         'failedFilter': True, **overlaps(0.95)},
    ]  # fmt: skip
    row = _row(invented.position(_LOSS, **{'oneKg': onekg, 'gnomad-cnv': cnv}))
    assert (row['onekg_sv_id'], row['onekg_sv_af'], row['onekg_sv_ac'], row['onekg_sv_an']) == ('esv2', 0.02, 100, 5008)
    assert (row['onekg_sv_overlap'], row['onekg_sv_distance']) == (1.0, None)
    assert (row['gnomad_cnv_id'], row['gnomad_cnv_sf'], row['gnomad_cnv_failed_filter']) == ('c2', 0.00002, True)
    assert row['gnomad_cnv_overlap'] == 0.95
    near = {'begin': 6021, 'end': 6020, 'id': 'esv3', 'variantType': 'insertion', 'allAf': 0.1}
    insertion = _row(invented.position(_INSERTION, oneKg=[near]))
    assert (insertion['onekg_sv_id'], insertion['onekg_sv_overlap'], insertion['onekg_sv_distance']) == (
        'esv3',
        None,
        20,
    )


_SINGLE = 'criteria provided, single submitter'
_PANEL = 'reviewed by expert panel'


def test_the_clinvar_records_of_a_structural_record_are_ordered_by_overlap_before_review() -> None:
    record, overlaps = invented.clinvar_entry, invented.overlaps
    entries = [
        record('VCV-panel', 500, 5000, 'copy_number_loss', 'Pathogenic', _PANEL, '2098-01-01', **overlaps(0.6)),
        record('VCV-later', 1001, 5000, 'deletion', 'Benign', _SINGLE, '2099-06-01', **overlaps(1.0)),
        record('VCV-b', 1001, 5000, 'deletion', 'Likely benign', _SINGLE, '2099-01-01', **overlaps(1.0)),
        record('VCV-a', 1001, 5000, 'deletion', 'Benign', _SINGLE, '2099-01-01', **overlaps(1.0)),
        record('VCV-unclassified', 1001, 5000, 'deletion', 'x', None, '', **overlaps(1.0)),
    ]
    row = _row(invented.position(_DELETION, clinvar=entries))
    assert [r['accession'] for r in row['clinvar_records']] == [
        'VCV-later',
        'VCV-a',
        'VCV-b',
        'VCV-unclassified',
        'VCV-panel',
    ]
    assert (row['clinvar_accession'], row['clinvar_classification'], row['clinvar_review_status']) == (
        'VCV-later',
        'Benign',
        _SINGLE,
    )
    assert (row['clinvar_stars'], row['clinvar_last_evaluated'], row['clinvar_overlap']) == (1, '2099-06-01', 1.0)
    assert row['clinvar_conditions'] == ['An invented condition']
    panel = row['clinvar_records'][-1]
    assert (panel['pos'], panel['end_pos'], panel['stars'], panel['variant_type']) == (499, 5000, 3, 'copy_number_loss')
    assert (panel['reciprocal_overlap'], panel['annotation_overlap']) == (0.6, 0.6)
    assert _row(invented.position(_DELETION, clinvar=entries[::-1]))['clinvar_records'] == row['clinvar_records']


def test_a_small_variant_has_the_alleles_own_gnomad_values_and_clinvar_records() -> None:
    genomes = {'allAf': 0.25, 'allAc': 5, 'allAn': 20, 'allHc': 1, 'failedFilter': True}
    exomes = {'allAf': 0.125, 'allAc': 1, 'allAn': 8, 'allHc': 0}
    record = invented.clinvar_entry
    own = record('VCV-own', 50001, 50003, 'deletion', 'Benign', _SINGLE, '2099-01-01', isAlleleSpecific=True)
    better = record('VCV-better', 50001, 50003, 'deletion', 'Pathogenic', _PANEL, '2098-01-01', isAlleleSpecific=True)
    other = record('VCV-other', 50001, 50001, 'SNV', 'Pathogenic', _PANEL, '2099-01-01')
    variant = {
        'gnomad': genomes,
        'gnomad-exome': exomes,
        'clinvar': [own, other, better],
        'dbsnp': ['rs1'],
        'revel': {},
    }
    row = _row(invented.position(_SHORT, variant))
    assert row['annotated_as'] == 'small_variant'
    assert (row['gnomad_genomes_af'], row['gnomad_genomes_hom'], row['gnomad_genomes_failed_filter']) == (0.25, 1, True)
    assert (row['gnomad_exomes_af'], row['gnomad_exomes_an'], row['gnomad_exomes_failed_filter']) == (0.125, 8, False)
    # The records of the allele itself, the best reviewed first, as in the variants table; none has an overlap.
    assert [r['accession'] for r in row['clinvar_records']] == ['VCV-better', 'VCV-own']
    assert (row['clinvar_accession'], row.get('clinvar_overlap')) == ('VCV-better', None)
    assert row['clinvar_records'][0].get('reciprocal_overlap') is None
    assert row.get('gnomad_sv_id') is None
    pa.Table.from_pylist([row], schema=sv_annotation.SV)


def _transcripts() -> list[dict[str, Any]]:
    transcript = invented.transcript
    return [
        transcript('NM_02.1', 'RefSeq', 'GENEB', '202', 'intron_variant', 'feature_truncation', 'exon_variant',
                   isManeSelect=True, isCanonical=True, proteinId='NP_02.1', exons='2-3/9', introns='2/8',
                   cdnaPos='200-450/3000', cdsPos='100-350/2000', proteinPos='34-117/666',
                   hgvsc='NM_02.1:c.100_350del'),
        transcript('ENST02.1', 'Ensembl', 'GENEB', 'ENSG02', 'copy_number_decrease', 'transcript_ablation',
                   completeOverlap=True, isCanonical=True),
        transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'unidirectional_gene_fusion', 'intron_variant',
                   introns='1/4', bioType='lncRNA', geneFusions=[{'hgnc': 'GENEC'}]),
    ]  # fmt: skip


def test_each_transcript_gives_a_row_with_what_tells_how_the_variant_hits_it() -> None:
    _, rows = _rows(invented.position(_DELETION, {'transcripts': _transcripts()}))
    assert [(r['variant_type'], r['transcript_id']) for r in rows] == [
        ('deletion', t) for t in ('ENST01.1', 'ENST02.1', 'NM_02.1')
    ]
    fused, whole, coding = rows
    assert (coding['source'], coding['gene_symbol'], coding['gene_id'], coding['protein_id']) == (
        'RefSeq',
        'GENEB',
        '202',
        'NP_02.1',
    )
    assert (coding['exon'], coding['intron'], coding['cdna_pos'], coding['protein_pos']) == (
        '2-3/9',
        '2/8',
        '200-450/3000',
        '34-117/666',
    )
    assert (coding['is_mane_select'], coding['is_canonical'], coding['complete_overlap']) == (True, True, False)
    assert (whole['complete_overlap'], whole['exon'], whole['protein_pos'], whole['is_mane_select']) == (
        True,
        None,
        None,
        False,
    )
    assert (fused['biotype'], fused['intron'], fused['protein_id']) == ('lncRNA', '1/4', None)
    # A transcript has the HGNC id of its gene by the genes section, which gives GENEB none.
    assert (fused['hgnc_id'], whole['hgnc_id'], coding['hgnc_id']) == (invented.HGNC['GENEA'], None, None)
    # The terms of a transcript from most to least severe, whatever order the Annotator gave them in.
    assert coding['consequences'] == ['feature_truncation', 'exon_variant', 'intron_variant']
    assert whole['consequences'] == ['transcript_ablation', 'copy_number_decrease']
    assert fused['consequences'] == ['unidirectional_gene_fusion', 'intron_variant']
    # The table has no column of the fields that are left out: the coding position, the impact, the fusion, the HGVS.
    assert not {'cds_pos', 'impact', 'hgvsc'} & set(sv_annotation.SV_TRANSCRIPTS.names)
    assert all(set(r) == set(sv_annotation.TRANSCRIPT_COLUMNS.names) for r in rows)
    pa.Table.from_pylist(rows, schema=sv_annotation.TRANSCRIPT_COLUMNS)


def test_a_site_has_the_hgnc_ids_of_the_genes_of_its_transcripts() -> None:
    transcript = invented.transcript
    row, _ = _rows(invented.position(_DELETION, {'transcripts': _transcripts()}))
    assert row['hgnc_ids'] == [invented.HGNC['GENEA']]
    # Without a transcript, or with transcripts of genes without an id or without a symbol alone, it has none.
    unnamed = transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG03', 'intron_variant') | {'hgnc': None}
    for transcripts in ([], _transcripts()[:2], [unnamed]):
        row, rows = _rows(invented.position(_DELETION, {'transcripts': transcripts}))
        assert row['hgnc_ids'] is None
        assert all(r['hgnc_id'] is None for r in rows)
    assert sv_annotation.COLUMNS.field('hgnc_ids').equals(annotation_tables.VARIANTS.field('hgnc_ids'))


def test_the_transcripts_table_has_the_columns_it_shares_with_a_small_variants_under_one_name_and_type() -> None:
    small = {field.name: field for field in annotation_tables.TRANSCRIPTS}
    larger = {field.name: field for field in sv_annotation.TRANSCRIPT_COLUMNS}
    shared = [name for name in small if name in larger]
    assert {'transcript_id', 'gene_symbol', 'gene_id', 'hgnc_id', 'consequences', 'is_mane_select'} <= set(shared)
    assert all(small[name].equals(larger[name]) for name in shared)
    # They stand in the same order in both, and a structural variant's own fields are the others.
    assert [name for name in larger if name in small] == shared
    assert [name for name in larger if name not in small] == ['variant_type', 'complete_overlap']


def test_the_rows_do_not_depend_on_the_order_of_transcripts_entries_or_terms() -> None:
    entry, overlaps = invented.gnomad_entry, invented.overlaps
    entries = [entry(f'g-{i}', 1001, 5000 - i, 'deletion', 0.01 * i, **overlaps(1 - i / 10)) for i in range(5)]
    expected = _rows(invented.position(_DELETION, {'transcripts': _transcripts()}, gnomad=entries))
    for turn in range(1, 4):
        transcripts = copy.deepcopy(_transcripts())
        transcripts = transcripts[turn:] + transcripts[:turn]
        for held in transcripts:
            held['consequence'].reverse()
        turned = entries[turn:] + entries[:turn]
        assert _rows(invented.position(_DELETION, {'transcripts': transcripts}, gnomad=turned[::-1])) == expected


def test_an_entry_that_a_source_repeats_counts_once_and_a_position_needs_no_cytogenetic_band() -> None:
    entry = invented.gnomad_entry('g-1', 1001, 5000, 'deletion', 0.1, **invented.overlaps(1.0))
    position = invented.position(_DELETION, gnomad=[entry, dict(entry)], mitomap=[{'begin': 1}, {'begin': 1}])
    del position['cytogeneticBand']
    row = _row(position)
    assert [m['id'] for m in row['gnomad_sv_matches']] == ['g-1']


def test_the_sources_that_are_left_out_are_known_and_give_no_column() -> None:
    held = {
        'decipher': [{'begin': 1, 'end': 9000, 'numDeletions': 3}],
        'clingen': [{'begin': 1001, 'end': 5000, 'id': 'nsv1'}],
        'clingenDosageSensitivityMap': [{'begin': 1, 'end': 90000, 'haploinsufficiency': 'x'}],
        'mitomap': [{'begin': 1001, 'end': 5000, 'variantType': 'deletion'}],
    }
    row, _, on_position, _ = sv_annotation.position_rows(invented.position(_DELETION, **held), _BY_ID, invented.HGNC)
    assert on_position == set(held)
    own = {'site_id', 'chrom', 'pos', 'variant_type', 'annotated_as', 'annotated_pos', 'annotated_end_pos', 'hgnc_ids'}
    assert set(row) == own
    assert not any(name.startswith(('decipher', 'clingen', 'mitomap')) for name in sv_annotation.COLUMNS.names)


def test_a_record_of_two_alleles_keeps_each_alleles_matches_and_no_summary() -> None:
    entry, record, overlaps = invented.gnomad_entry, invented.clinvar_entry, invented.overlaps
    loss = _row(
        invented.position(
            _LOSS,
            gnomad=[entry('g-loss', 1001, 5000, 'deletion', 0.5, **overlaps(1.0))],
            clinvar=[
                record('VCV-loss', 1001, 5000, 'copy_number_loss', 'Benign', _SINGLE, '2099-01-01', **overlaps(1.0))
            ],
        )
    )
    gain = _row(invented.position(_GAIN, gnomad=[entry('g-gain', 1001, 5000, 'duplication', 0.01, **overlaps(1.0))]))
    merged = sv_annotation.record_of_alleles([loss, gain])
    assert [m['id'] for m in merged['gnomad_sv_matches']] == ['g-loss', 'g-gain']
    assert [r['accession'] for r in merged['clinvar_records']] == ['VCV-loss']
    assert (merged['annotated_as'], merged['annotated_pos'], merged['annotated_end_pos']) == ('structural', 1000, 5000)
    assert set(merged) == {'annotated_as', 'annotated_pos', 'annotated_end_pos', 'hgnc_ids', *sv_annotation.LISTS}
    # The genes are those of either allele's transcripts, each once.
    named = [loss | {'hgnc_ids': [90002, 90001]}, gain | {'hgnc_ids': [90001, 90003]}]
    assert sv_annotation.record_of_alleles(named)['hgnc_ids'] == [90001, 90002, 90003]
    assert merged['hgnc_ids'] is None
    # Alleles without a match leave no list, and alleles of other intervals no interval.
    unmatched = [_row(invented.position(site)) | {'gnomad_sv_matches': None} for site in (_LOSS, _DUPLICATION)]
    plain = sv_annotation.record_of_alleles(unmatched)
    assert (plain['gnomad_sv_matches'], plain['annotated_pos'], plain['annotated_end_pos']) == (None, None, None)
    assert 'clinvar_records' not in plain
    with pytest.raises(ValueError, match=r"annotated as \['small_variant', 'structural'\]"):
        sv_annotation.record_of_alleles([loss, _row(invented.position(_SHORT))])


def _with(site: sv_sites.Site, variant: dict[str, Any] | None = None, /, **position: object) -> dict[str, Any]:
    return invented.position(site, variant, **position)


def _without(site: sv_sites.Site, key: str) -> dict[str, Any]:
    position = invented.position(site)
    del position[key]
    return position


def _gnomad(**extra: object) -> list[dict[str, Any]]:
    return [invented.gnomad_entry('g-1', 1001, 5000, 'deletion', 0.1, **extra)]


_OVERLAP = invented.overlaps(1.0)
POSITIONS_REFUSED: dict[str, tuple[Callable[[], dict[str, Any]], str]] = {
    'an id of no site': (lambda: _with(_DELETION, id='0' * 32), 'is that of no site of the union'),
    'no id': (lambda: _without(_DELETION, 'id'), 'is that of no site of the union'),
    'another chromosome than the site': (lambda: _with(_DELETION, chromosome='chr2'), 'and its site'),
    'another position than the site': (lambda: _with(_DELETION, position=1001), 'and its site'),
    'another REF than the site': (lambda: _with(_DELETION, refAllele='N'), 'and its site'),
    'another allele than the site': (lambda: _with(_DELETION, altAlleles=['<DUP>']), 'and its site'),
    'two alleles': (lambda: _with(_DELETION, altAlleles=['<DEL>', '<DUP>']), 'and its site'),
    'two variants': (
        lambda: (p := _with(_DELETION)) | {'variants': p['variants'] * 2},
        'not one variant',
    ),
    'no variant': (lambda: _with(_DELETION) | {'variants': []}, 'not one variant'),
    'a record with an extent without svEnd': (lambda: _without(_DELETION, 'svEnd'), 'svEnd is None'),
    'a long sequence without svEnd': (lambda: _without(_LONG, 'svEnd'), 'svEnd is None'),
    'an svEnd that is not the END of the site': (lambda: _with(_DELETION, svEnd=4999), 'svEnd is 4999'),
    'a variant that ends before its svEnd': (lambda: _with(_DELETION, {'end': 1001}), 'the variant ends at 1001'),
    'a small variant with svEnd': (lambda: _with(_SHORT, svEnd=50003), 'annotated as small_variant'),
    'a breakend with svEnd': (lambda: _with(_BREAKEND, svEnd=8000), 'svEnd is 8000'),
    'a breakend annotated elsewhere': (lambda: _with(_BREAKEND, {'begin': 8001, 'end': 8001}), 'annotated from 8001'),
    'no variant type': (lambda: _with(_DELETION, {'variantType': None}), 'no variantType'),
    'a variant type that this module does not know': (
        lambda: _with(_DELETION, {'variantType': 'complex_structural_alteration'}),
        "'complex_structural_alteration' of a record of SVTYPE DEL",
    ),
    'a deletion annotated as an insertion': (lambda: _with(_DELETION, {'variantType': 'insertion'}), 'is none this'),
    'a copy-number loss annotated as a gain': (
        lambda: _with(_LOSS, {'variantType': 'copy_number_gain'}),
        "'copy_number_gain' of a record of SVTYPE CNV",
    ),
    'a short deletion annotated as a copy-number loss': (
        lambda: _with(_SHORT, {'variantType': 'copy_number_loss'}),
        'annotated as small_variant, is none',
    ),
    'a short duplication written as a sequence, annotated as a small variant': (
        lambda: invented.position(_SHORT_DUPLICATION),
        "'tandem_duplication' of a record of SVTYPE DUP",
    ),
    'a symbolic allele that is not annotated as structural': (
        lambda: (
            (p := _with(_DELETION))
            | {'variants': [{k: v for k, v in p['variants'][0].items() if k != 'isStructuralVariant'}]}
        ),
        'is not annotated as a structural variant',
    ),
    'a breakend that is not annotated as structural': (
        lambda: _with(_BREAKEND, {'isStructuralVariant': False}),
        'is not annotated as a structural variant',
    ),
    'a key of a position that this module does not know': (
        lambda: _with(_DELETION, newSource=[]),
        r"keys \['newSource'\] of a position",
    ),
    'samples on a position': (lambda: _with(_DELETION, samples=[{}]), r"keys \['samples'\] of a position"),
    'a key of a variant that this module does not know': (
        lambda: _with(_DELETION, {'newScore': 1}),
        r"keys \['newScore'\] of a variant",
    ),
    'a term outside the severity order': (
        lambda: _with(
            _DELETION, {'transcripts': [invented.transcript('T1', 'Ensembl', 'GENEA', 'E', 'invented_variant')]}
        ),
        'not in the severity order',
    ),
    'a transcript without terms': (
        lambda: _with(_DELETION, {'transcripts': [invented.transcript('T1', 'Ensembl', 'GENEA', 'E')]}),
        'without consequences',
    ),
    'a transcript twice': (
        lambda: _with(
            _DELETION, {'transcripts': [invented.transcript('T1', 'Ensembl', 'GENEA', 'E', 'intron_variant')] * 2}
        ),
        'the transcript T1 twice',
    ),
    'a transcript of a gene that the genes section lacks': (
        lambda: _with(
            _DELETION, {'transcripts': [invented.transcript('T1', 'Ensembl', 'GENEX', 'E', 'intron_variant')]}
        ),
        "the genes section lacks 'GENEX', the gene of the transcript T1",
    ),
    'a transcript without its id': (
        lambda: _with(_DELETION, {'transcripts': [{'source': 'Ensembl', 'consequence': ['intron_variant']}]}),
        'no transcript',
    ),
    'entries that are no list': (lambda: _with(_DELETION, gnomad={'allAf': 0.1}), 'gnomad: not a list of entries'),
    'a gnomAD entry without its overlaps on a deletion': (
        lambda: _with(_DELETION, gnomad=_gnomad()),
        'without its overlaps',
    ),
    'a gnomAD entry without its annotation overlap': (
        lambda: _with(_DELETION, gnomad=_gnomad(reciprocalOverlap=1.0)),
        'without its overlaps',
    ),
    'a gnomAD entry below the threshold': (
        lambda: _with(_DELETION, gnomad=_gnomad(**invented.overlaps(0.49999))),
        'reciprocal overlap of 0.49999, below 0.5',
    ),
    'a duplication matched to a deletion': (
        lambda: _with(_DELETION, gnomad=[invented.gnomad_entry('g-1', 1001, 5000, 'duplication', 0.1, **_OVERLAP)]),
        'an entry of type duplication on a record of type deletion',
    ),
    'a loss matched to a copy-number gain': (
        lambda: _with(
            _GAIN, oneKg=[{'begin': 1001, 'end': 5000, 'id': 'e', 'variantType': 'copy_number_loss', **_OVERLAP}]
        ),
        'oneKg: an entry of type copy_number_loss on a record of type copy_number_gain',
    ),
    'a ClinVar gain matched to a deletion': (
        lambda: _with(
            _DELETION,
            clinvar=[
                invented.clinvar_entry('V', 1001, 5000, 'copy_number_gain', 'Benign', _SINGLE, '2099-01-01', **_OVERLAP)
            ],
        ),
        'clinvar: an entry of type copy_number_gain on a record of type deletion',
    ),
    'a gnomAD CNV entry below the threshold': (
        lambda: _with(_LOSS, **{'gnomad-cnv': [{'begin': 1, 'end': 9, 'variantId': 'c', **invented.overlaps(0.2)}]}),
        'gnomad-cnv: an entry at a reciprocal overlap of 0.2',
    ),
    'an entry with an overlap on an insertion': (
        lambda: _with(_INSERTION, gnomad=_gnomad(**_OVERLAP)),
        'an entry with an overlap on a record of type insertion',
    ),
    'an entry with an overlap on a breakend': (
        lambda: _with(_BREAKEND, oneKg=[{'begin': 1, 'end': 9, 'id': 'e', **_OVERLAP}]),
        'an entry with an overlap on a record of type translocation',
    ),
    'a gnomAD CNV entry on an insertion': (
        lambda: _with(_INSERTION, **{'gnomad-cnv': [{'begin': 6001, 'end': 6000, 'variantId': 'c'}]}),
        'gnomad-cnv: an entry on a record without an extent',
    ),
    'a ClinVar entry without an overlap on an insertion': (
        lambda: _with(
            _INSERTION, clinvar=[invented.clinvar_entry('V', 6001, 6000, 'insertion', 'Benign', _SINGLE, '')]
        ),
        'clinvar: an entry without an overlap',
    ),
    'an entry matched by overlap on a small variant': (
        lambda: _with(_SHORT, gnomad=[invented.gnomad_entry('g-1', 50001, 50003, 'deletion', 0.1, **_OVERLAP)]),
        'on a record annotated as a small variant',
    ),
    "the allele's own gnomAD values on a structural record": (
        lambda: _with(_LONG, {'gnomad': {'allAf': 0.1}}),
        r"\['gnomad'\] of the allele itself",
    ),
    "the allele's own ClinVar records on a structural record": (
        lambda: _with(_DELETION, {'clinvar': []}),
        r"\['clinvar'\] of the allele itself",
    ),
    'a gnomAD entry without its id': (
        lambda: _with(_DELETION, gnomad=[{'begin': 1001, 'end': 5000, 'allAf': 0.1, **_OVERLAP}]),
        'no variantId',
    ),
    'a nearby entry without its place': (
        lambda: _with(_INSERTION, gnomad=[{'variantId': 'g', 'allAf': 0.1}]),
        'no begin',
    ),
}
# A record that the sites accept and that the Annotator is not known to treat as a small variant.
_SHORT_DUPLICATION = invented.site('chr1\t60000\t.\tA\tACGT\t.\t.\tSVTYPE=DUP;END=60000')
_BY_ID[_SHORT_DUPLICATION.site_id] = _SHORT_DUPLICATION


@pytest.mark.parametrize('case', POSITIONS_REFUSED)
def test_a_position_outside_what_the_annotator_is_known_to_write_is_refused(case: str) -> None:
    position, message = POSITIONS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        sv_annotation.position_rows(position(), _BY_ID, invented.HGNC)


def _union(tmp_path: pathlib.Path, sites: Sequence[sv_sites.Site] = _SITES) -> pathlib.Path:
    lines = sv_sites.header_lines(dict.fromkeys(('chr1', 'chr2'), 200_000_000), [site.alt for site in sites])
    with gzip.open(tmp_path / 'union.vcf.gz', 'wt') as vcf:
        vcf.write('\n'.join([*lines, *(site.line() for site in sites)]) + '\n')
    return tmp_path / 'union.vcf.gz'


def _document(*sources: str, **positions: Any) -> dict[str, Any]:  # noqa: ANN401 (keys of an invented position)
    entry = invented.gnomad_entry('g-1', 1001, 5000, 'deletion', 0.1, **_OVERLAP)
    transcripts = {'transcripts': [invented.transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant')]}
    held = [invented.position(site, transcripts) for site in _SITES]
    held[0] = invented.position(_DELETION, transcripts, gnomad=[entry], **positions)
    return {'header': invented.header(*sources), 'positions': held}


def _convert(
    tmp_path: pathlib.Path, document: dict[str, Any], sources: Sequence[str] = invented.ALL, closed: bool = True
) -> annotation_tables.AnnotationInfo:
    return sv_annotation.convert(
        invented.write_json(tmp_path / 'union.json.gz', document, closed),
        _union(tmp_path),
        invented.run_config(*sources),
        invented.BUNDLE,
        invented.CATALOG,
        tmp_path / 'tables',
    )


def test_the_union_converts_into_its_two_tables(tmp_path: pathlib.Path) -> None:
    info = _convert(tmp_path, _document())
    assert info.sources == tuple(sorted(invented.ALL))
    con = duckdb.connect()
    con.read_parquet(str(tmp_path / 'tables' / 'sv' / '*.parquet')).create_view('sv')
    con.read_parquet(str(tmp_path / 'tables' / 'sv_transcripts' / '*.parquet')).create_view('transcripts')
    assert {r[0] for r in con.sql('SELECT site_id FROM sv').fetchall()} == {site.site_id for site in _SITES}
    assert con.sql('SELECT count(*), count(DISTINCT site_id) FROM transcripts').fetchall() == [(len(_SITES),) * 2]
    assert con.sql('SELECT site_id FROM sv WHERE gnomad_sv_id IS NOT NULL').fetchall() == [(_DELETION.site_id,)]
    assert con.sql('SELECT * FROM sv').to_arrow_table().schema.names == sv_annotation.SV.names
    assert con.sql('SELECT * FROM transcripts').to_arrow_table().schema.names == sv_annotation.SV_TRANSCRIPTS.names


def test_the_command_writes_the_tables_and_what_annotated_them(tmp_path: pathlib.Path) -> None:
    invented.write_json(tmp_path / 'union.json.gz', _document())
    (tmp_path / 'run.json').write_text(json.dumps(invented.run_config(*invented.ALL)))
    (tmp_path / 'bundle.json').write_text(json.dumps(invented.BUNDLE))
    command = ['sv-annotation-tables', '--json', str(tmp_path / 'union.json.gz'), '--sites', str(_union(tmp_path))]
    command += ['--run-config', str(tmp_path / 'run.json'), '--bundle-config', str(tmp_path / 'bundle.json')]
    cli.main([*command, '--catalog', invented.CATALOG, '--output', str(tmp_path / 'tables')])
    recorded = annotation_tables.AnnotationInfo.model_validate_json(
        (tmp_path / 'tables' / 'annotation.json').read_text()
    )
    assert recorded == invented.info(*invented.ALL)
    assert sorted(p.name for p in (tmp_path / 'tables').iterdir()) == ['annotation.json', 'sv', 'sv_transcripts']


def _dropping(index: int) -> dict[str, Any]:
    document = _document()
    del document['positions'][index]
    return document


def _repeating(index: int) -> dict[str, Any]:
    document = _document()
    document['positions'].append(document['positions'][index])
    return document


_DBSNP = {'dbsnp': ['rs1']}
DOCUMENTS_REFUSED: dict[str, tuple[Callable[[], dict[str, Any]], Sequence[str], str]] = {
    'a document without one site': (lambda: _dropping(3), invented.ALL, 'sites of the union that the JSON lacks'),
    'a document with a site twice': (lambda: _repeating(2), invented.ALL, f'the site {_SITES[2].site_id} twice'),
    'a document with samples': (
        lambda: _document() | {'header': invented.header(samples=['S1'])},
        invented.ALL,
        'this annotation has samples',
    ),
    'ClinVar records on a position although the run did not name ClinVar': (
        lambda: _document(
            clinvar=[invented.clinvar_entry('V', 1001, 5000, 'deletion', 'Benign', _SINGLE, '2099-01-01', **_OVERLAP)]
        ),
        ('Ensembl', 'gnomad', 'oneKg'),
        r"position keys \['clinvar'\] although the run did not name \['clinvar'\]",
    ),
    'legacy ClinGen records although the run did not name them': (
        lambda: _document(clingen=[{'id': 'nsv1'}]),
        invented.ALL,
        r"position keys \['clingen'\] although the run did not name \['clingen \(legacy\)'\]",
    ),
    'the dosage map although the run did not name ClinGen': (
        lambda: _document(clingenDosageSensitivityMap=[{'begin': 1}]),
        (*invented.ALL, 'clingen (legacy)'),
        r"\['clingenDosageSensitivityMap'\] although the run did not name \['clingen'\]",
    ),
    'MITOMAP entries although the run did not name MITOMAP': (
        lambda: _document(mitomap=[{'begin': 1}]),
        invented.ALL,
        r"position keys \['mitomap'\] although the run did not name \['mitomap'\]",
    ),
    'dbSNP ids on a small variant although the run did not name dbSNP': (
        lambda: (
            _document() | {'positions': [*_document()['positions'][:-1], invented.position(_SHORT_INSERTION, _DBSNP)]}
        ),
        invented.ALL,
        r"variant keys \['dbsnp'\] although the run did not name \['dbSNP'\]",
    ),
    "the allele's ClinVar records although the run did not name ClinVar": (
        lambda: (
            _document()
            | {'positions': [*_document()['positions'][:-1], invented.position(_SHORT_INSERTION, {'clinvar': []})]}
        ),
        ('Ensembl', 'gnomad', 'oneKg'),
        r"variant keys \['clinvar'\] although the run did not name \['clinvar'\]",
    ),
    'a header entry of a dated source that the run did not name': (
        lambda: _document('ClinVar'),
        ('Ensembl', 'gnomad', 'oneKg'),
        r"the JSON header holds \['ClinVar'\], of 'clinvar', which the run did not name",
    ),
}


@pytest.mark.parametrize('case', DOCUMENTS_REFUSED)
def test_a_document_that_is_not_this_unions_annotation_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    document, sources, message = DOCUMENTS_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, document(), sources)


def test_the_named_dated_sources_may_be_in_the_document(tmp_path: pathlib.Path) -> None:
    document = _document(clingen=[{'id': 'nsv1'}], clingenDosageSensitivityMap=[{'begin': 1}], mitomap=[{'begin': 1}])
    document['positions'][-1] = invented.position(_SHORT_INSERTION, _DBSNP)
    info = _convert(tmp_path, document, (*invented.ALL, 'dbSNP', 'clingen', 'clingen (legacy)', 'mitomap'))
    assert {'clingen', 'clingen (legacy)', 'mitomap', 'dbSNP'} <= set(info.sources)


def test_a_document_that_is_cut_short_is_refused(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='the positions are not closed'):
        _convert(tmp_path, _document(), closed=False)


def test_the_severity_order_knows_the_terms_of_structural_variants_and_keeps_the_others_order() -> None:
    order = annotation_tables.CONSEQUENCES
    added = (
        'unidirectional_gene_fusion',
        'inframe_indel',
        'exon_variant',
        'five_prime_duplicated_transcript',
        'three_prime_duplicated_transcript',
        'copy_number_increase',
        'copy_number_decrease',
        'structural_variant',
    )
    assert set(added) <= set(order)
    assert len(set(order)) == len(order)
    # An in-frame replacement of bases ranks with the in-frame terms: below a frameshift, above a missense variant.
    assert annotation_tables.by_severity(['missense_variant', 'inframe_indel', 'frameshift_variant']) == [
        'frameshift_variant',
        'inframe_indel',
        'missense_variant',
    ]
    assert order.index('inframe_deletion') < order.index('inframe_indel') < order.index('missense_variant')
    # The terms that say only that a variant is structural rank below every term that says what it hits.
    assert annotation_tables.worst(['copy_number_decrease', 'intron_variant', 'structural_variant']) == 'intron_variant'
    assert annotation_tables.worst(['feature_truncation', 'exon_variant']) == 'feature_truncation'
