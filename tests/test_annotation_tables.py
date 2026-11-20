"""The conversion of an Annotator JSON into the union's tables, on invented positions in the Annotator's layout."""

from __future__ import annotations

import copy
import gzip
import json
import pathlib
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pysam
import pytest

from case_annotation import annotation_tables, cli

_CATALOG = '209901'
_BUNDLE = {
    'data': {
        'catalog': _CATALOG,
        'versionSelector': {
            'Ensembl': {'SmallVariant': '1'},
            'RefSeq': {'SmallVariant': '1'},
            'gnomad': {'SmallVariant': '1'},
            'dbSNP': {'SmallVariant': '1'},
            'clinvar': {'SmallVariant': '1', 'StructuralVariant': '1'},
            'omim': {'Gene': '1'},
            'FusionCatcher': {'Gene': '1'},
            'NewSource': {'SmallVariant': '1'},
        },
    }
}


def _run_config(*sources: str) -> dict[str, Any]:
    config = copy.deepcopy(_BUNDLE)
    config['data']['versionSelector'] = {s: _BUNDLE['data']['versionSelector'][s] for s in sources}
    return config


_ALL = ('Ensembl', 'RefSeq', 'gnomad', 'dbSNP', 'clinvar')


def _header(*names: str, samples: list[str] | None = None) -> dict[str, Any]:
    return {
        'annotator': 'Invented Annotator 1.0',
        'creationTime': '2099-01-01 00:00:00',
        'genomeAssembly': 'GRCh38',
        'schemaVersion': 6,
        'dataSources': [{'name': n, 'version': '1', 'releaseDate': '2099-01-01'} for n in names],
        **({'samples': samples} if samples is not None else {}),
    }


def _transcript(tid: str, source: str, gene: str, gene_id: str, *consequences: str, **extra: object) -> dict[str, Any]:
    return {
        'transcript': tid,
        'source': source,
        'bioType': 'mRNA',
        'geneId': gene_id,
        'hgnc': gene,
        'consequence': list(consequences),
        'impact': 'moderate',
        'hgvsc': f'{tid}:c.1A>G',
        **extra,
    }


def _promoterai(transcript: str, score: float, distance: int, **extra: object) -> dict[str, Any]:
    return {'transcriptId': transcript, 'score': score, 'distanceFromTss': distance, **extra}


def _gnomad(coverage: int, **groups: tuple[float, int, int, int]) -> dict[str, Any]:
    """A gnomAD entry: its coverage, and of each group named its frequency, count, number and homozygotes."""
    entry: dict[str, Any] = {'coverage': coverage}
    for group, values in groups.items():
        entry |= dict(zip((f'{group}Af', f'{group}Ac', f'{group}An', f'{group}Hc'), values, strict=True))
    return entry


def _clinvar_record(
    accession: str,
    status: str | None,
    evaluated: str,
    classification: str,
    updated: str = '2099-06-01',
    cited: tuple[str, ...] = (),
    submissions: tuple[tuple[str, ...], ...] = (),
) -> dict[str, Any]:
    """An invented ClinVar record of the allele C>T.

    `cited` are the PubMed ids of its germline classification, and `submissions` those of each of its submissions.
    """
    germline = {
        'classification': classification,
        'reviewStatus': status,
        'dateLastEvaluated': evaluated,
        'conditions': [{'traits': [{'name': {'value': f'condition of {accession}'}}]}],
        **({'pubMedIds': list(cited)} if cited else {}),
    }
    return {
        'accession': accession,
        'version': 1,
        'recordType': 'classified',
        'refAllele': 'C',
        'altAllele': 'T',
        'isAlleleSpecific': True,
        'dateLastUpdated': updated,
        'classifications': {'germlineClassification': germline} if status else {},
        'clinicalAssertions': [
            {'accession': f'SCV{index}'} | ({'pubMedIds': list(ids)} if ids else {})
            for index, ids in enumerate(submissions)
        ],
    }


def _other_allele(record: dict[str, Any], alt: str) -> dict[str, Any]:
    """The invented ClinVar record as the Annotator lists it on an allele that it is not of."""
    return {key: value for key, value in record.items() if key != 'isAlleleSpecific'} | {'altAllele': alt}


def _position(chrom: str, pos: int, ref: str, alt: str, **variant: object) -> dict[str, Any]:
    return {
        'chromosome': chrom,
        'position': pos,
        'refAllele': ref,
        'altAlleles': [alt],
        'variants': [{'vid': 'x', 'chromosome': chrom, 'begin': pos, 'refAllele': ref, 'altAllele': alt, **variant}],
    }


# The genes section of the invented documents: GENEA has an HGNC id, and GENEB, like many real genes, has none. The
# ids are invented, as the genes are.
_GENES: list[dict[str, Any]] = [
    {'name': 'GENEA', 'ensemblGeneId': 'ENSG01', 'ncbiGeneId': '1001', 'hgncId': 90001},
    {'name': 'GENEB', 'ensemblGeneId': 'ENSG02'},
]
# Each symbol of the invented positions' transcripts, as hgnc_ids reads a genes section that lists them all.
_HGNC: dict[str, int | None] = {'GENEA': 90001, 'GENEB': None, 'GENEP': None, 'G': None}


def _position_rows(position: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]], set[str]]:
    return annotation_tables.position_rows(position, _HGNC)


def _rich_position() -> dict[str, Any]:
    """An invented missense allele in an invented gene GENEA, with every source this module reads."""
    return _position(
        'chr2',
        100,
        'A',
        'G',
        variantType='SNV',
        hgvsg='NC_000002.12:g.100A>G',
        dbsnp=['rs1', 'rs2'],
        inLowComplexityRegion=True,
        gnomad=_gnomad(31, all=(0.001, 2, 2000, 0), afr=(0.004, 2, 500, 0), male=(0.002, 2, 1000, 1))
        | {'failedFilter': True},
        transcripts=[
            _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant'),
            _transcript(
                'ENST02.3',
                'Ensembl',
                'GENEA',
                'ENSG01',
                'missense_variant',
                isCanonical=True,
                aminoAcidConservation={'scores': [0.34]},
            ),
            _transcript(
                'NM_01.2', 'RefSeq', 'GENEA', '1001', 'missense_variant', isManeSelect=True, hgvsp='NP_01.2:p.(K1E)'
            ),
            _transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'upstream_gene_variant'),
        ],
        alphaMissense=[
            {
                'transcriptId': 'ENST02.2',
                'pathogenicity': 0.9,
                'classification': 'likely_pathogenic',
                'isIsoform': True,
            },
            {'transcriptId': 'ENST01.1', 'pathogenicity': 0.2, 'classification': 'likely_benign'},
        ],
        spliceAI=[
            {
                'geneId': 'ENSG01',
                'donorGainScore0': 0.1,
                'donorGainDistance0': 12,
                'acceptorLossScore0': 0.4,
                'acceptorLossDistance0': -3,
            }
        ],
        promoterAI=[
            _promoterai('ENST02.2', -0.3, -120, geneId='ENSG01', strand=1),
            _promoterai('ENST03.1', 0.2, 45, geneId='ENSG02', strand=-1),
        ],
        regulatoryRegions=[
            {'id': 'R02', 'type': 'promoter', 'consequence': ['regulatory_region_variant']},
            {'id': 'R01', 'type': 'enhancer', 'consequence': ['regulatory_region_variant']},
        ],
        **{
            'primateAI-3D': [
                {
                    'score': 0.7,
                    'scorePercentile': 0.81,
                    'classification': 'pathogenic',
                    'ensemblTranscriptId': 'ENST02.3',
                    'refSeqTranscriptId': 'NM_01.1',
                }
            ]
        },
        revel={'score': 0.5},
        dannScore=0.9,
        phylopScore=3.1,
        phyloPPrimateScore=0.4,
        gerpScore=4.2,
        clinvar=[
            {
                'accession': 'VCV000000002',
                'refAllele': 'A',
                'altAllele': 'G',
                'isAlleleSpecific': True,
                'dateLastUpdated': '2099-02-01',
                'classifications': {
                    'germlineClassification': {
                        'classification': 'Uncertain significance',
                        'reviewStatus': 'criteria provided, single submitter',
                        'dateLastEvaluated': '2099-01-01',
                        'conditions': [{'traits': [{'name': {'value': 'invented condition'}}]}],
                        'pubMedIds': ['990000002'],
                    }
                },
                'clinicalAssertions': [{'accession': 'SCV000000001', 'pubMedIds': ['990000001']}],
            },
            {
                'accession': 'VCV000000003',
                'refAllele': 'A',
                'altAllele': 'T',
                'dateLastUpdated': '2099-03-01',
                'classifications': {
                    'germlineClassification': {
                        'classification': 'Pathogenic',
                        'reviewStatus': 'reviewed by expert panel',
                        'dateLastEvaluated': '2098-06-01',
                        'conditions': [{'traits': [{'name': {'value': 'another invented condition'}}]}],
                    }
                },
            },
        ],
    )


def test_a_position_gives_its_allele_and_its_transcripts() -> None:
    row, transcripts, keys = _position_rows(_rich_position())
    assert row['variant_id'] == 'chr2-100-A-G'
    assert row['rsids'] == ['rs1', 'rs2']  # each id of its own, in the Annotator's order
    # Two missense transcripts: the MANE Select one is chosen, and its gene's Ensembl id comes from the Ensembl
    # transcripts.
    assert (row['worst_consequence'], row['worst_impact']) == ('missense_variant', 'moderate')
    first = row['genes'][0]
    assert (first['worst_transcript_id'], first['gene_symbol'], first['ensembl_gene_id']) == (
        'NM_01.2',
        'GENEA',
        'ENSG01',
    )
    assert (first['mane_select_transcript_id'], first['mane_select_hgvsp']) == ('NM_01.2', 'NP_01.2:p.(K1E)')
    assert (row['gnomad_genomes_af'], row['gnomad_genomes_failed_filter']) == (0.001, True)
    assert 'gnomad_exomes_af' not in row
    # The table holds the ids as a list of strings and a frequency as a 32-bit float.
    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    assert table.schema.field('rsids').type == pa.list_(pa.string())
    assert table.column('rsids').to_pylist() == [['rs1', 'rs2']]
    assert {table.schema.field(f'gnomad_{kind}_af').type for kind in ('genomes', 'exomes')} == {pa.float32()}
    assert (row['alphamissense_max'], row['spliceai_max'], row['primateai3d_max']) == (0.9, 0.4, 0.7)
    assert (row['phylop'], row['phylop_primate']) == (3.1, 0.4)
    # Only the allele-specific ClinVar record is this allele's.
    assert row['clinvar_accession'] == 'VCV000000002'
    assert (row['clinvar_stars'], row['clinvar_conditions']) == (1, ['invented condition'])
    by_id = {t['transcript_id']: t for t in transcripts}
    assert by_id['ENST02.3']['alphamissense_pathogenicity'] == 0.9  # matched without the version
    assert by_id['NM_01.2']['primateai3d_score'] == 0.7
    assert by_id['NM_01.2']['spliceai_max'] == 0.4  # through the gene symbol
    assert by_id['ENST03.1']['spliceai_max'] is None
    assert by_id['ENST02.3']['promoterai_score'] == -0.3  # matched without the version
    assert {'dbsnp', 'clinvar', 'transcripts'} <= keys


def _spliceai(row: dict[str, Any]) -> dict[str, tuple[float, int]]:
    """The SpliceAI scores a transcript row holds, each with its distance, by their kind."""
    kinds = ('acceptor_gain', 'acceptor_loss', 'donor_gain', 'donor_loss')
    held = {kind: (row[f'spliceai_{kind}'], row[f'spliceai_{kind}_distance']) for kind in kinds}
    # A distance is null with its score, and the highest of the scores is the row's `spliceai_max`.
    assert all((score is None) == (distance is None) for score, distance in held.values())
    scores = [score for score, _ in held.values() if score is not None]
    assert row['spliceai_max'] == (max(scores) if scores else None)
    return {kind: pair for kind, pair in held.items() if pair[0] is not None}


def test_a_transcript_has_each_spliceai_score_of_its_gene_with_that_scores_distance() -> None:
    _, transcripts, _ = _position_rows(_rich_position())
    by_id = {t['transcript_id']: _spliceai(t) for t in transcripts}
    # The entry is of the Ensembl gene, and the RefSeq transcript of the same symbol has its values too.
    of_genea = {'acceptor_loss': (0.4, -3), 'donor_gain': (0.1, 12)}
    assert by_id == {'ENST01.1': of_genea, 'ENST02.3': of_genea, 'NM_01.2': of_genea, 'ENST03.1': {}}
    table = pa.Table.from_pylist(transcripts, schema=annotation_tables.TRANSCRIPTS)
    assert {table.schema.field(f'spliceai_{kind}_distance').type for kind in of_genea} == {pa.int32()}
    assert table.column('spliceai_acceptor_loss_distance').to_pylist() == [-3, -3, -3, None]


def test_equal_spliceai_scores_of_one_kind_give_the_same_distance_in_either_order() -> None:
    entries = [
        {'geneId': 'ENSG01', 'donorLossScore0': 0.3, 'donorLossDistance0': -8},
        {'geneId': 'ENSG01', 'donorLossScore0': 0.3, 'donorLossDistance0': 20},
    ]
    transcripts = [_transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant')]
    held = []
    for spliceai in (entries, entries[::-1]):
        _, (row,), _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts, spliceAI=spliceai))
        held.append(_spliceai(row))
    assert held[0] == held[1] == {'donor_loss': (0.3, 20)}


def test_an_allele_has_its_lowest_and_its_highest_promoterai_score_and_a_transcript_its_own_entry() -> None:
    transcripts = [
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'upstream_gene_variant'),
        _transcript('ENST02.1', 'Ensembl', 'GENEA', 'ENSG01', 'upstream_gene_variant'),
        _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', 'upstream_gene_variant'),
    ]
    entries = [_promoterai('ENST01', -0.4, -250), _promoterai('ENST02', 0.1, 30), _promoterai('ENST09', 0.3, 1000)]
    row, rows, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts, promoterAI=entries))
    # Both ends of the scores, also of an entry whose transcript the allele does not list: a filter on one signed
    # summary would miss the scores of the other sign.
    assert (row['promoterai_min'], row['promoterai_max']) == (-0.4, 0.3)
    held = {r['transcript_id']: (r['promoterai_score'], r['promoterai_tss_distance']) for r in rows}
    assert held == {'ENST01.1': (-0.4, -250), 'ENST02.1': (0.1, 30), 'NM_01.1': (None, None)}
    table = pa.Table.from_pylist(rows, schema=annotation_tables.TRANSCRIPTS)
    assert table.schema.field('promoterai_tss_distance').type == pa.int32()
    assert table.column('promoterai_tss_distance').to_pylist() == [-250, 30, None]

    one, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', promoterAI=entries[:1]))
    assert (one['promoterai_min'], one['promoterai_max']) == (-0.4, -0.4)
    none, _, _ = _position_rows(_position('chr1', 5, 'C', 'T'))
    assert (none['promoterai_min'], none['promoterai_max']) == (None, None)


def test_an_alleles_regulatory_features_are_sorted_by_id_and_hold_the_id_and_the_type_alone() -> None:
    row, _, _ = _position_rows(_rich_position())
    assert row['regulatory_features'] == [{'id': 'R01', 'type': 'enhancer'}, {'id': 'R02', 'type': 'promoter'}]
    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    assert table.schema.field('regulatory_features').type == pa.list_(annotation_tables.REGULATORY_FEATURE)
    assert table.column('regulatory_features').to_pylist() == [row['regulatory_features']]


@pytest.mark.parametrize('variant', [{}, {'regulatoryRegions': []}], ids=['no entry', 'an empty one'])
def test_an_allele_without_regulatory_features_has_no_list_rather_than_an_empty_one(variant: dict[str, Any]) -> None:
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', **variant))
    assert row['regulatory_features'] is None
    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    assert table.column('regulatory_features').to_pylist() == [None]


def test_an_allele_without_a_primate_phylop_score_has_none() -> None:
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', phylopScore=1.5))
    assert (row['phylop'], row['phylop_primate']) == (1.5, None)


def _groups(row: dict[str, Any], kind: str) -> dict[str, tuple[float, int, int, int]]:
    """The groups that a row of `variants` holds of gnomAD's genomes or exomes, each with its four values."""
    held = row[f'gnomad_{kind}_groups']
    return {group: tuple(values.values()) for group, values in held.items() if values is not None}


def test_an_allele_has_gnomads_coverage_and_the_values_of_each_group_that_the_annotator_gives() -> None:
    exomes = _gnomad(52, all=(0.5, 4, 8, 1), nfe=(0.5, 2, 4, 1), remaining=(0, 0, 2, 0), female=(0.25, 1, 4, 0))
    position = _rich_position()
    position['variants'][0]['gnomad-exome'] = exomes
    row, _, _ = _position_rows(position)
    assert (row['gnomad_genomes_coverage'], row['gnomad_exomes_coverage']) == (31, 52)
    # A group has its four values under the names of the all-sample columns, and a group that the entry leaves out
    # has none. The all-sample values are no group: they stay in their own columns.
    assert _groups(row, 'genomes') == {'afr': (0.004, 2, 500, 0), 'male': (0.002, 2, 1000, 1)}
    assert _groups(row, 'exomes') == {'nfe': (0.5, 2, 4, 1), 'remaining': (0, 0, 2, 0), 'female': (0.25, 1, 4, 0)}
    assert row['gnomad_genomes_groups']['male'] == {'af': 0.002, 'ac': 2, 'an': 1000, 'hom': 1}
    assert (row['gnomad_genomes_af'], row['gnomad_exomes_af']) == (0.001, 0.5)

    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    fields = {kind: table.schema.field(f'gnomad_{kind}_groups').type for kind in ('genomes', 'exomes')}
    # The exomes have the groups of the genomes without `ami`.
    assert 'ami' in [f.name for f in fields['genomes']]
    assert [f.name for f in fields['exomes']] == [f.name for f in fields['genomes'] if f.name != 'ami']
    assert {f.type for kind in fields.values() for f in kind} == {annotation_tables.GNOMAD_GROUP}
    assert [(f.name, f.type) for f in annotation_tables.GNOMAD_GROUP] == [
        ('af', pa.float32()),
        ('ac', pa.int32()),
        ('an', pa.int32()),
        ('hom', pa.int32()),
    ]
    assert table.schema.field('gnomad_genomes_coverage').type == pa.int32()
    (written,) = table.select(['gnomad_exomes_groups']).to_pylist()
    assert written['gnomad_exomes_groups']['afr'] is None
    assert written['gnomad_exomes_groups']['female'] == {'af': 0.25, 'ac': 1, 'an': 4, 'hom': 0}


@pytest.mark.parametrize(('key', 'kind'), [('gnomad', 'genomes'), ('gnomad-exome', 'exomes')])
def test_every_field_of_gnomads_groups_is_a_group_under_the_name_of_its_keys(key: str, kind: str) -> None:
    names = [f.name for f in annotation_tables.VARIANTS.field(f'gnomad_{kind}_groups').type]
    # An entry with the four keys of every field's name is accepted and fills every field with its own values.
    given = {name: (index / 100, index, 100, 0) for index, name in enumerate(names)}
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', **{key: _gnomad(20, all=(0.1, 1, 10, 0), **given)}))
    assert _groups(row, kind) == given


def test_an_allele_that_gnomad_lacks_has_neither_coverage_nor_groups() -> None:
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', gnomad=_gnomad(18, all=(0.1, 1, 10, 0))))
    # In the genomes without a group, and not in the exomes at all.
    assert row['gnomad_genomes_coverage'] == 18
    assert _groups(row, 'genomes') == {}
    assert not {'gnomad_exomes_coverage', 'gnomad_exomes_groups'} & set(row)
    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    (written,) = table.select(['gnomad_exomes_coverage', 'gnomad_exomes_groups']).to_pylist()
    assert written == {'gnomad_exomes_coverage': None, 'gnomad_exomes_groups': None}


def test_an_allele_is_in_a_low_complexity_region_only_where_the_annotator_says_so() -> None:
    inside, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', inLowComplexityRegion=True))
    outside, _, _ = _position_rows(_position('chr1', 5, 'C', 'T'))
    # The Annotator writes the flag only where it is true, so an allele without it is outside, not unknown.
    assert (inside['in_low_complexity_region'], outside['in_low_complexity_region']) == (True, False)
    assert annotation_tables.VARIANTS.field('in_low_complexity_region').type == pa.bool_()


def test_a_primateai3d_score_has_the_percentile_and_the_class_of_its_own_entry() -> None:
    transcripts = [
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant'),
        _transcript('ENST02.1', 'Ensembl', 'GENEA', 'ENSG01', 'inframe_deletion'),
        _transcript('ENST03.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant'),
    ]
    entries = [
        {'score': 0.3, 'scorePercentile': 0.4, 'classification': 'benign', 'ensemblTranscriptId': 'ENST01.1'},
        {'score': 0.8, 'scorePercentile': 0.95, 'classification': 'pathogenic', 'ensemblTranscriptId': 'ENST01.1'},
        # An entry with a score alone, as the Annotator writes some.
        {'score': 0.6, 'ensemblTranscriptId': 'ENST02.1', 'indelType': 'deletion'},
    ]
    for primateai3d in (entries, entries[::-1]):
        row, rows, _ = _position_rows(
            _position('chr1', 5, 'C', 'T', transcripts=transcripts, **{'primateAI-3D': primateai3d})
        )
        assert (row['primateai3d_max'], row['primateai3d_max_percentile'], row['primateai3d_max_class']) == (
            0.8,
            0.95,
            'pathogenic',
        )
        held = {
            r['transcript_id']: (r['primateai3d_score'], r['primateai3d_percentile'], r['primateai3d_class'])
            for r in rows
        }
        assert held == {
            'ENST01.1': (0.8, 0.95, 'pathogenic'),
            'ENST02.1': (0.6, None, None),
            'ENST03.1': (None, None, None),
        }

    alone, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', **{'primateAI-3D': entries[2:]}))
    assert (alone['primateai3d_max'], alone['primateai3d_max_percentile'], alone['primateai3d_max_class']) == (
        0.6,
        None,
        None,
    )


def test_an_alphamissense_score_says_whether_it_is_of_an_isoform() -> None:
    _, transcripts, _ = _position_rows(_rich_position())
    held = {t['transcript_id']: t['alphamissense_is_isoform'] for t in transcripts}
    # True where the entry says so, false for an entry without the flag, and nothing for a transcript without an
    # entry.
    assert held == {'ENST02.3': True, 'ENST01.1': False, 'NM_01.2': None, 'ENST03.1': None}
    assert annotation_tables.TRANSCRIPTS.field('alphamissense_is_isoform').type == pa.bool_()


def test_a_transcript_has_the_amino_acid_conservation_scores_as_the_annotator_lists_them() -> None:
    transcripts = [
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant'),
        _transcript(
            'ENST02.1', 'Ensembl', 'GENEA', 'ENSG01', 'inframe_deletion', aminoAcidConservation={'scores': [0.9, 0.12]}
        ),
        _transcript('ENST03.1', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant', aminoAcidConservation={'scores': []}),
    ]
    _, rows, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    held = {r['transcript_id']: r['amino_acid_conservation'] for r in rows}
    # In the Annotator's order, and no list for a transcript without a score rather than an empty one.
    assert held == {'ENST01.1': None, 'ENST02.1': [0.9, 0.12], 'ENST03.1': None}
    table = pa.Table.from_pylist(rows, schema=annotation_tables.TRANSCRIPTS)
    assert table.schema.field('amino_acid_conservation').type == pa.list_(pa.float64())
    assert table.column('amino_acid_conservation').to_pylist() == [None, [0.9, 0.12], None]


_SCORED = {'geneId': 'ENSG01', 'donorGainScore0': 0.5, 'donorGainDistance0': 2}
_CONSERVED = _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant', aminoAcidConservation={})
_OWN = _clinvar_record('VCV1', 'criteria provided, single submitter', '2099-01-01', 'Benign')
_OTHER = _other_allele(_OWN, 'G')


def _without(entry: dict[str, Any], key: str) -> dict[str, Any]:
    return {name: value for name, value in entry.items() if name != key}


ENTRIES_REFUSED = {
    'a SpliceAI entry with a score of a second set': (
        {'spliceAI': [_SCORED | {'donorGainScore1': 0.2, 'donorGainDistance1': 9}]},
        r"chr1:5: a SpliceAI entry with the keys \['donorGainDistance1', 'donorGainScore1'\]",
    ),
    'a SpliceAI entry with a key that is no score': (
        {'spliceAI': [_SCORED | {'transcriptId': 'ENST01.1'}]},
        r"chr1:5: a SpliceAI entry with the keys \['transcriptId'\]",
    ),
    'a SpliceAI score without its distance': (
        {'spliceAI': [_SCORED | {'acceptorGainScore0': 0.2}]},
        'chr1:5: a SpliceAI entry with one of acceptorGainScore0 and acceptorGainDistance0 and not the other',
    ),
    'a SpliceAI distance without its score': (
        {'spliceAI': [_SCORED | {'acceptorLossDistance0': -4}]},
        'chr1:5: a SpliceAI entry with one of acceptorLossScore0 and acceptorLossDistance0 and not the other',
    ),
    'a PromoterAI entry without a score': (
        {'promoterAI': [{'transcriptId': 'ENST01', 'distanceFromTss': 30}]},
        r"chr1:5: a PromoterAI entry without \['score'\]",
    ),
    'a PromoterAI entry without its distance': (
        {'promoterAI': [{'transcriptId': 'ENST01', 'score': 0.1}]},
        r"chr1:5: a PromoterAI entry without \['distanceFromTss'\]",
    ),
    'a regulatory entry without an id': (
        {'regulatoryRegions': [{'type': 'promoter', 'consequence': ['regulatory_region_variant']}]},
        r"chr1:5: a regulatory entry without \['id'\]",
    ),
    'a regulatory entry without a type': (
        {'regulatoryRegions': [{'id': 'R01', 'consequence': ['regulatory_region_variant']}]},
        r"chr1:5: a regulatory entry without \['type'\]",
    ),
    'a gnomAD entry without its coverage': (
        {'gnomad': _without(_gnomad(31, all=(0.1, 1, 10, 0)), 'coverage')},
        r"chr1:5: a gnomad entry without \['coverage'\]",
    ),
    'a gnomAD exomes entry without its coverage': (
        {'gnomad-exome': _without(_gnomad(31, all=(0.1, 1, 10, 0)), 'coverage')},
        r"chr1:5: a gnomad-exome entry without \['coverage'\]",
    ),
    'a gnomAD entry with a group that the columns have no field for': (
        {'gnomad': _gnomad(31, all=(0.1, 1, 10, 0), oth=(0.5, 1, 2, 0))},
        r"chr1:5: a gnomad entry with the keys \['othAc', 'othAf', 'othAn', 'othHc'\]",
    ),
    'a gnomAD exomes entry with a group that only the genomes have': (
        {'gnomad-exome': _gnomad(31, all=(0.1, 1, 10, 0), ami=(0.5, 1, 2, 0))},
        r"chr1:5: a gnomad-exome entry with the keys \['amiAc', 'amiAf', 'amiAn', 'amiHc'\]",
    ),
    'a gnomAD group without one of its four values': (
        {'gnomad': _without(_gnomad(31, all=(0.1, 1, 10, 0), male=(0.2, 1, 5, 0)), 'maleHc')},
        "chr1:5: a gnomad entry with some of the four values of 'male' and not the others",
    ),
    'an amino-acid conservation entry without scores': (
        {'transcripts': [_CONSERVED]},
        r"chr1:5: a transcript aminoAcidConservation entry without \['scores'\]",
    ),
    'a ClinVar record of the allele without the date of its last update': (
        {'clinvar': [_without(_OWN, 'dateLastUpdated')]},
        r"chr1:5: a ClinVar entry without \['dateLastUpdated'\]",
    ),
    'a ClinVar record of another allele without its accession': (
        {'clinvar': [_without(_OTHER, 'accession')]},
        r"chr1:5: a ClinVar entry without \['accession'\]",
    ),
    'a ClinVar record of another allele without its reference allele': (
        {'clinvar': [_without(_OTHER, 'refAllele')]},
        r"chr1:5: a ClinVar entry without \['refAllele'\]",
    ),
    'a ClinVar record of another allele without its alternate allele': (
        {'clinvar': [_without(_OTHER, 'altAllele')]},
        r"chr1:5: a ClinVar entry without \['altAllele'\]",
    ),
}


@pytest.mark.parametrize('case', ENTRIES_REFUSED)
def test_an_entry_that_lacks_what_a_column_reads_or_holds_an_unknown_key_is_refused(case: str) -> None:
    variant, message = ENTRIES_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _position_rows(_position('chr1', 5, 'C', 'T', **variant))


def test_a_transcript_has_its_genes_hgnc_id_and_an_allele_the_ids_of_its_genes() -> None:
    hgnc = {'GENEA': 90001, 'GENEB': None, 'GENEC': 80002, 'GENED': 80002000}
    transcripts = [
        _transcript('ENST01.1', 'Ensembl', 'GENED', 'ENSG04', 'intron_variant'),
        _transcript('ENST02.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant'),
        _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', 'intron_variant'),
        _transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'intron_variant'),
        _transcript('ENST04.1', 'Ensembl', 'GENEC', 'ENSG03', 'intron_variant'),
    ]
    row, rows, _ = annotation_tables.position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts), hgnc)
    # A transcript's id is its symbol's, whichever gene model the transcript is of; a symbol without one has none.
    assert [r['hgnc_id'] for r in rows] == [80002000, 90001, 90001, None, 80002]
    # The allele lists each gene's id once, in ascending order, as 32-bit integers.
    assert row['hgnc_ids'] == [80002, 90001, 80002000]
    variants = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    assert variants.schema.field('hgnc_ids').type == pa.list_(pa.int32())
    assert variants.column('hgnc_ids').to_pylist() == [[80002, 90001, 80002000]]
    assert pa.Table.from_pylist(rows, schema=annotation_tables.TRANSCRIPTS).schema.field('hgnc_id').type == pa.int32()


def test_an_allele_whose_genes_have_no_hgnc_id_has_no_list_rather_than_an_empty_one() -> None:
    transcripts = [_transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'intron_variant')]
    row, rows, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    assert (row['hgnc_ids'], rows[0]['hgnc_id']) == (None, None)
    without_transcripts, _, _ = _position_rows(_position('chr1', 5, 'C', 'T'))
    assert pa.Table.from_pylist([row, without_transcripts], schema=annotation_tables.VARIANTS).column(
        'hgnc_ids'
    ).to_pylist() == [None, None]


def test_a_symbol_has_the_hgnc_id_that_its_records_of_the_genes_section_give() -> None:
    genes = [
        # A gene with one record, one whose two gene models have a record each, of which only one holds the id, a gene
        # without an id, and a symbol that many unnamed genes share.
        {'name': 'GENEA', 'ensemblGeneId': 'ENSG01', 'ncbiGeneId': '1001', 'hgncId': 90001},
        {'name': 'GENEC', 'ensemblGeneId': 'ENSG03'},
        {'name': 'GENEC', 'ncbiGeneId': '1003', 'hgncId': 90003},
        {'name': 'GENEC', 'ensemblGeneId': 'ENSG05', 'hgncId': 90003},
        {'name': 'GENEB', 'ensemblGeneId': 'ENSG02'},
        {'name': 'Y_RNA', 'ensemblGeneId': 'ENSG06'},
        {'name': 'Y_RNA', 'ensemblGeneId': 'ENSG07'},
    ]
    assert annotation_tables.hgnc_ids(genes) == {'GENEA': 90001, 'GENEC': 90003, 'GENEB': None, 'Y_RNA': None}


def test_a_transcript_whose_gene_the_genes_section_lacks_is_refused() -> None:
    transcripts = [_transcript('ENST09.1', 'Ensembl', 'GENEZ', 'ENSG09', 'intron_variant')]
    with pytest.raises(
        ValueError, match=r"chr1:5: the genes section lacks 'GENEZ', the gene of the transcript ENST09\.1"
    ):
        _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))


GENES_REFUSED = {
    'a symbol whose records give two ids': (
        [{'name': 'GENEA', 'hgncId': 90001}, {'name': 'GENEA', 'hgncId': 90002}],
        "'GENEA' two HGNC ids, 90001 and 90002",
    ),
    'an id written as a string': ([{'name': 'GENEA', 'hgncId': 'HGNC:90001'}], 'not a positive integer'),
    'an id of 0': ([{'name': 'GENEA', 'hgncId': 0}], 'not a positive integer'),
    'a boolean for an id': ([{'name': 'GENEA', 'hgncId': True}], 'not a positive integer'),
    'a record without a symbol': ([{'ensemblGeneId': 'ENSG01', 'hgncId': 90001}], 'names no gene'),
    'a record with an empty symbol': ([{'name': '', 'ensemblGeneId': 'ENSG01'}], 'names no gene'),
}


@pytest.mark.parametrize('case', GENES_REFUSED)
def test_a_genes_section_that_cannot_be_right_is_refused(case: str) -> None:
    genes, message = GENES_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        annotation_tables.hgnc_ids(genes)


def test_a_position_without_annotation_still_gives_its_allele() -> None:
    position = _position('chr1', 5, 'C', 'T')
    position['variants'] = []
    row, transcripts, keys = _position_rows(position)
    assert (row['variant_id'], transcripts, keys) == ('chr1-5-C-T', [], set())


@pytest.mark.parametrize('variant', [{}, {'dbsnp': []}], ids=['no dbSNP entry', 'an empty one'])
def test_an_allele_without_dbsnp_ids_has_no_list_rather_than_an_empty_one(variant: dict[str, Any]) -> None:
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', **variant))
    assert row['rsids'] is None
    assert pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS).column('rsids').to_pylist() == [None]


def test_a_transcript_listed_twice_for_one_allele_is_refused() -> None:
    transcripts = [
        _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', 'missense_variant'),
        _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', 'intron_variant'),
    ]
    with pytest.raises(ValueError, match=r'chr1:5: the transcript NM_01\.1 twice'):
        _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))


def test_an_unknown_consequence_fails_rather_than_being_ranked() -> None:
    position = _position('chr1', 5, 'C', 'T', transcripts=[_transcript('T1', 'Ensembl', 'G', 'E', 'invented_variant')])
    with pytest.raises(ValueError, match='not in the severity order'):
        _position_rows(position)


def test_a_transcripts_consequences_run_from_the_most_to_the_least_severe_whatever_the_annotators_order() -> None:
    listed = ['intron_variant', 'splice_region_variant', 'non_coding_transcript_variant']
    ordered = ['splice_region_variant', 'intron_variant', 'non_coding_transcript_variant']
    for consequences in (listed, listed[::-1], ordered):
        position = _position('chr1', 5, 'C', 'T', transcripts=[_transcript('T1', 'Ensembl', 'G', 'E', *consequences)])
        _, (transcript,), _ = _position_rows(position)
        assert transcript['consequences'] == ordered
        # The first is the transcript's most severe consequence, and the order is the one that `sv_transcripts` has.
        assert transcript['worst_consequence'] == transcript['consequences'][0]
        assert transcript['consequences'] == annotation_tables.by_severity(consequences)


def test_a_transcripts_biotype_is_the_annotators_own_word_for_it() -> None:
    transcripts = [
        _transcript('ENST01.1', 'Ensembl', 'G', 'E', 'intron_variant', bioType='lncRNA'),
        _transcript('NM_01.1', 'RefSeq', 'G', 'E', 'intron_variant'),
    ]
    row, rows, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    # Passed through unchanged, and a transcript's alone: an entry of `genes` has none.
    assert {r['transcript_id']: r['biotype'] for r in rows} == {'ENST01.1': 'lncRNA', 'NM_01.1': 'mRNA'}
    assert all('biotype' not in gene for gene in row['genes'])


def test_a_transcript_without_consequences_fails() -> None:
    position = _position('chr1', 5, 'C', 'T', transcripts=[_transcript('T1', 'Ensembl', 'G', 'E')])
    with pytest.raises(ValueError, match='without consequences'):
        _position_rows(position)


def test_a_position_of_two_alleles_fails() -> None:
    position = _position('chr1', 5, 'C', 'T')
    position['altAlleles'] = ['T', 'G']
    with pytest.raises(ValueError, match='2 alleles and 1 variants; the union holds one each'):
        _position_rows(position)


def test_a_position_of_two_variants_fails() -> None:
    position = _position('chr1', 5, 'C', 'T')
    position['variants'] *= 2
    with pytest.raises(ValueError, match='1 alleles and 2 variants; the union holds one each'):
        _position_rows(position)


def test_every_clinvar_record_of_the_allele_is_kept_and_the_best_reviewed_summarises() -> None:
    records = [
        _clinvar_record('VCV3', None, '2099-03-01', 'Pathogenic'),  # no germline classification: last
        _clinvar_record('VCV2', 'criteria provided, single submitter', '2099-01-01', 'Benign'),
        _clinvar_record('VCV1', 'criteria provided, single submitter', '2099-02-01', 'Likely pathogenic'),
        _clinvar_record('VCV4', 'reviewed by expert panel', '2098-01-01', 'Pathogenic'),
    ]
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', clinvar=records))
    assert [r['accession'] for r in row['clinvar_records']] == ['VCV4', 'VCV1', 'VCV2', 'VCV3']
    assert (row['clinvar_accession'], row['clinvar_classification'], row['clinvar_stars']) == ('VCV4', 'Pathogenic', 3)
    assert row['clinvar_records'][1]['conditions'] == ['condition of VCV1']
    assert row['clinvar_records'][3]['stars'] is None
    # The records survive the table's schema.
    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    assert len(table.column('clinvar_records')[0]) == 4


def test_the_summarised_clinvar_record_gives_its_last_update_and_the_pubmed_ids_it_cites() -> None:
    records = [
        _clinvar_record(
            'VCV2', 'criteria provided, single submitter', '2099-01-01', 'Benign', '2099-04-01', ('990001',)
        ),
        _clinvar_record(
            'VCV4',
            'reviewed by expert panel',
            '2098-01-01',
            'Pathogenic',
            '2099-05-01',
            cited=('990003', '990002'),
            submissions=(('990004', '990002'), (), ('990001',)),
        ),
    ]
    for entries in (records, records[::-1]):
        row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', clinvar=entries))
        # The record that the other summary columns are of, wherever it stands among the entries. Its ids are those
        # of its classification and of its submissions, each once.
        assert (row['clinvar_accession'], row['clinvar_last_updated']) == ('VCV4', '2099-05-01')
        assert row['clinvar_pubmed_ids'] == ['990001', '990002', '990003', '990004']
    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    assert table.schema.field('clinvar_pubmed_ids').type == pa.list_(pa.string())
    assert table.column('clinvar_pubmed_ids').to_pylist() == [row['clinvar_pubmed_ids']]

    # A record that cites nothing has no list rather than an empty one.
    uncited = _clinvar_record('VCV1', 'criteria provided, single submitter', '2099-01-01', 'Benign', submissions=((),))
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', clinvar=[uncited]))
    assert (row['clinvar_last_updated'], row['clinvar_pubmed_ids']) == ('2099-06-01', None)


def test_the_clinvar_records_of_other_alleles_at_the_position_are_listed_apart_from_the_alleles_own() -> None:
    own = _clinvar_record('VCV1', 'criteria provided, single submitter', '2099-02-01', 'Uncertain significance')
    others = [
        _other_allele(_clinvar_record('VCV7', 'criteria provided, single submitter', '2099-01-01', 'Benign'), 'G'),
        _other_allele(_clinvar_record('VCV8', 'reviewed by expert panel', '2098-01-01', 'Pathogenic'), 'A'),
    ]
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', clinvar=[others[0], own, others[1]]))
    # The allele's own columns are as without the other records, and the other records are ordered by their review.
    alone, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', clinvar=[own]))
    assert 'clinvar_other_alleles' not in alone
    assert {name: value for name, value in row.items() if name != 'clinvar_other_alleles'} == alone
    assert [r['accession'] for r in row['clinvar_records']] == ['VCV1']
    held = [
        (r['accession'], r['ref_allele'], r['alt_allele'], r['classification'], r['stars'], r['conditions'])
        for r in row['clinvar_other_alleles']
    ]
    assert held == [
        ('VCV8', 'C', 'A', 'Pathogenic', 3, ['condition of VCV8']),
        ('VCV7', 'C', 'G', 'Benign', 1, ['condition of VCV7']),
    ]
    table = pa.Table.from_pylist([row, alone], schema=annotation_tables.VARIANTS)
    assert table.schema.field('clinvar_other_alleles').type == pa.list_(annotation_tables.CLINVAR_OTHER_ALLELE)
    assert table.column('clinvar_other_alleles').to_pylist() == [row['clinvar_other_alleles'], None]
    # An entry of the list has every field of an entry of `clinvar_records`, under the same name and type.
    assert list(annotation_tables.CLINVAR_OTHER_ALLELE)[: len(annotation_tables.CLINVAR_RECORD)] == list(
        annotation_tables.CLINVAR_RECORD
    )

    # An allele that ClinVar has no record of, at a position where it has one of another allele.
    unrecorded, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', clinvar=others))
    assert [r['accession'] for r in unrecorded['clinvar_other_alleles']] == ['VCV8', 'VCV7']
    assert not {'clinvar_accession', 'clinvar_records', 'clinvar_last_updated'} & set(unrecorded)


def test_every_column_that_clinvar_fills_is_one_that_a_run_without_clinvar_leaves_out() -> None:
    with_clinvar, _, _ = _position_rows(_rich_position())
    without, _, _ = _position_rows(_without_clinvar(_rich_position()))
    filled = set(with_clinvar) - set(without)
    assert filled == set(annotation_tables.SOURCE_COLUMNS['clinvar'])
    assert {'clinvar_last_updated', 'clinvar_pubmed_ids', 'clinvar_other_alleles'} <= filled


def test_repeated_entries_of_one_transcript_or_gene_keep_the_highest_score() -> None:
    position = _position(
        'chr1',
        5,
        'C',
        'T',
        transcripts=[
            _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant'),
            _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', 'missense_variant'),
        ],
        alphaMissense=[
            {'transcriptId': 'ENST01.1', 'pathogenicity': 0.3},
            {'transcriptId': 'ENST01.2', 'pathogenicity': 0.8},
        ],
        spliceAI=[
            {'geneId': 'ENSG01', 'donorGainScore0': 0.5, 'donorGainDistance0': 2},
            {'geneId': 'ENSG01', 'acceptorGainScore0': 0.2, 'acceptorGainDistance0': -7},
            {'geneId': 'ENSG01', 'donorGainScore0': 0.3, 'donorGainDistance0': 40},
        ],
    )
    row, transcripts, _ = _position_rows(position)
    by_id = {t['transcript_id']: t for t in transcripts}
    assert by_id['ENST01.1']['alphamissense_pathogenicity'] == 0.8
    assert by_id['ENST01.1']['spliceai_max'] == 0.5
    assert by_id['NM_01.1']['spliceai_max'] == 0.5
    assert row['spliceai_max'] == 0.5
    # Each kind of score keeps its own highest, with the distance of that score and not of a later entry.
    assert _spliceai(by_id['ENST01.1']) == {'donor_gain': (0.5, 2), 'acceptor_gain': (0.2, -7)}
    assert _spliceai(by_id['NM_01.1']) == _spliceai(by_id['ENST01.1'])


def test_a_symbol_of_two_ensembl_genes_leaves_a_refseq_transcripts_gene_open() -> None:
    position = _position(
        'chrX',
        20000,
        'C',
        'T',
        transcripts=[
            _transcript('ENST01.1', 'Ensembl', 'GENEP', 'ENSG01', 'intron_variant'),
            _transcript('ENST02.1', 'Ensembl', 'GENEP', 'ENSG02', 'intron_variant'),
            _transcript('NM_01.1', 'RefSeq', 'GENEP', '1001', 'missense_variant'),
        ],
        spliceAI=[
            {'geneId': 'ENSG01', 'donorGainScore0': 0.1, 'donorGainDistance0': 5},
            {'geneId': 'ENSG01', 'acceptorGainScore0': 0.3, 'acceptorGainDistance0': 11},
            {'geneId': 'ENSG02', 'donorGainScore0': 0.6, 'donorGainDistance0': -40},
        ],
    )
    row, transcripts, _ = _position_rows(position)
    gene = row['genes'][0]
    assert (gene['worst_transcript_id'], gene['gene_symbol'], gene['ensembl_gene_id']) == ('NM_01.1', 'GENEP', None)
    refseq = next(t for t in transcripts if t['source'] == 'RefSeq')
    assert refseq['spliceai_max'] == 0.6
    # The RefSeq transcript has the highest score of each kind over both genes of its symbol.
    assert _spliceai(refseq) == {'donor_gain': (0.6, -40), 'acceptor_gain': (0.3, 11)}


# What no column of `variants` holds: each is of one gene or one transcript of an allele, and a field of an entry.
PER_GENE = (
    'gene_symbol',
    'ensembl_gene_id',
    'worst_transcript_id',
    'mane_select_transcript_id',
    'mane_select_consequence',
    'mane_select_hgvsc',
    'mane_select_hgvsp',
)


def test_no_column_of_variants_is_of_one_gene_or_one_transcript_of_the_allele() -> None:
    row, _, _ = _position_rows(_rich_position())
    assert not set(PER_GENE) & set(row)
    assert not set(PER_GENE) & set(annotation_tables.VARIANTS.names)
    assert set(PER_GENE) <= {field.name for field in annotation_tables.GENE}
    assert {'worst_consequence', 'worst_impact', 'hgnc_ids', 'genes'} <= set(annotation_tables.VARIANTS.names)


def _gene(
    symbol: str | None,
    ensembl_gene_id: str | None,
    hgnc_id: int | None,
    worst: tuple[str, str, str],
    mane: tuple[str, str, str, str | None] | None = None,
) -> dict[str, Any]:
    """An entry of `genes`.

    `worst` is the consequence, the impact and the id of the transcript chosen for the gene, and `mane` the id, the
    consequence and the two HGVS names of its MANE Select transcript.
    """
    worst_consequence, worst_impact, worst_transcript_id = worst
    mane_id, mane_consequence, mane_hgvsc, mane_hgvsp = mane or (None, None, None, None)
    return {
        'gene_symbol': symbol,
        'ensembl_gene_id': ensembl_gene_id,
        'hgnc_id': hgnc_id,
        'worst_consequence': worst_consequence,
        'worst_impact': worst_impact,
        'worst_transcript_id': worst_transcript_id,
        'mane_select_transcript_id': mane_id,
        'mane_select_consequence': mane_consequence,
        'mane_select_hgvsc': mane_hgvsc,
        'mane_select_hgvsp': mane_hgvsp,
    }


def test_each_gene_of_an_allele_has_its_own_mane_select_transcript() -> None:
    transcripts = [
        _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', 'intron_variant', impact='modifier', isManeSelect=True),
        _transcript(
            'NM_02.1', 'RefSeq', 'GENEB', '1002', 'missense_variant', isManeSelect=True, hgvsp='NP_02.1:p.(K1E)'
        ),
        _transcript('ENST02.1', 'Ensembl', 'GENEB', 'ENSG02', 'stop_gained', impact='high'),
    ]
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    # GENEB's most severe consequence is on another transcript than its MANE Select one, and each gene's entry names
    # its own MANE Select transcript with that transcript's consequence and HGVS names.
    assert row['genes'] == [
        _gene(
            'GENEB',
            'ENSG02',
            None,
            ('stop_gained', 'high', 'ENST02.1'),
            ('NM_02.1', 'missense_variant', 'NM_02.1:c.1A>G', 'NP_02.1:p.(K1E)'),
        ),
        _gene(
            'GENEA',
            None,
            90001,
            ('intron_variant', 'modifier', 'NM_01.1'),
            ('NM_01.1', 'intron_variant', 'NM_01.1:c.1A>G', None),
        ),
    ]


def test_a_gene_without_a_mane_select_transcript_has_none_of_its_fields() -> None:
    transcripts = [
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant', isCanonical=True),
        _transcript('NM_02.1', 'RefSeq', 'GENEB', '1002', 'intron_variant', impact='modifier', isManeSelect=True),
    ]
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    first, second = row['genes']
    assert (first['gene_symbol'], first['worst_transcript_id']) == ('GENEA', 'ENST01.1')
    assert {name: first[name] for name in PER_GENE if name.startswith('mane_select_')} == dict.fromkeys(
        (name for name in PER_GENE if name.startswith('mane_select_')), None
    )
    # The MANE Select transcript of the allele's other gene stays with that gene.
    assert (second['gene_symbol'], second['mane_select_transcript_id']) == ('GENEB', 'NM_02.1')
    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    assert table.column('genes').to_pylist() == [row['genes']]


def test_an_allele_in_two_genes_has_an_entry_for_each_with_that_genes_own_summary() -> None:
    # GENEA's most severe consequence is on a transcript other than a MANE Select one, and GENEB's is more severe. Of
    # GENEA's transcripts, one with a less severe consequence has a higher impact than the chosen one, and two are
    # flagged as MANE Select.
    transcripts = [
        _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', 'intron_variant', impact='modifier', isManeSelect=True),
        _transcript(
            'NM_00.1', 'RefSeq', 'GENEA', '1001', 'upstream_gene_variant', impact='modifier', isManeSelect=True
        ),
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant', impact='moderate', isCanonical=True),
        _transcript('ENST02.1', 'Ensembl', 'GENEA', 'ENSG01', 'splice_region_variant', impact='low'),
        _transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'missense_variant'),
    ]
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    assert (row['worst_consequence'], row['worst_impact'], row['hgnc_ids']) == ('missense_variant', 'moderate', [90001])
    assert row['genes'] == [
        _gene('GENEB', 'ENSG02', None, ('missense_variant', 'moderate', 'ENST03.1')),
        _gene(
            'GENEA',
            'ENSG01',
            90001,
            ('splice_region_variant', 'low', 'ENST02.1'),
            ('NM_00.1', 'upstream_gene_variant', 'NM_00.1:c.1A>G', None),
        ),
    ]
    table = pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS)
    assert table.schema.field('genes').type == pa.list_(annotation_tables.GENE)
    assert table.column('genes').to_pylist() == [row['genes']]


def test_a_gene_without_an_hgnc_id_has_an_entry_that_lacks_the_id() -> None:
    transcripts = [_transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'intron_variant', impact='modifier')]
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    assert row['hgnc_ids'] is None
    assert row['genes'] == [_gene('GENEB', 'ENSG02', None, ('intron_variant', 'modifier', 'ENST03.1'))]


GENE_MODELS = {
    'a gene with an HGNC id': 'GENEA',
    'a gene without one, by its symbol': 'GENEB',
}


@pytest.mark.parametrize('case', GENE_MODELS)
def test_the_transcripts_of_both_gene_models_of_one_gene_give_one_entry(case: str) -> None:
    symbol = GENE_MODELS[case]
    transcripts = [
        _transcript('ENST01.1', 'Ensembl', symbol, 'ENSG01', 'intron_variant', impact='modifier'),
        _transcript('NM_01.1', 'RefSeq', symbol, '1001', 'missense_variant', isManeSelect=True),
    ]
    row, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    # The entry is the RefSeq transcript's, and its gene id the Ensembl gene of the symbol.
    assert row['genes'] == [
        _gene(
            symbol,
            'ENSG01',
            _HGNC[symbol],
            ('missense_variant', 'moderate', 'NM_01.1'),
            ('NM_01.1', 'missense_variant', 'NM_01.1:c.1A>G', None),
        )
    ]


def test_an_hgnc_id_names_one_gene_whatever_its_symbols_and_a_transcript_without_a_symbol_its_own_gene() -> None:
    hgnc = {'GENEA': 90001, 'GENEA-OLD': 90001}
    transcripts = [
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant'),
        _transcript('NM_01.1', 'RefSeq', 'GENEA-OLD', '1001', 'intron_variant', impact='modifier'),
        _transcript('ENST08.1', 'Ensembl', 'unnamed', 'ENSG08', 'intron_variant', impact='modifier', hgnc=None),
        _transcript('ENST09.1', 'Ensembl', 'unnamed', 'ENSG09', 'intron_variant', impact='modifier', hgnc=None),
        _transcript('ENST09.2', 'Ensembl', 'unnamed', 'ENSG09', 'upstream_gene_variant', impact='modifier', hgnc=None),
    ]
    row, _, _ = annotation_tables.position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts), hgnc)
    assert row['genes'] == [
        _gene('GENEA', 'ENSG01', 90001, ('missense_variant', 'moderate', 'ENST01.1')),
        _gene(None, 'ENSG08', None, ('intron_variant', 'modifier', 'ENST08.1')),
        _gene(None, 'ENSG09', None, ('intron_variant', 'modifier', 'ENST09.1')),
    ]


# The fields of an entry of `genes` that `transcripts` has a column of the same name for.
SHARED_WITH_TRANSCRIPTS = sorted(
    {field.name for field in annotation_tables.GENE} & set(annotation_tables.TRANSCRIPTS.names)
)

# Invented alleles, each as the transcripts the Annotator lists for it.
ALLELES_IN_GENES = {
    'one gene of both gene models': [
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant', impact='modifier'),
        _transcript('ENST02.3', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant', isCanonical=True),
        _transcript('NM_01.2', 'RefSeq', 'GENEA', '1001', 'missense_variant', isManeSelect=True),
    ],
    'the worst consequence in a gene without a MANE Select transcript': [
        _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', '5_prime_UTR_variant', impact='modifier', isManeSelect=True),
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', '5_prime_UTR_variant', impact='modifier'),
        _transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'splice_donor_variant', impact='high'),
    ],
    'two genes with the same consequence, one of them on its MANE Select transcript': [
        _transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'missense_variant', isCanonical=True),
        _transcript('NM_01.1', 'RefSeq', 'GENEA', '1001', 'missense_variant', isManeSelect=True),
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant'),
    ],
    'two genes that only their transcript ids order': [
        _transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'intron_variant', impact='modifier'),
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant', impact='modifier'),
    ],
    'a symbol of two Ensembl genes, and a gene only RefSeq has': [
        _transcript('ENST01.1', 'Ensembl', 'GENEP', 'ENSG01', 'intron_variant', impact='modifier'),
        _transcript('ENST02.1', 'Ensembl', 'GENEP', 'ENSG02', 'intron_variant', impact='modifier'),
        _transcript('NM_01.1', 'RefSeq', 'GENEP', '1001', 'missense_variant', isManeSelect=True),
        _transcript('NR_02.1', 'RefSeq', 'GENEB', '1002', 'upstream_gene_variant', impact='modifier'),
    ],
    'a transcript without a symbol': [
        _transcript('ENST08.1', 'Ensembl', 'unnamed', 'ENSG08', 'stop_gained', impact='high', hgnc=None),
        _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'downstream_gene_variant', impact='modifier'),
    ],
}


@pytest.mark.parametrize('case', ALLELES_IN_GENES)
def test_the_first_gene_of_an_allele_is_the_one_of_its_worst_consequence(case: str) -> None:
    transcripts = ALLELES_IN_GENES[case]
    row, rows, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts))
    genes = row['genes']
    assert (genes[0]['worst_consequence'], genes[0]['worst_impact']) == (row['worst_consequence'], row['worst_impact'])
    for gene in genes:
        own = [r for r in rows if (r['gene_symbol'], r['hgnc_id']) == (gene['gene_symbol'], gene['hgnc_id'])]
        # The entry's transcript is one of the gene's own rows of transcripts, and the entry says what that row says.
        chosen = [r for r in own if r['transcript_id'] == gene['worst_transcript_id']]
        assert [(r['worst_consequence'], r['impact']) for r in chosen] == [
            (gene['worst_consequence'], gene['worst_impact'])
        ]
        # A field that an entry shares by name with `transcripts` holds what that row holds: one name means one thing.
        assert [{name: r[name] for name in SHARED_WITH_TRANSCRIPTS} for r in chosen] == [
            {name: gene[name] for name in SHARED_WITH_TRANSCRIPTS}
        ]
        # So is its MANE Select transcript, and a gene without such a row has none of its fields.
        mane = [r for r in own if r['is_mane_select'] and r['transcript_id'] == gene['mane_select_transcript_id']]
        named = (gene['mane_select_consequence'], gene['mane_select_hgvsc'], gene['mane_select_hgvsp'])
        assert [(r['worst_consequence'], r['hgvsc'], r['hgvsp']) for r in mane] == (
            [named] if any(r['is_mane_select'] for r in own) else []
        )
        assert mane or (gene['mane_select_transcript_id'], *named) == (None, None, None, None)
    # The genes follow from most to least severe, each gene once, and they are those of the allele's transcripts.
    ranks = [annotation_tables.CONSEQUENCES.index(g['worst_consequence']) for g in genes]
    assert ranks == sorted(ranks)
    symbols = [g['gene_symbol'] for g in genes if g['gene_symbol'] is not None]
    assert len(set(symbols)) == len(symbols)
    assert {g['gene_symbol'] for g in genes} == {r['gene_symbol'] for r in rows}
    assert sorted(g['hgnc_id'] for g in genes if g['hgnc_id'] is not None) == (row['hgnc_ids'] or [])
    # A gene's consequence is the most severe of its own transcripts.
    for gene in genes:
        own = [r['worst_consequence'] for r in rows if r['gene_symbol'] == gene['gene_symbol']]
        assert gene['worst_consequence'] == annotation_tables.worst(own)
    # The order in which the Annotator lists the transcripts changes neither the entries nor their order.
    backwards, _, _ = _position_rows(_position('chr1', 5, 'C', 'T', transcripts=transcripts[::-1]))
    assert backwards['genes'] == genes


def test_an_allele_without_transcripts_has_no_genes_rather_than_an_empty_list() -> None:
    row, rows, _ = _position_rows(_position('chr1', 5, 'C', 'T'))
    assert rows == []
    assert 'genes' not in row
    assert pa.Table.from_pylist([row], schema=annotation_tables.VARIANTS).column('genes').to_pylist() == [None]


def _write_json(
    path: pathlib.Path,
    header: dict[str, Any],
    positions: list[dict[str, Any]],
    *,
    whole: bool = True,
    genes: list[dict[str, Any]] | None = _GENES,
) -> None:
    """An Annotator JSON: whole, or cut short after its positions; with a genes section, or (None) without one."""
    lines = [f'{{"header":{json.dumps(header)},"positions":[']
    lines += [json.dumps(p) + (',' if i < len(positions) - 1 else '') for i, p in enumerate(positions)]
    if whole and genes is not None:
        lines += ['],"genes":[', *(json.dumps(g) + (',' if i < len(genes) - 1 else '') for i, g in enumerate(genes))]
    if whole:
        lines += [']}']
    with gzip.open(path, 'wt') as f:
        f.write('\n'.join(lines) + '\n')


def _write_sites(path: pathlib.Path, alleles: list[tuple[str, int, str, str]]) -> None:
    header = pysam.VariantHeader()
    for chrom in ('chr1', 'chr2', 'chr10'):
        header.contigs.add(chrom, length=1000)
    with pysam.VariantFile(str(path), 'wz', header=header) as vcf:
        for chrom, pos, ref, alt in alleles:
            vcf.write(vcf.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt)))


_ALLELES = [('chr1', 5, 'C', 'T'), ('chr2', 100, 'A', 'G'), ('chr10', 7, 'G', 'GA')]


def _positions() -> list[dict[str, Any]]:
    return [_position('chr1', 5, 'C', 'T'), _rich_position(), _position('chr10', 7, 'G', 'GA')]


def _convert(
    tmp_path: pathlib.Path,
    positions: list[dict[str, Any]] | None = None,
    alleles: list[tuple[str, int, str, str]] | None = None,
    header: dict[str, Any] | None = None,
    sources: tuple[str, ...] = _ALL,
    *,
    whole: bool = True,
    genes: list[dict[str, Any]] | None = _GENES,
) -> annotation_tables.AnnotationInfo:
    _write_json(
        tmp_path / 'union.json.gz', header or _header('ClinVar'), positions or _positions(), whole=whole, genes=genes
    )
    _write_sites(tmp_path / 'union.vcf.gz', alleles or _ALLELES)
    return annotation_tables.convert(
        tmp_path / 'union.json.gz',
        tmp_path / 'union.vcf.gz',
        _run_config(*sources),
        _BUNDLE,
        _CATALOG,
        tmp_path / 'tables',
        workers=2,
        chunk_size=1,
    )


def test_the_union_converts_into_its_two_tables(tmp_path: pathlib.Path) -> None:
    info = _convert(tmp_path)
    con = duckdb.connect()
    variants = con.read_parquet(str(tmp_path / 'tables' / 'variants' / '*.parquet')).select('variant_id').order('1')
    assert [r[0] for r in variants.fetchall()] == ['chr1-5-C-T', 'chr10-7-G-GA', 'chr2-100-A-G']
    transcripts = con.read_parquet(str(tmp_path / 'tables' / 'transcripts' / '*.parquet'))
    assert len(transcripts.fetchall()) == 4
    assert info.sources == tuple(sorted(_ALL))
    assert info.data_sources[0].name == 'ClinVar'
    omitted = {o.source: o.reason for o in info.omitted}
    assert omitted == {'omim': 'dated', 'FusionCatcher': 'not_applicable', 'NewSource': 'unplaced'}


def test_the_unions_tables_hold_the_hgnc_ids_that_the_genes_section_gives(tmp_path: pathlib.Path) -> None:
    # The genes section ends the document, after the positions whose rows hold its ids.
    _convert(tmp_path)
    con = duckdb.connect()
    transcripts = con.read_parquet(str(tmp_path / 'tables' / 'transcripts' / '*.parquet'))
    ids = dict(transcripts.select('transcript_id, hgnc_id').fetchall())
    assert ids == {'ENST01.1': 90001, 'ENST02.3': 90001, 'NM_01.2': 90001, 'ENST03.1': None}
    variants = con.read_parquet(str(tmp_path / 'tables' / 'variants' / '*.parquet'))
    assert dict(variants.select('variant_id, hgnc_ids').fetchall()) == {
        'chr1-5-C-T': None,
        'chr2-100-A-G': [90001],
        'chr10-7-G-GA': None,
    }


def test_each_part_stores_pos_as_differences_and_every_other_column_with_a_dictionary(tmp_path: pathlib.Path) -> None:
    _convert(tmp_path)
    for table, declared in annotation_tables._ENCODINGS.items():
        assert declared['pos'] == 'DELTA_BINARY_PACKED'
        written = 0
        for part in sorted((tmp_path / 'tables' / table).glob('*.parquet')):
            metadata = pq.read_metadata(part)
            # A position without transcripts gives a part of transcripts without rows, which has no pages to encode.
            if not metadata.num_rows:
                continue
            written += 1
            for index in range(metadata.num_columns):
                chunk = metadata.row_group(0).column(index)
                if chunk.path_in_schema in declared:
                    assert not chunk.has_dictionary_page
                    # RLE is the encoding of a chunk's levels.
                    assert set(chunk.encodings) - {'RLE'} == {declared[chunk.path_in_schema]}
                    assert chunk.statistics.has_min_max
                elif chunk.physical_type != 'BOOLEAN':
                    assert chunk.has_dictionary_page, chunk.path_in_schema
        assert written


def test_the_writer_options_name_a_nested_column_by_its_path_and_refuse_a_column_the_schema_lacks() -> None:
    schema = pa.schema([pa.field('pos', pa.int32()), pa.field('rsids', pa.list_(pa.string()))])
    options = annotation_tables.column_options(schema, {'rsids.list.element': 'PLAIN'})
    assert options == {'use_dictionary': ['pos'], 'column_encoding': {'rsids.list.element': 'PLAIN'}}
    with pytest.raises(ValueError, match=r"no columns \['rsids'\] to write without a dictionary"):
        annotation_tables.column_options(schema, {'rsids': 'PLAIN'})


def test_the_command_writes_the_tables_and_what_annotated_them(tmp_path: pathlib.Path) -> None:
    _write_json(tmp_path / 'union.json.gz', _header('ClinVar'), _positions())
    _write_sites(tmp_path / 'union.vcf.gz', _ALLELES)
    (tmp_path / 'run.json').write_text(json.dumps(_run_config(*_ALL)))
    (tmp_path / 'bundle.json').write_text(json.dumps(_BUNDLE))
    args = ['--json', str(tmp_path / 'union.json.gz'), '--sites', str(tmp_path / 'union.vcf.gz')]
    args += ['--run-config', str(tmp_path / 'run.json'), '--bundle-config', str(tmp_path / 'bundle.json')]
    cli.main(['annotation-tables', *args, '--catalog', _CATALOG, '--workers', '1', '--output', str(tmp_path / 'out')])
    info = annotation_tables.AnnotationInfo.model_validate_json((tmp_path / 'out' / 'annotation.json').read_text())
    assert info.catalog == _CATALOG


def _without_clinvar(position: dict[str, Any]) -> dict[str, Any]:
    position['variants'][0].pop('clinvar', None)
    return position


def test_the_genes_are_read_from_their_section_among_others() -> None:
    lines = [
        b'{"chromosome":"chr1"},\n',
        b'],"other":[\n',
        b'{"name":"NOTAGENE","hgncId":90009}\n',
        b'],"genes":[\n',
        b'{"name":"GENEA","hgncId":90001},\r\n',
        b'{"name":"GENEB"}\n',
        b'],"another":[\n',
        b'{"name":"NOTAGENE","hgncId":90009}\n',
        b']}\n',
    ]
    assert annotation_tables.read_genes(iter(lines)) == [{'name': 'GENEA', 'hgncId': 90001}, {'name': 'GENEB'}]


def test_a_document_without_a_genes_section_converts_if_no_transcript_names_a_gene(tmp_path: pathlib.Path) -> None:
    _convert(tmp_path, positions=[_position(*allele) for allele in _ALLELES], genes=None)
    con = duckdb.connect()
    variants = con.read_parquet(str(tmp_path / 'tables' / 'variants' / '*.parquet'))
    assert sorted(variants.select('variant_id, hgnc_ids').fetchall()) == [
        ('chr1-5-C-T', None),
        ('chr10-7-G-GA', None),
        ('chr2-100-A-G', None),
    ]


REFUSED = {
    'a union allele the JSON lacks': (
        {'positions': _positions()[:2]},
        'union alleles the JSON lacks',
    ),
    'a JSON allele the union lacks': (
        {'alleles': _ALLELES[:2]},
        'JSON alleles the union lacks',
    ),
    'a repeated union allele': (
        {'alleles': [*_ALLELES, _ALLELES[0]]},
        'the union repeats alleles',
    ),
    'a repeated JSON allele': (
        {'positions': [*_positions(), _position('chr1', 5, 'C', 'T')]},
        'the JSON repeats alleles',
    ),
    'an annotation with samples': (
        {'header': _header(samples=['S1'])},
        'has samples',
    ),
    'a document cut short': (
        {'whole': False},
        'cut short',
    ),
    'a document without a genes section, whose transcripts name genes': (
        {'genes': None},
        "the genes section lacks 'GENE[AB]'",
    ),
    'a genes section that gives a symbol two HGNC ids': (
        {'genes': [*_GENES, {'name': 'GENEA', 'ncbiGeneId': '1001', 'hgncId': 90002}]},
        'two HGNC ids',
    ),
    "a genes section without a transcript's gene": (
        {'genes': _GENES[:1]},
        "the genes section lacks 'GENEB'",
    ),
    'ClinVar records in a run without ClinVar': (
        {'sources': ('Ensembl', 'RefSeq', 'gnomad', 'dbSNP'), 'header': _header()},
        "'clinvar' records",
    ),
    'dbSNP ids in a run without dbSNP': (
        {'sources': ('Ensembl', 'RefSeq', 'gnomad', 'clinvar')},
        "'dbsnp' records although the run did not name 'dbSNP'",
    ),
    "a header holding OMIM's data in a run without OMIM": (
        {'header': _header('ClinVar', 'OMIM')},
        "of 'omim'",
    ),
    "ClinVar's header entry in a run without ClinVar": (
        {
            'sources': ('Ensembl', 'RefSeq', 'gnomad', 'dbSNP'),
            'positions': [_without_clinvar(p) for p in _positions()],
            'header': _header('ClinVar'),
        },
        "of 'clinvar'",
    ),
}


@pytest.mark.parametrize('case', REFUSED)
def test_a_json_that_is_not_this_unions_annotation_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    changes, message = REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _convert(tmp_path, **changes)  # type: ignore[arg-type]


def test_a_run_config_of_another_catalog_is_refused() -> None:
    config = _run_config('Ensembl')
    config['data']['catalog'] = '209812'
    with pytest.raises(ValueError, match='the run config is of catalog 209812'):
        annotation_tables.annotation_info(_header(), config, _BUNDLE, _CATALOG)


def test_a_bundle_of_another_catalog_is_refused() -> None:
    bundle = copy.deepcopy(_BUNDLE)
    bundle['data']['catalog'] = '209812'
    with pytest.raises(ValueError, match="the bundle's run config is of catalog 209812"):
        annotation_tables.annotation_info(_header(), _run_config('Ensembl'), bundle, _CATALOG)


def test_a_run_config_naming_a_source_the_bundle_lacks_is_refused() -> None:
    config = _run_config('Ensembl')
    config['data']['versionSelector']['cosmic'] = {'SmallVariant': '1'}
    with pytest.raises(ValueError, match='bundle lacks'):
        annotation_tables.annotation_info(_header(), config, _BUNDLE, _CATALOG)


def test_a_run_config_pinning_another_version_than_the_bundle_is_refused() -> None:
    config = _run_config('Ensembl')
    config['data']['versionSelector']['Ensembl'] = {'SmallVariant': '2'}
    with pytest.raises(ValueError, match='the bundle'):
        annotation_tables.annotation_info(_header(), config, _BUNDLE, _CATALOG)


def test_a_union_record_of_two_alleles_is_refused(tmp_path: pathlib.Path) -> None:
    header = pysam.VariantHeader()
    header.contigs.add('chr1', length=1000)
    with pysam.VariantFile(str(tmp_path / 'u.vcf.gz'), 'wz', header=header) as vcf:
        vcf.write(vcf.new_record(contig='chr1', start=4, alleles=('C', 'T', 'G')))
    with pytest.raises(ValueError, match='one alternate allele'):
        annotation_tables.read_union_sites(tmp_path / 'u.vcf.gz')


def test_a_first_line_that_is_not_an_annotator_header_is_refused(tmp_path: pathlib.Path) -> None:
    with gzip.open(tmp_path / 'union.json.gz', 'wt') as f:
        f.write('{"positions":[\n')
    _write_sites(tmp_path / 'union.vcf.gz', _ALLELES)
    with pytest.raises(ValueError, match='not the header'):
        annotation_tables.convert(
            tmp_path / 'union.json.gz',
            tmp_path / 'union.vcf.gz',
            _run_config(*_ALL),
            _BUNDLE,
            _CATALOG,
            tmp_path / 'tables',
            workers=1,
        )


def test_a_document_without_its_closing_line_is_refused(tmp_path: pathlib.Path) -> None:
    _write_json(tmp_path / 'union.json.gz', _header('ClinVar'), _positions(), whole=False)
    with gzip.open(tmp_path / 'union.json.gz', 'at') as f:
        f.write('],"genes":[\n{"name":"GENEA"}\n')
    _write_sites(tmp_path / 'union.vcf.gz', _ALLELES)
    with pytest.raises(ValueError, match='does not end'):
        annotation_tables.convert(
            tmp_path / 'union.json.gz',
            tmp_path / 'union.vcf.gz',
            _run_config(*_ALL),
            _BUNDLE,
            _CATALOG,
            tmp_path / 'tables',
            workers=1,
        )


def test_a_union_with_samples_is_refused(tmp_path: pathlib.Path) -> None:
    header = pysam.VariantHeader()
    header.contigs.add('chr1', length=1000)
    header.add_meta('FORMAT', items=[('ID', 'GT'), ('Number', '1'), ('Type', 'String'), ('Description', 'x')])
    header.add_sample('S1')
    with pysam.VariantFile(str(tmp_path / 'u.vcf.gz'), 'wz', header=header) as vcf:
        vcf.write(vcf.new_record(contig='chr1', start=4, alleles=('C', 'T')))
    with pytest.raises(ValueError, match='has samples'):
        annotation_tables.read_union_sites(tmp_path / 'u.vcf.gz')
