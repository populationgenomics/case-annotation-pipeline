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
        gnomad={'allAf': 0.001, 'allAc': 2, 'allAn': 2000, 'allHc': 0, 'failedFilter': True},
        transcripts=[
            _transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant'),
            _transcript('ENST02.3', 'Ensembl', 'GENEA', 'ENSG01', 'missense_variant', isCanonical=True),
            _transcript(
                'NM_01.2', 'RefSeq', 'GENEA', '1001', 'missense_variant', isManeSelect=True, hgvsp='NP_01.2:p.(K1E)'
            ),
            _transcript('ENST03.1', 'Ensembl', 'GENEB', 'ENSG02', 'upstream_gene_variant'),
        ],
        alphaMissense=[
            {'transcriptId': 'ENST02.2', 'pathogenicity': 0.9, 'classification': 'likely_pathogenic'},
            {'transcriptId': 'ENST01.1', 'pathogenicity': 0.2, 'classification': 'likely_benign'},
        ],
        spliceAI=[{'geneId': 'ENSG01', 'donorGainScore0': 0.1, 'acceptorLossScore0': 0.4}],
        **{'primateAI-3D': [{'score': 0.7, 'ensemblTranscriptId': 'ENST02.3', 'refSeqTranscriptId': 'NM_01.1'}]},
        revel={'score': 0.5},
        dannScore=0.9,
        phylopScore=3.1,
        gerpScore=4.2,
        clinvar=[
            {
                'accession': 'VCV000000002',
                'isAlleleSpecific': True,
                'classifications': {
                    'germlineClassification': {
                        'classification': 'Uncertain significance',
                        'reviewStatus': 'criteria provided, single submitter',
                        'dateLastEvaluated': '2099-01-01',
                        'conditions': [{'traits': [{'name': {'value': 'invented condition'}}]}],
                    }
                },
            },
            {'accession': 'VCV000000003', 'refAllele': 'A', 'altAllele': 'T'},
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
    # Only the allele-specific ClinVar record is this allele's.
    assert row['clinvar_accession'] == 'VCV000000002'
    assert (row['clinvar_stars'], row['clinvar_conditions']) == (1, ['invented condition'])
    by_id = {t['transcript_id']: t for t in transcripts}
    assert by_id['ENST02.3']['alphamissense_pathogenicity'] == 0.9  # matched without the version
    assert by_id['NM_01.2']['primateai3d_score'] == 0.7
    assert by_id['NM_01.2']['spliceai_max'] == 0.4  # through the gene symbol
    assert by_id['ENST03.1']['spliceai_max'] is None
    assert {'dbsnp', 'clinvar', 'transcripts'} <= keys


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


def _clinvar_record(accession: str, status: str | None, evaluated: str, classification: str) -> dict[str, Any]:
    germline = {
        'classification': classification,
        'reviewStatus': status,
        'dateLastEvaluated': evaluated,
        'conditions': [{'traits': [{'name': {'value': f'condition of {accession}'}}]}],
    }
    return {
        'accession': accession,
        'version': 1,
        'recordType': 'classified',
        'isAlleleSpecific': True,
        'classifications': {'germlineClassification': germline} if status else {},
    }


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
        spliceAI=[{'geneId': 'ENSG01', 'donorGainScore0': 0.5}, {'geneId': 'ENSG01', 'acceptorGainScore0': 0.2}],
    )
    _, transcripts, _ = _position_rows(position)
    by_id = {t['transcript_id']: t for t in transcripts}
    assert by_id['ENST01.1']['alphamissense_pathogenicity'] == 0.8
    assert by_id['ENST01.1']['spliceai_max'] == 0.5
    assert by_id['NM_01.1']['spliceai_max'] == 0.5


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
        spliceAI=[{'geneId': 'ENSG01', 'donorGainScore0': 0.1}, {'geneId': 'ENSG02', 'donorGainScore0': 0.6}],
    )
    row, transcripts, _ = _position_rows(position)
    gene = row['genes'][0]
    assert (gene['worst_transcript_id'], gene['gene_symbol'], gene['ensembl_gene_id']) == ('NM_01.1', 'GENEP', None)
    assert next(t for t in transcripts if t['source'] == 'RefSeq')['spliceai_max'] == 0.6


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
