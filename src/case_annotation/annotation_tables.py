"""Illumina's annotation of a sites union as Parquet tables (docs/annotation.md, the family tables).

The DRAGEN Annotator writes one JSON document: a header line, one line per VCF record (a position), then a section of
genes and one of samples. The union is split into one allele per record, so each position line holds one allele. This
module reads the document line by line, without its `.jsi` index (whose offsets are wrong after a non-ASCII character),
and writes two tables of the union: `variants`, one row per allele, and `transcripts`, one row per allele and
transcript. The genes section holds a record per gene, of which the tables take the gene's HGNC id alone: its
gene-level knowledge, such as OMIM and ClinGen, stays out of them. The section follows the positions, and the rows of
both tables need its ids, so the document is read twice: once for the section, then for the positions.

The workers each convert a chunk of position lines into a Parquet file of each table, so the tables are folders of
parts in no particular order; the family tables sort them. A part holds its lines in the union's order, so its
positions rise, and each table stores `pos` as the differences between neighbours instead of with a dictionary.
"""

from __future__ import annotations

import concurrent.futures
import gzip
import json
import pathlib
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Any, Literal

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam

# Sequence Ontology consequences from most to least severe: Ensembl's order
# (https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html), with the terms Illumina's Annotator
# also writes placed among them. A term outside this list fails the conversion rather than being ranked by a guess.
# The Annotator's own terms and their places: inframe_indel, which it gives an in-frame replacement of bases, beside
# the other in-frame terms; the gene fusion of a structural variant after Ensembl's terms for a whole transcript, as
# Illumina rates it high; exon_variant beside the exon term of a non-coding transcript; and the terms that say only
# that a transcript is partly duplicated, that a copy number changes or that the variant is structural, which Illumina
# rates as modifiers or not at all, at the end before the term for any variant.
CONSEQUENCES = (
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'feature_elongation',
    'feature_truncation',
    'unidirectional_gene_fusion',
    'inframe_insertion',
    'inframe_deletion',
    'inframe_indel',
    'missense_variant',
    'protein_altering_variant',
    'splice_donor_5th_base_variant',
    'splice_region_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'coding_transcript_variant',
    'transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'regulatory_region_variant',
    'intergenic_variant',
    'five_prime_duplicated_transcript',
    'three_prime_duplicated_transcript',
    'copy_number_increase',
    'copy_number_decrease',
    'structural_variant',
    'sequence_variant',
)
_RANK = {c: i for i, c in enumerate(CONSEQUENCES)}

# ClinVar's review statuses and the gold stars ClinVar shows for each. A status outside the list keeps its text and has
# no stars.
_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, multiple submitters': 2,
    'criteria provided, conflicting classifications': 1,
    'criteria provided, conflicting interpretations': 1,
    'criteria provided, single submitter': 1,
    'no assertion criteria provided': 0,
    'no assertion provided': 0,
    'no classification provided': 0,
    'no classification for the individual variant': 0,
    'no classifications from unflagged records': 0,
}

# The four SpliceAI delta scores by the name each has in a column, each with the Annotator's keys of the score and of
# its distance. Illumina's filtered release omits a score below its threshold, and the score's distance with it.
_SPLICEAI: Mapping[str, tuple[str, str]] = {
    'acceptor_gain': ('acceptorGainScore0', 'acceptorGainDistance0'),
    'acceptor_loss': ('acceptorLossScore0', 'acceptorLossDistance0'),
    'donor_gain': ('donorGainScore0', 'donorGainDistance0'),
    'donor_loss': ('donorLossScore0', 'donorLossDistance0'),
}
# Every key of a SpliceAI entry. Illumina's SpliceAI file holds a second set of scores, which the Annotator leaves
# out. An entry with a key outside these is refused, so that a second set is never read as the first.
_SPLICEAI_KEYS = frozenset({'geneId'}) | {key for keys in _SPLICEAI.values() for key in keys}

# gnomAD's ancestry groups and then its two sexes, under the names the Annotator gives them, by the key that each of
# gnomAD's two sets has on a variant. The exomes have no `ami` group.
_GNOMAD_GROUPS: Mapping[str, tuple[str, ...]] = {
    'gnomad': ('afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'sas', 'remaining', 'female', 'male'),
    'gnomad-exome': ('afr', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'sas', 'remaining', 'female', 'male'),
}
# The four values the Annotator gives of a group, by the name each has in a column, with the ending of its key:
# `afrAf` is the frequency in `afr`.
_GNOMAD_VALUES: Mapping[str, str] = {'af': 'Af', 'ac': 'Ac', 'an': 'An', 'hom': 'Hc'}
# Every key of a gnomAD entry, by the entry's key on a variant. An entry with a key outside these is refused, so that a
# group the columns have no field for is never left out of them without a word.
_GNOMAD_KEYS: Mapping[str, frozenset[str]] = {
    key: frozenset({'coverage', 'failedFilter'})
    | {f'{group}{ending}' for group in ('all', *groups) for ending in _GNOMAD_VALUES.values()}
    for key, groups in _GNOMAD_GROUPS.items()
}


def gnomad_fields(prefix: str) -> list[pa.Field]:
    # Illumina gives a frequency to six decimal places, which a 32-bit float holds.
    return [
        pa.field(f'{prefix}_af', pa.float32()),
        pa.field(f'{prefix}_ac', pa.int32()),
        pa.field(f'{prefix}_an', pa.int32()),
        pa.field(f'{prefix}_hom', pa.int32()),
        pa.field(f'{prefix}_failed_filter', pa.bool_()),
    ]


# What gnomAD counted in one ancestry group or one sex: the allele's frequency, its count, the number of alleles
# called and the Annotator's count of homozygotes.
GNOMAD_GROUP = pa.struct(
    [
        pa.field('af', pa.float32()),
        pa.field('ac', pa.int32()),
        pa.field('an', pa.int32()),
        pa.field('hom', pa.int32()),
    ]
)


def _gnomad_detail_fields(prefix: str, key: str) -> list[pa.Field]:
    """The columns of a gnomAD entry beyond its all-sample values: its coverage, and a struct with a field per group."""
    return [
        pa.field(f'{prefix}_coverage', pa.int32()),
        pa.field(f'{prefix}_groups', pa.struct([pa.field(group, GNOMAD_GROUP) for group in _GNOMAD_GROUPS[key]])),
    ]


_STRING = pa.string()
_FLOAT = pa.float64()
# An HGNC id as a number: HGNC:1100 is 1100.
_HGNC_ID = pa.int32()

# One ClinVar record (a VCV) of the allele. An allele can have several: a record of the allele alone and one of a
# haplotype or compound that holds it, or a duplicate submission.
CLINVAR_RECORD = pa.struct(
    [
        pa.field('accession', _STRING),
        pa.field('version', _STRING),
        pa.field('record_type', _STRING),
        pa.field('variant_type', _STRING),
        pa.field('classification', _STRING),
        pa.field('review_status', _STRING),
        pa.field('stars', pa.int8()),
        pa.field('last_evaluated', _STRING),
        pa.field('conditions', pa.list_(_STRING)),
    ]
)

# One ClinVar record of another allele at the allele's position: the fields of a record of the allele itself, and after
# them the two alleles that the record is of, as the Annotator writes them. The Annotator trims them, so they are not
# spelt as the `ref` and `alt` of `variants` are.
CLINVAR_OTHER_ALLELE = pa.struct(
    [
        *CLINVAR_RECORD,
        pa.field('ref_allele', _STRING, nullable=False),
        pa.field('alt_allele', _STRING, nullable=False),
    ]
)

# One gene of the allele (_gene_entry): the gene, then what the transcript chosen among the gene's own says of it,
# then the gene's MANE Select transcript. Whatever is of one gene or one transcript of an allele stands here and in
# `transcripts`: a column of `variants` beside `genes` describes the allele as a whole. The gene's id is named
# `ensembl_gene_id` because it is Ensembl's whatever the chosen transcript's source, while `gene_id` in `transcripts`
# is the id that the transcript's own source gives the gene.
GENE = pa.struct(
    [
        pa.field('gene_symbol', _STRING),
        pa.field('ensembl_gene_id', _STRING),
        pa.field('hgnc_id', _HGNC_ID),
        pa.field('worst_consequence', _STRING, nullable=False),
        pa.field('worst_impact', _STRING),
        pa.field('worst_transcript_id', _STRING, nullable=False),
        pa.field('mane_select_transcript_id', _STRING),
        pa.field('mane_select_consequence', _STRING),
        pa.field('mane_select_hgvsc', _STRING),
        pa.field('mane_select_hgvsp', _STRING),
    ]
)

# One regulatory feature that the allele overlaps: its id and its type, both as the Annotator gives them.
REGULATORY_FEATURE = pa.struct(
    [
        pa.field('id', _STRING, nullable=False),
        pa.field('type', _STRING, nullable=False),
    ]
)

VARIANTS = pa.schema(
    [
        pa.field('chrom', _STRING, nullable=False),
        pa.field('pos', pa.int32(), nullable=False),
        pa.field('ref', _STRING, nullable=False),
        pa.field('alt', _STRING, nullable=False),
        pa.field('variant_id', _STRING, nullable=False),
        pa.field('variant_type', _STRING),
        pa.field('hgvsg', _STRING),
        pa.field('rsids', pa.list_(_STRING)),
        pa.field('worst_consequence', _STRING),
        pa.field('worst_impact', _STRING),
        pa.field('hgnc_ids', pa.list_(_HGNC_ID)),
        pa.field('genes', pa.list_(GENE)),
        pa.field('regulatory_features', pa.list_(REGULATORY_FEATURE)),
        pa.field('in_low_complexity_region', pa.bool_()),
        *gnomad_fields('gnomad_genomes'),
        *_gnomad_detail_fields('gnomad_genomes', 'gnomad'),
        *gnomad_fields('gnomad_exomes'),
        *_gnomad_detail_fields('gnomad_exomes', 'gnomad-exome'),
        pa.field('alphamissense_max', _FLOAT),
        pa.field('alphamissense_max_class', _STRING),
        pa.field('revel', _FLOAT),
        pa.field('spliceai_max', _FLOAT),
        pa.field('primateai3d_max', _FLOAT),
        pa.field('primateai3d_max_percentile', _FLOAT),
        pa.field('primateai3d_max_class', _STRING),
        pa.field('promoterai_min', _FLOAT),
        pa.field('promoterai_max', _FLOAT),
        pa.field('dann', _FLOAT),
        pa.field('phylop', _FLOAT),
        pa.field('phylop_primate', _FLOAT),
        pa.field('gerp', _FLOAT),
        pa.field('clinvar_accession', _STRING),
        pa.field('clinvar_classification', _STRING),
        pa.field('clinvar_review_status', _STRING),
        pa.field('clinvar_stars', pa.int8()),
        pa.field('clinvar_conditions', pa.list_(_STRING)),
        pa.field('clinvar_last_evaluated', _STRING),
        pa.field('clinvar_last_updated', _STRING),
        pa.field('clinvar_pubmed_ids', pa.list_(_STRING)),
        pa.field('clinvar_records', pa.list_(CLINVAR_RECORD)),
        pa.field('clinvar_other_alleles', pa.list_(CLINVAR_OTHER_ALLELE)),
    ]
)

TRANSCRIPTS = pa.schema(
    [
        pa.field('chrom', _STRING, nullable=False),
        pa.field('pos', pa.int32(), nullable=False),
        pa.field('variant_id', _STRING, nullable=False),
        pa.field('transcript_id', _STRING, nullable=False),
        pa.field('source', _STRING, nullable=False),
        pa.field('biotype', _STRING),
        pa.field('gene_symbol', _STRING),
        pa.field('gene_id', _STRING),
        pa.field('hgnc_id', _HGNC_ID),
        pa.field('consequences', pa.list_(_STRING), nullable=False),
        pa.field('worst_consequence', _STRING, nullable=False),
        pa.field('impact', _STRING),
        pa.field('hgvsc', _STRING),
        pa.field('hgvsp', _STRING),
        pa.field('protein_id', _STRING),
        pa.field('is_mane_select', pa.bool_(), nullable=False),
        pa.field('is_canonical', pa.bool_(), nullable=False),
        pa.field('exon', _STRING),
        pa.field('intron', _STRING),
        pa.field('cdna_pos', _STRING),
        pa.field('cds_pos', _STRING),
        pa.field('protein_pos', _STRING),
        pa.field('codons', _STRING),
        pa.field('amino_acids', _STRING),
        pa.field('amino_acid_conservation', pa.list_(_FLOAT)),
        pa.field('alphamissense_pathogenicity', _FLOAT),
        pa.field('alphamissense_class', _STRING),
        pa.field('alphamissense_is_isoform', pa.bool_()),
        pa.field('primateai3d_score', _FLOAT),
        pa.field('primateai3d_percentile', _FLOAT),
        pa.field('primateai3d_class', _STRING),
        pa.field('spliceai_max', _FLOAT),
        *(
            field
            for kind in _SPLICEAI
            for field in (
                pa.field(f'spliceai_{kind}', _FLOAT),
                pa.field(f'spliceai_{kind}_distance', pa.int32()),
            )
        ),
        pa.field('promoterai_score', _FLOAT),
        pa.field('promoterai_tss_distance', pa.int32()),
    ]
)

# The columns each of the union's tables writes without a dictionary, each by its path in the file, and the encoding
# of its values. pos rises along a part and rarely repeats, so a dictionary would hold nearly every position once, and
# DELTA_BINARY_PACKED stores the small differences between neighbours.
_ENCODINGS: Mapping[str, Mapping[str, str]] = {
    'variants': {'pos': 'DELTA_BINARY_PACKED'},
    'transcripts': {'pos': 'DELTA_BINARY_PACKED'},
}


def column_options(schema: pa.Schema, encodings: Mapping[str, str]) -> dict[str, Any]:
    """The Parquet writer's options that give these columns their encodings, and every other column a dictionary.

    pyarrow takes the columns that get a dictionary, each by its path in the file, and an encoding only for a column
    that gets none. An empty file of the schema gives the paths as the writer names them, a nested column's included.
    A column that the schema lacks fails, since a misspelt name would otherwise leave its column's dictionary on.
    """
    empty = pa.BufferOutputStream()
    pq.ParquetWriter(empty, schema).close()
    written = pq.ParquetFile(pa.BufferReader(empty.getvalue())).schema
    paths = [written.column(i).path for i in range(len(written.names))]
    if unknown := sorted(set(encodings) - set(paths)):
        raise ValueError(f'no columns {unknown} to write without a dictionary, among {paths}')
    return {
        'use_dictionary': [path for path in paths if path not in encodings],
        'column_encoding': dict(encodings),
    }


# The columns of the variants table that hold a dated source's data, by the source's name in a run config. A family's
# variants table leaves out the columns of a source its run did not name, so that a missing source reads as missing
# rather than as an empty result.
SOURCE_COLUMNS: Mapping[str, tuple[str, ...]] = {
    'dbSNP': ('rsids',),
    'clinvar': tuple(f.name for f in VARIANTS if f.name.startswith('clinvar_')),
}
# The key each such source has on a variant of the JSON: a document holding it although the run did not name the
# source is refused.
_SOURCE_KEYS = {'dbSNP': 'dbsnp', 'clinvar': 'clinvar'}

# The dated sources of docs/annotation.md, and the names their data has among the JSON header's data sources. dbSNP
# shares its header name with globalAllele, a current source, and is checked in the records instead.
DATED = frozenset({'clinvar', 'clinvar-preview', 'omim', 'clingen', 'clingen (legacy)', 'mitomap', 'dbSNP'})
_DATED_HEADER_NAMES = {
    'clinvar': ('ClinVar',),
    'clinvar-preview': ('ClinVarPreview',),
    'omim': ('OMIM',),
    'clingen': ('ClinGen disease validity curations', 'ClinGen Dosage Sensitivity Map'),
    'clingen (legacy)': ('ClinGen',),
    'mitomap': ('MITOMAP', 'MITOMAP_SV'),
}
# The current sources of docs/annotation.md, and those that do not apply to germline small variants.
CURRENT = frozenset(
    {
        'Ensembl',
        'RefSeq',
        'HGNC',
        'GenomeAssembly',
        'gnomad',
        'gnomad-exome',
        'TOPMed',
        'oneKg',
        'ABRaOM',
        'GME',
        'globalAllele',
        'DECIPHER',
        'alphaMissense',
        'spliceAI',
        'primateAI-3D',
        'PrimateAI',
        'PromoterAI',
        'REVEL',
        'DANN',
        'Gerp',
        'phylopScore',
        'MultiZ100Way',
    }
)
NOT_APPLICABLE = frozenset({'FusionCatcher', 'UcscCpgIsland'})

type OmissionReason = Literal['dated', 'not_applicable', 'current', 'unplaced']


class DataSource(pydantic.BaseModel, frozen=True, extra='forbid'):
    """One data source as the JSON header records it."""

    name: str
    version: str
    release_date: str | None


class SourceVersion(pydantic.BaseModel, frozen=True, extra='forbid'):
    """One version the run config pins: a source, an annotation type and its version."""

    source: str
    annotation_type: str
    version: str


class OmittedSource(pydantic.BaseModel, frozen=True, extra='forbid'):
    """A source of the downloaded bundle that the run did not name, and why.

    The reason is `dated` for a dated source, which the deployment left out;
    `not_applicable` for one that does not apply to germline small variants; `current` for a current source the run did
    not name; and `unplaced` for a source docs/annotation.md does not place, which counts as dated.
    """

    source: str
    reason: OmissionReason


class AnnotationInfo(pydantic.BaseModel, frozen=True, extra='forbid'):
    """What annotated the union.

    The Annotator, the catalog, the sources the run named and the versions they pinned, the data sources the JSON
    header records, and the bundle's sources the run left out.
    """

    annotator: str
    json_schema_version: int
    assembly: str
    catalog: str
    sources: tuple[str, ...]
    source_versions: tuple[SourceVersion, ...]
    data_sources: tuple[DataSource, ...]
    omitted: tuple[OmittedSource, ...]


def variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f'{chrom}-{pos}-{ref}-{alt}'


def record_id(record: pysam.VariantRecord, where: str) -> str:
    """The variant id of a VCF record of one allele; a record with another number of alternate alleles fails."""
    if record.ref is None or record.alts is None or len(record.alts) != 1:
        raise ValueError(f'{where}: {record.chrom}:{record.pos} does not hold one reference and one alternate allele')
    return variant_id(record.chrom, record.pos, record.ref, record.alts[0])


def by_severity(consequences: Iterable[str]) -> list[str]:
    """The consequences of a transcript from most to least severe; an unknown term fails, and so does none at all."""
    ranked = []
    for c in consequences:
        if c not in _RANK:
            raise ValueError(f'consequence {c!r} is not in the severity order')
        ranked.append((_RANK[c], c))
    if not ranked:
        raise ValueError('a transcript without consequences')
    return [c for _, c in sorted(ranked)]


def worst(consequences: Iterable[str]) -> str:
    """The most severe of the consequences; an unknown term fails."""
    return by_severity(consequences)[0]


def _unversioned(identifier: str) -> str:
    return identifier.split('.', 1)[0]


def _max_by(entries: Iterable[Mapping[str, Any]], key: str) -> Mapping[str, Any] | None:
    scored = [e for e in entries if e.get(key) is not None]
    return max(scored, key=lambda e: e[key]) if scored else None


def _require(entry: Mapping[str, Any], keys: Sequence[str], what: str, where: str) -> None:
    if missing := [key for key in keys if entry.get(key) is None]:
        raise ValueError(f'{where}: a {what} entry without {missing}')


def _spliceai_values(entry: Mapping[str, Any], where: str) -> dict[str, tuple[float, int]]:
    """The scores of a SpliceAI entry by their kind (_SPLICEAI), each with its distance.

    Fails on a key outside _SPLICEAI_KEYS, and on a score without its distance or a distance without its score.
    """
    if unknown := sorted(set(entry) - _SPLICEAI_KEYS):
        raise ValueError(f'{where}: a SpliceAI entry with the keys {unknown}, which the conversion does not know')
    values: dict[str, tuple[float, int]] = {}
    for kind, (score_key, distance_key) in _SPLICEAI.items():
        score, distance = entry.get(score_key), entry.get(distance_key)
        if (score is None) != (distance is None):
            raise ValueError(f'{where}: a SpliceAI entry with one of {score_key} and {distance_key} and not the other')
        if score is not None and distance is not None:
            values[kind] = (score, distance)
    return values


def _keep_highest(held: dict[str, tuple[float, int]], values: Mapping[str, tuple[float, int]]) -> None:
    """Raises each kind's held SpliceAI score and distance to the given ones where those score higher.

    Of two equal scores the one with the greater distance stays, so that the order of the entries decides nothing.
    """
    for kind, pair in values.items():
        held[kind] = max(held.get(kind, pair), pair)


def _spliceai_columns(values: Mapping[str, tuple[float, int]]) -> dict[str, float | int | None]:
    """The SpliceAI columns of a transcript row: the highest of the scores, then each score with its distance."""
    columns: dict[str, float | int | None] = {
        'spliceai_max': max((score for score, _ in values.values()), default=None)
    }
    for kind in _SPLICEAI:
        score, distance = values.get(kind, (None, None))
        columns |= {f'spliceai_{kind}': score, f'spliceai_{kind}_distance': distance}
    return columns


def _regulatory_features(entries: Iterable[Mapping[str, Any]], where: str) -> list[dict[str, str]] | None:
    """The regulatory features of an allele (REGULATORY_FEATURE) by id, and no list for an allele without one.

    Fails on an entry without `id` or `type`.
    """
    features = []
    for entry in entries:
        _require(entry, ('id', 'type'), 'regulatory', where)
        features.append({'id': entry['id'], 'type': entry['type']})
    return sorted(features, key=lambda feature: (feature['id'], feature['type'])) or None


def gnomad_values(prefix: str, data: Mapping[str, Any] | None) -> dict[str, Any]:
    if data is None:
        return {}
    return {
        f'{prefix}_af': data.get('allAf'),
        f'{prefix}_ac': data.get('allAc'),
        f'{prefix}_an': data.get('allAn'),
        f'{prefix}_hom': data.get('allHc'),
        f'{prefix}_failed_filter': bool(data.get('failedFilter', False)),
    }


def _gnomad_details(prefix: str, key: str, data: Mapping[str, Any] | None, where: str) -> dict[str, Any]:
    """The coverage and the groups of a gnomAD entry (_gnomad_detail_fields), and nothing for an allele without one.

    A group that the entry leaves out has no struct. Fails on an entry without `coverage`, on a key outside
    _GNOMAD_KEYS, and on a group with some of its four values and not the others.
    """
    if data is None:
        return {}
    _require(data, ('coverage',), key, where)
    if unknown := sorted(set(data) - _GNOMAD_KEYS[key]):
        raise ValueError(f'{where}: a {key} entry with the keys {unknown}, which the conversion does not know')
    groups: dict[str, dict[str, Any] | None] = {}
    for group in _GNOMAD_GROUPS[key]:
        values = {name: data.get(f'{group}{ending}') for name, ending in _GNOMAD_VALUES.items()}
        given = [value is not None for value in values.values()]
        if any(given) and not all(given):
            raise ValueError(f'{where}: a {key} entry with some of the four values of {group!r} and not the others')
        groups[group] = values if all(given) else None
    return {f'{prefix}_coverage': data['coverage'], f'{prefix}_groups': groups}


def _conservation(transcript: Mapping[str, Any], where: str) -> list[float] | None:
    """The amino-acid conservation scores of a transcript as the Annotator lists them, and None without any.

    Fails on an `aminoAcidConservation` entry without `scores`.
    """
    entry = transcript.get('aminoAcidConservation')
    if entry is None:
        return None
    _require(entry, ('scores',), 'transcript aminoAcidConservation', where)
    return list(entry['scores']) or None


def clinvar_record(record: Mapping[str, Any]) -> dict[str, Any]:
    germline = record.get('classifications', {}).get('germlineClassification', {})
    status = germline.get('reviewStatus')
    return {
        'accession': record.get('accession'),
        'version': None if record.get('version') is None else str(record['version']),
        'record_type': record.get('recordType'),
        'variant_type': record.get('variantType'),
        'classification': germline.get('classification'),
        'review_status': status,
        'stars': _STARS.get(status) if status is not None else None,
        'last_evaluated': germline.get('dateLastEvaluated'),
        'conditions': [
            trait['name']['value']
            for condition in germline.get('conditions', ())
            for trait in condition.get('traits', ())
            if trait.get('name', {}).get('value')
        ]
        or None,
    }


def by_review(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """ClinVar records in the order that puts the best reviewed first.

    By their gold stars (most first, a record without a germline classification last), then by the date their
    classification was last evaluated (latest first), then by accession.
    """
    ordered = sorted(records, key=lambda r: r['accession'] or '')
    ordered.sort(key=lambda r: r['last_evaluated'] or '', reverse=True)
    ordered.sort(key=lambda r: -1 if r['stars'] is None else r['stars'], reverse=True)
    return ordered


def clinvar_columns(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The ClinVar columns of ordered records: all of them, and a summary from the first; nothing without a record."""
    if not records:
        return {}
    top = records[0]
    return {
        'clinvar_accession': top['accession'],
        'clinvar_classification': top['classification'],
        'clinvar_review_status': top['review_status'],
        'clinvar_stars': top['stars'],
        'clinvar_conditions': top['conditions'],
        'clinvar_last_evaluated': top['last_evaluated'],
        'clinvar_records': list(records),
    }


def allele_clinvar(entries: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The allele's ClinVar records, all of them, and a summary from the first.

    Entries for other alleles at the position are not allele-specific and are left out. The records are ordered by
    their review (by_review), and the summary columns are the first record's.
    """
    return clinvar_columns(by_review([clinvar_record(e) for e in entries if e.get('isAlleleSpecific')]))


def _pubmed_ids(entry: Mapping[str, Any]) -> list[str] | None:
    """The PubMed ids that a ClinVar record cites, each once and sorted, and None for a record that cites none.

    They are those of the record's germline classification and those of its submissions (`clinicalAssertions`).
    """
    germline = entry.get('classifications', {}).get('germlineClassification', {})
    cited = set(germline.get('pubMedIds', ()))
    for assertion in entry.get('clinicalAssertions', ()):
        cited.update(assertion.get('pubMedIds', ()))
    return sorted(cited) or None


def _other_allele(entry: Mapping[str, Any], where: str) -> dict[str, Any]:
    """A ClinVar record of another allele (CLINVAR_OTHER_ALLELE); fails without its accession or one of its alleles."""
    _require(entry, ('accession', 'refAllele', 'altAllele'), 'ClinVar', where)
    return clinvar_record(entry) | {'ref_allele': entry['refAllele'], 'alt_allele': entry['altAllele']}


def variant_clinvar(entries: Sequence[Mapping[str, Any]], where: str) -> dict[str, Any]:
    """The ClinVar columns of an allele of `variants`; nothing for an allele without an entry.

    The allele's own records give what allele_clinvar gives. The first of them, which the summary columns are of, also
    gives the date ClinVar last updated it and the PubMed ids it cites (_pubmed_ids). The entries that are not the
    allele's own are the records of other alleles at its position, ordered as the allele's own are (by_review). Fails
    on a record of the allele without `dateLastUpdated`, and as _other_allele does.
    """
    columns = allele_clinvar(entries)
    own = [entry for entry in entries if entry.get('isAlleleSpecific')]
    for entry in own:
        _require(entry, ('dateLastUpdated',), 'ClinVar', where)
    if own:
        # Each record with what the two columns hold of it, in the order that allele_clinvar gave the records alone.
        first = by_review(
            [
                clinvar_record(entry) | {'last_updated': entry['dateLastUpdated'], 'pubmed_ids': _pubmed_ids(entry)}
                for entry in own
            ]
        )[0]
        columns |= {'clinvar_last_updated': first['last_updated'], 'clinvar_pubmed_ids': first['pubmed_ids']}
    others = by_review([_other_allele(entry, where) for entry in entries if not entry.get('isAlleleSpecific')])
    if others:
        columns['clinvar_other_alleles'] = others
    return columns


def _highest(entries: Iterable[Mapping[str, Any]], score: str, *id_keys: str) -> dict[str, Mapping[str, Any]]:
    """A source's entries by unversioned transcript or gene id; where one id has several, the highest scoring."""
    best: dict[str, Mapping[str, Any]] = {}
    for e in entries:
        for k in id_keys:
            if not e.get(k):
                continue
            key = _unversioned(e[k])
            held = best.get(key)
            if held is None or (e.get(score) is not None and (held.get(score) is None or e[score] > held[score])):
                best[key] = e
    return best


def hgnc_ids(genes: Iterable[Mapping[str, Any]]) -> dict[str, int | None]:
    """Each gene symbol of the genes section with its HGNC id, or None for a symbol without one.

    A record names its gene by `name`, the symbol that a transcript of the gene holds as `hgnc`, and gives the gene's
    HGNC id as `hgncId` where HGNC names the gene; many genes, such as those Ensembl has not named, have none. Several
    records can hold one symbol, such as a record of Ensembl's gene and one of NCBI's that the Annotator did not join,
    of which one may lack the id: the symbol has the id that any of them gives. Fails on a record without a symbol, on
    an id that is not a positive integer, and on a symbol whose records give two ids.
    """
    ids: dict[str, int | None] = {}
    for gene in genes:
        name, hgnc = gene.get('name'), gene.get('hgncId')
        if not isinstance(name, str) or not name:
            raise ValueError(f'a record of the genes section names no gene: {dict(gene)}')
        if hgnc is not None and (isinstance(hgnc, bool) or not isinstance(hgnc, int) or hgnc <= 0):
            raise ValueError(f'the genes section gives {name!r} the HGNC id {hgnc!r}, which is not a positive integer')
        held = ids.setdefault(name, hgnc)
        if held is None:
            ids[name] = hgnc
        elif hgnc is not None and hgnc != held:
            raise ValueError(f'the genes section gives {name!r} two HGNC ids, {held} and {hgnc}')
    return ids


def _preference(row: Mapping[str, Any]) -> tuple[int, bool, bool, bool, str]:
    """Where a transcript row stands among those of its allele, or of one gene of it: the first one summarises them.

    The most severe consequence comes first; among transcripts with the same, the MANE Select one, then a canonical
    one, then Ensembl's, then the first by id.
    """
    return (
        _RANK[row['worst_consequence']],
        not row['is_mane_select'],
        not row['is_canonical'],
        row['source'] != 'Ensembl',
        row['transcript_id'],
    )


def _gene_key(row: Mapping[str, Any]) -> tuple[str, int | str] | None:
    """What tells the gene of a transcript row from the other genes of its allele.

    The HGNC id, which the transcripts of Ensembl's and of RefSeq's model of a gene share, while each has the gene id
    of its own model. A gene without an HGNC id is told by its symbol, which the two models share as well, and a
    transcript without a symbol by the gene id of its own model. The transcripts that name no gene at all count as of
    one gene.
    """
    for name in ('hgnc_id', 'gene_symbol', 'gene_id'):
        if row[name] is not None:
            return name, row[name]
    return None


def _gene_entry(rows: Sequence[Mapping[str, Any]], ensembl_gene_id: str | None) -> dict[str, Any]:
    """One entry of `genes` (GENE), from the transcript rows of one gene of an allele in their order (_preference).

    The first row is the transcript chosen for the gene: it gives the entry its symbol, HGNC id, consequence and
    impact, and `ensembl_gene_id` is that row's Ensembl gene. The gene's MANE Select transcript gives the four fields
    named after it, and of several such transcripts the first by id does.
    """
    chosen = rows[0]
    mane = min((r for r in rows if r['is_mane_select']), key=lambda r: r['transcript_id'], default=None)
    of_mane: Mapping[str, Any] = {} if mane is None else mane
    return {
        'gene_symbol': chosen['gene_symbol'],
        'ensembl_gene_id': ensembl_gene_id,
        'hgnc_id': chosen['hgnc_id'],
        'worst_consequence': chosen['worst_consequence'],
        'worst_impact': chosen['impact'],
        'worst_transcript_id': chosen['transcript_id'],
        'mane_select_transcript_id': of_mane.get('transcript_id'),
        'mane_select_consequence': of_mane.get('worst_consequence'),
        'mane_select_hgvsc': of_mane.get('hgvsc'),
        'mane_select_hgvsp': of_mane.get('hgvsp'),
    }


def position_rows(
    position: Mapping[str, Any], hgnc: Mapping[str, int | None]
) -> tuple[dict[str, Any], list[dict[str, Any]], set[str]]:
    """One position line's row of `variants`, its rows of `transcripts`, and the keys of the sources on its allele.

    `hgnc` holds each gene symbol of the document's genes section with its HGNC id (hgnc_ids). The section lists the
    gene of every transcript, so a transcript whose symbol it lacks fails. So does a transcript listed twice: with the
    allele it is the key of a row of `transcripts`, which an entry of `genes` names its transcripts by.

    A SpliceAI score belongs to a gene, and a transcript row holds the highest score of each kind among the entries of
    its genes, with that score's distance. A PromoterAI entry belongs to a transcript, and the allele's
    `promoterai_min` and `promoterai_max` are over all its entries. The entries of both and of the regulatory features
    are held to the keys this function reads (_spliceai_values, _require), and so are the entries of gnomAD
    (_gnomad_details), of ClinVar (variant_clinvar) and of a transcript's amino-acid conservation (_conservation).

    The Annotator writes a flag only where it is true, so an allele without `inLowComplexityRegion` is outside those
    regions, and an AlphaMissense entry without `isIsoform` is not of an isoform.
    """
    chrom, pos, ref = position['chromosome'], position['position'], position['refAllele']
    alts = position['altAlleles']
    variants = position.get('variants', [])
    where = f'{chrom}:{pos}'
    if len(alts) != 1 or len(variants) > 1:
        raise ValueError(f'{where}: {len(alts)} alleles and {len(variants)} variants; the union holds one each')
    vid = variant_id(chrom, pos, ref, alts[0])
    row: dict[str, Any] = {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alts[0], 'variant_id': vid}
    if not variants:
        return row, [], set()
    v = variants[0]
    keys = set(v)

    alphamissense = _highest(v.get('alphaMissense', ()), 'pathogenicity', 'transcriptId')
    primateai3d = _highest(v.get('primateAI-3D', ()), 'score', 'ensemblTranscriptId', 'refSeqTranscriptId')
    for e in v.get('promoterAI', ()):
        _require(e, ('score', 'distanceFromTss'), 'PromoterAI', where)
    promoterai = _highest(v.get('promoterAI', ()), 'score', 'transcriptId')
    promoterai_scores = [e['score'] for e in v.get('promoterAI', ())]
    spliceai: dict[str, dict[str, tuple[float, int]]] = {}
    for e in v.get('spliceAI', ()):
        values = _spliceai_values(e, where)
        if e.get('geneId'):
            _keep_highest(spliceai.setdefault(e['geneId'], {}), values)
    transcripts = v.get('transcripts', [])
    # A RefSeq transcript's gene id is NCBI's; its symbol leads to the Ensembl genes of the same allele's transcripts.
    # A symbol can name more than one Ensembl gene (on X's and Y's pseudo-autosomal regions, for one).
    ensembl_genes: dict[str, set[str]] = {}
    for t in transcripts:
        if t.get('source') == 'Ensembl' and t.get('hgnc') and t.get('geneId'):
            ensembl_genes.setdefault(t['hgnc'], set()).add(t['geneId'])

    def genes_of(t: Mapping[str, Any]) -> set[str]:
        if t.get('source') == 'Ensembl':
            return {t['geneId']} if t.get('geneId') else set()
        return ensembl_genes.get(t.get('hgnc') or '', set())

    def ensembl_gene(t: Mapping[str, Any]) -> str | None:
        """The transcript's Ensembl gene; None if its symbol names none, or more than one."""
        genes = genes_of(t)
        return next(iter(genes)) if len(genes) == 1 else None

    def gene_spliceai(t: Mapping[str, Any]) -> dict[str, tuple[float, int]]:
        best: dict[str, tuple[float, int]] = {}
        for g in genes_of(t):
            _keep_highest(best, spliceai.get(g, {}))
        return best

    transcript_rows = []
    listed: set[str] = set()
    for t in transcripts:
        tid = t['transcript']
        if tid in listed:
            raise ValueError(f'{where}: the transcript {tid} twice')
        listed.add(tid)
        symbol = t.get('hgnc')
        if symbol is not None and symbol not in hgnc:
            raise ValueError(f'{where}: the genes section lacks {symbol!r}, the gene of the transcript {tid}')
        consequences = by_severity(t['consequence'])
        am = alphamissense.get(_unversioned(tid))
        p3d = primateai3d.get(_unversioned(tid))
        promoter = promoterai.get(_unversioned(tid))
        transcript_rows.append(
            {
                'chrom': chrom,
                'pos': pos,
                'variant_id': vid,
                'transcript_id': tid,
                'source': t['source'],
                'biotype': t.get('bioType'),
                'gene_symbol': symbol,
                'gene_id': t.get('geneId'),
                'hgnc_id': None if symbol is None else hgnc[symbol],
                'consequences': consequences,
                'worst_consequence': consequences[0],
                'impact': t.get('impact'),
                'hgvsc': t.get('hgvsc'),
                'hgvsp': t.get('hgvsp'),
                'protein_id': t.get('proteinId'),
                'is_mane_select': bool(t.get('isManeSelect', False)),
                'is_canonical': bool(t.get('isCanonical', False)),
                'exon': t.get('exons'),
                'intron': t.get('introns'),
                'cdna_pos': t.get('cdnaPos'),
                'cds_pos': t.get('cdsPos'),
                'protein_pos': t.get('proteinPos'),
                'codons': t.get('codons'),
                'amino_acids': t.get('aminoAcids'),
                'amino_acid_conservation': _conservation(t, where),
                'alphamissense_pathogenicity': am.get('pathogenicity') if am else None,
                'alphamissense_class': am.get('classification') if am else None,
                'alphamissense_is_isoform': bool(am.get('isIsoform', False)) if am else None,
                'primateai3d_score': p3d.get('score') if p3d else None,
                'primateai3d_percentile': p3d.get('scorePercentile') if p3d else None,
                'primateai3d_class': p3d.get('classification') if p3d else None,
                **_spliceai_columns(gene_spliceai(t)),
                'promoterai_score': promoter['score'] if promoter else None,
                'promoterai_tss_distance': promoter['distanceFromTss'] if promoter else None,
            }
        )

    if transcript_rows:
        # Each row with the Annotator's transcript it was made from, in the rows' order of preference. The first row
        # gives the allele its consequence and impact, and the first row of each gene that gene's entry, so the genes
        # stand in the order of their first rows and the allele's consequence is that of its first gene.
        ranked = sorted(zip(transcript_rows, transcripts, strict=True), key=lambda pair: _preference(pair[0]))
        chosen = ranked[0][0]
        by_gene: dict[tuple[str, int | str] | None, list[tuple[dict[str, Any], Mapping[str, Any]]]] = {}
        for pair in ranked:
            by_gene.setdefault(_gene_key(pair[0]), []).append(pair)
        row |= {
            'worst_consequence': chosen['worst_consequence'],
            'worst_impact': chosen['impact'],
            # Every gene of the allele's transcripts that has an id, and no list for an allele without one.
            'hgnc_ids': sorted({r['hgnc_id'] for r in transcript_rows if r['hgnc_id'] is not None}) or None,
            'genes': [
                _gene_entry([r for r, _ in members], ensembl_gene(members[0][1])) for members in by_gene.values()
            ],
        }
    top_am = _max_by(v.get('alphaMissense', ()), 'pathogenicity')
    top_p3d = _max_by(v.get('primateAI-3D', ()), 'score')
    splice_scores = [score for values in spliceai.values() for score, _ in values.values()]
    row |= {
        'variant_type': v.get('variantType'),
        'hgvsg': v.get('hgvsg'),
        'rsids': v.get('dbsnp') or None,
        'regulatory_features': _regulatory_features(v.get('regulatoryRegions', ()), where),
        'in_low_complexity_region': bool(v.get('inLowComplexityRegion', False)),
        **gnomad_values('gnomad_genomes', v.get('gnomad')),
        **_gnomad_details('gnomad_genomes', 'gnomad', v.get('gnomad'), where),
        **gnomad_values('gnomad_exomes', v.get('gnomad-exome')),
        **_gnomad_details('gnomad_exomes', 'gnomad-exome', v.get('gnomad-exome'), where),
        'alphamissense_max': top_am['pathogenicity'] if top_am else None,
        'alphamissense_max_class': top_am.get('classification') if top_am else None,
        'revel': v.get('revel', {}).get('score'),
        'spliceai_max': max(splice_scores) if splice_scores else None,
        'primateai3d_max': top_p3d['score'] if top_p3d else None,
        'primateai3d_max_percentile': top_p3d.get('scorePercentile') if top_p3d else None,
        'primateai3d_max_class': top_p3d.get('classification') if top_p3d else None,
        'promoterai_min': min(promoterai_scores, default=None),
        'promoterai_max': max(promoterai_scores, default=None),
        'dann': v.get('dannScore'),
        'phylop': v.get('phylopScore'),
        'phylop_primate': v.get('phyloPPrimateScore'),
        'gerp': v.get('gerpScore'),
        **variant_clinvar(v.get('clinvar', ()), where),
    }
    return row, transcript_rows, keys


def _convert_chunk(lines: list[bytes], index: int, output: pathlib.Path, hgnc: Mapping[str, int | None]) -> set[str]:
    """Converts a chunk of position lines into one part of each table; returns the keys of the sources seen."""
    variants, transcripts, keys = [], [], set()
    for line in lines:
        row, rows, seen = position_rows(json.loads(line.rstrip(b',\r\n')), hgnc)
        variants.append(row)
        transcripts.extend(rows)
        keys |= seen
    for table, schema, rows in (('variants', VARIANTS, variants), ('transcripts', TRANSCRIPTS, transcripts)):
        pq.write_table(
            pa.Table.from_pylist(rows, schema=schema),
            output / table / f'part-{index:06d}.parquet',
            compression='zstd',
            **column_options(schema, _ENCODINGS[table]),
        )
    return keys


_HEADER_END = b',"positions":['


def read_header(line: bytes) -> dict[str, Any]:
    text = line.rstrip(b'\r\n')
    if not (text.startswith(b'{"header":') and text.endswith(_HEADER_END)):
        raise ValueError('the first line is not the header of an Annotator JSON')
    return json.loads(text[: -len(_HEADER_END)] + b'}')['header']


_POSITION_START = b'{"chromosome"'
_GENES_START = b'],"genes":['


def position_chunks(lines: Iterator[bytes], size: int) -> Iterator[list[bytes]]:
    """The position lines in chunks, up to the line that closes the positions."""
    chunk: list[bytes] = []
    for line in lines:
        if not line.startswith(_POSITION_START):
            break
        chunk.append(line)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def read_genes(lines: Iterator[bytes]) -> list[dict[str, Any]]:
    """The records of the genes section, one a line, from the lines that follow the document's header.

    A section starts on the line that closes the one before it, and a document without a genes section has no
    records. Fails unless the positions are closed and the document ends as a whole one does.
    """
    closing = next((line for line in lines if not line.startswith(_POSITION_START)), b'')
    if not closing.startswith(b']'):
        raise ValueError('the positions are not closed: the document is cut short or malformed')
    in_genes = closing.startswith(_GENES_START)
    genes = []
    last = closing
    for line in lines:
        last = line
        if line.startswith(b']'):
            in_genes = line.startswith(_GENES_START)
        elif in_genes and line.startswith(b'{'):
            genes.append(json.loads(line.rstrip(b',\r\n')))
    if not last.rstrip(b'\r\n').endswith(b']}'):
        raise ValueError('the document does not end as a whole Annotator JSON does')
    return genes


def read_union_sites(path: pathlib.Path) -> pa.Table:
    """The union's alleles, each by its variant id; a record with other than one alternate allele fails."""
    with pysam.VariantFile(str(path)) as vcf:
        if len(vcf.header.samples):
            raise ValueError(f'{path}: the union is sites-only, and this VCF has samples')
        ids = [record_id(record, str(path)) for record in vcf]
    return pa.table({'variant_id': pa.array(ids, _STRING)})


def check_alleles(sites: pa.Table, annotated: pathlib.Path, name: str) -> None:
    """Fails unless the alleles of the Parquet files are exactly the union's, each once: they annotate this union.

    `annotated` is a Parquet file or a pattern of several with a `variant_id` column, and `name` what the messages
    call it.
    """
    con = duckdb.connect()
    con.register('sites', sites)
    con.read_parquet(str(annotated)).select('variant_id').create_view('annotated')
    checks = {
        'the union repeats alleles': 'SELECT variant_id FROM sites GROUP BY 1 HAVING count(*) > 1',
        f'the {name} repeats alleles': 'SELECT variant_id FROM annotated GROUP BY 1 HAVING count(*) > 1',
        f'union alleles the {name} lacks': 'SELECT variant_id FROM sites ANTI JOIN annotated USING (variant_id)',
        f'{name} alleles the union lacks': 'SELECT variant_id FROM annotated ANTI JOIN sites USING (variant_id)',
    }
    for problem, query in checks.items():
        found = [r[0] for r in con.sql(f'{query} ORDER BY 1 LIMIT 5').fetchall()]
        if found:
            raise ValueError(f'{problem}, such as {found}')


def _omission(source: str) -> OmissionReason:
    if source in DATED:
        return 'dated'
    if source in NOT_APPLICABLE:
        return 'not_applicable'
    if source in CURRENT:
        return 'current'
    return 'unplaced'


def annotation_info(
    header: Mapping[str, Any], run_config: Mapping[str, Any], bundle_config: Mapping[str, Any], catalog: str
) -> AnnotationInfo:
    """What annotated the union, from the JSON header, the run config and the bundle's run config.

    Fails if either config is of another catalog than the one named, if the run config names a source the bundle lacks
    or pins another version of one, or if the header holds a dated source that the run did not name.
    """
    for name, config in (('run config', run_config), ("bundle's run config", bundle_config)):
        if config['data']['catalog'] != catalog:
            raise ValueError(f'the {name} is of catalog {config["data"]["catalog"]}, not {catalog}')
    selected: Mapping[str, Mapping[str, str]] = run_config['data']['versionSelector']
    bundle: Mapping[str, Mapping[str, str]] = bundle_config['data']['versionSelector']
    if extra := sorted(set(selected) - set(bundle)):
        raise ValueError(f'the run config names sources the bundle lacks: {extra}')
    for source, types in selected.items():
        if types != bundle[source]:
            raise ValueError(f'{source}: the run config pins {types}, the bundle {bundle[source]}')
    data_sources = tuple(
        DataSource(name=d['name'], version=str(d['version']), release_date=d.get('releaseDate'))
        for d in header['dataSources']
    )
    in_header = {d.name for d in data_sources}
    for source, names in _DATED_HEADER_NAMES.items():
        if source not in selected and (found := sorted(in_header & set(names))):
            raise ValueError(f'the JSON header holds {found}, of {source!r}, which the run did not name')
    return AnnotationInfo(
        annotator=header['annotator'],
        json_schema_version=header['schemaVersion'],
        assembly=header['genomeAssembly'],
        catalog=catalog,
        sources=tuple(sorted(selected)),
        source_versions=tuple(
            SourceVersion(source=s, annotation_type=t, version=str(version))
            for s, types in sorted(selected.items())
            for t, version in sorted(types.items())
        ),
        data_sources=data_sources,
        omitted=tuple(OmittedSource(source=s, reason=_omission(s)) for s in sorted(set(bundle) - set(selected))),
    )


def convert(
    json_path: pathlib.Path,
    sites_path: pathlib.Path,
    run_config: Mapping[str, Any],
    bundle_config: Mapping[str, Any],
    catalog: str,
    output: pathlib.Path,
    workers: int,
    chunk_size: int = 10_000,
) -> AnnotationInfo:
    """Writes `variants/` and `transcripts/` under the output folder from the union's Annotator JSON.

    Checks them against the union's sites, and returns what annotated the union.
    """
    for table in ('variants', 'transcripts'):
        (output / table).mkdir(parents=True)
    keys: set[str] = set()
    with gzip.open(json_path, 'rb') as document:
        lines = iter(document)
        header = read_header(next(lines))
        if header.get('samples'):
            raise ValueError(f'{json_path}: the union is sites-only, and this annotation has samples')
        info = annotation_info(header, run_config, bundle_config, catalog)
        hgnc = hgnc_ids(read_genes(lines))
    with gzip.open(json_path, 'rb') as document:
        lines = iter(document)
        next(lines)
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            pending: set[concurrent.futures.Future[set[str]]] = set()
            for index, chunk in enumerate(position_chunks(lines, chunk_size)):
                if len(pending) >= 2 * workers:
                    done, pending = concurrent.futures.wait(pending, return_when=concurrent.futures.FIRST_COMPLETED)
                    for future in done:
                        keys |= future.result()
                pending.add(pool.submit(_convert_chunk, chunk, index, output, hgnc))
            for future in concurrent.futures.as_completed(pending):
                keys |= future.result()
    for source, key in _SOURCE_KEYS.items():
        if key in keys and source not in info.sources:
            raise ValueError(f'the JSON holds {key!r} records although the run did not name {source!r}')
    check_alleles(read_union_sites(sites_path), output / 'variants' / '*.parquet', 'JSON')
    return info
