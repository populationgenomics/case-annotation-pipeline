"""Illumina's annotation of a group's structural variants as Parquet tables (docs/annotation.md, the family tables).

The Annotator reads the union of the group's sites (sv_sites), one allele per record, and writes one position line per
record, which echoes the record's ID. This module reads the document line by line, as annotation_tables reads the one
of the small variants, and writes two tables of the union, each a folder of parts:

- `sv`, one row per site: how the Annotator treated it, the interval it annotated, its genes and those that one of
  its breakpoints falls inside, and what it matched in gnomAD's structural variants and copy-number variants, in the
  1000 Genomes structural variants and in ClinVar;
- `sv_transcripts`, one row per site and transcript the Annotator lists for it, with the partners of each gene fusion
  that the Annotator reports for the transcript.

The Annotator annotates a record in one of two ways. A symbolic allele, a breakend and a sequence of more than 50
bases are structural: its sources are matched by overlap and sit on the position. A shorter sequence is a small
variant: its sources are those of the allele itself and sit on the variant, as for the family's small variants. The
Annotator takes the type of a structural record from its SVTYPE and that of a small variant from its alleles alone.

How it matches by overlap is not documented, and this module holds it to what its output shows. A record with an
extent (a deletion, a duplication, an inversion, a copy-number segment) gets an entry of gnomAD, 1000 Genomes or
ClinVar only at a reciprocal overlap of 0.5 or more and never one of the other direction, a loss for a gain. An
insertion or a breakend gets entries near its position, without an overlap. An entry that breaks this fails the
conversion, because the columns would no longer mean what docs/family-tables.md says. So does a key the module does
not know, which can be a source nobody has placed. Nothing here depends on the order of a record's transcripts or
entries, which differs between two runs over the same input.
"""

from __future__ import annotations

import gzip
import json
import pathlib
from collections.abc import Iterable, Mapping, Sequence
from typing import Any, Literal

import pyarrow as pa
import pyarrow.parquet as pq

from case_annotation import annotation_tables, sv_sites

type AnnotatedAs = Literal['structural', 'small_variant']

# The reciprocal overlap below which the Annotator reports no entry of a source it matches by overlap.
RECIPROCAL_OVERLAP = 0.5

_STRING = pa.string()
_INT = pa.int32()
# The Annotator gives a frequency to six decimal places and an overlap to five, which a 32-bit float holds.
_FLOAT = pa.float32()
# An HGNC id as a number, as in the tables of the small variants.
_HGNC_ID = pa.int32()

# One entry of gnomAD's structural variants that a site matched. pos and end_pos give its interval as a table of
# structural variants gives a record's: the bases concerned run from pos + 1 to end_pos.
GNOMAD_SV_MATCH = pa.struct(
    [
        pa.field('id', _STRING),
        pa.field('variant_type', _STRING),
        pa.field('pos', _INT),
        pa.field('end_pos', _INT),
        pa.field('af', _FLOAT),
        pa.field('ac', _INT),
        pa.field('an', _INT),
        pa.field('hom', _INT),
        pa.field('failed_filter', pa.bool_()),
        pa.field('reciprocal_overlap', _FLOAT),
        pa.field('annotation_overlap', _FLOAT),
    ]
)
# One ClinVar record that a site matched: the fields of a small variant's record, then the record's interval and the
# overlaps, which are empty for a record of the allele itself.
CLINVAR_SV_RECORD = pa.struct(
    [
        *annotation_tables.CLINVAR_RECORD,
        pa.field('pos', _INT),
        pa.field('end_pos', _INT),
        pa.field('reciprocal_overlap', _FLOAT),
        pa.field('annotation_overlap', _FLOAT),
    ]
)

# One partner of a gene fusion that the Annotator reports for a transcript, as the Annotator gives it: the partner
# transcript and its gene, and the exon or the intron of the partner that holds the junction.
FUSION_PARTNER = pa.struct(
    [
        pa.field('transcript_id', _STRING),
        pa.field('gene_symbol', _STRING),
        pa.field('hgnc_id', _HGNC_ID),
        pa.field('exon', _STRING),
        pa.field('intron', _STRING),
        pa.field('in_frame', pa.bool_()),
        pa.field('imprecise', pa.bool_()),
        pa.field('directionality', _STRING),
        pa.field('hgvsr', _STRING),
    ]
)

# The columns the annotation adds to a family's sv table.
COLUMNS = pa.schema(
    [
        pa.field('annotated_as', _STRING, nullable=False),
        pa.field('annotated_pos', _INT),
        pa.field('annotated_end_pos', _INT),
        pa.field('hgnc_ids', pa.list_(_HGNC_ID)),
        pa.field('breakpoint_hgnc_ids', pa.list_(_HGNC_ID)),
        pa.field('gnomad_sv_id', _STRING),
        pa.field('gnomad_sv_af', _FLOAT),
        pa.field('gnomad_sv_ac', _INT),
        pa.field('gnomad_sv_an', _INT),
        pa.field('gnomad_sv_hom', _INT),
        pa.field('gnomad_sv_failed_filter', pa.bool_()),
        pa.field('gnomad_sv_overlap', _FLOAT),
        pa.field('gnomad_sv_distance', _INT),
        pa.field('gnomad_sv_max_af', _FLOAT),
        pa.field('gnomad_sv_matches', pa.list_(GNOMAD_SV_MATCH)),
        pa.field('onekg_sv_id', _STRING),
        pa.field('onekg_sv_af', _FLOAT),
        pa.field('onekg_sv_ac', _INT),
        pa.field('onekg_sv_an', _INT),
        pa.field('onekg_sv_overlap', _FLOAT),
        pa.field('onekg_sv_distance', _INT),
        pa.field('gnomad_cnv_id', _STRING),
        pa.field('gnomad_cnv_sf', _FLOAT),
        pa.field('gnomad_cnv_failed_filter', pa.bool_()),
        pa.field('gnomad_cnv_overlap', _FLOAT),
        *annotation_tables.gnomad_fields('gnomad_genomes'),
        *annotation_tables.gnomad_fields('gnomad_exomes'),
        pa.field('clinvar_accession', _STRING),
        pa.field('clinvar_classification', _STRING),
        pa.field('clinvar_review_status', _STRING),
        pa.field('clinvar_stars', pa.int8()),
        pa.field('clinvar_conditions', pa.list_(_STRING)),
        pa.field('clinvar_last_evaluated', _STRING),
        pa.field('clinvar_overlap', _FLOAT),
        pa.field('clinvar_records', pa.list_(CLINVAR_SV_RECORD)),
    ]
)
# The columns that hold the matches of a site in full. A family's record of two alleles keeps these alone, each
# holding the matches of both alleles: a summary would be that of one allele, and the record's members carry either.
LISTS = ('gnomad_sv_matches', 'clinvar_records')

# The union's sv table: a site, the Annotator's type of its allele, and the columns above.
SV = pa.schema(
    [
        pa.field('site_id', _STRING, nullable=False),
        pa.field('chrom', _STRING, nullable=False),
        pa.field('pos', _INT, nullable=False),
        pa.field('variant_type', _STRING, nullable=False),
        *COLUMNS,
    ]
)
# One transcript that the Annotator lists for a site: an entry of a record's annotation (position_rows). The union's
# sv_transcripts table holds each entry under its site, and a family's table under its record's place and ID.
TRANSCRIPT_COLUMNS = pa.schema(
    [
        pa.field('variant_type', _STRING, nullable=False),
        pa.field('transcript_id', _STRING, nullable=False),
        pa.field('source', _STRING, nullable=False),
        pa.field('biotype', _STRING),
        pa.field('gene_symbol', _STRING),
        pa.field('gene_id', _STRING),
        pa.field('hgnc_id', _HGNC_ID),
        pa.field('consequences', pa.list_(_STRING), nullable=False),
        pa.field('protein_id', _STRING),
        pa.field('is_mane_select', pa.bool_(), nullable=False),
        pa.field('is_canonical', pa.bool_(), nullable=False),
        pa.field('complete_overlap', pa.bool_(), nullable=False),
        pa.field('exon', _STRING),
        pa.field('intron', _STRING),
        pa.field('cdna_pos', _STRING),
        pa.field('protein_pos', _STRING),
        pa.field('fusion_partners', pa.list_(FUSION_PARTNER)),
    ]
)
SV_TRANSCRIPTS = pa.schema([pa.field('site_id', _STRING, nullable=False), *TRANSCRIPT_COLUMNS])

# The columns of a family's sv table that hold a dated source's data, by the source's name in a run config. The table
# leaves out the columns of a source its run did not name, as the variants table does.
SOURCE_COLUMNS: Mapping[str, tuple[str, ...]] = {
    'clinvar': tuple(name for name in COLUMNS.names if name.startswith('clinvar_')),
}

# The keys of a position and of a variant that are no source's.
_POSITION_KEYS = frozenset(
    {
        'chromosome', 'position', 'id', 'refAllele', 'altAlleles', 'vcfInfo', 'cytogeneticBand', 'variants',
        'svEnd', 'svLength', 'ciPos', 'ciEnd',
    }
)  # fmt: skip
_VARIANT_KEYS = frozenset(
    {
        'vid', 'chromosome', 'begin', 'end', 'refAllele', 'altAllele', 'variantType', 'isStructuralVariant',
        'hgvsg', 'inLowComplexityRegion', 'simpleNomenclature', 'transcripts', 'regulatoryRegions',
    }
)  # fmt: skip
# The sources the Annotator matches to a structural record by overlap, by their key on a position. The first four are
# read; the others are known and left out (docs/annotation.md says why).
_POSITION_SOURCES = frozenset(
    {'gnomad', 'gnomad-cnv', 'oneKg', 'clinvar', 'clingen', 'clingenDosageSensitivityMap', 'decipher', 'mitomap'}
)
# The sources of an allele itself, which the Annotator gives a record it treats as a small variant, by their key on a
# variant. gnomAD's and ClinVar's are read, as for the family's small variants; the scores and the other frequencies
# are in the variants table for a record that names its row there.
_VARIANT_SOURCES = frozenset(
    {
        'gnomad', 'gnomad-exome', 'clinvar', 'dbsnp', 'mitomap', 'oneKg', 'topmed', 'abraom', 'gmeVariome',
        'globalAllele', 'alphaMissense', 'spliceAI', 'primateAI-3D', 'primateAI', 'promoterAI', 'revel',
        'dannScore', 'gerpScore', 'phylopScore', 'phyloPPrimateScore',
    }
)  # fmt: skip
# The keys among them that hold a dated source's data (docs/annotation.md, dated sources), each with the source that a
# run config names it by: a document holding one although the run did not name the source is refused.
_DATED_KEYS: Mapping[str, Mapping[str, str]] = {
    'position': {
        'clinvar': 'clinvar',
        'clingen': 'clingen (legacy)',
        'clingenDosageSensitivityMap': 'clingen',
        'mitomap': 'mitomap',
    },
    'variant': {'clinvar': 'clinvar', 'dbsnp': 'dbSNP', 'mitomap': 'mitomap'},
}
# The sources among them that this module reads for a small variant, and those it reads by overlap.
_ALLELE_SOURCES = ('gnomad', 'gnomad-exome', 'clinvar')
_OVERLAP_SOURCES = ('gnomad', 'gnomad-cnv', 'oneKg', 'clinvar')

# The Annotator's variant types: those of a record with an extent, which it matches by overlap, with the direction of
# each; those of a record at one position, which it matches by distance; and those it gives a small variant.
_DIRECTION: Mapping[str, str | None] = {
    'deletion': 'loss',
    'copy_number_loss': 'loss',
    'duplication': 'gain',
    'tandem_duplication': 'gain',
    'copy_number_gain': 'gain',
    'inversion': None,
}
_POINT_TYPES = frozenset({'insertion', 'translocation'})
# The types a record can have as a structural variant, by the record's SVTYPE and, for a copy-number segment, its
# allele.
_STRUCTURAL_TYPES: Mapping[tuple[str, str | None], frozenset[str]] = {
    ('DEL', None): frozenset({'deletion'}),
    ('INS', None): frozenset({'insertion'}),
    ('DUP', None): frozenset({'tandem_duplication', 'duplication'}),
    ('INV', None): frozenset({'inversion'}),
    ('BND', None): frozenset({'translocation'}),
    ('CNV', '<DEL>'): frozenset({'copy_number_loss'}),
    ('CNV', '<DUP>'): frozenset({'copy_number_gain'}),
}
# The types a record can have as a small variant, whatever its SVTYPE. The Annotator trims the bases that the two
# alleles share and types what is left, without reading SVTYPE. A caller's deletion that inserts as many bases as it
# deletes is therefore an MNV, or an SNV if one base differs.
_SMALL_TYPES = frozenset({'SNV', 'MNV', 'insertion', 'deletion', 'delins'})
# The terms of a transcript that no breakpoint of the record falls inside: one that the record takes or copies whole,
# which a deletion written as a sequence has without completeOverlap, and one that the record lies beside.
_NO_BREAKPOINT_TERMS = frozenset(
    {'transcript_ablation', 'transcript_amplification', 'upstream_gene_variant', 'downstream_gene_variant'}
)


def _need(held: Mapping[str, Any], key: str, where: str) -> Any:  # noqa: ANN401 (a value of the Annotator's JSON)
    if held.get(key) is None:
        raise ValueError(f'{where}: no {key}')
    return held[key]


def _entries(held: object, where: str) -> list[Mapping[str, Any]]:
    """The entries of a source, each once: the Annotator repeats an identical entry now and then."""
    if not isinstance(held, list) or not all(isinstance(entry, dict) for entry in held):
        raise ValueError(f'{where}: not a list of entries')
    seen: set[str] = set()
    entries = []
    for entry in held:
        text = json.dumps(entry, sort_keys=True)
        if text not in seen:
            seen.add(text)
            entries.append(entry)
    return entries


def _variant_type(site: sv_sites.Site, variant: Mapping[str, Any], where: str) -> tuple[str, AnnotatedAs]:
    """The Annotator's type of a site's allele and how it treated the site; a type the site cannot have fails.

    A structural record can have the types of its SVTYPE, and a small variant any type of a small variant: its type
    says what its alleles are, which a caller can have filed under another SVTYPE.
    """
    variant_type = _need(variant, 'variantType', where)
    sv_type = site.sv_type
    symbolic = site.alt.startswith('<') or sv_sites.is_breakend(site.alt)
    if variant.get('isStructuralVariant', False):
        allowed = _STRUCTURAL_TYPES.get((sv_type, site.alt if sv_type == 'CNV' else None), frozenset())
        annotated_as: AnnotatedAs = 'structural'
    elif symbolic:
        raise ValueError(f'{where}: the allele {site.alt} is not annotated as a structural variant')
    else:
        allowed = _SMALL_TYPES
        annotated_as = 'small_variant'
    if variant_type not in allowed:
        raise ValueError(
            f'{where}: the variant type {variant_type!r} of a record of SVTYPE {sv_type} with the allele '
            f'{site.alt[:40]!r}, annotated as {annotated_as}, is none this module knows for it'
        )
    return variant_type, annotated_as


def _check_site(position: Mapping[str, Any], site: sv_sites.Site, where: str) -> None:
    """Fails unless the position is the site's record, of one variant, and holds no key this module does not know."""
    echoed = (
        position.get('chromosome'),
        position.get('position'),
        position.get('refAllele'),
        position.get('altAlleles'),
    )
    if echoed != (site.chrom, site.pos, site.ref, [site.alt]):
        raise ValueError(f'{where}: the position is {echoed}, and its site {site.line()}')
    if unknown := sorted(position.keys() - _POSITION_KEYS - _POSITION_SOURCES):
        raise ValueError(f'{where}: the keys {unknown} of a position, which this module does not know')
    variants = position.get('variants')
    if not isinstance(variants, list) or len(variants) != 1:
        raise ValueError(f'{where}: not one variant, and the union holds one allele for each record')
    if unknown := sorted(variants[0].keys() - _VARIANT_KEYS - _VARIANT_SOURCES):
        raise ValueError(f'{where}: the keys {unknown} of a variant, which this module does not know')


def _interval(
    position: Mapping[str, Any], site: sv_sites.Site, annotated_as: AnnotatedAs, where: str
) -> dict[str, int]:
    """The interval the Annotator annotated, as pos and end_pos give a record's.

    A structural record's end is its site's END, which the position repeats as svEnd, whatever the Annotator did to
    the start: it trims and left-aligns a sequence allele. Without that END it annotates a single base, which is why
    a structural record whose end is another fails. A breakend is the one base at its position, and a small variant
    has the interval of its trimmed allele and no svEnd.
    """
    variant = position['variants'][0]
    begin, end = _need(variant, 'begin', where), _need(variant, 'end', where)
    structural_end = site.end if annotated_as == 'structural' else None
    if position.get('svEnd') != structural_end:
        raise ValueError(
            f'{where}: svEnd is {position.get("svEnd")}, and a record annotated as {annotated_as} whose site has '
            f'the END {site.end} has {structural_end}'
        )
    if sv_sites.is_breakend(site.alt):
        if (begin, end) != (site.pos, site.pos):
            raise ValueError(f'{where}: a breakend at {site.pos} annotated from {begin} to {end}')
        return {'annotated_pos': site.pos, 'annotated_end_pos': site.pos}
    if annotated_as == 'structural' and end != site.end:
        raise ValueError(f'{where}: the variant ends at {end}, and its site has the END {site.end}')
    return {'annotated_pos': begin - 1, 'annotated_end_pos': end}


def _checked(
    entries: Sequence[Mapping[str, Any]], source: str, variant_type: str, annotated_as: AnnotatedAs, where: str
) -> bool:
    """Whether the entries of a source matched by overlap have overlaps; fails unless they are matched as documented.

    A record with an extent has entries at a reciprocal overlap of RECIPROCAL_OVERLAP or more, none of them of the
    other direction. An insertion or a breakend has entries without an overlap. A small variant has none at all.
    """
    if not entries:
        return False
    where = f'{where}: {source}'
    if annotated_as != 'structural':
        raise ValueError(f'{where}: an entry matched by overlap, on a record annotated as a small variant')
    if variant_type in _POINT_TYPES:
        if any('reciprocalOverlap' in entry for entry in entries):
            raise ValueError(
                f'{where}: an entry with an overlap on a record of type {variant_type}, which has no extent'
            )
        return False
    direction = _DIRECTION[variant_type]
    for entry in entries:
        overlap = entry.get('reciprocalOverlap')
        if not isinstance(overlap, int | float) or not isinstance(entry.get('annotationOverlap'), int | float):
            raise ValueError(f'{where}: an entry without its overlaps on a record of type {variant_type}')
        if overlap < RECIPROCAL_OVERLAP:
            raise ValueError(f'{where}: an entry at a reciprocal overlap of {overlap}, below {RECIPROCAL_OVERLAP}')
        other = _DIRECTION.get(entry.get('variantType', ''))
        if direction is not None and other is not None and other != direction:
            raise ValueError(f'{where}: an entry of type {entry["variantType"]} on a record of type {variant_type}')
    return True


def _distance(entry: Mapping[str, Any], site: sv_sites.Site, where: str) -> int:
    """The bases between a site's position and an entry that was matched without an overlap.

    Both are places between two bases or short intervals, read as pos and end_pos give them: an insertion entry at
    the site's own position is no bases away.
    """
    return max(_need(entry, 'begin', where) - 1 - site.pos, site.pos - _need(entry, 'end', where), 0)


def _ranked(
    entries: Sequence[Mapping[str, Any]],
    overlapping: bool,
    site: sv_sites.Site,
    frequency: str,
    identifier: str,
    where: str,
) -> list[Mapping[str, Any]]:
    """The entries from the best match on.

    The best match of a record with an extent has the highest reciprocal overlap, and that of an insertion or a
    breakend is the nearest. Among equals the more frequent comes first, then the one whose identifier sorts first.
    """

    def key(entry: Mapping[str, Any]) -> tuple[float, float, str]:
        closeness = -entry['reciprocalOverlap'] if overlapping else _distance(entry, site, where)
        held = entry.get(frequency)
        return closeness, -(held if held is not None else -1.0), str(entry.get(identifier) or '')

    return sorted(entries, key=key)


def _gnomad_sv(
    entries: Sequence[Mapping[str, Any]], overlapping: bool, site: sv_sites.Site, where: str
) -> dict[str, Any]:
    ranked = _ranked(entries, overlapping, site, 'allAf', 'variantId', where)
    if not ranked:
        return {}
    matches = [
        {
            'id': _need(entry, 'variantId', where),
            'variant_type': entry.get('variantType'),
            'pos': _need(entry, 'begin', where) - 1,
            'end_pos': _need(entry, 'end', where),
            'af': entry.get('allAf'),
            'ac': entry.get('allAc'),
            'an': entry.get('allAn'),
            'hom': entry.get('allHc'),
            'failed_filter': bool(entry.get('failedFilter', False)),
            'reciprocal_overlap': entry.get('reciprocalOverlap'),
            'annotation_overlap': entry.get('annotationOverlap'),
        }
        for entry in ranked
    ]
    best = matches[0]
    frequencies = [match['af'] for match in matches if match['af'] is not None]
    return {
        'gnomad_sv_id': best['id'],
        'gnomad_sv_af': best['af'],
        'gnomad_sv_ac': best['ac'],
        'gnomad_sv_an': best['an'],
        'gnomad_sv_hom': best['hom'],
        'gnomad_sv_failed_filter': best['failed_filter'],
        'gnomad_sv_overlap': best['reciprocal_overlap'],
        'gnomad_sv_distance': None if overlapping else _distance(ranked[0], site, where),
        'gnomad_sv_max_af': max(frequencies) if frequencies else None,
        'gnomad_sv_matches': matches,
    }


def _onekg_sv(
    entries: Sequence[Mapping[str, Any]], overlapping: bool, site: sv_sites.Site, where: str
) -> dict[str, Any]:
    ranked = _ranked(entries, overlapping, site, 'allAf', 'id', where)
    if not ranked:
        return {}
    best = ranked[0]
    return {
        'onekg_sv_id': _need(best, 'id', where),
        'onekg_sv_af': best.get('allAf'),
        'onekg_sv_ac': best.get('allAc'),
        'onekg_sv_an': best.get('allAn'),
        'onekg_sv_overlap': best.get('reciprocalOverlap'),
        'onekg_sv_distance': None if overlapping else _distance(best, site, where),
    }


def _gnomad_cnv(
    entries: Sequence[Mapping[str, Any]], overlapping: bool, site: sv_sites.Site, where: str
) -> dict[str, Any]:
    if entries and not overlapping:
        raise ValueError(f'{where}: gnomad-cnv: an entry on a record without an extent')
    ranked = _ranked(entries, overlapping, site, 'allSf', 'variantId', where)
    if not ranked:
        return {}
    best = ranked[0]
    return {
        'gnomad_cnv_id': _need(best, 'variantId', where),
        'gnomad_cnv_sf': best.get('allSf'),
        'gnomad_cnv_failed_filter': bool(best.get('failedFilter', False)),
        'gnomad_cnv_overlap': best['reciprocalOverlap'],
    }


def _clinvar(entries: Sequence[Mapping[str, Any]], where: str) -> dict[str, Any]:
    """The ClinVar records a structural record matched, all of them, and a summary from the first.

    The records are those of other variants that overlap this one, so the best-overlapping record comes first, and
    among equals the best reviewed (annotation_tables.by_review).
    """
    records = [
        annotation_tables.clinvar_record(entry)
        | {
            'pos': _need(entry, 'begin', where) - 1,
            'end_pos': _need(entry, 'end', where),
            'reciprocal_overlap': entry['reciprocalOverlap'],
            'annotation_overlap': entry['annotationOverlap'],
        }
        for entry in entries
    ]
    records = sorted(annotation_tables.by_review(records), key=lambda r: -r['reciprocal_overlap'])
    columns = annotation_tables.clinvar_columns(records)
    return columns | ({'clinvar_overlap': records[0]['reciprocal_overlap']} if records else {})


def _fusion_partners(
    transcript: Mapping[str, Any], hgnc: Mapping[str, int | None], where: str
) -> list[dict[str, Any]] | None:
    """The partners of the gene fusions that the Annotator reports for a transcript, by partner transcript.

    Each is copied as the Annotator gives it (FUSION_PARTNER), with the HGNC id of its gene by the genes section. The
    Annotator gives the exon or the intron of the partner that holds the junction. A partner without inFrame or
    isImprecise has the flag false, as a transcript without isCanonical has. An entry without a key that the column
    reads fails, and so does a partner whose gene symbol the genes section lacks.
    """
    partners = []
    for entry in _entries(transcript.get('geneFusions', []), f'{where}: geneFusions'):
        partner_id = _need(entry, 'transcript', f'{where}: a fusion partner')
        partner = f'{where}: the fusion partner {partner_id}'
        symbol = _need(entry, 'hgnc', partner)
        if symbol not in hgnc:
            raise ValueError(f'{partner}: the genes section lacks {symbol!r}, its gene')
        if entry.get('exon') is None and entry.get('intron') is None:
            raise ValueError(f'{partner}: neither exon nor intron')
        partners.append(
            {
                'transcript_id': partner_id,
                'gene_symbol': symbol,
                'hgnc_id': hgnc[symbol],
                'exon': entry.get('exon'),
                'intron': entry.get('intron'),
                'in_frame': bool(entry.get('inFrame', False)),
                'imprecise': bool(entry.get('isImprecise', False)),
                'directionality': _need(entry, 'directionality', partner),
                'hgvsr': _need(entry, 'hgvsr', partner),
            }
        )
    return sorted(partners, key=lambda partner: (partner['transcript_id'], partner['hgvsr'])) or None


def _holds_breakpoint(transcript: Mapping[str, Any]) -> bool:
    """Whether a breakpoint of the record falls inside the transcript of an entry (TRANSCRIPT_COLUMNS).

    The Annotator lists a transcript for a record in one of three ways. The record covers the transcript whole: the
    entry has completeOverlap, or the term of an ablation or an amplification, which a deletion written as a sequence
    has without the flag. The record lies beside the transcript: the entry has an upstream or a downstream term. Or
    the record overlaps a part of the transcript, which is every other entry, and then an end of the record's interval
    lies inside the transcript. An insertion and a breakend are annotated at one position, so each of their entries
    that is not beside the transcript is of a transcript that holds that position.
    """
    return not (transcript['complete_overlap'] or _NO_BREAKPOINT_TERMS & set(transcript['consequences']))


def _hgnc_ids(transcripts: Iterable[Mapping[str, Any]]) -> list[int] | None:
    """The HGNC ids of the genes of the transcript entries given, each once; None where no gene of theirs has one."""
    return sorted({t['hgnc_id'] for t in transcripts if t['hgnc_id'] is not None}) or None


def _transcripts(
    variant_type: str, variant: Mapping[str, Any], hgnc: Mapping[str, int | None], where: str
) -> list[dict[str, Any]]:
    """A site's transcript entries, by transcript.

    A transcript listed twice fails, since it is an entry's key, and so does one whose gene symbol the genes section
    lacks, as for a small variant's transcript (annotation_tables.position_rows).
    """
    rows: dict[str, dict[str, Any]] = {}
    for transcript in variant.get('transcripts', ()):
        transcript_id = _need(transcript, 'transcript', where)
        if transcript_id in rows:
            raise ValueError(f'{where}: the transcript {transcript_id} twice')
        symbol = transcript.get('hgnc')
        if symbol is not None and symbol not in hgnc:
            raise ValueError(f'{where}: the genes section lacks {symbol!r}, the gene of the transcript {transcript_id}')
        rows[transcript_id] = {
            'variant_type': variant_type,
            'transcript_id': transcript_id,
            'source': _need(transcript, 'source', where),
            'biotype': transcript.get('bioType'),
            'gene_symbol': symbol,
            'gene_id': transcript.get('geneId'),
            'hgnc_id': None if symbol is None else hgnc[symbol],
            'consequences': annotation_tables.by_severity(transcript.get('consequence', ())),
            'protein_id': transcript.get('proteinId'),
            'is_mane_select': bool(transcript.get('isManeSelect', False)),
            'is_canonical': bool(transcript.get('isCanonical', False)),
            'complete_overlap': bool(transcript.get('completeOverlap', False)),
            'exon': transcript.get('exons'),
            'intron': transcript.get('introns'),
            'cdna_pos': transcript.get('cdnaPos'),
            'protein_pos': transcript.get('proteinPos'),
            'fusion_partners': _fusion_partners(transcript, hgnc, f'{where}: the transcript {transcript_id}'),
        }
    return [rows[transcript_id] for transcript_id in sorted(rows)]


def position_rows(
    position: Mapping[str, Any], sites: Mapping[str, sv_sites.Site], hgnc: Mapping[str, int | None]
) -> tuple[dict[str, Any], list[dict[str, Any]], set[str], set[str]]:
    """One position line as its site's annotation, and the source keys on the position and on its variant.

    The annotation is the site's row of `sv` and its transcript entries (TRANSCRIPT_COLUMNS), whatever tables they
    are then laid out in. `hgnc` holds each gene symbol of the document's genes section with its HGNC id
    (annotation_tables.hgnc_ids).
    """
    site_id = position.get('id')
    if not isinstance(site_id, str) or site_id not in sites:
        raise ValueError(
            f'{position.get("chromosome")}:{position.get("position")}: the id {site_id!r} is that of no site of the '
            'union'
        )
    site = sites[site_id]
    where = f'{site.chrom}:{site.pos} {site_id}'
    _check_site(position, site, where)
    variant = position['variants'][0]
    variant_type, annotated_as = _variant_type(site, variant, where)
    row: dict[str, Any] = {
        'site_id': site_id,
        'chrom': site.chrom,
        'pos': site.pos,
        'variant_type': variant_type,
        'annotated_as': annotated_as,
        **_interval(position, site, annotated_as, where),
    }
    matched = {source: _entries(position.get(source, []), f'{where}: {source}') for source in _OVERLAP_SOURCES}
    overlapping = {
        source: _checked(entries, source, variant_type, annotated_as, where) for source, entries in matched.items()
    }
    row |= _gnomad_sv(matched['gnomad'], overlapping['gnomad'], site, where)
    row |= _onekg_sv(matched['oneKg'], overlapping['oneKg'], site, where)
    row |= _gnomad_cnv(matched['gnomad-cnv'], overlapping['gnomad-cnv'], site, where)
    own = [source for source in _ALLELE_SOURCES if source in variant]
    if annotated_as == 'structural':
        if own:
            raise ValueError(f'{where}: {own} of the allele itself, on a record annotated as a structural variant')
        if matched['clinvar'] and not overlapping['clinvar']:
            raise ValueError(f'{where}: clinvar: an entry without an overlap, which this module cannot order')
        row |= _clinvar(matched['clinvar'], where)
    else:
        row |= annotation_tables.gnomad_values('gnomad_genomes', variant.get('gnomad'))
        row |= annotation_tables.gnomad_values('gnomad_exomes', variant.get('gnomad-exome'))
        row |= annotation_tables.allele_clinvar(_entries(variant.get('clinvar', []), f'{where}: clinvar'))
    transcripts = _transcripts(variant_type, variant, hgnc, where)
    row['hgnc_ids'] = _hgnc_ids(transcripts)
    row['breakpoint_hgnc_ids'] = _hgnc_ids(t for t in transcripts if _holds_breakpoint(t))
    return (
        row,
        transcripts,
        set(position.keys() & _POSITION_SOURCES),
        set(variant.keys() & _VARIANT_SOURCES),
    )


def record_of_alleles(alleles: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The columns of a family's record of several alleles, from those of its alleles' sites in the record's order.

    The lists hold the matches of every allele, one allele's after the other's, each match with its type. The
    interval is the alleles' where they share one, and the HGNC ids are those of every allele's transcripts, which
    are the record's: all of them in hgnc_ids, and in breakpoint_hgnc_ids those that a breakpoint of an allele falls
    inside. Every other column is empty: a best match or a frequency would be that of one allele, and a member of the
    family may carry either. Alleles that are not all annotated in one way fail, since no value of annotated_as would
    be true of the record.
    """
    kinds = {allele['annotated_as'] for allele in alleles}
    if len(kinds) != 1:
        raise ValueError(f'the alleles of one record are annotated as {sorted(kinds)}')
    row: dict[str, Any] = {'annotated_as': kinds.pop()}
    for name in ('annotated_pos', 'annotated_end_pos'):
        held = {allele[name] for allele in alleles}
        row[name] = held.pop() if len(held) == 1 else None
    for name in ('hgnc_ids', 'breakpoint_hgnc_ids'):
        row[name] = sorted({hgnc_id for allele in alleles for hgnc_id in allele[name] or ()}) or None
    for name in LISTS:
        if any(name in allele for allele in alleles):
            row[name] = [match for allele in alleles for match in allele.get(name) or ()] or None
    return row


def _write_part(rows: Iterable[Mapping[str, Any]], schema: pa.Schema, folder: pathlib.Path, index: int) -> None:
    pq.write_table(
        pa.Table.from_pylist(list(rows), schema=schema), folder / f'part-{index:06d}.parquet', compression='zstd'
    )


def convert(
    json_path: pathlib.Path,
    sites_path: pathlib.Path,
    run_config: Mapping[str, Any],
    bundle_config: Mapping[str, Any],
    catalog: str,
    output: pathlib.Path,
    chunk_size: int = 10_000,
) -> annotation_tables.AnnotationInfo:
    """Writes `sv/` and `sv_transcripts/` under the output folder from the Annotator JSON of a union of sites.

    Fails unless the JSON's positions are exactly the union's sites, each once, and returns what annotated them. A
    dated source's key in the JSON fails if the run did not name the source: the Annotator must not have read it.
    """
    sites = sv_sites.read_union(sites_path)
    for table in ('sv', 'sv_transcripts'):
        (output / table).mkdir(parents=True)
    seen: set[str] = set()
    position_sources: set[str] = set()
    variant_sources: set[str] = set()
    with gzip.open(json_path, 'rb') as document:
        lines = iter(document)
        header = annotation_tables.read_header(next(lines))
        if header.get('samples'):
            raise ValueError(f'{json_path}: the union is sites-only, and this annotation has samples')
        info = annotation_tables.annotation_info(header, run_config, bundle_config, catalog)
        hgnc = annotation_tables.hgnc_ids(annotation_tables.read_genes(lines))
    with gzip.open(json_path, 'rb') as document:
        lines = iter(document)
        next(lines)
        for index, chunk in enumerate(annotation_tables.position_chunks(lines, chunk_size)):
            records, transcripts = [], []
            for line in chunk:
                row, entries, on_position, on_variant = position_rows(json.loads(line.rstrip(b',\r\n')), sites, hgnc)
                if row['site_id'] in seen:
                    raise ValueError(f'{json_path}: the site {row["site_id"]} twice')
                seen.add(row['site_id'])
                records.append(row)
                transcripts.extend({'site_id': row['site_id']} | entry for entry in entries)
                position_sources |= on_position
                variant_sources |= on_variant
            _write_part(records, SV, output / 'sv', index)
            _write_part(transcripts, SV_TRANSCRIPTS, output / 'sv_transcripts', index)
    if missing := sorted(sites.keys() - seen)[:5]:
        raise ValueError(f'{json_path}: sites of the union that the JSON lacks, such as {missing}')
    for level, keys in (('position', position_sources), ('variant', variant_sources)):
        dated = _DATED_KEYS[level]
        if unnamed := sorted(key for key in keys if key in dated and dated[key] not in info.sources):
            raise ValueError(
                f'{json_path}: the JSON holds the {level} keys {unnamed} although the run did not name '
                f'{sorted({dated[key] for key in unnamed})}'
            )
    return info
