"""An invented annotation of structural variants, in the layout Illumina's Annotator writes one in.

Every record, id, frequency and ClinVar record here is invented. The helpers turn the invented callers' records of
structural_vcfs into what the family tables read: a family's sites as the workflow cuts them from a caller's VCF, the
union of a group's sites with its map, the Annotator's JSON of that union, and its tables.
"""

from __future__ import annotations

import copy
import gzip
import json
import pathlib
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import structural_vcfs

from case_annotation import annotation_tables, family_tables, sv_annotation, sv_sites

CATALOG = '209901'
BUNDLE: dict[str, Any] = {
    'data': {
        'catalog': CATALOG,
        'versionSelector': {
            'Ensembl': {'GeneModels': '1'},
            'gnomad': {'SmallVariant': '1', 'StructuralVariant': '1'},
            'oneKg': {'StructuralVariant': '1'},
            'dbSNP': {'SmallVariant': '1'},
            'clinvar': {'SmallVariant': '1', 'StructuralVariant': '1'},
            'clingen': {'StructuralVariant': '1'},
            'clingen (legacy)': {'StructuralVariant': '1'},
            'mitomap': {'StructuralVariant': '1'},
        },
    }
}
ALL = ('Ensembl', 'gnomad', 'oneKg', 'clinvar')
_CONTIGS = ('chr1', 'chr2', 'chr10', 'chrX', 'chrY', 'chrM')


def run_config(*sources: str) -> dict[str, Any]:
    config = copy.deepcopy(BUNDLE)
    config['data']['versionSelector'] = {s: BUNDLE['data']['versionSelector'][s] for s in sources}
    return config


def header(*names: str, samples: list[str] | None = None) -> dict[str, Any]:
    return {
        'annotator': 'Invented Annotator 1.0',
        'creationTime': '2099-01-01 00:00:00',
        'genomeAssembly': 'GRCh38',
        'schemaVersion': 6,
        'dataSources': [{'name': n, 'version': '1', 'releaseDate': '2099-01-01'} for n in names],
        **({'samples': samples} if samples is not None else {}),
    }


def site_line(record: Sequence[str]) -> str:
    """A caller's record as the workflow cuts it for a family's sites: no samples, QUAL or FILTER, INFO cut down."""
    chrom, pos, sv_id, ref, alt, _qual, _filter, info = record[:8]
    kept = [item for item in info.split(';') if item.partition('=')[0] in sv_sites.KEPT]
    return '\t'.join([chrom, pos, sv_id, ref, alt, '.', '.', ';'.join(kept) or '.'])


def write_sites(
    path: pathlib.Path, records: Sequence[Sequence[str]], contigs: Sequence[str] = _CONTIGS
) -> pathlib.Path:
    """A family's sites of one caller: its records with an alternate allele, on the contigs given."""
    lines = [
        '##fileformat=VCFv4.2',
        *(f'##contig=<ID={contig},length=200000000>' for contig in contigs),
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']),
        *(site_line(record) for record in records if record[0] in contigs and record[4] != '.'),
    ]
    with gzip.open(path, 'wt') as vcf:
        vcf.write('\n'.join(lines) + '\n')
    return path


def site(line: str) -> sv_sites.Site:
    """The site of a line of a union or of a family's sites that holds one allele."""
    (found,) = sv_sites.allele_sites(line.split('\t'), 'an invented site')
    return found


def is_structural(found: sv_sites.Site) -> bool:
    """Whether the Annotator treats the site as a structural variant: a symbolic allele, a breakend, a long sequence."""
    return found.alt.startswith('<') or sv_sites.is_breakend(found.alt) or max(len(found.ref), len(found.alt)) > 50


def variant_type(found: sv_sites.Site) -> str:
    if found.sv_type == 'CNV':
        return {'<DEL>': 'copy_number_loss', '<DUP>': 'copy_number_gain'}[found.alt]
    types = {
        'DEL': 'deletion',
        'INS': 'insertion',
        'DUP': 'tandem_duplication',
        'INV': 'inversion',
        'BND': 'translocation',
    }
    return types.get(found.sv_type, 'sequence_alteration')


def transcript(
    transcript_id: str, source: str, gene: str, gene_id: str, *terms: str, **extra: object
) -> dict[str, Any]:
    return {
        'transcript': transcript_id,
        'source': source,
        'bioType': 'mRNA',
        'geneId': gene_id,
        'hgnc': gene,
        'consequence': list(terms),
        'impact': 'high',
        **extra,
    }


def position(found: sv_sites.Site, variant: Mapping[str, Any] | None = None, /, **sources: object) -> dict[str, Any]:
    """The position line of a site, with the sources given on it and the variant's keys given on its variant."""
    structural = is_structural(found)
    breakend = sv_sites.is_breakend(found.alt)
    held: dict[str, Any] = {
        'chromosome': found.chrom,
        'position': found.pos,
        'id': found.site_id,
        'refAllele': found.ref,
        'altAlleles': [found.alt],
        'vcfInfo': {'SVTYPE': found.sv_type},
        'cytogeneticBand': '1p1',
        **({'svEnd': found.end} if structural and not breakend else {}),
        **sources,
    }
    held['variants'] = [
        {
            'vid': 'x',
            'chromosome': found.chrom,
            'begin': found.pos if breakend else found.pos + 1,
            'end': found.pos if breakend else found.end,
            'refAllele': found.ref,
            'altAllele': found.alt,
            'variantType': variant_type(found),
            **({'isStructuralVariant': True} if structural else {}),
            **(variant or {}),
        }
    ]
    return held


# The genes section of an invented annotation: HGNC names GENEA and GENEC, and GENEB is a gene without an HGNC id.
GENES: Sequence[Mapping[str, Any]] = (
    {'name': 'GENEA', 'hgncId': 90001},
    {'name': 'GENEB'},
    {'name': 'GENEC', 'hgncId': 90002},
)
# Each of those genes' symbols with its HGNC id, as annotation_tables.hgnc_ids reads the section.
HGNC: Mapping[str, int | None] = {'GENEA': 90001, 'GENEB': None, 'GENEC': 90002}


def write_json(path: pathlib.Path, held: Mapping[str, Any], closed: bool = True) -> pathlib.Path:
    """An Annotator JSON of the header and positions given, with GENES; without its closing lines if not closed."""
    with gzip.open(path, 'wt') as f:
        f.write('{"header":' + json.dumps(held['header']) + ',"positions":[\n')
        f.write(',\n'.join(json.dumps(p) for p in held['positions']) + '\n')
        if closed:
            f.write('],"genes":[\n' + ',\n'.join(json.dumps(gene) for gene in GENES) + '\n]}\n')
    return path


def gnomad_entry(variant_id: str, begin: int, end: int, kind: str, af: float | None, **extra: object) -> dict[str, Any]:
    """An entry of gnomAD's structural variants; give reciprocalOverlap and annotationOverlap for a match by overlap."""
    return {
        'chromosome': '1',
        'begin': begin,
        'end': end,
        'variantId': variant_id,
        'variantType': kind,
        'allAf': af,
        'allAc': 5,
        'allAn': 100,
        'allHc': 1,
        **extra,
    }


def clinvar_entry(
    accession: str,
    begin: int,
    end: int,
    kind: str,
    classification: str,
    status: str | None,
    evaluated: str,
    **extra: object,
) -> dict[str, Any]:
    """A ClinVar record matched to a structural record; give its overlaps for a match by overlap."""
    germline = {
        'classification': classification,
        'reviewStatus': status,
        'dateLastEvaluated': evaluated,
        'conditions': [{'traits': [{'name': {'value': 'An invented condition'}}]}],
    }
    return {
        'chromosome': '1',
        'begin': begin,
        'end': end,
        'accession': accession,
        'version': '1',
        'recordType': 'classified',
        'variantType': kind,
        'classifications': {'germlineClassification': germline} if status is not None else {},
        **extra,
    }


def overlaps(reciprocal: float, annotation: float | None = None) -> dict[str, float]:
    return {'reciprocalOverlap': reciprocal, 'annotationOverlap': reciprocal if annotation is None else annotation}


def annotate(
    folder: pathlib.Path,
    families: Mapping[str, Sequence[Sequence[Sequence[str]]]],
    annotated: Callable[[str, sv_sites.Site], dict[str, Any]] | None = None,
    sources: Sequence[str] = ALL,
) -> family_tables.SvAnnotation:
    """The annotation of a group's structural variants, from each family's records of each of its callers.

    `annotated` gives the position line of a site from the ID of one of its records and the site; by default every
    site has the plain line of position(), with one transcript.
    """
    folder.mkdir()
    sites = {
        family: [
            write_sites(folder / f'{family}.{index}.sites.vcf.gz', records) for index, records in enumerate(callers)
        ]
        for family, callers in families.items()
    }
    sv_sites.union(sites, folder / 'union.vcf', folder / 'map.parquet')
    with (folder / 'union.vcf').open() as plain, gzip.open(folder / 'union.vcf.gz', 'wt') as compressed:
        compressed.write(plain.read())
    union = sv_sites.read_union(folder / 'union.vcf.gz')
    record_of = {
        found.site_id: sv_id
        for paths in sites.values()
        for path in paths
        for sv_id, alleles in sv_sites.read_family(path)[1]
        for found in alleles
    }
    plain_line = annotated or (
        lambda _sv_id, found: position(
            found, {'transcripts': [transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant')]}
        )
    )
    document = {
        'header': header(),
        'positions': [plain_line(record_of[site_id], found) for site_id, found in union.items()],
    }
    write_json(folder / 'union.json.gz', document)
    info = sv_annotation.convert(
        folder / 'union.json.gz', folder / 'union.vcf.gz', run_config(*sources), BUNDLE, CATALOG, folder / 'tables'
    )
    (folder / 'tables' / 'annotation.json').write_text(info.model_dump_json())
    return family_tables.SvAnnotation(tables=folder / 'tables', site_map=folder / 'map.parquet')


def info(*sources: str) -> annotation_tables.AnnotationInfo:
    """What annotated a union whose run named the sources, as the tables of small variants record it."""
    return annotation_tables.annotation_info(header(), run_config(*sources), BUNDLE, CATALOG)


# The callers' records of structural_vcfs that have a row in sv, as a family's two callers give them.
CALLERS = (structural_vcfs.SV_RECORDS, structural_vcfs.CNV_RECORDS)
