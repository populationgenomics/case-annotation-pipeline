"""A family's tables, from an invented trio's VCF and an invented union's tables and gnomAD joint values."""

from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.resources
import json
import pathlib
import struct
from collections.abc import Mapping, Sequence
from typing import Any

import duckdb
import invented_sv_annotation
import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam
import pytest
import structural_vcfs

import case_annotation
from case_annotation import (
    annotation_tables,
    cli,
    family_tables,
    gnomad_joint,
    pedigree,
    repeat_tables,
    sv_annotation,
    sv_sites,
    sv_tables,
)

_PED = 'fam-1\tS1\t0\t0\t1\t1\nfam-1\tS2\t0\t0\t2\t1\nfam-1\tS3\tS1\tS2\t2\t2\n'
_CONTIGS = ('chr1', 'chr2', 'chr10', 'chrX', 'chrM', 'chr1_KI270706v1_random')

# Per record: chrom, pos, ref, alt, and each sample's (GT, AD, DP, GQ, PL), father, mother, child. In the VCF's order,
# which is natural chromosome order; a sort by name would put chr10 before chr2 and chrM before chrX.
type Call = tuple[tuple[int | None, ...], tuple[int, int], int, int, tuple[int, ...]]
_REF: Call = ((0, 0), (30, 0), 30, 60, (0, 60, 600))
_HET: Call = ((0, 1), (15, 15), 30, 60, (300, 0, 300))
_RECORDS: list[tuple[str, int, str, str, tuple[Call, Call, Call]]] = [
    # Low quality in the child.
    ('chr1', 5, 'C', 'T', (_REF, _REF, ((0, 1), (4, 2), 6, 10, (10, 0, 40)))),
    # Inherited from the father.
    ('chr2', 100, 'A', 'G', (_HET, _REF, _HET)),
    ('chr10', 7, 'G', 'GA', (_REF, _REF, _HET)),
    # Outside X's pseudo-autosomal regions: the father haploid, the mother and the daughter diploid.
    ('chrX', 5_000_000, 'T', 'C', (((0,), (30, 0), 30, 60, (0, 90)), _REF, _HET)),
    # Homoplasmic in the child.
    (
        'chrM',
        100,
        'A',
        'G',
        (((0,), (30, 0), 30, 60, (0, 90)), ((0,), (30, 0), 30, 60, (0, 90)), ((1,), (0, 30), 30, 60, (90, 0))),
    ),
]
_ALLELES = [annotation_tables.variant_id(c, p, r, a) for c, p, r, a, _ in _RECORDS]


def _write_vcf(path: pathlib.Path, records: Sequence[tuple[str, int, str, str, tuple[Call, ...]]]) -> None:
    header = pysam.VariantHeader()
    for contig in _CONTIGS:
        header.contigs.add(contig, length=200_000_000)
    header.add_meta('FILTER', items=[('ID', 'LowGQ'), ('Description', 'invented')])
    for key, number, kind in (('GT', '1', 'String'), ('AD', 'R', 'Integer'), ('DP', '1', 'Integer')):
        header.add_meta('FORMAT', items=[('ID', key), ('Number', number), ('Type', kind), ('Description', 'x')])
    for key, number, kind in (
        ('GQ', '1', 'Integer'),
        ('PL', 'G', 'Integer'),
        ('DN', '1', 'String'),
        ('DQ', '1', 'Float'),
    ):
        header.add_meta('FORMAT', items=[('ID', key), ('Number', number), ('Type', kind), ('Description', 'x')])
    for sample in ('S1', 'S2', 'S3'):
        header.add_sample(sample)
    with pysam.VariantFile(str(path), 'wz', header=header) as vcf:
        for chrom, pos, ref, alt, calls in records:
            record = vcf.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt), qual=50.0, filter='PASS')
            for sample, (gt, ad, dp, gq, pl) in zip(('S1', 'S2', 'S3'), calls, strict=True):
                record.samples[sample]['GT'] = gt
                record.samples[sample]['AD'] = ad
                record.samples[sample]['DP'] = dp
                record.samples[sample]['GQ'] = gq
                record.samples[sample]['PL'] = pl
            record.samples['S3']['DN'] = 'DeNovo'
            for parent in ('S1', 'S2'):
                record.samples[parent]['DN'] = '.'
            vcf.write(record)


def _info(*sources: str) -> annotation_tables.AnnotationInfo:
    return annotation_tables.AnnotationInfo(
        annotator='Invented Annotator 1.0',
        json_schema_version=6,
        assembly='GRCh38',
        catalog='209901',
        sources=tuple(sorted(sources)),
        source_versions=(),
        data_sources=(),
        omitted=(),
    )


def _float32(value: float) -> float:
    """The 32-bit float nearest to the value, as Python reads it from a 32-bit column."""
    return struct.unpack('f', struct.pack('f', value))[0]


# The dbSNP ids of the union's alleles, in the alleles' order: two on the first, none on the second, one on the others.
_RSIDS: tuple[list[str] | None, ...] = (['rs1', 'rs2'], None, ['rs3'], ['rs4'], ['rs5'], ['rs6'])
# The invented HGNC ids of the genes of the union's alleles, in the alleles' order: one gene, two, none, then one each.
_HGNC_IDS: tuple[list[int] | None, ...] = ([90001], [90001, 90002], None, [90002], [90001], [90003])


def _gene(symbol: str, hgnc_id: int | None) -> dict[str, Any]:
    """An entry of `genes` (annotation_tables.GENE) as the union's table holds it.

    A gene with an HGNC id has a MANE Select transcript, and a gene without one has none.
    """
    mane = None if hgnc_id is None else f'NM-{symbol}.1'
    return {
        'gene_symbol': symbol,
        'ensembl_gene_id': None,
        'hgnc_id': hgnc_id,
        'worst_consequence': 'intron_variant',
        'worst_impact': 'modifier',
        'worst_transcript_id': f'ENST-{symbol}.1',
        'mane_select_transcript_id': mane,
        'mane_select_consequence': mane and 'intron_variant',
        'mane_select_hgvsc': mane and f'{mane}:c.1+5A>G',
        'mane_select_hgvsp': None,
    }


# The genes of the union's alleles, in the alleles' order: those of the HGNC ids, between the second allele's two a
# gene without an id, and for the third allele no list.
_GENES: tuple[list[dict[str, Any]] | None, ...] = (
    [_gene('GENEA', 90001)],
    [_gene('GENEB', 90002), _gene('GENEX', None), _gene('GENEA', 90001)],
    None,
    [_gene('GENEB', 90002)],
    [_gene('GENEA', 90001)],
    [_gene('GENEC', 90003)],
)
# gnomAD's genomes and exomes values of every allele of the union.
_PER_TYPE = {
    f'gnomad_{kind}_{value}': held
    for kind in ('genomes', 'exomes')
    for value, held in (('af', 0.25), ('ac', 5), ('an', 20), ('hom', 1), ('failed_filter', False))
}


def _write_union(path: pathlib.Path, alleles: Sequence[str], info: annotation_tables.AnnotationInfo) -> None:
    """Union tables holding the alleles and an allele of another family, each in a gene, with one transcript."""
    (path / 'variants').mkdir(parents=True)
    (path / 'transcripts').mkdir()
    variants, transcripts = [], []
    for vid, rsids, hgnc_ids, genes in zip([*alleles, 'chr2-500-C-A'], _RSIDS, _HGNC_IDS, _GENES, strict=True):
        chrom, pos, ref, alt = vid.split('-')
        variants.append(
            {
                'chrom': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'variant_id': vid,
                'hgnc_ids': hgnc_ids,
                'genes': genes,
                'rsids': rsids,
                'clinvar_accession': 'VCV1',
                **_PER_TYPE,
            }
        )
        transcripts.append(
            {
                'chrom': chrom,
                'pos': int(pos),
                'variant_id': vid,
                'transcript_id': 'ENST01.1',
                'source': 'Ensembl',
                'hgnc_id': hgnc_ids[0] if hgnc_ids else None,
                'consequences': ['intron_variant'],
                'worst_consequence': 'intron_variant',
                'is_mane_select': False,
                'is_canonical': True,
            }
        )
    # Out of order across parts, as the workers write them.
    for index, (v, t) in enumerate(((variants[3:], transcripts[3:]), (variants[:3], transcripts[:3]))):
        pq.write_table(
            pa.Table.from_pylist(v, schema=annotation_tables.VARIANTS), path / 'variants' / f'p{index}.parquet'
        )
        pq.write_table(
            pa.Table.from_pylist(t, schema=annotation_tables.TRANSCRIPTS), path / 'transcripts' / f'p{index}.parquet'
        )
    (path / 'annotation.json').write_text(info.model_dump_json())


_JOINT_INFO = gnomad_joint.GnomadJointInfo(
    archive='invented.zip',
    bytes=1,
    sha256='0' * 64,
    echtvar_version='echtvar 0.0.0',
    fields=(
        gnomad_joint.ArchiveField(
            field='AC_joint', alias='AC', ftype='Integer', multiplier=1, missing_value=-1, missing_string='MISSING'
        ),
    ),
)
# gnomAD's joint values of the union's alleles as (AF, AC, homozygotes, FAF95): the first allele is one gnomAD lacks,
# the second one it counted 0 times, and the other family's allele has values of its own.
_JOINT_OTHER = (0.25, 400_000, 50_000, _float32(0.24))
_JOINT: dict[str, tuple[float, int, int, float] | None] = {
    _ALLELES[0]: None,
    _ALLELES[1]: (0.0, 0, 0, 0.0),
    _ALLELES[2]: (_float32(0.0000015), 2, 0, _float32(0.0000005)),
    _ALLELES[3]: (_float32(0.001), 1600, 2, _float32(0.0009)),
    _ALLELES[4]: None,
    'chr2-500-C-A': _JOINT_OTHER,
}


def _write_joint(path: pathlib.Path, values: Mapping[str, tuple[float, int, int, float] | None]) -> None:
    path.mkdir()
    rows = [
        {'variant_id': vid} | ({} if joint is None else dict(zip(gnomad_joint.COLUMNS, joint, strict=True)))
        for vid, joint in values.items()
    ]
    pq.write_table(pa.Table.from_pylist(rows, schema=gnomad_joint.TABLE), path / gnomad_joint.TABLE_NAME)
    (path / gnomad_joint.INFO_NAME).write_text(_JOINT_INFO.model_dump_json())


# The tables of a family's alleles, which every family has.
_ALLELE_TABLES = ('variants', 'transcripts')
# A family whose sheet row names no VCF of its larger variants.
_NO_STRUCTURAL = family_tables.Structural(sv_vcf=None, cnv_vcf=None, repeat_vcfs=())


def _inputs(
    tmp_path: pathlib.Path,
    records: Sequence[tuple[str, int, str, str, tuple[Call, ...]]] = _RECORDS,
    ped: str = _PED,
    info: annotation_tables.AnnotationInfo | None = None,
    joint: Mapping[str, tuple[float, int, int, float] | None] = _JOINT,
    proband: str | None = 'S3',
    structural: family_tables.Structural | None = None,
) -> dict[str, Any]:
    _write_vcf(tmp_path / 'family.vcf.gz', records)
    (tmp_path / 'family.ped').write_text(ped)
    _write_union(tmp_path / 'union', _ALLELES, info or _info('Ensembl', 'dbSNP', 'clinvar'))
    _write_joint(tmp_path / 'joint', joint)
    return {
        'family': 'fam-1',
        'vcf_path': tmp_path / 'family.vcf.gz',
        'ped_path': tmp_path / 'family.ped',
        'annotation': tmp_path / 'union',
        'joint': tmp_path / 'joint',
        'output': tmp_path / 'out',
        'proband': proband,
        'structural': structural or _NO_STRUCTURAL,
    }


def _rows(path: pathlib.Path, query: str = 'SELECT * FROM t') -> list[dict[str, Any]]:
    con = duckdb.connect()
    con.read_parquet(str(path)).create_view('t')
    relation = con.sql(query)
    return [dict(zip(relation.columns, r, strict=True)) for r in relation.fetchall()]


def test_a_familys_tables_hold_its_alleles_sorted_in_chromosome_order(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**_inputs(tmp_path))
    out = tmp_path / 'out'
    variants = _rows(out / 'variants.parquet')
    assert [v['variant_id'] for v in variants] == _ALLELES  # the other family's allele is left out
    assert (variants[0]['qual'], variants[0]['filters']) == (50.0, ['PASS'])
    assert [t['variant_id'] for t in _rows(out / 'transcripts.parquet')] == _ALLELES
    assert manifest.alleles == len(_ALLELES)
    assert {t.table: t.rows for t in manifest.tables} == {'variants': len(_ALLELES), 'transcripts': len(_ALLELES)}
    # The folder holds the two tables, the manifest and the notes on the tables, and each table as the manifest
    # records it.
    assert sorted(f.name for f in out.iterdir()) == [
        'manifest.json',
        'tables.md',
        'transcripts.parquet',
        'variants.parquet',
    ]
    assert {t.table: (t.path, t.bytes) for t in manifest.tables} == {
        table: (f'{table}.parquet', (out / f'{table}.parquet').stat().st_size) for table in ('variants', 'transcripts')
    }
    on_disk = json.loads((out / 'manifest.json').read_text())
    assert on_disk['family'] == 'fam-1'
    assert {m['sample_id']: m['role'] for m in on_disk['members']} == {'S1': 'father', 'S2': 'mother', 'S3': 'proband'}
    assert on_disk['schema_version'] == family_tables.SCHEMA_VERSION


def test_each_members_column_holds_its_own_call_with_its_likelihoods(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**_inputs(tmp_path))
    path = tmp_path / 'out' / 'variants.parquet'
    schema = pq.read_schema(path)
    # The members' columns follow the allele's, at the top level of the table and not inside a struct of their own.
    assert schema.names == [*family_tables.ALLELE.names, 'S1', 'S2', 'S3']
    assert all(schema.field(s).type.equals(family_tables.CALL, check_metadata=False) for s in ('S1', 'S2', 'S3'))
    rows = {c['variant_id']: c for c in _rows(path)}
    for (chrom, pos, ref, alt, calls), sample in (
        (record, sample) for record in _RECORDS for sample in ('S1', 'S2', 'S3')
    ):
        gt, (ad_ref, ad_alt), dp, gq, pl = calls[('S1', 'S2', 'S3').index(sample)]
        call = rows[annotation_tables.variant_id(chrom, pos, ref, alt)][sample]
        values = (call['dp'], call['gq'], call['ad_ref'], call['ad_alt'], call['pl'])
        assert values == (dp, gq, ad_ref, ad_alt, list(pl))
        assert call['gt'] == '/'.join(str(a) for a in gt)
    # pysam writes the haploid call's likelihoods padded with a missing value to the diploid calls' length.
    assert len(rows['chrX-5000000-T-C']['S1']['pl']) == 2
    assert rows['chr10-7-G-GA']['S3']['dn'] == 'DeNovo'
    assert rows['chr10-7-G-GA']['S1']['dn'] is None
    assert rows['chr2-100-A-G']['S1']['allele_balance'] == 0.5


def test_the_call_columns_follow_the_vcfs_sample_order(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    vcf = _text_vcf(tmp_path / 'call.vcf', _GT_AD, ('0/1:5,5', '0/0:9,0', '0/0:9,0'), 'GT:AD', ('S3', 'S1', 'S2'))
    family_tables.assemble(**{**inputs, 'vcf_path': vcf})
    assert pq.read_schema(tmp_path / 'out' / 'variants.parquet').names[-3:] == ['S3', 'S1', 'S2']


def test_the_variants_table_is_written_with_its_declared_schema(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**_inputs(tmp_path))
    written = pq.read_schema(tmp_path / 'out' / 'variants.parquet')
    assert written.equals(family_tables.variants_schema(('S1', 'S2', 'S3')), check_metadata=False)
    # The file itself says which columns always have a value: the allele, and each member's column and genotype.
    # Whether the member carries the allele is none of them, because a genotype can leave it open.
    parquet = pq.read_metadata(tmp_path / 'out' / 'variants.parquet').schema
    required = {
        parquet.column(i).path for i in range(len(parquet.names)) if parquet.column(i).max_definition_level == 0
    }
    members = {f'{sample}.gt' for sample in ('S1', 'S2', 'S3')}
    assert required == {'chrom', 'pos', 'ref', 'alt', 'variant_id'} | members


def test_the_declared_schema_leaves_out_the_columns_a_run_omits(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**_inputs(tmp_path, info=_info('Ensembl')))
    written = pq.read_schema(tmp_path / 'out' / 'variants.parquet')
    declared = family_tables.variants_schema(('S1', 'S2', 'S3'), manifest.omitted_columns)
    assert written.equals(declared, check_metadata=False)
    assert set(family_tables.ALLELE.names) - set(written.names) == set(manifest.omitted_columns)


def test_an_alleles_hgnc_ids_are_a_list_of_integers_and_a_lookup_by_one_finds_it(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**_inputs(tmp_path))
    out = tmp_path / 'out'
    assert pq.read_schema(out / 'variants.parquet').field('hgnc_ids').type == pa.list_(pa.int32())
    assert pq.read_schema(out / 'transcripts.parquet').field('hgnc_id').type == pa.int32()
    rows = {v['variant_id']: v['hgnc_ids'] for v in _rows(out / 'variants.parquet')}
    # An allele in no gene with an id has no list rather than an empty one.
    assert rows == dict(zip(_ALLELES, _HGNC_IDS[: len(_ALLELES)], strict=True))
    # An id is found on an allele that holds it beside another.
    found = _rows(out / 'variants.parquet', 'SELECT variant_id FROM t WHERE list_contains(hgnc_ids, 90002)')
    assert sorted(f['variant_id'] for f in found) == sorted([_ALLELES[1], _ALLELES[3]])
    by_transcript = _rows(out / 'transcripts.parquet', 'SELECT variant_id FROM t WHERE hgnc_id = 90002')
    assert by_transcript == [{'variant_id': _ALLELES[3]}]


def test_an_alleles_genes_are_a_list_of_structs_and_unnest_gives_a_row_for_each_gene(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**_inputs(tmp_path))
    path = tmp_path / 'out' / 'variants.parquet'
    schema = pq.read_schema(path)
    assert schema.field('genes').type == pa.list_(annotation_tables.GENE)
    # What is of one gene or one transcript of an allele is a field of an entry, and no column beside `genes`.
    of_a_gene = {field.name for field in annotation_tables.GENE} - {'worst_consequence', 'worst_impact', 'hgnc_id'}
    assert not of_a_gene & set(schema.names)
    # Each allele has its entries in the union's order, and an allele without genes no list rather than an empty one.
    rows = {v['variant_id']: v['genes'] for v in _rows(path)}
    assert rows == dict(zip(_ALLELES, _GENES[: len(_ALLELES)], strict=True))
    # A count by gene counts an allele under each of its genes, also under one without an HGNC id.
    counted = _rows(
        path,
        'SELECT gene.gene_symbol AS gene_symbol, count(*) AS alleles '
        'FROM (SELECT unnest(genes) AS gene FROM t) GROUP BY 1',
    )
    held = [g['gene_symbol'] for genes in _GENES[: len(_ALLELES)] for g in genes or ()]
    assert {c['gene_symbol']: c['alleles'] for c in counted} == {symbol: held.count(symbol) for symbol in set(held)}
    assert {'GENEA', 'GENEX'} <= set(held)


def test_a_call_without_likelihoods_has_none(tmp_path: pathlib.Path) -> None:
    formats = _GT_AD + '\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="x">'
    inputs = _inputs(tmp_path)
    vcf = _text_vcf(tmp_path / 'call.vcf', formats, ('0/0:9,0:.', './.:0,0:.', '0/1:5,5:50,0,50'), 'GT:AD:PL')
    family_tables.assemble(**{**inputs, 'vcf_path': vcf})
    rows = _rows(tmp_path / 'out' / 'variants.parquet', 'SELECT S1.pl AS s1, S2.pl AS s2, S3.pl AS s3 FROM t')
    assert rows == [{'s1': None, 's2': None, 's3': [50, 0, 50]}]


def test_each_allele_has_gnomads_joint_values_and_one_gnomad_lacks_has_none(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**_inputs(tmp_path))
    path = tmp_path / 'out' / 'variants.parquet'
    schema = pq.read_schema(path)
    assert [(c, str(schema.field(c).type)) for c in gnomad_joint.COLUMNS] == [
        ('gnomad_joint_af', 'float'),
        ('gnomad_joint_ac', 'int32'),
        ('gnomad_joint_hom', 'int32'),
        ('gnomad_joint_faf95_max', 'float'),
    ]
    rows = {v['variant_id']: tuple(v[c] for c in gnomad_joint.COLUMNS) for v in _rows(path)}
    # The family's alleles alone, each with its own values: none for an allele gnomAD lacks, 0 where gnomAD has 0. A
    # frequency is the union's 32-bit float, bit for bit.
    assert rows == {vid: (None, None, None, None) if _JOINT[vid] is None else _JOINT[vid] for vid in _ALLELES}
    assert rows[_ALLELES[3]][0] != 0.001
    assert manifest.gnomad_joint == _JOINT_INFO
    assert json.loads((tmp_path / 'out' / 'manifest.json').read_text())['gnomad_joint']['archive'] == 'invented.zip'
    assert json.loads(schema.metadata[b'gnomad_joint']) == _JOINT_INFO.model_dump(mode='json')
    # Illumina's separate genomes and exomes columns stay, their frequencies as 32-bit floats too, and the joint
    # columns are no dated source's.
    assert {str(schema.field(c).type) for c in ('gnomad_genomes_af', 'gnomad_exomes_af')} == {'float'}
    assert not set(gnomad_joint.COLUMNS) & set(manifest.omitted_columns)


def test_an_alleles_dbsnp_ids_are_a_list_and_an_allele_without_any_has_none(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**_inputs(tmp_path))
    path = tmp_path / 'out' / 'variants.parquet'
    kind = pq.read_schema(path).field('rsids').type
    assert pa.types.is_list(kind)
    assert kind.value_type == pa.string()
    rows = {v['variant_id']: v['rsids'] for v in _rows(path, 'SELECT variant_id, rsids FROM t')}
    assert rows == dict(zip(_ALLELES, _RSIDS[: len(_ALLELES)], strict=True))
    # An id is found on an allele that holds it beside another.
    assert _rows(path, "SELECT variant_id FROM t WHERE list_contains(rsids, 'rs2')") == [{'variant_id': _ALLELES[0]}]


def test_the_columns_of_a_dated_source_the_run_did_not_name_are_left_out(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**_inputs(tmp_path, info=_info('Ensembl')))
    columns = pq.read_schema(tmp_path / 'out' / 'variants.parquet').names
    assert 'rsids' not in columns
    assert not any(c.startswith('clinvar_') for c in columns)
    assert 'gnomad_genomes_af' in columns
    assert 'gnomad_joint_af' in columns
    assert 'rsids' in manifest.omitted_columns


def test_the_tables_carry_their_metadata_row_groups_and_bloom_filters(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(family_tables, 'ROW_GROUP_ROWS', 2)
    family_tables.assemble(**_inputs(tmp_path))
    path = tmp_path / 'out' / 'variants.parquet'
    metadata = pq.read_metadata(path)
    assert [metadata.row_group(i).num_rows for i in range(metadata.num_row_groups)] == [2, 2, 1]
    kv = pq.read_schema(path).metadata
    assert kv[b'table'] == b'variants'
    assert kv[b'schema_version'] == family_tables.SCHEMA_VERSION.encode()
    assert json.loads(kv[b'annotation'])['catalog'] == '209901'
    con = duckdb.connect()
    query = 'SELECT DISTINCT path_in_schema FROM parquet_metadata(?) WHERE bloom_filter_length > 0 ORDER BY 1'
    blooms = con.execute(query, [str(path)]).fetchall()
    # hgnc_ids is a list, and its Bloom filter is on the list's values.
    assert [b[0] for b in blooms] == ['hgnc_ids, list, element', 'variant_id']
    blooms = con.execute(query, [str(tmp_path / 'out' / 'transcripts.parquet')]).fetchall()
    assert [b[0] for b in blooms] == ['gene_id', 'gene_symbol', 'hgnc_id', 'transcript_id', 'variant_id']


def test_the_bloom_filter_of_a_column_the_run_omits_is_left_out_with_it(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bloom = {**family_tables._BLOOM, 'variants': {'variant_id': 10, 'rsids.list.element': 10}}
    monkeypatch.setattr(family_tables, '_BLOOM', bloom)
    family_tables.assemble(**_inputs(tmp_path, info=_info('Ensembl')))
    con = duckdb.connect()
    query = 'SELECT DISTINCT path_in_schema FROM parquet_metadata(?) WHERE bloom_filter_length > 0 ORDER BY 1'
    assert con.execute(query, [str(tmp_path / 'out' / 'variants.parquet')]).fetchall() == [('variant_id',)]


def _value_encodings(path: pathlib.Path) -> dict[str, str]:
    """How each column of a Parquet file holds its values, by its path: DICTIONARY, or the encoding without one.

    A column chunk lists RLE for its levels beside the encoding of its values, and a chunk with a dictionary lists the
    dictionary's own encodings. A boolean column, which never has a dictionary, is left out.
    """
    metadata = pq.read_metadata(path)
    found: dict[str, set[str]] = {}
    for group in range(metadata.num_row_groups):
        for index in range(metadata.num_columns):
            chunk = metadata.row_group(group).column(index)
            if chunk.physical_type == 'BOOLEAN':
                continue
            values = {'DICTIONARY'} if chunk.has_dictionary_page else set(chunk.encodings) - {'RLE'}
            found.setdefault(chunk.path_in_schema, set()).update(values)
    assert all(len(values) == 1 for values in found.values()), found
    return {column: values.pop() for column, values in found.items()}


def test_each_column_named_for_an_encoding_is_written_with_it_and_every_other_column_with_a_dictionary(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(family_tables, 'ROW_GROUP_ROWS', 2)
    family_tables.assemble(**_inputs(tmp_path))
    out = tmp_path / 'out'
    for table in _ALLELE_TABLES:
        declared = family_tables._ENCODINGS[table]
        written = _value_encodings(out / f'{table}.parquet')
        # The fixture names every dated source, so each declared column is in its table, a nested one included.
        assert declared.items() <= written.items()
        assert {encoding for column, encoding in written.items() if column not in declared} == {'DICTIONARY'}
        # pos is sorted within a chromosome in every table, and stored as differences.
        assert written['pos'] == 'DELTA_BINARY_PACKED'
        # A column without a dictionary is compressed like the others.
        metadata = pq.read_metadata(out / f'{table}.parquet')
        assert {metadata.row_group(0).column(i).compression for i in range(metadata.num_columns)} == {'ZSTD'}
    assert {'rsids.list.element', 'genes.list.element.mane_select_hgvsc'} <= family_tables._ENCODINGS['variants'].keys()


def test_a_column_named_for_an_encoding_that_the_table_lacks_is_refused(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    encodings = {**family_tables._ENCODINGS, 'transcripts': {'gnomad_joint_af': 'PLAIN'}}
    monkeypatch.setattr(family_tables, '_ENCODINGS', encodings)
    with pytest.raises(ValueError, match=r"no columns \['gnomad_joint_af'\] to write without a dictionary"):
        family_tables.assemble(**_inputs(tmp_path))


def test_the_encodings_of_an_omitted_column_are_left_out() -> None:
    assert 'rsids.list.element' in family_tables._encodings('variants', ())
    without = family_tables._encodings('variants', ('rsids',))
    assert 'rsids.list.element' not in without
    assert without.keys() == family_tables._ENCODINGS['variants'].keys() - {'rsids.list.element'}


def test_the_position_and_the_variant_id_still_pick_out_a_row_group_in_every_table(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(family_tables, 'ROW_GROUP_ROWS', 2)
    family_tables.assemble(**_inputs(tmp_path))
    con = duckdb.connect()
    probe = 'SELECT bloom_filter_excludes FROM parquet_bloom_probe(?, ?, ?) ORDER BY row_group_id'
    for table in _ALLELE_TABLES:
        path = tmp_path / 'out' / f'{table}.parquet'
        rows = _rows(path, 'SELECT pos, variant_id FROM t')
        metadata = pq.read_metadata(path)
        groups = range(metadata.num_row_groups)
        assert len(groups) == 3
        index = {metadata.schema.column(i).path: i for i in range(metadata.num_columns)}
        for group in groups:
            held = rows[group * 2 : group * 2 + 2]
            # A row group's statistics give the range of its values, which a filter on either column skips it by.
            for column in ('pos', 'variant_id'):
                chunk = metadata.row_group(group).column(index[column])
                values = [row[column] for row in held]
                assert (chunk.statistics.min, chunk.statistics.max) == (min(values), max(values))
                assert chunk.has_column_index
                assert chunk.has_offset_index
            # The Bloom filter of variant_id rules out every row group but the one that holds the id.
            excluded = con.execute(probe, [str(path), 'variant_id', held[0]['variant_id']]).fetchall()
            assert [e[0] for e in excluded] == [other != group for other in groups]
        assert all(e[0] for e in con.execute(probe, [str(path), 'variant_id', 'chr1-1-A-C']).fetchall())
        # A lookup by either column returns the row.
        wanted = (rows[-1]['pos'], rows[-1]['variant_id'])
        lookup = 'SELECT pos, variant_id FROM read_parquet(?) WHERE '
        assert con.execute(lookup + 'variant_id = ?', [str(path), wanted[1]]).fetchall() == [wanted]
        assert wanted in con.execute(lookup + 'pos = ?', [str(path), wanted[0]]).fetchall()


def test_the_command_refuses_another_schema_version(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    with pytest.raises(SystemExit, match=r'not 0\.0\.0'):
        cli.main([*_tables_command(inputs, version='0.0.0'), '--proband', 'S3'])
    cli.main([*_tables_command(inputs), '--proband', 'S3'])
    assert (inputs['output'] / 'manifest.json').exists()


def test_the_command_refuses_a_workflow_of_another_release_of_the_pipeline(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    with pytest.raises(SystemExit, match=r'and the workflow of 0\.0\.0'):
        cli.main([*_tables_command(inputs, pipeline='0.0.0'), '--proband', 'S3'])
    # Nothing is written for it, and the workflow's own release writes a manifest that gives it.
    assert not inputs['output'].exists()
    cli.main([*_tables_command(inputs), '--proband', 'S3'])
    written = json.loads((inputs['output'] / 'manifest.json').read_text())
    assert written['pipeline_version'] == family_tables.PIPELINE_VERSION


def _tables_command(
    inputs: Mapping[str, Any],
    output: pathlib.Path | None = None,
    version: str = family_tables.SCHEMA_VERSION,
    pipeline: str = family_tables.PIPELINE_VERSION,
) -> list[str]:
    """The family-tables command of the inputs, without its proband."""
    args = ['family-tables', '--family', 'fam-1', '--ped', str(inputs['ped_path']), '--vcf', str(inputs['vcf_path'])]
    args += ['--annotation', str(inputs['annotation']), '--gnomad-joint', str(inputs['joint'])]
    args += ['--output', str(output or inputs['output'])]
    return [*args, '--schema-version', version, '--pipeline-version', pipeline]


def _check_command(inputs: Mapping[str, Any], family: str = 'fam-1') -> list[str]:
    """The family-check command of the inputs, without its proband."""
    return ['family-check', '--family', family, '--ped', str(inputs['ped_path']), '--vcf', str(inputs['vcf_path'])]


def test_the_family_check_refuses_a_ped_file_of_another_family(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    cli.main([*_check_command(inputs), '--proband', 'S3'])
    with pytest.raises(ValueError, match='of family'):
        cli.main([*_check_command(inputs, family='fam-2'), '--proband', 'S3'])


def _roles(output: pathlib.Path) -> dict[str, str | None]:
    return {m['sample_id']: m['role'] for m in json.loads((output / 'manifest.json').read_text())['members']}


def test_the_commands_take_the_proband_as_stated_and_none_where_none_is_designated(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    cli.main([*_check_command(inputs), '--no-proband'])
    # The mother is named, although the PED file has the daughter as the only affected member.
    cli.main([*_check_command(inputs), '--proband', 'S2'])
    cli.main([*_tables_command(inputs), '--proband', 'S2'])
    assert _roles(inputs['output']) == {'S1': None, 'S2': 'proband', 'S3': None}
    cli.main([*_tables_command(inputs, output=tmp_path / 'none'), '--no-proband'])
    assert _roles(tmp_path / 'none') == {'S1': None, 'S2': None, 'S3': None}


@pytest.mark.parametrize('stated', [[], ['--proband', 'S3', '--no-proband']], ids=['neither', 'both'])
def test_the_commands_refuse_to_run_unless_the_proband_is_stated_one_way(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], stated: list[str]
) -> None:
    inputs = _inputs(tmp_path)
    for command in (_check_command(inputs), _tables_command(inputs)):
        with pytest.raises(SystemExit):
            cli.main([*command, *stated])
        assert '--proband' in capsys.readouterr().err
    assert not inputs['output'].exists()


# A brother of the daughter who is in the PED file and was not sequenced.
_PED_WITH_BROTHER = _PED + 'fam-1\tS4\tS1\tS2\t1\t2\n'
PROBANDS_REFUSED = {
    'a proband the PED file lacks': (_PED, 'S9', "the proband 'S9' of 'fam-1' is not an individual of its PED file"),
    'a proband the VCF lacks': (_PED_WITH_BROTHER, 'S4', "the proband 'S4' of 'fam-1' is not a sample of its VCF"),
    'an empty id': (_PED, '', "the proband '' of 'fam-1' is not an individual of its PED file"),
}


@pytest.mark.parametrize('case', PROBANDS_REFUSED)
def test_a_proband_that_is_not_a_sequenced_member_of_the_family_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    ped, proband, message = PROBANDS_REFUSED[case]
    inputs = _inputs(tmp_path, ped=ped, proband=proband)
    with pytest.raises(ValueError, match=message):
        family_tables.check_inputs('fam-1', inputs['ped_path'], inputs['vcf_path'], proband, _NO_STRUCTURAL)
    with pytest.raises(ValueError, match=message):
        cli.main([*_check_command(inputs), '--proband', proband])
    with pytest.raises(ValueError, match=message):
        family_tables.assemble(**inputs)
    assert not inputs['output'].exists()


_NOBODY: tuple[Call, Call, Call] = (_REF, _REF, _REF)
REFUSED = {
    'an allele the union lacks': (
        {'records': [*_RECORDS, ('chr2', 900, 'T', 'A', (_HET, _REF, _HET))]},
        "alleles the union's annotation lacks",
    ),
    "an allele the union's gnomAD joint values lack": (
        {'joint': {vid: _JOINT[vid] for vid in _ALLELES[1:]}},
        "alleles the union's gnomAD joint values lack, such as \\['chr1-5-C-T'\\]",
    ),
    'an allele off the primary chromosomes': (
        {'records': [*_RECORDS, ('chr1_KI270706v1_random', 9, 'T', 'A', (_HET, _REF, _HET))]},
        'primary chromosome',
    ),
    'an allele nobody carries': ({'records': [*_RECORDS, ('chr2', 100, 'A', 'G', _NOBODY)]}, 'no member carries'),
    'a repeated allele': ({'records': [*_RECORDS, _RECORDS[0]]}, 'repeats alleles'),
    'a PED file of another family': ({'ped': _PED.replace('fam-1', 'fam-2')}, 'of family'),
    'a VCF sample the PED file lacks': (
        {'ped': 'fam-1\tS2\t0\t0\t2\t1\nfam-1\tS3\t0\tS2\t2\t2\n'},
        'lacks the VCF samples',
    ),
}


@pytest.mark.parametrize('case', REFUSED)
def test_an_input_outside_the_contract_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    changes, message = REFUSED[case]
    with pytest.raises(ValueError, match=message):
        family_tables.assemble(**_inputs(tmp_path, **changes))


def test_a_record_of_two_alleles_is_refused(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    header = pysam.VariantFile(str(inputs['vcf_path'])).header
    with pysam.VariantFile(str(tmp_path / 'two.vcf.gz'), 'wz', header=header) as vcf:
        record = vcf.new_record(contig='chr1', start=4, alleles=('C', 'T', 'G'))
        record.samples['S3']['GT'] = (0, 1)
        vcf.write(record)
    with pytest.raises(ValueError, match='one alternate allele'):
        family_tables.assemble(**{**inputs, 'vcf_path': tmp_path / 'two.vcf.gz'})


def _text_vcf(
    path: pathlib.Path,
    formats: str,
    sample_values: Sequence[str],
    format_keys: str,
    samples: Sequence[str] = ('S1', 'S2', 'S3'),
) -> pathlib.Path:
    """A one-record VCF written as text, for values a pysam writer would refuse or convert."""
    lines = [
        '##fileformat=VCFv4.2',
        '##contig=<ID=chr1,length=1000>',
        formats,
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *samples]),
        '\t'.join(['chr1', '5', '.', 'C', 'T', '50', 'PASS', '.', format_keys, *sample_values]),
    ]
    path.write_text('\n'.join(lines) + '\n')
    return path


_GT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">'
_GT_AD = _GT + '\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="x">'
CALLS_REFUSED = {
    'allele depths for more than the two alleles': (_GT_AD, ('0/0:9,0', '0/0:9,0', '0/1:5,5,1'), 'GT:AD', 'AD'),
    'allele depths that are not integers': (
        _GT + '\n##FORMAT=<ID=AD,Number=R,Type=Float,Description="x">',
        ('0/0:9,0', '0/0:9,0', '0/1:5.5,5'),
        'GT:AD',
        'not a depth of the reference and one of the allele',
    ),
    'one depth for the two alleles': (
        _GT + '\n##FORMAT=<ID=AD,Number=1,Type=Integer,Description="x">',
        ('0/0:9', '0/0:9', '0/1:5'),
        'GT:AD',
        'not a depth of the reference and one of the allele',
    ),
    'a depth that is not an integer': (
        _GT_AD + '\n##FORMAT=<ID=DP,Number=1,Type=Float,Description="x">',
        ('0/0:9,0:9.5', '0/0:9,0:9', '0/1:5,5:10'),
        'GT:AD:DP',
        'not an integer',
    ),
    'a quality that is not a number': (
        _GT_AD + '\n##FORMAT=<ID=DQ,Number=1,Type=String,Description="x">',
        ('0/0:9,0:.', '0/0:9,0:.', '0/1:5,5:high'),
        'GT:AD:DQ',
        'not a number',
    ),
    'a sample filter that is not a string': (
        _GT_AD + '\n##FORMAT=<ID=FT,Number=1,Type=Integer,Description="x">',
        ('0/0:9,0:1', '0/0:9,0:1', '0/1:5,5:1'),
        'GT:AD:FT',
        'not a string',
    ),
    "likelihoods of a diploid call that are a haploid one's": (
        _GT_AD + '\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="x">',
        ('0/0:9,0:0,90', '0/0:9,0:0,90,900', '0/1:5,5:50,0,50'),
        'GT:AD:PL',
        'a call of 2 alleles has 3',
    ),
    "likelihoods of a haploid call that are a diploid one's": (
        _GT_AD + '\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="x">',
        ('0:9,0:0,90,900', '0/0:9,0:0,90,900', '0/1:5,5:50,0,50'),
        'GT:AD:PL',
        'a call of 1 alleles has 2',
    ),
    'a likelihood missing beside others': (
        _GT_AD + '\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="x">',
        ('0/0:9,0:0,.,900', '0/0:9,0:0,90,900', '0/1:5,5:50,0,50'),
        'GT:AD:PL',
        'each present',
    ),
    'likelihoods that are not a list': (
        _GT_AD + '\n##FORMAT=<ID=PL,Number=1,Type=Integer,Description="x">',
        ('0/0:9,0:7', '0/0:9,0:7', '0/1:5,5:7'),
        'GT:AD:PL',
        'not a list of likelihoods',
    ),
    'likelihoods that are not integers': (
        _GT_AD + '\n##FORMAT=<ID=PL,Number=G,Type=Float,Description="x">',
        ('0/0:9,0:0,90.5,900', '0/0:9,0:0,90,900', '0/1:5,5:50,0,50'),
        'GT:AD:PL',
        'integer likelihoods',
    ),
    'no genotypes': (
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="x">',
        ('9,0', '9,0', '5,5'),
        'AD',
        'no GT',
    ),
}


@pytest.mark.parametrize('case', CALLS_REFUSED)
def test_a_call_outside_the_contract_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    formats, values, keys, message = CALLS_REFUSED[case]
    inputs = _inputs(tmp_path)
    vcf = _text_vcf(tmp_path / 'call.vcf', formats, values, keys)
    with pytest.raises(ValueError, match=message):
        family_tables.assemble(**{**inputs, 'vcf_path': vcf})


def test_a_call_with_an_allele_missing_says_only_what_its_called_alleles_decide(tmp_path: pathlib.Path) -> None:
    # ./1 holds a copy of the allele whatever the missing allele is. 0/. and ./. say neither that the member carries the
    # allele nor that it does not, since a missing allele may be the allele. None of the three says on how many copies
    # the allele lies, so none has a zygosity.
    inputs = _inputs(tmp_path)
    vcf = _text_vcf(tmp_path / 'call.vcf', _GT_AD, ('./1:10,10', '0/.:20,0', './.:0,0'), 'GT:AD')
    family_tables.assemble(**{**inputs, 'vcf_path': vcf})
    (row,) = _rows(tmp_path / 'out' / 'variants.parquet')
    assert [(row[s]['gt'], row[s]['carries'], row[s]['zygosity']) for s in ('S1', 'S2', 'S3')] == [
        ('./1', True, None),
        ('0/.', None, None),
        ('./.', None, None),
    ]


def test_a_calls_filters_are_a_list_of_names_as_the_records_are(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    formats = _GT_AD + '\n##FORMAT=<ID=FT,Number=1,Type=String,Description="x">'
    vcf = _text_vcf(
        tmp_path / 'call.vcf', formats, ('0/1:10,10:PASS', '0/0:20,0:.', '0/1:5,5:LowDepth;LowGQ'), 'GT:AD:FT'
    )
    family_tables.assemble(**{**inputs, 'vcf_path': vcf})
    path = tmp_path / 'out' / 'variants.parquet'
    schema = pq.read_schema(path)
    assert schema.field('S1').type.field('ft').type == schema.field('filters').type == pa.list_(pa.string())
    (row,) = _rows(path)
    # A call that passed holds the one name PASS, as the record does, and a call without filters has no list.
    assert row['filters'] == ['PASS']
    assert [row[s]['ft'] for s in ('S1', 'S2', 'S3')] == [['PASS'], None, ['LowDepth', 'LowGQ']]
    # One condition finds a filter among several.
    found = _rows(
        path, "SELECT variant_id FROM t WHERE list_contains(S3.ft, 'LowGQ') AND NOT list_contains(S1.ft, 'LowGQ')"
    )
    assert len(found) == 1


_UNCALLED: Call = ((None, None), (0, 0), 0, 0, (0, 0, 0))


def test_a_filter_for_a_parent_who_does_not_carry_an_allele_leaves_out_a_parent_without_a_call(
    tmp_path: pathlib.Path,
) -> None:
    # The father has no call of the insertion that the daughter carries. At chr1:5 both parents are called for the
    # reference.
    records = [
        (chrom, pos, ref, alt, (_UNCALLED, *calls[1:]) if (chrom, pos) == ('chr10', 7) else calls)
        for chrom, pos, ref, alt, calls in _RECORDS
    ]
    family_tables.assemble(**_inputs(tmp_path, records=records))
    path = tmp_path / 'out' / 'variants.parquet'
    (insertion,) = _rows(path, "SELECT S1 AS father FROM t WHERE variant_id = 'chr10-7-G-GA'")
    father = insertion['father']
    assert (father['gt'], father['carries'], father['zygosity']) == ('./.', None, None)

    def in_neither_parent(father: str) -> set[str]:
        query = f'SELECT variant_id FROM t WHERE S3.carries AND {father} AND NOT S2.carries'  # noqa: S608 (the conditions of this test)
        return {row['variant_id'] for row in _rows(path, query)}

    called_without = in_neither_parent('NOT S1.carries')
    assert 'chr1-5-C-T' in called_without
    assert 'chr10-7-G-GA' not in called_without
    # A parent who is called and lacks the allele is the one that hom_ref names.
    assert called_without == in_neither_parent("S1.zygosity = 'hom_ref'")
    # A reader who also wants the alleles that the father is not shown to carry asks for that.
    assert in_neither_parent('S1.carries IS NOT TRUE') - called_without == {'chr10-7-G-GA'}


def test_a_member_of_unknown_phenotype_is_null_in_the_manifest(tmp_path: pathlib.Path) -> None:
    # The parents' phenotypes are unknown, written as 0 and as -9, and the daughter has an unaffected sister who was
    # not sequenced.
    ped = 'fam-1\tS1\t0\t0\t1\t0\nfam-1\tS2\t0\t0\t2\t-9\nfam-1\tS3\tS1\tS2\t2\t2\nfam-1\tS4\tS1\tS2\t2\t1\n'
    manifest = family_tables.assemble(**_inputs(tmp_path, ped=ped))
    written = json.loads((tmp_path / 'out' / 'manifest.json').read_text())
    assert {m['sample_id']: m['affected'] for m in written['members']} == {
        'S1': None,
        'S2': None,
        'S3': True,
        'S4': False,
    }
    # The roles are the sheet's and the PED file's parents, whatever the phenotypes.
    assert {m['sample_id']: m['role'] for m in written['members']} == {
        'S1': 'father',
        'S2': 'mother',
        'S3': 'proband',
        'S4': None,
    }
    assert family_tables.FamilyManifest.model_validate(written) == manifest
    # A manifest states each member's phenotype, also where it is unknown: the key is not left out for a null.
    without = [{key: value for key, value in m.items() if key != 'affected'} for m in written['members']]
    with pytest.raises(pydantic.ValidationError, match='affected\n  Field required'):
        family_tables.FamilyManifest.model_validate({**written, 'members': without})


def test_a_member_carries_an_allele_that_its_genotype_holds_a_copy_of(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**_inputs(tmp_path))
    path = tmp_path / 'out' / 'variants.parquet'
    kind = pq.read_schema(path).field('S1').type
    assert [kind.field(i).name for i in range(3)] == ['gt', 'carries', 'zygosity']
    assert kind.field('carries').type == pa.bool_()
    assert kind.field('carries').nullable
    rows = {v['variant_id']: v for v in _rows(path)}
    for (chrom, pos, ref, alt, calls), sample in (
        (record, sample) for record in _RECORDS for sample in ('S1', 'S2', 'S3')
    ):
        gt = calls[('S1', 'S2', 'S3').index(sample)][0]
        assert rows[annotation_tables.variant_id(chrom, pos, ref, alt)][sample]['carries'] == (1 in gt)
    # A haploid call of the allele carries it, as the child's on chrM, and one of the reference does not.
    assert [rows['chrM-100-A-G'][s]['carries'] for s in ('S1', 'S2', 'S3')] == [False, False, True]
    carriers = _rows(path, 'SELECT variant_id FROM t WHERE S1.carries')
    assert [c['variant_id'] for c in carriers] == ['chr2-100-A-G']


# None of these families brings a VCF of its larger variants: a column of sv or of repeats is refused all the same.
COLUMN_CLASHES = {
    'a sample named as an allele column, in another case': (('S1', 'S2', 'POS'), "column name 'pos'"),
    'a sample named as an annotation column': (('S1', 'S2', 'revel'), "column name 'revel'"),
    "a sample named as a column of the family's record": (('S1', 'S2', 'Filters'), "column name 'filters'"),
    'two samples whose names differ only in case': (('S1', 's1', 'S3'), "column name 'S1'"),
    'a sample named as a column of sv, in another case': (
        ('S1', 'S2', 'SV_ID'),
        "column name 'sv_id' in the sv table",
    ),
    'a sample named as a column of repeats': (
        ('S1', 'S2', 'repeat_unit'),
        "column name 'repeat_unit' in the repeats table",
    ),
}


@pytest.mark.parametrize('case', COLUMN_CLASHES)
def test_a_sample_that_cannot_name_its_own_column_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    samples, message = COLUMN_CLASHES[case]
    inputs = _inputs(tmp_path)
    ped = ''.join(f'fam-1\t{s}\t0\t0\t1\t1\n' for s in samples)
    (tmp_path / 'clash.ped').write_text(ped)
    vcf = _text_vcf(tmp_path / 'call.vcf', _GT_AD, ('0/1:5,5', '0/0:9,0', '0/0:9,0'), 'GT:AD', samples)
    with pytest.raises(ValueError, match=message):
        family_tables.check_inputs('fam-1', tmp_path / 'clash.ped', vcf, None, _NO_STRUCTURAL)
    with pytest.raises(ValueError, match=message):
        family_tables.assemble(**{**inputs, 'vcf_path': vcf, 'ped_path': tmp_path / 'clash.ped', 'proband': None})


def test_no_sample_can_take_the_name_of_a_column_of_a_table_that_holds_the_members() -> None:
    # Also a column that a run without its source leaves out: a name is refused whatever the run's sources.
    assert {'chrom', 'hgnc_ids', 'rsids', 'clinvar_records', 'gnomad_joint_af', 'qual'} <= set(
        family_tables.ALLELE.names
    )
    assert family_tables.ROW_COLUMNS.keys() == {'variants', 'sv', 'repeats'}
    first: dict[str, str] = {}
    for table, columns in family_tables.ROW_COLUMNS.items():
        for column in columns.names:
            # A column that several tables have is reported for the first of them.
            held = first.setdefault(column, table)
            for sample in (column, column.upper()):
                with pytest.raises(ValueError, match=f'column name {column!r} in the {held} table'):
                    family_tables.check_column_names(['S1', sample])
    assert {first[column] for column in ('caller', 'small_variant_id', 'repeat_unit')} == {'sv', 'repeats'}
    family_tables.check_column_names(['S1', 'S2', 'S3'])


def test_the_reserved_names_are_the_columns_before_the_members_in_each_table(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path)})
    for table, columns in family_tables.ROW_COLUMNS.items():
        written = pq.read_schema(tmp_path / 'out' / f'{table}.parquet')
        assert written.names == [*columns.names, 'S1', 'S2', 'S3'], table
        assert all(pa.types.is_struct(written.field(sample).type) for sample in ('S1', 'S2', 'S3')), table


def test_a_call_of_a_larger_variant_has_the_fields_it_shares_with_a_small_variants_under_one_name_and_type() -> None:
    small = {field.name: field for field in family_tables.CALL}
    larger = {field.name: field for field in sv_tables.CALL}
    shared = [name for name in small if name in larger]
    assert shared == ['gt', 'carries', 'zygosity', 'gq', 'pl', 'ft', 'dn', 'dq']
    assert all(small[name].equals(larger[name]) for name in shared)
    # They stand in the same order in both, the fields of one kind of variant between them.
    assert [name for name in larger if name in small] == shared
    assert list(larger)[:3] == list(small)[:3] == ['gt', 'carries', 'zygosity']
    # A repeat's genotype is written as a call's is; the table has a row for every repeat, so nothing is carried.
    repeat = {field.name: field for field in repeat_tables.CALL}
    assert repeat['gt'].equals(small['gt'])
    assert not {'carries', 'zygosity'} & repeat.keys()


# Two affected children of unaffected parents: the sheet's proband decides the roles, whoever is affected.
_SIBLINGS = _PED_WITH_BROTHER
_EVERYONE = {'S1', 'S2', 'S3', 'S4'}
ROLES = {
    'the child of a trio': (_PED, {'S1', 'S2', 'S3'}, 'S3', {'S1': 'father', 'S2': 'mother', 'S3': 'proband'}),
    'a parent not sequenced': (_PED, {'S2', 'S3'}, 'S3', {'S2': 'mother', 'S3': 'proband'}),
    'the first of two affected children': (
        _SIBLINGS,
        _EVERYONE,
        'S3',
        {'S1': 'father', 'S2': 'mother', 'S3': 'proband'},
    ),
    'the second of two affected children': (
        _SIBLINGS,
        _EVERYONE,
        'S4',
        {'S1': 'father', 'S2': 'mother', 'S4': 'proband'},
    ),
    'an unaffected founder': (_SIBLINGS, _EVERYONE, 'S1', {'S1': 'proband'}),
    'no proband, with an affected child': (_PED, {'S1', 'S2', 'S3'}, None, {}),
    'no proband, with two affected children': (_SIBLINGS, _EVERYONE, None, {}),
}


@pytest.mark.parametrize('case', ROLES)
def test_the_roles_are_the_named_proband_and_its_sequenced_parents(tmp_path: pathlib.Path, case: str) -> None:
    text, sequenced, proband, roles = ROLES[case]
    (tmp_path / 'family.ped').write_text(text)
    assert family_tables.member_roles(pedigree.read(tmp_path / 'family.ped'), sequenced, proband) == roles


# The tables of the larger variants, from the invented VCFs of structural_vcfs.

_SV_META = (*structural_vcfs.SV_INFO, *structural_vcfs.SV_FORMAT)
_CNV_META = (*structural_vcfs.CNV_INFO, *structural_vcfs.CNV_FORMAT)
# The sv table of a family whose run named every dated source: the record's columns, the annotation's, the members'.
_SV_SCHEMA = family_tables.sv_schema(('S1', 'S2', 'S3'))
# The rows of sv in the table's order: by chromosome, position and end, both callers' records together.
_SV_ROWS = (
    'inv:DEL:1',
    'inv:LOSS:1',
    'inv:DEL:2',
    'inv:BND:1:0',
    'inv:GAIN:1',
    'inv:INS:1',
    'inv:BND:1:1',
    'inv:CNV:1',
    'inv:LOSS:2',
    'inv:DEL:3',
    'inv:GAIN:2',
    'inv:GAIN:3',
    'inv:DUP:TANDEM:1',
    'inv:LOSS:3',
    'inv:GAIN:4',
)


# The sources of the run that annotated the family's alleles in _inputs, which annotates its structural variants too.
_SOURCES = ('Ensembl', 'dbSNP', 'clinvar')
# Another family of the group, with records of its own.
_OTHER_FAMILY = {'fam-2': invented_sv_annotation.CALLERS}


def _larger(
    tmp_path: pathlib.Path,
    sv: Sequence[structural_vcfs.Record] | None = structural_vcfs.SV_RECORDS,
    cnv: Sequence[structural_vcfs.Record] | None = structural_vcfs.CNV_RECORDS,
    repeats: bool = True,
    alleles: Mapping[str, tuple[str, str, str, str]] | None = structural_vcfs.SV_ALLELES,
    samples: Sequence[str] = structural_vcfs.SAMPLES,
    annotated: Any = None,  # noqa: ANN401 (what invented_sv_annotation.annotate takes)
    sources: Sequence[str] = _SOURCES,
    group: Mapping[str, Sequence[Sequence[structural_vcfs.Record]]] | None = None,
) -> dict[str, Any]:
    """The arguments of assemble that name a family's larger variants; None leaves a VCF out.

    A family with a caller's VCF also gets the annotation of its group's structural variants, which `annotated` gives
    each site's position line of and whose run named `sources`. The group is the family alone, or the families of
    `group`, whose records then stand for this family's too.
    """
    folder = tmp_path / 'larger'
    folder.mkdir()
    structural = family_tables.Structural(
        sv_vcf=None if sv is None else structural_vcfs.write(folder / 'sv.vcf', _SV_META, sv, samples),
        cnv_vcf=None if cnv is None else structural_vcfs.write(folder / 'cnv.vcf', _CNV_META, cnv, samples),
        repeat_vcfs=structural_vcfs.write_repeats(folder) if repeats else (),
    )
    sv_alleles = None if alleles is None else structural_vcfs.write_sv_alleles(folder / 'sv-alleles.vcf', alleles)
    callers = [records for records in (sv, cnv) if records is not None]
    sv_annotated = None
    if callers:
        families = {'fam-1': callers} if group is None else group
        sv_annotated = invented_sv_annotation.annotate(folder / 'annotation', families, annotated, sources)
    return {'structural': structural, 'sv_alleles': sv_alleles, 'sv_annotated': sv_annotated}


def test_a_family_with_the_vcfs_of_its_larger_variants_gets_their_tables(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path)})
    out = tmp_path / 'out'
    assert {t.table: t.rows for t in manifest.tables} == {
        'variants': len(_ALLELES),
        'transcripts': len(_ALLELES),
        'sv': len(_SV_ROWS),
        # One invented transcript for each allele: the segment that one member loses and another gains has two.
        'sv_transcripts': len(_SV_ROWS) + 1,
        'repeats': len(structural_vcfs.REPEAT_IDS),
    }
    assert manifest.sv_callers == ('sv', 'cnv')
    assert manifest.alleles == len(_ALLELES)
    assert all((out / t.path).stat().st_size == t.bytes for t in manifest.tables)
    # Both callers' records in one order, and nothing off the primary chromosomes.
    assert [r['sv_id'] for r in _rows(out / 'sv.parquet')] == list(_SV_ROWS)
    assert [r['repeat_id'] for r in _rows(out / 'repeats.parquet')] == list(structural_vcfs.REPEAT_IDS)
    callers = {r['sv_id']: r['caller'] for r in _rows(out / 'sv.parquet')}
    assert {sv_id for sv_id, caller in callers.items() if caller == 'cnv'} == set(structural_vcfs.CNV_CALLS)
    # The tables are written with their declared schemas; sv ends with a call column per member.
    assert pq.read_schema(out / 'sv.parquet').equals(_SV_SCHEMA, check_metadata=False)
    assert pq.read_schema(out / 'sv_transcripts.parquet').equals(
        sv_tables.transcripts_schema(sv_annotation.TRANSCRIPT_COLUMNS), check_metadata=False
    )
    assert pq.read_schema(out / 'repeats.parquet').equals(
        repeat_tables.schema(('S1', 'S2', 'S3')), check_metadata=False
    )
    # The alleles' tables are as without the larger variants.
    assert [v['variant_id'] for v in _rows(out / 'variants.parquet')] == _ALLELES
    assert json.loads((out / 'manifest.json').read_text())['sv_callers'] == ['sv', 'cnv']


def test_a_family_without_such_vcfs_has_no_table_of_them(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**_inputs(tmp_path))
    assert [t.table for t in manifest.tables] == list(_ALLELE_TABLES)
    assert manifest.sv_callers == ()
    assert (manifest.sv_annotation, manifest.sv_omitted_columns) == (None, ())
    assert sorted(f.name for f in (tmp_path / 'out').iterdir()) == [
        'manifest.json',
        'tables.md',
        'transcripts.parquet',
        'variants.parquet',
    ]


@pytest.mark.parametrize('larger', [False, True])
def test_a_familys_folder_holds_the_packages_notes_on_the_tables_as_its_manifest_names_them(
    tmp_path: pathlib.Path, larger: bool
) -> None:
    manifest = family_tables.assemble(**{**_inputs(tmp_path), **(_larger(tmp_path) if larger else {})})
    packaged = importlib.resources.files(case_annotation).joinpath(family_tables.NOTES).read_bytes()
    written = (tmp_path / 'out' / manifest.notes.path).read_bytes()
    # Byte for byte the file of the package that wrote the tables, so it describes the tables beside it.
    assert written == packaged
    assert manifest.notes == family_tables.NotesFile(
        path='tables.md', bytes=len(packaged), sha256=hashlib.sha256(packaged).hexdigest()
    )
    recorded = json.loads((tmp_path / 'out' / 'manifest.json').read_text())['notes']
    assert recorded == {'path': 'tables.md', 'bytes': len(written), 'sha256': hashlib.sha256(written).hexdigest()}


def test_the_manifest_gives_the_release_of_the_pipeline_that_wrote_the_folder(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**_inputs(tmp_path))
    written = json.loads((tmp_path / 'out' / 'manifest.json').read_text())
    # The version of the installed package, which tests/test_package.py holds to the workflow's and the contract's.
    installed = importlib.metadata.version('case-annotation-tools')
    assert written['pipeline_version'] == manifest.pipeline_version == installed
    without = {name: value for name, value in written.items() if name != 'pipeline_version'}
    with pytest.raises(pydantic.ValidationError, match='pipeline_version\n  Field required'):
        family_tables.FamilyManifest.model_validate(without)


def test_the_manifest_of_a_family_without_such_vcfs_states_that_it_has_none(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**_inputs(tmp_path))
    written = json.loads((tmp_path / 'out' / 'manifest.json').read_text())
    assert written.keys() == family_tables.FamilyManifest.model_fields.keys()
    assert (written['sv_callers'], written['sv_annotation'], written['sv_omitted_columns']) == ([], None, [])
    assert family_tables.FamilyManifest.model_validate(written) == manifest


@pytest.mark.parametrize('key', family_tables.FamilyManifest.model_fields)
def test_a_manifest_without_one_of_its_keys_is_refused(tmp_path: pathlib.Path, key: str) -> None:
    # Also a key whose value says that the family has nothing for it, as each key of the larger variants does here.
    family_tables.assemble(**_inputs(tmp_path))
    written = json.loads((tmp_path / 'out' / 'manifest.json').read_text())
    without = {name: value for name, value in written.items() if name != key}
    with pytest.raises(pydantic.ValidationError, match=f'{key}\n  Field required'):
        family_tables.FamilyManifest.model_validate(without)


def test_a_manifest_with_a_key_the_model_does_not_know_is_refused(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**_inputs(tmp_path))
    written = json.loads((tmp_path / 'out' / 'manifest.json').read_text())
    with pytest.raises(pydantic.ValidationError, match='sv_caller\n  Extra inputs are not permitted'):
        family_tables.FamilyManifest.model_validate({**written, 'sv_caller': ['sv']})


_SV_TABLES = ('sv', 'sv_transcripts')
TABLES_BY_INPUT = {
    'a structural-variant VCF alone': ({'cnv': None, 'repeats': False}, _SV_TABLES, ('sv',)),
    'a copy-number VCF alone': ({'sv': None, 'alleles': None, 'repeats': False}, _SV_TABLES, ('cnv',)),
    'repeat VCFs alone': ({'sv': None, 'alleles': None, 'cnv': None}, ('repeats',), ()),
}


@pytest.mark.parametrize('case', TABLES_BY_INPUT)
def test_a_family_gets_the_tables_of_the_vcfs_it_brings(tmp_path: pathlib.Path, case: str) -> None:
    given, tables, callers = TABLES_BY_INPUT[case]
    manifest = family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path, **given)})
    assert [t.table for t in manifest.tables] == [*_ALLELE_TABLES, *tables]
    assert manifest.sv_callers == callers
    if callers:
        assert {r['caller'] for r in _rows(tmp_path / 'out' / 'sv.parquet')} == set(callers)


def test_a_structural_variant_that_is_an_allele_of_variants_names_it(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path)})
    out = tmp_path / 'out'
    linked = {r['sv_id']: r['small_variant_id'] for r in _rows(out / 'sv.parquet') if r['small_variant_id']}
    # The insertion is the family's allele. The deletion's normalised allele is none of them, and a symbolic allele,
    # a breakend and a copy-number segment have no sequence to compare.
    assert linked == {'inv:INS:1': 'chr10-7-G-GA'}
    joined = 'SELECT sv_id, v.rsids FROM t JOIN read_parquet(?) v ON t.small_variant_id = v.variant_id'
    con = duckdb.connect()
    con.read_parquet(str(out / 'sv.parquet')).create_view('t')
    assert con.execute(joined, [str(out / 'variants.parquet')]).fetchall() == [('inv:INS:1', ['rs3'])]


def test_an_allele_that_the_familys_variants_lack_links_nothing(tmp_path: pathlib.Path) -> None:
    # chr2-500-C-A is another family's allele, which the union's tables hold and this family's variants do not.
    alleles = {**structural_vcfs.SV_ALLELES, 'inv:DEL:2': ('chr2', '500', 'C', 'A')}
    family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path, alleles=alleles)})
    rows = {r['sv_id']: r['small_variant_id'] for r in _rows(tmp_path / 'out' / 'sv.parquet')}
    assert rows['inv:DEL:2'] is None


def test_the_member_columns_of_the_larger_variants_follow_the_familys_vcf(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    vcf = _text_vcf(tmp_path / 'call.vcf', _GT_AD, ('0/1:5,5', '0/0:9,0', '0/0:9,0'), 'GT:AD', ('S3', 'S1', 'S2'))
    family_tables.assemble(**{**inputs, **_larger(tmp_path), 'vcf_path': vcf})
    out = tmp_path / 'out'
    sv = pq.read_schema(out / 'sv.parquet')
    assert sv.names == [*family_tables.ROW_COLUMNS['sv'].names, 'S3', 'S1', 'S2']
    assert all(sv.field(s).type.equals(sv_tables.CALL, check_metadata=False) for s in ('S1', 'S2', 'S3'))
    repeats = pq.read_schema(out / 'repeats.parquet')
    assert repeats.names == [*repeat_tables.LOCUS.names, 'S3', 'S1', 'S2']
    assert all(repeats.field(s).type.equals(repeat_tables.CALL, check_metadata=False) for s in ('S1', 'S2', 'S3'))
    # Each column holds its own member's call, whatever the order of the files' samples.
    deletion = _rows(
        out / 'sv.parquet', "SELECT S1.gt AS s1, S2.gt AS s2, S3.gt AS s3 FROM t WHERE sv_id = 'inv:DEL:1'"
    )
    assert deletion == [{'s1': '0/1', 's2': '0/0', 's3': '1/1'}]
    sizes = _rows(out / 'repeats.parquet', "SELECT S1.repcn AS s1, S3.repcn AS s3 FROM t WHERE repeat_id = 'RPTA'")
    assert sizes == [{'s1': [10, 12], 's3': [12, 45]}]


def test_the_calls_of_both_callers_say_who_carries_a_record(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path)})
    rows = {r['sv_id']: r for r in _rows(tmp_path / 'out' / 'sv.parquet')}
    for sv_id, expected in structural_vcfs.CNV_CALLS.items():
        assert tuple((rows[sv_id][s]['carries'], rows[sv_id][s]['zygosity']) for s in ('S1', 'S2', 'S3')) == expected
    # A deletion that the father carries on his one X, read from either caller's record of it.
    on_x = 'SELECT sv_id FROM t WHERE chrom = ? AND S1.zygosity = ? AND S1.carries ORDER BY 1'
    con = duckdb.connect()
    con.read_parquet(str(tmp_path / 'out' / 'sv.parquet')).create_view('t')
    assert con.execute(on_x, ['chrX', 'hom_alt']).fetchall() == [('inv:DEL:3',), ('inv:LOSS:2',)]
    # A gain is carried with no zygosity, so the two fields differ there: carries finds it and zygosity does not.
    gains = 'SELECT sv_id FROM t WHERE S1.carries AND S1.zygosity IS NULL ORDER BY 1'
    assert {g[0] for g in con.sql(gains).fetchall()} == {
        sv_id for sv_id, calls in structural_vcfs.CNV_CALLS.items() if calls[0] == (True, None)
    }
    assert con.sql(gains).fetchall()
    # A member without a genotype in a structural-variant record, or without a copy number in a segment, is not said to
    # lack the record: a filter for the records that the daughter does not carry returns none of them.
    uncalled = {r[0] for r in con.sql('SELECT sv_id FROM t WHERE S3.carries IS NULL').fetchall()}
    assert uncalled == {'inv:DEL:2', 'inv:DUP:TANDEM:1', 'inv:LOSS:3', 'inv:GAIN:4'}
    not_carried = {r[0] for r in con.sql('SELECT sv_id FROM t WHERE NOT S3.carries').fetchall()}
    assert not_carried
    assert not not_carried & uncalled
    assert not con.sql('SELECT sv_id FROM t WHERE S3.carries IS NULL AND S3.zygosity IS NOT NULL').fetchall()


def test_the_tables_of_the_larger_variants_carry_their_metadata_filters_and_encodings(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(family_tables, 'ROW_GROUP_ROWS', 2)
    # With matches, so that the annotation's columns of sv hold values to encode.
    family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path, annotated=_matched)})
    con = duckdb.connect()
    blooms = 'SELECT DISTINCT path_in_schema FROM parquet_metadata(?) WHERE bloom_filter_length > 0 ORDER BY 1'
    probe = 'SELECT bloom_filter_excludes FROM parquet_bloom_probe(?, ?, ?) ORDER BY row_group_id'
    transcript_keys = ['gene_id', 'gene_symbol', 'hgnc_id', 'sv_id', 'transcript_id']
    for table, keys, rows in (
        ('sv', ['sv_id'], 15),
        ('sv_transcripts', transcript_keys, 16),
        ('repeats', ['repeat_id'], 3),
    ):
        key = 'repeat_id' if table == 'repeats' else 'sv_id'
        path = tmp_path / 'out' / f'{table}.parquet'
        metadata = pq.read_metadata(path)
        groups = [metadata.row_group(i).num_rows for i in range(metadata.num_row_groups)]
        assert (sum(groups), set(groups[:-1])) == (rows, {2})
        # The two tables that hold the annotation of the group's structural variants carry its provenance beside their
        # own version and family. Nothing annotated the repeats.
        kv = pq.read_schema(path).metadata
        provenance = {b'annotation'} if table in ('sv', 'sv_transcripts') else set()
        assert set(kv) == {b'schema_version', b'family', b'table'} | provenance
        if provenance:
            recorded = annotation_tables.AnnotationInfo.model_validate_json(kv[b'annotation'])
            assert recorded == invented_sv_annotation.info(*_SOURCES)
        assert (kv[b'table'], kv[b'family']) == (table.encode(), b'fam-1')
        assert kv[b'schema_version'] == family_tables.SCHEMA_VERSION.encode()
        assert [b[0] for b in con.execute(blooms, [str(path)]).fetchall()] == keys
        # The key's Bloom filter rules out every row group for an id the table lacks, and not the one that holds an id.
        held = _rows(path)[0][key]
        assert all(e[0] for e in con.execute(probe, [str(path), key, 'inv:none']).fetchall())
        assert con.execute(probe, [str(path), key, held]).fetchall()[0] == (False,)
        index = {metadata.schema.column(i).path: i for i in range(metadata.num_columns)}
        chunks = [metadata.row_group(0).column(i) for i in range(metadata.num_columns)]
        assert {chunk.compression for chunk in chunks} == {'ZSTD'}
        assert chunks[index['pos']].has_column_index
        assert chunks[index['pos']].has_offset_index
        # Each column named for an encoding is written with it, in every row group that holds a value of it.
        for column, encoding in family_tables._ENCODINGS[table].items():
            written = {
                e
                for group in range(metadata.num_row_groups)
                for e in metadata.row_group(group).column(index[column]).encodings
            }
            assert not metadata.row_group(0).column(index[column]).has_dictionary_page
            assert written - {'RLE'} == {encoding}, (table, column)
        assert {'pos', 'end_pos'} <= family_tables._ENCODINGS[table].keys()


def test_a_callers_vcf_without_a_record_gives_an_empty_table(tmp_path: pathlib.Path) -> None:
    # A copy-number VCF of segments at the reference's copy number alone, and a structural-variant VCF of no record.
    # The group's annotation is that of another family's records: this family has no site in it.
    reference_only = [r for r in structural_vcfs.CNV_RECORDS if r[4] == '.']
    larger = _larger(tmp_path, sv=[], cnv=reference_only, alleles={}, repeats=False, group=_OTHER_FAMILY)
    manifest = family_tables.assemble(**{**_inputs(tmp_path), **larger})
    assert {t.table: t.rows for t in manifest.tables if t.table in _SV_TABLES} == dict.fromkeys(_SV_TABLES, 0)
    assert pq.read_schema(tmp_path / 'out' / 'sv.parquet').names == _SV_SCHEMA.names
    assert manifest.sv_callers == ('sv', 'cnv')


_OFF_PRIMARY = tuple(('1', *record[1:]) for record in structural_vcfs.SV_RECORDS[:1])
_UNPREFIXED_REPEAT = ('1', *structural_vcfs.repeat_records('S1')[0][1:])


def _unprefixed(tmp_path: pathlib.Path) -> dict[str, Any]:
    """Larger variants on contigs named 1, 2, ...: none is on a primary chromosome."""
    larger = _larger(tmp_path, cnv=None, repeats=False, alleles={})
    lines = larger['structural'].sv_vcf.read_text().replace('##contig=<ID=chr1,', '##contig=<ID=1,').splitlines()
    kept = [line for line in lines if line.startswith('#')] + ['\t'.join(_OFF_PRIMARY[0])]
    larger['structural'].sv_vcf.write_text('\n'.join(kept) + '\n')
    return larger


def _unprefixed_repeats(tmp_path: pathlib.Path) -> dict[str, Any]:
    larger = _larger(tmp_path, sv=None, cnv=None, alleles=None)
    for path in larger['structural'].repeat_vcfs:
        lines = path.read_text().replace('##contig=<ID=chr1,', '##contig=<ID=1,').splitlines()
        sample = lines[-1].split('\t')[-1]
        kept = [line for line in lines if line.startswith('#')] + ['\t'.join((*_UNPREFIXED_REPEAT[:-1], sample))]
        path.write_text('\n'.join(kept) + '\n')
    return larger


def _shared_id(tmp_path: pathlib.Path) -> dict[str, Any]:
    """A copy-number segment under the ID of a structural variant."""
    segment = structural_vcfs.cnv_record('chr1', 20_000, 'inv:DEL:1', '<DEL>', (('0/1', '1'),) * 3)
    # The union of a group refuses such a family too, so the annotation here is that of another family's records.
    return _larger(tmp_path, cnv=[segment], repeats=False, group=_OTHER_FAMILY)


def _alleles_text(tmp_path: pathlib.Path, lines: Sequence[str]) -> dict[str, Any]:
    """The larger variants with the normalised sequences replaced by the lines given."""
    larger = _larger(tmp_path, cnv=None, repeats=False)
    header = [line for line in larger['sv_alleles'].read_text().splitlines() if line.startswith('#')]
    larger['sv_alleles'].write_text('\n'.join([*header, *lines]) + '\n')
    return larger


_INSERTION = 'chr10\t7\tinv:INS:1\tG\tGA\t.\t.\t.'
_DELETION = 'chr2\t298\tinv:DEL:2\tTGAC\tT\t.\t.\t.'
LARGER_REFUSED: dict[str, tuple[Any, str]] = {
    'a structural-variant VCF without its normalised sequences': (
        lambda tmp_path: _larger(tmp_path, alleles=None),
        'are given together',
    ),
    'normalised sequences without a structural-variant VCF': (
        lambda tmp_path: _larger(tmp_path, sv=None),
        'are given together',
    ),
    'normalised sequences that lack a record': (
        lambda tmp_path: _alleles_text(tmp_path, [_INSERTION]),
        r"lacks sequences of the family's structural variants, such as \['inv:DEL:2'\]",
    ),
    'normalised sequences of a record that is no sequence': (
        lambda tmp_path: _alleles_text(tmp_path, [_INSERTION, _DELETION, _DELETION.replace('inv:DEL:2', 'inv:DEL:1')]),
        r"holds records that are no sequence of the family's structural variants: \['inv:DEL:1'\]",
    ),
    "normalised sequences of another family's record": (
        lambda tmp_path: _alleles_text(tmp_path, [_INSERTION, _DELETION, _DELETION.replace('inv:DEL:2', 'inv:DEL:8')]),
        r"holds records that are no sequence of the family's structural variants: \['inv:DEL:8'\]",
    ),
    'normalised sequences with a record twice': (
        lambda tmp_path: _alleles_text(tmp_path, [_INSERTION, _DELETION, _DELETION]),
        'has no ID, or the ID of another record',
    ),
    'normalised sequences with a record without its ID': (
        lambda tmp_path: _alleles_text(tmp_path, [_INSERTION, _DELETION.replace('inv:DEL:2', '.')]),
        'has no ID, or the ID of another record',
    ),
    'normalised sequences with a record of two alleles': (
        lambda tmp_path: _alleles_text(tmp_path, [_INSERTION, _DELETION.replace('\tT\t', '\tT,TG\t')]),
        'one alternate allele',
    ),
    'two records of one ID': (_shared_id, r"give several records one ID, such as \['inv:DEL:1'\]"),
    'structural variants on contigs without the chr prefix': (_unprefixed, 'sv.vcf: no record on the primary'),
    'repeats on contigs without the chr prefix': (_unprefixed_repeats, 'repeats.vcf: no record on the primary'),
    "a structural-variant VCF of another family's samples": (
        lambda tmp_path: _larger(tmp_path, samples=('S1', 'S2', 'S9')),
        r"sv.vcf: its samples are \['S1', 'S2', 'S9'\]",
    ),
    'a copy-number VCF without one member': (
        lambda tmp_path: _larger(tmp_path, sv=None, alleles=None, cnv=[], samples=('S1', 'S2'), group=_OTHER_FAMILY),
        r"cnv.vcf: its samples are \['S1', 'S2'\]",
    ),
}


@pytest.mark.parametrize('case', LARGER_REFUSED)
def test_larger_variants_outside_the_contract_are_refused(tmp_path: pathlib.Path, case: str) -> None:
    larger, message = LARGER_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        family_tables.assemble(**{**_inputs(tmp_path), **larger(tmp_path)})


def test_repeat_vcfs_of_some_members_are_refused_before_anything_is_written(tmp_path: pathlib.Path) -> None:
    inputs = _inputs(tmp_path)
    larger = _larger(tmp_path)
    some = family_tables.Structural(sv_vcf=None, cnv_vcf=None, repeat_vcfs=larger['structural'].repeat_vcfs[:2])
    message = r"no repeat VCF of the samples \['S3'\]"
    with pytest.raises(ValueError, match=message):
        family_tables.check_inputs('fam-1', inputs['ped_path'], inputs['vcf_path'], 'S3', some)
    with pytest.raises(ValueError, match=message):
        family_tables.assemble(**{**inputs, 'structural': some})
    assert not inputs['output'].exists()


def _larger_options(larger: Mapping[str, Any]) -> list[str]:
    structural: family_tables.Structural = larger['structural']
    options = ['--sv-vcf', str(structural.sv_vcf), '--cnv-vcf', str(structural.cnv_vcf)]
    return [*options, *(part for vcf in structural.repeat_vcfs for part in ('--repeat-vcf', str(vcf)))]


def _annotation_options(larger: Mapping[str, Any]) -> list[str]:
    annotated: family_tables.SvAnnotation = larger['sv_annotated']
    return ['--sv-annotation', str(annotated.tables), '--sv-map', str(annotated.site_map)]


def test_the_commands_take_the_vcfs_of_the_larger_variants(tmp_path: pathlib.Path) -> None:
    inputs, larger = _inputs(tmp_path), _larger(tmp_path)
    options = _larger_options(larger)
    cli.main([*_check_command(inputs), '--proband', 'S3', *options])
    # Without one member's repeat VCF the check refuses the family, and so does the command that writes its tables.
    with pytest.raises(ValueError, match='no repeat VCF of the samples'):
        cli.main([*_check_command(inputs), '--proband', 'S3', *options[:-2]])
    annotation = _annotation_options(larger)
    with pytest.raises(ValueError, match='are given together'):
        cli.main([*_tables_command(inputs), '--proband', 'S3', *options, *annotation])
    # The annotation of the group's structural variants comes as its tables and its map, both or neither.
    command = [*_tables_command(inputs), '--proband', 'S3', *options, '--sv-alleles', str(larger['sv_alleles'])]
    with pytest.raises(SystemExit, match='--sv-annotation and --sv-map are given together'):
        cli.main([*command, *annotation[:2]])
    with pytest.raises(ValueError, match="need the annotation of the group's sites"):
        cli.main(command)
    assert not (inputs['output'] / 'manifest.json').exists()
    cli.main(
        [*_tables_command(inputs, output=tmp_path / 'written'), *command[len(_tables_command(inputs)) :], *annotation]
    )
    manifest = json.loads((tmp_path / 'written' / 'manifest.json').read_text())
    assert [t['table'] for t in manifest['tables']] == [*_ALLELE_TABLES, *_SV_TABLES, 'repeats']
    assert manifest['schema_version'] == family_tables.SCHEMA_VERSION


# The annotation of the group's structural variants in a family's tables.

_GENE = invented_sv_annotation.transcript('ENST01.1', 'Ensembl', 'GENEA', 'ENSG01', 'intron_variant')


def _matched(sv_id: str, site: sv_sites.Site) -> dict[str, Any]:
    """An invented annotation: matches for the symbolic deletion and for each allele of the two-allele segment."""
    entry, record, overlaps = (
        invented_sv_annotation.gnomad_entry,
        invented_sv_annotation.clinvar_entry,
        invented_sv_annotation.overlaps,
    )
    single = 'criteria provided, single submitter'
    sources: dict[str, Any] = {}
    transcripts = [_GENE]
    if sv_id == 'inv:DEL:1':
        sources['gnomad'] = [entry('inv-gnomad-1', 1001, 5000, 'deletion', 0.25, **overlaps(1.0))]
        sources['clinvar'] = [
            record('VCV9', 1001, 4000, 'deletion', 'Pathogenic', single, '2099-01-01', **overlaps(0.75))
        ]
    if sv_id == 'inv:CNV:1':
        loss = site.alt == '<DEL>'
        kind = 'copy_number_loss' if loss else 'copy_number_gain'
        sources['gnomad'] = [
            entry(f'inv-gnomad-{kind}', 5001, 15000, kind.replace('copy_number_', ''), 0.5, **overlaps(1.0))
        ]
        sources['clinvar'] = [
            record('VCV1' if loss else 'VCV2', 5001, 15000, kind, 'Benign', single, '2099-01-01', **overlaps(1.0))
        ]
        transcripts = [
            _GENE | {'transcript': 'ENST09.1'} if loss else _GENE | {'transcript': 'ENST08.1', 'hgnc': 'GENEC'}
        ]
    for source in sources.values():
        for held in source:
            held['variantType'] = {'loss': 'copy_number_loss', 'gain': 'copy_number_gain'}.get(
                held['variantType'], held['variantType']
            )
    return invented_sv_annotation.position(site, {'transcripts': transcripts}, **sources)


def test_each_record_has_the_annotation_of_its_site(tmp_path: pathlib.Path) -> None:
    manifest = family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path, annotated=_matched)})
    rows = {r['sv_id']: r for r in _rows(tmp_path / 'out' / 'sv.parquet')}
    deletion = rows['inv:DEL:1']
    assert (deletion['annotated_as'], deletion['annotated_pos'], deletion['annotated_end_pos']) == (
        'structural',
        1000,
        5000,
    )
    assert (deletion['gnomad_sv_id'], deletion['gnomad_sv_af'], deletion['gnomad_sv_overlap']) == (
        'inv-gnomad-1',
        0.25,
        1.0,
    )
    assert (deletion['clinvar_accession'], deletion['clinvar_classification'], deletion['clinvar_overlap']) == (
        'VCV9',
        'Pathogenic',
        0.75,
    )
    assert [(r['accession'], r['pos'], r['end_pos']) for r in deletion['clinvar_records']] == [('VCV9', 1000, 4000)]
    # The genes of the record's transcripts by their HGNC ids, on the record and on each transcript's row.
    gene = invented_sv_annotation.HGNC['GENEA']
    assert deletion['hgnc_ids'] == [gene]
    transcripts = [r for r in _rows(tmp_path / 'out' / 'sv_transcripts.parquet') if r['sv_id'] == 'inv:DEL:1']
    assert [(r['gene_symbol'], r['hgnc_id']) for r in transcripts] == [('GENEA', gene)]
    # A record without a match has the columns, empty; the short insertion is annotated as a small variant.
    assert rows['inv:DEL:3']['gnomad_sv_id'] is None
    assert rows['inv:DEL:3']['gnomad_sv_matches'] is None
    assert rows['inv:INS:1']['annotated_as'] == 'small_variant'
    assert {r['annotated_as'] for r in rows.values()} == {'structural', 'small_variant'}
    assert manifest.sv_annotation == invented_sv_annotation.info(*_SOURCES)
    assert manifest.sv_omitted_columns == ()


def test_a_record_of_two_alleles_keeps_the_matches_of_both_and_no_summary(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path, annotated=_matched)})
    out = tmp_path / 'out'
    segment = next(r for r in _rows(out / 'sv.parquet') if r['sv_id'] == 'inv:CNV:1')
    assert segment['alt'] == '<DEL>,<DUP>'
    # The loss's matches, then the gain's, each with its type. A best match would be that of one of the two alleles.
    assert [(m['id'], m['variant_type']) for m in segment['gnomad_sv_matches']] == [
        ('inv-gnomad-copy_number_loss', 'copy_number_loss'),
        ('inv-gnomad-copy_number_gain', 'copy_number_gain'),
    ]
    assert [(r['accession'], r['variant_type']) for r in segment['clinvar_records']] == [
        ('VCV1', 'copy_number_loss'),
        ('VCV2', 'copy_number_gain'),
    ]
    summaries = [name for name in sv_annotation.COLUMNS.names if name.startswith(('gnomad_', 'onekg_', 'clinvar_'))]
    assert {name for name in summaries if segment[name] is not None} == set(sv_annotation.LISTS)
    assert (segment['annotated_as'], segment['annotated_pos'], segment['annotated_end_pos']) == (
        'structural',
        5000,
        15000,
    )
    # The genes of both alleles' transcripts, and each allele's own transcripts under its own type.
    assert segment['hgnc_ids'] == [invented_sv_annotation.HGNC[gene] for gene in ('GENEA', 'GENEC')]
    transcripts = [
        (r['variant_type'], r['transcript_id'], r['chrom'], r['pos'], r['end_pos'])
        for r in _rows(out / 'sv_transcripts.parquet')
        if r['sv_id'] == 'inv:CNV:1'
    ]
    assert transcripts == [
        ('copy_number_gain', 'ENST08.1', 'chr10', 5000, 15000),
        ('copy_number_loss', 'ENST09.1', 'chr10', 5000, 15000),
    ]


def test_the_transcripts_of_the_records_stand_in_the_records_order(tmp_path: pathlib.Path) -> None:
    family_tables.assemble(**{**_inputs(tmp_path), **_larger(tmp_path)})
    rows = _rows(tmp_path / 'out' / 'sv_transcripts.parquet')
    seen = list(dict.fromkeys(r['sv_id'] for r in rows))
    assert seen == list(_SV_ROWS)
    keys = [(r['sv_id'], r['variant_type'], r['transcript_id']) for r in rows]
    assert len(set(keys)) == len(keys)
    places = {r['sv_id']: (r['chrom'], r['pos'], r['end_pos']) for r in _rows(tmp_path / 'out' / 'sv.parquet')}
    assert all((r['chrom'], r['pos'], r['end_pos']) == places[r['sv_id']] for r in rows)


def test_a_family_reads_its_own_records_of_a_union_that_families_share(tmp_path: pathlib.Path) -> None:
    # The other family holds the same records and one of its own; a site that both hold is annotated once.
    own = (*structural_vcfs.SV_RECORDS, ('chr2', '900', 'inv:DEL:7', 'T', '<DEL>', '.', '.', 'SVTYPE=DEL;END=990'))
    group = {'fam-1': invented_sv_annotation.CALLERS, 'fam-2': (own, structural_vcfs.CNV_RECORDS)}
    larger = _larger(tmp_path, group=group)
    union = sv_sites.read_union(tmp_path / 'larger' / 'annotation' / 'union.vcf.gz')
    assert len(union) == len(_SV_ROWS) + 2  # each allele of the family's records, and the other family's deletion
    manifest = family_tables.assemble(**{**_inputs(tmp_path), **larger})
    assert {t.table: t.rows for t in manifest.tables if t.table in _SV_TABLES} == {
        'sv': len(_SV_ROWS),
        'sv_transcripts': len(_SV_ROWS) + 1,
    }
    assert 'inv:DEL:7' not in {r['sv_id'] for r in _rows(tmp_path / 'out' / 'sv_transcripts.parquet')}


def test_the_sv_columns_of_a_dated_source_the_run_did_not_name_are_left_out(tmp_path: pathlib.Path) -> None:
    sources = ('Ensembl', 'dbSNP')
    inputs = _inputs(tmp_path, info=_info(*sources))
    manifest = family_tables.assemble(**{**inputs, **_larger(tmp_path, sources=sources)})
    names = pq.read_schema(tmp_path / 'out' / 'sv.parquet').names
    clinvar = [name for name in sv_annotation.COLUMNS.names if name.startswith('clinvar_')]
    assert clinvar
    assert names == [name for name in _SV_SCHEMA.names if name not in clinvar]
    assert manifest.sv_omitted_columns == tuple(clinvar)
    assert 'gnomad_sv_matches' in names


def test_the_tables_of_the_alleles_are_the_same_bytes_with_and_without_larger_variants(tmp_path: pathlib.Path) -> None:
    (tmp_path / 'with').mkdir()
    (tmp_path / 'without').mkdir()
    with_larger = family_tables.assemble(
        **{**_inputs(tmp_path / 'with'), **_larger(tmp_path / 'with', annotated=_matched)}
    )
    without = family_tables.assemble(**_inputs(tmp_path / 'without'))
    hashes = [{t.table: t.sha256 for t in m.tables if t.table in _ALLELE_TABLES} for m in (with_larger, without)]
    assert hashes[0] == hashes[1]
    assert len(hashes[0]) == len(_ALLELE_TABLES)


def _without_a_site(tmp_path: pathlib.Path) -> dict[str, Any]:
    """The group's annotation with the tables' row of one site of the family removed."""
    larger = _larger(tmp_path, repeats=False)
    folder = larger['sv_annotated'].tables / 'sv'
    (part,) = sorted(folder.iterdir())
    table = pq.read_table(part)
    pq.write_table(table.slice(1), part)
    return larger


def _other_sources(tmp_path: pathlib.Path) -> dict[str, Any]:
    return _larger(tmp_path, repeats=False, sources=('Ensembl', 'dbSNP'))


ANNOTATION_REFUSED: dict[str, tuple[Any, str]] = {
    "records without the annotation of the group's structural variants": (
        lambda tmp_path: _larger(tmp_path, repeats=False) | {'sv_annotated': None},
        "need the annotation of the group's sites",
    ),
    "the annotation for a family without a caller's VCF": (
        lambda tmp_path: _larger(tmp_path, repeats=False) | {'structural': _NO_STRUCTURAL, 'sv_alleles': None},
        "for a family without a caller's VCF",
    ),
    "a map of the group's sites without one of the family's records": (
        lambda tmp_path: _larger(
            tmp_path, repeats=False, group={'fam-1': (structural_vcfs.SV_RECORDS[1:], structural_vcfs.CNV_RECORDS)}
        ),
        r"records that the map of the group's sites lacks, such as \['inv:DEL:1'\]",
    ),
    "a map of the group's sites with a record that is not the family's": (
        lambda tmp_path: _larger(tmp_path, cnv=None, repeats=False, group={'fam-1': invented_sv_annotation.CALLERS}),
        r"records of the map of the group's sites that the family's callers lack, such as \['inv:CNV:1', 'inv:GAIN:1'",
    ),
    'an annotation without the site of a record': (_without_a_site, "records whose site the group's annotation lacks"),
    "an annotation of other sources than the alleles'": (
        _other_sources,
        r"annotated with catalog 209901 and the sources \['Ensembl', 'dbSNP'\]",
    ),
}


@pytest.mark.parametrize('case', ANNOTATION_REFUSED)
def test_an_annotation_that_is_not_the_familys_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    larger, message = ANNOTATION_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        family_tables.assemble(**{**_inputs(tmp_path), **larger(tmp_path)})
