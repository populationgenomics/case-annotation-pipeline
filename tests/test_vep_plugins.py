"""The table of VEP's two plugins, from invented documents in the layout of VEP's JSON.

Every allele, transcript and value here is invented. The texts of UTRannotator's records are written as the plugin
writes them: the kinds joined with `&`, and each record as `kind_key:value` pairs joined with commas.
"""

from __future__ import annotations

import gzip
import json
import pathlib
import re
from collections.abc import Mapping, Sequence
from typing import Any

import pyarrow.parquet as pq
import pytest

from case_annotation import cli, contract, gene_models, vep_plugins

REPOSITORY = pathlib.Path(__file__).parent.parent
MODULES = (REPOSITORY / 'modules' / 'family_tables.nf').read_text('utf-8')
WORKFLOW = (REPOSITORY / 'family_tables.nf').read_text('utf-8')

# One record of each kind, with every key that the plugin writes for the kind.
_RECORDS: Mapping[str, Mapping[str, str]] = {
    'uAUG_gained': {
        'CapDistanceToStart': '40',
        'DistanceToCDS': '82',
        'DistanceToStop': '27',
        'KozakContext': 'GCCATGG',
        'KozakStrength': 'Strong',
        'type': 'uORF',
    },
    'uAUG_lost': {
        'CapDistanceToStart': '12',
        'DistanceToCDS': '110',
        'DistanceToStop': 'NA',
        'Evidence': 'True',
        'KozakContext': 'TTTATGC',
        'KozakStrength': 'Weak',
        'type': 'OutOfFrame_oORF',
    },
    'uSTOP_gained': {
        'Evidence': 'False',
        'KozakContext': 'ACCATGA',
        'KozakStrength': 'Moderate',
        'newSTOPDistanceToCDS': '31',
        'ref_StartDistanceToCDS': '70',
        'ref_type': 'InFrame_oORF',
    },
    'uSTOP_lost': {
        'AltStop': 'False',
        'AltStopDistanceToCDS': 'NA',
        'Evidence': 'True',
        'FrameWithCDS': 'outOfFrame',
        'KozakContext': '-',
        'KozakStrength': '',
    },
    'uFrameShift': {
        'Evidence': 'False',
        'KozakContext': 'GTTATGT',
        'KozakStrength': 'Moderate',
        'alt_type': 'OutOfFrame_oORF',
        'alt_type_length': 'NA',
        'ref_StartDistanceToCDS': '95',
        'ref_type': 'uORF',
        'ref_type_length': '36',
    },
}
_COUNTS = {'existing_uorfs': 2, 'existing_outofframe_oorfs': 1, 'existing_inframe_oorfs': 0}


def _text(kind: str, **changed: str | None) -> str:
    """One record of a kind as the plugin writes it, with the values given in place of the invented ones."""
    values = {**_RECORDS[kind], **changed}
    return ','.join(f'{kind}_{key}:{value}' for key, value in values.items() if value is not None)


def _utr(*records: tuple[str, str], counts: Mapping[str, Any] = _COUNTS) -> dict[str, Any]:
    """The keys that UTRannotator adds to a transcript of VEP's JSON, for records as (kind, text)."""
    written: dict[str, Any] = dict(counts)
    if records:
        written['five_prime_utr_variant_consequence'] = '&'.join(kind for kind, _text in records)
        written['five_prime_utr_variant_annotation'] = '&'.join(text for _kind, text in records)
    return written


def _transcript(transcript_id: str, terms: Sequence[str] = ('5_prime_UTR_variant',), **keys: object) -> dict[str, Any]:
    return {'transcript_id': transcript_id, 'consequence_terms': list(terms), 'impact': 'MODIFIER', **keys}


# ---- one transcript ----


@pytest.mark.parametrize('kind', _RECORDS)
def test_a_record_of_each_kind_has_the_plugins_values_under_its_fields(kind: str) -> None:
    (record,) = vep_plugins.uorf_records(kind, _text(kind))
    assert record['consequence'] == kind
    assert set(record) - {'consequence'} == {vep_plugins._KEYS[key][0] for key in vep_plugins.KINDS[kind]}
    # A field has the type of its column, and a value that the plugin leaves open is null.
    types = {field.name: str(field.type) for field in vep_plugins.UORF}
    expected = {'int32': int, 'bool': bool, 'string': str}
    assert all(value is None or type(value) is expected[types[field]] for field, value in record.items())
    written = _RECORDS[kind]
    for key, value in written.items():
        field = vep_plugins._KEYS[key][0]
        if value in ('NA', '-', ''):
            assert record[field] is None
        elif value in ('True', 'False'):
            assert record[field] is (value == 'True')
        else:
            assert str(record[field]) == value


def test_every_key_of_a_kind_has_a_field_and_every_field_a_key() -> None:
    keys = set().union(*vep_plugins.KINDS.values())
    assert keys == set(vep_plugins._KEYS)
    assert {field for field, _read in vep_plugins._KEYS.values()} | {'consequence'} == set(vep_plugins.UORF.names)
    assert set(_RECORDS) == set(vep_plugins.KINDS)
    assert all(set(_RECORDS[kind]) == vep_plugins.KINDS[kind] for kind in _RECORDS)


# Each key of a kind that holds the type of a reading frame.
_TYPE_KEYS = [
    ('uAUG_gained', 'type'),
    ('uAUG_lost', 'type'),
    ('uSTOP_gained', 'ref_type'),
    ('uFrameShift', 'ref_type'),
    ('uFrameShift', 'alt_type'),
]


def test_the_in_frame_type_has_one_spelling_where_the_plugin_writes_two() -> None:
    # The plugin writes the in-frame type with a small i for a gained start and with a capital for every other kind.
    (gained,) = vep_plugins.uorf_records('uAUG_gained', _text('uAUG_gained', type='inFrame_oORF'))
    (lost,) = vep_plugins.uorf_records('uAUG_lost', _text('uAUG_lost', type='InFrame_oORF'))
    assert gained['type'] == lost['type'] == 'InFrame_oORF'


@pytest.mark.parametrize(('kind', 'key'), _TYPE_KEYS)
def test_no_two_labels_of_a_type_differ_in_their_case_alone(kind: str, key: str) -> None:
    labels = {
        written: vep_plugins.uorf_records(kind, _text(kind, **{key: written}))[0][key] for written in vep_plugins._TYPES
    }
    # A label is a spelling that the plugin writes, changed in its case at most.
    assert all(label in labels and label.lower() == written.lower() for written, label in labels.items())
    assert len({label.lower() for label in labels.values()}) == len(set(labels.values()))


def test_every_key_that_holds_a_type_is_read_as_one() -> None:
    typed = {key for key, (_field, read) in vep_plugins._KEYS.items() if read is vep_plugins._type}
    assert typed == {key for key in vep_plugins._KEYS if key == 'type' or key.endswith('_type')}
    assert {(kind, key) for kind, keys in vep_plugins.KINDS.items() for key in keys & typed} == set(_TYPE_KEYS)


@pytest.mark.parametrize(('kind', 'key'), _TYPE_KEYS)
@pytest.mark.parametrize('written', ['inframe_oORF', 'Inframe_oORF', 'INFRAME_OORF', 'outOfFrame_oORF', 'UORF', ''])
def test_a_type_in_a_spelling_that_the_plugin_is_not_known_to_write_is_refused(
    kind: str, key: str, written: str
) -> None:
    with pytest.raises(ValueError, match=rf"{kind}_{key}: '{written}' is none of"):
        vep_plugins.uorf_records(kind, _text(kind, **{key: written}))


def test_several_kinds_and_several_records_of_one_kind_each_have_an_entry() -> None:
    near, far = _text('uFrameShift', ref_StartDistanceToCDS='20'), _text('uFrameShift', ref_StartDistanceToCDS='300')
    records = [('uFrameShift', near), ('uAUG_gained', _text('uAUG_gained')), ('uFrameShift', far)]
    columns = vep_plugins.transcript_columns(_transcript('ENST09.1', **_utr(*records)))
    assert columns is not None
    assert columns['utr5_consequences'] == ['uAUG_gained', 'uFrameShift']
    starts = ('distance_to_cds', 'ref_start_distance_to_cds')
    listed = [(r['consequence'], *(r.get(start) for start in starts)) for r in columns['utr5_uorfs']]
    # By kind, and within a kind from the 5' end.
    assert listed == [('uAUG_gained', 82, None), ('uFrameShift', None, 300), ('uFrameShift', None, 20)]
    # The order of the plugin's text, which follows Perl's hashes, does not reach the columns.
    again = vep_plugins.transcript_columns(_transcript('ENST09.1', **_utr(*reversed(records))))
    assert again == columns


def test_a_transcript_that_the_plugin_evaluated_without_an_effect_has_its_counts_alone() -> None:
    columns = vep_plugins.transcript_columns(_transcript('ENST09.1', **_utr()))
    assert columns == {
        'utr5_existing_uorfs': 2,
        'utr5_existing_oorfs_out_of_frame': 1,
        'utr5_existing_oorfs_in_frame': 0,
        'nmd_escaping': None,
    }


def test_a_transcript_that_no_plugin_evaluated_has_no_row() -> None:
    assert vep_plugins.transcript_columns(_transcript('ENST09.1', ('intron_variant',))) is None
    assert (
        vep_plugins.transcript_columns(_transcript('ENST09.1', ('missense_variant',), cds_start=4, cds_end=4)) is None
    )


@pytest.mark.parametrize('term', sorted(vep_plugins._NMD_TERMS))
def test_the_nmd_flag_is_true_where_written_false_where_evaluated_and_null_elsewhere(term: str) -> None:
    def escaping(**keys: object) -> bool | None:
        columns = vep_plugins.transcript_columns(_transcript('NM_09.1', (term, 'splice_region_variant'), **keys))
        return None if columns is None else columns['nmd_escaping']

    assert escaping(cds_start=40, cds_end=40, nmd='NMD_escaping_variant') is True
    # The plugin evaluates a variant of these terms that has a start and an end in the coding sequence.
    assert escaping(cds_start=400, cds_end=401) is False
    assert escaping() is None
    assert escaping(cds_start=400) is None
    assert escaping(cds_end=400) is None


def test_a_row_can_hold_both_plugins_values() -> None:
    transcript = _transcript('ENST09.1', ('5_prime_UTR_variant', 'frameshift_variant'), cds_start=1, cds_end=3)
    columns = vep_plugins.transcript_columns(transcript | _utr(('uAUG_lost', _text('uAUG_lost'))))
    assert columns is not None
    assert (columns['utr5_consequences'], columns['nmd_escaping']) == (['uAUG_lost'], False)


@pytest.mark.parametrize(
    ('consequence', 'annotation', 'message'),
    [
        ('uAUG_gained&uAUG_lost', _text('uAUG_gained'), '2 kinds in .* and 1 records'),
        ('uORF_gained', _text('uAUG_gained'), "'uORF_gained' is no kind of effect"),
        ('-', '-', "'-' is no kind of effect"),
        ('uAUG_gained', _text('uAUG_gained', Evidence='True'), 'a uAUG_gained record has the keys'),
        ('uAUG_lost', _text('uAUG_lost', Evidence=None), 'a uAUG_lost record has the keys'),
        ('uAUG_gained', _text('uAUG_lost', Evidence=None), 'is no key and value of a uAUG_gained record'),
        ('uAUG_gained', _text('uAUG_gained') + ',uAUG_gained_type:uORF', 'or repeats one'),
        ('uAUG_gained', _text('uAUG_gained').replace('type:uORF', 'type'), 'is no key and value'),
        ('uAUG_gained', _text('uAUG_gained', DistanceToCDS='near'), "uAUG_gained_DistanceToCDS: 'near' is not a whole"),
        ('uAUG_gained', _text('uAUG_gained', DistanceToCDS='NA'), "uAUG_gained_DistanceToCDS: 'NA' is not a whole"),
        ('uAUG_gained', _text('uAUG_gained', CapDistanceToStart='4.5'), "'4.5' is not a whole number"),
        ('uAUG_lost', _text('uAUG_lost', Evidence='NA'), "uAUG_lost_Evidence: 'NA' is neither True nor False"),
        ('uSTOP_lost', _text('uSTOP_lost', AltStop='yes'), "uSTOP_lost_AltStop: 'yes' is neither True nor False"),
        ('uAUG_gained', _text('uAUG_gained', type='oORF'), "uAUG_gained_type: 'oORF' is none of"),
        ('uFrameShift', _text('uFrameShift', alt_type='-'), "uFrameShift_alt_type: '-' is none of"),
        ('uSTOP_gained', _text('uSTOP_gained', KozakStrength='strong'), "'strong' is none of"),
        ('uSTOP_lost', _text('uSTOP_lost', FrameWithCDS='InFrame'), "uSTOP_lost_FrameWithCDS: 'InFrame' is none of"),
        ('uAUG_gained', _text('uAUG_gained', KozakContext=''), 'uAUG_gained_KozakContext: an empty Kozak context'),
    ],
)
def test_a_record_that_the_plugin_is_not_known_to_write_is_refused(
    consequence: str, annotation: str, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        vep_plugins.uorf_records(consequence, annotation)


@pytest.mark.parametrize(
    ('keys', 'message'),
    [
        ({'existing_uorfs': 2}, r"UTRannotator wrote \['existing_uorfs'\] without"),
        (
            {'five_prime_utr_variant_consequence': 'uAUG_gained', 'five_prime_utr_variant_annotation': 'x'},
            'UTRannotator wrote .* without',
        ),
        (_COUNTS | {'existing_uorfs': '2'}, "existing_uorfs is '2', not a count"),
        (_COUNTS | {'existing_inframe_oorfs': True}, 'existing_inframe_oorfs is True, not a count'),
        (_COUNTS | {'existing_outofframe_oorfs': -1}, 'existing_outofframe_oorfs is -1, not a count'),
        (_COUNTS | {'five_prime_utr_variant_consequence': 'uAUG_gained'}, 'not two texts'),
        (_COUNTS | {'five_prime_utr_variant_annotation': _text('uAUG_gained')}, 'not two texts'),
        (_COUNTS | {'nmd': 'NMD_variant'}, "nmd is 'NMD_variant'"),
    ],
)
def test_keys_that_a_plugin_is_not_known_to_write_are_refused(keys: Mapping[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        vep_plugins.transcript_columns(_transcript('ENST09.1', **keys))


def test_a_transcript_without_its_consequence_terms_is_refused() -> None:
    with pytest.raises(ValueError, match='consequence_terms is None'):
        vep_plugins.transcript_columns({'transcript_id': 'ENST09.1'})


@pytest.mark.parametrize('transcript_id', ['ENST00000000009.4', 'NM_000009.12', 'NR_09.1', 'XM_9.1', 'XR_09.3'])
def test_a_transcript_id_of_either_source_with_its_version_is_one(transcript_id: str) -> None:
    assert vep_plugins.TRANSCRIPT_ID.fullmatch(transcript_id)


@pytest.mark.parametrize('transcript_id', ['ENST00000000009', 'NM_000009', 'rna-INV-GTT27-1', 'id-INVG7', 'NP_09.1'])
def test_an_id_without_a_version_or_of_another_kind_is_no_transcript_id(transcript_id: str) -> None:
    assert not vep_plugins.TRANSCRIPT_ID.fullmatch(transcript_id)


# ---- a run ----

_ALLELES = (('chr1', 500, 'G', 'A'), ('chr1', 720, 'T', 'TA'), ('chr2', 90, 'CAT', 'C'), ('chr2', 300, 'A', 'G'))
_IDS = tuple(f'{chrom}-{pos}-{ref}-{alt}' for chrom, pos, ref, alt in _ALLELES)


def _prepared() -> gene_models.GeneModels:
    def file(source: gene_models.Source, transcript: str) -> gene_models.GeneModelFile:
        return gene_models.GeneModelFile(
            source=source,
            file=gene_models.FILES[source],
            bytes=1000,
            sha256='c' * 64,
            prepared_from=gene_models.PublicFile(
                name=f'invented.{source}.gff.gz', bytes=100, sha256='a' * 64, stated=(f'#!invented {source}',)
            ),
            controls=(
                gene_models.Control(
                    id=f'control-UTRannotator-{source}',
                    plugin='UTRannotator',
                    transcript_id=transcript,
                    chrom='chr1',
                    pos=100,
                    ref='C',
                    alt='CATG',
                ),
                gene_models.Control(
                    id=f'control-NMD-{source}',
                    plugin='NMD',
                    transcript_id=transcript,
                    chrom='chr1',
                    pos=200,
                    ref='T',
                    alt='TA',
                ),
            ),
        )

    return gene_models.GeneModels(
        gene_models=(file('Ensembl', 'ENST01.1'), file('RefSeq', 'NM_01.1')),
        assembly_report=gene_models.PublicFile(name='report.txt', bytes=10, sha256='b' * 64, stated=()),
    )


def _controls(prepared: gene_models.GeneModels) -> list[gene_models.Control]:
    return [control for file in prepared.gene_models for control in file.controls]


def _given(tmp_path: pathlib.Path, records: Sequence[tuple[str, int, str, str, str]] | None = None) -> pathlib.Path:
    """The VCF that VEP read: the invented alleles under their variant ids and the controls, or the records given."""
    if records is None:
        controls = [(c.chrom, c.pos, c.id, c.ref, c.alt) for c in _controls(_prepared())]
        alleles = [(chrom, pos, f'{chrom}-{pos}-{ref}-{alt}', ref, alt) for chrom, pos, ref, alt in _ALLELES]
        records = sorted(controls + alleles)
    lines = ['##fileformat=VCFv4.2', '##contig=<ID=chr1,length=1000>', '##contig=<ID=chr2,length=1000>']
    lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
    lines += [f'{chrom}\t{pos}\t{record_id}\t{ref}\t{alt}\t.\t.\t.' for chrom, pos, record_id, ref, alt in records]
    path = tmp_path / 'given.vcf'
    path.write_text('\n'.join(lines) + '\n')
    return path


def _control_documents(source: gene_models.Source) -> dict[str, list[dict[str, Any]]]:
    """VEP's transcripts of the four controls in the run of one source: each control is annotated on its own file."""
    documents: dict[str, list[dict[str, Any]]] = {}
    for control in _controls(_prepared()):
        if not control.id.endswith(source):
            documents[control.id] = []
        elif control.plugin == 'UTRannotator':
            gained = ('uAUG_gained', _text('uAUG_gained'))
            documents[control.id] = [_transcript(control.transcript_id, **_utr(gained))]
        else:
            keys = {'cds_start': 15, 'cds_end': 16, 'nmd': 'NMD_escaping_variant'}
            documents[control.id] = [_transcript(control.transcript_id, ('frameshift_variant',), **keys)]
    return documents


# What VEP lists for each invented allele with each gene model file: an allele in a 5' UTR with two kinds of effect on
# one transcript and none on another, one that ends a reading frame, one that neither plugin evaluates and one between
# genes. RefSeq's file also has a feature whose id is no transcript id of `transcripts`, and its row is kept.
_LISTED: Mapping[gene_models.Source, Mapping[str, list[dict[str, Any]]]] = {
    'Ensembl': {
        _IDS[0]: [
            _transcript('ENST05.2', **_utr(('uSTOP_lost', _text('uSTOP_lost')), ('uAUG_gained', _text('uAUG_gained')))),
            _transcript('ENST04.1', **_utr()),
            _transcript('ENST03.1', ('upstream_gene_variant',)),
        ],
        _IDS[1]: [_transcript('ENST05.2', ('frameshift_variant',), cds_start=300, cds_end=301)],
        _IDS[2]: [_transcript('ENST08.1', ('intron_variant',))],
        _IDS[3]: [],
    },
    'RefSeq': {
        _IDS[0]: [
            _transcript('NM_05.3', **_utr(('uAUG_gained', _text('uAUG_gained')))),
            _transcript('rna-INV-GTT27-1', **_utr()),
        ],
        _IDS[1]: [_transcript('NM_05.3', ('stop_gained',), cds_start=300, cds_end=300, nmd='NMD_escaping_variant')],
        _IDS[2]: [],
        _IDS[3]: [],
    },
}


def _run(sha256: str = 'c' * 64) -> dict[str, Any]:
    def recorded(name: str) -> dict[str, Any]:
        return {'name': name, 'bytes': 1000, 'sha256': sha256}

    return {
        'vep_version': '116.3',
        'utrannotator_commit': 'd' * 40,
        'uorf_evidence': recorded('invented-uorfs.txt'),
        'nmd_plugin': recorded('NMD.pm'),
        'options': ['--json', '--safe'],
        'environment': ['PERL_HASH_SEED=0'],
        'gene_models': [{'source': s, 'file': recorded(gene_models.FILES[s])} for s in gene_models.FILES],
    }


def _vep(
    tmp_path: pathlib.Path,
    listed: Mapping[gene_models.Source, Mapping[str, list[dict[str, Any]]]] = _LISTED,
    run: Mapping[str, Any] | None = None,
    reverse: bool = False,
) -> pathlib.Path:
    """The folder of the step that ran VEP: its record, and one document per record of each source."""
    folder = tmp_path / 'vep'
    folder.mkdir()
    (folder / vep_plugins.RUN_NAME).write_text(json.dumps(_run() if run is None else run))
    for source, name in vep_plugins.OUTPUTS.items():
        documents = _control_documents(source) | dict(listed[source])
        lines = [
            json.dumps(
                {'id': record_id, 'input': 'invented', 'most_severe_consequence': 'invented'}
                | ({'transcript_consequences': transcripts} if transcripts else {'intergenic_consequences': [{}]})
            )
            for record_id, transcripts in documents.items()
        ]
        with gzip.open(folder / name, 'wt', encoding='utf-8') as f:
            f.write(''.join(line + '\n' for line in (reversed(lines) if reverse else lines)))
    return folder


def _convert(tmp_path: pathlib.Path, **vep: Any) -> vep_plugins.VepPluginsInfo:  # noqa: ANN401 (what _vep takes)
    return vep_plugins.convert(_given(tmp_path), _vep(tmp_path, **vep), _prepared(), tmp_path / 'out')


def test_the_table_holds_a_row_for_each_allele_and_transcript_that_a_plugin_evaluated(tmp_path: pathlib.Path) -> None:
    _convert(tmp_path)
    table = pq.read_table(tmp_path / 'out' / vep_plugins.TABLE_NAME)
    assert table.schema.equals(vep_plugins.TABLE)
    rows = table.to_pylist()
    keys = [(r['variant_id'], r['transcript_id'], r['source']) for r in rows]
    assert keys == [
        (_IDS[0], 'ENST04.1', 'Ensembl'),
        (_IDS[0], 'ENST05.2', 'Ensembl'),
        (_IDS[0], 'NM_05.3', 'RefSeq'),
        (_IDS[0], 'rna-INV-GTT27-1', 'RefSeq'),
        (_IDS[1], 'ENST05.2', 'Ensembl'),
        (_IDS[1], 'NM_05.3', 'RefSeq'),
    ]
    by_key = {key[:2]: row for key, row in zip(keys, rows, strict=True)}
    both = by_key[(_IDS[0], 'ENST05.2')]
    assert both['utr5_consequences'] == ['uAUG_gained', 'uSTOP_lost']
    assert [r['consequence'] for r in both['utr5_uorfs']] == ['uAUG_gained', 'uSTOP_lost']
    assert (both['utr5_existing_uorfs'], both['nmd_escaping']) == (2, None)
    # Evaluated, with no effect: the counts, and no list.
    none = by_key[(_IDS[0], 'ENST04.1')]
    assert (none['utr5_consequences'], none['utr5_uorfs'], none['utr5_existing_oorfs_out_of_frame']) == (None, None, 1)
    assert [by_key[(_IDS[1], t)]['nmd_escaping'] for t in ('ENST05.2', 'NM_05.3')] == [False, True]
    assert by_key[(_IDS[1], 'NM_05.3')]['utr5_existing_uorfs'] is None
    # No row of a control or of an allele that no plugin evaluated.
    assert {r['variant_id'] for r in rows} == {_IDS[0], _IDS[1]}
    # A feature of RefSeq's file whose id is no transcript id has its row, which no row of `transcripts` matches.
    assert [vep_plugins.TRANSCRIPT_ID.fullmatch(key[1]) is None for key in keys] == [
        False,
        False,
        False,
        True,
        False,
        False,
    ]


def test_two_runs_whose_documents_come_in_another_order_write_the_same_table(tmp_path: pathlib.Path) -> None:
    first, second = tmp_path / 'first', tmp_path / 'second'
    first.mkdir()
    second.mkdir()
    _convert(first)
    _convert(second, reverse=True)
    name = pathlib.Path('out') / vep_plugins.TABLE_NAME
    assert (first / name).read_bytes() == (second / name).read_bytes()


def test_the_record_holds_the_run_and_each_gene_model_file(tmp_path: pathlib.Path) -> None:
    info = _convert(tmp_path)
    assert info == vep_plugins.read_info(tmp_path / 'out')
    run = _run()
    assert (info.vep_version, info.utrannotator_commit) == (run['vep_version'], run['utrannotator_commit'])
    assert (list(info.options), list(info.environment)) == (run['options'], run['environment'])
    assert info.uorf_evidence.name == 'invented-uorfs.txt'
    described = {file.source: file for file in _prepared().gene_models}
    assert [m.source for m in info.gene_models] == list(gene_models.FILES)
    for model in info.gene_models:
        assert model.file.name == gene_models.FILES[model.source]
        assert model.prepared_from == described[model.source].prepared_from


def test_the_record_is_the_same_whatever_alleles_the_union_holds(tmp_path: pathlib.Path) -> None:
    # A family's tables carry the record, so it holds nothing that changes with the other families of a group.
    first, second = tmp_path / 'first', tmp_path / 'second'
    first.mkdir()
    second.mkdir()
    _convert(first)
    fewer = {source: {**listed, _IDS[0]: []} for source, listed in _LISTED.items()}
    _convert(second, listed=fewer)
    name = pathlib.Path('out') / vep_plugins.INFO_NAME
    assert (first / name).read_bytes() == (second / name).read_bytes()
    assert pq.read_metadata(first / 'out' / vep_plugins.TABLE_NAME).num_rows == 6
    assert pq.read_metadata(second / 'out' / vep_plugins.TABLE_NAME).num_rows == 2


def test_the_command_writes_the_table_and_the_record(tmp_path: pathlib.Path) -> None:
    record = tmp_path / gene_models.RECORD_NAME
    record.write_text(_prepared().model_dump_json())
    arguments = ['--given', str(_given(tmp_path)), '--vep', str(_vep(tmp_path)), '--gene-models', str(record)]
    cli.main(['vep-plugins', *arguments, '--output', str(tmp_path / 'out')])
    assert {path.name for path in (tmp_path / 'out').iterdir()} == {vep_plugins.TABLE_NAME, vep_plugins.INFO_NAME}


def _replaced(source: gene_models.Source, record_id: str, transcripts: list[dict[str, Any]]) -> dict[str, Any]:
    return {**_LISTED, source: {**_LISTED[source], record_id: transcripts}}


@pytest.mark.parametrize(
    ('source', 'control', 'transcripts', 'message'),
    [
        # The plugin did not run: VEP lists the transcript in its 5' UTR, without the plugin's keys.
        ('Ensembl', 'control-UTRannotator-Ensembl', [_transcript('ENST01.1')], 'lacks what UTRannotator writes'),
        # It ran and found another kind of effect than the gained start.
        (
            'Ensembl',
            'control-UTRannotator-Ensembl',
            [_transcript('ENST01.1', **_utr(('uAUG_lost', _text('uAUG_lost'))))],
            'lacks what UTRannotator writes',
        ),
        (
            'RefSeq',
            'control-NMD-RefSeq',
            [_transcript('NM_01.1', ('frameshift_variant',), cds_start=15, cds_end=16)],
            'lacks what NMD writes',
        ),
        # VEP did not read the gene models: the control's transcript is not among those it lists.
        (
            'RefSeq',
            'control-NMD-RefSeq',
            [_transcript('NM_02.1', ('intron_variant',))],
            'lists the transcript NM_01.1 0',
        ),
        ('Ensembl', 'control-UTRannotator-Ensembl', [], 'lists the transcript ENST01.1 0 times'),
    ],
)
def test_a_control_without_what_its_plugin_writes_stops_the_run(
    tmp_path: pathlib.Path, source: gene_models.Source, control: str, transcripts: list[dict[str, Any]], message: str
) -> None:
    with pytest.raises(ValueError, match=f'{message}.*' if 'lists' in message else f'{control} {message}'):
        _convert(tmp_path, listed=_replaced(source, control, transcripts))
    assert not (tmp_path / 'out').exists()


def test_a_document_of_a_record_that_vep_was_not_given_is_refused(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='a document of chr9-1-A-C, which VEP was not given'):
        _convert(tmp_path, listed=_replaced('Ensembl', 'chr9-1-A-C', []))


def test_a_record_without_a_document_is_refused(tmp_path: pathlib.Path) -> None:
    without = {**_LISTED, 'RefSeq': {k: v for k, v in _LISTED['RefSeq'].items() if k != _IDS[3]}}
    with pytest.raises(ValueError, match=rf'refseq\.json\.gz: no document of records .* such as \[{_IDS[3]!r}\]'):
        _convert(tmp_path, listed=without)


def test_two_documents_of_one_record_are_refused(tmp_path: pathlib.Path) -> None:
    folder = _vep(tmp_path)
    output = folder / vep_plugins.OUTPUTS['Ensembl']
    with gzip.open(output, 'rt') as f:
        lines = f.readlines()
    with gzip.open(output, 'wt') as f:
        f.writelines([*lines, lines[-1]])
    with pytest.raises(ValueError, match='two documents of'):
        vep_plugins.convert(_given(tmp_path), folder, _prepared(), tmp_path / 'out')


def test_a_line_that_is_no_document_is_refused(tmp_path: pathlib.Path) -> None:
    folder = _vep(tmp_path)
    with gzip.open(folder / vep_plugins.OUTPUTS['Ensembl'], 'at') as f:
        f.write('{"id": "cut sho\n')
    with pytest.raises(ValueError, match=r'ensembl\.json\.gz: line 9 is no JSON document'):
        vep_plugins.convert(_given(tmp_path), folder, _prepared(), tmp_path / 'out')


def test_a_document_without_an_id_is_refused(tmp_path: pathlib.Path) -> None:
    folder = _vep(tmp_path)
    with gzip.open(folder / vep_plugins.OUTPUTS['RefSeq'], 'at') as f:
        f.write('{"input": "invented"}\n')
    with pytest.raises(ValueError, match=r'refseq\.json\.gz: line 9 is no document of a record with an id'):
        vep_plugins.convert(_given(tmp_path), folder, _prepared(), tmp_path / 'out')


def test_a_transcript_that_a_plugin_evaluated_twice_for_an_allele_is_refused(tmp_path: pathlib.Path) -> None:
    twice = [_transcript('ENST04.1', **_utr()), _transcript('ENST04.1', **_utr())]
    with pytest.raises(ValueError, match=f'{_IDS[3]}: a plugin evaluated the transcript ENST04.1 twice'):
        _convert(tmp_path, listed=_replaced('Ensembl', _IDS[3], twice))


def test_a_row_of_one_allele_and_transcript_from_both_sources_is_refused(tmp_path: pathlib.Path) -> None:
    shared = [_transcript('ENST05.2', **_utr())]
    with pytest.raises(ValueError, match=f'two sources give a row for {_IDS[0]} and ENST05.2'):
        _convert(tmp_path, listed=_replaced('RefSeq', _IDS[0], shared))


def test_malformed_keys_of_a_plugin_are_reported_with_their_allele_and_transcript(tmp_path: pathlib.Path) -> None:
    malformed = [_transcript('NM_05.3', **_utr(('uAUG_gained', _text('uAUG_gained', type='oORF'))))]
    with pytest.raises(ValueError, match=rf"refseq\.json\.gz: {_IDS[2]}, NM_05\.3: uAUG_gained_type: 'oORF'"):
        _convert(tmp_path, listed=_replaced('RefSeq', _IDS[2], malformed))


def test_a_run_that_read_other_gene_models_than_the_record_names_is_refused(tmp_path: pathlib.Path) -> None:
    one = _run() | {'gene_models': _run()['gene_models'][:1]}
    with pytest.raises(ValueError, match=re.escape("VEP read the gene models of ['Ensembl'], and the record names")):
        _convert(tmp_path, run=one)


def test_a_run_that_read_a_source_twice_is_refused(tmp_path: pathlib.Path) -> None:
    models = _run()['gene_models']
    with pytest.raises(ValueError, match='VEP read the gene models of'):
        _convert(tmp_path, run=_run() | {'gene_models': [*models, models[0]]})


@pytest.mark.parametrize(
    'changed', [{'name': 'other.gff3.gz'}, {'bytes': 999}, {'sha256': 'e' * 64}], ids=['name', 'size', 'checksum']
)
def test_a_run_that_read_another_file_than_the_prepared_one_is_refused(
    tmp_path: pathlib.Path, changed: Mapping[str, Any]
) -> None:
    models = _run()['gene_models']
    other = [{'source': 'Ensembl', 'file': models[0]['file'] | changed}, models[1]]
    with pytest.raises(ValueError, match=r'for Ensembl, which is not the file ensembl\.gff3\.gz of 1000 bytes'):
        _convert(tmp_path, run=_run() | {'gene_models': other})


@pytest.mark.parametrize(
    'run',
    [
        _run(sha256='not a checksum'),
        _run() | {'utrannotator_commit': '3d9d466'},
        _run() | {'vep_version': ''},
        _run() | {'unknown': 1},
        {k: v for k, v in _run().items() if k != 'options'},
    ],
    ids=['a checksum', 'a short commit', 'no version', 'an unknown key', 'no options'],
)
def test_a_record_of_the_run_that_the_step_does_not_write_is_refused(tmp_path: pathlib.Path, run: Any) -> None:  # noqa: ANN401 (a JSON value)
    with pytest.raises(ValueError, match='VepRun'):
        _convert(tmp_path, run=run)


def _records(changed: Mapping[str, tuple[str, int, str, str, str] | None]) -> list[tuple[str, int, str, str, str]]:
    controls = {c.id: (c.chrom, c.pos, c.id, c.ref, c.alt) for c in _controls(_prepared())}
    alleles = {
        f'{chrom}-{pos}-{ref}-{alt}': (chrom, pos, f'{chrom}-{pos}-{ref}-{alt}', ref, alt)
        for chrom, pos, ref, alt in _ALLELES
    }
    return sorted(record for record in ({**controls, **alleles} | dict(changed)).values() if record is not None)


@pytest.mark.parametrize(
    ('changed', 'message'),
    [
        ({_IDS[0]: ('chr1', 500, 'chr1-500-G-T', 'G', 'A')}, 'chr1:500 has the ID chr1-500-G-T, which is not'),
        ({_IDS[0]: ('chr1', 500, '.', 'G', 'A')}, 'chr1:500 has no ID'),
        ({'again': ('chr1', 720, _IDS[1], 'T', 'TA')}, 'chr1:720 has no ID, or the ID of another record'),
        ({'control-NMD-RefSeq': None}, r"the positive controls \['control-NMD-RefSeq'\] are not among"),
    ],
    ids=['an id of another allele', 'no id', 'a repeated id', 'a missing control'],
)
def test_a_vcf_that_the_step_before_did_not_write_is_refused(
    tmp_path: pathlib.Path, changed: Mapping[str, tuple[str, int, str, str, str] | None], message: str
) -> None:
    given = _given(tmp_path, _records(changed))
    with pytest.raises(ValueError, match=message):
        vep_plugins.read_given(given, _controls(_prepared()))


def test_the_alleles_that_vep_read_are_the_records_without_the_controls(tmp_path: pathlib.Path) -> None:
    assert vep_plugins.read_given(_given(tmp_path), _controls(_prepared())) == sorted(_IDS)


# ---- the steps that give the command its inputs ----


def test_the_workflow_names_the_files_as_the_tools_do() -> None:
    # VEP_PLUGINS reads each prepared file by its name and writes its outputs and its record under theirs.
    models = ' '.join(
        f'{source}:{gene_models.FILES[source]}:{vep_plugins.OUTPUTS[source].removesuffix(".json.gz")}'
        for source in gene_models.FILES
    )
    assert f'for model in {models}; do' in MODULES
    assert '--output_file vep/\\${name}.json.gz' in MODULES
    assert f'--synonyms "\\${{inputs}}/gene-models/{gene_models.SYNONYMS_NAME}"' in MODULES
    for name in (gene_models.EXONS_NAME, gene_models.CONTROLS_NAME, gene_models.RECORD_NAME):
        assert f"settings.gene_models.resolve('{name}')" in WORKFLOW
    # The workflow looks for every file of the prepared folder before a task runs.
    (listed,) = re.findall(r'\ndef geneModelFiles\(\) \{\n    return \[([^\]]+)\]\n\}', WORKFLOW)
    prepared = [*gene_models.FILES.values(), *(f'{name}.tbi' for name in gene_models.FILES.values())]
    prepared += [gene_models.SYNONYMS_NAME, gene_models.EXONS_NAME, gene_models.CONTROLS_NAME, gene_models.RECORD_NAME]
    assert sorted(name.strip().strip("'") for name in listed.split(',')) == sorted(prepared)
    assert 'geneModelFiles().each { name ->\n        if (!settings.gene_models.resolve(name).exists()) {' in WORKFLOW
    stub = f'touch vep_plugins/{vep_plugins.TABLE_NAME} vep_plugins/{vep_plugins.INFO_NAME}'
    assert stub in MODULES
    # The record's keys are the model's.
    (record,) = re.findall(r"\n    printf '(\{.*\})\\\\n' ", MODULES)
    assert re.findall(r'"(\w+)": ', record) == list(vep_plugins.VepRun.model_fields)
    assert f'> vep/{vep_plugins.RUN_NAME}\n' in MODULES


def test_the_step_fixes_perls_hash_order_and_fails_for_a_plugin_that_does_not_load() -> None:
    (script,) = re.findall(r'\nprocess VEP_PLUGINS \{\n.*?\n    script:\n(.*?)\n    stub:\n', MODULES, flags=re.DOTALL)
    assert 'export PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0' in script
    (options,) = re.findall(r'^    options=\((.*)\)$', script, flags=re.MULTILINE)
    assert {'--json', '--safe', '--transcript_version'} <= set(options.split())
    assert '--offline' not in script
    assert '--dir_plugins "\\${UTRANNOTATOR_DIR}"' in script
    assert "Plugin '[^']+' (went wrong|did not return)|Failed to (compile|instantiate) plugin" in script


def test_the_step_fails_for_a_run_in_which_vep_did_not_set_its_offline_mode() -> None:
    (script,) = re.findall(r'\nprocess VEP_PLUGINS \{\n.*?\n    script:\n(.*?)\n    stub:\n', MODULES, flags=re.DOTALL)
    assert "if ! grep -q 'Setting offline mode' vep/\\${name}.warnings.txt; then" in script
    assert '--database' not in script
    assert '--cache' not in script


def test_vep_reads_the_reference_and_the_gene_models_from_copies_on_the_tasks_disk() -> None:
    (process,) = re.findall(r'\nprocess VEP_PLUGINS \{\n(.*?)\n    stub:\n', MODULES, flags=re.DOTALL)
    directives, script = process.split('\n    script:\n')
    # Both are copied before VEP starts, with what a staged link points to.
    copied = script.index('cp -L ${reference} ${reference_index} "\\${inputs}/"\n')
    assert copied < script.index('cp -rL gene-models "\\${inputs}/gene-models"\n') < script.index('vep --input_file')
    # The folder of the copies has a fixed absolute path, since VEP's warnings quote it.
    (folder,) = re.findall(r'^    inputs=(\S+)$', script, flags=re.MULTILINE)
    assert re.fullmatch(r'(/[\w-]+)+', folder)
    # Every file that VEP reads many times is a copy: the command names neither staged path.
    (command,) = re.findall(r'\n        (vep --input_file .*?"\\\$\{options\[@\]\}")\n', script, flags=re.DOTALL)
    given = dict(re.findall(r'(--gff|--fasta|--synonyms) (\S+)', command))
    assert given == {
        '--gff': '"\\${inputs}/gene-models/\\${file}"',
        '--fasta': '"\\${inputs}/${reference.name}"',
        '--synonyms': f'"\\${{inputs}}/gene-models/{gene_models.SYNONYMS_NAME}"',
    }
    # The record's checksum is of the copy that VEP read.
    assert 'recorded "\\${inputs}/gene-models/\\${file}"' in script
    assert 'recorded gene-models/' not in script
    # A run that stalls ends.
    assert len(re.findall(r'^    time \d+\.h$', directives, flags=re.MULTILINE)) == 1


def test_every_column_of_the_table_has_a_type_that_the_contract_names() -> None:
    described = [contract._column(field) for field in vep_plugins.COLUMNS]
    assert [column['name'] for column in described] == vep_plugins.COLUMNS.names
    assert set(vep_plugins.TABLE.names) - set(vep_plugins.COLUMNS.names) == {'variant_id', 'transcript_id', 'source'}
