"""What two plugins of Ensembl VEP say of a union's alleles, as a Parquet table (docs/annotation.md).

The Annotator gives a variant in a 5' untranslated region the consequence `5_prime_UTR_variant` and no more. The
UTRannotator plugin of Zhang, Wakeling, Ware and Whiffin (https://github.com/ImperialCardioGenetics/UTRannotator, MIT)
says which upstream open reading frames such a variant creates, removes or changes, and the NMD plugin of Ensembl's
collection says whether a variant that ends a reading frame early is predicted to escape nonsense-mediated decay. VEP
runs once for each gene model file (gene_models), on the union's alleles that overlap an exon, and writes one JSON
document per allele. This module turns those documents into one row per allele and transcript that a plugin evaluated,
keyed as `transcripts` is, with a record of what ran beside it.

The plugins are run as they are, and nothing here classifies a variant. A value is the plugin's, under a name in snake
case and with the type that the value has: UTRannotator writes its records as text, a number or `NA` for a distance and
`True` or `False` for a flag. The one thing read beside the plugins' own keys is how VEP placed the variant: the NMD
plugin writes a key only where it predicts an escape, so a row that it evaluated and did not flag is told by the
consequence terms and coding positions that make the plugin evaluate a row (_NMD_TERMS). No other value of VEP's is
kept.

A run is refused unless its positive controls come back annotated (gene_models.Control), which is how a plugin that
did not run is told from one that found nothing.
"""

from __future__ import annotations

import gzip
import itertools
import json
import pathlib
import re
from collections.abc import Callable, Iterator, Mapping, Sequence
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam

from case_annotation import annotation_tables, gene_models

TABLE_NAME = 'vep_plugins.parquet'
INFO_NAME = 'vep_plugins.json'
# What the step that runs VEP writes into its folder: the record of its run, and VEP's output for each source.
RUN_NAME = 'vep.json'
OUTPUTS: Mapping[gene_models.Source, str] = {'Ensembl': 'ensembl.json.gz', 'RefSeq': 'refseq.json.gz'}

_STRING = pa.string()
_INT = pa.int32()
_BOOL = pa.bool_()

# One record of UTRannotator: one upstream start codon that the allele gains or loses, or whose reading frame it
# changes. `consequence` is the kind of effect, and every other field is a key of the plugin's record of that kind,
# which a record of another kind does not have.
UORF = pa.struct(
    [
        pa.field('consequence', _STRING, nullable=False),
        pa.field('type', _STRING),
        pa.field('ref_type', _STRING),
        pa.field('alt_type', _STRING),
        pa.field('kozak_context', _STRING),
        pa.field('kozak_strength', _STRING),
        pa.field('cap_distance_to_start', _INT),
        pa.field('distance_to_cds', _INT),
        pa.field('distance_to_stop', _INT),
        pa.field('ref_start_distance_to_cds', _INT),
        pa.field('new_stop_distance_to_cds', _INT),
        pa.field('alt_stop', _BOOL),
        pa.field('alt_stop_distance_to_cds', _INT),
        pa.field('frame_with_cds', _STRING),
        pa.field('ref_type_length', _INT),
        pa.field('alt_type_length', _INT),
        pa.field('evidence', _BOOL),
    ]
)
# The columns that the table gives a row of `transcripts`, in their order.
COLUMNS = pa.schema(
    [
        pa.field('utr5_consequences', pa.list_(_STRING)),
        pa.field('utr5_uorfs', pa.list_(UORF)),
        pa.field('utr5_existing_uorfs', _INT),
        pa.field('utr5_existing_oorfs_out_of_frame', _INT),
        pa.field('utr5_existing_oorfs_in_frame', _INT),
        pa.field('nmd_escaping', _BOOL),
    ]
)
TABLE = pa.schema(
    [
        pa.field('variant_id', _STRING, nullable=False),
        pa.field('transcript_id', _STRING, nullable=False),
        pa.field('source', _STRING, nullable=False),
        *COLUMNS,
    ]
)

# A transcript id as `transcripts` holds one: Ensembl's or RefSeq's, with its version. NCBI's file also has a few
# features with other ids, which begin `rna-` or `id-`, and no row of `transcripts` can match those.
TRANSCRIPT_ID = re.compile(r'(ENST|NM_|NR_|XM_|XR_)\d+\.\d+')

_TYPES = frozenset({'uORF', 'inFrame_oORF', 'InFrame_oORF', 'OutOfFrame_oORF'})
_STRENGTHS = frozenset({'Weak', 'Moderate', 'Strong'})
_FRAMES = frozenset({'inFrame', 'outOfFrame'})


def _whole(value: str) -> int:
    if not re.fullmatch(r'-?\d+', value):
        raise ValueError(f'{value!r} is not a whole number')
    return int(value)


def _whole_or_na(value: str) -> int | None:
    return None if value == 'NA' else _whole(value)


def _flag(value: str) -> bool:
    if value not in ('True', 'False'):
        raise ValueError(f'{value!r} is neither True nor False')
    return value == 'True'


def _one_of(known: frozenset[str]) -> Callable[[str], str]:
    def read(value: str) -> str:
        if value not in known:
            raise ValueError(f'{value!r} is none of {sorted(known)}')
        return value

    return read


def _context(value: str) -> str | None:
    """A Kozak context, which the plugin writes as `-` for a start codon within three bases of an end of the UTR."""
    if value == '-':
        return None
    if not value:
        raise ValueError('an empty Kozak context')
    return value


def _strength(value: str) -> str | None:
    """A Kozak strength, which the plugin leaves empty where it gives no context."""
    return _one_of(_STRENGTHS)(value) if value else None


# Each key of a record of UTRannotator, after the kind that prefixes it: the field of UORF that holds its value, and
# how the plugin's text is read. `NA` is read only where the plugin's documentation says that it writes it.
_KEYS: Mapping[str, tuple[str, Callable[[str], Any]]] = {
    'type': ('type', _one_of(_TYPES)),
    'ref_type': ('ref_type', _one_of(_TYPES)),
    'alt_type': ('alt_type', _one_of(_TYPES)),
    'KozakContext': ('kozak_context', _context),
    'KozakStrength': ('kozak_strength', _strength),
    'CapDistanceToStart': ('cap_distance_to_start', _whole),
    'DistanceToCDS': ('distance_to_cds', _whole),
    'DistanceToStop': ('distance_to_stop', _whole_or_na),
    'ref_StartDistanceToCDS': ('ref_start_distance_to_cds', _whole),
    'newSTOPDistanceToCDS': ('new_stop_distance_to_cds', _whole),
    'AltStop': ('alt_stop', _flag),
    'AltStopDistanceToCDS': ('alt_stop_distance_to_cds', _whole_or_na),
    'FrameWithCDS': ('frame_with_cds', _one_of(_FRAMES)),
    'ref_type_length': ('ref_type_length', _whole_or_na),
    'alt_type_length': ('alt_type_length', _whole_or_na),
    'Evidence': ('evidence', _flag),
}
# The kinds of effect, and the keys that a record of each kind has: all of them, and no other.
KINDS: Mapping[str, frozenset[str]] = {
    'uAUG_gained': frozenset(
        {'type', 'KozakContext', 'KozakStrength', 'CapDistanceToStart', 'DistanceToCDS', 'DistanceToStop'}
    ),
    'uAUG_lost': frozenset(
        {'type', 'KozakContext', 'KozakStrength', 'CapDistanceToStart', 'DistanceToCDS', 'DistanceToStop', 'Evidence'}
    ),
    'uSTOP_gained': frozenset(
        {'ref_type', 'KozakContext', 'KozakStrength', 'ref_StartDistanceToCDS', 'newSTOPDistanceToCDS', 'Evidence'}
    ),
    'uSTOP_lost': frozenset(
        {'AltStop', 'AltStopDistanceToCDS', 'FrameWithCDS', 'KozakContext', 'KozakStrength', 'Evidence'}
    ),
    'uFrameShift': frozenset(
        {
            'ref_type',
            'ref_type_length',
            'alt_type',
            'alt_type_length',
            'ref_StartDistanceToCDS',
            'KozakContext',
            'KozakStrength',
            'Evidence',
        }
    ),
}

# The keys of a transcript in VEP's JSON that UTRannotator writes: the kinds of effect and their records, each as one
# text, and the three counts of the reading frames that the transcript's 5' UTR already has.
_CONSEQUENCE = 'five_prime_utr_variant_consequence'
_ANNOTATION = 'five_prime_utr_variant_annotation'
_COUNTS: Mapping[str, str] = {
    'existing_uorfs': 'utr5_existing_uorfs',
    'existing_outofframe_oorfs': 'utr5_existing_oorfs_out_of_frame',
    'existing_inframe_oorfs': 'utr5_existing_oorfs_in_frame',
}
# The NMD plugin's key and the one value it writes, and what makes the plugin evaluate a row: one of these consequence
# terms, on a variant with a start and an end in the coding sequence (the plugin's own two conditions).
_NMD = 'nmd'
_NMD_ESCAPING = 'NMD_escaping_variant'
_NMD_TERMS = frozenset({'stop_gained', 'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant'})


class FileRecord(pydantic.BaseModel, frozen=True, extra='forbid'):
    """A file as a run read it."""

    name: str
    bytes: int
    sha256: str = pydantic.Field(pattern=r'^[0-9a-f]{64}$')


class GeneModelRead(pydantic.BaseModel, frozen=True, extra='forbid'):
    source: gene_models.Source
    file: FileRecord


class VepRun(pydantic.BaseModel, frozen=True, extra='forbid'):
    """What the step that runs VEP records of its run (RUN_NAME)."""

    vep_version: str = pydantic.Field(min_length=1)
    # The commit of https://github.com/ImperialCardioGenetics/UTRannotator that the image holds the plugin at.
    utrannotator_commit: str = pydantic.Field(pattern=r'^[0-9a-f]{40}$')
    # The plugin's file of upstream reading frames that were seen translated, which gives a record its `evidence`.
    uorf_evidence: FileRecord
    nmd_plugin: FileRecord  # NMD.pm as the image holds it
    # VEP's options as the step gave them, without the paths of a run, and the variables of its environment.
    options: tuple[str, ...]
    environment: tuple[str, ...]
    gene_models: tuple[GeneModelRead, ...]


class GeneModelAnnotated(pydantic.BaseModel, frozen=True, extra='forbid'):
    """One gene model file of a run: the file, the public file it was prepared from, and what VEP gave with it."""

    source: gene_models.Source
    file: FileRecord
    prepared_from: gene_models.PublicFile
    rows: int  # the table's rows of this source
    # The rows that a plugin evaluated on a transcript whose id is no transcript id of `transcripts` (TRANSCRIPT_ID),
    # which the table leaves out.
    skipped_transcript_ids: int


class VepPluginsInfo(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The VEP run that the columns of the two plugins come from (INFO_NAME)."""

    vep_version: str
    utrannotator_commit: str
    uorf_evidence: FileRecord
    nmd_plugin: FileRecord
    options: tuple[str, ...]
    environment: tuple[str, ...]
    gene_models: tuple[GeneModelAnnotated, ...]
    alleles: int  # the union's alleles that overlap an exon of a gene model file, which VEP read


def uorf_records(consequence: str, annotation: str) -> list[dict[str, Any]]:
    """UTRannotator's records of one allele and transcript, each as an entry of UORF, in a fixed order.

    The plugin joins the kinds of its records with `&`, and the records in the same order with `&` as well. A record
    is its keys and values as `key:value`, joined with commas, and each key begins with the record's kind.

    Raises:
        ValueError: if the two texts hold different numbers of records, a kind is none of KINDS, a record's keys are
            not exactly its kind's, or a value is not what its key holds.
    """
    kinds, texts = consequence.split('&'), annotation.split('&')
    if len(kinds) != len(texts):
        raise ValueError(f'{len(kinds)} kinds in {consequence!r} and {len(texts)} records in {annotation!r}')
    records = []
    for kind, text in zip(kinds, texts, strict=True):
        if kind not in KINDS:
            raise ValueError(f'{kind!r} is no kind of effect that UTRannotator is known to write')
        values: dict[str, str] = {}
        for pair in text.split(','):
            prefixed, separator, value = pair.partition(':')
            key = prefixed.removeprefix(f'{kind}_')
            if not separator or key == prefixed or key in values:
                raise ValueError(f'{pair!r} is no key and value of a {kind} record, or repeats one: {text!r}')
            values[key] = value
        if values.keys() != KINDS[kind]:
            raise ValueError(f'a {kind} record has the keys {sorted(KINDS[kind])}, and this one {sorted(values)}')
        record: dict[str, Any] = {'consequence': kind}
        for key, value in values.items():
            field, read = _KEYS[key]
            try:
                record[field] = read(value)
            except ValueError as malformed:
                raise ValueError(f'{kind}_{key}: {malformed}') from malformed
        records.append(record)
    return sorted(records, key=_record_order)


def _record_order(record: Mapping[str, Any]) -> tuple[str, int, str]:
    """Where a record stands among those of its row: by kind, then from the 5' end by its start codon, then by text."""
    distances = [record.get('distance_to_cds'), record.get('ref_start_distance_to_cds')]
    start = next((d for d in distances if d is not None), None)
    return (record['consequence'], 0 if start is None else -start, json.dumps(record, sort_keys=True))


def _utr5(transcript: Mapping[str, Any]) -> dict[str, Any] | None:
    """The UTRannotator columns of one transcript of VEP's JSON, or None where the plugin did not evaluate it.

    The plugin evaluates a transcript whose 5' UTR the allele lies in, and writes its three counts for each such
    transcript, with or without an effect.
    """
    written = [key for key in (*_COUNTS, _CONSEQUENCE, _ANNOTATION) if key in transcript]
    if not written:
        return None
    if missing := [key for key in _COUNTS if key not in transcript]:
        raise ValueError(f'UTRannotator wrote {written} without {missing}')
    columns: dict[str, Any] = {}
    for key, column in _COUNTS.items():
        count = transcript[key]
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(f'{key} is {count!r}, not a count')
        columns[column] = count
    consequence, annotation = transcript.get(_CONSEQUENCE), transcript.get(_ANNOTATION)
    if consequence is None and annotation is None:
        return columns
    if not isinstance(consequence, str) or not isinstance(annotation, str):
        raise ValueError(f'{_CONSEQUENCE} is {consequence!r} and {_ANNOTATION} {annotation!r}, not two texts')
    records = uorf_records(consequence, annotation)
    return columns | {'utr5_consequences': sorted({r['consequence'] for r in records}), 'utr5_uorfs': records}


def _nmd(transcript: Mapping[str, Any]) -> bool | None:
    """Whether the NMD plugin predicts an escape on one transcript of VEP's JSON; None where it did not evaluate it."""
    flag = transcript.get(_NMD)
    if flag is not None:
        if flag != _NMD_ESCAPING:
            raise ValueError(f'{_NMD} is {flag!r}, and the NMD plugin is known to write {_NMD_ESCAPING} alone')
        return True
    terms = transcript.get('consequence_terms')
    if not isinstance(terms, list):
        raise ValueError(f'consequence_terms is {terms!r}, not a list of terms')
    evaluated = bool(_NMD_TERMS & set(terms)) and 'cds_start' in transcript and 'cds_end' in transcript
    return False if evaluated else None


def transcript_columns(transcript: Mapping[str, Any]) -> dict[str, Any] | None:
    """The plugins' columns of one transcript of VEP's JSON, or None where neither plugin evaluated it."""
    utr5, nmd = _utr5(transcript), _nmd(transcript)
    if utr5 is None and nmd is None:
        return None
    return (utr5 or {}) | {'nmd_escaping': nmd}


def read_given(path: pathlib.Path, controls: Sequence[gene_models.Control]) -> list[str]:
    """The ids of the records that VEP was given: the union's alleles that overlap an exon, without the controls.

    Raises:
        ValueError: if a record that is no control has an ID that is not the variant id of its allele, an ID is
            repeated, or a control is not among the records.
    """
    expected = {control.id for control in controls}
    seen: set[str] = set()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            where = f'{path}: {record.chrom}:{record.pos}'
            if record.id is None or record.id in seen:
                raise ValueError(f'{where} has no ID, or the ID of another record')
            if record.id not in expected and record.id != annotation_tables.record_id(record, str(path)):
                raise ValueError(f'{where} has the ID {record.id}, which is not the variant id of its allele')
            seen.add(record.id)
    if missing := sorted(expected - seen):
        raise ValueError(f'{path}: the positive controls {missing} are not among its records')
    return sorted(seen - expected)


def _documents(path: pathlib.Path) -> Iterator[dict[str, Any]]:
    """VEP's JSON documents, one a line and one for each record it was given."""
    with gzip.open(path, 'rt', encoding='utf-8') as lines:
        for number, line in enumerate(lines, start=1):
            try:
                document = json.loads(line)
            except json.JSONDecodeError as malformed:
                raise ValueError(f'{path}: line {number} is no JSON document') from malformed
            if not isinstance(document, dict) or not isinstance(document.get('id'), str):
                raise ValueError(f'{path}: line {number} is no document of a record with an id')
            yield document


def check_control(control: gene_models.Control, document: Mapping[str, Any], where: str) -> None:
    """Fails unless VEP's document of a positive control holds what the control's plugin must write for it.

    A plugin that did not load, a gene model file that VEP could not read or a reference under other names leaves the
    control without its keys, where an allele of the union would only have empty columns.
    """
    listed = [t for t in document.get('transcript_consequences', ()) if t['transcript_id'] == control.transcript_id]
    if len(listed) != 1:
        raise ValueError(
            f'{where}: VEP lists the transcript {control.transcript_id} {len(listed)} times for the positive control '
            f'{control.id}, so it did not read the gene models, the reference or the chromosome synonyms as prepared'
        )
    columns = transcript_columns(listed[0]) or {}
    passed = (
        'uAUG_gained' in (columns.get('utr5_consequences') or ())
        if control.plugin == 'UTRannotator'
        else columns.get('nmd_escaping') is True
    )
    if not passed:
        raise ValueError(
            f'{where}: the positive control {control.id} lacks what {control.plugin} writes for it on '
            f'{control.transcript_id}, so the plugin did not run: {dict(listed[0])}'
        )


def _source_rows(
    path: pathlib.Path,
    source: gene_models.Source,
    given: frozenset[str],
    controls: Sequence[gene_models.Control],
    passed_over: frozenset[str],
) -> tuple[list[dict[str, Any]], int]:
    """The rows of one VEP output, and how many more it would hold of transcripts whose id is no transcript id.

    `given` are the ids of the union's alleles that VEP read, and `controls` the positive controls of the output's own
    gene model file. Each run also reads the controls of the other file, whose ids are `passed_over`: this file's
    transcripts are not theirs.

    Raises:
        ValueError: if the documents are not exactly those of the records given, each once, a control fails
            (check_control), a plugin's keys are malformed, or a transcript is listed twice for an allele.
    """
    own = {control.id: control for control in controls}
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    skipped = 0
    for document in _documents(path):
        record_id = document['id']
        if record_id in seen:
            raise ValueError(f'{path}: two documents of {record_id}')
        seen.add(record_id)
        if record_id in own:
            check_control(own[record_id], document, str(path))
            continue
        if record_id in passed_over:
            continue
        if record_id not in given:
            raise ValueError(f'{path}: a document of {record_id}, which VEP was not given')
        listed: set[str] = set()
        for transcript in document.get('transcript_consequences', ()):
            try:
                columns = transcript_columns(transcript)
            except ValueError as malformed:
                raise ValueError(f'{path}: {record_id}, {transcript.get("transcript_id")}: {malformed}') from malformed
            if columns is None:
                continue
            transcript_id = transcript['transcript_id']
            if not TRANSCRIPT_ID.fullmatch(transcript_id):
                skipped += 1
                continue
            if transcript_id in listed:
                raise ValueError(f'{path}: {record_id}: a plugin evaluated the transcript {transcript_id} twice')
            listed.add(transcript_id)
            rows.append({'variant_id': record_id, 'transcript_id': transcript_id, 'source': source} | columns)
    if missing := sorted((given | own.keys() | passed_over) - seen)[:5]:
        raise ValueError(f'{path}: no document of records that VEP was given, such as {missing}')
    return rows, skipped


def convert(
    given: pathlib.Path, vep: pathlib.Path, prepared: gene_models.GeneModels, output: pathlib.Path
) -> VepPluginsInfo:
    """Writes the plugins' table of a union and the record of the run into the output folder.

    `given` is the VCF that VEP read: the union's alleles that overlap an exon, each under its variant id, and the
    positive controls. `vep` is the folder of the step that ran VEP, with its record and one output for each gene
    model file, and `prepared` the record of the gene models that the run read. The table's rows are sorted by variant
    id and transcript id, so that two runs of the same inputs write the same table.

    Raises:
        ValueError: if an input is not what the step before wrote (read_given, _source_rows), the run read other gene
            model files than the record names, or two sources give a row for the same allele and transcript.
    """
    run = VepRun.model_validate_json((vep / RUN_NAME).read_text())
    files: dict[gene_models.Source, FileRecord] = {read.source: read.file for read in run.gene_models}
    described: dict[gene_models.Source, gene_models.GeneModelFile] = {f.source: f for f in prepared.gene_models}
    if files.keys() != described.keys() or len(run.gene_models) != len(files):
        raise ValueError(f'VEP read the gene models of {sorted(files)}, and the record names {sorted(described)}')
    for source, file in described.items():
        if files[source].name != file.file:
            raise ValueError(f'VEP read {files[source].name} for {source}, and the record names {file.file}')
    controls = [control for file in prepared.gene_models for control in file.controls]
    alleles = read_given(given, controls)
    rows: list[dict[str, Any]] = []
    annotated = []
    for source, file in described.items():
        others = frozenset(control.id for control in controls) - {control.id for control in file.controls}
        found, skipped = _source_rows(vep / OUTPUTS[source], source, frozenset(alleles), file.controls, others)
        rows.extend(found)
        annotated.append(
            GeneModelAnnotated(
                source=source,
                file=files[source],
                prepared_from=file.prepared_from,
                rows=len(found),
                skipped_transcript_ids=skipped,
            )
        )
    rows.sort(key=lambda row: (row['variant_id'], row['transcript_id']))
    for row, following in itertools.pairwise(rows):
        if (row['variant_id'], row['transcript_id']) == (following['variant_id'], following['transcript_id']):
            raise ValueError(f'two sources give a row for {row["variant_id"]} and {row["transcript_id"]}')
    output.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=TABLE), output / TABLE_NAME, compression='zstd')
    info = VepPluginsInfo(
        vep_version=run.vep_version,
        utrannotator_commit=run.utrannotator_commit,
        uorf_evidence=run.uorf_evidence,
        nmd_plugin=run.nmd_plugin,
        options=run.options,
        environment=run.environment,
        gene_models=tuple(annotated),
        alleles=len(alleles),
    )
    (output / INFO_NAME).write_text(info.model_dump_json(indent=2) + '\n')
    return info


def read_info(folder: pathlib.Path) -> VepPluginsInfo:
    """The run's record in a folder that `convert` wrote."""
    return VepPluginsInfo.model_validate_json((folder / INFO_NAME).read_text())
