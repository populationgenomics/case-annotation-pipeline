"""What the rules of two plugins of Ensembl VEP say of a union's alleles, as a Parquet table (docs/annotation.md).

The Annotator gives a variant in a 5' untranslated region the consequence `5_prime_UTR_variant` and no more. The
UTRannotator plugin says which upstream open reading frames such a variant creates, removes or changes, and the NMD
plugin of Ensembl's collection says whether a variant that ends a reading frame early is predicted to escape
nonsense-mediated decay. This module does not run VEP. It applies the plugins' rules, restated in Python (utr5, nmd),
to the transcripts of the prepared gene model files (gene_models, transcript_models), for the alleles and transcripts
that VEP would give to each plugin (plugin_gates), and writes one row per allele and transcript that a rule
evaluated, keyed as `transcripts` is, with a record of what the values come from beside it. The record holds nothing
of the union's own, such as a count of its alleles, so that it is the same for every union that the same image and
gene models annotate.

A row has the plugins' values under names in snake case, each with the type that the value has. Beside them it has
three things that no plugin gives: the Kozak context and strength of the transcript's own start codon, the same two
with the allele where the allele changes them, and the reason for a predicted escape, which only a later edition of
the NMD plugin writes.

The alleles are those of the union that overlap an exon, with the positive controls of the prepared gene models
among them (gene_models.Control). A run is refused unless each control comes back with the value that its rule gives
it, which is how gene models, a reference or chromosome names that do not fit one another are told from a union
without such variants.
"""

from __future__ import annotations

import concurrent.futures
import dataclasses
import hashlib
import itertools
import json
import multiprocessing
import pathlib
import re
import time
from collections.abc import Mapping, Sequence
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pydantic
import pysam

from case_annotation import annotation_tables, gene_models, nmd, plugin_gates, transcript_models, utr5

TABLE_NAME = 'vep_plugins.parquet'
INFO_NAME = 'vep_plugins.json'
# A process writes a line to the log for every so many alleles of a chromosome that it has read with a gene model file.
PROGRESS_EVERY = 20_000

# The release of Ensembl VEP whose reading of a gene model file and whose consequence terms the step follows
# (transcript_models, plugin_gates): the one that computed these columns before the rules were restated.
VEP_RELEASE = '116.3'

_STRING = pa.string()
_INT = pa.int32()
_BOOL = pa.bool_()

# One record of UTRannotator's rules: one upstream start codon that the allele gains or loses, or whose reading frame
# it changes. `consequence` is the kind of effect, and every other field is a key of the plugin's record of that
# kind, which a record of another kind does not have.
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
        pa.field('utr5_cds_kozak_context', _STRING),
        pa.field('utr5_cds_kozak_strength', _STRING),
        pa.field('utr5_cds_kozak_context_alt', _STRING),
        pa.field('utr5_cds_kozak_strength_alt', _STRING),
        pa.field('nmd_escape_reason', _STRING),
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
# features with other ids, which begin `rna-` or `id-`. The table keeps their rows, and no row of `transcripts` can
# match one.
TRANSCRIPT_ID = re.compile(r'(ENST|NM_|NR_|XM_|XR_)\d+\.\d+')

# The kinds of effect, and the fields of UORF that a record of each kind has beside `consequence`: all of them, and
# no other.
KINDS: Mapping[str, frozenset[str]] = {
    'uAUG_gained': frozenset(
        {'type', 'kozak_context', 'kozak_strength', 'cap_distance_to_start', 'distance_to_cds', 'distance_to_stop'}
    ),
    'uAUG_lost': frozenset(
        {
            'type',
            'kozak_context',
            'kozak_strength',
            'cap_distance_to_start',
            'distance_to_cds',
            'distance_to_stop',
            'evidence',
        }
    ),
    'uSTOP_gained': frozenset(
        {
            'ref_type',
            'kozak_context',
            'kozak_strength',
            'ref_start_distance_to_cds',
            'new_stop_distance_to_cds',
            'evidence',
        }
    ),
    'uSTOP_lost': frozenset(
        {'alt_stop', 'alt_stop_distance_to_cds', 'frame_with_cds', 'kozak_context', 'kozak_strength', 'evidence'}
    ),
    'uFrameShift': frozenset(
        {
            'ref_type',
            'ref_type_length',
            'alt_type',
            'alt_type_length',
            'ref_start_distance_to_cds',
            'kozak_context',
            'kozak_strength',
            'evidence',
        }
    ),
}

_ALLELE = re.compile(r'[ACGT]+')
_START_CODON = 'ATG'


class FileRecord(pydantic.BaseModel, frozen=True, extra='forbid'):
    """A file as a run read it."""

    name: str
    bytes: int
    sha256: str = pydantic.Field(pattern=r'^[0-9a-f]{64}$')


class GeneModelAnnotated(pydantic.BaseModel, frozen=True, extra='forbid'):
    """One gene model file of a run: the file as the run read it, and the public file that it was prepared from."""

    source: gene_models.Source
    file: FileRecord
    prepared_from: gene_models.PublicFile


class VepPluginsInfo(pydantic.BaseModel, frozen=True, extra='forbid'):
    """What the columns of the two plugins' rules come from (INFO_NAME)."""

    # The commit of https://github.com/ImperialCardioGenetics/UTRannotator that the rules of the utr5_ columns were
    # restated from (utr5), and the commit of https://github.com/Ensembl/VEP_plugins for the rule of the nmd_ columns
    # (nmd).
    utrannotator_commit: str = pydantic.Field(pattern=r'^[0-9a-f]{40}$')
    nmd_plugin_commit: str = pydantic.Field(pattern=r'^[0-9a-f]{40}$')
    # The release of Ensembl VEP that the step follows in reading the gene models and in choosing the alleles that a
    # rule evaluates.
    vep_release: str = pydantic.Field(min_length=1)
    # The plugin's file of upstream reading frames that were seen translated, which gives a record its `evidence`.
    uorf_evidence: FileRecord
    gene_models: tuple[GeneModelAnnotated, ...]


def _recorded(path: pathlib.Path) -> FileRecord:
    with path.open('rb') as f:
        return FileRecord(
            name=path.name, bytes=path.stat().st_size, sha256=hashlib.file_digest(f, 'sha256').hexdigest()
        )


def _record_order(record: Mapping[str, Any]) -> tuple[str, int, str]:
    """Where a record stands among those of its row: by kind, then from the 5' end by its start codon, then by text.

    A `uSTOP_lost` record has no distance of its start codon, so those of a row stand by their text alone.
    """
    distances = [record.get('distance_to_cds'), record.get('ref_start_distance_to_cds')]
    start = next((d for d in distances if d is not None), None)
    return (record['consequence'], 0 if start is None else -start, json.dumps(record, sort_keys=True))


def utr_change(transcript: transcript_models.Transcript, allele: plugin_gates.Allele, utr: str) -> utr5.Change | None:
    """An allele as a change of a transcript's 5' UTR, or None for one that the rules do not read.

    `utr` is the UTR's sequence. The plugin reads an allele whose bases all lie in the 5' UTR, and an insertion whose
    following base on the transcript does. Two things are not the plugin's.

    - The plugin reads a deletion that begins in one exon of the UTR and ends in another as if the intron between them
      were not there. That is not read here.
    - An insertion directly before the first coding base has no following base in the UTR, and the plugin does not
      read it. Here it is read where it repeats the UTR's last base at its end, because the same change can then be
      written inside the UTR, where the plugin reads it.
    """
    coding = transcript.coding
    if coding is None:
        return None

    def place(position: int) -> int | None:
        cdna = transcript.cdna(position)
        return cdna - 1 if cdna is not None and cdna < coding.cdna_start else None

    removed = len(allele.ref)
    low, high = allele.start, allele.start + removed - 1
    inserted = allele.alt if transcript.strand == 1 else transcript_models.reverse_complement(allele.alt)
    first = place(low if transcript.strand == 1 else high)
    last = first if not removed else place(high if transcript.strand == 1 else low)
    if first is None and not removed:
        # The two bases beside the insertion, on the transcript: the UTR's last base and the first coding base.
        before, behind = (allele.end, allele.start) if transcript.strand == 1 else (allele.start, allele.end)
        at_the_end = place(before) == len(utr) - 1 and transcript.cdna(behind) == coding.cdna_start
        if at_the_end and utr and inserted[-1:] == utr[-1]:
            return utr5.Change(position=len(utr), removed=0, inserted=inserted)
    if first is None or last is None or (removed and last - first != removed - 1):
        return None
    return utr5.Change(position=first, removed=removed, inserted=inserted)


class Annotated:
    """A transcript with what the rules read of it, each read once."""

    def __init__(self, transcript: transcript_models.Transcript, reference: pysam.FastaFile) -> None:
        if transcript.coding is None:
            raise ValueError(f'{transcript.transcript_id} has no coding sequence')
        self.transcript = transcript
        self.sequences = transcript_models.sequences(transcript, reference)
        self._coding = transcript.coding
        self._frames: utr5.ReadingFrames | None = None
        self._existing: dict[str, Any] | None = None

    @property
    def frames(self) -> utr5.ReadingFrames:
        if self._frames is None:
            self._frames = utr5.ReadingFrames(self.sequences.utr5, self.sequences.coding)
        return self._frames

    @property
    def own_start(self) -> int | None:
        """Where the transcript's start codon begins in its spliced sequence, counted from 0.

        None for a transcript whose coding sequence does not begin with a whole `ATG`.
        """
        start = self._coding.cdna_start - 1
        whole = self._coding.phase == 0 and self.sequences.spliced[start : start + 3] == _START_CODON
        return start if whole else None

    @property
    def of_utr(self) -> dict[str, Any]:
        """The columns that every allele in the transcript's 5' UTR has: the counts and the own start codon."""
        if self._existing is None:
            context, strength = (None, None)
            if self.own_start is not None:
                context, strength = utr5.kozak(self.sequences.spliced, self.own_start)
            self._existing = self.frames.existing() | {
                'utr5_cds_kozak_context': context,
                'utr5_cds_kozak_strength': strength,
            }
        return self._existing


def transcript_columns(
    annotated: Annotated, allele: plugin_gates.Allele, evidence: frozenset[int]
) -> dict[str, Any] | None:
    """The columns of one allele on one transcript, or None where no rule evaluates the two.

    `evidence` holds the places on the transcript's chromosome of the start codons that were seen translated.
    """
    transcript = annotated.transcript
    placed = plugin_gates.Placed(transcript, allele, annotated.sequences)
    if not placed.within_transcript:
        return None
    columns: dict[str, Any] = {}
    if placed.five_prime_utr_variant:
        columns |= annotated.of_utr
        change = utr_change(transcript, allele, annotated.sequences.utr5)
        if change is not None:
            records = utr5.effects(annotated.frames, change, lambda start: transcript.genomic(start + 1) in evidence)
            if records:
                columns['utr5_consequences'] = sorted({record['consequence'] for record in records})
                columns['utr5_uorfs'] = sorted(records, key=_record_order)
    own_start = annotated.own_start
    if own_start is not None and len(placed.segments) == 1 and placed.segments[0].mapped:
        piece = placed.segments[0]
        changed = utr5.own_start_change(
            annotated.sequences.spliced,
            own_start,
            piece.start - 1,
            len(allele.ref),
            allele.alt if transcript.strand == 1 else transcript_models.reverse_complement(allele.alt),
        )
        if changed is not None:
            columns |= {key: annotated.of_utr[key] for key in ('utr5_cds_kozak_context', 'utr5_cds_kozak_strength')}
            columns |= {'utr5_cds_kozak_context_alt': changed[0], 'utr5_cds_kozak_strength_alt': changed[1]}
    ends_frame = placed.stop_gained or placed.frameshift_variant or placed.splice_site_variant
    if ends_frame and placed.cds_start is not None and placed.cds_end is not None:
        reason = nmd.escape_reason(transcript, allele.end, placed.cds_end)
        columns |= {'nmd_escaping': reason is not None, 'nmd_escape_reason': reason}
    return columns or None


@dataclasses.dataclass(frozen=True)
class _Task:
    """One chromosome of a run."""

    chrom: str
    names: Mapping[str, str]
    sites: pathlib.Path
    models: Mapping[gene_models.Source, pathlib.Path]
    reference: pathlib.Path
    evidence: frozenset[int]


def _chromosome(task: _Task) -> tuple[pa.Table, int]:
    """The rows of one chromosome's alleles, for the transcripts of every gene model file, and the count of the alleles.

    It writes a line to the log for every PROGRESS_EVERY alleles that it has read with a gene model file, and one when
    it has read them all, so that a run shows in its log where it is and how fast it goes.
    """
    rows: list[dict[str, Any]] = []
    started = time.monotonic()
    with pysam.FastaFile(str(task.reference)) as reference, pysam.VariantFile(str(task.sites)) as sites:
        if task.chrom not in sites.header.contigs:
            return TABLE.empty_table(), 0
        alleles = []
        for record in sites.fetch(task.chrom):
            if record.ref is None or record.alts is None or len(record.alts) != 1 or record.id is None:
                raise ValueError(f'{task.sites}: {record.chrom}:{record.pos} is no record of one allele under its id')
            if _ALLELE.fullmatch(record.ref) and _ALLELE.fullmatch(record.alts[0]):
                alleles.append((record.id, plugin_gates.Allele.of_record(record.pos, record.ref, record.alts[0])))
        if not alleles:
            return TABLE.empty_table(), 0
        for source, path in task.models.items():
            transcripts = transcript_models.read(path, task.chrom, task.names)
            overlapping = transcript_models.Overlapping(transcripts)
            annotated: dict[int, Annotated] = {}
            for done, (variant_id, allele) in enumerate(alleles, start=1):
                for transcript in overlapping.of(allele.low, allele.high):
                    held = annotated.get(id(transcript))
                    if held is None:
                        held = annotated[id(transcript)] = Annotated(transcript, reference)
                    columns = transcript_columns(held, allele, task.evidence)
                    if columns is not None:
                        key = {'variant_id': variant_id, 'transcript_id': transcript.transcript_id, 'source': source}
                        rows.append(key | columns)
                if done % PROGRESS_EVERY == 0 or done == len(alleles):
                    seconds = time.monotonic() - started
                    print(f'{task.chrom}, {source}: {done} of {len(alleles)} alleles after {seconds:.0f} s', flush=True)
    return pa.Table.from_pylist(rows, schema=TABLE), len(alleles)


def read_given(path: pathlib.Path, controls: Sequence[gene_models.Control]) -> list[str]:
    """The ids of the records that the step was given: the union's alleles that overlap an exon, without the controls.

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


def check_control(control: gene_models.Control, rows: Sequence[Mapping[str, Any]], where: str) -> None:
    """Fails unless the rows of a positive control hold what the control's rule must give it on its transcript.

    Gene models that were not read as prepared, a reference with other sequences or chromosome names that do not
    meet leave the control without its value, where an allele of the union would only have empty columns.
    """
    listed = [row for row in rows if row['transcript_id'] == control.transcript_id]
    if len(listed) != 1:
        raise ValueError(
            f'{where}: the positive control {control.id} has {len(listed)} rows on the transcript '
            f'{control.transcript_id}, so the gene models, the reference or the chromosome synonyms are not as prepared'
        )
    passed = (
        'uAUG_gained' in (listed[0].get('utr5_consequences') or ())
        if control.plugin == 'UTRannotator'
        else listed[0].get('nmd_escaping') is True
    )
    if not passed:
        raise ValueError(
            f'{where}: the positive control {control.id} lacks what the rules of {control.plugin} give it on '
            f'{control.transcript_id}: {dict(listed[0])}'
        )


def annotate(
    given: pathlib.Path,
    prepared: pathlib.Path,
    reference: pathlib.Path,
    uorf_evidence: pathlib.Path,
    output: pathlib.Path,
    workers: int = 1,
) -> VepPluginsInfo:
    """Writes the table of a union's alleles and the record of what its values come from into the output folder.

    `given` is the union's alleles that overlap an exon, each under its variant id, with the positive controls: a
    VCF compressed with bgzip, with its index beside it. `prepared` is the folder of the prepared gene models,
    `reference` the FASTA that the alleles are called against, with its `.fai`, and `uorf_evidence` the plugin's
    file of the upstream reading frames that were seen translated. The chromosomes are shared among `workers`
    processes. The table's rows are sorted by variant id and transcript id, so that two runs of the same inputs write
    the same table. The log gets a line for each chromosome that is done, with the alleles read and the seconds since
    the start, beside the lines of each process (_chromosome).

    Raises:
        ValueError: if an input is not what its role needs (read_given, transcript_models.read, utr5.uorf_evidence),
            a prepared gene model file is not the one that the folder's record describes, a control fails
            (check_control), or the gene model files give a transcript twice, in one file or in both.
    """
    record = gene_models.read_record(prepared)
    files: dict[gene_models.Source, FileRecord] = {}
    for described in record.gene_models:
        files[described.source] = read = _recorded(prepared / described.file)
        if (read.bytes, read.sha256) != (described.bytes, described.sha256):
            raise ValueError(
                f'{prepared / described.file} has {read.bytes} bytes and the SHA-256 {read.sha256}, which is not the '
                f'file that {gene_models.RECORD_NAME} describes'
            )
    controls = [control for file in record.gene_models for control in file.controls]
    alleles = read_given(given, controls)
    with uorf_evidence.open(encoding='utf-8') as lines:
        evidence = utr5.uorf_evidence(lines)
    names = transcript_models.sequence_names(prepared / gene_models.SYNONYMS_NAME)
    tasks = [
        _Task(
            chrom=chrom,
            names=names[chrom],
            sites=given,
            models={described.source: prepared / described.file for described in record.gene_models},
            reference=reference,
            evidence=frozenset(position for listed, position in evidence if listed == chrom),
        )
        for chrom in gene_models.CHROMOSOMES
    ]
    started = time.monotonic()
    parts: list[pa.Table] = []
    counted = 0

    def finished(part: tuple[pa.Table, int]) -> None:
        nonlocal counted
        parts.append(part[0])
        counted += part[1]
        seconds = time.monotonic() - started
        print(f'{len(parts)} of {len(tasks)} chromosomes done, {counted} alleles, after {seconds:.0f} s', flush=True)

    if workers > 1:
        context = multiprocessing.get_context('spawn')
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
            waiting = [pool.submit(_chromosome, task) for task in tasks]
            try:
                for done in concurrent.futures.as_completed(waiting):
                    finished(done.result())
            finally:
                # After a failure the chromosomes that have not begun are not read.
                for future in waiting:
                    future.cancel()
    else:
        for task in tasks:
            finished(_chromosome(task))
    table = pa.concat_tables(parts).sort_by([('variant_id', 'ascending'), ('transcript_id', 'ascending')])
    control_ids = {control.id for control in controls}
    of_controls = [variant_id in control_ids for variant_id in table['variant_id'].to_pylist()]
    control_rows = table.filter(pa.array(of_controls)).to_pylist()
    for control in controls:
        source = next(file.source for file in record.gene_models if control in file.controls)
        own = [row for row in control_rows if row['variant_id'] == control.id and row['source'] == source]
        check_control(control, own, str(given))
    table = table.filter(pa.array([not of_control for of_control in of_controls]))
    keys = zip(table['variant_id'].to_pylist(), table['transcript_id'].to_pylist(), strict=True)
    for row, following in itertools.pairwise(keys):
        if row == following:
            raise ValueError(f'two rows for {row[0]} and {row[1]}: the gene model files give the transcript twice')
    output.mkdir(parents=True)
    pq.write_table(table, output / TABLE_NAME, compression='zstd')
    info = VepPluginsInfo(
        utrannotator_commit=utr5.SOURCE_COMMIT,
        nmd_plugin_commit=nmd.SOURCE_COMMIT,
        vep_release=VEP_RELEASE,
        uorf_evidence=_recorded(uorf_evidence),
        gene_models=tuple(
            GeneModelAnnotated(source=file.source, file=files[file.source], prepared_from=file.prepared_from)
            for file in record.gene_models
        ),
    )
    (output / INFO_NAME).write_text(info.model_dump_json(indent=2) + '\n')
    print(f'{len(alleles)} alleles read, {table.num_rows} rows of an allele and transcript that a rule evaluated')
    return info


def read_info(folder: pathlib.Path) -> VepPluginsInfo:
    """The record in a folder that `annotate` wrote."""
    return VepPluginsInfo.model_validate_json((folder / INFO_NAME).read_text())
