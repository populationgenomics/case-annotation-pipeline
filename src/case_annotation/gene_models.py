"""The gene model files that Ensembl VEP reads without its cache, prepared from the public files (docs/annotation.md).

VEP builds its transcripts from a GFF file that is sorted, compressed with bgzip and indexed with tabix, one file per
run. `prepare` makes such a file from Ensembl's GFF3 and from NCBI's RefSeq GFF as they are published, and writes
beside them what the step that runs VEP needs of both:

- the chromosome synonyms, from NCBI's assembly report: Ensembl's file names a chromosome `1`, NCBI's `NC_000001.11`,
  and the pipeline's VCFs and reference `chr1`;
- the regions of the exons of either file, under the names of the VCFs, to which the step cuts a union's alleles
  before VEP reads them;
- the positive controls, two invented alleles for each file, which the step gives to VEP with a union's alleles and
  looks for in VEP's output before it drops them (vep_plugins);
- a record of the public files that the folder was made from.

Nothing is filtered: a prepared file holds every feature line of the public file. Its comment lines are left out, as
in the preparation that VEP's documentation gives, and the record keeps the lines that give its release.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import itertools
import pathlib
import re
from collections.abc import Iterator, Mapping, Sequence
from typing import Literal

import pydantic
import pysam

RECORD_NAME = 'gene-models.json'
EXONS_NAME = 'exons.bed'
SYNONYMS_NAME = 'chromosome-synonyms.txt'
CONTROLS_NAME = 'controls.vcf'

# The gene model files by the source of their transcripts, as the Annotator names it in `transcripts.source`, each
# with the name of its prepared file.
type Source = Literal['Ensembl', 'RefSeq']
FILES: Mapping[Source, str] = {'Ensembl': 'ensembl.gff3.gz', 'RefSeq': 'refseq.gff.gz'}

# The chromosomes of the tables (family_tables.PRIMARY), which the assembly report must name.
CHROMOSOMES = tuple(f'chr{c}' for c in (*range(1, 23), 'X', 'Y', 'M'))

# A control is an insertion a few bases into its feature, clear of the feature's ends: the 15th base of the first exon
# for the upstream start codon, and the 15th base of the coding sequence for the frameshift, which lies in the first
# 100 coding bases whatever the transcript.
_OFFSET = 14
_MARGIN = 30
# An insertion that reads ATG on the transcript's strand, which is the forward one.
_START_CODON = 'ATG'
# An id that ends in a version, to which VEP adds none.
_VERSIONED = re.compile(r'\.\d+$')

type Plugin = Literal['UTRannotator', 'NMD']


class PublicFile(pydantic.BaseModel, frozen=True, extra='forbid'):
    """A public file as the preparation read it."""

    name: str
    bytes: int
    sha256: str
    # The lines of a GFF file's head that begin with `#!`, which state its genome build and its release. None for the
    # assembly report.
    stated: tuple[str, ...]


class Control(pydantic.BaseModel, frozen=True, extra='forbid'):
    """An invented allele whose annotation is known from the gene model: a positive control of one plugin.

    The allele is an insertion on a transcript of the forward strand. For UTRannotator it puts ATG into the
    transcript's 5' UTR, which gains an upstream start codon. For NMD it puts one base into the first 100 coding bases,
    a frameshift that the plugin's rules say escapes nonsense-mediated decay.
    """

    id: str
    plugin: Plugin
    transcript_id: str  # with its version, as VEP writes it
    chrom: str
    pos: int
    ref: str
    alt: str


class GeneModelFile(pydantic.BaseModel, frozen=True, extra='forbid'):
    source: Source
    # The prepared file in the folder, which has its .tbi beside it: its name, its size and its SHA-256, by which a
    # run's record is held to the folder that it read.
    file: str
    bytes: int
    sha256: str
    prepared_from: PublicFile
    controls: tuple[Control, ...]


class GeneModels(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The record of a folder of prepared gene models (RECORD_NAME)."""

    gene_models: tuple[GeneModelFile, ...]
    assembly_report: PublicFile


@dataclasses.dataclass(frozen=True)
class Chromosome:
    """One chromosome of the assembly under its three names."""

    name: str  # as in Ensembl's files: 1, X, MT
    refseq: str  # NC_000001.11
    ucsc: str  # chr1, as in the pipeline's VCFs and reference

    @property
    def names(self) -> tuple[str, str, str]:
        return (self.name, self.ucsc, self.refseq)


@dataclasses.dataclass
class _Transcript:
    """An mRNA of the forward strand with what a control needs of it."""

    transcript_id: str
    exons: list[tuple[int, int]]
    cds: list[tuple[int, int]]


def _sha256(path: pathlib.Path) -> str:
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def _public(path: pathlib.Path, stated: Sequence[str] = ()) -> PublicFile:
    return PublicFile(name=path.name, bytes=path.stat().st_size, sha256=_sha256(path), stated=tuple(stated))


def read_assembly_report(path: pathlib.Path) -> tuple[Chromosome, ...]:
    """The assembled chromosomes of NCBI's assembly report, in the order of CHROMOSOMES.

    A line of the report holds a sequence's name, its role, and among its other columns its RefSeq accession (the
    seventh) and its UCSC name (the tenth).

    Raises:
        ValueError: if a line has fewer columns, or the assembled molecules are not exactly the chromosomes of the
            tables, each once.
    """
    found: dict[str, Chromosome] = {}
    for number, line in enumerate(path.read_text('utf-8').splitlines(), start=1):
        if line.startswith('#') or not line:
            continue
        fields = line.split('\t')
        if len(fields) < 10:
            raise ValueError(f'{path}: line {number} has {len(fields)} columns, and an assembly report has 10')
        if fields[1] != 'assembled-molecule':
            continue
        chromosome = Chromosome(name=fields[0], refseq=fields[6], ucsc=fields[9])
        if chromosome.ucsc in found:
            raise ValueError(f'{path}: {chromosome.ucsc} is listed twice')
        found[chromosome.ucsc] = chromosome
    if set(found) != set(CHROMOSOMES):
        missing, extra = sorted(set(CHROMOSOMES) - set(found)), sorted(set(found) - set(CHROMOSOMES))
        raise ValueError(f'{path}: its assembled molecules lack {missing} and hold {extra} beside chr1 to chrM')
    return tuple(found[name] for name in CHROMOSOMES)


def synonyms(chromosomes: Sequence[Chromosome]) -> str:
    """VEP's chromosome synonyms file: every pair among a chromosome's three names, since VEP does not chain them."""
    pairs = (pair for chromosome in chromosomes for pair in itertools.combinations(chromosome.names, 2))
    return ''.join(f'{left}\t{right}\n' for left, right in pairs)


def _attributes(column: str) -> dict[str, str]:
    """The attributes of a GFF3 line, each under its tag."""
    return dict(pair.split('=', 1) for pair in column.split(';') if '=' in pair)


def _blocks(lines: Iterator[str], where: str) -> Iterator[tuple[str, list[tuple[int, int, str]]]]:
    """A GFF file's feature lines in blocks of one sequence each, in the file's order: start, end and the line.

    Raises:
        ValueError: if a feature line has not nine columns or no whole numbers as its start and end, or the lines of
            one sequence are not together, which tabix needs of the sorted file.
    """
    seen: set[str] = set()
    current: str | None = None
    block: list[tuple[int, int, str]] = []
    for number, line in enumerate(lines, start=1):
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) != 9 or not (fields[3].isdigit() and fields[4].isdigit()):
            raise ValueError(f'{where}: line {number} is no GFF feature line of nine columns with a start and an end')
        if fields[0] != current:
            if current is not None:
                yield current, block
            if fields[0] in seen:
                raise ValueError(f'{where}: line {number}: the lines of {fields[0]} are not together')
            current, block = fields[0], []
            seen.add(current)
        block.append((int(fields[3]), int(fields[4]), line))
    if current is not None:
        yield current, block


def _transcripts(block: Sequence[tuple[int, int, str]]) -> list[_Transcript]:
    """The mRNAs of the forward strand in a block, in the block's order, each with its exons and coding segments.

    VEP names a transcript by the attribute `transcript_id` and, with --transcript_version, adds the attribute
    `version` to an id that does not end in one. A transcript that the file marks as partial, as an exception to its
    sequence or as without a start is left out, since a control needs a model that VEP reads plainly.
    """
    found: dict[str, _Transcript] = {}
    for start, end, line in block:
        fields = line.rstrip('\n').split('\t')
        kind, strand = fields[2], fields[6]
        if strand != '+' or kind not in ('mRNA', 'exon', 'CDS'):
            continue
        attributes = _attributes(fields[8])
        if kind == 'mRNA':
            flagged = {'partial', 'exception'} & attributes.keys() or '_NF' in attributes.get('tag', '')
            if 'ID' in attributes and 'transcript_id' in attributes and not flagged:
                name = attributes['transcript_id']
                if 'version' in attributes and not _VERSIONED.search(name):
                    name = f'{name}.{attributes["version"]}'
                found[attributes['ID']] = _Transcript(name, [], [])
        elif (parent := found.get(attributes.get('Parent', ''))) is not None:
            (parent.exons if kind == 'exon' else parent.cds).append((start, end))
    return list(found.values())


def _controls(
    source: Source, transcript: _Transcript, chromosome: Chromosome, reference: pysam.FastaFile
) -> tuple[Control, Control] | None:
    """The two controls on a transcript, or None for one that has no room for them.

    The first exon must hold _MARGIN bases of 5' UTR and the first coding segment _MARGIN bases, so that each
    insertion lies well inside its feature, and the reference must give a base at both places.
    """
    if not transcript.exons or not transcript.cds:
        return None
    if chromosome.ucsc not in reference.references:
        raise ValueError(f'the reference has no sequence named {chromosome.ucsc}; its sequences are to be chr1, ...')
    first_exon, first_cds = min(transcript.exons), min(transcript.cds)
    utr_bases = min(first_exon[1], first_cds[0] - 1) - first_exon[0] + 1
    if utr_bases < _MARGIN or first_cds[1] - first_cds[0] + 1 < _MARGIN:
        return None
    places: tuple[tuple[Plugin, int, str], ...] = (
        ('UTRannotator', first_exon[0] + _OFFSET, _START_CODON),
        ('NMD', first_cds[0] + _OFFSET, 'A'),
    )
    made = []
    for plugin, pos, inserted in places:
        base = reference.fetch(chromosome.ucsc, pos - 1, pos).upper()
        if base not in ('A', 'C', 'G', 'T'):
            return None
        made.append(
            Control(
                id=f'control-{plugin}-{source}',
                plugin=plugin,
                transcript_id=transcript.transcript_id,
                chrom=chromosome.ucsc,
                pos=pos,
                ref=base,
                alt=base + inserted,
            )
        )
    return made[0], made[1]


def _merged(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Half-open intervals without overlaps, in order: bcftools would otherwise read a record once for each."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _prepare_file(
    source: Source,
    public: pathlib.Path,
    by_name: Mapping[str, Chromosome],
    reference: pysam.FastaFile,
    exons: dict[str, list[tuple[int, int]]],
    output: pathlib.Path,
) -> GeneModelFile:
    """Writes one prepared file with its index, adds its exons to `exons`, and returns its record."""
    stated: dict[str, None] = {}
    controls: tuple[Control, Control] | None = None
    plain = output / f'{FILES[source]}.sorted'

    def noting_what_is_stated(lines: Iterator[str]) -> Iterator[str]:
        for line in lines:
            if line.startswith('#!'):
                stated[line.rstrip('\n')] = None
            yield line

    with gzip.open(public, 'rt', encoding='utf-8') as lines, plain.open('w', encoding='utf-8') as out:
        for sequence, block in _blocks(noting_what_is_stated(lines), str(public)):
            chromosome = by_name.get(sequence)
            if chromosome is not None:
                # One base beyond each end, so that an insertion written at the base before an exon is kept.
                exons.setdefault(chromosome.ucsc, []).extend(
                    (max(start - 2, 0), end + 1) for start, end, line in block if line.split('\t', 3)[2] == 'exon'
                )
                if controls is None:
                    made = (_controls(source, t, chromosome, reference) for t in _transcripts(block))
                    controls = next((pair for pair in made if pair is not None), None)
            # VEP's own preparation sorts by sequence, start and end. The sequences stay in the file's order, which
            # tabix reads as well.
            out.writelines(line for _start, _end, line in sorted(block, key=lambda feature: feature[:2]))
    if controls is None:
        raise ValueError(f'{public}: no mRNA of the forward strand on a chromosome has room for the positive controls')
    prepared = output / FILES[source]
    pysam.tabix_compress(str(plain), str(prepared), force=True)
    pysam.tabix_index(str(prepared), preset='gff', force=True)
    plain.unlink()
    return GeneModelFile(
        source=source,
        file=prepared.name,
        bytes=prepared.stat().st_size,
        sha256=_sha256(prepared),
        prepared_from=_public(public, tuple(stated)),
        controls=controls,
    )


def controls_vcf(files: Sequence[GeneModelFile]) -> str:
    """The controls of every file as the records of a VCF, each under its id."""
    records = [f'{c.chrom}\t{c.pos}\t{c.id}\t{c.ref}\t{c.alt}\t.\t.\t.\n' for file in files for c in file.controls]
    return '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' + ''.join(records)


def prepare(
    public: Mapping[Source, pathlib.Path], assembly_report: pathlib.Path, reference: pathlib.Path, output: pathlib.Path
) -> GeneModels:
    """Writes the folder of prepared gene models from the public files, and returns its record.

    `public` gives Ensembl's GFF3 and NCBI's RefSeq GFF, gzipped as published, and `reference` the FASTA that the
    pipeline's VCFs are called against, with its `.fai`, whose sequences are named `chr1`, `chr2` and so on. The
    reference gives each control the base that it is written after.

    Raises:
        ValueError: if a file is not what its role needs (read_assembly_report, _blocks), a file has no transcript
            with room for the controls on a chromosome that the report names, or the reference lacks the chromosome
            of a control.
    """
    if set(public) != set(FILES):
        raise ValueError(f'the gene model files are {sorted(FILES)}, and {sorted(public)} were given')
    chromosomes = read_assembly_report(assembly_report)
    by_name = {name: chromosome for chromosome in chromosomes for name in (chromosome.name, chromosome.refseq)}
    output.mkdir(parents=True)
    exons: dict[str, list[tuple[int, int]]] = {}
    with pysam.FastaFile(str(reference)) as fasta:
        files = [_prepare_file(source, public[source], by_name, fasta, exons, output) for source in FILES]
    (output / SYNONYMS_NAME).write_text(synonyms(chromosomes))
    (output / EXONS_NAME).write_text(
        ''.join(
            f'{chromosome.ucsc}\t{start}\t{end}\n'
            for chromosome in chromosomes
            for start, end in _merged(exons.get(chromosome.ucsc, []))
        )
    )
    (output / CONTROLS_NAME).write_text(controls_vcf(files))
    record = GeneModels(gene_models=tuple(files), assembly_report=_public(assembly_report))
    (output / RECORD_NAME).write_text(record.model_dump_json(indent=2) + '\n')
    return record


def read_record(folder: pathlib.Path) -> GeneModels:
    """The record of a folder that `prepare` wrote."""
    return GeneModels.model_validate_json((folder / RECORD_NAME).read_text())
