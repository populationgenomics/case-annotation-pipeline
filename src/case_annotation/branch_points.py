"""The catalogue of intronic branch points on GRCh38, prepared from the published datasets (docs/annotation.md).

Zhang et al. (PNAS 2022, doi:10.1073/pnas.2211194119) published their catalogue of human branch points as two text
files on GRCh37. Dataset S1 lists every branch point with the sources that support it. Dataset S2 gives, for each
branch point that lies 3 to 100 bases before the acceptor of an intron of their gene models, one row for each such
intron, with the distance to the acceptor, the branch point's rank in the intron and the class of its motif.

`prepare` lifts the branch points to GRCh38 with a chain file and writes the one file that the tables read
(joined_values): a row for each row of Dataset S2 whose branch point lifts cleanly, and a row with the sources alone
for each branch point that only Dataset S1 lists and that lifts cleanly. No branch point is selected by a rule of this
package. A record of the datasets, of the chain file and of the counts is written beside the file, and the file's first
lines name the three files by their SHA-256, so that a catalogue says how it was lifted.

A branch point is one base, and the tables match an allele to the three bases at positions -2 to 0 of it, on its
strand. Those three bases are therefore what must lift: each is lifted by itself, and the three must come out next to
one another, in the same order along the strand, on the chromosome that they were on. A chain file can turn a stretch
around, so a branch point can come out on the other strand, and it is kept there.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import importlib.metadata
import pathlib
import re
from collections.abc import Iterator, Mapping, Sequence
from typing import Literal

import pydantic
import pyliftover
import pysam

from case_annotation import sv_sites

FILE_NAME = 'branch-points.tsv.gz'
RECORD_NAME = 'branch-points.json'

# The columns of the prepared file that are of a row of Dataset S2, and empty in the row of a branch point that only
# Dataset S1 lists.
INTRON_COLUMNS = (
    'gene',
    'transcript_introns',
    'intron_type',
    'intron_length',
    'acceptor_offset',
    'rank',
    'total',
    'consensus',
)
# The columns of the prepared file. `name` is the paper's name of the branch point with its GRCh38 position, and
# `grch37_name` the name as the datasets write it. Every row fills the last two, so that no line ends in empty fields,
# which a tool that trims the ends of lines would remove.
COLUMNS = ('chrom', 'pos', 'strand', 'name', *INTRON_COLUMNS, 'sources', 'grch37_name')
# A line of the prepared file before its header starts with this.
PREAMBLE = '##'
# The classes of a branch point's motif and the kinds of intron, as Dataset S2 writes them.
CONSENSUS = ('1:YTNAY', '2:YTNA', '3:TNA', '4:YNA', '0:none')
INTRON_TYPES = ('major', 'minor')
# A branch point's name: its kind (m for a branch point with an experimental and a predicted source, e for one with
# experimental sources alone, c for one with predicted sources alone), its chromosome, position and strand, and its
# base on that strand. A source is named after its kind, e or c.
NAME = r'([mec])BP_(chr[0-9A-Za-z]+)_([1-9][0-9]*)_([+-])_([ACGT])'
SOURCE = r'[ec]BP_[A-Za-z0-9]+'
# A transcript's intron in the paper's gene models: Ensembl's transcript id without its version and the intron's number.
TRANSCRIPT_INTRON = r'ENSTR?[0-9]+_IVS[1-9][0-9]*'

_S1_HEADER = ('CHROM', 'BP_POS', 'STRAND', 'BP_NAME', 'SOURCE')
_S2_HEADER = (
    'CHROM',
    'POS',
    'STRAND',
    'BP_NAME',
    'SOURCE_COUNT',
    'SOURCE_LIST',
    'GENE',
    'TRANSCRIPT_INTRON',
    'INTRON_LENGTH',
    'INTRON_TYPE',
    'DIST_3SS',
    'RANK',
    'CONSENSUS',
    'ENERGY',
    'BP_MAF',
    'BP_GERP',
    'BP_PHYLOP',
    'BP-2_MAF',
    'BP-2_GERP',
    'BP-2_PHYLOP',
)
_NAME = re.compile(NAME)
_SOURCES = re.compile(f'{SOURCE}(,{SOURCE})*')
_TRANSCRIPT_INTRONS = re.compile(f'{TRANSCRIPT_INTRON}(,{TRANSCRIPT_INTRON})*')
_RANK = re.compile(r'#([1-9][0-9]*)/([1-9][0-9]*)')
_WHOLE = re.compile(r'[1-9][0-9]*')
_KINDS: Mapping[str, frozenset[str]] = {'m': frozenset('ec'), 'e': frozenset('e'), 'c': frozenset('c')}
_CHROM_RANK = {chrom: rank for rank, chrom in enumerate(sv_sites.PRIMARY)}
# The positions of a branch point's motif that an allele is matched to, along its strand.
_WINDOW = (-2, -1, 0)

# What became of a branch point in the lift. The first two are written.
type Outcome = Literal['clean', 'strand_flipped', 'unmapped', 'several_places', 'other_chromosome', 'split_window']


class InputFile(pydantic.BaseModel, frozen=True, extra='forbid'):
    """A file as the preparation read it."""

    name: str
    bytes: int
    sha256: str


class Lifted(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The branch points of Dataset S1 by what the chain file made of them."""

    # The three bases came out together on the same chromosome and strand.
    clean: int
    # They came out together on the same chromosome and the other strand, where the branch point is written.
    strand_flipped: int
    # The chain file has no place for the branch point's base.
    unmapped: int
    # It has several.
    several_places: int
    # Its place is on another sequence, such as an alternate contig.
    other_chromosome: int
    # One of the two bases before it has no place, several, or one that is not next to it.
    split_window: int


class Catalogue(pydantic.BaseModel, frozen=True, extra='forbid'):
    """The record of a prepared catalogue (RECORD_NAME)."""

    # The prepared file: its name, its size and its SHA-256, by which a family's manifest names it.
    file: str
    bytes: int
    sha256: str
    dataset_s1: InputFile
    dataset_s2: InputFile
    chain: InputFile
    pyliftover: str
    # The branch points of Dataset S1, and those among them that Dataset S2 has a row for.
    branch_points: int
    branch_points_with_introns: int
    lifted: Lifted
    # What the prepared file holds: its branch points, those among them with a row of Dataset S2, and its rows.
    written_branch_points: int
    written_branch_points_with_introns: int
    written_rows: int


@dataclasses.dataclass(frozen=True)
class _BranchPoint:
    chrom: str
    pos: int
    strand: str
    name: str
    kind: str
    base: str
    sources: str


@dataclasses.dataclass(frozen=True)
class _Intron:
    """What a row of Dataset S2 says of a branch point in one intron, as the prepared file writes INTRON_COLUMNS."""

    gene: str
    transcript_introns: str
    intron_type: str
    intron_length: str
    acceptor_offset: str
    rank: str
    total: str
    consensus: str


def _gzipped(path: pathlib.Path) -> bool:
    with path.open('rb') as f:
        return f.read(2) == b'\x1f\x8b'


def _input(path: pathlib.Path) -> InputFile:
    with path.open('rb') as f:
        return InputFile(name=path.name, bytes=path.stat().st_size, sha256=hashlib.file_digest(f, 'sha256').hexdigest())


def _rows(path: pathlib.Path, header: Sequence[str], what: str) -> Iterator[tuple[int, list[str]]]:
    """The rows of a tab-separated dataset after its header, gzipped or not, each with its line's number.

    Raises:
        ValueError: if the file cannot be read as text, its header is not the dataset's, or a row has another number
            of fields.
    """
    try:
        with gzip.open(path, 'rt', encoding='utf-8') if _gzipped(path) else path.open(encoding='utf-8') as lines:
            found = lines.readline().rstrip('\r\n').split('\t')
            if found != list(header):
                raise ValueError(f'{path}: its header is {found}, and that of {what} is {list(header)}')
            for number, line in enumerate(lines, start=2):
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != len(header):
                    raise ValueError(f'{path}: line {number} has {len(fields)} fields, and {what} has {len(header)}')
                yield number, fields
    except (OSError, EOFError, UnicodeDecodeError) as unreadable:
        raise ValueError(f'{path}: not {what}: {unreadable}') from unreadable


def _named(where: str, chrom: str, pos: str, strand: str, name: str) -> re.Match[str]:
    """A branch point's name, taken apart.

    Raises:
        ValueError: if the name is not of the paper's form, the row's chromosome, position or strand is not the
            name's, or the chromosome is none of the tables.
    """
    match = _NAME.fullmatch(name)
    if match is None or (match[2], match[3], match[4]) != (chrom, pos, strand):
        raise ValueError(
            f'{where}: the name {name!r} is not a kind, {chrom}, {pos}, {strand} and a base, as in mBP_chr1_100_+_A'
        )
    if chrom not in _CHROM_RANK:
        raise ValueError(f'{where}: {chrom} is none of chr1 to chr22, chrX, chrY and chrM')
    return match


def read_s1(path: pathlib.Path) -> dict[str, _BranchPoint]:
    """The branch points of Dataset S1 by name, in the file's order.

    Raises:
        ValueError: if the file is not Dataset S1 (_rows, _named), a row's sources are no list of names such as
            eBP_Mercer or are not of the kinds that the name's first letter states, or two rows hold one branch point
            or one position of a strand.
    """
    found: dict[str, _BranchPoint] = {}
    places: set[tuple[str, int, str]] = set()
    for number, (chrom, pos, strand, name, sources) in _rows(path, _S1_HEADER, 'Dataset S1'):
        where = f'{path}: line {number}'
        match = _named(where, chrom, pos, strand, name)
        if _SOURCES.fullmatch(sources) is None:
            raise ValueError(f'{where}: the sources {sources!r} are no names such as eBP_Mercer with commas between')
        if {source[0] for source in sources.split(',')} != _KINDS[match[1]]:
            raise ValueError(f'{where}: the name {name!r} states other kinds of source than {sources!r}')
        place = (chrom, int(pos), strand)
        if name in found or place in places:
            raise ValueError(f'{where}: a second row of {chrom}:{pos} on the {strand} strand')
        places.add(place)
        found[name] = _BranchPoint(chrom, int(pos), strand, name, match[1], match[5], sources)
    return found


def read_s2(path: pathlib.Path, branch_points: Mapping[str, _BranchPoint]) -> dict[str, list[_Intron]]:
    """The rows of Dataset S2 by the name of their branch point, in the file's order.

    Raises:
        ValueError: if the file is not Dataset S2 (_rows, _named), a row is of a branch point that Dataset S1 lacks or
            gives it other sources or another number of them, a row's gene is empty, its introns are no list of labels
            such as ENST00000000001_IVS2, its intron's length is no whole number, its intron's type, its distance to
            the acceptor, its rank or its consensus is not what the dataset writes, or two rows hold one branch point
            in the same introns.
    """
    found: dict[str, list[_Intron]] = {}
    seen: set[tuple[str, str]] = set()
    for number, fields in _rows(path, _S2_HEADER, 'Dataset S2'):
        chrom, pos, strand, name, count, sources, gene, introns, length, kind, distance, rank, consensus = fields[:13]
        where = f'{path}: line {number}'
        _named(where, chrom, pos, strand, name)
        if name not in branch_points:
            raise ValueError(f'{where}: Dataset S1 has no branch point {name}')
        if sources != branch_points[name].sources or count != str(len(sources.split(','))):
            raise ValueError(
                f'{where}: {count} sources, {sources!r}, and Dataset S1 gives {name} {branch_points[name].sources!r}'
            )
        ranked = _RANK.fullmatch(rank)
        problems = {
            'a gene': not gene,
            'introns such as ENST00000000001_IVS2 with commas between': _TRANSCRIPT_INTRONS.fullmatch(introns) is None,
            "a whole number as its intron's length": _WHOLE.fullmatch(length) is None,
            f'one of {", ".join(INTRON_TYPES)} as its type of intron': kind not in INTRON_TYPES,
            'a negative whole number as its distance to the acceptor': _WHOLE.fullmatch(distance[1:]) is None
            or distance[0] != '-',
            'a rank such as #1/2 that is at most its total': ranked is None or int(ranked[1]) > int(ranked[2]),
            f'one of {", ".join(CONSENSUS)} as its consensus': consensus not in CONSENSUS,
        }
        lacking = [what for what, lacks in problems.items() if lacks]
        if ranked is None or lacking:
            raise ValueError(f'{where}: the row of {name} is without {"; ".join(lacking)}')
        if (name, introns) in seen:
            raise ValueError(f'{where}: a second row of {name} in {introns}')
        seen.add((name, introns))
        found.setdefault(name, []).append(
            _Intron(gene, introns, kind, length, distance, ranked[1], ranked[2], consensus)
        )
    return found


class _Chain:
    """A chain file, which takes a position of one assembly to its places in another."""

    def __init__(self, path: pathlib.Path) -> None:
        """Reads the chain file, gzipped or not.

        Raises:
            ValueError: if the file is no chain file.
        """
        self._path = path
        try:
            with gzip.open(path, 'rb') if _gzipped(path) else path.open('rb') as f:
                self._chains = pyliftover.LiftOver(f)
        except Exception as malformed:
            raise ValueError(f'{path}: not a chain file: {malformed}') from malformed
        self._names: dict[str, str] = {}

    def places(self, chrom: str, pos: int, strand: str) -> list[tuple[str, int, str]]:
        """Every place of a 1-based position and its strand in the other assembly: sequence, position and strand.

        The chain file names a chromosome as UCSC's files do, chr1, or as Ensembl's do, 1, and a place is given under
        the first of the two.

        Raises:
            ValueError: if the chain file has no chain from the chromosome under either name.
        """
        names = (chrom, chrom.removeprefix('chr'))
        for name in (self._names[chrom],) if chrom in self._names else names:
            found = self._chains.convert_coordinate(name, pos - 1, strand)
            if found is not None:
                self._names[chrom] = name
                return [(to if to.startswith('chr') else f'chr{to}', at + 1, on) for to, at, on, _score in found]
        raise ValueError(f'{self._path}: no chain from {" or ".join(names)}, which the datasets hold')

    def lift(self, branch_point: _BranchPoint) -> tuple[Outcome, tuple[str, int, str] | None]:
        """What becomes of a branch point, and its place where it is written."""
        step = 1 if branch_point.strand == '+' else -1
        chrom, strand = branch_point.chrom, branch_point.strand
        own, *before = (self.places(chrom, branch_point.pos + offset * step, strand) for offset in reversed(_WINDOW))
        if not own:
            return 'unmapped', None
        if len(own) > 1:
            return 'several_places', None
        ((to_chrom, to_pos, to_strand),) = own
        if to_chrom != chrom:
            return 'other_chromosome', None
        to_step = 1 if to_strand == '+' else -1
        if any(places != [(to_chrom, to_pos - away * to_step, to_strand)] for away, places in enumerate(before, 1)):
            return 'split_window', None
        return ('clean' if to_strand == strand else 'strand_flipped'), (to_chrom, to_pos, to_strand)


def prepare(dataset_s1: pathlib.Path, dataset_s2: pathlib.Path, chain: pathlib.Path, output: pathlib.Path) -> Catalogue:
    """Writes the folder of the prepared catalogue from the two datasets and a chain file, and returns its record.

    The datasets are the text files of the paper's supplement, gzipped or not, and the chain file one from GRCh37 to
    GRCh38, gzipped or not, whose chromosomes are named as UCSC or as Ensembl names them. The prepared file is
    compressed with bgzip and its rows are sorted by chromosome, position and strand, so the same three files give the
    same bytes.

    Raises:
        ValueError: if a file is not what its role needs (read_s1, read_s2, _Chain), or the chain file has no chain
            from a chromosome of the datasets.
    """
    branch_points = read_s1(dataset_s1)
    introns = read_s2(dataset_s2, branch_points)
    lifter = _Chain(chain)
    counts: dict[str, int] = dict.fromkeys(Lifted.model_fields, 0)
    rows: list[tuple[int, int, str, str, str, str]] = []
    written = with_introns = 0
    for name, branch_point in branch_points.items():
        outcome, place = lifter.lift(branch_point)
        counts[outcome] += 1
        if place is None:
            continue
        chrom, pos, strand = place
        placed = [chrom, str(pos), strand, f'{branch_point.kind}BP_{chrom}_{pos}_{strand}_{branch_point.base}']
        written += 1
        with_introns += name in introns
        for intron in introns.get(name) or [None]:
            described = [''] * len(INTRON_COLUMNS) if intron is None else list(dataclasses.astuple(intron))
            key = '' if intron is None else intron.transcript_introns
            line = '\t'.join([*placed, *described, branch_point.sources, name])
            rows.append((_CHROM_RANK[chrom], pos, strand, name, key, line))
    # + sorts before - as a character, and a branch point's rows by their introns.
    rows.sort(key=lambda row: row[:5])
    inputs = {'dataset_s1': _input(dataset_s1), 'dataset_s2': _input(dataset_s2), 'chain': _input(chain)}
    version = importlib.metadata.version('pyliftover')
    output.mkdir(parents=True)
    plain = output / f'{FILE_NAME}.plain'
    with plain.open('w', encoding='utf-8') as out:
        out.write(
            f'{PREAMBLE}branch points of Zhang et al., PNAS 2022 (doi:10.1073/pnas.2211194119), Datasets S1 and S2, '
            'lifted from GRCh37 to GRCh38\n'
        )
        out.writelines(f'{PREAMBLE}{role}: {file.name} sha256 {file.sha256}\n' for role, file in inputs.items())
        out.write(f'{PREAMBLE}pyliftover: {version}\n')
        out.write('\t'.join(COLUMNS) + '\n')
        out.writelines(f'{row[5]}\n' for row in rows)
    prepared = output / FILE_NAME
    pysam.tabix_compress(str(plain), str(prepared), force=True)
    plain.unlink()
    written_file = _input(prepared)
    record = Catalogue(
        file=written_file.name,
        bytes=written_file.bytes,
        sha256=written_file.sha256,
        **inputs,
        pyliftover=version,
        branch_points=len(branch_points),
        branch_points_with_introns=len(introns),
        lifted=Lifted(**counts),
        written_branch_points=written,
        written_branch_points_with_introns=with_introns,
        written_rows=len(rows),
    )
    (output / RECORD_NAME).write_text(record.model_dump_json(indent=2) + '\n')
    return record


def read_record(folder: pathlib.Path) -> Catalogue:
    """The record of a folder that `prepare` wrote."""
    return Catalogue.model_validate_json((folder / RECORD_NAME).read_text())
