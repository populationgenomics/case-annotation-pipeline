"""echtvar archives of AlphaGenome's precomputed scores: what builds one and checks it (docs/annotation.md).

DeepMind publishes AlphaGenome's scores of every single-base change as three tab-separated tables, each sorted by
position, with the columns `#CHROM POS REF ALT` and then one column per value. `echtvar encode` reads a VCF, so `to_vcf`
writes a table's rows as VCF records with one INFO field per value column, under the column's name. The encode
configuration of each table is a file of this package (`config`), which gives each field its alias in the archive.
`check` compares what `echtvar anno` returns from an archive with the rows that the archive was built from.

The conversion passes a value on as the table writes it, so the archive holds what echtvar makes of that text.
"""

from __future__ import annotations

import dataclasses
import functools
import importlib.resources
import itertools
import json
import pathlib
import re
from collections.abc import Iterable, Iterator, Mapping
from typing import Any, TextIO

import pysam

from case_annotation import sv_sites

# The notice that DeepMind's terms give as an example for output under its Output Terms of Use. The encode
# configurations of the two tables that are not for permissive use hold it in every field's description, which echtvar
# embeds in the archive and writes into the header of a VCF that it annotates.
NOTICE = (
    'By using this information, you agree to AlphaGenome Output Terms of Use found at '
    'https://deepmind.google.com/science/alphagenome/output-terms'
)
# The stand-in of a float field with which echtvar writes a missing value (`.`) for an allele that the archive lacks:
# the bits of htslib's missing float. Any other stand-in is a number, which a score could be.
FLOAT_MISSING = 0x7F800001
# echtvar stores a float as the whole part of its product with the field's multiplier, in 32 bits.
_STORED_BELOW = 2**31

_KEY = ('#CHROM', 'POS', 'REF', 'ALT')
# A value as the tables write it: a decimal number, in scientific notation for some small values.
_NUMBER = re.compile(r'-?[0-9]+(?:\.[0-9]+)?(?:e[-+]?[0-9]+)?')
_ZERO = r'-?0+(?:\.0+)?(?:e[-+]?[0-9]+)?'


@dataclasses.dataclass(frozen=True)
class Table:
    """One of DeepMind's tables, by the name of the archive that is built from it.

    `columns` are the value columns of its header, after the four of the allele and in the header's order. The columns
    in `zero` hold 0 in every row and are left out of the archive: a row with another value there fails.
    """

    name: str
    title: str
    columns: tuple[str, ...]
    zero: frozenset[str] = frozenset()

    @property
    def kept(self) -> tuple[str, ...]:
        return tuple(column for column in self.columns if column not in self.zero)


_BREAKDOWN = (
    'MERGED_SPLICING',
    'MAX_ABS_ATAC',
    'MAX_ABS_CONTACT_MAPS',
    'MAX_ABS_DNASE',
    'MAX_ABS_CHIP_TF',
    'MAX_ABS_CHIP_HISTONE',
    'MAX_ABS_CAGE',
    'MAX_ABS_PROCAP',
    'MAX_ABS_RNA_SEQ',
    'MAX_ABS_POLYADENYLATION',
    'ALPHAMISSENSE',
    'CACTUS_241_WAY',
    'PROTEIN_TERMINATION',
    'START_LOST',
    'STOP_LOST',
    'PHASTCONS_470_WAY',
)
# The breakdown's two columns for insertions and deletions, which a table of single-base changes holds as 0.
_INDEL = ('IS_INSERTION', 'IS_DELETION')

TABLES: Mapping[str, Table] = {
    table.name: table
    for table in (
        Table('alphagenome_avi', 'AVI SNV scores', ('raw_score', 'PHRED')),
        Table(
            'alphagenome_avi_breakdown',
            'AVI SNV feature importance scores',
            (*_BREAKDOWN, *_INDEL),
            frozenset(_INDEL),
        ),
        Table('alphagenome_splicing', 'SNV merged splicing scores', ('alphagenome_splicing',)),
    )
}


def config_text(table: Table) -> str:
    """The table's encode configuration as `echtvar encode` takes it, which is a file of this package."""
    return importlib.resources.files('case_annotation').joinpath('echtvar', f'{table.name}.json').read_text('utf-8')


def config(table: Table) -> dict[str, dict[str, Any]]:
    """The fields of the table's encode configuration, by the column that each encodes."""
    return {field['field']: field for field in json.loads(config_text(table))}


def _storable(multiplier: int) -> str:
    """A pattern for the numbers that a field of this multiplier stores: few enough digits before the point.

    A value of at most that many whole digits, written without a positive exponent, has a product with the multiplier
    that fits in the 32 bits echtvar stores it in. echtvar would store a larger one as another number.
    """
    digits = len(str((_STORED_BELOW - 1) // multiplier)) - 1
    if digits < 1:
        raise ValueError(f'the multiplier {multiplier} leaves no value that echtvar can store')
    return rf'-?[0-9]{{1,{digits}}}(?:\.[0-9]+)?(?:e-[0-9]+)?'


@functools.cache
def _row_pattern(table: Table) -> re.Pattern[str]:
    """What a row of the table matches, as one pattern: its groups are the allele's four cells and the kept values."""
    fields = config(table)
    values = ''.join(
        rf'\t(?:{_ZERO})' if column in table.zero else rf'\t({_storable(fields[column]["multiplier"])})'
        for column in table.columns
    )
    return re.compile(rf'([^\t\n]+)\t([1-9][0-9]*)\t([ACGT])\t([ACGT]){values}\n?')


def _problem(table: Table, line: str) -> str:
    """Why a line is not a row of the table."""
    cells = line.removesuffix('\n').split('\t')
    if len(cells) != len(_KEY) + len(table.columns):
        return f'{len(cells)} cells, and the header has {len(_KEY) + len(table.columns)}'
    chrom, pos, ref, alt, *values = cells
    if not chrom or not re.fullmatch(r'[1-9][0-9]*', pos):
        return f'no chromosome and position in {chrom!r} and {pos!r}'
    if not re.fullmatch(r'[ACGT]\t[ACGT]', f'{ref}\t{alt}'):
        return f'the alleles {ref!r} and {alt!r} are not two single bases'
    fields = config(table)
    for column, value in zip(table.columns, values, strict=True):
        if not _NUMBER.fullmatch(value):
            return f'{column} is {value!r}, not a number'
        if column in table.zero:
            if not re.fullmatch(_ZERO, value):
                return f'{column} is {value}, not 0'
        elif not re.fullmatch(_storable(fields[column]['multiplier']), value):
            return f'{column} is {value}, which echtvar cannot store with the multiplier {fields[column]["multiplier"]}'
    return 'not a row of the table'


def rows(table: Table, lines: Iterable[str]) -> Iterator[tuple[str, ...]]:
    """Each row of a table's text as its chromosome, position, REF, ALT and kept values, all as the table writes them.

    The first line is the table's header.

    Raises:
        ValueError: if the header is not the table's, a row is malformed, or a row holds another value than 0 in a
            column that the archive leaves out.
    """
    lines = iter(lines)
    header = next(lines, None)
    expected = '\t'.join((*_KEY, *table.columns))
    if header is None or header.removesuffix('\n') != expected:
        found = 'no line' if header is None else repr(header.removesuffix('\n'))
        raise ValueError(f'the header of {table.name} is {expected!r}, and the table starts with {found}')
    pattern = _row_pattern(table)
    for number, line in enumerate(lines, start=1):
        found = pattern.fullmatch(line)
        if found is None:
            raise ValueError(f'row {number} of {table.name}: {_problem(table, line)}: {line.removesuffix(chr(10))!r}')
        yield found.groups()


def read_contigs(reference_index: pathlib.Path) -> dict[str, int]:
    """The primary chromosomes of a reference with their lengths, from its `.fai`, in the tables' order of them.

    Raises:
        ValueError: if the index holds none, as that of a reference with contigs named 1, 2, ... does.
    """
    lengths = {}
    for line in reference_index.read_text('utf-8').splitlines():
        name, length, *_ = line.split('\t')
        lengths[name] = int(length)
    contigs = {name: lengths[name] for name in sv_sites.PRIMARY if name in lengths}
    if not contigs:
        raise ValueError(f'{reference_index}: no primary chromosome; are its contigs named chr1, chr2, ...?')
    return contigs


def to_vcf(table: Table, lines: Iterable[str], contigs: Mapping[str, int], out: TextIO) -> int:
    """Writes a table's rows as the VCF that `echtvar encode` reads with the table's configuration; returns the rows.

    Each kept column is an INFO field of one float under the column's name, with the description that the
    configuration gives its field where it gives one. Every record holds every field, so that a value that echtvar
    returns as missing says that the archive lacks the allele. The header declares the contigs given, which are those of
    every file that one archive is encoded from, since echtvar takes them from the first file.

    Raises:
        ValueError: as `rows` does, and if a row lies on another contig than the ones given or the table has no row.
    """
    fields = config(table)
    # The first row is read before anything is written, so that a table with another header writes nothing.
    cells = rows(table, lines)
    first = next(cells, None)
    if first is None:
        raise ValueError(f'{table.name}: the table has a header and no row')
    out.write('##fileformat=VCFv4.2\n')
    out.write(f"##source=AlphaGenome's {table.title}\n")
    out.writelines(f'##contig=<ID={name},length={length}>\n' for name, length in contigs.items())
    for column in table.kept:
        description = fields[column].get('description', f"{column} of AlphaGenome's {table.title}")
        out.write(f'##INFO=<ID={column},Number=1,Type=Float,Description="{description}">\n')
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    record = '%s\t%s\t.\t%s\t%s\t.\t.\t' + ';'.join(f'{column}=%s' for column in table.kept) + '\n'
    written = 0
    chromosome = None
    for row in itertools.chain((first,), cells):
        if row[0] != chromosome:
            chromosome = row[0]
            if chromosome not in contigs:
                raise ValueError(f'row {written + 1} of {table.name}: {chromosome} is none of {", ".join(contigs)}')
        out.write(record % row)
        written += 1
    return written


def indexed_lines(tsv: pathlib.Path, chromosome: str | None) -> Iterator[str]:
    """The header and the rows of a bgzipped table with a tabix index beside it: one chromosome's rows, or all.

    Raises:
        ValueError: if the table has no row on the chromosome.
    """
    with pysam.TabixFile(str(tsv)) as table:
        if chromosome is not None and chromosome not in table.contigs:
            raise ValueError(f'{tsv}: no row on {chromosome}; it has rows on {", ".join(table.contigs)}')
        yield from table.header
        yield from table.fetch(chromosome) if chromosome is not None else table.fetch()


def _tolerance(value: float, multiplier: int) -> float:
    """How far a value that echtvar returns may lie from the one encoded.

    echtvar cuts the product of the value and the multiplier to a whole number, which loses up to one step of
    1/multiplier, and it multiplies, stores and divides as 32-bit floats, each of which is exact to 24 bits.
    """
    return 1 / multiplier + abs(value) * 2**-22


def check(table: Table, lines: Iterable[str], annotated: pathlib.Path) -> int:
    """Compares rows of a table with what `echtvar anno` returned for them from the table's archive; returns the rows.

    `lines` is the text of the rows with the table's header, and `annotated` a VCF of the same alleles in the same
    order, as `to_vcf` writes them, after echtvar annotated it with the archive.

    Raises:
        ValueError: if the alleles differ, or the archive returned no value or another value than a row's for a field.
    """
    fields = config(table)
    problems: list[str] = []
    compared = 0
    with pysam.VariantFile(str(annotated)) as vcf:
        for column in table.kept:
            if fields[column]['alias'] not in vcf.header.info:
                raise ValueError(f'{annotated}: no INFO field {fields[column]["alias"]}, which the archive adds')
        for cells, record in itertools.zip_longest(rows(table, lines), vcf):
            if cells is None or record is None:
                more = 'the annotated VCF' if cells is None else 'the table'
                raise ValueError(f'{annotated}: {more} holds alleles after the {compared} that both hold')
            chrom, pos, ref, alt, *values = cells
            allele = f'{chrom}-{pos}-{ref}-{alt}'
            if (record.chrom, str(record.pos), record.ref, record.alts) != (chrom, pos, ref, (alt,)):
                raise ValueError(f'{annotated}: record {compared + 1} is not {allele}, the row at its place')
            for column, text in zip(table.kept, values, strict=True):
                field = fields[column]
                returned: Any = record.info.get(field['alias'])
                if isinstance(returned, tuple):
                    (returned,) = returned
                if returned is None:
                    problems.append(f'{allele}: no {field["alias"]}, and {column} is {text}')
                elif abs(returned - float(text)) > _tolerance(float(text), field['multiplier']):
                    problems.append(f'{allele}: {field["alias"]} is {returned!r}, and {column} is {text}')
            compared += 1
            if len(problems) >= 5:
                break
    if problems:
        raise ValueError(f'{annotated}: the archive does not return the rows of {table.name}: {"; ".join(problems)}')
    if not compared:
        raise ValueError(f'{table.name}: the table has a header and no row')
    return compared
