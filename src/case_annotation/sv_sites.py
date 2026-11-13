"""The sites of a group's structural variants and copy-number segments, for one Annotator run (docs/annotation.md).

Each family brings the records of its callers without samples, QUAL or FILTER, under the callers' record IDs, and with
INFO cut down to the fields in KEPT (the workflow's SV_SITES cuts them with bcftools). The union of a group holds one
record per allele:

- a record of two alleles, such as a copy-number segment that one member loses and another gains, becomes one record
  for each allele, because the Annotator matches its sources against the first allele of a record only;
- a record that several families hold alike is written once.

A site's id is derived from its content, so that it is the same in every family and in every run, and names nothing of
a family. The map says which site each allele of each family's records is, and the union's VCF is what the Annotator
reads: it echoes a record's ID, which is how its annotation is found again (sv_annotation).

The module reads and writes the records as text. INFO END stays as the caller wrote it: the Annotator takes a
record's extent from it, and annotates a deletion without it as a single base.
"""

from __future__ import annotations

import dataclasses
import gzip
import hashlib
import pathlib
import re
from collections.abc import Iterator, Mapping, Sequence

import pyarrow as pa
import pyarrow.parquet as pq

PRIMARY = tuple(f'chr{c}' for c in (*range(1, 23), 'X', 'Y', 'M'))
_CHROM_RANK = {c: i for i, c in enumerate(PRIMARY)}

# The INFO fields of a site, in the order they are written: the record's type, extent and length, and how uncertain
# its breakpoints are. A breakend's MATEID and EVENT name records of one family's VCF, and are left out of a record
# that several families share.
KEPT = ('END', 'SVTYPE', 'SVLEN', 'CIPOS', 'CIEND', 'IMPRECISE')
_HEADER = {
    'END': ('1', 'Integer', 'End position of the variant described in this record'),
    'SVTYPE': ('1', 'String', 'Type of structural variant'),
    'SVLEN': ('.', 'Integer', 'Difference in length between REF and ALT alleles'),
    'CIPOS': ('2', 'Integer', 'Confidence interval around POS'),
    'CIEND': ('2', 'Integer', 'Confidence interval around END'),
    'IMPRECISE': ('0', 'Flag', 'Imprecise structural variation'),
}

# The types of record the Annotator is known to annotate, as INFO SVTYPE names them, and the alleles of a copy-number
# segment. A record of another type is refused here, before the Annotator runs, since the conversion of its output
# would not know what the Annotator made of it.
SV_TYPES = ('DEL', 'INS', 'DUP', 'INV', 'BND', 'CNV')
SEGMENT_ALLELES = ('<DEL>', '<DUP>')

# A header line that declares a contig with its length, which the union's header repeats for the primary chromosomes.
_CONTIG = re.compile(r'##contig=<(?:.*,)?ID=(?P<id>[^,>]+),(?:.*,)?length=(?P<length>[0-9]+)[,>].*')

MAP = pa.schema(
    [
        pa.field('family', pa.string(), nullable=False),
        pa.field('sv_id', pa.string(), nullable=False),
        pa.field('alt', pa.string(), nullable=False),
        pa.field('site_id', pa.string(), nullable=False),
    ]
)


def is_breakend(alt: str) -> bool:
    """Whether an alternate allele is a breakend, which names its mate's place in brackets and has no end."""
    return '[' in alt or ']' in alt


@dataclasses.dataclass(frozen=True)
class Site:
    """One allele of one record: what the Annotator is given of it."""

    chrom: str
    pos: int
    ref: str
    alt: str
    info: tuple[tuple[str, str | None], ...]  # the fields of KEPT that the record has, in KEPT's order

    @property
    def info_text(self) -> str:
        return ';'.join(key if value is None else f'{key}={value}' for key, value in self.info)

    @property
    def site_id(self) -> str:
        """The site's id: a digest of its content, which is all that the Annotator reads of it."""
        content = '\t'.join([self.chrom, str(self.pos), self.ref, self.alt, self.info_text])
        return hashlib.sha256(content.encode()).hexdigest()[:32]

    @property
    def sv_type(self) -> str:
        """INFO SVTYPE, which every site has."""
        return str(dict(self.info)['SVTYPE'])

    @property
    def end(self) -> int | None:
        """INFO END; None for a breakend, which has none."""
        held = dict(self.info).get('END')
        return None if held is None else int(held)

    def line(self) -> str:
        return '\t'.join([self.chrom, str(self.pos), self.site_id, self.ref, self.alt, '.', '.', self.info_text])


def _info(text: str, where: str) -> dict[str, str | None]:
    fields: dict[str, str | None] = {}
    for item in [] if text == '.' else text.split(';'):
        key, equals, value = item.partition('=')
        if key not in KEPT:
            raise ValueError(f'{where}: INFO {key}, and the sites of a family hold only {", ".join(KEPT)}')
        if key in fields:
            raise ValueError(f'{where}: INFO {key} twice')
        fields[key] = value if equals else None
    return fields


def allele_sites(columns: Sequence[str], where: str) -> list[Site]:
    """The sites of one record of a family's sites, one for each of its alleles.

    The record has the eight fixed columns of a VCF and no samples. A record of several alleles gives each allele its
    own SVLEN, and the other fields to all of them. Every record that is not a breakend needs INFO END.
    """
    if len(columns) != 8:
        raise ValueError(f"{where}: {len(columns)} columns, and a family's sites have the eight without samples")
    chrom, pos, _, ref, alts, qual, filters, info_text = columns
    if chrom not in _CHROM_RANK:
        raise ValueError(f"{where}: not on a primary chromosome; a family's sites hold only those")
    if not pos.isdigit():
        raise ValueError(f'{where}: the position {pos!r} is not a number')
    if (qual, filters) != ('.', '.'):
        raise ValueError(f"{where}: QUAL {qual} and FILTER {filters}, and a family's sites have neither")
    info = _info(info_text, where)
    if info.get('SVTYPE') is None:
        raise ValueError(f'{where}: no SVTYPE, the type of the structural variant')
    if info['SVTYPE'] not in SV_TYPES:
        raise ValueError(f'{where}: SVTYPE {info["SVTYPE"]!r} is none of {", ".join(SV_TYPES)}')
    if 'IMPRECISE' in info and info['IMPRECISE'] is not None:
        raise ValueError(f'{where}: IMPRECISE is a flag, and has the value {info["IMPRECISE"]!r}')
    alleles = alts.split(',')
    if alts == '.' or len(set(alleles)) != len(alleles) or not all(alleles):
        raise ValueError(f'{where}: the alternate alleles {alts!r} are not one or several different alleles')
    lengths = None if info.get('SVLEN') is None else str(info['SVLEN']).split(',')
    if lengths is not None and len(lengths) != len(alleles):
        raise ValueError(f'{where}: SVLEN {info["SVLEN"]} gives {len(lengths)} lengths for {len(alleles)} alleles')
    sites = []
    if 'END' in info and not (str(info['END']).isdigit() and int(str(info['END'])) >= int(pos)):
        raise ValueError(f"{where}: END {info['END']!r} is not a position at or after the record's, {pos}")
    if info['SVTYPE'] == 'CNV' and not set(alleles) <= set(SEGMENT_ALLELES):
        raise ValueError(
            f'{where}: the alleles of a copy-number segment are {" and ".join(SEGMENT_ALLELES)}, not {alts}'
        )
    for index, alt in enumerate(alleles):
        if is_breakend(alt) and 'END' in info:
            raise ValueError(f'{where}: the breakend {alt!r} has END {info["END"]}, and a breakend has no end')
        if not is_breakend(alt) and 'END' not in info:
            raise ValueError(
                f'{where}: no END, which every record but a breakend needs as its caller wrote it: the Annotator '
                'takes how far the variant reaches from it'
            )
        own = dict(info) if lengths is None else info | {'SVLEN': lengths[index]}
        sites.append(
            Site(chrom=chrom, pos=int(pos), ref=ref, alt=alt, info=tuple((key, own[key]) for key in KEPT if key in own))
        )
    return sites


def _lines(path: pathlib.Path) -> Iterator[tuple[int, str]]:
    """The lines of a gzipped VCF, each with its number."""
    with gzip.open(path, 'rb') as vcf:
        for number, line in enumerate(vcf, start=1):
            yield number, line.decode().rstrip('\n')


def read_family(path: pathlib.Path) -> tuple[dict[str, int], list[tuple[str, list[Site]]]]:
    """The primary contigs a family's sites declare with their lengths, and each record's ID with its alleles' sites."""
    contigs: dict[str, int] = {}
    records: list[tuple[str, list[Site]]] = []
    for number, line in _lines(path):
        if (declared := _CONTIG.fullmatch(line)) is not None and declared['id'] in _CHROM_RANK:
            contigs[declared['id']] = int(declared['length'])
        if line.startswith('#'):
            continue
        columns = line.split('\t')
        where = f'{path}: line {number}'
        if len(columns) < 3 or columns[2] == '.':
            raise ValueError(f'{where}: a record without an ID, which is its key in the tables')
        records.append((columns[2], allele_sites(columns, where)))
    return contigs, records


def header_lines(contigs: Mapping[str, int], alts: Sequence[str]) -> list[str]:
    """The header of the union's VCF: the primary contigs given, the fields of KEPT and the symbolic alleles."""
    symbols = sorted({alt[1:-1] for alt in alts if alt.startswith('<') and alt.endswith('>')})
    return [
        '##fileformat=VCFv4.2',
        *(f'##contig=<ID={contig},length={contigs[contig]}>' for contig in PRIMARY if contig in contigs),
        *(
            f'##INFO=<ID={key},Number={number},Type={kind},Description="{description}">'
            for key, (number, kind, description) in _HEADER.items()
        ),
        *(f'##ALT=<ID={symbol},Description="{symbol}">' for symbol in symbols),
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']),
    ]


@dataclasses.dataclass(frozen=True)
class Union:
    """What a union holds: its sites, how many of them have INFO END, and the alleles of the families' records."""

    sites: int
    with_end: int
    alleles: int


def union(families: Mapping[str, Sequence[pathlib.Path]], vcf_path: pathlib.Path, map_path: pathlib.Path) -> Union:
    """Writes the union of the families' sites as a VCF, uncompressed, and the map from their records to its sites.

    Each family gives the files of its sites, one for each of its callers. The VCF is sorted by chromosome (1 to 22,
    X, Y, M), position and id. The map has one row per family, record and allele. A record's ID that a family's files
    hold twice fails, since it is the key of the family's tables.
    """
    sites: dict[str, Site] = {}
    contigs: dict[str, int] = {}
    rows: list[dict[str, str]] = []
    for family, paths in families.items():
        ids: set[str] = set()
        for path in paths:
            declared, records = read_family(path)
            for sv_id, alleles in records:
                if sv_id in ids:
                    raise ValueError(f"{path}: the ID {sv_id} is that of several records of the family's callers")
                ids.add(sv_id)
                for site in alleles:
                    sites.setdefault(site.site_id, site)
                    rows.append({'family': family, 'sv_id': sv_id, 'alt': site.alt, 'site_id': site.site_id})
            for contig, length in declared.items():
                if contigs.setdefault(contig, length) != length:
                    raise ValueError(f'{path}: {contig} is {length} bases long, and {contigs[contig]} in another file')
    ordered = sorted(sites.values(), key=lambda site: (_CHROM_RANK[site.chrom], site.pos, site.site_id))
    if missing := sorted({site.chrom for site in ordered} - contigs.keys(), key=_CHROM_RANK.__getitem__):
        raise ValueError(f"the families' sites declare no contig {missing}, which the union's header needs")
    with vcf_path.open('w') as vcf:
        for line in header_lines(contigs, [site.alt for site in ordered]):
            vcf.write(line + '\n')
        for site in ordered:
            vcf.write(site.line() + '\n')
    pq.write_table(pa.Table.from_pylist(rows, schema=MAP), map_path, compression='zstd')
    return Union(sites=len(ordered), with_end=sum(site.end is not None for site in ordered), alleles=len(rows))


def read_union(path: pathlib.Path) -> dict[str, Site]:
    """The sites of a union's VCF by id; a record whose ID is not that of its content fails, and so does one twice."""
    sites: dict[str, Site] = {}
    for number, line in _lines(path):
        if line.startswith('#'):
            continue
        columns = line.split('\t')
        where = f'{path}: line {number}'
        found = allele_sites(columns, where)
        if len(found) != 1:
            raise ValueError(f'{where}: {len(found)} alleles, and the union holds one for each record')
        (site,) = found
        if columns[2] != site.site_id or site.site_id in sites:
            raise ValueError(f"{where}: the ID {columns[2]} is not the id of the record's content, or is there twice")
        sites[site.site_id] = site
    return sites
