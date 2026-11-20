"""Invented output of `bcftools roh` and of `upd regions`, and the folder that roh_upd.convert writes from it.

Every run, region, frequency and version here is invented. The lines are in the layout that the two tools write.
"""

from __future__ import annotations

import json
import pathlib
import zipfile
from collections.abc import Mapping, Sequence
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from case_annotation import echtvar_values, roh_upd

# The fields of an archive of gnomAD's joint release, as echtvar embeds them.
CONFIG: list[dict[str, Any]] = [
    {'field': 'AC_joint', 'alias': 'gnomad_AC', 'missing_value': -2147483648, 'multiplier': 1, 'ftype': 'Integer'},
    {'field': 'AF_joint', 'alias': 'gnomad_AF', 'missing_value': 0, 'multiplier': 2_000_000, 'ftype': 'Float'},
    {
        'field': 'nhomalt_joint',
        'alias': 'gnomad_HomAlt',
        'missing_value': -2147483648,
        'multiplier': 1,
        'ftype': 'Integer',
    },
    {
        'field': 'fafmax_faf95_max_joint',
        'alias': 'gnomad_faf95max',
        'missing_value': 0,
        'multiplier': 2_000_000,
        'ftype': 'Float',
    },
]
# What the step passes of the tools that ran.
RAN: Mapping[str, Any] = {
    'bcftools_version': '0.0',
    'echtvar_version': '0.0.0',
    'upd': ('0.0', '0' * 40),
    'preparation': ('bcftools view --invented', 'echtvar anno'),
    'roh_command': 'bcftools roh --AF-tag gnomad_AF -G30 -O r',
    'upd_command': 'upd --af-tag gnomad_AF regions',
}
HEADER = (
    '# This file was produced by: bcftools roh(0.0+htslib-0.0)\n'
    '# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality\n'
)


def archive(path: pathlib.Path, **changed: object) -> pathlib.Path:
    """An archive whose field of the frequency has these values in place of CONFIG's."""
    config = [field | changed if field['field'] == 'AF_joint' else field for field in CONFIG]
    with zipfile.ZipFile(path, 'w') as z:
        z.writestr('echtvar/config.json', json.dumps([field | {'missing_string': 'MISSING'} for field in config]))
    return path


def run(member: str, chrom: str, pos: int, end_pos: int, markers: int = 40, quality: float = 62.5) -> str:
    """A line of `bcftools roh -O r`."""
    return f'RG\t{member}\t{chrom}\t{pos}\t{end_pos}\t{end_pos - pos + 1}\t{markers}\t{quality}'


def region(
    chrom: str,
    pos: int,
    end_pos: int,
    origin: str = 'MATERNAL',
    outer: tuple[int, int] | None = None,
    sites: int = 12,
    kind: str = 'HETERODISOMY',
    **changed: object,
) -> str:
    """A line of `upd regions` for the region from pos to end_pos; `changed` replaces values of its fourth column."""
    outer_pos, outer_end_pos = outer or (pos - 50, end_pos + 50)
    values: dict[str, object] = {
        'ORIGIN': origin,
        'TYPE': kind,
        'LOW_SIZE': end_pos - pos,
        'INF_SITES': sites,
        'SNPS': sites * 9,
        'HET_HOM': f'{sites * 4}/{sites * 5}',
        'OPP_SITES': 0,
        'START_LOW': outer_pos,
        'END_HIGH': outer_end_pos,
        'HIGH_SIZE': outer_end_pos - outer_pos,
    } | changed
    return f'{chrom}\t{pos - 1}\t{end_pos}\t' + ';'.join(f'{key}={value}' for key, value in values.items())


def folder(
    work: pathlib.Path,
    family: str,
    ped: pathlib.Path,
    vcf: pathlib.Path,
    runs: Sequence[str] = (),
    regions: Mapping[str, Sequence[str]] | None = None,
    sites: int = 5,
) -> pathlib.Path:
    """The folder that roh_upd.convert writes for a family from these lines of the tools, under `work`.

    `regions` holds the lines of `upd regions` by child, and by default none for each child that the family's pedigree
    and the VCF's samples give.
    """
    work.mkdir()
    (work / 'roh.txt').write_text(HEADER + ''.join(f'{line}\n' for line in runs))
    if regions is None:
        planned = roh_upd.plan(family, ped, vcf, archive(work / 'planned.zip'))
        regions = {trio.child: () for trio in planned.trios}
    beds = {}
    for child, lines in regions.items():
        beds[child] = work / f'{child}.bed'
        beds[child].write_text(''.join(f'{line}\n' for line in lines))
    roh_upd.convert(
        family, ped, vcf, archive(work / 'invented.zip'), work / 'roh.txt', beds, work / 'roh_upd', sites=sites, **RAN
    )
    return work / 'roh_upd'


def written(
    folder: pathlib.Path,
    samples: Sequence[str],
    children: Sequence[roh_upd.Trio],
    runs: Sequence[Mapping[str, Any]] = (),
    regions: Sequence[Mapping[str, Any]] = (),
) -> pathlib.Path:
    """A folder as roh_upd.convert writes one, with these rows, for a family whose inputs nothing has checked."""
    folder.mkdir()
    pq.write_table(pa.Table.from_pylist(list(runs), schema=roh_upd.ROH), folder / roh_upd.ROH_NAME)
    if children:
        pq.write_table(pa.Table.from_pylist(list(regions), schema=roh_upd.UPD), folder / roh_upd.UPD_NAME)
    frequency = next(field for field in CONFIG if field['field'] == 'AF_joint')
    info = roh_upd.RohUpdInfo(
        bcftools_version=RAN['bcftools_version'],
        upd_version=RAN['upd'][0],
        upd_commit=RAN['upd'][1],
        preparation=RAN['preparation'],
        roh_command=RAN['roh_command'],
        upd_command=RAN['upd_command'],
        frequency_archive=roh_upd.FrequencyArchive(
            archive='invented.zip',
            bytes=100,
            sha256='a' * 64,
            echtvar_version=RAN['echtvar_version'],
            field=echtvar_values.ArchiveField(**frequency, missing_string='MISSING'),
        ),
        sites=5,
        roh_members=tuple(samples),
        upd_children=tuple(children),
    )
    (folder / roh_upd.INFO_NAME).write_text(info.model_dump_json(indent=2) + '\n')
    return folder
