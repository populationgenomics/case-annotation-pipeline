"""Invented VCFs of a structural-variant caller, a copy-number caller and a repeat genotyper, written as text.

Every record, id and repeat here is invented, in the layout DRAGEN writes each kind of file in. The family is the trio
of test_family_tables: the father S1, who is called haploid on X and Y, the mother S2 and their daughter S3.
"""

from __future__ import annotations

import pathlib
from collections.abc import Mapping, Sequence

SAMPLES = ('S1', 'S2', 'S3')
_CONTIGS = ('chr1', 'chr2', 'chr10', 'chrX', 'chrY', 'chr1_KI270706v1_random')

# One record: the eight fixed columns, the FORMAT keys, and each sample's values.
type Record = tuple[str, ...]


def _meta(kind: str, key: str, number: str, type_: str) -> str:
    return f'##{kind}=<ID={key},Number={number},Type={type_},Description="invented">'


SV_INFO = (
    _meta('INFO', 'SVTYPE', '1', 'String'),
    _meta('INFO', 'SVLEN', '.', 'Integer'),
    _meta('INFO', 'END', '1', 'Integer'),
    _meta('INFO', 'IMPRECISE', '0', 'Flag'),
    _meta('INFO', 'CIPOS', '2', 'Integer'),
    _meta('INFO', 'CIEND', '2', 'Integer'),
    _meta('INFO', 'MATEID', '.', 'String'),
    _meta('INFO', 'EVENT', '1', 'String'),
)
SV_FORMAT = (
    _meta('FORMAT', 'GT', '1', 'String'),
    _meta('FORMAT', 'FT', '1', 'String'),
    _meta('FORMAT', 'GQ', '1', 'Integer'),
    _meta('FORMAT', 'PL', 'G', 'Integer'),
    _meta('FORMAT', 'PR', 'R', 'Integer'),
    _meta('FORMAT', 'SR', 'R', 'Integer'),
    _meta('FORMAT', 'VF', 'R', 'Integer'),
    _meta('FORMAT', 'DQ', '1', 'Float'),
    _meta('FORMAT', 'DN', '1', 'String'),
)
CNV_INFO = (
    _meta('INFO', 'SVTYPE', '1', 'String'),
    _meta('INFO', 'SVLEN', '.', 'Integer'),
    _meta('INFO', 'END', '1', 'Integer'),
    _meta('INFO', 'REFLEN', '1', 'Integer'),
)
CNV_FORMAT = (
    _meta('FORMAT', 'GT', '1', 'String'),
    _meta('FORMAT', 'SM', '1', 'Float'),
    _meta('FORMAT', 'CN', '1', 'Integer'),
    _meta('FORMAT', 'BC', '1', 'Integer'),
    _meta('FORMAT', 'PE', '2', 'Integer'),
    _meta('FORMAT', 'QS', '1', 'Integer'),
    _meta('FORMAT', 'FT', '1', 'String'),
    _meta('FORMAT', 'DQ', '1', 'Float'),
    _meta('FORMAT', 'DN', '1', 'String'),
)
REPEAT_INFO = (
    _meta('INFO', 'END', '1', 'Integer'),
    _meta('INFO', 'REF', '1', 'Integer'),
    _meta('INFO', 'REPID', '1', 'String'),
    _meta('INFO', 'RU', '1', 'String'),
    _meta('INFO', 'VARID', '1', 'String'),
)
REPEAT_FORMAT = (
    _meta('FORMAT', 'GT', '1', 'String'),
    _meta('FORMAT', 'SO', '1', 'String'),
    _meta('FORMAT', 'REPCN', '1', 'String'),
    _meta('FORMAT', 'REPCI', '1', 'String'),
    _meta('FORMAT', 'ADSP', '1', 'String'),
    _meta('FORMAT', 'ADFL', '1', 'String'),
    _meta('FORMAT', 'ADIR', '1', 'String'),
    _meta('FORMAT', 'LC', '1', 'Float'),
)


def write(
    path: pathlib.Path, meta: Sequence[str], records: Sequence[Record], samples: Sequence[str] = SAMPLES
) -> pathlib.Path:
    """A VCF of the records under the header lines given, as text."""
    lines = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=MinQUAL,Description="invented">',
        '##FILTER=<ID=SampleFT,Description="invented">',
        *(f'##contig=<ID={contig},length=200000000>' for contig in _CONTIGS),
        *meta,
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *samples]),
        *('\t'.join(record) for record in records),
    ]
    path.write_text('\n'.join(lines) + '\n')
    return path


_SV_KEYS = 'GT:FT:GQ:PL:PR:SR:VF:DQ:DN'
_SV_REF = '0/0:PASS:50:0,50,500:15,0:12,0:20,0:.:.'
_SV_HET = '0/1:PASS:40:90,0,80:10,5:8,4:12,6:.:.'
# In the file's order, which is not the tables': chr10 stands before chr2, and chr1's unplaced contig between them.
SV_RECORDS: tuple[Record, ...] = (
    # A deletion written as a symbolic allele, imprecise, heterozygous in the father and homozygous in the daughter.
    (
        'chr1',
        '1000',
        'inv:DEL:1',
        'T',
        '<DEL>',
        '80',
        'PASS',
        'SVTYPE=DEL;END=5000;SVLEN=-4000;IMPRECISE;CIPOS=-10,10;CIEND=-20,20',
        _SV_KEYS,
        _SV_HET,
        _SV_REF,
        '1/1:LowGQ:30:400,40,0:0,9:0,7:1,11:3.5:DeNovo',
    ),
    # An insertion written as a sequence: the allele chr10-7-G-GA of the family's small variants.
    (
        'chr10',
        '7',
        'inv:INS:1',
        'G',
        'GA',
        '90',
        'PASS',
        'SVTYPE=INS;END=7;SVLEN=1',
        _SV_KEYS,
        _SV_REF,
        _SV_REF,
        _SV_HET,
    ),
    # A breakend and its mate, parts of one event.
    (
        'chr10',
        '900',
        'inv:BND:1:1',
        'T',
        '[chr2:500[T',
        '70',
        'PASS',
        'SVTYPE=BND;MATEID=inv:BND:1:0;EVENT=inv:event:1',
        _SV_KEYS,
        _SV_REF,
        _SV_HET,
        _SV_HET,
    ),
    ('chr1_KI270706v1_random', '50', 'inv:DEL:9', 'T', '<DEL>', '80', 'PASS', 'SVTYPE=DEL;END=90;SVLEN=-40', _SV_KEYS)
    + (_SV_HET,) * 3,
    # A deletion written as a sequence that is no allele of the family's small variants, with two filters.
    (
        'chr2',
        '300',
        'inv:DEL:2',
        'ACGT',
        'A',
        '20',
        'MinQUAL;SampleFT',
        'SVTYPE=DEL;END=303;SVLEN=-3',
        'GT:FT:GQ',
        '0/1:PASS:9',
        '0/0:PASS:8',
        './.:.:.',
    ),
    (
        'chr2',
        '500',
        'inv:BND:1:0',
        'G',
        'G]chr10:900]',
        '70',
        'PASS',
        'SVTYPE=BND;MATEID=inv:BND:1:1;EVENT=inv:event:1;IMPRECISE;CIPOS=-5,5',
        _SV_KEYS,
        _SV_REF,
        _SV_HET,
        _SV_HET,
    ),
    # Outside X's pseudo-autosomal regions the father is called haploid, with two likelihoods.
    (
        'chrX',
        '5000100',
        'inv:DEL:3',
        'C',
        '<DEL>',
        '60',
        'PASS',
        'SVTYPE=DEL;END=5000900;SVLEN=-800',
        _SV_KEYS,
        '1:PASS:35:80,0:0,6:0,5:0,8:.:.',
        _SV_HET,
        _SV_REF,
    ),
    # On Y the father is called haploid and the mother and the daughter have no call.
    (
        'chrY',
        '3000000',
        'inv:DUP:TANDEM:1',
        'A',
        '<DUP:TANDEM>',
        '55',
        'PASS',
        'SVTYPE=DUP;END=3002000;SVLEN=2000',
        _SV_KEYS,
        '0:PASS:45:0,70:9,0:7,0:11,0:.:.',
        '.:.:.:.:.:.:.:.:.',
        '.:.:.:.:.:.:.:.:.',
    ),
)
# The tables' order of the records above, without the one off the primary chromosomes.
SV_IDS = ('inv:DEL:1', 'inv:DEL:2', 'inv:BND:1:0', 'inv:INS:1', 'inv:BND:1:1', 'inv:DEL:3', 'inv:DUP:TANDEM:1')
# The records written as a sequence, each with the allele a normalisation gives it: the insertion's own, and the
# deletion's two bases to the left.
SV_ALLELES = {'inv:INS:1': ('chr10', '7', 'G', 'GA'), 'inv:DEL:2': ('chr2', '298', 'TGAC', 'T')}


def write_sv_alleles(path: pathlib.Path, alleles: Mapping[str, tuple[str, str, str, str]] = SV_ALLELES) -> pathlib.Path:
    """The normalised sequences of a structural-variant VCF: each allele under its record's ID, without samples."""
    lines = [
        '##fileformat=VCFv4.2',
        *(f'##contig=<ID={contig},length=200000000>' for contig in _CONTIGS),
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']),
        *('\t'.join([chrom, pos, sv_id, ref, alt, '.', '.', '.']) for sv_id, (chrom, pos, ref, alt) in alleles.items()),
    ]
    path.write_text('\n'.join(lines) + '\n')
    return path


_CNV_KEYS = 'GT:SM:CN:BC:PE:QS:FT:DN:DQ'


def cnv_call(gt: str, cn: str) -> str:
    """A member's call of a segment with this genotype and copy number, and invented values for the rest.

    A member without a copy number has no value at all, as DRAGEN writes a female on Y.
    """
    if cn == '.':
        return ':'.join([gt, *('.' for _ in _CNV_KEYS.split(':')[1:])])
    return f'{gt}:0.75:{cn}:12:3,2:40:PASS:.:.'


def cnv_record(
    chrom: str, pos: int, sv_id: str, alt: str, calls: Sequence[tuple[str, str]], length: int = 10_000
) -> Record:
    """A segment of the length given, with each member's genotype and copy number."""
    end = pos + length
    info = f'END={end};REFLEN={length}'
    if alt != '.':
        lengths = ','.join(str(-length if allele == '<DEL>' else length) for allele in alt.split(','))
        info = f'SVTYPE=CNV;SVLEN={lengths};{info}'
    return (chrom, str(pos), sv_id, 'N', alt, '.', 'PASS', info, _CNV_KEYS, *(cnv_call(gt, cn) for gt, cn in calls))


_NEUTRAL = ('./.', '2')
_NO_CALL = ('.', '.')
# Each segment with the father's, the mother's and the daughter's genotype and copy number, as DRAGEN writes them.
CNV_RECORDS: tuple[Record, ...] = (
    # A segment at the reference's copy number, which has no alternate allele and no row.
    cnv_record('chr1', 10_000, 'inv:REF:1', '.', (_NEUTRAL, _NEUTRAL, _NEUTRAL)),
    cnv_record('chr1', 20_000, 'inv:LOSS:1', '<DEL>', (('0/1', '1'), _NEUTRAL, ('1/1', '0'))),
    cnv_record('chr2', 40_000, 'inv:GAIN:1', '<DUP>', (('./1', '3'), _NEUTRAL, ('./1', '4')), length=50_000),
    # A loss in the father and a gain in the mother.
    cnv_record('chr10', 5_000, 'inv:CNV:1', '<DEL>,<DUP>', (('0/1', '1'), ('./2', '3'), _NEUTRAL)),
    # Outside X's pseudo-autosomal regions the father has one copy, and is called haploid.
    cnv_record('chrX', 5_000_000, 'inv:LOSS:2', '<DEL>', (('1', '0'), ('0/1', '1'), _NEUTRAL)),
    cnv_record('chrX', 6_000_000, 'inv:GAIN:2', '<DUP>', (('0', '1'), ('./1', '3'), _NEUTRAL)),
    cnv_record('chrX', 7_000_000, 'inv:GAIN:3', '<DUP>', (('1', '2'), _NEUTRAL, _NEUTRAL)),
    # On Y the mother and the daughter have no copy number.
    cnv_record('chrY', 3_000_000, 'inv:LOSS:3', '<DEL>', (('1', '0'), _NO_CALL, _NO_CALL)),
    cnv_record('chrY', 4_000_000, 'inv:GAIN:4', '<DUP>', (('1', '3'), _NO_CALL, _NO_CALL)),
)
# Per segment with a row, each member's (carries, zygosity).
_LOSS_HET, _LOSS_HOM, _GAIN = (True, 'het'), (True, 'hom_alt'), (True, None)
_AT_EXPECTED, _WITHOUT = (False, 'hom_ref'), (None, None)
CNV_CALLS: dict[str, tuple[tuple[bool | None, str | None], ...]] = {
    'inv:LOSS:1': (_LOSS_HET, _AT_EXPECTED, _LOSS_HOM),
    'inv:GAIN:1': (_GAIN, _AT_EXPECTED, _GAIN),
    'inv:CNV:1': (_LOSS_HET, _GAIN, _AT_EXPECTED),
    'inv:LOSS:2': (_LOSS_HOM, _LOSS_HET, _AT_EXPECTED),
    'inv:GAIN:2': (_AT_EXPECTED, _GAIN, _AT_EXPECTED),
    'inv:GAIN:3': (_GAIN, _AT_EXPECTED, _AT_EXPECTED),
    'inv:LOSS:3': (_LOSS_HOM, _WITHOUT, _WITHOUT),
    'inv:GAIN:4': (_GAIN, _WITHOUT, _WITHOUT),
}

_REPEAT_KEYS = 'GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC'
# The invented catalogue: each repeat's place, id, unit and reference size.
_REPEATS = (
    ('chr2', 700, 'RPTB_GCC', 'GCC', 7),
    ('chr1', 2_000, 'RPTA', 'CAG', 10),
    ('chrX', 9_000, 'RPTX', 'CGG', 20),
)
REPEAT_IDS = ('RPTA', 'RPTB_GCC', 'RPTX')
# Each member's genotype of each repeat: its alternate alleles and its values. The daughter's longer allele of RPTA is
# sized from reads within the repeat, the father is called haploid on X, and the mother has no call of RPTB_GCC.
REPEAT_CALLS: dict[str, dict[str, tuple[str, str]]] = {
    'S1': {
        'RPTB_GCC': ('.', '0/0:SPANNING/SPANNING:7/7:7-7/7-7:9/9:4/4:0/0:31.5'),
        'RPTA': ('<STR12>', '0/1:SPANNING/SPANNING:10/12:10-10/12-12:5/6:7/8:0/0:30.5'),
        'RPTX': ('<STR23>', '1:SPANNING:23:23-23:9:11:0:15.25'),
    },
    'S2': {
        'RPTB_GCC': ('.', './.:./.:./.:./.:./.:./.:./.:2.5'),
        'RPTA': ('<STR11>', '0/1:SPANNING/SPANNING:10/11:10-10/11-11:6/6:8/8:0/0:29.0'),
        'RPTX': ('<STR23>,<STR30>', '1/2:SPANNING/FLANKING:23/30:23-23/28-34:8/0:10/5:0/0:33.0'),
    },
    'S3': {
        'RPTB_GCC': ('.', '0/0:SPANNING/SPANNING:7/7:7-7/7-7:8/8:5/5:0/0:30.0'),
        'RPTA': ('<STR12>,<STR45>', '1/2:SPANNING/INREPEAT:12/45:12-12/40-58:5/0:7/9:0/4:28.0'),
        'RPTX': ('<STR23>', '1/1:SPANNING/SPANNING:23/23:23-23/23-23:7/7:9/9:0/0:32.0'),
    },
}


def repeat_records(sample: str) -> list[Record]:
    """A member's records of the invented catalogue."""
    records = []
    for chrom, pos, repeat_id, unit, ref_copies in _REPEATS:
        alt, values = REPEAT_CALLS[sample][repeat_id]
        info = f'END={pos + ref_copies * len(unit)};REF={ref_copies};REPID={repeat_id};RU={unit};VARID={repeat_id}'
        records.append((chrom, str(pos), '.', 'C', alt, '.', 'PASS', info, _REPEAT_KEYS, values))
    return records


def write_repeats(folder: pathlib.Path, samples: Sequence[str] = SAMPLES) -> tuple[pathlib.Path, ...]:
    """Each member's repeat VCF in the folder."""
    return tuple(
        write(folder / f'{sample}.repeats.vcf', (*REPEAT_INFO, *REPEAT_FORMAT), repeat_records(sample), (sample,))
        for sample in samples
    )
