"""The preparation of the catalogue of branch points, from invented datasets and an invented chain file.

Every branch point, source, gene and transcript here is invented, and so is the chain file: it describes no assembly.
The datasets have the headers and the layout of Datasets S1 and S2 of Zhang et al. (PNAS 2022), and the chain file the
format of UCSC's. Its chains take stretches of an invented chr1 of 10,000 bases to another place on it, to its other
strand, apart from one another, to another sequence or to two places, and leave the rest of it out.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import pathlib
from collections.abc import Sequence

import pyarrow.parquet as pq
import pytest

from case_annotation import branch_points, cli, joined_values

S1_HEADER = 'CHROM\tBP_POS\tSTRAND\tBP_NAME\tSOURCE'
S2_HEADER = (
    'CHROM\tPOS\tSTRAND\tBP_NAME\tSOURCE_COUNT\tSOURCE_LIST\tGENE\tTRANSCRIPT_INTRON\tINTRON_LENGTH\tINTRON_TYPE\t'
    'DIST_3SS\tRANK\tCONSENSUS\tENERGY\tBP_MAF\tBP_GERP\tBP_PHYLOP\tBP-2_MAF\tBP-2_GERP\tBP-2_PHYLOP'
)
BOTH = 'eBP_Invented,cBP_Invented'
PREDICTED = 'cBP_Invented'

# A chain's header gives its score, then the name, size, strand, start and end of the stretch it is from and of the
# stretch it goes to, and its id. Each line after it is the size of a block that both stretches hold, with the bases
# that follow it in the one and in the other alone; the last line is the last block. Positions are 0-based, and those
# of a stretch on the - strand count from the sequence's end.
CHAINS = (
    # 1,000 bases that lie 500 bases further on.
    ('chain 1000 chr1 10000 + 1000 2000 chr1 10000 + 1500 2500 1', '1000'),
    # 1,000 bases that are turned around where they lie.
    ('chain 900 chr1 10000 + 3000 4000 chr1 10000 - 6000 7000 2', '1000'),
    # Two blocks that are next to one another before the lift and ten bases apart after it.
    ('chain 800 chr1 10000 + 5000 5200 chr1 10000 + 5000 5210 3', '100 0 10', '100'),
    # 100 bases that go to another sequence.
    ('chain 700 chr1 10000 + 7000 7100 chr1_invented_alt 500 + 0 100 4', '100'),
    # Two chains that hold the same 50 bases.
    ('chain 600 chr1 10000 + 8000 8100 chr1 10000 + 8000 8100 5', '100'),
    ('chain 500 chr1 10000 + 8050 8150 chr1 10000 + 9050 9150 6', '100'),
    ('chain 400 chr2 10000 + 0 1000 chr2 10000 + 0 1000 7', '1000'),
)


def _chain(path: pathlib.Path, chains: Sequence[Sequence[str]] = CHAINS, ensembl: bool = False) -> pathlib.Path:
    """A chain file, gzipped where the path's name says so, with its sequences named as Ensembl names them or not."""
    text = '\n\n'.join('\n'.join(chain) for chain in chains) + '\n'
    if ensembl:
        text = text.replace(' chr', ' ')
    path.write_bytes(gzip.compress(text.encode()) if path.suffix == '.gz' else text.encode())
    return path


def _s1_row(chrom: str, pos: int, strand: str, kind: str, base: str, sources: str) -> str:
    return '\t'.join([chrom, str(pos), strand, f'{kind}BP_{chrom}_{pos}_{strand}_{base}', sources])


def _s2_row(
    s1_row: str,
    introns: str,
    distance: int,
    rank: str,
    consensus: str = '1:YTNAY',
    gene: str = 'INVENTED1',
    length: int = 900,
    kind: str = 'major',
) -> str:
    """A row of Dataset S2 for a branch point of Dataset S1, with invented values in the columns that are not read."""
    chrom, pos, strand, name, sources = s1_row.split('\t')
    count = str(len(sources.split(',')))
    read = [chrom, pos, strand, name, count, sources, gene, introns, str(length), kind, str(distance), rank, consensus]
    return '\t'.join([*read, '-1.50', '.', '0.5', '0.25', '-3.2', '1.5', '-0.125'])


# The branch points of the invented Dataset S1, out of the order of their positions as a public file need not be in it.
ON_CHR2 = _s1_row('chr2', 301, '-', 'c', 'T', PREDICTED)
FLIPPED = _s1_row('chr1', 3201, '-', 'e', 'A', 'eBP_Invented')
MOVED = _s1_row('chr1', 1101, '+', 'm', 'A', BOTH)
NEIGHBOUR = _s1_row('chr1', 1103, '+', 'c', 'C', PREDICTED)
SPLIT = _s1_row('chr1', 5101, '+', 'c', 'A', PREDICTED)
ELSEWHERE = _s1_row('chr1', 7051, '+', 'c', 'A', PREDICTED)
TWICE = _s1_row('chr1', 8061, '+', 'c', 'A', PREDICTED)
NOWHERE = _s1_row('chr1', 501, '+', 'c', 'A', PREDICTED)
S1_ROWS = (ON_CHR2, FLIPPED, MOVED, NEIGHBOUR, SPLIT, ELSEWHERE, TWICE, NOWHERE)
# The invented Dataset S2: one branch point in two introns, the one that is turned around, and one that goes to
# another sequence.
S2_ROWS = (
    _s2_row(MOVED, 'ENST02_IVS1,ENST03_IVS1', -30, '#1/1', '2:YTNA', length=4000),
    _s2_row(MOVED, 'ENST01_IVS2', -24, '#1/2'),
    _s2_row(FLIPPED, 'ENST04_IVS7', -19, '#2/3', '0:none', 'INVENTED2', 85, 'minor'),
    _s2_row(ELSEWHERE, 'ENST05_IVS1', -40, '#1/1'),
)
# What the prepared file holds of them: chromosome, position, strand and name, then for a row of Dataset S2 the gene,
# the introns, the kind and length of the intron, the distance, the rank, the total and the consensus, and last the
# sources and the name in the datasets.
NOTHING = ('',) * 8
NEAR = ('INVENTED1', 'ENST01_IVS2', 'major', '900', '-24', '1', '2', '1:YTNAY')
FAR = ('INVENTED1', 'ENST02_IVS1,ENST03_IVS1', 'major', '4000', '-30', '1', '1', '2:YTNA')
MINOR = ('INVENTED2', 'ENST04_IVS7', 'minor', '85', '-19', '2', '3', '0:none')
WRITTEN = (
    ('chr1', '1601', '+', 'mBP_chr1_1601_+_A', *NEAR, BOTH, 'mBP_chr1_1101_+_A'),
    ('chr1', '1601', '+', 'mBP_chr1_1601_+_A', *FAR, BOTH, 'mBP_chr1_1101_+_A'),
    ('chr1', '1603', '+', 'cBP_chr1_1603_+_C', *NOTHING, PREDICTED, 'cBP_chr1_1103_+_C'),
    # The base at 3201 of the - strand is the 201st of the stretch, which is its 800th once it is turned around.
    ('chr1', '3800', '+', 'eBP_chr1_3800_+_A', *MINOR, 'eBP_Invented', 'eBP_chr1_3201_-_A'),
    ('chr2', '301', '-', 'cBP_chr2_301_-_T', *NOTHING, PREDICTED, 'cBP_chr2_301_-_T'),
)
LIFTED = {
    'clean': 3,
    'strand_flipped': 1,
    'unmapped': 1,
    'several_places': 1,
    'other_chromosome': 1,
    'split_window': 1,
}


def _dataset(path: pathlib.Path, header: str, rows: Sequence[str]) -> pathlib.Path:
    text = ('\n'.join([header, *rows]) + '\n').encode()
    path.write_bytes(gzip.compress(text) if path.suffix == '.gz' else text)
    return path


def _prepare(
    tmp_path: pathlib.Path,
    s1: Sequence[str] = S1_ROWS,
    s2: Sequence[str] = S2_ROWS,
    s1_header: str = S1_HEADER,
    s2_header: str = S2_HEADER,
    chain: pathlib.Path | None = None,
    names: tuple[str, str] = ('s1.invented.txt', 's2.invented.txt'),
) -> branch_points.Catalogue:
    tmp_path.mkdir(exist_ok=True)
    return branch_points.prepare(
        _dataset(tmp_path / names[0], s1_header, s1),
        _dataset(tmp_path / names[1], s2_header, s2),
        chain or _chain(tmp_path / 'invented.chain.gz'),
        tmp_path / 'out',
    )


def _lines(folder: pathlib.Path) -> tuple[list[str], list[tuple[str, ...]]]:
    """The prepared file of a folder: the lines before its header, and its rows after it."""
    with gzip.open(folder / branch_points.FILE_NAME, 'rt') as f:
        lines = f.read().splitlines()
    before = [line for line in lines if line.startswith(branch_points.PREAMBLE)]
    assert lines[len(before)].split('\t') == list(branch_points.COLUMNS)
    return before, [tuple(line.split('\t')) for line in lines[len(before) + 1 :]]


def test_the_file_holds_every_branch_point_that_lifts_with_the_two_bases_before_it(tmp_path: pathlib.Path) -> None:
    record = _prepare(tmp_path)
    _before, rows = _lines(tmp_path / 'out')
    # A row for each row of Dataset S2 and one for a branch point that only Dataset S1 lists, sorted by chromosome,
    # position and strand, and a branch point's rows by their introns.
    assert rows == list(WRITTEN)
    # Every row fills its last two columns, so no line ends in an empty field, which a tool that trims the ends of
    # lines would remove.
    assert all(row[-2] and row[-1] for row in rows)
    assert record.lifted.model_dump() == LIFTED
    assert (record.branch_points, record.branch_points_with_introns) == (len(S1_ROWS), 3)
    assert (record.written_branch_points, record.written_branch_points_with_introns, record.written_rows) == (4, 2, 5)
    # Every branch point of Dataset S1 is counted once, and the two kinds that lift are the ones written.
    assert sum(LIFTED.values()) == record.branch_points
    assert record.lifted.clean + record.lifted.strand_flipped == record.written_branch_points


def test_a_chain_file_under_ensembls_names_gives_the_same_rows(tmp_path: pathlib.Path) -> None:
    _prepare(tmp_path / 'ucsc')
    _prepare(tmp_path / 'ensembl', chain=_chain(tmp_path / 'invented.ensembl.chain', ensembl=True))
    assert _lines(tmp_path / 'ensembl' / 'out')[1] == _lines(tmp_path / 'ucsc' / 'out')[1] == list(WRITTEN)


def test_the_same_three_files_give_the_same_bytes(tmp_path: pathlib.Path) -> None:
    first, second = _prepare(tmp_path / 'first'), _prepare(tmp_path / 'second')
    assert first == second
    written = (tmp_path / 'first' / 'out' / branch_points.FILE_NAME).read_bytes()
    assert written == (tmp_path / 'second' / 'out' / branch_points.FILE_NAME).read_bytes()
    assert (first.file, first.bytes, first.sha256) == (
        branch_points.FILE_NAME,
        len(written),
        hashlib.sha256(written).hexdigest(),
    )
    # bgzip ends a file with an empty block, which a reader takes as its end.
    assert written.endswith(bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000'))


def test_the_file_and_the_record_name_the_three_files_by_their_checksums(tmp_path: pathlib.Path) -> None:
    record = _prepare(tmp_path, names=('s1.invented.txt.gz', 's2.invented.txt.gz'))
    files = {
        'dataset_s1': tmp_path / 's1.invented.txt.gz',
        'dataset_s2': tmp_path / 's2.invented.txt.gz',
        'chain': tmp_path / 'invented.chain.gz',
    }
    before, rows = _lines(tmp_path / 'out')
    # The datasets were gzipped here, and give the rows that they give as text.
    assert rows == list(WRITTEN)
    assert 'doi:10.1073/pnas.2211194119' in before[0]
    for line, (role, path) in zip(before[1:4], files.items(), strict=True):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert line == f'{branch_points.PREAMBLE}{role}: {path.name} sha256 {digest}'
        assert getattr(record, role) == branch_points.InputFile(
            name=path.name, bytes=path.stat().st_size, sha256=digest
        )
    assert before[4] == f'{branch_points.PREAMBLE}pyliftover: {record.pyliftover}'
    assert len(before) == 5
    assert branch_points.read_record(tmp_path / 'out') == record
    assert {path.name for path in (tmp_path / 'out').iterdir()} == {branch_points.FILE_NAME, branch_points.RECORD_NAME}


def test_the_command_writes_the_folder(tmp_path: pathlib.Path) -> None:
    s1 = _dataset(tmp_path / 's1.txt', S1_HEADER, S1_ROWS)
    s2 = _dataset(tmp_path / 's2.txt', S2_HEADER, S2_ROWS)
    chain = _chain(tmp_path / 'invented.chain')
    arguments = ['--dataset-s1', str(s1), '--dataset-s2', str(s2), '--chain', str(chain)]
    cli.main(['branch-points', *arguments, '--output', str(tmp_path / 'out')])
    record = json.loads((tmp_path / 'out' / branch_points.RECORD_NAME).read_text())
    assert [record[role]['name'] for role in ('dataset_s1', 'dataset_s2', 'chain')] == [
        's1.txt',
        's2.txt',
        'invented.chain',
    ]
    assert _lines(tmp_path / 'out')[1] == list(WRITTEN)


def test_the_tables_read_the_prepared_file(tmp_path: pathlib.Path) -> None:
    _prepare(tmp_path)
    # An allele at the base of the branch point that moved, which is the base two before its neighbour, and one at
    # the place that the branch point had before the lift.
    alleles = [('chr1', 1101, 'A', 'G'), ('chr1', 1601, 'A', 'G')]
    sites = tmp_path / 'union.vcf'
    sites.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        + ''.join(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n' for chrom, pos, ref, alt in alleles)
    )
    joined_values.convert(sites, {'branch_points': tmp_path / 'out' / branch_points.FILE_NAME}, {}, tmp_path / 'joined')
    rows = pq.read_table(tmp_path / 'joined' / joined_values.TABLE_NAME).to_pylist()
    assert rows[0] == {'variant_id': 'chr1-1101-A-G', 'branch_points': None}
    found = [(entry['name'], entry['hit'], entry['acceptor_offset']) for entry in rows[1]['branch_points']]
    assert found == [
        ('mBP_chr1_1601_+_A', [0], -24),
        ('mBP_chr1_1601_+_A', [0], -30),
        ('cBP_chr1_1603_+_C', [-2], None),
    ]
    # The record of the run names the prepared file by the checksum that the preparation recorded.
    (joined,) = joined_values.read_info(tmp_path / 'joined').files
    assert joined.sha256 == branch_points.read_record(tmp_path / 'out').sha256


def _edited(row: str, column: int, value: str) -> str:
    fields = row.split('\t')
    fields[column] = value
    return '\t'.join(fields)


S1_REFUSED: dict[str, tuple[dict[str, object], str]] = {
    'another header': (
        {'s1_header': S1_HEADER.replace('BP_POS', 'POS')},
        r"its header is \['CHROM', 'POS', .*and that of Dataset S1 is \['CHROM', 'BP_POS'",
    ),
    'a row of four fields': (
        {'s1': [*S1_ROWS, MOVED.rsplit('\t', 1)[0]]},
        'line 10 has 4 fields, and Dataset S1 has 5',
    ),
    'a position that is not the one in the name': (
        {'s1': [_edited(MOVED, 1, '1102')]},
        r"line 2: the name 'mBP_chr1_1101_\+_A' is not a kind, chr1, 1102, \+ and a base",
    ),
    'a chromosome that is not the one in the name': (
        {'s1': [_edited(MOVED, 0, 'chr2')]},
        r"the name 'mBP_chr1_1101_\+_A' is not a kind, chr2, 1101, \+ and a base",
    ),
    'a strand that is not the one in the name': (
        {'s1': [_edited(MOVED, 2, '-')]},
        r"the name 'mBP_chr1_1101_\+_A' is not a kind, chr1, 1101, - and a base",
    ),
    'a position that is no number': (
        {'s1': [_edited(_edited(MOVED, 1, 'x'), 3, 'mBP_chr1_x_+_A')]},
        r"the name 'mBP_chr1_x_\+_A' is not a kind, chr1, x, \+ and a base",
    ),
    'a name without a base': (
        {'s1': [_edited(MOVED, 3, 'mBP_chr1_1101_+_N')]},
        r"the name 'mBP_chr1_1101_\+_N' is not a kind",
    ),
    'a chromosome that the tables lack': (
        {'s1': [_s1_row('chrUn', 100, '+', 'c', 'A', PREDICTED)]},
        'line 2: chrUn is none of chr1 to chr22, chrX, chrY and chrM',
    ),
    'sources that are no names of sources': (
        {'s1': [_edited(MOVED, 4, 'eBP_Invented;cBP_Invented')]},
        r"the sources 'eBP_Invented;cBP_Invented' are no names such as eBP_Mercer",
    ),
    'no source': ({'s1': [_edited(NEIGHBOUR, 4, '')]}, "the sources '' are no names"),
    'a predicted source of a branch point that is named as experimental': (
        {'s1': [_edited(FLIPPED, 4, BOTH)]},
        r"the name 'eBP_chr1_3201_-_A' states other kinds of source than 'eBP_Invented,cBP_Invented'",
    ),
    'a branch point twice': ({'s1': [MOVED, NEIGHBOUR, MOVED]}, r'line 4: a second row of chr1:1101 on the \+ strand'),
    'a position of a strand under two names': (
        {'s1': [MOVED, _s1_row('chr1', 1101, '+', 'c', 'A', PREDICTED)]},
        r'line 3: a second row of chr1:1101 on the \+ strand',
    ),
}


@pytest.mark.parametrize('case', S1_REFUSED)
def test_a_file_that_is_not_dataset_s1_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    written, message = S1_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _prepare(tmp_path, **{'s2': (), **written})  # type: ignore[arg-type]
    # Nothing is written for a dataset that is refused.
    assert not (tmp_path / 'out').exists()


FIRST = S2_ROWS[1]
S2_REFUSED: dict[str, tuple[dict[str, object], str]] = {
    'another header': (
        {'s2_header': S2_HEADER.replace('DIST_3SS', 'DISTANCE')},
        r"its header is \['CHROM', 'POS', .*'DISTANCE'.*and that of Dataset S2 is",
    ),
    'a row of nineteen fields': (
        {'s2': [FIRST.rsplit('\t', 1)[0]]},
        'line 2 has 19 fields, and Dataset S2 has 20',
    ),
    'a branch point that Dataset S1 lacks': (
        {'s2': [_s2_row(_s1_row('chr1', 1201, '+', 'c', 'A', PREDICTED), 'ENST01_IVS2', -24, '#1/1')]},
        r'line 2: Dataset S1 has no branch point cBP_chr1_1201_\+_A',
    ),
    'a name that is not of its position': (
        {'s2': [_edited(FIRST, 1, '1102')]},
        r"the name 'mBP_chr1_1101_\+_A' is not a kind, chr1, 1102",
    ),
    'other sources than Dataset S1 gives': (
        {'s2': [_edited(_edited(FIRST, 5, PREDICTED), 4, '1')]},
        r"1 sources, 'cBP_Invented', and Dataset S1 gives mBP_chr1_1101_\+_A 'eBP_Invented,cBP_Invented'",
    ),
    'another number of sources than it lists': (
        {'s2': [_edited(FIRST, 4, '3')]},
        "3 sources, 'eBP_Invented,cBP_Invented', and Dataset S1 gives",
    ),
    'no gene': ({'s2': [_edited(FIRST, 6, '')]}, r'the row of mBP_chr1_1101_\+_A is without a gene$'),
    'introns that are no labels': (
        {'s2': [_edited(FIRST, 7, 'ENST01_IVS2;ENST02_IVS1')]},
        'is without introns such as ENST00000000001_IVS2 with commas between',
    ),
    'an intron without a number': (
        {'s2': [_edited(FIRST, 7, 'ENST01_IVS')]},
        'is without introns such as ENST00000000001_IVS2',
    ),
    'a length that is no number': (
        {'s2': [_edited(FIRST, 8, '9e2')]},
        "is without a whole number as its intron's length",
    ),
    'another type of intron': (
        {'s2': [_edited(FIRST, 9, 'Major')]},
        'is without one of major, minor as its type of intron',
    ),
    'a distance that is not negative': (
        {'s2': [_edited(FIRST, 10, '24')]},
        'is without a negative whole number as its distance to the acceptor',
    ),
    'no distance': (
        {'s2': [_edited(FIRST, 10, '')]},
        'is without a negative whole number as its distance to the acceptor',
    ),
    'a rank without a total': ({'s2': [_edited(FIRST, 11, '1')]}, 'is without a rank such as #1/2'),
    'a rank above its total': (
        {'s2': [_edited(FIRST, 11, '#3/2')]},
        'is without a rank such as #1/2 that is at most its total',
    ),
    'another consensus': (
        {'s2': [_edited(FIRST, 12, 'YTNAY')]},
        'is without one of 1:YTNAY, 2:YTNA, 3:TNA, 4:YNA, 0:none as its consensus',
    ),
    'two things at once': (
        {'s2': [_edited(_edited(FIRST, 6, ''), 12, '5:NA')]},
        'is without a gene; one of 1:YTNAY',
    ),
    'a branch point twice in the same introns': (
        {'s2': [FIRST, _edited(FIRST, 10, '-25')]},
        r'line 3: a second row of mBP_chr1_1101_\+_A in ENST01_IVS2',
    ),
}


@pytest.mark.parametrize('case', S2_REFUSED)
def test_a_file_that_is_not_dataset_s2_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    written, message = S2_REFUSED[case]
    with pytest.raises(ValueError, match=message):
        _prepare(tmp_path, **written)  # type: ignore[arg-type]
    assert not (tmp_path / 'out').exists()


def test_a_dataset_that_is_cut_off_is_refused(tmp_path: pathlib.Path) -> None:
    rows = [_s1_row('chr1', 1000 + number, '+', 'c', 'A', PREDICTED) for number in range(500)]
    whole = _dataset(tmp_path / 'whole.txt.gz', S1_HEADER, rows)
    cut = tmp_path / 's1.cut.txt.gz'
    cut.write_bytes(whole.read_bytes()[: whole.stat().st_size // 2])
    s2 = _dataset(tmp_path / 's2.txt', S2_HEADER, ())
    with pytest.raises(ValueError, match=r's1\.cut\.txt\.gz: not Dataset S1'):
        branch_points.prepare(cut, s2, _chain(tmp_path / 'invented.chain'), tmp_path / 'out')


def test_a_chain_file_without_a_chromosome_of_the_datasets_is_refused(tmp_path: pathlib.Path) -> None:
    chain = _chain(tmp_path / 'chr1-alone.chain', [chain for chain in CHAINS if ' chr2 ' not in chain[0]])
    with pytest.raises(ValueError, match=r'chr1-alone\.chain: no chain from chr2 or 2, which the datasets hold'):
        _prepare(tmp_path, chain=chain)
    assert not (tmp_path / 'out').exists()


CHAIN_REFUSED = {
    'a chain from a stretch of the - strand': ('chain 1000 chr1 10000 - 1000 2000 chr1 10000 + 1500 2500 1', '1000'),
    'a size that is no number': ('chain 1000 chr1 big + 1000 2000 chr1 10000 + 1500 2500 1', '1000'),
    'a header that ends early': ('chain 1000 chr1 10000 +', '1000'),
    'a last block of two numbers': ('chain 1000 chr1 10000 + 1000 2000 chr1 10000 + 1500 2500 1', '500 500'),
}


@pytest.mark.parametrize('case', CHAIN_REFUSED)
def test_a_file_that_is_no_chain_file_is_refused(tmp_path: pathlib.Path, case: str) -> None:
    chain = _chain(tmp_path / 'malformed.chain', [CHAIN_REFUSED[case]])
    with pytest.raises(ValueError, match=r'malformed\.chain: not a chain file'):
        _prepare(tmp_path, chain=chain)


def test_a_folder_that_exists_is_not_written_into(tmp_path: pathlib.Path) -> None:
    (tmp_path / 'out').mkdir()
    with pytest.raises(FileExistsError):
        _prepare(tmp_path)
    assert not list((tmp_path / 'out').iterdir())
