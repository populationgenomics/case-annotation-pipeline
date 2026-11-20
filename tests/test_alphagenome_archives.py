"""What builds and checks an echtvar archive of AlphaGenome's scores, on invented rows of DeepMind's three tables."""

from __future__ import annotations

import io
import json
import pathlib
import re
from collections.abc import Mapping, Sequence

import pysam
import pytest

from case_annotation import alphagenome_archives, cli

AVI = alphagenome_archives.TABLES['alphagenome_avi']
BREAKDOWN = alphagenome_archives.TABLES['alphagenome_avi_breakdown']
SPLICING = alphagenome_archives.TABLES['alphagenome_splicing']
_CONTIGS = {'chr1': 1000, 'chr2': 1000, 'chrX': 1000}

_AVI_ROWS = (
    'chr1\t5\tC\tA\t-0.04231\t1.20713',
    'chr1\t5\tC\tG\t0.0\t0.0',
    'chr1\t5\tC\tT\t2.71306\t48.90215',
    'chrX\t7\tG\tA\t-1.137\t97.3',
)
# A row of the breakdown: sixteen contributions, one in scientific notation and one negative, and the two columns of
# the insertions and deletions.
_CONTRIBUTIONS = '0.2817\t0.001293\t5.318e-05\t0.0\t-0.0004126' + '\t0.0125' * 11
_BREAKDOWN_ROWS = (f'chr1\t5\tC\tA\t{_CONTRIBUTIONS}\t0.0\t0.0', f'chr2\t9\tT\tC\t{_CONTRIBUTIONS}\t0\t0.0e0')


def _lines(table: alphagenome_archives.Table, rows: Sequence[str]) -> list[str]:
    return [line + '\n' for line in ('\t'.join(('#CHROM', 'POS', 'REF', 'ALT', *table.columns)), *rows)]


def _vcf(table: alphagenome_archives.Table, rows: Sequence[str], contigs: Mapping[str, int] = _CONTIGS) -> str:
    out = io.StringIO()
    alphagenome_archives.to_vcf(table, _lines(table, rows), contigs, out)
    return out.getvalue()


def _read(tmp_path: pathlib.Path, text: str) -> tuple[pysam.VariantHeader, list[tuple[str, dict[str, object]]]]:
    """A VCF's header and its records, each as its variant id and its INFO fields, as htslib reads the text."""
    path = tmp_path / 'rows.vcf'
    path.write_text(text)
    with pysam.VariantFile(str(path)) as vcf:
        records = [
            (f'{r.chrom}-{r.pos}-{r.ref}-{"/".join(r.alts or ())}', {key: r.info[key] for key in r.info}) for r in vcf
        ]
        return vcf.header, records


# ---- the tables and their configurations ----


@pytest.mark.parametrize('name', alphagenome_archives.TABLES)
def test_each_configuration_encodes_the_kept_columns_of_its_table(name: str) -> None:
    table = alphagenome_archives.TABLES[name]
    encoded = json.loads(alphagenome_archives.config_text(table))
    assert [field['field'] for field in encoded] == list(table.kept)
    assert set(table.zero) < set(table.columns)
    for field in encoded:
        # A score can be negative, and only a missing float is no score: any number as the stand-in could be one.
        assert field['zigzag'] is True
        assert field['missing_value'] == alphagenome_archives.FLOAT_MISSING
        # echtvar reads the field's type from the VCF, which declares each as a float.
        assert 'ftype' not in field


def test_no_two_fields_of_the_archives_share_an_alias() -> None:
    # One `echtvar anno` call reads every archive, and a later archive's field would replace one of the same alias.
    aliases = [
        field['alias']
        for table in alphagenome_archives.TABLES.values()
        for field in alphagenome_archives.config(table).values()
    ]
    assert len(set(aliases)) == len(aliases)


def test_the_tables_that_are_not_for_permissive_use_carry_the_notice_in_every_field() -> None:
    for table in (BREAKDOWN, SPLICING):
        for field in alphagenome_archives.config(table).values():
            assert field['description'].endswith(alphagenome_archives.NOTICE)
    assert all('description' not in field for field in alphagenome_archives.config(AVI).values())


# ---- the conversion ----


def test_each_row_becomes_a_record_with_a_field_of_each_kept_column(tmp_path: pathlib.Path) -> None:
    header, records = _read(tmp_path, _vcf(AVI, _AVI_ROWS))
    assert [variant for variant, _info in records] == ['chr1-5-C-A', 'chr1-5-C-G', 'chr1-5-C-T', 'chrX-7-G-A']
    assert [set(info) for _variant, info in records] == [set(AVI.kept)] * len(_AVI_ROWS)
    assert records[0][1] == {'raw_score': pytest.approx(-0.04231), 'PHRED': pytest.approx(1.20713)}
    # A value of 0 is written as a value, since a record without the field would read as an allele without a score.
    assert records[1][1] == {'raw_score': 0.0, 'PHRED': 0.0}
    assert {name: (header.info[name].number, header.info[name].type) for name in AVI.kept} == dict.fromkeys(
        AVI.kept, (1, 'Float')
    )
    assert {name: contig.length for name, contig in header.contigs.items()} == _CONTIGS


def test_a_value_is_passed_on_as_the_table_writes_it() -> None:
    text = _vcf(BREAKDOWN, _BREAKDOWN_ROWS)
    record = next(line for line in text.splitlines() if line.startswith('chr1\t'))
    assert 'MERGED_SPLICING=0.2817;MAX_ABS_ATAC=0.001293;MAX_ABS_CONTACT_MAPS=5.318e-05;' in record
    assert 'MAX_ABS_CHIP_TF=-0.0004126;' in record


def test_the_columns_that_hold_zero_are_left_out(tmp_path: pathlib.Path) -> None:
    header, records = _read(tmp_path, _vcf(BREAKDOWN, _BREAKDOWN_ROWS))
    assert len(BREAKDOWN.kept) == len(BREAKDOWN.columns) - len(BREAKDOWN.zero) == 16
    assert all(set(info) == set(BREAKDOWN.kept) for _variant, info in records)
    assert not set(BREAKDOWN.zero) & set(header.info)


def test_the_vcf_of_a_table_with_a_notice_carries_it_in_every_field(tmp_path: pathlib.Path) -> None:
    header, _records = _read(tmp_path, _vcf(SPLICING, ['chr1\t5\tC\tA\t0.002873']))
    description = header.info['alphagenome_splicing'].description
    assert description is not None
    assert description.endswith(alphagenome_archives.NOTICE)
    header, _records = _read(tmp_path, _vcf(AVI, _AVI_ROWS))
    assert 'Terms of Use' not in str(header)


def _with(row: str, index: int, value: str) -> str:
    cells = row.split('\t')
    cells[index] = value
    return '\t'.join(cells)


ROWS_REFUSED = {
    'a row with a cell too few': (AVI, 'chr1\t5\tC\tA\t0.5', '5 cells, and the header has 6'),
    'a row with a cell too many': (AVI, 'chr1\t5\tC\tA\t0.5\t1.5\t2.5', '7 cells, and the header has 6'),
    'an empty line': (AVI, '', '1 cells'),
    'a value that is no number': (AVI, 'chr1\t5\tC\tA\tNA\t1.5', "raw_score is 'NA', not a number"),
    'an empty value': (AVI, 'chr1\t5\tC\tA\t0.5\t', "PHRED is '', not a number"),
    'a line that ends with a carriage return': (AVI, 'chr1\t5\tC\tA\t0.5\t1.5\r', 'PHRED is'),
    'an insertion': (AVI, 'chr1\t5\tC\tCA\t0.5\t1.5', "the alleles 'C' and 'CA' are not two single bases"),
    'a base in lower case': (AVI, 'chr1\t5\tc\tA\t0.5\t1.5', "the alleles 'c' and 'A' are not two single bases"),
    'a position that is no number': (AVI, 'chr1\tfive\tC\tA\t0.5\t1.5', 'no chromosome and position'),
    'a position of 0': (AVI, 'chr1\t0\tC\tA\t0.5\t1.5', 'no chromosome and position'),
    'a chromosome the reference lacks': (AVI, 'chr3\t5\tC\tA\t0.5\t1.5', 'chr3 is none of chr1, chr2, chrX'),
    'a value too large for the archive': (
        AVI,
        'chr1\t5\tC\tA\t0.5\t2148.5',
        'PHRED is 2148.5, which echtvar cannot store with the multiplier 1000000',
    ),
    'a value with a positive exponent': (
        SPLICING,
        'chr1\t5\tC\tA\t1.5e+03',
        r'alphagenome_splicing is 1\.5e\+03, which echtvar cannot store',
    ),
    'a contribution of an insertion': (
        BREAKDOWN,
        _with(_BREAKDOWN_ROWS[0], -2, '0.0125'),
        'IS_INSERTION is 0.0125, not 0',
    ),
    'a contribution of a deletion': (
        BREAKDOWN,
        _with(_BREAKDOWN_ROWS[0], -1, '-1e-07'),
        'IS_DELETION is -1e-07, not 0',
    ),
}


@pytest.mark.parametrize('case', ROWS_REFUSED)
def test_a_row_that_is_not_one_of_the_table_is_refused(case: str) -> None:
    table, row, message = ROWS_REFUSED[case]
    good = _AVI_ROWS[0] if table is AVI else _BREAKDOWN_ROWS[0] if table is BREAKDOWN else 'chr1\t4\tC\tA\t0.5'
    with pytest.raises(ValueError, match=f'row 2 of {table.name}: .*{message}'):
        _vcf(table, [good, row])


HEADERS_REFUSED = {
    'a header of another table': ('\t'.join(('#CHROM', 'POS', 'REF', 'ALT', *SPLICING.columns)) + '\n', 'starts with'),
    'a header with a column more': ('#CHROM\tPOS\tREF\tALT\traw_score\tPHRED\tCADD\n', 'starts with'),
    'a header with its columns in another order': ('#CHROM\tPOS\tREF\tALT\tPHRED\traw_score\n', 'starts with'),
    'a header without its number sign': ('CHROM\tPOS\tREF\tALT\traw_score\tPHRED\n', 'starts with'),
    'a table that starts with a row': (_AVI_ROWS[0] + '\n', 'starts with'),
}


@pytest.mark.parametrize('case', HEADERS_REFUSED)
def test_a_header_that_is_not_the_tables_is_refused(case: str) -> None:
    header, message = HEADERS_REFUSED[case]
    with pytest.raises(ValueError, match=f'the header of alphagenome_avi is .*{message}'):
        alphagenome_archives.to_vcf(AVI, [header, _AVI_ROWS[0] + '\n'], _CONTIGS, io.StringIO())


def test_a_table_without_a_line_or_without_a_row_is_refused() -> None:
    with pytest.raises(ValueError, match='the table starts with no line'):
        alphagenome_archives.to_vcf(AVI, [], _CONTIGS, io.StringIO())
    with pytest.raises(ValueError, match='a header and no row'):
        _vcf(AVI, [])


def test_a_table_with_another_header_writes_nothing() -> None:
    out = io.StringIO()
    with pytest.raises(ValueError, match='the header of alphagenome_avi'):
        alphagenome_archives.to_vcf(AVI, _lines(SPLICING, ['chr1\t5\tC\tA\t0.5']), _CONTIGS, out)
    assert out.getvalue() == ''


def test_a_field_stores_the_values_whose_product_with_its_multiplier_fits_in_32_bits() -> None:
    millionths = re.compile(alphagenome_archives._storable(1_000_000))
    assert millionths.fullmatch('999.99999')
    assert millionths.fullmatch('-2.5e-07')
    assert not millionths.fullmatch('1000.0')
    with pytest.raises(ValueError, match='leaves no value that echtvar can store'):
        alphagenome_archives._storable(1_000_000_000)


def test_the_contigs_are_the_primary_chromosomes_of_the_reference(tmp_path: pathlib.Path) -> None:
    index = tmp_path / 'reference.fa.fai'
    index.write_text('chr2\t900\t0\t60\t61\nchr1\t1000\t9\t60\t61\nchr1_KI270706v1_random\t70\t9\t60\t61\n')
    assert list(alphagenome_archives.read_contigs(index).items()) == [('chr1', 1000), ('chr2', 900)]
    index.write_text('1\t1000\t9\t60\t61\n2\t900\t9\t60\t61\n')
    with pytest.raises(ValueError, match='no primary chromosome'):
        alphagenome_archives.read_contigs(index)


def _indexed(tmp_path: pathlib.Path, table: alphagenome_archives.Table, rows: Sequence[str]) -> pathlib.Path:
    """The table as DeepMind publishes it: bgzipped, with a tabix index beside it."""
    path = tmp_path / f'{table.name}.tsv'
    path.write_text(''.join(_lines(table, rows)))
    return pathlib.Path(pysam.tabix_index(str(path), seq_col=0, start_col=1, end_col=1, meta_char='#', force=True))


def test_an_indexed_table_gives_the_rows_of_one_chromosome_or_of_all(tmp_path: pathlib.Path) -> None:
    tsv = _indexed(tmp_path, AVI, _AVI_ROWS)
    whole = _lines(AVI, _AVI_ROWS)
    assert [line.rstrip('\n') for line in alphagenome_archives.indexed_lines(tsv, None)] == [
        line.rstrip('\n') for line in whole
    ]
    one = io.StringIO()
    assert alphagenome_archives.to_vcf(AVI, alphagenome_archives.indexed_lines(tsv, 'chrX'), _CONTIGS, one) == 1
    assert one.getvalue() == _vcf(AVI, _AVI_ROWS[3:])
    with pytest.raises(ValueError, match='no row on chr2; it has rows on chr1, chrX'):
        list(alphagenome_archives.indexed_lines(tsv, 'chr2'))


# ---- the commands ----


def _reference_index(tmp_path: pathlib.Path) -> pathlib.Path:
    index = tmp_path / 'reference.fa.fai'
    index.write_text(''.join(f'{name}\t{length}\t0\t60\t61\n' for name, length in _CONTIGS.items()))
    return index


def test_the_command_converts_a_table_from_standard_input_or_from_an_indexed_file(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    arguments = ['alphagenome-vcf', '--table', AVI.name, '--reference-index', str(_reference_index(tmp_path))]
    monkeypatch.setattr('sys.stdin', io.StringIO(''.join(_lines(AVI, _AVI_ROWS))))
    cli.main(arguments)
    assert capsys.readouterr().out == _vcf(AVI, _AVI_ROWS)
    cli.main([*arguments, '--tsv', str(_indexed(tmp_path, AVI, _AVI_ROWS)), '--chromosome', 'chr1'])
    assert capsys.readouterr().out == _vcf(AVI, _AVI_ROWS[:3])
    with pytest.raises(SystemExit, match='--chromosome picks the rows of an indexed table'):
        cli.main([*arguments, '--chromosome', 'chr1'])


@pytest.mark.parametrize('name', alphagenome_archives.TABLES)
def test_the_command_writes_the_configuration_of_a_table(name: str, capsys: pytest.CaptureFixture[str]) -> None:
    cli.main(['alphagenome-config', '--table', name])
    written = capsys.readouterr().out
    assert written == alphagenome_archives.config_text(alphagenome_archives.TABLES[name])
    assert [field['field'] for field in json.loads(written)] == list(alphagenome_archives.TABLES[name].kept)


def test_the_commands_refuse_a_table_that_is_none_of_the_three() -> None:
    with pytest.raises(SystemExit):
        cli.main(['alphagenome-config', '--table', 'alphagenome'])


# ---- the check of an archive ----

_ALIASES = ('avi_raw', 'avi_phred')
# What echtvar returns for the invented rows from an archive of them: each value cut to a whole number of millionths,
# which is the value or one step below it, as a 32-bit float.
_RETURNED = (
    'avi_raw=-0.042309;avi_phred=1.207129',
    'avi_raw=0;avi_phred=0',
    'avi_raw=2.71306;avi_phred=48.902149',
    'avi_raw=-1.137;avi_phred=97.300003',
)


def _annotated(
    tmp_path: pathlib.Path,
    rows: Sequence[str] = _AVI_ROWS,
    returned: Sequence[str] = _RETURNED,
    aliases: Sequence[str] = _ALIASES,
) -> pathlib.Path:
    """The rows' alleles as echtvar annotates them: without the fields of the conversion, with the archive's."""
    lines = ['##fileformat=VCFv4.2', *(f'##contig=<ID={name},length={length}>' for name, length in _CONTIGS.items())]
    lines += [f'##INFO=<ID={alias},Number=1,Type=Float,Description="added by echtvar">' for alias in aliases]
    lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
    for row, info in zip(rows, returned, strict=True):
        chrom, pos, ref, alt, *_ = row.split('\t')
        lines.append(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}')
    path = tmp_path / 'annotated.vcf'
    path.write_text('\n'.join(lines) + '\n')
    return path


def test_an_archive_that_returns_the_rows_passes_the_check(tmp_path: pathlib.Path) -> None:
    assert alphagenome_archives.check(AVI, _lines(AVI, _AVI_ROWS), _annotated(tmp_path)) == len(_AVI_ROWS)


def test_the_command_checks_the_rows_on_standard_input(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr('sys.stdin', io.StringIO(''.join(_lines(AVI, _AVI_ROWS))))
    cli.main(['alphagenome-check', '--table', AVI.name, '--annotated', str(_annotated(tmp_path))])
    assert f'the {len(_AVI_ROWS)} rows of alphagenome_avi' in capsys.readouterr().out


def _returned(index: int, info: str) -> list[str]:
    return [info if at == index else kept for at, kept in enumerate(_RETURNED)]


ARCHIVES_FAILING = {
    'a value two steps below the row': (
        {'returned': _returned(0, 'avi_raw=-0.042313;avi_phred=1.207129')},
        r'chr1-5-C-A: avi_raw is -0\.04231\d+, and raw_score is -0\.04231',
    ),
    'a large value further off than a 32-bit float can account for': (
        {'returned': _returned(2, 'avi_raw=2.71306;avi_phred=48.90213')},
        r'chr1-5-C-T: avi_phred is 48\.9021\d+, and PHRED is 48\.90215',
    ),
    'an allele that the archive lacks': (
        {'returned': _returned(1, 'avi_raw=.;avi_phred=.')},
        'chr1-5-C-G: no avi_raw, and raw_score is 0.0; chr1-5-C-G: no avi_phred, and PHRED is 0.0',
    ),
    'the two fields exchanged': (
        {'returned': _returned(0, 'avi_raw=1.207129;avi_phred=-0.042309')},
        'chr1-5-C-A: avi_raw is 1.20',
    ),
    'an annotated VCF without a field of the archive': (
        {'aliases': _ALIASES[:1], 'returned': ['avi_raw=0'] * len(_AVI_ROWS)},
        'no INFO field avi_phred, which the archive adds',
    ),
    'fewer records than rows': (
        {'rows': _AVI_ROWS[:3], 'returned': _RETURNED[:3]},
        'the table holds alleles after the 3 that both hold',
    ),
    'records in another order than the rows': (
        {'rows': [_AVI_ROWS[1], _AVI_ROWS[0], *_AVI_ROWS[2:]], 'returned': _RETURNED},
        'record 1 is not chr1-5-C-A, the row at its place',
    ),
}


@pytest.mark.parametrize('case', ARCHIVES_FAILING)
def test_an_archive_that_does_not_return_the_rows_fails_the_check(tmp_path: pathlib.Path, case: str) -> None:
    changes, message = ARCHIVES_FAILING[case]
    with pytest.raises(ValueError, match=message):
        alphagenome_archives.check(AVI, _lines(AVI, _AVI_ROWS), _annotated(tmp_path, **changes))


def test_more_records_than_rows_fail_the_check(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='the annotated VCF holds alleles after the 3 that both hold'):
        alphagenome_archives.check(AVI, _lines(AVI, _AVI_ROWS[:3]), _annotated(tmp_path))


def test_a_check_of_no_row_fails(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='a header and no row'):
        alphagenome_archives.check(AVI, _lines(AVI, []), _annotated(tmp_path, rows=[], returned=[]))
