"""The preparation of the gene models that the transcript rules read, on the workflow test's invented gene models.

tests/workflows/family_tables/gene-models/ holds one invented gene in the layout of Ensembl's GFF3 and in that of NCBI's
RefSeq GFF, over the invented reference beside it, and an invented assembly report.
"""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import pathlib
import re
from collections.abc import Callable, Sequence

import pysam
import pytest

from case_annotation import cli, family_tables, gene_models

FIXTURE = pathlib.Path(__file__).parent / 'workflows' / 'family_tables'
REFERENCE = FIXTURE / 'reference.fa'
REPORT = FIXTURE / 'gene-models' / 'assembly-report.txt'
ENSEMBL = (FIXTURE / 'gene-models' / 'ensembl.gff3').read_text('utf-8')
REFSEQ = (FIXTURE / 'gene-models' / 'refseq.gff').read_text('utf-8')


def _gzipped(path: pathlib.Path, text: str) -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='utf-8') as f:
        f.write(text)
    return path


def _prepare(
    tmp_path: pathlib.Path,
    ensembl: str = ENSEMBL,
    refseq: str = REFSEQ,
    report: pathlib.Path = REPORT,
    reference: pathlib.Path = REFERENCE,
) -> gene_models.GeneModels:
    public: dict[gene_models.Source, pathlib.Path] = {
        'Ensembl': _gzipped(tmp_path / 'invented.gff3.gz', ensembl),
        'RefSeq': _gzipped(tmp_path / 'invented.gff.gz', refseq),
    }
    return gene_models.prepare(public, report, reference, tmp_path / 'gene-models')


def _features(text: str) -> list[list[str]]:
    return [line.split('\t') for line in text.splitlines() if line and not line.startswith('#')]


def _edited(text: str, edit: Callable[[list[str]], list[str] | None]) -> str:
    """A GFF text with each feature line as `edit` returns it, and without a line for which it returns None."""
    lines = []
    for line in text.splitlines():
        if line.startswith('#') or not line:
            lines.append(line)
        elif (fields := edit(line.split('\t'))) is not None:
            lines.append('\t'.join(fields))
    return '\n'.join(lines) + '\n'


def test_the_chromosomes_are_the_ones_of_the_tables() -> None:
    assert gene_models.CHROMOSOMES == family_tables.PRIMARY


@pytest.mark.parametrize(('source', 'text'), [('Ensembl', ENSEMBL), ('RefSeq', REFSEQ)])
def test_a_prepared_file_holds_every_feature_line_sorted_within_its_sequence(
    tmp_path: pathlib.Path, source: gene_models.Source, text: str
) -> None:
    # The lines of a second sequence follow, out of order within it, as a public file has them.
    second = '2' if source == 'Ensembl' else 'INV_000002.1'
    more = ''.join(
        '\t'.join([second, 'invented', 'exon', start, end, '.', '+', '.', 'Parent=none']) + '\n'
        for start, end in (('50', '90'), ('10', '40'), ('10', '20'))
    )
    _prepare(tmp_path, **{source.lower(): text + more})  # type: ignore[arg-type]
    prepared = tmp_path / 'gene-models' / gene_models.FILES[source]
    with gzip.open(prepared, 'rt') as f:
        written = [line.rstrip('\n').split('\t') for line in f]
    assert sorted(written) == sorted(_features(text + more))
    by_sequence: dict[str, list[tuple[int, int]]] = {}
    for fields in written:
        by_sequence.setdefault(fields[0], []).append((int(fields[3]), int(fields[4])))
    assert len(by_sequence) == 2
    assert all(places == sorted(places) for places in by_sequence.values())
    # tabix reads the file by region, under the sequence's name in the file.
    with pysam.TabixFile(str(prepared)) as indexed:
        assert sorted(indexed.contigs) == sorted(by_sequence)
        assert len(list(indexed.fetch(second, 0, 30))) == 2


def test_the_synonyms_pair_every_two_names_of_each_chromosome(tmp_path: pathlib.Path) -> None:
    _prepare(tmp_path)
    pairs = [
        tuple(line.split('\t'))
        for line in (tmp_path / 'gene-models' / gene_models.SYNONYMS_NAME).read_text().splitlines()
    ]
    chromosomes = gene_models.read_assembly_report(REPORT)
    assert [c.ucsc for c in chromosomes] == list(gene_models.CHROMOSOMES)
    assert len(pairs) == len(set(pairs)) == 3 * len(chromosomes)
    for chromosome in chromosomes:
        name, ucsc, refseq = chromosome.names
        assert {(name, ucsc), (name, refseq), (ucsc, refseq)} <= set(pairs)
    # A sequence that is no assembled chromosome has no synonym.
    assert not any('chrUn_test' in pair for pair in pairs)


def test_the_regions_hold_every_exon_of_both_files_with_a_base_on_each_side(tmp_path: pathlib.Path) -> None:
    # An exon of RefSeq's alone, which overlaps one of both files, and a feature that is no exon.
    refseq = REFSEQ + ''.join(
        '\t'.join(['INV_000001.1', 'invented', kind, start, end, '.', '+', '.', 'Parent=none']) + '\n'
        for kind, start, end in (('exon', '250', '300'), ('exon', '20', '30'), ('region', '590', '600'))
    )
    _prepare(tmp_path, refseq=refseq)
    regions = [
        (chrom, int(start), int(end))
        for chrom, start, end in (
            line.split('\t') for line in (tmp_path / 'gene-models' / gene_models.EXONS_NAME).read_text().splitlines()
        )
    ]
    assert {chrom for chrom, _start, _end in regions} == {'chr1'}
    # In order and apart, so that bcftools reads no record for two of them.
    assert regions == sorted(regions)
    assert all(before[2] < after[1] for before, after in itertools.pairwise(regions))
    exons = [(int(f[3]), int(f[4])) for text in (ENSEMBL, refseq) for f in _features(text) if f[2] == 'exon']
    assert len(exons) == 10
    for start, end in exons:
        # BED counts from 0 and leaves the end out: the base before the exon is start - 2, the one after it end.
        assert any(first <= start - 2 and end + 1 <= last for _chrom, first, last in regions), (start, end)
    covered = sum(end - start for _chrom, start, end in regions)
    assert covered == len({p for start, end in exons for p in range(start - 1, end + 2)})


def _first_exon_and_cds(text: str) -> tuple[tuple[int, int], tuple[int, int]]:
    features = _features(text)
    exons = sorted((int(f[3]), int(f[4])) for f in features if f[2] == 'exon')
    cds = sorted((int(f[3]), int(f[4])) for f in features if f[2] == 'CDS')
    return exons[0], cds[0]


@pytest.mark.parametrize(
    ('source', 'text', 'transcript'), [('Ensembl', ENSEMBL, 'ENST01.1'), ('RefSeq', REFSEQ, 'NM_01.1')]
)
def test_the_controls_are_insertions_into_the_utr_and_the_coding_start_of_a_transcript(
    tmp_path: pathlib.Path, source: gene_models.Source, text: str, transcript: str
) -> None:
    record = _prepare(tmp_path)
    (file,) = [f for f in record.gene_models if f.source == source]
    utr, nmd = file.controls
    assert (utr.plugin, nmd.plugin) == ('UTRannotator', 'NMD')
    # The id that VEP writes with --transcript_version: the id with the version that Ensembl's file gives beside it,
    # and RefSeq's id, which ends in its version.
    assert utr.transcript_id == nmd.transcript_id == transcript
    first_exon, first_cds = _first_exon_and_cds(text)
    with pysam.FastaFile(str(REFERENCE)) as reference:
        for control in (utr, nmd):
            assert control.chrom == 'chr1'
            assert control.ref == reference.fetch('chr1', control.pos - 1, control.pos)
            assert control.alt.startswith(control.ref)
    # ATG after a base of the first exon, which lies wholly before the coding sequence.
    assert first_exon[0] < utr.pos < min(first_exon[1], first_cds[0] - 1)
    assert utr.alt == utr.ref + 'ATG'
    # One base within the first 100 coding bases.
    assert first_cds[0] < nmd.pos < min(first_cds[1], first_cds[0] + 99)
    assert len(nmd.alt) == len(nmd.ref) + 1 == 2
    assert len({c.id for f in record.gene_models for c in f.controls}) == 4


def test_the_controls_file_holds_each_control_as_a_record_under_its_id(tmp_path: pathlib.Path) -> None:
    record = _prepare(tmp_path)
    with pysam.VariantFile(str(tmp_path / 'gene-models' / gene_models.CONTROLS_NAME)) as vcf:
        written = [(r.id, r.chrom, r.pos, r.ref, r.alts) for r in vcf]
    controls = [c for f in record.gene_models for c in f.controls]
    assert written == [(c.id, c.chrom, c.pos, c.ref, (c.alt,)) for c in controls]
    # No id is a variant id, which begins with its chromosome.
    assert not any(c.id.startswith('chr') for c in controls)


def test_the_record_names_each_public_file_by_its_checksum_and_what_it_states(tmp_path: pathlib.Path) -> None:
    record = _prepare(tmp_path)
    assert record == gene_models.read_record(tmp_path / 'gene-models')
    assert [f.source for f in record.gene_models] == list(gene_models.FILES)
    given = (('invented.gff3.gz', ENSEMBL), ('invented.gff.gz', REFSEQ))
    for file, (public, text) in zip(record.gene_models, given, strict=True):
        content = (tmp_path / public).read_bytes()
        assert file.file == gene_models.FILES[file.source]
        prepared = (tmp_path / 'gene-models' / file.file).read_bytes()
        assert (file.bytes, file.sha256) == (len(prepared), hashlib.sha256(prepared).hexdigest())
        assert (file.prepared_from.name, file.prepared_from.bytes) == (public, len(content))
        assert file.prepared_from.sha256 == hashlib.sha256(content).hexdigest()
        assert list(file.prepared_from.stated) == [line for line in text.splitlines() if line.startswith('#!')]
        assert (tmp_path / 'gene-models' / f'{file.file}.tbi').is_file()
    report = REPORT.read_bytes()
    assert record.assembly_report == gene_models.PublicFile(
        name=REPORT.name, bytes=len(report), sha256=hashlib.sha256(report).hexdigest(), stated=()
    )
    # Nothing is left of the unsorted text.
    names = {path.name for path in (tmp_path / 'gene-models').iterdir()}
    assert names == {
        *gene_models.FILES.values(),
        *(f'{name}.tbi' for name in gene_models.FILES.values()),
        gene_models.EXONS_NAME,
        gene_models.SYNONYMS_NAME,
        gene_models.CONTROLS_NAME,
        gene_models.RECORD_NAME,
    }


def test_the_command_writes_the_folder(tmp_path: pathlib.Path) -> None:
    ensembl, refseq = _gzipped(tmp_path / 'e.gff3.gz', ENSEMBL), _gzipped(tmp_path / 'r.gff.gz', REFSEQ)
    arguments = ['--ensembl', str(ensembl), '--refseq', str(refseq), '--assembly-report', str(REPORT)]
    cli.main(['gene-models', *arguments, '--reference', str(REFERENCE), '--output', str(tmp_path / 'out')])
    record = json.loads((tmp_path / 'out' / gene_models.RECORD_NAME).read_text())
    assert [f['prepared_from']['name'] for f in record['gene_models']] == ['e.gff3.gz', 'r.gff.gz']


def test_a_control_skips_a_transcript_that_the_file_flags_or_that_has_no_room(tmp_path: pathlib.Path) -> None:
    # Four transcripts before the invented one, each unfit in one way: marked partial, without a start by its tag, on
    # the reverse strand, and with a first exon of 20 bases.
    def transcript(name: str, strand: str, exon: tuple[int, int], extra: str) -> str:
        rows = (
            ('mRNA', 5, 60, f'ID=transcript:{name};transcript_id={name};version=1{extra}'),
            ('exon', *exon, f'Parent=transcript:{name}'),
            ('CDS', 41, 60, f'Parent=transcript:{name}'),
        )
        return ''.join(
            '\t'.join(['1', 'invented', kind, str(start), str(end), '.', strand, '.', attributes]) + '\n'
            for kind, start, end, attributes in rows
        )

    unfit = (
        transcript('ENST91', '+', (5, 60), ';partial=true')
        + transcript('ENST92', '+', (5, 60), ';tag=cds_start_NF')
        + transcript('ENST93', '-', (5, 60), '')
        + transcript('ENST94', '+', (21, 60), '')
    )
    head, body = ENSEMBL.split('\n1\t', 1)
    record = _prepare(tmp_path, ensembl=f'{head}\n{unfit}1\t{body}')
    assert {c.transcript_id for c in record.gene_models[0].controls} == {'ENST01.1'}


def test_a_transcript_id_that_ends_in_a_version_gets_no_second_one(tmp_path: pathlib.Path) -> None:
    versioned = ENSEMBL.replace('transcript_id=ENST01;version=1', 'transcript_id=ENST01.7;version=1')
    record = _prepare(tmp_path, ensembl=versioned)
    assert {c.transcript_id for c in record.gene_models[0].controls} == {'ENST01.7'}


def _report(tmp_path: pathlib.Path, edit: Callable[[list[str]], Sequence[str]]) -> pathlib.Path:
    path = tmp_path / 'report.txt'
    path.write_text(''.join(line + '\n' for row in REPORT.read_text().splitlines() for line in edit(row.split('\t'))))
    return path


def test_an_assembly_report_with_a_short_line_is_refused(tmp_path: pathlib.Path) -> None:
    short = _report(tmp_path, lambda row: ['\t'.join(row[:9] if row[0] == '2' else row)])
    with pytest.raises(ValueError, match='line 5 has 9 columns'):
        gene_models.read_assembly_report(short)


def test_an_assembly_report_that_lists_a_chromosome_twice_is_refused(tmp_path: pathlib.Path) -> None:
    twice = _report(tmp_path, lambda row: ['\t'.join(row)] * (2 if row[0] == 'X' else 1))
    with pytest.raises(ValueError, match='chrX is listed twice'):
        gene_models.read_assembly_report(twice)


def test_an_assembly_report_without_a_chromosome_of_the_tables_is_refused(tmp_path: pathlib.Path) -> None:
    without = _report(tmp_path, lambda row: [] if row[0] == 'MT' else ['\t'.join(row)])
    with pytest.raises(ValueError, match=re.escape("lack ['chrM'] and hold []")):
        gene_models.read_assembly_report(without)


def test_an_assembly_report_with_another_assembled_molecule_is_refused(tmp_path: pathlib.Path) -> None:
    def assembled(row: list[str]) -> list[str]:
        return ['\t'.join(['assembled-molecule' if field == 'unplaced-scaffold' else field for field in row])]

    with pytest.raises(ValueError, match=re.escape("lack [] and hold ['chrUn_test']")):
        gene_models.read_assembly_report(_report(tmp_path, assembled))


def test_a_feature_line_without_nine_columns_is_refused(tmp_path: pathlib.Path) -> None:
    short = _edited(ENSEMBL, lambda fields: fields[:8] if fields[2] == 'gene' else fields)
    with pytest.raises(ValueError, match='line 5 is no GFF feature line'):
        _prepare(tmp_path, ensembl=short)


def test_a_feature_line_whose_start_is_no_number_is_refused(tmp_path: pathlib.Path) -> None:
    unplaced = _edited(REFSEQ, lambda fields: [*fields[:3], '.', *fields[4:]] if fields[2] == 'gene' else fields)
    with pytest.raises(ValueError, match='line 7 is no GFF feature line'):
        _prepare(tmp_path, refseq=unplaced)


def test_a_file_whose_sequences_are_not_together_is_refused(tmp_path: pathlib.Path) -> None:
    other = '\t'.join(['2', 'invented', 'exon', '10', '20', '.', '+', '.', 'Parent=none']) + '\n'
    again = '\t'.join(['1', 'invented', 'exon', '10', '20', '.', '+', '.', 'Parent=none']) + '\n'
    with pytest.raises(ValueError, match='the lines of 1 are not together'):
        _prepare(tmp_path, ensembl=ENSEMBL + other + again)


def test_a_file_without_a_transcript_for_the_controls_is_refused(tmp_path: pathlib.Path) -> None:
    reverse = _edited(ENSEMBL, lambda fields: [*fields[:6], '-', *fields[7:]])
    with pytest.raises(ValueError, match='no mRNA of the forward strand on a chromosome has room'):
        _prepare(tmp_path, ensembl=reverse)


def test_a_file_under_other_sequence_names_is_refused(tmp_path: pathlib.Path) -> None:
    # A file whose sequence names are not in the assembly report has no transcript on a chromosome, so the
    # preparation finds no candidate for a control and stops.
    with pytest.raises(ValueError, match=r'invented\.gff\.gz: no mRNA of the forward strand on a chromosome'):
        _prepare(tmp_path, refseq=ENSEMBL.replace('\n1\t', '\nsequence-1\t'))


def test_a_reference_under_other_sequence_names_is_refused(tmp_path: pathlib.Path) -> None:
    renamed = tmp_path / 'renamed.fa'
    renamed.write_text(REFERENCE.read_text().replace('>chr', '>'))
    pysam.faidx(str(renamed))
    with pytest.raises(ValueError, match='the reference has no sequence named chr1'):
        _prepare(tmp_path, reference=renamed)


def test_a_control_is_not_written_after_a_base_that_the_reference_does_not_give(tmp_path: pathlib.Path) -> None:
    lines = REFERENCE.read_text().splitlines()
    # The reference's first sequence with N at the place of the UTR's control, 84.
    lines[2] = lines[2][:23] + 'N' + lines[2][24:]
    masked = tmp_path / 'masked.fa'
    masked.write_text('\n'.join(lines) + '\n')
    pysam.faidx(str(masked))
    with pysam.FastaFile(str(masked)) as reference:
        assert reference.fetch('chr1', 83, 84) == 'N'
    with pytest.raises(ValueError, match='no mRNA of the forward strand on a chromosome has room'):
        _prepare(tmp_path, reference=masked)


def test_the_two_files_are_both_given(tmp_path: pathlib.Path) -> None:
    only = {'Ensembl': _gzipped(tmp_path / 'invented.gff3.gz', ENSEMBL)}
    with pytest.raises(ValueError, match=re.escape("the gene model files are ['Ensembl', 'RefSeq']")):
        gene_models.prepare(only, REPORT, REFERENCE, tmp_path / 'gene-models')  # type: ignore[arg-type]
