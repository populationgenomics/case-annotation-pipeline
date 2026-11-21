"""The `case-annotation` command line: one subcommand per pipeline step."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from collections.abc import Callable, Sequence

from case_annotation import (
    alphagenome_archives,
    annotation_tables,
    branch_points,
    contract,
    echtvar_values,
    family_tables,
    gene_models,
    joined_values,
    roh_upd,
    strchive,
    sv_annotation,
    sv_sites,
    transcript_rules,
)


def _proband_arguments(parser: argparse.ArgumentParser) -> None:
    """The family's proband as the run's sheet states it: one of the two is required, and nothing derives a proband."""
    stated = parser.add_mutually_exclusive_group(required=True)
    stated.add_argument(
        '--proband',
        help=(
            "the sample id of the family's proband, which must be an individual of the PED file and a sample of the "
            'VCF; its phenotype is not checked'
        ),
    )
    stated.add_argument(
        '--no-proband',
        action='store_true',
        help='no proband is designated in the family, so no member has a role in the manifest',
    )


def _structural_arguments(parser: argparse.ArgumentParser) -> None:
    """The VCFs of a family's larger variants as the run's sheet names them, each of them optional."""
    path = pathlib.Path
    parser.add_argument(
        '--sv-vcf',
        type=path,
        help="the VCF of the family's structural-variant caller, whose samples are the family's, for the sv table",
    )
    parser.add_argument(
        '--cnv-vcf',
        type=path,
        help="the VCF of the family's copy-number caller, whose samples are the family's, for the sv table",
    )
    parser.add_argument(
        '--repeat-vcf',
        type=path,
        action='append',
        default=[],
        help="one member's repeat VCF, for repeats; given once for every sample of the family's VCF, or not at all",
    )


def _structural(args: argparse.Namespace) -> family_tables.Structural:
    return family_tables.Structural(sv_vcf=args.sv_vcf, cnv_vcf=args.cnv_vcf, repeat_vcfs=tuple(args.repeat_vcf))


def _family_check_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument('--family', required=True, help="the family's id, which its PED file must hold")
    parser.add_argument('--ped', type=path, required=True, help="the family's PED file")
    parser.add_argument(
        '--vcf',
        type=path,
        required=True,
        help=(
            "the family's VCF, whose samples the PED file names and whose sample ids name the members' columns of the "
            'variants table, and of sv and repeats'
        ),
    )
    _proband_arguments(parser)
    _structural_arguments(parser)
    _release_argument(parser)


def _release_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--pipeline-version',
        required=True,
        help=(
            "the version of the pipeline whose workflow runs this command, which a family's manifest records as the "
            "release that wrote its folder; another than this package's fails"
        ),
    )


def _require_release(given: str) -> None:
    """Fails unless the workflow that runs a command is of this package's release of the pipeline.

    A run takes its tools from an image, which a params file names apart from the workflow. The first command of a
    run checks it, before anything is annotated, and so does the command that writes the release into a manifest.
    """
    if given != family_tables.PIPELINE_VERSION:
        raise SystemExit(
            f'these tools are of version {family_tables.PIPELINE_VERSION} of the pipeline, and the workflow of '
            f"{given}: a run takes the tools image of the workflow's release"
        )


def _family_check(args: argparse.Namespace) -> None:
    _require_release(args.pipeline_version)
    family_tables.check_inputs(args.family, args.ped, args.vcf, args.proband, _structural(args))


def _annotation_tables_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument('--json', type=path, required=True, help="the Annotator's JSON of the union, gzipped")
    parser.add_argument('--sites', type=path, required=True, help='the union, a sites-only VCF the JSON annotates')
    parser.add_argument('--run-config', type=path, required=True, help='the run config the Annotator ran with')
    parser.add_argument(
        '--bundle-config', type=path, required=True, help="the bundle's run config, runConfigs/bundle.json"
    )
    parser.add_argument('--catalog', required=True, help="the catalog's id")
    parser.add_argument('--workers', type=int, required=True, help='processes converting the JSON')
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help='the folder to write variants/, transcripts/ and annotation.json into; it must not exist',
    )


def _annotation_tables(args: argparse.Namespace) -> None:
    info = annotation_tables.convert(
        json_path=args.json,
        sites_path=args.sites,
        run_config=json.loads(args.run_config.read_text()),
        bundle_config=json.loads(args.bundle_config.read_text()),
        catalog=args.catalog,
        output=args.output,
        workers=args.workers,
    )
    (args.output / 'annotation.json').write_text(info.model_dump_json(indent=2) + '\n')


def _sv_union_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument(
        '--sites',
        nargs=2,
        action='append',
        required=True,
        metavar=('FAMILY', 'VCF'),
        help=(
            "a family's id and the sites of one of its callers: the caller's records on the primary chromosomes, "
            f'bgzipped, without samples, QUAL or FILTER, and with INFO cut down to {", ".join(sv_sites.KEPT)}; given '
            "once for each caller's VCF of each family"
        ),
    )
    parser.add_argument('--output-vcf', type=path, required=True, help="the union's VCF to write, uncompressed")
    parser.add_argument(
        '--output-map', type=path, required=True, help="the Parquet map from the families' records to the sites"
    )


def _sv_union(args: argparse.Namespace) -> None:
    families: dict[str, list[pathlib.Path]] = {}
    for family, vcf in args.sites:
        families.setdefault(family, []).append(pathlib.Path(vcf))
    union = sv_sites.union(families, args.output_vcf, args.output_map)
    if not union.sites:
        raise SystemExit("the families' callers give no record on the primary chromosomes, so nothing to annotate")
    print(f'{union.sites} sites of {union.alleles} alleles, {union.with_end} with END')


def _sv_annotation_tables_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument('--json', type=path, required=True, help="the Annotator's JSON of the union, gzipped")
    parser.add_argument(
        '--sites', type=path, required=True, help="the union of the group's sites, the VCF the JSON annotates, gzipped"
    )
    parser.add_argument('--run-config', type=path, required=True, help='the run config the Annotator ran with')
    parser.add_argument(
        '--bundle-config', type=path, required=True, help="the bundle's run config, runConfigs/bundle.json"
    )
    parser.add_argument('--catalog', required=True, help="the catalog's id")
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help='the folder to write sv/, sv_transcripts/ and annotation.json into; it must not exist',
    )


def _sv_annotation_tables(args: argparse.Namespace) -> None:
    info = sv_annotation.convert(
        json_path=args.json,
        sites_path=args.sites,
        run_config=json.loads(args.run_config.read_text()),
        bundle_config=json.loads(args.bundle_config.read_text()),
        catalog=args.catalog,
        output=args.output,
    )
    (args.output / 'annotation.json').write_text(info.model_dump_json(indent=2) + '\n')


def _echtvar_values_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument(
        '--annotated',
        type=path,
        required=True,
        help='the union as one `echtvar anno` call annotated it with every archive given',
    )
    parser.add_argument('--sites', type=path, required=True, help='the union, a sites-only VCF')
    parser.add_argument(
        '--archive',
        nargs=2,
        action='append',
        required=True,
        metavar=('SOURCE', 'ZIP'),
        help=(
            f'an echtvar archive that annotated it, after what it holds: one of {", ".join(echtvar_values.SOURCES)}; '
            f'given once for {echtvar_values.REQUIRED.name}, which every run reads, and once for each other archive '
            'of the run'
        ),
    )
    parser.add_argument('--echtvar-version', required=True, help='the version of the echtvar that annotated it')
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help=f'the folder to write {echtvar_values.TABLE_NAME} and {echtvar_values.INFO_NAME} into; it must not exist',
    )


def _echtvar_values(args: argparse.Namespace) -> None:
    archives: dict[str, pathlib.Path] = {}
    for source, archive in args.archive:
        if source in archives:
            raise SystemExit(f'--archive is given twice for {source}')
        archives[source] = pathlib.Path(archive)
    echtvar_values.convert(args.annotated, args.sites, archives, args.echtvar_version, args.output)


def _joined_values_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument('--sites', type=path, required=True, help='the union, a sites-only VCF')
    parser.add_argument(
        '--file',
        nargs=2,
        action='append',
        default=[],
        metavar=('SOURCE', 'FILE'),
        help=(
            f'a public file to join to the alleles, after what it holds: one of {", ".join(joined_values.SOURCES)}; '
            'given once for each file of the run, and not at all in a run without one'
        ),
    )
    parser.add_argument(
        '--normalised',
        nargs=2,
        action='append',
        default=[],
        metavar=('SOURCE', 'VCF'),
        help=(
            'with --file of a source that is a VCF of alleles: the file as bcftools norm rewrote it against the '
            "union's reference with the options that normalise a family, which is what is read; such sources are "
            f'{", ".join(name for name, source in joined_values.SOURCES.items() if source.normalised)}'
        ),
    )
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help=(
            f'the folder to write {joined_values.TABLE_NAME}, {joined_values.TRANSCRIPTS_NAME} and '
            f'{joined_values.INFO_NAME} into; it must not exist'
        ),
    )
    parser.add_argument('--scratch', type=path, help='a local folder for temporary files, rather than the default')


def _by_source(option: str, pairs: Sequence[Sequence[str]]) -> dict[str, pathlib.Path]:
    """The files of an option that is given once for a source, by the source's name."""
    files: dict[str, pathlib.Path] = {}
    for source, file in pairs:
        if source in files:
            raise SystemExit(f'{option} is given twice for {source}')
        files[source] = pathlib.Path(file)
    return files


def _joined_values(args: argparse.Namespace) -> None:
    joined_values.convert(
        args.sites,
        _by_source('--file', args.file),
        _by_source('--normalised', args.normalised),
        args.output,
        scratch=args.scratch,
    )


def _gene_models_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument(
        '--ensembl', type=path, required=True, help="Ensembl's GFF3 of the gene models, gzipped as published"
    )
    parser.add_argument(
        '--refseq', type=path, required=True, help="NCBI's RefSeq GFF of the gene models, gzipped as published"
    )
    parser.add_argument(
        '--assembly-report',
        type=path,
        required=True,
        help="NCBI's assembly report of the genome build, which names each chromosome as both files do",
    )
    parser.add_argument(
        '--reference',
        type=path,
        required=True,
        help='the reference the VCFs are called against, with its .fai beside it, whose sequences are chr1, chr2, ...',
    )
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help=(
            f'the folder to write the prepared files into, with {gene_models.EXONS_NAME}, '
            f'{gene_models.SYNONYMS_NAME}, {gene_models.CONTROLS_NAME} and {gene_models.RECORD_NAME}; it must not exist'
        ),
    )


def _gene_models(args: argparse.Namespace) -> None:
    gene_models.prepare(
        {'Ensembl': args.ensembl, 'RefSeq': args.refseq}, args.assembly_report, args.reference, args.output
    )


def _branch_points_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument(
        '--dataset-s1',
        type=path,
        required=True,
        help="Dataset S1 of the paper's supplement, the text file of every branch point on GRCh37, gzipped or not",
    )
    parser.add_argument(
        '--dataset-s2',
        type=path,
        required=True,
        help="Dataset S2 of the paper's supplement, the text file of the branch points by intron, gzipped or not",
    )
    parser.add_argument(
        '--chain',
        type=path,
        required=True,
        help=(
            'a chain file from GRCh37 to GRCh38, gzipped or not, whose chromosomes are named as UCSC names them '
            '(chr1) or as Ensembl does (1)'
        ),
    )
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help=(f'the folder to write {branch_points.FILE_NAME} and {branch_points.RECORD_NAME} into; it must not exist'),
    )


def _branch_points(args: argparse.Namespace) -> None:
    branch_points.prepare(args.dataset_s1, args.dataset_s2, args.chain, args.output)


def _transcript_rules_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument(
        '--given',
        type=path,
        required=True,
        help=(
            "the union's alleles that a rule can evaluate, each under its variant id, and the positive controls: a "
            'VCF compressed with bgzip, with its index beside it'
        ),
    )
    parser.add_argument(
        '--gene-models',
        type=path,
        required=True,
        help=f'the folder of the prepared gene models, with {gene_models.RECORD_NAME} (case-annotation gene-models)',
    )
    parser.add_argument(
        '--reference',
        type=path,
        required=True,
        help='the reference FASTA that the alleles are called against, with its .fai beside it',
    )
    parser.add_argument(
        '--uorf-evidence',
        type=path,
        required=True,
        help=(
            "UTRannotator's file of the upstream reading frames that were seen translated, "
            'uORF_5UTR_GRCh38_PUBLIC.txt of its repository'
        ),
    )
    parser.add_argument(
        '--branch-points',
        type=path,
        help=(
            f'the prepared catalogue of branch points, {branch_points.FILE_NAME} (case-annotation branch-points), '
            'with which the AG dinucleotides that an allele gains before an acceptor are found as well'
        ),
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=1,
        help='the number of processes that share the chromosomes (default: 1)',
    )
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help=(
            f'the folder to write {transcript_rules.TABLE_NAME} and {transcript_rules.INFO_NAME} into; '
            'it must not exist'
        ),
    )


def _transcript_rules(args: argparse.Namespace) -> None:
    if args.workers < 1:
        raise ValueError(f'--workers is {args.workers}, and at least one process is needed')
    transcript_rules.annotate(
        args.given,
        args.gene_models,
        args.reference,
        args.uorf_evidence,
        args.output,
        workers=args.workers,
        branch_points=args.branch_points,
    )


def _family_tables_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument('--family', required=True, help="the family's id, which its PED file must hold")
    parser.add_argument('--ped', type=path, required=True, help="the family's PED file")
    parser.add_argument('--vcf', type=path, required=True, help="the family's normalised VCF, with its index beside it")
    _proband_arguments(parser)
    _structural_arguments(parser)
    parser.add_argument(
        '--sv-alleles',
        type=path,
        help=(
            'with --sv-vcf: its records that are written as a sequence, on the primary chromosomes, normalised as the '
            "family's VCF is and each under its ID, which link a structural variant to the same allele in variants"
        ),
    )
    parser.add_argument(
        '--sv-annotation',
        type=path,
        help=(
            "with --sv-vcf or --cnv-vcf: sv-annotation-tables' output folder of the group's sites, which annotates "
            "the family's records of both callers; given with --sv-map"
        ),
    )
    parser.add_argument(
        '--sv-map', type=path, help="with --sv-annotation: sv-union's map from the families' records to the sites"
    )
    parser.add_argument('--annotation', type=path, required=True, help="annotation-tables' output folder of the union")
    parser.add_argument('--echtvar-values', type=path, required=True, help="echtvar-values' output folder of the union")
    parser.add_argument('--joined-values', type=path, required=True, help="joined-values' output folder of the union")
    parser.add_argument(
        '--transcript-rules', type=path, required=True, help="transcript-rules' output folder of the union"
    )
    parser.add_argument(
        '--roh-upd',
        type=path,
        required=True,
        help="roh-upd's output folder of the family, which holds its runs of homozygosity and regions of disomy",
    )
    parser.add_argument(
        '--strchive',
        type=path,
        help=(
            "the loci of one release of STRchive, a JSON file, whose size ranges go into the family's repeats table; "
            'given with --strchive-release and --strchive-release-date, or not at all'
        ),
    )
    parser.add_argument(
        '--strchive-release',
        help='with --strchive: the release that the file is of, such as v2.28.0, which the file does not state',
    )
    parser.add_argument(
        '--strchive-release-date',
        help='with --strchive: the day on which that release was published, as YYYY-MM-DD',
    )
    parser.add_argument(
        '--schema-version',
        required=True,
        help="the tables' layout version the caller expects, which names its output folder; another fails",
    )
    _release_argument(parser)
    parser.add_argument(
        '--output', type=path, required=True, help='the folder to write the tables into; it must not exist'
    )
    parser.add_argument('--scratch', type=path, help='a local folder for temporary files, rather than the default')


def _family_tables(args: argparse.Namespace) -> None:
    if args.schema_version != family_tables.SCHEMA_VERSION:
        raise SystemExit(f'the tables are of version {family_tables.SCHEMA_VERSION}, not {args.schema_version}')
    _require_release(args.pipeline_version)
    if (args.sv_annotation is None) != (args.sv_map is None):
        raise SystemExit('--sv-annotation and --sv-map are given together')
    sv_annotated = None
    if args.sv_annotation is not None:
        sv_annotated = family_tables.SvAnnotation(tables=args.sv_annotation, site_map=args.sv_map)
    release = (args.strchive, args.strchive_release, args.strchive_release_date)
    if any(part is None for part in release) and any(part is not None for part in release):
        raise SystemExit('--strchive, --strchive-release and --strchive-release-date are given together')
    ranges = None if args.strchive is None else strchive.stated(*release)
    family_tables.assemble(
        args.family,
        args.vcf,
        args.ped,
        args.annotation,
        args.echtvar_values,
        args.joined_values,
        args.transcript_rules,
        args.output,
        proband=args.proband,
        structural=_structural(args),
        segments=args.roh_upd,
        sv_alleles=args.sv_alleles,
        sv_annotated=sv_annotated,
        ranges=ranges,
        scratch=args.scratch,
    )


def _roh_upd_family_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument('--family', required=True, help="the family's id, which its PED file must hold")
    parser.add_argument('--ped', type=path, required=True, help="the family's PED file")
    parser.add_argument(
        '--archive',
        type=path,
        required=True,
        help=f"the echtvar archive of {echtvar_values.REQUIRED.name}, which gives each site gnomAD's frequency",
    )


def _roh_upd_plan_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    _roh_upd_family_arguments(parser)
    parser.add_argument('--vcf', type=path, required=True, help="the family's VCF, as the run's sheet names it")
    parser.add_argument(
        '--trios',
        type=path,
        required=True,
        help=(
            'the file to write the children that upd assesses into: a line for each child whose parents are both '
            'samples of the VCF, with the child, its mother and its father, separated by tabs'
        ),
    )


def _roh_upd_plan(args: argparse.Namespace) -> None:
    planned = roh_upd.plan(args.family, args.ped, args.vcf, args.archive)
    args.trios.write_text(''.join(f'{trio.child}\t{trio.mother}\t{trio.father}\n' for trio in planned.trios))
    print(planned.frequency.alias)


def _roh_upd_sites_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--prepared', type=pathlib.Path, required=True, help='the VCF that the step prepared for both tools, bgzipped'
    )
    parser.add_argument(
        '--af-tag', required=True, help='the INFO field of the frequency, as upd is given it and roh-upd-plan prints it'
    )
    parser.add_argument('--min-af', type=float, required=True, help="upd's lowest frequency of a site, as upd is given")
    parser.add_argument('--min-gq', type=int, required=True, help="upd's lowest GQ of a call, as upd is given it")


def _roh_upd_sites(args: argparse.Namespace) -> None:
    kept = roh_upd.keeps_a_site(args.prepared, args.af_tag, args.min_af, args.min_gq)
    print('kept' if kept else 'none')


def _roh_upd_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    _roh_upd_family_arguments(parser)
    parser.add_argument(
        '--prepared',
        type=path,
        required=True,
        help='the VCF that both tools read: the one the step prepared and echtvar annotated with the archive',
    )
    parser.add_argument('--sites', type=int, required=True, help='the number of records of the prepared VCF')
    parser.add_argument('--roh', type=path, required=True, help='the output of `bcftools roh -O r` on the prepared VCF')
    parser.add_argument(
        '--regions',
        nargs=2,
        action='append',
        default=[],
        metavar=('CHILD', 'BED'),
        help=(
            "a child's sample and the output of `upd regions` for it and its parents; given once for each child "
            'that roh-upd-plan lists and upd ran for'
        ),
    )
    parser.add_argument(
        '--not-assessed',
        action='append',
        default=[],
        metavar='CHILD',
        help=(
            'a child that roh-upd-plan lists and upd was not run for, because roh-upd-sites found no site that upd '
            'keeps; given once for each such child'
        ),
    )
    parser.add_argument('--bcftools-version', required=True, help='the version of the bcftools that ran')
    parser.add_argument('--echtvar-version', required=True, help='the version of the echtvar that annotated the VCF')
    parser.add_argument(
        '--preparation',
        action='append',
        required=True,
        help='a command that prepared the VCF, without its files; given once for each, in the order they ran',
    )
    parser.add_argument('--roh-command', required=True, help='the bcftools roh command that ran, without its files')
    parser.add_argument(
        '--upd-command', required=True, help='the upd command that ran, without its files and its three samples'
    )
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help=(
            f'the folder to write {roh_upd.ROH_NAME}, {roh_upd.UPD_NAME} for a family with a child of two samples, and '
            f'{roh_upd.INFO_NAME} into; it must not exist'
        ),
    )


def _roh_upd(args: argparse.Namespace) -> None:
    roh_upd.convert(
        args.family,
        args.ped,
        args.prepared,
        args.archive,
        args.roh,
        _by_source('--regions', args.regions),
        args.output,
        not_assessed=args.not_assessed,
        sites=args.sites,
        bcftools_version=args.bcftools_version,
        echtvar_version=args.echtvar_version,
        upd=roh_upd.installed_upd(),
        preparation=args.preparation,
        roh_command=args.roh_command,
        upd_command=args.upd_command,
    )


def _alphagenome_table_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--table',
        required=True,
        choices=tuple(alphagenome_archives.TABLES),
        help="which of DeepMind's tables, by the name of the archive that is built from it",
    )


def _alphagenome_vcf_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    _alphagenome_table_argument(parser)
    parser.add_argument(
        '--reference-index',
        type=path,
        required=True,
        help="the .fai of the reference, whose primary chromosomes the VCF's header declares with their lengths",
    )
    parser.add_argument(
        '--tsv',
        type=path,
        help=(
            'the table as DeepMind publishes it, bgzipped and with its .tbi beside it; without the option the '
            "table's text is read from standard input, header first"
        ),
    )
    parser.add_argument(
        '--chromosome', help='with --tsv: the one chromosome to write the rows of, such as chr1, rather than all'
    )


def _alphagenome_vcf(args: argparse.Namespace) -> None:
    if args.chromosome is not None and args.tsv is None:
        raise SystemExit('--chromosome picks the rows of an indexed table, which --tsv gives')
    table = alphagenome_archives.TABLES[args.table]
    lines = sys.stdin if args.tsv is None else alphagenome_archives.indexed_lines(args.tsv, args.chromosome)
    alphagenome_archives.to_vcf(table, lines, alphagenome_archives.read_contigs(args.reference_index), sys.stdout)


def _alphagenome_config(args: argparse.Namespace) -> None:
    sys.stdout.write(alphagenome_archives.config_text(alphagenome_archives.TABLES[args.table]))


def _alphagenome_check_arguments(parser: argparse.ArgumentParser) -> None:
    _alphagenome_table_argument(parser)
    parser.add_argument(
        '--annotated',
        type=pathlib.Path,
        required=True,
        help=(
            'the alleles of the rows on standard input as alphagenome-vcf writes them, in the same order, after '
            "`echtvar anno` annotated them with the table's archive"
        ),
    )


def _alphagenome_check(args: argparse.Namespace) -> None:
    compared = alphagenome_archives.check(alphagenome_archives.TABLES[args.table], sys.stdin, args.annotated)
    print(f'the archive returns the {compared} rows of {args.table} that were compared')


def _contract_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--out',
        type=pathlib.Path,
        required=True,
        help="the folder to write the contract's files into, which replace those of the same names",
    )


def _contract(args: argparse.Namespace) -> None:
    contract.write(args.out)


# Each subcommand: its help text, a function that adds its arguments to its parser, and the function that runs it.
_SUBCOMMANDS: dict[str, tuple[str, Callable[[argparse.ArgumentParser], None], Callable[[argparse.Namespace], None]]] = {
    'family-check': (
        (
            "checks a family's PED file against its id and its VCF's samples, that each sample id can name a column "
            "of the variants table (it is no allele column's name and no other sample's id, whatever the case of its "
            'letters, and no column of the sv or the repeats table either), that the proband the run names for the '
            'family, if any, is an individual of the PED file and a sample of the VCF, and that the VCFs of its '
            "structural variants, copy numbers and repeats hold the family's samples, before any annotation runs"
        ),
        _family_check_arguments,
        _family_check,
    ),
    'annotation-tables': (
        (
            "converts the Annotator's JSON of a sites union into its variants and transcripts tables, each gene with "
            'its HGNC id from the genes section'
        ),
        _annotation_tables_arguments,
        _annotation_tables,
    ),
    'sv-union': (
        (
            "writes the union of a group's sites of structural variants and copy-number segments, one record per "
            "allele under an id of its content, and the map from each family's records to them"
        ),
        _sv_union_arguments,
        _sv_union,
    ),
    'sv-annotation-tables': (
        "converts the Annotator's JSON of a union of such sites into its sv and sv_transcripts tables",
        _sv_annotation_tables_arguments,
        _sv_annotation_tables,
    ),
    'echtvar-values': (
        (
            "writes the values of a sites union from its echtvar annotation, which are gnomAD's joint values and, "
            "where the run has their archives, AlphaGenome's scores, with each archive's record"
        ),
        _echtvar_values_arguments,
        _echtvar_values,
    ),
    'joined-values': (
        (
            'writes the values of a sites union from the public files that a run joins to its alleles by allele or by '
            "position, as a table by allele and one by allele and transcript, with each file's record: gnomAD's "
            "mitochondrial frequencies, MitoTIP's scores, GENCODE's polyadenylation signals, APARENT2's effects on "
            "polyadenylation, gnomAD's genomic constraint by window, regions of missense constraint and MPC scores, "
            'and the prepared catalogue of branch points; fails for a file that is not what its source publishes'
        ),
        _joined_values_arguments,
        _joined_values,
    ),
    'roh-upd-plan': (
        (
            'says what bcftools roh and upd are given for a family: prints the INFO field that holds the frequency '
            'of an allele once echtvar has annotated a VCF with the archive, and lists the children whose parents '
            'are both sequenced; fails for an archive that writes another stand-in than 0 for the frequency of an '
            'allele it lacks, and for a family with such a child whose VCF has no GQ'
        ),
        _roh_upd_plan_arguments,
        _roh_upd_plan,
    ),
    'roh-upd-sites': (
        (
            'says whether upd would keep a site of the prepared VCF, by the rule of its code: prints `kept` where a '
            'record that is not the last has a frequency of at least --min-af and a GQ of at least --min-gq in every '
            'sample, and `none` otherwise, where upd would fail'
        ),
        _roh_upd_sites_arguments,
        _roh_upd_sites,
    ),
    'roh-upd': (
        (
            "writes a family's runs of homozygosity and regions of uniparental disomy as tables, from the output of "
            'bcftools roh and of upd regions on the VCF that the step prepared, with a record of what ran; fails on '
            'output that is not what the tools write'
        ),
        _roh_upd_arguments,
        _roh_upd,
    ),
    'gene-models': (
        (
            "prepares the gene models that the transcript rules read from Ensembl's and NCBI's public files: each file "
            "sorted, bgzipped and indexed, the chromosome synonyms, the exons' regions and the positive controls"
        ),
        _gene_models_arguments,
        _gene_models,
    ),
    'branch-points': (
        (
            'prepares the catalogue of branch points that a run joins to its alleles, from Datasets S1 and S2 of '
            'Zhang et al. (PNAS 2022) and a chain file: every branch point that lifts to GRCh38 with the two bases '
            'before it, sorted and bgzipped, with a record of the three files and of the counts; fails for a dataset '
            'that is not what the paper published'
        ),
        _branch_points_arguments,
        _branch_points,
    ),
    'transcript-rules': (
        (
            "writes what the rules of VEP's UTRannotator and NMD plugins say of a sites union as a table by allele "
            'and transcript, from the prepared gene models and the reference and without running VEP, after checking '
            'the positive controls, with the record of what the values come from'
        ),
        _transcript_rules_arguments,
        _transcript_rules,
    ),
    'family-tables': (
        (
            "writes a family's variants table, with a column of each member's calls, and its transcripts table from "
            "its VCF and the union's tables, with what the transcript rules say of each transcript, its sv, "
            'sv_transcripts and repeats tables from the VCFs of its structural variants, copy numbers and repeats '
            'where it has '
            "them, with the annotation of the group's "
            "structural variants and, where the run gives a release of STRchive, STRchive's size ranges of the "
            'repeats, its roh table and, for a family with a child whose parents are both sequenced, its upd table '
            "from roh-upd's output, its manifest with the roles of the proband the run names and of its parents, and "
            'a copy of the notes on the tables (tables.md)'
        ),
        _family_tables_arguments,
        _family_tables,
    ),
    'alphagenome-vcf': (
        (
            "writes one of AlphaGenome's tables of scores as the VCF that `echtvar encode` reads with the table's "
            "configuration, on standard output; fails on a header that is not the table's, on a malformed row and "
            'on a row with another value than 0 in a column that the archive leaves out'
        ),
        _alphagenome_vcf_arguments,
        _alphagenome_vcf,
    ),
    'alphagenome-config': (
        "writes the configuration that `echtvar encode` takes for one of AlphaGenome's tables, on standard output",
        _alphagenome_table_argument,
        _alphagenome_config,
    ),
    'alphagenome-check': (
        (
            "compares rows of one of AlphaGenome's tables, read from standard input with the table's header, with "
            "what `echtvar anno` returned for their alleles from the table's archive; fails on a value that differs"
        ),
        _alphagenome_check_arguments,
        _alphagenome_check,
    ),
    'contract': (
        (
            "writes the pipeline's published contract (contract/): the JSON Schema of a family's manifest, the "
            "families sheet's columns, the workflows' params and published folders, the tables' columns and their "
            'version, and an invented sheet and manifest'
        ),
        _contract_arguments,
        _contract,
    ),
}


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog='case-annotation', description=__doc__)
    subcommands = root.add_subparsers(dest='subcommand', required=True)
    for name, (help_text, add_arguments, run) in _SUBCOMMANDS.items():
        subcommand = subcommands.add_parser(name, help=help_text, description=help_text)
        add_arguments(subcommand)
        subcommand.set_defaults(run=run)
    return root


def main(argv: Sequence[str] | None = None) -> None:
    args = parser().parse_args(argv)
    args.run(args)
