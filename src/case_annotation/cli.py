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
    contract,
    family_tables,
    gnomad_joint,
    sv_annotation,
    sv_sites,
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


def _gnomad_joint_arguments(parser: argparse.ArgumentParser) -> None:
    path = pathlib.Path
    parser.add_argument(
        '--annotated', type=path, required=True, help='the union as `echtvar anno` annotated it with the archive'
    )
    parser.add_argument('--sites', type=path, required=True, help='the union, a sites-only VCF')
    parser.add_argument(
        '--archive', type=path, required=True, help="the echtvar archive of gnomAD's joint release that annotated it"
    )
    parser.add_argument('--echtvar-version', required=True, help='the version of the echtvar that annotated it')
    parser.add_argument(
        '--output',
        type=path,
        required=True,
        help=f'the folder to write {gnomad_joint.TABLE_NAME} and {gnomad_joint.INFO_NAME} into; it must not exist',
    )


def _gnomad_joint(args: argparse.Namespace) -> None:
    gnomad_joint.convert(args.annotated, args.sites, args.archive, args.echtvar_version, args.output)


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
    parser.add_argument('--gnomad-joint', type=path, required=True, help="gnomad-joint's output folder of the union")
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
    family_tables.assemble(
        args.family,
        args.vcf,
        args.ped,
        args.annotation,
        args.gnomad_joint,
        args.output,
        proband=args.proband,
        structural=_structural(args),
        sv_alleles=args.sv_alleles,
        sv_annotated=sv_annotated,
        scratch=args.scratch,
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
    'gnomad-joint': (
        "writes gnomAD's joint values of a sites union from its echtvar annotation, with the archive's record",
        _gnomad_joint_arguments,
        _gnomad_joint,
    ),
    'family-tables': (
        (
            "writes a family's variants table, with a column of each member's calls, and its transcripts table from "
            "its VCF and the union's tables, its sv, sv_transcripts and repeats tables from the VCFs of its "
            "structural variants, copy numbers and repeats where it has them, with the annotation of the group's "
            'structural variants, its manifest with the roles of the proband the run names and of its parents, and '
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
