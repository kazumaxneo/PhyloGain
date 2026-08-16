from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .project import InputError, build_project, validate_inputs
from .server import serve_project


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        prog="phylogain",
        description="Build interactive branch-level gene-family gain/loss maps from OrthoFinder or PIRATE results.",
    )
    root.add_argument("--version", action="version", version=__version__)
    commands = root.add_subparsers(dest="command", required=True)

    validate = commands.add_parser("validate", help="Validate comparative-genomics and phenotype inputs")
    _input_arguments(validate)
    validate.add_argument("--json", action="store_true", help="Print machine-readable JSON")

    build = commands.add_parser("build", help="Infer events and build an interactive map")
    _input_arguments(build)
    build.add_argument("--output", required=True, help="New or empty output directory")
    build.add_argument("--gain-cost", type=float, default=1.0)
    build.add_argument("--loss-cost", type=float, default=1.0)
    build.add_argument("--root-state", choices=["auto", "0", "1"], default="auto")
    build.add_argument("--presence-threshold", type=int, default=1)
    build.add_argument(
        "--no-members",
        action="store_true",
        help="Skip indexing gene-family member IDs for a smaller output",
    )
    build.add_argument(
        "--annotate",
        choices=["eggnog"],
        help="Run optional functional annotation (currently: eggnog)",
    )
    build.add_argument(
        "--annotations",
        help="Import an existing eggNOG-mapper .emapper.annotations file",
    )
    build.add_argument(
        "--proteomes",
        help="Protein FASTA directory used by OrthoFinder (required with --annotate eggnog)",
    )
    build.add_argument("--eggnog-emapper", default="emapper.py")
    build.add_argument("--eggnog-data-dir")
    build.add_argument(
        "--go-obo",
        help="Official go-basic.obo file used to add readable names to GO identifiers",
    )
    build.add_argument(
        "--fetch-kegg-names",
        action="store_true",
        help="Fetch and cache official KEGG names at build time (academic use only)",
    )
    build.add_argument("--annotation-cpu", type=int, default=1)

    serve = commands.add_parser("serve", help="Open a generated map in a local web server")
    serve.add_argument("project", help="Generated Gene Gain/Loss Viewer directory")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8000)
    serve.add_argument("--no-open", action="store_true", help="Do not open a browser")
    return root


def _input_arguments(command: argparse.ArgumentParser) -> None:
    source = command.add_mutually_exclusive_group(required=True)
    source.add_argument("--orthofinder", help="OrthoFinder results directory (recommended)")
    source.add_argument("--pirate", help="PIRATE output directory")
    command.add_argument(
        "--species-tree",
        help="Optional rooted Newick tree replacing the input tool's default tree",
    )
    command.add_argument(
        "--gtdb-taxonomy",
        help="Optional GTDB or GTDB-Tk taxonomy TSV for rank-based clade collapsing",
    )
    command.add_argument(
        "--genome-metadata",
        help=(
            "Optional TSV containing a species ID column and genome_size_bp "
            "or genome_size_mb for the tree-side bar graph"
        ),
    )
    command.add_argument(
        "--tip-metadata",
        help=(
            "Optional TSV containing a species ID column and strain_name, "
            "organism_name, display_name, or another supported tip-label column"
        ),
    )
    command.add_argument("--phenotypes", help="Optional wide TSV with species_id in the first column")
    command.add_argument(
        "--phenotype",
        action="append",
        help="Phenotype column to analyze; repeat to select several (default: all)",
    )


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    try:
        if args.command == "validate":
            report = validate_inputs(
                args.orthofinder,
                args.phenotypes,
                args.phenotype,
                args.species_tree,
                args.gtdb_taxonomy,
                args.pirate,
                args.genome_metadata,
                args.tip_metadata,
            )
            if args.json:
                print(json.dumps(report, ensure_ascii=False, indent=2))
            else:
                _print_validation(report)
            if not report["ok"]:
                raise SystemExit(2)
        elif args.command == "build":
            result = build_project(
                orthofinder=args.orthofinder,
                output=args.output,
                phenotype_path=args.phenotypes,
                selected_phenotypes=args.phenotype,
                gain_cost=args.gain_cost,
                loss_cost=args.loss_cost,
                root_state=args.root_state,
                presence_threshold=args.presence_threshold,
                include_members=not args.no_members,
                species_tree_path=args.species_tree,
                gtdb_taxonomy_path=args.gtdb_taxonomy,
                annotate=args.annotate,
                annotation_path=args.annotations,
                go_obo_path=args.go_obo,
                fetch_kegg_names=args.fetch_kegg_names,
                proteomes=args.proteomes,
                eggnog_emapper=args.eggnog_emapper,
                eggnog_data_dir=args.eggnog_data_dir,
                annotation_cpu=args.annotation_cpu,
                progress=lambda message: print(message, flush=True),
                pirate=args.pirate,
                genome_metadata_path=args.genome_metadata,
                tip_metadata_path=args.tip_metadata,
            )
            print(
                f"Built {result['output']} ({result['species']} species, "
                f"{result['families']} families, {result['branches']} branches)."
            )
            print(f"Open it with: phylogain serve \"{result['output']}\"")
            for warning in result["warnings"]:
                print(f"WARNING: {warning}", file=sys.stderr)
        elif args.command == "serve":
            serve_project(args.project, args.host, args.port, not args.no_open)
    except (InputError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


def _print_validation(report: dict[str, object]) -> None:
    print(f"Input format:       {str(report.get('input_format', '')).upper()}")
    if "tree_tips" in report:
        print(f"Species tree tips:  {report['tree_tips']}")
        print(f"Gene-count species: {report['count_species']}")
        print(f"Matched species IDs: {report['matched_species']}")
    phenotypes = report.get("phenotypes") or []
    if phenotypes:
        print(f"Phenotypes:         {', '.join(phenotypes)}")
    taxonomy = report.get("taxonomy") or {}
    if taxonomy:
        print(f"GTDB taxonomy:      {taxonomy['mapped_species']} species mapped")
    genome_metadata = report.get("genome_metadata") or {}
    if genome_metadata:
        print(
            "Genome sizes:       "
            f"{genome_metadata['mapped_species']} species mapped"
        )
    tip_metadata = report.get("tip_metadata") or {}
    if tip_metadata:
        print(
            "Tip labels:         "
            f"{tip_metadata['mapped_species']} species mapped"
        )
    for warning in report["warnings"]:
        print(f"WARNING: {warning}")
    for error in report["errors"]:
        print(f"ERROR: {error}")
    print("Validation passed." if report["ok"] else "Validation failed.")
