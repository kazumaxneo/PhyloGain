from __future__ import annotations

import csv
import json
import platform
import shutil
import sqlite3
import sys
import zlib
from collections import Counter
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path

from . import __version__
from .analysis import (
    PIRATE_METADATA_COLUMNS,
    Reconstruction,
    find_orthofinder_files,
    find_pirate_files,
    iter_gene_counts,
    iter_pirate_counts,
    read_count_header,
    read_pirate_header,
    read_phenotypes,
    split_pirate_loci,
)
from .annotation import (
    fetch_official_kegg_names,
    import_go_term_names,
    import_eggnog_annotations,
    run_eggnog_mapper,
    write_representative_fasta,
)
from .tree import as_project_nodes, leaf_labels, parse_newick, preorder
from .taxonomy import read_gtdb_taxonomy
from .metadata import read_genome_sizes


class InputError(ValueError):
    pass


def validate_inputs(
    orthofinder: str | Path | None,
    phenotype_path: str | Path | None = None,
    selected_phenotypes: list[str] | None = None,
    species_tree_path: str | Path | None = None,
    gtdb_taxonomy_path: str | Path | None = None,
    pirate: str | Path | None = None,
    genome_metadata_path: str | Path | None = None,
) -> dict[str, object]:
    input_format, found = _resolve_input_files(orthofinder, pirate)
    errors: list[str] = []
    warnings: list[str] = []
    tree_path = Path(species_tree_path) if species_tree_path else found["tree"]
    if tree_path is None or not tree_path.is_file():
        location = (
            "Species_Tree/" if input_format == "orthofinder" else "the PIRATE directory"
        )
        errors.append(f"No tree found in {location}; provide --species-tree")
    if found["counts"] is None:
        expected = (
            "Orthogroups/Orthogroups.GeneCount.tsv"
            if input_format == "orthofinder"
            else "PIRATE.gene_families.tsv"
        )
        errors.append(f"{expected} was not found")
    if errors:
        return {
            "ok": False,
            "errors": errors,
            "warnings": warnings,
            "input_format": input_format,
        }

    root = parse_newick(tree_path)
    tips = leaf_labels(root)
    if not all(tips):
        errors.append("Every species-tree tip must have a label")
    duplicate_tips = sorted(name for name, count in Counter(tips).items() if count > 1)
    if duplicate_tips:
        errors.append(f"Duplicate species-tree tips: {', '.join(duplicate_tips[:10])}")

    count_species, _ = (
        read_count_header(found["counts"])
        if input_format == "orthofinder"
        else read_pirate_header(found["counts"])
    )
    tree_only = sorted(set(tips) - set(count_species))
    counts_only = sorted(set(count_species) - set(tips))
    if tree_only or counts_only:
        errors.append(
            f"Tree/count species mismatch: {len(tree_only)} tree-only and {len(counts_only)} count-only"
        )

    phenotype_names: list[str] = []
    phenotype_species: set[str] = set()
    if phenotype_path:
        phenotype_names, phenotype_data = read_phenotypes(
            phenotype_path, selected_phenotypes
        )
        for values in phenotype_data.values():
            phenotype_species.update(values)
        phenotype_only = sorted(phenotype_species - set(tips))
        if phenotype_only:
            errors.append(
                f"Phenotype TSV contains {len(phenotype_only)} species absent from the tree"
            )
        absent = sorted(set(tips) - phenotype_species)
        if absent:
            warnings.append(
                f"Phenotype TSV omits {len(absent)} tree species; they will be treated as unknown"
            )

    if found["members"] is None:
        member_name = (
            "Orthogroups/Orthogroups.tsv"
            if input_format == "orthofinder"
            else "PIRATE.gene_families.tsv"
        )
        warnings.append(f"{member_name} was not found; gene IDs cannot be shown in the HTML")
    if input_format == "pirate" and not species_tree_path:
        warnings.append(
            "Using PIRATE binary_presence_absence.nwk, a gene-content tree; "
            "a rooted external species tree is recommended for gain/loss inference"
        )

    taxonomy_report: dict[str, object] | None = None
    if gtdb_taxonomy_path:
        taxonomy_file = Path(gtdb_taxonomy_path)
        if not taxonomy_file.is_file():
            errors.append(f"GTDB taxonomy TSV was not found: {taxonomy_file}")
        else:
            _, taxonomy_report = read_gtdb_taxonomy(taxonomy_file, tips)
            if taxonomy_report["mapped_species"] == 0:
                errors.append("GTDB taxonomy TSV did not match any species-tree tips")
            elif taxonomy_report["unmapped_species"]:
                warnings.append(
                    "GTDB taxonomy TSV omits or cannot match "
                    f"{taxonomy_report['unmapped_species']} tree species"
                )
            if taxonomy_report["ambiguous_rows"]:
                warnings.append(
                    f"GTDB taxonomy TSV has {taxonomy_report['ambiguous_rows']} ambiguous ID rows"
                )

    genome_metadata_report: dict[str, object] | None = None
    if genome_metadata_path:
        genome_metadata_file = Path(genome_metadata_path)
        if not genome_metadata_file.is_file():
            errors.append(f"Genome metadata TSV was not found: {genome_metadata_file}")
        else:
            _, genome_metadata_report = read_genome_sizes(genome_metadata_file, tips)
            if genome_metadata_report["mapped_species"] == 0:
                errors.append("Genome metadata TSV did not match any species-tree tips")
            elif genome_metadata_report["unmapped_species"]:
                warnings.append(
                    "Genome metadata TSV omits or cannot match "
                    f"{genome_metadata_report['unmapped_species']} tree species"
                )
            for key, label in (
                ("ambiguous_rows", "ambiguous ID"),
                ("invalid_rows", "invalid genome-size"),
                ("duplicate_rows", "duplicate species"),
            ):
                if genome_metadata_report[key]:
                    warnings.append(
                        f"Genome metadata TSV has {genome_metadata_report[key]} {label} rows"
                    )

    return {
        "ok": not errors,
        "errors": errors,
        "warnings": warnings,
        "tree_tips": len(tips),
        "count_species": len(count_species),
        "matched_species": len(set(tips) & set(count_species)),
        "phenotypes": phenotype_names,
        "tree": str(tree_path),
        "counts": str(found["counts"]),
        "members": str(found["members"]) if found["members"] else None,
        "taxonomy": taxonomy_report,
        "genome_metadata": genome_metadata_report,
        "input_format": input_format,
    }


def _resolve_input_files(
    orthofinder: str | Path | None,
    pirate: str | Path | None,
) -> tuple[str, dict[str, Path | None]]:
    if bool(orthofinder) == bool(pirate):
        raise InputError("Specify exactly one of --orthofinder or --pirate")
    if orthofinder:
        return "orthofinder", find_orthofinder_files(orthofinder)
    return "pirate", find_pirate_files(pirate)


def build_project(
    orthofinder: str | Path | None,
    output: str | Path,
    phenotype_path: str | Path | None = None,
    selected_phenotypes: list[str] | None = None,
    gain_cost: float = 1.0,
    loss_cost: float = 1.0,
    root_state: str | int = "auto",
    presence_threshold: int = 1,
    include_members: bool = True,
    species_tree_path: str | Path | None = None,
    gtdb_taxonomy_path: str | Path | None = None,
    annotate: str | None = None,
    annotation_path: str | Path | None = None,
    go_obo_path: str | Path | None = None,
    fetch_kegg_names: bool = False,
    proteomes: str | Path | None = None,
    eggnog_emapper: str = "emapper.py",
    eggnog_data_dir: str | Path | None = None,
    annotation_cpu: int = 1,
    progress=None,
    pirate: str | Path | None = None,
    genome_metadata_path: str | Path | None = None,
) -> dict[str, object]:
    progress = progress or (lambda message: None)
    report = validate_inputs(
        orthofinder,
        phenotype_path,
        selected_phenotypes,
        species_tree_path,
        gtdb_taxonomy_path,
        pirate,
        genome_metadata_path,
    )
    if not report["ok"]:
        raise InputError("; ".join(report["errors"]))
    input_format, found = _resolve_input_files(orthofinder, pirate)
    if annotate and annotation_path:
        raise InputError("Use either --annotate or --annotations, not both")
    if annotate == "eggnog" and input_format == "orthofinder" and not proteomes:
        raise InputError("--proteomes is required with --annotate eggnog")
    if (annotate or annotation_path) and not found["members"]:
        raise InputError("A gene-family membership table is required for functional annotation")
    if annotate == "eggnog" and input_format == "pirate" and not found.get("representatives"):
        raise InputError("representative_sequences.faa is required with PIRATE --annotate eggnog")
    if annotation_path and not Path(annotation_path).is_file():
        raise InputError(f"Annotation file was not found: {annotation_path}")
    if go_obo_path and not Path(go_obo_path).is_file():
        raise InputError(f"GO ontology file was not found: {go_obo_path}")
    if go_obo_path and not (annotate or annotation_path):
        raise InputError("--go-obo requires --annotate or --annotations")
    if fetch_kegg_names and not (annotate or annotation_path):
        raise InputError("--fetch-kegg-names requires --annotate or --annotations")
    output_path = Path(output).resolve()
    if output_path.exists() and any(output_path.iterdir()):
        raise InputError(f"Output directory is not empty: {output_path}")
    output_path.mkdir(parents=True, exist_ok=True)

    tree_path = Path(report["tree"])
    root = parse_newick(tree_path)
    nodes = preorder(root)
    species = leaf_labels(root)
    count_species, _ = (
        read_count_header(found["counts"])
        if input_format == "orthofinder"
        else read_pirate_header(found["counts"])
    )
    count_iterator = (
        iter_gene_counts if input_format == "orthofinder" else iter_pirate_counts
    )
    reconstruction = Reconstruction(root, gain_cost=gain_cost, loss_cost=loss_cost)
    taxonomy_data: dict[str, dict[str, str]] = {}
    taxonomy_report: dict[str, object] | None = None
    if gtdb_taxonomy_path:
        taxonomy_data, taxonomy_report = read_gtdb_taxonomy(gtdb_taxonomy_path, species)
        _write_taxonomy(output_path / "gtdb_taxonomy.tsv", species, taxonomy_data)
    genome_sizes: dict[str, int] = {}
    if genome_metadata_path:
        genome_sizes, _ = read_genome_sizes(genome_metadata_path, species)
        _write_genome_metadata(
            output_path / "genome_metadata.tsv", species, genome_sizes
        )

    database_path = output_path / "species_map.sqlite"
    connection = sqlite3.connect(database_path)
    connection.execute("PRAGMA journal_mode=MEMORY")
    connection.execute("PRAGMA synchronous=OFF")
    _create_schema(connection)

    branch_rows = [
        (
            node.branch_id,
            node.parent.node_id if node.parent else None,
            node.node_id,
            node.label if node.is_leaf else "",
            node.depth,
            node.length,
        )
        for node in nodes
        if node.parent is not None
    ]
    connection.executemany(
        "INSERT INTO branches(branch_id,parent_node,child_node,tip_label,depth,branch_length) VALUES(?,?,?,?,?,?)",
        branch_rows,
    )

    branch_counts = {row[0]: [0, 0] for row in branch_rows}
    event_tsv = output_path / "gene_gain_loss.tsv"
    progress("Inferring gene-family gains and losses")
    family_count = 0
    with event_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["branch_id", "family_id", "event", "parent_state", "child_state"]
        )
        event_batch: list[tuple[str, str, str, int, int]] = []
        family_batch: list[tuple[str]] = []
        for family_id, counts in count_iterator(found["counts"]):
            family_count += 1
            states = {
                species_id: int(count >= presence_threshold)
                for species_id, count in zip(count_species, counts)
            }
            transitions, _, _ = reconstruction.infer(states, root_state=root_state)
            family_batch.append((family_id,))
            for transition in transitions:
                row = (
                    transition.branch_id,
                    family_id,
                    transition.event,
                    transition.parent_state,
                    transition.child_state,
                )
                event_batch.append(row)
                writer.writerow(row)
                branch_counts[transition.branch_id][0 if transition.event == "gain" else 1] += 1
            if len(family_batch) >= 1000:
                connection.executemany(
                    "INSERT OR IGNORE INTO families(family_id) VALUES(?)", family_batch
                )
                connection.executemany(
                    "INSERT INTO events(branch_id,family_id,event,parent_state,child_state) VALUES(?,?,?,?,?)",
                    event_batch,
                )
                family_batch.clear()
                event_batch.clear()
                if family_count % 10000 == 0:
                    progress(f"  processed {family_count:,} gene families")
        if family_batch:
            connection.executemany(
                "INSERT OR IGNORE INTO families(family_id) VALUES(?)", family_batch
            )
        if event_batch:
            connection.executemany(
                "INSERT INTO events(branch_id,family_id,event,parent_state,child_state) VALUES(?,?,?,?,?)",
                event_batch,
            )

    connection.executemany(
        "UPDATE branches SET gain_count=?, loss_count=? WHERE branch_id=?",
        [(counts[0], counts[1], branch) for branch, counts in branch_counts.items()],
    )
    _write_branches(output_path / "branches.tsv", branch_rows, branch_counts)

    phenotype_names: list[str] = []
    if phenotype_path:
        progress("Inferring phenotype transitions")
        phenotype_names, phenotype_data = read_phenotypes(
            phenotype_path, selected_phenotypes
        )
        with (output_path / "phenotype_gain_loss.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["phenotype_id", "branch_id", "event", "parent_state", "child_state"]
            )
            phenotype_events = []
            for phenotype_id in phenotype_names:
                transitions, _, _ = reconstruction.infer(
                    phenotype_data[phenotype_id], root_state=root_state
                )
                for transition in transitions:
                    row = (
                        phenotype_id,
                        transition.branch_id,
                        transition.event,
                        transition.parent_state,
                        transition.child_state,
                    )
                    phenotype_events.append(row)
                    writer.writerow(row)
            connection.executemany(
                "INSERT INTO phenotype_events(phenotype_id,branch_id,event,parent_state,child_state) VALUES(?,?,?,?,?)",
                phenotype_events,
            )

    if include_members and found["members"]:
        progress("Indexing gene IDs for interactive lookup")
        if input_format == "orthofinder":
            _import_members(connection, found["members"], progress)
        else:
            _import_pirate_members(connection, found["members"], progress)

    annotation_report: dict[str, int] | None = None
    imported_annotation: Path | None = None
    if annotate == "eggnog":
        annotation_dir = output_path / "annotations"
        if input_format == "orthofinder":
            representative_path = annotation_dir / "orthogroup_representatives.faa"
            progress("Selecting one representative protein per orthogroup")
            representative_count = write_representative_fasta(
                found["members"], proteomes, representative_path
            )
        else:
            representative_path = found["representatives"]
            with representative_path.open("r", encoding="utf-8") as handle:
                representative_count = sum(1 for line in handle if line.startswith(">"))
            progress("Using PIRATE representative_sequences.faa for annotation")
        progress(f"Running eggNOG-mapper for {representative_count:,} representatives")
        imported_annotation = run_eggnog_mapper(
            representative_path,
            annotation_dir,
            emapper=eggnog_emapper,
            data_dir=eggnog_data_dir,
            cpu=annotation_cpu,
        )
    elif annotation_path:
        annotation_dir = output_path / "annotations"
        annotation_dir.mkdir(parents=True, exist_ok=True)
        imported_annotation = annotation_dir / "eggnog.emapper.annotations"
        shutil.copyfile(annotation_path, imported_annotation)
    if imported_annotation:
        progress("Importing functional annotations")
        annotation_report = import_eggnog_annotations(connection, imported_annotation)
        if fetch_kegg_names:
            progress("Fetching official KEGG term names")
            kegg_report = fetch_official_kegg_names(
                connection,
                output_path / "annotations" / "kegg_term_names.tsv",
            )
            annotation_report["named_kegg_terms"] = kegg_report["total"]
        if go_obo_path:
            progress("Importing official Gene Ontology term names")
            go_names = import_go_term_names(connection, go_obo_path)
            annotation_report["named_go_terms"] = go_names
            go_copy = output_path / "annotations" / "go-basic.obo"
            if Path(go_obo_path).resolve() != go_copy.resolve():
                shutil.copyfile(go_obo_path, go_copy)
        _write_functional_annotations(connection, output_path / "functional_annotations.tsv")

    connection.executescript(
        """
        CREATE INDEX idx_events_branch ON events(branch_id, event);
        CREATE INDEX idx_events_family ON events(family_id, event);
        CREATE INDEX idx_family_terms_family ON family_terms(family_id);
        CREATE INDEX idx_family_terms_term ON family_terms(source, term_id, family_id);
        CREATE INDEX idx_phenotype_branch ON phenotype_events(phenotype_id, branch_id, event);
        """
    )
    progress("Building candidate rankings")
    _build_candidates(connection)
    _write_candidates(connection, output_path / "candidate_genes.tsv")
    connection.executescript(
        """
        CREATE INDEX idx_candidates_phenotype ON candidates(phenotype_id, score DESC);
        ANALYZE;
        """
    )
    connection.commit()
    connection.close()

    project = {
        "format_version": 1,
        "tool": "PhyloGain",
        "tool_version": __version__,
        "title": "Gene Gain/Loss Viewer",
        "input_format": input_format,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "species_count": len(species),
        "family_count": family_count,
        "phenotypes": phenotype_names,
        "nodes": as_project_nodes(root),
        "taxonomy": {
            "source": str(Path(gtdb_taxonomy_path).resolve()) if gtdb_taxonomy_path else None,
            "ranks": taxonomy_report["ranks"] if taxonomy_report else [],
            "mapped_species": len(taxonomy_data),
            "species": taxonomy_data,
        },
        "metadata": {
            "source": str(Path(genome_metadata_path).resolve())
            if genome_metadata_path
            else None,
            "genome_size_bp": genome_sizes,
        },
        "settings": {
            "gain_cost": gain_cost,
            "loss_cost": loss_cost,
            "root_state": root_state,
            "presence_threshold": presence_threshold,
            "gene_members_indexed": bool(include_members and found["members"]),
            "annotation_engine": "eggnog" if imported_annotation else None,
            "annotated_families": annotation_report["annotated_families"] if annotation_report else 0,
            "named_go_terms": annotation_report.get("named_go_terms", 0) if annotation_report else 0,
            "named_kegg_terms": annotation_report.get("named_kegg_terms", 0) if annotation_report else 0,
        },
    }
    (output_path / "project.json").write_text(
        json.dumps(project, ensure_ascii=False, separators=(",", ":")), encoding="utf-8"
    )
    metadata = {
        "tool_version": __version__,
        "python": sys.version,
        "platform": platform.platform(),
        "input_format": input_format,
        "orthofinder_directory": str(Path(orthofinder).resolve()) if orthofinder else None,
        "pirate_directory": str(Path(pirate).resolve()) if pirate else None,
        "tree_file": str(tree_path),
        "gene_count_file": str(found["counts"]),
        "gene_members_file": str(found["members"]) if found["members"] else None,
        "phenotype_file": str(Path(phenotype_path).resolve()) if phenotype_path else None,
        "gtdb_taxonomy_file": str(Path(gtdb_taxonomy_path).resolve()) if gtdb_taxonomy_path else None,
        "genome_metadata_file": str(Path(genome_metadata_path).resolve())
        if genome_metadata_path
        else None,
        "annotation_file": str(imported_annotation) if imported_annotation else None,
        "go_ontology_file": str(Path(go_obo_path).resolve()) if go_obo_path else None,
        "kegg_names_source": "KEGG REST API" if fetch_kegg_names else None,
        "proteomes_directory": str(Path(proteomes).resolve()) if proteomes else None,
        "warnings": report["warnings"],
    }
    (output_path / "run_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    template = files("species_innovation_map.templates").joinpath("index.html")
    shutil.copyfile(template, output_path / "index.html")
    return {
        "output": str(output_path),
        "species": len(species),
        "families": family_count,
        "branches": len(branch_rows),
        "phenotypes": phenotype_names,
        "taxonomy_species": len(taxonomy_data),
        "genome_metadata_species": len(genome_sizes),
        "annotated_families": annotation_report["annotated_families"] if annotation_report else 0,
        "warnings": report["warnings"],
    }


def _write_genome_metadata(
    path: Path, species: list[str], genome_sizes: dict[str, int]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species_id", "genome_size_bp", "genome_size_mb"])
        for species_id in species:
            if species_id not in genome_sizes:
                continue
            size_bp = genome_sizes[species_id]
            writer.writerow([species_id, size_bp, f"{size_bp / 1_000_000:.6f}"])


def _create_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE TABLE branches(
          branch_id TEXT PRIMARY KEY,
          parent_node TEXT,
          child_node TEXT NOT NULL,
          tip_label TEXT,
          depth INTEGER NOT NULL,
          branch_length REAL,
          gain_count INTEGER NOT NULL DEFAULT 0,
          loss_count INTEGER NOT NULL DEFAULT 0
        );
        CREATE TABLE families(
          family_id TEXT PRIMARY KEY,
          members_json BLOB,
          preferred_name TEXT NOT NULL DEFAULT '',
          description TEXT NOT NULL DEFAULT ''
        );
        CREATE TABLE events(
          branch_id TEXT NOT NULL,
          family_id TEXT NOT NULL,
          event TEXT NOT NULL,
          parent_state INTEGER NOT NULL,
          child_state INTEGER NOT NULL
        );
        CREATE TABLE phenotype_events(
          phenotype_id TEXT NOT NULL,
          branch_id TEXT NOT NULL,
          event TEXT NOT NULL,
          parent_state INTEGER NOT NULL,
          child_state INTEGER NOT NULL
        );
        CREATE TABLE candidates(
          phenotype_id TEXT NOT NULL,
          family_id TEXT NOT NULL,
          score INTEGER NOT NULL,
          coincident_gains INTEGER NOT NULL,
          phenotype_gains INTEGER NOT NULL,
          family_gains INTEGER NOT NULL,
          PRIMARY KEY(phenotype_id, family_id)
        );
        CREATE TABLE annotation_terms(
          source TEXT NOT NULL,
          term_id TEXT NOT NULL,
          term_name TEXT NOT NULL,
          family_count INTEGER NOT NULL DEFAULT 0,
          PRIMARY KEY(source, term_id)
        );
        CREATE TABLE family_terms(
          family_id TEXT NOT NULL,
          source TEXT NOT NULL,
          term_id TEXT NOT NULL,
          PRIMARY KEY(family_id, source, term_id)
        );
        """
    )


def _write_taxonomy(
    path: Path,
    species: list[str],
    taxonomy: dict[str, dict[str, str]],
) -> None:
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["species_id", *ranks])
        for species_id in species:
            values = taxonomy.get(species_id, {})
            writer.writerow([species_id, *(values.get(rank, "") for rank in ranks)])


def _write_branches(path: Path, rows, counts) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "branch_id",
                "parent_node",
                "child_node",
                "tip_label",
                "depth",
                "branch_length",
                "gain_count",
                "loss_count",
            ]
        )
        for row in rows:
            writer.writerow([*row, *counts[row[0]]])


def _import_members(connection: sqlite3.Connection, path: Path, progress) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            return
        species = header[1:]
        batch = []
        for row_number, row in enumerate(reader, start=2):
            if not row:
                continue
            members: dict[str, list[str]] = {}
            for species_id, cell in zip(species, row[1:]):
                genes = [gene.strip() for gene in cell.split(",") if gene.strip()]
                if genes:
                    members[species_id] = genes
            encoded = json.dumps(
                members, ensure_ascii=False, separators=(",", ":")
            ).encode("utf-8")
            batch.append((sqlite3.Binary(zlib.compress(encoded, level=6)), row[0]))
            if len(batch) >= 500:
                connection.executemany(
                    "UPDATE families SET members_json=? WHERE family_id=?", batch
                )
                batch.clear()
            if row_number % 10000 == 0:
                progress(f"  indexed {row_number - 1:,} gene families")
        if batch:
            connection.executemany(
                "UPDATE families SET members_json=? WHERE family_id=?", batch
            )


def _import_pirate_members(connection: sqlite3.Connection, path: Path, progress) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        species = header[PIRATE_METADATA_COLUMNS:]
        batch = []
        for row_number, row in enumerate(reader, start=2):
            if not row:
                continue
            members: dict[str, list[str]] = {}
            for species_id, cell in zip(species, row[PIRATE_METADATA_COLUMNS:]):
                genes = split_pirate_loci(cell)
                if genes:
                    members[species_id] = genes
            encoded = json.dumps(
                members, ensure_ascii=False, separators=(",", ":")
            ).encode("utf-8")
            preferred_name = row[2].strip()
            description = row[3].strip()
            batch.append(
                (
                    sqlite3.Binary(zlib.compress(encoded, level=6)),
                    "" if preferred_name in {"-", "NA"} else preferred_name,
                    "" if description in {"-", "NA"} else description,
                    row[1],
                )
            )
            if len(batch) >= 500:
                connection.executemany(
                    "UPDATE families SET members_json=?,preferred_name=?,description=? "
                    "WHERE family_id=?",
                    batch,
                )
                batch.clear()
            if row_number % 10000 == 0:
                progress(f"  indexed {row_number - 1:,} PIRATE gene families")
        if batch:
            connection.executemany(
                "UPDATE families SET members_json=?,preferred_name=?,description=? "
                "WHERE family_id=?",
                batch,
            )


def _build_candidates(connection: sqlite3.Connection) -> None:
    connection.execute(
        """
        INSERT INTO candidates(
          phenotype_id, family_id, score, coincident_gains, phenotype_gains, family_gains
        )
        SELECT
          pe.phenotype_id,
          e.family_id,
          COUNT(*) AS score,
          COUNT(*) AS coincident_gains,
          (SELECT COUNT(*) FROM phenotype_events p2
             WHERE p2.phenotype_id=pe.phenotype_id AND p2.event='gain'),
          (SELECT COUNT(*) FROM events e2
             WHERE e2.family_id=e.family_id AND e2.event='gain')
        FROM phenotype_events pe
        JOIN events e ON e.branch_id=pe.branch_id AND e.event='gain'
        WHERE pe.event='gain'
        GROUP BY pe.phenotype_id, e.family_id
        """
    )


def _write_candidates(connection: sqlite3.Connection, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "phenotype_id",
                "family_id",
                "score",
                "coincident_gains",
                "phenotype_gains",
                "family_gains",
            ]
        )
        writer.writerows(
            connection.execute(
                "SELECT phenotype_id,family_id,score,coincident_gains,phenotype_gains,family_gains FROM candidates ORDER BY phenotype_id,score DESC,family_id"
            )
        )


def _write_functional_annotations(connection: sqlite3.Connection, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["family_id", "preferred_name", "description", "source", "term_id", "term_name"]
        )
        writer.writerows(
            connection.execute(
                """
                SELECT f.family_id,f.preferred_name,f.description,
                       ft.source,ft.term_id,at.term_name
                FROM families f
                LEFT JOIN family_terms ft ON ft.family_id=f.family_id
                LEFT JOIN annotation_terms at ON at.source=ft.source AND at.term_id=ft.term_id
                ORDER BY f.family_id,ft.source,ft.term_id
                """
            )
        )
