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
    Reconstruction,
    find_orthofinder_files,
    iter_gene_counts,
    read_count_header,
    read_phenotypes,
)
from .tree import as_project_nodes, leaf_labels, parse_newick, preorder
from .taxonomy import read_gtdb_taxonomy


class InputError(ValueError):
    pass


def validate_inputs(
    orthofinder: str | Path,
    phenotype_path: str | Path | None = None,
    selected_phenotypes: list[str] | None = None,
    species_tree_path: str | Path | None = None,
    gtdb_taxonomy_path: str | Path | None = None,
) -> dict[str, object]:
    found = find_orthofinder_files(orthofinder)
    errors: list[str] = []
    warnings: list[str] = []
    tree_path = Path(species_tree_path) if species_tree_path else found["tree"]
    if tree_path is None or not tree_path.is_file():
        errors.append("No rooted species tree found in Species_Tree/")
    if found["counts"] is None:
        errors.append("Orthogroups/Orthogroups.GeneCount.tsv was not found")
    if errors:
        return {"ok": False, "errors": errors, "warnings": warnings}

    root = parse_newick(tree_path)
    tips = leaf_labels(root)
    if not all(tips):
        errors.append("Every species-tree tip must have a label")
    duplicate_tips = sorted(name for name, count in Counter(tips).items() if count > 1)
    if duplicate_tips:
        errors.append(f"Duplicate species-tree tips: {', '.join(duplicate_tips[:10])}")

    count_species, _ = read_count_header(found["counts"])
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
        warnings.append(
            "Orthogroups/Orthogroups.tsv was not found; gene IDs cannot be shown in the HTML"
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
    }


def build_project(
    orthofinder: str | Path,
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
    progress=None,
) -> dict[str, object]:
    progress = progress or (lambda message: None)
    report = validate_inputs(
        orthofinder,
        phenotype_path,
        selected_phenotypes,
        species_tree_path,
        gtdb_taxonomy_path,
    )
    if not report["ok"]:
        raise InputError("; ".join(report["errors"]))
    found = find_orthofinder_files(orthofinder)
    output_path = Path(output).resolve()
    if output_path.exists() and any(output_path.iterdir()):
        raise InputError(f"Output directory is not empty: {output_path}")
    output_path.mkdir(parents=True, exist_ok=True)

    tree_path = Path(report["tree"])
    root = parse_newick(tree_path)
    nodes = preorder(root)
    species = leaf_labels(root)
    count_species, _ = read_count_header(found["counts"])
    reconstruction = Reconstruction(root, gain_cost=gain_cost, loss_cost=loss_cost)
    taxonomy_data: dict[str, dict[str, str]] = {}
    taxonomy_report: dict[str, object] | None = None
    if gtdb_taxonomy_path:
        taxonomy_data, taxonomy_report = read_gtdb_taxonomy(gtdb_taxonomy_path, species)
        _write_taxonomy(output_path / "gtdb_taxonomy.tsv", species, taxonomy_data)

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
        for family_id, counts in iter_gene_counts(found["counts"]):
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
        _import_members(connection, found["members"], progress)

    connection.executescript(
        """
        CREATE INDEX idx_events_branch ON events(branch_id, event);
        CREATE INDEX idx_events_family ON events(family_id, event);
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
        "tool": "species-innovation-map",
        "tool_version": __version__,
        "title": "Species Innovation Map",
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
        "settings": {
            "gain_cost": gain_cost,
            "loss_cost": loss_cost,
            "root_state": root_state,
            "presence_threshold": presence_threshold,
            "gene_members_indexed": bool(include_members and found["members"]),
        },
    }
    (output_path / "project.json").write_text(
        json.dumps(project, ensure_ascii=False, separators=(",", ":")), encoding="utf-8"
    )
    metadata = {
        "tool_version": __version__,
        "python": sys.version,
        "platform": platform.platform(),
        "orthofinder_directory": str(Path(orthofinder).resolve()),
        "tree_file": str(tree_path),
        "gene_count_file": str(found["counts"]),
        "gene_members_file": str(found["members"]) if found["members"] else None,
        "phenotype_file": str(Path(phenotype_path).resolve()) if phenotype_path else None,
        "gtdb_taxonomy_file": str(Path(gtdb_taxonomy_path).resolve()) if gtdb_taxonomy_path else None,
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
        "warnings": report["warnings"],
    }


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
        CREATE TABLE families(family_id TEXT PRIMARY KEY, members_json BLOB);
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
