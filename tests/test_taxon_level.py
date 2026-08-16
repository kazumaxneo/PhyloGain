from __future__ import annotations

import json
import sqlite3
import tempfile
import unittest
from contextlib import closing
from pathlib import Path

from species_innovation_map.aggregation import (
    aggregate_family_counts,
    classify_occupancy,
)
from species_innovation_map.analysis import Reconstruction
from species_innovation_map.models import TaxonGroup
from species_innovation_map.project import build_project
from species_innovation_map.tree import _Parser, assign_ids, leaf_labels
from species_innovation_map.tree_collapse import collapse_tree_to_taxa


HERE = Path(__file__).parent
ORTHOFINDER = HERE / "fixtures" / "orthofinder"
GTDB = HERE / "fixtures" / "gtdb_taxonomy.tsv"
GENOME_METADATA = HERE / "fixtures" / "genome_metadata.tsv"


def tree(text: str):
    root = _Parser(text).parse()
    assign_ids(root)
    return root


class TaxonAggregationTests(unittest.TestCase):
    def test_occupancy_boundary_is_present(self):
        group = TaxonGroup("g__A", "genus", tuple(f"A{i}" for i in range(10)))
        states = aggregate_family_counts(
            "OG1",
            list(group.member_genomes),
            [1] * 9 + [0],
            {group.taxon_id: group},
            present_threshold=0.9,
            absent_threshold=0.1,
        )
        self.assertEqual(states["g__A"].occupancy, 0.9)
        self.assertEqual(states["g__A"].observed_state, "present")

    def test_polymorphic_is_ambiguous_for_sankoff(self):
        self.assertEqual(classify_occupancy(0.5, 0.9, 0.1), "polymorphic")
        root = tree("(g__A,g__B);")
        reconstruction = Reconstruction(root)
        transitions, _, _, inferred = reconstruction.infer_states(
            {"g__A": frozenset({0, 1}), "g__B": frozenset({1})}
        )
        self.assertEqual(transitions, [])
        self.assertEqual(inferred["g__A"], 1)

    def test_confidence_mode_uses_sample_size(self):
        self.assertEqual(
            classify_occupancy(
                1.0, 0.9, 0.1, present_count=3, total_count=3,
                state_method="confidence", confidence=0.95,
            ),
            "polymorphic",
        )
        self.assertEqual(
            classify_occupancy(
                1.0, 0.9, 0.1, present_count=100, total_count=100,
                state_method="confidence", confidence=0.95,
            ),
            "present",
        )


class TaxonTreeTests(unittest.TestCase):
    def test_monophyletic_taxa_collapse(self):
        root = tree("((A1,A2),(B1,B2));")
        groups = {
            "g__A": TaxonGroup("g__A", "genus", ("A1", "A2")),
            "g__B": TaxonGroup("g__B", "genus", ("B1", "B2")),
        }
        result = collapse_tree_to_taxa(root, groups, min_genomes_per_taxon=2)
        self.assertEqual(leaf_labels(result.root), ["g__A", "g__B"])
        self.assertEqual(result.excluded_non_monophyletic, ())

    def test_non_monophyletic_taxon_is_not_forced(self):
        root = tree("(((A1,B1),(A2,B2)),(C1,C2));")
        groups = {
            "g__A": TaxonGroup("g__A", "genus", ("A1", "A2")),
            "g__C": TaxonGroup("g__C", "genus", ("C1", "C2")),
        }
        result = collapse_tree_to_taxa(root, groups, min_genomes_per_taxon=2)
        self.assertEqual(result.excluded_non_monophyletic, ("g__A",))
        self.assertNotIn("g__A", leaf_labels(result.root))

    def test_taxon_project_and_genome_backward_compatibility(self):
        with tempfile.TemporaryDirectory() as directory:
            genome_output = Path(directory) / "genome"
            build_project(ORTHOFINDER, genome_output)
            self.assertFalse((genome_output / "taxon_occupancy.tsv").exists())
            with closing(sqlite3.connect(genome_output / "species_map.sqlite")) as connection:
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM taxon_gene_states").fetchone()[0],
                    0,
                )

            taxon_output = Path(directory) / "taxon"
            build_project(
                ORTHOFINDER,
                taxon_output,
                gtdb_taxonomy_path=GTDB,
                taxon_rank="genus",
                min_genomes_per_taxon=1,
            )
            project = json.loads((taxon_output / "project.json").read_text(encoding="utf-8"))
            self.assertTrue(project["taxon_analysis"]["enabled"])
            self.assertEqual(project["taxon_analysis"]["rank"], "genus")
            self.assertEqual(project["species_count"], 3)
            with closing(sqlite3.connect(taxon_output / "species_map.sqlite")) as connection:
                row = connection.execute(
                    "SELECT occupancy,observed_state FROM taxon_gene_states "
                    "WHERE taxon_id='g__Genus_A' AND family_id='OG0000004'"
                ).fetchone()
            self.assertEqual(row, (0.5, "polymorphic"))

    def test_taxon_genome_size_is_mean_with_standard_deviation(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "taxon"
            build_project(
                ORTHOFINDER,
                output,
                gtdb_taxonomy_path=GTDB,
                genome_metadata_path=GENOME_METADATA,
                taxon_rank="genus",
                min_genomes_per_taxon=1,
            )
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            summary = project["metadata"]["genome_size_summary"]["g__Genus_A"]
            self.assertEqual(project["metadata"]["genome_size_bp"]["g__Genus_A"], 4_500_000)
            self.assertEqual(summary["mean_bp"], 4_500_000)
            self.assertEqual(summary["n"], 2)
            self.assertEqual(summary["sd_bp"], 500_000)
            html = (output / "index.html").read_text(encoding="utf-8")
            self.assertIn('class: "genome-size-error"', html)

    def test_taxon_bootstrap_writes_event_support(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "bootstrap"
            build_project(
                ORTHOFINDER,
                output,
                gtdb_taxonomy_path=GTDB,
                taxon_rank="genus",
                min_genomes_per_taxon=1,
                bootstrap_replicates=3,
                bootstrap_seed=7,
            )
            support = output / "event_bootstrap.tsv"
            self.assertTrue(support.exists())
            lines = support.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[0].split("\t")[-2:], ["replicates", "support"])
            self.assertTrue(any("\t3\t" in line for line in lines[1:]))
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            self.assertEqual(project["bootstrap"]["replicates"], 3)
            html = "Mean genome size: ${megabases.toFixed(2)} " + chr(177) + "}"
            self.assertIn('Mean genome size: ${megabases.toFixed(2)} ±', html)


if __name__ == "__main__":
    unittest.main()
