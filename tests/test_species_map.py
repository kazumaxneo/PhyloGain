from __future__ import annotations

import json
import sqlite3
import tempfile
import unittest
import zlib
from contextlib import closing
from pathlib import Path

from species_innovation_map.analysis import Reconstruction
from species_innovation_map.project import build_project, validate_inputs
from species_innovation_map.tree import leaf_labels, parse_newick


HERE = Path(__file__).parent
FIXTURE = HERE / "fixtures" / "orthofinder"
PHENOTYPES = HERE / "fixtures" / "phenotypes.tsv"


class TreeTests(unittest.TestCase):
    def test_parse_tree(self):
        root = parse_newick(FIXTURE / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt")
        self.assertEqual(leaf_labels(root), ["species_A", "species_B", "species_C", "species_D"])
        self.assertEqual(len({node.branch_id for node in root.children}), 2)

    def test_sankoff_finds_single_clade_gain(self):
        root = parse_newick(FIXTURE / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt")
        result = Reconstruction(root).infer(
            {"species_A": 1, "species_B": 1, "species_C": 0, "species_D": 0}
        )
        transitions, root_state, score = result
        self.assertEqual(root_state, 0)
        self.assertEqual(score, 1)
        self.assertEqual([transition.event for transition in transitions], ["gain"])


class ProjectTests(unittest.TestCase):
    def test_validate_fixture(self):
        report = validate_inputs(FIXTURE, PHENOTYPES, ["nitrogen_fixation"])
        self.assertTrue(report["ok"])
        self.assertEqual(report["tree_tips"], 4)
        self.assertEqual(report["phenotypes"], ["nitrogen_fixation"])

    def test_build_project(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "map"
            result = build_project(
                FIXTURE,
                output,
                phenotype_path=PHENOTYPES,
                selected_phenotypes=["nitrogen_fixation"],
            )
            self.assertEqual(result["families"], 4)
            for name in (
                "index.html",
                "project.json",
                "species_map.sqlite",
                "branches.tsv",
                "gene_gain_loss.tsv",
                "phenotype_gain_loss.tsv",
                "candidate_genes.tsv",
            ):
                self.assertTrue((output / name).is_file(), name)
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            self.assertEqual(project["species_count"], 4)
            with closing(sqlite3.connect(output / "species_map.sqlite")) as connection:
                candidate = connection.execute(
                    "SELECT score FROM candidates WHERE phenotype_id=? AND family_id=?",
                    ("nitrogen_fixation", "OG0000001"),
                ).fetchone()
                self.assertEqual(candidate, (1,))
                members = connection.execute(
                    "SELECT members_json FROM families WHERE family_id='OG0000001'"
                ).fetchone()[0]
                self.assertIn("gene_A1", zlib.decompress(members).decode("utf-8"))


if __name__ == "__main__":
    unittest.main()
