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
from species_innovation_map.taxonomy import parse_gtdb_taxonomy, read_gtdb_taxonomy


HERE = Path(__file__).parent
FIXTURE = HERE / "fixtures" / "orthofinder"
PHENOTYPES = HERE / "fixtures" / "phenotypes.tsv"
GTDB_TAXONOMY = HERE / "fixtures" / "gtdb_taxonomy.tsv"


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

    def test_build_project_with_gtdb_taxonomy(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "map"
            result = build_project(FIXTURE, output, gtdb_taxonomy_path=GTDB_TAXONOMY)
            self.assertEqual(result["taxonomy_species"], 4)
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            self.assertEqual(project["taxonomy"]["mapped_species"], 4)
            self.assertEqual(project["taxonomy"]["species"]["species_A"]["genus"], "Genus_A")
            self.assertIn("family", project["taxonomy"]["ranks"])
            self.assertTrue((output / "gtdb_taxonomy.tsv").is_file())


class TaxonomyTests(unittest.TestCase):
    def test_parse_and_match_taxonomy(self):
        parsed = parse_gtdb_taxonomy("d__Bacteria;f__Nostocaceae;g__Nostoc;s__")
        self.assertEqual(parsed, {"domain": "Bacteria", "family": "Nostocaceae", "genus": "Nostoc"})
        taxonomy, report = read_gtdb_taxonomy(
            GTDB_TAXONOMY, ["species_A", "species_B", "species_C", "species_D"]
        )
        self.assertEqual(report["mapped_species"], 4)
        self.assertEqual(taxonomy["species_D"]["family"], "Family_B")

    def test_match_gtdbtk_summary_accession(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gtdbtk.bac120.summary.tsv"
            path.write_text(
                "user_genome\tclassification\n"
                "GCA_123456789\td__Bacteria;p__Cyanobacteriota;f__Nostocaceae;g__Nostoc\n",
                encoding="utf-8",
            )
            taxonomy, report = read_gtdb_taxonomy(
                path, ["Nostocaceae__GCA_123456789.1_GTDB"]
            )
            self.assertEqual(report["mapped_species"], 1)
            self.assertEqual(
                taxonomy["Nostocaceae__GCA_123456789.1_GTDB"]["genus"], "Nostoc"
            )


if __name__ == "__main__":
    unittest.main()
