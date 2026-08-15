from __future__ import annotations

import json
import sqlite3
import tempfile
import unittest
import zlib
from contextlib import closing
from pathlib import Path
from unittest.mock import patch

from species_innovation_map.analysis import Reconstruction
from species_innovation_map.annotation import (
    annotation_sources,
    branch_enrichment,
    write_representative_fasta,
)
from species_innovation_map.project import build_project, validate_inputs
from species_innovation_map.tree import leaf_labels, parse_newick
from species_innovation_map.taxonomy import parse_gtdb_taxonomy, read_gtdb_taxonomy


HERE = Path(__file__).parent
FIXTURE = HERE / "fixtures" / "orthofinder"
PHENOTYPES = HERE / "fixtures" / "phenotypes.tsv"
GTDB_TAXONOMY = HERE / "fixtures" / "gtdb_taxonomy.tsv"
EGGNOG_ANNOTATIONS = HERE / "fixtures" / "eggnog.emapper.annotations"
PROTEOMES = HERE / "fixtures" / "proteomes"


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
            html = (output / "index.html").read_text(encoding="utf-8")
            self.assertIn('`${node.rankValue} (1)`', html)
            self.assertIn(".taxonomy-tip-label", html)
            self.assertIn('tipClass = rank && node.rankValue ? "taxonomy-tip-label"', html)

    def test_build_project_with_existing_eggnog_annotations(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "map"
            result = build_project(FIXTURE, output, annotation_path=EGGNOG_ANNOTATIONS)
            self.assertEqual(result["annotated_families"], 4)
            with closing(sqlite3.connect(output / "species_map.sqlite")) as connection:
                connection.row_factory = sqlite3.Row
                family = connection.execute(
                    "SELECT preferred_name,description FROM families WHERE family_id='OG0000001'"
                ).fetchone()
                self.assertEqual(family["preferred_name"], "nifH")
                terms = connection.execute(
                    "SELECT COUNT(*) FROM family_terms WHERE family_id='OG0000001'"
                ).fetchone()[0]
                self.assertGreaterEqual(terms, 5)
                branch = connection.execute(
                    "SELECT branch_id FROM events WHERE family_id='OG0000001' AND event='gain'"
                ).fetchone()[0]
                enrichment = branch_enrichment(connection, branch, "gain", min_overlap=1)
                self.assertTrue(enrichment["results"])
                self.assertIsNone(enrichment["source"])
                source = enrichment["results"][0]["source"]
                source_enrichment = branch_enrichment(
                    connection, branch, "gain", min_overlap=1, source=source
                )
                self.assertEqual(source_enrichment["source"], source)
                self.assertTrue(
                    all(row["source"] == source for row in source_enrichment["results"])
                )
                self.assertIn(source, [item["source"] for item in annotation_sources(connection)])
                with self.assertRaisesRegex(ValueError, "Unknown annotation database"):
                    branch_enrichment(
                        connection, branch, "gain", min_overlap=1, source="Not a database"
                    )
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            self.assertEqual(project["settings"]["annotation_engine"], "eggnog")
            self.assertTrue((output / "annotations" / "eggnog.emapper.annotations").is_file())
            self.assertTrue((output / "functional_annotations.tsv").is_file())
            html = (output / "index.html").read_text(encoding="utf-8")
            self.assertIn("Gained-family enrichment", html)
            self.assertIn("/api/enrichment?branch=", html)
            self.assertIn("Functional annotation", html)
            self.assertIn("FDR corrected within database", html)
            self.assertIn("/api/annotation-sources", html)
            self.assertIn("chart-axis", html)
            self.assertIn("Background", html)
            self.assertIn("bar length = hits", html)
            self.assertIn("createFdrLegend", html)
            self.assertIn("red = lower", html)
            self.assertIn("width: min(var(--detail-panel-width, 760px)", html)
            self.assertIn("bar.style.backgroundColor = fdrColor", html)
            self.assertIn("detailResizeHandle", html)
            self.assertIn("enableDetailResize", html)
            self.assertIn("background: var(--detail-panel)", html)
            self.assertIn('--detail-panel: #202428', html)
            self.assertIn('id="detailBack"', html)
            self.assertIn('configureDetailBack("Back"', html)
            self.assertIn("font-size: 9px", html)
            self.assertNotIn('textContent = "Inferred events"', html)
            self.assertNotIn('textContent = "Genes by species"', html)
            self.assertIn('assemblyNode.textContent = `Assembly: ${assembly}`', html)
            self.assertIn(".member-list { display: grid; gap: 0; margin-top: 32px; }", html)
            self.assertIn("text-decoration: underline", html)
            self.assertIn('content: "・"', html)
            self.assertIn('annotation.className = "function-annotation"', html)
            self.assertIn('functionTitle.className = "section-title functional-title"', html)
            self.assertIn("font-weight: 400", html)
            self.assertIn('classList.add("family-title")', html)
            self.assertIn('detailSubtitle").hidden = true', html)

    def test_build_project_runs_optional_eggnog_annotation(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "map"
            with patch(
                "species_innovation_map.project.run_eggnog_mapper",
                return_value=EGGNOG_ANNOTATIONS,
            ) as runner:
                result = build_project(
                    FIXTURE,
                    output,
                    annotate="eggnog",
                    proteomes=PROTEOMES,
                    annotation_cpu=3,
                )
            self.assertEqual(result["annotated_families"], 4)
            self.assertTrue((output / "annotations" / "orthogroup_representatives.faa").is_file())
            self.assertEqual(runner.call_args.kwargs["cpu"], 3)

    def test_write_orthogroup_representatives(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "representatives.faa"
            count = write_representative_fasta(
                FIXTURE / "Orthogroups" / "Orthogroups.tsv", PROTEOMES, fasta
            )
            self.assertEqual(count, 4)
            text = fasta.read_text(encoding="utf-8")
            self.assertIn(">OG0000001\n", text)
            self.assertNotIn(">gene_A1", text)


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
