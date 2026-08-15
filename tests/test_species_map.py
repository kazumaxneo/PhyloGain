from __future__ import annotations

import json
import sqlite3
import tempfile
import unittest
import zlib
from io import BytesIO
from contextlib import closing
from pathlib import Path
from unittest.mock import patch

from species_innovation_map.analysis import Reconstruction
from species_innovation_map.annotation import (
    annotation_sources,
    branch_enrichment,
    fetch_official_kegg_names,
    write_representative_fasta,
)
from species_innovation_map.project import build_project, validate_inputs
from species_innovation_map.tree import as_project_nodes, leaf_labels, parse_newick
from species_innovation_map.taxonomy import parse_gtdb_taxonomy, read_gtdb_taxonomy


HERE = Path(__file__).parent
FIXTURE = HERE / "fixtures" / "orthofinder"
PIRATE_FIXTURE = HERE / "fixtures" / "pirate"
PHENOTYPES = HERE / "fixtures" / "phenotypes.tsv"
GTDB_TAXONOMY = HERE / "fixtures" / "gtdb_taxonomy.tsv"
EGGNOG_ANNOTATIONS = HERE / "fixtures" / "eggnog.emapper.annotations"
PROTEOMES = HERE / "fixtures" / "proteomes"


class TreeTests(unittest.TestCase):
    def test_parse_tree(self):
        root = parse_newick(FIXTURE / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt")
        self.assertEqual(leaf_labels(root), ["species_A", "species_B", "species_C", "species_D"])
        self.assertEqual(len({node.branch_id for node in root.children}), 2)

    def test_parse_and_normalize_branch_support(self):
        with tempfile.TemporaryDirectory() as directory:
            tree_path = Path(directory) / "supported_tree.nwk"
            tree_path.write_text("((A:1,B:1)95:1,C:1)0.875;", encoding="utf-8")
            root = parse_newick(tree_path)
            self.assertEqual(root.support, 0.875)
            self.assertEqual(root.children[0].support, 0.95)
            self.assertIsNone(root.children[1].support)
            project_nodes = as_project_nodes(root)
            supported = {node["label"]: node["support"] for node in project_nodes}
            self.assertEqual(supported["95"], 0.95)

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
    def test_validate_requires_exactly_one_input_format(self):
        with self.assertRaisesRegex(ValueError, "exactly one"):
            validate_inputs(None)
        with self.assertRaisesRegex(ValueError, "exactly one"):
            validate_inputs(FIXTURE, pirate=PIRATE_FIXTURE)

    def test_build_project_from_pirate(self):
        report = validate_inputs(None, pirate=PIRATE_FIXTURE)
        self.assertTrue(report["ok"])
        self.assertEqual(report["input_format"], "pirate")
        self.assertTrue(any("gene-content tree" in warning for warning in report["warnings"]))
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "map"
            result = build_project(None, output, pirate=PIRATE_FIXTURE)
            self.assertEqual(result["families"], 4)
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            self.assertEqual(project["input_format"], "pirate")
            metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["input_format"], "pirate")
            self.assertIsNone(metadata["orthofinder_directory"])
            self.assertEqual(Path(metadata["pirate_directory"]), PIRATE_FIXTURE.resolve())
            with closing(sqlite3.connect(output / "species_map.sqlite")) as connection:
                family = connection.execute(
                    "SELECT preferred_name,description,members_json FROM families WHERE family_id='g000002'"
                ).fetchone()
                self.assertEqual(family[:2], ("trbI", "Bacterial conjugation TrbI-like protein"))
                members = json.loads(zlib.decompress(family[2]).decode("utf-8"))
                self.assertEqual(members["species_A"], ["gene_A2", "gene_A3"])
            html = (output / "index.html").read_text(encoding="utf-8")
            self.assertIn('Gene Gain/Loss Viewer · ${inputName}', html)

    def test_pirate_eggnog_alleles_map_back_to_gene_family(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            annotation = root / "pirate.emapper.annotations"
            annotation.write_text(
                "#query\tPreferred_name\tDescription\tGOs\tKEGG_ko\tKEGG_Pathway\t"
                "KEGG_Module\tKEGG_Reaction\tCOG_category\tPFAMs\n"
                "g000001_1\tnifH\tnitrogenase iron protein\tGO:0009399\tko:K02588\t"
                "map00910\tM00175\tR00001\tC\tFer4_NifH\n",
                encoding="utf-8",
            )
            output = root / "map"
            result = build_project(
                None,
                output,
                pirate=PIRATE_FIXTURE,
                annotation_path=annotation,
            )
            self.assertEqual(result["annotated_families"], 1)
            with closing(sqlite3.connect(output / "species_map.sqlite")) as connection:
                terms = connection.execute(
                    "SELECT COUNT(*) FROM family_terms WHERE family_id='g000001'"
                ).fetchone()[0]
                self.assertGreater(terms, 0)

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
                all_enrichment = branch_enrichment(
                    connection, branch, "gain", limit=100, min_overlap=1, source="All"
                )
                self.assertEqual(all_enrichment["source"], "All")
                self.assertGreater(len({row["source"] for row in all_enrichment["results"]}), 1)
                ko_enrichment = branch_enrichment(
                    connection, branch, "gain", min_overlap=1, source="KEGG KO"
                )
                self.assertEqual(ko_enrichment["results"][0]["term_id"], "K02588")
                self.assertEqual(
                    ko_enrichment["results"][0]["term_name"], "nitrogenase iron protein"
                )
                with self.assertRaisesRegex(ValueError, "Unknown annotation database"):
                    branch_enrichment(
                        connection, branch, "gain", min_overlap=1, source="Not a database"
                    )
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            self.assertEqual(project["settings"]["annotation_engine"], "eggnog")
            self.assertTrue((output / "annotations" / "eggnog.emapper.annotations").is_file())
            self.assertTrue((output / "functional_annotations.tsv").is_file())
            html = (output / "index.html").read_text(encoding="utf-8")
            self.assertIn("Functional enrichment of gained gene families", html)
            self.assertIn("Functional enrichment of lost gene families", html)
            self.assertIn("Gene Gain/Loss Viewer", html)
            self.assertNotIn("Species Innovation Map", html)
            self.assertIn('<div class="summary" aria-label="Dataset summary">', html)
            self.assertNotIn('<section class="summary">', html)
            self.assertIn('.card span::after { content: ":"; }', html)
            self.assertIn("flex: 1 1 auto", html)
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
            self.assertIn('allOption.value = "All"', html)
            self.assertIn('["dot", "Dot plot"]', html)
            self.assertIn('dot.className = "chart-dot"', html)
            self.assertIn('x = fold enrichment · dot size = hits', html)
            self.assertIn('state.enrichmentPlotType', html)
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
            self.assertIn('eventLabel.className = "event-kind"', html)
            self.assertIn("grid-template-columns: 34px 92px minmax(0, 1fr)", html)
            self.assertIn("text-overflow: ellipsis", html)
            self.assertIn("enrichmentTermLabel(row)", html)
            self.assertIn('enrichment.source === "KEGG reaction"', html)
            self.assertIn("--term-column: minmax(280px, 38%)", html)
            self.assertIn("--term-column: minmax(340px, 48%)", html)
            self.assertIn('id="gainMin"', html)
            self.assertIn('id="gainMax"', html)
            self.assertIn('id="lossMin"', html)
            self.assertIn('id="lossMax"', html)
            self.assertIn("function branchMatchesFilter(branch)", html)
            self.assertIn('branchValueMatchesFilter(branch, "gain")', html)
            self.assertIn('branchValueMatchesFilter(branch, "loss")', html)
            self.assertIn('class: gainClass', html)
            self.assertIn('class: lossClass', html)
            self.assertIn(".range-filter output { grid-column: 2; grid-row: 1; color: var(--ink); font: 11px ui-monospace", html)
            self.assertIn('"data-event": "gain"', html)
            self.assertIn('"data-event": "loss"', html)
            self.assertIn("async function selectBranch(branchId, eventKind = null)", html)
            self.assertIn("const eventKinds = eventKind ? [eventKind]", html)
            self.assertIn('eventKind === "gain" ? "Gained gene families"', html)
            self.assertIn("selectBranch(pendingBranch, pendingEvent)", html)
            self.assertIn("manuallyCollapsed: new Set()", html)
            self.assertIn("flippedNodes: new Set()", html)
            self.assertIn('"data-collapse-node": node.id', html)
            self.assertIn('"data-flip-node": node.id', html)
            self.assertIn('svgEl("rect", hitAttributes)', html)
            self.assertIn("Collapsed clade", html)
            self.assertIn("state.flippedNodes.has(node.id)", html)
            self.assertIn('addSvgTooltip(connectorHit, "Flip clades")', html)
            self.assertIn('node.collapsed ? "Expand clade" : "Collapse clade"', html)
            self.assertIn('addSvgTooltip(marker, "Expand clade")', html)
            self.assertIn('id="treeActionHint"', html)
            self.assertIn("node.dataset.actionHint = text", html)
            self.assertIn(".tree-action-hint { position: absolute", html)
            self.assertNotIn('const title = svgEl("title")', html)
            self.assertIn('id="treeOrderSelect"', html)
            self.assertIn('id="supportCollapseSelect"', html)
            self.assertIn('value="original">Original', html)
            self.assertIn('value="0">All', html)
            for threshold in range(1, 11):
                self.assertIn(f'&lt;{threshold / 10:.2f}</option>', html)
            self.assertIn("state.treeOrder !== \"original\"", html)
            self.assertIn("node.supportCollapsed = Boolean", html)
            self.assertIn("function expandCollapsedNode(nodeId)", html)
            self.assertIn("font-size: 15px", html)
            self.assertIn('row.source === "COG category"', html)
            self.assertNotIn('id="hideUnmatched"', html)
            self.assertIn("filtered-branch-label", html)
            self.assertIn("normalizedBranchSupport(node)", html)
            self.assertIn('class: "support-label"', html)
            self.assertIn("support.toFixed(3)", html)
            self.assertIn('id="treeZoom"', html)
            self.assertIn("function enableTreeZoom()", html)
            self.assertIn('transform: `scale(${zoom})`', html)
            self.assertIn("centerX * state.treeZoom", html)
            self.assertIn('id="treeLayoutSelect"', html)
            self.assertIn('value="rectangular">Rectangular', html)
            self.assertIn('value="circular">Circular', html)
            self.assertNotIn('value="unrooted">Unrooted', html)
            self.assertNotIn('state.treeLayout === "unrooted"', html)
            self.assertIn('id="horizontalSpacing" type="range"', html)
            self.assertIn('id="verticalSpacing" type="range"', html)
            self.assertIn('id="tipFontSize" type="range"', html)
            self.assertIn('transform: `rotate(${labelAngle}', html)
            self.assertIn("function enableTreeLayoutControls()", html)
            self.assertIn('id="downloadSvg"', html)
            self.assertIn("function downloadCurrentTreeSvg()", html)
            self.assertIn('image/svg+xml;charset=utf-8', html)
            self.assertIn('gene-gain-loss-tree-${state.treeLayout}.svg', html)
            self.assertIn('Download plot', html)
            self.assertIn('Download table', html)
            self.assertIn('function downloadEnrichmentPlot(', html)
            self.assertIn('legendDirection.textContent = "red = lower"', html)
            self.assertIn('legendTitle.textContent = "FDR"', html)
            self.assertIn('axisTitle.textContent = plotType === "dot" ? "Fold enrichment" : "Count"', html)
            self.assertIn('function downloadEnrichmentTable(', html)
            self.assertIn('text/tab-separated-values;charset=utf-8', html)
            self.assertIn('state.treeLayout !== "rectangular"', html)
            self.assertIn('class: "radial-branch-hit"', html)
            self.assertIn("25 * state.verticalSpacing", html)
            self.assertIn("110 * state.horizontalSpacing", html)
            self.assertIn('<div class="header-controls" aria-label="Tree controls">', html)
            self.assertNotIn('<div class="toolbar">', html)
            self.assertIn('<div class="tree-tools">', html)
            self.assertIn('id="filterReset"', html)
            self.assertIn(".gain-filter { grid-column: 1; grid-row: 1; }", html)
            self.assertIn(".loss-filter { grid-column: 2; grid-row: 1; }", html)
            self.assertIn('<span class="range-filter-name">Filter (gain)</span>', html)
            self.assertIn('<span class="range-filter-name">Filter (loss)</span>', html)
            self.assertIn('<label for="horizontalSpacing">Horizontal</label>', html)
            self.assertIn('<label for="verticalSpacing">Vertical</label>', html)
            self.assertIn('<label for="cladeLabelSelect">GTDB labels</label>', html)
            self.assertIn('option.value = "family,genus"', html)
            self.assertIn('function cladeTaxonomyLabels(node)', html)
            self.assertIn('visibleNodes.forEach(node => appendCladeTaxonomyLabel(node, true))', html)
            self.assertIn('visibleNodes.forEach(node => appendCladeTaxonomyLabel(node))', html)
            self.assertIn('html[data-theme="black"] .range-filter-name', html)
            self.assertIn("grid-template-rows: 15px 18px", html)
            self.assertIn("min-width: 235px", html)
            self.assertIn('<div class="filter-status">', html)
            self.assertIn("border: 0; border-radius: 0", html)
            self.assertIn("officialTermUrl(source, termId)", html)
            self.assertIn("https://www.kegg.jp/entry/", html)
            self.assertIn("https://www.kegg.jp/module/", html)
            self.assertIn("https://www.kegg.jp/pathway/", html)
            self.assertIn("https://amigo.geneontology.org/amigo/term/", html)
            self.assertIn('link.target = "_blank"', html)

    def test_import_official_go_names(self):
        with tempfile.TemporaryDirectory() as directory:
            obo = Path(directory) / "go-basic.obo"
            obo.write_text(
                "format-version: 1.2\n\n"
                "[Term]\n"
                "id: GO:0009399\n"
                "name: nitrogen fixation\n\n"
                "[Term]\n"
                "id: GO:0003677\n"
                "name: DNA binding\n",
                encoding="utf-8",
            )
            output = Path(directory) / "map"
            build_project(
                FIXTURE,
                output,
                annotation_path=EGGNOG_ANNOTATIONS,
                go_obo_path=obo,
            )
            with closing(sqlite3.connect(output / "species_map.sqlite")) as connection:
                name = connection.execute(
                    "SELECT term_name FROM annotation_terms WHERE source='GO' AND term_id='GO:0009399'"
                ).fetchone()[0]
                self.assertEqual(name, "nitrogen fixation")
            self.assertTrue((output / "annotations" / "go-basic.obo").is_file())
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            self.assertEqual(project["settings"]["named_go_terms"], 2)

    def test_fetch_official_kegg_names_once_and_cache(self):
        responses = {
            "ko": "K02588\tnifH; nitrogenase iron protein\n",
            "pathway": "map00910\tNitrogen metabolism\n",
            "module": "M00175\tNitrogen fixation, nitrogen => ammonia\n",
            "reaction": "R00001\tpolyphosphate polyphosphohydrolase\n",
        }

        def fake_urlopen(request, timeout):
            database = request.full_url.rsplit("/", 1)[-1]
            return BytesIO(responses[database].encode("utf-8"))

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "map"
            with patch("species_innovation_map.annotation.urlopen", side_effect=fake_urlopen):
                with patch("species_innovation_map.annotation.time.sleep"):
                    build_project(
                        FIXTURE,
                        output,
                        annotation_path=EGGNOG_ANNOTATIONS,
                        fetch_kegg_names=True,
                    )
            with closing(sqlite3.connect(output / "species_map.sqlite")) as connection:
                name = connection.execute(
                    "SELECT term_name FROM annotation_terms WHERE source='KEGG KO' AND term_id='K02588'"
                ).fetchone()[0]
                self.assertEqual(name, "nifH; nitrogenase iron protein")
            cache = output / "annotations" / "kegg_term_names.tsv"
            self.assertIn("KEGG KO\tK02588\tnifH", cache.read_text(encoding="utf-8"))
            project = json.loads((output / "project.json").read_text(encoding="utf-8"))
            self.assertEqual(project["settings"]["named_kegg_terms"], 4)

    def test_retired_kegg_modules_are_identified(self):
        responses = {
            "ko": "",
            "pathway": "",
            "module": "M00175\tNitrogen fixation, nitrogen => ammonia\n",
            "reaction": "",
        }

        def fake_urlopen(request, timeout):
            database = request.full_url.rsplit("/", 1)[-1]
            return BytesIO(responses[database].encode("utf-8"))

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "map"
            build_project(FIXTURE, output, annotation_path=EGGNOG_ANNOTATIONS)
            with closing(sqlite3.connect(output / "species_map.sqlite")) as connection:
                connection.execute(
                    "INSERT INTO annotation_terms(source,term_id,term_name,family_count) VALUES('KEGG module','M00207','M00207',1)"
                )
                with patch("species_innovation_map.annotation.urlopen", side_effect=fake_urlopen):
                    with patch("species_innovation_map.annotation.time.sleep"):
                        report = fetch_official_kegg_names(connection)
                label = connection.execute(
                    "SELECT term_name FROM annotation_terms WHERE source='KEGG module' AND term_id='M00207'"
                ).fetchone()[0]
                self.assertIn("Retired KEGG module", label)
                self.assertEqual(report["retired_modules"], 1)

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
    def test_read_normalized_rank_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "taxonomy.tsv"
            path.write_text(
                "species_id\tdomain\tfamily\tgenus\tspecies\n"
                "species_A\tBacteria\tFamily_A\tGenus_A\tGenus_A species_A\n",
                encoding="utf-8",
            )
            mapped, report = read_gtdb_taxonomy(path, ["species_A"])
            self.assertEqual(mapped["species_A"]["genus"], "Genus_A")
            self.assertEqual(report["mapped_species"], 1)

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
