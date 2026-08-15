# Changelog

## 0.4.0 — 2026-08-15

- Add a database selector for KEGG pathway/module/KO, GO, Pfam, COG, and KEGG reaction enrichment.
- Correct multiple testing separately within the selected annotation database.
- Redesign enrichment charts with gene-family count on the axis and FDR encoded by a labelled red-to-blue color scale.
- Add background frequencies to the exact enrichment tables.
- Make the branch-detail panel width draggable and persistent.
- Give the detail panel a subtly distinct background in every color theme.
- Add compact orthogroup member lists and a back button to return to the originating branch or candidate list.
- Remove the redundant inferred-event branch list from individual orthogroup views.
- Label orthogroup member groups as assemblies instead of species and remove the redundant member-list heading.
- Separate assembly members from annotations with compact, underlined, theme-neutral assembly labels.
- Prefix each functional annotation item with a middle dot and increase spacing before assembly members.
- Use colon-separated, regular-weight assembly labels and emphasize functional annotation headings and items.
- Unify functional-annotation typography, simplify orthogroup headers, and emphasize the Back button.
- Replace boxed gene-family event buttons with compact single-line gain/loss, orthogroup, and annotation rows.
- Fetch and cache official KEGG KO, pathway, module, and reaction names at build time.
- Link KEGG and GO term identifiers in branch and orthogroup panels to their official database pages.
- Mark legacy module IDs missing from the current KEGG release instead of leaving unexplained bare IDs.
- Give KEGG reaction charts a wider label column and a more compact plotting area for long reaction descriptions.
- Widen term-label columns across all enrichment databases while keeping KEGG reaction labels widest.
- Add independent Gain and Loss range filters that fade only unmatched count labels while preserving every branch, node, and tip.
- Omit single-letter COG category codes from chart labels while retaining them in the result table.
- Stack Gain and Loss filter bars vertically, enlarge their labels and values, and simplify Reset to an unboxed text control.
- Rename the public interface from Species Innovation Map to Gene Gain/Loss Viewer.
- Move the four dataset statistics into the title bar to preserve vertical space on laptop screens.
- Clarify the gained- and lost-family functional-enrichment headings in the branch detail panel.
- Add `--go-obo` support and display official GO term names, including obsolete terms retained by older annotations.

## 0.3.1 — 2026-08-15

- Widen the branch-detail panel for enrichment results.
- Add ranked horizontal FDR charts above the exact enrichment tables.
- Expand enrichment tables with source, hit count, fold enrichment, and FDR columns.

## 0.3.0 — 2026-08-15

- Add optional eggNOG-mapper annotation with `--annotate eggnog`.
- Import precomputed eggNOG results with `--annotations`.
- Calculate branch-specific Gain/Loss enrichment with Fisher's exact test and BH FDR.
- Show enriched terms, preferred gene names, and functional descriptions in the branch panel.
- Export `functional_annotations.tsv` for reuse outside the viewer.

## 0.2.1 — 2026-08-15

- Show taxonomy-assigned singleton tips as `Taxon (1)` in collapsed rank views.
- Keep genome IDs only for tips without an assignment at the selected rank.

## 0.2.0 — 2026-08-15

- Accept optional GTDB/GTDB-Tk taxonomy TSV files.
- Collapse monophyletic clades interactively by phylum, class, order, family, or genus.
- Expand individual collapsed groups without leaving the selected rank view.
- Accept an optional rooted Newick tree in place of the OrthoFinder species tree.
- Add selectable interface themes, a sliding branch-details panel, and drag-to-pan navigation.
- Simplify branch labels and Gain/Loss summaries for readability.

## 0.1.0 — 2026-08-15

- Read rooted species trees and orthogroup tables directly from OrthoFinder results.
- Infer branch-level gene-family gains and losses with weighted Sankoff parsimony.
- Accept optional multi-phenotype TSV input and infer phenotype transitions.
- Rank gene families gained on phenotype-gain branches.
- Generate an interactive species-tree viewer backed by a compressed SQLite index.
- Export branch, gene-family, phenotype, and candidate tables as TSV.
