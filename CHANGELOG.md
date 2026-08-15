# Changelog

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
