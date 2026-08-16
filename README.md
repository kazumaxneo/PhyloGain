# Gene Gain/Loss Viewer

**Gene Gain/Loss Viewer** turns an [OrthoFinder](https://github.com/davidemms/OrthoFinder) or [PIRATE](https://github.com/SionBayliss/PIRATE) result directory into an interactive species tree annotated with gene-family gains and losses. OrthoFinder is the recommended and primary input.

Each branch shows `+gains / −losses`. The viewer can switch these numbers to gain/loss pie charts: slice proportions show the relative gain/loss composition, while circle size follows the total event count with a capped threefold diameter range. Click a branch event to inspect its orthogroups, then click an orthogroup to see member gene IDs. An optional phenotype table marks inferred phenotype transitions and ranks gene families gained on the same branches.

> Alpha software. The current release uses weighted Sankoff parsimony on orthogroup presence/absence. Treat inferred events as hypotheses, especially for incomplete or contaminated genomes.

## Install

Python 3.10 or later is required.

```bash
pip install git+https://github.com/kazumaxneo/species-innovation-map.git
```

For development:

```bash
git clone https://github.com/kazumaxneo/species-innovation-map.git
cd species-innovation-map
pip install -e .
```

## Quick start

Build a map directly from an OrthoFinder results directory:

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --output species_innovation_map
```

Open it:

```bash
species-map serve species_innovation_map
```

The input options are mutually exclusive: specify exactly one of `--orthofinder` or `--pirate`.

For PIRATE, provide the PIRATE output directory and preferably a rooted external species tree. The tree tips must match the genome columns in `PIRATE.gene_families.tsv`.

```bash
species-map build \
  --pirate pyrites_ANI95_PIRATE_output \
  --species-tree rooted_core_gene_tree.nwk \
  --output pirate_gain_loss_map
```

If `--species-tree` is omitted, the program can use PIRATE's `binary_presence_absence.nwk`, but reports a warning because that tree is inferred from the same gene-content matrix being mapped. PIRATE consensus gene names, products, copy counts, duplicated loci, and parenthesized fission loci are retained in the interactive lookup.

The output directory contains:

- `index.html` — interactive viewer
- `branches.tsv` — gain/loss totals for every branch
- `gene_gain_loss.tsv` — inferred event for every affected orthogroup
- `species_map.sqlite` — indexed data used by the viewer
- `project.json` and `run_metadata.json` — tree, settings, and provenance

## Add genome-size metadata

Genome size is not contained in OrthoFinder or PIRATE output, so supply it as an optional tab-separated file. Use a tree-tip identifier column (`species_id`, `genome_id`, `user_genome`, or `assembly`) and either `genome_size_bp` or `genome_size_mb`.

```tsv
species_id	genome_size_bp
species_A	4000000
species_B	5000000
```

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --genome-metadata genome_metadata.tsv \
  --output map_with_genome_sizes
```

The Basic tab then provides `Genome size: Off / Bar graph`. The optional horizontal bars are aligned with tips in the rectangular layout and are included in the downloaded tree SVG. A collapsed clade shows the mean of its descendant genomes with available values. The normalized values are also written to `genome_metadata.tsv` in the output directory.

## Add phenotypes

Create a tab-separated file with `species_id` followed by one or more phenotype columns. Species IDs must match the OrthoFinder tree tips.

```tsv
species_id	nitrogen_fixation	heterocyst
species_A	+	+
species_B	+	+
species_C	-	?
species_D	-	-
```

Accepted states are `+`, `-`, `?`, `1`, `0`, `present`, `absent`, and `unknown`.

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --phenotypes phenotypes.tsv \
  --phenotype nitrogen_fixation \
  --output nitrogen_fixation_map

species-map serve nitrogen_fixation_map
```

This additionally creates:

- `phenotype_gain_loss.tsv` — inferred phenotype transitions
- `candidate_genes.tsv` — gene families gained on phenotype-gain branches

Repeat `--phenotype` to select several columns. Omit it to analyze all phenotype columns.

## Add functional enrichment with eggNOG-mapper

Functional annotation is optional. Supply the protein FASTA directory used by OrthoFinder and request eggNOG-mapper during the build:

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --annotate eggnog \
  --proteomes proteomes \
  --eggnog-data-dir /path/to/eggnog_data \
  --go-obo /path/to/go-basic.obo \
  --fetch-kegg-names \
  --annotation-cpu 16 \
  --output annotated_map

species-map serve annotated_map
```

The executable defaults to `emapper.py`. Use `--eggnog-emapper /path/to/emapper.py` when it is not on `PATH`. Gene Gain/Loss Viewer selects the first listed protein from each orthogroup as its representative, runs eggNOG-mapper once, and stores GO, KEGG, COG category, and Pfam assignments. Pass the official Gene Ontology `go-basic.obo` file with `--go-obo` to display GO names alongside GO identifiers. `--fetch-kegg-names` retrieves official KO, pathway, module, and reaction names from the KEGG REST API once during the build and stores them in SQLite and `annotations/kegg_term_names.tsv`; viewing the map then requires no KEGG API calls. The KEGG REST API is limited to academic use by academic users.

With PIRATE input, `--annotate eggnog` uses PIRATE's own `representative_sequences.faa`; `--proteomes` is not needed. Because PIRATE may retain several alleles for one family, all returned terms are merged into the corresponding PIRATE `gene_family`.

```bash
species-map build \
  --pirate pyrites_ANI95_PIRATE_output \
  --species-tree rooted_core_gene_tree.nwk \
  --annotate eggnog \
  --eggnog-data-dir /path/to/eggnog_data \
  --go-obo /path/to/go-basic.obo \
  --fetch-kegg-names \
  --output pirate_annotated_map
```

If eggNOG-mapper was run separately, import its standard output without repeating the search:

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --annotations eggnog.emapper.annotations \
  --output annotated_map
```

Clicking a branch calculates gained- and lost-family enrichment on demand. The database selector separates KEGG pathway, KEGG module, KEGG KO, GO, Pfam, COG category, and KEGG reaction results. The test uses all orthogroups in the map as the background, a one-sided Fisher exact test, and Benjamini-Hochberg correction within the selected database. Results with FDR <= 0.05 appear as a horizontal bar chart and exact table above the gene-family list. Bar length is the number of foreground gene families and bar color encodes FDR on a logarithmic red-to-blue scale, with lower FDR shown in red. The table reports foreground hits, background frequency, fold enrichment, and FDR. The output also includes `functional_annotations.tsv` and the original eggNOG-mapper file under `annotations/`.

Select `All` to rank significant terms from every available annotation database together while retaining separate Benjamini-Hochberg correction within each database. The plot selector switches between a bar plot (bar length = hits) and a dot plot (x = fold enrichment, dot size = hits); color represents FDR in both views.

## Collapse the tree by GTDB rank

Add a GTDB taxonomy TSV to enable interactive collapse controls for phylum, class, order, family, and genus. Standard GTDB-Tk summary files are accepted directly (`user_genome` plus `classification`). A generic TSV may instead use `species_id` or `genome_id` plus `gtdb_taxonomy`.

```tsv
species_id	gtdb_taxonomy
species_A	d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__Cyanobacteriales;f__Nostocaceae;g__Nostoc;s__Nostoc sp.
```

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --gtdb-taxonomy gtdbtk.bac120.summary.tsv \
  --output gtdb_collapsible_map
```

The viewer collapses each maximal monophyletic clade sharing the selected rank. A repeated taxon label means that taxon is non-monophyletic in the supplied species tree. Click a collapsed label to expand only that group. Gain/Loss values remain the events on the branch entering the collapsed clade; hidden descendant events are not added to that number.

The `GTDB labels` control can independently show family labels, genus labels, or both inside the uncollapsed tree. A label is placed only at the crown node of a monophyletic group with at least two sampled genomes; single-genome taxa remain represented only by their normal tip names. Select `Off` to hide all internal GTDB labels. The GTDB legend lists one colored taxon per line. Double-click GTDB or tip labels to rename them locally; `Restore defaults` in the Output tab clears renamed labels and manually moved GTDB label positions.

To use a different rooted topology, including a pruned and relabeled GTDB tree whose tips exactly match the OrthoFinder species IDs:

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --species-tree gtdb_pruned_rooted.nwk \
  --gtdb-taxonomy gtdbtk.bac120.summary.tsv \
  --output gtdb_tree_map
```

## OrthoFinder inputs

The program reads only standard OrthoFinder outputs:

```text
Results_xxx/
├── Species_Tree/
│   ├── SpeciesTree_rooted_node_labels.txt
│   └── SpeciesTree_rooted.txt                 # fallback
└── Orthogroups/
    ├── Orthogroups.GeneCount.tsv              # required
    └── Orthogroups.tsv                        # optional but recommended
```

Run input checks before a long analysis:

```bash
species-map validate \
  --orthofinder Results_Jul02 \
  --phenotypes phenotypes.tsv
```

## Main options

```text
--gain-cost FLOAT          Cost of a 0→1 transition (default: 1)
--loss-cost FLOAT          Cost of a 1→0 transition (default: 1)
--root-state auto|0|1      Root family/phenotype state (default: auto)
--presence-threshold INT   Copies required for presence (default: 1)
--orthofinder DIR          OrthoFinder results directory (recommended)
--pirate DIR               PIRATE output directory (exclusive with --orthofinder)
--no-members               Skip gene-family member IDs for a smaller output
--species-tree FILE        Rooted Newick tree replacing the input tool's default tree
--gtdb-taxonomy FILE       GTDB/GTDB-Tk taxonomy TSV for rank collapsing
--genome-metadata FILE     Optional TSV with genome size in bp or Mb
--annotate eggnog          Run eggNOG-mapper for orthogroup representatives
--annotations FILE         Import an existing .emapper.annotations file
--proteomes DIR            Protein FASTA directory used by OrthoFinder
--eggnog-emapper FILE      eggNOG-mapper executable (default: emapper.py)
--eggnog-data-dir DIR      Existing eggNOG database directory
--go-obo FILE              Official go-basic.obo file for readable GO term names
--fetch-kegg-names         Fetch and cache official KEGG names during the build
--annotation-cpu INT       CPUs used by eggNOG-mapper (default: 1)
```

Numeric internal-node support labels in the selected Newick tree are shown below
their branches on a 0–1 scale. Values supplied as percentages (for example, `95`)
are normalized to `0.950`; branches without support labels remain unlabelled.
The floating control window is divided into Basic, Filters, and Output tabs. Basic
contains analysis, tree, and display settings; Filters contains the independent
Gain/Loss ranges; Output contains reset and tree-SVG controls. The selected tab is
remembered in the browser.
The Zoom slider scales the complete tree while preserving the current viewport center.
Use H−/H+ to change the distance between successive tree depths and V−/V+ to
change the spacing between displayed tips. These controls are especially useful
after collapsing a large clade. Tip labels adjusts label text from 6 to 16 px.
Layout switches between Rectangular and Circular views without changing the
topology, branch events, or collapsed-node state. Circular tip names follow their
terminal-branch angles and reverse orientation on the left side for readability.
Branches switches between a default Cladogram with equal depth intervals and a
Phylogram whose horizontal or radial distances reflect non-negative branch lengths
from the selected Newick tree. Gain/Loss labels remain centered on their branches.
Download SVG saves the complete currently rendered tree, including its layout,
collapsed clades, Gain/Loss labels, support values, filters, and active color theme.
Each significant enrichment section also provides Download plot (SVG) and
Download table (UTF-8 TSV) for the selected branch, event type, and database.
Click a horizontal internal branch to collapse its descendant clade, and click the
collapsed marker to expand it again. Clicking a vertical connector flips the display
order of its child clades without changing the topology or inferred events.
Subtle hover tooltips identify these actions as Collapse clade, Expand clade, or Flip clades.
The Order control ladderizes every node by descendant-tip count in ascending or
descending order. The Support control can collapse branches below a selected
0–1 bootstrap threshold in 0.10 increments; its default All setting leaves the
complete tree visible.

Example with gains penalized relative to losses:

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --gain-cost 2 \
  --loss-cost 1 \
  --output weighted_map
```

## Interpretation and limitations

- A missing gene in an incomplete genome can be misidentified as a loss.
- Presence/absence parsimony does not explicitly model horizontal gene transfer.
- Candidate scores count coincident inferred gains; they are not statistical significance values.
- If a phenotype is defined from particular genes, exclude those known markers when interpreting novel candidates to avoid circular validation.
- The current version does not use gene-tree/species-tree reconciliation or copy-number evolution.

## License

MIT
