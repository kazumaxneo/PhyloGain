# Species Innovation Map

**Species Innovation Map** turns an [OrthoFinder](https://github.com/davidemms/OrthoFinder) result directory into an interactive species tree annotated with gene-family gains and losses.

Each branch shows `+gains / −losses`. Click a branch to inspect its orthogroups, then click an orthogroup to see member gene IDs. An optional phenotype table marks inferred phenotype transitions and ranks gene families gained on the same branches.

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

The output directory contains:

- `index.html` — interactive viewer
- `branches.tsv` — gain/loss totals for every branch
- `gene_gain_loss.tsv` — inferred event for every affected orthogroup
- `species_map.sqlite` — indexed data used by the viewer
- `project.json` and `run_metadata.json` — tree, settings, and provenance

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

The executable defaults to `emapper.py`. Use `--eggnog-emapper /path/to/emapper.py` when it is not on `PATH`. Species Innovation Map selects the first listed protein from each orthogroup as its representative, runs eggNOG-mapper once, and stores GO, KEGG, COG category, and Pfam assignments. Pass the official Gene Ontology `go-basic.obo` file with `--go-obo` to display GO names alongside GO identifiers. `--fetch-kegg-names` retrieves official KO, pathway, module, and reaction names from the KEGG REST API once during the build and stores them in SQLite and `annotations/kegg_term_names.tsv`; viewing the map then requires no KEGG API calls. The KEGG REST API is limited to academic use by academic users.

If eggNOG-mapper was run separately, import its standard output without repeating the search:

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --annotations eggnog.emapper.annotations \
  --output annotated_map
```

Clicking a branch calculates gained- and lost-family enrichment on demand. The database selector separates KEGG pathway, KEGG module, KEGG KO, GO, Pfam, COG category, and KEGG reaction results. The test uses all orthogroups in the map as the background, a one-sided Fisher exact test, and Benjamini-Hochberg correction within the selected database. Results with FDR <= 0.05 appear as a horizontal bar chart and exact table above the gene-family list. Bar length is the number of foreground gene families and bar color encodes FDR on a logarithmic red-to-blue scale, with lower FDR shown in red. The table reports foreground hits, background frequency, fold enrichment, and FDR. The output also includes `functional_annotations.tsv` and the original eggNOG-mapper file under `annotations/`.

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
--no-members               Skip Orthogroups.tsv to make a smaller output
--species-tree FILE        Rooted Newick tree replacing the OrthoFinder tree
--gtdb-taxonomy FILE       GTDB/GTDB-Tk taxonomy TSV for rank collapsing
--annotate eggnog          Run eggNOG-mapper for orthogroup representatives
--annotations FILE         Import an existing .emapper.annotations file
--proteomes DIR            Protein FASTA directory used by OrthoFinder
--eggnog-emapper FILE      eggNOG-mapper executable (default: emapper.py)
--eggnog-data-dir DIR      Existing eggNOG database directory
--go-obo FILE              Official go-basic.obo file for readable GO term names
--fetch-kegg-names         Fetch and cache official KEGG names during the build
--annotation-cpu INT       CPUs used by eggNOG-mapper (default: 1)
```

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
