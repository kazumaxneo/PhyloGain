# PhyloGain

PhyloGain is a visualization tool that uses **OrthoFinder results to infer and display gene family gains and losses across a phylogenetic tree**. It counts gain and loss events along each branch, making it easy to identify lineages with major changes in gene content. Users can also inspect the orthogroups associated with each gain or loss event and view their member gene IDs.

## Install

```bash
pip install git+https://github.com/kazumaxneo/PhyloGain.git
```

## Quick start

The standard workflow is **OrthoFinder → eggNOG-mapper → GTDB-Tk → PhyloGain**.

### 1. Run OrthoFinder

Place the protein FASTA files in one directory and run OrthoFinder.

```bash
orthofinder -f proteomes -t 16 -a 16
```

The directory created under `proteomes/OrthoFinder/` is used as the PhyloGain input.

### 2. Run eggNOG-mapper

Combine the same protein FASTA files used by OrthoFinder, then annotate them with eggNOG-mapper. Gene IDs must remain identical to those in the OrthoFinder results.

```bash
cat proteomes/*.faa > all_proteins.faa
emapper.py -i all_proteins.faa --itype proteins --output eggnog --cpu 16
curl -L https://purl.obolibrary.org/obo/go/go-basic.obo -o go-basic.obo
```

### 3. Run GTDB-Tk

Run GTDB-Tk on the corresponding genome assemblies.

```bash
gtdbtk classify_wf --genome_dir genomes --out_dir gtdbtk_output --extension fna --cpus 16
```

For bacterial genomes, use `gtdbtk_output/gtdbtk.bac120.summary.tsv` in the next step.

### 4. Run PhyloGain

```bash
phylogain build --orthofinder proteomes/OrthoFinder/Results_MmmDD --annotations eggnog.emapper.annotations --gtdb-taxonomy gtdbtk_output/gtdbtk.bac120.summary.tsv --go-obo go-basic.obo --fetch-kegg-names --output phylogain_output
phylogain serve phylogain_output
```

Replace `Results_MmmDD` with the actual OrthoFinder result directory. The viewer opens in a web browser.

## Add phenotype metadata

Phenotypes can be supplied as a tab-separated file. The first column must be `species_id`, and its values must match the tip names in the phylogenetic tree.

```text
species_id	nitrogen_fixation	heterocyst
species_A	+	+
species_B	+	-
species_C	-	-
species_D	?	?
```

Accepted states are `+`, `-`, `?`, `1`, `0`, `present`, `absent`, and `unknown`.

Add `--phenotypes` to the PhyloGain command. Use `--phenotype` to select a column; omit it to analyze every phenotype column.

```bash
phylogain build --orthofinder proteomes/OrthoFinder/Results_MmmDD --annotations eggnog.emapper.annotations --gtdb-taxonomy gtdbtk_output/gtdbtk.bac120.summary.tsv --go-obo go-basic.obo --fetch-kegg-names --phenotypes phenotypes.tsv --phenotype nitrogen_fixation --output phylogain_output
phylogain serve phylogain_output
```
