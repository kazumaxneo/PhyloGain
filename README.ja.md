# PhyloGain

OrthoFinderの結果から、species treeの各branchにgene-familyの獲得数・欠失数を表示する対話型HTMLを作成します。

## インストール

```bash
pip install git+https://github.com/kazumaxneo/PhyloGain.git
```

主要コマンドは `phylogain` です。従来の `species-map` も互換エイリアスとして引き続き利用できます。

## 基本的な実行方法

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --output species_innovation_map

species-map serve species_innovation_map
```

branchをクリックすると、そのbranchでGain/Lossしたorthogroupが表示されます。さらにorthogroupをクリックすると、`Orthogroups.tsv`に含まれる遺伝子IDを確認できます。

## phenotypeを追加する

`phenotypes.tsv`：

```tsv
species_id	nitrogen_fixation
species_A	+
species_B	+
species_C	-
species_D	?
```

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --phenotypes phenotypes.tsv \
  --phenotype nitrogen_fixation \
  --output nitrogen_fixation_map

species-map serve nitrogen_fixation_map
```

phenotypeの獲得branchと同じbranchでGainしたgene familyが`candidate_genes.tsv`に出力されます。

## 分類群レベルのpangenome解析

`--taxon-rank`を指定すると、GTDB taxonomyに従ってゲノム単位のpresence/absenceをspecies、genus、family、order単位のoccupancyへ集約してからGain/Lossを推定します。未指定時と`genome`指定時は従来どおりのゲノム単位解析です。

```bash
phylogain build \
  --orthofinder Results_Jul02 \
  --gtdb-taxonomy gtdbtk.bac120.summary.tsv \
  --taxon-rank genus \
  --present-threshold 0.90 \
  --absent-threshold 0.10 \
  --min-genomes-per-taxon 3 \
  --output phylogain_genus
```

absence閾値とpresence閾値の中間は`polymorphic`として保持し、Sankoff解析には曖昧状態`{0,1}`として渡します。入力系統樹上で単系統となる分類群だけをcollapseして解析し、非単系統群は警告とともに除外します。genus以上の広い比較ではOrthoFinderを推奨します。

## 主な出力

- `index.html`：対話型Gene Gain/Loss Viewer
- `branches.tsv`：branchごとのGain/Loss数
- `gene_gain_loss.tsv`：各gene familyの推定Gain/Loss
- `phenotype_gain_loss.tsv`：phenotypeの推定変化
- `candidate_genes.tsv`：phenotypeと同時にGainした候補
- `species_map.sqlite`：HTML表示用データベース
- `taxon_occupancy.tsv`：分類群×gene familyのoccupancyと観測・推定状態（`--taxon-rank`指定時）
- `taxon_tree.nwk`：単系統群をcollapseした分類群レベル系統樹（`--taxon-rank`指定時）

## 注意

現在のバージョンはpresence/absenceに対するSankoff parsimonyを使用します。不完全なMAGでは、未検出遺伝子が偽のLossとして推定される可能性があります。

## GTDB分類で系統樹を折り畳む

GTDB-Tkのsummary TSV、または`species_id`と`gtdb_taxonomy`列を持つTSVを指定すると、Phylum、Class、Order、Family、Genus単位で単系統クレードを折り畳めます。

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --gtdb-taxonomy gtdbtk.bac120.summary.tsv \
  --output gtdb_collapsible_map
```

折り畳み名をクリックすると、その分類群だけ展開されます。同じ分類群名が複数表示される場合、その分類群は入力系統樹上で単系統ではありません。折り畳み表示のGain/Lossは、そのクレードへ入るbranchの値であり、内部の全イベントの合計ではありません。

OrthoFinder treeの代わりに、tip名をOrthoFinder species IDへ揃えたrooted Newickを指定することもできます。

```bash
species-map build \
  --orthofinder Results_Jul02 \
  --species-tree gtdb_pruned_rooted.nwk \
  --gtdb-taxonomy gtdbtk.bac120.summary.tsv \
  --output gtdb_tree_map
```
