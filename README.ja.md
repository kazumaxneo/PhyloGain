# Species Innovation Map

OrthoFinderの結果から、species treeの各branchにgene-familyの獲得数・欠失数を表示する対話型HTMLを作成します。

## インストール

```bash
pip install git+https://github.com/kazumaxneo/species-innovation-map.git
```

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

## 主な出力

- `index.html`：対話型Species Innovation Map
- `branches.tsv`：branchごとのGain/Loss数
- `gene_gain_loss.tsv`：各gene familyの推定Gain/Loss
- `phenotype_gain_loss.tsv`：phenotypeの推定変化
- `candidate_genes.tsv`：phenotypeと同時にGainした候補
- `species_map.sqlite`：HTML表示用データベース

## 注意

現在のバージョンはpresence/absenceに対するSankoff parsimonyを使用します。不完全なMAGでは、未検出遺伝子が偽のLossとして推定される可能性があります。
