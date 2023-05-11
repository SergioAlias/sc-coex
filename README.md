# sc-net-analysis

Code used for my Bioinformatics and Computational Biology MSc thesis: Inferencia de redes moleculares y caracterización de relaciones genotipo-fenotipo a partir de datos de single-cell RNA-seq. You can ask me for a copy: sergioaliaseg[at]gmail[dot]com (please note that the document is written in Spanish).

### Abstract

Clinical signs are useful for describing the spectrum of human pathologies. Integrating these phenotypes with single-cell RNA-seq data allows the identification of potential relationships between phenotypes and the specific cell types causing them. In this Master's thesis, single-cell transcriptomic data were used to infer these relationships.Differential expression and co-expression of genes associated with abnormal phenotypes were computed as possible measures to identify the cell types associated with a phenotype. The results obtained by these metrics were compared with the set of phenotype-cell type associations described in the literature. The results showed that co-expression, not so much differential expression, emerges as a measure that allows, on the one hand, to identify significant relationships already present in the literature; and on the other hand, to point out cell types potentially associated with abnormal phenotypes that have not yet been extensively described in the literature.

### Scripts

The following list is exahustive and covers every script on the repo. Feel free to `ctrl+F` to the specific scripts you are interested in.

---

- `HPA/extract-dataset.sh`: Bash script for dividing the HPA dataset according to the tissues. It is used from the command line as follows:

```bash
./extract-dataset.sh tissue last-cluster
```

Being `tissue` the tissue name (e.g. liver) and last-cluster the number of the las cluster of the tissue (17). The output is a TSV file with the HPA dataset for that specific tissue.

---

- `HPA/results/COTAN-cluster-script.R`: R script for performing the Gene Pair Analysis step of COTAN. It is used inside `run_cluster_script.sh` as follows:

```bash
Rscript COTAN-cluster-script.R tissue-cluster
```

Being tissue-cluster the combination of the tissue name and the cluster name (e.g. liver-0). It reads the RDS file with the `scCOTAN` object already preprocessed and outputs another RDS file with the co-expression matrix.

---

- `HPA/results/run_cluster_script.sh`: Bash script for running `COTAN-cluster-script.R` for each cluster of a given tissue. It is used from the command line as follows:

```bash
./run_cluster_script.sh tissue last-cluster
```

Being `tissue` the tissue name (e.g. liver) and last-cluster the number of the las cluster of the tissue (17). The output is an RDS file with a co-expression matrix for each loop iteration.

---

- `annotations/annotation-files.py`: Python script for integrating all the annotations generated in the MSc Thesis. It is used form the command line as follows:

```bash
python3 annotation-files.py
```

It takes Uberon child terms, Uberon-HPO relationships, HPO-genes relationships and the HPO OBO file (for HPO children terms) and generates two TSV files: the first one has the structure tissue **`\t`** HPO code **`\t`** HPO name **`\t`** num genes; and the second one has the structure HPO code **`\t`** gene code **`\t`** gene name.

---