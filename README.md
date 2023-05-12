# sc-net-analysis

Code used for my Bioinformatics and Computational Biology MSc Thesis: Inferencia de redes moleculares y caracterización de relaciones genotipo-fenotipo a partir de datos de single-cell RNA-seq. You can ask me for a copy: sergioaliaseg[at]gmail[dot]com (please note that the document is written in Spanish).

### Abstract

Clinical signs are useful for describing the spectrum of human pathologies. Integrating these phenotypes with single-cell RNA-seq data allows the identification of potential relationships between phenotypes and the specific cell types causing them. In this Master's thesis, single-cell transcriptomic data were used to infer these relationships. Differential expression and co-expression of genes associated with abnormal phenotypes were computed as possible measures to identify the cell types associated with a phenotype. The results obtained by these metrics were compared with the set of phenotype-cell type associations described in the literature. The results showed that co-expression, not so much differential expression, emerges as a measure that allows, on the one hand, to identify significant relationships already present in the literature; and on the other hand, to point out cell types potentially associated with abnormal phenotypes that have not yet been extensively described in the literature.

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

It takes Uberon child terms, Uberon-HPO relationships, HPO-genes relationships and the HPO OBO file (for HPO children terms) and generates two TSV files: the first one has the structure tissue **\t** HPO code **\t** HPO name **\t** num genes; and the second one has the structure HPO code **\t** gene code **\t** gene name.

---

- `annotations/specific-annotation-file.py`: Python script for filtering the first output file of `annotation-files.py` to keep only HPO terms with 20 or more associated genes that do not have any children term with 20 or more associated genes. It is used from the command line as follows:

```bash
python3 specific-annotation-file.py
```

It takes output file 1 from `annotation-files.py` and the HPO OBO file, and generates a filtered TSV file with the same format.

---

- `annotations/HPO/ENSEMBL-to-Entrez-id.py`: Python script for obtaining Entrez ID- Ensembl ID relationships for HPA genes. It is used from the comand line as follows:

```bash
python3 ENSEMBL-to-Entrez-id.py
```

It takes HPA genes Ensembl IDs and Ensembl-Entrez relationships retrieved from BioMart, and generates two JSON files: the first file with the Ensembl IDs as keys and the Entrez IDs as the content of the entries, and the second file the other way around.

---

- `annotations/HPO/HPO-UBERON.py`: Python script for generating HPO-Uberon relationships. It is used from the command line as follows:

```bash
python3 HPO-UBERON.py
```

It takes the HPO OWL file and generates two JSON files: the first file with the HPO terms as keys and the Uberon terms as the content of the entries, and the second file the other way around.

---

- `annotations/HPO/entrez-HPO.py`: Python script for generating HPO-genes relationships. It is used from the command line as follows:

```bash
python3 entrez-HPO.py
```

It takes Ensembl-Entrez relationships and HPO annotation file and generates a JSON file with the HPO terms as keys and the gene Entrez ID, Ensembl ID and regular name as the content of the entries.

---

- `annotations/HPO/hpo-id-to-name.py`: Python script for having an easy access to HPO terms names. It is used from the command line as follows:

```bash
python3 hpo-is-to-name.py
```

It takes the HPO OBO file and generates the file hpo-id-to-name, with the format HPO code **\t** HPO name.

---

- `annotations/UBERON/children-terms.py`: Python script for obtaining children terms of Uberon terms (or any ontology). It is used inside `run-children-terms.sh` as follows:

```bash
python3 children-terms.py OBO tissue
```

Being `OBO` the ontology OBO file (Uberon in this case) and `tissue` the term whose childs we want to retrieve (in this case the name of an Uberon term). The output is printed direclty to the console and it has the structure UBERON term **\t** term name **\t** UBERON subterm **\t** subterm name **\t** relationship.

---

- `annotations/UBERON/run-children-terms.sh`: Bash script for running `children-terms.py` for all tissues in HPA and redirecting the output to actual files. It is used from the command line as follows:

```bash
./run-children-terms.sh OBO tissue-list outdir
```

Being `OBO` the Uberon OBO file, `tissue-list` the list of tissue names and `outdir` the output directory. It generates a txt file with the output of `run-children-terms.sh` for each loop iteration.

---

- `annotations/comention/get_comentions.R`: R script for annotating K-S test results with co-mention information. It is used from the command line as follows:

```bash
./get_comentions.R
```

We need to comment/uncomment lines depending of if we are using fold-change results or co-expression results. It takes K-S test output TSV files and co-mention file, and generates a new K-S test file with additional columns for co-mention information.

---

- `annotations/comention/get_conf_matrix.R`: R script for generating the confusion matrices. It is used from the command line as follows:

```bash
./get_conf_matrix.R
```

We need to change the `metric` and `pval_thr` variables according to our needs. It generates TSV files with the confusion matrices for each condition.

---

- `annotations/metrics/check_TP.R`: R script for checking for each celltype-HPO pair if those pairs are the same for high TP and low TP. It is used from the command line as follows:

```bash
./check_TP.R
```

It takes the K-S test result files with the co-mention information added. It generates a TSV file with the shared TP pairs.

---

- `fold_change/compute_FC_ks.R`: R script for obtaining gene FC distributions associated to HPOs and making K-S tests. It is used from the command line as follows:

```bash
./compute_FC_ks.R
```

It takes fold-change results, filtered annotations and HPO-genes relationships. The output is a TSV file for each tissue and fold-change type that contains K-S p-values for extreme. high and low values tests.

---

- `fold_change/get_FC.R`: R script for obtaining the logFC and logFC-st of genes on the pooled RNA-seq HPA dataset. It is used from the command line as follows:

```bash
./get_FC.R
```

It takes the pooled RNA-seq HPA dataset and it generates a TSV file with the logFC and logFC-st values.

---

- `automatic-COTAN-script.R`: R script for performing the pre-processing step of COTAN. It is used from the command line as follows:

```bash
Rscript automatic-COTAN-script.R
```

Once executed, the script will ask for some input (tissue name, number of the last cluster and dataset used). It takes the HPA dataset for the tissue and cluster selected and generates an RDS file with the pre-processed `scCOTAN` R object.

---

- `coex-to-tsv.R`: R script for converting RDS files with co-expression matrices to TSV files. It is used inside `workflow_coex.sh` as follows:

```bash
Rscript coex-to-tsv.R --tissue tissue --n_cluster n_cluster --dataset dataset
```

Being `tissue` the tissue name (e.g. liver), `n_cluster` the number of clusters (e.g. 17) and `dataset` the dataset name (HPA). It takes an RDS file with a co-expression matrix and it outputs a TSV file containing that matrix.

---

- `plot_coex.R`: R script for performing Kolmogórov-Smirnov tests of HPO-related genes co-expression in a tissue (despite the "plot" in the script name, this script does not plot anything since its name is due to historical reasons). It is used inside `workflow_coex.sh` as follows:

```bash
Rscript plot_coex.R --tissue tissue --max_cluster max_cluster --hpo hpo --dataset dataset
```

Being `tissue` the tissue name (e.g. liver), `max_cluster` the number of clusters (e.g. 17), `hpo` the HPO code without the colon (e.g. HP0006561) and `dataset` the dataset name (HPA). It takes a file with the cluster annotation for the tissue, HPO code-HPO name relationships, a list of genes related to the HPO term and the co-expression matrix and the output is a TSV file with the p-values of all K-S tests performed.

---

- `plot_ks.R`: R script for plotting distributions and the corresponding K-S test p-value. It is used from the command line as follows:

```bash
./plot_ks.R
```
We need to change the `metric`, `test`, `hpo_term`, `tissue` and `cluster` variables according to our needs. It takes the cluster annotation file for the tissue, K-S results for the tissue and measure used (differencial expression or COEX), a list of genes related to the HPO term and the co-expression matrix. It generates a PDF file with the plot.

---

- `plot_percent_HPO.R`: R script for plotting the percentage of cell-type-related co-mentions per HPO term. It is used from the command line as follows:

```bash
./plot_percent_HPO.R
```

It takes the co-mention file and outputs a PDF file with the plot.

---

- `subset-genes.py`: Python script for retrieving all genes related to an HPO term in a cluster. It is used inside `workflow_coex.sh` as follows:

```bash
python3 subset-genes.py hpo tissue n_clusters dataset
```

Being `hpo` the HPO term without the colon (e.g. HP0006561), `tissue` the tissue name (e.g. liver), `n_clusters` the number of clusters (e.g. 17) and `dataset` the dataset name (HPA). It takes the HPO-genes relationships, the Entrez-Ensembl relationships and the TSV file with the co-expression matrix. It generates a text file with all genes associated with the HPO term that are in the cluster.

---

- `workflow_coex.sh`: Bash script for performing some steps of the COEX workflow (from RDS matrices to K-S tests results). It is used from the command line as follows:

```bash
./workflow_coex.sh tissue
```

Being `tissue` the tissue name (e.g. liver). There are other optional parameters we normally don't need to use (e.g. argument 2 can be a text file with HPO names, defaults to all; argument 3 can be another dataset, defaults to HPA). The output is the output of scripts `coex-to-tsv.R`, `subset-genes.py` and `plot_coex.R`.