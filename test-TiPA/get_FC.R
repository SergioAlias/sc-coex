#!/usr/bin/env Rscript

# Sergio Alías, 20230417
# Last modified 20230417

# Script for obtaining the FC of genes on the pooled RNAseq HPA dataset

### Libs ###

library(data.table)

### Main script ###

RNAseq <- fread("../../rna_single_cell_type_tissue.tsv") # Read data

RNAseq[, sum_pTPM := sum(pTPM), by = .(Gene, Tissue)] # Create a new column with the sum of pTPM for each gene and tissue

RNAseq[, n_clusters := uniqueN(Cluster), by = .(Gene, Tissue)] # Create a new column with the count of clusters for each gene and tissue

RNAseq[, log2FC := log2(pTPM / ((sum_pTPM - pTPM) / (n_clusters - 1)))] # Calculate log2 fold-change for each row

RNAseq[, c("sum_pTPM", "n_clusters") := NULL] # Remove intermediate columns
