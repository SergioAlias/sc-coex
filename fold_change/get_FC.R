#!/usr/bin/env Rscript

# Sergio Alías, 20230417
# Last modified 20230419

# Script for obtaining the FC of genes on the pooled RNAseq HPA dataset

### Libs ###

library(data.table)

### Main script ###

RNAseq <- fread("../rna_single_cell_type_tissue.tsv") # Read data

colnames(RNAseq)[c(2, 5)] <- c("Gene_name", "Cell_type") # Change some colnames with spaces

RNAseq[, sum_pTPM := sum(pTPM), by = .(Gene, Tissue)] # Create a new column with the sum of pTPM for each gene and tissue

RNAseq[, n_clusters := uniqueN(Cluster), by = .(Gene, Tissue)] # Create a new column with the count of clusters for each gene and tissue

RNAseq[, log2FC := log2(pTPM / ((sum_pTPM - pTPM) / (n_clusters - 1)))] # Calculate log2 fold-change for each row

# RNAseq[, c("sum_pTPM", "n_clusters") := NULL] # Remove intermediate columns





# Compute the total pTPM and number of clusters for each Gene, Tissue, and Cell_type combination
RNAseq[, `:=` (sum_pTPM_per_ct = sum(pTPM), n_clusters_per_ct = .N), by = .(Gene, Tissue, Cell_type)]

# Create a new column with the sum of pTPM for each Gene and Tissue combination, excluding the current row's Cell_type value
RNAseq[, sum_pTPM_single_tissue := sum_pTPM - sum_pTPM_per_ct]

# Create a new column with the count of clusters for each Gene and Tissue combination, excluding the current row's Cell_type value
RNAseq[, n_clusters_single_tissue := n_clusters - n_clusters_per_ct]

# Create a new column with the log2 fold-change for each row, relative to other clusters with different Cell_type values
RNAseq[, log2FC_single_tissue := log2(pTPM / (sum_pTPM_single_tissue / n_clusters_single_tissue))]

# Remove the temporary columns
RNAseq[, c("sum_pTPM", "n_clusters", "sum_pTPM_per_ct", "n_clusters_per_ct", "sum_pTPM_single_tissue", "n_clusters_single_tissue") := NULL]

fwrite(RNAseq,
       file = "computed_FC.tsv",
       sep = "\t")
