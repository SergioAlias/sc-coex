#!/usr/bin/env Rscript

# Sergio Alías, 20230411
# Last modified 20230411

# Script for calculating genes fold change on single-cell count matrices using Seurat

# Imports

library(data.table)
library(Seurat)

## Args

args = commandArgs(trailingOnly=TRUE)

# Main script

tissue <- args[1]
mclus <- args[2]

data <- fread("../rna_single_cell_read_count.tsv") # Read HPA dataset
ann <- data[, c(1:3)]
counts <- data[, -c(1:3)]
rm(data)

### Filter the tissue (e.g. Liver)

counts <- counts[c(1:nrow(ann))[ann[,1] == tissue],]
ann <- ann[c(1:nrow(ann))[ann[,1] == tissue],]

### Prepare count matrix for Seurat

raw.counts <- t(counts)
rm(counts)
colnames(raw.counts) <- ann[, Cell]

### Create Seurat object

seu <- CreateSeuratObject(counts = raw.counts, project = "testing")

### Create vectors of cell groups [MAKE A FOR LOOP]

cells.1 <- which(ann[,3] == i)
cells.2 <- which(ann[,3] != i)

### Calculate fold change

fc <- FoldChange(object = seu, ident.1 = cells.1)