#!/usr/bin/env Rscript

# Sergio Alías, 20230321
# Last modified 20230324

# Script for runnig DE analysis on single-cell count matrices using R package DEsingle

## Imports

library(DEsingle)
library(data.table)

## Test data

# data(TestData)
# results <- DEsingle(counts = counts, group = group)
# results.classified <- DEtype(results = results, threshold = 0.05)
# results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]

## Main script

tissue <- "Liver"

data <- fread("../rna_single_cell_read_count.tsv") # Read HPA dataset
ann <- data[, c(1:3)]
counts <- data[, -c(1:3)]
rm(data)

# c(1:nrow(ann))[ann[,1] == "Liver" & ann[,3] == 0] # For checking

### Filter the tissue (e.g. Liver)

counts <- counts[c(1:nrow(ann))[ann[,1] == tissue],]
ann <- ann[c(1:nrow(ann))[ann[,1] == tissue],]


### Preparing the input for DEsingle

# counts is a count matrix or an sce object
# The rows of the matrix are genes and columns are cells

counts <- t(counts)

# group is a vector of factor which specifies the two groups in the matrix to be compared, corresponding to the columns in counts

group <- numeric(length = nrow(ann))
group[which(ann[,3] == 0)] <- 1
group[which(ann[,3] != 0)] <- 2
group <- as.factor(group)


results <- DEsingle(counts = counts, group = group)
results.classified <- DEtype(results = results, threshold = 0.05)
results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]