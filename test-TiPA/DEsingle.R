#!/usr/bin/env Rscript

# Sergio Alías, 20230321
# Last modified 20230323

# Script for runnig DE analysis on single-cell count matrices using R package DEsingle

## Imports

library(DEsingle)
library(data.table)

## Main script

data(TestData)
results <- DEsingle(counts = counts, group = group)
results.classified <- DEtype(results = results, threshold = 0.05)
results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]


data <- fread("TFM/rna_single_cell_read_count.tsv")
ann <- data[, c(1:3)]
counts <- data[, -c(1:3)]
rm(data)
c(1:nrow(ann))[ann[,1] == "Adipose" & ann[,3] == 0]