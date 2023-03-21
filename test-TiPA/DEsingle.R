#!/usr/bin/env Rscript

# Sergio Alías, 20230321
# Last modified 20230321

# Script for runnig DE analysis on single-cell count matrices using R package DEsingle

## Imports

library(DEsingle)

## Main script

data(TestData)
results <- DEsingle(counts = counts, group = group)
results.classified <- DEtype(results = results, threshold = 0.05)
