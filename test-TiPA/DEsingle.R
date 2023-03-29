#!/usr/bin/env Rscript

# Sergio Alías, 20230321
# Last modified 20230329

# Script for runnig DE analysis on single-cell count matrices using R package DEsingle

## Imports

library(DEsingle)
library(data.table)
library(BiocParallel)

## Args

args = commandArgs(trailingOnly=TRUE)

## Test data

# data(TestData)
# results <- DEsingle(counts = counts, group = group)
# results.classified <- DEtype(results = results, threshold = 0.05)
# results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]

## Functions

writeToLog <-
function(msg,
         log = logfile
         )
{
  write(msg, log, append = TRUE)
}


## Main script

logfile <- paste0("test", ".log")
file.create(logfile)

tissue <- args[1]
mclus <- args[2]

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

for (i in c(0:mclus)){

  message(paste0("Running for cluster ", i, "...\n"))
  
  group <- numeric(length = nrow(ann))
  group[which(ann[,3] == i)] <- 1
  group[which(ann[,3] != i)] <- 2
  group <- as.factor(group)

  ### Parallelization

  param <- SnowParam(workers = 40, type = "SOCK", progressbar = TRUE)
  register(param)
  before <- proc.time()
  #####
  results <- DEsingle(counts = counts, group = group, parallel = TRUE, BPPARAM = param)
  #####
  after <- proc.time()
  time_paral <- after - before 

  writeToLog(paste0("\nWith parallelization (40 cores, cluster ", i, "):\n"))
  writeToLog(time_paral)



  results.classified <- DEtype(results = results, threshold = 0.05)
  results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]

  # we are interested in DEg and DEa

  results.DEga <- results.sig[results.sig$Type != "DEs", ]
  DEgenes <- rownames(results.DEga)

  saveRDS(DEgenes,
          file = paste0("DEresults/", tissue, "/DEgenes-vs-tissue-", tissue, "-", cluster, ".RDS")
          )

}