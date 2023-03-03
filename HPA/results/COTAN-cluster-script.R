#!/usr/bin/env Rscript
# Sergio Alías, 20220405
# Last modified 20220614

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) stop("One argument must be supplied (tissue)", call.=FALSE)

suppressMessages(library(COTAN))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(ggrepel))
# library(factoextra)
suppressMessages(library(Rtsne))
suppressMessages(library(utils))
suppressMessages(library(plotly))
suppressMessages(library(tidyverse))
suppressMessages(library(htmlwidgets))
suppressMessages(library(MASS))
suppressMessages(library(dendextend))

# Main script

tissue <- args[1]

paste('Executing for tissue ', tissue, '...', sep = '')

obj <- readRDS(paste('HPA.', tissue, '.input.cotan.RDS', sep = ''))

obj <- cotan_analysis(obj,cores = 2)

obj <- get.coex(obj)

coex <- extract.coex(object = obj)

saveRDS(coex, file = paste('HPA.', tissue, '.matrix.cotan.RDS', sep = ''))

paste('COTAN object saved --> HPA.', tissue, '.matrix.cotan.RDS', sep = '')

print("Finished!")
