#!/usr/bin/env Rscript

# Sergio Alías, 20230421
# Last modified 20230421

# Script for scatterplot genes of specific examples we find interesting


### Libs ###

library(data.table)
library(jsonlite)
library(ggplot2)

### Main script ###

measure <- "fcst" # coex, fc or fcst
tissue <- "adipose-tissue"
hpo <- "HP0100578"
hpo_semicolon <- "HP:0100578"
clus_a <- "13"
cond_a <- "high"
clus_b <- "10"
cond_b <- "low"
x_label <- paste0("c-", clus_a, "-", cond_a, "-", measure)
y_label <- paste0("c-", clus_b, "-", cond_b, "-", measure)

if (measure != "coex"){
  data <- fread("../../test-TiPA/computed_FC.tsv")
  data <- data[!is.na(data$log2FC)]
  data <- data[data$log2FC != -Inf]
  data <- data[data$log2FC != Inf]
  data <- data[!is.na(data$log2FC_single_tissue)]
  data <- data[data$log2FC_single_tissue != -Inf]
  data <- data[data$log2FC_single_tissue != Inf]
  data[Tissue == "adipose tissue", Tissue := gsub("adipose tissue", "adipose-tissue", Tissue)]
  data[Tissue == "bone marrow", Tissue := gsub("bone marrow", "bone-marrow", Tissue)]
  data[Tissue == "heart muscle", Tissue := gsub("heart muscle", "heart", Tissue)]
  data[Tissue == "lymph node", Tissue := gsub("lymph node", "lymph-node", Tissue)]
  data[Tissue == "prostate", Tissue := gsub("prostate", "prostate-gland", Tissue)]
  data[Tissue == "skeletal muscle", Tissue := gsub("skeletal muscle", "skeletal-muscle-organ", Tissue)]
  data[Tissue == "skin", Tissue := gsub("skin", "skin-of-body", Tissue)]
  data[Tissue == "small intestine", Tissue := gsub("small intestine", "small-intestine", Tissue)]
  data[Tissue == "pbmc", Tissue := gsub("pbmc", "blood", Tissue)]
}


hpo_genes <- fromJSON("../HPO/hpo_to_genes_HPA.json")[[hpo_semicolon]]$EnsemblID


data <- data[Tissue == tissue]
data <-data[Gene %in% hpo_genes]
data <-data[Cluster %in% c(paste0("c-", clus_a), paste0("c-", clus_b))]

pivoted_data <- dcast(data, Gene ~ Cluster, value.var = "log2FC_single_tissue")

clean_data <- na.omit(pivoted_data)

colnames(clean_data)[c(2, 3)] <- c(paste0("c", clus_b), paste0("c", clus_a))


pdf("scatter_test.pdf")

ggplot(clean_data, aes(x = c13, y = c10)) +
  geom_point() +
  labs(x = x_label, y = y_label)

dev.off()

