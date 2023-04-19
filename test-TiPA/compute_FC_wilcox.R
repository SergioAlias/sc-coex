#!/usr/bin/env Rscript

# Sergio Alías, 20230419
# Last modified 20230419

# Script for obtaining gene FC distributions associated to HPOs and making wilcox test against random distributions

# We need:
#    - FC results -> computed_FC.tsv
#    - Tissue-HPO -> TFM/annotations/filtered_tis_hpo_ngenes_HPA.tsv
#    - HPO-ENSEMBL -> TFM/annotations/HPO/hpo_to_genes_HPA.json


### Libs ###

library(data.table)
library(jsonlite)

### Main script ###

# Input 

FC <- fread("computed_FC.tsv")
tis_hpo <- fread("../annotations/filtered_tis_hpo_ngenes_HPA.tsv")
hpo_genes <- fromJSON("../annotations/HPO/hpo_to_genes_HPA.json")

# Input filtering

colnames(tis_hpo)[1] <- "Tissue"

FC <- FC[!is.na(FC$log2FC)]
FC <- FC[FC$log2FC != -Inf]
FC <- FC[FC$log2FC != Inf]
FC <- FC[!is.na(FC$log2FC_single_tissue)]
FC <- FC[FC$log2FC_single_tissue != -Inf]
FC <- FC[FC$log2FC_single_tissue != Inf]
FC[Tissue == "adipose tissue", Tissue := gsub("adipose tissue", "adipose-tissue", Tissue)]
FC[Tissue == "bone marrow", Tissue := gsub("bone marrow", "bone-marrow", Tissue)]
FC[Tissue == "heart muscle", Tissue := gsub("heart muscle", "heart", Tissue)]
FC[Tissue == "lymph node", Tissue := gsub("lymph node", "lymph-node", Tissue)]
FC[Tissue == "prostate", Tissue := gsub("prostate", "prostate-gland", Tissue)]
FC[Tissue == "skeletal muscle", Tissue := gsub("skeletal muscle", "skeletal-muscle-organ", Tissue)]
FC[Tissue == "skin", Tissue := gsub("skin", "skin-of-body", Tissue)]
FC[Tissue == "small intestine", Tissue := gsub("small intestine", "small-intestine", Tissue)]
FC[Tissue == "pbmc", Tissue := gsub("pbmc", "blood", Tissue)]


# Wilcoxon tests

seed <- 10
rnd_samples <- 1000
tis <- unique(FC$Tissue)

for (i in seq_along(tis)){
  tissue <- tis[i]
  tis_FC <- FC[Tissue == tissue]
  tis_hpo_reduced <- tis_hpo[Tissue == tissue]
  hpocodes <- tis_hpo_reduced$hpo_code
  hponames <- tis_hpo_reduced$hpo_name
  message("Running for tissue: ", tissue, " (", length(hponames), "HPO terms)")
  for (h in seq_along(hpocodes)){
    hp_code <- hpocodes[h]
    hp_name <- hponames[h]
    genes <- hpo_genes[[hp_code]]$EnsemblID
    message("-> ", hp_code, " ", hp_name)
    clus <- unique(tis_FC$Cluster)
    for (j in seq_along(clus)){
      cluster <- clus[j]
      cluster_name <- unique(tis_FC[Cluster == cluster, Cell_type])
      message("    - Cluster ", cluster, ": ", cluster_name, "...")
      tis_clus_FC <- tis_FC[Cluster == cluster]
      fc_hpo_genes <- tis_clus_FC[Gene %in% genes]$log2FC
      fcst_hpo_genes <- tis_clus_FC[Gene %in% genes]$log2FC_single_tissue
      rnd_sample_size <- length(fc_hpo_genes)
      # creating objects for storing distributions
      fc_df <- data.frame(Subset = rep('HPOgenes', length(fc_hpo_genes)), FC = fc_hpo_genes)
      fcst_df <- data.frame(Subset = rep('HPOgenes', length(fcst_hpo_genes)), FC = fcst_hpo_genes)
      df_list_fc <- list('HPOgenes' = fc_df)
      df_list_fcst <- list('HPOgenes' = fcst_df)
      rm(fc_df)
      rm(fcst_df)
      # sampling 1000 random FC distributions
      set.seed(seed)
      for (r in c(1:rnd_samples)) {
        rnd_genes <- sample(1:nrow(tis_clus_FC), rnd_sample_size)
        rnd_fc <- tis_clus_FC[rnd_genes, log2FC]
        rnd_fcst <- tis_clus_FC[rnd_genes, log2FC_single_tissue]
        cname <- paste0('random_', r)
        to_add_fc <- data.frame(Subset = rep(cname, length(rnd_fc)), FC = rnd_fc)
        to_add_fcst <- data.frame(Subset = rep(cname, length(rnd_fcst)), FC = rnd_fcst)
        df_list_fc[[cname]] <- to_add_fc
        df_list_fcst[[cname]] <- to_add_fcst
      }
      # wilcoxon tests
      
    }
  }
}