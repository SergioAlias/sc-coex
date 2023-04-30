#!/usr/bin/env Rscript

# Sergio Alías, 20230430
# Last modified 20230430

# Script for checking for each celltype-HPO pair if those pairs are the same for high TP and low TP (only coex results)

# We need:
#    - K-S results with comention -> TFM/results_w_comention/[tissue]/[tissue]_HPA_coex_results_w_comention.tsv // [tissue]_HPA_05pval_coex_results_w_comention.tsv


### Packages

library(data.table)


### Input

urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"
urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"

tissues <- c("blood", "liver", "lung", "pancreas", "spleen", "stomach", "testis")

for (t in seq_along(tissues)){
  
  tis <- tissues[t]
  res_w_com <- fread(file.path(urales_home,
                               "TFM/results_w_comention",
                               tis,
                               paste0(tis, "_HPA_coex_results_w_comention.tsv")))
  
}