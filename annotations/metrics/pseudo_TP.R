#!/usr/bin/env Rscript

# Sergio Alías, 20230608
# Last modified 20230608

# Script for counting TP for 0,05 p-val threshold but not for 10^-3

# We need:
#    - K-S results with comention (0,05 p-val threshold) -> [tissue]_HPA_05pval_coex_results_w_comention.tsv

### Packages

library(data.table)

### Input

# urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"
urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"

tissues <- c("blood", "liver", "lung", "pancreas", "spleen", "stomach", "testis", "skin")

results <- data.table(tissue = character(),
                      TP05_high = numeric(),
                      pseudo_TP05_high = numeric(),
                      perc_TP05_high = numeric(),
                      TP05_low = numeric(),
                      pseudo_TP05_low = numeric(),
                      perc_TP05_low = numeric(),
                      TP05_extreme = numeric(),
                      pseudo_TP05_extreme = numeric(),
                      perc_TP05_extreme = numeric(),
                      TPf05_high = numeric(),
                      pseudo_TPf05_high = numeric(),
                      perc_TPf05_high = numeric(),
                      TPf05_low = numeric(),
                      pseudo_TPf05_low = numeric(),
                      perc_TPf05_low = numeric(),
                      TPf05_extreme = numeric(),
                      pseudo_TPf05_extreme = numeric(),
                      perc_TPf05_extreme = numeric()
)

pval_thr <- 0.001
pval_thr_2 <- 0.05


for (t in seq_along(tissues)){
  
  res_05pval <- fread(file.path(urales_home,
                                "TFM/results_w_comention",
                                tis,
                                paste0(tis, "_HPA_05pval_coex_results_w_comention.tsv")))
  fil_high_05pval <- res_05pval[, .SD[which.min(high_wil_pval)], by = .(hpo, annotation)]
  fil_low_05pval <- res_05pval[, .SD[which.min(low_wil_pval)], by = .(hpo, annotation)]
  fil_extreme_05pval <- res_05pval[, .SD[which.min(wilcoxon.pval)], by = .(hpo, annotation)]
  
  ###
  
  high_05pval <- res_05pval[high_wil_pval <= pval_thr & coment.pval > pval_thr & coment.pval <= pval_thr_2]
  low_05pval <- res_05pval[low_wil_pval <= pval_thr & coment.pval > pval_thr & coment.pval <= pval_thr_2]
  extreme_05pval <- res_05pval[wilcoxon.pval <= pval_thr & coment.pval > pval_thr & coment.pval <= pval_thr_2]
  
  pseudo_high_05pval <- res_05pval[high_wil_pval <= pval_thr & coment.pval > pval_thr & coment.pval <= pval_thr_2]
  pseudo_low_05pval <- res_05pval[low_wil_pval <= pval_thr & coment.pval > pval_thr & coment.pval <= pval_thr_2]
  pseudo_extreme_05pval <- res_05pval[wilcoxon.pval <= pval_thr & coment.pval > pval_thr & coment.pval <= pval_thr_2]
  
}