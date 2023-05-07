#!/usr/bin/env Rscript

# Sergio Alías, 20230430
# Last modified 20230501

# Script for checking for each celltype-HPO pair if those pairs are the same for high TP and low TP (only coex results)

# We need:
#    - K-S results with comention -> TFM/results_w_comention/[tissue]/[tissue]_HPA_coex_results_w_comention.tsv // [tissue]_HPA_05pval_coex_results_w_comention.tsv


### Packages

library(data.table)


### Input

urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"
#urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"

tissues <- c("blood", "liver", "lung", "pancreas", "spleen", "stomach", "testis", "skin")

results <- data.table(tissue = character(),
                      TP_high = numeric(),
                      TP_low = numeric(),
                      common = numeric(),
                      TPf_high = numeric(),
                      TPf_low = numeric(),
                      common_f = numeric(),
                      TP05_high = numeric(),
                      TP05_low = numeric(),
                      common_05 = numeric(),
                      TPf05_high = numeric(),
                      TPf05_low = numeric(),
                      common_f_05 = numeric()
                      )

pval_thr <- 0.001
pval_thr_2 <- 0.05

for (t in seq_along(tissues)){
  
  tis <- tissues[t]
  res_w_com <- fread(file.path(urales_home,
                               "TFM/results_w_comention",
                               tis,
                               paste0(tis, "_HPA_coex_results_w_comention.tsv")))
  fil_high <- res_w_com[, .SD[which.min(high_wil_pval)], by = .(hpo, annotation)]
  fil_low <- res_w_com[, .SD[which.min(low_wil_pval)], by = .(hpo, annotation)]
  
  res_05pval <- fread(file.path(urales_home,
                               "TFM/results_w_comention",
                               tis,
                               paste0(tis, "_HPA_05pval_coex_results_w_comention.tsv")))
  fil_high_05pval <- res_05pval[, .SD[which.min(high_wil_pval)], by = .(hpo, annotation)]
  fil_low_05pval <- res_05pval[, .SD[which.min(low_wil_pval)], by = .(hpo, annotation)]
  
  ###
  
  high <- res_w_com[high_wil_pval <= pval_thr & coment.pval <= pval_thr]
  low <- res_w_com[low_wil_pval <= pval_thr & coment.pval <= pval_thr]
  
  high_05pval <- res_05pval[high_wil_pval <= pval_thr_2 & coment.pval <= pval_thr_2]
  low_05pval <- res_05pval[low_wil_pval <= pval_thr_2 & coment.pval <= pval_thr_2]
  
  high_f <- fil_high[high_wil_pval <= pval_thr & coment.pval <= pval_thr]
  low_f <- fil_low[low_wil_pval <= pval_thr & coment.pval <= pval_thr]
  
  high_f_05pval <- fil_high_05pval[high_wil_pval <= pval_thr_2 & coment.pval <= pval_thr_2]
  low_f_05pval <- fil_low_05pval[low_wil_pval <= pval_thr_2 & coment.pval <= pval_thr_2]
  
  merged <- merge(high, low, by = c("hpo", "tissue"))
  
  merged_05pval <- merge(high_05pval, low_05pval, by = c("hpo", "tissue"))
  
  merged_f <- merge(high_f, low_f, by = c("hpo", "tissue"))
  
  merged_f_05pval <- merge(high_f_05pval, low_f_05pval, by = c("hpo", "tissue"))
  
  ###
  
  to.add <- data.table(tissue = tis,
                       TP_high = nrow(high),
                       TP_low = nrow(low),
                       common = nrow(merged),
                       TPf_high = nrow(high_f),
                       TPf_low = nrow(low_f),
                       common_f = nrow(merged_f),
                       TP05_high = nrow(high_05pval),
                       TP05_low = nrow(low_05pval),
                       common_05 = nrow(merged_05pval),
                       TPf05_high = nrow(high_f_05pval),
                       TPf05_low = nrow(low_f_05pval),
                       common_f_05 = nrow(merged_f_05pval)
                       )
  
  results <- rbind(results, to.add)
  
}

fwrite(results,
       file.path(urales_home,
                 "TFM/results_w_comention/shared_TP.tsv"),
       sep = "\t")