#!/usr/bin/env Rscript

# Sergio Alías, 20240221
# Last modified 20230222

# Script for applying FDR to p-values (grouped by HPO)

library(data.table)
library(dplyr)

urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias/TFM"
urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias/TFM"
urales_home <- "/run/user/1001/gvfs/sftp:host=urales,user=salias/home/salias/TFM"

for (metric in c("fold_change/fc", "fold_change/fcst", "coex-analysis/HPA")){

if (metric == "coex-analysis/HPA"){
  metric2 <- "HPA"
  tissues <- c("adipose-tissue",
               "blood",
               "brain",
               "breast",
               "bronchus",
               "colon",
               "endometrium",
               "esophagus",
               "eye",
               "heart-muscle",
               "kidney",
               "liver",
               "lung",
               "ovary",
               "pancreas",
               "prostate-gland",
               "skeletal-muscle-organ",
               "skin",
               "small-intestine",
               "spleen",
               "stomach",
               "testis")
} else {
  if (metric == "fold_change/fc"){
    metric2 <- "fc"
  } else if (metric == "fold_change/fcst"){
    metric2 <- "fcst"
  }
  tissues <- c("adipose-tissue",
               "blood",
               "brain",
               "breast",
               "bronchus",
               "colon",
               "endometrium",
               "esophagus",
               "eye",
               "heart",
               "kidney",
               "liver",
               "lung",
               "ovary",
               "pancreas",
               "prostate-gland",
               "skeletal-muscle-organ",
               "skin-of-body",
               "small-intestine",
               "spleen",
               "stomach",
               "testis") 
}
for (tissue in tissues){

dt <- fread(file.path(urales_home,
                      metric,
                      tissue,
                      paste0("wil_results_",
                             metric2,
                             "_",
                             tissue,
                             ".tsv")))


dt <- dt %>%
  group_by(hpo) %>%
  mutate(wilcoxon.pval = p.adjust(wilcoxon.pval, method = "fdr"),
         high_wil_pval = p.adjust(high_wil_pval, method = "fdr"),
         low_wil_pval = p.adjust(low_wil_pval, method = "fdr"))

fwrite(dt,
       file = file.path(urales_home,
                         metric,
                         tissue,
                         paste0("wil_results_",
                                metric2,
                                "_",
                                tissue,
                                "_corrected.tsv")),
       sep = "\t")

}
}
