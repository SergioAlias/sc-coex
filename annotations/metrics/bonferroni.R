#!/usr/bin/env Rscript

# Sergio Alías, 20230324
# Last modified 20230424

# Script for applying Bonferroni correction to p-values

library(data.table)

urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias/TFM"
tissues <- c("blood", "liver", "lung", "pancreas", "spleen", "stomach", "testis")
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
             "lymph-node",
             "ovary",
             "pancreas",
             "placenta",
             "prostate-gland",
             "rectum",
             "skeletal-muscle-organ",
             "skin-of-body",
             "small-intestine",
             "spleen",
             "stomach",
             "testis")

for (metric in c("fc", "fcst")){

for (tissue in tissues){

dt <- fread(file.path(urales_home,
                      "fold_change",
                      metric,
                      tissue,
                      paste0("wil_results_",
                             metric,
                             "_",
                             tissue,
                             ".tsv")))


cols_to_correct <- c("wilcoxon.pval", "high_wil_pval", "low_wil_pval")
dt[, (cols_to_correct) := lapply(.SD, function(x) pmin(1, x * nrow(dt))), .SDcols = cols_to_correct]

cols_to_remove <- c("wicoxon.sign", "high_wil_sign", "low_wil_sign") # , "effsize", "effsize.sign")
dt <- dt[, !cols_to_remove, with = FALSE]


fwrite(dt,
       file = file.path(urales_home,
                         "fold_change",
                         metric,
                         tissue,
                         paste0("wil_results_",
                                metric,
                                "_",
                                tissue,
                                "_corrected.tsv")),
       sep = "\t")

}
}