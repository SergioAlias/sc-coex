#!/usr/bin/env Rscript

# Sergio Alías, 20240605
# Last modified 20240605

# Script for keeping only phenotypes with at least one significative comention (remove blind predictions)

### Packages

library(data.table)
library(dplyr)

### Input specs

urales_home <- "/run/user/1001/gvfs/sftp:host=urales,user=salias/home/salias"
#urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"
#urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"

metric <- "coex" # coex, fc, fcst
pval_thr <- 0.001 # 0.001

if (metric == "coex"){
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

comentions <- fread(file.path(urales_home,
                              "TFM/annotations/comention/ALL_ACT-HPO_0.05_add_s.tsv"))

colnames(comentions) <- c("ACT.id",
                          "ACT.name",
                          "HPO.id",
                          "HPO.name",
                          "str.sim",
                          "mentions.ACT",
                          "mentions.HPO",
                          "comentions",
                          "ratio",
                          "coment.pval")

for (i in seq_along(tissues)){
  
  tissue <- tissues[i]
  
  if (metric == "coex"){
    results <- fread(file.path(urales_home,
                               "TFM/coex-analysis/HPA",
                               tissue,
                               paste0("wil_results_HPA_",
                                      tissue,
                                      "_corrected.tsv")))
  } else {
    results <- fread(file.path(urales_home,
                               "TFM/fold_change",
                               metric,
                               tissue,
                               paste0("wil_results_",
                                      metric,
                                      "_",
                                      tissue,
                                      "_corrected.tsv")))
  }
  
  unique_strings <- tolower(unique(results[[2]]))
  
  strings_to_discard <- c()
  
  # Iterar sobre cada valor único
  for (string in unique_strings) {
    # Filtrar las filas de comentions donde ACT.name es igual a la cadena actual
    relevant_rows <- comentions[HPO.name == string]
    
    # Verificar si todos los valores en coment.pval son mayores o iguales a 0.001
    if (nrow(relevant_rows) == 0 || all(relevant_rows$coment.pval >= pval_thr)) {
      # Si es así, agregar la cadena a strings_to_discard
      strings_to_discard <- c(strings_to_discard, string)
    }
  }
  
  # Filtrar las filas de results para descartar aquellas con valores en strings_to_discard
  filtered_results <- results[!(tolower(results[[2]]) %in% strings_to_discard)]
  
  if (metric == "coex"){
    fwrite(filtered_results, sep = "\t", file = file.path(urales_home,
                               "TFM/coex-analysis/HPA",
                               tissue,
                               paste0("wil_results_HPA_",
                                      tissue,
                                      "_no_blind.tsv")))
  } else {
    fwrite(filtered_results, sep = "\t", file = file.path(urales_home,
                               "TFM/fold_change",
                               metric,
                               tissue,
                               paste0("wil_results_",
                                      metric,
                                      "_",
                                      tissue,
                                      "_no_blind.tsv")))
  }
  
}

