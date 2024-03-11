#!/usr/bin/env Rscript

# Sergio Alías, 20230426
# Last modified 20240226

# Script for getting confusion matrices that explain coex / fc / fcst results (per celltype, NOT cluster)

### Packages

library(data.table)
library(dplyr)

### Input specs

urales_home <- "/run/user/1001/gvfs/sftp:host=urales,user=salias/home/salias"
#urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"
#urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"

metric <- "fc" # coex, fc, fcst
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

truth_tables_filt <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric())
truth_tables <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric())

truth_tables_high_filt <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric())
truth_tables_high <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric())

truth_tables_low_filt <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric())
truth_tables_low <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric())

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
  
  if (!file.exists(file.path(urales_home,
                             "TFM/results_w_comention",
                             tissue))) {
    dir.create(file.path(urales_home,
                         "TFM/results_w_comention",
                         tissue))
  }
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
  
  annotations <- fread(file.path(urales_home,
                                 "/TFM/annotations/cluster-annotation/HPA",
                                 paste0(tissue,
                                        "-cluster-annotation")),
                       header = FALSE)
  strings <- unique(annotations[[3]])
  
  result_list <- list()
  
  for (string in strings) {
    matching_rows <- annotations[grepl(string, V3)]
    matching_numbers <- as.numeric(sub("^c-", "", as.character(matching_rows[[2]])))
    result_list[[string]] <- matching_numbers
  }
  
  if (metric == "coex"){
    HPO.comentions <- comentions[HPO.id %in% unique(gsub("HP", "HP:", results$hpo))] # coex
  } else {
    HPO.comentions <- comentions[HPO.id %in% unique(results$hpo)] # fc, fcst
  }
  
  for (pval_col in c("wilcoxon.pval", "high_wil_pval", "low_wil_pval")){
  
  ks_values <- dcast(results, tissue ~ hpo, value.var = pval_col)
  ks_values[, c("text", "num") := tstrsplit(tissue, "(?<=\\D)(?=\\d)", perl = TRUE)] # para asegurarnos de que están bien ordenadas las filas 
  ks_values[, num := as.numeric(num)]
  ks_values <- ks_values[order(num)]
  ks_values[, c("text", "num") := NULL]
  rownames(ks_values) <- ks_values$tissue
  ks_values$tissue <- annotations$V3
  
  coment_values <- matrix(1, nrow = nrow(ks_values), ncol = ncol(ks_values))
  rownames(coment_values) <- rownames(ks_values)
  colnames(coment_values) <- colnames(ks_values)
  coment_values <- as.data.frame(coment_values)
  coment_values$tissue <- annotations$V3
  
  for (r in seq_along(rownames(coment_values))){
    for (c in 2:ncol(coment_values)){
      if (metric == "coex"){
        possible_row <- HPO.comentions[ACT.name == coment_values$tissue[r] & HPO.id == gsub("HP", "HP:", colnames(coment_values)[c]),]
      } else {
        possible_row <- HPO.comentions[ACT.name == coment_values$tissue[r] & HPO.id == colnames(coment_values)[c],]
      }
      if (nrow(possible_row) != 0){
        coment_values[r, c] <- possible_row$coment.pval
      }
    }
  }
  
  
  # Keep ties 
  #ks_values <- as.data.frame(ks_values)
  #ks_values_filtered <- ks_values %>%
   # group_by(tissue) %>%
    #mutate(across(starts_with("HP"), ~ ifelse(.x == min(.x), .x, NA))) %>%
    #ungroup()
  #ks_values_filtered <- as.data.frame(ks_values_filtered)
  
  # Only one minimum
  ks_values <- as.data.frame(ks_values)
  ks_values_filtered <- ks_values %>%
    group_by(tissue) %>%
    mutate(across(starts_with("HP"), ~ ifelse(rank(.x, ties.method = "first") == 1, .x, NA))) %>%
    ungroup()
  ks_values_filtered <- as.data.frame(ks_values_filtered)
  
  
  TP <- 0
  TN <- 0
  FP <- 0
  FN <- 0
  
  TP_f <- 0
  TN_f <- 0
  FP_f <- 0
  FN_f <- 0
  
  for (r in seq_along(rownames(coment_values))){
    for (c in 2:ncol(coment_values)){
      sign_ks <- ks_values[r, c] <= pval_thr
      sign_coment <- coment_values[r, c] <= pval_thr
      sign_ks_f <- "nope"
      if(!is.na(ks_values_filtered[r, c])){
        sign_ks_f <- ks_values_filtered[r, c] <= pval_thr
      }
      if (sign_ks & sign_coment){
        TP <- TP + 1
      } else if (!sign_ks & !sign_coment){
        TN <- TN + 1
      } else if (sign_ks & !sign_coment){
        FP <- FP + 1
      } else if (!sign_ks & sign_coment){
        FN <- FN + 1
      }
      if (sign_ks_f != "nope"){
      if (sign_ks & sign_coment){
        TP_f <- TP_f + 1
      } else if (!sign_ks & !sign_coment){
        TN_f <- TN_f + 1
      } else if (sign_ks & !sign_coment){
        FP_f <- FP_f + 1
      } else if (!sign_ks & sign_coment){
        FN_f <- FN_f + 1
      }
      }
    }
  }
  
  if (pval_col == "wilcoxon.pval"){
    to.add <- data.frame(TP = TP, TN = TN, FP = FP, FN = FN)
    truth_tables <- rbind(truth_tables, to.add)
    rownames(truth_tables)[nrow(truth_tables)] <- tissue
    
    to.add <- data.frame(TP = TP_f, TN = TN_f, FP = FP_f, FN = FN_f)
    truth_tables_filt <- rbind(truth_tables_filt, to.add)
    rownames(truth_tables_filt)[nrow(truth_tables_filt)] <- tissue
  } else if (pval_col == "high_wil_pval"){
    to.add <- data.frame(TP = TP, TN = TN, FP = FP, FN = FN)
    truth_tables_high <- rbind(truth_tables_high, to.add)
    rownames(truth_tables_high)[nrow(truth_tables_high)] <- tissue
    
    to.add <- data.frame(TP = TP_f, TN = TN_f, FP = FP_f, FN = FN_f)
    truth_tables_high_filt <- rbind(truth_tables_high_filt, to.add)
    rownames(truth_tables_high_filt)[nrow(truth_tables_high_filt)] <- tissue
  } else if (pval_col == "low_wil_pval"){
    to.add <- data.frame(TP = TP, TN = TN, FP = FP, FN = FN)
    truth_tables_low <- rbind(truth_tables_low, to.add)
    rownames(truth_tables_low)[nrow(truth_tables_low)] <- tissue
    
    to.add <- data.frame(TP = TP_f, TN = TN_f, FP = FP_f, FN = FN_f)
    truth_tables_low_filt <- rbind(truth_tables_low_filt, to.add)
    rownames(truth_tables_low_filt)[nrow(truth_tables_low_filt)] <- tissue
  }

  }
}


sums <- colSums(truth_tables)
truth_tables <- rbind(truth_tables, sums)
sums <- colSums(truth_tables_filt)
truth_tables_filt <- rbind(truth_tables_filt, sums)
sums <- colSums(truth_tables_high)
truth_tables_high <- rbind(truth_tables_high, sums)
sums <- colSums(truth_tables_high_filt)
truth_tables_high_filt <- rbind(truth_tables_high_filt, sums)
sums <- colSums(truth_tables_low)
truth_tables_low <- rbind(truth_tables_low, sums)
sums <- colSums(truth_tables_low_filt)
truth_tables_low_filt <- rbind(truth_tables_low_filt, sums)
  
rownames(truth_tables)[nrow(truth_tables)] <- "all"
rownames(truth_tables_filt)[nrow(truth_tables_filt)] <- "all"
rownames(truth_tables_high)[nrow(truth_tables_high)] <- "all"
rownames(truth_tables_high_filt)[nrow(truth_tables_high_filt)] <- "all"
rownames(truth_tables_low)[nrow(truth_tables_low)] <- "all"
rownames(truth_tables_low_filt)[nrow(truth_tables_low_filt)] <- "all"




write.table(truth_tables,
            file = file.path(urales_home,
                             "TFM/results_w_comention",
                             paste0(metric,
                                    "_extreme_confusion_table_corr_byct.tsv")),
            sep = "\t",
            quote = FALSE,
            col.names = NA)


write.table(truth_tables_filt,
            file = file.path(urales_home,
                             "TFM/results_w_comention",
                             paste0(metric,
                                    "_extreme_confusion_table_filtered_corr_byct.tsv")),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

#-----#

write.table(truth_tables_high,
            file = file.path(urales_home,
                             "TFM/results_w_comention",
                             paste0(metric,
                                    "_high_confusion_table_corr_byct.tsv")),
            sep = "\t",
            quote = FALSE,
            col.names = NA)


write.table(truth_tables_high_filt,
            file = file.path(urales_home,
                             "TFM/results_w_comention",
                             paste0(metric,
                                    "_high_confusion_table_filtered_corr_byct.tsv")),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

#-----#

write.table(truth_tables_low,
            file = file.path(urales_home,
                             "TFM/results_w_comention",
                             paste0(metric,
                                    "_low_confusion_table_corr_byct.tsv")),
            sep = "\t",
            quote = FALSE,
            col.names = NA)


write.table(truth_tables_low_filt,
            file = file.path(urales_home,
                             "TFM/results_w_comention",
                             paste0(metric,
                                    "_low_confusion_table_filtered_corr_byct.tsv")),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

