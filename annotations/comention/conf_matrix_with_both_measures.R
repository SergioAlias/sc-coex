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

pval_thr <- 0.001 # 0.001


tissues_coex <- c("adipose-tissue",
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

tissues_fc <- c("adipose-tissue",
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

truth_tables_filt <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric())

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


for (i in seq_along(tissues_coex)){
  
  tissue_coex <- tissues_coex[i]
  tissue_fc <- tissues_fc[i]
  
  if (!file.exists(file.path(urales_home,
                             "TFM/results_w_comention",
                             tissue_coex))) {
    dir.create(file.path(urales_home,
                         "TFM/results_w_comention",
                         tissue_coex))
  }

  results_coex <- fread(file.path(urales_home,
                                  "TFM/res_finales/coex",
                                  paste0("wil_results_HPA_",
                                         tissue_coex,
                                         "_no_blind.tsv")))

  results_fc <- fread(file.path(urales_home,
                                "TFM/res_finales/fc",
                                paste0("wil_results_fc_",
                                       tissue_fc,
                                       "_no_blind.tsv")))
  
  annotations <- fread(file.path(urales_home,
                                 "/TFM/annotations/cluster-annotation/HPA",
                                 paste0(tissue_coex,
                                        "-cluster-annotation")),
                       header = FALSE)
  strings <- unique(annotations[[3]])
  
  result_list <- list()
  
  for (string in strings) {
    matching_rows <- annotations[grepl(string, V3)]
    matching_numbers <- as.numeric(sub("^c-", "", as.character(matching_rows[[2]])))
    result_list[[string]] <- matching_numbers
  }
  

  HPO.comentions <- comentions[HPO.id %in% unique(gsub("HP", "HP:", results_coex$hpo))]

  
# COMIENZO BUCLE PVAL COLS
  
  pval_col_coex <- "wilcoxon.pval"
  pval_col_fc <- "high_wil_pval"
  
  ks_values_coex <- dcast(results_coex, tissue ~ hpo, value.var = pval_col_coex)
  ks_values_fc <- dcast(results_fc, tissue ~ hpo, value.var = pval_col_fc)
  
  ks_values_coex[, c("text", "num") := tstrsplit(tissue, "(?<=\\D)(?=\\d)", perl = TRUE)] # para asegurarnos de que están bien ordenadas las filas 
  ks_values_fc[, c("text", "num") := tstrsplit(tissue, "(?<=\\D)(?=\\d)", perl = TRUE)] # para asegurarnos de que están bien ordenadas las filas 
  
  ks_values_coex[, num := as.numeric(num)]
  ks_values_fc[, num := as.numeric(num)]
  
  ks_values_coex <- ks_values_coex[order(num)]
  ks_values_fc <- ks_values_fc[order(num)]
  
  ks_values_coex[, c("text", "num") := NULL]
  ks_values_fc[, c("text", "num") := NULL]
  
  rownames(ks_values_coex) <- ks_values_coex$tissue
  rownames(ks_values_fc) <- ks_values_fc$tissue
  
  ks_values_coex$tissue <- annotations$V3
  ks_values_fc$tissue <- annotations$V3
  
  
  coment_values <- matrix(1, nrow = nrow(ks_values_coex), ncol = ncol(ks_values_coex))
  rownames(coment_values) <- rownames(ks_values_coex)
  colnames(coment_values) <- colnames(ks_values_coex)
  coment_values <- as.data.frame(coment_values)
  coment_values$tissue <- annotations$V3
  
  for (r in seq_along(rownames(coment_values))){
    for (c in 2:ncol(coment_values)){
        possible_row <- HPO.comentions[ACT.name == coment_values$tissue[r] & HPO.id == gsub("HP", "HP:", colnames(coment_values)[c]),]
      if (nrow(possible_row) != 0){
        coment_values[r, c] <- possible_row$coment.pval
      }
    }
  }
  
  # Only one minimum
  ks_values_coex <- as.data.frame(ks_values_coex)
  ks_values_fc <- as.data.frame(ks_values_fc)
  
  ks_values_filtered_coex <- ks_values_coex %>%
    group_by(tissue) %>%
    mutate(across(starts_with("HP"), ~ ifelse(rank(.x, ties.method = "first") == 1, .x, NA))) %>%
    ungroup()
  
  ks_values_filtered_fc <- ks_values_fc %>%
    group_by(tissue) %>%
    mutate(across(starts_with("HP"), ~ ifelse(rank(.x, ties.method = "first") == 1, .x, NA))) %>%
    ungroup()
  
  ks_values_filtered_coex <- as.data.frame(ks_values_filtered_coex)
  ks_values_filtered_fc <- as.data.frame(ks_values_filtered_fc)
  
  colnames(ks_values_filtered_fc) <- gsub(":", "", colnames(ks_values_filtered_fc))
  
  coment_values <- coment_values %>% distinct(tissue, .keep_all = TRUE)
  
  ks_values_filtered_fc <- ks_values_filtered_fc %>%
    group_by(tissue) %>%
    summarise(across(everything(), ~ na.omit(.)[1]))
  
  ks_values_filtered_coex <- ks_values_filtered_coex %>%
    group_by(tissue) %>%
    summarise(across(everything(), ~ na.omit(.)[1]))
  
  sort_dataframe <- function(df) {
    df %>%
      arrange(tissue) %>%         # Sort rows by 'tissue' column
      select(order(colnames(df))) # Order columns alphabetically
    df %>%
      select(tissue, everything())
  }
  
  rownames(coment_values) <- NULL
  
  ks_values_filtered_fc <- sort_dataframe(ks_values_filtered_fc)
  ks_values_filtered_coex <- sort_dataframe(ks_values_filtered_coex)
  coment_values <- sort_dataframe(coment_values)
  coment_values <- coment_values[match(ks_values_filtered_coex$tissue, coment_values$tissue), ]
  
  TP_f <- 0
  TN_f <- 0
  FP_f <- 0
  FN_f <- 0
  NE_f <- 0
  
  for (r in seq_along(rownames(coment_values))){
    for (c in 2:ncol(coment_values)){
      sign_coment <- coment_values[r, c] < pval_thr
      sign_ks_f_coex <- ks_values_filtered_coex[r, c] < pval_thr
      sign_ks_f_fc <- ks_values_filtered_fc[r, c] < pval_thr
      if (all(!sign_ks_f_coex, !sign_ks_f_fc, sign_coment)){
        NE_f <- NE_f + 1
      }
      if (all(sign_ks_f_coex, sign_ks_f_fc, sign_coment)){
        TP_f <- TP_f + 1
      } else if (all(sign_ks_f_coex, sign_ks_f_fc,  !sign_coment)){
        FP_f <- FP_f + 1
      } else if (!sign_coment){
        TN_f <- TN_f + 1
      } else {
        FN_f <- FN_f + 1
      }
    }
  }
    
    to.add <- data.frame(TP = TP_f, TN = TN_f, FP = FP_f, FN = FN_f)
    truth_tables_filt <- rbind(truth_tables_filt, to.add)
    rownames(truth_tables_filt)[nrow(truth_tables_filt)] <- tissue_coex

  }



sums <- colSums(truth_tables_filt)
truth_tables_filt <- rbind(truth_tables_filt, sums)
rownames(truth_tables_filt)[nrow(truth_tables_filt)] <- "all"

write.table(truth_tables_filt,
            file = file.path(urales_home,
                             "TFM/results_w_comention",
                             "both_confusion_table_filtered_no_blind.tsv"),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

