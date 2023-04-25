#!/usr/bin/env Rscript

# Sergio Alías, 20230416
# Last modified 20230424

# Script for getting metrics that explains confusion table results
# We need:
#    - Confusion table -> TFM/results_w_comention/truth_table.tsv ot truth_table_filtered.tsv


### Packages

library(data.table)


### Input

urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias/TFM"

file_names <- list.files(file.path(urales_home, "results_w_comention"), pattern = "\\.tsv$", full.names = FALSE)

file_names_full <- list.files(file.path(urales_home, "results_w_comention"), pattern = "\\.tsv$", full.names = TRUE)


for (i in seq_along(file_names_full)){

filename_full <- file_names_full[i]
filename <- file_names[i]  
  
table <- fread(filename_full)

colnames(table)[1] <- "tissue"

results <- data.table(tissue = character(),
                      TP = numeric(),
                      TN = numeric(),
                      FP = numeric(),
                      FN = numeric(),
                      TPR = numeric(),
                      TNR = numeric(),
                      PPV = numeric(),
                      # NPV = numeric(),
                      # FNR = numeric(),
                      # FPR = numeric(),
                      # FDR = numeric(),
                      # FOR = numeric(),
                      # LRplus = numeric(),
                      # LRminus = numeric(),
                      # PT = numeric(),
                      # CSI = numeric(),
                      # Prev = numeric(),
                      # ACC = numeric(),
                      F1 = numeric(),
                      MCC = numeric())#,
                      # FM = numeric(),
                      # BM = numeric(),
                      # MK = numeric(),
                      # DOR = numeric())

for (i in c(1:nrow(table))){
  tissue <- table[i, tissue]
  TP <- table[i, TP]
  TN <- table[i, TN]
  FP <- table[i, FP]
  FN <- table[i, FN]
  row.add <- data.table(tissue = tissue,
                        TP = TP,
                        TN = TN,
                        FP = FP,
                        FN = FN,
                        TPR = TP/(TP+FN),
                        TNR = TN/(TN+FP),
                        PPV = TP/(TP+FP),
                        # NPV = TN/(TN+FN),
                        # FNR = FN/(FN+TP),
                        # FPR = FP/(FP+TN),
                        # FDR = FP/(FP+TP),
                        # FOR = FN/(FN+TN),
                        # LRplus = (TP/(TP+FN))/(FP/(FP+TN)),
                        # LRminus = (FN/(FN+TP))/(TN/(TN+FP)),
                        # PT = sqrt(FP/(FP+TN))/(sqrt(TP/(TP+FN))+sqrt(FP/(FP+TN))),
                        # CSI = TP/(TP+FN+FP),
                        # Prev = (TP+FN)/(TP+TN+FP+FN),
                        # ACC = (TP+TN)/(TP+TN+FP+FN),
                        F1 = (2*TP)/((2*TP)+FP+FN),
                        MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))#,
                        # FM = sqrt((TP/(TP+FP))*(TP/(TP+FN))),
                        # BM = (TP/(TP+FN))+(TN/(TN+FP))-1,
                        # MK = (TP/(TP+FP))+(TN/(TN+FN))-1,
                        # DOR = ((TP/(TP+FN))/(FP/(FP+TN)))/((FN/(FN+TP))/(TN/(TN+FP))))
  results <- rbind(results, row.add)
  
}


fwrite(results,
       file = gsub("confusion_table", "metrics", filename_full),
       sep = "\t")

}