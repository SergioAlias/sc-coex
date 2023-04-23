#!/usr/bin/env Rscript

# Sergio Alías, 20230310
# Last modified 20230420

# Script for getting commentions that explain coex results
# We need:
#   - Coex analysis results -> wil_results-[HPA/TabulaSapiens]-[tissue].tsv
#   - Comention files -> ALL_CL-HPO_add_s.tsv /// ALL_ACT-HPO_add_s.tsv


### Packages

library(data.table)


### Input files setup

tissues <- c("blood", "liver", "lung", "pancreas", "spleen", "stomach", "testis")
dataset <- "HPA"
urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"
urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"


for (i in seq_along(tissues)){

tissue <- tissues[i]
  
results <- fread(file.path(urales_home,
                           "TFM/coex-analysis",
                           dataset,
                           tissue,
                           paste0("wil_results_",
                                  dataset,
                                  "_",
                                  tissue,
                                  ".tsv")
                           )
                 )

annotations <- fread(file.path(urales_home,
                               "/TFM/annotations/cluster-annotation",
                               dataset,
                               paste0(tissue,
                                      "-cluster-annotation")),
                     header = FALSE)

#results <- results[wicoxon.sign != "ns"] # For keeping only significant wilcoxon p-values 
#results <- results[wilcoxon.pval <= 10**-3] # For keeping only significant wilcoxon p-values 

comentions <- fread(file.path(urales_home,
                              "TFM/annotations/comention/ALL_ACT-HPO_add_s.tsv"))

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

### separate cluster numbers by cell type

strings <- unique(annotations[[3]])

result_list <- list()

for (string in strings) {
  matching_rows <- annotations[grepl(string, V3)]
  matching_numbers <- as.numeric(sub("^c-", "", as.character(matching_rows[[2]])))
  result_list[[string]] <- matching_numbers
}


### Filtering comention file by HPO terms

HPO.comentions <- comentions[HPO.id %in% unique(gsub("HP", "HP:", results$hpo))]

rm(comentions)


### Check if the comention is between the HPO and the desired cell type


sign.results <- data.table()

for (i in 1:nrow(results)){
  has.comention <- FALSE
  for (j in 1:nrow(HPO.comentions)){
    if (gsub("HP", "HP:", results[i, hpo]) == HPO.comentions[j, HPO.id]){
      if (results[i, annotation] == HPO.comentions[j, ACT.name]){
        row.to.add <- cbind(results[i], HPO.comentions[j, list(str.sim, comentions, ratio, coment.pval)])
        sign.results <- rbind(sign.results, row.to.add)
        has.comention <- TRUE
      }
    }
  }
  if (!has.comention){
    row.to.add <- cbind(results[i], data.table(str.sim = 1,
                                               comentions = 1,
                                               ratio = 1,
                                               coment.pval = 1))
    sign.results <- rbind(sign.results, row.to.add)
  }
}


fwrite(sign.results, file = paste0(file.path(urales_home,
                                             "TFM/results_w_comention",
                                             tissue,
                                             paste0(tissue,
                                                    "_",
                                                    dataset,
                                                    "_results_w_comention_ALL.tsv"))),
       sep = "\t")
#}

###########    PARA QUEDARSE SOLO CON EL MÍNIMO DE CADA PAR HPO-CELLTYPE ########

truth_tables_filt <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric()) # solo la primera vez

only.lowest.results <- sign.results[, .SD[which.min(wilcoxon.pval)], by = .(hpo, annotation)]

pval_thr <- 0.001

# Calculate TP, TN, FP, FN
TP <- only.lowest.results[wilcoxon.pval <= pval_thr & coment.pval <= pval_thr, .N]
TN <- only.lowest.results[wilcoxon.pval > pval_thr & coment.pval > pval_thr, .N]
FP <- only.lowest.results[wilcoxon.pval <= pval_thr & coment.pval > pval_thr, .N]
FN <- only.lowest.results[wilcoxon.pval > pval_thr & coment.pval <= pval_thr, .N]

# Print the counts
cat("TP:", TP, "\n")
cat("TN:", TN, "\n")
cat("FP:", FP, "\n")
cat("FN:", FN, "\n")


to.add <- data.frame(TP = TP, TN = TN, FP = FP, FN = FN)
truth_tables_filt <- rbind(truth_tables_filt, to.add)
rownames(truth_tables_filt)[nrow(truth_tables_filt)] <- tissue





#################################################################################



### PARA LOS LOG10

# min_nonzero_coment_pval <- min(sign.results[coment.pval != 0][["coment.pval"]])
# 
# sign.results[, `:=`(log10.wilcoxon.pval = -log10(wilcoxon.pval + (min_nonzero_coment_pval/10)),
#                     log10.coment.pval = -log10(coment.pval + (min_nonzero_coment_pval/10)))]
# 
# fwrite(sign.results, file = paste0(file.path("..",
#                                              "..",
#                                              "results_w_comention",
#                                              tissue,
#                                              paste0(tissue,
#                                                     "_",
#                                                     dataset,
#                                                     "_results_w_comention.tsv"))),
#        sep = "\t")



### For loop considering we want results per cell type
# for (i in seq_along(result_list)){
#   
#   clusters_to_keep <- result_list[[i]]
#   celltype_name <- names(result_list)[i]
#   
#   # filter the data.table to keep only the rows with matching numbers
#   results_filtered <- sign.results[sapply(strsplit(tissue, "-"), function(x) {
#     num <- as.numeric(x[2])
#     num %in% clusters_to_keep
#   })]
#   
#   min_nonzero_coment_pval <- min(results_filtered[coment.pval != 0][["coment.pval"]])
#   
#   results_filtered[, `:=`(log10.wilcoxon.pval = -log10(wilcoxon.pval + (min_nonzero_coment_pval/10)),
#                       log10.coment.pval = -log10(coment.pval + (min_nonzero_coment_pval/10)))]
#   
#   fwrite(results_filtered, file = paste0(file.path("..",
#                                                    "..",
#                                                    "results_w_comention",
#                                                    tissue,
#                                                    paste0(tissue,
#                                                           "_",
#                                                           gsub(" ", "-", celltype_name),
#                                                           "_",
#                                                           dataset,
#                                                           "_results_w_comention.tsv"))),
#          sep = "\t")
#   
#   
# }

### Truth table creation

#truth_tables <- data.frame(TP = numeric(), TN = numeric(), FP = numeric(), FN = numeric())


#### TP, TN, FP, FN

pval_thr <- 0.001

# Calculate TP, TN, FP, FN
TP <- sign.results[wilcoxon.pval <= pval_thr & coment.pval <= pval_thr, .N]
TN <- sign.results[wilcoxon.pval > pval_thr & coment.pval > pval_thr, .N]
FP <- sign.results[wilcoxon.pval <= pval_thr & coment.pval > pval_thr, .N]
FN <- sign.results[wilcoxon.pval > pval_thr & coment.pval <= pval_thr, .N]

# Print the counts
cat("TP:", TP, "\n")
cat("TN:", TN, "\n")
cat("FP:", FP, "\n")
cat("FN:", FN, "\n")


to.add <- data.frame(TP = TP, TN = TN, FP = FP, FN = FN)
truth_tables <- rbind(truth_tables, to.add)
rownames(truth_tables)[nrow(truth_tables)] <- tissue


}


# res_table <- data.frame(
#   Negative = c(TN, FN, TN + FN),
#   Positive = c(TP, FP, TP + FP),
#   Total = c(TN + TP, FN + FP, TN + FN + TP + FP)
# )
# 
# # Set row names
# rownames(res_table) <- c("True", "False", "Total")
# 
# 
# # Export as TSV file
# write.table(result, file = "output.tsv", sep = "\t", quote = FALSE)

write.table(truth_tables_filt,
            file = file.path(urales_home,
                             "TFM/results_w_comention/confusion_table_filtered.tsv"),
            sep = "\t",
            quote = FALSE,
            col.names = NA)

write.table(truth_tables,
            file = file.path(urales_home,
                             "TFM/results_w_comention/confusion_table.tsv"),
            sep = "\t",
            quote = FALSE,
            col.names = NA)