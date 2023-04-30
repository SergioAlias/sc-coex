#!/usr/bin/env Rscript

# Sergio Alías, 20230310
# Last modified 20230428

# Script for getting commentions that explain coex / fc / fcst results
# We need:
#   - Coex analysis results -> wil_results-[HPA/TabulaSapiens]-[tissue].tsv
#   - Comention files -> ALL_CL-HPO_add_s.tsv /// ALL_ACT-HPO_add_s.tsv


### Packages

library(data.table)


### Input files setup

tissues <- c("blood", "liver", "lung", "pancreas", "spleen", "stomach", "testis")
# tissues <- c("adipose-tissue",
#              "blood",
#              "brain",
#              "breast",
#              "bronchus",
#              "colon",
#              "endometrium",
#              "esophagus",
#              "eye",
#              "heart",
#              "kidney",
#              "liver",
#              "lung",
#              "lymph-node",
#              "ovary",
#              "pancreas",
#              #"placenta",
#              "prostate-gland",
#              #"rectum",
#              "skeletal-muscle-organ",
#              "skin-of-body",
#              "small-intestine",
#              "spleen",
#              "stomach",
#              "testis")
dataset <- "HPA"
urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"
# urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"


for (t in seq_along(tissues)){

tissue <- tissues[t]

if (!file.exists(file.path(urales_home,
                           "TFM/results_w_comention",
                           tissue))) {
  dir.create(file.path(urales_home,
                       "TFM/results_w_comention",
                       tissue))
}


results <- fread(file.path(urales_home,
                           "TFM/coex-analysis/HPA",
                           tissue,
                           paste0("wil_results_HPA_",
                                  tissue,
                                  #"_corrected.tsv"))
                                  ".tsv"))
                 )

annotations <- fread(file.path(urales_home,
                               "/TFM/annotations/cluster-annotation",
                               dataset,
                               paste0(tissue,
                                      "-cluster-annotation")),
                     header = FALSE)

#-----#


#results <- results[wicoxon.sign != "ns"] # For keeping only significant wilcoxon p-values 
#results <- results[wilcoxon.pval <= 10**-3] # For keeping only significant wilcoxon p-values 


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

### separate cluster numbers by cell type

strings <- unique(annotations[[3]])

result_list <- list()

for (string in strings) {
  matching_rows <- annotations[grepl(string, V3)]
  matching_numbers <- as.numeric(sub("^c-", "", as.character(matching_rows[[2]])))
  result_list[[string]] <- matching_numbers
}


### Filtering comention file by HPO terms

HPO.comentions <- comentions[HPO.id %in% unique(gsub("HP", "HP:", results$hpo))] # coex
# HPO.comentions <- comentions[HPO.id %in% unique(results$hpo)] # fc, fcst

rm(comentions)


### Check if the comention is between the HPO and the desired cell type

message("Generating file")

sign.results <- data.table()

for (i in 1:nrow(results)){
  message("Row ", i, "/", nrow(results), " ", t, "...")
  has.comention <- FALSE
  for (j in 1:nrow(HPO.comentions)){
    if (gsub("HP", "HP:", results[i, hpo]) == HPO.comentions[j, HPO.id]){ # coex
    # if (results[i, hpo] == HPO.comentions[j, HPO.id]){ # fc, fcst
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
                                                    "_05pval_coex_results_w_comention.tsv"))),
       sep = "\t")


#-----#

message("Done")

}