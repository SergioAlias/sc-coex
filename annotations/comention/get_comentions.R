#!/usr/bin/env Rscript

# Sergio Alías, 20230310
# Last modified 20230315

# Script for getting commentions that explain coex results
# We need:
#   - Coex analysis results -> wil_results-[HPA/TabulaSapiens]-[tissue].tsv
#   - Comention files -> ALL_CL-HPO_add_s.tsv /// ALL_ACT-HPO_add_s.tsv


### Packages

library(data.table)


### Input files setup

tissue <- "blood"
dataset <- "HPA"

results <- fread(file.path("../../coex-analysis/",
                           dataset,
                           tissue,
                           paste0("wil_results_",
                                  dataset,
                                  "_",
                                  tissue,
                                  ".tsv")
                           )
                 )

results <- results[wicoxon.sign != "ns"]


comentions <- fread("ALL_ACT-HPO_add_s.tsv")

colnames(comentions) <- c("ACT.id",
                          "ACT.name",
                          "HPO.id",
                          "HPO.name",
                          "str.sim",
                          "mentions.ACT",
                          "mentions.HPO",
                          "comentions",
                          "ratio",
                          "p.val")

### Filtering comention file by HPO terms

HPO.comentions <- comentions[HPO.id %in% unique(gsub("HP", "HP:", results$hpo))]

rm(comentions)

### Check if the comention is between the HPO and the desired cell type


sign.results <- data.table()

for (i in 1:nrow(results)){
  for (j in 1:nrow(HPO.comentions)){
    if (gsub("HP", "HP:", results[i, hpo]) == HPO.comentions[j, HPO.id]){
      if (results[i, annotation] == HPO.comentions[j, ACT.name]){
        row.to.add <- cbind(results[i], HPO.comentions[j, list(str.sim, comentions, ratio, p.val)])
        sign.results <- rbind(sign.results, row.to.add)
      }
    }
  }
}


fwrite(sign.results, file = paste0(file.path("..", "..", "results_w_comention", paste0(tissue, "_", dataset, "_results_w_comention.tsv"))), sep = "\t")
