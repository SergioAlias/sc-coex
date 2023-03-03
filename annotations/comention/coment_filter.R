#!/usr/bin/env Rscript

# Sergio Alías, 20221209
# Last modified 20230303

library(data.table)

coment_hpo_cl <- fread("ALL_CL-HPO_add_s.tsv",
                       col.names = c("CL.ID", "CL.name", "HPO.ID", "HPO.name",
                                     "str.sim", "n.papers.CL", "n.papers.HPO",
                                     "n.papers.both", "comention", "pval"))
coment_hpo_act <- fread("ALL_ACT-HPO_add_s.tsv",
                        col.names = c("ACT.ID", "ACT.name", "HPO.ID", "HPO.name",
                                      "str.sim", "n.papers.ACT", "n.papers.HPO",
                                      "n.papers.both", "comention", "pval"))

results <- fread("../../coex-analysis/HPA/liver/wil_results_liver.tsv")
results <- results[wicoxon.sign != "ns"]
results.hpo <- gsub("HP", "HP:", results$hpo)


filtered_cl <- results[results.hpo %in% coment_hpo_cl$HPO.ID,]

filtered_act <- results[results.hpo %in% coment_hpo_act$HPO.ID,]

if (nrow(filtered_cl) == 0 & nrow(filtered_act) == 0){
  message("No comention found")
} else {message("Comentions found!")
}
