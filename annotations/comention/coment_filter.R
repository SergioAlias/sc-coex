#!/usr/bin/env Rscript

# Sergio Alías, 20221209
# Last modified 20221209

library(data.table)

coment_hpo_cl <- fread("../annotations/comention/ALL_CL-HPO_add_s.tsv",
                       col.names = c("CL.ID", "CL.name", "HPO.ID", "HPO.name",
                                     "str.sim", "n.papers.CL", "n.papers.HPO",
                                     "n.papers.both", "comention", "pval"))
coment_hpo_act <- fread("../annotations/comention/ALL_ACT-HPO_add_s.tsv",
                        col.names = c("ACT.ID", "ACT.name", "HPO.ID", "HPO.name",
                                      "str.sim", "n.papers.ACT", "n.papers.HPO",
                                      "n.papers.both", "comention", "pval"))

results <- fread("../coex-analysis/HPA/pancreas/wil_results_pancreas.tsv")
results <- results[wicoxon.sign != "ns"]

filtered_cl <- results[hpo %in% coment_hpo_cl$HPO.ID]

filtered_act <- results[hpo %in% coment_hpo_act$HPO.ID]

if (nrow(filtered_cl) == 0 & nrow(filtered_act) == 0){
  message("No comention found")
} else {
  
}
