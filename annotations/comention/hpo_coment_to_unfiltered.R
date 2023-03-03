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
unfiltered_results <- fread("tis_hpo_ngenes_HPA.tsv")

filtered_act <- unfiltered_results[hpo_code %in% coment_hpo_act$HPO.ID]
