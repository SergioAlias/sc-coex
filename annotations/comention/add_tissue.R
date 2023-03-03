#!/usr/bin/env Rscript

# Sergio Alías, 20230303
# Last modified 20230303

# Script for adding tissue row (and dataset rows)

library(data.table)

coment_hpo_cl <- fread("ALL_CL-HPO_add_s.tsv",
                       col.names = c("CL.ID", "CL.name", "HPO.ID", "HPO.name",
                                     "str.sim", "n.papers.CL", "n.papers.HPO",
                                     "n.papers.both", "comention", "pval")
                       )
coment_hpo_act <- fread("ALL_ACT-HPO_add_s.tsv",
                        col.names = c("ACT.ID", "ACT.name", "HPO.ID", "HPO.name",
                                      "str.sim", "n.papers.ACT", "n.papers.HPO",
                                      "n.papers.both", "comention", "pval")
                        )

act.HPA    <- readLines("HPA-CellTypes.txt")
act.TabSap <- readLines("missing-terms-TabSap.txt")