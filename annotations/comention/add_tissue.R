#!/usr/bin/env Rscript

# Sergio Alías, 20230303
# Last modified 20230307

# Script for adding tissue row (and dataset rows)

# CL includes Cell Line Ontology terms (almost all included in TabulaSapiens)
# ACT includes the 76 HPA cluster names and 4 Tabula Sapien cluster names
# So for HPA - use the ACT discarding those 4 terms from TabSap (it seems that they have no comentions so consider that done)
# And for TabSap - use both the CL and the 4 terms of the ACT (again - they don't have comentions so just use CL)

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

act.HPA    <- readLines("HPA-CellTypes.txt")        # cell types as in HPA (not following CL nomenclature)
act.TabSap <- readLines("missing-terms-TabSap.txt") # some terms that TabulaSapiens people just made up

HPA.clus.folder <- list.files("../cluster-annotation/HPA/")


HPA.comentions <- c()

for (i in seq_along(HPA.clus.folder)){
  clus.data <- fread(file.path("..", "cluster-annotation", "HPA", HPA.clus.folder[i]),
                     header = FALSE)
  name <- unlist(strsplit(HPA.clus.folder[i], "-cluster-"))[1]
  n.coment <- nrow(coment_hpo_act[coment_hpo_act$ACT.name %in% unique(clus.data$V3)])
  HPA.comentions[name] <- n.coment
}

HPA.coment.dt <- data.table(tissue = names(HPA.comentions),
                           n.comentions = HPA.comentions)
HPA.coment.dt <- HPA.coment.dt[order(-n.comentions),]


fwrite("")