#!/usr/bin/env Rscript

# Sergio Alías, 20230430
# Last modified 20230430

# Script for getting number of cell-types comentioned per HPO per tissue
# We need:
#    - Annotation file -> TFM/annotations/filtered_tis_hpo_ngenes_HPA.tsv
#    - Comention files -> ALL_ACT-HPO_add_s.tsv /// ALL_ACT-HPO_0.05_add_s.tsv

### Packages

library(data.table)


### Input

urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"
urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"

tissues <- c("adipose-tissue",
             "blood",
             "brain",
             "breast",
             "bronchus",
             "colon",
             "endometrium",
             "esophagus",
             "eye",
             "heart",
             "kidney",
             "liver",
             "lung",
             "lymph-node",
             "ovary",
             "pancreas",
             #"placenta",
             "prostate-gland",
             #"rectum",
             "skeletal-muscle-organ",
             "skin-of-body",
             "small-intestine",
             "spleen",
             "stomach",
             "testis")

comentions <- fread(file.path(urales_home,
                              "TFM/annotations/comention/ALL_ACT-HPO_add_s.tsv"))

comentions_05pval <- fread(file.path(urales_home,
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

colnames(comentions_05pval) <- c("ACT.id",
                          "ACT.name",
                          "HPO.id",
                          "HPO.name",
                          "str.sim",
                          "mentions.ACT",
                          "mentions.HPO",
                          "comentions",
                          "ratio",
                          "coment.pval")

annotations <- fread(file.path(urales_home,
                               "TFM/annotations/filtered_tis_hpo_ngenes_HPA.tsv"))

results <- data.table(tissue = character(),
                      n_celltypes = numeric(),
                      hpo.code = character(),
                      hpo.name = character(),
                      total = numeric(),
                      intersect = numeric(),
                      total.05pval = numeric(),
                      intersect.05pval = numeric()
                      )

for (t in seq_along(tissues)){
  
  tis <- tissues[t]
  tis_celltypes <- unique(fread(file.path(urales_home,
                                "TFM/annotations/cluster-annotation/HPA",
                                paste0(tis, "-cluster-annotation")),
                     header = FALSE)$V3)
  tis_ann <- annotations[tissue == tis]
  
  for (h in 1:nrow(tis_ann)){
  hp.code <- tis_ann[h, hpo_code]
  hp.name <- tis_ann[h, hpo_name]
  
  cm.total <- nrow(comentions[HPO.id == hp.code])
  cm.total.05pval <- nrow(comentions_05pval[HPO.id == hp.code])
  inter <- nrow(comentions[HPO.id == hp.code & ACT.name %in% tis_celltypes])
  inter.05pval <- nrow(comentions_05pval[HPO.id == hp.code & ACT.name %in% tis_celltypes])
  
  results <- rbind(results, data.table(tissue = tis,
                                       n_celltypes = length(tis_celltypes),
                                       hpo.code = hp.code,
                                       hpo.name = hp.name,
                                       total = cm.total,
                                       intersect = inter,
                                       total.05pval = cm.total.05pval,
                                       intersect.05pval = inter.05pval))
  }
}  