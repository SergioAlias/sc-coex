#!/usr/bin/env Rscript
# Sergio Alías, 20220411
# Last modified 20220711

# Script para hacer subsetting de los tissues del dataset Tabula Sapiens por cell types (clusters)

# Genera los clusters separados en ficheros .tsv,
# un fichero con el nº de clusters por tissue (si existe añade el tissue a este)
# y un fichero con la anotación para ese tissue

message("\n#####################")
message("## TabulaSapiens.R ##")
message("#####################\n")

if (!require("zellkonverter", quietly = TRUE))
  BiocManager::install("zellkonverter")

suppressMessages(library(SingleCellExperiment))
suppressMessages(library(optparse))


### Command lien arguments ###

option_list <- list(
  make_option(c("-t", "--tissue"), type = "character",
              help="Tissue name")
)

opt <- parse_args(OptionParser(option_list = option_list))



### Funciones ###

#' Dataset division by cell type
#'
#' Extracts cells assigned to a cell type from a count matrix and store it in a tsv file (mainly intended for the Tabula Sapiens dataset)
#' @param cell_type The cell type we want to subset (e.g. hepatocyte)
#' @param clust_num The "numeric name" of the cluster (e.g. liver-0)
#' @param sce The SingleCellExperiment object
#' @param rm_smartseq TRUE if we don't want SmartSeq data (Tabula Sapiens have 10X and SmartSeq)
#' @param assay The assay where raw counts are stored. Default is "raw_counts" as in the Tabula Sapiens dataset
#'
#' @return It does not return anything, but it creates a tsv file with the reduced count matrix and a file with the cluster annotation
#'
#' @examples
#' sce <- zellkonverter::readH5AD("datasets/TS_Liver.h5ad")
#' divideByCellType(sce, "macrophage")
divideByCellType <- function(cell_type, clust_num, sce, rm_smartseq = TRUE, assay = "raw_counts"){
  # Count matrix reduction
  if (rm_smartseq){
    message(paste0("Removing smartseq2 data and keeping only 10X data -> ", clust_num, " (", cell_type, ")"))
    sce <- sce[, colData(sce)$method != "smartseq2"]
  }
  raw_counts <- assay(sce, assay)
  message(paste0("Reducing count matrix for ", clust_num, " (", cell_type, ")"))
  reduced_counts <- raw_counts[, rownames(colData(sce)[colData(sce)$cell_ontology_class == cell_type,])]
  dir <- paste0("datasets/", tissue)
  dir.create(dir, showWarnings = FALSE)
  write.table(as.matrix(reduced_counts), file = paste0(dir, "/TabulaSapiens-", clust_num, "-dataset.tsv"), col.names = NA, sep = "\t")
  # Annotation file creation
  message(paste0("Annotating ", clust_num, " (", cell_type, ")\n"))
  tis_clus_ann_file <- paste0("../annotations/cluster-annotation/TabulaSapiens/", tissue, "-cluster-annotation")
  line <- paste(tissue, clust_num, cell_type, sep = "\t")
  write(line, file = tis_clus_ann_file, append = file.exists(tis_clus_ann_file))
}



### Main script ###

tissue <- opt$tissue

message(paste0("Tissue: ", tissue, "\n"))

first_mayus_tissue <- paste0(toupper(substring(tissue, 1, 1)), substring(tissue, 2))
sce <- zellkonverter::readH5AD(paste0("datasets/", tissue, "/TS_", first_mayus_tissue, ".h5ad"))
message("sce object loaded!\n")

# Cell types ordered by cluster size
cell_types <- table(colData(sce)$cell_ontology_class)[order(table(colData(sce)$cell_ontology_class), decreasing = TRUE)]
names <- names(cell_types)
cell_types <- paste0(tissue, "-", 0:(length(cell_types)-1))
names(cell_types) <- names

# Adding nº of clusters to TabulaSapiens-ncluste
nclus_file <- "../annotations/TabulaSapiens-nclusters"
line <- paste(tissue, length(cell_types), sep = "\t")
write(line, file = nclus_file, append = file.exists(nclus_file))

# Applying function to each cell type

for (i in 1:(length(cell_types))){
  divideByCellType(names(cell_types[i]), unname(cell_types[i]), sce)
}

message(paste0("Nº of clusters for ", tissue, " saved in annotations/TabulaSapiens-nclusters"))
message(paste0("Cluster annotations for ", tissue, " saved in annotations/cluster-annotation/TabulaSapiens/", tissue, "-cluster-annotation"))
message(paste0("Count matrices for ", tissue, " saved in TabulaSapiens/datasets/", tissue))

message("Finished! :D")










