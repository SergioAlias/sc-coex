#!/usr/bin/env Rscript

# Sergio Alías, 20230131
# Last modified 20230201

# Script for obtaining the input files for TiPA analysis

### Libs ###
library(edgeR) # edgeR works on a table of read counts, with rows corresponding to genes and columns to independent libraries
# library(limma) Already an edgeR dependency

### Functions ###

getEntrezIDs <-
  function(from = "ENSEMBL")
  {
    suppressMessages(require(org.Hs.eg.db))
    ensembldb::mapIds(org.Hs.eg.db,
                      keys = ensembldb::keys(org.Hs.eg.db,
                                             keytype = from),
                      keytype = from,
                      column = "ENTREZID")
  }

### Main script ###

entrezIDs <- getEntrezIDs()

files <- list.files("~/COTAN/datasets/liver/")
names <- paste0(unlist(strsplit(files, "-dataset.tsv")), "-all")

for (i in seq_along(files)){
  
tissue_name <- names[i]
message(paste0("Running ", tissue_name, "..."))

HPAraw <- read.table(paste0('~/COTAN/datasets/liver/', files[i]), header = TRUE, sep = "\t")
HPAraw <- as.data.frame(t(as.matrix(HPAraw)))
HPAraw <- HPAraw[4:nrow(HPAraw), ]
row.names <- rownames(HPAraw)
HPAraw <- data.frame(apply(HPAraw, 2, function(x) as.numeric(as.character(x))))
rownames(HPAraw) <- row.names

### Filtering ###

#HPAraw <- HPAraw[rowSums(HPAraw) > 10, ]

### Norm counts ###

HPAraw <- DGEList(HPAraw)

HPAraw <- calcNormFactors(HPAraw, method="TMM")

cpms <- cpm.DGEList(HPAraw, log=TRUE)

### Transform ###

HPA <- voom(cpms)

### Linear model ###

HPAlm <- lmFit(HPA)
HPAlm <- eBayes(HPAlm)

### LogFC export ### 

top.table <- topTable(HPAlm, n = Inf)

top.table <- top.table[rownames(top.table) %in% names(entrezIDs),]

output <- data.frame(gene_name = entrezIDs[rownames(top.table)],
                     gene_ENS = rownames(top.table),
                     logFC = top.table$logFC)

write.csv(output, file = paste0("~/TFM/test-TiPA/TiPA-main/tissues_FC_directory/", tissue_name, ".csv"), row.names = FALSE, quote = FALSE)

}

message("Finished! :D")