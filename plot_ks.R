#!/usr/bin/env Rscript

# Sergio Alías, 20230428
# Last modified 20230428

# plot_ks.R

# Distribution plots with associated ks p-value


# Libs

library(data.table)
library(ggplot2)


# Main script

dataset <- "HPA" # dataset name
metric <- "fcst" # coex, fc or fcst
test <- "high" # extreme, high, low
hpo_term <- "HP:0001629" # hpo code ("HP" for coex, "HP:" for fc / fcst)
hpo_code_plain <- hpo_term
tissue <- "heart" # tissue
cluster <- "0" # cluster
tissue_clus <- paste0(tissue, "-", cluster)

ann_file <- paste0('annotations/cluster-annotation/', dataset, '/', tissue, '-cluster-annotation')
annotations <- read.table(file = ann_file, sep = '\t', header = FALSE)
cluster_ann <- annotations[strtoi(cluster)+1, 3]

if (metric == "coex"){
  
# p-value del K-S y anotación HPO

ks_pval <- fread(paste0("coex-analysis/", dataset, "/", tissue, "/wil_results_", dataset, "_", tissue, ".tsv"))

hpo_name <- ks_pval[hpo == hpo_term & tissue == tissue_clus]$hpo.name 

if (test == "extreme"){
  ks_pval <- ks_pval[hpo == hpo_term & tissue == tissue_clus][,wilcoxon.pval]
} else if (test == "high"){
  ks_pval <- ks_pval[hpo == hpo_term & tissue == tissue_clus][,high_wil_pval]
} else if (test == "low"){
  ks_pval <- ks_pval[hpo == hpo_term & tissue == tissue_clus][,low_wil_pval]
}

ks_pval <- as.numeric(ks_pval)

# Carga de los gene names
gene_file <- paste0("coex-analysis/", dataset, "/", tissue, "/", hpo_term, '/', hpo_term, '_subset_genes_', tissue, '-', cluster, '.txt')
genes <- readLines(gene_file)

# Carga y subsetting de la matrix de COEX

coex <- readRDS(paste0(dataset, '/results/', tissue, '/', dataset, '.', tissue, '-', cluster, '.matrix.cotan.RDS'))

sub_coex <- coex[genes, genes] 
coex_values <- sub_coex[lower.tri(sub_coex, diag = FALSE)]
coex_df <- data.frame(Subset = rep('HPO-related genes', length(coex_values)), COEX = coex_values)

sub_background <- coex[-which(rownames(coex) %in% genes), -which(colnames(coex) %in% genes)]
back_values <- sub_background[lower.tri(sub_background, diag = FALSE)]
back_df <- data.frame(Subset = rep('Background', length(back_values)), COEX = back_values)

df <- rbind(coex_df, back_df)

if (test == "extreme"){
  df$COEX <- abs(df$COEX)
}

plot <- ggplot(df, aes(COEX, fill = Subset, colour = Subset)) + geom_density(alpha = 0.2) +
  labs(title=paste0(dataset, ' ', tissue, '-', cluster, ' (', cluster_ann, ')\n', paste0(substring(hpo_term, 1, 2), ":", substring(hpo_term, 3)), ' (', hpo_name, ')'), y = 'Frecuencia') +
  theme(plot.title = element_text(hjust=0.5),
        legend.title= element_blank())

ylim <- layer_scales(plot)$y$range$range[2]
xlim <- layer_scales(plot)$x$range$range[2]

plot <- plot +
  annotate("label", x = xlim*0.7, y = ylim, label = paste0('K-S test p-val: \n', signif(ks_pval, 5)))


} else { # fc, fcst
  
  library(jsonlite)
  
  FC <- fread("fold_change/computed_FC.tsv")
  hpo_genes <- fromJSON("annotations/HPO/hpo_to_genes_HPA.json")
  hpo_code_plain <- paste0(strsplit(hpo_term, ":")[[1]][1], strsplit(hpo_term, ":")[[1]][2])
  genes <- hpo_genes[[hpo_term]]$EnsemblID
  
  FC <- FC[!is.na(FC$log2FC)]
  FC <- FC[FC$log2FC != -Inf]
  FC <- FC[FC$log2FC != Inf]
  FC <- FC[!is.na(FC$log2FC_single_tissue)]
  FC <- FC[FC$log2FC_single_tissue != -Inf]
  FC <- FC[FC$log2FC_single_tissue != Inf]
  FC[Tissue == "adipose tissue", Tissue := gsub("adipose tissue", "adipose-tissue", Tissue)]
  FC[Tissue == "bone marrow", Tissue := gsub("bone marrow", "bone-marrow", Tissue)]
  FC[Tissue == "heart muscle", Tissue := gsub("heart muscle", "heart", Tissue)]
  FC[Tissue == "lymph node", Tissue := gsub("lymph node", "lymph-node", Tissue)]
  FC[Tissue == "prostate", Tissue := gsub("prostate", "prostate-gland", Tissue)]
  FC[Tissue == "skeletal muscle", Tissue := gsub("skeletal muscle", "skeletal-muscle-organ", Tissue)]
  FC[Tissue == "skin", Tissue := gsub("skin", "skin-of-body", Tissue)]
  FC[Tissue == "small intestine", Tissue := gsub("small intestine", "small-intestine", Tissue)]
  FC[Tissue == "pbmc", Tissue := gsub("pbmc", "blood", Tissue)]
  
  FC <- FC[Tissue == tissue & Cluster == paste0("c-", cluster)]
  
  if (metric == "fc"){
    hpo_genes <- FC[Gene %in% genes]$log2FC
    back_genes <- FC[!(Gene %in% genes)]$log2FC # los background
    x_axis_title <- "FC"
  } else if (metric == "fcst"){
    hpo_genes <- FC[Gene %in% genes]$log2FC_single_tissue 
    back_genes <- FC[!(Gene %in% genes)]$log2FC_single_tissue # otra vez los background pero de fcst
    x_axis_title <- "FC-st"
  }

  df <- data.frame(Subset = rep('HPO-related genes', length(hpo_genes)), FC = hpo_genes)
  df <- rbind(df, data.frame(Subset = rep('Background', length(back_genes)), FC = back_genes))
  
  if (test == "extreme"){
    df$FC <- abs(df$FC)
  }
  
  # p-value del K-S y anotación HPO
  
  ks_pval <- fread(paste0("fold_change/", metric, "/", tissue, "/wil_results_", metric, "_", tissue, ".tsv"))
  
  hpo_name <- ks_pval[hpo == hpo_term & tissue == paste0("c-", cluster)]$hpo.name 
  
  if (test == "extreme"){
    ks_pval <- ks_pval[hpo == hpo_term & tissue == paste0("c-", cluster)][,wilcoxon.pval]
  } else if (test == "high"){
    ks_pval <- ks_pval[hpo == hpo_term & tissue == paste0("c-", cluster)][,high_wil_pval]
  } else if (test == "low"){
    ks_pval <- ks_pval[hpo == hpo_term & tissue == paste0("c-", cluster)][,low_wil_pval]
  }
  
  ks_pval <- as.numeric(ks_pval)
  
  plot <- ggplot(df, aes(FC, fill = Subset, colour = Subset)) + geom_density(alpha = 0.2) +
    labs(title=paste0(dataset, ' ', tissue, '-', cluster, ' (', cluster_ann, ')\n', hpo_term, ' (', hpo_name, ')'),
         y = 'Frecuencia',
         x = x_axis_title) +
    theme(plot.title = element_text(hjust=0.5),
          legend.title= element_blank())
  
  ylim <- layer_scales(plot)$y$range$range[2]
  xlim <- layer_scales(plot)$x$range$range[2]
  
  plot <- plot +
    annotate("label", x = xlim*0.7, y = ylim, label = paste0('K-S test p-val: \n', signif(ks_pval, 5)))
  
}

pdf(paste0('plots/', metric, '_', test, '_', hpo_code_plain, '_', tissue, '_', cluster, '.pdf'))

plot

invisible(dev.off())

message("Plot saved! -> ", paste0(metric, '_', test, '_', hpo_code_plain, '_', tissue, '_', cluster, '.pdf'))