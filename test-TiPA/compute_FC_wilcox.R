#!/usr/bin/env Rscript

# Sergio Alías, 20230419
# Last modified 20230420

# Script for obtaining gene FC distributions associated to HPOs and making wilcox test against random distributions

# We need:
#    - FC results -> computed_FC.tsv
#    - Tissue-HPO -> TFM/annotations/filtered_tis_hpo_ngenes_HPA.tsv
#    - HPO-ENSEMBL -> TFM/annotations/HPO/hpo_to_genes_HPA.json


### Libs ###

library(data.table) # better than data.frame
library(jsonlite) # for reading JSON files
library(tidyverse) # for the %>%
library(ggplot2)
library(rstatix)
library(ggpubr)

### Main script ###

# Input 

FC <- fread("computed_FC.tsv")
tis_hpo <- fread("../annotations/filtered_tis_hpo_ngenes_HPA.tsv")
hpo_genes <- fromJSON("../annotations/HPO/hpo_to_genes_HPA.json")

# Input filtering

colnames(tis_hpo)[1] <- "Tissue"

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


# Wilcoxon tests

seed <- 10
rnd_samples <- 1000
tis <- unique(FC$Tissue)

pval_cols <- c('hpo',
               'hpo-name',
               'tissue-name',
               'tissue',
               'annotation',
               'wilcoxon-pval',
               'wicoxon-sign',
               'high_wil_pval',
               'high_wil_sign',
               'low_wil_pval',
               'low_wil_sign')



for (i in seq_along(tis)){
  tissue <- tis[i]
  tis_FC <- FC[Tissue == tissue]
  tis_hpo_reduced <- tis_hpo[Tissue == tissue]
  hpocodes <- tis_hpo_reduced$hpo_code
  hponames <- tis_hpo_reduced$hpo_name
  message("Running for tissue: ", tissue, " (", length(hponames), " HPO terms)")
  dir.create(paste0('fc/', tissue))
  dir.create(paste0('fcst/', tissue))
  for (h in seq_along(hpocodes)){
    pval_results_fc <- data.frame(matrix(ncol=length(pval_cols),nrow=0, dimnames=list(NULL, pval_cols)))
    pval_results_fcst <- data.frame(matrix(ncol=length(pval_cols),nrow=0, dimnames=list(NULL, pval_cols)))
    hp_code <- hpocodes[h]
    hpo_code_plain <- paste0(strsplit(hp_code, ":")[[1]][1], strsplit(hp_code, ":")[[1]][2])
    dir.create(paste0('fc/', tissue, "/", hpo_code_plain))
    dir.create(paste0('fcst/', tissue, "/", hpo_code_plain))
    hp_name <- hponames[h]
    genes <- hpo_genes[[hp_code]]$EnsemblID
    message("-> ", hp_code, " - ", hp_name)
    clus <- unique(tis_FC$Cluster)
    hpo_plots_fc <- vector(mode='list', length = length(clus))
    hpo_plots_high_fc <- vector(mode='list', length = length(clus))
    hpo_plots_low_fc <- vector(mode='list', length = length(clus))
    
    hpo_plots_fcst <- vector(mode='list', length = length(clus))
    hpo_plots_high_fcst <- vector(mode='list', length = length(clus))
    hpo_plots_low_fcst <- vector(mode='list', length = length(clus))
    
    for (j in seq_along(clus)){
      cluster <- clus[j]
      cluster_name <- unique(tis_FC[Cluster == cluster, Cell_type])
      message("    - Cluster ", cluster, ": ", cluster_name, "...")
      tis_clus_FC <- tis_FC[Cluster == cluster]
      fc_hpo_genes <- tis_clus_FC[Gene %in% genes]$log2FC
      fcst_hpo_genes <- tis_clus_FC[Gene %in% genes]$log2FC_single_tissue
      rnd_sample_size <- length(fc_hpo_genes)
      # creating objects for storing distributions
      fc_df <- data.frame(Subset = rep('HPOgenes', length(fc_hpo_genes)), FC = fc_hpo_genes)
      fcst_df <- data.frame(Subset = rep('HPOgenes', length(fcst_hpo_genes)), FC = fcst_hpo_genes)
      df_list_fc <- list('HPOgenes' = fc_df)
      df_list_fcst <- list('HPOgenes' = fcst_df)
      rm(fc_df)
      rm(fcst_df)
      # sampling 1000 random FC distributions
      set.seed(seed)
      for (r in c(1:rnd_samples)) {
        rnd_genes <- sample(1:nrow(tis_clus_FC), rnd_sample_size)
        rnd_fc <- tis_clus_FC[rnd_genes, log2FC]
        rnd_fcst <- tis_clus_FC[rnd_genes, log2FC_single_tissue]
        cname <- paste0('random_', r)
        to_add_fc <- data.frame(Subset = rep(cname, length(rnd_fc)), FC = rnd_fc)
        to_add_fcst <- data.frame(Subset = rep(cname, length(rnd_fcst)), FC = rnd_fcst)
        df_list_fc[[cname]] <- to_add_fc
        df_list_fcst[[cname]] <- to_add_fcst
      }
      # wilcoxon tests
      rnd_1000_distr_fc <- c()
      rnd_1000_distr_fcst <- c()
      rnd_1000_distr_fc_abs <- c()
      rnd_1000_distr_fcst_abs <- c()
      for (i in 2:length(df_list_fc)) {
        rnd_1000_distr_fc <- c(rnd_1000_distr_fc, df_list_fc[[i]]$FC)
        rnd_1000_distr_fcst <- c(rnd_1000_distr_fcst, df_list_fcst[[i]]$FC)
        rnd_1000_distr_fc_abs <- c(rnd_1000_distr_fc_abs, df_list_fc[[i]]$FC %>% abs())
        rnd_1000_distr_fcst_abs <- c(rnd_1000_distr_fcst_abs, df_list_fcst[[i]]$FC %>% abs())
      }
      message('    - Performing 3 wilcox_tests x 2 FC values')
      
      hpo_df_fc_abs <- df_list_fc[[1]]
      hpo_df_fc_abs$FC <- hpo_df_fc_abs$FC %>% abs()
      wil_test_1000_fc_abs <- hpo_df_fc_abs
      wil_test_1000_fc <- df_list_fc[[1]]
      
      hpo_df_fcst_abs <- df_list_fcst[[1]]
      hpo_df_fcst_abs$FC <- hpo_df_fcst_abs$FC %>% abs()
      wil_test_1000_fcst_abs <- hpo_df_fcst_abs
      wil_test_1000_fcst <- df_list_fcst[[1]]
      
      to_add_fc_abs <- data.frame(Subset = rep('Random', length(rnd_1000_distr_fc_abs)), FC = rnd_1000_distr_fc_abs)
      to_add_fc <- data.frame(Subset = rep('Random', length(rnd_1000_distr_fc)), FC = rnd_1000_distr_fc)
      wil_test_1000_fc_abs <- rbind(wil_test_1000_fc_abs, to_add_fc_abs)
      wil_test_1000_fc <- rbind(wil_test_1000_fc, to_add_fc)
      
      to_add_fcst_abs <- data.frame(Subset = rep('Random', length(rnd_1000_distr_fcst_abs)), FC = rnd_1000_distr_fcst_abs)
      to_add_fcst <- data.frame(Subset = rep('Random', length(rnd_1000_distr_fcst)), FC = rnd_1000_distr_fcst)
      wil_test_1000_fcst_abs <- rbind(wil_test_1000_fcst_abs, to_add_fcst_abs)
      wil_test_1000_fcst <- rbind(wil_test_1000_fcst, to_add_fcst)
      
      stat.test.fc <- wil_test_1000_fc_abs %>% 
        wilcox_test(FC ~ Subset, alternative = "greater") %>%
        add_significance()
      stat.test.fc
      
      stat.test.high.fc <- wil_test_1000_fc %>% 
        wilcox_test(FC ~ Subset, alternative = "greater") %>%
        add_significance()
      stat.test.high.fc
      
      stat.test.low.fc <- wil_test_1000_fc %>% 
        wilcox_test(FC ~ Subset, alternative = "less") %>%
        add_significance()
      stat.test.low.fc
      
      
      
      stat.test.fcst <- wil_test_1000_fcst_abs %>% 
        wilcox_test(FC ~ Subset, alternative = "greater") %>%
        add_significance()
      stat.test.fcst
      
      stat.test.high.fcst <- wil_test_1000_fcst %>% 
        wilcox_test(FC ~ Subset, alternative = "greater") %>%
        add_significance()
      stat.test.high.fcst
      
      stat.test.low.fcst <- wil_test_1000_fcst %>% 
        wilcox_test(FC ~ Subset, alternative = "less") %>%
        add_significance()
      stat.test.low.fcst
      
      # plots
      
      plot.fc <- ggplot(wil_test_1000_fc_abs, aes(FC, fill = Subset, colour = Subset)) + geom_density(alpha = 0.2) +
        labs(title=paste0(tissue, ' ', cluster, ' (', cluster_name, ') for\n', hp_code, ' (', hp_name, ')'), y = 'frecuency') +
        theme(plot.title = element_text(hjust=0.5))
      
      ylim <- layer_scales(plot.fc)$y$range$range[2]
      xlim <- layer_scales(plot.fc)$x$range$range[2]
      
      plot.fc <- plot.fc +
        annotate("label", x = xlim*0.7, y = ylim, label = paste0('wilcox_test p-val: \n', stat.test.fc$p, ' (', stat.test.fc$p.signif, ')'))
      
      plot.no.abs.fc <- ggplot(wil_test_1000_fc, aes(FC, fill = Subset, colour = Subset)) + geom_density(alpha = 0.2) +
        labs(title=paste0(tissue, ' ', cluster, ' (', cluster_name, ') for\n', hp_code, ' (', hp_name, ')'), y = 'frecuency') +
        theme(plot.title = element_text(hjust=0.5))
      
      ylim <- layer_scales(plot.no.abs.fc)$y$range$range[2]
      xlim <- layer_scales(plot.no.abs.fc)$x$range$range[2]
      
      plot.high.fc <- plot.no.abs.fc +
        annotate("label", x = xlim*0.7, y = ylim, label = paste0('wilcox_test p-val: \n', stat.test.high.fc$p, ' (', stat.test.high.fc$p.signif, ')'))
      
      plot.low.fc <- plot.no.abs.fc +
        annotate("label", x = xlim*0.7, y = ylim, label = paste0('wilcox_test p-val: \n', stat.test.low.fc$p, ' (', stat.test.low.fc$p.signif, ')'))
      
      
      
      plot.fcst <- ggplot(wil_test_1000_fcst_abs, aes(FC, fill = Subset, colour = Subset)) + geom_density(alpha = 0.2) +
        labs(title=paste0(tissue, ' ', cluster, ' (', cluster_name, ') for\n', hp_code, ' (', hp_name, ')'), y = 'frecuency') +
        theme(plot.title = element_text(hjust=0.5))
      
      ylim <- layer_scales(plot.fcst)$y$range$range[2]
      xlim <- layer_scales(plot.fcst)$x$range$range[2]
      
      plot.fcst <- plot.fcst +
        annotate("label", x = xlim*0.7, y = ylim, label = paste0('wilcox_test p-val: \n', stat.test.fcst$p, ' (', stat.test.fcst$p.signif, ')'))
      
      plot.no.abs.fcst <- ggplot(wil_test_1000_fcst, aes(FC, fill = Subset, colour = Subset)) + geom_density(alpha = 0.2) +
        labs(title=paste0(tissue, ' ', cluster, ' (', cluster_name, ') for\n', hp_code, ' (', hp_name, ')'), y = 'frecuency') +
        theme(plot.title = element_text(hjust=0.5))
      
      ylim <- layer_scales(plot.no.abs.fcst)$y$range$range[2]
      xlim <- layer_scales(plot.no.abs.fcst)$x$range$range[2]
      
      plot.high.fcst <- plot.no.abs.fcst +
        annotate("label", x = xlim*0.7, y = ylim, label = paste0('wilcox_test p-val: \n', stat.test.high.fcst$p, ' (', stat.test.high.fcst$p.signif, ')'))
      
      plot.low.fcst <- plot.no.abs.fcst +
        annotate("label", x = xlim*0.7, y = ylim, label = paste0('wilcox_test p-val: \n', stat.test.low.fcst$p, ' (', stat.test.low.fcst$p.signif, ')'))
      
      # output
      
      row_to_add_fc <- data.frame('hpo' = hp_code,
                                  'hpo-name' = hp_name,
                                  'tissue-name' = tissue,
                                  'tissue' = cluster,
                                  'annotation' = cluster_name,
                                  'wilcoxon-pval' = stat.test.fc$p,
                                  'wicoxon-sign' = stat.test.fc$p.signif,
                                  'high_wil_pval' = stat.test.high.fc$p,
                                  'high_wil_sign' = stat.test.high.fc$p.signif,
                                  'low_wil_pval' = stat.test.low.fc$p,
                                  'low_wil_sign' = stat.test.low.fc$p.signif)
      
      row_to_add_fcst <- data.frame('hpo' = hp_code,
                                  'hpo-name' = hp_name,
                                  'tissue-name' = tissue,
                                  'tissue' = cluster,
                                  'annotation' = cluster_name,
                                  'wilcoxon-pval' = stat.test.fcst$p,
                                  'wicoxon-sign' = stat.test.fcst$p.signif,
                                  'high_wil_pval' = stat.test.high.fcst$p,
                                  'high_wil_sign' = stat.test.high.fcst$p.signif,
                                  'low_wil_pval' = stat.test.low.fcst$p,
                                  'low_wil_sign' = stat.test.low.fcst$p.signif)
      
      pval_results_fc <- rbind(pval_results_fc, row_to_add_fc)
      pval_results_fcst <- rbind(pval_results_fcst, row_to_add_fcst)
      
      hpo_plots_fc[[j]] <- plot.fc
      hpo_plots_high_fc[[j]] <- plot.high.fc
      hpo_plots_low_fc[[j]] <- plot.low.fc
      hpo_plots_fcst[[j]] <- plot.fcst
      hpo_plots_high_fcst[[j]] <- plot.high.fcst
      hpo_plots_low_fcst[[j]] <- plot.low.fcst
      
    }
  
    outfile <- paste0('fc/', tissue, '/wil_results_fc_', tissue, '.tsv')
    
    write.table(pval_results_fc,
                file = outfile,
                quote = FALSE,
                sep = '\t',
                col.names = !file.exists(outfile),
                row.names = FALSE,
                append = file.exists(outfile)
    )
    
    message(paste0("wilcoxon fc results saved in ", outfile))
    
    outfile <- paste0('fcst/', tissue, '/wil_results_fcst_', tissue, '.tsv')
    
    write.table(pval_results_fcst,
                file = outfile,
                quote = FALSE,
                sep = '\t',
                col.names = !file.exists(outfile),
                row.names = FALSE,
                append = file.exists(outfile)
    )
    
    message(paste0("wilcoxon fcst results saved in ", outfile))
    
    #### Saving plots
    
    
    pdf(paste0('fc/', tissue, "/", hpo_code_plain, '/wil_fc_plots_extreme_', tissue, '_', hpo_code_plain, '.pdf'))
    
    for (j in 1:length(hpo_plots_fc)){
      print(hpo_plots_fc[[j]])
    }
    
    invisible(dev.off())
    
    message(paste0("Plots for fc extreme values saved in ", paste0('fc/', tissue, "/", hpo_code_plain, '/wil_fc_plots_extreme_', tissue, '_', hpo_code_plain, '.pdf')))
    
    ##### High values
    
    pdf(paste0('fc/', tissue, "/", hpo_code_plain, '/wil_fc_plots_high_', tissue, '_', hpo_code_plain, '.pdf'))
    
    for (j in 1:length(hpo_plots_high_fc)){
      print(hpo_plots_high_fc[[j]])
    }
    
    invisible(dev.off())
    
    message(paste0("Plots for fc high values saved in ", paste0('fc/', tissue, "/", hpo_code_plain, '/wil_fc_plots_high_', tissue, '_', hpo_code_plain, '.pdf')))
    
    ##### Low values
    
    pdf(paste0('fc/', tissue, "/", hpo_code_plain, '/wil_fc_plots_low_', tissue, '_', hpo_code_plain, '.pdf'))
    
    for (j in 1:length(hpo_plots_low_fc)){
      print(hpo_plots_low_fc[[j]])
    }
    
    invisible(dev.off())
    
    message(paste0("Plots for fc low values saved in ", paste0('fc/', tissue, "/", hpo_code_plain, '/wil_fc_plots_low_', tissue, '_', hpo_code_plain, '.pdf')))
    
    
    #### fcst ####
    
    
    pdf(paste0('fcst/', tissue, "/", hpo_code_plain, '/wil_fcst_plots_extreme_', tissue, '_', hpo_code_plain, '.pdf'))
    
    for (j in 1:length(hpo_plots_fcst)){
      print(hpo_plots_fcst[[j]])
    }
    
    invisible(dev.off())
    
    message(paste0("Plots for fcst extreme values saved in ", paste0('fcst/', tissue, "/", hpo_code_plain, '/wil_fcst_plots_extreme_', tissue, '_', hpo_code_plain, '.pdf')))
    
    ##### High values
    
    pdf(paste0('fcst/', tissue, "/", hpo_code_plain, '/wil_fcst_plots_high_', tissue, '_', hpo_code_plain, '.pdf'))
    
    for (j in 1:length(hpo_plots_high_fcst)){
      print(hpo_plots_high_fcst[[j]])
    }
    
    invisible(dev.off())
    
    message(paste0("Plots for fcst high values saved in ", paste0('fcst/', tissue, "/", hpo_code_plain, '/wil_fcst_plots_high_', tissue, '_', hpo_code_plain, '.pdf')))
    
    ##### Low values
    
    pdf(paste0('fcst/', tissue, "/", hpo_code_plain, '/wil_fcst_plots_low_', tissue, '_', hpo_code_plain, '.pdf'))
    
    for (j in 1:length(hpo_plots_low_fcst)){
      print(hpo_plots_low_fcst[[j]])
    }
    
    invisible(dev.off())
    
    message(paste0("Plots for fcst low values saved in ", paste0('fcst/', tissue, "/", hpo_code_plain, '/wil_fcst_plots_low_', tissue, '_', hpo_code_plain, '.pdf')))
    
  
  }
}