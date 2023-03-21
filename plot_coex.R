# Sergio Alías, 20220614
# Last modified 20230321

# plot_coex.R

# Wilcoxon test, p-values and plots of HPO-related genes in a tissue

# Call it from the bash for multiple HPO/tissue analysis (or use workflow_coex.sh directly)


suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(rstatix))
suppressMessages(library(ggpubr))
suppressMessages(library(optparse))
#library(plyr)
#library(coin)
#library(rjson)


############################
## Command line arguments ##
############################

option_list <- list(
  make_option(c("-t", "--tissue"), type = "character",
              help = "Tissue name. Example: liver"),
  make_option(c("-c", "--max_cluster"), type = "integer",
              help = "Number of clusters of the tissue. Example for liver: 17"),
  make_option(c("-o", "--hpo"), type = "character",
              help = "HPO term (without the ':'). Example: HP0006561"),
  make_option(c("-r", "--random_samples"), type = "integer", default = 1000,
              help = "Number of random samples [default = %default]"),
  make_option(c("-d", "--dataset"), type = "character",
              help = "Dataset used")
)

opt <- parse_args(OptionParser(option_list = option_list))

###############################
## Variables and input parse ##
###############################

dataset <- opt$dataset
tissue <- opt$tissue
max_cluster <- opt$max_cluster
ann_file <- paste0('annotations/cluster-annotation/', dataset, '/', tissue, '-cluster-annotation')
annotations <- read.table(file = ann_file, sep = '\t', header = FALSE)
hpo <- opt$hpo
hponamesdata <- read.table('annotations/HPO/hpo-id-to-name', header = FALSE, sep = '\t', quote = "")
hpo_name <- hponamesdata[hponamesdata[1] == hpo][2]
seed <- 10
rnd_samples <- opt$random_samples
include_diag <- FALSE

wdir <- paste0('coex-analysis/', dataset, '/', tissue)
old_wdir <- getwd()
setwd(wdir)

pval_cols <- c('hpo',
               'hpo-name',
               'tissue',
               'annotation',
               'pval',
               'wilcoxon-pval',
               'wicoxon-sign',
               'effsize',
               'effsize-sign',
               'mean_coex')

pval_results <- data.frame(matrix(ncol=length(pval_cols),nrow=0, dimnames=list(NULL, pval_cols)))


###############
## Functions ##
###############

get_pval_wil <- function(hpo, cluster){
  cluster_ann <- annotations[strtoi(cluster)+1, 3]
  
  # Carga de los gene names
  gene_file <- paste0(hpo, '/', hpo, '_subset_genes_', tissue, '-', cluster, '.txt')
  genes <- readLines(gene_file)
  
  # Carga y subsetting de la matrix de COEX
  setwd(old_wdir)
  if(dataset != "HPA"){
    coex <- readRDS(paste0(dataset, '/results/', tissue, '/', tissue, '-', cluster, '/', dataset, '.', tissue, '-', cluster, '.matrix.cotan.RDS'))
  }else{
    coex <- readRDS(paste0(dataset, '/results/', tissue, '/', dataset, '.', tissue, '-', cluster, '.matrix.cotan.RDS'))
  }
  setwd(wdir)
  sub_coex <- coex[genes, genes]  

  # Extracción de valores únicos de la matriz (es simétrica)
  coex_values <- sub_coex[lower.tri(sub_coex, diag = include_diag)]
  
  # Creación de la lista de dataframes
  coex_df <- data.frame(Subset = rep('HPOgenes', length(coex_values)), COEX = coex_values)
  
  df_list <- list('HPOgenes' = coex_df)
  
  rm(coex_df)
  
  # 1000 réplicas de genes al azar
  set.seed(seed)
  for (i in c(1:rnd_samples)) {
    rnd_genes <- sample(1:nrow(coex), length(genes))
    rnd_sub_coex <- coex[rnd_genes, rnd_genes]
    rnd_values <- rnd_sub_coex[lower.tri(rnd_sub_coex, diag = include_diag)]
    cname <- paste0('random_', i)
    to_add <- data.frame(Subset = rep(cname, length(rnd_values)), COEX = rnd_values)
    df_list[[cname]] <- to_add
    # coex_df <- rbind(coex_df, to_add)
  }
  rm(i)

  # p-value de las medias
  rnd_means <- length(df_list) %>% numeric() # Pre-allocating the vector for efficiency purposes
  for (i in 1:length(df_list)) {
    rnd_means[i] <- df_list[[i]]$COEX %>% abs() %>% mean()
  }
  rm(i)
  
  times_gt <- 0
  for (i in 2:length(rnd_means)) {
    if (rnd_means[i] > rnd_means[1]){
      times_gt <- times_gt + 1
    }
  }
  rm(i)
  pvalue <- times_gt / rnd_samples
  message(paste0('Cluster ', cluster, ' random mean was greater than the true mean ', times_gt, ' times (p-value: ', pvalue, ")"))
  
  # wilcoxon test de las distribuciones
  rnd_1000_distr <- c()
  
  for (i in 2:length(df_list)) {
    rnd_1000_distr <- c(rnd_1000_distr, df_list[[i]]$COEX %>% abs())
  }
  message(paste0('Performing wilcox_test (', hpo, ' ', tissue, '-', cluster, ')... ', i-1, ' random samples\n'))
  
  hpo_df <- df_list[[1]]
  hpo_df$COEX <- hpo_df$COEX %>% abs()
  wil_test_1000 <- hpo_df
  to_add <- data.frame(Subset = rep('Random', length(rnd_1000_distr)), COEX = rnd_1000_distr)
  wil_test_1000 <- rbind(wil_test_1000, to_add)
  
  stat.test <- wil_test_1000 %>% 
    wilcox_test(COEX ~ Subset, alternative = "greater") %>%
    add_significance()
  stat.test
  
  effsize <- wil_test_1000 %>% wilcox_effsize(COEX ~ Subset, alternative = "greater")
  
  # plot
  plot <- ggplot(wil_test_1000, aes(COEX, fill = Subset, colour = Subset)) + geom_density(alpha = 0.2) +
    labs(title=paste0(dataset, ' ', tissue, '-', cluster, ' (', cluster_ann, ') for\n', paste0(substring(hpo, 1, 2), ":", substring(hpo, 3)), ' (', hpo_name, ')'), y = 'frecuency') +
    theme(plot.title = element_text(hjust=0.5))

  ylim <- layer_scales(plot)$y$range$range[2]
  xlim <- layer_scales(plot)$x$range$range[2]
  
  plot <- plot +
    annotate("label", x = xlim*0.7, y = ylim, label = paste0('wilcox_test p-val: \n', stat.test$p, ' (', stat.test$p.signif, ')'))
  
  # output row to the dataframe
  row_to_add <- data.frame('hpo' = hpo,
                           'hpo-name' = hpo_name,
                           'tissue' = paste0(tissue, '-', cluster),
                           'annotation' = cluster_ann,
                           'pval' = pvalue,
                           'wilcoxon-pval' = stat.test$p,
                           'wicoxon-sign' = stat.test$p.signif,
                           'effsize' = unname(effsize$effsize),
                           'effsize-sign' = effsize$magnitude,
                           'mean_coex' = df_list[['HPOgenes']]$COEX %>% abs() %>% mean())
  return(list(row_to_add, plot))
}


#################
## Main script ##
#################

### Execution

message("\n#################")
message("## plot_coex.R ##")
message("#################\n")
message(paste0("Starting analysis for ", opt$hpo, " in ", dataset, ' ', opt$tissue, " (", opt$max_cluster + 1, " clusters, ", opt$random_samples, " random samples)\n"))

hpo_plots <- vector(mode='list', length = max_cluster+1)
for (j in (0:max_cluster)){
  results <- get_pval_wil(hpo, j)
   pval_results <- rbind(pval_results, results[[1]])
  hpo_plots[[j+1]] <- results[[2]]
}

rm(j)


### Output

outfile <- paste0('wil_results_', dataset, '_', tissue, '.tsv')

write.table(pval_results,
            file = outfile,
            quote = FALSE,
            sep = '\t',
            col.names = !file.exists(outfile),
            row.names = FALSE,
            append = file.exists(outfile)
            )

message(paste0("p-val and wilcoxon results saved in ", outfile))

pdf(paste0(hpo, '/wil_plots_', dataset, '_', tissue, '_', hpo, '.pdf'))

for (j in 1:length(hpo_plots)){
  print(hpo_plots[[j]])
}

invisible(dev.off())

message(paste0("Plots saved in ", paste0(hpo, '/wil_plots_', dataset, '_', tissue, '_', hpo, '.pdf')))
