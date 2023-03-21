# Sergio Alías, 20220622
# Last modified 20230321

library(optparse)

option_list <- list(
  make_option(c("-t", "--tissue"), type = "character",
              help="Tissue name. Example: liver"),
  make_option(c("-c", "--n_cluster"), type = "integer",
              help="Number of clusters. Example: 17"),
  make_option(c("-d", "--dataset"), type = "character",
              help="Dataset name. Example: HPA")
)

opt <- parse_args(OptionParser(option_list = option_list))

tissue <- opt$tissue
max_cluster <- opt$n_cluster - 1

dataset <- opt$dataset
setwd(paste0('./', dataset, '/results'))

message("\n###################")
message("## coex-to-tsv.R ##")
message("###################\n")

for(i in 0:max_cluster) {
  tryCatch({
      cluster <- i
      message(paste0('Generating TSV from ', tissue, '-', cluster, ' COEX matrix...\n'))
      if(opt$dataset != "HPA"){
        coex <- readRDS(paste0(tissue, '/', tissue, '-', cluster, '/', dataset, '.', tissue, '-', cluster, '.matrix.cotan.RDS'))
      }else{
        coex <- readRDS(paste0(tissue, '/', dataset, '.', tissue, '-', cluster, '.matrix.cotan.RDS'))
      }
      write.table(colnames(coex), paste0(tissue, '/genes-', tissue, '-', cluster, '.tsv'), sep = '\t')
  }, warning = function(w) {
      message(paste0(tissue, '-', cluster, ' does not seem to exist. Please make sure this makes sense'))
  }, error = function(e) {
      message(paste0(tissue, '-', cluster, ' does not seem to exist. Please make sure this makes sense'))
  })
}
