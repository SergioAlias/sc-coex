# Sergio Alías, 20221116
# Last modified 20221118
# Adapted from "Guided_tutorial.Rmd"

message("\n##############################")
message("## automatic-COTAN-script.R ##")
message("##############################\n")

suppressMessages(library(COTAN))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(ggrepel))
suppressMessages(library(factoextra))
suppressMessages(library(Rtsne))
suppressMessages(library(utils))
suppressMessages(library(plotly))
suppressMessages(library(tidyverse))
suppressMessages(library(htmlwidgets))
suppressMessages(library(MASS))
suppressMessages(library(dendextend))

#setwd("~/COTAN/results")
#setwd("~/TFM/TabulaSapiens/results")

# Colours and stuff
mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                 axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                 axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                 axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))
#--------------------#

writeToLog <-
function(msg,
         log = logfile
         )
{
  write(msg, log, append = TRUE)
}

#--------------------#

input <-
function(msg)
{
    cat(msg)
    readLines(con = "stdin", n = 1)
}

#--------------------#

tmpPlot <-
function(plot,
         outfile = "tmp.pdf"
         )
{
    pdf(outfile)
    print(plot)
    dev.off()
}


# Tissue

tissue_name <- input("Tissue name (e.g. liver): ")
n_clus <- input("Number of the last cluster (e.g. 16): ")
dataset_name <- input("Dataset (e.g. TabulaSapiens): ")

setwd(paste0("~/TFM/", dataset_name, "/results"))

logfile <- paste0(tissue_name, "/COTAN-", dataset_name, "-", tissue_name, ".log")
file.create(logfile)

writeToLog(paste0("COTAN para ", dataset_name, " ", tissue_name, ":"))

for (i in 0:n_clus){
  tryCatch({
  cluster <- i
  tissue <- paste0(tissue_name, '-', cluster)
  message("####################################")
  message(paste0("COTAN for ", dataset_name, " ", tissue))
  message("####################################")
  
  # Dataset
  
  # HPAraw <- read.table(paste('../datasets/', tissue_name, '/', dataset_name, "-", tissue, '-dataset.tsv', sep = ''), header = TRUE, sep = "\t") # PARA HPA
  HPAraw <- read.table(paste('../datasets/', tissue_name, '/', dataset_name, "-", tissue, '-dataset.tsv', sep = ''), row.names = 1, header = TRUE, sep = "\t") # PARA TABULA SAPIENS
  # HPAraw <- as.data.frame(t(as.matrix(HPAraw))) # PARA HPA
  
  ### CAMBIAR SI SE USA HPA O TABULASAPIENS ###
  # HPAraw <- HPAraw[4:nrow(HPAraw), ] # Para HPA
  # HPAraw <- HPAraw[1:10000, ] # para reducir el dataset
  
  #data("ERCCraw", package = "COTAN")
  #ERCCraw = as.data.frame(ERCCraw)
  #rownames(ERCCraw) = ERCCraw$V1
  #ERCCraw = ERCCraw[,2:ncol(ERCCraw)]
  #ERCCraw[1:5,1:5]
  
  
  # Outdir
  
  # out_path <- '/home/salias/COTAN/results/' # para HPA
  # out_path <- './'
  
  out_dir <- paste0(tissue_name, '/', tissue)
  
  
  # Initializing COTAN object with the row count table and the metadata for the experiment.
  obj <- new("scCOTAN",raw = HPAraw)
  #obj = initRaw(obj,GEO="GSM2861514" ,sc.method="Drop_seq",cond = "mouse cortex E17.5")
  obj <- initRaw(obj,GEO="several studies", sc.method="several methods", cond = paste0(dataset_name, ' ', tissue, ' dataset'))
  
  
  # Now we can start the cleaning. Analysis requires and starts
  # from a matrix of raw UMI counts after removing possible cell
  # doublets or multiplets and low quality or dying cells
  # (with too high mtRNA percentage, easily done with Seurat or
  # other tools).
  # If we do not want to consider the mitochondrial genes we
  # can remove them before starting the analysis.
  
  # genes_to_rem = get.genes(obj)[grep('^mt', get.genes(obj))] 
  # cells_to_rem = names(get.cell.size(obj)[which(get.cell.size(obj) == 0)])
  # obj = drop.genes.cells(obj,genes_to_rem,cells_to_rem )
  
  
  # We want also to define a prefix to identify the sample.
  t <- tissue
  
  print(paste("Condition ",t,sep = ""))
  #--------------------------------------
  n_cells = length(get.cell.size(object = obj))
  print(paste("n cells", n_cells, sep = " "))
  
  writeToLog(paste0("\t- ", tissue, ": ", n_cells, " células"))
  
  n_it <- 1
  
  
  ## Data cleaning
  # First we create the directory to store all information
  # regarding the data cleaning.
  if(!file.exists(tissue_name)){
    dir.create(file.path(tissue_name))
  }
  
  if(!file.exists(out_dir)){
    dir.create(file.path(out_dir))
  }
  
  if(!file.exists(paste(out_dir,"cleaning", sep = ""))){   
    dir.create(file.path(out_dir, "cleaning"))
  }
  
  
  ttm <- clean(obj)
  
  obj <- ttm$object
  
  tmpPlot(ttm$pca.cell.2)
  
  message("Please check temp.pdf and decide")
  
  while (grepl(input("Cleaning? (y/n): "), "Yy")){
    writeToLog(paste0("\t\t- Cleaning ", n_it))
    message(paste0("Performing cleaning nº ", n_it))
    pdf(file.path(out_dir,"cleaning",paste(t,"_",n_it,"_plots_before_cells_exlusion.pdf", sep = "")))
    print(ttm$pca.cell.2)
    print(ggplot(ttm$D, aes(x=n,y=means)) + geom_point() +
      geom_text_repel(data=subset(ttm$D, n > (max(ttm$D$n)- 15) ), aes(n,means,label=rownames(ttm$D[ttm$D$n > (max(ttm$D$n)- 15),])),
                      nudge_y      = 0.05,
                      nudge_x      = 0.05,
                      direction    = "x",
                      angle        = 90,
                      vjust        = 0,
                      segment.size = 0.2)+
      ggtitle("B cell group genes mean expression")+my_theme +
      theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 5,hjust = 0.02 )))
    dev.off()
    
    if (length(ttm$cl1) < length(ttm$cl2)) {
      to_rem = ttm$cl1
    }else{
      to_rem = ttm$cl2
    }
    n_it = n_it+1
    obj = drop.genes.cells(object = obj,genes = c(),cells = to_rem)
    
    gc()
    
    ttm <- clean(obj)
    #ttm = clean.sqrt(obj, cells)
    obj <- ttm$object
    
    tmpPlot(ttm$pca.cell.2)

  }
  
  # Run this only in the last iteration, instead the previous code,
  # when B cells group has not to be removed
  pdf(file.path(out_dir,"cleaning",paste(t,"_",n_it,"_plots_before_cells_exlusion.pdf", sep = "")))
  print(ttm$pca.cell.2)
  print(ggplot(ttm$D, aes(x=n,y=means)) + geom_point() +
    geom_text_repel(data=subset(ttm$D, n > (max(ttm$D$n)- 15) ), aes(n,means,label=rownames(ttm$D[ttm$D$n > (max(ttm$D$n)- 15),])),
                    nudge_y      = 0.05,
                    nudge_x      = 0.05,
                    direction    = "x",
                    angle        = 90,
                    vjust        = 0,
                    segment.size = 0.2)+
    ggtitle(label = "B cell group genes mean expression", subtitle = " - B group NOT removed -")+my_theme +
    theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 10,hjust = 0.02 ),
          plot.subtitle = element_text(color = "darkred",vjust = - 15,hjust = 0.01 )))
  
  dev.off()
  
  
  # To color the pca based on nu_j (so the cells' efficiency)
  nu_est <- round(get.nu(object = obj), digits = 7)
  
  plot.nu <- ggplot(ttm$pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))
  
  plot.nu = plot.nu + geom_point(size = 1,alpha= 0.8)+
    scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" ,
                          midpoint = log(mean(nu_est)),name = "ln(nu)")+
    ggtitle("Cells PCA coloured by cells efficiency") +
    my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                      legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                      legend.text = element_text(color = "#3C5488FF", size = 11),
                      legend.key.width = unit(2, "mm"),
                      legend.position="right")
  
  pdf(file.path(out_dir,"cleaning",paste(t,"_plots_PCA_efficiency_colored.pdf", sep = "")))
  print(plot.nu)
  dev.off()
  
  tmpPlot(plot.nu)
  
  message("Please look at the plot and decide")
  
  option <- 100
  while (option < 1 | option > 4){
    option <- as.integer(input("Opinion? \n1. No polarizado\n2. No parece muy polarizado\n3. Parece algo polarizado\n4. Bastante polarizado\nChoose number: "))
  }
  
  opinion_opt <- c("No polarizado",
                   "No parece muy polarizado",
                   "Parece algo polarizado",
                   "Bastante polarizado")
  
  writeToLog(paste0("\t\t- ", opinion_opt[option]))
  
  nu_df <- data.frame("nu"= sort(get.nu(obj)), "n"=c(1:length(get.nu(obj))))
  
  ggplot(nu_df, aes(x = n, y=nu)) + 
    geom_point(colour = "#8491B4B2", size=1)+
    my_theme #+ ylim(0,1) + xlim(0,70)
  
  
  # We can zoom on the smallest values and, if we detect
  # a clear elbow, we can decide to remove the cells.
  yset = 0.5 #threshold to remove low UDE cells
  plot.ude <- ggplot(nu_df, aes(x = n, y=nu)) + 
    geom_point(colour = "#8491B4B2", size=1) + 
    my_theme + ylim(0,1) + xlim(0,400) +
    geom_hline(yintercept=yset, linetype="dashed", color = "darkred") +
    annotate(geom="text", x=200, y=0.25, 
             label=paste("to remove cells with nu < ",yset,sep = " "), 
             color="darkred", size=4.5)
  
  tmpPlot(plot.ude)
  
  while (grepl(input("Would you like to change the threshold? (y/n): "), "Yy")){
    yset <- as.numeric(input("New threshold: "))
    plot.ude <- ggplot(nu_df, aes(x = n, y=nu)) + 
      geom_point(colour = "#8491B4B2", size=1) + 
      my_theme + ylim(0,1) + xlim(0,400) +
      geom_hline(yintercept=yset, linetype="dashed", color = "darkred") +
      annotate(geom="text", x=200, y=0.25, 
               label=paste("to remove cells with nu < ",yset,sep = " "), 
               color="darkred", size=4.5)
    tmpPlot(plot.ude)
  }
  
  if (grepl(input("Remove cells? (y/n): "), "Yy")){
    writeToLog(paste0("\t\t- Threshold a ", yset))
    pdf(file.path(out_dir,"cleaning",paste(t,"_plots_efficiency.pdf", sep = "")))
    print(plot.ude)
    dev.off()
    
    plot.ude
    # We also save the defined threshold in the metadata
    # and re run the estimation
    obj <- add.row.to.meta(obj,c("Threshold low UDE cells:",yset)) 
    
    to_rem <- rownames(nu_df[which(nu_df$nu < yset),])
    
    obj <- drop.genes.cells(object = obj, genes = c(),cells =  to_rem)
    
    
    # Repeat the estimation after the cells are removed
    ttm <- clean(obj)
    obj <- ttm$object
    ttm$pca.cell.2
  } else {
    writeToLog("\t\t- No elimino células")
    # In case we do not want to remove anything, we can run:
    pdf(file.path(out_dir,"cleaning",paste(t,"_plots_efficiency.pdf", sep = "")))
    print(ggplot(nu_df, aes(x = n, y=nu)) + geom_point(colour = "#8491B4B2", size=1) +my_theme + #xlim(0,100)+
      annotate(geom="text", x=50, y=0.25, label="nothing to remove ", color="darkred"))
    dev.off()
  }
  
  # Just to check again, we plot the final efficiency colored PCA
  nu_est <- round(get.nu(object = obj), digits = 7)
  plot.nu <- ggplot(ttm$pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))
  plot.nu = plot.nu + geom_point(size = 2,alpha= 0.8)+
    scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" ,
                          midpoint = log(mean(nu_est)),name = "ln(nu)")+
    ggtitle("Cells PCA coloured by cells efficiency: last") +
    my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                      legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                      legend.text = element_text(color = "#3C5488FF", size = 11),
                      legend.key.width = unit(2, "mm"),
                      legend.position="right")
  
  pdf(file.path(out_dir,"cleaning",paste(t,"_plots_PCA_efficiency_colored_FINAL.pdf", sep = "")))
  print(plot.nu)
  dev.off()
  
  
  tmpPlot(plot.nu)
  
  input("Take your time to check the final plot (press any): ")
  
  saveRDS(obj, file = paste('./', tissue_name, '/', tissue, '/', dataset_name, '.', tissue, '.input.cotan.RDS', sep = ''))
  
  
  rm(list=setdiff(ls(), c("mycolours",
                          "my_theme",
                          "writeToLog",
                          "input",
                          "tmpPlot",
                          "tissue_name",
                          "n_clus",
                          "dataset_name",
                          "logfile"
                         )))
  },
  error = function(e)
          {
          message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
          message(paste0("Error during ", dataset_name, " ", tissue, ": ", e))
          message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
          writeToLog("\t\t- ### Error durante el preprocessing ###")
          }
  )
}
message("All finished! :D")
