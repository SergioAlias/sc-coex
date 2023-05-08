#!/usr/bin/env Rscript

# Sergio Alías, 20230508
# Last modified 20230508

# plot_ks.R

# HPO percent of co-mentions plot


# Libs

library(data.table)
library(ggplot2)
library(patchwork)

### Input

urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"
urales_home <- "/run/user/1013/gvfs/sftp:host=urales,user=salias/home/salias"

data <- fread(file.path(urales_home,
                        "TFM/results_w_comention/coment_per_HPO.tsv"))

# Main script

setorder(data, tissue, -percent)
data[, order := 1:.N]
data$tissue <- gsub("-", " ", paste0(toupper(substring(data$tissue, 1, 1)), substring(data$tissue, 2)))
tejidos <- unique(data$tissue)

data1 <- data[tissue %in% tejidos[1:8],]
data_endo <- data[tissue == "Endometrium",]
data2 <- data[tissue %in% tejidos[9:16],]
data_ln <- data[tissue == "Lymph node",]
data3 <- data[tissue %in% tejidos[17:23],]
data_pg <- data[tissue == "Prostate gland",]

plot1 <- ggplot(data1, aes(x = tissue, y = percent)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = "#99ffcc") +
         geom_boxplot(data = data_endo, width = 0.1, fill = "white", outlier.shape = NA) +
         labs(x = NULL, y = NULL) +
         theme_bw()

plot2 <- ggplot(data2, aes(x = tissue, y = percent)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = "#99ffcc") +
  geom_boxplot(data = data_ln, width = 0.1, fill = "white", outlier.shape = NA) +
  labs( x = NULL, y = NULL) +
  theme_bw()

plot3 <- ggplot(data3, aes(x = tissue, y = percent)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = "#99ffcc") +
  geom_boxplot(data = data_pg, width = 0.1, fill = "white", outlier.shape = NA) +
  labs( x = "Tejido", y = NULL) +
  theme_bw() +
  scale_x_discrete(labels=c("Prostate", "Skeletal muscle", "Skin", "Small intestine", "Spleen", "Stomach", "Testis"))


plot <- (plot1 / plot2 / plot3) + plot_annotation(title = "Porcentaje de co-mención de tipos celulares",
                                                 theme = theme(plot.title = element_text(hjust = 0.5)))


pdf('plots/percent_HPO.pdf')

plot

invisible(dev.off())

message("Plot saved! -> plots/percent_HPO.pdf")
