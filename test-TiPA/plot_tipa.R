#!/usr/bin/env Rscript

# Sergio Alías, 20230207
# Last modified 20230207

### Libs ###

library(ggplot2)
library(patchwork)


### TiPA scores exploration ###

tipa <- read.csv("~/TFM/test-TiPA/TiPA-main/tipa_scores_matrix_HPA_liver_all.csv")

forbidden <- c("cell",
               "negative",
               "positive",
               "regulation",
               "involved",
               "process",
               "pathway",
               "to",
               "of",
               "in",
               "by"
)

cols <- colnames(tipa)[-1]

for (col in cols){

tipa.c <- tipa[order(tipa[[col]], decreasing = TRUE),][c("X", col)] # sorting by TiPA score

tipa.c <- tipa.c[!is.na(tipa.c[[col]]), ]

tipa.c <- tipa.c[1:(nrow(tipa.c)*0.01), ] # keeping top 1% scores


words <- unlist(strsplit(tipa.c$X, " ", fixed = TRUE))

words <- words[!(words %in% forbidden)]

counts <- as.data.frame(table(words))

counts <- counts[order(counts$Freq, decreasing = TRUE),]

top.8.counts <- counts[c(1:8),]

plot <- ggplot(top.8.counts, aes(x = Freq, y = factor(words, levels = rev(words)))) +
        geom_bar(stat="identity", fill = "firebrick") +
        theme_bw() +
        theme(axis.title.y = element_blank()) +
        xlab("Count") +
        ggtitle(col)

if (col == cols[1]){all.plots <- plot}else{all.plots <- all.plots + plot}

}


#############################

pdf("~/TFM/test-TiPA/HPA-liver-all.pdf", width = 16, height = 10)
all.plots
dev.off()

