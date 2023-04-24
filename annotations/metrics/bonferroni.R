#!/usr/bin/env Rscript

# Sergio Alías, 20230324
# Last modified 20230424

# Script for applying Bonferroni correction to p-values

library(data.table)

urales_home <- "/run/user/1000/gvfs/sftp:host=urales/home/salias"

# read the tsv file as a data.table
dt <- fread(file.path(urales_home,
                      "filename.tsv"))

# specify the column(s) to apply Bonferroni correction
cols_to_correct <- c("column1", "column2")

# apply Bonferroni correction on the specified columns
dt[, (cols_to_correct) := lapply(.SD, function(x) pmin(1, x * length(.SD))), .SDcols = cols_to_correct]

# save the corrected data.table as a tsv with the same name but adding "_corrected" before the .tsv
write.table(dt, file = gsub(".tsv", "_corrected.tsv", "filename.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)