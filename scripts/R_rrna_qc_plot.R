#!/usr/bin/env Rscript
# R_rrna_qc_plot.R — Generate interactive stacked barplot of rRNA mapping rates
# Usage: Rscript R_rrna_qc_plot.R <output_dir>
# Expects <output_dir>/rRNA_mapping_summary_R1.tsv to exist

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

summary_file <- file.path(output_dir, "rRNA_mapping_summary_R1.tsv")

if (!file.exists(summary_file)) {
  cat("ERROR: Summary file not found:", summary_file, "\n")
  quit(status = 1)
}

suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmlwidgets))

rRnaData <- read.delim(summary_file, stringsAsFactors = FALSE)

p <- plot_ly(rRnaData,
             x = ~Sample,
             y = ~Sense_Pct,
             type = "bar",
             name = "sense") %>%
  add_trace(y = ~Antisense_Pct, name = "antisense") %>%
  layout(title = "rRNA Silva Mapping - Bowtie2, end-to-end Alignment",
         yaxis = list(title = "rRNA-Mapping-Rate in %"),
         xaxis = list(title = "", tickangle = -45),
         barmode = 'stack',
         margin = list(b = 100))

output_html <- file.path(output_dir, "rRNA_mapping_barplot.html")
saveWidget(p, file = output_html, selfcontained = TRUE)

