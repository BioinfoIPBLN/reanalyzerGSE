#!/usr/bin/env Rscript
# render_enrichment_report.R
# Usage: render_enrichment_report.R <dge_dir> <project_name> <organism>
#
# Renders the functional enrichment RMarkdown report into a self-contained HTML
# file inside the DGE directory. Will fail loudly if rmarkdown or pandoc is
# not available.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: render_enrichment_report.R <dge_dir> <project_name> <organism>")
}

dge_dir      <- args[1]
project_name <- args[2]
organism     <- args[3]

# Fail loudly if rmarkdown is not installed
if (!requireNamespace("rmarkdown", quietly = TRUE))
  stop("ERROR: 'rmarkdown' package is not installed. Install it with: install.packages('rmarkdown')")

# Fail loudly if pandoc is not available
if (!rmarkdown::pandoc_available())
  stop("ERROR: pandoc is not available. Install pandoc or ensure it is on your PATH.")

# Locate the template (same directory as this script)
script_dir <- dirname(sub("^--file=", "", commandArgs()[grep("--file=", commandArgs())]))
if (length(script_dir) == 0 || script_dir == "") script_dir <- "."
template <- file.path(script_dir, "R_enrichment_report.Rmd")

if (!file.exists(template))
  stop(paste0("ERROR: RMarkdown template not found at: ", template))

output_file <- file.path(dge_dir, "functional_enrichment_report.html")

cat(sprintf("Rendering functional enrichment report...\n  DGE dir:  %s\n  Project:  %s\n  Organism: %s\n  Output:   %s\n",
            dge_dir, project_name, organism, output_file))

rmarkdown::render(
  input       = template,
  output_file = output_file,
  params      = list(
    dge_dir      = normalizePath(dge_dir),
    project_name = project_name,
    organism     = organism
  ),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

cat(sprintf("\nDone! Report saved to: %s\n", output_file))
