#!/usr/bin/env Rscript
# render_enrichment_report.R
# Usage: render_enrichment_report.R <dge_dir> <project_name> <organism> [full|minimal] [cp_universe] [folder|embedded]
#
# Renders the functional enrichment RMarkdown report into the DGE directory. The full
# report keeps its plot previews and HTML libraries in functional_enrichment_report_files/
# (report_images 'folder', the default) or embeds everything in one HTML ('embedded');
# the minimal report (tables only) is always one self-contained HTML. Will fail loudly if
# rmarkdown or pandoc is not available.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: render_enrichment_report.R <dge_dir> <project_name> <organism> [full|minimal] [cp_universe] [folder|embedded]")
}

dge_dir      <- args[1]
project_name <- args[2]
organism     <- args[3]
mode         <- if (length(args) >= 4) args[4] else "full"
cp_universe  <- if (length(args) >= 5) args[5] else ""
report_images <- if (length(args) >= 6 && nzchar(args[6])) args[6] else "folder"
if (!report_images %in% c("folder", "embedded"))
  stop(paste0("ERROR: unknown report_images value '", report_images, "', use folder or embedded"))
if (!mode %in% c("full", "minimal"))
  stop(paste0("ERROR: unknown report mode '", mode, "', use full or minimal"))

# Locate the template (same directory as this script)
script_dir <- dirname(sub("^--file=", "", commandArgs()[grep("--file=", commandArgs())]))
if (length(script_dir) == 0 || script_dir == "") script_dir <- "."
template <- file.path(script_dir, "R_enrichment_report.Rmd")

if (!file.exists(template))
  stop(paste0("ERROR: RMarkdown template not found at: ", template))

# Copy template to the final output directory to prevent read-only filesystem errors
tmp_template <- file.path(dge_dir, basename(template))
file.copy(template, tmp_template, overwrite = TRUE)

output_file <- file.path(dge_dir, if (mode == "minimal") "functional_enrichment_report_minimal.html" else "functional_enrichment_report.html")
files_dir <- "functional_enrichment_report_files"
img_dir <- if (mode == "full" && report_images == "folder") file.path(normalizePath(dge_dir), files_dir, "images") else ""
if (nzchar(img_dir)) {
  unlink(dirname(img_dir), recursive = TRUE)
  dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("DGE dir:  %s\n  Project:  %s\n  Organism: %s\n  Output:   %s\n  Images:   %s\n",
            dge_dir, project_name, organism, output_file, if (nzchar(img_dir)) img_dir else "embedded in the HTML"))

rmarkdown::render(
  input       = tmp_template,
  output_file = output_file,
  params      = list(
    dge_dir      = normalizePath(dge_dir),
    project_name = project_name,
    organism     = organism,
    minimal      = mode == "minimal",
    cp_universe  = cp_universe,
    img_dir      = img_dir
  ),
  output_options = if (nzchar(img_dir)) list(self_contained = FALSE, lib_dir = files_dir) else NULL,
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

cat(sprintf("\nDone! Report saved to: %s\n", output_file))
