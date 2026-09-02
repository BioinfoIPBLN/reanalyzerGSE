#!/usr/bin/env Rscript
### Thin orthologr wrapper: pairwise ortholog detection between two protein FASTA files.
### Usage: R_orthologs_human.R <query.faa> <subject.faa> <ortho_detection> <cores> <out.tsv> <work_dir> [eval] [sensitivity_mode]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  cat("Usage: R_orthologs_human.R <query.faa> <subject.faa> <ortho_detection> <cores> <out.tsv> <work_dir> [eval] [sensitivity_mode]\n")
  quit(save = "no", status = 2)
}
query_file <- normalizePath(args[1], mustWork = TRUE)
subject_file <- normalizePath(args[2], mustWork = TRUE)
ortho_detection <- args[3]
comp_cores <- suppressWarnings(as.integer(args[4]))
if (is.na(comp_cores) || comp_cores < 1) comp_cores <- 1
out_tsv <- args[5]
work_dir <- args[6]
eval_thr <- if (length(args) >= 7 && nzchar(args[7])) args[7] else "1E-5"
sens_mode <- if (length(args) >= 8 && nzchar(args[8])) args[8] else "fast"

valid <- c("DIAMOND_RBH", "DIAMOND_BH", "RBH", "BH")
if (!ortho_detection %in% valid) {
  cat("ERROR: ortho_detection '", ortho_detection, "' is not supported (expected one of ",
      paste(valid, collapse = ", "), ").\n", sep = "")
  quit(save = "no", status = 2)
}

if (!requireNamespace("orthologr", quietly = TRUE)) {
  cat("ERROR: the R package 'orthologr' is not installed, so ortholog detection cannot run.\n")
  cat("Install it with: remotes::install_github('drostlab/orthologr')\n")
  quit(save = "no", status = 3)
}

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
old_wd <- setwd(work_dir)
on.exit(setwd(old_wd), add = TRUE)

cat("orthologr::orthologs(ortho_detection = '", ortho_detection, "', comp_cores = ", comp_cores,
    ", eval = '", eval_thr, "')\n", sep = "")
cat("  query  : ", query_file, "\n  subject: ", subject_file, "\n", sep = "")

res <- try(orthologr::orthologs(
  query_file      = query_file,
  subject_files   = subject_file,
  seq_type        = "protein",
  format          = "fasta",
  ortho_detection = ortho_detection,
  eval            = eval_thr,
  sensitivity_mode = sens_mode,
  comp_cores      = comp_cores,
  clean_folders   = FALSE,
  quiet           = FALSE
), silent = FALSE)

if (inherits(res, "try-error") || is.null(res) || nrow(as.data.frame(res)) == 0) {
  cat("ERROR: orthologr returned no ortholog pairs.\n")
  quit(save = "no", status = 4)
}

df <- as.data.frame(res, stringsAsFactors = FALSE)
nms <- names(df)
pick <- function(patterns, fallback = NA_integer_) {
  for (p in patterns) {
    hit <- grep(p, nms, ignore.case = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  fallback
}
q_i <- pick(c("^query_id$", "^query"), 1L)
s_i <- pick(c("^subject_id$", "^subject"), 2L)
id_i <- pick(c("^perc_identity$", "perc_ident", "^pident$", "identity"))
ev_i <- pick(c("^evalue$", "e_value", "^eval$"))
bs_i <- pick(c("^bit_score$", "bitscore", "bit_sc"))

out <- data.frame(
  query_id      = as.character(df[[q_i]]),
  subject_id    = as.character(df[[s_i]]),
  perc_identity = if (!is.na(id_i)) as.character(df[[id_i]]) else "",
  evalue        = if (!is.na(ev_i)) as.character(df[[ev_i]]) else "",
  bit_score     = if (!is.na(bs_i)) as.character(df[[bs_i]]) else "",
  stringsAsFactors = FALSE
)
out <- out[nzchar(out$query_id) & nzchar(out$subject_id), , drop = FALSE]
out <- out[!duplicated(paste(out$query_id, out$subject_id, sep = "\t")), , drop = FALSE]

dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)
write.table(out, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote ", out_tsv, " (", nrow(out), " pairs; orthologr columns used: ",
    paste(nms[c(q_i, s_i)], collapse = ", "), ")\n", sep = "")
quit(save = "no", status = 0)
