#!/usr/bin/env Rscript
.rgse_scripts_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.rgse_scripts_dir, "R_qs_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
leave_out_arg <- trimws(args[2])
cores <- suppressWarnings(as.integer(args[3]))
if (is.na(cores) || cores < 1) cores <- 1L

suppressMessages(library("edgeR", quiet = T, warn.conflicts = F))
suppressMessages(library("parallel", quiet = T, warn.conflicts = F))

cat(paste0("\n\nLeave-out differential expression runs (", leave_out_arg, ")...\n")); print(paste0("Current date: ", date()))

dge_dir <- file.path(output_dir, "DGE")
lo_dir <- file.path(dge_dir, "leave_out")
unlink(lo_dir, recursive = TRUE)
dir.create(lo_dir, recursive = TRUE, showWarnings = FALSE)
info_file <- file.path(lo_dir, "leave_out_info.txt")
finish <- function(status, message, samples = character(0), runs = 0, code = 0) {
  writeLines(c(paste0("status=", status), paste0("mode=", if (leave_out_arg == "auto") "auto" else "samples"),
               paste0("samples=", paste(samples, collapse = ",")), paste0("runs=", runs),
               paste0("message=", message)), info_file)
  cat(paste0("\n", message, "\n"))
  quit(save = "no", status = code)
}

comp_ids <- sort(as.integer(sub("^DGE_analysis_comp([0-9]+)\\.qs2$", "\\1",
                               list.files(dge_dir, pattern = "^DGE_analysis_comp[0-9]+\\.qs2$"))))
if (length(comp_ids) == 0) finish("none", "No differential expression results were found, so no leave-out runs were made.")

gene_counts <- NULL
comps <- list()
for (k in comp_ids) {
  e <- new.env()
  rgse_load_env(file.path(dge_dir, paste0("DGE_analysis_comp", k, ".qs2")), envir = e)
  if (is.null(e$de_count_cols) || is.null(e$de_keep)) {
    finish("error", paste0("DGE_analysis_comp", k, ".qs2 was written by an older version of the pipeline; re-run step 4 to use the leave-out runs."), code = 1)
  }
  if (is.null(gene_counts)) {
    gene_counts <- e$gene_counts
    filter <- e$filter
    DESeq2_compare <- e$DESeq2_compare
  }
  prim <- read.delim(file.path(dge_dir, paste0("DGE_analysis_comp", k, ".txt")), check.names = FALSE, stringsAsFactors = FALSE)
  lfc_name <- colnames(prim)[3]
  contrast <- paste(sapply(strsplit(sub("^logFC", "", lfc_name), "__VS__", fixed = TRUE)[[1]], function(p) sub("^__", "", p)), collapse = " vs ")
  comps[[length(comps) + 1]] <- list(
    id = k, comp = e$list_combinations[[e$i]], samples = e$de_count_cols[e$de_keep], condition = as.character(e$condition),
    time = if (!is.null(e$covariab) && e$covariab != "none") e$Time else NULL, covariab_format = e$covariab_format,
    diff_soft = e$diff_soft, filter_option = e$filter_option, contrast = contrast,
    lfc = setNames(prim[[3]], prim$Gene_ID), fdr = setNames(prim$FDR, prim$Gene_ID))
  rm(e); invisible(gc())
}
if (any(sapply(comps, function(cp) cp$diff_soft == "DESeq2"))) suppressMessages(library("DESeq2", quiet = T, warn.conflicts = F))
gene_cols <- c(grep("Gene_ID", colnames(gene_counts)), grep("Length", colnames(gene_counts)))

all_samples <- unique(unlist(lapply(comps, `[[`, "samples")))
all_ids <- rgse_sample_id(all_samples)
if (leave_out_arg == "auto") {
  so_file <- file.path(output_dir, "QC_and_others", "sample_outliers", "sample_outlier_summary_norm.tsv")
  if (!file.exists(so_file)) finish("error", "leave_out_samples is 'auto' but the sample outlier check table (QC_and_others/sample_outliers/sample_outlier_summary_norm.tsv) was not found.", code = 1)
  so <- read.delim(so_file, check.names = FALSE, stringsAsFactors = FALSE)
  flagged <- so$Sample[so$Flag %in% c("outlier", "check")]
  if (length(flagged) == 0) finish("none", "leave_out_samples is 'auto' and the sample outlier check flagged no sample, so no leave-out runs were made.")
  requested <- all_samples[all_ids %in% flagged | all_samples %in% flagged]
  if (length(requested) == 0) finish("none", paste0("The sample outlier check flagged ", paste(flagged, collapse = ", "), ", but none of them is in the differential expression analyses, so no leave-out runs were made."))
} else {
  wanted <- unique(trimws(unlist(strsplit(leave_out_arg, ","))))
  wanted <- wanted[wanted != ""]
  hit <- sapply(wanted, function(w) { j <- which(all_samples == w | all_ids == w); if (length(j) == 1) j else NA_integer_ })
  if (anyNA(hit)) {
    finish("error", paste0("leave_out_samples: not found among the samples of the differential expression analyses: ",
                           paste(wanted[is.na(hit)], collapse = ", "), ". Samples removed with pattern_to_remove cannot be left out. Available samples: ",
                           paste(all_ids, collapse = ", ")), code = 1)
  }
  requested <- all_samples[unique(hit)]
}
requested_ids <- rgse_sample_id(requested)
n_runs <- 2^length(requested) - 1
if (n_runs > 63) {
  finish("error", paste0("leave_out_samples: ", length(requested), " samples give ", n_runs, " runs; the limit is 63 runs (6 samples). Samples: ",
                         paste(requested_ids, collapse = ", ")), samples = requested_ids, code = 1)
}
subsets <- unlist(lapply(seq_along(requested), function(m) combn(seq_along(requested), m, simplify = FALSE)), recursive = FALSE)
cat(paste0("Leaving out every combination of: ", paste(requested_ids, collapse = ", "), " (", length(subsets), " runs x ", length(comps), " comparisons, ", min(cores, length(subsets)), " in parallel)\n"))

run_comp <- function(cp, drop) {
  keep <- !(cp$samples %in% drop)
  cond <- cp$condition[keep]
  grp <- if (all(!startsWith(cond, "__"))) paste0("__", cond) else cond
  n <- c(sum(grp == cp$comp[1]), sum(grp == cp$comp[2]))
  if (any(n < 2)) return(list(status = paste0("skipped (n<2 in ", sub("^__", "", cp$comp[which(n < 2)[1]]), ")"), n = n))
  obj <- DGEList(counts = gene_counts[, cp$samples[keep]], group = cond, genes = gene_counts[, gene_cols])
  obj <- filter(filter = cp$filter_option, obj)
  obj <- normLibSizes(obj)
  if (is.null(cp$time)) {
    obj <- estimateDisp(obj, robust = TRUE)
    if (is.na(obj$common.dispersion)) obj$common.dispersion <- 0.4^2
  } else {
    Treat <- obj$samples$group
    Time <- cp$time[keep]
    if (is.factor(Time)) Time <- droplevels(Time)
    design <- model.matrix(~0 + Treat + Time)
    rownames(design) <- colnames(obj)
    obj <- estimateDisp(obj, design, robust = TRUE)
  }
  if (all(!startsWith(as.character(obj$samples$group), "__"))) obj$samples$group <- as.factor(paste0("__", obj$samples$group))
  if (cp$diff_soft == "DESeq2") {
    tab <- DESeq2_compare(comp = cp$comp, object = obj, covariab = if (is.null(cp$time)) "none" else paste(as.character(Time), collapse = ","),
                          covariab_format = cp$covariab_format)$table
  } else if (is.null(cp$time)) {
    tab <- topTags(exactTest(obj, pair = rev(cp$comp)), n = nrow(obj), adjust.method = "BH", sort.by = "PValue")$table
  } else {
    fit <- glmQLFit(obj, design, robust = TRUE)
    contrast <- rep(0, ncol(design))
    idx <- rev(match(paste0("Treat", sub("__", "", cp$comp)), colnames(design)))
    if (length(idx) != 2 || anyNA(idx)) stop("the contrast could not be matched to the design")
    contrast[idx] <- c(-1, 1)
    tab <- topTags(glmQLFTest(fit, contrast = contrast), n = nrow(obj), adjust.method = "BH", sort.by = "PValue")$table
  }
  list(status = "ok", n = n, lfc = setNames(tab[[3]], tab$Gene_ID), fdr = setNames(tab$FDR, tab$Gene_ID))
}

safe_cor <- function(a, b) if (length(a) >= 3 && sd(a) > 0 && sd(b) > 0) round(cor(a, b), 4) else NA_real_

run_subset <- function(s) {
  drop <- requested[s]
  run_name <- paste0("without_", paste(gsub("[^A-Za-z0-9._-]", "_", requested_ids[s]), collapse = "+"))
  run_dir <- file.path(lo_dir, run_name)
  dir.create(run_dir, showWarnings = FALSE)
  rows <- list(); keep_fdr <- list()
  for (cp in comps) {
    res <- tryCatch({ out <- NULL; invisible(capture.output(out <- run_comp(cp, drop))); out },
                    error = function(e) list(status = paste0("failed: ", conditionMessage(e)), n = c(NA, NA)))
    grp_all <- if (all(!startsWith(cp$condition, "__"))) paste0("__", cp$condition) else cp$condition
    in_comp <- any(grp_all[cp$samples %in% drop] %in% cp$comp)
    row <- data.frame(Run = run_name, Left_out = paste(requested_ids[s], collapse = ","), n_left_out = length(s),
                      Comparison = paste0("comp", cp$id), Contrast = cp$contrast, In_comparison = if (in_comp) "yes" else "no",
                      Samples = paste(res$n, collapse = " / "), Status = res$status,
                      DEGs = NA, Up = NA, Down = NA, Main_DEGs = NA, Kept = NA, Lost = NA, Gained = NA, Sign_flips = NA,
                      r_logFC_all = NA, r_logFC_main_DEGs = NA, stringsAsFactors = FALSE)
    if (res$status == "ok") {
      deg_p <- names(cp$fdr)[!is.na(cp$fdr) & cp$fdr < 0.05]
      deg_n <- names(res$fdr)[!is.na(res$fdr) & res$fdr < 0.05]
      lost <- setdiff(deg_p, deg_n); gained <- setdiff(deg_n, deg_p)
      common <- intersect(names(cp$lfc), names(res$lfc))
      common_deg <- intersect(deg_p, names(res$lfc))
      row$DEGs <- length(deg_n); row$Up <- sum(res$lfc[deg_n] > 0); row$Down <- sum(res$lfc[deg_n] < 0)
      row$Main_DEGs <- length(deg_p); row$Kept <- length(intersect(deg_p, deg_n))
      row$Lost <- length(lost); row$Gained <- length(gained)
      row$Sign_flips <- sum(sign(cp$lfc[common_deg]) != sign(res$lfc[common_deg]))
      row$r_logFC_all <- safe_cor(cp$lfc[common], res$lfc[common])
      row$r_logFC_main_DEGs <- safe_cor(cp$lfc[common_deg], res$lfc[common_deg])
      if (length(lost) + length(gained) > 0) {
        genes <- c(lost, gained)
        ch <- data.frame(Gene_ID = genes, Change = rep(c("lost", "gained"), c(length(lost), length(gained))),
                         logFC_main = round(unname(cp$lfc[genes]), 4), FDR_main = signif(unname(cp$fdr[genes]), 4),
                         logFC_leave_out = round(unname(res$lfc[genes]), 4), FDR_leave_out = signif(unname(res$fdr[genes]), 4))
        ch <- ch[order(ch$Change != "lost", ifelse(ch$Change == "lost", ch$FDR_main, ch$FDR_leave_out)), ]
        write.table(ch, file = file.path(run_dir, paste0("comp", cp$id, "_changes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      }
      keep_fdr[[as.character(cp$id)]] <- unname(res$fdr[deg_p])
    }
    rows[[length(rows) + 1]] <- row
  }
  if (length(list.files(run_dir)) == 0) unlink(run_dir, recursive = TRUE)
  cat(paste0("  ", run_name, ": done\n"))
  list(rows = do.call(rbind, rows), keep_fdr = keep_fdr)
}

results <- mclapply(subsets, function(s) tryCatch(run_subset(s), error = function(e) conditionMessage(e)),
                    mc.cores = min(cores, length(subsets)), mc.preschedule = FALSE)
failed <- !sapply(results, is.list)
if (any(failed)) cat(paste0("\nWARNING: ", sum(failed), " leave-out run(s) failed: ", paste(unique(unlist(results[failed])), collapse = " | "), "\n"))
results <- results[!failed]
if (length(results) == 0) finish("error", "Every leave-out run failed; see DGE/leave_out.log.", samples = requested_ids, runs = length(subsets), code = 1)

summary_tab <- do.call(rbind, lapply(results, `[[`, "rows"))
summary_tab <- summary_tab[order(as.integer(sub("^comp", "", summary_tab$Comparison)), summary_tab$n_left_out, summary_tab$Left_out), ]
write.table(summary_tab, file = file.path(lo_dir, "leave_out_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

stable_rows <- list()
for (cp in comps) {
  deg_p <- names(cp$fdr)[!is.na(cp$fdr) & cp$fdr < 0.05]
  runs_fdr <- lapply(results, function(r) r$keep_fdr[[as.character(cp$id)]])
  runs_fdr <- runs_fdr[!sapply(runs_fdr, is.null)]
  worst <- rep(NA_real_, length(deg_p))
  if (length(runs_fdr) > 0 && length(deg_p) > 0) {
    m <- do.call(cbind, runs_fdr)
    m[is.na(m)] <- 1
    worst <- apply(m, 1, max)
  }
  stable <- deg_p[!is.na(worst) & worst < 0.05]
  if (length(stable) > 0) {
    st <- data.frame(Gene_ID = stable, logFC_main = round(unname(cp$lfc[stable]), 4), FDR_main = signif(unname(cp$fdr[stable]), 4),
                     worst_FDR_leave_out = signif(worst[!is.na(worst) & worst < 0.05], 4))
    st <- st[order(st$FDR_main), ]
    write.table(st, file = file.path(lo_dir, paste0("stable_DEGs_comp", cp$id, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  stable_rows[[length(stable_rows) + 1]] <- data.frame(Comparison = paste0("comp", cp$id), Contrast = cp$contrast, Main_DEGs = length(deg_p),
                                                       Runs_tested = length(runs_fdr), Stable_DEGs = if (length(runs_fdr) > 0) length(stable) else NA,
                                                       stringsAsFactors = FALSE)
}
write.table(do.call(rbind, stable_rows), file = file.path(lo_dir, "stable_DEGs_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

writeLines(c(
  "Leave-out differential expression runs. For comparison only: the main differential expression results are not changed.", "",
  paste0("Samples left out, in every non-empty combination (", length(subsets), " runs): ", paste(requested_ids, collapse = ", "),
         if (leave_out_arg == "auto") " (flagged 'outlier' or 'check' by the sample outlier check)" else ""),
  "Each run repeats every comparison with the same engine, filter, normalisation, covariates and contrast as the main run, on the remaining samples.", "",
  "leave_out_summary.tsv, one row per run and comparison:",
  "- In_comparison: 'yes' if a left-out sample belongs to one of the two conditions compared; otherwise the comparison changes only through the gene filter, the normalisation and the dispersion estimate.",
  "- Samples: samples left in the first / second condition of the comparison. Status: 'ok', 'skipped (n<2 in X)' when a condition would keep fewer than 2 samples, or 'failed: <error>'.",
  "- DEGs, Up, Down: genes at FDR < 0.05 in the leave-out run. Main_DEGs: the same in the main run. Kept, Lost, Gained: main-run DEGs still DEGs, no longer DEGs (including genes removed by the filter), and new DEGs.",
  "- Sign_flips: main-run DEGs whose logFC changes sign. r_logFC_all: Pearson correlation of the logFC with the main run over all genes tested in both; r_logFC_main_DEGs: the same over the main-run DEGs, the more informative one.",
  "without_<samples>/comp<N>_changes.tsv: the lost and gained genes, with their logFC and FDR in the main and the leave-out run.",
  "stable_DEGs_summary.tsv and stable_DEGs_comp<N>.tsv: main-run DEGs that stay at FDR < 0.05 in every run where the comparison could be tested, with their worst FDR across those runs."
), file.path(lo_dir, "README_leave_out.txt"))

finish("done", paste0("Leave-out runs done: ", length(results), " of ", length(subsets), " runs, results in ", lo_dir),
       samples = requested_ids, runs = length(subsets))
