#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
organism <- args[2]
cores <- as.numeric(args[3])
padjustmethod <- args[4]
file <- args[5]
tables <- if (length(args) >= 6 && nzchar(args[6])) args[6] else "validated"
minGS <- if (length(args) >= 7 && nzchar(args[7])) as.numeric(args[7]) else 10
maxGS <- if (length(args) >= 8 && nzchar(args[8])) as.numeric(args[8]) else 500

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "R_gene_id_helpers.R"))
.org_info <- rgse_organism_info(organism)
organism_cp <- if (is.null(.org_info)) gsub("_", " ", organism) else .org_info$name
org <- if (is.null(.org_info)) NA_character_ else .org_info$multimir
if (is.na(org)) {
  cat(paste0("multiMiR holds miRNA targets for human, mouse and rat only, so no target enrichment is done for ", organism_cp, "\n"))
  quit(save = "no", status = 0)
}
.org_status <- rgse_orgdb_check(.org_info)
if (.org_status != "ok") {
  cat(rgse_organism_note(organism, .org_info, .org_status), "\n")
  quit(save = "no", status = 0)
}
suppressMessages({
  library(clusterProfiler, quiet = TRUE, warn.conflicts = FALSE)
  library(enrichplot, quiet = TRUE, warn.conflicts = FALSE)
  library(ggplot2, quiet = TRUE, warn.conflicts = FALSE)
  library(multiMiR, quiet = TRUE, warn.conflicts = FALSE)
  library(data.table, quiet = TRUE, warn.conflicts = FALSE)
})
orgDB <- rgse_orgdb_attach(.org_info)
react_org <- .org_info$reactome

a <- read.table(file.path(path, file), header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
fc_col <- grep("^logFC", colnames(a), value = TRUE)[1]
lfc <- suppressWarnings(as.numeric(a[[fc_col]]))
fdr <- suppressWarnings(as.numeric(a$FDR))
ids <- as.character(a$Gene_ID)
ids_q <- ifelse(grepl("^[a-z]{3}-", ids), ids, paste0(org, "-", ids))
sets <- list(fdr_05 = fdr < 0.05, fdr_01 = fdr < 0.01,
             fdr_05_logpos = fdr < 0.05 & lfc > 0, fdr_05_logneg = fdr < 0.05 & lfc < 0,
             fdr_01_logpos = fdr < 0.01 & lfc > 0, fdr_01_logneg = fdr < 0.01 & lfc < 0)
sets <- lapply(sets, function(k) unique(ids_q[which(k)]))
sets <- sets[lengths(sets) > 0]
name_internal <- sub("\\.txt$", "", file)
path2 <- file.path(path, paste0(name_internal, "_funct_enrich_clusterProfiler"))
for (d in c("go_figs", "kegg_paths_snapshots", "reactome_paths_snapshots")) dir.create(file.path(path2, d), recursive = TRUE, showWarnings = FALSE)
setwd(path2)
cat(paste0("Processing ", file, ": ", length(ids), " miRNAs tested, DE sets: ", paste(names(sets), lengths(sets), sep = "=", collapse = ", "), "\n"))

fetch_targets <- function(mirnas) {
  mirnas <- unique(mirnas[!is.na(mirnas) & nzchar(mirnas)])
  if (length(mirnas) == 0) return(NULL)
  out <- list()
  for (chunk in split(mirnas, ceiling(seq_along(mirnas) / 150))) {
    res <- NULL
    for (attempt in 1:3) {
      res <- tryCatch(suppressMessages(multiMiR::get_multimir(org = org, mirna = chunk, table = tables, summary = FALSE,
                                                             predicted.cutoff.type = "p", predicted.cutoff = 20, use.tibble = FALSE)),
                      error = function(e) { message("multiMiR query failed (attempt ", attempt, "): ", conditionMessage(e)); NULL })
      if (!is.null(res)) break
      Sys.sleep(20 * attempt)
    }
    if (!is.null(res) && nrow(res@data) > 0) out[[length(out) + 1]] <- as.data.frame(res@data, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(NULL)
  d <- as.data.frame(data.table::rbindlist(out, fill = TRUE))
  d <- d[!is.na(d$target_symbol) & nzchar(d$target_symbol), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  keep <- intersect(c("mature_mirna_id", "mature_mirna_acc", "target_symbol", "target_entrez", "target_ensembl", "database", "type", "support_type", "experiment", "pubmed_id", "score"), colnames(d))
  unique(d[, keep, drop = FALSE])
}

universe_file <- file.path(path, paste0("mirna_universe_targets_", tables, ".tsv"))
universe <- NULL
if (file.exists(universe_file)) universe <- tryCatch(as.data.frame(data.table::fread(universe_file)), error = function(e) NULL)
if (is.null(universe) || nrow(universe) == 0) {
  cat(paste0("Collecting the ", tables, " targets of the ", length(unique(ids_q)), " tested miRNAs (background)...\n"))
  universe <- fetch_targets(ids_q)
  if (is.null(universe)) {
    writeLines("multiMiR returned no targets for the tested miRNAs", "mirna_targets_err.txt")
    stop("multiMiR returned no targets for the tested miRNAs of ", file)
  }
  tmp <- paste0(universe_file, ".tmp", Sys.getpid())
  data.table::fwrite(universe, tmp, sep = "\t")
  file.rename(tmp, universe_file)
}
uni_sym <- unique(universe$target_symbol)
uni_ent <- unique(as.character(universe$target_entrez[!is.na(universe$target_entrez) & nzchar(universe$target_entrez)]))
writeLines(c(paste0("mode=", tables, " multiMiR targets of the tested miRNAs"), paste0("genes=", length(uni_sym)), paste0("entrez=", length(uni_ent)),
             paste0("mirnas=", length(unique(ids_q)))), "universe.txt")

add_mirna_columns <- function(df, targets, gene_col = "geneID") {
  if (is.null(df) || nrow(df) == 0) return(df)
  tmap <- split(targets$mature_mirna_id, targets$target_symbol)
  genes <- strsplit(as.character(df[[gene_col]]), "/", fixed = TRUE)
  df$miRNAs <- vapply(genes, function(g) paste(sort(unique(unlist(tmap[g], use.names = FALSE))), collapse = ","), character(1))
  df$miRNA_count <- vapply(strsplit(df$miRNAs, ",", fixed = TRUE), function(x) sum(nzchar(x)), integer(1))
  df
}
save_plots <- function(obj, stem, dir) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0 || !any(obj@result$p.adjust < 0.05, na.rm = TRUE)) return(invisible(NULL))
  try(suppressMessages(ggsave(barplot(obj, showCategory = 20), filename = file.path(dir, paste0(stem, "_barplot.pdf")), width = 15, height = 15)), silent = TRUE)
  try(suppressMessages(ggsave(dotplot(obj, showCategory = 20), filename = file.path(dir, paste0(stem, "_dotplot.pdf")), width = 15, height = 15)), silent = TRUE)
}
i <- padjustmethod
for (geneset in names(sets)) {
  cat(paste0("\n", geneset, ": ", length(sets[[geneset]]), " DE miRNAs\n"))
  targets <- universe[universe$mature_mirna_id %in% sets[[geneset]], , drop = FALSE]
  if (nrow(targets) == 0) {
    writeLines(paste0("none of the ", length(sets[[geneset]]), " DE miRNAs has ", tables, " targets in multiMiR"), paste0("mirna_targets_", geneset, "_err.txt"))
    next
  }
  data.table::fwrite(targets[order(targets$mature_mirna_id, targets$target_symbol), ], paste0("mirna_targets_", geneset, ".tsv"), sep = "\t")
  genes_sym <- unique(targets$target_symbol)
  genes_ent <- unique(as.character(targets$target_entrez[!is.na(targets$target_entrez) & nzchar(targets$target_entrez)]))
  cat(paste0("  ", nrow(targets), " miRNA-target pairs, ", length(genes_sym), " target genes (background: ", length(uni_sym), " targets of all tested miRNAs)\n"))
  if (length(genes_sym) < 5) {
    writeLines(paste0("only ", length(genes_sym), " target genes, too few for an over-representation test"), paste0("mirna_targets_", geneset, "_err.txt"))
    next
  }
  for (ont in c("BP", "MF", "CC")) {
    tryCatch({
      ego <- suppressMessages(enrichGO(gene = genes_sym, universe = uni_sym, OrgDb = orgDB, keyType = rgse_symbol_keytype(.org_info), ont = ont,
                                       pAdjustMethod = i, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = minGS, maxGSSize = maxGS, readable = TRUE))
      df <- add_mirna_columns(as.data.frame(ego), targets)
      write.table(df, file = paste0("GO_overrepresentation_test_", ont, "_", i, "_", geneset, ".txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
      save_plots(ego, paste0("GO_overrepresentation_test_", ont, "_", i, "_", geneset), "go_figs")
    }, error = function(e) writeLines(as.character(e), paste0("GO_overrepresentation_test_", ont, "_", i, "_", geneset, "_err.txt")))
  }
  tryCatch({
    kk <- suppressMessages(enrichKEGG(gene = genes_ent, organism = org, universe = uni_ent, pAdjustMethod = i, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = minGS, maxGSSize = maxGS))
    kk <- suppressMessages(setReadable(kk, OrgDb = orgDB, keyType = "ENTREZID"))
    df <- as.data.frame(kk)
    df$Gene_ID <- df$geneID
    df <- add_mirna_columns(df, targets)
    write.table(df, file = paste0("KEGG_enrich_", i, "_", geneset, ".txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    save_plots(kk, paste0("KEGG_overrepresentation_test", i, "_", geneset), "kegg_paths_snapshots")
  }, error = function(e) writeLines(as.character(e), paste0("KEGG_enrich_", i, "_", geneset, "_err.txt")))
  tryCatch({
    rp <- suppressMessages(ReactomePA::enrichPathway(gene = genes_ent, organism = react_org, universe = uni_ent, pAdjustMethod = i, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = minGS, maxGSSize = maxGS, readable = TRUE))
    df <- as.data.frame(rp)
    df$Gene_ID <- df$geneID
    df <- add_mirna_columns(df, targets)
    write.table(df, file = paste0("REACT_", i, "_", geneset, ".txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    save_plots(rp, paste0("REACT_", i, "_", geneset), "reactome_paths_snapshots")
  }, error = function(e) writeLines(as.character(e), paste0("REACT_", i, "_", geneset, "_err.txt")))
}
writeLines(c("Functional enrichment of a small RNA-seq comparison.",
             paste0("The DE miRNAs of each set (fdr_05, fdr_05_logpos, ...) were mapped to their ", tables, " target genes with multiMiR (mirna_targets_<set>.tsv),"),
             "and the target genes were tested for GO, KEGG and Reactome over-representation against the targets of every miRNA tested in the comparison (universe.txt).",
             "The miRNAs and miRNA_count columns of each table list the DE miRNAs whose targets fall in the term."), "README_mirna_targets.txt")
system(paste0("find ", shQuote(path2), " -type f -exec awk 'NR>1{exit 1}' {} \\; -exec rm -f {} \\;"))
cat(paste0("\nDONE with ", file, ". Current date: ", date(), "\n"))
