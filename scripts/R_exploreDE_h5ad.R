#!/usr/bin/env Rscript
.rgse_scripts_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.rgse_scripts_dir, "R_qs_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
final_dir <- args[1]
samples_info <- args[2]
out_file <- args[3]
organism_arg <- if (length(args) >= 4) trimws(args[4]) else ""
annotation_label <- if (length(args) >= 5) trimws(args[5]) else ""

suppressPackageStartupMessages(library(exploreDE))
suppressPackageStartupMessages(library(clusterProfiler))

cat(paste0("\n\nBuilding the exploreDE file ", out_file, "...\n")); print(paste0("Current date: ", date()))

dge_dir <- file.path(final_dir, "DGE")
read_tab <- function(f) as.data.frame(data.table::fread(f, colClasses = c(Gene_ID = "character"), check.names = FALSE), check.names = FALSE)
match_ids <- function(x, table) {
  m <- match(x, table)
  miss <- is.na(m)
  if (any(miss)) m[miss] <- match(toupper(x[miss]), toupper(table))
  m
}
as_layer <- function(df, rows, cols, label) {
  i <- match_ids(rows, df$Gene_ID)
  if (anyNA(i)) stop(sum(is.na(i)), " genes of CPM_counts_genes.txt are missing from ", label)
  if (!all(cols %in% colnames(df))) stop("samples missing from ", label, ": ", paste(setdiff(cols, colnames(df)), collapse = ", "))
  m <- as.matrix(df[i, cols, drop = FALSE])
  dimnames(m) <- list(rows, cols)
  m
}

### Counts: the genes kept by the expression filter (the ones tested), raw, TMM-normalised CPM and TPM
cpm <- read_tab(file.path(final_dir, "CPM_counts_genes.txt"))
genes <- cpm$Gene_ID
samples <- setdiff(colnames(cpm), "Gene_ID")
layers <- list(
  counts = as_layer(read_tab(file.path(final_dir, "Raw_counts_genes.txt")), genes, samples, "Raw_counts_genes.txt"),
  normalised = as_layer(cpm, genes, samples, "CPM_counts_genes.txt"),
  tpm = as_layer(read_tab(file.path(final_dir, "TPM_counts_genes.txt")), genes, samples, "TPM_counts_genes.txt")
)
cat(paste0("Genes: ", length(genes), ", samples: ", length(samples), "\n"))

### Samples
si <- read.delim(samples_info, header = FALSE, colClasses = "character")
si_ids <- NULL
for (k in intersect(c("V1", "V2"), colnames(si))) {
  cand <- rgse_sample_id(si[[k]])
  if (all(samples %in% cand)) { si_ids <- cand; break }
}
if (is.null(si_ids)) stop("the count columns (", paste(samples, collapse = ", "), ") could not be matched to ", samples_info)
cond <- si[[if (ncol(si) >= 3) "V3" else "V2"]]
coldata <- data.frame(condition = cond[match(samples, si_ids)], row.names = samples, stringsAsFactors = FALSE)
batch_file <- file.path(dirname(samples_info), "batch_vector.txt")
if (file.exists(batch_file)) {
  batch <- trimws(readLines(batch_file))
  batch <- batch[nzchar(batch)]
  if (length(batch) == nrow(si)) coldata$batch <- batch[match(samples, si_ids)]
}
coldata$condition <- factor(coldata$condition, levels = intersect(unique(cond), coldata$condition))
ord <- order(coldata$condition)
coldata <- coldata[ord, , drop = FALSE]
layers <- lapply(layers, function(m) m[, rownames(coldata), drop = FALSE])

### Gene annotation: names and descriptions from the annotation file, GO terms per ontology
rowdata <- data.frame(gene_id = genes, gene_name = genes, row.names = genes, check.names = FALSE)
merged <- list.files(dge_dir, pattern = "^DGE_analysis_comp[0-9]+_merged_CPM\\.txt$", full.names = TRUE)
if (length(merged) > 0) {
  ann <- do.call(rbind, lapply(merged, function(f) {
    d <- read_tab(f)
    nm <- intersect(c("gene_name", "gene", "Name", "gene_symbol", "symbol"), colnames(d))[1]
    ds <- intersect(c("description", "product", "gene_description", "Note", "note"), colnames(d))[1]
    data.frame(Gene_ID = d$Gene_ID,
               name = if (is.na(nm)) NA_character_ else as.character(d[[nm]]),
               description = if (is.na(ds)) NA_character_ else as.character(d[[ds]]),
               stringsAsFactors = FALSE)
  }))
  ann <- ann[!duplicated(ann$Gene_ID), , drop = FALSE]
  i <- match_ids(genes, ann$Gene_ID)
  nm <- ann$name[i]
  ok <- !is.na(nm) & nzchar(nm)
  rowdata$gene_name[ok] <- nm[ok]
  if (any(!is.na(ann$description[i]))) rowdata$description <- ann$description[i]
}

go_by_gene <- NULL
org_key <- gsub("[ ]+", "_", organism_arg)
orgdb_pkg <- unname(c(Homo_sapiens = "org.Hs.eg.db", Mus_musculus = "org.Mm.eg.db", Rattus_norvegicus = "org.Rn.eg.db")[org_key])
if (!is.na(orgdb_pkg) && requireNamespace(orgdb_pkg, quietly = TRUE)) {
  db <- getExportedValue(orgdb_pkg, orgdb_pkg)
  ens <- grepl("^ENS[A-Z]*G[0-9]+", genes)
  keytype <- if (mean(ens) > 0.5) "ENSEMBL" else "SYMBOL"
  ids <- if (keytype == "ENSEMBL") sub("\\.[0-9]+$", "", genes) else genes
  db_keys <- AnnotationDbi::keys(db, keytype = keytype)
  k <- db_keys[match_ids(ids, db_keys)]
  res <- suppressMessages(AnnotationDbi::select(db, keys = unique(k[!is.na(k)]), keytype = keytype,
                                                 columns = setdiff(c("SYMBOL", "GENENAME", "GO", "ONTOLOGY"), keytype)))
  info <- res[!duplicated(res[[keytype]]), , drop = FALSE]
  j <- match(k, info[[keytype]])
  if (keytype == "ENSEMBL") {
    sym <- info$SYMBOL[j]
    rowdata$gene_name[!is.na(sym)] <- sym[!is.na(sym)]
  }
  desc <- info$GENENAME[j]
  if (is.null(rowdata$description)) rowdata$description <- NA_character_
  rowdata$description[!is.na(desc)] <- desc[!is.na(desc)]
  res <- res[!is.na(res$GO), , drop = FALSE]
  go_by_gene <- data.frame(gene = genes[match(res[[keytype]], k)], go = res$GO, ontology = res$ONTOLOGY, stringsAsFactors = FALSE)
} else {
  go_file <- list.files(dge_dir, pattern = "\\.automatically_extracted_GO_terms\\.txt$", full.names = TRUE)[1]
  if (!is.na(go_file)) {
    g <- as.data.frame(data.table::fread(go_file, header = FALSE, fill = TRUE, colClasses = "character"))[, 1:2]
    g <- g[g[[1]] != "source_id", , drop = FALSE]
    sp <- strsplit(g[[2]], "[,; ]+")
    long <- data.frame(gene = rep(genes[match_ids(gsub(":.*", "", g[[1]]), genes)], lengths(sp)), go = trimws(unlist(sp)), stringsAsFactors = FALSE)
    long <- unique(long[!is.na(long$gene) & grepl("^GO:[0-9]{7}$", long$go), , drop = FALSE])
    onto <- if (nrow(long) > 0) tryCatch(suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(long$go), keytype = "GOID", columns = "ONTOLOGY")),
                                         error = function(e) { cat(paste0("WARNING: GO terms not split by ontology (", conditionMessage(e), "), the file has no GO columns\n")); NULL })
    if (!is.null(onto)) {
      long$ontology <- onto$ONTOLOGY[match(long$go, onto$GOID)]
      go_by_gene <- long
    }
  }
}
if (!is.null(go_by_gene)) {
  go_by_gene <- unique(go_by_gene[!is.na(go_by_gene$gene) & go_by_gene$ontology %in% c("BP", "MF", "CC"), , drop = FALSE])
  for (ont in c("BP", "MF", "CC")) {
    s <- go_by_gene[go_by_gene$ontology == ont, , drop = FALSE]
    if (nrow(s) == 0) next
    per_gene <- vapply(split(s$go, s$gene), paste, character(1), collapse = "; ")
    rowdata[[paste("GO", ont)]] <- unname(per_gene[genes])
  }
}
cat(paste0("Genes with a name: ", sum(rowdata$gene_name != rowdata$gene_id), "; with GO BP/MF/CC: ",
           paste(vapply(c("GO BP", "GO MF", "GO CC"), function(x) if (x %in% colnames(rowdata)) sum(!is.na(rowdata[[x]])) else 0L, integer(1)), collapse = "/"), "\n"))

### Differential expression and the clusterProfiler objects of each comparison
recase <- function(x, ids) {
  if (inherits(x, c("enrichResult", "gseaResult"))) {
    up <- setNames(ids, toupper(ids))
    fix <- function(v) {
      hit <- !(v %in% ids) & toupper(v) %in% names(up)
      v[hit] <- up[toupper(v[hit])]
      v
    }
    col <- if (inherits(x, "gseaResult")) "core_enrichment" else "geneID"
    if (col %in% colnames(x@result)) {
      x@result[[col]] <- vapply(strsplit(as.character(x@result[[col]]), "/", fixed = TRUE), function(g) paste(fix(g), collapse = "/"), character(1))
    }
    x@geneSets <- lapply(x@geneSets, fix)
    if (inherits(x, "gseaResult")) {
      names(x@geneList) <- fix(names(x@geneList))
    } else {
      x@gene <- fix(x@gene)
      x@universe <- fix(x@universe)
    }
    return(x)
  }
  if (is.list(x)) return(lapply(x, recase, ids = ids))
  x
}

dge_files <- list.files(dge_dir, pattern = "^DGE_analysis_comp[0-9]+\\.txt$")
dge_files <- dge_files[order(as.integer(gsub("\\D", "", dge_files)))]
de <- list(); de_columns <- list(); pathways <- list()
for (f in dge_files) {
  d <- read_tab(file.path(dge_dir, f))
  fc <- grep("^logFC", colnames(d), value = TRUE)[1]
  if (is.na(fc)) next
  grp <- strsplit(sub("^logFC", "", fc), "__VS__", fixed = TRUE)[[1]]
  contrast <- make.unique(c(names(de), paste(grp, collapse = "_vs_")))[length(de) + 1]
  i <- match_ids(d$Gene_ID, genes)
  keep <- !is.na(i) & !duplicated(i)
  d <- d[keep, , drop = FALSE]
  rownames(d) <- genes[i[keep]]
  de[[contrast]] <- d
  de_columns[[contrast]] <- list(effect_log2fc = fc, score_pvalue = "PValue", score_fdr = "FDR")
  obj <- file.path(dge_dir, "clusterProfiler_objects", sub("\\.txt$", ".qs2", f))
  if (file.exists(obj)) pathways[[contrast]] <- recase(qs2::qs_read(obj), genes)
  cat(paste0(sub("\\.txt$", "", f), " -> ", contrast, ": ", nrow(d), " genes, ",
             if (file.exists(obj)) paste(names(pathways[[contrast]]), collapse = " + ") else "no clusterProfiler objects", "\n"))
}
if (length(de) == 0) cat("No DGE_analysis_comp*.txt table found, the file holds the counts only.\n")

### Write and check
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
build_explore_h5ad(
  x = layers,
  out = out_file,
  omics = "rnaseq",
  coldata = coldata,
  rowdata = rowdata,
  rowdata_columns = list(id_col = "gene_id", symbol_col = "gene_name"),
  de = de,
  de_columns = de_columns,
  pathways = pathways,
  enrichment_padj_max = list(ora = 0.05, gsea = 0.05),
  factors = colnames(coldata),
  organism = gsub("_", " ", org_key),
  ref_build = if (nzchar(annotation_label)) paste0(org_key, "/", annotation_label) else org_key
)

eds <- new_explore_dataset(out_file)
for (cn in eds$contrasts()) {
  d <- eds$de(cn)
  src <- de[[cn]]
  if (all(is.na(d$effect_log2fc)) || (!is.null(src) && any(!is.na(src$FDR)) && all(is.na(d$score_fdr)))) {
    cat(paste0("WARNING: ", cn, ": the log2 fold change or the FDR did not come through to the exploreDE file\n"))
  }
  cat(sprintf("%s: FDR<0.05 up %d / down %d, ORA terms %d, GSEA terms %d\n", cn,
              sum(d$score_fdr < 0.05 & d$effect_log2fc > 0, na.rm = TRUE), sum(d$score_fdr < 0.05 & d$effect_log2fc < 0, na.rm = TRUE),
              NROW(eds$enrichment_table(cn, "ora")), NROW(eds$enrichment_table(cn, "gsea"))))
}
cat(paste0("Written ", out_file, " (", round(file.size(out_file) / 1e6, 1), " MB)\n"))
