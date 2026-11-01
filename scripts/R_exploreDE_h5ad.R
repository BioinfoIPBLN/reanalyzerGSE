#!/usr/bin/env Rscript
.rgse_scripts_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.rgse_scripts_dir, "R_qs_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
final_dir <- args[1]
samples_info <- args[2]
out_prefix <- sub("\\.h5ad$", "", args[3])
organism_arg <- if (length(args) >= 4) trimws(args[4]) else ""
annotation_file <- if (length(args) >= 5) trimws(args[5]) else ""
feature_key <- if (length(args) >= 6) trimws(args[6]) else ""
annotation_label <- basename(annotation_file)

suppressPackageStartupMessages(library(exploreDE))
suppressPackageStartupMessages(library(clusterProfiler))

cat(paste0("\n\nBuilding the exploreDE file(s) ", out_prefix, "*.h5ad...\n")); print(paste0("Current date: ", date()))

dge_dir <- file.path(final_dir, "DGE")
read_tab <- function(f) as.data.frame(data.table::fread(f, colClasses = c(Gene_ID = "character"), check.names = FALSE), check.names = FALSE)
match_ids <- function(x, table) {
  m <- match(x, table)
  miss <- is.na(m)
  if (any(miss)) m[miss] <- match(toupper(x[miss]), toupper(table))
  m
}
annotation_descriptions <- function(f, key) {
  if (!nzchar(f) || !file.exists(f) || !nzchar(key)) return(NULL)
  con <- gzfile(f)
  x <- readLines(con, warn = FALSE)
  close(con)
  a <- sub("^([^\t]*\t){8}", "", x[!startsWith(x, "#") & nzchar(x)])
  if (!any(grepl("product|description|Note", a))) return(NULL)
  gff3 <- !any(grepl("\"", head(a, 1000), fixed = TRUE))
  attr1 <- function(k) {
    pat <- if (gff3) paste0("(?:^|;)\\s*", k, "=([^;]*)") else paste0("(?:^|;)\\s*", k, "\\s+\"([^\"]*)\"")
    m <- regexpr(pat, a, perl = TRUE)
    v <- rep(NA_character_, length(a))
    v[m > 0] <- sub(pat, "\\1", regmatches(a, m), perl = TRUE)
    enc <- gff3 & !is.na(v) & grepl("%", v, fixed = TRUE)
    v[enc] <- vapply(v[enc], function(z) tryCatch(suppressWarnings(utils::URLdecode(z)), error = function(e) z), character(1), USE.NAMES = FALSE)
    v
  }
  keys <- attr1(key)
  desc <- attr1("product")
  for (k in c("description", "Note")) {
    miss <- is.na(desc) | !nzchar(desc)
    if (any(miss)) desc[miss] <- attr1(k)[miss]
  }
  if (gff3) {
    ids <- attr1("ID")
    parent <- sub(",.*", "", attr1("Parent"))
    for (pass in 1:3) {
      known <- !is.na(ids) & !is.na(keys)
      id2key <- setNames(keys[known], ids[known])
      fill <- is.na(keys) & !is.na(parent) & parent %in% names(id2key)
      if (!any(fill)) break
      keys[fill] <- id2key[parent[fill]]
    }
  }
  keep <- !is.na(keys) & !is.na(desc) & nzchar(desc)
  keep[keep] <- !duplicated(keys[keep])
  setNames(desc[keep], keys[keep])
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

layers$log2 <- log2(layers$normalised + 1)

### Samples, with one factor per design file (the comparisons of each design use its groups)
si <- read.delim(samples_info, header = FALSE, colClasses = "character")
si_ids <- NULL
for (k in intersect(c("V1", "V2"), colnames(si))) {
  cand <- rgse_sample_id(si[[k]])
  if (all(samples %in% cand)) { si_ids <- cand; break }
}
if (is.null(si_ids)) stop("the count columns (", paste(samples, collapse = ", "), ") could not be matched to ", samples_info)
rows <- match(samples, si_ids)
designs <- list()
for (z in list.files(dirname(samples_info), pattern = "design_possible_full", recursive = TRUE, full.names = TRUE)) {
  v <- tryCatch(as.character(read.table(z, header = FALSE, blank.lines.skip = FALSE)$V1), error = function(e) NULL)
  if (length(v) != nrow(si) || length(unique(v[rows])) < 2 || any(vapply(designs, identical, logical(1), v))) next
  designs[[length(designs) + 1]] <- v
}
if (length(designs) == 0) designs <- list(si[[if (ncol(si) >= 3) "V3" else "V2"]])
names(designs) <- c("condition", paste0("condition_", seq_along(designs)[-1]))
coldata <- data.frame(lapply(designs, function(v) factor(v[rows], levels = intersect(unique(v), v[rows]))), row.names = samples)
batch_file <- file.path(dirname(samples_info), "batch_vector.txt")
if (file.exists(batch_file)) {
  batch <- trimws(unlist(strsplit(readLines(batch_file), ",")))
  batch <- batch[nzchar(batch)]
  if (length(batch) == nrow(si)) coldata$batch <- batch[rows]
}
cat(paste0("Sample groups: ", paste(names(designs), collapse = ", "), if ("batch" %in% colnames(coldata)) ", batch" else "", "\n"))

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
  ad <- tryCatch(annotation_descriptions(annotation_file, feature_key), error = function(e) {
    cat(paste0("WARNING: no gene descriptions read from ", annotation_file, " (", conditionMessage(e), ")\n"))
    NULL
  })
  if (length(ad) > 0) {
    i <- match_ids(genes, names(ad))
    if (is.null(rowdata$description)) rowdata$description <- NA_character_
    fill <- (is.na(rowdata$description) | !nzchar(rowdata$description)) & !is.na(i)
    rowdata$description[fill] <- unname(ad[i[fill]])
  }
  go_file <- list.files(dge_dir, pattern = "\\.automatically_extracted_GO_terms\\.txt$", full.names = TRUE)[1]
  if (!is.na(go_file)) {
    g <- as.data.frame(data.table::fread(go_file, header = FALSE, fill = TRUE, colClasses = "character"))[, 1:2]
    g <- g[g[[1]] != "source_id", , drop = FALSE]
    sid <- gsub(":.*", "", g[[1]])
    hit <- genes[match_ids(sid, genes)]
    if (all(is.na(hit))) {
      u <- unique(sid)
      inside <- vapply(u, function(x) { h <- genes[grep(x, genes, fixed = TRUE)]; if (length(h) > 0) h[1] else NA_character_ }, character(1))
      hit <- unname(inside[sid])
    }
    sp <- strsplit(g[[2]], "[,; ]+")
    long <- data.frame(gene = rep(hit, lengths(sp)), go = trimws(unlist(sp)), stringsAsFactors = FALSE)
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
    rowdata[[paste0("go_", tolower(ont))]] <- unname(per_gene[genes])
  }
}
cat(paste0("Genes with a name: ", sum(rowdata$gene_name != rowdata$gene_id), "; with a description: ",
           if (is.null(rowdata$description)) 0L else sum(!is.na(rowdata$description) & nzchar(rowdata$description)), "; with GO BP/MF/CC: ",
           paste(vapply(c("go_bp", "go_mf", "go_cc"), function(x) if (x %in% colnames(rowdata)) sum(!is.na(rowdata[[x]])) else 0L, integer(1)), collapse = "/"), "\n"))

### Differential expression and the clusterProfiler objects of each comparison
recase <- function(x, ids) {
  if (inherits(x, c("enrichResult", "gseaResult"))) {
    up <- setNames(ids, toupper(ids))
    fix <- function(v) {
      hit <- !(v %in% ids) & toupper(v) %in% names(up)
      v[hit] <- up[toupper(v[hit])]
      v
    }
    gsea <- inherits(x, "gseaResult")
    col <- if (gsea) "core_enrichment" else "geneID"
    if (col %in% colnames(x@result)) {
      x@result[[col]] <- vapply(strsplit(as.character(x@result[[col]]), "/", fixed = TRUE), function(g) paste(fix(g), collapse = "/"), character(1))
    }
    if (gsea) {
      kept <- intersect(names(x@geneSets), x@result$ID)
      x@geneSets[kept] <- lapply(x@geneSets[kept], fix)
      names(x@geneList) <- fix(names(x@geneList))
    }
    return(x)
  }
  if (is.list(x)) return(lapply(x, recase, ids = ids))
  x
}

dge_files <- list.files(dge_dir, pattern = "^DGE_analysis_comp[0-9]+\\.txt$")
dge_files <- dge_files[order(as.integer(gsub("\\D", "", dge_files)))]
groups <- lapply(designs, function(v) unique(v[rows]))
de <- list(); de_columns <- list(); pathways <- list(); design_of <- character(0)
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
  k <- names(groups)[vapply(groups, function(l) length(grp) == 2 && all(grp %in% l), logical(1))][1]
  if (is.na(k)) {
    cat(paste0("NOTE: the groups of ", contrast, " are not the groups of any design file, so it goes with the first one\n"))
    k <- names(groups)[1]
  }
  design_of[contrast] <- k
  obj <- file.path(dge_dir, "clusterProfiler_objects", sub("\\.txt$", ".qs2", f))
  if (file.exists(obj)) {
    pw <- tryCatch(recase(qs2::qs_read(obj), genes), error = function(e) {
      cat(paste0("WARNING: ", obj, " could not be read (", conditionMessage(e), "), so ", contrast, " has no pathway results\n"))
      NULL
    })
    if (!is.null(pw)) pathways[[contrast]] <- pw
  }
  cat(paste0(sub("\\.txt$", "", f), " -> ", contrast, " (", k, "): ", nrow(d), " genes, ",
             if (is.null(pathways[[contrast]])) "no clusterProfiler objects" else paste(names(pathways[[contrast]]), collapse = " + "), "\n"))
}
if (length(de) == 0) cat("No DGE_analysis_comp*.txt table found, the file holds the counts only.\n")

### Write and check, one file per design that has comparisons
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
used <- if (length(de) > 0) unique(design_of) else names(designs)[1]
for (k in used) {
  cn <- names(de)[design_of[names(de)] == k]
  cd <- coldata[, c(k, setdiff(colnames(coldata), k)), drop = FALSE]
  cd <- cd[order(cd[[k]]), , drop = FALSE]
  out_file <- if (length(used) == 1) paste0(out_prefix, ".h5ad") else paste0(out_prefix, "_", k, ".h5ad")
  build_explore_h5ad(
    x = lapply(layers, function(m) m[, rownames(cd), drop = FALSE]),
    out = out_file,
    omics = "rnaseq",
    coldata = cd,
    rowdata = rowdata,
    rowdata_columns = list(id_col = "gene_id", symbol_col = "gene_name"),
    de = de[cn],
    de_columns = de_columns[cn],
    pathways = pathways[intersect(cn, names(pathways))],
    enrichment_padj_max = list(ora = 0.05, gsea = 0.05),
    factors = colnames(cd),
    organism = gsub("_", " ", org_key),
    ref_build = if (nzchar(annotation_label)) paste0(org_key, "/", annotation_label) else org_key
  )

  eds <- new_explore_dataset(out_file)
  for (contrast in eds$contrasts()) {
    d <- eds$de(contrast)
    src <- de[[contrast]]
    if (all(is.na(d$effect_log2fc)) || (!is.null(src) && any(!is.na(src$FDR)) && all(is.na(d$score_fdr)))) {
      cat(paste0("WARNING: ", contrast, ": the log2 fold change or the FDR did not come through to the exploreDE file\n"))
    }
    cat(sprintf("%s: FDR<0.05 up %d / down %d, ORA terms %d, GSEA terms %d\n", contrast,
                sum(d$score_fdr < 0.05 & d$effect_log2fc > 0, na.rm = TRUE), sum(d$score_fdr < 0.05 & d$effect_log2fc < 0, na.rm = TRUE),
                NROW(eds$enrichment_table(contrast, "ora")), NROW(eds$enrichment_table(contrast, "gsea"))))
  }
  cat(paste0("Written ", out_file, " (grouped by ", k, ", ", round(file.size(out_file) / 1e6, 1), " MB)\n"))
}
