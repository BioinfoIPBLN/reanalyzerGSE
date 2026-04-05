#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

input_raw_counts <- args[1]
input_tpm_counts <- args[2]
input_metadata <- args[3]
input_list_comp <- args[4]
input_path_deg_results <- args[5]
pattern_deg_results <- args[6]
bakk_name <- args[7]

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(annotables))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))

# Source shared ENSEMBL helper
script_dir <- dirname(sub("^--file=", "", commandArgs()[grep("--file=", commandArgs())]))
ensembl_helper <- file.path(script_dir, "R_ensembl_to_symbol.R")
if (file.exists(ensembl_helper)) source(ensembl_helper)

cat("\n========================================")
cat("\n  prepare_SE.R - With Pathway Analysis")
cat("\n========================================\n")

# === Counts ====
cat("\n[1/6] Reading raw counts...")
raw_counts <- read_tsv(input_raw_counts, show_col_types = FALSE) %>% 
  data.frame(check.names = F)
if ("Length" %in% colnames(raw_counts)){
  raw_counts <- raw_counts %>% dplyr::select(-Length)
}
if ("Gene_ID" %in% colnames(raw_counts)){
  raw_counts <- raw_counts %>% column_to_rownames("Gene_ID")
}

cat("\n[2/6] Reading TPM counts...")
tpm_counts <- data.table::fread(input_tpm_counts) %>% 
  data.frame(check.names = F) %>% 
  column_to_rownames("V1")

cat("\n      Sanity check: identical colnames raw counts and tpm counts: ", 
    identical(colnames(raw_counts), colnames(tpm_counts)))
cat("\n      Sanity check: identical genes in raw counts and tpm counts: ", 
    identical(rownames(raw_counts), rownames(tpm_counts)))
cat("\n      Subsetting genes in raw counts to the ones in tpm counts...")
raw_counts <- raw_counts[match(rownames(tpm_counts), rownames(raw_counts)),]
rownames(raw_counts) <- toupper(rownames(raw_counts))
rownames(tpm_counts) <- toupper(rownames(tpm_counts))
cat("\n      Sanity check: identical genes in raw counts and tpm counts: ", 
    identical(rownames(raw_counts), rownames(tpm_counts)))

# Detect if gene IDs are ENSEMBL format
gene_ids_are_ensembl <- exists("is_ensembl_id") && is_ensembl_id(rownames(tpm_counts))
if (gene_ids_are_ensembl) {
  cat("\n      Detected ENSEMBL gene IDs (e.g.,", head(rownames(tpm_counts), 2), ")")
}


# === rowData (DE results) ===
cat("\n\n[3/6] Reading DE results and comparisons...")
list_comps <- read_tsv(input_list_comp, col_names = F, show_col_types = FALSE)
deg_results_files <- list.files(input_path_deg_results, pattern_deg_results, full.names = T)

# Parse comparison names and filter out empty/blank ones
comparisons <- lapply(1:length(deg_results_files), function(i) {
  paste(grep(paste0("Comparison number ", i), list_comps$X1, value = TRUE) %>% gsub(".*__", "", .), collapse = "__")
}) %>% unlist
# Filter out empty comparison names (from trailing empty entries in list_comp.txt)
valid_comps <- nchar(trimws(comparisons)) > 0
if (sum(!valid_comps) > 0) {
  cat("\n      Filtering out", sum(!valid_comps), "empty comparison entries...")
}
comparisons <- comparisons[valid_comps]
deg_results_files <- deg_results_files[valid_comps]
cat("\n      Found", length(comparisons), "valid comparisons:", paste(comparisons, collapse = ", "))

de_results <- read_tsv(deg_results_files[1], show_col_types = FALSE)
check <- identical(sort(rownames(tpm_counts)), sort(toupper(de_results$Gene_ID)))
cat("\n      Sanity check: identical genes in tpm counts and DE table: ", check)
if(!check) {
  cat("\n      WARNING: The genes that differ will be discarded...")
}

# Build DE results list with annotations
# For ENSEMBL IDs: match via ensgene (stripping versions) instead of symbol
de_results_list <- lapply(1:length(deg_results_files), function(i) {
  de_results <- read_tsv(deg_results_files[i], show_col_types = FALSE)
  colnames(de_results) <- c("gene_name", "gene_length", "log2Ratio", "logCPM", "pValue", "fdr")
  de_results$gene_name <- toupper(de_results$gene_name)
  de_results <- de_results[match(rownames(tpm_counts), de_results$gene_name),]
  de_results <- de_results[complete.cases(de_results),]
  
  if (gene_ids_are_ensembl) {
    # Match ENSEMBL IDs (stripping version) against grch38$ensgene
    ids_stripped <- toupper(strip_ensembl_version(de_results$gene_name))
    grch38_stripped <- toupper(grch38$ensgene)
    match_idx <- match(ids_stripped, grch38_stripped)
    de_results$gene_id <- grch38$ensgene[match_idx]
    de_results$description <- grch38$description[match_idx]
    de_results$entrez_id <- grch38$entrez[match_idx]
  } else {
    # Standard symbol-based matching
    de_results$gene_id <- grch38$ensgene[match(de_results$gene_name, grch38$symbol)]
    de_results$description <- grch38$description[match(de_results$gene_name, grch38$symbol)]
    de_results$entrez_id <- grch38$entrez[match(de_results$gene_name, grch38$symbol)]
  }
  
  de_results
}) %>% setNames(paste0("contrast_", comparisons))

# Use common genes across all DE result tables for rowData alignment
common_genes <- Reduce(intersect, lapply(de_results_list, function(x) x$gene_name))
cat("\n      Genes common across all contrasts:", length(common_genes))

# Subset DE results to common genes
de_results_list <- lapply(de_results_list, function(x) {
  x[match(common_genes, x$gene_name), ]
})

de_results_DF <- DataFrame(row.names = de_results_list[[1]]$gene_name)
de_results_DF@listData <- de_results_list


# === colData (sample metadata) ===
cat("\n\n[4/6] Reading sample metadata...")
metadata <- read_tsv(input_metadata, col_names = F, show_col_types = FALSE)
colnames(metadata) <- c("Name","Tissue [Factor]","Batch [Factor]")
# Try matching with and without suffixes (handle _STAR.bam, _condition, etc.)
sample_names <- colnames(raw_counts)
# Try direct match first
direct_match <- match(sample_names, metadata$Name)
if (all(!is.na(direct_match))) {
  metadata <- metadata[direct_match, ]
} else {
  # Try fuzzy match: strip common suffixes from sample names
  sample_cores <- sub("(_STAR\\.bam|_hisat2\\.bam|_[^_]+)$", "", sample_names)
  meta_cores <- sub("(_STAR\\.bam|_hisat2\\.bam|_[^_]+)$", "", metadata$Name)
  fuzzy_match <- match(sample_cores, meta_cores)
  if (all(!is.na(fuzzy_match))) {
    metadata <- metadata[fuzzy_match, ]
  } else {
    # Last resort: keep order as-is, pad if needed
    metadata <- metadata[1:length(sample_names), ]
  }
}
matched_count <- sum(!is.na(match(sample_names, metadata$Name)))
cat("\n      Samples matched: ", matched_count, "/", nrow(metadata))


# Subset counts to common genes
raw_counts <- raw_counts[rownames(raw_counts) %in% common_genes, , drop = FALSE]
tpm_counts <- tpm_counts[rownames(tpm_counts) %in% common_genes, , drop = FALSE]
# Ensure identical row ordering across assays and rowData
raw_counts <- raw_counts[match(common_genes, rownames(raw_counts)), , drop = FALSE]
tpm_counts <- tpm_counts[match(common_genes, rownames(tpm_counts)), , drop = FALSE]
cat("\n      Final gene count for SE:", nrow(raw_counts))


# === Pathway Analysis ===
cat("\n\n[5/6] Computing pathway analyses for each contrast...")
cat("\n      Using FDR < 0.05 for DEG selection")
cat("\n      Universe: genes in TPM counts (n =", nrow(tpm_counts), ")\n")

# Get universe (all genes in TPM counts with Entrez IDs) - must be CHARACTER for clusterProfiler
universe_genes <- de_results_list[[1]]$gene_name
universe_entrez <- as.character(de_results_list[[1]]$entrez_id)
universe_entrez <- universe_entrez[!is.na(universe_entrez) & universe_entrez != "NA"]
cat("      Genes with Entrez IDs:", length(universe_entrez), "\n")

# Parameters for pathway analysis
param <- list(
  runGO = TRUE,
  refBuild = "Homo_sapiens/GENCODE/GRCh38",
  fdrThreshORA = 0.05,
  fdrThreshGSEA = 0.05,
  pValThreshGO = 0.05,
  log2RatioThreshGO = 0  # No LFC threshold, just FDR
)

# Store enrichment results per contrast
enrichResultList <- list()
enrichInputList <- list()

for (contrast_name in names(de_results_list)) {
  cat("\n      Processing:", contrast_name)
  
  de_data <- de_results_list[[contrast_name]]
  
  # Select DEGs based on FDR < 0.05
  sig_genes <- de_data[de_data$fdr < 0.05 & !is.na(de_data$fdr), ]
  up_genes <- sig_genes[sig_genes$log2Ratio > 0, ]
  down_genes <- sig_genes[sig_genes$log2Ratio < 0, ]
  
  cat(" (UP:", nrow(up_genes), ", DOWN:", nrow(down_genes), ")")
  
  # Get Entrez IDs (must be CHARACTER for clusterProfiler)
  up_entrez <- as.character(up_genes$entrez_id[!is.na(up_genes$entrez_id)])
  up_entrez <- up_entrez[up_entrez != "NA"]
  down_entrez <- as.character(down_genes$entrez_id[!is.na(down_genes$entrez_id)])
  down_entrez <- down_entrez[down_entrez != "NA"]
  both_entrez <- c(up_entrez, down_entrez)
  
  # Store enrichment input
  enrichInputList[[contrast_name]] <- list(
    selections = list(
      upGenes = up_genes$gene_id[!is.na(up_genes$gene_id)],
      downGenes = down_genes$gene_id[!is.na(down_genes$gene_id)],
      bothGenes = sig_genes$gene_id[!is.na(sig_genes$gene_id)]
    ),
    seqAnno = de_data,
    log2Ratio = setNames(de_data$log2Ratio, de_data$gene_id)
  )
  
  # Initialize ORA results structure
  ora_results <- list(
    BP = list(upGenes = NULL, downGenes = NULL, bothGenes = NULL),
    MF = list(upGenes = NULL, downGenes = NULL, bothGenes = NULL),
    CC = list(upGenes = NULL, downGenes = NULL, bothGenes = NULL)
  )
  
  # Initialize GSEA results structure
  gsea_results <- list(BP = NULL, MF = NULL, CC = NULL)
  
  # Run ORA for each GO category and gene set
  for (ont in c("BP", "MF", "CC")) {
    # ORA for upGenes
    if (length(up_entrez) >= 5) {
      tryCatch({
        ora_results[[ont]]$upGenes <- enrichGO(
          gene = up_entrez,
          universe = universe_entrez,
          OrgDb = org.Hs.eg.db,
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = 1,  # Keep all for visualization
          qvalueCutoff = 1,
          readable = TRUE
        )
        # Add geneName column for compatibility
        if (!is.null(ora_results[[ont]]$upGenes) && nrow(ora_results[[ont]]$upGenes@result) > 0) {
          ora_results[[ont]]$upGenes@result$geneName <- ora_results[[ont]]$upGenes@result$geneID
        }
      }, error = function(e) cat(" [ORA-UP-", ont, " failed]", sep = ""))
    }
    
    # ORA for downGenes
    if (length(down_entrez) >= 5) {
      tryCatch({
        ora_results[[ont]]$downGenes <- enrichGO(
          gene = down_entrez,
          universe = universe_entrez,
          OrgDb = org.Hs.eg.db,
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          readable = TRUE
        )
        if (!is.null(ora_results[[ont]]$downGenes) && nrow(ora_results[[ont]]$downGenes@result) > 0) {
          ora_results[[ont]]$downGenes@result$geneName <- ora_results[[ont]]$downGenes@result$geneID
        }
      }, error = function(e) cat(" [ORA-DOWN-", ont, " failed]", sep = ""))
    }
    
    # ORA for bothGenes
    if (length(both_entrez) >= 5) {
      tryCatch({
        ora_results[[ont]]$bothGenes <- enrichGO(
          gene = both_entrez,
          universe = universe_entrez,
          OrgDb = org.Hs.eg.db,
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          readable = TRUE
        )
        if (!is.null(ora_results[[ont]]$bothGenes) && nrow(ora_results[[ont]]$bothGenes@result) > 0) {
          ora_results[[ont]]$bothGenes@result$geneName <- ora_results[[ont]]$bothGenes@result$geneID
        }
      }, error = function(e) cat(" [ORA-BOTH-", ont, " failed]", sep = ""))
    }
    
    # GSEA - requires ranked gene list with character Entrez IDs as names
    # Must remove duplicates - keep the gene with highest absolute log2FC per Entrez ID
    gene_list <- de_data$log2Ratio
    names(gene_list) <- as.character(de_data$entrez_id)
    gene_list <- gene_list[!is.na(names(gene_list)) & names(gene_list) != "NA"]
    
    # Remove duplicates: for each Entrez ID, keep the value with max absolute log2FC
    if (any(duplicated(names(gene_list)))) {
      gene_df <- data.frame(entrez = names(gene_list), lfc = gene_list, abs_lfc = abs(gene_list))
      gene_df <- gene_df[order(gene_df$abs_lfc, decreasing = TRUE), ]
      gene_df <- gene_df[!duplicated(gene_df$entrez), ]
      gene_list <- setNames(gene_df$lfc, gene_df$entrez)
    }
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    if (length(gene_list) >= 10) {
      tryCatch({
        gsea_results[[ont]] <- gseGO(
          geneList = gene_list,
          OrgDb = org.Hs.eg.db,
          ont = ont,
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          verbose = FALSE
        )
        # Add geneName column for compatibility
        if (!is.null(gsea_results[[ont]]) && nrow(gsea_results[[ont]]@result) > 0) {
          gsea_results[[ont]]@result$geneName <- gsea_results[[ont]]@result$core_enrichment
        }
      }, error = function(e) cat(" [GSEA-", ont, " failed: ", conditionMessage(e), "]", sep = ""))
    }
  }
  
  enrichResultList[[contrast_name]] <- list(
    ora = ora_results,
    gsea = gsea_results
  )
}


# === Build GO annotation columns for GO Tile plots ===
cat("\n\n      Adding GO annotations to rowData...")
# Get GO terms for each gene
gene_symbols <- de_results_list[[1]]$gene_name

if (gene_ids_are_ensembl) {
  # For ENSEMBL IDs: use entrez_id directly from the already-mapped de_results
  entrez_ids <- de_results_list[[1]]$entrez_id
} else {
  entrez_ids <- grch38$entrez[match(gene_symbols, grch38$symbol)]
}

# Map to GO terms  
valid_entrez <- as.character(entrez_ids[!is.na(entrez_ids)])
valid_entrez <- valid_entrez[valid_entrez != "NA" & nchar(valid_entrez) > 0]

go_bp <- tryCatch({
  if (length(valid_entrez) == 0) list() else {
    res <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = valid_entrez,
                                  columns = "GO",
                                  keytype = "ENTREZID")
    res <- res[res$ONTOLOGY == "BP", ]
    split(res$GO, res$ENTREZID)
  }
}, error = function(e) list())

go_mf <- tryCatch({
  if (length(valid_entrez) == 0) list() else {
    res <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = valid_entrez,
                                  columns = "GO",
                                  keytype = "ENTREZID")
    res <- res[res$ONTOLOGY == "MF", ]
    split(res$GO, res$ENTREZID)
  }
}, error = function(e) list())

go_cc <- tryCatch({
  if (length(valid_entrez) == 0) list() else {
    res <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = valid_entrez,
                                  columns = "GO",
                                  keytype = "ENTREZID")
    res <- res[res$ONTOLOGY == "CC", ]
    split(res$GO, res$ENTREZID)
  }
}, error = function(e) list())

# Add GO columns to each DE result
for (contrast_name in names(de_results_list)) {
  de_data <- de_results_list[[contrast_name]]
  de_data$`GO BP` <- sapply(as.character(de_data$entrez_id), function(e) {
    if (is.na(e) || !e %in% names(go_bp)) return(NA_character_)
    paste(go_bp[[e]], collapse = "; ")
  })
  de_data$`GO MF` <- sapply(as.character(de_data$entrez_id), function(e) {
    if (is.na(e) || !e %in% names(go_mf)) return(NA_character_)
    paste(go_mf[[e]], collapse = "; ")
  })
  de_data$`GO CC` <- sapply(as.character(de_data$entrez_id), function(e) {
    if (is.na(e) || !e %in% names(go_cc)) return(NA_character_)
    paste(go_cc[[e]], collapse = "; ")
  })
  de_results_list[[contrast_name]] <- de_data
  
  # Also update enrichInput seqAnno
  enrichInputList[[contrast_name]]$seqAnno <- de_data
}

# Update de_results_DF with new data
de_results_DF <- DataFrame(row.names = de_results_list[[1]]$gene_name)
de_results_DF@listData <- de_results_list


# === Make SE ===
cat("\n\n[6/6] Building SummarizedExperiment object...")

# Final alignment check
stopifnot("Assay/rowData dimension mismatch" = 
  nrow(raw_counts) == nrow(de_results_DF) && 
  nrow(tpm_counts) == nrow(de_results_DF))
stopifnot("Assay/rowData row names mismatch" = 
  identical(rownames(raw_counts), rownames(de_results_DF)))

se <- SummarizedExperiment(
  assays = SimpleList("counts" = as.matrix(raw_counts), "TPM" = as.matrix(tpm_counts)),
  rowData = de_results_DF,
  colData = DataFrame(metadata, check.names = F)
)

# Add pathway analysis to metadata
metadata(se)$param <- param
metadata(se)$enrichResultList <- enrichResultList
metadata(se)$enrichInputList <- enrichInputList

cat("\n      SE object created successfully!")
cat("\n      Contrasts:", paste(names(enrichResultList), collapse = ", "))
cat("\n      Assays:", paste(names(assays(se)), collapse = ", "))
cat("\n      Pathway data: enrichResultList, enrichInputList in metadata(se)")


# === Save ===
output_dir <- dirname(sub("^--file=", "", commandArgs()[grep("--file=", commandArgs())]))
qs2::qs_save(se, file.path(output_dir, paste0("deResults_", bakk_name, ".qs2")))
qs2::qs_save(se, file.path(output_dir, "deResults.qs2"))

cat("\n\n========================================")
cat("\n  Output saved to:")
cat("\n    ", file.path(output_dir, paste0("deResults_", bakk_name, ".qs2")))
cat("\n    ", file.path(output_dir, "deResults.qs2"))
cat("\n========================================\n\n")
