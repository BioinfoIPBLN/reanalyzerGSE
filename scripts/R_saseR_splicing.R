#!/usr/bin/env Rscript
#
# R_saseR_splicing.R — Differential splicing analysis using saseR with adapted offsets
#
# Usage: R_saseR_splicing.R <bam_dir> <annotation_file> <output_dir> <library_layout>
#        <strand> <samples_info> <design_file> <cores> <pattern_to_remove> <comparisons>
#
# References:
#   Segers et al. (2023) — saseR: Scalable Aberrant Splicing and Expression Retrieval
#   https://bioconductor.org/packages/release/bioc/html/saseR.html

args <- commandArgs(trailingOnly = TRUE)
bam_dir           <- args[1]
annotation_file   <- args[2]
output_dir        <- args[3]
library_layout    <- args[4]   # "SINGLE" or "PAIRED"
strand_info       <- args[5]   # "yes", "no", "reverse"
samples_info_file <- args[6]
design_file       <- args[7]
cores             <- as.integer(args[8])
pattern_to_remove <- args[9]   # "none" or regex
comparisons_arg   <- args[10]  # "no" or "AvsB,CvsD"

cat(paste0("\n\n========================================\n"))
cat(paste0("saseR splicing analysis\n"))
cat(paste0("========================================\n"))
cat(paste0("BAM directory: ", bam_dir, "\n"))
cat(paste0("Annotation: ", annotation_file, "\n"))
cat(paste0("Output: ", output_dir, "\n"))
cat(paste0("Library layout: ", library_layout, "\n"))
cat(paste0("Strand: ", strand_info, "\n"))
cat(paste0("Cores: ", cores, "\n"))
cat(paste0("Pattern to remove: ", pattern_to_remove, "\n"))
cat(paste0("Comparisons: ", comparisons_arg, "\n"))
cat(paste0("Current date: ", date(), "\n\n"))

# ── Load libraries ──
suppressPackageStartupMessages({
  library(saseR)
  library(edgeR)
  library(GenomicAlignments)
  library(GenomicFeatures)
  library(txdbmaker)
  library(BiocParallel)
})

# ── 1. Build targets data.frame from BAM files and sample info ──
cat("Building targets from BAM files and sample metadata...\n")

bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)
bam_files <- bam_files[file.exists(paste0(bam_files, ".bai")) |
                         file.exists(sub("\\.bam$", ".bai", bam_files))]
if (length(bam_files) == 0) {
  stop("No indexed BAM files found in ", bam_dir)
}
cat(paste0("Found ", length(bam_files), " indexed BAM files\n"))

# Read samples info to get condition mapping
samples_info <- read.delim(samples_info_file, header = FALSE, stringsAsFactors = FALSE)
# samples_info columns: V1=SRR/filename, V2=sample_condition, V3=condition

# Read design (condition per sample, one per line)
conditions <- readLines(design_file)
conditions <- trimws(conditions)

# Match BAM file names to sample info
bam_basenames <- basename(bam_files)
# Extract sample prefix from BAM name (remove aligner suffix like _hisat2.bam, _star.bam)
bam_sample_names <- sub("_(hisat2|star|STAR|Aligned\\.sortedByCoord\\.out)\\.bam$", "", bam_basenames)

# Build targets
targets <- data.frame(
  bam = bam_files,
  sample_name = bam_sample_names,
  stringsAsFactors = FALSE
)

# Match conditions: try matching bam_sample_names to V2 column of samples_info
targets$condition <- NA
for (i in seq_len(nrow(targets))) {
  # Try matching with sample name from samples_info V2
  match_idx <- which(sapply(samples_info$V2, function(x) {
    grepl(x, targets$sample_name[i], fixed = TRUE) ||
      grepl(targets$sample_name[i], x, fixed = TRUE)
  }))
  if (length(match_idx) > 0) {
    targets$condition[i] <- samples_info$V3[match_idx[1]]
  }
}

# Fallback: if no match, use conditions vector in order
if (any(is.na(targets$condition))) {
  cat("WARNING: Could not match all BAM files to conditions via samples_info.\n")
  cat("Falling back to design file order...\n")
  if (length(conditions) == nrow(targets)) {
    targets$condition <- conditions
  } else {
    stop("Number of conditions (", length(conditions), ") != number of BAM files (",
         nrow(targets), "). Cannot assign conditions.")
  }
}

# Apply pattern_to_remove filter
if (!is.null(pattern_to_remove) && pattern_to_remove != "none" && pattern_to_remove != "") {
  keep <- !grepl(pattern_to_remove, targets$sample_name)
  cat(paste0("Removing ", sum(!keep), " samples matching pattern '", pattern_to_remove, "'\n"))
  targets <- targets[keep, ]
  if (nrow(targets) < 2) stop("Too few samples after filtering")
}

targets$condition <- factor(targets$condition)
rownames(targets) <- targets$sample_name
cat("\nTargets:\n")
print(targets[, c("sample_name", "condition")])
cat("\n")

# ── 2. Create TxDb from annotation ──
cat("Creating TxDb from annotation file...\n")
genomeTxDb <- txdbmaker::makeTxDbFromGFF(annotation_file)
cat("TxDb created successfully.\n")

# ── 3. Determine strandedness for counting ──
ignore_strand <- TRUE
if (!is.null(strand_info) && strand_info %in% c("yes", "reverse")) {
  ignore_strand <- FALSE
}
is_single_end <- (toupper(library_layout) == "SINGLE")
cat(paste0("Counting mode: singleEnd=", is_single_end, ", ignore.strand=", ignore_strand, "\n"))

# ── 4. Differential splicing with adapted offsets (Section 3 of saseR vignette) ──
cat("\n── Differential splicing analysis with adapted offsets ──\n")

tryCatch({
  # Extract exonic parts from annotation
  cat("Extracting exonic parts from annotation...\n")
  flattenedAnnotation <- GenomicFeatures::exonicParts(genomeTxDb,
                                                       linked.to.single.gene.only = TRUE)
  cat(paste0("Exonic parts extracted: ", length(flattenedAnnotation), " bins\n"))

  # Count overlaps with BAM files
  cat("Counting overlaps with BAM files (this may take a while)...\n")
  bpparam <- MulticoreParam(workers = min(cores, nrow(targets)))
  se <- GenomicAlignments::summarizeOverlaps(
    flattenedAnnotation,
    BamFileList(targets$bam),
    singleEnd = is_single_end,
    ignore.strand = ignore_strand,
    BPPARAM = bpparam
  )
  colnames(se) <- targets$sample_name
  colData(se)$condition <- targets$condition
  cat(paste0("Counting done: ", nrow(se), " bins x ", ncol(se), " samples\n"))

  # Save raw bin counts
  raw_counts <- assays(se)$counts
  write.table(data.frame(BinID = rownames(raw_counts), raw_counts, check.names = FALSE),
              file = file.path(output_dir, "saseR_bin_counts_raw.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Compute per-gene offsets (sum of bin counts per gene)
  counts_mat <- assays(se)$counts
  gene_ids <- rowData(se)$gene_id
  offsetsGene <- aggregate(counts_mat, by = list(gene = gene_ids), FUN = sum)
  offsets_mat <- offsetsGene[match(gene_ids, offsetsGene$gene), ]
  offsets_mat$gene <- NULL
  offsets_mat <- as.matrix(offsets_mat)
  rownames(offsets_mat) <- rownames(counts_mat)

  # Handle zero offsets
  zero_idx <- offsets_mat == 0
  counts_mat[zero_idx] <- 1L
  offsets_mat[zero_idx] <- 1

  # Determine comparisons
  unique_conditions <- levels(targets$condition)
  cat(paste0("Conditions found: ", paste(unique_conditions, collapse = ", "), "\n"))

  if (comparisons_arg != "no" && comparisons_arg != "") {
    # User-specified comparisons
    comp_list <- strsplit(comparisons_arg, ",")[[1]]
    comp_pairs <- lapply(comp_list, function(x) {
      parts <- unlist(strsplit(x, "vs|VS|__VS__"))
      if (length(parts) == 2) return(trimws(parts))
      return(NULL)
    })
    comp_pairs <- Filter(Negate(is.null), comp_pairs)
  } else {
    # All pairwise comparisons
    if (length(unique_conditions) >= 2) {
      comp_pairs <- combn(unique_conditions, 2, simplify = FALSE)
    } else {
      comp_pairs <- list()
      cat("WARNING: Only one condition found. Skipping differential splicing.\n")
    }
  }

  # Run edgeR with adapted offsets for each comparison
  for (comp_idx in seq_along(comp_pairs)) {
    cond1 <- comp_pairs[[comp_idx]][1]
    cond2 <- comp_pairs[[comp_idx]][2]
    cat(paste0("\n  Comparison ", comp_idx, ": ", cond1, " vs ", cond2, "\n"))

    # Subset samples for this comparison
    keep_samples <- targets$condition %in% c(cond1, cond2)
    if (sum(keep_samples) < 2) {
      cat("    Skipping: fewer than 2 samples\n")
      next
    }

    counts_sub <- counts_mat[, keep_samples, drop = FALSE]
    offsets_sub <- offsets_mat[, keep_samples, drop = FALSE]
    cond_sub <- factor(targets$condition[keep_samples], levels = c(cond2, cond1))

    # Filter low-count bins
    keep_bins <- rowSums(counts_sub > 0) >= 2
    counts_sub <- counts_sub[keep_bins, , drop = FALSE]
    offsets_sub <- offsets_sub[keep_bins, , drop = FALSE]
    cat(paste0("    Bins after filtering: ", nrow(counts_sub), "\n"))

    # edgeR with adapted offsets
    DGE <- DGEList(counts = counts_sub)
    DGE$offset <- log(offsets_sub)
    design <- model.matrix(~ cond_sub)
    DGE <- estimateDisp(DGE, design = design)
    fitDGE <- glmFit(DGE, design = design)
    results_DGE <- glmLRT(fitDGE, coef = 2)

    # Extract results
    tt <- topTags(results_DGE, n = Inf, sort.by = "PValue")$table
    tt$BinID <- rownames(tt)
    tt$gene_id <- gene_ids[match(tt$BinID, rownames(counts_mat))]

    # Reorder columns
    tt <- tt[, c("BinID", "gene_id", "logFC", "logCPM", "LR", "PValue", "FDR")]

    # Write results
    out_file <- file.path(output_dir,
                          paste0("saseR_diff_splicing_", cond1, "_vs_", cond2, ".txt"))
    write.table(tt, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(paste0("    Results written: ", out_file, "\n"))
    cat(paste0("    Significant bins (FDR < 0.05): ", sum(tt$FDR < 0.05, na.rm = TRUE), "\n"))

    # Write summary of significant bins per gene
    sig <- tt[tt$FDR < 0.05 & !is.na(tt$FDR), ]
    if (nrow(sig) > 0) {
      gene_summary <- aggregate(BinID ~ gene_id, data = sig, FUN = length)
      colnames(gene_summary) <- c("gene_id", "n_significant_bins")
      gene_summary <- gene_summary[order(-gene_summary$n_significant_bins), ]
      out_summary <- file.path(output_dir,
                                paste0("saseR_diff_splicing_", cond1, "_vs_", cond2, "_gene_summary.txt"))
      write.table(gene_summary, file = out_summary, sep = "\t", quote = FALSE, row.names = FALSE)
      cat(paste0("    Genes with significant differential splicing: ", nrow(gene_summary), "\n"))
    }
  }

  cat("\nDifferential splicing analysis completed.\n")

}, error = function(e) {
  cat(paste0("\nERROR in differential splicing analysis: ", e$message, "\n"))
  cat("Continuing with aberrant splicing analysis...\n")
})


# ── 5. Aberrant splicing analysis (Section 2 of saseR vignette) ──
cat("\n── Aberrant splicing analysis (outlier detection) ──\n")

tryCatch({
  # Extract features using binGenome
  cat("Extracting features with binGenome (this may take a while)...\n")
  features <- saseR::binGenome(genomeTxDb, logTo = NULL)

  # Build ASpli targets (all samples treated equally for aberrant analysis)
  aspli_targets <- data.frame(
    row.names = targets$sample_name,
    bam = targets$bam,
    f1 = rep("A", nrow(targets)),
    stringsAsFactors = FALSE
  )

  # Count with ASpli
  lib_type <- ifelse(is_single_end, "SE", "PE")
  cat(paste0("Counting with BamtoAspliCounts (libType=", lib_type, ")...\n"))
  bpparam <- MulticoreParam(workers = min(cores, nrow(targets)))
  ASpliSE <- saseR::BamtoAspliCounts(
    features = features,
    targets = aspli_targets,
    minReadLength = 50,
    libType = lib_type,
    BPPARAM = bpparam
  )

  # Convert to SummarizedExperiment objects
  SEbins <- saseR::convertASpli(ASpliSE, type = "bin")
  SEjunctions <- saseR::convertASpli(ASpliSE, type = "junction")

  # Set design (intercept only for aberrant analysis)
  metadata(SEbins)$design <- ~1
  metadata(SEjunctions)$design <- ~1

  # Calculate offsets
  cat("Calculating offsets for bins...\n")
  SEbins <- saseR::calculateOffsets(SEbins, method = "AS", aggregation = "locus")
  cat("Calculating offsets for junctions...\n")
  SEjunctions <- saseR::calculateOffsets(SEjunctions, method = "AS",
                                          aggregation = "symbol",
                                          mergeGeneASpli = TRUE)

  # Filter low-count features
  filterbins <- edgeR::filterByExpr(SEbins)
  SEbins <- SEbins[filterbins, ]
  filterjunctions <- edgeR::filterByExpr(SEjunctions)
  SEjunctions <- SEjunctions[filterjunctions, ]
  cat(paste0("After filtering: ", nrow(SEbins), " bins, ", nrow(SEjunctions), " junctions\n"))

  # Find optimal encoding dimension
  cat("Finding optimal encoding dimension...\n")
  SEbins <- saseR::saseRfindEncodingDim(SEbins, method = "GD")
  SEjunctions <- saseR::saseRfindEncodingDim(SEjunctions, method = "GD")

  # Fit model
  cat("Fitting saseR model for bins...\n")
  SEbins <- saseR::saseRfit(SEbins, analysis = "AS", padjust = "BH", fit = "fast")
  cat("Fitting saseR model for junctions...\n")
  SEjunctions <- saseR::saseRfit(SEjunctions, analysis = "AS", padjust = "BH", fit = "fast")

  # Extract results: for each sample, report significant aberrant bins/junctions
  cat("\nExtracting aberrant splicing results...\n")

  # Bins results
  if ("pValueAdjust" %in% assayNames(SEbins)) {
    pvals_bins <- assays(SEbins)$pValueAdjust
    for (s in seq_len(ncol(SEbins))) {
      sample_name <- colnames(SEbins)[s]
      sig_idx <- which(pvals_bins[, s] < 0.05)
      if (length(sig_idx) > 0) {
        result_df <- data.frame(
          BinID = rownames(SEbins)[sig_idx],
          pValueAdjust = pvals_bins[sig_idx, s],
          stringsAsFactors = FALSE
        )
        if ("locus" %in% colnames(rowData(SEbins))) {
          result_df$gene <- rowData(SEbins)$locus[sig_idx]
        }
        result_df <- result_df[order(result_df$pValueAdjust), ]
        out_file <- file.path(output_dir,
                               paste0("saseR_aberrant_bins_", sample_name, ".txt"))
        write.table(result_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  }

  # Junctions results
  if ("pValueAdjust" %in% assayNames(SEjunctions)) {
    pvals_junctions <- assays(SEjunctions)$pValueAdjust
    for (s in seq_len(ncol(SEjunctions))) {
      sample_name <- colnames(SEjunctions)[s]
      sig_idx <- which(pvals_junctions[, s] < 0.05)
      if (length(sig_idx) > 0) {
        result_df <- data.frame(
          JunctionID = rownames(SEjunctions)[sig_idx],
          pValueAdjust = pvals_junctions[sig_idx, s],
          stringsAsFactors = FALSE
        )
        if ("locus" %in% colnames(rowData(SEjunctions))) {
          result_df$gene <- rowData(SEjunctions)$locus[sig_idx]
        }
        result_df <- result_df[order(result_df$pValueAdjust), ]
        out_file <- file.path(output_dir,
                               paste0("saseR_aberrant_junctions_", sample_name, ".txt"))
        write.table(result_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  }

  # Summary table: number of aberrant bins and junctions per sample
  summary_df <- data.frame(
    sample = colnames(SEbins),
    condition = targets$condition[match(colnames(SEbins), targets$sample_name)],
    n_aberrant_bins = colSums(pvals_bins < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  if ("pValueAdjust" %in% assayNames(SEjunctions)) {
    junc_counts <- colSums(pvals_junctions < 0.05, na.rm = TRUE)
    summary_df$n_aberrant_junctions <- junc_counts[match(summary_df$sample,
                                                          colnames(SEjunctions))]
  }
  write.table(summary_df, file = file.path(output_dir, "saseR_aberrant_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nAberrant splicing summary:\n")
  print(summary_df)

  cat("\nAberrant splicing analysis completed.\n")

}, error = function(e) {
  cat(paste0("\nERROR in aberrant splicing analysis: ", e$message, "\n"))
  cat("This analysis is optional. The differential splicing results (if any) are still available.\n")
})

cat(paste0("\n\nsaseR splicing analysis completed. Current date: ", date(), "\n"))
