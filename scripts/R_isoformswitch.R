#!/usr/bin/env Rscript
#
# R_isoformswitch.R — Isoform switch / alternative splicing analysis using IsoformSwitchAnalyzeR
#
# Usage: R_isoformswitch.R <quant_dir> <annotation_file> <transcript_fasta> <output_dir>
#        <samples_info> <design_file> <cores> <pattern_to_remove> <comparisons>
#        <aligner>
#
# Reference:
#   Vitting-Seerup & Sandelin (2019) — IsoformSwitchAnalyzeR
#   https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html

args <- commandArgs(trailingOnly = TRUE)
quant_dir          <- args[1]   # Kallisto results dir or Salmon quant dir
annotation_file    <- args[2]   # GTF annotation
transcript_fasta   <- args[3]   # Transcript FASTA (cDNA)
output_dir         <- args[4]
samples_info_file  <- args[5]
design_file        <- args[6]
cores              <- as.integer(args[7])
pattern_to_remove  <- args[8]   # "none" or regex
comparisons_arg    <- args[9]   # "no" or "AvsB,CvsD"
aligner            <- args[10]  # "star", "hisat2", or "kallisto"

cat(paste0("\n\n========================================\n"))
cat(paste0("IsoformSwitchAnalyzeR analysis\n"))
cat(paste0("========================================\n"))
cat(paste0("Quantification dir: ", quant_dir, "\n"))
cat(paste0("Annotation: ", annotation_file, "\n"))
cat(paste0("Transcript FASTA: ", transcript_fasta, "\n"))
cat(paste0("Output: ", output_dir, "\n"))
cat(paste0("Aligner: ", aligner, "\n"))
cat(paste0("Cores: ", cores, "\n"))
cat(paste0("Pattern to remove: ", pattern_to_remove, "\n"))
cat(paste0("Comparisons: ", comparisons_arg, "\n"))
cat(paste0("Current date: ", date(), "\n\n"))

# ── Load libraries ──
suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
})

# ── 1. Read sample info and build design matrix ──
cat("Building design matrix from sample metadata...\n")

samples_info <- read.delim(samples_info_file, header = FALSE, stringsAsFactors = FALSE)
# V1=SRR/filename, V2=sample_name_condition, V3=condition

conditions <- readLines(design_file)
conditions <- trimws(conditions)

# ── 2. Import quantification data ──
cat("Importing quantification data...\n")

if (aligner == "kallisto") {
  # Look for Kallisto abundance.tsv files
  abundance_files <- list.files(quant_dir, pattern = "abundance.tsv",
                                 recursive = TRUE, full.names = TRUE)
  if (length(abundance_files) == 0) {
    stop("No abundance.tsv files found in ", quant_dir, " for Kallisto import")
  }
  cat(paste0("Found ", length(abundance_files), " Kallisto abundance files\n"))

  # Import using IsoformSwitchAnalyzeR's function
  kallistoQuant <- importIsoformExpression(sampleVector = abundance_files)

  # Build sample names from file paths
  sample_names <- basename(dirname(abundance_files))

  # Build design matrix
  myDesign <- data.frame(
    sampleID = sample_names,
    condition = character(length(sample_names)),
    stringsAsFactors = FALSE
  )

  # Match conditions: try matching sample names to V2 of samples_info
  for (i in seq_along(sample_names)) {
    match_idx <- which(sapply(samples_info$V2, function(x) {
      grepl(x, sample_names[i], fixed = TRUE) ||
        grepl(sample_names[i], x, fixed = TRUE)
    }))
    if (length(match_idx) > 0) {
      myDesign$condition[i] <- samples_info$V3[match_idx[1]]
    }
  }

  # Fallback to design file order
  if (any(myDesign$condition == "")) {
    cat("WARNING: Could not match all samples to conditions via samples_info.\n")
    if (length(conditions) == nrow(myDesign)) {
      myDesign$condition <- conditions
    } else {
      stop("Number of conditions != number of samples")
    }
  }

} else {
  # For STAR/hisat2: look for Salmon quant results or featureCounts at transcript level
  # Check if Salmon quant was run (e.g. during strandness prediction)
  salmon_dirs <- list.dirs(quant_dir, recursive = FALSE)
  salmon_quant_files <- list.files(quant_dir, pattern = "quant.sf",
                                    recursive = TRUE, full.names = TRUE)

  if (length(salmon_quant_files) > 0) {
    cat(paste0("Found ", length(salmon_quant_files), " Salmon quant.sf files\n"))
    kallistoQuant <- importIsoformExpression(sampleVector = salmon_quant_files)

    sample_names <- basename(dirname(salmon_quant_files))
    myDesign <- data.frame(
      sampleID = sample_names,
      condition = character(length(sample_names)),
      stringsAsFactors = FALSE
    )
    for (i in seq_along(sample_names)) {
      match_idx <- which(sapply(samples_info$V2, function(x) {
        grepl(x, sample_names[i], fixed = TRUE) ||
          grepl(sample_names[i], x, fixed = TRUE)
      }))
      if (length(match_idx) > 0) {
        myDesign$condition[i] <- samples_info$V3[match_idx[1]]
      }
    }
    if (any(myDesign$condition == "")) {
      if (length(conditions) == nrow(myDesign)) {
        myDesign$condition <- conditions
      } else {
        stop("Number of conditions != number of samples")
      }
    }
  } else {
    stop("IsoformSwitchAnalyzeR requires transcript-level quantification (Kallisto or Salmon).\n",
         "No abundance.tsv (Kallisto) or quant.sf (Salmon) files found in: ", quant_dir, "\n",
         "Please re-run the pipeline with aligner='kallisto', or provide Salmon quantification.\n",
         "Alternatively, run Salmon on the BAM files and point to the Salmon output directory.")
  }
}

# Apply pattern_to_remove filter
if (!is.null(pattern_to_remove) && pattern_to_remove != "none" && pattern_to_remove != "") {
  keep <- !grepl(pattern_to_remove, myDesign$sampleID)
  cat(paste0("Removing ", sum(!keep), " samples matching pattern '", pattern_to_remove, "'\n"))
  myDesign <- myDesign[keep, ]
  # Also filter the count/abundance matrices
  keep_cols <- c(TRUE, colnames(kallistoQuant$counts)[-1] %in% myDesign$sampleID)
  kallistoQuant$counts <- kallistoQuant$counts[, keep_cols, drop = FALSE]
  kallistoQuant$abundance <- kallistoQuant$abundance[, keep_cols, drop = FALSE]
}

cat("\nDesign matrix:\n")
print(myDesign)
cat("\n")

# ── 3. Determine comparisons ──
unique_conditions <- unique(myDesign$condition)
cat(paste0("Conditions found: ", paste(unique_conditions, collapse = ", "), "\n"))

if (comparisons_arg != "no" && comparisons_arg != "") {
  comp_list <- strsplit(comparisons_arg, ",")[[1]]
  myComparisons <- data.frame(
    condition_1 = character(0),
    condition_2 = character(0),
    stringsAsFactors = FALSE
  )
  for (comp in comp_list) {
    parts <- unlist(strsplit(comp, "vs|VS|__VS__"))
    if (length(parts) == 2) {
      myComparisons <- rbind(myComparisons, data.frame(
        condition_1 = trimws(parts[1]),
        condition_2 = trimws(parts[2]),
        stringsAsFactors = FALSE
      ))
    }
  }
} else {
  # All pairwise comparisons
  if (length(unique_conditions) >= 2) {
    pairs <- combn(unique_conditions, 2)
    myComparisons <- data.frame(
      condition_1 = pairs[1, ],
      condition_2 = pairs[2, ],
      stringsAsFactors = FALSE
    )
  } else {
    stop("Only one condition found. IsoformSwitchAnalyzeR requires at least 2 conditions.")
  }
}

cat("\nComparisons to make:\n")
print(myComparisons)
cat("\n")

# ── 4. Create switchAnalyzeRlist ──
cat("Creating switchAnalyzeRlist (importing data)...\n")

# Check if transcript_fasta exists
if (!file.exists(transcript_fasta)) {
  cat(paste0("WARNING: Transcript FASTA not found: ", transcript_fasta, "\n"))
  cat("Attempting to proceed without it (some consequence analyses may be limited)...\n")
  transcript_fasta <- NULL
}

aSwitchList <- importRdata(
  isoformCountMatrix = kallistoQuant$counts,
  isoformRepExpression = kallistoQuant$abundance,
  designMatrix = myDesign,
  comparisonsToMake = myComparisons,
  isoformExonAnnoation = annotation_file,
  isoformNtFasta = transcript_fasta,
  showProgress = TRUE
)

cat("\nswitchAnalyzeRlist summary:\n")
print(summary(aSwitchList))
cat("\n")

# ── 5. Pre-filtering ──
cat("Pre-filtering low-expression genes and isoforms...\n")
aSwitchList <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)
cat(paste0("After filtering: ", nrow(aSwitchList$isoformFeatures), " isoform features\n"))

# ── 6. Test for isoform switches using DEXSeq ──
cat("\nTesting for isoform switches using DEXSeq...\n")
cat("(This may take a while depending on dataset size)\n")

alpha <- 0.05
dIFcutoff <- 0.10

aSwitchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchList,
  reduceToSwitchingGenes = FALSE,
  alpha = alpha,
  dIFcutoff = dIFcutoff,
  onlySigIsoforms = FALSE
)

# Extract and save switch summary
switch_summary <- extractSwitchSummary(aSwitchList)
cat("\nSwitch summary:\n")
print(switch_summary)
write.table(switch_summary, file = file.path(output_dir, "isoformswitch_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Extract top switches
top_switches <- extractTopSwitches(
  aSwitchList,
  filterForConsequences = FALSE,
  n = 50,
  inEachComparison = TRUE
)
write.table(top_switches, file = file.path(output_dir, "isoformswitch_top_switches.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 7. Analyze alternative splicing ──
cat("\nAnalyzing alternative splicing...\n")
aSwitchList <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = aSwitchList,
  onlySwitchingGenes = FALSE,
  quiet = TRUE
)

# Save splicing analysis results
if (!is.null(aSwitchList$AlternativeSplicingAnalysis)) {
  write.table(aSwitchList$AlternativeSplicingAnalysis,
              file = file.path(output_dir, "isoformswitch_alternative_splicing.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Save splicing summary plot
  tryCatch({
    pdf(file.path(output_dir, "isoformswitch_splicing_summary.pdf"), width = 10, height = 7)
    extractSplicingSummary(aSwitchList, plotGenes = FALSE)
    dev.off()
    cat("Splicing summary plot saved.\n")
  }, error = function(e) cat(paste0("Could not generate splicing summary plot: ", e$message, "\n")))

  # Splicing enrichment
  tryCatch({
    pdf(file.path(output_dir, "isoformswitch_splicing_enrichment.pdf"), width = 10, height = 7)
    extractSplicingEnrichment(aSwitchList, splicingToAnalyze = 'all')
    dev.off()
    cat("Splicing enrichment plot saved.\n")
  }, error = function(e) cat(paste0("Could not generate splicing enrichment plot: ", e$message, "\n")))
}

# ── 8. Analyze switch consequences (intron retention at minimum) ──
cat("\nAnalyzing switch consequences...\n")
tryCatch({
  # Start with intron_retention (always available)
  consequencesOfInterest <- c("intron_retention")

  # Add more consequences if ORF info is available
  if (!is.null(aSwitchList$orfAnalysis) && any(!is.na(aSwitchList$orfAnalysis$orf_origin))) {
    consequencesOfInterest <- c(consequencesOfInterest,
                                 "ORF_seq_similarity",
                                 "NMD_status")
  }

  aSwitchList <- analyzeSwitchConsequences(
    aSwitchList,
    consequencesToAnalyze = consequencesOfInterest,
    dIFcutoff = 0.10
  )

  # Summary with consequences
  switch_summary_consequences <- extractSwitchSummary(aSwitchList,
                                                       filterForConsequences = TRUE)
  write.table(switch_summary_consequences,
              file = file.path(output_dir, "isoformswitch_summary_with_consequences.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nSwitch summary (with consequences filter):\n")
  print(switch_summary_consequences)

}, error = function(e) {
  cat(paste0("Could not analyze consequences: ", e$message, "\n"))
  cat("This is expected if external tools (CPC2, Pfam, etc.) were not run.\n")
})

# ── 9. Generate volcano plots ──
cat("\nGenerating volcano and summary plots...\n")
tryCatch({
  library(ggplot2)

  datToPlot <- aSwitchList$isoformFeatures
  minPValToPlot <- 1E-50
  datToPlot$isoform_switch_q_value <- ifelse(
    datToPlot$isoform_switch_q_value < minPValToPlot,
    minPValToPlot,
    datToPlot$isoform_switch_q_value
  )

  p1 <- ggplot(data = datToPlot, aes(x = dIF, y = -log10(isoform_switch_q_value))) +
    geom_point(aes(color = abs(dIF) > dIFcutoff & isoform_switch_q_value < alpha), size = 1) +
    geom_hline(yintercept = -log10(alpha), linetype = 'dashed') +
    geom_vline(xintercept = c(-dIFcutoff, dIFcutoff), linetype = 'dashed') +
    facet_wrap(~ condition_1 + condition_2) +
    scale_color_manual('Significant\nIsoform Switch', values = c('grey40', 'red')) +
    labs(x = 'dIF', y = '-Log10 (Isoform Switch Q Value)') +
    theme_bw()

  ggsave(file.path(output_dir, "isoformswitch_volcano.pdf"), p1, width = 10, height = 7)

  p2 <- ggplot(data = aSwitchList$isoformFeatures,
               aes(x = gene_log2_fold_change, y = dIF)) +
    geom_point(aes(color = abs(dIF) > dIFcutoff & isoform_switch_q_value < alpha), size = 1) +
    facet_wrap(~ condition_1 + condition_2) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    scale_color_manual('Significant\nIsoform Switch', values = c('grey40', 'red')) +
    labs(x = 'Gene log2 fold change', y = 'dIF') +
    theme_bw()

  ggsave(file.path(output_dir, "isoformswitch_genelfc_vs_dIF.pdf"), p2, width = 10, height = 7)
  cat("Volcano and dIF plots saved.\n")
}, error = function(e) cat(paste0("Could not generate plots: ", e$message, "\n")))

# ── 10. Generate switch plots for top genes ──
cat("\nGenerating switch plots for top genes...\n")
tryCatch({
  sig_genes <- unique(top_switches$gene_name)
  n_to_plot <- min(20, length(sig_genes))
  if (n_to_plot > 0) {
    pdf(file.path(output_dir, "isoformswitch_top_gene_plots.pdf"), width = 12, height = 8)
    for (g in sig_genes[1:n_to_plot]) {
      for (comp_i in seq_len(nrow(myComparisons))) {
        tryCatch({
          switchPlot(
            aSwitchList,
            gene = g,
            condition1 = myComparisons$condition_1[comp_i],
            condition2 = myComparisons$condition_2[comp_i],
            plotTopology = FALSE
          )
        }, error = function(e) NULL)
      }
    }
    dev.off()
    cat(paste0("Switch plots saved for top ", n_to_plot, " genes.\n"))
  }
}, error = function(e) cat(paste0("Could not generate switch plots: ", e$message, "\n")))

# ── 11. Save full isoform features table ──
cat("\nSaving full isoform features table...\n")
write.table(aSwitchList$isoformFeatures,
            file = file.path(output_dir, "isoformswitch_all_isoform_features.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Save the full object for potential downstream use
tryCatch({
  saveRDS(aSwitchList, file = file.path(output_dir, "switchAnalyzeRlist.rds"))
  cat("Full switchAnalyzeRlist saved as RDS.\n")
}, error = function(e) cat(paste0("Could not save RDS: ", e$message, "\n")))

# ── 12. Overlap between comparisons ──
if (nrow(myComparisons) > 1) {
  tryCatch({
    pdf(file.path(output_dir, "isoformswitch_overlap.pdf"), width = 8, height = 6)
    extractSwitchOverlap(aSwitchList, filterForConsequences = FALSE, plotIsoforms = FALSE)
    dev.off()
    cat("Overlap plot saved.\n")
  }, error = function(e) cat(paste0("Could not generate overlap plot: ", e$message, "\n")))
}

cat(paste0("\n\nIsoformSwitchAnalyzeR analysis completed. Current date: ", date(), "\n"))
