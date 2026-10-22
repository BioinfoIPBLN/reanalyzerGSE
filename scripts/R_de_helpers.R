rgse_de_settings <- function() {
  lfc <- suppressWarnings(as.numeric(Sys.getenv("RGSE_LFC_THRESHOLD", "0")))
  if (is.na(lfc) || lfc < 0) lfc <- 0
  list(edger_test = tolower(Sys.getenv("RGSE_EDGER_TEST", "ql")),
       lfc = lfc,
       shrink = tolower(Sys.getenv("RGSE_LFC_SHRINK", "no")))
}

rgse_parse_covariates <- function(covariab, covariab_format) {
  if (is.null(covariab) || is.na(covariab) || covariab %in% c("none", "")) return(NULL)
  vals <- lapply(strsplit(as.character(covariab), ";", fixed = TRUE)[[1]], function(v) trimws(strsplit(v, ",", fixed = TRUE)[[1]]))
  fmts <- rep_len(trimws(strsplit(as.character(covariab_format), ";", fixed = TRUE)[[1]]), length(vals))
  n <- lengths(vals)
  if (length(unique(n)) != 1) stop("Every covariable given in -C must have one value per sample; got ", paste(n, collapse = ", "), " values")
  covs <- as.data.frame(lapply(seq_along(vals), function(k) if (fmts[k] == "fact") factor(vals[[k]]) else as.numeric(vals[[k]])))
  colnames(covs) <- paste0("Cov", seq_along(vals))
  if (any(sapply(covs, function(x) is.numeric(x) && anyNA(x)))) stop("A covariable declared 'num' in -Cf has non-numeric values: ", covariab)
  covs
}

rgse_subset_covariates <- function(covs, keep) {
  if (is.null(covs)) return(NULL)
  if (length(keep) == nrow(covs)) covs <- covs[keep, , drop = FALSE]
  covs
}

rgse_covariates_for_batch_removal <- function(covs) {
  if (is.null(covs)) return(NULL)
  covs <- droplevels(covs)
  m <- model.matrix(~ ., data = covs)
  m[, colnames(m) != "(Intercept)", drop = FALSE]
}

rgse_covariate_design <- function(Treat, covs, label = "") {
  Treat <- droplevels(as.factor(Treat))
  notes <- character(0)
  if (!is.null(covs)) {
    if (nrow(covs) != length(Treat)) stop("The covariables have ", nrow(covs), " values but the design has ", length(Treat), " samples", if (nzchar(label)) paste0(" (", label, ")"))
    covs <- droplevels(covs)
    const <- sapply(covs, function(x) length(unique(x)) < 2)
    if (any(const)) {
      notes <- c(notes, paste0("constant in these samples, not used: ", paste(colnames(covs)[const], collapse = ", ")))
      covs <- covs[, !const, drop = FALSE]
    }
  }
  build <- function(cv) {
    if (is.null(cv) || ncol(cv) == 0) return(model.matrix(~0 + Treat))
    model.matrix(as.formula(paste("~0 + Treat +", paste(colnames(cv), collapse = " + "))), data = cbind(data.frame(Treat = Treat), cv))
  }
  design <- build(covs)
  while (!is.null(covs) && ncol(covs) > 0 && qr(design)$rank < ncol(design)) {
    notes <- c(notes, paste0("confounded with the condition or with another covariable, not used: ", colnames(covs)[ncol(covs)]))
    covs <- covs[, -ncol(covs), drop = FALSE]
    design <- build(covs)
  }
  if (!is.null(covs) && ncol(covs) == 0) covs <- NULL
  list(design = design, covs = covs, notes = notes)
}

rgse_contrast <- function(design, comp) {
  idx <- match(paste0("Treat", sub("^__", "", comp)), colnames(design))
  if (length(idx) != 2 || anyNA(idx)) {
    idx <- match(paste0("Treat", comp), colnames(design))
  }
  if (length(idx) != 2 || anyNA(idx)) stop("The comparison ", paste(comp, collapse = " vs "), " could not be matched to the design columns: ", paste(colnames(design), collapse = ", "))
  contrast <- rep(0, ncol(design))
  contrast[idx] <- c(1, -1)
  contrast
}

rgse_de_method_label <- function(engine, settings, covs) {
  cov_txt <- if (is.null(covs)) "" else paste0(" + ", paste(colnames(covs), collapse = " + "))
  lfc_txt <- if (settings$lfc > 0) paste0(", H0: |log2FC| <= ", settings$lfc) else ""
  switch(engine,
         edgeR = if (settings$edger_test == "exact" && is.null(covs) && settings$lfc == 0) "edgeR exactTest (~ condition)"
                 else paste0("edgeR quasi-likelihood ", if (settings$lfc > 0) "glmTreat" else "F-test", " (~0 + condition", cov_txt, lfc_txt, ")"),
         DESeq2 = paste0("DESeq2 Wald test (~", sub("^ \\+ ", " ", paste0(cov_txt, " + condition")), lfc_txt,
                         if (settings$shrink == "yes") paste0(", shrunkLFC by lfcShrink type ", if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal") else "", ")"),
         limma = paste0("limma-voom moderated ", if (settings$lfc > 0) "treat" else "t-test", " (~0 + condition", cov_txt, lfc_txt, ", robust eBayes)"),
         wilcoxon = "Wilcoxon rank-sum test on TMM-normalised CPM (covariables and the fold-change threshold do not apply)",
         engine)
}

rgse_de_test <- function(y, comp, covs = NULL, engine = "edgeR", settings = rgse_de_settings()) {
  groups <- as.character(y$samples$group)
  cd <- rgse_covariate_design(groups, covs)
  covs <- cd$covs
  design <- cd$design
  rownames(design) <- colnames(y)
  out_cols <- function(logFC, logCPM, stat, stat_name, p, fdr) {
    tab <- data.frame(Gene_ID = y$genes$Gene_ID, Length = y$genes$Length, logFC = logFC, logCPM = logCPM, stringsAsFactors = FALSE)
    if (!is.null(stat)) tab[[stat_name]] <- stat
    tab$PValue <- p
    tab$FDR <- fdr
    rownames(tab) <- rownames(y)
    tab
  }
  if (engine == "edgeR") {
    if (settings$edger_test == "exact" && is.null(covs) && settings$lfc == 0) {
      if (is.null(y$common.dispersion)) y <- edgeR::estimateDisp(y, robust = TRUE)
      if (is.na(y$common.dispersion)) y$common.dispersion <- 0.4^2
      tab <- edgeR::topTags(edgeR::exactTest(y, pair = rev(comp)), n = nrow(y), adjust.method = "BH", sort.by = "PValue")$table
    } else {
      if (ncol(y) - ncol(design) < 1) stop("no residual degrees of freedom for the quasi-likelihood fit")
      fit <- edgeR::glmQLFit(y, design, robust = TRUE)
      contrast <- rgse_contrast(design, comp)
      res <- if (settings$lfc > 0) edgeR::glmTreat(fit, contrast = contrast, lfc = settings$lfc) else edgeR::glmQLFTest(fit, contrast = contrast)
      tab <- edgeR::topTags(res, n = nrow(y), adjust.method = "BH", sort.by = "PValue")$table
      tab$unshrunk.logFC <- NULL
    }
  } else if (engine == "DESeq2") {
    counts_mat <- round(as.matrix(y$counts)); mode(counts_mat) <- "integer"
    coldata <- data.frame(condition = factor(sub("^__", "", groups)), row.names = colnames(counts_mat))
    if (!is.null(covs)) coldata <- cbind(coldata, covs)
    form <- as.formula(paste("~", paste(c(colnames(covs), "condition"), collapse = " + ")))
    dds <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = form), quiet = TRUE)
    ctr <- c("condition", sub("^__", "", comp[1]), sub("^__", "", comp[2]))
    res <- DESeq2::results(dds, contrast = ctr, lfcThreshold = settings$lfc, altHypothesis = "greaterAbs")
    tab <- out_cols(res$log2FoldChange, edgeR::aveLogCPM(y), res$stat, "stat", res$pvalue, res$padj)
    tab$PValue[is.na(tab$PValue)] <- 1
    tab$FDR[is.na(tab$FDR)] <- 1
    if (settings$shrink == "yes") {
      shr_type <- if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal"
      shr <- suppressMessages(DESeq2::lfcShrink(dds, contrast = ctr, res = res, type = shr_type, quiet = TRUE))
      tab$shrunkLFC <- shr$log2FoldChange
    }
    tab <- tab[order(tab$PValue), ]
  } else if (engine == "limma") {
    if (ncol(y) - ncol(design) < 1) stop("no residual degrees of freedom for limma")
    v <- limma::voom(y, design)
    fit <- limma::contrasts.fit(limma::lmFit(v, design), rgse_contrast(design, comp))
    if (settings$lfc > 0) {
      fit <- limma::treat(fit, lfc = settings$lfc, robust = TRUE)
      tt <- limma::topTreat(fit, coef = 1, number = Inf, sort.by = "none")
    } else {
      fit <- limma::eBayes(fit, robust = TRUE)
      tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
    }
    tab <- out_cols(tt$logFC, tt$AveExpr, tt$t, "t", tt$P.Value, tt$adj.P.Val)
    tab <- tab[order(tab$PValue), ]
  } else if (engine == "wilcoxon") {
    g1 <- groups == comp[1]; g2 <- groups == comp[2]
    if (min(sum(g1), sum(g2)) < 8) cat(paste0("WARNING: the Wilcoxon rank-sum test with ", sum(g1), " vs ", sum(g2), " samples has little power and coarse p-values; it is meant for about 8 or more samples per condition\n"))
    cpm_m <- edgeR::cpm(y)
    lcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
    w <- t(apply(cpm_m, 1, function(x) {
      r <- suppressWarnings(wilcox.test(x[g1], x[g2]))
      c(r$statistic, r$p.value)
    }))
    p <- w[, 2]; p[is.na(p)] <- 1
    tab <- out_cols(rowMeans(lcpm[, g1, drop = FALSE]) - rowMeans(lcpm[, g2, drop = FALSE]), edgeR::aveLogCPM(y), w[, 1], "W", p, p.adjust(p, "BH"))
    tab <- tab[order(tab$PValue), ]
  } else {
    stop("Unknown differential expression engine: ", engine)
  }
  list(table = tab, method = rgse_de_method_label(engine, settings, covs), notes = cd$notes)
}
