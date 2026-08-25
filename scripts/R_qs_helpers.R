#!/usr/bin/env Rscript

local({
  n <- suppressWarnings(as.integer(Sys.getenv("RGSE_QS_THREADS", unset = "")))
  if (length(n) != 1 || is.na(n) || n < 1L) n <- 4L
  qs2::qopt("nthreads", as.integer(min(n, 8L)))
})

rgse_save_env <- function(file, envir = .GlobalEnv) {
  qs2::qs_save(as.list(envir, all.names = TRUE), file)
  invisible(file)
}

rgse_load_env <- function(file, envir = .GlobalEnv) {
  list2env(qs2::qs_read(file), envir = envir)
  invisible(envir)
}

# Canonical sample id from any per-sample artifact (BAM, flagstat, samtools stats,
# featureCounts .tab, count matrix colname, targets entry). Set mate_suffix=TRUE only
# for per-mate artifacts (fastq/fastqc file names), where a trailing _1/_2 identifies
# the read mate rather than the sample: every aligner/counting artifact is named from
# the mate-less prefix, so stripping it there eats real sample names (e.g. YP0_1).
rgse_sample_id <- function(x, mate_suffix = FALSE) {
  x <- basename(as.character(x))
  x <- sub("_categ(_2)?$", "", x)
  x <- sub("_readcount\\.tab(\\.summary)?$", "", x)
  x <- sub("_(hisat2|HISAT2|star|STAR)([._].*)?$", "", x)
  x <- sub("_(flagstat|stats)(\\..*)?$", "", x)
  x <- sub("_fastqc(\\..*)?$", "", x)
  x <- sub("\\.bam$", "", x)
  x <- sub("\\.(fastq|fq)(\\.gz)?$", "", x)
  if (mate_suffix) x <- sub("_R?[12]$", "", x)
  x
}

# Warn and keep the un-normalized names for entries that collided, rather than
# disambiguating them (make.unique) and silently mislabelling samples.
rgse_resolve_collisions <- function(clean, raw, what = "sample names") {
  if (!anyDuplicated(clean)) return(clean)
  dups <- unique(clean[duplicated(clean)])
  bad <- clean %in% dups
  cat("\nWARNING: normalization of ", what, " collapsed distinct samples into: ",
      paste(dups, collapse = ", "),
      "\n  Keeping the un-normalized name for: ",
      paste(raw[bad], collapse = ", "), "\n", sep = "")
  clean[bad] <- raw[bad]
  if (anyDuplicated(clean)) {
    stop("The ", what, " are still not unique after falling back to the raw names: ",
         paste(unique(clean[duplicated(clean)]), collapse = ", "))
  }
  clean
}
