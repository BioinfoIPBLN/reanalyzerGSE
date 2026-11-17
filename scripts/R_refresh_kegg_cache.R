#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "R_net_config.R"))
.rgse_net_configure(verbose = FALSE)
if (!requireNamespace("rappdirs", quietly = TRUE)) quit(save = "no", status = 0)

files <- c("kegg_category.rda", "kegg_species.rda")
urls <- c("https://yulab-smu.top/clusterProfiler", "https://raw.githubusercontent.com/YuLab-SMU/clusterProfiler/gh-pages")
cache_dir <- rappdirs::user_data_dir("clusterProfiler")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
md5_of <- function(f) unname(tolower(tools::md5sum(f)))

md5 <- NULL
for (url in urls) {
  md5 <- tryCatch(suppressWarnings(utils::read.table(paste0(url, "/md5.txt"), sep = "", colClasses = "character", fill = TRUE)), error = function(e) NULL)
  if (!is.null(md5) && ncol(md5) >= 2) break
  md5 <- NULL
}
if (is.null(md5)) {
  cat("The clusterProfiler KEGG tables could not be checked online. Keeping the cached copies in ", cache_dir, "\n", sep = "")
  quit(save = "no", status = 0)
}

for (f in files) {
  remote <- tolower(md5[md5[[1]] == f, 2])
  if (length(remote) != 1) next
  dest <- file.path(cache_dir, f)
  if (file.exists(dest) && identical(md5_of(dest), remote)) next
  tmp <- tempfile(paste0(f, "."), tmpdir = cache_dir)
  done <- FALSE
  for (url in urls) {
    ok <- tryCatch(utils::download.file(paste0(url, "/", f), tmp, mode = "wb", quiet = TRUE) == 0, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && identical(md5_of(tmp), remote) && file.rename(tmp, dest)) {
      done <- TRUE
      break
    }
  }
  unlink(tmp)
  if (done) {
    cat("Updated the clusterProfiler KEGG table ", dest, "\n", sep = "")
  } else {
    cat("Could not update the clusterProfiler KEGG table ", dest, ". Keeping the cached copy, if any\n", sep = "")
  }
}
