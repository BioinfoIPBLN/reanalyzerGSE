#!/usr/bin/env Rscript

local({
  n <- suppressWarnings(as.integer(Sys.getenv("RGSE_QS_THREADS", unset = "")))
  if (length(n) != 1 || is.na(n) || n < 1L) n <- 4L
  qs2::qopt(nthreads = as.integer(min(n, 8L)))
})

rgse_save_env <- function(file, envir = .GlobalEnv) {
  qs2::qs_save(as.list(envir, all.names = TRUE), file)
  invisible(file)
}

rgse_load_env <- function(file, envir = .GlobalEnv) {
  list2env(qs2::qs_read(file), envir = envir)
  invisible(envir)
}
