#!/usr/bin/env Rscript
################################################################################
#### Shared helpers: resolve the organism annotation and canonicalise gene IDs
################################################################################

.rgse_orgdb_ref <- local({
  cached <- NULL
  resolved <- FALSE
  function() {
    if (resolved) return(cached)
    resolved <<- TRUE
    pkg <- NULL
    if (exists("orgDB", envir = globalenv(), inherits = FALSE)) {
      v <- get("orgDB", envir = globalenv())
      if (is.character(v) && length(v) == 1L && nzchar(v)) pkg <- v
    }
    if (is.null(pkg) && exists("organism", envir = globalenv(), inherits = FALSE)) {
      o <- get("organism", envir = globalenv())
      if (is.character(o) && length(o) == 1L) {
        pkg <- if (grepl("Mus", o, ignore.case = TRUE)) "org.Mm.eg.db"
               else if (grepl("Homo", o, ignore.case = TRUE)) "org.Hs.eg.db"
               else if (grepl("Rattus", o, ignore.case = TRUE)) "org.Rn.eg.db"
               else NULL
      }
    }
    if (is.null(pkg) || !requireNamespace(pkg, quietly = TRUE)) return(cached)
    cached <<- tryCatch(get(pkg, envir = asNamespace(pkg)), error = function(e) NULL)
    cached
  }
})

.rgse_symbol_lut <- local({
  cached <- NULL
  resolved <- FALSE
  function() {
    if (resolved) return(cached)
    resolved <<- TRUE
    db <- .rgse_orgdb_ref()
    if (is.null(db)) return(cached)
    cached <<- tryCatch({
      syms <- AnnotationDbi::keys(db, keytype = "SYMBOL")
      list(official = syms, lut = stats::setNames(syms, tolower(syms)))
    }, error = function(e) NULL)
    cached
  }
})

canonicalise_gene_ids <- function(ids, fallback = NULL) {
  if (length(ids) == 0) return(ids)
  if (is.null(fallback)) fallback <- ids
  ref <- .rgse_symbol_lut()
  if (is.null(ref)) return(fallback)
  out <- ids
  todo <- !(ids %in% ref$official)
  if (any(todo)) {
    hit <- ref$lut[tolower(ids[todo])]
    out[todo] <- ifelse(is.na(hit), fallback[todo], unname(hit))
  }
  changed <- !is.na(out) & !is.na(ids) & out != ids
  if (any(changed)) {
    collide <- changed & (out %in% out[duplicated(out)])
    if (any(collide)) out[collide] <- ids[collide]
  }
  out
}
