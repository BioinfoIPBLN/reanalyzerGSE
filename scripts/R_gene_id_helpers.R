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
      lut <- stats::setNames(syms, tolower(syms))
      if ("ALIAS" %in% AnnotationDbi::keytypes(db)) {
        ali <- AnnotationDbi::keys(db, keytype = "ALIAS")
        ali <- ali[!tolower(ali) %in% names(lut)]
        if (length(ali) > 0) {
          amap <- suppressWarnings(suppressMessages(AnnotationDbi::select(
            db, keys = ali, columns = "SYMBOL", keytype = "ALIAS")))
          amap <- amap[!is.na(amap$SYMBOL) & !duplicated(tolower(amap$ALIAS)), , drop = FALSE]
          lut <- c(lut, stats::setNames(amap$SYMBOL, tolower(amap$ALIAS)))
        }
      }
      list(official = syms, lut = lut[!duplicated(names(lut))])
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
  out
}
