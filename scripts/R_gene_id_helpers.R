#!/usr/bin/env Rscript
################################################################################
#### Shared helpers: resolve the organism annotation and canonicalise gene IDs
################################################################################

.rgse_organisms <- data.frame(
  name     = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Drosophila melanogaster", "Caenorhabditis elegans",
               "Saccharomyces cerevisiae", "Arabidopsis thaliana", "Gallus gallus", "Sus scrofa", "Bos taurus", "Canis lupus familiaris",
               "Macaca mulatta", "Pan troglodytes", "Xenopus laevis", "Escherichia coli", "Anopheles gambiae", "Plasmodium falciparum"),
  orgdb    = c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dr.eg.db", "org.Dm.eg.db", "org.Ce.eg.db",
               "org.Sc.sgd.db", "org.At.tair.db", "org.Gg.eg.db", "org.Ss.eg.db", "org.Bt.eg.db", "org.Cf.eg.db",
               "org.Mmu.eg.db", "org.Pt.eg.db", "org.Xl.eg.db", "org.EcK12.eg.db", "org.Ag.eg.db", "org.Pf.plasmo.db"),
  kegg     = c("hsa", "mmu", "rno", "dre", "dme", "cel", "sce", "ath", "gga", "ssc", "bta", "cfa", "mcc", "ptr", "xla", "eco", "aga", "pfa"),
  reactome = c("human", "mouse", "rat", "zebrafish", "fly", "celegans", "yeast", "arabidopsis", "chicken", "pig", "bovine", "canine",
               "rhesus", "chimp", "xenopus", "ecolik12", "anopheles", "malaria"),
  panther  = c(9606, 10090, 10116, 7955, 7227, 6239, 559292, 3702, 9031, 9823, 9913, 9615, 9544, 9598, 8364, 83333, 7165, 36329),
  enrichr  = c("Enrichr", "Enrichr", "Enrichr", "FishEnrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", NA, "Enrichr", "Enrichr",
               "Enrichr", "Enrichr", "Enrichr", "Enrichr", NA, NA, NA, NA),
  multimir = c("hsa", "mmu", "rno", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  symbol_keytype = c(rep("SYMBOL", 6), "GENENAME", rep("SYMBOL", 11)),
  wikipathways = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Drosophila melanogaster", "Caenorhabditis elegans",
               "Saccharomyces cerevisiae", "Arabidopsis thaliana", "Gallus gallus", "Sus scrofa", "Bos taurus", "Canis familiaris",
               NA, "Pan troglodytes", NA, NA, "Anopheles gambiae", NA),
  graphite = c("hsapiens", "mmusculus", "rnorvegicus", "drerio", "dmelanogaster", "celegans", "scerevisiae", NA, "ggallus", "sscrofa",
               "btaurus", "cfamiliaris", NA, NA, NA, NA, NA, NA),
  msigdbr  = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Drosophila melanogaster", "Caenorhabditis elegans",
               "Saccharomyces cerevisiae", NA, "Gallus gallus", "Sus scrofa", "Bos taurus", "Canis lupus familiaris",
               "Macaca mulatta", "Pan troglodytes", "Xenopus tropicalis", NA, NA, NA),
  aliases  = c("human", "mouse", "rat", "zebrafish", "fruit fly|drosophila", "worm|nematode", "yeast|budding yeast", "arabidopsis|thale cress",
               "chicken", "pig|swine", "cow|cattle|bovine", "dog|canis familiaris", "rhesus|rhesus macaque|macaque", "chimpanzee|chimp",
               "frog|xenopus", "e. coli|escherichia coli k-12", "mosquito", "malaria|plasmodium falciparum 3d7"),
  stringsAsFactors = FALSE)

rgse_organism_info <- function(organism) {
  if (is.null(organism) || length(organism) != 1L || is.na(organism)) return(NULL)
  o <- tolower(trimws(gsub("[_[:space:]]+", " ", organism)))
  if (!nzchar(o)) return(NULL)
  tab <- .rgse_organisms
  nm <- tolower(tab$name)
  i <- match(o, nm); how <- "exact"
  if (is.na(i)) { i <- which(startsWith(o, nm))[1]; how <- "prefix" }
  if (is.na(i)) { i <- which(vapply(strsplit(tab$aliases, "|", fixed = TRUE), function(a) o %in% a, logical(1)))[1]; how <- "alias" }
  if (is.na(i)) { j <- which(sub(" .*", "", nm) == sub(" .*", "", o)); if (length(j) == 1L) { i <- j; how <- "genus" } }
  if (is.na(i)) return(NULL)
  info <- as.list(tab[i, , drop = FALSE])
  info$match <- how
  info$input <- organism
  info
}

rgse_symbol_keytype <- function(info) {
  if (is.null(info) || !suppressPackageStartupMessages(suppressMessages(requireNamespace(info$orgdb, quietly = TRUE)))) return(NA_character_)
  db <- tryCatch(get(info$orgdb, envir = asNamespace(info$orgdb)), error = function(e) NULL)
  if (is.null(db)) return(NA_character_)
  kt <- suppressMessages(AnnotationDbi::keytypes(db))
  cand <- unique(c(info$symbol_keytype, "SYMBOL", "GENENAME", "ORF", "TAIR"))
  hit <- cand[cand %in% kt]
  if (length(hit)) hit[1] else NA_character_
}

rgse_orgdb_check <- function(info, need = NULL) {
  if (is.null(info)) return("unknown")
  if (!suppressPackageStartupMessages(suppressMessages(requireNamespace(info$orgdb, quietly = TRUE)))) return("not_installed")
  db <- tryCatch(get(info$orgdb, envir = asNamespace(info$orgdb)), error = function(e) NULL)
  if (is.null(db)) return("not_installed")
  if (is.null(need)) { if (is.na(rgse_symbol_keytype(info))) return("no_symbol") }
  else if (!all(need %in% suppressMessages(AnnotationDbi::keytypes(db)))) return("no_symbol")
  "ok"
}

rgse_orgdb_attach <- function(info) {
  suppressMessages(library(info$orgdb, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
  if (!identical(info$match, "exact")) cat(paste0("Organism '", info$input, "' mapped to ", info$name, " (", info$orgdb, ")\n"))
  info$orgdb
}

rgse_organism_note <- function(organism, info, status) {
  switch(status,
    unknown = paste0("Organism '", organism, "' has no annotation package in reanalyzerGSE (supported: ", paste(.rgse_organisms$name, collapse = ", "), ")."),
    not_installed = paste0("The annotation package ", info$orgdb, " for ", info$name, " is not installed in this container: add it to the ORGDBS build argument of the image, or bind an R library that contains it."),
    no_symbol = paste0(info$orgdb, " for ", info$name, " lacks the identifier type this step needs."),
    "")
}

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
        info <- rgse_organism_info(o)
        pkg <- if (is.null(info)) NULL else info$orgdb
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
      kt <- suppressMessages(AnnotationDbi::keytypes(db))
      kt <- if ("SYMBOL" %in% kt) "SYMBOL" else if ("GENENAME" %in% kt) "GENENAME" else stop("no symbol keytype")
      syms <- AnnotationDbi::keys(db, keytype = kt)
      list(official = syms, lut = stats::setNames(syms, tolower(syms)))
    }, error = function(e) NULL)
    cached
  }
})

canonicalise_gene_ids <- function(ids, fallback = NULL) {
  if (length(ids) == 0) return(ids)
  if (identical(Sys.getenv("RGSE_GENE_ID_CANONICALISE", "yes"), "no")) return(if (is.null(fallback)) ids else fallback)
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
