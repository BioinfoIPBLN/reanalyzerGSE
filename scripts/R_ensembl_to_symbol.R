#!/usr/bin/env Rscript
# R_ensembl_to_symbol.R — Shared helper for ENSEMBL ID detection and conversion
# Source this file from other R scripts that need ENSEMBL → Symbol mapping.

#' Detect if gene IDs are ENSEMBL format (ENSG, ENST, ENSMUSG, etc.)
#' Returns TRUE if the majority (>50%) of non-NA IDs match the pattern.
is_ensembl_id <- function(ids) {
  ids <- ids[!is.na(ids) & nchar(ids) > 0]
  if (length(ids) == 0) return(FALSE)
  # Match both gene (G) and transcript (T) ENSEMBL IDs, with optional version suffix
  pattern <- "^ENS[A-Z]*[GT]\\d+(\\.\\d+)?$"
  mean(grepl(pattern, ids, ignore.case = TRUE)) > 0.5
}

#' Strip ENSEMBL version suffix (.N) from IDs
#' E.g., "ENSG00000172399.7" -> "ENSG00000172399"
strip_ensembl_version <- function(ids) {
  sub("\\.\\d+$", "", ids)
}

#' Build ENSEMBL-to-Symbol mapping from a GTF/GFF annotation file.
#' Extracts gene_id and gene_name from the "gene" feature rows of the GTF.
#' Returns a named character vector: names = ensembl_id (no version), values = gene_name
build_gtf_ensembl_map <- function(gtf_file) {
  if (is.null(gtf_file) || !file.exists(gtf_file)) return(NULL)
  
  gtf <- tryCatch(
    data.table::fread(gtf_file, header = FALSE, sep = "\t",
                      quote = "", comment.char = "#", fill = TRUE,
                      select = c(3L, 9L)),
    error = function(e) NULL
  )
  if (is.null(gtf) || nrow(gtf) == 0) return(NULL)
  
  # Prefer "gene" feature type for one-per-gene mapping; fall back to "exon" or all
  gene_rows <- gtf[gtf$V3 == "gene", ]
  if (nrow(gene_rows) == 0) gene_rows <- gtf
  v9 <- gene_rows$V9
  
  # Extract gene_id
  gene_id <- rep(NA_character_, length(v9))
  rx_id <- regexpr('gene_id\\s+"([^"]*)"', v9, perl = TRUE)
  has_id <- rx_id > 0
  if (any(has_id)) {
    gene_id[has_id] <- sub('.*gene_id\\s+"([^"]*)".*', "\\1", v9[has_id])
  } else {
    # GFF3 fallback
    rx_id2 <- regexpr("(?:^|;)\\s*gene_id=([^;]*)", v9, perl = TRUE)
    has_id2 <- rx_id2 > 0
    if (any(has_id2)) {
      gene_id[has_id2] <- sub(".*gene_id=([^;]*).*", "\\1", v9[has_id2])
    }
  }
  
  # Extract gene_name
  gene_name <- rep(NA_character_, length(v9))
  rx_name <- regexpr('gene_name\\s+"([^"]*)"', v9, perl = TRUE)
  has_name <- rx_name > 0
  if (any(has_name)) {
    gene_name[has_name] <- sub('.*gene_name\\s+"([^"]*)".*', "\\1", v9[has_name])
  } else {
    # GFF3 fallback: gene_name= or Name=
    rx_name2 <- regexpr("(?:^|;)\\s*(?:gene_name|Name)=([^;]*)", v9, perl = TRUE)
    has_name2 <- rx_name2 > 0
    if (any(has_name2)) {
      gene_name[has_name2] <- sub(".*(?:gene_name|Name)=([^;]*).*", "\\1", v9[has_name2], perl = TRUE)
    }
  }
  
  # Build mapping: strip version from gene_id, deduplicate
  valid <- !is.na(gene_id) & !is.na(gene_name) & nchar(gene_name) > 0
  if (sum(valid) == 0) return(NULL)
  
  mapping <- setNames(gene_name[valid], strip_ensembl_version(gene_id[valid]))
  mapping <- mapping[!duplicated(names(mapping))]
  return(mapping)
}

#' Convert ENSEMBL IDs to gene symbols.
#' Strategy: 1) GTF mapping (preferred), 2) org.*.eg.db ENSEMBL keytype (fallback)
#' Returns a character vector of the same length as `ids`, with symbols where possible.
#' Unmapped IDs are returned as-is.
ensembl_to_symbol <- function(ids, gtf_file = NULL, orgdb_name = NULL) {
  if (length(ids) == 0) return(ids)
  
  original_ids <- ids
  ids_stripped <- strip_ensembl_version(ids)
  result <- ids  # Default: return original if no mapping found
  
  # Strategy 1: GTF-based mapping
  gtf_map <- build_gtf_ensembl_map(gtf_file)
  if (!is.null(gtf_map)) {
    # Case-insensitive lookup
    mapped <- gtf_map[toupper(ids_stripped)]
    found <- !is.na(mapped)
    if (any(found)) {
      result[found] <- unname(mapped[found])
      cat(paste0("  ENSEMBL→Symbol via GTF: mapped ", sum(found), "/", length(ids), " IDs\n"))
    }
    # If all mapped, we're done
    if (all(found)) return(result)
  }
  
  # Strategy 2: org.*.eg.db fallback for remaining unmapped
  if (!is.null(orgdb_name)) {
    unmapped_idx <- which(result == original_ids)  # Still unmapped
    if (length(unmapped_idx) > 0) {
      tryCatch({
        suppressMessages(library(orgdb_name, character.only = TRUE, quietly = TRUE))
        orgdb <- eval(parse(text = orgdb_name))
        db_result <- suppressMessages(
          AnnotationDbi::select(orgdb,
                                keys = unique(ids_stripped[unmapped_idx]),
                                columns = "SYMBOL",
                                keytype = "ENSEMBL")
        )
        if (!is.null(db_result) && nrow(db_result) > 0) {
          db_map <- setNames(db_result$SYMBOL, db_result$ENSEMBL)
          db_map <- db_map[!is.na(db_map) & !duplicated(names(db_map))]
          db_mapped <- db_map[ids_stripped[unmapped_idx]]
          db_found <- !is.na(db_mapped)
          if (any(db_found)) {
            result[unmapped_idx[db_found]] <- unname(db_mapped[db_found])
            cat(paste0("  ENSEMBL→Symbol via ", orgdb_name, ": mapped ", sum(db_found), 
                        " additional IDs\n"))
          }
        }
      }, error = function(e) {
        cat(paste0("  Warning: org.db fallback failed: ", e$message, "\n"))
      })
    }
  }
  
  return(result)
}
