#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
cores <- as.integer(args[2])
if (is.na(cores) || cores < 1L) cores <- 1L
pattern_to_avoid <- as.character(args[3])

MAX_CHARS <- 32767L   # Excel's hard per-cell character limit

# Cheap, vectorised check: does any cell exceed Excel's limit and therefore need
# the (expensive) row-splitting below? Only text columns can, so skip the rest.
needs_split <- function(df) {
  for (j in seq_along(df)) {
    x <- df[[j]]
    if (is.character(x) || is.factor(x)) {
      if (any(nchar(as.character(x)) > MAX_CHARS, na.rm = TRUE)) return(TRUE)
    }
  }
  FALSE
}

split_large_cells <- function(df, max_chars = MAX_CHARS) {
  # Create an empty dataframe with the same structure as the original
  new_df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
  colnames(new_df) <- colnames(df)

  for (i in 1:nrow(df)) {
    row_data <- df[i, ]
    split_rows <- list(row_data)
    
    for (j in 1:ncol(df)) {
      cell_content <- as.character(row_data[[j]])
      if (!is.na(cell_content) && nchar(cell_content) > max_chars) {
        # Split the cell content into chunks
        chunks <- strsplit(cell_content, paste0("(?<=.{" , max_chars , "})"), perl = TRUE)[[1]]
        split_rows[[1]][[j]] <- chunks[1]
        for (k in 2:length(chunks)) {
          # If there aren't enough rows in split_rows, add one
          if (length(split_rows) < k) {
            split_rows[[k]] <- row_data
            split_rows[[k]][-j] <- NA  # Set other columns to NA
          }
          split_rows[[k]][[j]] <- chunks[k]
        }
      }
    }
    new_df <- rbind(new_df, do.call(rbind, split_rows))
  }

  return(new_df)
}

suppressMessages(library(parallel,quiet = T,warn.conflicts = F))
print("Writing xlsx files...")


files=grep(pattern_to_avoid,list.files(path=path,pattern = "\\.txt$", full.names = TRUE,recursive=T),val=T,invert=T)

process_file <- function(file){
  data <- data.table::fread(file,sep="\t")
  # Check if the table has more than one column
  if (ncol(data) > 1) {
    data <- as.data.frame(data)
    data[is.na(data) | data == ""] <- "-"
    # Only pay for the O(n^2) row-splitting when a cell actually overflows Excel's
    # limit; count/annotation matrices never do, so they skip it entirely.
    data_final <- if (needs_split(data)) split_large_cells(data) else data
    writexl::write_xlsx(as.data.frame(data_final), paste0(sub("\\.txt$", "", file), ".xlsx"))
    #unlink(file)
  }
}

mclapply(
    mc.cores = cores,
    X = files,
    FUN = process_file
)
