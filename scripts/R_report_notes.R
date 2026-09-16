.rgse_capping_note <- function(dest_dir, analysis, subject, kept, total, reason) {
  if (length(dest_dir) != 1 || !nzchar(dest_dir) || !dir.exists(dest_dir)) return(invisible(NULL))
  line <- paste0(paste(analysis, subject, kept, total, reason, sep = "\t"), "\n")
  tryCatch(cat(line, file = file.path(dest_dir, "analysis_capping_notes.txt"), append = TRUE),
           error = function(e) invisible(NULL))
  invisible(NULL)
}
