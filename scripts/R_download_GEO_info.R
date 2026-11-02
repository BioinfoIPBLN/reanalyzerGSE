#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GEO_ID <- args[1]
path <- args[2]

GEO_ID_path <- unlist(lapply(strsplit(GEO_ID,","),function(x){paste(sort(x),collapse="_")}))
cat(paste0("\nI'm in: /",path,"/",GEO_ID_path,"/reads_study_info\n\n"))
dir.create(paste(path,GEO_ID_path,"reads_study_info",sep="/"), showWarnings = FALSE, recursive = TRUE); setwd(paste(path,GEO_ID_path,"reads_study_info",sep="/"))

retry_delays <- c(0, 10, 30)

http_status <- function(url) {
  code <- suppressWarnings(system2("curl", c("-s", "-o", "/dev/null", "-w", shQuote("%{http_code}"), "--max-time", "60", shQuote(url)), stdout = TRUE, stderr = FALSE))
  suppressWarnings(as.integer(tail(code, 1)))
}

wget_with_retry <- function(dir, url) {
  for (attempt in seq_along(retry_delays)) {
    if (retry_delays[attempt] > 0) Sys.sleep(retry_delays[attempt])
    if (system(paste0("cd ", dir, " && wget -nv ", url)) == 0) return(TRUE)
    cat(paste0("Download attempt ", attempt, " failed for: ", url, "\n"))
  }
  return(FALSE)
}

stop_download <- function(GEO_ID, url) {
  cat(paste0("\nERROR: the GEO metadata of ", GEO_ID, " could not be downloaded from NCBI.\n",
             "If NCBI answers 403, its file server is refusing this network for a while, usually after too many requests: wait and rerun later. This prints 200 once access is back:\n",
             "  curl -s -o /dev/null -w '%{http_code}\\n' ", url, "\n",
             "If it answers 404, check the GEO ID: a private or very recent series has no files there yet.\n\n"))
  quit(save = "no", status = 1)
}

for (i in unlist(strsplit(GEO_ID_path,"_"))){
	GEO_ID <- i
	nnn_id <- paste0(substr(GEO_ID,1,nchar(GEO_ID)-3),"nnn")
	series_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", nnn_id, "/", GEO_ID, "/")
	target_dir <- paste(path,GEO_ID_path,"reads_study_info",sep="/")
	success <- FALSE
	refused <- FALSE
	for (attempt in seq_along(retry_delays)) {
		if (retry_delays[attempt] > 0) Sys.sleep(retry_delays[attempt])
		code <- http_status(series_url)
		refused <- !is.na(code) && code %in% c(403, 404)
		if (refused) {
			cat(paste0("NCBI answered HTTP ", code, " for ", series_url, " (attempt ", attempt, ")\n"))
			next
		}
		tryCatch({
			a <- GEOquery::getGEO(GEO_ID, destdir = target_dir)
			GEOquery::getGEOSuppFiles(GEO_ID, baseDir = target_dir)
			success <- TRUE
		}, error=function(e){
			cat(paste0("GEOquery attempt ", attempt, " failed for ", GEO_ID, ": ", e$message, "\n"))
		})
		if (success) break
	}
	if (refused) stop_download(GEO_ID, series_url)
	if (!success) {
		cat(paste0("GEOquery failed after ", length(retry_delays), " attempts for ", GEO_ID, ". Falling back to direct wget with retries...\n"))
		system(paste0("cd ", target_dir, " && rm -rf *"))
		wget_with_retry(target_dir, paste0(series_url, "matrix/", GEO_ID, "_series_matrix.txt.gz"))
		wget_with_retry(target_dir, paste0(series_url, "soft/", GEO_ID, "_family.soft.gz"))
		suppl_dir <- paste(path,GEO_ID_path,"reads_study_info",GEO_ID,sep="/")
		system(paste0("mkdir -p ", suppl_dir))
		if (wget_with_retry(suppl_dir, paste0(series_url, "suppl/"))) {
			system(paste0('cd ', suppl_dir, ' && for f in $(grep href index.html 2>/dev/null | grep -v www. | grep -v // | grep -v nnn | cut -d ">" -f 2 | cut -d "<" -f 1); do wget -nv ', series_url, 'suppl/$f; done'))
		}
	}
	if (length(list.files(target_dir, pattern = paste0("^", GEO_ID, ".*_series_matrix\\.txt\\.gz$"))) == 0) stop_download(GEO_ID, series_url)
}

# All the info... for example:
# a[[1]]@experimentData@other$title
# a[[1]]@phenoData@data$title
# a[[1]]@phenoData@data$characteristics_ch1.3
# a[[1]]@phenoData@data$relation.1
cat("\nDownloading from GEO done\n\n")

