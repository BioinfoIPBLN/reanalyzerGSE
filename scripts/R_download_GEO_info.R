#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GEO_ID <- args[1]
path <- args[2]

GEO_ID_path <- unlist(lapply(strsplit(GEO_ID,","),function(x){paste(sort(x),collapse="_")}))
cat(paste0("\nI'm in: /",path,"/",GEO_ID_path,"/reads_study_info\n\n"))
dir.create(paste(path,GEO_ID_path,"reads_study_info",sep="/"), showWarnings = FALSE, recursive = TRUE); setwd(paste(path,GEO_ID_path,"reads_study_info",sep="/"))

wget_with_retry <- function(url_cmd) {
  delays <- c(0, 10, 20, 30)
  for (attempt in 1:4) {
    if (delays[attempt] > 0) Sys.sleep(delays[attempt])
    res <- system(url_cmd)
    if (res == 0) return(TRUE)
    cat(paste0("Download attempt ", attempt, " failed for: ", url_cmd, ". Retrying...\n"))
  }
  return(FALSE)
}

for (i in unlist(strsplit(GEO_ID_path,"_"))){
	GEO_ID <- i
	success <- FALSE
	delays <- c(0, 10, 20, 30)
	for (attempt in 1:4) {
		if (delays[attempt] > 0) Sys.sleep(delays[attempt])
		tryCatch({
			a <- GEOquery::getGEO(GEO_ID, destdir = paste(path,GEO_ID_path,"reads_study_info",sep="/"))
			GEOquery::getGEOSuppFiles(GEO_ID, baseDir = paste(path,GEO_ID_path,"reads_study_info",sep="/"))
			success <- TRUE
		}, error=function(e){
			cat(paste0("GEOquery attempt ", attempt, " failed for ", GEO_ID, ": ", e$message, "\n"))
		})
		if (success) break
	}
	if (!success) {
		cat(paste0("GEOquery failed after 4 attempts for ", GEO_ID, ". Falling back to direct wget with retries...\n"))
		target_dir <- paste(path,GEO_ID_path,"reads_study_info",sep="/")
		system(paste0("cd ", target_dir, " && rm -rf *"))
		nnn_id <- paste0(substr(GEO_ID,1,nchar(GEO_ID)-3),"nnn")
		wget_with_retry(paste0("cd ", target_dir, " && wget https://ftp.ncbi.nlm.nih.gov/geo/series/", nnn_id, "/", GEO_ID, "/matrix/", GEO_ID, "_series_matrix.txt.gz"))
		wget_with_retry(paste0("cd ", target_dir, " && wget https://ftp.ncbi.nlm.nih.gov/geo/series/", nnn_id, "/", GEO_ID, "/soft/", GEO_ID, "_family.soft.gz"))
		suppl_dir <- paste(path,GEO_ID_path,"reads_study_info",GEO_ID,sep="/")
		system(paste0("mkdir -p ", suppl_dir))
		wget_with_retry(paste0("cd ", suppl_dir, " && wget https://ftp.ncbi.nlm.nih.gov/geo/series/", nnn_id, "/", GEO_ID, "/suppl/"))
		system(paste0('cd ', suppl_dir, ' && for f in $(grep href index.html 2>/dev/null | grep -v www. | grep -v // | grep -v nnn | cut -d ">" -f 2 | cut -d "<" -f 1); do wget https://ftp.ncbi.nlm.nih.gov/geo/series/', nnn_id, '/', GEO_ID, '/suppl/$f; done'))
	}
}

# All the info... for example:
# a[[1]]@experimentData@other$title
# a[[1]]@phenoData@data$title
# a[[1]]@phenoData@data$characteristics_ch1.3
# a[[1]]@phenoData@data$relation.1
cat("\nDownloading from GEO done\n\n")

