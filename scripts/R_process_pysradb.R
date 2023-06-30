#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GEO_ID <- args[1]
path <- args[2]
suppressMessages(library(dplyr,quiet = T,warn.conflicts = F))


# Check if there are same headers if there are multiple pysradb files and rbind if yes:
setwd(paste0(path,"/",GEO_ID,"/GEO_info"))

a <- lapply(as.list(list.files(pattern="pysradb")),
            function(x){as.data.frame(data.table::fread(x))})
a <- lapply(a,function(x){x[,gtools::mixedsort(colnames(x))]})
b <- lapply(a,colnames)
if (length(unique(b)) == 1){
	d <- bind_rows(a, .id = "column_label")
	d <- as.data.frame(d[gtools::mixedorder(d[,names(which.min(nchar(grep("^GSM",d[1,],val=T))))[1]]),]) # Ensure the order of the metadata is following the GSMXXXX ids and the SRR ids, because this is the order in the series_matrix tables from which I get more info
	not_all_na <- function(x) any(!is.na(x))
	e <- d %>% select(where(not_all_na))
	# Function to detect that a column is duplicated of other, despite different colnames:
	dup_cols <- sapply(2:ncol(e), function(j) {
	  any(sapply(1:(j - 1), function(i) {
	    isTRUE(all.equal(e[, i], e[, j]))
	  }))
	})
	columns_to_discard <- unique(c("column_label",grep("AWS|ena_fastq|http|url",colnames(e),val=T),names(e)[dup_cols]))
	f <- e[,!(colnames(e) %in% columns_to_discard)]
	
	dir.create(paste0(path,"/",GEO_ID,"/final_results_reanalysis"),showWarnings = FALSE)
	write.table(f,
	    file=paste0(path,"/",GEO_ID,"/final_results_reanalysis/sample_info.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
} else {
	cat("\nWe are dealing with multiple pysradb files and GSEXXX ids here... the merging of the metadata has not worked. Please proceed with caution, and double check all the files, naming, ordering, etc. You probably MUST look manually into this\n\n")
}


