#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GEO_ID <- args[1]
path <- args[2]

GEO_ID_path <- unlist(lapply(strsplit(GEO_ID,","),function(x){paste(sort(x),collapse="_")}))
cat(paste0("\nI'm in: /",path,"/",GEO_ID_path,"/GEO_info\n\n"))
dir.create(paste(path,GEO_ID_path,"GEO_info",sep="/"), showWarnings = FALSE, recursive = TRUE); setwd(paste(path,GEO_ID_path,"GEO_info",sep="/"))

for (i in unlist(strsplit(GEO_ID_path,"_"))){
	GEO_ID <- i
	tryCatch({
		a <- GEOquery::getGEO(GEO_ID,destdir = paste(path,GEO_ID_path,"GEO_info",sep="/"))
		GEOquery::getGEOSuppFiles(GEO_ID, baseDir = paste(path,GEO_ID_path,"GEO_info",sep="/"))
	}, error=function(e){
		system(paste0("cd ",paste(path,GEO_ID_path,"GEO_info",sep="/")," && rm -rf *"))
		system(paste0("cd ",paste(path,GEO_ID_path,"GEO_info",sep="/")," && wget https://ftp.ncbi.nlm.nih.gov/geo/series/",paste0(substr(GEO_ID,1,nchar(GEO_ID)-3),"nnn"),"/",GEO_ID,"/matrix/",GEO_ID,"_series_matrix.txt.gz"))
		system(paste0("cd ",paste(path,GEO_ID_path,"GEO_info",sep="/")," && wget https://ftp.ncbi.nlm.nih.gov/geo/series/",paste0(substr(GEO_ID,1,nchar(GEO_ID)-3),"nnn"),"/",GEO_ID,"/soft/",GEO_ID,"_family.soft.gz"))
		system(paste0("mkdir -p ",paste(path,GEO_ID_path,"GEO_info",GEO_ID,sep="/")," && cd ",paste(path,GEO_ID_path,"GEO_info",GEO_ID,sep="/")," && wget ftp.ncbi.nlm.nih.gov/geo/series/",paste0(substr(GEO_ID,1,nchar(GEO_ID)-3),"nnn"),"/",GEO_ID,"/suppl/"))
		system(paste0('cd ',paste(path,GEO_ID_path,'GEO_info',GEO_ID,sep='/'),' && for f in $(grep href index.html | grep -v www. | grep -v // | grep -v nnn |  grep href index.html | grep -v www. | grep -v // | grep -v nnn | cut -d ">" -f 2 | cut -d "<" -f 1); do wget https://ftp.ncbi.nlm.nih.gov/geo/series/',paste0(substr(GEO_ID,1,nchar(GEO_ID)-3),'nnn'),'/',GEO_ID,'/suppl/$f; done'))
	})
}

# All the info... for example:
# a[[1]]@experimentData@other$title
# a[[1]]@phenoData@data$title
# a[[1]]@phenoData@data$characteristics_ch1.3
# a[[1]]@phenoData@data$relation.1
cat("\nDownloading from GEO done\n\n")

