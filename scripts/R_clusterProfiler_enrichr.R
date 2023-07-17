#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
annotation <- args[1]
expression <- args[2]
path <- args[3]
cores <- args[4]

suppressMessages(library(parallel,quiet = T,warn.conflicts = F))
suppressMessages(library(clusterProfiler,quiet = T,warn.conflicts = F))

print(paste0("Current date: ",date()))
print("Performing clusterProfiler's functional enrichment analyses from data in the annotation...")

setwd(path)
sink(paste0(path,"/clusterProfiler_enrichr_funct_enrichment.log"));setwd(path)

expression_table <- data.table::fread(expression,head=T)

for (f in list.files(pattern = "^DGE_analysis_comp.*.txt")){
  a <- data.table::fread(f,head=T,fill=T)
  b <- a[a$FDR < 0.05,1]
  d <- a[a$FDR < 0.05 & a$logFC>0,1]
  e <- a[a$FDR < 0.05 & a$logFC<0,1]
  if (dim(b)[1]!=0){
    write.table(b,file=paste0(gsub(".txt","",f),"_fdr_05.txt"),col.names = F,row.names = F,quote = F,sep="\t")  
    write.table(b$Gene_ID,file=paste0(gsub(".txt","",f),"_fdr_05_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
  }
  if (dim(d)[1]!=0){
    write.table(d,file=paste0(gsub(".txt","",f),"_fdr_05_logpos.txt"),col.names = F,row.names = F,quote = F,sep="\t")
    write.table(d$Gene_ID,file=paste0(gsub(".txt","",f),"_fdr_05_logpos_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
  }
  if (dim(e)[1]!=0){
    write.table(e,file=paste0(gsub(".txt","",f),"_fdr_05_logneg.txt"),col.names = F,row.names = F,quote = F,sep="\t")
    write.table(e$Gene_ID,file=paste0(gsub(".txt","",f),"_fdr_05_logneg_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
  }
  b <- a[a$PValue < 0.05,1]
  d <- a[a$PValue < 0.05 & a$logFC>0,1]
  e <- a[a$PValue < 0.05 & a$logFC<0,1]
  if (dim(b)[1]!=0){
    write.table(b,file=paste0(gsub(".txt","",f),"_pval_05.txt"),col.names = F,row.names = F,quote = F,sep="\t")  
    write.table(b$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05.txt_Gene_IDs"),col.names = F,row.names = F,quote = F,sep="\n")  
  }
  if (dim(d)[1]!=0){
    write.table(d,file=paste0(gsub(".txt","",f),"_pval_05_logpos.txt"),col.names = F,row.names = F,quote = F,sep="\t")
    write.table(d$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05_logpos_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
  }
  if (dim(e)[1]!=0){
    write.table(e,file=paste0(gsub(".txt","",f),"_pval_05_logneg.txt"),col.names = F,row.names = F,quote = F,sep="\t")
    write.table(e$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05_logneg_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
  }
}

process_file <- function(file){
  elements_interest <- data.table::fread(file,head=F)$V1
  a <- data.table::fread(annotation,fill=T)[,1:2]
  colnames(a)[2] <- "Computed GO Process IDs"
  b <- GOfuncR::get_names(a$"Computed GO Process IDs")
  a$"Computed GO Processes" <- b$go_name
  a$Type <- b$root_node
  a_2 <- a[a$source_id %in% expression_table$Gene_ID,]


  # If user have GO annotation data (in data.frame format with first column of gene ID and second column of GO ID), they can use enricher() and gseGO() functions to perform over-representation test and gene set enrichment analysis.
  # If genes are annotated by direction annotation, it should also annotated by its ancestor GO nodes (indirect annation). If user only has direct annotation, they can pass their annotation to buildGOmap function, which will infer indirection annotation and generate a data.frame that suitable for both enricher() and gseGO()

  # https://rdrr.io/github/GuangchuangYu/clusterProfiler/man/enricher.html
  # https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
  # https://rdrr.io/bioc/clusterProfiler/man/buildGOmap.html


  ### GO BP
  b <- a[a$Type=="biological_process",]
  background_1 <- b[,c("Computed GO Process IDs","source_id")]
  background_2 <- b[,c("Computed GO Process IDs","Computed GO Processes")]
  # Move the terms that are separated by commas to new individual rows and use as instructed buildGOmap to get the ancestor IDs.
  s <- strsplit(as.character(background_1$"Computed GO Process IDs"), split = ";")
  background_1 <- data.frame(BP_IDs = unlist(s), Isoforms = rep(background_1$source_id, sapply(s, length)))
  background_1 <- background_1[background_1$BP_IDs!="N/A",]
  background_2 <- data.frame(BP_IDs = unlist(strsplit(as.character(background_2$`Computed GO Process IDs`),split=";")),BP_Names = unlist(strsplit(as.character(background_2$`Computed GO Processes`),split=";")))
  background_2 <- background_2[background_2$BP_IDs!="N/A",]
  bp <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1, TERM2NAME = background_2)
  background_1 <- clusterProfiler::buildGOmap(background_1)
  bp_2 <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1)
  write.table(bp,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_BP_1.txt"),col.names = F,row.names = F,quote = F,sep="\t")
  write.table(bp_2,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_BP_2.txt"),col.names = F,row.names = F,quote = F,sep="\t")

  # With the universe of expressed genes instead of all annotated genes:
  b <- a_2[a_2$Type=="biological_process",]
  background_1 <- b[,c("Computed GO Process IDs","source_id")]
  background_2 <- b[,c("Computed GO Process IDs","Computed GO Processes")]
  # Move the terms that are separated by commas to new individual rows and use as instructed buildGOmap to get the ancestor IDs.
  s <- strsplit(as.character(background_1$"Computed GO Process IDs"), split = ";")
  background_1 <- data.frame(BP_IDs = unlist(s), Isoforms = rep(background_1$source_id, sapply(s, length)))
  background_1 <- background_1[background_1$BP_IDs!="N/A",]
  background_2 <- data.frame(BP_IDs = unlist(strsplit(as.character(background_2$`Computed GO Process IDs`),split=";")),BP_Names = unlist(strsplit(as.character(background_2$`Computed GO Processes`),split=";")))
  background_2 <- background_2[background_2$BP_IDs!="N/A",]
  bp <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1, TERM2NAME = background_2)
  background_1 <- clusterProfiler::buildGOmap(background_1)
  bp_2 <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1)
  write.table(bp,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_BP_expr_1.txt"),col.names = F,row.names = F,quote = F,sep="\t")
  write.table(bp_2,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_BP_expr_2.txt"),col.names = F,row.names = F,quote = F,sep="\t")


  ### GO MF
  b <- a[a$Type=="molecular_function",]
  background_1 <- b[,c("Computed GO Function IDs","source_id")]
  background_2 <- b[,c("Computed GO Function IDs","Computed GO Functions")]
  # Move the terms that are separated by commas to new individual rows and use as instructed buildGOmap to get the ancestor IDs.
  s <- strsplit(as.character(background_1$"Computed GO Function IDs"), split = ";")
  background_1 <- data.frame(MF_IDs = unlist(s), Isoforms = rep(background_1$source_id, sapply(s, length)))
  background_1 <- background_1[background_1$MF_IDs!="N/A",]
  background_2 <- data.frame(MF_IDs = unlist(strsplit(as.character(background_2$`Computed GO Function IDs`),split=";")),MF_Names = unlist(strsplit(as.character(background_2$`Computed GO Functions`),split=";")))
  background_2 <- background_2[background_2$MF_IDs!="N/A",]
  MF <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1, TERM2NAME = background_2)
  background_1 <- clusterProfiler::buildGOmap(background_1)
  MF_2 <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1)
  write.table(MF,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_MF_1.txt"),col.names = F,row.names = F,quote = F,sep="\t")
  write.table(MF_2,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_MF_2.txt"),col.names = F,row.names = F,quote = F,sep="\t")

  # With the universe of expressed genes instead of all annotated genes:
  b <- a_2[a_2$Type=="molecular_function",]
  background_1 <- b[,c("Computed GO Function IDs","source_id")]
  background_2 <- b[,c("Computed GO Function IDs","Computed GO Functions")]
  # Move the terms that are separated by commas to new individual rows and use as instructed buildGOmap to get the ancestor IDs.
  s <- strsplit(as.character(background_1$"Computed GO Function IDs"), split = ";")
  background_1 <- data.frame(MF_IDs = unlist(s), Isoforms = rep(background_1$source_id, sapply(s, length)))
  background_1 <- background_1[background_1$MF_IDs!="N/A",]
  background_2 <- data.frame(MF_IDs = unlist(strsplit(as.character(background_2$`Computed GO Function IDs`),split=";")),MF_Names = unlist(strsplit(as.character(background_2$`Computed GO Functions`),split=";")))
  background_2 <- background_2[background_2$MF_IDs!="N/A",]
  MF <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1, TERM2NAME = background_2)
  background_1 <- clusterProfiler::buildGOmap(background_1)
  MF_2 <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1)
  write.table(MF,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_MF_expr_1.txt"),col.names = F,row.names = F,quote = F,sep="\t")
  write.table(MF_2,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_MF_expr_2.txt"),col.names = F,row.names = F,quote = F,sep="\t")


  ### GO CC
  b <- a[a$Type=="cellular_component",]
  background_1 <- b[,c("Computed GO Component IDs","source_id")]
  background_2 <- b[,c("Computed GO Component IDs","Computed GO Components")]
  # Move the terms that are separated by commas to new individual rows and use as instructed buildGOmap to get the ancestor IDs.
  s <- strsplit(as.character(background_1$"Computed GO Component IDs"), split = ";")
  background_1 <- data.frame(CC_IDs = unlist(s), Isoforms = rep(background_1$source_id, sapply(s, length)))
  background_1 <- background_1[background_1$CC_IDs!="N/A",]
  background_2 <- data.frame(CC_IDs = unlist(strsplit(as.character(background_2$`Computed GO Component IDs`),split=";")),CC_Names = unlist(strsplit(as.character(background_2$`Computed GO Components`),split=";")))
  background_2 <- background_2[background_2$CC_IDs!="N/A",]
  CC <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1, TERM2NAME = background_2)
  background_1 <- clusterProfiler::buildGOmap(background_1)
  CC_2 <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1)
  write.table(CC,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_CC_1.txt"),col.names = F,row.names = F,quote = F,sep="\t")
  write.table(CC_2,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_CC_2.txt"),col.names = F,row.names = F,quote = F,sep="\t")

  # With the universe of expressed genes instead of all annotated genes:
  b <- a_2[a_2$Type=="cellular_component",]
  background_1 <- b[,c("Computed GO Component IDs","source_id")]
  background_2 <- b[,c("Computed GO Component IDs","Computed GO Components")]
  # Move the terms that are separated by commas to new individual rows and use as instructed buildGOmap to get the ancestor IDs.
  s <- strsplit(as.character(background_1$"Computed GO Component IDs"), split = ";")
  background_1 <- data.frame(CC_IDs = unlist(s), Isoforms = rep(background_1$source_id, sapply(s, length)))
  background_1 <- background_1[background_1$CC_IDs!="N/A",]
  background_2 <- data.frame(CC_IDs = unlist(strsplit(as.character(background_2$`Computed GO Component IDs`),split=";")),CC_Names = unlist(strsplit(as.character(background_2$`Computed GO Components`),split=";")))
  background_2 <- background_2[background_2$CC_IDs!="N/A",]
  CC <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1, TERM2NAME = background_2)
  background_1 <- clusterProfiler::buildGOmap(background_1)
  CC_2 <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1)
  write.table(CC,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_CC_expr_1.txt"),col.names = F,row.names = F,quote = F,sep="\t")
  write.table(CC_2,file=paste0(gsub(".txt","",file),"_clusterProfiler_enrichr_CC_expr_2.txt"),col.names = F,row.names = F,quote = F,sep="\t")
}


mclapply(
    mc.cores = cores,
    X = list.files(path = path, pattern = "_Gene_IDs\\.txt$"),
    FUN = process_file
)

setwd(path);sink()
print("ALL DONE clusterProfiler from data in the annotation...")
print(paste0("Current date: ",date()))
