#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
annotation <- args[1]      # GO annotation file (can be empty string "")
expression <- args[2]
path <- args[3]
pattern_search <- args[4]
kegg_annotation <- if (length(args) >= 5 && nchar(args[5]) > 0) args[5] else NULL
run_ora <- length(args) < 6 || args[6] != "no"   # "no" skips the over-representation analyses
minGS <- if (length(args) >= 7) as.numeric(args[7]) else 10
maxGS <- if (length(args) >= 8) as.numeric(args[8]) else 500

.rgse_scripts_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.rgse_scripts_dir, "R_qs_helpers.R"))
source(file.path(.rgse_scripts_dir, "R_net_config.R"))
.rgse_net_configure()

suppressMessages(library(clusterProfiler,quiet = T,warn.conflicts = F))
suppressMessages(library(aPEAR,quiet = T,warn.conflicts = F))
suppressMessages(library(ggplot2,quiet = T,warn.conflicts = F))

# GOfuncR is used *only* to turn GO IDs into human-readable names (go_name) and to
# classify them by ontology (root_node: biological_process / molecular_function /
# cellular_component); it is not needed for the enrichment itself. Depending on it
# hard-crashed the whole script for non-model organisms when the package was not
# installed. Resolve names via GOfuncR when present, otherwise via GO.db (which is
# organism-agnostic and already a clusterProfiler dependency), and only as a last
# resort fall back to the raw IDs so enrichment still runs. Returns a data.frame
# with go_id / go_name / root_node columns, one row per input id and in input order.
get_go_names <- function(ids) {
  ids <- as.character(ids)
  if (requireNamespace("GOfuncR", quietly = TRUE)) {
    return(GOfuncR::get_names(ids))
  }
  onto <- c(BP = "biological_process", MF = "molecular_function", CC = "cellular_component")
  if (requireNamespace("GO.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    message("NOTE: GOfuncR not installed; resolving GO term names and ontology via GO.db.")
    uid <- unique(ids[!is.na(ids) & grepl("^GO:", ids)])
    nm <- setNames(rep(NA_character_, length(uid)), uid)
    rt <- setNames(rep(NA_character_, length(uid)), uid)
    if (length(uid)) {
      res <- tryCatch(suppressWarnings(suppressMessages(
        AnnotationDbi::select(GO.db::GO.db, keys = uid,
                              columns = c("TERM", "ONTOLOGY"), keytype = "GOID"))),
        error = function(e) NULL)
      if (!is.null(res)) {
        res <- res[!duplicated(res$GOID), ]
        nm[res$GOID] <- res$TERM
        rt[res$GOID] <- onto[res$ONTOLOGY]
      }
    }
    return(data.frame(go_id = ids, go_name = unname(nm[ids]),
                      root_node = unname(rt[ids]), stringsAsFactors = FALSE))
  }
  message("WARNING: neither GOfuncR nor GO.db is installed; GO term names/ontology are unavailable, so GO enrichment will be empty. Install either package to enable it.")
  data.frame(go_id = ids, go_name = ids, root_node = NA_character_, stringsAsFactors = FALSE)
}

print(paste0("Current date: ",date()))
print(paste0("Performing clusterProfiler's functional enrichment analyses from data in the provided annotation file...",annotation))

setwd(path)
# sink(paste0(path,"/clusterProfiler_enrichr_funct_enrichment.log"));setwd(path)

expression_table <- as.data.frame(data.table::fread(paste0(expression),head=T))
expression_table$Gene_ID <- gsub(":.*","",expression_table$Gene_ID)

dge_files <- list.files(pattern = pattern_search)
for (f in if (run_ora) dge_files else character(0)){
  a <- as.data.frame(data.table::fread(paste0(f),head=T,fill=T))
  a$Gene_ID <- gsub(":.*","",a$Gene_ID)
  b <- a[a$FDR < 0.05,]
  d <- a[a$FDR < 0.05 & a$logFC>0,]
  e <- a[a$FDR < 0.05 & a$logFC<0,]
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
  b <- a[a$PValue < 0.05,]
  d <- a[a$PValue < 0.05 & a$logFC>0,]
  e <- a[a$PValue < 0.05 & a$logFC<0,]
  if (dim(b)[1]!=0){
    write.table(b,file=paste0(gsub(".txt","",f),"_pval_05.txt"),col.names = F,row.names = F,quote = F,sep="\t")  
    write.table(b$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")  
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

for (file in if (run_ora) list.files(path = path, pattern = "_Gene_IDs\\.txt$") else character(0)){
  print(paste0("Processing ",file," ..."))
  elements_interest <- data.table::fread(paste0(file),head=F)$V1
  elements_interest <- toupper(elements_interest)
  expression_table$Gene_ID <- toupper(expression_table$Gene_ID)

  # --- GO enrichment (only if GO annotation provided) ---
  has_go <- nchar(annotation) > 0 && file.exists(annotation)
  if (has_go) {
  a <- as.data.frame(data.table::fread(paste0(annotation),head=T,fill=T)[,1:2])
  a$source_id <- gsub(":.*","",a$source_id)
  colnames(a)[2] <- "Computed_GO_Process_IDs"
  b <- get_go_names(a$Computed_GO_Process_IDs)
  a$"Computed GO Processes" <- b$go_name
  a$Type <- b$root_node
  
  a$source_id <- toupper(a$source_id)
  
  # Match annotation source_id (e.g. ABLU38_00060) to full gene IDs in expression table (e.g. GCA_..._cds_ABLU38_00060_1) if no exact overlap
  expr_genes <- toupper(expression_table$Gene_ID)
  if (length(intersect(a$source_id, expr_genes)) == 0) {
    unique_sids <- unique(a$source_id)
    mapped <- vapply(unique_sids, function(sid) {
      hit <- expr_genes[grep(sid, expr_genes, fixed = TRUE)]
      if (length(hit) > 0) hit[1] else sid
    }, character(1))
    sid_to_full <- setNames(mapped, unique_sids)
    a$source_id <- unname(sid_to_full[a$source_id])
  }
  
  a_2 <- a[a$source_id %in% expr_genes,]


  # If user have GO annotation data (in data.frame format with first column of gene ID and second column of GO ID), they can use enricher() and gseGO() functions to perform over-representation test and gene set enrichment analysis.
  # If genes are annotated by direction annotation, it should also annotated by its ancestor GO nodes (indirect annation). If user only has direct annotation, they can pass their annotation to buildGOmap function, which will infer indirection annotation and generate a data.frame that suitable for both enricher() and gseGO()

  # https://rdrr.io/github/GuangchuangYu/clusterProfiler/man/enricher.html
  # https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
  # https://rdrr.io/bioc/clusterProfiler/man/buildGOmap.html


  ### GO BP
  try({
    b <- a[a$Type=="biological_process",]
    background_1 <- b[,c("Computed_GO_Process_IDs","source_id")]
    background_2 <- b[,c("Computed_GO_Process_IDs","Computed GO Processes")]
    # Move the terms that are separated by commas to new individual rows and use as instructed buildGOmap to get the ancestor IDs.
    s <- strsplit(as.character(background_1$Computed_GO_Process_IDs), split = ";")
    background_1 <- data.frame(BP_IDs = unlist(s), Isoforms = rep(background_1$source_id, sapply(s, length)))
    background_1 <- background_1[background_1$BP_IDs!="N/A",]
    background_2 <- data.frame(BP_IDs = unlist(strsplit(as.character(background_2$Computed_GO_Process_IDs),split=";")),BP_Names = unlist(strsplit(as.character(background_2$`Computed GO Processes`),split=";")))
    background_2 <- background_2[background_2$BP_IDs!="N/A",]
    background_1 <- background_1[complete.cases(background_1),]
    background_2 <- background_2[complete.cases(background_2),]
    bp <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1, TERM2NAME = background_2)
    background_1 <- clusterProfiler::buildGOmap(background_1)
    bp_2 <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1)
    write_1 <- as.data.frame(bp@result)
    write_2 <- as.data.frame(bp_2@result)
    write_1$Description <- get_go_names(write_1$ID)$go_name; write_2$Description <- get_go_names(write_2$ID)$go_name
    write.table(write_1,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_1.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(write_2,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_2.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_1, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_1_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_1_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_1)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_1_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_1_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_2, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_2_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_2_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_2)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_2_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_2_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  },silent=T)

  # With the universe of expressed genes instead of all annotated genes:
  try({
    b <- a_2[a_2$Type=="biological_process",]
    background_1 <- b[,c("Computed_GO_Process_IDs","source_id")]
    background_2 <- b[,c("Computed_GO_Process_IDs","Computed GO Processes")]
    # Move the terms that are separated by commas to new individual rows and use as instructed buildGOmap to get the ancestor IDs.
    s <- strsplit(as.character(background_1$Computed_GO_Process_IDs), split = ";")
    background_1 <- data.frame(BP_IDs = unlist(s), Isoforms = rep(background_1$source_id, sapply(s, length)))
    background_1 <- background_1[background_1$BP_IDs!="N/A",]
    background_2 <- data.frame(BP_IDs = unlist(strsplit(as.character(background_2$Computed_GO_Process_IDs),split=";")),BP_Names = unlist(strsplit(as.character(background_2$`Computed GO Processes`),split=";")))
    background_2 <- background_2[background_2$BP_IDs!="N/A",]
    bp <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1, TERM2NAME = background_2)
    background_1 <- clusterProfiler::buildGOmap(background_1)
    bp_2 <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05,TERM2GENE=background_1)
    write_1 <- as.data.frame(bp@result)
    write_2 <- as.data.frame(bp_2@result)
    write_1$Description <- get_go_names(write_1$ID)$go_name; write_2$Description <- get_go_names(write_2$ID)$go_name
    write.table(write_1,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_1.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(write_2,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_2.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_1, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_1_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_1_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_1)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_1_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_1_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_2, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_2_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_2_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_2)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_2_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_BP_expr_2_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  },silent=T)


  ### GO MF
  try({
    b <- a[a$Type=="molecular_function",]
    colnames(b)[2:3] <- c("Computed GO Function IDs","Computed GO Functions")
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
    write_1 <- as.data.frame(MF@result)
    write_2 <- as.data.frame(MF_2@result)
    write_1$Description <- get_go_names(write_1$ID)$go_name; write_2$Description <- get_go_names(write_2$ID)$go_name
    write.table(write_1,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_1.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(write_2,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_2.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_1, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_1_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_1_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_1)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_1_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_1_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_2, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_2_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_2_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_2)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_2_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_2_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  },silent=T)
  
  # With the universe of expressed genes instead of all annotated genes:
  try({
    b <- a_2[a_2$Type=="molecular_function",]
    colnames(b)[2:3] <- c("Computed GO Function IDs","Computed GO Functions")
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
    write_1 <- as.data.frame(MF@result)
    write_2 <- as.data.frame(MF_2@result)
    write_1$Description <- get_go_names(write_1$ID)$go_name; write_2$Description <- get_go_names(write_2$ID)$go_name
    write.table(write_1,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_1.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(write_2,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_2.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_1, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_1_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_1_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_1)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_1_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_1_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_2, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_2_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_2_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_2)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_2_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_MF_expr_2_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  },silent=T)


  ### GO CC
  try({
    b <- a[a$Type=="cellular_component",]
    colnames(b)[2:3] <- c("Computed GO Component IDs","Computed GO Components")
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
    write_1 <- as.data.frame(CC@result)
    write_2 <- as.data.frame(CC_2@result)
    write_1$Description <- get_go_names(write_1$ID)$go_name; write_2$Description <- get_go_names(write_2$ID)$go_name
    write.table(write_1,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_1.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(write_2,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_2.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_1, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_1_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_1_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_1)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_1_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_1_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_2, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_2_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_2_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_2)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_2_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_2_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  },silent=T)

  # With the universe of expressed genes instead of all annotated genes:
  try({
    b <- a_2[a_2$Type=="cellular_component",]
    colnames(b)[2:3] <- c("Computed GO Component IDs","Computed GO Components")
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
    write_1 <- as.data.frame(CC@result)
    write_2 <- as.data.frame(CC_2@result)
    write_1$Description <- get_go_names(write_1$ID)$go_name; write_2$Description <- get_go_names(write_2$ID)$go_name
    write.table(write_1,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_1.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(write_2,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_2.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_1, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_1_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_1_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_1)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_1_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_1_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    p <- enrichmentNetwork(write_2, repelLabels = TRUE, drawEllipses = TRUE)
    ggsave(p, filename = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_2_aPEAR.pdf"),width=15, height=15)
    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_2_aPEAR.html"),selfcontained = TRUE))
    clusters <- findPathClusters(write_2)
    write.table(clusters$clusters,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_2_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
    write.table(clusters$similarity,file=paste0(gsub(".txt","",file),"_functional_clusterProfiler_enrichr_CC_expr_2_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  },silent=T)

  } # end if (has_go)

  # --- KEGG Orthology enrichment (only if KEGG annotation provided) ---
  if (!is.null(kegg_annotation) && file.exists(kegg_annotation)) {
    print(paste0("  Running KEGG enrichment for ", file, " ..."))
    kegg_df <- as.data.frame(data.table::fread(kegg_annotation, head=T, fill=T)[,1:2])
    colnames(kegg_df) <- c("source_id", "KEGG_KO_ID")
    kegg_df$source_id <- toupper(gsub(":.*", "", kegg_df$source_id))
    kegg_df <- kegg_df[kegg_df$KEGG_KO_ID != "" & !is.na(kegg_df$KEGG_KO_ID),]
    kegg_df_expr <- kegg_df[kegg_df$source_id %in% expression_table$Gene_ID,]

    # Whole annotation background
    try({
      bg <- data.frame(KEGG_ID = kegg_df$KEGG_KO_ID, Gene = kegg_df$source_id)
      bg <- bg[complete.cases(bg),]
      kegg_res <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05, TERM2GENE = bg)
      write_k <- as.data.frame(kegg_res@result)
      write.table(write_k, file = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_1.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
      p <- enrichmentNetwork(write_k, repelLabels = TRUE, drawEllipses = TRUE)
      ggsave(p, filename = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_1_aPEAR.pdf"), width = 15, height = 15)
      suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p), file = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_1_aPEAR.html"), selfcontained = TRUE))
      clusters <- findPathClusters(write_k)
      write.table(clusters$clusters, file = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_1_aPEAR_clusters.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
      write.table(clusters$similarity, file = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_1_aPEAR_similarity.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    }, silent = T)

    # Expressed genes background
    try({
      bg_expr <- data.frame(KEGG_ID = kegg_df_expr$KEGG_KO_ID, Gene = kegg_df_expr$source_id)
      bg_expr <- bg_expr[complete.cases(bg_expr),]
      kegg_res <- clusterProfiler::enricher(elements_interest, pvalueCutoff = 0.05, TERM2GENE = bg_expr)
      write_k <- as.data.frame(kegg_res@result)
      write.table(write_k, file = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_expr_1.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
      p <- enrichmentNetwork(write_k, repelLabels = TRUE, drawEllipses = TRUE)
      ggsave(p, filename = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_expr_1_aPEAR.pdf"), width = 15, height = 15)
      suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p), file = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_expr_1_aPEAR.html"), selfcontained = TRUE))
      clusters <- findPathClusters(write_k)
      write.table(clusters$clusters, file = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_expr_1_aPEAR_clusters.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
      write.table(clusters$similarity, file = paste0(gsub(".txt", "", file), "_functional_clusterProfiler_enrichr_KEGG_expr_1_aPEAR_similarity.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    }, silent = T)
  } # end KEGG
}

# --- GSEA on all the tested genes, ranked by sign(logFC) * -log10(PValue) ---
gsea_sets <- list()
if (nchar(annotation) > 0 && file.exists(annotation)) {
  go <- as.data.frame(data.table::fread(paste0(annotation),head=T,fill=T)[,1:2])
  colnames(go) <- c("source_id","GO_ID")
  go$source_id <- toupper(gsub(":.*","",go$source_id))
  go <- go[!is.na(go$GO_ID) & grepl("^GO:",go$GO_ID),]
  go$Type <- get_go_names(go$GO_ID)$root_node
  for (ont in c("BP","MF","CC")) {
    t2g <- go[go$Type %in% c(BP="biological_process",MF="molecular_function",CC="cellular_component")[[ont]], c("GO_ID","source_id")]
    if (nrow(t2g) > 0) gsea_sets[[paste0("GO_",ont)]] <- clusterProfiler::buildGOmap(t2g)
  }
}
if (!is.null(kegg_annotation) && file.exists(kegg_annotation)) {
  kegg <- as.data.frame(data.table::fread(kegg_annotation, head=T, fill=T)[,1:2])
  colnames(kegg) <- c("source_id","KEGG_ID")
  kegg$source_id <- toupper(gsub(":.*","",kegg$source_id))
  kegg <- kegg[!is.na(kegg$KEGG_ID) & kegg$KEGG_ID != "", c("KEGG_ID","source_id")]
  if (nrow(kegg) > 0) gsea_sets[["KEGG"]] <- kegg
}
match_source_ids <- function(sids, genes) {
  if (length(intersect(sids, genes)) > 0) return(sids)
  u <- unique(sids)
  hit <- vapply(u, function(sid) { h <- genes[grep(sid, genes, fixed = TRUE)]; if (length(h) > 0) h[1] else sid }, character(1))
  unname(setNames(hit, u)[sids])
}
for (f in if (length(gsea_sets) > 0) dge_files else character(0)) {
  d <- as.data.frame(data.table::fread(paste0(f),head=T,fill=T))
  fc <- suppressWarnings(as.numeric(d[[grep("^logFC", colnames(d), value = TRUE)[1]]]))
  pv <- suppressWarnings(as.numeric(d$PValue))
  ids <- toupper(gsub(":.*","",d$Gene_ID))
  keep <- !is.na(fc) & !is.na(pv) & nzchar(ids) & !duplicated(ids)
  ranked <- sort(setNames(sign(fc[keep]) * -log10(pmax(pv[keep], .Machine$double.xmin)), ids[keep]), decreasing = TRUE)
  print(paste0("GSEA of ",f," on ",length(ranked)," ranked genes..."))
  for (set_name in names(gsea_sets)) {
    t2g <- gsea_sets[[set_name]]
    t2g[[2]] <- match_source_ids(t2g[[2]], names(ranked))
    t2n <- if (startsWith(set_name, "GO_")) data.frame(term = unique(t2g[[1]]), name = get_go_names(unique(t2g[[1]]))$go_name) else NA
    tryCatch({
      set.seed(123)
      res <- clusterProfiler::GSEA(ranked, TERM2GENE = t2g, TERM2NAME = t2n, minGSSize = minGS, maxGSSize = maxGS, pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, verbose = FALSE)
      if (is.null(res) || nrow(res@result) == 0) {
        print(paste0("  ",set_name,": no gene set with ",minGS,"-",maxGS," ranked genes"))
      } else {
        write.table(as.data.frame(res@result),file=paste0(gsub(".txt","",f),"_functional_clusterProfiler_GSEA_",set_name,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
        print(paste0("  ",set_name,": ",nrow(res@result)," gene sets, ",sum(res@result$p.adjust < 0.05)," with adjusted p < 0.05"))
      }
    }, error = function(e) print(paste0("  WARNING: GSEA ",set_name," failed for ",f,": ",conditionMessage(e))))
  }
}

dir.create(file.path(getwd(), "enrichr_clusterProfiler")); rgse_save_env(file.path(getwd(), "enrichr_clusterProfiler/globalenvir.qs2"))
# Move each file to the new directory
for (file in list.files(path = getwd(), pattern = "_functional_clusterProfiler_(enrichr|GSEA)_", full.names = TRUE)) {
  file.rename(from = file, to = file.path(file.path(getwd(), "enrichr_clusterProfiler"), basename(file)))
}

#setwd(path);sink()
print("ALL DONE clusterProfiler from data in the annotation...")
print(paste0("Current date: ",date()))
