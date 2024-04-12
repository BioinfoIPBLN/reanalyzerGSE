#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
organism <- args[2]
cores <- as.numeric(args[3])
padjustmethod <- args[4]
all_analyses <- args[5] # if not provided, "no"
cluster_enrich <- args[6] # if not provided, "no"
pattern_search <- args[7]


### Preparing:
print(paste0("Current date: ",date()))
print("Performing clusterProfiler's functional enrichment analyses...")
print("The majority of statistical analyses are going to be tested, and the most common plots performed. Even so, please do comprehensively check the whole vignette/manual, as for example there are dozens of possible plots. Links include:")
print("https://yulab-smu.top/biomedical-knowledge-mining-book/index.html"); print("https://bioc.ism.ac.jp/packages/3.7/bioc/vignettes/enrichplot/inst/doc/enrichplot.html") 
suppressMessages(library(clusterProfiler,quiet = T,warn.conflicts = F))
suppressMessages(library(AnnotationDbi,quiet = T,warn.conflicts = F))
suppressMessages(library(ReactomePA,quiet = T,warn.conflicts = F))
suppressMessages(library(DOSE,quiet = T,warn.conflicts = F))
suppressMessages(library(pathview,quiet = T,warn.conflicts = F))
suppressMessages(library(enrichplot,quiet = T,warn.conflicts = F))
suppressMessages(library(parallel,quiet = T,warn.conflicts = F))
suppressMessages(library(ggplot2,quiet = T,warn.conflicts = F))
suppressMessages(library(GOxploreR,quiet = T,warn.conflicts = F))
suppressMessages(library(aPEAR,quiet = T,warn.conflicts = F))

print(paste0("Trying to use as many as ",cores," cores, but if many subsets of genes/comparisons, please expect a lenghty process of at least a few hours"))
organism_cp <- gsub("_"," ",organism)
if(organism_cp=="Homo sapiens"){
  organism_cp_react <- "human"  
  orgDB <- "org.Hs.eg.db"
  suppressMessages(library("org.Hs.eg.db",quiet = T,warn.conflicts = F))
} else if(organism_cp=="Mus musculus"){
  organism_cp_react <- "mouse"
  orgDB <- "org.Mm.eg.db"
  suppressMessages(library("org.Mm.eg.db",quiet = T,warn.conflicts = F))
}
org <- search_kegg_organism(organism_cp, by='scientific_name')[1,1]
entrez_ids_keys <- select(eval(parse(text=orgDB)), keys=keys(eval(parse(text=orgDB)),keytype="ENTREZID"), columns=c("SYMBOL"), keytype="ENTREZID")
entrez_ids_keys$Custom_ID <- paste(entrez_ids_keys$SYMBOL,entrez_ids_keys$ENTREZID,sep="_")

### Create objects of interest to iterate later:
print(paste0("Creating subsets of genes of interest..."))
files <- list.files(path=path,pattern=pattern_search)
for (f in files){
  a <- read.table(paste0(path,"/",f),head=T)
  #readlist_cpm_fdr_05 <- a$logCP[a$FDR<0.05]
  #names(readlist_cpm_fdr_05) <- a$Gene_ID[a$FDR<0.05]
  #assign(paste0(f,"_","readlist_cpm_fdr_05"), readlist_cpm_fdr_05,envir = .GlobalEnv)
  #readlist_cpm_fdr_01 <- a$logCP[a$FDR<0.01]
  #names(readlist_cpm_fdr_01) <- a$Gene_ID[a$FDR<0.01]
  #assign(paste0(f,"_","readlist_cpm_fdr_01"), readlist_cpm_fdr_01,envir = .GlobalEnv)
  #readlist_cpm_pval_05 <- a$logCP[a$PValue<0.05]
  #names(readlist_cpm_pval_05) <- a$Gene_ID[a$PValue<0.05]
  #assign(paste0(f,"_","readlist_cpm_pval_05"), readlist_cpm_pval_05,envir = .GlobalEnv)
  #readlist_cpm_pval_01 <- a$logCP[a$PValue<0.01]
  #names(readlist_cpm_pval_01) <- a$Gene_ID[a$PValue<0.01]
  #assign(paste0(f,"_","readlist_cpm_pval_01"), readlist_cpm_pval_01,envir = .GlobalEnv)
  readlist_fc_fdr_05 <- a$logFC[a$FDR<0.05]
  names(readlist_fc_fdr_05) <- a$Gene_ID[a$FDR<0.05]
  readlist_fc_fdr_01 <- a$logFC[a$FDR<0.01]
  names(readlist_fc_fdr_01) <- a$Gene_ID[a$FDR<0.01]
  #readlist_fc_pval_05 <- a$logFC[a$PValue<0.05]
  #names(readlist_fc_pval_05) <- a$Gene_ID[a$PValue<0.05]
  #readlist_fc_pval_01 <- a$logFC[a$PValue<0.01]
  #names(readlist_fc_pval_01) <- a$Gene_ID[a$PValue<0.01]
  readlist_fc_backg <- a$logFC
  names(readlist_fc_backg) <- a$Gene_ID

  genes_of_interest <- list(fdr_05=a$Gene_ID[a$FDR<0.05],
                            fdr_01=a$Gene_ID[a$FDR<0.01],
                            #pval_05=a$Gene_ID[a$PValue<0.05],
                            #pval_01=a$Gene_ID[a$PValue<0.01],
                            fdr_05_fc_pos=a$Gene_ID[a$FDR<0.05 & a$logFC>0],
                            fdr_01_fc_pos=a$Gene_ID[a$FDR<0.01 & a$logFC>0],
                            #pval_05_fc_pos=a$Gene_ID[a$PValue<0.05 & a$logFC>0],
                            #pval_01_fc_pos=a$Gene_ID[a$PValue<0.01 & a$logFC>0],
                            fdr_05_fc_neg=a$Gene_ID[a$FDR<0.05 & a$logFC<0],
                            fdr_01_fc_neg=a$Gene_ID[a$FDR<0.01 & a$logFC<0],
                            #pval_05_fc_neg=a$Gene_ID[a$PValue<0.05 & a$logFC<0],
                            #pval_01_fc_neg=a$Gene_ID[a$PValue<0.01 & a$logFC<0],
                            readlist_fc_fdr_05=readlist_fc_fdr_05,
                            readlist_fc_fdr_01=readlist_fc_fdr_01,
                            #readlist_fc_pval_05=readlist_fc_pval_05,
                            #readlist_fc_pval_01=readlist_fc_pval_01,
                            readlist_fc_backg=readlist_fc_backg,
                            genes_backg=a$Gene_ID
                            )
  genes_of_interest <- genes_of_interest[unname(unlist(lapply(genes_of_interest,length)))!=0]
  assign(paste0(f,"_","genes_of_interest"), genes_of_interest,envir = .GlobalEnv) 
}
print(paste0("Done!"))
# grep("^readlist_|^entrez|^genes_of",ls(),val=T)
rm(genes_of_interest);rm(readlist_fc_fdr_01);rm(readlist_fc_fdr_05)

# Deduce the naming convention in the orgDB package:
check_naming <- function(names) {
  names <- grep("[[:punct:]]|orf|p43|p45",names,val=T,invert=T)  
  print(paste0("Looking at the annotation of ",length(names)," genes..."))
  all_upper <- grepl("^[A-Z0-9]+$", names)
  all_lower <- grepl("^[a-z0-9]+$", names)
  
  # Checks if the first letter is uppercase followed by lowercase letters or numbers, for names with at least one letter
  first_upper_rest_lower <- sapply(names, function(name) {
    if (grepl("[A-Za-z]", name)) { # Check if the name contains at least one letter
      # Extract the first letter and the rest of the string separately
      first_letter <- substr(name, regexpr("[A-Za-z]", name), regexpr("[A-Za-z]", name))
      rest <- substr(name, regexpr("[A-Za-z]", name) + 1, nchar(name))
      # Check if the first letter is uppercase and the rest of the string is lowercase or numeric (ignoring leading numbers)
      return(grepl("^[A-Z]$", first_letter) && grepl("^[a-z0-9]*$", rest))
    } else {
      return(TRUE) # If the name doesn't contain letters, it trivially satisfies the condition
    }
  })

  if(!all(all_upper) && !all(all_lower) && !all(first_upper_rest_lower)){
    pattern_results <- which.max(c(sum(all_upper),sum(all_lower),sum(first_upper_rest_lower)))
    num_genes <- c(sum(all_upper),sum(all_lower),sum(first_upper_rest_lower))[pattern_results]
    pattern_result_final <- c("all_upper","all_lower","first_upper_rest_lower")[pattern_results]
    print(paste0("Identified pattern for gene naming is ",pattern_result_final, ", accounting for ",num_genes," genes"))
  } else {
    pattern_result_final <- c("all_upper","all_lower","first_upper_rest_lower")[c(all_upper,all_lower,first_upper_rest_lower)]
    print(paste0("Identified pattern for gene naming is ",pattern_result_final, ", accounting for all annotated genes"))
  }
  
  return(pattern_result_final)

  
}
convert_ids <- function(ids,mode) {
  if(mode=="all_upper"){
    ids2 <- toupper(ids)
  } else if (mode=="all_lower") {
    ids2 <- tolower(ids)
  } else if (mode=="first_upper_rest_lower") {
    ids2 <- stringr::str_to_title(ids)
  }
  return(ids2)
}
mode <- check_naming(keys(eval(parse(text=orgDB)), keytype = "SYMBOL"))


process_file <- function(file){
  invisible(biomaRt::biomartCacheClear())
  # file="DGE_analysis_comp1.txt"
  print(paste0("Processing ",file, "_Current date: ",date()))

  basename <- gsub(".txt","",file)
  name <- paste0(file,"_genes_of_interest")
  
  a <- read.table(paste0(path,"/",file),head=T)
  a$logFC_sense <- a$logFC>0
  a$logFC_sense[a$logFC_sense] <- "POS"; a$logFC_sense[a$logFC_sense=="FALSE"] <- "NEG"
  path2 <- paste0(path,"/",basename,"_funct_enrich_clusterProfiler")
  genes_of_interest <- eval(parse(text=name))
  
  # FUNCTION INTERNAL:
  process_file_within <- function(name_internal){
    invisible(biomaRt::biomartCacheClear())
    print(paste0("Processing ",file,"_",name_internal))
    geneset <- name_internal
    i <- padjustmethod
    entrez_ids <- entrez_ids_keys$ENTREZID[entrez_ids_keys$SYMBOL %in% convert_ids(genes_of_interest[[geneset]],mode)]
    entrez_ids_backg <- entrez_ids_keys$ENTREZID[entrez_ids_keys$SYMBOL %in% convert_ids(genes_of_interest[["genes_backg"]],mode)]
    dir.create(paste0(path2,"/kegg_paths_snapshots"),recursive=T,showWarnings=F); dir.create(paste0(path2,"/reactome_paths_snapshots"),recursive=T,showWarnings=F);dir.create(paste0(path2,"/go_figs/"),recursive=T,showWarnings=F)
    setwd(path2)
    suppressMessages(library(orgDB, character.only = TRUE,quiet = T,warn.conflicts = F)) # Crucial apparently, so the functions using the orgDB object can be done in parallel

    if(all_analyses=="yes"){
      print(paste0("Processing ",file,"_",name_internal,"... Gene classification based on GO distribution at a specific level (2-6)"))
      ###### 1. Gene classification based on GO distribution at a specific level:
        for (levelgo in 2:6){
          # print(paste0("groupGO_level_",levelgo))
          tryCatch({
            ggo <- c(); Sys.sleep(2)
            ggo <- groupGO(gene     = convert_ids(genes_of_interest[[geneset]],mode),
                           OrgDb    = orgDB,
                           keyType  = "SYMBOL",
                           ont      = "BP",
                           level    = levelgo,
                           readable = TRUE)
            # ggo@result
            # ggo@result$Count
            write.table(ggo@result,file=paste0("GO_description_level_",levelgo,"_",geneset,"_groupGO_BP.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          }, error = function(e) {
            writeLines(as.character(e), paste0("GO_description_level_",levelgo,"_",geneset,"_groupGO_BP_err.txt"))
          })
          tryCatch({
            ggo <- c(); Sys.sleep(2)
            ggo <- groupGO(gene     = convert_ids(genes_of_interest[[geneset]],mode),
                           OrgDb    = orgDB,
                           keyType  = "SYMBOL",
                           ont      = "MF",
                           level    = levelgo,
                           readable = TRUE)
            # ggo@result
            # ggo@result$Count
            write.table(ggo@result,file=paste0("GO_description_level_",levelgo,"_",geneset,"_groupGO_MF.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          }, error = function(e) {
            writeLines(as.character(e), paste0("GO_description_level_",levelgo,"_",geneset,"_groupGO_MF_err.txt"))
          })
          tryCatch({
            ggo <- c(); Sys.sleep(2)
            ggo <- groupGO(gene     = convert_ids(genes_of_interest[[geneset]],mode),
                         OrgDb    = orgDB,
                         keyType  = "SYMBOL",
                         ont      = "CC",
                         level    = levelgo,
                         readable = TRUE)
            # ggo@result
            # ggo@result$Count
            write.table(ggo@result,file=paste0("GO_description_level_",levelgo,"_",geneset,"_groupGO_CC.txt"),col.names = T,row.names = F,quote = F,sep="\t")
        }, error = function(e) {
            writeLines(as.character(e), paste0("GO_description_level_",levelgo,"_",geneset,"_groupGO_CC_err.txt"))
          })
        }
        tryCatch({
          b <- c(); Sys.sleep(2)
          b <- Gene2GOTermAndLevel_ON(genes = entrez_ids_keys$ENTREZID[entrez_ids_keys$SYMBOL %in% convert_ids(genes_of_interest[[geneset]],mode)], organism = organism_cp, domain = "BP")
          b$GO_Description <- Term(b$"GO ID"); b$Gene_ID <- sapply(b$"Entrezgene ID",function(x){paste(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x],collapse=",")})
          # Assuming your data is in a data frame named 'df'
          d <- do.call(rbind, lapply(split(b, b$"GO ID"), function(group) {data.frame(GO_ID = unique(group$"GO ID"),GO_Description = unique(group$GO_Description),
                                                                                      Level = mean(group$Level), Gene_ID = paste(group$Gene_ID, collapse = ","),ENTREZ_ID = paste(group$"Entrezgene ID", collapse = ","))}))
          write.table(d,file=paste0("GO_description_all_",geneset,"_BP.txt"),col.names = T,row.names = F,quote = F,sep="\t")
        }, error = function(e) {
            writeLines(as.character(e), paste0("GO_description_all_",geneset,"_BP_err.txt"))
        })
        tryCatch({
          b <- c(); Sys.sleep(2)
          b <- Gene2GOTermAndLevel_ON(genes = entrez_ids_keys$ENTREZID[entrez_ids_keys$SYMBOL %in% convert_ids(genes_of_interest[[geneset]],mode)], organism = organism_cp, domain = "MF")
          b$GO_Description <- Term(b$"GO ID"); b$Gene_ID <- sapply(b$"Entrezgene ID",function(x){paste(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x],collapse=",")})
          # Assuming your data is in a data frame named 'df'
          d <- do.call(rbind, lapply(split(b, b$"GO ID"), function(group) {data.frame(GO_ID = unique(group$"GO ID"),GO_Description = unique(group$GO_Description),
                                                                                      Level = mean(group$Level), Gene_ID = paste(group$Gene_ID, collapse = ","),ENTREZ_ID = paste(group$"Entrezgene ID", collapse = ","))}))
          write.table(d,file=paste0("GO_description_all_",geneset,"_MF.txt"),col.names = T,row.names = F,quote = F,sep="\t")
        }, error = function(e) {
            writeLines(as.character(e), paste0("GO_description_all_",geneset,"_MF_err.txt"))
        })
        tryCatch({
          b <- c(); Sys.sleep(2)
          b <- Gene2GOTermAndLevel_ON(genes = entrez_ids_keys$ENTREZID[entrez_ids_keys$SYMBOL %in% convert_ids(genes_of_interest[[geneset]],mode)], organism = organism_cp, domain = "CC")
          b$GO_Description <- Term(b$"GO ID"); b$Gene_ID <- sapply(b$"Entrezgene ID",function(x){paste(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x],collapse=",")})
          # Assuming your data is in a data frame named 'df'
          d <- do.call(rbind, lapply(split(b, b$"GO ID"), function(group) {data.frame(GO_ID = unique(group$"GO ID"),GO_Description = unique(group$GO_Description),
                                                                                      Level = mean(group$Level), Gene_ID = paste(group$Gene_ID, collapse = ","),ENTREZ_ID = paste(group$"Entrezgene ID", collapse = ","))}))
          write.table(d,file=paste0("GO_description_all_",geneset,"_CC.txt"),col.names = T,row.names = F,quote = F,sep="\t")
        }, error = function(e) {
            writeLines(as.character(e), paste0("GO_description_all_",geneset,"_CC_err.txt"))
        })
      
    }
    ###### 2. GO over-representation analyses:
      print(paste0("Processing ",file,"_",name_internal,"... GO over-representation analyses"))
      tryCatch({
                  ego <- c(); Sys.sleep(2)
                  ego <- enrichGO(gene          = convert_ids(genes_of_interest[[geneset]],mode),
                                # universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))), # It's the same than default
                                universe = convert_ids(genes_of_interest[["genes_backg"]],mode),
                                OrgDb         = orgDB,
                                keyType  = "SYMBOL",
                                ont           = "ALL",
                                pAdjustMethod = i,
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
                  # head(ego)
                  ego_df <- as.data.frame(ego)
                  ego_df$summary_LogFC <- unlist(lapply(lapply(ego_df$geneID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                  ego_df$summary_LogFC_2 <- unlist(lapply(strsplit(ego_df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); ego_df$summary_LogFC <- unlist(lapply(strsplit(ego_df$geneID,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                  ego_df$summary_POS <- unlist(lapply(strsplit(ego_df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                  ego_df$summary_NEG <- unlist(lapply(strsplit(ego_df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                  write.table(ego_df,file=paste0("GO_overrepresentation_test_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  if(cluster_enrich=="yes"){
                    p <- enrichmentNetwork(ego@result, repelLabels = TRUE, drawEllipses = TRUE)
                    ggsave(p, filename = paste0("GO_overrepresentation_test_",i,"_",geneset,"_aPEAR.pdf"),width=30, height=30)
                    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_overrepresentation_test_",i,"_",geneset,"_aPEAR.html"),selfcontained = TRUE))
                    clusters <- findPathClusters(ego@result)
                    write.table(clusters$clusters,file=paste0("GO_overrepresentation_test_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    write.table(clusters$similarity,file=paste0("GO_overrepresentation_test_",i,"_",geneset,"_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  }
      }, error = function(e) {
          writeLines(as.character(e), paste0("GO_overrepresentation_test_",i,"_",geneset,"_err.txt"))
      })
      tryCatch({
                  ego <- c(); Sys.sleep(2)
                  ego <- enrichGO(gene          = convert_ids(genes_of_interest[[geneset]],mode),
                                # universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                                universe = convert_ids(genes_of_interest[["genes_backg"]],mode),
                                OrgDb         = orgDB,
                                keyType  = "SYMBOL",
                                ont           = "BP",
                                pAdjustMethod = i,
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
                  ego_df <- as.data.frame(ego)
                  ego_df$summary_LogFC <- unlist(lapply(lapply(ego_df$geneID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                  ego_df$summary_LogFC_2 <- unlist(lapply(strsplit(ego_df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); ego_df$summary_LogFC <- unlist(lapply(strsplit(ego_df$geneID,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                  ego_df$summary_POS <- unlist(lapply(strsplit(ego_df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                  ego_df$summary_NEG <- unlist(lapply(strsplit(ego_df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                  write.table(ego_df,file=paste0("GO_overrepresentation_test_BP_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  suppressMessages(ggsave(goplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,".pdf"),width=30, height=30))
                  suppressMessages(ggsave(barplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_barplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(dotplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_dotplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(cnetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_cnetplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(heatplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_heatplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(treeplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_treeplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(emapplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_emapplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(upsetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_upsetplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(pmcplot(ego$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_pmcplot.pdf"),width=30, height=30))
                  if(cluster_enrich=="yes"){          
                    p <- enrichmentNetwork(ego@result, repelLabels = TRUE, drawEllipses = TRUE)
                    ggsave(p, filename = paste0("GO_overrepresentation_test_BP_",i,"_",geneset,"_aPEAR.pdf"),width=30, height=30)
                    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_overrepresentation_test_BP_",i,"_",geneset,"_aPEAR.html"),selfcontained = TRUE))
                    clusters <- findPathClusters(ego@result)
                    write.table(clusters$clusters,file=paste0("GO_overrepresentation_test_BP_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    write.table(clusters$similarity,file=paste0("GO_overrepresentation_test_BP_",i,"_",geneset,"_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")                  
                  }
      }, error = function(e) {
          writeLines(as.character(e), paste0("GO_overrepresentation_test_BP_",i,"_",geneset,"_err.txt"))
      })
      tryCatch({
                  ego <- c(); Sys.sleep(2)
                  ego <- enrichGO(gene          = convert_ids(genes_of_interest[[geneset]],mode),
                                # universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                                universe = convert_ids(genes_of_interest[["genes_backg"]],mode),
                                OrgDb         = orgDB,
                                keyType  = "SYMBOL",
                                ont           = "MF",
                                pAdjustMethod = i,
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
                  ego_df <- as.data.frame(ego)
                  ego_df <- as.data.frame(ego)
                  ego_df$summary_LogFC <- unlist(lapply(lapply(ego_df$geneID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                  ego_df$summary_LogFC_2 <- unlist(lapply(strsplit(ego_df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); ego_df$summary_LogFC <- unlist(lapply(strsplit(ego_df$geneID,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                  ego_df$summary_POS <- unlist(lapply(strsplit(ego_df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                  ego_df$summary_NEG <- unlist(lapply(strsplit(ego_df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                  write.table(ego_df,file=paste0("GO_overrepresentation_test_MF_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  suppressMessages(ggsave(goplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,".pdf"),width=30, height=30))
                  suppressMessages(ggsave(barplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_barplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(dotplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_dotplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(cnetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_cnetplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(heatplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_heatplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(treeplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_treeplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(emapplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_emapplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(upsetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_upsetplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(pmcplot(ego$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_pmcplot.pdf"),width=30, height=30))
                  if(cluster_enrich=="yes"){
                    p <- enrichmentNetwork(ego@result, repelLabels = TRUE, drawEllipses = TRUE)
                    ggsave(p, filename = paste0("GO_overrepresentation_test_MF_",i,"_",geneset,"_aPEAR.pdf"),width=30, height=30)
                    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_overrepresentation_test_MF_",i,"_",geneset,"_aPEAR.html"),selfcontained = TRUE))
                    clusters <- findPathClusters(ego@result)
                    write.table(clusters$clusters,file=paste0("GO_overrepresentation_test_MF_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    write.table(clusters$similarity,file=paste0("GO_overrepresentation_test_MF_",i,"_",geneset,"_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")                  
                  }
      }, error = function(e) {
          writeLines(as.character(e), paste0("GO_overrepresentation_test_MF_",i,"_",geneset,"_err.txt"))
      })
      tryCatch({
                  ego <- c(); Sys.sleep(2)
                  ego <- enrichGO(gene          = convert_ids(genes_of_interest[[geneset]],mode),
                              #universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                              universe = convert_ids(genes_of_interest[["genes_backg"]],mode),
                              OrgDb         = orgDB,
                              keyType  = "SYMBOL",
                              ont           = "CC",
                              pAdjustMethod = i,
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 1,
                              readable      = TRUE)
                  ego_df <- as.data.frame(ego)
                  ego_df$summary_LogFC <- unlist(lapply(lapply(ego_df$geneID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                  ego_df$summary_LogFC_2 <- unlist(lapply(strsplit(ego_df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); ego_df$summary_LogFC <- unlist(lapply(strsplit(ego_df$geneID,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                  ego_df$summary_POS <- unlist(lapply(strsplit(ego_df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                  ego_df$summary_NEG <- unlist(lapply(strsplit(ego_df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                  write.table(ego_df,file=paste0("GO_overrepresentation_test_CC_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  suppressMessages(ggsave(goplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,".pdf"),width=30, height=30))
                  suppressMessages(ggsave(barplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_barplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(dotplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_dotplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(cnetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_cnetplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(heatplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_heatplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(treeplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_treeplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(emapplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_emapplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(upsetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_upsetplot.pdf"),width=30, height=30))
                  suppressMessages(ggsave(pmcplot(ego$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_pmcplot.pdf"),width=30, height=30))
                  if(cluster_enrich=="yes"){
                    p <- enrichmentNetwork(ego@result, repelLabels = TRUE, drawEllipses = TRUE)
                    ggsave(p, filename = paste0("GO_overrepresentation_test_CC_",i,"_",geneset,"_aPEAR.pdf"),width=30, height=30)
                    suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_overrepresentation_test_CC_",i,"_",geneset,"_aPEAR.html"),selfcontained = TRUE))
                    clusters <- findPathClusters(ego@result)
                    write.table(clusters$clusters,file=paste0("GO_overrepresentation_test_CC_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    write.table(clusters$similarity,file=paste0("GO_overrepresentation_test_CC_",i,"_",geneset,"_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")                  
                  }
      }, error = function(e) {
          writeLines(as.character(e), paste0("GO_overrepresentation_test_CC_",i,"_",geneset,"_err.txt"))
      })
    

    if(all_analyses=="yes"){
      print(paste0("Processing ",file,"_",name_internal,"... Gene Set Enrichment Analysis of Gene Ontology"))
        ###### 3. Gene Set Enrichment Analysis of Gene Ontology:        
          if(geneset=="fdr_05" || geneset=="fdr_01"){
              f <- paste0("readlist_fc_",geneset)
              b <- sort(unlist(genes_of_interest[f]),decreasing=T); names(b) <- convert_ids(gsub(paste0(f,"."),"",names(b)),mode)
              tryCatch({
                        gse_enrich <- c(); Sys.sleep(2)
                        gse_enrich <- suppressMessages(gseGO(geneList= b,
                                           OrgDb        = orgDB,
                                           ont          = "ALL",
                                           keyType  = "SYMBOL",
                                           pvalueCutoff = 1,
                                           pAdjustMethod = i,
                                           verbose      = FALSE,
                                           by="fgsea"))
                        df <- gse_enrich@result
                        df$summary_LogFC <- unlist(lapply(lapply(df$core_enrichment, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                        df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                        df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                        df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                        write.table(df,file=paste0("GO_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        if(cluster_enrich=="yes"){
                          p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                          ggsave(p, filename = paste0("GO_GSEA_",f,"_",i,"_fgsea_aPEAR.pdf"),width=30, height=30)
                          suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_GSEA_",f,"_",i,"_fgsea_aPEAR.html"),selfcontained = TRUE))
                          clusters <- findPathClusters(df)
                          write.table(clusters$clusters,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                          write.table(clusters$similarity,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        }
              }, error = function(e) {
                        writeLines(as.character(e), paste0("GO_GSEA_",f,"_",i,"_fgsea_err.txt"))
              })
              tryCatch({
                        gse_enrich <- c(); Sys.sleep(2)
                        gse_enrich <- suppressMessages(gseGO(geneList=b,
                                         OrgDb        = orgDB,
                                         ont          = "BP",
                                         keyType  = "SYMBOL",
                                         pvalueCutoff = 1,
                                         pAdjustMethod = i,
                                         verbose      = FALSE,
                                         by="fgsea"))
                        df <- gse_enrich@result
                        df$summary_LogFC <- unlist(lapply(lapply(df$core_enrichment, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                        df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                        df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                        df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                        write.table(df,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_BP.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        suppressMessages(ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP.pdf"),width=30, height=30))
                        suppressMessages(ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_dotplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_cnetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_heatplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_treeplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_emapplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_upsetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_pmcplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_ridgeplot.pdf"),width=30, height=30))
                        if(cluster_enrich=="yes"){
                          p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                          ggsave(p, filename = paste0("GO_GSEA_",f,"_",i,"_fgsea_BP_aPEAR.pdf"),width=30, height=30)
                          suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_GSEA_",f,"_",i,"_fgsea_BP_aPEAR.html"),selfcontained = TRUE))
                          clusters <- findPathClusters(df)
                          write.table(clusters$clusters,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_BP_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                          write.table(clusters$similarity,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_BP_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        }
              }, error = function(e) {
                        writeLines(as.character(e), paste0("GO_GSEA_",f,"_",i,"_fgsea_BP_err.txt"))
              })
              tryCatch({
                        gse_enrich <- c(); Sys.sleep(2)
                        gse_enrich <- suppressMessages(gseGO(geneList=b,
                                         OrgDb        = orgDB,
                                         ont          = "MF",
                                         keyType  = "SYMBOL",
                                         pvalueCutoff = 1,
                                         pAdjustMethod = i,
                                         verbose      = FALSE,
                                         by="fgsea"))
                        df <- gse_enrich@result
                        df$summary_LogFC <- unlist(lapply(lapply(df$core_enrichment, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                        df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                        df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                        df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                        write.table(df,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_MF.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        suppressMessages(ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF.pdf"),width=30, height=30))
                        suppressMessages(ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_dotplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_cnetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_heatplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_treeplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_emapplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_upsetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_pmcplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_ridgeplot.pdf"),width=30, height=30))
                        if(cluster_enrich=="yes"){
                          p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                          ggsave(p, filename = paste0("GO_GSEA_",f,"_",i,"_fgsea_MF_aPEAR.pdf"),width=30, height=30)
                          suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_GSEA_",f,"_",i,"_fgsea_MF_aPEAR.html"),selfcontained = TRUE))
                          clusters <- findPathClusters(df)
                          write.table(clusters$clusters,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_MF_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                          write.table(clusters$similarity,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_MF_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        }
              }, error = function(e) {
                        writeLines(as.character(e), paste0("GO_GSEA_",f,"_",i,"_fgsea_MF_err.txt"))
              })
              tryCatch({
                        gse_enrich <- c(); Sys.sleep(2)
                        gse_enrich <- suppressMessages(gseGO(geneList=b,
                                         OrgDb        = orgDB,
                                         ont          = "CC",
                                         keyType  = "SYMBOL",
                                         pvalueCutoff = 1,
                                         pAdjustMethod = i,
                                         verbose      = FALSE,
                                         by="fgsea"))
                        df <- gse_enrich@result
                        df$summary_LogFC <- unlist(lapply(lapply(df$core_enrichment, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                        df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                        df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                        df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                        write.table(df,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_CC.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        suppressMessages(ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC.pdf"),width=30, height=30))
                        suppressMessages(ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_dotplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_cnetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_heatplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_treeplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_emapplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_upsetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_pmcplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_ridgeplot.pdf"),width=30, height=30))
                        if(cluster_enrich=="yes"){
                          p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                          ggsave(p, filename = paste0("GO_GSEA_",f,"_",i,"_fgsea_CC_aPEAR.pdf"),width=30, height=30)
                          suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_GSEA_",f,"_",i,"_fgsea_CC_aPEAR.html"),selfcontained = TRUE))
                          clusters <- findPathClusters(df)
                          write.table(clusters$clusters,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_CC_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                          write.table(clusters$similarity,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_CC_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        }
              }, error = function(e) {
                        writeLines(as.character(e), paste0("GO_GSEA_",f,"_",i,"_fgsea_CC_err.txt"))
              })
              tryCatch({
                        gse_enrich <- c(); Sys.sleep(2)
                        gse_enrich <- suppressMessages(gseGO(geneList=b,
                                         OrgDb        = orgDB,
                                         ont          = "ALL",
                                         keyType  = "SYMBOL",
                                         pvalueCutoff = 1,
                                         pAdjustMethod = i,
                                         verbose      = FALSE,
                                         by="DOSE",
                                         nPerm = 100))
                        df <- gse_enrich@result
                        df$summary_LogFC <- unlist(lapply(lapply(df$core_enrichment, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                        df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                        df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                        df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                        write.table(df,file=paste0("GO_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        if(cluster_enrich=="yes"){
                          p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                          ggsave(p, filename = paste0("GO_GSEA_",f,"_",i,"_DOSE_aPEAR.pdf"),width=30, height=30)
                          suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_GSEA_",f,"_",i,"_DOSE_aPEAR.html"),selfcontained = TRUE))
                          clusters <- findPathClusters(df)
                          write.table(clusters$clusters,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                          write.table(clusters$similarity,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        }
              }, error = function(e) {
                        writeLines(as.character(e), paste0("GO_GSEA_",f,"_",i,"_DOSE_err.txt"))
              })
              tryCatch({
                        gse_enrich <- c(); Sys.sleep(2)
                        gse_enrich <- suppressMessages(gseGO(geneList=b,
                                        OrgDb        = orgDB,
                                        ont          = "BP",
                                        keyType  = "SYMBOL",
                                        pvalueCutoff = 1,
                                        pAdjustMethod = padjustmethod,
                                        verbose      = FALSE,
                                        by="DOSE",
                                        nPerm = 100))
                        df <- gse_enrich@result
                        df$summary_LogFC <- unlist(lapply(lapply(df$core_enrichment, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                        df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                        df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                        df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                        write.table(df,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_BP.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        suppressMessages(ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP.pdf"),width=30, height=30))
                        suppressMessages(ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_dotplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_cnetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_heatplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_treeplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_emapplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_upsetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_pmcplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_ridgeplot.pdf"),width=30, height=30))
                        if(cluster_enrich=="yes"){
                          p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                          ggsave(p, filename = paste0("GO_GSEA_",f,"_",i,"_DOSE_BP_aPEAR.pdf"),width=30, height=30)
                          suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_GSEA_",f,"_",i,"_DOSE_BP_aPEAR.html"),selfcontained = TRUE))
                          clusters <- findPathClusters(df)
                          write.table(clusters$clusters,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_BP_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                          write.table(clusters$similarity,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_BP_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        }
              }, error = function(e) {
                        writeLines(as.character(e), paste0("GO_GSEA_",f,"_",i,"_DOSE_BP_err.txt"))
              })
              tryCatch({
                        gse_enrich <- c(); Sys.sleep(2)
                        gse_enrich <- suppressMessages(gseGO(geneList=b,
                                        OrgDb        = orgDB,
                                        ont          = "MF",
                                        keyType  = "SYMBOL",
                                        pvalueCutoff = 1,
                                        pAdjustMethod = i,
                                        verbose      = FALSE,
                                        by="DOSE",
                                        nPerm = 100))
                        df <- gse_enrich@result
                        df$summary_LogFC <- unlist(lapply(lapply(df$core_enrichment, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                        df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                        df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                        df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                        write.table(df,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_MF.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        suppressMessages(ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF.pdf"),width=30, height=30))
                        suppressMessages(ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_dotplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_cnetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_heatplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_treeplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_emapplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_upsetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_pmcplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_ridgeplot.pdf"),width=30, height=30))
                        if(cluster_enrich=="yes"){
                          p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                          ggsave(p, filename = paste0("GO_GSEA_",f,"_",i,"_DOSE_MF_aPEAR.pdf"),width=30, height=30)
                          suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_GSEA_",f,"_",i,"_DOSE_MF_aPEAR.html"),selfcontained = TRUE))
                          clusters <- findPathClusters(df)
                          write.table(clusters$clusters,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_MF_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                          write.table(clusters$similarity,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_MF_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        }
              }, error = function(e) {
                        writeLines(as.character(e), paste0("GO_GSEA_",f,"_",i,"_DOSE_MF_err.txt"))
              })
              tryCatch({
                        gse_enrich <- c(); Sys.sleep(2)
                        gse_enrich <- suppressMessages(gseGO(geneList=b,
                                        OrgDb        = orgDB,
                                        ont          = "CC",
                                        keyType  = "SYMBOL",
                                        pvalueCutoff = 1,
                                        pAdjustMethod = i,
                                        verbose      = FALSE,
                                        by="DOSE",
                                        nPerm = 100))
                        df <- gse_enrich@result
                        df$summary_LogFC <- unlist(lapply(lapply(df$core_enrichment, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,"/"))]}),function(y){paste(y,collapse=",")}))
                        df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                        df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                        df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                        write.table(df,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_CC.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        suppressMessages(ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC.pdf"),width=30, height=30))
                        suppressMessages(ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_dotplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_cnetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_heatplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_treeplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_emapplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_upsetplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_pmcplot.pdf"),width=30, height=30))
                        suppressMessages(ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_ridgeplot.pdf"),width=30, height=30))
                        if(cluster_enrich=="yes"){
                          p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                          ggsave(p, filename = paste0("GO_GSEA_",f,"_",i,"_DOSE_CC_aPEAR.pdf"),width=30, height=30)
                          suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("GO_GSEA_",f,"_",i,"_DOSE_CC_aPEAR.html"),selfcontained = TRUE))
                          clusters <- findPathClusters(df)
                          write.table(clusters$clusters,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_CC_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                          write.table(clusters$similarity,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_CC_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        }
              }, error = function(e) {
                        writeLines(as.character(e), paste0("GO_GSEA_",f,"_",i,"_DOSE_CC_err.txt"))
              })
          }  
        
        
        ###### 4. KEGG over-representation:    
          print(paste0("Processing ",file,"_",name_internal,"... KEGG over-representation"))
          tryCatch({
              kk <- suppressMessages(enrichKEGG(gene= entrez_ids,
                               organism     = org,
                               universe = entrez_ids_backg,
                               pAdjustMethod = i,
                               pvalueCutoff = 1,
                               qvalueCutoff = 1))
              suppressMessages(ggsave(barplot(kk, showCategory=20), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_barplot.pdf"),width=30, height=30))
              suppressMessages(ggsave(dotplot(kk, showCategory=20), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_dotplot.pdf"),width=30, height=30))
              suppressMessages(ggsave(cnetplot(kk), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_cnetplot.pdf"),width=30, height=30))
              suppressMessages(ggsave(heatplot(kk, showCategory=20), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_heatplot.pdf"),width=30, height=30))
              suppressMessages(ggsave(emapplot(pairwise_termsim(kk)), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_emapplot.pdf"),width=30, height=30))
              suppressMessages(ggsave(upsetplot(kk), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_upsetplot.pdf"),width=30, height=30))
              kk_write <- kk@result
              kk_write$Gene_ID <-  unlist(lapply(strsplit(kk_write$geneID,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
              kk_write$Custom_id <- unlist(lapply(strsplit(kk_write$geneID,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
              kk_write$summary_LogFC <- unlist(lapply(lapply(kk_write$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
              kk_write$summary_LogFC_2 <- unlist(lapply(strsplit(kk_write$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); kk_write$summary_LogFC <- unlist(lapply(strsplit(kk_write$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
              kk_write$summary_up <- unlist(lapply(strsplit(kk_write$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
              kk_write$summary_down <- unlist(lapply(strsplit(kk_write$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
              write.table(kk_write,file=paste0("KEGG_enrich_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
              if(cluster_enrich=="yes"){
                p <- enrichmentNetwork(kk_write, repelLabels = TRUE, drawEllipses = TRUE)
                ggsave(p, filename = paste0("KEGG_enrich_",i,"_",geneset,"_aPEAR.pdf"),width=30, height=30)
                suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("KEGG_enrich_",i,"_",geneset,"_aPEAR.html"),selfcontained = TRUE))
                clusters <- findPathClusters(kk_write)
                write.table(clusters$clusters,file=paste0("KEGG_enrich_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                write.table(clusters$similarity,file=paste0("KEGG_enrich_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              }
          }, error = function(e) {
              writeLines(as.character(e), paste0("KEGG_enrich_",i,"_",geneset,"_err.txt"))
          })
        
        
        ###### 5. KEGG pathways visualization:
              print(paste0("Processing ",file,"_",name_internal,"... KEGG pathways visualization"))
              paths_list <- unique(unlist(lapply(list.files(path = getwd(), pattern = "^KEGG", full.names = TRUE), function(x){
                data <- read.delim(x)
                data_filtered <- data$ID[data$pvalue < 0.05]
              })))
              for (f in paths_list){
                tryCatch({
                  suppressMessages(pathview(gene.data=kk@result,pathway.id=f,species=org,kegg.dir=paste0(getwd(),"/kegg_paths_snapshots")))
                }, error = function(e) {
                  writeLines(as.character(e), paste0(getwd(),"/kegg_paths_snapshots/err.txt"))
                })
              }
              invisible(file.remove(list.files(path=getwd(),pattern="*.pathview.png",full.names = T)))
                
        ###### 6. Gene Set Enrichment Analysis of KEGG:
          if(geneset=="fdr_05" || geneset=="fdr_01"){
              print(paste0("Processing ",file,"_",name_internal,"... Gene Set Enrichment Analysis of KEGG"))
              f <- paste0("readlist_fc_",geneset)
              b <- sort(unlist(genes_of_interest[f]),decreasing=T); names(b) <- convert_ids(gsub(paste0(f,"."),"",names(b)),mode)
              b <- b[names(b) %in% entrez_ids_keys$SYMBOL]; names(b) <- entrez_ids_keys$ENTREZID[match(names(b),entrez_ids_keys$SYMBOL)]
              tryCatch({
                      gse_enrich <- c(); Sys.sleep(2)
                      gse_enrich <- suppressMessages(gseKEGG(geneList= b,
                                                             organism = org,
                                                             pvalueCutoff = 1,
                                                             pAdjustMethod = i,
                                                             verbose      = FALSE,
                                                             by="fgsea"))
                      df <- gse_enrich@result
                      df$Gene_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                      df$Custom_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                      df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                      df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                      df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                      df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                      write.table(df,file=paste0("KEGG_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      for (k in gse_enrich@result$ID[gse_enrich@result$pvalue < 0.05]){
                        suppressMessages(pathview(gene.data=b,
                                                  pathway.id=k,species=org))
                        invisible(file.remove(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),invert=T,val=T)))
                        invisible(file.rename(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),val=T),paste0(getwd(),"/kegg_paths_snapshots/",k,"_","KEGG_GSEA_",f,"_",i,"_fgsea.png")))
                      }
                      if(cluster_enrich=="yes"){
                        p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                        ggsave(p, filename = paste0("KEGG_GSEA_",f,"_",i,"_fgsea_aPEAR.pdf"),width=30, height=30)
                        suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("KEGG_GSEA_",f,"_",i,"_fgsea_aPEAR.html"),selfcontained = TRUE))
                        clusters <- findPathClusters(df)
                        write.table(clusters$clusters,file=paste0("KEGG_GSEA_",f,"_",i,"_fgsea_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        write.table(clusters$similarity,file=paste0("KEGG_GSEA_",f,"_",i,"_fgsea_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      }
              }, error = function(e) {
                  writeLines(as.character(e), paste0("KEGG_GSEA_",f,"_",i,"_fgsea_err.txt"))
              })
              tryCatch({
                      gse_enrich <- c(); Sys.sleep(2)
                      gse_enrich <- suppressMessages(gseKEGG(geneList= b,
                                                             organism = org,
                                                             pvalueCutoff = 1,
                                                             pAdjustMethod = i,
                                                             verbose      = FALSE,
                                                             by="DOSE"))
                      df <- gse_enrich@result
                      df$Gene_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                      df$Custom_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                      df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                      df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                      df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                      df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                      write.table(df,file=paste0("KEGG_DOSE_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      for (k in gse_enrich@result$ID[gse_enrich@result$pvalue < 0.05]){
                        suppressMessages(pathview(gene.data=b,
                                                  pathway.id=k,species=org))
                        invisible(file.remove(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),invert=T,val=T)))
                        invisible(file.rename(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),val=T),paste0(getwd(),"/kegg_paths_snapshots/",k,"_","KEGG_DOSE_",f,"_",i,"_DOSE.png")))
                      }
                      if(cluster_enrich=="yes"){
                        p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                        ggsave(p, filename = paste0("KEGG_GSEA_",f,"_",i,"_DOSE_aPEAR.pdf"),width=30, height=30)
                        suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("KEGG_GSEA_",f,"_",i,"_DOSE_aPEAR.html"),selfcontained = TRUE))
                        clusters <- findPathClusters(df)
                        write.table(clusters$clusters,file=paste0("KEGG_GSEA_",f,"_",i,"_DOSE_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        write.table(clusters$similarity,file=paste0("KEGG_GSEA_",f,"_",i,"_DOSE_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      }
              }, error = function(e) {
                  writeLines(as.character(e), paste0("KEGG_DOSE_",f,"_",i,"_DOSE_err.txt"))
              })
          }
        
        
        ###### 7. KEGG Module over-representation:
          print(paste0("Processing ",file,"_",name_internal,"... KEGG Module over-representation"))
          tryCatch({
                gse_enrich <- c(); Sys.sleep(2)
                gse_enrich <- suppressMessages(enrichMKEGG(gene= entrez_ids,
                                                           universe = entrez_ids_backg,
                                                           organism = org,
                                                           pvalueCutoff = 1,
                                                           qvalueCutoff = 1,
                                                           pAdjustMethod = i))
                kk_write <- gse_enrich@result
                kk_write$Gene_ID <-  unlist(lapply(strsplit(kk_write$geneID,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                kk_write$Custom_id <- unlist(lapply(strsplit(kk_write$geneID,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                kk_write$summary_LogFC <- unlist(lapply(lapply(kk_write$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                kk_write$summary_LogFC_2 <- unlist(lapply(strsplit(kk_write$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); kk_write$summary_LogFC <- unlist(lapply(strsplit(kk_write$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                kk_write$summary_up <- unlist(lapply(strsplit(kk_write$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                kk_write$summary_down <- unlist(lapply(strsplit(kk_write$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                write.table(kk_write,file=paste0("MKEGG_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
                if(cluster_enrich=="yes"){
                  p <- enrichmentNetwork(kk_write, repelLabels = TRUE, drawEllipses = TRUE)
                  ggsave(p, filename = paste0("MKEGG_",i,"_",geneset,"_aPEAR.pdf"),width=30, height=30)
                  suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("MKEGG_",i,"_",geneset,"_aPEAR.html"),selfcontained = TRUE))
                  clusters <- findPathClusters(kk_write)
                  write.table(clusters$clusters,file=paste0("MKEGG_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  write.table(clusters$similarity,file=paste0("MKEGG_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                }
          }, error = function(e) {
                  writeLines(as.character(e), paste0("MKEGG_",i,"_",geneset,"_err.txt"))
          })
        
        
        ###### 8. Gene Set Enrichment Analysis of KEGG modules:      
          if(geneset=="fdr_05" || geneset=="fdr_01"){
                print(paste0("Processing ",file,"_",name_internal,"... Gene Set Enrichment Analysis of KEGG modules"))
                f <- paste0("readlist_fc_",geneset)
                b <- sort(unlist(genes_of_interest[f]),decreasing=T); names(b) <- convert_ids(gsub(paste0(f,"."),"",names(b)),mode)
                b <- b[names(b) %in% entrez_ids_keys$SYMBOL]; names(b) <- entrez_ids_keys$ENTREZID[match(names(b),entrez_ids_keys$SYMBOL)]
                tryCatch({
                      gse_enrich <- c(); Sys.sleep(2)
                      gse_enrich <- suppressMessages(gseMKEGG(geneList= b,
                                                              organism = org,
                                                              pvalueCutoff = 1,
                                                              pAdjustMethod = i,
                                                              verbose      = FALSE,
                                                              by="fgsea"))
                      df <- gse_enrich@result
                      df$Gene_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                      df$Custom_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                      df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                      df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                      df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                      df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                      write.table(df,file=paste0("MKEGG_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      if(cluster_enrich=="yes"){
                        p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                        ggsave(p, filename = paste0("MKEGG_GSEA_",f,"_",i,"_fgsea_aPEAR.pdf"),width=30, height=30)
                        suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("MKEGG_GSEA_",f,"_",i,"_fgsea_aPEAR.html"),selfcontained = TRUE))
                        clusters <- findPathClusters(df)
                        write.table(clusters$clusters,file=paste0("MKEGG_GSEA_",f,"_",i,"_fgsea_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        write.table(clusters$similarity,file=paste0("MKEGG_GSEA_",f,"_",i,"_fgsea_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      }
                }, error = function(e) {
                  writeLines(as.character(e), paste0("MKEGG_GSEA_",f,"_",i,"_fgsea_err.txt"))
                })
                tryCatch({
                      gse_enrich <- c(); Sys.sleep(2)
                      gse_enrich <- suppressMessages(gseMKEGG(geneList= b,
                                                              organism = org,
                                                              pvalueCutoff = 1,
                                                              pAdjustMethod = i,
                                                              verbose      = FALSE,
                                                              by="DOSE",
                                                              nPerm = 100))
                      df <- gse_enrich@result
                      df$Gene_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                      df$Custom_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                      df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                      df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                      df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                      df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                      write.table(df,file=paste0("MKEGG_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      if(cluster_enrich=="yes"){
                        p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                        ggsave(p, filename = paste0("MKEGG_GSEA_",f,"_",i,"_DOSE_aPEAR.pdf"),width=30, height=30)
                        suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("MKEGG_GSEA_",f,"_",i,"_DOSE_aPEAR.html"),selfcontained = TRUE))
                        clusters <- findPathClusters(df)
                        write.table(clusters$clusters,file=paste0("MKEGG_GSEA_",f,"_",i,"_DOSE_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                        write.table(clusters$similarity,file=paste0("MKEGG_GSEA_",f,"_",i,"_DOSE_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      }
                }, error = function(e) {
                  writeLines(as.character(e), paste0("MKEGG_GSEA_",f,"_",i,"_DOSE_err.txt"))
                })
          }    
        
        
        ###### 9. WikiPathways over-representation:
          print(paste0("Processing ",file,"_",name_internal,"... WikiPathways over-representation"))
          tryCatch({
            gse_enrich <- c(); Sys.sleep(2)
            gse_enrich <- suppressMessages(enrichWP(gene= entrez_ids,
                                                       organism = organism_cp,
                                                       universe = entrez_ids_backg,
                                                       pvalueCutoff = 1,
                                                       qvalueCutoff = 1,
                                                       pAdjustMethod = i))
            df <- gse_enrich@result
            df$Gene_ID <- unlist(lapply(strsplit(df$geneID,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
            df$Custom_ID <- unlist(lapply(strsplit(df$geneID,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
            df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
            df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
            df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
            df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
            write.table(df,file=paste0("WP_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
            if(cluster_enrich=="yes"){
              p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
              ggsave(p, filename = paste0("WP_",i,"_",geneset,"_aPEAR.pdf"),width=30, height=30)
              suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("WP_",i,"_",geneset,"_aPEAR.html"),selfcontained = TRUE))
              clusters <- findPathClusters(df)
              write.table(clusters$clusters,file=paste0("WP_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              write.table(clusters$similarity,file=paste0("WP_",i,"_",geneset,"_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            }
          }, error = function(e) {
                  writeLines(as.character(e), paste0("WP_",i,"_",geneset,"_err.txt"))
          })            
        
        
        ###### 10. Gene set enrichment analyses of WikiPathways: 
          if(geneset=="fdr_05" || geneset=="fdr_01"){
            print(paste0("Processing ",file,"_",name_internal,"... Gene set enrichment analyses of WikiPathways"))
            f <- paste0("readlist_fc_",geneset)
            b <- sort(unlist(genes_of_interest[f]),decreasing=T); names(b) <- convert_ids(gsub(paste0(f,"."),"",names(b)),mode)
            b <- b[names(b) %in% entrez_ids_keys$SYMBOL]; names(b) <- entrez_ids_keys$ENTREZID[match(names(b),entrez_ids_keys$SYMBOL)]
            tryCatch({
                    gse_enrich <- c(); Sys.sleep(2)
                    gse_enrich <- suppressMessages(gseWP(geneList= b,
                                                            organism = organism_cp,
                                                            pvalueCutoff = 1,
                                                            pAdjustMethod = i,
                                                            verbose      = FALSE,
                                                            by="fgsea"))
                    df <- gse_enrich@result
                    df$Gene_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                    df$Custom_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                    df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                    df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                    df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                    df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                    write.table(df,file=paste0("WP_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    if(cluster_enrich=="yes"){
                      p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                      ggsave(p, filename = paste0("WP_GSEA_",i,"_",geneset,"_fgsea_aPEAR.pdf"),width=30, height=30)
                      suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("WP_GSEA",i,"_",geneset,"_fgsea_aPEAR.html"),selfcontained = TRUE))
                      clusters <- findPathClusters(df)
                      write.table(clusters$clusters,file=paste0("WP_GSEA_",i,"_",geneset,"_fgsea_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      write.table(clusters$similarity,file=paste0("WP_GSEA_",i,"_",geneset,"_fgsea_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }
            }, error = function(e) {
                  writeLines(as.character(e), paste0("WP_GSEA_",f,"_",i,"_fgsea_err.txt"))
            })
            tryCatch({
                    gse_enrich <- c(); Sys.sleep(2)
                    gse_enrich <- suppressMessages(gseWP(geneList= b,
                                                            organism = organism_cp,
                                                            pvalueCutoff = 1,
                                                            pAdjustMethod = i,
                                                            verbose      = FALSE,
                                                            by="DOSE",
                                                            nPerm = 100))
                    df <- gse_enrich@result
                    df$Gene_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                    df$Custom_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                    df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                    df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                    df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                    df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                    write.table(df,file=paste0("WP_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    if(cluster_enrich=="yes"){
                      p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                      ggsave(p, filename = paste0("WP_GSEA_",i,"_",geneset,"_DOSE_aPEAR.pdf"),width=30, height=30)
                      suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("WP_GSEA",i,"_",geneset,"_DOSE_aPEAR.html"),selfcontained = TRUE))
                      clusters <- findPathClusters(df)
                      write.table(clusters$clusters,file=paste0("WP_GSEA_",i,"_",geneset,"_DOSE_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      write.table(clusters$similarity,file=paste0("WP_GSEA_",i,"_",geneset,"_DOSE_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }
            }, error = function(e) {
                  writeLines(as.character(e), paste0("WP_GSEA_",f,"_",i,"_DOSE_err.txt"))
            })
          }
        
        
        ###### 11. Reactome over-representation:      
          print(paste0("Processing ",file,"_",name_internal,"... Reactome over-representation"))
          tryCatch({
            gse_enrich <- c(); Sys.sleep(2)
            gse_enrich <- suppressMessages(enrichPathway(gene= entrez_ids,
                                                    organism = organism_cp_react,
                                                    universe = entrez_ids_backg,
                                                    pvalueCutoff = 1,
                                                    qvalueCutoff = 1,
                                                    pAdjustMethod = i))
            df <- gse_enrich@result
            df$Gene_ID <- unlist(lapply(strsplit(df$geneID,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
            df$Custom_ID <- unlist(lapply(strsplit(df$geneID,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
            df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
            df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
            df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
            df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
            write.table(df,file=paste0("REACT_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
            if(cluster_enrich=="yes"){
              p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
              ggsave(p, filename = paste0("REACT_",i,"_",geneset,"_aPEAR.pdf"),width=30, height=30)
              suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("REACT_",i,"_",geneset,"_aPEAR.html"),selfcontained = TRUE))
              clusters <- findPathClusters(df)
              write.table(clusters$clusters,file=paste0("REACT_",i,"_",geneset,"_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              write.table(clusters$similarity,file=paste0("REACT_",i,"_",geneset,"_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            }
          }, error = function(e) {
                  writeLines(as.character(e), paste0("REACT_",i,"_",geneset,"_err.txt"))
          })
        
        
        ###### 12. Gene set enrichment analyses of Reactome:
          if(geneset=="fdr_05" || geneset=="fdr_01"){
                  print(paste0("Processing ",file,"_",name_internal,"... Gene set enrichment analyses of Reactome"))
                  f <- paste0("readlist_fc_",geneset)
                  b <- sort(unlist(genes_of_interest[f]),decreasing=T); names(b) <- convert_ids(gsub(paste0(f,"."),"",names(b)),mode)
                  b <- b[names(b) %in% entrez_ids_keys$SYMBOL]; names(b) <- entrez_ids_keys$ENTREZID[match(names(b),entrez_ids_keys$SYMBOL)]
                  tryCatch({
                    gse_enrich <- c(); Sys.sleep(2)
                    gse_enrich <- suppressMessages(gsePathway(geneList=b,
                                                         organism = organism_cp_react,
                                                         pvalueCutoff = 1,
                                                         pAdjustMethod = i,
                                                         verbose      = FALSE,
                                                         by="fgsea"))
                    df <- gse_enrich@result
                    df$Gene_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                    df$Custom_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                    df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                    df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                    df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                    df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                    write.table(df,file=paste0("REACT_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    if(cluster_enrich=="yes"){
                      p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                      ggsave(p, filename = paste0("REACT_GSEA_",i,"_",geneset,"_fgsea_aPEAR.pdf"),width=30, height=30)
                      suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("REACT_GSEA",i,"_",geneset,"_fgsea_aPEAR.html"),selfcontained = TRUE))
                      clusters <- findPathClusters(df)
                      write.table(clusters$clusters,file=paste0("REACT_GSEA_",i,"_",geneset,"_fgsea_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      write.table(clusters$similarity,file=paste0("REACT_GSEA_",i,"_",geneset,"_fgsea_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }
                  }, error = function(e) {
                    writeLines(as.character(e), paste0("REACT_GSEA_",f,"_",i,"_fgsea_err.txt"))
                  })
                  tryCatch({
                    gse_enrich <- c(); Sys.sleep(2)
                    gse_enrich <- suppressMessages(gsePathway(geneList=b,
                                                         organism = organism_cp_react,
                                                         pvalueCutoff = 1,
                                                         pAdjustMethod = i,
                                                         verbose      = FALSE,
                                                         by="DOSE",
                                                         nPerm = 100))
                    df <- gse_enrich@result
                    df$Gene_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$SYMBOL[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                    df$Custom_ID <- unlist(lapply(strsplit(df$core_enrichment,"/"),function(x){paste(unique(entrez_ids_keys$Custom_ID[entrez_ids_keys$ENTREZID %in% x]),collapse=",")}))
                    df$summary_LogFC <- unlist(lapply(lapply(df$Gene_ID, function(x){a$logFC_sense[a$Gene_ID %in% unlist(strsplit(x,","))]}),function(y){paste(y,collapse=",")}))
                    df$summary_LogFC_2 <- unlist(lapply(strsplit(df$summary_LogFC,",",),function(x){paste(names(table(x)),table(x),sep="_",collapse=",")})); df$summary_LogFC <- unlist(lapply(strsplit(df$Gene_ID,","),function(x){paste(x,a$logFC_sense[a$Gene_ID %in% x],sep="_",collapse="/")}))
                    df$summary_POS <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_POS","",grep("_POS",x,val=T)))),collapse=",")}))
                    df$summary_NEG <- unlist(lapply(strsplit(df$summary_LogFC,"/"),function(x){paste(sort(unique(gsub("_NEG","",grep("_NEG",x,val=T)))),collapse=",")}))
                    write.table(df,file=paste0("REACT_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    if(cluster_enrich=="yes"){
                      p <- enrichmentNetwork(df, repelLabels = TRUE, drawEllipses = TRUE)
                      ggsave(p, filename = paste0("REACT_GSEA_",i,"_",geneset,"_DOSE_aPEAR.pdf"),width=30, height=30)
                      suppressWarnings(htmlwidgets::saveWidget(widget = plotly::ggplotly(p),file = paste0("REACT_GSEA",i,"_",geneset,"_DOSE_aPEAR.html"),selfcontained = TRUE))
                      clusters <- findPathClusters(df)
                      write.table(clusters$clusters,file=paste0("REACT_GSEA_",i,"_",geneset,"_DOSE_aPEAR_clusters.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                      write.table(clusters$similarity,file=paste0("REACT_GSEA_",i,"_",geneset,"_DOSE_aPEAR_similarity.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }
                  }, error = function(e) {
                    writeLines(as.character(e), paste0("REACT_GSEA_",f,"_",i,"_DOSE_err.txt"))
                  })
          }        
        
        ###### 13. Reactome pathways visualization:
          print(paste0("Processing ",file,"_",name_internal,"... Reactome pathways visualization"))
          paths_list <- unique(unlist(lapply(list.files(path = getwd(), pattern = "^REACT", full.names = TRUE), function(x){
              data <- read.delim(x)
              data_filtered <- data$Description[data$pvalue < 0.05]
          })))
          for (f in paths_list){
            tryCatch({
              suppressMessages(ggsave(viewPathway(f,readable = TRUE,organism = organism_cp_react),filename = paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome.pdf"),width=30, height=30))
            }, error = function(e) {
              writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
            })
          }        
            
        ###### 14. Over-representation analyses for human databases (DO, NCG and DGN):
          if(organism_cp=="Homo sapiens"){
                  print(paste0("Processing ",file,"_",name_internal,"... Over-representation analyses for human databases (DO, NCG and DGN)"))
                  tryCatch({
                    ego <- c(); Sys.sleep(2)
                    ego <- enrichDO(gene          = entrez_ids,
                                    ont           = "DO",
                                    universe = entrez_ids_backg,
                                    pAdjustMethod = i,
                                    pvalueCutoff  = 1,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)
                    # head(ego)
                    write.table(ego,file=paste0("DO_overrepresentation_test_",i,"_",geneset,"_DO.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  }, error = function(e) {
                    writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                  })
                  tryCatch({
                    ego <- c(); Sys.sleep(2)
                    ego <- enrichDO(gene          = entrez_ids,
                                    ont           = "DOLite",
                                    universe = entrez_ids_backg,
                                    pAdjustMethod = i,
                                    pvalueCutoff  = 1,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)
                    # head(ego)
                    write.table(ego,file=paste0("DO_overrepresentation_test_",i,"_",geneset,"_DOLite.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  }, error = function(e) {
                    writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                  })
                  tryCatch({
                      ego <- c(); Sys.sleep(2)  
                      ego <- enrichNCG(gene          = entrez_ids,
                                       pAdjustMethod = i,
                                       universe = entrez_ids_backg,
                                       pvalueCutoff  = 1,
                                       qvalueCutoff  = 1,
                                       readable      = TRUE)
                      # head(ego)
                      write.table(ego,file=paste0("NGC_overrepresentation_test_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  }, error = function(e) {
                    writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                  })
                  tryCatch({
                      ego <- c(); Sys.sleep(2)  
                      ego <- enrichDGN(gene          = entrez_ids,
                                       pAdjustMethod = i,
                                       universe = entrez_ids_backg,
                                       pvalueCutoff  = 1,
                                       qvalueCutoff  = 1,
                                       readable      = TRUE)
                      # head(ego)
                      write.table(ego,file=paste0("DGN_overrepresentation_test_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
                  }, error = function(e) {
                    writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                  })
                       
                
                  if(geneset=="fdr_05" || geneset=="fdr_01"){
                    f <- paste0("readlist_fc_",geneset)
                    b <- sort(unlist(genes_of_interest[f]),decreasing=T); names(b) <- convert_ids(gsub(paste0(f,"."),"",names(b)),mode)
                    b <- b[names(b) %in% entrez_ids_keys$SYMBOL]; names(b) <- entrez_ids_keys$ENTREZID[match(names(b),entrez_ids_keys$SYMBOL)]
                    tryCatch({
                          gse_enrich <- c(); Sys.sleep(2)
                          gse_enrich <- suppressMessages(gseDO(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                                    pvalueCutoff = 1,
                                                                    pAdjustMethod = i,
                                                                    verbose      = FALSE,
                                                                    by="fgsea"))
                          write.table(gse_enrich@result,file=paste0("DO_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }, error = function(e) {
                      writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                    })
                    tryCatch({
                          gse_enrich <- c(); Sys.sleep(2)  
                          gse_enrich <- suppressMessages(gseDO(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                                    pvalueCutoff = 1,
                                                                    pAdjustMethod = i,
                                                                    verbose      = FALSE,
                                                                    by="DOSE",
                                                                    nPerm = 100))
                          write.table(gse_enrich@result,file=paste0("DO_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }, error = function(e) {
                      writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                    })
                    
                    tryCatch({
                          gse_enrich <- c(); Sys.sleep(2)  
                          gse_enrich <- suppressMessages(gseNCG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                               pvalueCutoff = 1,
                                                               pAdjustMethod = i,
                                                               verbose      = FALSE,
                                                               by="fgsea"))
                          write.table(gse_enrich@result,file=paste0("NCG_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }, error = function(e) {
                      writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                    })
                    tryCatch({
                          gse_enrich <- c(); Sys.sleep(2)  
                          gse_enrich <- suppressMessages(gseNCG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                               pvalueCutoff = 1,
                                                               pAdjustMethod = i,
                                                               verbose      = FALSE,
                                                               by="DOSE",
                                                               nPerm = 100))
                          write.table(gse_enrich@result,file=paste0("NCG_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }, error = function(e) {
                      writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                    })
                    
                    tryCatch({
                          gse_enrich <- c(); Sys.sleep(2)  
                          gse_enrich <- suppressMessages(gseDGN(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                                pvalueCutoff = 1,
                                                                pAdjustMethod = i,
                                                                verbose      = FALSE,
                                                                by="fgsea"))
                          write.table(gse_enrich@result,file=paste0("DGN_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }, error = function(e) {
                      writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                    })
                    tryCatch({
                          gse_enrich <- c(); Sys.sleep(2)  
                          gse_enrich <- suppressMessages(gseDGN(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                                pvalueCutoff = 1,
                                                                pAdjustMethod = i,
                                                                verbose      = FALSE,
                                                                by="DOSE",
                                                                nPerm = 100))
                          write.table(gse_enrich@result,file=paste0("DGN_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
                    }, error = function(e) {
                      writeLines(as.character(e), paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome_err.pdf"))
                    })
                  }
          }
        }
    }
    if(cores==1){cores_2=1} else {cores_2=6}
    mclapply(
      mc.cores = cores_2,
      X = grep("_backg",grep("readlist_",names(genes_of_interest),invert=T,val=T),invert=T,val=T),
      FUN = process_file_within
    )     

    print(paste0("Processing ",file,"... Tidying up..."))
    ### Remove empty files and tar everything:
    system(paste0("find ", path2, " -type f -exec awk 'NR>1{exit 1}' {} \\; -exec rm -f {} \\;"))      
    print(paste0("DONE with ",file, "_Current date: ",date()))
  }

### Parallel processing:
print("List of files:");cat(files, sep = "\n")
print(paste0("Processing... Many figures are going to be showed, but please double check the parameters used (default) and the possibilities within https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html. it is likely that some of the plots should be redo manually to optimize them..."))
print(paste0("Processing... Remember that to visualize particular pathways of interest in KEGG, the functions pathway() and browseKEGG can be used, https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html#visualize-enriched-kegg-pathways"))
print(paste0("Processing... Here all the pathways that may be of interest in KEGG (i..e pval < 0.05) are shown for all alternatives"))
print(paste0("Processing... Remember that to visualize particular pathways of interest in Reactome, the function viewPathway() can be used, https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html#pathway-visualization"))
print(paste0("Processing... Here all the pathways that may be of interest in Reactome (i..e pval < 0.05) are shown, but without coloring because it would change for each case... please redo if required adding the parameter 'foldChange=' in the function"))

if(cores==1){cores_3=1} else {cores_3=round(as.numeric(cores)/6)}
mclapply(
    mc.cores = cores_3,
    X = files,
    FUN = process_file
)

# Tidying...
setwd(path)
print("Tidying...")
removeEmptyDirs <- function(directory) {
  # List all directories
  dirs <- list.dirs(directory, recursive = TRUE)
  dirs <- dirs[grep("funct",dirs)]
  
  # Check each directory
  for (fold in dirs) {
    # If the directory is empty
    if (length(dir(fold)) < 3 && fold!=path) {
      # Remove the directory
      invisible(unlink(fold, recursive = TRUE))
      print(paste0("Removing ",fold))
    }
  }
}
removeEmptyDirs(path)

save.image(file.path(path,"funct_enrich_clusterProfiler_globalenvir.RData"))
print("ALL DONE clusterProfiler")
print(paste0("Current date: ",date()))
