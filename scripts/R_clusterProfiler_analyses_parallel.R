#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
organism <- args[2]
cores <- args[3]

### Preparing:
print("Performing clusterProfilers functional enrichment analyses...")
print("The majority of statistical analyses are going to be tested, and the most common plots performed. Even so, please do comprehensively check the whole vignette/manual, as for example there are dozens of possible plots. Links include:")
print("https://yulab-smu.top/biomedical-knowledge-mining-book/index.html"); print("https://bioc.ism.ac.jp/packages/3.7/bioc/vignettes/enrichplot/inst/doc/enrichplot.html") 
suppressMessages(library(clusterProfiler,quiet = T,warn.conflicts = F))
suppressMessages(library(AnnotationDbi,quiet = T,warn.conflicts = F))
suppressMessages(library(ReactomePA,quiet = T,warn.conflicts = F))
suppressMessages(library(DOSE,quiet = T,warn.conflicts = F))
suppressMessages(library(pathview,quiet = T,warn.conflicts = F))
suppressMessages(library(enrichplot,quiet = T,warn.conflicts = F))
suppressMessages(library(parallel,quiet = T,warn.conflicts = F))



process_file <- function(file){
  a <- read.table(paste0(path,"/",file),head=T)
  readlist_cpm_fdr_05 <- a$logCP[a$FDR<0.05]
  names(readlist_cpm_fdr_05) <- a$Gene_ID[a$FDR<0.05]
  readlist_cpm_fdr_01 <- a$logCP[a$FDR<0.01]
  names(readlist_cpm_fdr_01) <- a$Gene_ID[a$FDR<0.01]
  readlist_cpm_pval_05 <- a$logCP[a$PValue<0.05]
  names(readlist_cpm_pval_05) <- a$Gene_ID[a$PValue<0.05]
  readlist_cpm_pval_01 <- a$logCP[a$PValue<0.01]
  names(readlist_cpm_pval_01) <- a$Gene_ID[a$PValue<0.01]
  readlist_fc_fdr_05 <- a$logFC[a$FDR<0.05]
  names(readlist_fc_fdr_05) <- a$Gene_ID[a$FDR<0.05]
  readlist_fc_fdr_01 <- a$logFC[a$FDR<0.01]
  names(readlist_fc_fdr_01) <- a$Gene_ID[a$FDR<0.01]
  readlist_fc_pval_05 <- a$logFC[a$PValue<0.05]
  names(readlist_fc_pval_05) <- a$Gene_ID[a$PValue<0.05]
  readlist_fc_pval_01 <- a$logFC[a$PValue<0.01]
  names(readlist_fc_pval_01) <- a$Gene_ID[a$PValue<0.01]

  # organism <- organism
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
  genes_of_interest <- list(fdr_05=a$Gene_ID[a$FDR<0.05],
                            fdr_01=a$Gene_ID[a$FDR<0.01],
                            pval_05=a$Gene_ID[a$PValue<0.05],
                            pval_01=a$Gene_ID[a$PValue<0.01],
                            fdr_05_fc_pos=a$Gene_ID[a$FDR<0.05 & a$logFC>0],
                            fdr_01_fc_pos=a$Gene_ID[a$FDR<0.01 & a$logFC>0],
                            pval_05_fc_pos=a$Gene_ID[a$PValue<0.05 & a$logFC>0],
                            pval_01_fc_pos=a$Gene_ID[a$PValue<0.01 & a$logFC>0],
                            fdr_05_fc_neg=a$Gene_ID[a$FDR<0.05 & a$logFC<0],
                            fdr_01_fc_neg=a$Gene_ID[a$FDR<0.01 & a$logFC<0],
                            pval_05_fc_neg=a$Gene_ID[a$PValue<0.05 & a$logFC<0],
                            pval_01_fc_neg=a$Gene_ID[a$PValue<0.01 & a$logFC<0])
  genes_of_interest <- genes_of_interest[unname(unlist(lapply(genes_of_interest,length)))!=0]

  for (f in grep("_entreznames",grep("^readlist_",ls(),val=T),invert = T,val=T)){
    # print(f)
    vec <- get(f)
    names(vec) <- as.character(na.omit(unique(suppressMessages(select(eval(parse(text=orgDB)), keys=names(vec), columns=c("ENTREZID"), keytype="SYMBOL"))$ENTREZID)))
    if(length(vec)>1){
      assign(paste0(f,"_entreznames"), vec,envir = .GlobalEnv)
    }
  }

  dir.create(paste0(path,"/",tools::file_path_sans_ext(file),"_funct_enrich_clusterProfiler"),showWarnings=F);setwd(paste0(path,"/",tools::file_path_sans_ext(file),"_funct_enrich_clusterProfiler"))
  dir.create(paste0(getwd(),"/go_figs/"),showWarnings=F)
    print(paste0("Processing ",file,"... Gene classification based on GO distribution at a specific level"))
    ### Gene classification based on GO distribution at a specific level:
      for (geneset in names(genes_of_interest)){
        for (levelgo in 1:5){
          # print(paste0("groupGO_level_",levelgo))
          try({
            ggo <- groupGO(gene     = genes_of_interest[[geneset]],
                           OrgDb    = orgDB,
                           keyType  = "SYMBOL",
                           ont      = "BP",
                           level    = levelgo,
                           readable = TRUE)
            # ggo@result
            # ggo@result$Count
            write.table(ggo@result,file=paste0("level_",levelgo,"_",geneset,"_groupGO_BP.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ggo <- groupGO(gene     = genes_of_interest[[geneset]],
                           OrgDb    = orgDB,
                           keyType  = "SYMBOL",
                           ont      = "MF",
                           level    = levelgo,
                           readable = TRUE)
            # ggo@result
            # ggo@result$Count
            write.table(ggo@result,file=paste0("level_",levelgo,"_",geneset,"_groupGO_MF.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ggo <- groupGO(gene     = genes_of_interest[[geneset]],
                         OrgDb    = orgDB,
                         keyType  = "SYMBOL",
                         ont      = "CC",
                         level    = levelgo,
                         readable = TRUE)
            # ggo@result
            # ggo@result$Count
            write.table(ggo@result,file=paste0("level_",levelgo,"_",geneset,"_groupGO_CC.txt"),col.names = T,row.names = F,quote = F,sep="\t")
        },silent=T)
        }
      }

      print(paste0("Processing ",file,"... GO over-representation analyses"))
      ### GO over-representation analyses:
      for (geneset in names(genes_of_interest)){
        for (i in p.adjust.methods){
          try({
            ego <- enrichGO(gene          = genes_of_interest[[geneset]],
                          universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                          OrgDb         = orgDB,
                          keyType  = "SYMBOL",
                          ont           = "ALL",
                          pAdjustMethod = i,
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("GO_overrepresentation_test_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ego <- enrichGO(gene          = genes_of_interest[[geneset]],
                          universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                          OrgDb         = orgDB,
                          keyType  = "SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = i,
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
            write.table(ego,file=paste0("GO_overrepresentation_test_BP_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
            suppressMessages(ggplot2::ggsave(goplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,".pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(barplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_barplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(dotplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_dotplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(cnetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_cnetplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(heatplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_heatplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_treeplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_emapplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(upsetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_upsetplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(pmcplot(ego$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_BP_",i,"_",geneset,"_pmcplot.pdf"),width=30, height=30))
          },silent=T)
          try({
            ego <- enrichGO(gene          = genes_of_interest[[geneset]],
                          universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                          OrgDb         = orgDB,
                          keyType  = "SYMBOL",
                          ont           = "MF",
                          pAdjustMethod = i,
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
            write.table(ego,file=paste0("GO_overrepresentation_test_MF_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
            suppressMessages(ggplot2::ggsave(goplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,".pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(barplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_barplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(dotplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_dotplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(cnetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_cnetplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(heatplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_heatplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_treeplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_emapplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(upsetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_upsetplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(pmcplot(ego$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_MF_",i,"_",geneset,"_pmcplot.pdf"),width=30, height=30))
          },silent=T)
          try({  
            ego <- enrichGO(gene          = genes_of_interest[[geneset]],
                        universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                        OrgDb         = orgDB,
                        keyType  = "SYMBOL",
                        ont           = "CC",
                        pAdjustMethod = i,
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)
            write.table(ego,file=paste0("GO_overrepresentation_test_CC_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
            suppressMessages(ggplot2::ggsave(goplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,".pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(barplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_barplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(dotplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_dotplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(cnetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_cnetplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(heatplot(ego, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_heatplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_treeplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(ego)), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_emapplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(upsetplot(ego), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_upsetplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(pmcplot(ego$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_overrepresentation_test_CC_",i,"_",geneset,"_pmcplot.pdf"),width=30, height=30))
          },silent=T)
          try({
            ego <- enrichGO(gene          = genes_of_interest[[geneset]],
                        universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                        OrgDb         = orgDB,
                        keyType  = "SYMBOL",
                        ont           = "ALL",
                        pAdjustMethod = i,
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 1,
                        readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("GO_overrepresentation_test_",i,"_",geneset,"_no_qval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ego <- enrichGO(gene          = genes_of_interest[[geneset]],
                        universe      = unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL"))))),
                        OrgDb         = orgDB,
                        keyType  = "SYMBOL",
                        ont           = "ALL",
                        pAdjustMethod = i,
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("GO_overrepresentation_test_",i,"_",geneset,"_no_pval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
      }
      }

      print(paste0("Processing ",file,"... Gene Set Enrichment Analysis of Gene Ontology"))
      ### Gene Set Enrichment Analysis of Gene Ontology:
      for (f in grep("entrez",grep("^readlist_",ls(),val=T),invert = T,val=T)){
          for (i in p.adjust.methods){
            try({  
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                                 OrgDb        = orgDB,
                                 ont          = "ALL",
                                 keyType  = "SYMBOL",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = i,
                                 verbose      = FALSE,
                                 by="fgsea"))
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({  
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                               OrgDb        = orgDB,
                               ont          = "BP",
                               keyType  = "SYMBOL",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = i,
                               verbose      = FALSE,
                               by="fgsea"))
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_BP.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              suppressMessages(ggplot2::ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_dotplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_cnetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_heatplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_treeplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_emapplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_upsetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_pmcplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_BP_ridgeplot.pdf"),width=30, height=30))
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                               OrgDb        = orgDB,
                               ont          = "MF",
                               keyType  = "SYMBOL",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = i,
                               verbose      = FALSE,
                               by="fgsea"))
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_MF.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              suppressMessages(ggplot2::ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_dotplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_cnetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_heatplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_treeplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_emapplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_upsetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_pmcplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_MF_ridgeplot.pdf"),width=30, height=30))
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                               OrgDb        = orgDB,
                               ont          = "CC",
                               keyType  = "SYMBOL",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = i,
                               verbose      = FALSE,
                               by="fgsea"))
          
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_CC.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              suppressMessages(ggplot2::ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_dotplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_cnetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_heatplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_treeplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_emapplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_upsetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_pmcplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_fgsea_CC_ridgeplot.pdf"),width=30, height=30))
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                               OrgDb        = orgDB,
                               ont          = "ALL",
                               keyType  = "SYMBOL",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = i,
                               verbose      = FALSE,
                               by="DOSE",
                               nPerm = 100))
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                              OrgDb        = orgDB,
                              ont          = "BP",
                              keyType  = "SYMBOL",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = i,
                              verbose      = FALSE,
                              by="DOSE",
                              nPerm = 100))
          
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_BP.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              suppressMessages(ggplot2::ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_dotplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_cnetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_heatplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_treeplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_emapplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_upsetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_pmcplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_BP_ridgeplot.pdf"),width=30, height=30))
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                              OrgDb        = orgDB,
                              ont          = "MF",
                              keyType  = "SYMBOL",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = i,
                              verbose      = FALSE,
                              by="DOSE",
                              nPerm = 100))
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_MF.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              suppressMessages(ggplot2::ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_dotplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_cnetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_heatplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_treeplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_emapplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_upsetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_pmcplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_MF_ridgeplot.pdf"),width=30, height=30))
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                              OrgDb        = orgDB,
                              ont          = "CC",
                              keyType  = "SYMBOL",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = i,
                              verbose      = FALSE,
                              by="DOSE",
                              nPerm = 100))
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_CC.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              suppressMessages(ggplot2::ggsave(goplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(dotplot(gse_enrich, showCategory=20), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_dotplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(cnetplot(gse_enrich,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_cnetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(heatplot(gse_enrich, showCategory=20,foldChange=sort(eval(parse(text=f)),decreasing = T)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_heatplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(treeplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_treeplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(gse_enrich)), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_emapplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(upsetplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_upsetplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(pmcplot(gse_enrich@result$Description[1:10], 2010:paste0("20",unlist(lapply(strsplit(date(),"20"),function(x){x[2]})))), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_pmcplot.pdf"),width=30, height=30))
              suppressMessages(ggplot2::ggsave(ridgeplot(gse_enrich), filename = paste0(getwd(),"/go_figs/","GO_GSEA_",f,"_",i,"_DOSE_CC_ridgeplot.pdf"),width=30, height=30))
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                               OrgDb        = orgDB,
                               ont          = "ALL",
                               keyType  = "SYMBOL",
                               pvalueCutoff = 1,
                               pAdjustMethod = i,
                               verbose      = FALSE,
                               by="fgsea"))
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_fgsea_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)    
            try({
              gse_enrich <- suppressMessages(gseGO(geneList= sort(eval(parse(text=f)),decreasing = T),
                               OrgDb        = orgDB,
                               ont          = "ALL",
                               keyType  = "SYMBOL",
                               pvalueCutoff = 1,
                               pAdjustMethod = i,
                               verbose      = FALSE,
                               by="DOSE",
                               nPerm = 100))
              write.table(gse_enrich@result,file=paste0("GO_GSEA_",f,"_",i,"_DOSE_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
      }
      }

      print(paste0("Processing ",file,"... KEGG over-representation"))
      dir.create(paste0(getwd(),"/kegg_paths_snapshots"),showWarnings=F)
      ### KEGG over-representation:
      for (geneset in names(genes_of_interest)){
        entrez_ids <- as.character(na.omit(unique(suppressMessages(select(eval(parse(text=orgDB)), keys=genes_of_interest[[geneset]], columns=c("ENTREZID"), keytype="SYMBOL"))$ENTREZID)))
        for (i in p.adjust.methods){
          try({
            kk <- suppressMessages(enrichKEGG(gene= entrez_ids,
                             organism     = org,
                             pAdjustMethod = i,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05))
            suppressMessages(ggplot2::ggsave(barplot(kk, showCategory=20), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_barplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(dotplot(kk, showCategory=20), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_dotplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(cnetplot(kk), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_cnetplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(heatplot(kk, showCategory=20), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_heatplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(emapplot(pairwise_termsim(kk)), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_emapplot.pdf"),width=30, height=30))
            suppressMessages(ggplot2::ggsave(upsetplot(kk), filename = paste0(getwd(),"/kegg_paths_snapshots/","KEGG_overrepresentation_test",i,"_",geneset,"_upsetplot.pdf"),width=30, height=30))
            kk_write <- kk@result; kk_write$Gene_ID <- unname(unlist(lapply(sapply(kk_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(kk_write,file=paste0("KEGG_enrich_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            kk <- suppressMessages(enrichKEGG(gene= entrez_ids,
                                              organism     = org,
                                              pAdjustMethod = i,
                                              pvalueCutoff = 0.05,
                                              qvalueCutoff = 1))
            kk_write <- kk@result; kk_write$Gene_ID <- unname(unlist(lapply(sapply(kk_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(kk_write,file=paste0("KEGG_enrich_",i,"_",geneset,"_noqval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            kk <- suppressMessages(enrichKEGG(gene= entrez_ids,
                                              organism     = org,
                                              pAdjustMethod = i,
                                              pvalueCutoff = 1,
                                              qvalueCutoff = 1))
            kk_write <- kk@result; kk_write$Gene_ID <- unname(unlist(lapply(sapply(kk_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(kk_write,file=paste0("KEGG_enrich_",i,"_",geneset,"_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
        }
      }

      print(paste0("Processing ",file,"... KEGG pathways visualization"))
      ### KEGG pathways visualization:
      paths_list <- unique(unlist(lapply(list.files(path = getwd(), pattern = "^KEGG", full.names = TRUE), function(x){
        data <- read.delim(x)
        data_filtered <- data$ID[data$pvalue < 0.05]
      })))
      for (f in paths_list){try({suppressMessages(pathview(gene.data=kk@result,pathway.id=f,species=org,kegg.dir=paste0(getwd(),"/kegg_paths_snapshots")))},silent=T)}
      invisible(file.remove(list.files(path=getwd(),pattern="*.pathview.png",full.names = T)))

      print(paste0("Processing ",file,"... Gene Set Enrichment Analysis of KEGG"))
      ### Gene Set Enrichment Analysis of KEGG:
      for (f in grep("_entreznames$",ls(),val=T)){
          for (i in p.adjust.methods){
            try({
              gse_enrich <- suppressMessages(gseKEGG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                     organism = org,
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = i,
                                                     verbose      = FALSE,
                                                     by="fgsea"))
              write.table(gse_enrich@result,file=paste0("KEGG_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              for (k in gse_enrich@result$ID[gse_enrich@result$pvalue < 0.05]){
                suppressMessages(pathview(gene.data=sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                          pathway.id=k,species=org))
                invisible(file.remove(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),invert=T,val=T)))
                invisible(file.rename(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),val=T),paste0(getwd(),"/kegg_paths_snapshots/",k,"_","KEGG_GSEA_",f,"_",i,"_fgsea.png")))
              }
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseKEGG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                     organism = org,
                                                     pvalueCutoff = 0.05,
                                                     pAdjustMethod = i,
                                                     verbose      = FALSE,
                                                     by="DOSE",
                                                     nPerm = 100))
              write.table(gse_enrich@result,file=paste0("KEGG_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              for (k in gse_enrich@result$ID[gse_enrich@result$pvalue < 0.05]){
                suppressMessages(pathview(gene.data=sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                          pathway.id=k,species=org))
                invisible(file.remove(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),invert=T,val=T)))
                invisible(file.rename(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),val=T),paste0(getwd(),"/kegg_paths_snapshots/",k,"_","KEGG_GSEA_",f,"_",i,"_fgsea.png")))
              }
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseKEGG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                     organism = org,
                                                     pvalueCutoff = 1,
                                                     pAdjustMethod = i,
                                                     verbose      = FALSE,
                                                     by="fgsea"))
              write.table(gse_enrich@result,file=paste0("KEGG_GSEA_",f,"_",i,"_fgsea_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              for (k in gse_enrich@result$ID[gse_enrich@result$pvalue < 0.05]){
                suppressMessages(pathview(gene.data=sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                          pathway.id=k,species=org))
                invisible(file.remove(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),invert=T,val=T)))
                invisible(file.rename(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),val=T),paste0(getwd(),"/kegg_paths_snapshots/",k,"_","KEGG_GSEA_",f,"_",i,"_fgsea.png")))
              }
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseKEGG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                     organism = org,
                                                     pvalueCutoff = 1,
                                                     pAdjustMethod = i,
                                                     verbose      = FALSE,
                                                     by="DOSE",
                                                     nPerm = 100))
              write.table(gse_enrich@result,file=paste0("KEGG_GSEA_",f,"_",i,"_DOSE_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
              for (k in gse_enrich@result$ID[gse_enrich@result$pvalue < 0.05]){
                suppressMessages(pathview(gene.data=sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                          pathway.id=k,species=org))
                invisible(file.remove(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),invert=T,val=T)))
                invisible(file.rename(grep("pathview",list.files(path=getwd(),pattern=paste0(k,".*"),full.names = T),val=T),paste0(getwd(),"/kegg_paths_snapshots/",k,"_","KEGG_GSEA_",f,"_",i,"_fgsea.png")))
              }
            },silent=T)
          }
      }

      print(paste0("Processing ",file,"... KEGG Module over-representation"))
      ### KEGG Module over-representation:
      for (geneset in names(genes_of_interest)){
        entrez_ids <- as.character(na.omit(unique(suppressMessages(select(eval(parse(text=orgDB)), keys=genes_of_interest[[geneset]], columns=c("ENTREZID"), keytype="SYMBOL"))$ENTREZID)))
        for (i in p.adjust.methods){
            try({
              gse_enrich <- suppressMessages(enrichMKEGG(gene= entrez_ids,
                                                         organism = org,
                                                         pvalueCutoff = 0.05,
                                                         pAdjustMethod = i))
              gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
              write.table(gse_enrich_write,file=paste0("MKEGG_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(enrichMKEGG(gene= entrez_ids,
                                                         organism = org,
                                                         pvalueCutoff = 0.05,
                                                         qvalueCutoff = 0.05,
                                                         pAdjustMethod = i))
              gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
              write.table(gse_enrich_write,file=paste0("MKEGG_",i,"_",geneset,"_qval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(enrichMKEGG(gene= entrez_ids,
                                                         organism = org,
                                                         pvalueCutoff = 0.05,
                                                         qvalueCutoff = 1,
                                                         pAdjustMethod = i))
              gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
              write.table(gse_enrich_write,file=paste0("MKEGG_",i,"_",geneset,"_noqval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(enrichMKEGG(gene= entrez_ids,
                                                         organism = org,
                                                         pvalueCutoff = 1,
                                                         qvalueCutoff = 1,
                                                         pAdjustMethod = i))
              gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
              write.table(gse_enrich_write,file=paste0("MKEGG_",i,"_",geneset,"_no_pval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
        }
      }

      print(paste0("Processing ",file,"... Gene Set Enrichment Analysis of KEGG modules"))
      ### Gene Set Enrichment Analysis of KEGG modules:
      for (f in grep("_entreznames$",ls(),val=T)){
        for (i in p.adjust.methods){
          try({
            gse_enrich <- suppressMessages(gseMKEGG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    organism = org,
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="fgsea"))
            write.table(gse_enrich@result,file=paste0("MKEGG_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gseMKEGG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    organism = org,
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="DOSE",
                                                    nPerm = 100))
            write.table(gse_enrich@result,file=paste0("MKEGG_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gseMKEGG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    organism = org,
                                                    pvalueCutoff = 1,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="fgsea"))
            write.table(gse_enrich@result,file=paste0("MKEGG_GSEA_",f,"_",i,"_fgsea_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gseMKEGG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    organism = org,
                                                    pvalueCutoff = 1,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="DOSE",
                                                    nPerm = 100))
            write.table(gse_enrich@result,file=paste0("MKEGG_GSEA_",f,"_",i,"_DOSE_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
        }
      }

      print(paste0("Processing ",file,"... WikiPathways over-representation"))
      ### WikiPathways over-representation:
      for (geneset in names(genes_of_interest)){
        entrez_ids <- as.character(na.omit(unique(suppressMessages(select(eval(parse(text=orgDB)), keys=genes_of_interest[[geneset]], columns=c("ENTREZID"), keytype="SYMBOL"))$ENTREZID)))
        for (i in p.adjust.methods){
          try({
            gse_enrich <- suppressMessages(enrichWP(gene= entrez_ids,
                                                       organism = organism_cp,
                                                       pvalueCutoff = 0.05,
                                                       pAdjustMethod = i))
            gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(gse_enrich_write,file=paste0("WP_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(enrichWP(gene= entrez_ids,
                                                       organism = organism_cp,
                                                       pvalueCutoff = 0.05,
                                                       qvalueCutoff = 0.05,
                                                       pAdjustMethod = i))
            gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(gse_enrich_write,file=paste0("WP_",i,"_",geneset,"_qval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(enrichWP(gene= entrez_ids,
                                                       organism = organism_cp,
                                                       pvalueCutoff = 0.05,
                                                       qvalueCutoff = 1,
                                                       pAdjustMethod = i))
            gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(gse_enrich_write,file=paste0("WP_",i,"_",geneset,"_noqval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(enrichWP(gene= entrez_ids,
                                                       organism = organism_cp,
                                                       pvalueCutoff = 1,
                                                       qvalueCutoff = 1,
                                                       pAdjustMethod = i))
            gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(gse_enrich_write,file=paste0("WP_",i,"_",geneset,"_no_pval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
        }
      }

      print(paste0("Processing ",file,"... Gene set enrichment analyses of WikiPathways"))
      ### Gene set enrichment analyses of WikiPathways: 
      for (f in grep("_entreznames$",ls(),val=T)){
        for (i in p.adjust.methods){
          try({
            gse_enrich <- suppressMessages(gseWP(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    organism = organism_cp,
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="fgsea"))
            write.table(gse_enrich@result,file=paste0("WP_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gseWP(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    organism = organism_cp,
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="DOSE",
                                                    nPerm = 100))
            write.table(gse_enrich@result,file=paste0("WP_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gseWP(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    organism = organism_cp,
                                                    pvalueCutoff = 1,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="fgsea"))
            write.table(gse_enrich@result,file=paste0("WP_GSEA_",f,"_",i,"_fgsea_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gseWP(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    organism = organism_cp,
                                                    pvalueCutoff = 1,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="DOSE",
                                                    nPerm = 100))
            write.table(gse_enrich@result,file=paste0("WP_GSEA_",f,"_",i,"_DOSE_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
        }
      }

      print(paste0("Processing ",file,"... Reactome over-representation"))
      ### Reactome over-representation:
      for (geneset in names(genes_of_interest)){
        entrez_ids <- as.character(na.omit(unique(suppressMessages(select(eval(parse(text=orgDB)), keys=genes_of_interest[[geneset]], columns=c("ENTREZID"), keytype="SYMBOL"))$ENTREZID)))
        for (i in p.adjust.methods){
          try({
            gse_enrich <- suppressMessages(enrichPathway(gene= entrez_ids,
                                                    organism = organism_cp_react,
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = i))
            gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(gse_enrich_write,file=paste0("REACT_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(enrichPathway(gene= entrez_ids,
                                                    organism = organism_cp_react,
                                                    pvalueCutoff = 0.05,
                                                    qvalueCutoff = 0.05,
                                                    pAdjustMethod = i))
            gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(gse_enrich_write,file=paste0("REACT_",i,"_",geneset,"_qval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(enrichPathway(gene= entrez_ids,
                                                    organism = organism_cp_react,
                                                    pvalueCutoff = 0.05,
                                                    qvalueCutoff = 1,
                                                    pAdjustMethod = i))
            gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(gse_enrich_write,file=paste0("REACT_",i,"_",geneset,"_noqval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(enrichPathway(gene= entrez_ids,
                                                    organism = organism_cp_react,
                                                    pvalueCutoff = 1,
                                                    qvalueCutoff = 1,
                                                    pAdjustMethod = i))
            gse_enrich_write <- gse_enrich@result; gse_enrich_write$Gene_ID <- unname(unlist(lapply(sapply(gse_enrich_write$geneID,function(x){strsplit(x,"/")}),function(y){paste(suppressMessages(select(eval(parse(text=orgDB)), keys=y, columns=c("SYMBOL"), keytype="ENTREZID"))$SYMBOL,collapse=",")})))
            write.table(gse_enrich_write,file=paste0("REACT_",i,"_",geneset,"_no_pval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
        }
      }

      print(paste0("Processing ",file,"... Gene set enrichment analyses of Reactome"))
      ### Gene set enrichment analyses of Reactome: 
      for (f in grep("_entreznames$",ls(),val=T)){
        for (i in p.adjust.methods){
          try({
            gse_enrich <- suppressMessages(gsePathway(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                 organism = organism_cp_react,
                                                 pvalueCutoff = 0.05,
                                                 pAdjustMethod = i,
                                                 verbose      = FALSE,
                                                 by="fgsea"))
            write.table(gse_enrich@result,file=paste0("REACT_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gsePathway(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                 organism = organism_cp_react,
                                                 pvalueCutoff = 0.05,
                                                 pAdjustMethod = i,
                                                 verbose      = FALSE,
                                                 by="DOSE",
                                                 nPerm = 100))
            write.table(gse_enrich@result,file=paste0("REACT_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gsePathway(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                 organism = organism_cp_react,
                                                 pvalueCutoff = 1,
                                                 pAdjustMethod = i,
                                                 verbose      = FALSE,
                                                 by="fgsea"))
            write.table(gse_enrich@result,file=paste0("REACT_GSEA_",f,"_",i,"_fgsea_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            gse_enrich <- suppressMessages(gsePathway(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                 organism = organism_cp_react,
                                                 pvalueCutoff = 1,
                                                 pAdjustMethod = i,
                                                 verbose      = FALSE,
                                                 by="DOSE",
                                                 nPerm = 100))
            write.table(gse_enrich@result,file=paste0("REACT_GSEA_",f,"_",i,"_DOSE_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
        }
      }

      print(paste0("Processing ",file,"... Reactome pathways visualization"))
      ### Reactome pathways visualization:
      paths_list <- unique(unlist(lapply(list.files(path = getwd(), pattern = "^REACT", full.names = TRUE), function(x){
        data <- read.delim(x)
        data_filtered <- data$Description[data$pvalue < 0.05]
      })))
      dir.create(paste0(getwd(),"/reactome_paths_snapshots"),showWarnings=F)
      for (f in paths_list){try({suppressMessages(ggplot2::ggsave(viewPathway(f,readable = TRUE,organism = organism_cp_react),filename = paste0(getwd(),"/reactome_paths_snapshots/",gsub(" ","_",substr(f,0,40)),"_Reactome.pdf"),width=30, height=30))},silent=T)}

      print(paste0("Processing ",file,"... Over-representation analyses for human databases (DO, NCG and DGN)"))
      ### Over-representation analyses for human databases (DO, NCG and DGN):
      if(organism_cp=="Homo sapiens"){
        for (geneset in names(genes_of_interest)){
        entrez_ids <- as.character(na.omit(unique(suppressMessages(select(eval(parse(text=orgDB)), keys=genes_of_interest[[geneset]], columns=c("ENTREZID"), keytype="SYMBOL"))$ENTREZID)))
        for (i in p.adjust.methods){
          try({
            ego <- enrichDO(gene          = entrez_ids,
                            ont           = "DO",
                            pAdjustMethod = i,
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("DO_overrepresentation_test_",i,"_",geneset,"_DO.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ego <- enrichDO(gene          = entrez_ids,
                            ont           = "DOLite",
                            pAdjustMethod = i,
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("DO_overrepresentation_test_",i,"_",geneset,"_DOLite.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ego <- enrichDO(gene          = entrez_ids,
                            ont           = "DO",
                            pAdjustMethod = i,
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 1,
                            readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("DO_overrepresentation_test_",i,"_",geneset,"_DO_noqval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ego <- enrichDO(gene          = entrez_ids,
                            ont           = "DOLite",
                            pAdjustMethod = i,
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 1,
                            readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("DO_overrepresentation_test_",i,"_",geneset,"_DOLite_noqval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ego <- enrichDO(gene          = entrez_ids,
                            ont           = "DO",
                            pAdjustMethod = i,
                            pvalueCutoff  = 1,
                            qvalueCutoff  = 1,
                            readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("DO_overrepresentation_test_",i,"_",geneset,"_DO_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
          try({
            ego <- enrichDO(gene          = entrez_ids,
                            ont           = "DOLite",
                            pAdjustMethod = i,
                            pvalueCutoff  = 1,
                            qvalueCutoff  = 1,
                            readable      = TRUE)
            # head(ego)
            write.table(ego,file=paste0("DO_overrepresentation_test_",i,"_",geneset,"_DOLite_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
          },silent=T)
        }
        }
        
        for (geneset in names(genes_of_interest)){
          entrez_ids <- as.character(na.omit(unique(suppressMessages(select(eval(parse(text=orgDB)), keys=genes_of_interest[[geneset]], columns=c("ENTREZID"), keytype="SYMBOL"))$ENTREZID)))
          for (i in p.adjust.methods){
            try({
              ego <- enrichNCG(gene          = entrez_ids,
                              pAdjustMethod = i,
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)
              # head(ego)
              write.table(ego,file=paste0("NGC_overrepresentation_test_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              ego <- enrichNCG(gene          = entrez_ids,
                               pAdjustMethod = i,
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 1,
                               readable      = TRUE)
              # head(ego)
              write.table(ego,file=paste0("NGC_overrepresentation_test_",i,"_",geneset,"_noqval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              ego <- enrichNCG(gene          = entrez_ids,
                               pAdjustMethod = i,
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1,
                               readable      = TRUE)
              # head(ego)
              write.table(ego,file=paste0("NGC_overrepresentation_test_",i,"_",geneset,"_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            
          }
        }
        
        for (geneset in names(genes_of_interest)){
          entrez_ids <- as.character(na.omit(unique(suppressMessages(select(eval(parse(text=orgDB)), keys=genes_of_interest[[geneset]], columns=c("ENTREZID"), keytype="SYMBOL"))$ENTREZID)))
          for (i in p.adjust.methods){
            try({
              ego <- enrichDGN(gene          = entrez_ids,
                               pAdjustMethod = i,
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE)
              # head(ego)
              write.table(ego,file=paste0("DGN_overrepresentation_test_",i,"_",geneset,".txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              ego <- enrichDGN(gene          = entrez_ids,
                               pAdjustMethod = i,
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 1,
                               readable      = TRUE)
              # head(ego)
              write.table(ego,file=paste0("DGN_overrepresentation_test_",i,"_",geneset,"_noqval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              ego <- enrichDGN(gene          = entrez_ids,
                               pAdjustMethod = i,
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1,
                               readable      = TRUE)
              # head(ego)
              write.table(ego,file=paste0("DGN_overrepresentation_test_",i,"_",geneset,"_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            
          }
        }
        
        for (f in grep("_entreznames$",ls(),val=T)){
          for (i in p.adjust.methods){
            try({
              gse_enrich <- suppressMessages(gseDO(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                        pvalueCutoff = 0.05,
                                                        pAdjustMethod = i,
                                                        verbose      = FALSE,
                                                        by="fgsea"))
              write.table(gse_enrich@result,file=paste0("DO_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseDO(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                        pvalueCutoff = 0.05,
                                                        pAdjustMethod = i,
                                                        verbose      = FALSE,
                                                        by="DOSE",
                                                        nPerm = 100))
              write.table(gse_enrich@result,file=paste0("DO_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseDO(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                        pvalueCutoff = 1,
                                                        pAdjustMethod = i,
                                                        verbose      = FALSE,
                                                        by="fgsea"))
              write.table(gse_enrich@result,file=paste0("DO_GSEA_",f,"_",i,"_fgsea_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseDO(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                        pvalueCutoff = 1,
                                                        pAdjustMethod = i,
                                                        verbose      = FALSE,
                                                        by="DOSE",
                                                        nPerm = 100))
              write.table(gse_enrich@result,file=paste0("DO_GSEA_",f,"_",i,"_DOSE_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
          }
        }
        
        for (f in grep("_entreznames$",ls(),val=T)){
          for (i in p.adjust.methods){
            try({
              gse_enrich <- suppressMessages(gseNCG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                   pvalueCutoff = 0.05,
                                                   pAdjustMethod = i,
                                                   verbose      = FALSE,
                                                   by="fgsea"))
              write.table(gse_enrich@result,file=paste0("NCG_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseNCG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                   pvalueCutoff = 0.05,
                                                   pAdjustMethod = i,
                                                   verbose      = FALSE,
                                                   by="DOSE",
                                                   nPerm = 100))
              write.table(gse_enrich@result,file=paste0("NCG_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseNCG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                   pvalueCutoff = 1,
                                                   pAdjustMethod = i,
                                                   verbose      = FALSE,
                                                   by="fgsea"))
              write.table(gse_enrich@result,file=paste0("NCG_GSEA_",f,"_",i,"_fgsea_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseNCG(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                   pvalueCutoff = 1,
                                                   pAdjustMethod = i,
                                                   verbose      = FALSE,
                                                   by="DOSE",
                                                   nPerm = 100))
              write.table(gse_enrich@result,file=paste0("NCG_GSEA_",f,"_",i,"_DOSE_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
          }
        }
        
        for (f in grep("_entreznames$",ls(),val=T)){
          for (i in p.adjust.methods){
            try({
              gse_enrich <- suppressMessages(gseDGN(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="fgsea"))
              write.table(gse_enrich@result,file=paste0("DGN_GSEA_",f,"_",i,"_fgsea.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseDGN(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    pvalueCutoff = 0.05,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="DOSE",
                                                    nPerm = 100))
              write.table(gse_enrich@result,file=paste0("DGN_GSEA_",f,"_",i,"_DOSE.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseDGN(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    pvalueCutoff = 1,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="fgsea"))
              write.table(gse_enrich@result,file=paste0("DGN_GSEA_",f,"_",i,"_fgsea_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
            try({
              gse_enrich <- suppressMessages(gseDGN(geneList= sort(eval(parse(text=f)),decreasing=T)[!is.na(names(sort(eval(parse(text=f)),decreasing=T)))],
                                                    pvalueCutoff = 1,
                                                    pAdjustMethod = i,
                                                    verbose      = FALSE,
                                                    by="DOSE",
                                                    nPerm = 100))
              write.table(gse_enrich@result,file=paste0("DGN_GSEA_",f,"_",i,"_DOSE_nopval.txt"),col.names = T,row.names = F,quote = F,sep="\t")
            },silent=T)
          }
        }
      }

      print(paste0("Processing ",file,"... Tidying up..."))
      ### Remove empty files and tar everything:
      system(paste0("find ", getwd(), " -type f -exec awk 'NR>1{exit 1}' {} \\; -exec rm -f {} \\;"))
      tar(paste0(path,"/clusterProfiler_out.tar.gz"), compression = "gzip")
      invisible(file.remove(list.files(getwd(), include.dirs = F, full.names = T, recursive = T)))
      invisible(file.rename(paste0(path,"/clusterProfiler_out.tar.gz"),paste0(getwd(),"/clusterProfiler_out.tar.gz")))
}

### Parallel processing:
files <- list.files(path=path,pattern="^DGE_analysis_comp[0-9]+.txt$")
print(paste0("Processing ",paste(files,collapse=","),"... Many figures are going to be showed, but please double check the parameters used (default) and the possibilities within https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html. it is likely that some of the plots should be redo manually to optimize them..."))
print(paste0("Processing ",paste(files,collapse=","),"... Remember that to visualize particular pathways of interest in KEGG, the functions pathway() and browseKEGG can be used, https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html#visualize-enriched-kegg-pathways"))
print(paste0("Processing ",paste(files,collapse=","),"... Here all the pathways that may be of interest in KEGG (i..e pval < 0.05) are shown for all alternatives"))
print(paste0("Processing ",paste(files,collapse=","),"... Remember that to visualize particular pathways of interest in Reactome, the function viewPathway() can be used, https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html#pathway-visualization"))
print(paste0("Processing ",paste(files,collapse=","),"... Here all the pathways that may be of interest in Reactome (i..e pval < 0.05) are shown, but without coloring because it would change for each case... please redo if required adding the parameter 'foldChange=' in the function"))
mclapply(
    mc.cores = cores,
    X = files,
    FUN = process_file
)
print("ALL DONE")