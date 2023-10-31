#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
table_enrich <- args[1]
table_dge <- args[2]
organism <- args[3]
rev_thres <- as.numeric(args[4])
suppressMessages(library(data.table,quiet = T,warn.conflicts = F))

a <- data.table::fread(table_enrich); b <- data.table::fread(table_dge)

if(!("geneID" %in% colnames(a))){
  colnames(a)[dim(a)[2]] <- "geneID"
}

### Add up and down:

a$geneID <- gsub(";","/",a$geneID)
b$Gene_ID <- gsub("\\..*","",b$Gene_ID)
b_Gene_ID_up <- b$Gene_ID[b$logFC > 0]

a$Gene_ID_num <- unlist(lapply(strsplit(a$geneID,"/"),length))
a$Gene_ID_up <- gsub("/-","",unlist(lapply(strsplit(a$geneID,"/"),function(x){paste(c(x[toupper(x) %in% toupper(b_Gene_ID_up)],"-"),collapse="/")})))
a$Gene_ID_up_num <- unlist(lapply(strsplit(gsub("-","",a$Gene_ID_up),"/"),length))
a$Gene_ID_down <- gsub("/-","",unlist(lapply(strsplit(a$geneID,"/"),function(x){paste(c(x[!(toupper(x) %in% toupper(b_Gene_ID_up))],"-"),collapse="/")})))
a$Gene_ID_down_num <- unlist(lapply(strsplit(gsub("-","",a$Gene_ID_down),"/"),length))


### Apply Revigo:
suppressMessages(library(rrvgo,quiet = T,warn.conflicts = F))
organism_cp <- gsub("_"," ",organism)
if(organism_cp=="Homo sapiens"){
  orgDB <- "org.Hs.eg.db"
} else if(organism_cp=="Mus musculus"){
  orgDB <- "org.Mm.eg.db"  
}

if (exists("orgDB")){
  a <- data.table::fread(table_enrich)
  if(grepl("Cellular",table_enrich) || grepl("CC",table_enrich)){
    ontology="CC"
  } else if(grepl("Molecular",table_enrich) || grepl("MF",table_enrich)){
    ontology="MF"
  } else if(grepl("Biological",table_enrich) || grepl("BP",table_enrich)){
    ontology="BP"
  }

  if (grepl("/GO_",table_enrich)){  
    if ("Term" %in% colnames(a)){
      ids=gsub(")","",gsub(".*GO:","GO:",a$Term))
      simMatrix <- calculateSimMatrix(ids,
                                      orgdb=orgDB,
                                      ont=ontology,
                                      method="Rel")
      scores <- setNames(-log10(a[,grep("^pval|^p.val|^P.valu",colnames(a),val=T)]), ids)
      reducedTerms <- reduceSimMatrix(simMatrix,
                                      scores,
                                      threshold=rev_thres,
                                      orgdb=orgDB)

      pdf(paste0(table_enrich,"_revigo.pdf"),paper="a4")
      heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

      scatterPlot(simMatrix, reducedTerms)
      treemapPlot(reducedTerms)
      wordcloudPlot(reducedTerms, min.freq=1, colors="black")
      dev.off()
    }
  }
}