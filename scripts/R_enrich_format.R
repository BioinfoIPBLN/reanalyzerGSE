#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
table_enrich <- args[1]
table_dge <- args[2]
organism <- args[3]
rev_thres <- as.numeric(args[4])
suppressMessages(library(data.table,quiet = T,warn.conflicts = F))

a <- data.table::fread(table_enrich)

### Add up and down:
if ("geneID" %in% colnames(a) && !("Gene_ID_up" %in% colnames(a)) && file.exists(table_dge)) {
  b <- data.table::fread(table_dge)
  fc_col <- grep("^logFC", colnames(b), value = TRUE)[1]
  if ("Gene_ID" %in% colnames(b) && !is.na(fc_col)) {
    ids <- toupper(as.character(b$Gene_ID))
    fc <- suppressWarnings(as.numeric(b[[fc_col]]))
    ids_unversioned <- sub("\\.[0-9]+$", "", ids)
    up <- unique(c(ids[which(fc > 0)], ids_unversioned[which(fc > 0)]))
    down <- unique(c(ids[which(fc < 0)], ids_unversioned[which(fc < 0)]))
    genes <- lapply(strsplit(gsub(";", "/", ifelse(is.na(a$geneID), "", as.character(a$geneID))), "/"), function(x) x[nzchar(x)])
    pick <- function(set) vapply(genes, function(x) { s <- paste(x[toupper(x) %in% set], collapse = "/"); if (nzchar(s)) s else "-" }, character(1))
    count <- function(v) vapply(strsplit(v, "/"), function(x) sum(nzchar(x) & x != "-"), integer(1))
    a$Gene_ID_num <- vapply(genes, length, integer(1))
    a$Gene_ID_up <- pick(up)
    a$Gene_ID_up_num <- count(a$Gene_ID_up)
    a$Gene_ID_down <- pick(down)
    a$Gene_ID_down_num <- count(a$Gene_ID_down)
    write.table(a, file = table_enrich, quote = F, row.names = F, col.names = T, sep = "\t")
  }
}

### Apply Revigo:
organism_cp <- gsub("_"," ",organism)
orgDB <- switch(organism_cp, "Homo sapiens" = "org.Hs.eg.db", "Mus musculus" = "org.Mm.eg.db", NULL)
table_base <- basename(table_enrich)

try({
  if (!is.null(orgDB) && grepl("^GO_", table_base, ignore.case = TRUE) && !file.exists(paste0(table_enrich, "_revigo.pdf"))) {
    suppressMessages(library(rrvgo,quiet = T,warn.conflicts = F))
    ontology <- if (grepl("Cellular|_CC", table_base, ignore.case = TRUE)) "CC" else
                if (grepl("Molecular|_MF", table_base, ignore.case = TRUE)) "MF" else
                if (grepl("Biological|_BP", table_base, ignore.case = TRUE)) "BP" else NA
    a <- as.data.frame(data.table::fread(table_enrich))
    pval_col <- grep("^pval|^p.val|^P.valu", colnames(a), value = TRUE)[1]
    if (!is.na(ontology) && "Term" %in% colnames(a) && !is.na(pval_col)) {
      print(paste0("Applying REVIGO with threshold ",rev_thres," for ",organism," and table ", table_enrich," ..."))
      ids <- gsub(")","",gsub(".*GO:","GO:",a$Term))
      simMatrix <- calculateSimMatrix(ids,
                                      orgdb=orgDB,
                                      ont=ontology,
                                      method="Rel")
      scores <- setNames(-log10(a[[pval_col]]), ids)
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
},silent=T)
