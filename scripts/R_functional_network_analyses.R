#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
expr <- args[2]
pattern_search <- args[3]
organism <- as.numeric(args[4])
wgcna_mode <- if (length(args) >= 5) tolower(args[5]) else "all"  # "all" (canonical) or "degs"
organism_name <- if (length(args) >= 6) args[6] else ""
samples_info <- if (length(args) >= 7) args[7] else ""   # reads_study_info/samples_info.txt

# Condition per sample, taken from samples_info.txt (col1 = sample, col3 = condition).
# Falls back to stripping "_Rep..." off the sample name when the file is absent or a
# sample is missing from it, which is what this script did for every sample before.
sample_conditions <- function(sample_names) {
    fallback <- gsub("_Rep.*", "", sample_names)
    if (!nzchar(samples_info) || !file.exists(samples_info)) return(fallback)
    si <- tryCatch(as.data.frame(data.table::fread(samples_info, header = FALSE)),
                   error = function(e) NULL)
    if (is.null(si) || ncol(si) < 3) return(fallback)
    cond <- as.character(si[[3]])[match(sample_names, as.character(si[[1]]))]
    if (all(is.na(cond))) return(fallback)
    cond[is.na(cond)] <- fallback[is.na(cond)]
    cond
}

# `expr` is the RAW count table. Co-expression needs a homoscedastic matrix, so the
# counts are variance-stabilised with DESeq2 before WGCNA and BONOBO see them.
# vst() needs enough genes to fit dispersions on; small or sparse matrices fall back
# to the exact transformation, and then to log2(CPM+1) so the step never dies here.
load_expression <- function(counts_file) {
    raw <- as.data.frame(data.table::fread(counts_file))
    m <- as.matrix(raw[, grep("Gene_ID|Length", colnames(raw), invert = TRUE), drop = FALSE])
    rownames(m) <- raw$Gene_ID
    m <- round(m); m[is.na(m)] <- 0; mode(m) <- "integer"
    trans <- tryCatch({
        v <- DESeq2::vst(m, blind = TRUE)
        cat("Expression for network analyses: DESeq2::vst on raw counts\n"); v
    }, error = function(e) {
        cat(paste0("  vst() not applicable (", conditionMessage(e),
                   "); trying varianceStabilizingTransformation()\n"))
        tryCatch({
            v <- DESeq2::varianceStabilizingTransformation(m, blind = TRUE)
            cat("Expression for network analyses: DESeq2::varianceStabilizingTransformation\n"); v
        }, error = function(e2) {
            cat(paste0("  VST unavailable (", conditionMessage(e2),
                       "); falling back to log2(CPM+1)\n"))
            log2(edgeR::cpm(m) + 1)
        })
    })
    data.frame(Gene_ID = rownames(trans), trans, check.names = FALSE)
}

.expr_cache <- NULL
get_expression <- function() {
    if (is.null(.expr_cache)) .expr_cache <<- load_expression(expr)
    .expr_cache
}

suppressMessages(library(WGCNA, quiet = T, warn.conflicts = F))
suppressMessages(library(tidyr, quiet = T, warn.conflicts = F))
suppressMessages(library(STRINGdb, quiet = T, warn.conflicts = F))
suppressMessages(library(ggplot2, quiet = T, warn.conflicts = F))
suppressMessages(library(clusterProfiler, quiet = T, warn.conflicts = F))
suppressMessages(library(AnnotationDbi, quiet = T, warn.conflicts = F))

# Organism-aware OrgDb selection
orgdb <- NULL
kegg_org <- NULL
species_label <- NULL
if (grepl("Mus", organism_name, ignore.case = TRUE)) {
  suppressMessages(library(org.Mm.eg.db, quiet = T, warn.conflicts = F))
  orgdb <- org.Mm.eg.db; kegg_org <- "mmu"; species_label <- "Mus musculus"
  cat("Organism detected: Mouse (Mus musculus)\n")
} else if (grepl("Homo", organism_name, ignore.case = TRUE)) {
  suppressMessages(library(org.Hs.eg.db, quiet = T, warn.conflicts = F))
  orgdb <- org.Hs.eg.db; kegg_org <- "hsa"; species_label <- "Homo sapiens"
  cat("Organism detected: Human (Homo sapiens)\n")
} else {
  cat(paste0("Non-model organism ('", organism_name, "'): GSEA and module enrichment will be skipped\n"))
}

# Create output base directory under DGE/network_analyses/
out_base <- paste0(path, "/network_analyses/")
system(paste("mkdir -p", out_base))


################################################################################
#### Helper: map gene IDs to Entrez IDs
################################################################################
map_to_entrez <- function(gene_ids, orgdb) {
    if (is.null(orgdb)) return(NULL)
    # Try ENSEMBL first (strip version suffix)
    ids_clean <- gsub("\\..*", "", toupper(gene_ids))
    mapped <- tryCatch({
        res <- suppressMessages(AnnotationDbi::select(orgdb,
                                   keys = ids_clean,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "ENSEMBL"))
        res <- res[!is.na(res$ENTREZID), ]
        res <- res[!duplicated(res$ENSEMBL), ]
        res
    }, error = function(e) NULL)
    
    # If ENSEMBL mapping failed or yielded nothing, try SYMBOL
    if (is.null(mapped) || nrow(mapped) == 0) {
        mapped <- tryCatch({
            res <- suppressMessages(AnnotationDbi::select(orgdb,
                                       keys = gene_ids,
                                       columns = c("ENTREZID", "ENSEMBL"),
                                       keytype = "SYMBOL"))
            res <- res[!is.na(res$ENTREZID), ]
            res <- res[!duplicated(res$SYMBOL), ]
            colnames(res)[colnames(res) == "SYMBOL"] <- "input_id"
            res
        }, error = function(e) NULL)
    } else {
        colnames(mapped)[colnames(mapped) == "ENSEMBL"] <- "input_id"
    }
    return(mapped)
}


################################################################################
#### Helper: run GSEA
################################################################################
run_gsea <- function(deg_df, label, gsea_base_dir, orgdb, kegg_org, species_label) {
    if (is.null(orgdb)) {
        cat("  Skipping GSEA: no OrgDb available for this organism\n")
        return(invisible(NULL))
    }
    if (nrow(deg_df) == 0) {
        cat("  Skipping GSEA: no genes in DEG table\n")
        return(invisible(NULL))
    }
    
    # Map gene IDs to Entrez
    gene_ids <- deg_df$Gene_ID
    entrez_map <- map_to_entrez(gene_ids, orgdb)
    if (is.null(entrez_map) || nrow(entrez_map) == 0) {
        cat("  Skipping GSEA: could not map gene IDs to Entrez\n")
        return(invisible(NULL))
    }
    
    # Merge with DEG data
    deg_mapped <- merge(deg_df, entrez_map, by.x = "Gene_ID", by.y = "input_id")
    deg_mapped <- deg_mapped[!is.na(deg_mapped$ENTREZID) & !is.na(deg_mapped$logFC), ]
    
    if (nrow(deg_mapped) < 10) {
        cat("  Skipping GSEA: too few mapped genes\n")
        return(invisible(NULL))
    }
    
    # Helper to save GSEA results
    save_gsea_results <- function(gsea_obj, db_name, out_dir, metric_label) {
        if (is.null(gsea_obj)) return()
        res_df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
        if (is.null(res_df) || nrow(res_df) == 0) {
            cat(paste0("    No pathways found for: ", db_name, "\n"))
            return()
        }
        write.table(res_df, file = paste0(out_dir, "GSEA_", db_name, "_", metric_label, ".txt"),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        
        sig_df <- res_df[res_df$p.adjust < 0.05, ]
        if (nrow(sig_df) > 0) {
            tryCatch({
                gsea_sig <- gsea_obj
                gsea_sig@result <- sig_df
                p <- enrichplot::dotplot(gsea_sig, showCategory = min(20, nrow(sig_df)),
                                        title = paste0(db_name, " - ", metric_label))
                ggsave(paste0(out_dir, "GSEA_", db_name, "_", metric_label, ".png"),
                       p, width = 10, height = 8)
                cat(paste0("    ", db_name, ": ", nrow(sig_df), " significant pathways (dotplot saved)\n"))
            }, error = function(e) {
                cat(paste0("    ", db_name, ": ", nrow(sig_df), " significant pathways (dotplot failed: ", e$message, ")\n"))
            })
        } else {
            cat(paste0("    ", db_name, ": table saved, no significant pathways (p.adjust < 0.05)\n"))
        }
    }
    
    # Run GSEA with two ranking metrics
    for (metric in c("log2FC", "FDR")) {
        cat(paste0("\n  GSEA [", label, "] ranking by ", metric, "...\n"))
        
        if (metric == "FDR") {
            deg_mapped$rank_val <- -log10(deg_mapped$FDR + 1e-10) * sign(deg_mapped$logFC)
        } else {
            deg_mapped$rank_val <- deg_mapped$logFC
        }
        
        # Build ranked list: deduplicate Entrez IDs, keep highest |rank_val|
        ranked_df <- deg_mapped[order(-abs(deg_mapped$rank_val)), ]
        ranked_df <- ranked_df[!duplicated(ranked_df$ENTREZID), ]
        ranked <- sort(setNames(ranked_df$rank_val, ranked_df$ENTREZID), decreasing = TRUE)
        ranked <- ranked[!is.na(ranked)]
        
        if (length(ranked) < 10) {
            cat("    Too few ranked genes, skipping this metric\n")
            next
        }
        
        metric_label <- paste0(label, "_", metric)
        out_dir <- paste0(gsea_base_dir, metric, "/")
        system(paste("mkdir -p", out_dir))
        
        # GO Biological Process
        gsea_go_bp <- tryCatch(gseGO(geneList = ranked, OrgDb = orgdb, ont = "BP",
                                      keyType = "ENTREZID", minGSSize = 15, maxGSSize = 500,
                                      pvalueCutoff = 1, eps = 0, verbose = FALSE),
                               error = function(e) { cat(paste0("    GO BP error: ", e$message, "\n")); NULL })
        save_gsea_results(gsea_go_bp, "GO_BP", out_dir, metric_label)
        
        # GO Molecular Function
        gsea_go_mf <- tryCatch(gseGO(geneList = ranked, OrgDb = orgdb, ont = "MF",
                                      keyType = "ENTREZID", minGSSize = 15, maxGSSize = 500,
                                      pvalueCutoff = 1, eps = 0, verbose = FALSE),
                               error = function(e) { cat(paste0("    GO MF error: ", e$message, "\n")); NULL })
        save_gsea_results(gsea_go_mf, "GO_MF", out_dir, metric_label)
        
        # GO Cellular Component
        gsea_go_cc <- tryCatch(gseGO(geneList = ranked, OrgDb = orgdb, ont = "CC",
                                      keyType = "ENTREZID", minGSSize = 15, maxGSSize = 500,
                                      pvalueCutoff = 1, eps = 0, verbose = FALSE),
                               error = function(e) { cat(paste0("    GO CC error: ", e$message, "\n")); NULL })
        save_gsea_results(gsea_go_cc, "GO_CC", out_dir, metric_label)
        
        # KEGG
        gsea_kegg <- tryCatch(gseKEGG(geneList = ranked, organism = kegg_org,
                                       minGSSize = 15, maxGSSize = 500,
                                       pvalueCutoff = 1, eps = 0, verbose = FALSE),
                              error = function(e) { cat(paste0("    KEGG error: ", e$message, "\n")); NULL })
        save_gsea_results(gsea_kegg, "KEGG", out_dir, metric_label)
        
        # Reactome (via ReactomePA)
        gsea_reactome <- tryCatch({
            suppressMessages(library(ReactomePA, quiet = T, warn.conflicts = F))
            ReactomePA::gsePathway(geneList = ranked, organism = species_label,
                                   minGSSize = 15, maxGSSize = 500,
                                   pvalueCutoff = 1, eps = 0, verbose = FALSE)
        }, error = function(e) { cat(paste0("    Reactome error: ", e$message, "\n")); NULL })
        save_gsea_results(gsea_reactome, "Reactome", out_dir, metric_label)
        
        # MSigDB Canonical Pathways (C2:CP) via msigdbr
        gsea_msigdb_cp <- tryCatch({
            suppressMessages(library(msigdbr, quiet = T, warn.conflicts = F))
            cp_gs <- msigdbr(species = species_label, collection = "C2", subcollection = "CP")
            gene_col <- if ("ncbi_gene" %in% colnames(cp_gs)) "ncbi_gene" else "entrez_gene"
            term2gene_cp <- cp_gs[, c("gs_name", gene_col)]
            colnames(term2gene_cp)[colnames(term2gene_cp) == gene_col] <- "entrez_gene"
            term2gene_cp$entrez_gene <- as.character(term2gene_cp$entrez_gene)
            GSEA(geneList = ranked, TERM2GENE = term2gene_cp,
                 minGSSize = 15, maxGSSize = 500,
                 pvalueCutoff = 1, eps = 0, verbose = FALSE)
        }, error = function(e) { cat(paste0("    MSigDB CP error: ", e$message, "\n")); NULL })
        save_gsea_results(gsea_msigdb_cp, "MSigDB_CP", out_dir, metric_label)

        # MSigDB Reactome (C2:CP:REACTOME) via msigdbr
        gsea_msigdb_reactome <- tryCatch({
            suppressMessages(library(msigdbr, quiet = T, warn.conflicts = F))
            reactome_gs <- msigdbr(species = species_label, collection = "C2", subcollection = "CP:REACTOME")
            gene_col <- if ("ncbi_gene" %in% colnames(reactome_gs)) "ncbi_gene" else "entrez_gene"
            term2gene_reactome <- reactome_gs[, c("gs_name", gene_col)]
            colnames(term2gene_reactome)[colnames(term2gene_reactome) == gene_col] <- "entrez_gene"
            term2gene_reactome$entrez_gene <- as.character(term2gene_reactome$entrez_gene)
            GSEA(geneList = ranked, TERM2GENE = term2gene_reactome,
                 minGSSize = 15, maxGSSize = 500,
                 pvalueCutoff = 1, eps = 0, verbose = FALSE)
        }, error = function(e) { cat(paste0("    MSigDB Reactome error: ", e$message, "\n")); NULL })
        save_gsea_results(gsea_msigdb_reactome, "MSigDB_Reactome", out_dir, metric_label)
    }
    
    cat(paste0("  GSEA complete for: ", label, "\n"))
    return(invisible(NULL))
}


################################################################################
#### Helper: run WGCNA analysis
################################################################################
run_wgcna <- function(nexpr_mat, new_path, label, orgdb = NULL, all_gene_ids = NULL, minSize = 30, MEDissThres = 0.25) {
    nexpr_mat <- nexpr_mat[rowSums(nexpr_mat, na.rm = TRUE) > 0, , drop = FALSE]
    if (nrow(nexpr_mat) < minSize) {
        cat(paste0("Skipping WGCNA for ", label, ": only ", nrow(nexpr_mat),
                   " genes (minimum ", minSize, " required)\n"))
        return(invisible(NULL))
    }

    allowWGCNAThreads()
    system(paste("mkdir -p", new_path))

    t_exprs <- t(nexpr_mat)
    gsg <- goodSamplesGenes(t_exprs, verbose = 0)
    if (!gsg$allOK) {
        cat(paste0("WGCNA QC for ", label, ": dropping ", sum(!gsg$goodGenes), " gene(s) and ",
                   sum(!gsg$goodSamples), " sample(s) with zero variance or too many missing values\n"))
        t_exprs <- t_exprs[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    }
    if (ncol(t_exprs) < minSize) {
        cat(paste0("Skipping WGCNA for ", label, ": only ", ncol(t_exprs),
                   " genes left after QC (minimum ", minSize, " required)\n"))
        return(invisible(NULL))
    }

    ## Choose a set of soft-thresholding powers
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold(t_exprs, powerVector = powers)
    if (is.na(sft$powerEstimate)){
      power = sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
    }  else {power = sft$powerEstimate}

    pindex=which(sft$fitIndices$Power==power)

    # Plot results
    pdf(file =  paste(new_path,"WGCNA_plots_",label,".pdf",sep=""),paper="a4")
    par(mfrow = c(1,2))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[-pindex,1], -sign(sft$fitIndices[-pindex,3])*sft$fitIndices[-pindex,2],
         labels=powers[-pindex],col="black", cex=0.7)
    text(sft$fitIndices[pindex,1], -sign(sft$fitIndices[pindex,3])*sft$fitIndices[pindex,2],
         labels=power,col="red")
    abline(h=0.90,col="blue")

    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[-pindex,1], sft$fitIndices[-pindex,5], labels=powers[-pindex],col="black")
    text(sft$fitIndices[pindex,1], sft$fitIndices[pindex,5], labels=power,col="red")
    par(mfrow = c(1,1))
    dev.off()

    ## Correlation network construction and module identification (one-step)
    temp_cor <- cor
    cor <- WGCNA::cor

    sink(file = paste0(new_path,label,"_blockwiseModules.out"))
    netwk <- blockwiseModules(t_exprs,
                              corType = "pearson",
                              # == Adjacency Function ==
                              power = power,
                              networkType = "signed",
                              # == Tree and Block Options ==
                              deepSplit = 2,
                              pamRespectsDendro = F,
                              minModuleSize = minSize,
                              # == Module Adjustments ==
                              reassignThreshold = 0,
                              mergeCutHeight = MEDissThres,
                              # == Output Options
                              numericLabels = F,
                              verbose = 3)
    sink()

    cor <- temp_cor

    ## Plot module merge results
    MEList = moduleEigengenes(t_exprs, colors = netwk$unmergedColors)
    MEs = MEList$eigengenes
    colnames(MEs)=gsub("ME","",colnames(MEs))
    MEDiss = 1-cor(MEs);
    METree = hclust(as.dist(MEDiss), method = "average");

    pdf(file =  paste(new_path,"hclust_plots_",label,".pdf",sep=""),paper="a4")
    plot(METree, main = paste0("Module eigengenes merge result (",
                               length(unique(netwk$unmergedColors))," -> ",
                               length(unique(netwk$colors)),")"),
         xlab = "", sub = "")
    abline(h=MEDissThres, col = "red")
    dev.off()

    # Plot the dendrogram and the module colors underneath
    for (i in 1:length(netwk$dendrograms)) {
      pdf(file =  paste(new_path,"dendogram_plots_",i,"_",label,".pdf",sep=""),paper="a4")
      plotDendroAndColors(netwk$dendrograms[[i]],
                          netwk$colors[netwk$blockGenes[[i]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05 )
      dev.off()
    }

    # Get Module Eigengenes per cluster
    MEs0 = orderMEs(netwk$MEs)
    module_order = gsub("ME","",names(MEs0))
    MEs0$treatment = sample_conditions(rownames(t_exprs))

    mME = MEs0 %>%
      pivot_longer(-treatment) %>%
      dplyr::mutate(name = gsub("ME", "", name), name = factor(name, levels = module_order))

    mt.plot <-ggplot(mME, aes(x=treatment, y=name, fill=value)) +
      geom_tile() +
      theme_bw() +
      scale_fill_gradient2(
        low = "blue",
        high = "red",
        mid = "white",
        midpoint = 0,
        limit = c(-1,1)) +
      theme(axis.text.x = element_text(angle=90, size = 8, vjust = 0.5, hjust = 1)) +
      labs(title = "Module-type Relationships", y = "Modules", x="Condition", fill="corr")
    
    pdf(file =  paste(new_path,"WGCNA_module_type_",label,".pdf",sep=""),paper="a4")
    plot(mt.plot)
    dev.off()

    png(paste(new_path,"WGCNA_module_type_",label,".png",sep=""))
    plot(mt.plot)
    dev.off()

    ## Module-trait significance (p-values)
    tryCatch({
        nSamples <- nrow(t_exprs)
        MEs_ordered <- orderMEs(netwk$MEs)
        conditions <- sample_conditions(rownames(t_exprs))
        unique_conds <- unique(conditions)
        trait_mat <- matrix(0, nrow = nSamples, ncol = length(unique_conds))
        colnames(trait_mat) <- unique_conds
        for (cond in unique_conds) { trait_mat[conditions == cond, cond] <- 1 }
        moduleTraitCor <- cor(MEs_ordered, trait_mat, use = "p")
        moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
        write.table(data.frame(Module = rownames(moduleTraitCor), moduleTraitCor, check.names = FALSE),
                    file = paste0(new_path, "WGCNA_module_trait_cor_", label, ".txt"),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        write.table(data.frame(Module = rownames(moduleTraitPvalue), moduleTraitPvalue, check.names = FALSE),
                    file = paste0(new_path, "WGCNA_module_trait_pvalue_", label, ".txt"),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }, error = function(e) {
        cat(paste0("Module-trait significance could not be computed: ", e$message, "\n"))
    })

    # Save module assignments
    module_df = data.frame(
      gene_id = names(netwk$colors),
      clusters = netwk$colors)
    colnames(module_df)[1] <- "Gene_ID"
    write.table(module_df[order(module_df$clusters),], file = paste(new_path,"WGCNA_all_modules_",label,".txt",sep=""),sep="\t",row.names = F,col.names=T,quote=F)
    for (cluster in unique(module_df$clusters)){
      module_df_tmp <- module_df[module_df$clusters==cluster,]
      write.table(module_df_tmp[,1], file = paste(new_path,"WGCNA_modules_",label,"_cluster",cluster,"_Gene_IDs.txt",sep=""),sep="\n",row.names = F, col.names=F,quote=F)
    }

    ## Hub gene identification
    tryCatch({
        adj <- adjacency(t_exprs, power = power, type = "signed")
        connectivity <- intramodularConnectivity(adj, netwk$colors)
        connectivity$Gene_ID <- names(netwk$colors)
        connectivity$module <- netwk$colors
        connectivity <- connectivity[order(-connectivity$kWithin), ]
        write.table(connectivity,
                    file = paste0(new_path, "WGCNA_hub_genes_connectivity_", label, ".txt"),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        hub_genes <- do.call(rbind, lapply(unique(netwk$colors), function(mod) {
            mod_conn <- connectivity[connectivity$module == mod, ]
            head(mod_conn[order(-mod_conn$kWithin), c("Gene_ID", "module", "kTotal", "kWithin", "kOut", "kDiff")], 10)
        }))
        write.table(hub_genes,
                    file = paste0(new_path, "WGCNA_top10_hub_genes_", label, ".txt"),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        cat(paste0("Hub gene analysis complete: ", nrow(hub_genes), " top hub genes across ",
                   length(unique(netwk$colors)), " modules\n"))
    }, error = function(e) {
        cat(paste0("Hub gene identification could not be computed: ", e$message, "\n"))
    })

    ## Module functional enrichment (compareCluster GOBP)
    if (!is.null(orgdb)) {
        tryCatch({
            module_genes <- data.frame(
                gene_id = names(netwk$colors),
                module = netwk$colors,
                stringsAsFactors = FALSE
            )
            # Exclude grey (unassigned) module
            module_genes <- module_genes[module_genes$module != "grey", ]
            
            if (nrow(module_genes) > 0) {
                # Map to Entrez IDs
                entrez_map <- map_to_entrez(module_genes$gene_id, orgdb)
                if (!is.null(entrez_map) && nrow(entrez_map) > 0) {
                    module_mapped <- merge(module_genes, entrez_map, by.x = "gene_id", by.y = "input_id")
                    module_mapped <- module_mapped[!is.na(module_mapped$ENTREZID), ]
                    
                    if (nrow(module_mapped) > 0 && length(unique(module_mapped$module)) > 1) {
                        cat("  Running comparative GOBP enrichment on WGCNA modules...\n")
                        
                        # Universe: all expressed genes mapped to Entrez
                        universe_ids <- if (!is.null(all_gene_ids)) all_gene_ids else names(netwk$colors)
                        universe_map <- map_to_entrez(universe_ids, orgdb)
                        universe_entrez <- if (!is.null(universe_map)) unique(universe_map$ENTREZID) else NULL
                        
                        comp_go <- compareCluster(
                            ENTREZID ~ module,
                            data = module_mapped,
                            fun = "enrichGO",
                            OrgDb = orgdb,
                            ont = "BP",
                            universe = universe_entrez,
                            pvalueCutoff = 1
                        )
                        
                        if (!is.null(comp_go) && nrow(as.data.frame(comp_go)) > 0) {
                            res_df <- as.data.frame(comp_go)
                            write.table(res_df, file = paste0(new_path, "WGCNA_module_enrichment_GOBP_", label, ".txt"),
                                        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
                            
                            sig_df <- res_df[res_df$p.adjust < 0.05, ]
                            if (nrow(sig_df) > 0) {
                                comp_go_sig <- comp_go
                                comp_go_sig@compareClusterResult <- sig_df
                                p_comp <- enrichplot::dotplot(comp_go_sig, showCategory = 5,
                                            title = paste("WGCNA Module GOBP Enrichment -", label))
                                ggsave(paste0(new_path, "WGCNA_module_enrichment_GOBP_", label, ".png"),
                                       p_comp, width = 12, height = 8)
                                cat(paste0("  Module enrichment: ", nrow(sig_df), " significant GOBP terms (dotplot saved)\n"))
                            } else {
                                cat("  Module enrichment: table saved, no significant terms (p.adjust < 0.05)\n")
                            }
                        } else {
                            cat("  No GOBP terms found for WGCNA modules\n")
                        }
                    } else {
                        cat("  Skipping module enrichment: not enough non-grey modules with mapped genes\n")
                    }
                } else {
                    cat("  Skipping module enrichment: could not map module genes to Entrez IDs\n")
                }
            } else {
                cat("  Skipping module enrichment: all genes assigned to grey module\n")
            }
        }, error = function(e) {
            cat(paste0("  Module enrichment error: ", e$message, "\n"))
        })
    }

    cat(paste0("WGCNA complete: ", length(unique(netwk$colors)), " modules detected from ",
               nrow(nexpr_mat), " genes\n\n"))
    return(invisible(netwk))
}


################################################################################
#### Helper: run STRINGdb (with capping + chunking)
################################################################################
run_stringdb <- function(deg_df, label, new_path, organism_taxid, orgdb = NULL) {
    tryCatch({
        cat(paste0("\n  STRINGdb [", label, "]...\n"))

        if (nrow(deg_df) == 0) {
            cat("    Skipping STRINGdb: no genes\n")
            return(invisible(NULL))
        }

        system(paste("mkdir -p", new_path))

        if(!("Gene_ID" %in% colnames(deg_df))){
            cat("WARNING: 'Gene_ID' column not found. Trying to amend by renaming column 1\n")
            colnames(deg_df)[1] <- "Gene_ID"
        }

        # Map to gene symbols via OrgDb if available (for better STRING mapping)
        map_col <- "Gene_ID"
        if (!is.null(orgdb)) {
            sym_map <- tryCatch({
                ids_clean <- gsub("\\..*", "", toupper(deg_df$Gene_ID))
                res <- suppressMessages(AnnotationDbi::select(orgdb,
                                           keys = ids_clean,
                                           columns = "SYMBOL",
                                           keytype = "ENSEMBL"))
                res <- res[!is.na(res$SYMBOL), ]
                res <- res[!duplicated(res$ENSEMBL), ]
                res
            }, error = function(e) NULL)
            
            if (!is.null(sym_map) && nrow(sym_map) > 0) {
                deg_df$ENSEMBL_clean <- gsub("\\..*", "", toupper(deg_df$Gene_ID))
                deg_df <- merge(deg_df, sym_map, by.x = "ENSEMBL_clean", by.y = "ENSEMBL", all.x = TRUE)
                deg_df$SYMBOL[is.na(deg_df$SYMBOL)] <- deg_df$Gene_ID[is.na(deg_df$SYMBOL)]
                map_col <- "SYMBOL"
            }
        }

        string_db <- STRINGdb$new(version="12.0", species=organism_taxid,
                                  score_threshold=400, network_type="full", input_directory="")
        example1_mapped <- string_db$map(deg_df, map_col, removeUnmappedRows = TRUE)

        if (nrow(example1_mapped) == 0) {
            cat("    No genes mapped to STRING database\n")
            return(invisible(NULL))
        }

        hits <- example1_mapped$STRING_id

        # Network plot (capped at 400 genes for readability)
        plot_ids <- hits
        if (length(plot_ids) > 400) {
            cat(paste0("    Capping STRINGdb network plot to top 400 most significant genes (total mapped: ", length(plot_ids), ")\n"))
            plot_ids <- plot_ids[1:400]
        }

        tryCatch({
            pdf(file = paste(new_path, "String_network_", label, ".pdf", sep=""), paper="a4")
            string_db$plot_network(plot_ids)
            dev.off()
        }, error = function(e) {
            cat(paste0("    Network plot error: ", e$message, "\n"))
            if (dev.cur() > 1) dev.off()
        })

        # Colored network plot (logFC)
        col_fc <- grep("logFC", colnames(example1_mapped), val=T)
        if (length(col_fc) > 0) {
            tryCatch({
                example1_mapped_pval05 <- string_db$add_diff_exp_color(example1_mapped, logFcColStr=col_fc)
                payload_id <- string_db$post_payload(example1_mapped_pval05$STRING_id,
                                                     colors=example1_mapped_pval05$color)
                pdf(file = paste(new_path, "String_network_colored_", label, ".pdf", sep=""), paper="a4")
                string_db$plot_network(plot_ids, payload_id=payload_id)
                dev.off()
            }, error = function(e) {
                cat(paste0("    Colored network plot error: ", e$message, "\n"))
                if (dev.cur() > 1) dev.off()
            })
        }

        # Functional enrichment (chunked in blocks of 2000)
        chunk_size <- 2000
        num_chunks <- ceiling(length(hits) / chunk_size)
        enrichment_list <- list()
        for (i in 1:num_chunks) {
            start_idx <- (i - 1) * chunk_size + 1
            end_idx <- min(i * chunk_size, length(hits))
            chunk_ids <- hits[start_idx:end_idx]
            cat(sprintf("    Querying STRINGdb enrichment chunk %d/%d (%d genes)...\n", i, num_chunks, length(chunk_ids)))
            tryCatch({
                enrich_chunk <- string_db$get_enrichment(chunk_ids)
                if (!is.null(enrich_chunk) && nrow(enrich_chunk) > 0) {
                    colnames(enrich_chunk)[colnames(enrich_chunk) == "inputGenes"] <- "geneID"
                    enrichment_list[[i]] <- enrich_chunk
                }
            }, error = function(e) {
                cat(sprintf("    Enrichment chunk %d error: %s\n", i, e$message))
            })
        }
        if (length(enrichment_list) > 0) {
            enrichment_combined <- do.call(rbind, enrichment_list)
            enrichment_combined <- enrichment_combined[!duplicated(enrichment_combined), ]
            write.table(enrichment_combined,
                        file = paste0(new_path, "STRINGdb_functional_enrichment_", label, ".txt"),
                        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }

        # Interactions table (chunked in blocks of 2000)
        interactions_list <- list()
        for (i in 1:num_chunks) {
            start_idx <- (i - 1) * chunk_size + 1
            end_idx <- min(i * chunk_size, length(hits))
            chunk_ids <- hits[start_idx:end_idx]
            cat(sprintf("    Querying STRINGdb interactions chunk %d/%d (%d genes)...\n", i, num_chunks, length(chunk_ids)))
            tryCatch({
                interact_chunk <- string_db$get_interactions(chunk_ids)
                if (!is.null(interact_chunk) && nrow(interact_chunk) > 0) {
                    interactions_list[[i]] <- interact_chunk
                }
            }, error = function(e) {
                cat(sprintf("    Interactions chunk %d error: %s\n", i, e$message))
            })
        }
        if (length(interactions_list) > 0) {
            interactions_combined <- do.call(rbind, interactions_list)
            interactions_combined <- interactions_combined[!duplicated(interactions_combined), ]
            write.table(interactions_combined,
                        file = paste0(new_path, "STRINGdb_interactions_", label, ".txt"),
                        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
            cat(paste0("    Saved STRINGdb pairwise interactions table (", nrow(interactions_combined), " interactions)\n"))
        }

        # PPI enrichment test
        tryCatch({
            ppi_enrich <- string_db$get_ppi_enrichment(hits)
            write.table(ppi_enrich, file = paste0(new_path, "STRINGdb_PPI_enrichment_test_", label, ".txt"),
                        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
            cat(paste0("    PPI enrichment p-value: ", ppi_enrich$p_value, "\n"))
        }, error = function(e) {
            cat(paste0("    PPI enrichment test error: ", e$message, "\n"))
        })

        # Clusters
        tryCatch({
            clustersList <- string_db$get_clusters(example1_mapped$STRING_id)
            for(i in seq(1:length(clustersList))){
                pdf(file = paste(new_path, "String_network_clusters_", label, "_", i, ".pdf", sep=""), paper="a4")
                string_db$plot_network(clustersList[[i]])
                dev.off()
            }
            clustersList_2 <- lapply(clustersList, function(x) {
                example1_mapped$Gene_ID[example1_mapped$STRING_id %in% x]
            })
            clustersList_3 <- do.call(rbind, lapply(1:length(clustersList_2), function(i) {
                data.frame(Gene_ID = clustersList_2[[i]], clusters = i)
            }))
            write.table(clustersList_3, file = paste(new_path, "STRINGdb_all_modules_", label, ".txt", sep=""),
                        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
            for (cluster in unique(clustersList_3$clusters)){
                module_df_tmp <- clustersList_3[clustersList_3$clusters==cluster,]
                write.table(module_df_tmp[,1],
                            file = paste(new_path, "STRINGdb_modules_", label, "_cluster", cluster, "_Gene_IDs.txt", sep=""),
                            sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
            }
        }, error = function(e) {
            cat(paste0("    Cluster analysis error: ", e$message, "\n"))
        })

        cat(paste0("  STRINGdb complete for: ", label, "\n"))
    }, error = function(e) {
        cat(paste0("STRINGdb analysis failed (organism taxonid=", organism_taxid,
                   " may not be available in STRINGdb): ", e$message, "\n"))
    })
}


################################################################################
#### WGCNA on all expressed genes (canonical mode, run once)
################################################################################
if (wgcna_mode == "all") {
    cat("Running WGCNA on all expressed genes (canonical mode)...\n\n")
    setwd(path)
    a <- get_expression()
    nexpr_mat <- as.matrix(a[, -grep("Gene_ID", colnames(a))])
    rownames(nexpr_mat) <- a$Gene_ID
    new_path <- paste0(out_base, "WGCNA/")
    run_wgcna(nexpr_mat, new_path, label = "all_expressed_genes",
              orgdb = orgdb, all_gene_ids = a$Gene_ID)
}


################################################################################
#### Per-comparison analyses loop
################################################################################
for (genes_interest in list.files(pattern = pattern_search, path = path)){

#### WGCNA on DEGs only (alternative mode, per comparison)
    if (wgcna_mode == "degs") {
		cat(paste0("Processing WGCNA of DEGs from ",genes_interest,"...\n\n"))
		setwd(path)
		a <- get_expression()
		b <- as.data.frame(data.table::fread(genes_interest))
		if("FDR" %in% colnames(b)){b <- b[b$FDR < 0.05,]}
		a <- a[a$Gene_ID %in% b$Gene_ID,]
		nexpr_mat <- as.matrix(a[,-grep("Gene_ID",colnames(a))])
		rownames(nexpr_mat) <- a$Gene_ID
		label=sub("\\..*","",basename(genes_interest))
		new_path=paste0(out_base,"WGCNA/")
		run_wgcna(nexpr_mat, new_path, label = label,
		          orgdb = orgdb, all_gene_ids = get_expression()$Gene_ID)
    }

  
#### Computing STRINGdb (with direction-specific subsets):
    cat(paste0("\n\nProcessing STRINGdb from ",genes_interest,"...\n\n"))
    setwd(path)
    b <- as.data.frame(data.table::fread(genes_interest))

    # Identify FDR and logFC columns
    fdr_col <- NULL
    fc_col <- NULL
    if ("FDR" %in% colnames(b)) fdr_col <- "FDR"
    fc_candidates <- grep("logFC|log2Ratio|log2FoldChange", colnames(b), value = TRUE, ignore.case = TRUE)
    if (length(fc_candidates) > 0) fc_col <- fc_candidates[1]

    # Filter by FDR
    if (!is.null(fdr_col)) b <- b[b[[fdr_col]] < 0.05, ]

    label_base <- sub("\\..*", "", basename(genes_interest))

    # Prepare standardized columns for GSEA
    deg_for_gsea <- b
    if ("Gene_ID" %in% colnames(deg_for_gsea)) {
        if (!is.null(fc_col)) deg_for_gsea$logFC <- deg_for_gsea[[fc_col]]
        if (!is.null(fdr_col)) deg_for_gsea$FDR <- deg_for_gsea[[fdr_col]]
    }

    # Run GSEA (on all DEGs, uses full ranked list)
    if (!is.null(orgdb) && "Gene_ID" %in% colnames(deg_for_gsea) && !is.null(fc_col)) {
        # For GSEA, use the FULL unfiltered DEG table (not FDR-filtered) to get a proper ranking
        setwd(path)
        b_full <- as.data.frame(data.table::fread(genes_interest))
        b_full$logFC <- b_full[[fc_col]]
        if (!is.null(fdr_col)) b_full$FDR <- b_full[[fdr_col]] else b_full$FDR <- 1
        gsea_dir <- paste0(out_base, "GSEA/")
        system(paste("mkdir -p", gsea_dir))
        run_gsea(b_full, label_base, gsea_dir, orgdb, kegg_org, species_label)
    }

    # Direction-specific STRINGdb analyses
    for (direction in c("all", "up", "down")) {
        b_subset <- b
        dir_label <- paste0(label_base, "_", toupper(direction))

        if (direction == "up" && !is.null(fc_col)) {
            b_subset <- b[b[[fc_col]] > 0, ]
        } else if (direction == "down" && !is.null(fc_col)) {
            b_subset <- b[b[[fc_col]] < 0, ]
        }
        # If no logFC column, only run "all" direction
        if (direction != "all" && is.null(fc_col)) next

        if (nrow(b_subset) == 0) {
            cat(paste0("  No DEGs for direction: ", direction, ", skipping\n"))
            next
        }

        new_path <- paste0(out_base, "STRINGdb/", direction, "/")
        run_stringdb(b_subset, dir_label, new_path, organism, orgdb)
    }


#### Computing BONOBO:
    tryCatch({
		if (nzchar(Sys.which("netzoopy"))) {
			cat(paste0("Processing BONOBO of ",genes_interest,"...\n\n"))
			setwd(path)
			a <- get_expression()
			b <- as.data.frame(data.table::fread(genes_interest))
			if("FDR" %in% colnames(b)){b <- b[b$FDR < 0.05,]}
			a <- a[a$Gene_ID %in% b$Gene_ID,]
			nexpr_mat <- as.matrix(a[,-grep("Gene_ID",colnames(a))])
			rownames(nexpr_mat) <- a$Gene_ID
			bonobo_out <- paste0(out_base, "BONOBO/results_", tools::file_path_sans_ext(genes_interest), "/")
			dir.create(bonobo_out, recursive = TRUE, showWarnings = FALSE)
			bonobo_expr <- paste0(bonobo_out, "expression_vst_degs_no_header.txt")
			write.table(nexpr_mat, file = bonobo_expr,sep="\t",row.names = T, col.names=F,quote=F)
			system(paste0("netzoopy bonobo --expression_file ",bonobo_expr," --output_folder '", bonobo_out, "' --output_format '.csv' --sparsify --save_pvals"))
		} else {
			cat("netzoopy is not installed or not in PATH. Skipping BONOBO analysis.\n")
		}
    }, error = function(e) {
        cat(paste0("BONOBO analysis failed: ", e$message, "\n"))
    })
}

					  
save.image(file=paste0(out_base,"network_analyses_envir.RData"))
