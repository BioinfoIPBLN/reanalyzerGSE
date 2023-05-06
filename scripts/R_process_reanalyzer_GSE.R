#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
genes <- args[2]
filter_option <- args[3]
organism <- args[4]
enrichment_databases <- args[5]
targets_file <- args[6]


###### Load read counts, format, filter, start differential expression analyses, get RPKM, save...:
cat("\nProcessing counts and getting figures...\n")
a <- lapply(paste0(path,"/miARma_out/str_readcount_results/",list.files(paste0(path,"/miARma_out/str_readcount_results"),pattern = ".tab$")),
            function(x){data.table::fread(x)[,c(1,7)]})
b <- data.table::fread(paste0(path,"/miARma_out/Readcount_results/str-Size.tab"))
colnames(b)[1] <- "Geneid"
b$Length[is.na(b$Length)] <- 0
gene_counts <- Reduce(merge,a)
colnames(gene_counts)[2:dim(gene_counts)[2]] <- gsub(".*//star_results/|_nat_str_Aligned.sortedByCoord.out.bam","",colnames(gene_counts)[2:dim(gene_counts)[2]])
gene_counts <- as.data.frame(merge(gene_counts,b))
rownames(gene_counts) <- gene_counts$Geneid
gene_counts <- gene_counts[,-1]
gene_counts$Gene_ID <- stringr::str_to_title(rownames(gene_counts))
colnames(gene_counts) <- basename(colnames(gene_counts))
# Reorder so gene_counts columns follow the order of GSMXXXXXXX, or alfanumeric if GSM not present in the colnames:
if (length(grep("_GSM",colnames(gene_counts)))==0){
  gene_counts <- gene_counts[,gtools::mixedorder(colnames(gene_counts),decreasing=T)]
} else {
  idx <- colnames(gene_counts)[unlist(lapply(strsplit(colnames(gene_counts),"_|__"),function(x){any(startsWith(x,"GSM"))}))]
  gene_counts <- gene_counts[,c(idx[order(unlist(lapply(strsplit(idx,"_|__"),function(x){grep("^GSM",x,val=T)})))],"Gene_ID","Length")]
}
dir.create(paste0(path,"/final_results_reanalysis"),showWarnings=F)
write.table(gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)),grep("Gene_ID|Length",colnames(gene_counts),invert=T))],
            file=paste0(path,"/final_results_reanalysis/Raw_counts_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")

###### Batch effect correction if requested:
if (file.exists(paste0(path,"/GEO_info/batch_vector.txt")) && !file.exists(paste0(path,"/GEO_info/batch_biological_variables.txt"))){
  count_matrix <- as.matrix(gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)])
  batch <- data.table::fread(paste0(path,"/GEO_info/batch_vector.txt"),head=F,sep="*")$V1
  batch <- unlist(strsplit(batch,","))
  adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.numeric(batch), group=NULL)
  write.table(batch,file=paste0(path,"/GEO_info/batch_vector.txt"),quote = F,row.names = F, col.names = F,sep = "\n"); print("BATCH DONE_1")
  write.table(adjusted_counts,
              file=paste0(path,"/final_results_reanalysis/counts_adjusted_ComBat_seq.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
} else if (file.exists(paste0(path,"/GEO_info/batch_vector.txt")) && file.exists(paste0(path,"/GEO_info/batch_biological_variables.txt"))){
  count_matrix <- as.matrix(gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)])
  batch <- data.table::fread(paste0(path,"/GEO_info/batch_vector.txt"),head=F,sep="*")$V1
  batch <- unlist(strsplit(batch,","))
  biological_cov <- data.table::fread(paste0(path,"/GEO_info/batch_biological_variables.txt"),head=F,sep="*")$V1
  if (stringr::str_count(biological_cov," ") == 0){
    group <- unlist(strsplit(biological_cov,","))
    adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.numeric(batch), group=as.numeric(group))
    write.table(batch,file=paste0(path,"/GEO_info/batch_vector.txt"),quote = F,row.names = F, col.names = F,sep = "\n")
    write.table(group,file=paste0(path,"/GEO_info/batch_biological_variables.txt"),quote = F,row.names = F, col.names = F,sep = "\n"); print("BATCH DONE_2")
    write.table(adjusted_counts,
              file=paste0(path,"/final_results_reanalysis/counts_adjusted_ComBat_seq.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  } else {
    covar_mat <- matrix(rep(1,stringr::str_count(gsub(" .*","",biological_cov),",")+1))
    for (i in 1:(stringr::str_count(biological_cov," ") + 1)){
      covar_mat <- cbind(covar_mat,as.numeric(unlist(strsplit(strsplit(biological_cov," ")[[1]][i],","))))
    }
    covar_mat <- covar_mat[,-1]
    adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.numeric(batch), group=NULL, covar_mod=covar_mat)
    write.table(batch,file=paste0(path,"/GEO_info/batch_vector.txt"),quote = F,row.names = F, col.names = F,sep = "\n")
    write.table(covar_mat,file=paste0(path,"/GEO_info/batch_biological_variables.txt"),quote = F,row.names = F, col.names = F,sep = "\t"); print("BATCH DONE_3")
    write.table(adjusted_counts,
              file=paste0(path,"/final_results_reanalysis/Raw_counts_adjusted_ComBat_seq.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  }

}

pheno <- data.table::fread(paste0(path,"/GEO_info/samples_info.txt"),head=F)[,2:3]
colnames(pheno) <- c("sample","condition")
pheno <- pheno[gtools::mixedorder(pheno$sample,decreasing=T),]

layout <- read.table(list.files(pattern = "library_layout_info.txt", recursive = TRUE, full.names=T, path=path))$V1
# Make sure that in the case of paired studies I'm taking the correct groups:
if (file.exists(paste0(path,"/GEO_info/library_layout_info.txt")) & (dim(pheno)[1] > dim(gene_counts)[2]-2) & layout=="PAIRED"){
  pheno <- unique(pheno)
  print("Readjusting pheno in this paired study...")
}
# Also fix the designs if required:
for (z in list.files(pattern = "design_possible_full", recursive = TRUE, full.names=T, path=paste0(path,"/GEO_info"))){
  if (file.exists(paste0(path,"/GEO_info/library_layout_info.txt")) & (as.numeric(system(paste0("cat ",z," | wc -l"),intern=T)) > dim(gene_counts)[2]-2) & layout=="PAIRED"){
      tmp_condition <- read.table(z,head=F,blank.lines.skip=FALSE)$V1
      x=1:length(tmp_condition)
      tmp_condition <- tmp_condition[x%%2==1]
      write.table(tmp_condition,
            file=z,quote = F,row.names = F, col.names = F,sep = "\n")
      print("Fixing designs in this paired study...")
  }
}

###### Get edgeR object and normalized counts:
cat("Please double check the following lists are in the same order (automatically extracted and ordered column names of the counts vs rows of the pheno/targets data):\n")
print(colnames(gene_counts)[grep("Gene_ID|Length",colnames(gene_counts),invert=T)]); print(pheno)
edgeR_object <- edgeR::DGEList(counts=gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)],
                 group=pheno$condition,
                 genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])

###### Filter counts:
filter <- function(filter="standard",data,min_group=3){
  if(filter == "standard"){
    keep <- rowSums(edgeR::cpm(data)>1) >= min_group
    data <- data[keep,]
    data$samples$lib.size  <- colSums(data$counts)
  }
  else if(filter == "bin"){
    keep <- rowSums(data$counts)>0
    data <- data[keep,]
    data$samples$lib.size  <- colSums(data$counts)
  }
  else{
    stop("At the moment only bin/standard are supported")
  }
  return(data)
}
cat(paste0("Chosen filter is ",filter_option,"\n"))
cat(paste0("Number of genes before filter: ", nrow(edgeR_object)),"\n")
edgeR_object_prefilter <- edgeR_object
edgeR_object <- filter(filter=filter_option,edgeR_object) # Make sure of use bin to capture Cort and the lower expressed genes
cat(paste0("Number of genes after filter: ", nrow(edgeR_object)),"\n")
edgeR_object_norm <- edgeR::calcNormFactors(edgeR_object)
edgeR_object_norm <- edgeR::estimateCommonDisp(edgeR_object_norm, robust=TRUE)
if (is.na(edgeR_object_norm$common.dispersion)){
  edgeR_object_norm$common.dispersion <- 0.4 ^ 2
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
  cat("\nEstimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate...\n")
}

edgeR_object_norm <- edgeR::estimateTagwiseDisp(edgeR_object_norm)
gene_counts_rpkm <- as.data.frame(edgeR::rpkm(edgeR_object_norm,normalized.lib.sizes=TRUE))
colnames(gene_counts_rpkm) <- rownames(edgeR_object_norm$samples)
gene_counts_rpkm$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm))

###### Reorder so gene_counts columns follow the alfanumeric order, easier for the users in the written tables, although for the figures the script needs GSMXXXX-ordered:
gene_counts_rpkm_to_write <- gene_counts_rpkm[,c(grep("Gene_ID",colnames(gene_counts_rpkm)),grep("Gene_ID",colnames(gene_counts_rpkm),invert=T))]
gene_counts_rpkm_to_write <- gene_counts_rpkm_to_write[,c("Gene_ID",sort(colnames(gene_counts_rpkm_to_write)[-1]))]
write.table(gene_counts_rpkm_to_write,
            file=paste0(path,"/final_results_reanalysis/RPKM_counts_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
gene_counts_rpkm_log <- log2(gene_counts_rpkm_to_write[,grep("Gene_ID",colnames(gene_counts_rpkm_to_write),invert=T)] + 0.1)
gene_counts_rpkm_log <- cbind(gene_counts_rpkm_to_write$Gene_ID,gene_counts_rpkm_log)
colnames(gene_counts_rpkm_log)[1] <- "Gene_ID"
write.table(gene_counts_rpkm_log,
            file=paste0(path,"/final_results_reanalysis/RPKM_counts_genes_log2_0.1.txt"),quote = F,row.names = F, col.names = T,sep = "\t")

###### Similar to above, obtain RPKM counts but from the ComBat-Seq-adjusted counts instead of the raw counts
if (exists("adjusted_counts")){
  edgeR_object_combat <- edgeR::DGEList(counts=adjusted_counts,
                   group=pheno$condition,
             genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])
  cat(paste0("Number of genes combat before filter: ", nrow(edgeR_object_combat)),"\n")
  edgeR_object_prefilter_combat <- edgeR_object_combat
  edgeR_object_combat <- filter(filter=filter_option,edgeR_object_combat) # Make sure of use bin to capture Cort and the lower expressed genes
  cat(paste0("Number of genes combat after filter: ", nrow(edgeR_object_combat)),"\n")
  edgeR_object_norm_combat <- edgeR::calcNormFactors(edgeR_object_combat)
  edgeR_object_norm_combat <- edgeR::estimateCommonDisp(edgeR_object_norm_combat, robust=TRUE)
  if (is.na(edgeR_object_norm_combat$common.dispersion)){
    edgeR_object_norm_combat$common.dispersion <- 0.4 ^ 2
    # https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
    cat("Estimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate...")
  }
  edgeR_object_norm_combat <- edgeR::estimateTagwiseDisp(edgeR_object_norm_combat)
  gene_counts_rpkm_combat <- as.data.frame(edgeR::rpkm(edgeR_object_norm_combat,normalized.lib.sizes=TRUE))
  colnames(gene_counts_rpkm_combat) <- rownames(edgeR_object_norm_combat$samples)
  gene_counts_rpkm_combat$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm_combat))

  gene_counts_rpkm_combat_to_write <- gene_counts_rpkm_combat[,c(grep("Gene_ID",colnames(gene_counts_rpkm_combat)),grep("Gene_ID",colnames(gene_counts_rpkm_combat),invert=T))]
  gene_counts_rpkm_combat_to_write <- gene_counts_rpkm_combat_to_write[,c("Gene_ID",sort(colnames(gene_counts_rpkm_combat_to_write)[-1]))]
  write.table(gene_counts_rpkm_combat_to_write,
        file=paste0(path,"/final_results_reanalysis/RPKM_counts_ComBat_seq_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  gene_counts_rpkm_combat_log <- log2(gene_counts_rpkm_combat_to_write[,grep("Gene_ID",colnames(gene_counts_rpkm_combat_to_write),invert=T)] + 0.1)
  gene_counts_rpkm_combat_log <- cbind(gene_counts_rpkm_combat_to_write$Gene_ID,gene_counts_rpkm_combat_log)
  colnames(gene_counts_rpkm_combat_log)[1] <- "Gene_ID"
  write.table(gene_counts_rpkm_combat_log,
        file=paste0(path,"/final_results_reanalysis/RPKM_counts_ComBat_seq_genes_log2_0.1.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
}


###### Figure of the expr of certain genes of interest:
## Introduce in the violin plot statistics, loop through the different designs to get different coloring and grouping... etc
dir.create(paste0(path,"/final_results_reanalysis"),showWarnings = FALSE)
dir.create(paste0(path,"/final_results_reanalysis/violin"),showWarnings = FALSE)
dir.create(paste0(path,"/final_results_reanalysis/QC_and_others"),showWarnings = FALSE)
dir.create(paste0(path,"/final_results_reanalysis/violin_stats"),showWarnings = FALSE)
write.table(c("PLEASE NOTE:","If any of the violin plot do not display the expected statistics based on the title, it's because the test wasn't possible due to the data distribution or similar reasons, and therefore it's not the most suitable. If it is, then looking at it manually is required"),
      file=paste0(path,"/final_results_reanalysis/violin_stats/readme.txt"),quote = F,row.names = F, col.names = F,sep = "\n")
if (!exists("adjusted_counts")){
  gene_counts_rpkm_to_plot <- gene_counts_rpkm
} else {
  gene_counts_rpkm_to_plot <- gene_counts_rpkm_combat
}

### If the user has requested to analyze only some of the samples/column based on the GSMXXXXX, filter and restrict here.
for (i in unlist(strsplit(genes,","))){
  for (z in list.files(pattern = "design_possible_full", recursive = TRUE, full.names=T, path=paste0(path,"/GEO_info"))){
    if (length(unique(read.table(z,head=F,blank.lines.skip=T)$V1)) == 1){
      next
    }
    system(paste0("sed 's/^$/-/g' -i ",z)); print(i); print(z)
    if (i %in% gene_counts_rpkm_to_plot$Gene_ID){
      a <- gene_counts_rpkm_to_plot[gene_counts_rpkm_to_plot$Gene_ID==i,-grep("Gene_ID",colnames(gene_counts_rpkm_to_plot))]
      df <- data.frame(sample=names(a),
                       expr_RPKM=as.numeric(unname(a)),
                       condition=read.table(z,head=F,blank.lines.skip=FALSE)$V1)
    if (file.exists(paste0(path,"/GEO_info/gsm_manual_filter.txt"))){
        gsm_manual_filter <- data.table::fread(paste0(path,"/GEO_info/gsm_manual_filter.txt"),head=F,sep="*")$V1
        idxs_gsm_manual <- which(unlist(lapply(strsplit(colnames(gene_counts_rpkm),"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))
        gene_counts_rpkm_to_plot <- gene_counts_rpkm_to_plot[,
                                                      c(idxs_gsm_manual,which(colnames(gene_counts_rpkm)=="Gene_ID"))]
        a <- gene_counts_rpkm_to_plot[gene_counts_rpkm_to_plot$Gene_ID==i,-grep("Gene_ID",colnames(gene_counts_rpkm_to_plot))]                      
        df <- data.frame(sample=names(a),
                  expr_RPKM=as.numeric(unname(a)),
                  condition=read.table(z,head=F,blank.lines.skip=FALSE)$V1[idxs_gsm_manual])
      if (length(unique(read.table(z,head=F,blank.lines.skip=FALSE)$V1[idxs_gsm_manual])) == 1){
        next
      }
    }
      df$"Expr_RPKM_log2_01" <- round(log2(df$expr_RPKM + 0.1),2)
      lab_title <- tryCatch({
                tmp_lab_title <- data.table::fread(list.files(pattern = "series_matrix.txt.gz$", recursive = TRUE, full.names=T, path=path)[1],fill=T)
                lab_title <- paste(tmp_lab_title$V2[grep("Series_geo_accession",tmp_lab_title$V1)],tmp_lab_title$V2[grep("Series_title",tmp_lab_title$V1)],sep=": ")
               }, error=function(e){
                lab_title <- read.table(list.files(pattern = "study_title.txt$", recursive = TRUE, full.names=T, path=path))$V1
               })
      print(paste0("Title is ", lab_title))
      df$sample <- gsub("_t|m_Rep|_seq|_KO|_WT","",df$sample)
      suppressMessages(library(ggpubr,quiet = T,warn.conflicts = F))
      p <- ggbarplot(df, x = "sample", y = "Expr_RPKM_log2_01", color = "condition",
        add = "mean_se", label=T,lab.vjust = 4,
        position = ggplot2::position_dodge()) +
        theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",color="Condition", title=paste0(lab_title,"// GENE SHOWN: ", i))
      ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin/",i,"_barplot_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)

      df2 <- reshape::melt(df)
      df2 <- unique(df2[df2$variable=="Expr_RPKM_log2_01",c("condition","value","sample")])
      suppressMessages(library(dplyr,quiet = T,warn.conflicts = F))
      df2 <- df2 %>%
        add_count(condition, name = "condition_n")
        df2$sample2 <- as.numeric(1:dim(df2)[1])
      list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(df2$condition), unique(df2$condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T) # Used ***** only to make sure is a separator that's not going to be used in the conditions name
      df2$i <- 1
      labs_total <- aggregate(i~condition,df2,sum)
      labs_total$value <- NA
      labs_total$i <- paste0("Total points: ",labs_total$i)
      p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      geom_text(data=labs_total,aes(x=condition,y= max(df2$value) + 0.1,label=i)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0(lab_title,"// GENE SHOWN: ", i))
      ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin/",i,"_violin_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
      write.table(paste0("Samples_numbering:\n",paste0("Number_",1:length(df2$sample),": ",df2$sample,collapse="\n")),
                  file=paste0(path,"/final_results_reanalysis/label_samples.txt"),quote = F,row.names = F, col.names = F,sep = "\n")

      p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      # geom_signif(comparisons = list(unique(df2$condition)),
      #             map_signif_level=T, color="black",linetype="dashed",text=2.5) +
      stat_compare_means(method = "t.test", comparisons = list_combinations) + # Add pairwise comparisons p-value
      # stat_compare_means(method = "t.test", label.y = 10) + # Add global p-value
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0(lab_title,"// GENE SHOWN: ", i))
      ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_ttest_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
      p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      # geom_signif(comparisons = list(unique(df2$condition)),
      #             map_signif_level=T, color="black",linetype="dashed",text=2.5) +
      stat_compare_means(method = "wilcox.test",comparisons = list_combinations) + # Add pairwise comparisons p-value
      # stat_compare_means(method = "wilcox.test",label.y = 10) + # Add global p-value
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0(lab_title,"// GENE SHOWN: ", i))
      ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_wilcoxtest_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
      p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      # geom_signif(comparisons = list(unique(df2$condition)),
      #             map_signif_level=T, color="black",linetype="dashed",text=2.5) +
      stat_compare_means(method = "kruskal.test",comparisons = list_combinations) + # Add pairwise comparisons p-value
      # stat_compare_means(method = "kruskal.test",label.y = 10) + # Add global p-value
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0(lab_title,"// GENE SHOWN: ", i))
      ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_kruskaltest_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
      p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      # geom_signif(comparisons = list(unique(df2$condition)),
      #             map_signif_level=T, color="black",linetype="dashed",text=2.5) +
      stat_compare_means(method = "anova",comparisons = list_combinations) + # Add pairwise comparisons p-value
      # stat_compare_means(method = "anova",label.y = 10) + # Add global p-value
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0(lab_title,"// GENE SHOWN: ", i))
       ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_anova_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
cat("/n/n")
}
}
}

print("PDF barplot and violin plots done. Genes highlighted are:")
print(genes)


###### Attempt of Differential Gene Expression Analyses... modified from Bioinfo Unit to use here:
suppressMessages(library(edgeR,quiet = T,warn.conflicts = F))
suppressMessages(library(ggplot2,quiet = T,warn.conflicts = F))
suppressMessages(library(ggrepel,quiet = T,warn.conflicts = F))
Volcano<-function(data,file,main){
  options(ggrepel.max.overlaps = Inf)
  data<-data$table
  data$Gene<-rownames(data)
  data$threshold<-"nDEG"
  data[data$FDR <= 0.05,"threshold"]<-"DEG"
  #
  #data[data$FDR==0,"adj.P.Val"]<-1e-318
  #
  data$FDR[data$FDR==0] <- 1e-318

  ##Los 10 con mejor FC
  real_DE<-data[data$FDR<=0.05,]
  selected_FC<-head(real_DE[order(abs(real_DE$logFC),decreasing = T),], 20)
  max_value=max(abs(real_DE$logFC))
  p <- ggplot(data=data, aes(x=logFC, y=-log10(FDR), colour=threshold)) + coord_cartesian(xlim = c(-max_value, max_value )) +
    geom_point(alpha=0.4, size=1.75) +
    xlab("log2 fold change") + ylab("-log10 adj.P.Val") +
    geom_text_repel(data=selected_FC, aes(label=Gene),colour="black",size=3) + ggtitle(main) + theme(plot.title = element_text(hjust = 0.5)) #adding text for the top1 20 genes  
  ggsave(p, filename = file,width=30, height=30)
}
DGE <- function(comp,organism="mouse", myFDR=0.05, myFC=0){
  et <- exactTest(edgeR_object_norm,pair = comp)
  #Extracting the statistical data order by p-value
  top1 <- topTags(et, n=nrow(et), adjust.method="BH",sort.by="PValue")

  summary(decideTests(et))
  nrow(top1$table[top1$table$FDR<=myFDR & abs(top1$table$logFC)>= myFC,])
  print(comp); print("Top results:")
  print(head(top1$table,10)[,c(-1,-2)])
  myLabel1=paste(comp, collapse = '_vs_')
  myLabel2=paste(myLabel1, "FDR", myFDR, "FC", myFC, sep="_")
  # createExcel(top1$table,paste("DEG_",myLabel1,".xlsx", sep=""),organism = organism)
  # my_enrichment(top1$table,FA_label=myLabel2,cutoff=0.05,organism = organism, FDR=myFDR,FC=myFC)
  Volcano(top1,paste(path,"/final_results_reanalysis/DGE/Volcano_plot_",myLabel1,".pdf", sep=""),myLabel1)
  return(top1)
}

print("Attempting differential gene expression analyses between the conditions:")
dir.create(paste0(path,"/final_results_reanalysis/DGE"),showWarnings=F)
for (z in list.files(pattern = "design_possible_full", recursive = TRUE, full.names=T, path=paste0(path,"/GEO_info"))){
  if (length(unique(read.table(z,head=F,blank.lines.skip=T)$V1)) == 1){
     next
  }
  if (exists("gsm_manual_filter")){
    if (length(unique(read.table(z,head=F,blank.lines.skip=FALSE)$V1[idxs_gsm_manual])) == 1){
      next
    } else {
      condition <- read.table(z,head=F,blank.lines.skip=FALSE)$V1[idxs_gsm_manual]
    }
  } else {
    condition <- read.table(z,head=F,blank.lines.skip=FALSE)$V1
  }
  list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(condition), unique(condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T)
  list_combinations <- lapply(list_combinations,function(x){if (sum(!startsWith(x,"__"))==2){paste0("__",x)} else {x}})
  print(list_combinations)
  existing <- length(list.files(path=paste0(path,"/final_results_reanalysis/DGE/"),pattern=".RData"))
for (i in 1:length(list_combinations)){
  # print(i)
  if(sum(!startsWith(condition,"__")) == length(condition)){condition <- paste0("__",condition)}
  if (exists("gsm_manual_filter")){
    idxs_gsm_manual_2 <- which(unlist(lapply(strsplit(colnames(gene_counts),"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))
    idxs_gsm_manual_2_2 <- which(unlist(lapply(strsplit(pheno$sample,"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))

    edgeR_object <- edgeR::DGEList(counts=gene_counts[,idxs_gsm_manual_2],
                     group=pheno[idxs_gsm_manual_2_2,"condition"]$condition,
                     genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])
    edgeR_object_prefilter <- edgeR_object
    edgeR_object <- filter(filter=filter_option,edgeR_object) # Make sure of use bin to capture Cort and the lower expressed genes
    edgeR_object_norm <- edgeR::calcNormFactors(edgeR_object)
    edgeR_object_norm <- edgeR::estimateCommonDisp(edgeR_object_norm, robust=TRUE)
    if (is.na(edgeR_object_norm$common.dispersion)){
      edgeR_object_norm$common.dispersion <- 0.4 ^ 2
      # https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
      cat("Estimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate...")
    }
    edgeR_object_norm <- edgeR::estimateTagwiseDisp(edgeR_object_norm)
    gene_counts_rpkm <- as.data.frame(edgeR::rpkm(edgeR_object_norm,normalized.lib.sizes=TRUE))
    colnames(gene_counts_rpkm) <- rownames(edgeR_object_norm$samples)
    gene_counts_rpkm$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm))
  } else {
    edgeR_object_norm$samples$group <- as.factor(condition)
  }
  write.table(paste0("\nComparison number ",i+existing,": ",list_combinations[[i]]),
          file=paste0(path,"/final_results_reanalysis/DGE/list_comp.txt"),quote = F,row.names = F, col.names = F,sep = "\n", append=T)
  edgeR_results <- DGE(comp=list_combinations[[i]])
  colnames(edgeR_results$table)[3] <- paste0(colnames(edgeR_results$table)[3],list_combinations[[i]][1],"_VS_",list_combinations[[i]][2])
  save.image(file=paste0(path,"/final_results_reanalysis/DGE/DGE_analysis_comp",i+existing,".RData"))
  write.table(edgeR_results$table,
          file=paste0(path,"/final_results_reanalysis/DGE/DGE_analysis_comp",i+existing,".txt"),quote = F,row.names = F, col.names = T,sep = "\t")
}
}
print("Done")


###### Attempt of automatic functional enrichment analyses...:
print("Attempting automatic gene ontology enrichment analyses of the DGE results...")
suppressMessages(library(autoGO,quiet = T,warn.conflicts = F))
# choose_database()
### 03/2023:
  # [1] "Achilles_fitness_decrease"                          "Achilles_fitness_increase"                         
  # [3] "Aging_Perturbations_from_GEO_down"                  "Aging_Perturbations_from_GEO_up"                   
  # [5] "Allen_Brain_Atlas_10x_scRNA_2021"                   "Allen_Brain_Atlas_down"                            
  # [7] "Allen_Brain_Atlas_up"                               "ARCHS4_Cell-lines"                                 
  # [9] "ARCHS4_IDG_Coexp"                                   "ARCHS4_Kinases_Coexp"                              
  # [11] "ARCHS4_TFs_Coexp"                                   "ARCHS4_Tissues"                                    
  # [13] "Azimuth_Cell_Types_2021"                            "BioCarta_2013"                                     
  # [15] "BioCarta_2015"                                      "BioCarta_2016"                                     
  # [17] "BioPlanet_2019"                                     "BioPlex_2017"                                      
  # [19] "Cancer_Cell_Line_Encyclopedia"                      "CCLE_Proteomics_2020"                              
  # [21] "CellMarker_Augmented_2021"                          "ChEA_2013"                                         
  # [23] "ChEA_2015"                                          "ChEA_2016"                                         
  # [25] "ChEA_2022"                                          "Chromosome_Location"                               
  # [27] "Chromosome_Location_hg19"                           "ClinVar_2019"                                      
  # [29] "CORUM"                                              "COVID-19_Related_Gene_Sets"                        
  # [31] "COVID-19_Related_Gene_Sets_2021"                    "Data_Acquisition_Method_Most_Popular_Genes"        
  # [33] "dbGaP"                                              "DepMap_WG_CRISPR_Screens_Broad_CellLines_2019"     
  # [35] "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019"     "Descartes_Cell_Types_and_Tissue_2021"              
  # [37] "Diabetes_Perturbations_GEO_2022"                    "Disease_Perturbations_from_GEO_down"               
  # [39] "Disease_Perturbations_from_GEO_up"                  "Disease_Signatures_from_GEO_down_2014"             
  # [41] "Disease_Signatures_from_GEO_up_2014"                "DisGeNET"                                          
  # [43] "Drug_Perturbations_from_GEO_2014"                   "Drug_Perturbations_from_GEO_down"                  
  # [45] "Drug_Perturbations_from_GEO_up"                     "DrugMatrix"                                        
  # [47] "DSigDB"                                             "Elsevier_Pathway_Collection"                       
  # [49] "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"          "ENCODE_Histone_Modifications_2013"                 
  # [51] "ENCODE_Histone_Modifications_2015"                  "ENCODE_TF_ChIP-seq_2014"                           
  # [53] "ENCODE_TF_ChIP-seq_2015"                            "Enrichr_Libraries_Most_Popular_Genes"              
  # [55] "Enrichr_Submissions_TF-Gene_Coocurrence"            "Enrichr_Users_Contributed_Lists_2020"              
  # [57] "Epigenomics_Roadmap_HM_ChIP-seq"                    "ESCAPE"                                            
  # [59] "FANTOM6_lncRNA_KD_DEGs"                             "Gene_Perturbations_from_GEO_down"                  
  # [61] "Gene_Perturbations_from_GEO_up"                     "Genes_Associated_with_NIH_Grants"                  
  # [63] "GeneSigDB"                                          "Genome_Browser_PWMs"                               
  # [65] "GlyGen_Glycosylated_Proteins_2022"                  "GO_Biological_Process_2013"                        
  # [67] "GO_Biological_Process_2015"                         "GO_Biological_Process_2017"                        
  # [69] "GO_Biological_Process_2017b"                        "GO_Biological_Process_2018"                        
  # [71] "GO_Biological_Process_2021"                         "GO_Cellular_Component_2013"                        
  # [73] "GO_Cellular_Component_2015"                         "GO_Cellular_Component_2017"                        
  # [75] "GO_Cellular_Component_2017b"                        "GO_Cellular_Component_2018"                        
  # [77] "GO_Cellular_Component_2021"                         "GO_Molecular_Function_2013"                        
  # [79] "GO_Molecular_Function_2015"                         "GO_Molecular_Function_2017"                        
  # [81] "GO_Molecular_Function_2017b"                        "GO_Molecular_Function_2018"                        
  # [83] "GO_Molecular_Function_2021"                         "GTEx_Aging_Signatures_2021"                        
  # [85] "GTEx_Tissue_Expression_Down"                        "GTEx_Tissue_Expression_Up"                         
  # [87] "GWAS_Catalog_2019"                                  "HDSigDB_Human_2021"                                
  # [89] "HDSigDB_Mouse_2021"                                 "HMDB_Metabolites"                                  
  # [91] "HMS_LINCS_KinomeScan"                               "HomoloGene"                                        
  # [93] "HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression" "HuBMAP_ASCTplusB_augmented_2022"                   
  # [95] "Human_Gene_Atlas"                                   "Human_Phenotype_Ontology"                          
  # [97] "HumanCyc_2015"                                      "HumanCyc_2016"                                     
  # [99] "huMAP"                                              "IDG_Drug_Targets_2022"                             
  # [101] "InterPro_Domains_2019"                              "Jensen_COMPARTMENTS"                               
  # [103] "Jensen_DISEASES"                                    "Jensen_TISSUES"                                    
  # [105] "KEA_2013"                                           "KEA_2015"                                          
  # [107] "KEGG_2013"                                          "KEGG_2015"                                         
  # [109] "KEGG_2016"                                          "KEGG_2019_Human"                                   
  # [111] "KEGG_2019_Mouse"                                    "KEGG_2021_Human"                                   
  # [113] "Kinase_Perturbations_from_GEO_down"                 "Kinase_Perturbations_from_GEO_up"                  
  # [115] "KOMP2_Mouse_Phenotypes_2022"                        "L1000_Kinase_and_GPCR_Perturbations_down"          
  # [117] "L1000_Kinase_and_GPCR_Perturbations_up"             "Ligand_Perturbations_from_GEO_down"                
  # [119] "Ligand_Perturbations_from_GEO_up"                   "LINCS_L1000_Chem_Pert_Consensus_Sigs"              
  # [121] "LINCS_L1000_Chem_Pert_down"                         "LINCS_L1000_Chem_Pert_up"                          
  # [123] "LINCS_L1000_CRISPR_KO_Consensus_Sigs"               "LINCS_L1000_Ligand_Perturbations_down"             
  # [125] "LINCS_L1000_Ligand_Perturbations_up"                "lncHUB_lncRNA_Co-Expression"                       
  # [127] "MAGMA_Drugs_and_Diseases"                           "MCF7_Perturbations_from_GEO_down"                  
  # [129] "MCF7_Perturbations_from_GEO_up"                     "Metabolomics_Workbench_Metabolites_2022"           
  # [131] "MGI_Mammalian_Phenotype_2013"                       "MGI_Mammalian_Phenotype_2017"                      
  # [133] "MGI_Mammalian_Phenotype_Level_3"                    "MGI_Mammalian_Phenotype_Level_4"                   
  # [135] "MGI_Mammalian_Phenotype_Level_4_2019"               "MGI_Mammalian_Phenotype_Level_4_2021"              
  # [137] "Microbe_Perturbations_from_GEO_down"                "Microbe_Perturbations_from_GEO_up"                 
  # [139] "miRTarBase_2017"                                    "Mouse_Gene_Atlas"                                  
  # [141] "MSigDB_Computational"                               "MSigDB_Hallmark_2020"                              
  # [143] "MSigDB_Oncogenic_Signatures"                        "NCI-60_Cancer_Cell_Lines"                          
  # [145] "NCI-Nature_2015"                                    "NCI-Nature_2016"                                   
  # [147] "NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions"     "NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions"    
  # [149] "NIH_Funded_PIs_2017_Human_AutoRIF"                  "NIH_Funded_PIs_2017_Human_GeneRIF"                 
  # [151] "NURSA_Human_Endogenous_Complexome"                  "Old_CMAP_down"                                     
  # [153] "Old_CMAP_up"                                        "OMIM_Disease"                                      
  # [155] "OMIM_Expanded"                                      "Orphanet_Augmented_2021"                           
  # [157] "PanglaoDB_Augmented_2021"                           "Panther_2015"                                      
  # [159] "Panther_2016"                                       "Pfam_Domains_2019"                                 
  # [161] "Pfam_InterPro_Domains"                              "PFOCR_Pathways"                                    
  # [163] "PhenGenI_Association_2021"                          "PheWeb_2019"                                       
  # [165] "Phosphatase_Substrates_from_DEPOD"                  "PPI_Hub_Proteins"                                  
  # [167] "Proteomics_Drug_Atlas_2023"                         "ProteomicsDB_2020"                                 
  # [169] "Rare_Diseases_AutoRIF_ARCHS4_Predictions"           "Rare_Diseases_AutoRIF_Gene_Lists"                  
  # [171] "Rare_Diseases_GeneRIF_ARCHS4_Predictions"           "Rare_Diseases_GeneRIF_Gene_Lists"                  
  # [173] "Reactome_2013"                                      "Reactome_2015"                                     
  # [175] "Reactome_2016"                                      "Reactome_2022"                                     
  # [177] "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"  "RNAseq_Automatic_GEO_Signatures_Human_Down"        
  # [179] "RNAseq_Automatic_GEO_Signatures_Human_Up"           "RNAseq_Automatic_GEO_Signatures_Mouse_Down"        
  # [181] "RNAseq_Automatic_GEO_Signatures_Mouse_Up"           "Serine_Threonine_Kinome_Atlas_2023"                
  # [183] "SILAC_Phosphoproteomics"                            "SubCell_BarCode"                                   
  # [185] "SynGO_2022"                                         "SysMyo_Muscle_Gene_Sets"                           
  # [187] "Table_Mining_of_CRISPR_Studies"                     "Tabula_Muris"                                      
  # [189] "Tabula_Sapiens"                                     "TargetScan_microRNA"                               
  # [191] "TargetScan_microRNA_2017"                           "TF_Perturbations_Followed_by_Expression"           
  # [193] "TF-LOF_Expression_from_GEO"                         "TG_GATES_2020"                                     
  # [195] "Tissue_Protein_Expression_from_Human_Proteome_Map"  "Tissue_Protein_Expression_from_ProteomicsDB"       
  # [197] "Transcription_Factor_PPIs"                          "TRANSFAC_and_JASPAR_PWMs"                          
  # [199] "TRRUST_Transcription_Factors_2019"                  "UK_Biobank_GWAS_v1"                                
  # [201] "Virus_Perturbations_from_GEO_down"                  "Virus_Perturbations_from_GEO_up"                   
  # [203] "Virus-Host_PPI_P-HIPSTer_2020"                      "VirusMINT"                                         
  # [205] "WikiPathway_2021_Human"                             "WikiPathways_2013"                                 
  # [207] "WikiPathways_2015"                                  "WikiPathways_2016"                                 
  # [209] "WikiPathways_2019_Human"                            "WikiPathways_2019_Mouse"

if (!is.na(enrichment_databases)){
    if (length(enrichment_databases) > 0){
      enrichment_databases <- c("GO_Biological_Process_2021","GO_Molecular_Function_2021","GO_Cellular_Component_2021")}
}
enrichment_databases <- unlist(strsplit(enrichment_databases,","))
if (grepl("sapiens", organism, fixed=TRUE)){
  databases_autoGO <- unique(c(enrichment_databases,"WikiPathway_2021_Human","RNAseq_Automatic_GEO_Signatures_Human_Down","RNAseq_Automatic_GEO_Signatures_Human_Up","Reactome_2022","KEGG_2021_Human","HDSigDB_Human_2021"))
  databases_autoGO_print <- paste(databases_autoGO,collapse=",")
  print(paste0("The databases selected for autoGO are ",databases_autoGO_print,". Please double check autoGO::choose_database(), which has > 200 databases, in case you want to add any extra by using the pipeline argument..."))
} else if (grepl("musculus", organism, fixed=TRUE)){
  databases_autoGO <- unique(c(enrichment_databases,"WikiPathways_2019_Mouse","RNAseq_Automatic_GEO_Signatures_Mouse_Down","RNAseq_Automatic_GEO_Signatures_Mouse_Up","Reactome_2022","Mouse_Gene_Atlas","KEGG_2021_Mouse","HDSigDB_Mouse_2021"))
  databases_autoGO_print <- paste(databases_autoGO,collapse=",")
  print(paste0("The databases selected for autoGO are ",databases_autoGO_print,". Please double check autoGO::choose_database(), which has > 200 databases, in case you want to add any extra by using the pipeline argument..."))
} else {
  print(paste0("Your organism is ",organism,", and unfortunately the pipeline for automatic functional enrichment currently fully supports only mouse and human. We'll include non-model organisms eventually, but in the meantime please don't give up and double check if you can provide any extra database containing your organism from autoGO::choose_database(), which has > 200 databases. For now trying to run with the default databases..."))
  databases_autoGO <- enrichment_databases
  databases_autoGO_print <- paste(databases_autoGO,collapse=",")
  print(databases_autoGO_print)
}

setwd(paste0(path,"/final_results_reanalysis/DGE"))
for (f in list.files(pattern = "^DGE_analysis_comp.*.txt")){
  a <- data.table::fread(f,head=T,fill=T)
  b <- a[a$FDR < 0.05,1]
  d <- a[a$FDR < 0.05 & a$logFC>0,1]
  e <- a[a$FDR < 0.05 & a$logFC<0,1]
  if (dim(b)[1]!=0){
    write.table(b,file=paste0(gsub(".txt","",f),"_fdr_05.txt"),col.names = F,row.names = F,quote = F,sep="\t")  
  }
  if (dim(d)[1]!=0){
    write.table(d,file=paste0(gsub(".txt","",f),"_fdr_05_logpos.txt"),col.names = F,row.names = F,quote = F,sep="\t")
  }
  if (dim(e)[1]!=0){
    write.table(e,file=paste0(gsub(".txt","",f),"_fdr_05_logneg.txt"),col.names = F,row.names = F,quote = F,sep="\t")
  }  
}

setwd(paste0(path,"/final_results_reanalysis/DGE/"))
if (length(list.files(pattern="*_fdr_05.txt"))!=0){
  gene_lists_2 <- read_gene_lists(
    gene_lists_path = getwd(),
    which_list = c("everything"),
    from_autoGO = F,
    files_format = "_fdr_05.txt"
  )
  for (f in names(gene_lists_2)){
    print(paste0("Performing functional enrichment analyses for: ",basename(f)))
    setwd(paste0(path,"/final_results_reanalysis/DGE"))
    autoGO(gene_lists_2[f],databases_autoGO)  
    file.rename(paste0(getwd(),"/enrichment_tables/"), paste0(getwd(),"/",basename(f),"_funct_enrichment/"))
    setwd(paste0(getwd(),"/",basename(f),"_funct_enrichment/"))
    enrich_tables <- read_enrich_tables(
      enrich_table_path = getwd(),
      which_list = "everything",
      from_autoGO = F,
      files_format = ".tsv")
    for (i in 1:length(names(enrich_tables))){
      barplotGO(enrich_tables = enrich_tables[[i]],
                title = c(names(enrich_tables)[i],basename(f)),
                outfolder = getwd(),
                outfile = paste0(basename(f),"_",names(enrich_tables)[i],"_barplot.png"),
                from_autoGO = F)
      lolliGO(enrich_tables = enrich_tables[[i]],
              title = c(names(enrich_tables)[i],basename(f)),
              outfolder = getwd(),
              outfile = paste0(basename(f),"_",names(enrich_tables)[i],"_lolliplot.png"),
              from_autoGO = F)
    }
  }
} else {  
  print("No significant DEGs (FDR < 0.05) found")
}
setwd(paste0(path,"/final_results_reanalysis/DGE/"))
if (length(list.files(pattern="*_fdr_05_logneg.txt"))!=0){
  gene_lists_2 <- read_gene_lists(
    gene_lists_path = getwd(),
    which_list = c("everything"),
    from_autoGO = F,
    files_format = "_fdr_05_logneg.txt"
  )
  for (f in names(gene_lists_2)){
    print(paste0("Performing functional enrichment analyses for: ",basename(f)))
    setwd(paste0(path,"/final_results_reanalysis/DGE/"))
    autoGO(gene_lists_2[f],databases_autoGO)
    file.rename(paste0(getwd(),"/enrichment_tables/"), paste0(getwd(),"/",basename(f),"_funct_enrichment/"))
    setwd(paste0(getwd(),"/",basename(f),"_funct_enrichment/"))
    enrich_tables <- read_enrich_tables(
      enrich_table_path = getwd(),
      which_list = "everything",
      from_autoGO = F,
      files_format = ".tsv")
    for (i in 1:length(names(enrich_tables))){
      barplotGO(enrich_tables = enrich_tables[[i]],
                title = c(names(enrich_tables)[i],basename(f)),
                outfolder = getwd(),
                outfile = paste0(basename(f),"_",names(enrich_tables)[i],"_barplot.png"),
                from_autoGO = F)
      lolliGO(enrich_tables = enrich_tables[[i]],
              title = c(names(enrich_tables)[i],basename(f)),
              outfolder = getwd(),
              outfile = paste0(basename(f),"_",names(enrich_tables)[i],"_lolliplot.png"),
              from_autoGO = F)
    }
  }
} else {
  print("No significant DEGs (FDR < 0.05) with log2FC < 0 found")
}
setwd(paste0(path,"/final_results_reanalysis/DGE/"))
if (length(list.files(pattern="*_fdr_05_logpos.txt"))!=0){
  setwd(paste0(path,"/final_results_reanalysis/DGE/"))
  gene_lists_2 <- read_gene_lists(
    gene_lists_path = getwd(),
    which_list = c("everything"),
    from_autoGO = F,
    files_format = "_fdr_05_logpos.txt"
  )
  for (f in names(gene_lists_2)){
    print(paste0("Performing functional enrichment analyses for: ",basename(f)))
    setwd(paste0(path,"/final_results_reanalysis/DGE/"))
    autoGO(gene_lists_2[f],databases_autoGO)
    file.rename(paste0(getwd(),"/enrichment_tables/"), paste0(getwd(),"/",basename(f),"_funct_enrichment/"))
    setwd(paste0(getwd(),"/",basename(f),"_funct_enrichment/"))
    enrich_tables <- read_enrich_tables(
      enrich_table_path = getwd(),
      which_list = "everything",
      from_autoGO = F,
      files_format = ".tsv")
    for (i in 1:length(names(enrich_tables))){
      barplotGO(enrich_tables = enrich_tables[[i]],
                title = c(names(enrich_tables)[i],basename(f)),
                outfolder = getwd(),
                outfile = paste0(basename(f),"_",names(enrich_tables)[i],"_barplot.png"),
                from_autoGO = F)
      lolliGO(enrich_tables = enrich_tables[[i]],
              title = c(names(enrich_tables)[i],basename(f)),
              outfolder = getwd(),
              outfile = paste0(basename(f),"_",names(enrich_tables)[i],"_lolliplot.png"),
              from_autoGO = F)
    }
  }
} else {
  print("No significant DEGs (FDR < 0.05) with log2FC > 0 found")
}


###### QC PDF from Bioinfo and Laura:
cat("\nPerforming QC_PDF...\n")
suppressMessages(library("edgeR",quiet = T,warn.conflicts = F))
suppressMessages(library("RColorBrewer",quiet = T,warn.conflicts = F))
suppressMessages(library("genefilter",quiet = T,warn.conflicts = F))
suppressMessages(library("ggrepel",quiet = T,warn.conflicts = F))
suppressMessages(library("ggfortify",quiet = T,warn.conflicts = F))
suppressMessages(library("cluster",quiet = T,warn.conflicts = F))
suppressMessages(library("factoextra",quiet = T,warn.conflicts = F))
suppressMessages(library("ggplot2",quiet = T,warn.conflicts = F))
suppressMessages(library("M3C",quiet = T,warn.conflicts = F))
suppressMessages(library("AnnotationDbi",quiet = T,warn.conflicts = F))
# suppressMessages(library("xlsx",quiet = T,warn.conflicts = F))
suppressMessages(library("clusterProfiler",quiet = T,warn.conflicts = F))
suppressMessages(library("png",quiet = T,warn.conflicts = F))
suppressMessages(library("curl",quiet = T,warn.conflicts = F))
suppressMessages(library("corrplot",quiet = T,warn.conflicts = F))

label <- basename(path)

x_prefilter <- edgeR_object_prefilter
lcpm_prefilter <- cpm(x_prefilter, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...

x <- edgeR_object
cpm <- cpm(x) # This is normalized, altough not through edgeR, but the argument normalized.lib.sizes=TRUE by default in cpm...
lcpm <- cpm(x, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...
lcpm_no_log <- cpm(x, log=FALSE)

nsamples <- ncol(x)

if (!is.na(targets_file)){
    if (length(targets_file) > 0) {
        targets <- read.table(file = targets_file,header = F,stringsAsFactors=F)
    } else {
        targets <- read.table(file = paste0(path,"/GEO_info/samples_info.txt"),header = F,stringsAsFactors=F)
    }
} else {
    targets <- read.table(file = paste0(path,"/GEO_info/samples_info.txt"),header = F,stringsAsFactors=F)
}
colnames(targets) <- c("Filename","Name","Type")
targets$Filename <- gsub("_t|m_Rep|_seq|_KO|_WT","",paste(unique(targets$Filename,targets$Name),sep="_"))
targets <- targets[gtools::mixedorder(targets$Filename,decreasing=T),]
cat("Please double check the following is in the same order than the lists above (counts and pheno/targets):\n")
print(targets)

if (exists("gsm_manual_filter")){
    idxs_gsm_manual_3 <- which(unlist(lapply(strsplit(targets$Name,"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))    
    targets <- targets[idxs_gsm_manual_3,]
}

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
# c(L, M)
lcpm.cutoff <- log2(10/M + 2/L)

# Here the original x2 was normalized, so I used from above directly the edgeR object transformed the same way...
# x <- calcNormFactors(x)
# x$samples
# x <- estimateCommonDisp(x, robust=TRUE)
# x <- estimateTagwiseDisp(x)
# x2 <- x
x2 <- edgeR_object_norm
#x2$samples$norm.factors <- 1

lcpm2 <- cpm(x2, log=TRUE)
lcpm2_no_log <- cpm(x2, log=FALSE)

group <- x$samples$group
col.group <- as.factor(group)
# levels(col.group) <- RColorBrewer::brewer.pal(nlevels(col.group), "Set1")
# not enough colors sometimes, so I get random colors:
color = gplots::col2hex(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)])
levels(col.group) <- sample(color,nlevels(col.group))
col.group <- as.character(col.group)

# summary(lcpm)
# table(rowSums(x$counts==0)==6)


### The actual pdf file:
pdf(paste0(path,"/final_results_reanalysis/QC_and_others/",label,"_QC.pdf"),paper="A4")

### 1.1. Density rawcounts log2, cpm...:
col <- RColorBrewer::brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm_prefilter), col=col[1], lwd=2, las=2, main="", xlab="")
title(main="Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm_prefilter[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright",legend=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), text.col=col, bty = "n", cex = 0.5)
### 1.2. Density rawcounts log2, cpm...:
# x is the raw counts with the bin or standard filter:
plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
title(main="Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", legend=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), text.col=col, bty="n", cex = 0.5)
### 2.1. Unnormalised:
par(mfrow=c(1,2))
boxplot(lcpm, las=2, col=col.group, main="", names=targets$Type, cex.axis=0.4)
title(main="Unnormalized data",ylab="Log-cpm",xlab="sample_type")
### 2.2. Normalised:
boxplot(lcpm2, las=2, col=col.group, main="", names=targets$Type, cex.axis=0.4)
title(main="Normalized data",ylab="Log-cpm",xlab="sample_type")
### 3. Library size:
par(mfrow=c(1,1))
barplot(x$samples$lib.size,names.arg = gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name),las=2, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))))
### 4. Corrplot no log
tmp <- lcpm_no_log; colnames(tmp) <- gsub("_t|m_Rep|_seq|_KO|_WT","",colnames(tmp))
corrplot(cor(tmp,method="spearman"), method='number',type = 'upper')
corrplot(cor(tmp,method="spearman"), order='AOE')
### 5.1. MDS # Commented out because I've checked it's identical to the norm one
#z <- plotMDS(lcpm_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
#edge <- sd(z$x)
#plotMDS(lcpm_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
#title(main="MDS-PCoA Sample Names")
### 5.2. MDS_log # Commented out because I've checked it's almost identical to the norm one
# par(mfrow=c(1,1))
# z <- plotMDS(lcpm, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
# edge <- sd(z$x)
#cat(edge)
# plotMDS(lcpm, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
# title(main="MDS-PCoA log2 Sample Names")
### 5.3. MDS_norm
z <- plotMDS(lcpm2_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
edge <- sd(z$x)
plotMDS(lcpm2_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
title(main="MDS-PCoA Sample Names Norm")
### 5.4. MDS_log_norm
par(mfrow=c(1,1))
z <- plotMDS(lcpm2, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
edge <- sd(z$x)
#cat(edge)
plotMDS(lcpm2, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
title(main="MDS-PCoA log2 Sample Names Norm")
### 6. PCA por tipos
data_pca <- as.matrix(x)
data_pca <- as.data.frame(t(data_pca))
rownames(data_pca) <- targets$Filename
data_pca.PC = prcomp(data_pca)
data_pca$Type <- targets$Type
data_pca$Filename <- targets$Filename
data_pca$Name <- targets$Name
data_pca$Sex <- targets$Sex
data_pca$Age <- targets$Age
data_pca$VAS_Group <- targets$VAS_Group
data_pca$TypeII <- targets$TypeII
plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Type',xlim = c(-0.8,0.8),label.size=3,label.repel=T))
### 7. Heatmap 250 mots differential entities
rsd <- rowSds(as.matrix(x))
sel <- order(rsd, decreasing=TRUE)[1:250]
samplenames <- gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name)
heatmap(na.omit(as.matrix(x[sel,])),margins=c(10,8),main="Heatmap 250 most diff entities raw counts",cexRow=0.01,cexCol=0.5,labCol=samplenames)
### 8.1. Dendogram cluster raw  # Commented out because I've checked it's identical to the norm one
#par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
#pr.hc.c <- hclust(na.omit(dist(t(data))))
#plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of ", label, sep=""), labels=targets$Filemane, cex=0.5)
#Normalized clustering analysis plot
#pr.hc.c <- hclust(na.omit(dist(t(edgeR_object_norm$counts))))
#pr.hc.c <- hclust(na.omit(dist(t(cpm(x$counts,log=F)),method = "euclidean")))
#plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of raw counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)
### 8.2. Dendogram cluster raw_log # Commented out because I've checked it's identical to the norm one
#par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
#pr.hc.c <- hclust(na.omit(dist(t(cpm(x$counts,log=T)),method = "euclidean")))
#plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 raw counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)
### 8.3. Dendogram cluster raw norm
par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
pr.hc.c <- hclust(na.omit(dist(t(cpm(x2$counts,log=F)),method = "euclidean")))
plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of normalized counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)
### 8.3. Dendogram cluster raw norm
par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
pr.hc.c <- hclust(na.omit(dist(t(cpm(x2$counts,log=T)),method = "euclidean")))
plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 normalized counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)

#tSNE
#a <- tsne(x$counts,seed=100,labels=as.factor(targets$Type), perplex=perplex, legendtitle="Types",text=targets$Type ,dotsize=3, legendtextsize = 8) + ggtitle("Tsne") + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
#plot(a)
dev.off()

###### Add the figures using the counts corrected by ComBat-seq:
if (exists("adjusted_counts")){
  cat("\nQC_PDF ComBat-seq counts...\n")

  label <- basename(path)

  x_prefilter <- edgeR_object_prefilter_combat
  lcpm_prefilter <- cpm(x_prefilter, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...

  x <- edgeR_object_combat
  cpm <- cpm(x) # This is normalized, altough not through edgeR, but the argument normalized.lib.sizes=TRUE by default in cpm...
  lcpm <- cpm(x, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...
  lcpm_no_log <- cpm(x, log=FALSE)

  nsamples <- ncol(x)

  if (!is.na(targets_file)) {
    if (length(targets_file) > 0) {
        targets <- read.table(file = targets_file,header = F,stringsAsFactors=F)
    } else {
        targets <- read.table(file = paste0(path,"/GEO_info/samples_info.txt"),header = F,stringsAsFactors=F)
    }
  } else {
      targets <- read.table(file = paste0(path,"/GEO_info/samples_info.txt"),header = F,stringsAsFactors=F)
  }
  colnames(targets) <- c("Filename","Name","Type")
  targets$Filename <- gsub("_t|m_Rep|_seq|_KO|_WT","",paste(unique(targets$Filename,targets$Name),sep="_"))
  targets <- targets[gtools::mixedorder(targets$Filename,decreasing=T),]
  cat("Please double check the following is in the same order than the lists above (counts and pheno/targets):\n")
  print(targets)
  L <- mean(x$samples$lib.size) * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  # c(L, M)
  lcpm.cutoff <- log2(10/M + 2/L)

  # Here the original x2 was normalized, so I used from above directly the edgeR object transformed the same way...
  # x <- calcNormFactors(x)
  # x$samples
  # x <- estimateCommonDisp(x, robust=TRUE)
  # x <- estimateTagwiseDisp(x)
  # x2 <- x
  x2 <- edgeR_object_norm_combat
  #x2$samples$norm.factors <- 1

  lcpm2 <- cpm(x2, log=TRUE)
  lcpm2_no_log <- cpm(x2, log=FALSE)

  group <- x$samples$group
  col.group <- as.factor(group)
  # levels(col.group) <- RColorBrewer::brewer.pal(nlevels(col.group), "Set1")
  # not enough colors sometimes, so I get random colors:
  color = gplots::col2hex(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)])
  levels(col.group) <- sample(color,nlevels(col.group))
  col.group <- as.character(col.group)

  # summary(lcpm)
  # table(rowSums(x$counts==0)==6)


  ### The actual pdf file:
  pdf(paste0(path,"/final_results_reanalysis/QC_and_others/",label,"_QC_ComBat-seq.pdf"),paper="A4")

  ### 1.1. Density rawcounts log2, cpm...:
  col <- RColorBrewer::brewer.pal(nsamples, "Paired")
  par(mfrow=c(1,2))
  plot(density(lcpm_prefilter), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="Raw data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm_prefilter[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright",legend=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), text.col=col, bty = "n", cex = 0.5)
  ### 1.2. Density rawcounts log2, cpm...:
  # x is the raw counts with the bin or standard filter:
  plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="Filtered data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", legend=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), text.col=col, bty="n", cex = 0.5)
  ### 2.1. Unnormalised:
  par(mfrow=c(1,2))
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Type, cex.axis=0.4)
  title(main="Unnormalized data",ylab="Log-cpm",xlab="sample_type")
  ### 2.2. Normalised:
  boxplot(lcpm2, las=2, col=col.group, main="", names=targets$Type, cex.axis=0.4)
  title(main="Normalized data",ylab="Log-cpm",xlab="sample_type")
  ### 3. Library size:
  par(mfrow=c(1,1))
  barplot(x$samples$lib.size,names.arg = gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name),las=2, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))))
  ### 4. Corrplot no log
  tmp <- lcpm_no_log; colnames(tmp) <- gsub("_t|m_Rep|_seq|_KO|_WT","",colnames(tmp))
  corrplot(cor(tmp,method="spearman"), method='number',type = 'upper')
  corrplot(cor(tmp,method="spearman"), order='AOE')
  ### 5.1. MDS
  #z <- plotMDS(lcpm_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
  #edge <- sd(z$x)
  #plotMDS(lcpm_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  #title(main="MDS-PCoA Sample Names")
  ### 5.2. MDS_log
  #par(mfrow=c(1,1))
  #z <- plotMDS(lcpm, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
  #edge <- sd(z$x)
  #cat(edge)
  #plotMDS(lcpm, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  #title(main="MDS-PCoA log2 Sample Names")
  ### 5.3. MDS_norm
  z <- plotMDS(lcpm2_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
  edge <- sd(z$x)
  plotMDS(lcpm2_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="MDS-PCoA Sample Names Norm")
  ### 5.4. MDS_log_norm
  par(mfrow=c(1,1))
  z <- plotMDS(lcpm2, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
  edge <- sd(z$x)
  #cat(edge)
  plotMDS(lcpm2, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="MDS-PCoA log2 Sample Names Norm")
  ### 6. PCA por tipos
  data_pca <- as.matrix(x)
  data_pca <- as.data.frame(t(data_pca))
  rownames(data_pca) <- targets$Filename
  data_pca.PC = prcomp(data_pca)
  data_pca$Type <- targets$Type
  data_pca$Filename <- targets$Filename
  data_pca$Name <- targets$Name
  data_pca$Sex <- targets$Sex
  data_pca$Age <- targets$Age
  data_pca$VAS_Group <- targets$VAS_Group
  data_pca$TypeII <- targets$TypeII
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Type',xlim = c(-0.8,0.8),label.size=3,label.repel=T))
  ### 7. Heatmap 250 mots differential entities
  rsd <- rowSds(as.matrix(x))
  sel <- order(rsd, decreasing=TRUE)[1:250]
  samplenames <- gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name)
  heatmap(na.omit(as.matrix(x[sel,])),margins=c(10,8),main="Heatmap 250 most diff entities raw counts",cexRow=0.01,cexCol=0.5,labCol=samplenames)
  ### 8.1. Dendogram cluster raw
  #par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  #pr.hc.c <- hclust(na.omit(dist(t(cpm(x$counts,log=F)),method = "euclidean")))
  #plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of raw counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)
  ### 8.2. Dendogram cluster raw_log
  #par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  #pr.hc.c <- hclust(na.omit(dist(t(cpm(x$counts,log=T)),method = "euclidean")))
  #plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 raw counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)
  ### 8.3. Dendogram cluster raw norm
  par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  pr.hc.c <- hclust(na.omit(dist(t(cpm(x2$counts,log=F)),method = "euclidean")))
  plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of normalized counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)
  ### 8.3. Dendogram cluster raw norm
  par(mfrow=c(1,1), col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  pr.hc.c <- hclust(na.omit(dist(t(cpm(x2$counts,log=T)),method = "euclidean")))
  plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 normalized counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)

  #tSNE
  #a <- tsne(x$counts,seed=100,labels=as.factor(targets$Type), perplex=perplex, legendtitle="Types",text=targets$Type ,dotsize=3, legendtextsize = 8) + ggtitle("Tsne") + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  #plot(a)
  dev.off()
}


print("ALL DONE")
