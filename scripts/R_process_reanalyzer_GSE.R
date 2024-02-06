#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
input_dir <- args[2]
output_dir <- args[3]
genes <- args[4]
filter_option <- args[5]
organism <- args[6]
targets_file <- args[7]
diff_soft <- args[8]
covariab <- args[9] # if not provided, "none"
cdseq_exec <- args[10] # if not provided, "no"
restrict_comparisons <- args[11] # if not provided, "no"

###### Load read counts, format, filter, start differential expression analyses, get RPKM, save...:
  cat("\nProcessing counts and getting figures...\n")
  print(paste0(input_dir,". Current date: ",date()))
  a <- lapply(paste0(input_dir,"/str_readcount_results/",list.files(paste0(input_dir,"/str_readcount_results"),pattern = ".tab$")),
              function(x){data.table::fread(x)[,c(1,7)]})
  b <- data.table::fread(paste0(input_dir,"/Readcount_results/str-Size.tab"))
  colnames(b)[1] <- "Geneid"
  b$Length[is.na(b$Length)] <- 0
  gene_counts <- Reduce(merge,a)
  colnames(gene_counts)[2:dim(gene_counts)[2]] <- gsub(".*_results/|_nat.*","",colnames(gene_counts)[2:dim(gene_counts)[2]])
  gene_counts <- as.data.frame(merge(gene_counts,b))
  rownames(gene_counts) <- gene_counts$Geneid
  gene_counts <- gene_counts[,-1]
  gene_counts$Gene_ID <- stringr::str_to_title(rownames(gene_counts))
  colnames(gene_counts) <- basename(colnames(gene_counts))
  # Reorder so gene_counts columns follow the order of GSMXXXXXXX, or alfanumeric if GSM not present in the colnames:
  if (length(grep("_GSM[0-9]",colnames(gene_counts)))==0){
    gene_counts <- gene_counts[,c(order(grep("Length|Gene_ID",colnames(gene_counts),invert=T)),grep("Length|Gene_ID",colnames(gene_counts)))]
  } else {
    idx <- colnames(gene_counts)[unlist(lapply(strsplit(colnames(gene_counts),"_|__"),function(x){any(startsWith(x,"GSM"))}))]
    gene_counts <- gene_counts[,c(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"GSM"))}))],"_+"),function(x){grep("^GSM",x,val=T)})),decreasing=T)],"Gene_ID","Length")]
  }
  # This is crucial, also double checking reordering based on SRR* if present in the names:
  if (length(grep("_SRR[0-9]",colnames(gene_counts)))!=0){
    idx <- colnames(gene_counts)[unlist(lapply(strsplit(colnames(gene_counts),"_|__"),function(x){any(startsWith(x,"SRR"))}))]
    gene_counts <- gene_counts[,c(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"SRR"))}))],"_+"),function(x){grep("^SRR",x,val=T)})),decreasing=T)],"Gene_ID","Length")]
  }
  dir.create(paste0(output_dir,""),showWarnings=F)
  write.table(gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)),grep("Gene_ID|Length",colnames(gene_counts),invert=T))],
              file=paste0(output_dir,"/Raw_counts_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  # CDSeq preliminary deconvolution:
  if (cdseq_exec!="no"){
    suppressMessages(library(CDSeq,quiet = T,warn.conflicts = F))
    read_length <- data.table::fread(list.files(dirname(input_dir),pattern="*.ini",full.names=T)[1],fill=T)
    read_length <- as.numeric(gsub("read_length=","",grep("read_length",read_length$V1,val=T))) + 1 # in the miARma.ini the read_length is -1, so here +1 to restore
    gene_length_arg <- as.vector(gene_counts[,c("Length")]-read_length+1); gene_length_arg[gene_length_arg <= 1] <- 2
    result_cdseq <- CDSeq(bulk_data =  as.matrix(gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)]), 
                          cell_type_number = 2:30,
                          gene_length = gene_length_arg,
                          dilution_factor = 10,
                          mcmc_iterations = 1000,
                          verbose = T)
    save(result_cdseq,file=paste0(output_dir,"/result_cdseq.RData"))
    print(paste0("Number of cell types estimated by CDSeq: ",result_cdseq$estT))
    print("Estimated proportions for the number of cell types estimated by CDSeq: ")
    result_cdseq$estProp; cat("\n\n")
  }

###### Batch effect correction if requested:
  if (file.exists(paste0(path,"/GEO_info/batch_vector.txt")) && !file.exists(paste0(path,"/GEO_info/batch_biological_variables.txt"))){
    count_matrix <- as.matrix(gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)])
    batch <- data.table::fread(paste0(path,"/GEO_info/batch_vector.txt"),head=F,sep="*")$V1
    batch <- unlist(strsplit(batch,","))
    adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.numeric(batch), group=NULL)
    write.table(batch,file=paste0(path,"/GEO_info/batch_vector.txt"),quote = F,row.names = F, col.names = F,sep = "\n"); print("BATCH DONE_1")
    write.table(adjusted_counts,
                file=paste0(output_dir,"/counts_adjusted_ComBat_seq.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
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
                file=paste0(output_dir,"/counts_adjusted_ComBat_seq.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
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
                file=paste0(output_dir,"/Raw_counts_adjusted_ComBat_seq.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
    }

  }

  pheno <- as.data.frame(data.table::fread(paste0(path,"/GEO_info/samples_info.txt"),head=F))
  if(dim(pheno)[2]>2){
    pheno <- pheno[,2:3]
  }
  colnames(pheno) <- c("sample","condition")
  pheno <- pheno[gtools::mixedorder(pheno$sample,decreasing=T),]
  # Ensure reordering based on GSM* or SRR* if present in the names of pheno:
  if (length(grep("_GSM[0-9]",pheno$sample))!=0){
    idx <- pheno$sample[unlist(lapply(strsplit(pheno$sample,"_+"),function(x){any(startsWith(x,"GSM"))}))]
    pheno <- pheno[match(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"GSM"))}))],"_+"),function(x){grep("^GSM",x,val=T)})),decreasing=T)],pheno$sample),]
  }
  if (length(grep("_SRR[0-9]",pheno$sample))!=0){
    idx <- pheno$sample[unlist(lapply(strsplit(pheno$sample,"_+"),function(x){any(startsWith(x,"SRR"))}))]
    pheno <- pheno[match(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"SRR"))}))],"_+"),function(x){grep("^SRR",x,val=T)})),decreasing=T)],pheno$sample),]
  } else {
    # If there are not SRR nor GSM in the names (and then likely this has not been a download from database but local analyses with raw_data), ordering is also required... (not required mixed order and so, so the printed lists are in the same order)
    pheno <- pheno[order(pheno$sample),]
    # Should be redundant because done above, but make suer to order also gene_counts if this has not been a download from database but local analyses with raw_data (and then not SRR or GSM in sample name):
    gene_counts <- gene_counts[,c(order(colnames(gene_counts)[grep("Gene_ID|Length",colnames(gene_counts),invert=T)]),grep("Gene_ID|Length",colnames(gene_counts)))]
  }

  layout <- read.table(list.files(pattern = "library_layout_info.txt", recursive = TRUE, full.names=T, path=path))$V1
  # Make sure that in the case of paired studies I'm taking the correct groups:
  if (file.exists(paste0(path,"/GEO_info/library_layout_info.txt")) & (dim(pheno)[1] > dim(gene_counts)[2]-2) & layout=="PAIRED"){
    pheno$sample <- gsub("_SRR.*","",pheno$sample)
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
  suppressMessages(library(edgeR,quiet = T,warn.conflicts = F))
  cat("Please double check the following lists are in the same order (automatically extracted and ordered column names of the counts vs rows of the pheno/targets data):\n")  
  print(colnames(gene_counts)[grep("Gene_ID|Length",colnames(gene_counts),invert=T)]); print(pheno)
  edgeR_object <- DGEList(counts=gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)],
                   group=pheno$condition,
                   genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])

###### Filter counts:
  filter <- function(filter="standard",data,min_group=3){
    if(filter == "standard"){
      keep <- rowSums(cpm(data)>1) >= min_group
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
  edgeR_object_norm <- calcNormFactors(edgeR_object)

  if(covariab == "none"){
    edgeR_object_norm <- estimateCommonDisp(edgeR_object_norm, robust=TRUE)
    if (is.na(edgeR_object_norm$common.dispersion)){
      edgeR_object_norm$common.dispersion <- 0.4 ^ 2
      cat("\nEstimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate. Please consult the vignette for further information...\n")
    }
    edgeR_object_norm <- estimateTagwiseDisp(edgeR_object_norm)
  } else { 
    Treat <- edgeR_object$samples$group
    Time <- as.factor(unlist(strsplit(as.character(covariab),",")))
    design <- model.matrix(~0+Treat+Time)
    rownames(design) <- colnames(edgeR_object_norm)
    edgeR_object_norm <- estimateDisp(edgeR_object_norm, design, robust=TRUE) # Preparing for one covariable following edgeR vignette new methods Nov2023
  }

  gene_counts_rpkm <- as.data.frame(rpkm(edgeR_object_norm,normalized.lib.sizes=TRUE))
  colnames(gene_counts_rpkm) <- rownames(edgeR_object_norm$samples)
  gene_counts_rpkm$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm))

###### Reorder so gene_counts columns follow the alfanumeric order, easier for the users in the written tables, although for the figures the script needs GSMXXXX-ordered:
  gene_counts_rpkm_to_write <- gene_counts_rpkm[,c(grep("Gene_ID",colnames(gene_counts_rpkm)),grep("Gene_ID",colnames(gene_counts_rpkm),invert=T))]
  gene_counts_rpkm_to_write <- gene_counts_rpkm_to_write[,c("Gene_ID",sort(colnames(gene_counts_rpkm_to_write)[-1]))]
  # High/medium/low categ:    
  write.table(gene_counts_rpkm_to_write,
              file=paste0(output_dir,"/RPKM_counts_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  gene_counts_rpkm_to_write_categ <- gene_counts_rpkm_to_write
  for (col in colnames(gene_counts_rpkm_to_write_categ[,-1])){
    a <- Hmisc::cut2(gene_counts_rpkm_to_write_categ[,col],g=3); b <- as.character(a)
    b[b==levels(a)[1]] <- "Low"; b[b==levels(a)[2]] <- "Medium"; b[b==levels(a)[3]] <- "High"
    gene_counts_rpkm_to_write_categ[,paste0(col,"_categ")] <- b
    a <- Hmisc::cut2(gene_counts_rpkm_to_write_categ[,col],g=5); b <- as.character(a)
    b[b==levels(a)[1]] <- "Very_Low"; b[b==levels(a)[2]] <- "Low"; b[b==levels(a)[3]] <- "Medium"; b[b==levels(a)[4]] <- "High"; b[b==levels(a)[5]] <- "Very_High"
    gene_counts_rpkm_to_write_categ[,paste0(col,"_categ_2")] <- b
  }
  write.table(gene_counts_rpkm_to_write_categ,
              file=paste0(output_dir,"/RPKM_counts_genes_categ.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  gene_counts_rpkm_log <- log2(gene_counts_rpkm_to_write[,-1]+0.1)
  gene_counts_rpkm_log$Gene_ID <- gene_counts_rpkm_to_write$Gene_ID
  gene_counts_rpkm_log <- gene_counts_rpkm_log[,c(length(colnames(gene_counts_rpkm_log)),grep("Gene_ID",colnames(gene_counts_rpkm_log),invert=T))]
  write.table(gene_counts_rpkm_log,
              file=paste0(output_dir,"/RPKM_counts_genes_log2_0.1.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  gene_counts_rpkm_log_categ <- gene_counts_rpkm_log
  for (col in colnames(gene_counts_rpkm_log_categ[,-1])){
    a <- Hmisc::cut2(gene_counts_rpkm_log_categ[,col],g=3); b <- as.character(a)
    b[b==levels(a)[1]] <- "Low"; b[b==levels(a)[2]] <- "Medium"; b[b==levels(a)[3]] <- "High"
    gene_counts_rpkm_log_categ[,paste0(col,"_categ")] <- b
    a <- Hmisc::cut2(gene_counts_rpkm_log_categ[,col],g=5); b <- as.character(a)
    b[b==levels(a)[1]] <- "Very_Low"; b[b==levels(a)[2]] <- "Low"; b[b==levels(a)[3]] <- "Medium"; b[b==levels(a)[4]] <- "High"; b[b==levels(a)[5]] <- "Very_High"
    gene_counts_rpkm_log_categ[,paste0(col,"_categ_2")] <- b
  }
  write.table(gene_counts_rpkm_log_categ,
              file=paste0(output_dir,"/RPKM_counts_genes_log2_0.1_categ.txt"),quote = F,row.names = F, col.names = T,sep = "\t")

###### Similar to above, obtain RPKM counts but from the ComBat-Seq-adjusted counts instead of the raw counts
  if (exists("adjusted_counts")){
    edgeR_object_combat <- DGEList(counts=adjusted_counts,
                     group=pheno$condition,
               genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])
    cat(paste0("Number of genes combat before filter: ", nrow(edgeR_object_combat)),"\n")
    edgeR_object_prefilter_combat <- edgeR_object_combat
    edgeR_object_combat <- filter(filter=filter_option,edgeR_object_combat) # Make sure of use bin to capture Cort and the lower expressed genes
    cat(paste0("Number of genes combat after filter: ", nrow(edgeR_object_combat)),"\n")
    edgeR_object_norm_combat <- calcNormFactors(edgeR_object_combat)
    edgeR_object_norm_combat <- estimateCommonDisp(edgeR_object_norm_combat, robust=TRUE)
    if (is.na(edgeR_object_norm_combat$common.dispersion)){
      edgeR_object_norm_combat$common.dispersion <- 0.4 ^ 2
      # https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
      cat("Estimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate...")
    }
    edgeR_object_norm_combat <- estimateTagwiseDisp(edgeR_object_norm_combat)
    gene_counts_rpkm_combat <- as.data.frame(rpkm(edgeR_object_norm_combat,normalized.lib.sizes=TRUE))
    colnames(gene_counts_rpkm_combat) <- rownames(edgeR_object_norm_combat$samples)
    gene_counts_rpkm_combat$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm_combat))

    gene_counts_rpkm_combat_to_write <- gene_counts_rpkm_combat[,c(grep("Gene_ID",colnames(gene_counts_rpkm_combat)),grep("Gene_ID",colnames(gene_counts_rpkm_combat),invert=T))]
    gene_counts_rpkm_combat_to_write <- gene_counts_rpkm_combat_to_write[,c("Gene_ID",sort(colnames(gene_counts_rpkm_combat_to_write)[-1]))]
    write.table(gene_counts_rpkm_combat_to_write,
          file=paste0(output_dir,"/RPKM_counts_ComBat_seq_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
    gene_counts_rpkm_combat_log <- log2(gene_counts_rpkm_combat_to_write[,grep("Gene_ID",colnames(gene_counts_rpkm_combat_to_write),invert=T)] + 0.1)
    gene_counts_rpkm_combat_log <- cbind(gene_counts_rpkm_combat_to_write$Gene_ID,gene_counts_rpkm_combat_log)
    colnames(gene_counts_rpkm_combat_log)[1] <- "Gene_ID"
    write.table(gene_counts_rpkm_combat_log,
          file=paste0(output_dir,"/RPKM_counts_ComBat_seq_genes_log2_0.1.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  }
  print(paste0("Counts written. Current date: ",date()))

###### Figure of the expr of certain genes of interest:
  ## Introduce in the violin plot statistics, loop through the different designs to get different coloring and grouping... etc
  dir.create(paste0(output_dir,""),showWarnings = FALSE)
  dir.create(paste0(output_dir,"/QC_and_others"),showWarnings = FALSE)
  dir.create(paste0(output_dir,"/violin"),showWarnings = FALSE)
  write.table(c(paste0("PLEASE NOTE:","This folder contains barplots and violin plots for the expression of the of interest: ",genes,". In the filenames the statistical tests attempted are shown. If any of the violin plots do not display the expected statistics based on the title, it's because the test wasn't possible due to the data distribution or similar reasons, and therefore it's not the most suitable")),
        file=paste0(output_dir,"/violin/readme.txt"),quote = F,row.names = F, col.names = F,sep = "\n")
  write.table(c("PLEASE VISIT https://bioinfoipbln.shinyapps.io/expressionvisualizationapp/ to use the Expression Visualization App and obtain the same figures for any other gene of interest uploading the expression tables"),
        file=paste0(output_dir,"/violin/readme_external_ExpressionVisualizationApp.txt"),quote = F,row.names = F, col.names = F,sep = "\n")
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
        a <- a[,gtools::mixedsort(colnames(a))]
        df <- data.frame(sample=names(a),
                         expr_RPKM=as.numeric(unname(a)),
                         condition=gtools::mixedsort(read.table(z,head=F,blank.lines.skip=FALSE)$V1))
      if (file.exists(paste0(path,"/GEO_info/gsm_manual_filter.txt"))){
          gsm_manual_filter <- data.table::fread(paste0(path,"/GEO_info/gsm_manual_filter.txt"),head=F,sep="*")$V1
          idxs_gsm_manual <- which(unlist(lapply(strsplit(gtools:mixedsort(colnames(gene_counts_rpkm)),"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))
          gene_counts_rpkm_to_plot <- gene_counts_rpkm[,gtools:mixedsort(colnames(gene_counts_rpkm))]
          gene_counts_rpkm_to_plot <- gene_counts_rpkm_to_plot[,
                                                        c(idxs_gsm_manual,which(colnames(gene_counts_rpkm_to_plot)=="Gene_ID"))]
          a <- gene_counts_rpkm_to_plot[gene_counts_rpkm_to_plot$Gene_ID==i,-grep("Gene_ID",colnames(gene_counts_rpkm_to_plot))]
          a <- a[,gtools::mixedsort(colnames(a))]
          df <- data.frame(sample=names(a),
                    expr_RPKM=as.numeric(unname(a)),
                    condition=gtools::mixedsort(read.table(z,head=F,blank.lines.skip=FALSE)$V1)[idxs_gsm_manual])
        if (length(unique(read.table(z,head=F,blank.lines.skip=FALSE)$V1[idxs_gsm_manual])) == 1){
          next
        }
      }
      print("Please double check that the order betwen the columns is correct and association between conditions and samples in the rows is correct:"); print(df)
        df$"Expr_RPKM_log2_01" <- round(log2(df$expr_RPKM + 0.1),2)
        lab_title <- tryCatch({
                  tmp_lab_title <- data.table::fread(list.files(pattern = "series_matrix.txt.gz$", recursive = TRUE, full.names=T, path=path)[1],fill=T)
                  lab_title <- paste(tmp_lab_title$V2[grep("Series_geo_accession",tmp_lab_title$V1)],tmp_lab_title$V2[grep("Series_title",tmp_lab_title$V1)],sep=": ")
                 }, error=function(e){
                  lab_title <- read.table(list.files(pattern = "study_title.txt$", recursive = TRUE, full.names=T, path=path))$V1
                 })
        print(paste0("Title of the project is ", lab_title))
        df$sample <- gsub("_t|m_Rep|_seq|_KO|_WT","",df$sample)
        suppressMessages(library(ggpubr,quiet = T,warn.conflicts = F))
        suppressMessages(library(plotly,quiet = T,warn.conflicts = F))
        suppressMessages(library(dplyr,quiet = T,warn.conflicts = F))
        p <- ggbarplot(df, x = "sample", y = "Expr_RPKM_log2_01", color = "condition",
          add = "mean_se", label=T,lab.vjust = 4,
          position = ggplot2::position_dodge()) +
          theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
        labs(x="",color="Condition", title=paste0(lab_title,"// GENE SHOWN: ", i))
        ggsave(p, filename = paste0(output_dir,"/violin/",i,"_barplot_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
        suppressWarnings(htmlwidgets::saveWidget(widget = ggplotly(p),file = paste0(output_dir,"/violin/",i,"_barplot_",gsub(".txt","",basename(z)),".html"),selfcontained = TRUE))

        df2 <- reshape::melt(df)
        df2 <- unique(df2[df2$variable=="Expr_RPKM_log2_01",c("condition","value","sample")])        
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
        ggsave(p, filename = paste0(output_dir,"/violin/",i,"_violin_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
        suppressWarnings(htmlwidgets::saveWidget(widget = ggplotly(p),file = paste0(output_dir,"/violin/",i,"_violin_",gsub(".txt","",basename(z)),".html"),selfcontained = TRUE))
        write.table(paste0("Samples_numbering:\n",paste0("Number_",1:length(df2$sample),": ",df2$sample,collapse="\n")),
                    file=paste0(output_dir,"/violin/label_samples.txt"),quote = F,row.names = F, col.names = F,sep = "\n")

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
        ggsave(p, filename = paste0(output_dir,"/violin/",i,"_violin_ttest_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
        suppressWarnings(htmlwidgets::saveWidget(widget = ggplotly(p),file = paste0(output_dir,"/violin/",i,"_violin_ttest_",gsub(".txt","",basename(z)),".html"),selfcontained = TRUE))
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
        ggsave(p, filename = paste0(output_dir,"/violin/",i,"_violin_wilcoxtest_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
        suppressWarnings(htmlwidgets::saveWidget(widget = ggplotly(p),file = paste0(output_dir,"/violin/",i,"_violin_wilcoxtest_",gsub(".txt","",basename(z)),".html"),selfcontained = TRUE))
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
        ggsave(p, filename = paste0(output_dir,"/violin/",i,"_violin_kruskaltest_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
        suppressWarnings(htmlwidgets::saveWidget(widget = ggplotly(p),file = paste0(output_dir,"/violin/",i,"_violin_kruskaltest_",gsub(".txt","",basename(z)),".html"),selfcontained = TRUE))
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
        ggsave(p, filename = paste0(output_dir,"/violin/",i,"_violin_anova_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
        suppressWarnings(htmlwidgets::saveWidget(widget = ggplotly(p),file = paste0(output_dir,"/violin/",i,"_violin_anova_",gsub(".txt","",basename(z)),".html"),selfcontained = TRUE))  
  }
  }
  }

  print(paste0("PDF barplot and violin plots done if required. Keep in mind you can also use the Expression Visualization App at https://bioinfoipbln.shinyapps.io/expressionvisualizationapp/."))
  print(paste0("Current date: ",date()))
  print("Genes highlighted are:")
  print(genes)


###### Attempt of Differential Gene Expression Analyses... modified from Bioinfo Unit to use here:
  suppressMessages(library(edgeR,quiet = T,warn.conflicts = F))
  suppressMessages(library(ggplot2,quiet = T,warn.conflicts = F))
  suppressMessages(library(plotly,quiet = T,warn.conflicts = F))
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
    suppressWarnings(htmlwidgets::saveWidget(widget = ggplotly(p),file = gsub("pdf","html",file),selfcontained = TRUE))
  }
  DGE <- function(comp,organism="mouse", myFDR=0.05, myFC=0){
    et <- exactTest(edgeR_object_norm,pair = comp)
    #Extracting the statistical data order by p-value
    top1 <- topTags(et, n=nrow(et), adjust.method="BH",sort.by="PValue")

    print(summary(decideTests(et)))
    print(nrow(top1$table[top1$table$FDR<=myFDR & abs(top1$table$logFC)>= myFC,]))
    print(comp); print("Top results:")
    print(head(top1$table,10)[,c(-1,-2)])
    myLabel1=paste(comp, collapse = '_vs_')
    myLabel1=gsub("^_","",gsub("_+","_",gsub("[^[:alnum:]_]+", "_", myLabel1)))
    myLabel2=paste(myLabel1, "FDR", myFDR, "FC", myFC, sep="_")
    # createExcel(top1$table,paste("DEG_",myLabel1,".xlsx", sep=""),organism = organism)
    # my_enrichment(top1$table,FA_label=myLabel2,cutoff=0.05,organism = organism, FDR=myFDR,FC=myFC)
    Volcano(top1,paste(output_dir,"/DGE/Volcano_plot_",myLabel1,".pdf", sep=""),myLabel1)
    return(top1)
  }


  ## Create targets:
    x <- edgeR_object
    if (targets_file!="no"){
        print(paste0("Reading targets in ",targets_file))
        targets <- read.table(file = targets_file,header = F,stringsAsFactors=F)      
    } else {
        print("Reading targets in /GEO_info/samples_info.txt")
        targets <- read.table(file = paste0(path,"/GEO_info/samples_info.txt"),header = F,stringsAsFactors=F)
        targets[match(rownames(x$samples),targets$V1),]
    }
    colnames(targets) <- c("Filename","Name","Type")
    targets$Filename <- gsub("_t|m_Rep|_seq|_KO|_WT","",paste(unique(targets$Filename,targets$Name),sep="_"))
    # Again, reorder based on GSM/SRR, if any:
    if (length(grep("_GSM[0-9]",targets$Filename))!=0){
      idx <- targets$Filename[unlist(lapply(strsplit(targets$Filename,"_+"),function(x){any(startsWith(x,"GSM"))}))]
      targets <- targets[match(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"GSM"))}))],"_+"),function(x){grep("^GSM",x,val=T)})),decreasing=T)],targets$Filename),]
    }
    if (length(grep("_SRR[0-9]",targets$Filename))!=0){
      idx <- targets$Filename[unlist(lapply(strsplit(targets$Filename,"_+"),function(x){any(startsWith(x,"SRR"))}))]
      targets <- targets[match(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"SRR"))}))],"_+"),function(x){grep("^SRR",x,val=T)})),decreasing=T)],targets$Filename),]
    }
    targets$Name <- gsub("_SRR.*","",targets$Name); targets$Filename <- targets$Name; targets <- unique(targets); targets$Name <- 1:length(targets$Name)
    cat("Please double check the following is in the same order or showing the same and correct categories than the lists above (counts and pheno/targets), and use this to interpret the plots:\n")
    print(targets)

    if (exists("gsm_manual_filter")){
        idxs_gsm_manual_3 <- which(unlist(lapply(strsplit(targets$Name,"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))    
        targets <- targets[idxs_gsm_manual_3,]
    }


  ## Proceed with DGE analyses:
    print("Attempting differential gene expression analyses between the conditions:")
    dir.create(paste0(output_dir,"/DGE"),showWarnings=F)
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
      print("These are the combinations that will be analyzed in differential gene expression analyses, if you want to restrict these, please use the argument -Dec")
      if (restrict_comparisons!="no"){
        cat("Restriction requested. Reading the comma-separated list with the indexes of the combinations that you want to keep...\n")
        indexes_restrict_combinations <- as.numeric(unlist(strsplit(restrict_comparisons, ",")))
        list_combinations <- list_combinations[indexes_restrict_combinations]
        cat("Restricted to:\n")
        print(indexes_restrict_combinations)
        print(list_combinations)
      }
      existing <- length(list.files(path=paste0(output_dir,"/DGE/"),pattern=".RData$"))
    for (i in 1:length(list_combinations)){
      # print(i)    
      if (exists("gsm_manual_filter")){
        idxs_gsm_manual_2 <- which(unlist(lapply(strsplit(colnames(gene_counts),"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))
        idxs_gsm_manual_2_2 <- which(unlist(lapply(strsplit(pheno$sample,"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))

        edgeR_object <- DGEList(counts=gene_counts[,idxs_gsm_manual_2],
                         group=pheno[idxs_gsm_manual_2_2,"condition"]$condition,
                         genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])
        edgeR_object_prefilter <- edgeR_object
        edgeR_object <- filter(filter=filter_option,edgeR_object) # Make sure of use bin to capture Cort and the lower expressed genes
        edgeR_object_norm <- calcNormFactors(edgeR_object)
        edgeR_object_norm <- estimateCommonDisp(edgeR_object_norm, robust=TRUE)
        if (is.na(edgeR_object_norm$common.dispersion)){
          edgeR_object_norm$common.dispersion <- 0.4 ^ 2
          # https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
          cat("Estimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate...")
        }
        edgeR_object_norm <- estimateTagwiseDisp(edgeR_object_norm)
        gene_counts_rpkm <- as.data.frame(rpkm(edgeR_object_norm,normalized.lib.sizes=TRUE))
        colnames(gene_counts_rpkm) <- rownames(edgeR_object_norm$samples)
        gene_counts_rpkm$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm))
      } 
      if(sum(!startsWith(as.character(edgeR_object_norm$samples$group),"__")) == length(as.character(edgeR_object_norm$samples$group))){
        edgeR_object_norm$samples$group <- as.factor(paste0("__",edgeR_object_norm$samples$group))
      }          
      write.table(paste0("\nComparison number ",i+existing,": ",list_combinations[[i]]),
              file=paste0(output_dir,"/DGE/list_comp.txt"),quote = F,row.names = F, col.names = F,sep = "\n", append=T)      
      if(covariab == "none"){
        edgeR_results <- DGE(comp=list_combinations[[i]])
        colnames(edgeR_results$table)[3] <- paste0(colnames(edgeR_results$table)[3],list_combinations[[i]][1],"_VS_",list_combinations[[i]][2])
      } else {
        fit <- glmQLFit(edgeR_object_norm, design, robust=TRUE)
        contrast <- rep(0,dim(design)[2])
        contrast[grep(paste(gsub("__|_seq1|_seq2","",list_combinations[[i]]),collapse="|"),colnames(design))] <- c(-1,1)
        qlf <- glmQLFTest(fit,contrast=contrast)
        edgeR_results <- topTags(qlf,n=nrow(qlf),adjust.method="BH",sort.by="PValue")
        print(summary(decideTests(qlf)))
        print(nrow(edgeR_results$table[edgeR_results$table$FDR<=0.05 & abs(edgeR_results$table$logFC)>= 0,]))
        print(list_combinations[[i]]); print("Top results:")
        print(head(edgeR_results$table,10)[,c(-1,-2)])
        myLabel1=paste(list_combinations[[i]], collapse = '_vs_')
        myLabel1=gsub("^_","",gsub("_+","_",gsub("[^[:alnum:]_]+", "_", myLabel1)))
        Volcano(edgeR_results,paste(output_dir,"/DGE/Volcano_plot_",myLabel1,".pdf", sep=""),myLabel1)
      }
      save.image(file=paste0(output_dir,"/DGE/DGE_analysis_comp",i+existing,".RData"))
      write.table(edgeR_results$table,
              file=paste0(output_dir,"/DGE/DGE_analysis_comp",i+existing,".txt"),quote = F,row.names = F, col.names = T,sep = "\t")
    }
    }
    print("Done")
    print(paste0("Differential gene expression analyses done. Current date: ",date()))
  
  ## Computing house-keeping/hallmark genes:
    print("Obtaining house-keeping genes...")
    suppressMessages(library(NormqPCR,quiet = T,warn.conflicts = F))
    suppressMessages(library(limma,quiet = T,warn.conflicts = F))
    setwd(paste0(output_dir,"/DGE"))    
      tryCatch({
        RPKM <- gene_counts_rpkm_to_write[,-grep("Gene_ID",colnames(gene_counts_rpkm_to_write))]
        a <- data.frame(Type=gsub("_Rep.*","",colnames(RPKM)))
        rownames(a) <- colnames(RPKM)
        write.table(a,
                    file="temp_targets.txt",quote = F,row.names = T, col.names = T,sep = "\t")
  
        sink("HK_genes.log")
        target <- readTargets("temp_targets.txt")
        matriz_obj<-new("qPCRBatch", exprs=as.matrix(RPKM))
        
        #1
        pData(matriz_obj)<-data.frame(Name=colnames(RPKM),Type=target$Type)
        Class <- as.factor(pData(matriz_obj)[,"Type"])      
        HK_normPCR_normfinder <- selectHKs(matriz_obj,Symbols=featureNames(matriz_obj),method="NormFinder",group=Class,minNrHKs=10)      
        ranking_NormFinder <- data.frame(
          rank=c(1:10),
          Name=as.character(HK_normPCR_normfinder$ranking)[1:10],
          Rho=as.numeric(HK_normPCR_normfinder$rho)[1:10],
          AveExp=as.numeric(rowMeans(RPKM[as.character(HK_normPCR_normfinder$ranking)[1:10],])),
          MedianExp=as.numeric(rowMedians(as.matrix(RPKM[as.character(HK_normPCR_normfinder$ranking)[1:10],])))
        )
        #2
        matriz_rho <- stabMeasureRho(matriz_obj, group = Class)
        matriz_rho <- sort(matriz_rho)
        ranking_Rho <- data.frame(
          rank=c(1:10),
          Name=names(matriz_rho)[1:10],
          Rho=as.numeric(matriz_rho)[1:10],
          AveExp=as.numeric(rowMeans(RPKM[names(matriz_rho)[1:10],])),
          MedianExp=as.numeric(rowMedians(as.matrix(RPKM[names(matriz_rho)[1:10],])))
        )
        #3
        print("Top 10 hallmark/house-keeping genes according to NormFinder and Rho methods, respectively:")
        print(ranking_NormFinder) # 10 genes
        print(ranking_Rho) # 10 genes
        sink()      
        
        write.table(ranking_NormFinder,file=paste0("HK_genes_normfinder.txt"),row.names = F,sep="\t")
        write.table(ranking_Rho,file=paste0("HK_genes_rho.txt"),row.names = F,sep="\t")
  
        print("Hallmark/house-keeping genes NormFinder and Rho methods combined:")
        hallmarks_comb <- intersect(ranking_NormFinder$Name,ranking_Rho$Name)
        print(hallmarks_comb)
        write.table(hallmarks_comb,file=paste0("HK_genes_combined.txt"),row.names = F,sep="\t")
        
      #Make barplots:
      #for (i in hallmarks_comb){
      #a <- RPKM[rownames(RPKM)==i,]
      #df <- data.frame(sample=names(a),
            #expr_RPKM=as.numeric(unname(a)))
      #suppressMessages(library(ggpubr,quiet = T,warn.conflicts = F))
      #p <- ggbarplot(df, x = "sample", y = "expr_RPKM",
        #add = "mean_se", label=T,lab.vjust = 4,
        #position = ggplot2::position_dodge()) +
        #theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
        #labs(x="", title=paste0("GENE SHOWN: ", i))
      #ggsave(p, filename = paste0("KOvsWT12m_hallmark_bars_",i,".pdf"),width=30, height=30)
      #}
      # qpdf::pdf_combine(input = list.files(pattern="KOvsWT12m_hallmark_bars_"),output="KOvsWT12m_hallmark_bars.pdf")
      # file.remove(list.files(pattern="KOvsWT12m_hallmark_bars_"))
      file.remove(list.files(pattern="temp_targets.txt"))
    })
    


###### Performing Venn diagrams for the DEGs:
Venn_funct <- function(files){
  list_of_tables <- lapply(files, read.delim)
  group <- sub("\\..*$", "",sub("DGE_limma_timecourse_|DGE_analysis_","",basename(files))); col.group <- as.factor(group)
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)] # Get a list of non-gray colors
    color_rgb <- col2rgb(color); luminance <- 0.299*color_rgb[1,] + 0.587*color_rgb[2,] + 0.114*color_rgb[3,]
    # Filter out light colors based on a luminance threshold
    color <- color[luminance < 400] # You can adjust the threshold value as needed
    levels(col.group) <- sample(color, nlevels(col.group))
    col.group <- as.character(col.group)
    list_of_ids <- lapply(list_of_tables,function(y){y$Gene_ID[y$FDR<0.05]})
    names(list_of_ids) <- group

    if(length(list_of_ids) < 8){
      # A more complex Venn diagram with an independent file showing the intersect:
    tmp_file <- tempfile()
    sink(tmp_file)
    hey <- nVennR::plotVenn(list_of_ids,nCycles=14000,showPlot=F)
    nVennR::showSVG(nVennObj = hey, outFile=paste0(unique(dirname(files)),"/Venn_diagram_complete.svg"),opacity=0.2,borderWidth=0.5)
    sink(); unlink(tmp_file)
    hey <- nVennR::listVennRegions(hey);names(hey) <- gsub(".* \\(","(",names(hey));names(hey) <- gsub(", ","-",names(hey));names(hey) <- gsub("\\(|\\)","",names(hey))
    write.table(data.frame(Combination=names(hey),Shared_genes=unname(unlist(lapply(hey,function(x){paste(x,collapse=",")})))),file=paste0(unique(dirname(files)),"/Venn_diagram_complete.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  } else {
    print(paste0("Too many conditions (",length(list_of_ids),"). We decided to only allow up to 7 conditions for the more complex Venn Diagram by nVennR. Otherwise, it would take too much time. You can check all the Venn Diagrams with up to 4 comparisons or perform the more complex one externally..."))
  }

  # A more typical Venn diagram. If more than 4 sets, all possible iterations
  if(length(list_of_ids)>=4){
    list_of_combinations <- combn(1:length(list_of_ids), 4, simplify = FALSE)
  } else {
    list_of_combinations <- list(unique_combn=1:length(list_of_ids))
  }
  
  for (i in 1:length(list_of_combinations)){
    write.table(paste0("Comparison number ",i,": ",paste0(names(list_of_ids[list_of_combinations[[i]]]),collapse=" // ")),file=paste0(unique(dirname(files)),"/list_combn.txt"),col.names = F,row.names = F,quote = F,sep="\n",append=T)
    list_of_ids[list_of_combinations[[i]]]
    tmp_file <- tempfile();sink(tmp_file)
      VennDiagram::venn.diagram(
            x = list_of_ids[list_of_combinations[[i]]],
            category.names = names(list_of_ids)[list_of_combinations[[i]]],
            filename = paste0(unique(dirname(files)),"/Venn_diagram_combn_",i,".png"),
            output=TRUE,
            disable.logging=T,
            
            # Output features
            imagetype="png" ,
            height = 480 , 
            width = 480 , 
            resolution = 300,
            compression = "lzw",
            
            # Circles
            lwd = 2,
            lty = 'blank',
            fill = col.group[list_of_combinations[[i]]],
            
            # Numbers
            cex = .2,
            fontface = "bold",
            fontfamily = "sans",
            
            # Set names
            cat.cex = 0.2,
            cat.fontface = "bold",
            cat.default.pos = "outer",
            cat.fontfamily = "sans")
    sink(); unlink(tmp_file)
  }
}

if(length(list.files(path=paste0(output_dir,"/DGE"),full.names=T,pattern="^DGE_analysis_comp\\d+\\.txt$"))>1){
  print("Attempting to perform Venn diagrams for DGE analyses...")
  Venn_funct(list.files(path=paste0(output_dir,"/DGE"),full.names=T,pattern="^DGE_analysis_comp\\d+\\.txt$"))
}


save.image(paste0(output_dir,"/QC_and_others/globalenvir.RData"))
###### QC PDF from Bioinfo and Laura:
  cat("\nPerforming QC_PDF...\n");print(paste0("Current date: ",date()))
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
  suppressMessages(library("ggpubr",quiet = T,warn.conflicts = F))
  suppressMessages(library("ggpmisc",quiet = T,warn.conflicts = F))

  label <- basename(path)

  x_prefilter <- edgeR_object_prefilter
  lcpm_prefilter <- cpm(x_prefilter, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...

  x <- edgeR_object
  cpm <- cpm(x) # This is normalized, altough not through edgeR, but the argument normalized.lib.sizes=TRUE by default in cpm...
  lcpm <- cpm(x, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...
  lcpm_no_log <- cpm(x, log=FALSE)

  nsamples <- ncol(x)
  
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
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)] # Get a list of non-gray colors
  color_rgb <- col2rgb(color)
  luminance <- 0.299*color_rgb[1,] + 0.587*color_rgb[2,] + 0.114*color_rgb[3,]
  # Filter out light colors based on a luminance threshold
  color <- color[luminance < 400] # You can adjust the threshold value as needed
  levels(col.group) <- sample(color, nlevels(col.group))
  col.group <- as.character(col.group)

  # summary(lcpm)
  # table(rowSums(x$counts==0)==6)

  
  ### The actual pdf file:
  pdf(paste0(output_dir,"/QC_and_others/",label,"_QC.pdf"),paper="A4")
  ### 0 Reminder of the samples:
  ggplot() + theme_void(base_size=1) + coord_flip() +
    annotate(geom = "table",
                   x = 0,
                   y = 0,
                   size = 1,
                   label = list(as.data.frame(targets)))
  ### 1.1. Density rawcounts log2, cpm...:
  col <- RColorBrewer::brewer.pal(nsamples, "Paired")
  if(sum(duplicated(col))>0){
    col <- grDevices::rainbow(nsamples)
    print(paste0("Replacing the palette with ",nsamples," random colors..."))
  }
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
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Name, cex.axis=0.4)
  title(main="Unnormalized data",ylab="Log-cpm",xlab="sample_type")
  ### 2.2. Normalised:
  boxplot(lcpm2, las=2, col=col.group, main="", names=targets$Name, cex.axis=0.4)
  title(main="Normalized data",ylab="Log-cpm",xlab="sample_type")
  ### 3. Library size:
  par(mfrow=c(1,1))
  bar_mids <- barplot(x$samples$lib.size,names.arg = gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name),las=2, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))))
  # Loop over the bar midpoints and add the text on top of each bar
  for(i in 1:length(bar_mids)) {
    # The y position is slightly above the top of the bar
    y_pos <- x$samples$lib.size[i] + 0.02 * max(x$samples$lib.size)    
    # Add the text, centered on the bar midpoint
    text(bar_mids[i], y_pos, labels = x$samples$lib.size[i], cex = 0.8, pos = 3)
  }
  ### Figures with the number of reads
  reads <- c()
  files <- grep("_stats.txt",list.files(path=input_dir,full.names=T,recursive=T),val=T)
  for (f in files){reads <- c(reads,system(paste0("cat ",f," | grep '1st fragments' | sed 's,.*:\t,,g'"),intern=T))}
  bam_reads <- data.frame(names=gsub("_nat.*","",basename(files)),reads=as.numeric(reads))

  targets_to_merge <- as.data.frame(targets)
  for (i in 1:dim(targets_to_merge)[1]){
    targets_to_merge$Filename[i] <- gsub(paste0("_",targets_to_merge$Type[i]),"",targets_to_merge$Filename[i])
  }
  bam_reads_2 <- merge(bam_reads,targets_to_merge,by.x="names",by.y="Filename")
  bam_reads_2$color <- col.group
  bam_reads_2$Name <- as.character(bam_reads_2$Name)
  bar_plot <- ggbarplot(
    bam_reads_2, 
    x = "Name", 
    y = "reads", 
    fill = "Type",
    color = "color",
    stat = "identity"
  )
  bar_plot + 
    geom_text(aes(label = reads), vjust = -0.5, color = "black", size=3) + labs(title="raw_reads") + guides(color = "none")

  reads <- c()
  files <- grep("_1_fastqc.html",list.files(path=input_dir,full.names=T,recursive=T),val=T)
  for (f in files){reads <- c(reads,system(paste0("cat ",f," | sed 's,<td>,\\n,g;s,</td>,\\n,g' | grep -A2 'Total Sequences' | tail -1"),intern=T))}
  if(length(files)!=0){ # Control that sometimes if these are repeated runs, fastqc is not going to be executed    
    fastq_reads <- data.frame(names=gsub("_1_fastqc.*","",basename(files)),reads=as.numeric(reads))

    targets_to_merge <- as.data.frame(targets)
    for (i in 1:dim(targets_to_merge)[1]){
      targets_to_merge$Filename[i] <- gsub(paste0("_",targets_to_merge$Type[i]),"",targets_to_merge$Filename[i])
    }
    fastq_reads_2 <- merge(fastq_reads,targets_to_merge,by.x="names",by.y="Filename")
    fastq_reads_2$color <- col.group
    fastq_reads_2$Name <- as.character(fastq_reads_2$Name)
    bar_plot <- ggbarplot(
      fastq_reads_2, 
      x = "Name", 
      y = "reads", 
      fill = "Type",
      color = "color",
      stat = "identity"
    )
    perc <- c()
    for (i in 1:length(reads)){perc<-c(perc,round(bam_reads_2$reads[i]*100/fastq_reads_2$reads[i],2))}
    bar_plot + 
       geom_text(aes(label = paste0(reads," (bam/fastq: ",perc, " %)"), angle=45), vjust = -0.5, color = "black", size=2) + labs(title="bam_reads") + guides(color = "none")
  }
  ### 4. Corrplot no log
  tmp <- lcpm_no_log; colnames(tmp) <- gsub("_t|m_Rep|_seq|_KO|_WT","",colnames(tmp))
  colnames(tmp) <- targets$Name[match(colnames(tmp),targets$Name)]
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
  rownames(data_pca) <- targets$Name
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
  cat("\n\nRemember that batch effect correction/covariables have been only provided to Combat-Seq for visualization purposes, if you need to include covariables in the DGE model after checking the visualization, please rerun the main program using the argument -C\n\n")
  cat("\nQC_PDF ComBat-seq counts\n")
  label <- basename(path)

  x_prefilter <- edgeR_object_prefilter_combat
  lcpm_prefilter <- cpm(x_prefilter, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...

  x <- edgeR_object_combat
  cpm <- cpm(x) # This is normalized, altough not through edgeR, but the argument normalized.lib.sizes=TRUE by default in cpm...
  lcpm <- cpm(x, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...
  lcpm_no_log <- cpm(x, log=FALSE)

  nsamples <- ncol(x)

  
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


  suppressMessages(library(ggpmisc,quiet = T,warn.conflicts = F))
  ### The actual pdf file:
  pdf(paste0(output_dir,"/QC_and_others/",label,"QC_ComBat-seq.pdf"),paper="A4")
  ### 0 Reminder of the samples:
  ggplot() + theme_void(base_size=1) + coord_flip() +
    annotate(geom = "table",
                   x = 0,
                   y = 0,
                   size = 1,
                   label = list(as.data.frame(targets)))
  
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
  rownames(data_pca) <- targets$Name
  data_pca.PC = prcomp(data_pca)
  data_pca$Type <- targets$Type
  data_pca$Filename <- targets$Filename
  data_pca$Name <- targets$Name
  data_pca$Sex <- targets$Sex
  data_pca$Age <- targets$Age
  data_pca$VAS_Group <- targets$VAS_Group
  data_pca$TypeII <- targets$TypeII
  plot(autoplot(data_pca.PC,label=T,data=data_pca,colour='Type',xlim = c(-0.8,0.8),label.size=3,label.repel=T))
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

###### WIP add DESeq2 as a full alternative to edgeR, for now, generate and provide/write independently the counts if the user ask for it:
if (diff_soft=="DESeq2"){  
  sampleFiles <- grep("Gene_ID|Length",colnames(gene_counts),invert=T)
  sampleCondition <- pheno$condition
  sampleInf <- paste0(d, "_rep", ave(pheno$condition, pheno$condition, FUN = seq_along))
  sampleTable <- data.frame(sampleName = sampleFiles,
                            condition = sampleCondition,
                            infection = sampleInf)
  sampleTable_2 <- sampleTable[,-1]
  rownames(sampleTable_2) <- sampleTable$sampleName
  sampleTable_2
  
  a <- gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)]  
  rownames(a) <- gene_counts[,"Gene_ID"]
  filter_for_deseq2 <- function(filter="standard",data,min_group=3){
    if(filter == "standard"){
      keep <- rowSums(a>1) >= min_group
      data <- data[keep,]      
    }
    else if(filter == "bin"){
      keep <- rowSums(a)>0
      data <- data[keep,]      
    }
    else{
      stop("At the moment only bin/standard are supported")
    }
    return(data)
  }
  b <- filter_for_deseq2(filter=filter_option,a)
  deseq2_object <- DESeqDataSetFromMatrix(countData = b,
                                   colData = sampleTable_2,
                                   design = ~ infection + condition)
  deseq2_object

  # Normalization:
  deseq2_object <- estimateSizeFactors(deseq2_object)
  deseq2_object_counts <- counts(deseq2_object, normalized=TRUE) # These are the normalized counts
  dim(deseq2_object_counts) # 
  dim(a) # 
  
  # Differential DESeq analysis:
  deseq2_object <- DESeq(deseq2_object)
  for (i in 1:dim(combn(unique(d),2))[2]){
    res_deseq2_object <- results(deseq2_object,contrast=c("condition",combn(unique(d),2)[,i][1],combn(unique(d),2)[,i][2]))
    table <- as.data.frame(res_deseq2_object)
    table$Gene_ID <- rownames(table)
    write.table(table,
            file=paste0(output_dir,"/DGE/DGE_analysis_DESEQ2_comp",combn(unique(d),2)[,i][1],"vs",combn(unique(d),2)[,i][2],".txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  }
}
print("ALL DONE")

print(paste0("Current date: ",date()))
