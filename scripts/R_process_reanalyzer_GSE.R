#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
input_dir <- args[2]
output_dir <- args[3]
genes <- args[4] # If not provided, "none"
filter_option <- args[5] # For now, bin or standard
organism <- args[6]
targets_file <- args[7] # if not provided, "no"
diff_soft <- args[8] # if not provided, "edgeR"
batch_format <- args[9] # if not provided, "fact"
covariab <- args[10] # if not provided, "none"
covariab_format <- args[11] # if not provided, "num"
cdseq_exec <- args[12] # if not provided, "no"
restrict_comparisons <- args[13] # if not provided, "no"
full_analyses <- args[14] # if not provided, "yes"
venn_volcano <- args[15] # if not provided, "none"
pattern_to_remove <- args[16] # if not provided, "no"

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
    gene_counts <- gene_counts[,c(gtools::mixedorder(grep("Length|Gene_ID",colnames(gene_counts),invert=T)),grep("Length|Gene_ID",colnames(gene_counts)))]
  } else {
    idx <- colnames(gene_counts)[unlist(lapply(strsplit(colnames(gene_counts),"_|__"),function(x){any(startsWith(x,"GSM"))}))]
    gene_counts <- gene_counts[,c(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"GSM"))}))],"_+"),function(x){grep("^GSM",x,val=T)})))],"Gene_ID","Length")]
  }
  # This is crucial, also double checking reordering based on SRR* if present in the names:
  if (length(grep("_SRR[0-9]",colnames(gene_counts)))!=0){
    idx <- colnames(gene_counts)[unlist(lapply(strsplit(colnames(gene_counts),"_|__"),function(x){any(startsWith(x,"SRR"))}))]
    gene_counts <- gene_counts[,c(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"SRR"))}))],"_+"),function(x){grep("^SRR",x,val=T)})))],"Gene_ID","Length")]
  }
  dir.create(output_dir,showWarnings=F)
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

###### Batch effect correction/count adjustment if requested:
  # First check if covariable has been provided, and then that will be used with limma's removeBatchEffect() to perform an adjustment. 
  # If batch has also been provided, then it will be overwritten immediately below, as ComBat_seq is the preferred method
  if(covariab != "none"){
    count_matrix <- as.matrix(gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)])
    if(covariab_format == "fact"){
      cov <- as.factor(unlist(strsplit(as.character(covariab),",")))
      cov <- model.matrix(~ cov)
    } else if(covariab_format == "num"){
      cov <- as.numeric(unlist(strsplit(as.character(covariab),",")))      
    }    
    adjusted_counts <- limma::removeBatchEffect(x=count_matrix,covariates=cov)    
  }
  if (file.exists(paste0(path,"/GEO_info/batch_vector.txt")) && !file.exists(paste0(path,"/GEO_info/batch_biological_variables.txt"))){
    count_matrix <- as.matrix(gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)])
    batch <- data.table::fread(paste0(path,"/GEO_info/batch_vector.txt"),head=F,sep="*")$V1
    batch <- unlist(strsplit(batch,","))
    if(covariab_format == "fact"){
      adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.factor(batch), group=NULL)
    } else {
      adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.numeric(batch), group=NULL)
    }
    write.table(batch,file=paste0(path,"/GEO_info/batch_vector.txt"),quote = F,row.names = F, col.names = F,sep = "\n"); print("BATCH DONE_1")    
  } else if (file.exists(paste0(path,"/GEO_info/batch_vector.txt")) && file.exists(paste0(path,"/GEO_info/batch_biological_variables.txt"))){
    count_matrix <- as.matrix(gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)])
    batch <- data.table::fread(paste0(path,"/GEO_info/batch_vector.txt"),head=F,sep="*")$V1
    batch <- unlist(strsplit(batch,","))
    biological_cov <- data.table::fread(paste0(path,"/GEO_info/batch_biological_variables.txt"),head=F,sep="*")$V1
    if (stringr::str_count(biological_cov," ") == 0){
      group <- unlist(strsplit(biological_cov,","))
      if(covariab_format == "fact"){
        adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.factor(batch), group=as.factor(group))
      } else {
        adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.numeric(batch), group=as.numeric(group))
      }
      write.table(batch,file=paste0(path,"/GEO_info/batch_vector.txt"),quote = F,row.names = F, col.names = F,sep = "\n")
      write.table(group,file=paste0(path,"/GEO_info/batch_biological_variables.txt"),quote = F,row.names = F, col.names = F,sep = "\n"); print("BATCH DONE_2")      
    } else {
      covar_mat <- matrix(rep(1,stringr::str_count(gsub(" .*","",biological_cov),",")+1))
      for (i in 1:(stringr::str_count(biological_cov," ") + 1)){
        covar_mat <- cbind(covar_mat,as.numeric(unlist(strsplit(strsplit(biological_cov," ")[[1]][i],","))))
      }
      covar_mat <- covar_mat[,-1]
      if(covariab_format == "fact"){
        adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.factor(batch), group=NULL, covar_mod=covar_mat)
      } else {
        adjusted_counts <- sva::ComBat_seq(count_matrix, batch=as.numeric(batch), group=NULL, covar_mod=covar_mat)
      }
      write.table(batch,file=paste0(path,"/GEO_info/batch_vector.txt"),quote = F,row.names = F, col.names = F,sep = "\n")
      write.table(covar_mat,file=paste0(path,"/GEO_info/batch_biological_variables.txt"),quote = F,row.names = F, col.names = F,sep = "\t"); print("BATCH DONE_3")      
    }
  }
  quantile.normalize <- function(x, target){
    sort(target)[rank(x)]
  }
  # Before writing the counts, quantile normalization for the adjusted counts if negative values present:
  if (exists("adjusted_counts")){  
    if(any(adjusted_counts<0)){
        counts.quant <- matrix(0, nrow=nrow(count_matrix), ncol=ncol(count_matrix))
        for (i in 1:ncol(counts.quant)){
          counts.quant[,i] <- quantile.normalize(adjusted_counts[,i], count_matrix[,i])
        }
        adjusted_counts <- counts.quant
      }
    write.table(adjusted_counts,
                file=paste0(output_dir,"/counts_adjusted.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  }  

  pheno <- as.data.frame(data.table::fread(paste0(path,"/GEO_info/samples_info.txt"),head=F))
  if(dim(pheno)[2]>2){
    pheno <- pheno[,2:3]
  }
  colnames(pheno) <- c("sample","condition")
  pheno <- pheno[gtools::mixedorder(pheno$sample),]
  # Ensure reordering based on GSM* or SRR* if present in the names of pheno:
  if (length(grep("GSM[0-9]",pheno$sample))!=0){
    idx <- pheno$sample[unlist(lapply(strsplit(pheno$sample,"_+"),function(x){any(startsWith(x,"GSM"))}))]
    pheno <- pheno[match(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"GSM"))}))],"_+"),function(x){grep("^GSM",x,val=T)})))],pheno$sample),]
  }
  if (length(grep("SRR[0-9]",pheno$sample))!=0){
    idx <- pheno$sample[unlist(lapply(strsplit(pheno$sample,"_+"),function(x){any(startsWith(x,"SRR"))}))]
    pheno <- pheno[match(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"SRR"))}))],"_+"),function(x){grep("^SRR",x,val=T)})))],pheno$sample),]
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
  edgeR_object <- DGEList(counts=gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)],
                   group=pheno$condition,
                   genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])
  cat("\n\nPlease note that counts have been normalized and figures have been performed using the groups and samples:\n"); print(pheno$condition); print(grep("Gene_ID|Length",colnames(gene_counts),invert=T,val=T))

###### Filter counts:
  filter <- function(filter="standard",data,min_group=3){
    if(filter == "standard"){
      print("Applying standard filtering...")
      keep <- rowSums(cpm(data)>1) >= min_group
      data <- data[keep,]
      data$samples$lib.size <- colSums(data$counts)
    }
    else if(filter == "bin"){
      print("Applying bin filtering...")
      keep <- rowSums(data$counts)>0
      data <- data[keep,]
      data$samples$lib.size <- colSums(data$counts)
    }
    else if(filter == "filterbyexpr"){
      print("Applying filterByExpr by edgeR together with the grouping...")
      keep <- filterByExpr(y,group=data$samples$group)
      data <- data[keep,]
    }
    else{
      stop("At the moment only bin/standard/filterByExpr, are supported")
    }
    return(data)
  }
  cat(paste0("Chosen filter is ",filter_option,"\n"))
  cat(paste0("Number of genes before filter: ", nrow(edgeR_object)),"\n")
  edgeR_object_prefilter <- edgeR_object
  edgeR_object <- filter(filter=filter_option,edgeR_object) # Make sure of use bin to capture Cort and the lower expressed genes
  cat(paste0("Number of genes after filter: ", nrow(edgeR_object)),"\n")
  # edgeR_object_norm <- calcNormFactors(edgeR_object)
  edgeR_object_norm <- normLibSizes(edgeR_object)

  if(covariab == "none"){
    #edgeR_object_norm <- estimateCommonDisp(edgeR_object_norm, robust=TRUE)
    edgeR_object_norm <- estimateDisp(edgeR_object_norm, robust=TRUE)
    if (is.na(edgeR_object_norm$common.dispersion)){
      edgeR_object_norm$common.dispersion <- 0.4 ^ 2
      cat("\nEstimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate. Please consult the vignette for further information...\n")
    }
    #edgeR_object_norm <- estimateTagwiseDisp(edgeR_object_norm)
  } else { 
    if(filter_option == "filterbyexpr"){
      print("Applying filterByExpr by edgeR together with the covariables in the model...Overwriting previous application...")
      edgeR_object <- edgeR_object_prefilter
      Treat <- edgeR_object$samples$group
      if(covariab_format == "fact"){
        Time <- as.factor(unlist(strsplit(as.character(covariab),",")))
      } else if(covariab_format == "num"){
        Time <- as.numeric(unlist(strsplit(as.character(covariab),",")))
      }      
      design <- model.matrix(~0+Treat+Time)
      keep <- filterByExpr(edgeR_object,design)
      edgeR_object <- edgeR_object[keep,,keep.lib.sizes=FALSE]
      edgeR_object_norm <- normLibSizes(edgeR_object)
    } else {
      Treat <- edgeR_object$samples$group
      if(covariab_format == "fact"){
        Time <- as.factor(unlist(strsplit(as.character(covariab),",")))
      } else if(covariab_format == "num"){
        Time <- as.numeric(unlist(strsplit(as.character(covariab),",")))
      }      
      design <- model.matrix(~0+Treat+Time)
      rownames(design) <- colnames(edgeR_object_norm)
      edgeR_object_norm <- estimateDisp(edgeR_object_norm, design, robust=TRUE) # Preparing for one covariable following edgeR vignette new methods Nov2023
    }
  }

  gene_counts_rpkm <- as.data.frame(rpkm(edgeR_object_norm,normalized.lib.sizes=TRUE))
  colnames(gene_counts_rpkm) <- rownames(edgeR_object_norm$samples)
  gene_counts_rpkm$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm))

  # Function for saving tpms later:
  rpkm_to_tpm <- function(rpkm) {    
    rpkm_sum <- sum(rpkm)
    tpm <- (rpkm / rpkm_sum) * 1e6    
    return(tpm)
  }

###### Reorder so gene_counts columns follow the alfanumeric order, easier for the users in the written tables, although for the figures the script needs GSMXXXX-ordered:
  gene_counts_rpkm_to_write <- gene_counts_rpkm[,c(grep("Gene_ID",colnames(gene_counts_rpkm)),grep("Gene_ID",colnames(gene_counts_rpkm),invert=T))]
  gene_counts_rpkm_to_write <- gene_counts_rpkm_to_write[,c("Gene_ID",sort(colnames(gene_counts_rpkm_to_write)[-1]))]
  write.table(gene_counts_rpkm_to_write,
              file=paste0(output_dir,"/RPKM_counts_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
  tpm_counts <- rpkm_to_tpm(gene_counts_rpkm_to_write[,grep("Gene_ID",colnames(gene_counts_rpkm_to_write),invert=T)]); rownames(tpm_counts) <- gene_counts_rpkm_to_write$Gene_ID
  write.table(tpm_counts,
              file=paste0(output_dir,"/TPM_counts_genes.txt"),quote = F,row.names = T, col.names = T,sep = "\t")
  write.table(log2(tpm_counts+0.1),
              file=paste0(output_dir,"/TPM_counts_genes_log2_0.1.txt"),quote = F,row.names = T, col.names = T,sep = "\t")
  # High/medium/low categ:
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

  
###### Similar to above, obtain RPKM counts but from the adjusted counts instead of the raw counts (ComBat-Seq, edgeR...)
  if (exists("adjusted_counts")){
    edgeR_object_adjusted <- DGEList(counts=adjusted_counts,
                                   group=pheno$condition,
                                   genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])
    cat(paste0("Counts adjustment... Number of genes before filter: ", nrow(edgeR_object_adjusted)),"\n")
    edgeR_object_prefilter_adjusted <- edgeR_object_adjusted
    edgeR_object_adjusted <- filter(filter=filter_option,edgeR_object_adjusted) # Make sure of use bin to capture Cort and the lower expressed genes
    cat(paste0("Counts adjustment... Number of genes after filter: ", nrow(edgeR_object_adjusted)),"\n")
    #edgeR_object_norm_adjusted <- calcNormFactors(edgeR_object_adjusted)
    #edgeR_object_norm_adjusted <- estimateCommonDisp(edgeR_object_norm_adjusted, robust=TRUE)
    edgeR_object_norm_adjusted <- normLibSizes(edgeR_object_adjusted)
    edgeR_object_norm_adjusted <- estimateDisp(edgeR_object_norm_adjusted, robust=TRUE)
    if (is.na(edgeR_object_norm_adjusted$common.dispersion)){
      edgeR_object_norm_adjusted$common.dispersion <- 0.4 ^ 2
      # https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
      cat("Estimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate...")
    }
    #edgeR_object_norm_adjusted <- estimateTagwiseDisp(edgeR_object_norm_adjusted)
    gene_counts_rpkm_adjusted <- as.data.frame(rpkm(edgeR_object_norm_adjusted,normalized.lib.sizes=TRUE))
    colnames(gene_counts_rpkm_adjusted) <- rownames(edgeR_object_norm_adjusted$samples)
    gene_counts_rpkm_adjusted$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm_adjusted))

    gene_counts_rpkm_adjusted_to_write <- gene_counts_rpkm_adjusted[,c(grep("Gene_ID",colnames(gene_counts_rpkm_adjusted)),grep("Gene_ID",colnames(gene_counts_rpkm_adjusted),invert=T))]
    gene_counts_rpkm_adjusted_to_write <- gene_counts_rpkm_adjusted_to_write[,c("Gene_ID",sort(colnames(gene_counts_rpkm_adjusted_to_write)[-1]))]
    write.table(gene_counts_rpkm_adjusted_to_write,
          file=paste0(output_dir,"/RPKM_counts_adjusted_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
    gene_counts_rpkm_adjusted_log <- log2(gene_counts_rpkm_adjusted_to_write[,grep("Gene_ID",colnames(gene_counts_rpkm_adjusted_to_write),invert=T)] + 0.1)
    gene_counts_rpkm_adjusted_log <- cbind(gene_counts_rpkm_adjusted_to_write$Gene_ID,gene_counts_rpkm_adjusted_log)
    colnames(gene_counts_rpkm_adjusted_log)[1] <- "Gene_ID"
    write.table(gene_counts_rpkm_adjusted_log,
          file=paste0(output_dir,"/RPKM_counts_adjusted_genes_log2_0.1.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
    tpm_counts <- rpkm_to_tpm(gene_counts_rpkm_adjusted_to_write[,grep("Gene_ID",colnames(gene_counts_rpkm_to_write),invert=T)]); rownames(tpm_counts) <- gene_counts_rpkm_adjusted_to_write$Gene_ID
    write.table(tpm_counts,
                file=paste0(output_dir,"/TPM_counts_adjusted_genes.txt"),quote = F,row.names = T, col.names = T,sep = "\t")
    write.table(log2(tpm_counts+0.1),
                file=paste0(output_dir,"/TPM_counts_adjusted_genes_log2_0.1.txt"),quote = F,row.names = T, col.names = T,sep = "\t")
  }
  cat(paste0("\nCounts written. Current date: ",date()))

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
    gene_counts_rpkm_to_plot <- gene_counts_rpkm_adjusted
  }

  for (i in unlist(strsplit(genes,","))){
    for (z in list.files(pattern = "design_possible_full", recursive = TRUE, full.names=T, path=paste0(path,"/GEO_info"))){
      if (length(unique(read.table(z,head=F,blank.lines.skip=T)$V1)) == 1){
        next
      }
      system(paste0("sed -i -e 's/[^a-zA-Z0-9]/_/g' -e '/^$/s//_/' ",z))
      cat("\n\nPerforming figures for "); print(i); cat("\n\nDesign:"); print(z); print(read.table(z)$V1) # Deal with the designs, also dealing with special characters
      if (i %in% gene_counts_rpkm_to_plot$Gene_ID){
        a <- gene_counts_rpkm_to_plot[gene_counts_rpkm_to_plot$Gene_ID==i,-grep("Gene_ID",colnames(gene_counts_rpkm_to_plot))]
        a <- a[,gtools::mixedsort(colnames(a))]
        df <- data.frame(sample=names(a),
                         expr_RPKM=as.numeric(unname(a)),
                         condition=gtools::mixedsort(read.table(z,head=F,blank.lines.skip=FALSE)$V1))
      ### If the user has requested to analyze only some of the samples/column based on the GSMXXXXX, filter and restrict here.
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

  cat(paste0("\n\nPDF barplots and violin plots done if required. Keep in mind you can also use the Expression Visualization App at https://bioinfoipbln.shinyapps.io/expressionvisualizationapp/."))
  cat("\n\nGenes highlighted are:");print(genes);print(paste0("Current date: ",date()))


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
    colnames(data)[grep("logFC",colnames(data))] <- "logFC"
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
  DGE <- function(comp,object=edgeR_object_norm, organism="mouse", myFDR=0.05, myFC=0){
    et <- exactTest(object,pair = rev(comp))
    #Extracting the statistical data order by p-value
    top1 <- topTags(et, n=nrow(et), adjust.method="BH",sort.by="PValue")
    print(summary(decideTests(et)))
    print(nrow(top1$table[top1$table$FDR<=myFDR & abs(top1$table$logFC)>= myFC,]))
    print(comp); print("Top 10 results (each sense):"); top2 <- top1[top1$table$FDR<=myFDR,3:6]
    if(dim(top2)[1] > 20){
      print(as.data.frame(top2)[order(as.data.frame(top2)$logFC)[c(as.numeric(dim(top2)[1]):as.numeric(dim(top2)[1]-10),1:10)],])
    } else {
      print(as.data.frame(top2)[order(as.data.frame(top2)$logFC),])
    }   
    if (venn_volcano!="no"){
      myLabel1=gsub("^_","",gsub("_+","_",gsub("[^[:alnum:]_]+", "_", paste(comp, collapse = '_vs_'))))
      Volcano(top1,paste(output_dir,"/DGE/Volcano_plot_",myLabel1,".pdf", sep=""),myLabel1)
    }
    # myLabel2=paste(myLabel1, "FDR", myFDR, "FC", myFC, sep="_")
    # createExcel(top1$table,paste("DEG_",myLabel1,".xlsx", sep=""),organism = organism)
    # my_enrichment(top1$table,FA_label=myLabel2,cutoff=0.05,organism = organism, FDR=myFDR,FC=myFC)
    
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
    }
    colnames(targets) <- c("Filename","Name","Type")
    targets$Filename <- gsub("_t|m_Rep|_seq|_KO|_WT","",paste(unique(targets$Filename,targets$Name),sep="_"))    
    # Again, reorder based on GSM/SRR, if any:
    if (length(grep("_GSM[0-9]",targets$Filename))!=0){
      idx <- targets$Filename[unlist(lapply(strsplit(targets$Filename,"_+"),function(x){any(startsWith(x,"GSM"))}))]
      targets <- targets[match(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"GSM"))}))],"_+"),function(x){grep("^GSM",x,val=T)})))],targets$Filename),]
    }
    if (length(grep("SRR[0-9]",targets$Filename))!=0){
      idx <- targets$Filename[unlist(lapply(strsplit(targets$Filename,"_+"),function(x){any(startsWith(x,"SRR"))}))]
      targets <- targets[match(idx[gtools::mixedorder(unlist(lapply(strsplit(idx[unlist(lapply(strsplit(idx,"_|__"),function(x){any(startsWith(x,"SRR"))}))],"_+"),function(x){grep("^SRR",x,val=T)})))],targets$Filename),]
    }    
    bakktarget <- targets
    targets$Name <- gsub("_SRR.*","",targets$Name); targets$Filename <- targets$Name; targets <- unique(targets); targets$Name <- 1:length(targets$Name)

    ### Checking that the automatic ordering, which may be quite messy, has worked...
    cat("\n\n\nPLEASE double check that all the following is showing the correct order after automatically extracting and ordering data from the counts, pheno/targets data... This is CRUCIAL and an incorrect automatic ordering could INVALIDATE ALL RESULTS. If anything does not look like it should, please ask for support...\n")
    cat("The names containing numbers can be used to interpret the plots...\n")
    print(colnames(gene_counts)[grep("Gene_ID|Length",colnames(gene_counts),invert=T)])
    cat("\n\n"); print(pheno, row.names=F)
    cat("\n\n"); print(bakktarget, row.names=F)
    cat("\n\n"); print(targets, row.names=F)
    cat("\n\n\n")

    if (exists("gsm_manual_filter")){
        idxs_gsm_manual_3 <- which(unlist(lapply(strsplit(targets$Name,"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))    
        targets <- targets[idxs_gsm_manual_3,]
    }

    write.table(as.data.frame(targets),file=paste0(output_dir,"/QC_and_others/targets_info_file.txt"),quote = F,row.names = F, col.names = T,sep = "\t")

  ## Proceed with DGE analyses:
    print("Attempting differential gene expression analyses between the conditions:")
    dir.create(paste0(output_dir,"/DGE"),showWarnings=F)
    files_designs <- list.files(pattern = "design_possible_full", recursive = TRUE, full.names=T, path=paste0(path,"/GEO_info"))
    for (n in 1:length(files_designs)){
      z <- files_designs[n]
      
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

      # Create edgeR_objects for each grouping/design if required:
      edgeR_object_tmp <- DGEList(counts=gene_counts[,grep("Gene_ID|Length",colnames(gene_counts),invert=T)],
                               		group=condition,
                               		genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])
      edgeR_object_tmp <- filter(filter=filter_option,edgeR_object_tmp)
      edgeR_object_norm_temp <- normLibSizes(edgeR_object_tmp)
      # edgeR_object_norm_temp <- calcNormFactors(edgeR_object_tmp)      
      if(covariab == "none"){
          #edgeR_object_norm_temp <- estimateCommonDisp(edgeR_object_norm_temp, robust=TRUE)
          edgeR_object_norm_temp <- estimateDisp(edgeR_object_norm_temp, robust=TRUE)
          if (is.na(edgeR_object_norm_temp$common.dispersion)){
            edgeR_object_norm_temp$common.dispersion <- 0.4 ^ 2
            cat("\nEstimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate. Please consult the vignette for further information...\n")
          }
          #edgeR_object_norm_temp <- estimateTagwiseDisp(edgeR_object_norm_temp)
      } else { 
          Treat <- edgeR_object_tmp$samples$group
          if(covariab_format == "fact"){
            Time <- as.factor(unlist(strsplit(as.character(covariab),",")))
          } else if(covariab_format == "num"){
            Time <- as.numeric(unlist(strsplit(as.character(covariab),",")))
          }
          design <- model.matrix(~0+Treat+Time)
          rownames(design) <- colnames(edgeR_object_norm_temp)
          edgeR_object_norm_temp <- estimateDisp(edgeR_object_norm_temp, design, robust=TRUE) # Preparing for one covariable following edgeR vignette new methods Nov2023          
      }
      
      edgeR_object_norm_temp_to_process <- edgeR_object_norm_temp

      # Process each combination
      list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(condition), unique(condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T)
      list_combinations <- lapply(list_combinations,function(x){if (sum(!startsWith(x,"__"))==2){paste0("__",x)} else {x}})
      print("These are the combinations that will be analyzed in differential gene expression analyses, if you want to restrict or change these, please use the argument -Dec or -pR")
      if (restrict_comparisons!="no"){
        cat("You have manually provided a list. Reading the comma-separated list with the ordered comparisons that you want to perform...\n")
        cat("These are:\n")
        list_combinations <- strsplit(unlist(strsplit(restrict_comparisons,",")),"//")        
        print(list_combinations)
      }
      if (pattern_to_remove!="none"){
        list_combinations <- list_combinations[grep(pattern_to_remove,list_combinations,invert=T)]
      }
      existing <- length(list.files(path=paste0(output_dir,"/DGE/"),pattern=".RData$"))
      print(list_combinations)
    for (i in 1:length(list_combinations)){
      # print(i)    
      # This if should be updated and probably moved above if the filter is to be used properly
      if (exists("gsm_manual_filter")){
        idxs_gsm_manual_2 <- which(unlist(lapply(strsplit(colnames(gene_counts),"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))
        idxs_gsm_manual_2_2 <- which(unlist(lapply(strsplit(pheno$sample,"_"),function(x){any(x %in% unlist(strsplit(gsm_manual_filter,",")))})))

        edgeR_object <- DGEList(counts=gene_counts[,idxs_gsm_manual_2],
                         group=pheno[idxs_gsm_manual_2_2,"condition"]$condition,
                         genes=gene_counts[,c(grep("Gene_ID",colnames(gene_counts)),grep("Length",colnames(gene_counts)))])
        edgeR_object_prefilter <- edgeR_object
        edgeR_object <- filter(filter=filter_option,edgeR_object) # Make sure of use bin to capture Cort and the lower expressed genes
        #edgeR_object_norm <- calcNormFactors(edgeR_object)
        #edgeR_object_norm <- estimateCommonDisp(edgeR_object_norm, robust=TRUE)
        edgeR_object_norm <- normLibSizes(edgeR_object_norm)
        edgeR_object_norm <- estimateDisp(edgeR_object_norm, robust=TRUE)
        if (is.na(edgeR_object_norm$common.dispersion)){
          edgeR_object_norm$common.dispersion <- 0.4 ^ 2
          # https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
          cat("Estimating Dispersion... Errors or warnings? Likely because no replicates, addressing providing a fixed value for dispersion, but do not trust comparative analyses because it's likely not accurate...")
        }
        #edgeR_object_norm <- estimateTagwiseDisp(edgeR_object_norm)
        gene_counts_rpkm <- as.data.frame(rpkm(edgeR_object_norm,normalized.lib.sizes=TRUE))
        colnames(gene_counts_rpkm) <- rownames(edgeR_object_norm$samples)
        gene_counts_rpkm$Gene_ID <- stringr::str_to_title(rownames(gene_counts_rpkm))
      } 
      
      if(sum(!startsWith(as.character(edgeR_object_norm$samples$group),"__")) == length(as.character(edgeR_object_norm$samples$group))){
        edgeR_object_norm$samples$group <- as.factor(paste0("__",edgeR_object_norm$samples$group))
      }
      if(sum(!startsWith(as.character(edgeR_object_norm_temp_to_process$samples$group),"__")) == length(as.character(edgeR_object_norm_temp_to_process$samples$group))){
        edgeR_object_norm_temp_to_process$samples$group <- as.factor(paste0("__",edgeR_object_norm_temp_to_process$samples$group))
      }
      
      write.table(paste0("\nComparison number ",i+existing,": ",list_combinations[[i]]),
              file=paste0(output_dir,"/DGE/list_comp.txt"),quote = F,row.names = F, col.names = F,sep = "\n", append=T)

      # Perform DGE:
      if (full_analyses!="no"){
        if(covariab == "none"){
          edgeR_results <- DGE(comp=list_combinations[[i]],object=edgeR_object_norm_temp_to_process)
          colnames(edgeR_results$table)[3] <- paste0(colnames(edgeR_results$table)[3],paste(sub("__","",list_combinations[[i]]),collapse = "__VS__"))
        } else {
          fit <- glmQLFit(edgeR_object_norm_temp_to_process, design, robust=TRUE)
          contrast <- rep(0,dim(design)[2])
          idxs_design <- rev(c(grep(sub("__","",list_combinations[[i]][1]),colnames(design)),grep(sub("__","",list_combinations[[i]][2]),colnames(design))))
          if(length(idxs_design)!=2){cat(paste0("\n\nSomething is WRONG as one of your contrasts has been required to compare ",length(idxs_design)," conditions. Probably conflicting naming of biological conditions...\n\n"));stop("Exiting, please review the naming of the conditions...")}
          contrast[idxs_design] <- c(-1,1)
          qlf <- glmQLFTest(fit,contrast=contrast)
          edgeR_results <- topTags(qlf,n=nrow(qlf),adjust.method="BH",sort.by="PValue")
          print(summary(decideTests(qlf)))
          print(nrow(edgeR_results$table[edgeR_results$table$FDR<=0.05 & abs(edgeR_results$table$logFC)>= 0,]))
          print(list_combinations[[i]]); print("Top 10 results (each sense):")
          print(as.data.frame(edgeR_results)[order(as.data.frame(edgeR_results)$logFC)[c(as.numeric(dim(edgeR_results)[1]):as.numeric(dim(edgeR_results)[1]-10),1:10)],3:6])
          colnames(edgeR_results$table)[3] <- paste0(colnames(edgeR_results$table)[3],paste(sub("__","",list_combinations[[i]]),collapse = "__VS__"))
          if (venn_volcano!="no"){
            myLabel1=paste(list_combinations[[i]], collapse = '_vs_')
            myLabel1=gsub("^_","",gsub("_+","_",gsub("[^[:alnum:]_]+", "_", myLabel1)))
            Volcano(edgeR_results,paste(output_dir,"/DGE/Volcano_plot_",myLabel1,".pdf", sep=""),myLabel1)
          }
        }
        conflicts <- intersect(ls(envir = environment()), ls(envir = .GlobalEnv))
        if (length(conflicts) > 0) {warning("The following objects will be overwritten in the global environment: ", paste(conflicts, collapse = ", "))}
        list2env(as.list(environment()), envir = .GlobalEnv); save.image(file=paste0(output_dir,"/DGE/DGE_analysis_comp",i+existing,".RData"))
        write.table(edgeR_results$table,
                file=paste0(output_dir,"/DGE/DGE_analysis_comp",i+existing,".txt"),quote = F,row.names = F, col.names = T,sep = "\t")
      }
    }
    }
    print("Done")
    print(paste0("Differential gene expression analyses done. Current date: ",date()))
  
  ## Computing house-keeping/hallmark genes:
    if (full_analyses!="no"){
      print("Obtaining house-keeping genes...")
      suppressMessages(library(NormqPCR,quiet = T,warn.conflicts = F))
      suppressMessages(library(limma,quiet = T,warn.conflicts = F))
      setwd(paste0(output_dir,"/DGE"))
        tryCatch({
          RPKM <- gene_counts_rpkm_to_write[,-grep("Gene_ID",colnames(gene_counts_rpkm_to_write))]
          if (pattern_to_remove!="none"){
            RPKM <- RPKM[,grep(pattern_to_remove,colnames(RPKM),invert=T)]
          }
          a <- data.frame(Type=sub("_Rep.*","",sub("(.*)_.*", "\\1",colnames(RPKM))))
          rownames(a) <- colnames(RPKM)
          print("Obtaining targets from colum names removing everything after the last underline... This is the result, please double check as it may be the source of errors...")
          print(a)
          write.table(a,
                      file="temp_targets.txt",quote = F,row.names = T, col.names = T,sep = "\t")
    
          #sink("HK_genes.log")
          target <- readTargets("temp_targets.txt")
          matriz_obj<-new("qPCRBatch", exprs=as.matrix(RPKM))
          
          ### Get a number of HK genes: 10 by default
          hk_genes_number_input <- 10
          #max_number_degs <- max(unname(sapply(list.files(path=output_dir,pattern="DGE_analysis_comp*",recursive=T,full=T),function(x){length(unique(read.delim(x)$Gene_ID[read.delim(x)$FDR<0.05]))})))
          #for (hk_genes_number in c(hk_genes_number_input,100,max_number_degs)){
          for (hk_genes_number in c(hk_genes_number_input,100)){
            #1
            pData(matriz_obj)<-data.frame(Name=colnames(RPKM),Type=target$Type)
            Class <- as.factor(pData(matriz_obj)[,"Type"])
            HK_normPCR_normfinder <- selectHKs(matriz_obj,Symbols=featureNames(matriz_obj),method="NormFinder",group=Class,minNrHKs=hk_genes_number,trace=F)
            ranking_NormFinder <- data.frame(
              rank=c(1:hk_genes_number),
              Name=as.character(HK_normPCR_normfinder$ranking)[1:hk_genes_number],
              Rho=as.numeric(HK_normPCR_normfinder$rho)[1:hk_genes_number],
              AveExp=as.numeric(rowMeans(RPKM[as.character(HK_normPCR_normfinder$ranking)[1:hk_genes_number],])),
              MedianExp=as.numeric(rowMedians(as.matrix(RPKM[as.character(HK_normPCR_normfinder$ranking)[1:hk_genes_number],])))
            )
            #2
            matriz_rho <- stabMeasureRho(matriz_obj, group = Class)
            matriz_rho <- sort(matriz_rho)
            ranking_Rho <- data.frame(
              rank=c(1:hk_genes_number),
              Name=names(matriz_rho)[1:hk_genes_number],
              Rho=as.numeric(matriz_rho)[1:hk_genes_number],
              AveExp=as.numeric(rowMeans(RPKM[names(matriz_rho)[1:hk_genes_number],])),
              MedianExp=as.numeric(rowMedians(as.matrix(RPKM[names(matriz_rho)[1:hk_genes_number],])))
            )
            #3
            #print(paste0("Top ",hk_genes_number," hallmark/house-keeping genes according to NormFinder and Rho methods, respectively:"))
            #print(ranking_NormFinder)
            #print(ranking_Rho)          
            
            write.table(ranking_NormFinder,file=paste0("HK_genes_normfinder_",hk_genes_number,".txt"),row.names = F,sep="\t",quote = F)
            write.table(ranking_Rho,file=paste0("HK_genes_rho_",hk_genes_number,".txt"),row.names = F,sep="\t",quote = F)
      
            #print("Hallmark/house-keeping genes NormFinder and Rho methods combined:")
            hallmarks_comb <- intersect(ranking_NormFinder$Name,ranking_Rho$Name)
            #print(hallmarks_comb)
            write.table(hallmarks_comb,file=paste0("HK_genes_combined_",hk_genes_number,".txt"),row.names = F,sep="\t",quote = F,col.names=F)
            print(paste0("Computed ",hk_genes_number," HK genes!"))
          }
          
          #sink()
  
          
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
        }, error = function(e) {
            writeLines(as.character(e), paste0("Housekeeping_err.txt"))
        })
      }
    


###### Performing Venn diagrams for the DEGs:
if (venn_volcano!="no"){  
  Venn_funct <- function(files){
      list_of_tables <- lapply(files, read.delim)
      group <- sub("\\..*$", "",sub("DGE_limma_timecourse_|DGE_analysis_","",basename(files))); col.group <- as.factor(group)
      color = grDevices::colors()[grep('gr(a|e)y|white', grDevices::colors(), invert = TRUE)] # Get a list of non-gray or white colors
      contrast <- sapply(color,colorspace::contrast_ratio); contrast <- contrast[contrast>4] # Ensure a high contrast here and below (>4 on W3C standard)
      contrast2 <- unique(t(combn(unique(names(contrast)),2))[apply(t(combn(unique(names(contrast)),2)),1,function(x){colorspace::contrast_ratio(x[1],col2=x[2])}) > 4])
      levels(col.group) <- sample(contrast2, nlevels(col.group)); col.group <- as.character(col.group)
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
  
  if (full_analyses!="no"){
    if(length(list.files(path=paste0(output_dir,"/DGE"),full.names=T,pattern="^DGE_analysis_comp\\d+\\.txt$"))>1){
      print("Attempting to perform Venn diagrams for DGE analyses...")
      Venn_funct(list.files(path=paste0(output_dir,"/DGE"),full.names=T,pattern="^DGE_analysis_comp\\d+\\.txt$"))
    }
    
    try(system("tar cf venn_diagrams.tar Venn_diagram* 2>/dev/null; rm Venn_diagram* 2>/dev/null"))
  }
}
conflicts <- intersect(ls(envir = environment()), ls(envir = .GlobalEnv))
if (length(conflicts) > 0) {warning("The following objects will be overwritten in the global environment: ", paste(conflicts, collapse = ", "))}
list2env(as.list(environment()), envir = .GlobalEnv); save.image(paste0(output_dir,"/QC_and_others/globalenvir.RData"))

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
# print("ALL DONE")

# print(paste0("Current date: ",date()))
