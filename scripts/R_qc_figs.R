#!/usr/bin/env Rscript
load(file.path(commandArgs(trailingOnly=TRUE)[3],"QC_and_others/globalenvir.RData"))

args = commandArgs(trailingOnly=TRUE)
path <- args[1]
input_dir <- args[2]
output_dir <- args[3]
edgeR_object_prefilter <- eval(as.symbol(args[4]))
edgeR_object <- eval(as.symbol(args[5]))
edgeR_object_norm <- eval(as.symbol(args[6]))
pattern_to_remove <- args[7]
label <- basename(path)
label2 <- sub(".*_","",args[6])

cat(paste0("\nPerforming QC_PDF_",label,"_",label2,"...\n"));print(paste0("Current date: ",date()))

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
suppressMessages(library("clusterProfiler",quiet = T,warn.conflicts = F))
suppressMessages(library("png",quiet = T,warn.conflicts = F))
suppressMessages(library("curl",quiet = T,warn.conflicts = F))
suppressMessages(library("corrplot",quiet = T,warn.conflicts = F))
suppressMessages(library("ggpubr",quiet = T,warn.conflicts = F))
suppressMessages(library("ggpmisc",quiet = T,warn.conflicts = F))
suppressMessages(library("dplyr",quiet = T,warn.conflicts = F))
suppressMessages(library("ggdendro",quiet = T,warn.conflicts = F))

### Prepare counts:
  lcpm_prefilter <- cpm(edgeR_object_prefilter, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...
  if (pattern_to_remove!="none"){
    cat("\n\nRepeating QC figures removing the samples matching: ")
    print(pattern_to_remove)
    lcpm_prefilter <- cpm(edgeR_object_prefilter[,grep(pattern_to_remove,colnames(edgeR_object_prefilter),invert=T,val=T)],
                          log=TRUE)  
  }
  
  x <- edgeR_object
  if (pattern_to_remove!="none"){
    x <- edgeR_object[,grep(pattern_to_remove,colnames(edgeR_object),invert=T,val=T)]
  }
  cpm <- cpm(x) # This is normalized, altough not through edgeR, but the argument normalized.lib.sizes=TRUE by default in cpm...
  lcpm <- cpm(x, log=TRUE)  # This is log2 and normalized due to the argument normalized.lib.sizes=TRUE by default in cpm...
  lcpm_no_log <- cpm(x, log=FALSE)

  nsamples <- ncol(x)
  
  L <- mean(x$samples$lib.size) * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  # c(L, M)
  lcpm.cutoff <- log2(10/M + 2/L)

  x2 <- edgeR_object_norm
  if (pattern_to_remove!="none"){
    x2 <- edgeR_object_norm[,grep(pattern_to_remove,colnames(edgeR_object_norm),invert=T,val=T)]
  }
  #x2$samples$norm.factors <- 1
  lcpm2 <- cpm(x2, log=TRUE)
  lcpm2_no_log <- cpm(x2, log=FALSE)

  if (pattern_to_remove!="none"){
    targets <- targets[grep(pattern_to_remove,targets$Filename,invert=T),]
  }


### Prepare colors :
  group <- x$samples$group
  col.group <- as.factor(group)
  # levels(col.group) <- RColorBrewer::brewer.pal(nlevels(col.group), "Set1")
  # Not enough colors sometimes, so I get random colors:
  color = grDevices::colors()[grep('gr(a|e)y|white', grDevices::colors(), invert = TRUE)] # Get a list of non-gray or white colors
  contrast <- sapply(color,colorspace::contrast_ratio); contrast <- contrast[contrast>4] # Ensure a high contrast here and below (>4 on W3C standard)
  contrast2 <- unique(t(combn(unique(names(contrast)),2))[apply(t(combn(unique(names(contrast)),2)),1,function(x){colorspace::contrast_ratio(x[1],col2=x[2])}) > 4])
  levels(col.group) <- sample(contrast2, nlevels(col.group)); col.group <- as.character(col.group)

# summary(lcpm)
# table(rowSums(x$counts==0)==6)


### QC figures:
  pdf(paste0(output_dir,"/QC_and_others/",label,"_",label2,"_QC.pdf"),paper="A4")
  

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
  plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
  title(main="Filtered data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  legend("topright", legend=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), text.col=col, bty="n", cex = 0.5)
  
  ### 2.1. Boxplots non-normalised:
  par(mfrow=c(1,2))
  boxplot(lcpm, las=2, col=col.group, main="", names=targets$Name, cex.axis=0.4)
  title(main="Unnormalized data",ylab="Log-cpm",xlab="sample_type")
  
  ### 2.1. Boxplots normalised:
  boxplot(lcpm2, las=2, col=col.group, main="", names=targets$Name, cex.axis=0.4)
  title(main="Normalized data",ylab="Log-cpm",xlab="sample_type")
  
  ### 3. Library size figures:
  par(mfrow=c(1,1))
  bar_mids <- barplot(x$samples$lib.size,names.arg = gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name),las=2, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))))
  # Loop over the bar midpoints and add the text on top of each bar
  for(i in 1:length(bar_mids)) {
    # The y position is slightly above the top of the bar
    y_pos <- x$samples$lib.size[i] + 0.02 * max(x$samples$lib.size)    
    # Add the text, centered on the bar midpoint
    text(bar_mids[i], y_pos, labels = x$samples$lib.size[i], cex = 0.8, pos = 3)
  }
  
  if(length(grep("_skip_",list.files(path=input_dir,full.names=T,recursive=T)))==0){
    ### Figures with the number of reads
    reads <- c()
    files <- grep("_stats.txt",list.files(path=input_dir,full.names=T,recursive=T),val=T)
    for (f in files){reads <- c(reads,system(paste0("cat ",f," | grep '1st fragments' | sed 's,.*:\t,,g'"),intern=T))}
    bam_reads_2 <- data.frame(names=as.character(gsub("_nat.*","",basename(files))),reads=as.numeric(reads))
    # substr <- unlist(strsplit(bam_reads_2$names, "[\\W_]+"))
    # bam_reads_2$color <- substr[substr %in% names(table(substr)[table(substr)>1])]
    bam_reads_2$color <- bam_reads_2$names
  
    cat("Please check ordering and number of bam reads...\n"); print(bam_reads_2)  
  
    bar_plot <- ggbarplot(
      bam_reads_2, 
      x = "names", 
      y = "reads", 
      fill = "color",
      color = "color",
      stat = "identity"
    ) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    bar_plot + 
      geom_text(aes(label = reads), vjust = -0.5, color = "black", size=3) + labs(title="raw_reads") + guides(color = "none") + theme(legend.position = "none")
  
    reads <- c()
    files <- grep("_1_fastqc.html",list.files(path=input_dir,full.names=T,recursive=T),val=T)
    for (f in files){reads <- c(reads,system(paste0("cat ",f," | sed 's,<td>,\\n,g;s,</td>,\\n,g' | grep -A2 'Total Sequences' | tail -1"),intern=T))}
    
    if(length(files)!=0){ # Control that sometimes if these are repeated runs, fastqc is not going to be executed    
      fastq_reads_2 <- data.frame(names=as.character(gsub("_1_fastqc.*","",basename(files))),reads=as.numeric(reads))    
      # substr <- unlist(strsplit(fastq_reads_2$names, "[\\W_]+"))    
      # fastq_reads_2$color <- substr[substr %in% names(table(substr)[table(substr)>1])]
      fastq_reads_2$color <- fastq_reads_2$names
          
      cat("Please check ordering and number of fastq reads...\n"); print(fastq_reads_2)
      
      bar_plot <- ggbarplot(
        fastq_reads_2, 
        x = "names", 
        y = "reads", 
        fill = "color",
        color = "color",
        stat = "identity"
      ) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      perc <- c()
      for (i in 1:length(reads)){perc<-c(perc,round(bam_reads_2$reads[i]*100/fastq_reads_2$reads[i],2))}
      bar_plot + 
         geom_text(aes(label = paste0(reads," (bam/fastq: ",perc, " %)"), angle=45), vjust = -0.5, color = "black", size=2) + labs(title="bam_reads") + guides(color = "none") + theme(legend.position = "none")
    }
    reads_info <- fastq_reads_2[,1:2]
    reads_info$library_size <- x$samples$lib.size
    # Write the read numbers if there's need of rerunning or subsampling:
    write.table(reads_info,file=paste0(output_dir,"/QC_and_others/reads_numbers.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  }

  ### 4. Corrplot no log
  tmp <- lcpm_no_log; colnames(tmp) <- gsub("_t|m_Rep|_seq|_KO|_WT","",colnames(tmp))
  colnames(tmp) <- targets$Name[match(colnames(tmp),targets$Name)]
  corrplot(cor(tmp,method="spearman"), method='number',type = 'upper')
  corrplot(cor(tmp,method="spearman"), order='AOE')
  
  ### 5.1. MDS_norm
  z <- plotMDS(lcpm2_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise", plot=F)
  edge <- sd(z$x)
  plotMDS(lcpm2_no_log, labels=gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name), col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="MDS-PCoA Sample Names Norm")

  ### 5.2. MDS_log_norm
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
  plot(autoplot(data_pca.PC,title="PCA_over_edgeR",label=T,data=data_pca,colour='Type',xlim = c(-0.8,0.8),label.size=3,label.repel=T))
  
  ### 7. Heatmap 250 mots differential entities
  rsd <- rowSds(as.matrix(x))
  sel <- order(rsd, decreasing=TRUE)[1:250]
  samplenames <- gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name)
  heatmap(na.omit(as.matrix(x[sel,])),margins=c(10,8),main="Heatmap 250 most diff entities raw counts",cexRow=0.01,cexCol=0.5,labCol=samplenames)
  
  ### 8.1. Dendogram cluster raw norm
  par(mfrow=c(1,1), mar=c(8, 4, 4, 2),col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  pr.hc.c <- hclust(na.omit(dist(t(cpm(x2$counts,log=F)),method = "euclidean")))
  plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of normalized counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)
  
  ### 8.2. Dendogram cluster raw norm colored
  groups <- as.factor(x$samples$group)
  group_index <- as.numeric(groups[pr.hc.c$order])
  ggdendrogram(pr.hc.c, rotate = FALSE, theme_dendro=F) + theme_minimal() + theme(axis.text.x = element_text(face="bold", color=unique(col.group)[group_index])) + labs(title=paste("Hierarchical Clustering of normalized counts from samples of ", label, sep=""), y="Height",x="Samples") + ggpubr::rotate_x_text() + ylim(c(min(pr.hc.c$height),tail(sort(pr.hc.c$height),2)[1]))
  
  ### 8.3. Dendogram cluster raw norm log
  par(mfrow=c(1,1), mar=c(8, 4, 4, 2),col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  pr.hc.c <- hclust(na.omit(dist(t(cpm(x2$counts,log=T)),method = "euclidean")))
  plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 normalized counts from samples of ", label, sep=""), labels=targets$Filename, cex=0.5)
  
  ### 8.4. Dendogram cluster raw norm log colored
  groups <- as.factor(x$samples$group)
  group_index <- as.numeric(groups[pr.hc.c$order])
  ggdendrogram(pr.hc.c, rotate = FALSE, theme_dendro=F) + theme_minimal() + theme(axis.text.x = element_text(face="bold", color=unique(col.group)[group_index])) + labs(title=paste("Hierarchical Clustering of normalized counts from samples of ", label, sep=""), y="Height",x="Samples") + ggpubr::rotate_x_text() + ylim(c(min(pr.hc.c$height),tail(sort(pr.hc.c$height),2)[1]))

  dev.off()



