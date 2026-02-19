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

cat(paste0("\n\nPerforming QC_PDF_",label,"_",label2,"...\n"));print(paste0("Current date: ",date()))

suppressMessages(library("edgeR",quiet = T,warn.conflicts = F))
suppressMessages(library("RColorBrewer",quiet = T,warn.conflicts = F))
suppressMessages(library("genefilter",quiet = T,warn.conflicts = F))
suppressMessages(library("ggrepel",quiet = T,warn.conflicts = F))
suppressMessages(library("ggfortify",quiet = T,warn.conflicts = F))
suppressMessages(library("cluster",quiet = T,warn.conflicts = F))
suppressMessages(library("ggplot2",quiet = T,warn.conflicts = F))
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

  cat("\n\nSummary cpm log=TRUE per sample...")
  print(summary(lcpm))


### QC figures:
  pdf(paste0(output_dir,"/QC_and_others/",label,"_",label2,"_QC.pdf"),paper="A4")
  
  ### 0. Reminder of the samples:
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
  title(main="Unnormalized data (lcpm)",ylab="Log-cpm")
  
  ### 2.1. Boxplots normalised:
  boxplot(lcpm2, las=2, col=col.group, main="", names=targets$Name, cex.axis=0.4)
  title(main="Normalized data (lcpm2)",ylab="Log-cpm")
  
  ### 3. Library size and read counts figures:
  par(mfrow=c(1,1))
  bar_mids <- barplot(x$samples$lib.size,names.arg = gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Name),las=2, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))))
  # Loop over the bar midpoints and add the text on top of each bar
  for(i in 1:length(bar_mids)) {
    # The y position is slightly above the top of the bar
    y_pos <- x$samples$lib.size[i] + 0.02 * max(x$samples$lib.size)    
    # Add the text, centered on the bar midpoint
    text(bar_mids[i], y_pos, labels = x$samples$lib.size[i], cex = 0.8, pos = 3)
  }
  
  par(mfrow=c(1,1))
  bar_mids <- barplot(x$samples$lib.size,names.arg = gsub("_t|m_Rep|_seq|_KO|_WT","",targets$Filename),las=1, main="Library Size",col=col.group, ylim=range(pretty(c(0, x$samples$lib.size))))
  # Loop over the bar midpoints and add the text on top of each bar
  for(i in 1:length(bar_mids)) {
    # The y position is slightly above the top of the bar
    y_pos <- x$samples$lib.size[i] + 0.02 * max(x$samples$lib.size)    
    # Add the text, centered on the bar midpoint
    text(bar_mids[i], y_pos, labels = x$samples$lib.size[i], cex = 0.8, pos = 3)
  }
  
  
  
  bam_files_present <- length(list.files(path=input_dir, pattern="\\.bam$", recursive=TRUE)) > 0

  if(bam_files_present) {

  if(length(grep("_skip_",list.files(path=input_dir,full.names=T,recursive=T)))==0 & length(grep("_stats.txt",list.files(path=input_dir,full.names=T,recursive=T)))>0){
    ### Figures with the number of reads
    reads <- c()
    files <- grep("_flagstat.txt",list.files(path=input_dir,full.names=T,recursive=T),val=T)
    for (f in files){reads <- c(reads,sub(" .*","",read.delim(f)[9,]))}
    bam_reads_2 <- data.frame(names=sub("_hisat2.*|_STAR.*","",basename(files)),reads=as.numeric(reads))
    bam_reads_2$color <- col.group
  
    cat("Please check ordering and number of bam reads...\n"); print(bam_reads_2)  
  
    ## Barplot 1 no. of reads
    bar_plot <- ggbarplot(
      bam_reads_2, 
      x = "names", 
      y = "reads", 
      fill = "color",
      color = "color",
      stat = "identity"
    ) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    print(bar_plot + 
      geom_text(aes(label = reads), vjust = -0.5, color = "black", size=3) + labs(title="bam_read1") + guides(color = "none") + theme(legend.position = "none"))
  
  
    reads <- c()
    files <- grep("_1_fastqc.zip",list.files(path=input_dir,full.names=T,recursive=T),val=T)
    for (f in files){reads <- c(reads,read.delim(unz(f, file.path(sub(".zip","",(basename(f))),"fastqc_data.txt")))[6,2])}
    
    ## Barplot 2 no. of reads
    if(length(files)!=0){ # Control that sometimes if these are repeated runs, fastqc is not going to be executed    
      fastq_reads_2 <- data.frame(names=as.character(gsub("_1_fastqc.*","",basename(files))),reads=as.numeric(reads))    
      fastq_reads_2$color <- col.group
          
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
      print(bar_plot + 
        geom_text(aes(label = paste0(reads, "\n(bam/fastq: ", perc, " %)")), 
                  vjust = -0.3, color = "black", size = 2.5) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +  # 20% padding on top
        labs(title = "fastq_reads") + 
        guides(color = "none") + 
        theme(legend.position = "none"))
      
    }
    
    
    reads_info <- fastq_reads_2[,1:2]
    reads_info$library_size <- x$samples$lib.size
    # Write the read numbers if there's need of rerunning or subsampling:
    write.table(reads_info,file=paste0(output_dir,"/QC_and_others/reads_numbers.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  }

  ## Barplot 3 no. of alignments
  multiqc_dir <- file.path(input_dir,"../multiqc_out/multiqc_report_data/")
  aln_files <- grep("pe_plot.txt|star_alignment",list.files(multiqc_dir,full=T),val=T)
  if(length(aln_files) > 0) {
  aln <- read.table(aln_files[1], 
                    header = TRUE, sep = "\t", check.names = FALSE)
  
  aln_long <- aln %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "category", values_to = "count") %>%
    mutate(category = factor(category, levels = colnames(aln)[-1]))
 
   n_cats <- length(levels(aln_long$category))
   aln_colors <- setNames(
    colorRampPalette(c("#2166AC", "#74ADD1", "#ABD9E9", "#FEE090", "#FDAE61", "#F46D43", "#D9D9D9"))(n_cats),
    levels(aln_long$category)
  )
  
  # Total
  p1 <- ggplot(aln_long, aes(x = Sample, y = count / 1e6, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = ifelse(count > 0, round(count / 1e6, 2), "")), 
              position = position_stack(vjust = 0.5), size = 2.5, color = "black") +
    scale_fill_manual(values = aln_colors) +
    scale_y_continuous(labels = function(x) paste0(x, "M")) +
    labs(title = "Alignment categories total", x = "sample", y = "Count", fill = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color=col.group),
          legend.position = "right")
  
  # Proportions
  aln_prop <- aln_long %>%
    group_by(Sample) %>%
    mutate(prop = count / sum(count))
  
  p2 <- ggplot(aln_prop, aes(x = Sample, y = prop, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = ifelse(prop > 0.02, scales::percent(prop, accuracy = 1), "")), 
              position = position_stack(vjust = 0.5), size = 2.5, color = "black") +
    scale_fill_manual(values = aln_colors) +
    labs(title = "Alignment categories proportions", x = "sample", y = "Proportion", fill = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color=col.group),
          legend.position = "right")
  
  p1; p2
  } else {
    cat("\nNo alignment summary files found in MultiQC output, skipping alignment barplots...\n")
  }

  } else {
    cat("\nNo BAM files found in input_dir, skipping BAM-dependent QC plots (read counts, alignment barplots)...\n")
  }
  
  
  ### 4. Corrplot no log
  tmp <- lcpm_no_log; colnames(tmp) <- gsub("_hisat.*|_STAR.*","",colnames(tmp))
  corrplot(cor(tmp,method="spearman"), order='AOE',type = 'full',title="Spearman_correlation",tl.col = col.group,tl.srt = 45)
  corrplot(cor(tmp,method="spearman"), method='number',type = 'full', title="Spearman_correlation",tl.col = col.group,tl.srt = 45)
  corrplot(cor(tmp,method="pearson"), order='AOE',type = 'full',title="Pearson_correlation",tl.col = col.group,tl.srt = 45)
  corrplot(cor(tmp,method="pearson"), method='number',type = 'full', title="Pearson_correlation",tl.col = col.group,tl.srt = 45)
  
  ### 5.1. MDS_norm
  z <- plotMDS(lcpm2_no_log, labels=targets$Name, col=col.group, gene.selection = "pairwise", plot=F)
  edge <- sd(z$x)
  plotMDS(lcpm2_no_log, labels=targets$Name, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="MDS-PCoA Sample Names Norm")
  
  z <- plotMDS(lcpm2_no_log, labels=targets$Filename, col=col.group, gene.selection = "pairwise", plot=F)
  edge <- sd(z$x)
  plotMDS(lcpm2_no_log, labels=targets$Filename, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="MDS-PCoA Sample Names Norm")

  ### 5.2. MDS_log_norm
  z <- plotMDS(lcpm2, labels=targets$Name, col=col.group, gene.selection = "pairwise", plot=F)
  edge <- sd(z$x)
  plotMDS(lcpm2, labels=targets$Name, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="MDS-PCoA Sample Names log2 Norm")
  
  z <- plotMDS(lcpm2, labels=targets$Filename, col=col.group, gene.selection = "pairwise", plot=F)
  edge <- sd(z$x)
  plotMDS(lcpm2, labels=targets$Filename, col=col.group, gene.selection = "pairwise",xlim=c(min(z$x)-edge,max(z$x) + edge))
  title(main="MDS-PCoA Sample Names log2 Norm")
  
  ### 6. PCA por tipos
  data_pca <- as.matrix(x)
  data_pca <- as.data.frame(t(data_pca))
  rownames(data_pca) <- targets$Name
  data_pca.PC = prcomp(data_pca)
  data_pca$Type <- targets$Type
  data_pca$Filename <- targets$Filename
  data_pca$Name <- targets$Name
  autoplot(data_pca.PC,label=T,data=data_pca,colour=col.group,xlim = c(-0.8,0.8),label.size=3,label.repel=T) + theme_minimal() + ggtitle("PCA")
  
  rownames(data_pca) <- data_pca$Filename
  autoplot(data_pca.PC,label=T,data=data_pca,colour=col.group,xlim = c(-0.8,0.8),label.size=3,label.repel=T) + theme_minimal() + ggtitle("PCA")
  
  ### 7. Heatmap 250 mots differential entities
  rsd <- rowSds(as.matrix(x))
  sel <- order(rsd, decreasing=TRUE)[1:250]
  
  heatmap(na.omit(as.matrix(x[sel,])),margins=c(10,8),main="Heatmap 250 most diff entities raw counts",cexRow=0.01,cexCol=0.5,labCol=sub("_hisat2.*|_STAR.*","",rownames(x$samples)))
  
  ### 8.1. Dendogram cluster raw norm
  par(mfrow=c(1,1), mar=c(8, 4, 4, 2),col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  pr.hc.c <- hclust(na.omit(dist(t(cpm(x2$counts,log=F)),method = "euclidean")))
  pr.hc.c$labels <- sub("_hisat2.*|_STAR.*","",pr.hc.c$labels)
  plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of normalized counts from samples of ", label, sep=""), cex=0.5)
  
  ### 8.2. Dendogram cluster raw norm colored
  groups <- as.factor(x$samples$group)
  group_index <- as.numeric(groups[pr.hc.c$order])
  ggdendrogram(pr.hc.c, rotate = FALSE, theme_dendro=F) + theme_minimal() + theme(axis.text.x = element_text(face="bold", color=unique(col.group)[group_index])) + 
    labs(title=paste("Hierarchical Clustering of normalized counts from samples of ", label, sep=""), y="Height",x="Samples") + ggpubr::rotate_x_text() + ylim(c(min(pr.hc.c$height),tail(sort(pr.hc.c$height),2)[1]))
  
  ### 8.3. Dendogram cluster raw norm log
  par(mfrow=c(1,1), mar=c(8, 4, 4, 2),col.main="royalblue4", col.lab="royalblue4", col.axis="royalblue4", bg="white", fg="royalblue4", font=2, cex.axis=0.6, cex.main=0.8)
  pr.hc.c <- hclust(na.omit(dist(t(cpm(x2$counts,log=T)),method = "euclidean")))
  pr.hc.c$labels <- sub("_hisat2.*|_STAR.*","",pr.hc.c$labels)
  plot(pr.hc.c, xlab="Sample Distance",main=paste("Hierarchical Clustering of log2 normalized counts from samples of ", label, sep=""), cex=0.5)
  
  ### 8.4. Dendogram cluster raw norm log colored
  groups <- as.factor(x$samples$group)
  group_index <- as.numeric(groups[pr.hc.c$order])
  ggdendrogram(pr.hc.c, rotate = FALSE, theme_dendro=F) + theme_minimal() + theme(axis.text.x = element_text(face="bold", color=unique(col.group)[group_index])) +
    labs(title=paste("Hierarchical Clustering of normalized counts from samples of ", label, sep=""), y="Height",x="Samples") + ggpubr::rotate_x_text() + ylim(c(min(pr.hc.c$height),tail(sort(pr.hc.c$height),2)[1]))



  
dev.off()



