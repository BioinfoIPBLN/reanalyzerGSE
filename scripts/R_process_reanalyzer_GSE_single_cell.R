#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GEO_ID <- args[1]
path <- args[2]
genes <- args[3]

### I think normally the matrix of counts are not normalized, and I try to load the non-normalized counts. Then, I normalize:
### I also load the annotations:
setwd(paste(path,GEO_ID,"GEO_info",GEO_ID,sep="/"))
potential_counts <- grep("rpkm|norm",grep("expr|matr",list.files(),val=T,ignore.case=T),val=T,ignore.case=T,invert=T)
potential_annotation <- list.files()[!(list.files() %in% potential_counts)]

a <- as.data.frame(data.table::fread(potential_counts))
rownames(a) <- a$V1
a <- a[,-1]
b <- as.data.frame(data.table::fread(potential_annotation))
print("Good annotation loaded:")
identical(b[,1],colnames(a)) # T
print("If FALSE, PLEASE DON'T TRUST THESE RESULTS AND DOUBLE-CHECK MANUALLY")

print("I'M TRYING TO DO STUFF AUTOMATICALLY, BUT THE ANNOTATION PROVIDED BY THE USERS IS PROBABLY VERY USER-DEPENDENT, SO THIS IS PROBABLY GOING TO GIVE ERRORS AND QUICK MANUAL WORK TO PROCESS THE ANNOTATION TABLE WOULD BE REQUIRED")

b$Condition_full <- paste(b$Condition,b$Lesion,b$Celltypes,sep="___")
b$Condition_full <- factor(b$Condition_full,levels=sort(unique(b$Condition_full)))

# Remove the ones that are only one:
tmp <- names(table(b$Condition_full))[table(b$Condition_full) > 1]
b_2 <- b[b$Condition_full %in% tmp,]
b_2$Condition_full <- as.character(b_2$Condition_full)
a_2 <- a[,colnames(a) %in% b_2$Detected]

colnames(a_2) <- b_2$Condition_full

### Normalization
options(httr_oauth_cache="n")
a_norm <- scRecover::normalization(as.matrix(a_2))

### Get rowsums for each... as if it were bulk just summing
final_bulk_normalized_counts <- data.frame(Gene_ID=stringr::str_to_title(rownames(a_norm))) # Getting first cap letter in Gene_ID
for (i in names(table(colnames(a_2)))){
  final_bulk_normalized_counts <- cbind(final_bulk_normalized_counts,
                                        data.frame(i=rowSums(a_norm[,grep(i,colnames(a_norm))])))
  # print(i)
  colnames(final_bulk_normalized_counts)[dim(final_bulk_normalized_counts)[2]] <- i
}
# dim(final_bulk_normalized_counts)[2]



dir.create(paste0(path,"/",GEO_ID,"/final_results_reanalysis"),showWarnings = FALSE)
write.table(final_bulk_normalized_counts,
            file=paste0(path,"/",GEO_ID,"/final_results_reanalysis/Norm_counts_genes.txt"),quote = F,row.names = F, col.names = T,sep = "\t")
write.table("Please note that these counts come from matrix of single cell counts, normalized and manually sumed per condition to resemble bulk counts and make possible the figures",
            file=paste0(path,"/",GEO_ID,"/final_results_reanalysis/readme.txt"),quote = F,row.names = F, col.names = F,sep = "\n")







### Figures of the expr of certain genes of interest:
gene_counts_rpkm <- final_bulk_normalized_counts # Not sure these are rpkm, but I use this naming for the sake of copying the code...
conditions <- list()
for (i in 1:(stringr::str_count(colnames(gene_counts_rpkm)[2],"___")+1)){
conditions <- c(conditions,list(unlist(lapply(strsplit(colnames(gene_counts_rpkm)[-1],"___"),function(x){x[i]}))))
}

dir.create(paste0(path,"/",GEO_ID,"/final_results_reanalysis/violin_stats"),showWarnings = FALSE)
path <- paste0(path,"/",GEO_ID)
write.table(c("PLEASE NOTE:","If any of the violin plot do not display the expected statistics based on the title, it's because the test wasn't possible due to the data distribution or similar reasons, and therefore it's not the most suitable. If it is, then looking at it manually is required"),
	    file=paste0(path,"/final_results_reanalysis/violin_stats/readme.txt"),quote = F,row.names = F, col.names = F,sep = "\n")

for (i in unlist(strsplit(genes,","))){
  for (z in 1:length(conditions)){
	  # system(paste0("sed 's/^$/-/g' -i ",z)); print(i); print(z)
	  if (i %in% gene_counts_rpkm$Gene_ID){
		  a <- gene_counts_rpkm[gene_counts_rpkm$Gene_ID==i,-grep("Gene_ID",colnames(gene_counts_rpkm))]
		  df <- data.frame(sample=names(a),
					  expr_RPKM=as.numeric(unname(a)),
					  condition=conditions[[z]])
		  df$"Expr_RPKM (log2 + .1)" <- round(log2(df$expr_RPKM + 0.1),2)
		  a <- data.table::fread(list.files(pattern = "series_matrix.txt.gz$", recursive = TRUE, full.names=T, path=path),fill=T)

		  suppressMessages(library(ggpubr,quiet = T,warn.conflicts = F))
		  p <- ggbarplot(df, x = "sample", y = "Expr_RPKM (log2 + .1)", color = "condition",
				add = "mean_se", label=T,lab.vjust = 4,
				position = ggplot2::position_dodge()) +
				theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
			labs(x="",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/",i,"_barplot_",paste0("custom_design_",z),".pdf"),width=30, height=30)

		  df2 <- reshape::melt(df)
		  df2 <- unique(df2[df2$variable=="Expr_RPKM (log2 + .1)",c("condition","value","sample")])
		  suppressMessages(library(dplyr,quiet = T,warn.conflicts = F))
		  df2 <- df2 %>%
			  add_count(condition, name = "condition_n")
  		  df2$sample2 <- as.numeric(factor(df2$sample))
		  list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(df2$condition), unique(df2$condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T) # Used ***** only to make sure is a separator that's not going to be used in the conditions name
		  p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
			geom_violin(trim=T) +
			geom_boxplot(width=0.1) +
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			geom_point(data = dplyr::filter(df2, condition_n == 1)) +
			theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
			labs(x="",y="Expr_RPKM (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/",i,"_violin_",paste0("custom_design_",z),".pdf"),width=30, height=30)

		  df2$i <- 1
		  labs_total <- aggregate(i~condition,df2,sum)
		  labs_total$value <- NA
		  labs_total$i <- paste0("Total points: ",labs_total$i)
		  p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
			geom_violin(trim = T) +
			geom_boxplot(width=0.1) +
			geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
			ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
			geom_point(data = dplyr::filter(df2, condition_n == 1)) +
			geom_text(data=labs_total,aes(x=condition,y= max(df2$value) + 0.1,label=i)) +
			annotate("text", x=df2$condition[length(df2$condition)], y=max(df2$value) + 2.5, label=paste0("Samples_numbering:\n",paste0("Number_",1:length(df2$sample),": ",df2$sample,collapse="\n")), hjust=0, vjust=1) +
			theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
			labs(x="",y="Expr_marray (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/",i,"_violin_",paste0("custom_design_",z),"_label_samples.pdf"),width=30, height=30)
		  write.table(paste0("Samples_numbering:\n",paste0("Number_",1:length(df2$sample),": ",df2$sample,collapse="\n")),
		              file=paste0(path,"/final_results_reanalysis/label_samples.txt"),quote = F,row.names = F, col.names = F,sep = "\n")


		  p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
			geom_violin(trim=T) +
			# geom_signif(comparisons = list(unique(df2$condition)),
			#             map_signif_level=T, color="black",linetype="dashed",text=2.5) +
			stat_compare_means(method = "t.test", comparisons = list_combinations) + # Add pairwise comparisons p-value
			# stat_compare_means(method = "t.test", label.y = 10) + # Add global p-value
			geom_boxplot(width=0.1) +
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			geom_point(data = dplyr::filter(df2, condition_n == 1)) +
			theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
			labs(x="",y="Expr_RPKM (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_ttest_",paste0("custom_design_",z),".pdf"),width=30, height=30)
		  p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
			geom_violin(trim=T) +
			# geom_signif(comparisons = list(unique(df2$condition)),
			#             map_signif_level=T, color="black",linetype="dashed",text=2.5) +
			stat_compare_means(method = "wilcox.test",comparisons = list_combinations) + # Add pairwise comparisons p-value
			# stat_compare_means(method = "wilcox.test",label.y = 10) + # Add global p-value
			geom_boxplot(width=0.1) +
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			geom_point(data = dplyr::filter(df2, condition_n == 1)) +
			theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
			labs(x="",y="Expr_RPKM (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_wilcoxtest_",paste0("custom_design_",z),".pdf"),width=30, height=30)
		  p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
			geom_violin(trim=T) +
			# geom_signif(comparisons = list(unique(df2$condition)),
			#             map_signif_level=T, color="black",linetype="dashed",text=2.5) +
			stat_compare_means(method = "kruskal.test",comparisons = list_combinations) + # Add pairwise comparisons p-value
			# stat_compare_means(method = "kruskal.test",label.y = 10) + # Add global p-value
			geom_boxplot(width=0.1) +
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			geom_point(data = dplyr::filter(df2, condition_n == 1)) +
			theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
			labs(x="",y="Expr_RPKM (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_kruskaltest_",paste0("custom_design_",z),".pdf"),width=30, height=30)
		  p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
			geom_violin(trim=T) +
			# geom_signif(comparisons = list(unique(df2$condition)),
			#             map_signif_level=T, color="black",linetype="dashed",text=2.5) +
			stat_compare_means(method = "anova",comparisons = list_combinations) + # Add pairwise comparisons p-value
			# stat_compare_means(method = "anova",label.y = 10) + # Add global p-value
			geom_boxplot(width=0.1) +
			geom_jitter(shape=16, position=position_jitter(0.2)) +
			geom_point(data = dplyr::filter(df2, condition_n == 1)) +
			theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
			labs(x="",y="Expr_RPKM (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		   ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_anova_",paste0("custom_design_",z),".pdf"),width=30, height=30)
}
}
}

print("All done")
print("Genes highlighted in PDF plots are:")
print(genes)
