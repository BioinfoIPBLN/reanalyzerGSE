#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GEO_ID <- args[1]
path <- args[2]
genes <- args[3]
# GEO_ID="GSE38010"
# path="/mnt/lustre/scratch/nlsas/home/csic/bds/caa/out_dir/database_GSE"
# genes="Cort,Zfp990"

path=paste0(path,"/",GEO_ID)
dir.create(path, showWarnings = FALSE); setwd(path)

### Download info from GEO and from bioconductor, install the libraries required for the precise array, load and organize the data...:
a <- GEOquery::getGEO(GEO_ID)
# a <- GEOquery::getGEO(GEO_ID,destdir = path)
# GEOquery::getGEOSuppFiles(GEO_ID, baseDir = path)
phenodata <- a[[1]]@phenoData@data[grep("characteristic",names(a[[1]]@phenoData@data),value = T)]
phenodata <- phenodata[,apply(phenodata,2,function(x){length(unique(x)) != 1})]
if (class(phenodata) == "character"){
	design=unlist(lapply(strsplit(phenodata," "),function(x){tail(x,1)}))
} else {
	design=gsub(" |/","_",unname(apply(phenodata,1,function(x){paste(gsub(".*: ","",x),collapse="_")})))
}

sample_names <- paste(a[[1]]@phenoData@data$title,
                      design,
                      a[[1]]@phenoData@data$geo_accession,
                      sep="_")
sample_names <- gsub(" |/","_",sample_names)
write.table(cbind(design,sample_names),
            file=paste0(path,"/samples_info.txt"),quote = F,row.names = F, col.names = T,sep = "\t")

id <- unlist(lapply(strsplit(grep("array=",a[[1]]@featureData@varMetadata$Description,val=T),"array="),function(x){gsub("&.*","",x[2])}))
if (!is.null(id)){
	b <- read.table("https://bioconductor.org/packages/release/data/annotation/src/contrib/PACKAGES",fill=T)
	print("Installing and loading microarray R packages")
	for (i in grep(id,b$V2,ignore.case=T,val=T))
	{
	print(i)
	if (!require(i, quietly = TRUE))
	    BiocManager::install(i,update = F, ask = FALSE)
	}
	for (i in grep(id,b$V2,ignore.case=T,val=T))
	{
	print(i)
	suppressMessages(library(i,character.only = TRUE,quiet = T,warn.conflicts = F))
	}

	print("Analyzing platform..."); print(id)
	path2=paste(path,dirname(grep("_RAW",list.files(path=path,recursive=T), val=T)),sep="/")
	print(paste0("I'm in ", path2)); setwd(path2); system("tar -xvf *.tar")
	GSEXXXXX_allfiles <- affy::list.celfiles(path=path2,full.names=T); setwd(path)
}

soft_file <- read.table(grep("soft",list.files(path=path,recursive=T,full.names=T),val=T),fill=T)
if (length(unlist(apply(soft_file,1,function(x){grep("agilent",x,ignore.case=T)}))) > 0){
	path2=paste(path,dirname(grep("_RAW",list.files(path=path, recursive=T), val=T)),sep="/")
	print(paste0("I'm in ", path2)); setwd(path2); system("tar -xvf *.tar")
	files=list.files(pattern = "^GSM", path = path2,full.names = T)
	project <- limma::read.maimages(files,source="agilent", green.only=TRUE)
	id="agilent"; print(id); setwd(path)
}


### Info from the study:
soft <- GEOquery::getGEO(filename=grep("soft",list.files(path=path,recursive=T,full.names=T),val=T))
soft_meta <- soft@header
# https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/
# soft <- GEOquery::getGEO(filename=grep("soft",list.files(path=path,recursive=T,full.names=T),val=T))
print("GEO information:")
# GEOquery::Meta(soft)$description
# soft_meta <- soft@header
soft_meta$organism
soft_meta$title
soft_meta$description
soft_meta$manufacturer
soft_meta$status
soft_meta$web_link
print("Last update:")
soft_meta$last_update_date
soft_meta$geo_accession
soft_meta$platform
soft_meta$type

if (length(grep("HumanHT-12 v4|HumanHT-12 v4.0|HT12 V4",soft_meta$description,val=T,ignore.case=T)) > 0){
	path2=paste(path,dirname(grep("_RAW",list.files(path=path, recursive=T), val=T)),sep="/")
	print(paste0("I'm in ", path2)); setwd(path2); system("tar -xvf *.tar"); system("pigz -p 4 -d *.gz")
	id="HT12V4"; print(id)
}

### Analyzing the data:
dir.create(paste0(path,"/final_results_reanalysis"),showWarnings = FALSE)
dir.create(paste0(path,"/final_results_reanalysis/QC_and_others"),showWarnings = FALSE)
dir.create(paste0(path,"/final_results_reanalysis/violin_stats"),showWarnings = FALSE)
write.table(c("PLEASE NOTE:","If any of the violin plot do not display the expected statistics based on the title, it's because the test wasn't possible due to the data distribution or similar reasons, and therefore it's not the most suitable. If it is, then looking at it manually is required"),
	    file=paste0(path,"/final_results_reanalysis/violin_stats/readme.txt"),quote = F,row.names = F, col.names = F,sep = "\n")

if (id!="agilent" & id!="HT12V4"){
	# Reading cel files
	GSEXXXXX_raw <- affy::read.affybatch(filenames=GSEXXXXX_allfiles,verbose=T,compress=T)
	# Normalizing
	GSEXXXXX <- affy::rma(GSEXXXXX_raw) # If the required libraries for the precise array were not installed previously, this function would also install them
	# Gathering expression values
	GSEXXXXX_1 <- affy::exprs(GSEXXXXX) # Eset object
	save(GSEXXXXX_1,file=paste0(path,"/final_results_reanalysis/QC_and_others/final_eset_object.RData"))

	matriz_GSEXXXXX<-as.data.frame(GSEXXXXX_1)
	# Obtaining the relation gene - probe
	symbols <- as.data.frame(annotate::getSYMBOL(rownames(matriz_GSEXXXXX), sort(grep(".db",grep(id,b$V2,ignore.case=T,val=T),val=T))[1]))
	colnames(symbols) <- "gene"
	matriz_GSEXXXXX$gene <- symbols$gene
	# Merging values for same genes
	GSEXXXXX_agg <- aggregate(. ~ gene, data = matriz_GSEXXXXX, mean,na.rm=T)
	rownames(GSEXXXXX_agg) <- stringr::str_to_title(GSEXXXXX_agg$gene)
	GSEXXXXX_agg <- GSEXXXXX_agg[setdiff(names(GSEXXXXX_agg), "gene")]
	colnames(GSEXXXXX_agg) <- gsub(".CEL.*","",colnames(GSEXXXXX_agg))
	# head(GSEXXXXX_agg)
} else if (id=="agilent"){
	# Background correction
	project.bgc <- limma::backgroundCorrect(project, method="normexp")
	# Normalize the data with the 'quantile' method for 1-color
	project.NormData <- limma::normalizeBetweenArrays(project.bgc,method="quantile")
	# Adding genenames & colapsing
	GSEXXXXX_2 <- project.NormData$E # Eset object
	save(GSEXXXXX_2,file=paste0(path,"/final_results_reanalysis/QC_and_others/final_eset_object.RData"))

	GSEXXXXX <- as.data.frame(GSEXXXXX_2)
	GSEXXXXX$Probes <- project.NormData$genes$GeneName
	GSEXXXXX_agg <- aggregate(. ~ Probes, data = GSEXXXXX, mean,na.rm=T)
	rownames(GSEXXXXX_agg) <- stringr::str_to_title(GSEXXXXX_agg$Probes)
	GSEXXXXX_agg <- GSEXXXXX_agg[setdiff(names(GSEXXXXX_agg),"Probes")]

	# a <- unlist(lapply(strsplit(colnames(GSEXXXXX_agg),"_"),function(x){paste0("GSM",gsub(".*GSM","",grep("GSM",x,val=T)))}))
	b <- read.table(paste0(path,"/samples_info.txt"),header=T)
	# b$GSM <- unlist(lapply(strsplit(b$sample_names,"_"),function(x){paste0("GSM",gsub(".*GSM","",grep("GSM",x,val=T)))}))
	# identical(a,b$GSM) # T

	colnames(GSEXXXXX_agg) <- b$sample_names
	# head(GSEXXXXX_agg)
} else if (id=="HT12V4"){
	path_script <- funr::get_script_path()
	print("Sourcing analyzeBead.R from..."); print(path_script)
	source(paste0(path_script,"/analyzeBead.R"))
	phenodata <- as.data.frame(data.table::fread(paste0(path,"/GEO_info/phenodata_extracted.txt"),sep="\t",head=F,fill=T))
	phenodata_possible_designs <- phenodata[,names(which(sort(unlist(lapply(apply(phenodata,2,table),function(x){length(x)}))) <= 10))]
	names <- gsub("__","_",gsub("__","_",gsub("[[:punct:]]| ","_",paste(phenodata$V1,apply(phenodata_possible_designs,1,function(x){paste(x,collapse="_")}),sep="_"))))
	names <- make.unique(names, sep = '_')

	bead.results <- beadAnalyze(idats = grep(".idat",list.files(path=path2,full.names=T),val=T),
		                    names = names,
		                    condition = read.table(paste0(path,"/GEO_info/design_possible_full_1.txt"))$V1,
		                    ref.condition = unique(read.table(paste0(path,"/GEO_info/design_possible_full_1.txt"))$V1)[1])
	GSEXXXXX <- bead.results$expr_df[!is.na(bead.results$expr_df$SYMBOL),-c((length(colnames(bead.results$expr_df))-3):length(colnames(bead.results$expr_df)))]
	GSEXXXXX_agg <- aggregate(. ~ SYMBOL, data = GSEXXXXX, mean,na.rm=T)
	rownames(GSEXXXXX_agg) <- stringr::str_to_title(GSEXXXXX_agg$SYMBOL)
	GSEXXXXX_agg <- GSEXXXXX_agg[setdiff(names(GSEXXXXX),"SYMBOL")]
	# head(GSEXXXXX_agg)
	GSEXXXXX_3 <- bead.results$eset # Eset object
	save(GSEXXXXX_3,file=paste0(path,"/final_results_reanalysis/QC_and_others/final_eset_object.RData"))
}

write.table(GSEXXXXX_agg,
	    file=paste0(path,"/final_results_reanalysis/Expr_matrix_marray.txt"),quote = F,row.names = T, col.names = T,sep = "\t")
cat("\n"); print("Normalized expression written")

### Figure of the expr of certain genes of interest:
## Introduce in the violin plot statistics, loop through the different designs to get different coloring and grouping... etc
for (i in unlist(strsplit(genes,","))){
  for (z in list.files(pattern = "design_possible_full", recursive = TRUE, full.names=T, path=paste0(path,"/GEO_info"))){
	  system(paste0("sed 's/^$/-/g' -i ",z)); print(i); print(z)
	  if (i %in% rownames(GSEXXXXX_agg)){
		  a <- GSEXXXXX_agg[rownames(GSEXXXXX_agg)==i,]
		  df <- data.frame(sample=names(a),
					  expr_marray=as.numeric(unname(a)),
					  condition=read.table(z,head=F,blank.lines.skip=FALSE)$V1)
		  df$"Expr_marray (log2 + .1)" <- round(log2(df$expr_marray + 0.1),2)
		  a <- data.table::fread(list.files(pattern = "series_matrix.txt.gz$", recursive = TRUE, full.names=T, path=path),fill=T)

		  suppressMessages(library(ggpubr,quiet = T,warn.conflicts = F))
		  p <- ggbarplot(df, x = "sample", y = "Expr_marray (log2 + .1)", color = "condition",
				add = "mean_se", label=T,lab.vjust = 4,
				position = ggplot2::position_dodge()) +
				theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
			labs(x="",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/",i,"_barplot_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)

		  df2 <- reshape::melt(df)
		  df2 <- unique(df2[df2$variable=="Expr_marray (log2 + .1)",c("condition","value","sample")])
		  suppressMessages(library(dplyr,quiet = T,warn.conflicts = F))
		  df2 <- df2 %>%
			  add_count(condition, name = "condition_n")
		  df2$sample2 <- as.numeric(factor(df2$sample))
		  list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(df2$condition), unique(df2$condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T) # Used ***** only to make sure is a separator that's not going to be used in the conditions name
		  p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
			geom_violin(trim = T) +
			geom_boxplot(width=0.1) +
			geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
			ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
			geom_point(data = dplyr::filter(df2, condition_n == 1)) +
			theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
			labs(x="",y="Expr_marray (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/",i,"_violin_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)

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
			annotate("text", x=df2$condition[1], y=max(df2$value) + 2.5, label=paste0("Samples_numbering:\n",paste0("Number_",1:length(df2$sample),": ",df2$sample,collapse="\n")), hjust=0, vjust=1) +
			theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8)) +
			labs(x="",y="Expr_marray (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		  ggsave(p, filename = paste0(path,"/final_results_reanalysis/",i,"_violin_",gsub(".txt","",basename(z)),"_label_samples.pdf"),width=30, height=30)
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
			labs(x="",y="Expr_marray (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
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
			labs(x="",y="Expr_marray (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
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
			labs(x="",y="Expr_marray (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
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
			labs(x="",y="Expr_marray (log2 + .1)",color="Condition", title=paste0(paste(a$V2[grep("Series_geo_accession",a$V1)],a$V2[grep("Series_title",a$V1)],sep=": "),". GENE SHOWN: ", i))
		   ggsave(p, filename = paste0(path,"/final_results_reanalysis/violin_stats/",i,"_violin_anova_",gsub(".txt","",basename(z)),".pdf"),width=30, height=30)
}
}
}

print("All done")
print("Genes highlighted in PDF plots are:")
print(genes)

### QC from GEO
# suppressMessages(library(GEOquery,quiet = T,warn.conflicts = F))
suppressMessages(library(limma,quiet = T,warn.conflicts = F))
suppressMessages(library(umap,quiet = T,warn.conflicts = F))
suppressMessages(library(maptools,quiet = T,warn.conflicts = F))

pdf(paste0(path,"/final_results_reanalysis/QC_and_others/QC_box_whisker_plot.pdf"))
par(mar=c(7,4,2,1))
boxplot(GSEXXXXX_agg, boxwex=0.7, notch=T, outline=FALSE, las=2)
dev.off()
pdf(paste0(path,"/final_results_reanalysis/QC_and_others/QC_expr_distribution.pdf"))
par(mar=c(4,4,2,1))
plotDensities(GSEXXXXX_agg, legend=F)
dev.off()
tryCatch({
pdf(paste0(path,"/final_results_reanalysis/QC_and_others/QC_mean_variance_trend.pdf"))
ex <- na.omit(GSEXXXXX_agg) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend")
dev.off()
}, error=function(e){cat("Error in the plotSA function, double check manually if it's of interest")})
pdf(paste0(path,"/final_results_reanalysis/QC_and_others/QC_umap_plot_multidimensional_scaling.pdf"))
ex <- GSEXXXXX_agg[!duplicated(GSEXXXXX_agg), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=3", xlab="", ylab="", pch=20, cex=1.5)
# point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
dev.off()


### Differential expression analyses:
# Based on the analyzeBead.R code for differntiale xpression analyses (limma auto.DE function) I think I could extend this code to automatically perform differntial expression analyses also for microarrays ins this script, coming from the normalized eset objects above.
# The problem is that from the possible_design_i.txt files that I created in another step, I would need to take pairs of conditions, filter the samples, do it independently... it's possible, but I don't think it's worth it for now. I could just do it manually for the conditions that I'm required and when I'm required, and that's it.
