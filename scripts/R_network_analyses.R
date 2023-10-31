#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
expr <- args[2]
pattern_search <- args[3]
organism <- as.numeric(args[4])


suppressMessages(library(WGCNA,quiet = T,warn.conflicts = F))
suppressMessages(library(tidyr,quiet = T,warn.conflicts = F))
suppressMessages(library(STRINGdb,quiet = T,warn.conflicts = F))
suppressMessages(library(ggplot2,quiet = T,warn.conflicts = F))


for (genes_interest in list.files(pattern = pattern_search,path=path)){
	## Computing WGCNA:
		cat(paste0("Processing WGCNA of ",genes_interest,"...\n"))
		setwd(path)
		a <- as.data.frame(data.table::fread(expr))
		b <- as.data.frame(data.table::fread(genes_interest))
		a <- a[a$Gene_ID %in% b$Gene_ID,]
		nexpr_mat <- as.matrix(a[,-grep("Gene_ID",colnames(a))])
		rownames(nexpr_mat) <- a$Gene_ID

		allowWGCNAThreads()
		new_path=paste0(path,"/WGCNA/"); system(paste("mkdir -p ", new_path,sep=""))
		label=sub("\\..*","",basename(genes_interest))
		minSize=30
		MEDissThres=0.25

		t_exprs=t(nexpr_mat)

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
		                          deepSplit = 2,              # 0 - 4. simplified control over how sensitive module detection should be to module splitting
		                          pamRespectsDendro = F,
		                          # detectCutHeight = 0.75,
		                          minModuleSize = minSize,
		                          
		                          # == Module Adjustments ==
		                          reassignThreshold = 0,
		                          mergeCutHeight = MEDissThres,
		                          
		                          # == TOM == Archive the run results in TOM file (saves time)
		                          # saveTOMs = T,
		                          # saveTOMFileBase = paste0(label,"_blockwiseTOM"),
		                          
		                          # == Output Options
		                          numericLabels = F,
		                          verbose = 3)
		sink()

		cor <- temp_cor

		## Plot module merge results
		# Calculate dissimilarity of module eigengenes
		MEList = moduleEigengenes(t_exprs, colors = netwk$unmergedColors)
		MEs = MEList$eigengenes
		colnames(MEs)=gsub("ME","",colnames(MEs))
		MEDiss = 1-cor(MEs);
		# Cluster module eigengenes
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
		# Reorder modules so similar modules are next to each other
		module_order = gsub("ME","",names(MEs0))

		# Add treatment names
		MEs0$treatment = gsub("_Rep.*","",colnames(a)[-1])

		# tidy & plot data
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

		# Relate Module to Type Groups
		module_df = data.frame(
		  gene_id = names(netwk$colors),
		  clusters = netwk$colors)
		colnames(module_df)[1] <- "Gene_ID"
		write.table(module_df[order(module_df$clusters),], file = paste(new_path,"WGCNA_all_modules_",label,".txt",sep=""),sep="\t",row.names = F,col.names=T,quote=F)
		for (cluster in unique(module_df$clusters)){
		  module_df_tmp <- module_df[module_df$clusters==cluster,]
		  write.table(module_df_tmp[,1], file = paste(new_path,"WGCNA_modules_",label,"_cluster",cluster,"_Gene_IDs.txt",sep=""),sep="\n",row.names = F, col.names=F,quote=F)
		}

	  
	## Computing STRINGdb:    
	    cat(paste0("Processing STRINGdb with threshold FDR 0.05 of ",genes_interest,"...\n"))
	    setwd(path)
	    new_path=paste0(path,"/STRINGdb/"); system(paste("mkdir -p ", new_path,sep=""))
	    b <- as.data.frame(data.table::fread(genes_interest))
	    b <- b[b$FDR<0.05,]

	    if (length(b$Gene_ID)<2000){
		    label=sub("\\..*","",basename(genes_interest))
		    string_db <- STRINGdb$new( version="11.5", species=organism,score_threshold=200, network_type="full", input_directory="")
		    example1_mapped <- string_db$map( b, "Gene_ID", removeUnmappedRows = TRUE )
		    hits <- example1_mapped$STRING_id
		    pdf(file =  paste(new_path,"String_network_",label,".pdf",sep=""),paper="a4")
		    string_db$plot_network( hits )
		    dev.off()

		    # filter by p-value and add a color column
			# (i.e. green down-regulated gened and red for up-regulated genes)
			# example1_mapped_pval05 <- string_db$add_diff_exp_color( subset(example1_mapped, pvalue<0.05),logFcColStr="logFC" )
			col_fc <- grep("logFC",colnames(example1_mapped),val=T)
			example1_mapped_pval05 <- string_db$add_diff_exp_color( example1_mapped,logFcColStr=col_fc ) # No need to subset because this is going to use already FDR < 0.05
			# post payload information to the STRING server
			payload_id <- string_db$post_payload( example1_mapped_pval05$STRING_id,colors=example1_mapped_pval05$color )
			# display a STRING network png with the "halo"
			pdf(file =  paste(new_path,"String_network_colored",label,".pdf",sep=""),paper="a4")
			string_db$plot_network( hits, payload_id=payload_id )
			dev.off()
			
			enrichment <- string_db$get_enrichment( hits )
			colnames(enrichment)[7] <- "geneID"
			write.table(enrichment, file = paste0(new_path,"STRINGdb_functional_enrichment_",label,".txt"),sep="\t",row.names = F,col.names=T,quote=F)

			# string_db$set_background(backgroundV)

			# get clusters
			clustersList <- string_db$get_clusters(example1_mapped$STRING_id)
			# plot clusters		
			for(i in seq(1:length(clustersList))){
				pdf(file =  paste(new_path,"String_network_clusters",label,"_",i,".pdf",sep=""),paper="a4")
				string_db$plot_network(clustersList[[i]])
				dev.off()
			}
			# save_clusters:
			clustersList_2 <- lapply(clustersList,function(x){example1_mapped_pval05$Gene_ID[example1_mapped_pval05$STRING_id %in% x]})
			clustersList_3 <- do.call(rbind, lapply(1:length(clustersList_2), function(i){data.frame(Gene_ID = clustersList_2[[i]], clusters = i)}))
			write.table(clustersList_3, file = paste(new_path,"STRINGdb_all_modules_",label,".txt",sep=""),sep="\t",row.names = F,col.names=T,quote=F)
			for (cluster in unique(clustersList_3$clusters)){
			  module_df_tmp <- clustersList_3[clustersList_3$clusters==cluster,]
			  write.table(module_df_tmp[,1], file = paste(new_path,"STRINGdb_modules_",label,"_cluster",cluster,"_Gene_IDs.txt",sep=""),sep="\n",row.names = F, col.names=F,quote=F)
			}

		} else {
			print("STRINGdb cannot be executed for sets of more than 2000 genes... Skipping...")
		}
}
save.image(file=paste0(new_path,"network_analyses_envir.RData"))
