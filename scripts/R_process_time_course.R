#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
input <- args[2]
object <- args[3]
min.std <- as.numeric(args[4])
mestimate_value <- as.numeric(args[5])

suppressMessages(library(limma,quiet = T,warn.conflicts = F))
suppressMessages(library(Mfuzz,quiet = T,warn.conflicts = F))


### Function for performing Venn diagrams:
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
		suppressMessages(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
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
	try(system("tar cf venn_diagrams.tar Venn_diagram*; rm Venn_diagram*"))
}


###### Load read counts, format... etc:
setwd(paste0(path,"/DGE"))
tmpfile <- tempfile(); system(paste("tar -xzf", shQuote("allRData.tar.gz"), input, "-O >", tmpfile))
load(tmpfile); unlink(tmpfile); path <- dirname(getwd())

dir.create(paste0(path,"/time_course_analyses"),showWarnings = FALSE); setwd(paste0(path,"/time_course_analyses"))
eset <- eval(parse(text=object))

###### Limma:
	  print("Performing time course functional analyses with limma")
	  lev <- unique(gtools::mixedsort(targets$Type))
	  f <- factor(targets$Type, levels=lev)
	  design <- model.matrix(~0+f)
	  colnames(design) <- lev
	  fit <- lmFit(as.matrix(eset), design)

	  pairs_list <- split(combn(lev, 2), col(combn(lev, 2)))
	  for (i in 1:length(pairs_list)){
	    	print(paste0(paste(pairs_list[[i]],collapse="_vs_"),"..."))
		    contrast_formula <- paste(pairs_list[[i]][1], "-", pairs_list[[i]][2])
			cont <- makeContrasts(contrast_formula, levels=design)
			fit2 <- contrasts.fit(fit, cont)
			fit2 <- eBayes(fit2)
			df <- topTable(fit2, adjust="BH"); df$Gene_ID <- rownames(df)
			colnames(df)[c(4,5)] <- c("PValue","FDR")
			write.table(df[,c(7,3,1,2,4,5)],file=paste0("DGE_limma_timecourse_",pairs_list[[i]][1],"_vs_",pairs_list[[i]][2],".txt"),col.names = T,row.names = F,quote = F,sep="\t")
	  }

###### Create Venn Diagrams for the result of limma analyses in time-course: (for the Mufuzz does not make sense because there won't be shared genes amongst the clusters)
if(length(list.files(path=paste0(path,"/time_course_analyses"),full.names=T,pattern="^DGE_limma_timecourse_.*_vs_.*\\.txt$"))!=0){
	print("Attempting to perform Venn diagrams for time_course analyses...")
	Venn_funct(list.files(path=paste0(path,"/time_course_analyses"),full.names=T,pattern="^DGE_limma_timecourse_.*_vs_.*\\.txt$"))
}


###### Mfuzz:
	  # Media de las réplicas de los time-points:
	  groups_reps <- sapply(levels(eset$samples$group),function(x){rownames(eset$samples)[eset$samples$group==x]})
	  print("Detected replicates and groups are:"); print(groups_reps); print("Obtaining the mean expression per condition...");
	  if(class(groups_reps)[1]=="list"){
	  	df_mfuzz <- as.data.frame(lapply(groups_reps,function(x){rowMeans(as.matrix(eset)[,colnames(as.matrix(eset)) %in% as.character(x)])}))
	  } else if(class(groups_reps)[1]=="matrix"){
	  	df_mfuzz <- as.data.frame(apply(groups_reps,2,function(x){rowMeans(as.matrix(eset)[,colnames(as.matrix(eset)) %in% as.character(x)])}))
	  }

	  # I have to remove columns with 0 variance
	  print("Performing soft-clustering analyses with Mfuzz")
	  variances = apply(df_mfuzz,1,var)
	  if (length(which(variances==0))!=0){
		 df_mfuzz <- df_mfuzz[-unname(which(variances==0)),]
	  }
          
      # Convert my data to expressionSet object:
      require(Biobase)
      clust<-new("ExpressionSet", exprs=as.matrix(df_mfuzz))

      clust

      print("Filtering 25% NA"); clust.r <- filter.NA(clust, thres=0.25)

      #clust.f <- fill.NA(clust.r,mode="mean")
      clust.f <- clust.r

      # Filter based on std is advisable
      print(paste0("You can check the file 'mfuzz_std_plot.pdf' and specify a new threshold, but for now filtering based on a min std of ",min.std))
      
      pdf("mfuzz_std_plot.pdf")
      clust.std <- filter.std(clust.f,min.std=min.std)
      dev.off()

      clust.s <- standardise(clust.std)
      # Standardisation of the expression values of every gene/transcript/protein is carried out, so that the average expression value for each gene/transcript/protein is zero and the standard deviation of its expression profile is one.

      # Estimates:
      if(as.numeric(mestimate_value)==0){
      	mestimate_value <- mestimate(clust.s)
      } else {
      	mestimate_value <- as.numeric(mestimate_value)
      }
      print(paste0("You can provide another one, but but for now the predicted or provided fuzzifier m statistic was ",mestimate_value))

      # How to interpret this?. The recommend number of clusters is where we start to have repetitions and it ends going up the line.
      # An automatica approach to try and capture the number of clusters where repetitions starts. If no repetitions, then the max number of clusters by default in cselection (32) is used.

	  pdf("mfuzz_cselection_plot.pdf")
	  tmp  <- cselection(clust.s,m=mestimate_value,repeats=5,visu=TRUE)
	  dev.off()

	  colnames(tmp) <- gsub("c:","",colnames(tmp))
	  cluster_number <- as.numeric(names(which(unlist(lapply(apply(tmp,2,function(x){unique(x)}),function(y){length(y)!=1})))[1]))

	  if(is.na(cluster_number)){
	  	cluster_number_final <- 32
	  } else {
	  	cluster_number_final <- cluster_number
	  }

	  print(paste0("The number of clusters is ",cluster_number_final))

      # We may tune this a little bit... fuzziness is a very experimental value, so we can decrease the predicted one a little bit or play with it... for example Juan Tena in their publication got a 2.5 predicted and they used 1.7

      # Alternative estimates of the number of clusters according to the manual:
      #Dmin(clust.s,2.58,crange=seq(4,40,4),repeats=3,visu=TRUE)
      #Dmin(clust.s,1.7,crange=seq(4,40,4),repeats=3,visu=TRUE)
      # The minimum centroid distance can be used as cluster validity index. For an optimal cluster number,
      # we may see a ‘drop’ of minimum centroid distance wh plotted versus a range of cluster number and
      # a slower decrease of the minimum centroid distance for higher cluster number. However, it
      # should be used with care, as the determination remains difficult especially for short time series and
      # overlapping clusters. Alternatively, the function cselection can be used or functional enrichment
      # analysis (e.g. using Gene Ontology) can help to adjust the cluster number.

    
      # Final computation:
      cl <- mfuzz(clust.s,c=cluster_number_final,m=mestimate_value)
      
	  grid_dimensions <- function(n) {
	    root <- floor(sqrt(n))
	    
	    if (root * root == n) {
	      return(c(root, root))
	    } else {
	      return(c(root + 1, root + 1))
	    }
	  }	  
	  
	  grids <- grid_dimensions(cluster_number_final)
	  if(grids[1] >= 6){
	  	grids <- c(5,ceiling(cluster_number_final/5))
	  }
	  pdf("mfuzz_plots.pdf", width=8.27, height=11.69)
      mfuzz.plot2(clust.s,cl=cl,mfrow=grids,x11=FALSE,Xwidth=8,Xheight=11,time.labels=colnames(df_mfuzz),cex.axis=0.5)
      dev.off()

      # Yellow or green colored lines correspond to genes with low membership value; red and purple colored lines correspond to genes with high membership value

	  elements_per_cluster <- acore(clust.s,cl,min.acore=0.51)
	  names(elements_per_cluster) <- paste0("Cluster_n",1:length(elements_per_cluster))
      elements_per_cluster_name <- lapply(elements_per_cluster, function(x){x$NAME})

      elements_per_cluster_df <- do.call(rbind, lapply(names(elements_per_cluster), function(name) {
		  df <- elements_per_cluster[[name]]
		  if (nrow(df) > 0) {
		    df$slot_name <- name
		    return(df)
		  } else {
		    return(NULL)
		  }
		}))

	  # Filter out NULL elements (which may arise from empty dataframes)
	  elements_per_cluster_df <- elements_per_cluster_df[!sapply(elements_per_cluster_df, is.null), ]; colnames(elements_per_cluster_df)[1] <- "Gene_ID"
	  write.table(elements_per_cluster_df,file="mfuzz_elements_clusters.txt",quote = F,row.names = F, col.names = T,sep = "\t")
	  for (i in 1:length(elements_per_cluster_name)){
	  	write.table(elements_per_cluster_name[[i]],file=paste0("mfuzz_elements_cluster_",names(elements_per_cluster_name[i]),"_Gene_IDs.txt"),quote = F,row.names = F, col.names = F,sep = "\n")
	  }

	  cl_dims = unlist(lapply(elements_per_cluster, function(x){dim(x)[1]})); 
	  print(paste0(sum(cl_dims)," genes classified in clusters"))
      
      print("Soft clustering was performed, please check the files mfuzz_plots.pdf and mfuzz_elements_cluster.txt to determine if any further processing is advisable (e.g. manually merging clusters with very similar patterns, correlation or PCA distribution, in order to apply the preferred method of functional enrichment) ")

      O <- overlap(cl)

      pdf("mfuzz_clusters_PCA.pdf")
      Ptmp <- overlap.plot(cl,over=O,thres=0.05) # 4PCAs
	  dev.off()

      suppressMessages(library(corrplot,quiet = T,warn.conflicts = F))
      pdf("mfuzz_clusters_corrplot.pdf")
      corrplot(O, method="square",is.corr = FALSE)
      dev.off()
      pdf("mfuzz_clusters_corrplot_2.pdf")
      # To increase the scale:
      corrplot(O, method="square")
      dev.off()


save.image(file="time_course_envir.RData")
