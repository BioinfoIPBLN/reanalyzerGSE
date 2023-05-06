#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
GEO_ID <- args[1]
path <- args[2]

setwd(paste(path,GEO_ID,"GEO_info",sep="/")); dir.create(paste0(path,"/",GEO_ID,"/final_results_reanalysis"),showWarnings = FALSE)

GEO_ID_path <- GEO_ID

for (z in unlist(strsplit(GEO_ID_path,"_"))){
	if (length(list.files(pattern="pysradb"))==0){
		GEO_ID <- z
		tryCatch({
			geoquery_data <- GEOquery::getGEO(GEO_ID)
		}, error=function(e){
			gpl_id <- paste0("GPL",gsub('"','',noquote(system(paste0("zcat ",getwd(),"/",grep(GEO_ID,list.files(pattern='_series_matrix'),val=T), " | grep GPL | grep Series | sed 's,.*GPL,,g'"),intern=T))))
			geoquery_data <- GEOquery::getGEO(filename=grep(gpl_id,list.files(pattern="soft.gz"),val=T))
		})
		phenodata <- c()
		for (i in 1:length(geoquery_data)){
		  phenodata <- dplyr::bind_rows(phenodata,geoquery_data[[i]]@phenoData@data[grep("characteristic|relation.1|title",names(geoquery_data[[i]]@phenoData@data),value = T)])
		}
		phenodata <- phenodata[,apply(phenodata,2,function(x){length(unique(x)) != 1})]

		srx_ids <- gsub(".*term=","",as.character(unlist(phenodata[grep("relation",names(phenodata))])))
		write.table(srx_ids,
		            file=paste(path,GEO_ID_path,"GEO_info","srx_ids.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\n", append=TRUE)

		phenodata <- phenodata[,apply(phenodata,2,function(x){length(unique(x)) != dim(phenodata)[1]})]

		if (class(phenodata) == "character"){
			design <- gsub("[[:punct:]]| |treatment", "_",phenodata)
		} else if (class(phenodata) == "data.frame" && dim(phenodata)[2]>2){
			phenodata <- phenodata[order(rownames(phenodata)),]
			design <- gsub("[[:punct:]]| |treatment", "_",unname(apply(phenodata[,-grep("relation",names(phenodata))],1,function(x){paste(gsub(".*: ","",x),collapse="_")})))
		} else if (class(phenodata) == "data.frame" && dim(phenodata)[2]<=2){
			phenodata <- phenodata[order(rownames(phenodata)),]
			design <- gsub("[[:punct:]]| |treatment", "_",phenodata[,-grep("relation",names(phenodata))])
		}
		design <- gsub("____|___|__","_",gsub("^____|^___|^__|^_","",design))

		write.table(design,
								file=paste(path,GEO_ID_path,"GEO_info","design.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\n", append=TRUE)

	# Alternative code to iteratively deal with pysradb info... I comment here and write new code below to do outside the loop, because I already merged the info in previous scripts:
		#if (file.exists("sample_info_pysradb.txt")){
	#	if (length(list.files(pattern="pysradb"))!=0){
	#		a <- as.data.frame(data.table::fread(grep(GEO_ID,list.files(pattern="pysradb"),val=T)))
	#		a <- a[order(a[,names(which.min(nchar(grep("^GSM",a[1,],val=T))))[1]]),] # Ensure the order of the metadata is following the GSMXXXX ids and the SRR ids, because this is the order in the series_matrix tables from which I get more info
	#		write.table(a$run_accession,
	#		            file=paste(path,GEO_ID_path,"GEO_info","srr_ids.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\n", append=TRUE)
	#		if ("disease state" %in% colnames(a)){
	#			design <- gsub("____|___|__","_",gsub("[[:punct:]]| |treatment", "_",a$"disease state"))
	#		} else if ("treatment" %in% colnames(a)){
	#			design <- gsub("____|___|__","_",gsub("[[:punct:]]| |treatment", "_",a$treatment))
	#		}
	#		write.table(design,
	#		            file=paste(path,GEO_ID_path,"GEO_info","design.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\n", append=TRUE)
	#		sample_names <- paste(design,
	#		                  a[,names(which.min(nchar(grep("^GSM",a[1,],val=T))))[1]],
	#		                  sep="_")
	#		sample_names <- gsub("[[:punct:]]| ","_",sample_names) # Make sure no commas in the name and no instances of _1 or _2, important for miARma-seq later on
	#		sample_names <- gsub("_2","2",gsub("_1","1",sample_names))
	#		write.table(sample_names,
	#		            file=paste(path,GEO_ID_path,"GEO_info","sample_names.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\n", append=TRUE)
	#		cat("\nExtracting info from pysradb output...\n")
	#	} else {

		geoquery_data_df <- c()
		for (i in 1:length(geoquery_data)){
			geoquery_data_df <- dplyr::bind_rows(geoquery_data_df,geoquery_data[[i]]@phenoData@data)
		}
		geoquery_data_df <- geoquery_data_df[order(rownames(geoquery_data_df)),]
		sample_names_1 <- geoquery_data_df$title
		sample_names_3 <- geoquery_data_df$geo_accession
		sample_names <- paste(substr(sample_names_1,1,50),
				      design,
				      substr(sample_names_1,1,20),
					  rownames(geoquery_data_df),
				      sep="_")
		sample_names <- gsub("[[:punct:]]| ","_",sample_names) # Make sure no commas in the name and no instances of _1 or _2, important for miARma-seq later on
		sample_names <- gsub("_2","2",gsub("_1","1",sample_names))
		sample_names <- gsub("____|___|__","_",gsub("^____|^___|^__|^_","",sample_names))
		sample_names <- unlist(lapply(strsplit(sample_names,"_"),function(x){paste(unique(x),collapse="_")}))
		write.table(sample_names,
		            file=paste(path,GEO_ID_path,"GEO_info","sample_names.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\n", append=TRUE)
	  cat("\nAttempting to manually extract info from the series_matrix...\n")

		# Should be, but make sure phenodata is ordered when it's a data.frame by the GSM, because this is the order in the series_matrix tables from which I get more info
		tryCatch({
		if (!is.null(dim(phenodata)) & all(startsWith(rownames(phenodata),"GSM"))){
			phenodata <- phenodata[rownames(phenodata)[order(rownames(phenodata))],]
			cat("\nPhenodata_:\n")
			print(phenodata)
		}}, error=function(e){
			cat("\nPhenodata__:\n")
			print(phenodata)
		})

		# Write phenodata:
		write.table(phenodata,
			    file=paste(path,GEO_ID_path,"GEO_info","phenodata_extracted.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\t", append=TRUE)
		# Write possible designs:

		if (!is.null(dim(phenodata))){
			#phenodata_possible_designs <- phenodata[,names(which(sort(unlist(lapply(apply(phenodata,2,table),function(x){length(x)}))) <= 10))]
			phenodata_possible_designs <- phenodata[,names(which(apply(phenodata,2,function(x){length(table(x))}) <= 10 & apply(phenodata,2,function(x){length(table(x))}) > 1 & apply(phenodata,2,function(x){length(table(x))}) <= dim(phenodata)[1]))] # Keep in mind if over 10 conditions in large studies...
			# Add combinations of columns
			combinations <- c()
			for (i in 2:dim(phenodata_possible_designs)[2]){
				apply(combn(1:dim(phenodata_possible_designs)[2],i),2,function(x){
					combinations <<- c(combinations,paste(x,collapse="."))
				})
			}

			for (i in combinations){
				lapply(strsplit(i,".",fixed=T),function(x){
					phenodata_possible_designs <<- cbind(phenodata_possible_designs,apply(phenodata_possible_designs[,as.numeric(x)],1,function(y){
					gsub(" ","_",paste(gsub(".*: ","",y),collapse="."))
					}))
				})
			}

			num_des=length(grep("design_possible_full",list.files(),val=T))
			for (i in (num_des+1):(dim(phenodata_possible_designs)[2]+num_des)){
				idx=i
				if(length(grep("design_possible_full",list.files(),val=T))!=0){
					idx=i-num_des
				}
				design_possible <- gsub("[[:punct:]]| ","_",phenodata_possible_designs[,idx])
				write.table(unique(design_possible),
						file=paste0(path,"/",GEO_ID_path,"/GEO_info/design_possible",i,"_",GEO_ID,".txt"),quote = F,row.names = F, col.names = F,sep = "\n")
				write.table(design_possible,
						file=paste0(path,"/",GEO_ID_path,"/GEO_info/design_possible_full",i,"_",GEO_ID,".txt"),quote = F,row.names = F, col.names = F,sep = "\n")
				cat(paste0("\nDetected ",i,": "))
				print(table(design_possible)); cat("\n")
			}
		} else {
			i=length(grep("design_possible_full",list.files(),val=T))+1
			if (stringr::str_count(GEO_ID_path,"_")!=0){
				phenodata <- paste0(phenodata,"_",GEO_ID)
			}
			design_possible <- gsub("[[:punct:]]| ","_",names(table(phenodata)))
			write.table(design_possible,
						file=paste0(path,"/",GEO_ID_path,"/GEO_info/design_possible_",i,"_",GEO_ID,".txt"),quote = F,row.names = F, col.names = F,sep = "\n")
			write.table(paste0("__",gsub("[[:punct:]]| ","_",phenodata)),
		  				file=paste0(path,"/",GEO_ID_path,"/GEO_info/design_possible_full_",i,"_",GEO_ID,".txt"),quote = F,row.names = F, col.names = F,sep = "\n")
			cat(paste0("\nDetected ",i,": "))
			print(table(phenodata)); cat("\n")
		}
	# Copy info to final folder:
	system(paste0("for i in $(ls -d ",path,"/",GEO_ID_path,"/GEO_info/* | grep 'design_possible_full'); do echo -e '\n'$(basename $i) >> ",path,"/",GEO_ID_path,"/final_results_reanalysis/possible_designs_all.txt && cat $i | sort | uniq >> ",path,"/",GEO_ID_path,"/final_results_reanalysis/possible_designs_all.txt; done && sed -i '1{/^$/d}' ",path,"/",GEO_ID_path,"/final_results_reanalysis/possible_designs_all.txt"))
	system(paste0("paste ",path,"/",GEO_ID_path,"/GEO_info/sample_names.txt ",path,"/",GEO_ID_path,"/GEO_info/phenodata_extracted.txt | sed 's/ /_/g' > ",path,"/",GEO_ID_path,"/final_results_reanalysis/phenotypic_data_samples.txt"))
	}
}

# The code above is to deal iteratively with samples in case pysradb has not worked (non-RNA-seq studies for example). But ideally, I should be able to extract all that is required from the file from pysradb (in the final_results folder, I processed in previous scripts the pysradb files and created samples_info)
if (file.exists(paste0(path,"/",GEO_ID_path,"/final_results_reanalysis/sample_info.txt"))){
	# Get SRR ids
	info <- as.data.frame(data.table::fread(paste0(path,"/",GEO_ID_path,"/final_results_reanalysis/sample_info.txt")))
	write.table(unique(info$run_accession),
		    file=paste(path,GEO_ID_path,"GEO_info","srr_ids.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\n")
	# Get sample_names:
	cols_to_ignore <- unique(c("column_label",
                    grep("AWS",colnames(info),val=T),
                    grep("run_accession",grep("access",colnames(info),val=T),invert=T,val=T),
                    grep("run",grep("access",colnames(info),val=T),invert=T,val=T),
                    grep("url",colnames(info),val=T),
                    grep("instrument",colnames(info),val=T),
                    grep("library",colnames(info),val=T),
                    grep("organism",colnames(info),val=T),
                    grep("NCBI",colnames(info),val=T),
					grep("run_total",colnames(info),val=T),
                    "experiment_desc","experiment_title","study_title","total_size","total_spots","run_alias"
                    ))
	info_filt <- info[,!(colnames(info) %in% cols_to_ignore)]
	info_filt_2 <- info_filt[,apply(info_filt,2,function(x){length(table(x))!=1})]
	info_filt_2 <- as.data.frame(apply(info_filt_2,2,function(x){tolower(sub("([+-])|[[:punct:]]|RNA-Seq|Homo sapiens|Mus musculus","\\1",x))}))
	info_filt_2 <- as.data.frame(apply(info_filt_2,2,function(x){gsub(" ","-",x)}))
	info_filt_2 <- as.data.frame(apply(info_filt_2,2,function(x){gsub("--","",x)}))
	info_filt_2 <- as.data.frame(apply(info_filt_2,2,function(x){gsub("^gsm","GSM",x)}))
	info_filt_2 <- as.data.frame(apply(info_filt_2,2,function(x){gsub("^srr","SRR",x)}))
	info_filt_2 <- as.data.frame(apply(info_filt_2,2,function(x){gsub(")$","",x)}))
	sample_names <- apply(info_filt_2,1,function(x){paste(x,collapse="_")})
	sample_names <- unlist(lapply(strsplit(sample_names,"-|_"),function(x){paste(unique(x),collapse="_")}))
	sample_names <- substr(sample_names,1,150)
	sample_names <- sub("([+-])|[[:punct:]]","_",sample_names) # Make sure no commas in the name or spaces, and no instances of _1 or _2, important for miARma-seq later on
	sample_names <- gsub("_2","-2",gsub("_1","-1",sample_names))
	sample_names <- gsub(";","",sample_names)
	sample_names <- unique(gsub("__","-_",sample_names))
	
	write.table(sample_names,
		    file=paste(path,GEO_ID_path,"GEO_info","sample_names.txt",sep="/"),quote = F,row.names = F, col.names = F,sep = "\n")
	# Get possible_designs... columns in the sample info that were all the same were removed before, and now I remove the ones that are all different plus the GSM and the SRR because I don't want them in the designs
	# I also remove the ones that are just two, and one of them empty
	# And I replace the empty cells with "-"
	info_filt_2 <- unique(info_filt_2)
	info_filt_2 <- info_filt_2[,!(colnames(info_filt_2) %in% c("run_accession","experiment_alias"))]
	info_filt_2 <- info_filt_2[,unname(unlist(lapply(apply(info_filt_2,2,function(x){table(x)}),function(y){length(names(y)[names(y)!=""])!=1})))]
	info_filt_2[info_filt_2 == ""] <- "-"
	info_filt_design_final <- as.data.frame(info_filt_2[,apply(info_filt_2,2,function(x){length(table(x))!=dim(info)[1]})])
	# Get combinations between the columns if more than one and add to the possible designs:
	if(dim(info_filt_design_final)[2] > 1){
		info_filt_design_final <- cbind(info_filt_design_final,
										as.data.frame(lapply(combn(info_filt_design_final, 2, simplify=FALSE),function(x){apply(x,1,function(y){paste(y,collapse="_")})})))
	}

	# Write:
	for (i in 1:dim(info_filt_design_final)[2]){
		design_possible <- paste0("cond-",info_filt_design_final[,i])
		if (length(table(design_possible))!=dim(info_filt_design_final)[1]){
			write.table(unique(design_possible),
					file=paste0(path,"/",GEO_ID_path,"/GEO_info/design_possible_",i,"_",GEO_ID,".txt"),quote = F,row.names = F, col.names = F,sep = "\n")
			write.table(design_possible,
						file=paste0(path,"/",GEO_ID_path,"/GEO_info/design_possible_full_",i,"_",GEO_ID,".txt"),quote = F,row.names = F, col.names = F,sep = "\n")

			cat(paste0("\nDetected: ",i," "))
			print(table(design_possible)); cat("\n")
		}
	}

	# Copy info to final folder:
	system(paste0("for i in $(ls -d ",path,"/",GEO_ID_path,"/GEO_info/* | grep 'design_possible_full'); do echo -e '\n'$(basename $i) >> ",path,"/",GEO_ID_path,"/final_results_reanalysis/possible_designs_all.txt && cat $i | sort | uniq >> ",path,"/",GEO_ID_path,"/final_results_reanalysis/possible_designs_all.txt; done && sed -i '1{/^$/d}' ",path,"/",GEO_ID_path,"/final_results_reanalysis/possible_designs_all.txt"))
}
cat("\nExtracting info from GEO download done\n")
