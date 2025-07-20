#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
pattern_to_match <- args[2]
organism <- args[3]

suppressMessages(library("AnnotationHub",quiet = T,warn.conflicts = F))
suppressMessages(library("dplyr",quiet = T,warn.conflicts = F))

# Load annotation info:
organism_cp <- gsub("_"," ",organism)
if(organism_cp=="Homo sapiens"){
  organism_cp_react <- "human"  
  orgDB <- "org.Hs.eg.db"
  suppressMessages(library("org.Hs.eg.db",quiet = T,warn.conflicts = F))  
} else if(organism_cp=="Mus musculus"){
  organism_cp_react <- "mouse"
  orgDB <- "org.Mm.eg.db"
  suppressMessages(library("org.Mm.eg.db",quiet = T,warn.conflicts = F))
}


# Deduce the naming convention in the orgDB package:
check_naming <- function(names) {
  names <- grep("[[:punct:]]|orf|p43|p45",names,val=T,invert=T)  
  print(paste0("Looking at the annotation of ",length(names)," genes..."))
  all_upper <- grepl("^[A-Z0-9]+$", names)
  all_lower <- grepl("^[a-z0-9]+$", names)
  
  # Checks if the first letter is uppercase followed by lowercase letters or numbers, for names with at least one letter
  first_upper_rest_lower <- sapply(names, function(name) {
    if (grepl("[A-Za-z]", name)) { # Check if the name contains at least one letter
      # Extract the first letter and the rest of the string separately
      first_letter <- substr(name, regexpr("[A-Za-z]", name), regexpr("[A-Za-z]", name))
      rest <- substr(name, regexpr("[A-Za-z]", name) + 1, nchar(name))
      # Check if the first letter is uppercase and the rest of the string is lowercase or numeric (ignoring leading numbers)
      return(grepl("^[A-Z]$", first_letter) && grepl("^[a-z0-9]*$", rest))
    } else {
      return(TRUE) # If the name doesn't contain letters, it trivially satisfies the condition
    }
  })

  if(!all(all_upper) && !all(all_lower) && !all(first_upper_rest_lower)){
    pattern_results <- which.max(c(sum(all_upper),sum(all_lower),sum(first_upper_rest_lower)))
    num_genes <- c(sum(all_upper),sum(all_lower),sum(first_upper_rest_lower))[pattern_results]
    pattern_result_final <- c("all_upper","all_lower","first_upper_rest_lower")[pattern_results]
    print(paste0("Identified pattern for gene naming is ",pattern_result_final, ", accounting for ",num_genes," genes"))
  } else {
    pattern_result_final <- c("all_upper","all_lower","first_upper_rest_lower")[c(all(all_upper),all(all_lower),all(first_upper_rest_lower))]
    print(paste0("Identified pattern for gene naming is ",pattern_result_final, ", accounting for all annotated genes"))
  }
  
  return(pattern_result_final)

  
}
convert_ids <- function(ids,mode) {
  if(mode=="all_upper"){
    ids2 <- toupper(ids)
  } else if (mode=="all_lower") {
    ids2 <- tolower(ids)
  } else if (mode=="first_upper_rest_lower") {
    ids2 <- stringr::str_to_title(ids)
  }
  return(ids2)
}
mode <- check_naming(keys(eval(parse(text=orgDB)), keytype = "SYMBOL"))


if (exists("orgDB")){
  print(paste0("Loaded ",orgDB," annotation to add information to all list of genes in the current analyses..."))
  # To choose columns with the type of annotation:
  # columns(eval(parse(text=orgDB)))
  # All the symbols (take a lot of time)
  # unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL")))))

  ### Build a master table for all the genes with quantified normalized expression in the project:
  annot <- suppressMessages(suppressWarnings(AnnotationDbi::select(eval(parse(text=orgDB)),
    keys = convert_ids(read.delim(list.files(path=path,recursive=T, include.dirs=T, full.names=T,pattern="^RPKM_counts_genes.txt")[1])$Gene_ID,mode),
    columns = c("ALIAS","GENENAME","GOALL"),
    keytype = 'SYMBOL')))

  annot_summary <- as.data.frame(annot %>% 
    group_by(SYMBOL) %>% 
    summarise(across(everything(), list(~ paste(unique(.), collapse = ", ")), .names = "concatenated_{col}")))

  colnames(annot_summary) <- gsub("concatenated_","",colnames(annot_summary))
  annot_summary$GOALL_names <- unlist(lapply(strsplit(annot_summary$GOALL,", "),function(y){paste(GOfuncR::get_names(y)$go_name,collapse=" // ")}))
  a <- strsplit(annot_summary[,"GOALL"],", ")
  b <- strsplit(annot_summary[,"GOALL_names"]," // ")
  d <- strsplit(annot_summary[,"ONTOLOGYALL"],", ")
  annot_summary$GO_terms <- unlist(apply(annot_summary,1,function(x){suppressWarnings(paste(unname(mapply(paste,strsplit(x["GOALL"],", "),strsplit(x["GOALL_names"]," // "),sep="_")),collapse=" /// "))}))
  annot_summary_final <- annot_summary[,c(1,2,3,8)]

  ### Take the files and annotate them
  for (file in grep("annotation",list.files(path=path,recursive=T, include.dirs=T, full.names=T,pattern=pattern_to_match),invert=T,val=T)){
  	a <- as.data.frame(data.table::fread(file))
    if("Gene_ID" %in% colnames(a)){
      a$Gene_ID <- convert_ids(a$Gene_ID,mode)
      annot_summary_final$SYMBOL <- convert_ids(annot_summary_final$SYMBOL,mode)
      b <- merge(a,annot_summary_final,by.x="Gene_ID",by.y="SYMBOL",all.x=T)
  	  write.table(b,file=gsub(".txt$","_annotation.txt",file),col.names = T,row.names = F,quote = F,sep="\t")	
      print(paste0("Annotating ",basename(file)," ..."))
    }
  }
} else {
  print("For now the only organism supported for annotating gene IDs are human and mouse. There are plans to add more in the future, including both the already available orgDB and support via annotationForge to build custom databases...")
  print(paste0(organism," is currently not supported..."))
}
