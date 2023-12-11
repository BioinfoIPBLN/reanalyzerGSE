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

if (exists("orgDB")){
  print(paste0("Loaded ",orgDB," annotation to add information to all list of genes in the current analyses..."))
  # To choose columns with the type of annotation:
  # columns(eval(parse(text=orgDB)))
  # All the symbols (take a lot of time)
  # unname(as.character(eval(parse(text=paste0(gsub(".db","",orgDB),"SYMBOL")))))

  ### Build a master table for all the genes with quantified normalized expression in the project:
  annot <- suppressMessages(suppressWarnings(AnnotationDbi::select(eval(parse(text=orgDB)),
    keys = read.delim(list.files(path=path,recursive=T, include.dirs=T, full.names=T,pattern="^RPKM_count")[1])$Gene_ID,
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
      b <- merge(a,annot_summary_final,by.x="Gene_ID",by.y="SYMBOL",all.x=T)
  	  write.table(b,file=gsub(".txt$","_annotation.txt",file),col.names = T,row.names = F,quote = F,sep="\t")	
      print(paste0("Annotating ",basename(file)," ..."))
    }
  }
} else {
  print("For now the only organism supported for annotating gene IDs are human and mouse. There are plans to add more in the future, including the already available orgDB or support via annotationForge to build custom databases...")
  print(paste0(organism," is currently not supported..."))
}