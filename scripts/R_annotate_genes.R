#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
pattern_to_match <- args[2]
organism <- args[3]

suppressMessages(library("AnnotationHub",quiet = T,warn.conflicts = F))
suppressMessages(library("dplyr",quiet = T,warn.conflicts = F))

# Source shared ENSEMBL helper
script_dir <- dirname(sub("^--file=", "", commandArgs()[grep("--file=", commandArgs())]))
ensembl_helper <- file.path(script_dir, "R_ensembl_to_symbol.R")
if (file.exists(ensembl_helper)) source(ensembl_helper)

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
} else {
  print("For now the only organisms supported for further adding annotation and info to gene IDs are human and mouse. There are plans to add more in the future, including both the already available orgDB and support via annotationForge to build custom databases...")
  stop(paste0(organism," is currently not supported..."))
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

  ### Detect if gene IDs are ENSEMBL and handle accordingly
  rpkm_file <- list.files(path=path,recursive=T, include.dirs=T, full.names=T,pattern="^RPKM_counts_genes.txt")[1]
  gene_ids_raw <- read.delim(rpkm_file)$Gene_ID
  
  # Determine keytype and prepare gene IDs
  use_ensembl_keytype <- FALSE
  if (exists("is_ensembl_id") && is_ensembl_id(gene_ids_raw)) {
    use_ensembl_keytype <- TRUE
    print("Detected ENSEMBL gene IDs — using ENSEMBL keytype for annotation lookup...")
    # Strip version suffixes for AnnotationDbi lookup
    gene_ids_for_lookup <- strip_ensembl_version(gene_ids_raw)
    # Also try to get symbols via GTF for better downstream matching
    gtf_file <- Sys.getenv("ANNOTATION_FILE", unset = "")
    .ensembl_map_annot <- if (nchar(gtf_file) > 0 && file.exists(gtf_file)) {
      build_gtf_ensembl_map(gtf_file)
    } else NULL
  } else {
    gene_ids_for_lookup <- convert_ids(gene_ids_raw, mode)
  }
  
  ### Build a master table for all the genes with quantified normalized expression in the project:
  annot <- tryCatch({
    if (use_ensembl_keytype) {
      # Pre-filter to valid ENSEMBL keys in the org.db to avoid "none valid" error
      valid_ensembl_keys <- keys(eval(parse(text=orgDB)), keytype = "ENSEMBL")
      lookup_ids <- unique(gene_ids_for_lookup)
      lookup_ids <- lookup_ids[toupper(lookup_ids) %in% toupper(valid_ensembl_keys)]
      print(paste0("ENSEMBL IDs matching org.db: ", length(lookup_ids), " / ", length(unique(gene_ids_for_lookup))))
      if (length(lookup_ids) == 0) {
        print("No ENSEMBL IDs matched org.db keys. Skipping AnnotationDbi lookup.")
        NULL
      } else {
        suppressMessages(suppressWarnings(AnnotationDbi::select(eval(parse(text=orgDB)),
          keys = lookup_ids,
          columns = c("SYMBOL","ALIAS","GENENAME","GOALL"),
          keytype = 'ENSEMBL')))
      }
    } else {
      suppressMessages(suppressWarnings(AnnotationDbi::select(eval(parse(text=orgDB)),
        keys = gene_ids_for_lookup,
        columns = c("ALIAS","GENENAME","GOALL"),
        keytype = 'SYMBOL')))
    }
  }, error = function(e) {
    print(paste0("Warning: Annotation lookup failed: ", e$message))
    print("Continuing with empty annotation...")
    NULL
  })

  if (!is.null(annot) && nrow(annot) > 0) {
    # Determine the key column name (ENSEMBL or SYMBOL depending on keytype used)
    key_col <- if (use_ensembl_keytype) "ENSEMBL" else "SYMBOL"
    
    annot_summary <- as.data.frame(annot %>% 
      group_by(across(all_of(key_col))) %>% 
      summarise(across(everything(), list(~ paste(unique(.), collapse = ", ")), .names = "concatenated_{col}")))

    colnames(annot_summary) <- gsub("concatenated_","",colnames(annot_summary))
    
    # Ensure GOALL column exists before processing
    if ("GOALL" %in% colnames(annot_summary)) {
      annot_summary$GOALL_names <- unlist(lapply(strsplit(annot_summary$GOALL,", "),function(y){
        tryCatch(paste(GOfuncR::get_names(y)$go_name,collapse=" // "), error = function(e) NA_character_)
      }))
      annot_summary$GO_terms <- unlist(apply(annot_summary,1,function(x){
        suppressWarnings(paste(unname(mapply(paste,strsplit(x["GOALL"],", "),strsplit(x["GOALL_names"]," // "),sep="_")),collapse=" /// "))
      }))
    } else {
      annot_summary$GOALL_names <- NA_character_
      annot_summary$GO_terms <- NA_character_
    }
    
    # Select output columns
    if (use_ensembl_keytype) {
      # Include SYMBOL for ENSEMBL-keyed results
      out_cols <- intersect(c(key_col, "SYMBOL", "ALIAS", "GENENAME", "GO_terms"), colnames(annot_summary))
    } else {
      out_cols <- intersect(c(key_col, "ALIAS", "GENENAME", "GO_terms"), colnames(annot_summary))
    }
    annot_summary_final <- annot_summary[, out_cols, drop = FALSE]

    ### Take the files and annotate them
    for (file in grep("annotation",list.files(path=path,recursive=T, include.dirs=T, full.names=T,pattern=pattern_to_match),invert=T,val=T)){
      a <- as.data.frame(data.table::fread(file))
      if("Gene_ID" %in% colnames(a)){
        if (use_ensembl_keytype) {
          # Match by stripping ENSEMBL version from Gene_ID
          a$Gene_ID_stripped <- strip_ensembl_version(a$Gene_ID)
          b <- merge(a, annot_summary_final, by.x="Gene_ID_stripped", by.y=key_col, all.x=T)
          b$Gene_ID_stripped <- NULL  # Remove temp column
        } else {
          a$Gene_ID <- convert_ids(a$Gene_ID,mode)
          annot_summary_final[[key_col]] <- convert_ids(annot_summary_final[[key_col]],mode)
          b <- merge(a,annot_summary_final,by.x="Gene_ID",by.y=key_col,all.x=T)
        }
        write.table(b,file=gsub(".txt$","_annotation.txt",file),col.names = T,row.names = F,quote = F,sep="\t")	
        print(paste0("Annotating ",basename(file)," ..."))
      }
    }
  } else {
    print("No annotation data available. Skipping annotation step.")
  }
}
