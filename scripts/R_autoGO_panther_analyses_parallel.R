#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
organism <- args[2]
cores <- args[3]
enrichment_databases <- args[4]
pattern_search <- args[5]
padjustmethod <- args[6]

print("Attempting automatic gene ontology enrichment analyses by autoGO of the results...")
print(paste0("Current date: ",date()))
suppressMessages(library(parallel,quiet = T,warn.conflicts = F))
suppressMessages(library(autoGO,quiet = T,warn.conflicts = F))
suppressMessages(library(rbioapi,quiet = T,warn.conflicts = F))

a <- rba_panther_info(what="organisms")
org_panther <- as.numeric(a$taxon_id[grep(gsub(" ","_",organism),gsub(" ","_",a$long_name))])
methods <- rba_panther_info(what = "datasets")$id

# choose_database() # Check out enrichr webpage to get updated list
### 04/2024:
#  [1] "Genome_Browser_PWMs"                               
#  [2] "TRANSFAC_and_JASPAR_PWMs"                          
#  [3] "Transcription_Factor_PPIs"                         
#  [4] "ChEA_2013"                                         
#  [5] "Drug_Perturbations_from_GEO_2014"                  
#  [6] "ENCODE_TF_ChIP-seq_2014"                           
#  [7] "BioCarta_2013"                                     
#  [8] "Reactome_2013"                                     
#  [9] "WikiPathways_2013"                                 
# [10] "Disease_Signatures_from_GEO_up_2014"               
# [11] "KEGG_2013"                                         
# [12] "TF-LOF_Expression_from_GEO"                        
# [13] "TargetScan_microRNA"                               
# [14] "PPI_Hub_Proteins"                                  
# [15] "GO_Molecular_Function_2015"                        
# [16] "GeneSigDB"                                         
# [17] "Chromosome_Location"                               
# [18] "Human_Gene_Atlas"                                  
# [19] "Mouse_Gene_Atlas"                                  
# [20] "GO_Cellular_Component_2015"                        
# [21] "GO_Biological_Process_2015"                        
# [22] "Human_Phenotype_Ontology"                          
# [23] "Epigenomics_Roadmap_HM_ChIP-seq"                   
# [24] "KEA_2013"                                          
# [25] "NURSA_Human_Endogenous_Complexome"                 
# [26] "CORUM"                                             
# [27] "SILAC_Phosphoproteomics"                           
# [28] "MGI_Mammalian_Phenotype_Level_3"                   
# [29] "MGI_Mammalian_Phenotype_Level_4"                   
# [30] "Old_CMAP_up"                                       
# [31] "Old_CMAP_down"                                     
# [32] "OMIM_Disease"                                      
# [33] "OMIM_Expanded"                                     
# [34] "VirusMINT"                                         
# [35] "MSigDB_Computational"                              
# [36] "MSigDB_Oncogenic_Signatures"                       
# [37] "Disease_Signatures_from_GEO_down_2014"             
# [38] "Virus_Perturbations_from_GEO_up"                   
# [39] "Virus_Perturbations_from_GEO_down"                 
# [40] "Cancer_Cell_Line_Encyclopedia"                     
# [41] "NCI-60_Cancer_Cell_Lines"                          
# [42] "Tissue_Protein_Expression_from_ProteomicsDB"       
# [43] "Tissue_Protein_Expression_from_Human_Proteome_Map" 
# [44] "HMDB_Metabolites"                                  
# [45] "Pfam_InterPro_Domains"                             
# [46] "GO_Biological_Process_2013"                        
# [47] "GO_Cellular_Component_2013"                        
# [48] "GO_Molecular_Function_2013"                        
# [49] "Allen_Brain_Atlas_up"                              
# [50] "ENCODE_TF_ChIP-seq_2015"                           
# [51] "ENCODE_Histone_Modifications_2015"                 
# [52] "Phosphatase_Substrates_from_DEPOD"                 
# [53] "Allen_Brain_Atlas_down"                            
# [54] "ENCODE_Histone_Modifications_2013"                 
# [55] "Achilles_fitness_increase"                         
# [56] "Achilles_fitness_decrease"                         
# [57] "MGI_Mammalian_Phenotype_2013"                      
# [58] "BioCarta_2015"                                     
# [59] "HumanCyc_2015"                                     
# [60] "KEGG_2015"                                         
# [61] "NCI-Nature_2015"                                   
# [62] "Panther_2015"                                      
# [63] "WikiPathways_2015"                                 
# [64] "Reactome_2015"                                     
# [65] "ESCAPE"                                            
# [66] "HomoloGene"                                        
# [67] "Disease_Perturbations_from_GEO_down"               
# [68] "Disease_Perturbations_from_GEO_up"                 
# [69] "Drug_Perturbations_from_GEO_down"                  
# [70] "Genes_Associated_with_NIH_Grants"                  
# [71] "Drug_Perturbations_from_GEO_up"                    
# [72] "KEA_2015"                                          
# [73] "Gene_Perturbations_from_GEO_up"                    
# [74] "Gene_Perturbations_from_GEO_down"                  
# [75] "ChEA_2015"                                         
# [76] "dbGaP"                                             
# [77] "LINCS_L1000_Chem_Pert_up"                          
# [78] "LINCS_L1000_Chem_Pert_down"                        
# [79] "GTEx_Tissue_Expression_Down"                       
# [80] "GTEx_Tissue_Expression_Up"                         
# [81] "Ligand_Perturbations_from_GEO_down"                
# [82] "Aging_Perturbations_from_GEO_down"                 
# [83] "Aging_Perturbations_from_GEO_up"                   
# [84] "Ligand_Perturbations_from_GEO_up"                  
# [85] "MCF7_Perturbations_from_GEO_down"                  
# [86] "MCF7_Perturbations_from_GEO_up"                    
# [87] "Microbe_Perturbations_from_GEO_down"               
# [88] "Microbe_Perturbations_from_GEO_up"                 
# [89] "LINCS_L1000_Ligand_Perturbations_down"             
# [90] "LINCS_L1000_Ligand_Perturbations_up"               
# [91] "L1000_Kinase_and_GPCR_Perturbations_down"          
# [92] "L1000_Kinase_and_GPCR_Perturbations_up"            
# [93] "Reactome_2016"                                     
# [94] "KEGG_2016"                                         
# [95] "WikiPathways_2016"                                 
# [96] "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"         
# [97] "Kinase_Perturbations_from_GEO_down"                
# [98] "Kinase_Perturbations_from_GEO_up"                  
# [99] "BioCarta_2016"                                     
#[100] "HumanCyc_2016"                                     
#[101] "NCI-Nature_2016"                                   
#[102] "Panther_2016"                                      
#[103] "DrugMatrix"                                        
#[104] "ChEA_2016"                                         
#[105] "huMAP"                                             
#[106] "Jensen_TISSUES"                                    
#[107] "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO" 
#[108] "MGI_Mammalian_Phenotype_2017"                      
#[109] "Jensen_COMPARTMENTS"                               
#[110] "Jensen_DISEASES"                                   
#[111] "BioPlex_2017"                                      
#[112] "GO_Cellular_Component_2017"                        
#[113] "GO_Molecular_Function_2017"                        
#[114] "GO_Biological_Process_2017"                        
#[115] "GO_Cellular_Component_2017b"                       
#[116] "GO_Molecular_Function_2017b"                       
#[117] "GO_Biological_Process_2017b"                       
#[118] "ARCHS4_Tissues"                                    
#[119] "ARCHS4_Cell-lines"                                 
#[120] "ARCHS4_IDG_Coexp"                                  
#[121] "ARCHS4_Kinases_Coexp"                              
#[122] "ARCHS4_TFs_Coexp"                                  
#[123] "SysMyo_Muscle_Gene_Sets"                           
#[124] "miRTarBase_2017"                                   
#[125] "TargetScan_microRNA_2017"                          
#[126] "Enrichr_Libraries_Most_Popular_Genes"              
#[127] "Enrichr_Submissions_TF-Gene_Coocurrence"           
#[128] "Data_Acquisition_Method_Most_Popular_Genes"        
#[129] "DSigDB"                                            
#[130] "GO_Biological_Process_2018"                        
#[131] "GO_Cellular_Component_2018"                        
#[132] "GO_Molecular_Function_2018"                        
#[133] "TF_Perturbations_Followed_by_Expression"           
#[134] "Chromosome_Location_hg19"                          
#[135] "NIH_Funded_PIs_2017_Human_GeneRIF"                 
#[136] "NIH_Funded_PIs_2017_Human_AutoRIF"                 
#[137] "Rare_Diseases_AutoRIF_ARCHS4_Predictions"          
#[138] "Rare_Diseases_GeneRIF_ARCHS4_Predictions"          
#[139] "NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions"    
#[140] "NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions"    
#[141] "Rare_Diseases_GeneRIF_Gene_Lists"                  
#[142] "Rare_Diseases_AutoRIF_Gene_Lists"                  
#[143] "SubCell_BarCode"                                   
#[144] "GWAS_Catalog_2019"                                 
#[145] "WikiPathways_2019_Human"                           
#[146] "WikiPathways_2019_Mouse"                           
#[147] "TRRUST_Transcription_Factors_2019"                 
#[148] "KEGG_2019_Human"                                   
#[149] "KEGG_2019_Mouse"                                   
#[150] "InterPro_Domains_2019"                             
#[151] "Pfam_Domains_2019"                                 
#[152] "DepMap_WG_CRISPR_Screens_Broad_CellLines_2019"     
#[153] "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019"    
#[154] "MGI_Mammalian_Phenotype_Level_4_2019"              
#[155] "UK_Biobank_GWAS_v1"                                
#[156] "BioPlanet_2019"                                    
#[157] "ClinVar_2019"                                      
#[158] "PheWeb_2019"                                       
#[159] "DisGeNET"                                          
#[160] "HMS_LINCS_KinomeScan"                              
#[161] "CCLE_Proteomics_2020"                              
#[162] "ProteomicsDB_2020"                                 
#[163] "lncHUB_lncRNA_Co-Expression"                       
#[164] "Virus-Host_PPI_P-HIPSTer_2020"                     
#[165] "Elsevier_Pathway_Collection"                       
#[166] "Table_Mining_of_CRISPR_Studies"                    
#[167] "COVID-19_Related_Gene_Sets"                        
#[168] "MSigDB_Hallmark_2020"                              
#[169] "Enrichr_Users_Contributed_Lists_2020"              
#[170] "TG_GATES_2020"                                     
#[171] "Allen_Brain_Atlas_10x_scRNA_2021"                  
#[172] "Descartes_Cell_Types_and_Tissue_2021"              
#[173] "KEGG_2021_Human"                                   
#[174] "WikiPathway_2021_Human"                            
#[175] "HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression"
#[176] "GO_Biological_Process_2021"                        
#[177] "GO_Cellular_Component_2021"                        
#[178] "GO_Molecular_Function_2021"                        
#[179] "MGI_Mammalian_Phenotype_Level_4_2021"              
#[180] "CellMarker_Augmented_2021"                         
#[181] "Orphanet_Augmented_2021"                           
#[182] "COVID-19_Related_Gene_Sets_2021"                   
#[183] "PanglaoDB_Augmented_2021"                          
#[184] "Azimuth_Cell_Types_2021"                           
#[185] "PhenGenI_Association_2021"                         
#[186] "RNAseq_Automatic_GEO_Signatures_Human_Down"        
#[187] "RNAseq_Automatic_GEO_Signatures_Human_Up"          
#[188] "RNAseq_Automatic_GEO_Signatures_Mouse_Down"        
#[189] "RNAseq_Automatic_GEO_Signatures_Mouse_Up"          
#[190] "GTEx_Aging_Signatures_2021"                        
#[191] "HDSigDB_Human_2021"                                
#[192] "HDSigDB_Mouse_2021"                                
#[193] "HuBMAP_ASCTplusB_augmented_2022"                   
#[194] "FANTOM6_lncRNA_KD_DEGs"                            
#[195] "MAGMA_Drugs_and_Diseases"                          
#[196] "PFOCR_Pathways"                                    
#[197] "Tabula_Sapiens"                                    
#[198] "ChEA_2022"                                         
#[199] "Diabetes_Perturbations_GEO_2022"                   
#[200] "LINCS_L1000_Chem_Pert_Consensus_Sigs"              
#[201] "LINCS_L1000_CRISPR_KO_Consensus_Sigs"              
#[202] "Tabula_Muris"                                      
#[203] "Reactome_2022"                                     
#[204] "SynGO_2022"                                        
#[205] "GlyGen_Glycosylated_Proteins_2022"                 
#[206] "IDG_Drug_Targets_2022"                             
#[207] "KOMP2_Mouse_Phenotypes_2022"                       
#[208] "Metabolomics_Workbench_Metabolites_2022"           
#[209] "Proteomics_Drug_Atlas_2023"                        
#[210] "The_Kinase_Library_2023"                           
#[211] "GTEx_Tissues_V8_2023"                              
#[212] "GO_Biological_Process_2023"                        
#[213] "GO_Cellular_Component_2023"                        
#[214] "GO_Molecular_Function_2023"                        
#[215] "PFOCR_Pathways_2023"                               
#[216] "GWAS_Catalog_2023"                                 
#[217] "GeDiPNet_2023"                                     
#[218] "MAGNET_2023"                                       
#[219] "Azimuth_2023"                                      
#[220] "Rummagene_kinases"                                 
#[221] "Rummagene_signatures"                              
#[222] "Rummagene_transcription_factors"                   
#[223] "MoTrPAC_2023"                                      
#[224] "WikiPathway_2023_Human"                            
#[225] "SynGO_2024"#

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

# Decide autoGO/enrichr databases:
if (!is.na(enrichment_databases)){
    if (length(enrichment_databases) > 0){
      enrichment_databases <- c(enrichment_databases,"GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023")
    } else {
    enrichment_databases <- c("GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023")
    }
} else {
    enrichment_databases <- c("GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023")
}

enrichment_databases <- unique(unlist(strsplit(enrichment_databases,",")))

if (grepl("sapiens", organism, fixed=TRUE)){
  databases_autoGO <- unique(c(enrichment_databases,"WikiPathway_2023_Human","RNAseq_Automatic_GEO_Signatures_Human_Down","RNAseq_Automatic_GEO_Signatures_Human_Up","Reactome_2022","KEGG_2021_Human","HDSigDB_Human_2021"))
  databases_autoGO_print <- paste(databases_autoGO,collapse=",")
  print(paste0("The databases selected for autoGO are ",databases_autoGO_print,". Please double check autoGO::choose_database(), which has > 200 databases, in case you want to add any extra by using the pipeline argument..."))
  suppressMessages(library("org.Hs.eg.db",quiet = T,warn.conflicts = F))
  mode <- check_naming(keys(org.Hs.eg.db, keytype = "SYMBOL"))
} else if (grepl("musculus", organism, fixed=TRUE)){
  databases_autoGO <- unique(c(enrichment_databases,"WikiPathways_2019_Mouse","RNAseq_Automatic_GEO_Signatures_Mouse_Down","RNAseq_Automatic_GEO_Signatures_Mouse_Up","Reactome_2022","Mouse_Gene_Atlas","KEGG_2021_Mouse","HDSigDB_Mouse_2021"))
  databases_autoGO_print <- paste(databases_autoGO,collapse=",")
  print(paste0("The databases selected for autoGO are ",databases_autoGO_print,". Please double check autoGO::choose_database(), which has > 200 databases, in case you want to add any extra by using the pipeline argument..."))
  suppressMessages(library("org.Mm.eg.db",quiet = T,warn.conflicts = F))
  mode <- check_naming(keys(org.Mm.eg.db, keytype = "SYMBOL"))
} else {
  print(paste0("Your organism is ",organism,", and unfortunately the pipeline for automatic functional enrichment currently fully supports only mouse and human. We'll include non-model organisms eventually, but in the meantime please don't give up and double check if you can use autoGO manually with any of the rest of databases that may contain your organism from autoGO::choose_database(), which has > 200 databases"))
  stop("Exiting")
}


# Write sets of genes of interest
key_files <- data.frame()
for (f in grep("05|_annotation|Gene_ID",list.files(pattern = pattern_search,path=path,recursive=T, full.names=T),invert=T,val=T)){  
  try({
    a <- data.table::fread(f,head=T,fill=T)
    b <- a[a$FDR < 0.05,1]
    d <- a[a$FDR < 0.05 & a$logFC>0,1]
    e <- a[a$FDR < 0.05 & a$logFC<0,1]
    if (dim(b)[1]!=0){
      print(paste0("Writing ",basename(paste0(gsub(".txt","",f),"_fdr_05.txt"))," and ",basename(paste0(gsub(".txt","",f),"_fdr_05_Gene_IDs.txt")),"... ",dim(b)[1]," genes"))
      write.table(b,file=paste0(gsub(".txt","",f),"_fdr_05.txt"),col.names = F,row.names = F,quote = F,sep="\t")  
      write.table(convert_ids(b$Gene_ID,mode),file=paste0(gsub(".txt","",f),"_fdr_05_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    }
    if (dim(d)[1]!=0){
      print(paste0("Writing ",basename(paste0(gsub(".txt","",f),"_fdr_05_logpos.txt"))," and ",basename(paste0(gsub(".txt","",f),"_fdr_05_logpos_Gene_IDs.txt")),"... ",dim(d)[1]," genes"))
      write.table(d,file=paste0(gsub(".txt","",f),"_fdr_05_logpos.txt"),col.names = F,row.names = F,quote = F,sep="\t")
      write.table(convert_ids(d$Gene_ID,mode),file=paste0(gsub(".txt","",f),"_fdr_05_logpos_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    }
    if (dim(e)[1]!=0){
      print(paste0("Writing ",basename(paste0(gsub(".txt","",f),"_fdr_05_logneg.txt"))," and ",basename(paste0(gsub(".txt","",f),"_fdr_05_logneg_Gene_IDs.txt")),"... ",dim(e)[1]," genes"))
      write.table(e,file=paste0(gsub(".txt","",f),"_fdr_05_logneg.txt"),col.names = F,row.names = F,quote = F,sep="\t")
      write.table(convert_ids(e$Gene_ID,mode),file=paste0(gsub(".txt","",f),"_fdr_05_logneg_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    }
    key_files <- rbind(key_files,data.frame(new_files=c(paste0(gsub(".txt","",f),"_fdr_05_Gene_IDs.txt"),paste0(gsub(".txt","",f),"_fdr_05_logpos_Gene_IDs.txt"),paste0(gsub(".txt","",f),"_fdr_05_logneg_Gene_IDs.txt")),
                                            old_file=rep(f,3)))
    #b <- a[a$PValue < 0.05,1]
    #d <- a[a$PValue < 0.05 & a$logFC>0,1]
    #e <- a[a$PValue < 0.05 & a$logFC<0,1]
    #if (dim(b)[1]!=0){
      #write.table(b,file=paste0(gsub(".txt","",f),"_pval_05.txt"),col.names = F,row.names = F,quote = F,sep="\t")  
      #write.table(b$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")  
    #}
    #if (dim(d)[1]!=0){
      #write.table(d,file=paste0(gsub(".txt","",f),"_pval_05_logpos.txt"),col.names = F,row.names = F,quote = F,sep="\t")
      #write.table(d$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05_logpos_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    #}
    #if (dim(e)[1]!=0){
      #write.table(e,file=paste0(gsub(".txt","",f),"_pval_05_logneg.txt"),col.names = F,row.names = F,quote = F,sep="\t")
      #write.table(e$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05_logneg_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    #}
  }, silent = TRUE)
}

# Main function to apply in parallel
process_file <- function(file){
  # Copy file to a subfolder:
  while(dev.cur() > 1) dev.off()
  file2=sub("\\..+$", "", basename(file))
  path2=paste0(dirname(file),"/",file2,"_funct_enrichment/")
  dir.create(path2, showWarnings = FALSE);setwd(path2);invisible(file.copy(file,paste0(path2,basename(file))))
  print(paste0("Processing autoGO for ",file2," and ",length(read.table(file,head=F)$V1)," genes..."))
  
  # autoGO:
  tryCatch({
    autoGO(read_gene_lists(gene_lists_path=path2,which_list="everything",from_autoGO=F,files_format=basename(file)),
           databases_autoGO)
  }, error = function(e) {
          print("autoGO with errors"); print(e)
  })
  path3=paste0(path2,"/","enrichment_tables")
  if (dir.exists(path3) & length(list.files(path3)) > 0){
    setwd(path3)
    enrich_tables <- read_enrich_tables(
        enrich_table_path = path3,
        which_list = "everything",
        from_autoGO = F,
        files_format = ".tsv")
    tryCatch({
        for (i in 1:length(names(enrich_tables))){
          barplotGO(enrich_tables = enrich_tables[[i]],
                    title = c(names(enrich_tables)[i],file2),
                    outfolder = getwd(),
                    outfile = paste0(file2,"_",names(enrich_tables)[i],"_barplot.png"),
                    from_autoGO = F)
          lolliGO(enrich_tables = enrich_tables[[i]],
                  title = c(names(enrich_tables)[i],file2),
                  outfolder = getwd(),
                  outfile = paste0(file2,"_",names(enrich_tables)[i],"_lolliplot.png"),
                  from_autoGO = F)
          }
        }, error = function(e) {
          print("autoGO with errors"); print(e)
        })
  }
  invisible(file.rename(path3,sub("//enrichment_tables","_autoGO",path3)))

  # Panther:
  print(paste0("Processing Panther for ",file2," and ",length(read.table(file,head=F)$V1)," genes..."))
  setwd(path2)
  dataset <- read.table(paste0(file2,".txt"),head=F)$V1
  expr_back <- convert_ids(read.table(key_files$old_file[grep(file2,key_files$new_files)])$V1,mode)
  if (length(dataset) < 10000){
    for(annot_panther in methods){
      annot_panther2 <- annot_panther
      annot_panther2 <- gsub("ANNOT_TYPE_ID_PANTHER_","",gsub("GO:0003674","GO_MF",gsub("GO:0008150","GO_BP",gsub("GO:0005575","GO_CC",annot_panther2))))
      tryCatch({
        enriched_fisher_whole_back <- rba_panther_enrich(genes = dataset,
                                       organism = org_panther,
                                       annot_dataset = annot_panther,
                                       test_type = "FISHER",
                                       correction = padjustmethod,
                                       cutoff = 0.05)
        if (dim(enriched_fisher_whole_back$result)[1]!=0){
          write.table(enriched_fisher_whole_back$result,file=paste0("enriched_fisher_whole_back_",annot_panther2,"_",padjustmethod,"_panther.txt"),col.names = T,row.names = F,quote = F,sep="\t")      
        }
      }, error = function(e) {
          writeLines(as.character(e), paste0("enriched_fisher_whole_back_",annot_panther2,"_",padjustmethod,"_panther_err.txt"))
      })
      tryCatch({
        enriched_binom_whole_back <- rba_panther_enrich(genes = dataset,
                                       organism = org_panther,
                                       annot_dataset = annot_panther,
                                       test_type = "BINOMIAL",
                                       correction = padjustmethod,
                                       cutoff = 0.05)
        if (dim(enriched_binom_whole_back$result)[1]!=0){
          write.table(enriched_binom_whole_back$result,file=paste0("enriched_binom_whole_back_",annot_panther2,"_",padjustmethod,"_panther.txt"),col.names = T,row.names = F,quote = F,sep="\t")      
        }
      }, error = function(e) {
          writeLines(as.character(e),paste0("enriched_binom_whole_back_",annot_panther2,"_",padjustmethod,"_panther_err.txt"))
      })
      tryCatch({
        enriched_fisher_expr_back <- rba_panther_enrich(genes = dataset,
                                       organism = org_panther,
                                       annot_dataset = annot_panther,
                                       test_type = "FISHER",
                                       correction = padjustmethod,
                                       cutoff = 0.05,
                                       ref_genes = expr_back,
                                       ref_organism = org_panther)
        if (dim(enriched_fisher_expr_back$result)[1]!=0){
          write.table(enriched_fisher_expr_back$result,file=paste0("enriched_fisher_expr_back_",annot_panther2,"_",padjustmethod,"_panther.txt"),col.names = T,row.names = F,quote = F,sep="\t")      
        }
      }, error = function(e) {
          writeLines(as.character(e), paste0("enriched_fisher_expr_back_",annot_panther2,"_",padjustmethod,"_panther_err.txt"))
      })
      tryCatch({
        enriched_binom_expr_back <- rba_panther_enrich(genes = dataset,
                                       organism = org_panther,
                                       annot_dataset = annot_panther,
                                       test_type = "BINOMIAL",
                                       correction = padjustmethod,
                                       cutoff = 0.05,
                                       ref_genes = expr_back,
                                       ref_organism = org_panther)
        if (dim(enriched_binom_expr_back$result)[1]!=0){
          write.table(enriched_binom_expr_back$result,file=paste0("enriched_binom_expr_back_",annot_panther2,"_",padjustmethod,"_panther.txt"),col.names = T,row.names = F,quote = F,sep="\t")      
        }
      }, error = function(e) {
          writeLines(as.character(e), paste0("enriched_binom_expr_back_",annot_panther2,"_",padjustmethod,"_panther_err.txt"))
      })
    }
  } else {
    cat(file2, "... too many genes to apply PANTHER\n")
  } 
  setwd(path)
  invisible(file.rename(path2,sub("_funct_enrichment/","_funct_enrichment_panther",path2)))
}

final_files_list <- grep(tools::file_path_sans_ext(pattern_search),list.files(path = path, pattern = "_Gene_IDs\\.txt$",recursive=T,full=T),val=T)
mclapply(
    mc.cores = cores,
    X = final_files_list[which(!(duplicated(basename(final_files_list))))], # If you are at any time using this script outside the main pipeline be aware of these input files... there may be confounding...
    FUN = process_file
)

# Tidying...
#setwd(path)
#save.image("autoGO_panther_globalenvir.RData")
#print("Tidying...")
#removeEmptyDirs <- function(directory) {
  # List all directories
  #dirs <- list.dirs(directory, recursive = TRUE)
  #dirs <- dirs[grep("funct",dirs)]
  
  # Check each directory
  #for (fold in dirs) {
    # If the directory is empty
    #if (length(dir(fold)) < 3 && fold!=path) {
      # Remove the directory
      #invisible(unlink(fold, recursive = TRUE))
      #print(paste0("Removing ",fold))
    #}
  #}
#}
#removeEmptyDirs(path)
print("ALL DONE autoGO and Panther")
print(paste0("Current date: ",date()))


#files=list.files(path=getwd(),pattern="Gene_IDs.txt$")
#invisible(file.remove(files))
