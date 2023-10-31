#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
organism <- args[2]
cores <- args[3]
enrichment_databases <- args[4]
pattern_search <- args[5]

print("Attempting automatic gene ontology enrichment analyses by autoGO of the results...")
print(paste0("Current date: ",date()))
suppressMessages(library(parallel,quiet = T,warn.conflicts = F))
suppressMessages(library(autoGO,quiet = T,warn.conflicts = F))
# choose_database()
### 03/2023:
  # [1] "Achilles_fitness_decrease"                          "Achilles_fitness_increase"                         
  # [3] "Aging_Perturbations_from_GEO_down"                  "Aging_Perturbations_from_GEO_up"                   
  # [5] "Allen_Brain_Atlas_10x_scRNA_2021"                   "Allen_Brain_Atlas_down"                            
  # [7] "Allen_Brain_Atlas_up"                               "ARCHS4_Cell-lines"                                 
  # [9] "ARCHS4_IDG_Coexp"                                   "ARCHS4_Kinases_Coexp"                              
  # [11] "ARCHS4_TFs_Coexp"                                   "ARCHS4_Tissues"                                    
  # [13] "Azimuth_Cell_Types_2021"                            "BioCarta_2013"                                     
  # [15] "BioCarta_2015"                                      "BioCarta_2016"                                     
  # [17] "BioPlanet_2019"                                     "BioPlex_2017"                                      
  # [19] "Cancer_Cell_Line_Encyclopedia"                      "CCLE_Proteomics_2020"                              
  # [21] "CellMarker_Augmented_2021"                          "ChEA_2013"                                         
  # [23] "ChEA_2015"                                          "ChEA_2016"                                         
  # [25] "ChEA_2022"                                          "Chromosome_Location"                               
  # [27] "Chromosome_Location_hg19"                           "ClinVar_2019"                                      
  # [29] "CORUM"                                              "COVID-19_Related_Gene_Sets"                        
  # [31] "COVID-19_Related_Gene_Sets_2021"                    "Data_Acquisition_Method_Most_Popular_Genes"        
  # [33] "dbGaP"                                              "DepMap_WG_CRISPR_Screens_Broad_CellLines_2019"     
  # [35] "DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019"     "Descartes_Cell_Types_and_Tissue_2021"              
  # [37] "Diabetes_Perturbations_GEO_2022"                    "Disease_Perturbations_from_GEO_down"               
  # [39] "Disease_Perturbations_from_GEO_up"                  "Disease_Signatures_from_GEO_down_2014"             
  # [41] "Disease_Signatures_from_GEO_up_2014"                "DisGeNET"                                          
  # [43] "Drug_Perturbations_from_GEO_2014"                   "Drug_Perturbations_from_GEO_down"                  
  # [45] "Drug_Perturbations_from_GEO_up"                     "DrugMatrix"                                        
  # [47] "DSigDB"                                             "Elsevier_Pathway_Collection"                       
  # [49] "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"          "ENCODE_Histone_Modifications_2013"                 
  # [51] "ENCODE_Histone_Modifications_2015"                  "ENCODE_TF_ChIP-seq_2014"                           
  # [53] "ENCODE_TF_ChIP-seq_2015"                            "Enrichr_Libraries_Most_Popular_Genes"              
  # [55] "Enrichr_Submissions_TF-Gene_Coocurrence"            "Enrichr_Users_Contributed_Lists_2020"              
  # [57] "Epigenomics_Roadmap_HM_ChIP-seq"                    "ESCAPE"                                            
  # [59] "FANTOM6_lncRNA_KD_DEGs"                             "Gene_Perturbations_from_GEO_down"                  
  # [61] "Gene_Perturbations_from_GEO_up"                     "Genes_Associated_with_NIH_Grants"                  
  # [63] "GeneSigDB"                                          "Genome_Browser_PWMs"                               
  # [65] "GlyGen_Glycosylated_Proteins_2022"                  "GO_Biological_Process_2013"                        
  # [67] "GO_Biological_Process_2015"                         "GO_Biological_Process_2017"                        
  # [69] "GO_Biological_Process_2017b"                        "GO_Biological_Process_2018"                        
  # [71] "GO_Biological_Process_2021"                         "GO_Cellular_Component_2013"                        
  # [73] "GO_Cellular_Component_2015"                         "GO_Cellular_Component_2017"                        
  # [75] "GO_Cellular_Component_2017b"                        "GO_Cellular_Component_2018"                        
  # [77] "GO_Cellular_Component_2021"                         "GO_Molecular_Function_2013"                        
  # [79] "GO_Molecular_Function_2015"                         "GO_Molecular_Function_2017"                        
  # [81] "GO_Molecular_Function_2017b"                        "GO_Molecular_Function_2018"                        
  # [83] "GO_Molecular_Function_2021"                         "GTEx_Aging_Signatures_2021"                        
  # [85] "GTEx_Tissue_Expression_Down"                        "GTEx_Tissue_Expression_Up"                         
  # [87] "GWAS_Catalog_2019"                                  "HDSigDB_Human_2021"                                
  # [89] "HDSigDB_Mouse_2021"                                 "HMDB_Metabolites"                                  
  # [91] "HMS_LINCS_KinomeScan"                               "HomoloGene"                                        
  # [93] "HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression" "HuBMAP_ASCTplusB_augmented_2022"                   
  # [95] "Human_Gene_Atlas"                                   "Human_Phenotype_Ontology"                          
  # [97] "HumanCyc_2015"                                      "HumanCyc_2016"                                     
  # [99] "huMAP"                                              "IDG_Drug_Targets_2022"                             
  # [101] "InterPro_Domains_2019"                              "Jensen_COMPARTMENTS"                               
  # [103] "Jensen_DISEASES"                                    "Jensen_TISSUES"                                    
  # [105] "KEA_2013"                                           "KEA_2015"                                          
  # [107] "KEGG_2013"                                          "KEGG_2015"                                         
  # [109] "KEGG_2016"                                          "KEGG_2019_Human"                                   
  # [111] "KEGG_2019_Mouse"                                    "KEGG_2021_Human"                                   
  # [113] "Kinase_Perturbations_from_GEO_down"                 "Kinase_Perturbations_from_GEO_up"                  
  # [115] "KOMP2_Mouse_Phenotypes_2022"                        "L1000_Kinase_and_GPCR_Perturbations_down"          
  # [117] "L1000_Kinase_and_GPCR_Perturbations_up"             "Ligand_Perturbations_from_GEO_down"                
  # [119] "Ligand_Perturbations_from_GEO_up"                   "LINCS_L1000_Chem_Pert_Consensus_Sigs"              
  # [121] "LINCS_L1000_Chem_Pert_down"                         "LINCS_L1000_Chem_Pert_up"                          
  # [123] "LINCS_L1000_CRISPR_KO_Consensus_Sigs"               "LINCS_L1000_Ligand_Perturbations_down"             
  # [125] "LINCS_L1000_Ligand_Perturbations_up"                "lncHUB_lncRNA_Co-Expression"                       
  # [127] "MAGMA_Drugs_and_Diseases"                           "MCF7_Perturbations_from_GEO_down"                  
  # [129] "MCF7_Perturbations_from_GEO_up"                     "Metabolomics_Workbench_Metabolites_2022"           
  # [131] "MGI_Mammalian_Phenotype_2013"                       "MGI_Mammalian_Phenotype_2017"                      
  # [133] "MGI_Mammalian_Phenotype_Level_3"                    "MGI_Mammalian_Phenotype_Level_4"                   
  # [135] "MGI_Mammalian_Phenotype_Level_4_2019"               "MGI_Mammalian_Phenotype_Level_4_2021"              
  # [137] "Microbe_Perturbations_from_GEO_down"                "Microbe_Perturbations_from_GEO_up"                 
  # [139] "miRTarBase_2017"                                    "Mouse_Gene_Atlas"                                  
  # [141] "MSigDB_Computational"                               "MSigDB_Hallmark_2020"                              
  # [143] "MSigDB_Oncogenic_Signatures"                        "NCI-60_Cancer_Cell_Lines"                          
  # [145] "NCI-Nature_2015"                                    "NCI-Nature_2016"                                   
  # [147] "NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions"     "NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions"    
  # [149] "NIH_Funded_PIs_2017_Human_AutoRIF"                  "NIH_Funded_PIs_2017_Human_GeneRIF"                 
  # [151] "NURSA_Human_Endogenous_Complexome"                  "Old_CMAP_down"                                     
  # [153] "Old_CMAP_up"                                        "OMIM_Disease"                                      
  # [155] "OMIM_Expanded"                                      "Orphanet_Augmented_2021"                           
  # [157] "PanglaoDB_Augmented_2021"                           "Panther_2015"                                      
  # [159] "Panther_2016"                                       "Pfam_Domains_2019"                                 
  # [161] "Pfam_InterPro_Domains"                              "PFOCR_Pathways"                                    
  # [163] "PhenGenI_Association_2021"                          "PheWeb_2019"                                       
  # [165] "Phosphatase_Substrates_from_DEPOD"                  "PPI_Hub_Proteins"                                  
  # [167] "Proteomics_Drug_Atlas_2023"                         "ProteomicsDB_2020"                                 
  # [169] "Rare_Diseases_AutoRIF_ARCHS4_Predictions"           "Rare_Diseases_AutoRIF_Gene_Lists"                  
  # [171] "Rare_Diseases_GeneRIF_ARCHS4_Predictions"           "Rare_Diseases_GeneRIF_Gene_Lists"                  
  # [173] "Reactome_2013"                                      "Reactome_2015"                                     
  # [175] "Reactome_2016"                                      "Reactome_2022"                                     
  # [177] "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"  "RNAseq_Automatic_GEO_Signatures_Human_Down"        
  # [179] "RNAseq_Automatic_GEO_Signatures_Human_Up"           "RNAseq_Automatic_GEO_Signatures_Mouse_Down"        
  # [181] "RNAseq_Automatic_GEO_Signatures_Mouse_Up"           "Serine_Threonine_Kinome_Atlas_2023"                
  # [183] "SILAC_Phosphoproteomics"                            "SubCell_BarCode"                                   
  # [185] "SynGO_2022"                                         "SysMyo_Muscle_Gene_Sets"                           
  # [187] "Table_Mining_of_CRISPR_Studies"                     "Tabula_Muris"                                      
  # [189] "Tabula_Sapiens"                                     "TargetScan_microRNA"                               
  # [191] "TargetScan_microRNA_2017"                           "TF_Perturbations_Followed_by_Expression"           
  # [193] "TF-LOF_Expression_from_GEO"                         "TG_GATES_2020"                                     
  # [195] "Tissue_Protein_Expression_from_Human_Proteome_Map"  "Tissue_Protein_Expression_from_ProteomicsDB"       
  # [197] "Transcription_Factor_PPIs"                          "TRANSFAC_and_JASPAR_PWMs"                          
  # [199] "TRRUST_Transcription_Factors_2019"                  "UK_Biobank_GWAS_v1"                                
  # [201] "Virus_Perturbations_from_GEO_down"                  "Virus_Perturbations_from_GEO_up"                   
  # [203] "Virus-Host_PPI_P-HIPSTer_2020"                      "VirusMINT"                                         
  # [205] "WikiPathway_2021_Human"                             "WikiPathways_2013"                                 
  # [207] "WikiPathways_2015"                                  "WikiPathways_2016"                                 
  # [209] "WikiPathways_2019_Human"                            "WikiPathways_2019_Mouse"
if (!is.na(enrichment_databases)){
    if (length(enrichment_databases) > 0){
      enrichment_databases <- c("GO_Biological_Process_2021","GO_Molecular_Function_2021","GO_Cellular_Component_2021")}
}
enrichment_databases <- unlist(strsplit(enrichment_databases,","))
if (grepl("sapiens", organism, fixed=TRUE)){
  databases_autoGO <- unique(c(enrichment_databases,"WikiPathway_2021_Human","RNAseq_Automatic_GEO_Signatures_Human_Down","RNAseq_Automatic_GEO_Signatures_Human_Up","Reactome_2022","KEGG_2021_Human","HDSigDB_Human_2021"))
  databases_autoGO_print <- paste(databases_autoGO,collapse=",")
  print(paste0("The databases selected for autoGO are ",databases_autoGO_print,". Please double check autoGO::choose_database(), which has > 200 databases, in case you want to add any extra by using the pipeline argument..."))
} else if (grepl("musculus", organism, fixed=TRUE)){
  databases_autoGO <- unique(c(enrichment_databases,"WikiPathways_2019_Mouse","RNAseq_Automatic_GEO_Signatures_Mouse_Down","RNAseq_Automatic_GEO_Signatures_Mouse_Up","Reactome_2022","Mouse_Gene_Atlas","KEGG_2021_Mouse","HDSigDB_Mouse_2021"))
  databases_autoGO_print <- paste(databases_autoGO,collapse=",")
  print(paste0("The databases selected for autoGO are ",databases_autoGO_print,". Please double check autoGO::choose_database(), which has > 200 databases, in case you want to add any extra by using the pipeline argument..."))
} else {
  print(paste0("Your organism is ",organism,", and unfortunately the pipeline for automatic functional enrichment currently fully supports only mouse and human. We'll include non-model organisms eventually, but in the meantime please don't give up and double check if you can use autoGO manually with any of the rest of databases that may contain your organism from autoGO::choose_database(), which has > 200 databases"))
  stop("Exiting")
}

for (f in list.files(pattern = pattern_search,path=path,recursive=T)){
  try({
    a <- data.table::fread(f,head=T,fill=T)
    b <- a[a$FDR < 0.05,1]
    d <- a[a$FDR < 0.05 & a$logFC>0,1]
    e <- a[a$FDR < 0.05 & a$logFC<0,1]
    if (dim(b)[1]!=0){
      write.table(b,file=paste0(gsub(".txt","",f),"_fdr_05.txt"),col.names = F,row.names = F,quote = F,sep="\t")  
      write.table(b$Gene_ID,file=paste0(gsub(".txt","",f),"_fdr_05_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    }
    if (dim(d)[1]!=0){
      write.table(d,file=paste0(gsub(".txt","",f),"_fdr_05_logpos.txt"),col.names = F,row.names = F,quote = F,sep="\t")
      write.table(d$Gene_ID,file=paste0(gsub(".txt","",f),"_fdr_05_logpos_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    }
    if (dim(e)[1]!=0){
      write.table(e,file=paste0(gsub(".txt","",f),"_fdr_05_logneg.txt"),col.names = F,row.names = F,quote = F,sep="\t")
      write.table(e$Gene_ID,file=paste0(gsub(".txt","",f),"_fdr_05_logneg_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    }
    b <- a[a$PValue < 0.05,1]
    d <- a[a$PValue < 0.05 & a$logFC>0,1]
    e <- a[a$PValue < 0.05 & a$logFC<0,1]
    if (dim(b)[1]!=0){
      write.table(b,file=paste0(gsub(".txt","",f),"_pval_05.txt"),col.names = F,row.names = F,quote = F,sep="\t")  
      write.table(b$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05.txt_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")  
    }
    if (dim(d)[1]!=0){
      write.table(d,file=paste0(gsub(".txt","",f),"_pval_05_logpos.txt"),col.names = F,row.names = F,quote = F,sep="\t")
      write.table(d$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05_logpos_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    }
    if (dim(e)[1]!=0){
      write.table(e,file=paste0(gsub(".txt","",f),"_pval_05_logneg.txt"),col.names = F,row.names = F,quote = F,sep="\t")
      write.table(e$Gene_ID,file=paste0(gsub(".txt","",f),"_pval_05_logneg_Gene_IDs.txt"),col.names = F,row.names = F,quote = F,sep="\n")
    }
  }, silent = TRUE)
}

process_file <- function(file){
  file2=sub("\\..+$", "", basename(file))
  path2=paste0(dirname(file),"/",file2,"_funct_enrichment/")
  dir.create(path2, showWarnings = FALSE);setwd(path2);file.copy(file,paste0(path2,basename(file)))
  try({
    autoGO(read_gene_lists(gene_lists_path=path2,which_list="everything",from_autoGO=F,files_format=basename(file)),
           databases_autoGO)
  }, silent = TRUE)
  path3=paste0(path2,"/","enrichment_tables")
  if (dir.exists(path3) & length(list.files(path3)) > 0){
    setwd(path3)
    enrich_tables <- read_enrich_tables(
        enrich_table_path = path3,
        which_list = "everything",
        from_autoGO = F,
        files_format = ".tsv")
      try({
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
          }, silent = TRUE)
  }
  file.rename(path3,sub("//enrichment_tables","_autoGO",path3)); unlink(path2,recursive=T)
}

mclapply(
    mc.cores = cores,
    X = list.files(path = path, pattern = paste0(pattern_search,"|_Gene_IDs\\.txt$"),recursive=T,full=T),
    FUN = process_file
)

# Tidying...
setwd(path)
print("Tidying...")
removeEmptyDirs <- function(directory) {
  # List all directories
  dirs <- list.dirs(directory, recursive = TRUE)
  
  # Check each directory
  for (dir in dirs) {
    # If the directory is empty
    if (length(dir(dir)) < 3) {
      # Remove the directory
      invisible(unlink(dir, recursive = TRUE))
    }
  }
}
removeEmptyDirs(path)

#files=list.files(path=getwd(),pattern="Gene_IDs.txt$")
#invisible(file.remove(files))
