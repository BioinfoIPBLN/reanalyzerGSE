# reanalyzerGSE
reanalyzerGSE is a pipeline to assist with and streamline transcriptomic analyses of various datasets (i.e. microarrays, RNA-seq, scRNA-seq) by automatically reanalyzing raw data submitted to public databases like [GEO](https://www.ncbi.nlm.nih.gov/geo/), [ENA](https://www.ebi.ac.uk/ena/browser/home) or [SRA](https://www.ncbi.nlm.nih.gov/sra). Local data can also be provided by the user. The pipeline is based on several steps implementing standard tools and novel scripts. (i.e. data download, quality control, alignment to reference genome, quantification, differential gene expression analyses, functional enrichment analyses...)

## Installation
We suggest alternatives for installation. Please choose one of:

1) An Apptainer/Singularity container (~2.5 GB) is provided. You can either:

1.1) Use the .def file to create the .sif image by executing:
```
git clone https://github.com/BioinfoIPBLN/reanalyzerGSE
apptainer build reanalyzerGSE.sif reanalyzerGSE/external_software/installation/reanalyzerGSE.def | tee -a reanalyzerGSE.sif.build.log
```
1.2) Download the ready-to-use .sif image:
```
wget -q https://bit.ly/reana_apptainer -O reanalyzerGSE.sif
```


2) Another option is to use the folder 'external_software', which contain some of the required software (i.e. miARma-seq), and within the 'external_software/installation' folder a wrapper script installs and configures all dependencies (mainly through miniconda and pip, 'external_software/installation/install.sh'). To perform a conda-based installation and setup everything required to run reanalyzerGSE, please execute:
```
git clone https://github.com/BioinfoIPBLN/reanalyzerGSE
bash reanalyzerGSE/external_software/installation/install.sh 2>&1 | tee -a reanalyzerGSE/external_software/installation/install.sh.log # Check log to make sure that installation of all dependencies has been succesful
```

This should work if you already have miniconda3 installed, and also install miniconda3 within the reanalyzerGSE folder if not available or if you have kept it out of the PATH. Plese keep in mind that in the 'install.sh' script most of the versions of the tools installed by conda are frozen (by means of multiple '.yml' files corresponding to different environments), so please open an issue or try to install with conda if there are dependency-related problems or any software is not installed.


3) The less recommended option is to manually install the required software.
If you want to manually install the software, check out the list of required tools in the .def file (Apptainer/Singularity) or in the files '.yml' within the folder 'external_software/installation'. Please be aware that many scripts (bash, perl...) within the 'scripts' folder are also used, so you may need to manually change the interpreter in the corresponding statements (first line #!) to ensure that everything works in your system. You may also check out the 'install.sh' script to conform to other needs, such as making scripts executable ('chmod' command).


## Quick start / Minimal examples
Please go to test_data/README and follow instructions.


## reanalyzerGSE arguments
Please refer to the help page for futher details. Please be aware that a txt containing arguments can be provided to '-options'
```
reanalyzerGSE.sh -h
		-h | -help # Type this to get help
	        -options | Provide the file containing the parameters to be used. You can adapt the file 'manual_options.txt' provided in the scripts folder, alternative to manually input in the command line all the options...)
	        
	        #### Input/output: 
	        -i | -input # GEO_ID (GSEXXXXXX, separated by comma if more than one), or folder containing raw reads (please provide full absolute path, e.g. /path/folder_name/, containing only fastq.gz files and not folders, links or any other item), or almost any accession from ENA/SRA to download .fastq from (any of the ids with the prefixes PRJEB,PRJNA,PRJDB,ERP,DRP,SRP,SAMD,SAME,SAMN,ERS,DRS,SRS,ERX,DRX,SRX,ERR,DRR,SRR, please separated by commas if more than one id as input)
	        -iG | -input_GEO_reads # If you want to combine downloading metadata from GEO with reads from GEO or any database already downloaded, maybe from a previous attempt, please provide an absolute path
	        -n | -name # Name of the project/folder to create and store results
	        -o | -output_folder # Destination folder
	        -g | -genes # Genes to highlight their expression in plots (one or several, separated by comma and no space, none by default)
	        -TMP | -TMPDIR # Directory to export the environmental variable TMPDIR (by default or if left empty an internal folder of the output directory is used, or please enter 'system' to use system's default, or an absolute pathway that will be created if it does not exist)
	
	        #### Reference and databases:
	        -r | -reference_genome # Reference genome to be used (.fasta file or .gz, absolute pathway)
	        -ri | -reference_genome_index # If the reference genome to be used already has an index that would like to reuse, please provide full pathway here (by default the provided genome is indexed)
	        -a | -annotation # Reference annotation to be used (.gtf file, absolute pathway). If hisat2 is used, a gff file (make sure format is '.gff' and not '.gff3') is accepted (some QC steps like 'qualimap rnaseqqc' may be skipped though). You can provide a comma-separated list of the pathways to different annotation, and multiple/independent quantification/outputs from the same alignments will be generated.
	        -t | -transcripts # Reference transcripts to be used (.fasta cDNA file, absolute pathway, only used if '-s' argument not provided so salmon prediction of strandness is required)
	        -Dk | -kraken2_databases # Folder (absolute pathway) containing the database that should be used by Kraken2 (any input here, e.g. 'standard_eupathdb_48_kraken2_db', would activate the kraken2-based decontamination step)
	        -Ds | -sortmerna_databases # The database (absolute pathway) that should be used by SortMeRNA (any input here, e.g. '/path/to/rRNA_databases/smr_v4.3_sensitive_db.fasta', would activate the sortmerna-based rRNA removal step)
	        -Df | -databases_function # Manually provide a comma separated list of databases to be used in automatic functional enrichment analyses of DEGs (check out the R package autoGO::choose_database(), but the most popular GO terms are used by default)
	
	        #### Metadata and sample info:
	        -D | -design_custom_local # Specifying here the experimental design for the local dataset (by default an interactive prompt will ask for a comma-separated list of the same length than the number of samples, if you want to avoid that manual input please provide the list in this argument. If more than one design to provide, please input comma-separated list separed by a '/', without spaces. Please avoid naming that would match the same pattern in grep, e.g. XXXX and XXXX_TREAT)
	        -d | -design_custom # Manually specifying the experimental design for GEO download ('no' by default and if 'yes', please expect an interactive prompt after data download from GEO, and please enter the assignment to groups when asked in the terminal, with a comma-separated list of the same length than the number of samples)
	        -O | -organism # Specifying here the scientific name of the organism for the local dataset (by default an interactive prompt will ask for it, if you want to avoid that manual input please provide the full organism name in this argument, please use underline instead of space)
	        -Tx | -taxon_id # NCBI's taxon id of the organism
	        -R | -number_reads_to_subsample # Information and number of reads to subsample to the sequences before the analyses (none by default, a path to the 'reads_numbers.txt' file from a previous execution and a number of reads must be provided, separated with comma, and proportions will be computed, with all samples being scaled to approximately, +- 10% of that number)
	        -bv | -batch_vector # Comma-separated list of numbers for use as batch vector with Combat-seq
	        -bc | -batch_biological_covariable # Comma-separated list of numbers for use as batch vector of covariables of biological interest with Combat-seq
		-bf | -batch_format # Format of the provided batch variables ('num' for numeric/vector variables or 'fact' for factors, by default)
	        -C | -covariables # Please input a comma-separated list for the covariable that should be included in the limma model for removeBatchEffect or in the edgeR model for DGE (for now only one covariable allowed, for example an expected batch effect)
	        -Cf | -covariables_format # Format of the provided covariate ('num' by default for numeric covariables, or 'fact' for factors)
	        -T | -target # Protopical target file for attempts to differential gene expression analyses (containing filenames and covariates, automatically built if not provided)
	
	        #### Activate alternative modes:
	        -Dm | -debug_module # For debugging, step to remove the content of the corresponding folders and to resume a failed or incomplete run without repeating (one of 'step1', 'step1a', 'step1b', 'step1c', 'step2', 'step3a', 'step3b', 'step4', 'step5', 'step6', 'step7', 'step8', 'step9', or 'all' to execute everything, by default)
	        -q | -qc_raw_reads # Whether to perform quality control on the raw reads ('yes' by default, or 'no')
	        -fd | -full_differential_analyses # Whether to perform full differential enrichment analyses (for example including computation of DEGs or Venn diagrams, 'no' or 'yes', by default)
	        -fe | -functional_enrichment_analyses # Whether to perform functional enrichment analyses ('no' or 'yes', by default)
	        -cPf | -clusterProfiler_full # Whether to perform additional functional enrichment analyses with multiple databases using clusterProfiler, by default only ORA for GO BP, GO MF and GO CC, and KEGG and REACTOME enrichment, will be performed, as additional analyses may be slow if many significant DEGs or multiple number of comparisons ('yes' or 'no', by default)
	        -b | -batch # Batch effect present? (no by default, yes if correction through Combat-seq and model is to be performed, and info is going to be required in other arguments or prompts)
	        -B | -bed_mode # Whether to convert list of files to bed format so they can be visualized in genome browsers ('yes' or 'no', by default)
	        -Dc | -deconvolution # Whether to perform deconvolution of the bulk RNA-seq data by CDSeq ('yes', which may require few hours to complete, or 'no', by default)
	        -vv | -perform_volcano_venn # Whether to perform all Volcano plots and Venn diagrams, which may take a long time if many comparisons ('no' or 'yes', by default)
	        -aP | -aPEAR_execution # Whether to simplify pathway enrichment analysis results by detecting clusters of similar pathways and visualizing enrichment networks by aPEAR package, which may be slow ('yes' or 'no', by default)
	        -Ti | -tidy_tmp_files # Space-efficient run, with a last step removing raw reads if downloaded, converting bam to cram, removing tmp files... etc ('yes' or 'no', by default)
	        -Txls | -convert_tables_excel # Convert all tables in results from .txt format, without limitation of size to Excel's .xlsx format, with a limitation of 32,767 characters ('yes' or 'no', by default)
	        -Tc | -time_course # Whether to perform additional time-course analyses as a last step ('yes' or 'no', by default)
	        -Na | -network_analyses # Whether to perform network analyses, only for human or mouse analyses ('yes' or 'no', by default)
	        -apl | -auto_panther_log # Whether to perform additional autoGO and Panther analyses for DEGs separated by log2Fc positive or negative ('yes' or 'no', by default)
	
	        #### Processing parameters:
	        -s | -strand # Strandness of the library ('yes, 'no', 'reverse'). If not provided and '-t' used, this would be predicted by salmon. Please use this parameter if prediction not correct, see explanations in for example in bit.ly/strandness0 and bit.ly/strandness
	        -f | -filter # Threshold of gene counts to use ('bin' to capture the lower expressed genes, 'filterbyexpr' to use the edgeR solution, 'or 'standard', by default). Please provide a comma separated list with the filters to use at each quantification if multiple annotation are provided
	        -Of | -options_featureCounts_feat # The feature type to use to count in featureCounts (default 'exon')
	        -Os | -options_featureCounts_seq # The seqid type to use to count in featureCounts (default 'gene_name')
	        -A | -aligner # Aligner software to use ('hisat2' or 'star', by default)
	        -Des | -differential_expr_software # Software to be used in the differential expression analyses ('edgeR' by default, or 'DESeq2')
	        -fp | -fastp_mode # Whether to perform fastp analyses over the raw reads in default mode, except for adapter trimming and end trimming ('yes' or 'no', by default)
	        -fpa | -fastp_adapter # Whether to perform adapter trimming on the raw reads by fastp ('yes' or 'no', by default, to perform automatic trimming, or a path to a fasta file to perform trimming of its sequences)
	        -fpt | -fastp_trimming # Whether to trim the raw reads by fastp ('none' by default, if two numbers separated by comma, the indicated number of bases will be trimmed from the front and tail, respectively)
	        -std | -time_course_std # Standard deviation threshold to filter in time course analyses (numeric, 1 by default)
	        -fuzz | -time_course_fuzz # Fuziness value for the soft clustering approach (by default an estimate is automatically computed but manual testing is encouraged)
	
	        #### Filtering out samples/comparisons:
	        -G | -GSM_filter # GSM ids (one or several, separated by comma and no space) within the GSE entry to restrict the analysis to. An alternative to requesting a stop with -S to reorganize the downloaded files manually
	        -S | -stop # Manual stop so the automatically downloaded files can be manually modified ('yes' or 'no', by default)
	        -pR | -pattern_to_remove # A pattern to remove some matching samples from QC figures and DGE analyses (by default 'none')
	        -Dec | -differential_expr_comparisons # Whether to restrict the differential expression analyses to only some of the possible comparisons or reorder the 'treatment' and 'control' elements of the comparison ('no', by default, or a comma-separated list specifying separated by '\\' the elements in the comparison, which you could get from a preliminar previous run, e.g. 'A//B,C//D,D//A'...)
	
	        #### Functional enrichment/networking analyses
	        -cPm | -clusterProfiler_method # Method for adjusting p.value in clusterProfiler iterations (one of 'holm','hochberg','hommel','bonferroni','BH','BY,'none', or 'fdr', by default)
	        -cPu | -clusterProfiler_universe # Universe to use in clusterProfiler iterations ('all' or 'detected', by default)
	        -mGS | -clusterProfiler_minGSSize # minGSSize parameter to use in clusterProfiler iterations (a number, '10' by default)
	        -MGS | -clusterProfiler_maxGSSize # maxGSSize parameter to use in clusterProfiler iterations (a number, '500' by default)
	        -Pm | -panther_method # Method for adjusting p.value in panther analyses via rbioapi (one of 'NONE','BONFERRONI', or 'FDR', by default)
	        -rev | -revigo_threshold_similarity # Similarity threshold for Revigo summaries of GO terms (0-1, suggested values are 0.9, 0.7, 0.5, 0.4 for large, medium, small, and tiny levels of similarity, respectively, being default 0.7
	
	        #### Performance:
	        -p | -cores # Number of cores
	        -cR | -cores_reads_to_subsample # Cores to use in subsampling by seqtk (10 by default)
	        -pi | -cores_index # Number of cores for genome indexing in aligning step (by default, same than -p)
	        -M | -memory_max # Max RAM memory to be used by aligner or JAVA in bytes (by default 257698037760, or 240GB, used)
	        -P | -number_parallel # Number of files to be processed in parallel (10 by default)
	        -cG | -compression_level # Specify the compression level to gzip the downloaded fastq files from GEO (numeric '0' to '9', default '9')
	        -Ac | -aligner_index_cache # Whether to try and keep the genome index on the cache/loaded RAM so concurrent jobs do not have to reload it and can use it more easily ('no', which will empty cache at the end, or 'yes', by default)
	        -K | -Kraken2_fast # Kraken2 fast mode, consisting on copying the Kraken2 database to /dev/shm (RAM) so execution is faster ('yes' or 'no' by default)
```

Parameters are not positional. If you did not provide a required parameter, the pipeline may exit or use default values if possible (check the help page above, the log after execution, or the 'Arguments and variables' first section in the main script 'reanalyzerGSE.sh'). For example, if the argument '-s' is not provided, strandness will be predicted using Salmon and how_are_we_stranded_here and transcript sequences would be required, so the pipeline would exit if not provided with the argument '-t'. Similarly, reference genome and annotation are likely going to be required for the alignment and quantifying steps (arguments '-r' and '-a').

In general, from a folder containing raw sequences (.fastq.gz) or a GEO entry (GSEXXXXX) as input (argument '-i'), reanalyzerGSE is going to provide a standard transcriptomic analysis named with the argument '-n', in the folder provided by the argument '-o'. The argument '-p' provides the number of cores to be used when multithreading is possible, '-P' the number of files to be processed simultaneously when possible, leveraging GNU's parallel, and '-M' the maximum amount of RAM memory available (required mostly for the index generation and alignment step). If the reanalysis of a GEO entry is requested, all the samples in the database will be included by default. The parameter '-G' allows to restrict the analysis to some of them (providing GSMXXXXX ids), while '-S' interrupts the pipeline until the user has made manual changes to the downloaded files, and '-d' allows the user to input a custom experimental design, to be used instead of the automatically-detected one.

An improved version of miARma-seq has been included in reanalyzerGSE and used by default (subfolder 'external_software'). Amognst others, reanalyzerGSE also includes a module to perform batch correction (the design matrix must be provided in a prompt after using the flag '-b'), the possibility to use multiple thresholds when filtering out low expressed genes (argument '-f'), and a module to output plots from the differential gene expression analyses, highlighting the genes provided with the argument '-g'.

Please refer to the help page or open an issue for any further clarification.

## Output:
The output folder contains all results in a structure format. Quality control results, amongst others, are aggregate in a MultiQC HTML report. The results by miARma-seq (i.e. quantification and alignment) are container in the folder 'miARma_out'. The folder "final_results" contain all final results and tables, including a Sphynx HTML report guiding users through the output.


## Comments
Please cite this reference when using reanalyzerGSE for your publications:

> Ruiz, J. L., Terrón-Camero, L. C., Castillo-González, J., Fernández-Rengel, I., Delgado, M., Gonzalez-Rey, E., & Andrés-León, E. (2023). reanalyzerGSE: tackling the everlasting lack of reproducibility and reanalyses in transcriptomics. bioRxiv, 2023-07. https://doi.org/10.1101/2023.07.12.548663

```
@article{ruiz2023reanalyzergse,
  title={reanalyzerGSE: tackling the everlasting lack of reproducibility and reanalyses in transcriptomics},
  author={Ruiz, Jose L and Terr{\'o}n-Camero, Laura Carmen and Castillo-Gonz{\'a}lez, Julia and Fern{\'a}ndez-Rengel, Iv{\'a}n and Delgado, Mario and Gonzalez-Rey, Elena and Andr{\'e}s-Le{\'o}n, Eduardo},
  journal={bioRxiv},
  pages={2023--07},
  year={2023},
  publisher={Cold Spring Harbor Laboratory}
}
```


## Support
Please report any [issue](https://github.com/BioinfoIPBLN/reanalyzerGSE/issues) or [contact us](mailto:bioinformatica@ipb.csic.es?subject=[GitHub]%20Source%20reanalyzerGSE%20Support) to request support.
