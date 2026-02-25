#!/bin/bash
CURRENT_DIR_SCRIPTS=$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export PATH=$CURRENT_DIR_SCRIPTS:$PATH; echo -e "\n\nAdding to PATH the scripts folder...$CURRENT_DIR_SCRIPTS"

##### From the command line, get arguments by looping through an index...
index=0
for argument in $options; do

### Incrementing index
	index=`expr $index + 1`

### Gather the parameters, default values, exit if essential not provided...
	case $argument in
		-h*) echo "reanalyzerGSE v3.2.0 - usage: reanalyzerGSE.sh [options]
	        -h | -help # Type this to get help
	        -options | Provide a YAML configuration file containing the parameters to be used. You can adapt the file 'config_template.yaml' provided in the scripts folder. CLI arguments override values from the YAML file.)
	        
	        #### Input/output: 
	        -i | -input # GEO_ID (GSEXXXXXX, separated by comma if more than one), or folder containing raw reads (please provide full absolute path, e.g. /path/folder_name/, containing only fastq.gz files and not folders, links or any other item, and please rename samples with meaningful names if possible), or almost any accession from ENA/SRA to download .fastq from (any of the ids with the prefixes PRJEB,PRJNA,PRJDB,ERP,DRP,SRP,SAMD,SAME,SAMN,ERS,DRS,SRS,ERX,DRX,SRX,ERR,DRR,SRR, please separated by commas if more than one id as input)
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
	        -Dk | -kraken2_databases # Comma-separated list of Kraken2 database folders (e.g. '/path/to/core_nt,/path/to/gtdb'). Any input here activates the kraken2-based decontamination step. All DB+confidence combinations will be run
	        -Kc | -kraken2_confidence # Comma-separated confidence scores for Kraken2 classification (default '0', e.g. '0,0.20,0.50'). Each score is run for each database
	        -Ds | -sortmerna_databases # The database (absolute pathway) that should be used by SortMeRNA (any input here, e.g. '/path/to/rRNA_databases/smr_v4.3_sensitive_db.fasta', would activate the sortmerna-based rRNA removal step)
	        -Df | -databases_function # Manually provide a comma separated list of databases to be used in automatic functional enrichment analyses of DEGs (check out the R package autoGO::choose_database(), but the most popular GO terms are used by default)
	        -nrf | -non_reference_funct_enrichm # Pathway to a file containing functional annotation (GAF, GFF, or GTF) to be used for functional enrichment when the organism is not Human or Mouse. If provided, this overrides the default behavior of using the main annotation file.
	
	        #### Metadata and sample info:
	        -D | -design_custom_local # Specifying here the experimental design for the local dataset (by default an interactive prompt will ask for a comma-separated list of the same length than the number of samples, if you want to avoid that manual input please provide the list in this argument. If no design, please provide a dummy one, e.g. with every sample in the same group or a separate one. If more than one design to provide, please input comma-separated list separated by a '/', without spaces. Please avoid naming that would match the same pattern in grep, e.g. XXXX and XXXX_TREAT)
	        -d | -design_custom # Manually specifying the experimental design for GEO download ('no' by default and if 'yes', please expect an interactive prompt after data download from GEO, and please enter the assignment to groups when asked in the terminal, with a comma-separated list of the same length than the number of samples)
	        -O | -organism # Specifying here the scientific name of the organism for the local dataset (by default an interactive prompt will ask for it, if you want to avoid that manual input please provide the full organism name in this argument, please use underline instead of space)
	        -Tx | -taxon_id # NCBI's taxon id of the organism
	        -R | -number_reads_to_subsample # Subsample raw reads to reach an approximate target library size (none by default). Provide a comma-separated pair: path to the 'reads_numbers.txt' file from a previous execution, and target library size (i.e. the number of counted/mapped reads you want to reach, NOT the number of raw reads). The pipeline will use the raw-read and library-size columns in reads_numbers.txt to compute the proportional number of raw reads to keep for each sample so that the resulting library size is approximately +-10% of the target. If the sample already has fewer reads than the target, all its reads are kept.
	        -bv | -batch_vector # Comma-separated list of numbers for use as batch vector with Combat-seq
	        -bc | -batch_biological_covariable # Comma-separated list of numbers for use as batch vector of covariables of biological interest with Combat-seq
			-bf | -batch_format # Format of the provided batch variables ('num' for numeric/vector variables or 'fact' for factors, by default)
	        -C | -covariables # Please input a comma-separated list for the covariable that should be included in the limma model for removeBatchEffect or in the edgeR model for DGE (for now only one covariable allowed, for example an expected batch effect)
	        -Cf | -covariables_format # Format of the provided covariate ('num' by default for numeric covariables, or 'fact' for factors)
	        -T | -target # Protopical target file for attempts to differential gene expression analyses (containing filenames and covariates, automatically built if not provided)
	
	        #### Activate alternative modes:
	        -Dm | -debug_module # For debugging, step to remove the content of the corresponding folders and to resume a failed or incomplete run without repeating (one of 'step1', 'step1a', 'step1b', 'step1c', 'step1d', 'step2', 'step3a', 'step3b', 'step4', 'step5', 'step6', 'step7', 'step8', step9', or 'all' to execute everything, by default)
	        -q | -qc_raw_reads # Whether to perform quality control on the raw reads ('yes' by default, or 'no')
	        -fd | -full_differential_analyses # Whether to perform full differential enrichment analyses (for example including computation of DEGs or Venn diagrams, 'no' or 'yes', by default)
	        -fe | -functional_enrichment_analyses # Whether to perform functional enrichment analyses ('no' or 'yes', by default)
	        -cPf | -clusterProfiler_full # Whether to perform additional functional enrichment analyses with multiple databases using clusterProfiler, by default only ORA for GO BP, GO MF and GO CC, and KEGG and REACTOME enrichment, will be performed, as additional analyses may be slow if many significant DEGs or multiple number of comparisons ('yes' or 'no', by default)
	        -b | -batch # Batch effect present? (no by default, yes if correction through Combat-seq and model is to be performed, and info is going to be required in other arguments or prompts)
	        -B | -bed_mode # Whether to convert list of files to bed format so they can be visualized in genome browsers ('yes' or 'no', by default)
	        -Dc | -deconvolution # Deconvolution method for bulk RNA-seq data: 'no' (default), 'CDSeq' (unsupervised, may take hours), or 'BisqueRNA' (reference-based, requires -scM, -scP, and -bulkM)
	        -scM | -sc_count_matrix # Path to sc/snRNA count matrix for BisqueRNA deconvolution (tab-delimited, first column 'Gene', remaining columns are cells, values are counts)
	        -scP | -sc_phenotype # Path to sc/snRNA phenotype file for BisqueRNA (tab-delimited, columns: 'SubjectName' and 'cellType', rows = cell IDs matching count matrix columns)
	        -bulkM | -bulk_expression_matrix # Path to bulk expression matrix for BisqueRNA (tab-delimited, first column 'Gene', remaining columns are samples)
	        -vv | -perform_volcano_venn # Whether to perform all Volcano plots and Venn diagrams, which may take a long time if many comparisons ('no' or 'yes', by default)
	        -aP | -aPEAR_execution # Whether to simplify pathway enrichment analysis results by detecting clusters of similar pathways and visualizing enrichment networks by aPEAR package, which may be slow ('yes' or 'no', by default)
	        -Ti | -tidy_tmp_files # Space-efficient run, with a last step removing raw reads if downloaded, converting bam to cram, removing tmp files... etc ('yes' or 'no', by default)
	        -Txls | -convert_tables_excel # Convert all tables in results from .txt format, without limitation of size to Excel's .xlsx format, with a limitation of 32,767 characters ('yes' or 'no', by default)
	        -Tc | -time_course # Whether to perform additional time-course analyses as a last step ('yes' or 'no', by default)
	        -Na | -network_analyses # Whether to perform network analyses, only for human or mouse analyses ('yes' or 'no', by default)
	        -apl | -auto_panther_log # Whether to perform additional autoGO and Panther analyses for DEGs separated by log2Fc positive or negative ('yes' or 'no', by default)
	        -eDe | -exploreDE_se # Whether to generate a SummarizedExperiment object (.qs2) for the exploreDE Shiny app ('yes' or 'no', by default). Only works for Human (Homo_sapiens) analyses.
	
	        #### Processing parameters:
	        -s | -strand # Strandness of the library ('yes, 'no', 'reverse'). If not provided and '-t' used, this would be predicted by salmon. Please use this parameter if prediction not correct, see explanations in for example in bit.ly/strandness0 and bit.ly/strandness
	        -f | -filter # Threshold of gene counts to use ('bin' to capture the lower expressed genes, 'filterbyexpr' to use the edgeR solution, 'or 'standard', by default). Please provide a comma separated list with the filters to use at each quantification if multiple annotation are provided
	        -Of | -options_featureCounts_feat # The feature type to use to count in featureCounts (default 'exon')
	        -Os | -options_featureCounts_seq # The seqid type to use to count in featureCounts (default 'gene_name')
	        -A | -aligner # Aligner software to use ('hisat2' or 'star', by default)
	        input_filter_regex=""
	        input_filter_regex_exclude=""
        alignment_removal=""
        cores=8
	        -Des | -differential_expr_software # Software to be used in the differential expression analyses ('edgeR' by default, or 'DESeq2')
	        -fp | -fastp_mode # Whether to perform quality filtering on the raw reads by fastp ('yes' or 'no', by default)
	        -fpa | -fastp_adapter # Whether to perform adapter trimming on the raw reads by fastp ('yes' or 'no', by default, to perform automatic trimming, or a path to a fasta file to perform trimming of its sequences)
	        -fpt | -fastp_trimming # Whether to trim the raw reads by fastp ('none' by default, if two numbers separated by comma, the indicated number of bases will be trimmed from the front and tail, respectively)
	        -std | -time_course_std # Standard deviation threshold to filter in time course analyses (numeric, 1 by default)
	        -fuzz | -time_course_fuzz # Fuziness value for the soft clustering approach (by default an estimate is automatically computed but manual testing is encouraged)
	
	        #### Filtering out samples/comparisons:
	        -G | -GSM_filter # GSM ids (one or several, separated by comma and no space) within the GSE entry to restrict the analysis to. An alternative to requesting a stop with -S to reorganize the downloaded files manually
	        -S | -stop # Manual stop so the automatically downloaded files can be manually modified ('yes' or 'no', by default)
	        -pR | -pattern_to_remove # A pattern to exclude matching samples from downstream R processing only, i.e. QC figures and DGE analyses (by default 'none'). Unlike -regex/-regexExclude which filter raw reads before alignment, this option keeps all samples through alignment and counting but excludes matching ones at the R analysis stage. Useful for removing outlier samples without re-running the full pipeline (e.g. resume from -Dm step4)
	        -Dec | -differential_expr_comparisons # Whether to restrict the differential expression analyses to only some of the possible comparisons or reorder the 'treatment' and 'control' elements of the comparison ('no', by default, or a comma-separated list specifying separated by '\\' the elements in the comparison, which you could get from a preliminar previous run, e.g. 'A//B,C//D,D//A'...)
	
	        #### Functional enrichment/networking analyses
	        -cPm | -clusterProfiler_method # Method for adjusting p.value in clusterProfiler iterations (one of 'holm','hochberg','hommel','bonferroni','BH','BY,'none', or 'fdr', by default)
	        -cPu | -clusterProfiler_universe # Universe to use in clusterProfiler iterations ('all' or 'detected', by default)
	        -mGS | -clusterProfiler_minGSSize # minGSSize parameter to use in clusterProfiler iterations (a number, '10' by default)
	        -MGS | -clusterProfiler_maxGSSize # maxGSSize parameter to use in clusterProfiler iterations (a number, '500' by default)
	        -Pm | -panther_method # Method for adjusting p.value in panther analyses via rbioapi (one of 'NONE','BONFERRONI', or 'FDR', by default)
	        -rev | -revigo_threshold_similarity # Similarity threshold for Revigo summaries of GO terms (0-1, suggested values are 0.9, 0.7, 0.5, 0.4 for large, medium, small, and tiny levels of similarity, respectively, being default 0.7
	
	        #### BAM filtering (post-alignment):
	        -mQ | -bam_mapq_threshold # Min MAPQ for samtools view -q and featureCounts -Q (0 = use default)
	        -Fex | -bam_exclude_flags # samtools -F flags to exclude, e.g. '4' (unmapped), '256' (secondary), '2308' (combined)
	        -Freq | -bam_require_flags # samtools -f flags to require. Suggestion: use '2' for Paired-End (proper pair), leave empty or '4' for Single-End.
	        -Fdup | -bam_dedup # Duplicate removal: 'no' (default), 'samtools' (markdup -r), 'picard' (REMOVE_DUPLICATES), 'picard_optical' (REMOVE_SEQUENCING_DUPLICATES)
	
	        #### Performance:
	        -regex | -input_filter_regex # Regex to keep only matching input files in local mode, removing the rest (e.g. "Sample_A|Sample_B")
	        -regexExclude | -input_filter_regex_exclude # Regex to exclude matching input files in local mode, keeping the rest (e.g. "Sample_Bad|Sample_Outlier")
	        -Ar | -alignment_removal # Fasta file to map against and remove aligned reads (e.g. host genome)
	        -p | -cores # Number of cores
	        -cR | -cores_reads_to_subsample # Cores to use in subsampling by seqtk (10 by default)
	        -pi | -cores_index # Number of cores for genome indexing in aligning step (by default, same than -p)
	        -M | -memory_max # Max RAM memory to be used by aligner or JAVA in bytes (by default 257698037760, or 240GB, used)
	        -P | -number_parallel # Number of files to be processed in parallel (10 by default)
	        -cG | -compression_level # Specify the compression level to gzip the downloaded fastq files from GEO (numeric '0' to '9', default '9')
	        -Ac | -aligner_index_cache # Whether to try and keep the genome index on the cache/loaded RAM so concurrent jobs do not have to reload it and can use it more easily ('no', which will empty cache at the end, or 'yes', by default)
	        -K | -Kraken2_fast # DEPRECATED: daemon mode now replaces /dev/shm approach. This option is kept for backward compatibility but has no effect" && exit 1;;
		-options) options_file=${arguments[index]} ;;
		-i) input=${arguments[index]} ;;
		-n) name=${arguments[index]} ;;
		-o) output_folder=${arguments[index]} ;;
		-p) cores=${arguments[index]} ;;
        -regex | -input_filter_regex) input_filter_regex=${arguments[index]} ;;
        -regexExclude | -input_filter_regex_exclude) input_filter_regex_exclude=${arguments[index]} ;;
        -Ar | -alignment_removal) alignment_removal=${arguments[index]} ;;
		-pi) cores_index=${arguments[index]} ;;
		-M) memory_max=${arguments[index]} ;;
		-s) strand=${arguments[index]} ;;
		-r) reference_genome=${arguments[index]} ;;
		-ri) reference_genome_index=${arguments[index]} ;;
		-g) genes=${arguments[index]} ;;
		-G) GSM_filter=${arguments[index]} ;;
		-rev) revigo_threshold_similarity=${arguments[index]} ;;
		-a) annotation=${arguments[index]} ;;
		-f) filter=${arguments[index]} ;;
		-d) design_custom=${arguments[index]} ;;
		-D) design_custom_local=${arguments[index]} ;;
		-Df) databases_function=${arguments[index]} ;;
		-nrf) non_reference_funct_enrichm=${arguments[index]} ;;
		-b) batch=${arguments[index]} ;;
		-bv) batch_vector=${arguments[index]} ;;
		-bc) batch_biological_covariates=${arguments[index]} ;;
		-bf) batch_format=${arguments[index]} ;;
		-B) bed_mode=${arguments[index]} ;;
		-C) covariables=${arguments[index]} ;;
		-Cf) covariables_format=${arguments[index]} ;;
		-S) stop=${arguments[index]} ;;
		-P) number_parallel=${arguments[index]} ;;
		-R) number_reads_to_subsample=${arguments[index]} ;;
		-t) transcripts=${arguments[index]} ;;
		-T) target=${arguments[index]} ;;
		-A) aligner=${arguments[index]} ;;
		-A) aligner_index_cache=${arguments[index]} ;;
		-K) kraken2_fast=${arguments[index]}; echo "WARNING: -K/Kraken2_fast is deprecated; k2 daemon mode is used instead. This option has no effect." ;;
		-Dk) kraken2_databases=${arguments[index]} ;;
		-Kc | -kraken2_confidence) kraken2_confidence=${arguments[index]} ;;
		-Ds) sortmerna_databases=${arguments[index]} ;;
		-Des) differential_expr_soft=${arguments[index]} ;;
		-Dm) debug_module=${arguments[index]} ;;
		-Dec) differential_expr_comparisons=${arguments[index]} ;;
		-Dc) deconvolution=${arguments[index]} ;;
		-scM | -sc_count_matrix) sc_count_matrix=${arguments[index]} ;;
		-scP | -sc_phenotype) sc_phenotype=${arguments[index]} ;;
		-bulkM | -bulk_expression_matrix) bulk_expression_matrix=${arguments[index]} ;;
		-cPf) clusterProfiler_full=${arguments[index]} ;;
		-fe) functional_enrichment_analyses=${arguments[index]} ;;
		-fd) full_differential_analyses=${arguments[index]} ;;
		-vv) perform_volcano_venn=${arguments[index]} ;;
		-aP) aPEAR_execution=${arguments[index]} ;;
		-Of) optionsFeatureCounts_feat=${arguments[index]} ;;
		-O) organism=${arguments[index]} ;;
		-Os) optionsFeatureCounts_seq=${arguments[index]} ;;
		-iG) input_geo_reads=${arguments[index]} ;;
		-cG) compression_level=${arguments[index]} ;;
		-Ti) tidy_tmp_files=${arguments[index]} ;;
		-Txls) convert_tables_excel=${arguments[index]} ;;
		-Tx) taxonid=${arguments[index]} ;;
		-Tc) time_course=${arguments[index]} ;;
		-std) time_course_std=${arguments[index]} ;;
		-fuzz) time_course_fuzz=${arguments[index]} ;;
		-TMP) TMPDIR_arg=${arguments[index]} ;;
		-q) qc_raw_reads=${arguments[index]} ;;
		-cPm) clusterProfiler_method=${arguments[index]} ;;
		-cPu) clusterProfiler_universe=${arguments[index]} ;;
		-mGS) clusterProfiler_minGSSize=${arguments[index]} ;;
		-MGS) clusterProfiler_maxGSSize=${arguments[index]} ;;
		-Pm) panther_method=${arguments[index]} ;;
		-Na) network_analyses=${arguments[index]} ;;
		-fp) fastp_mode=${arguments[index]} ;;
		-fpa) fastp_adapter=${arguments[index]} ;;
		-fpt) fastp_trimming=${arguments[index]} ;;
		-cR) cores_reads_to_subsample=${arguments[index]} ;;
		-pR) pattern_to_remove=${arguments[index]} ;;
		-apl) auto_panther_log=${arguments[index]} ;;
		-mQ) bam_mapq_threshold=${arguments[index]} ;;
		-Fex) bam_exclude_flags=${arguments[index]} ;;
		-Freq) bam_require_flags=${arguments[index]} ;;
		-Fdup) bam_dedup=${arguments[index]} ;;
		-nrf) non_reference_funct_enrichm=${arguments[index]} ;;
		-eDe) exploreDE_se=${arguments[index]} ;;
	esac
done

##### From the YAML configuration file...
if [ ! -z "$options_file" ]; then
	if ! command -v yq &> /dev/null; then
		echo "Error: yq is required to parse YAML config files but was not found in PATH. Please install yq."; exit 1
	fi
	if [ ! -f "$options_file" ]; then
		echo "Error: Config file '$options_file' not found."; exit 1
	fi
	echo -e "\nLoading configuration from YAML file: $options_file\n"
	while IFS='=' read -r key val; do
		# Skip empty keys or null values
		[ -z "$key" ] && continue
		[ "$val" = "null" ] && continue
		# Only set if not already defined by CLI arguments (CLI takes priority)
		if [ -z "${!key}" ]; then
			export "$key=$val"
		fi
	done < <(yq -r 'to_entries[] | select(.value != null) | "\(.key)=\(.value)"' "$options_file")
fi

##### Deal with defaults or with the user not providing some...
if { [ -z "$input" ] && [ -z "$input_geo_reads" ]; } || [ -z "$output_folder" ] || [ -z "$cores" ] || [ -z "$reference_genome" ] || [ -z "$annotation" ]; then
	echo "Please check usage with 'reanalyzer_GSE_RNA_seq.sh -h'. At least one required argument has not been provided..."; exit 1
fi

echo -e "\n\nArguments:"
echo -e "\ninput=$input\n"
echo -e "\noutput_folder=$output_folder\n"
echo -e "\ncores=$cores\n"
if [ -z "$name" ]; then
	if [[ $input == G* ]]; then
		arrIN=(${input//,/ }); name=$(for a in "${arrIN[@]}"; do echo "$a"; done | sort | tr '\n' '_' | sed 's,_$,,g')
	elif [[ $input == /* ]]; then
		name=$(basename $input)
	else
		echo "Please double check you have provided either a GEO ID (GSEXXXXXX) or a full pathway to input folder (/path/folder_name/) as the argument '-i'"; exit 1
	fi
else
	if [[ $input == G* ]]; then
		arrIN=(${input//,/ }); name=$(for a in "${arrIN[@]}"; do echo "$a"; done | sort | tr '\n' '_' | sed 's,_$,,g')
		echo -e "\nOverwriting name because you are downloading from GEO...\nname=$name"
	fi
fi
echo -e "\nname=$name\n"
if [ -z "$number_parallel" ]; then
	number_parallel=10
fi
echo -e "\nnumber_parallel=$number_parallel\n"
echo "Please make sure that number_parallel above is a number lower than the total number of samples/sequences that are going to be analyzed in the run..."
echo -e "\nreference_genome=$reference_genome\n"
echo -e "\nannotation=$annotation\n"
if [ -z "$transcripts" ] && [ -z "$strand" ]; then
	echo "Please check usage with 'reanalyzer_GSE_RNA_seq.sh -h'. At least one of '-t' or '-s' are required..."; exit 1
elif [ ! -z "$strand" ]; then
	echo -e "\nstrand=$strand\n"
elif [ ! -z "$transcripts" ]; then
	echo -e "\ntranscripts=$transcripts\n"
fi
if [ -z "$strand" ]; then
	echo -e "\nsalmon is going to be used to predict strandness using the transcript sequences in '-t'. If not correct, please provide it with '-s'\n"
fi
if [ -z "$genes" ]; then
	genes="none"
else
	echo -e "\ngenes=$genes\n"
fi
if [ ! -z "$GSM_filter" ]; then
	echo -e "\nGSM_filter=$GSM_filter\n"
fi
if [ ! -z "$number_reads_to_subsample" ]; then
	echo -e "\number_reads_to_subsample=$number_reads_to_subsample\n"
fi
if [ -z "$filter" ]; then
	filter=standard
fi
echo -e "\nfilter=$filter\n"
if [ -z "$batch" ]; then
	batch="no"
fi
echo -e "\nbatch=$batch\n"
if [ -z "$covariables" ]; then
	covariables="none"
fi
if [ -z "$batch_format" ]; then
	batch_format="fact"
fi
echo -e "\ncovariables_format=$covariables_format\n"
echo -e "\ncovariables=$covariables\n"
if [ -z "$covariables_format" ]; then
	covariables_format="num"
fi
echo -e "\ncovariables_format=$covariables_format\n"
if [ -z "$deconvolution" ]; then
	deconvolution="no"
fi
# Backward compat: treat 'yes' as 'CDSeq'
if [ "$deconvolution" == "yes" ]; then
	deconvolution="CDSeq"
fi
if [ -z "$sc_count_matrix" ]; then
	sc_count_matrix="none"
fi
if [ -z "$sc_phenotype" ]; then
	sc_phenotype="none"
fi
if [ -z "$bulk_expression_matrix" ]; then
	bulk_expression_matrix="none"
fi
if [ -z "$design_custom" ]; then
	design_custom="no"
fi
echo -e "\ndesign_custom=$design_custom\n"
if [ -z "$stop" ]; then
	stop="no"
fi
echo -e "\nstop=$stop\n"
if [ -z "$memory_max" ]; then
	memory_max=257698037760
fi
echo -e "\nmemory_max=$memory_max\n"
if [ -z "$miarma_path" ]; then
	miarma_path=$(dirname $CURRENT_DIR_SCRIPTS)/external_software/miARma-seq
fi
if [ -d "$output_folder/$name" ]; then
	echo -e "Please note that $output_folder/$name already exists... reanalyzerGSE is going to attempt a new run or resume running, but you may want to remove the folder, change the destination folder with '-o' or '-n', use downloaded raw data from an external software... etc. Sleeping for a while to give you time to exit if you want, and then continuing..."
	secs=$((1 * 15))
	while [ $secs -gt 0 ]; do
		echo -ne "$secs\033[0K\r"
		sleep 1
		: $((secs--))
	done
fi
if [ -z "$databases_function" ]; then
	databases_function="GO_Biological_Process_2023,GO_Molecular_Function_2023,GO_Cellular_Component_2023"
fi
if [ -z "$aligner" ]; then
	aligner="star"
fi
if [ -z "$aligner_index_cache" ]; then
	aligner_index_cache="yes"
fi
if [ -z "$cores_index" ]; then
	cores_index=$cores
fi
if [ -z "$clusterProfiler_full" ]; then
	clusterProfiler_full="no"
fi
if [ -z "$functional_enrichment_analyses" ]; then
	functional_enrichment_analyses="yes"
fi
if [ -z "$full_differential_analyses" ]; then
	full_differential_analyses="yes"
fi
if [ -z "$perform_volcano_venn" ]; then
	perform_volcano_venn="yes"
fi
if [ -z "$aPEAR_execution" ]; then
	aPEAR_execution="no"
fi
if [ -z "$time_course" ]; then
	time_course="no"
fi
if [ -z "$convert_tables_excel" ]; then
	convert_tables_excel="no"
fi
if [ -z "$time_course_fuzz" ]; then
	time_course_fuzz=0
fi
if [ -z "$time_course_std" ]; then
	time_course_std=1
fi
if [ -z "$clusterProfiler_method" ]; then
	clusterProfiler_method="fdr"
fi
if [ -z "$panther_method" ]; then
	panther_method="FDR"
fi
if [ -z "$differential_expr_soft" ]; then
	differential_expr_soft="edgeR"
fi
if [ -z "$differential_expr_comparisons" ]; then
	differential_expr_comparisons="no"
fi
if [ -z "$target" ]; then
	target="no"
fi
if [ -z "$tidy_tmp_files" ]; then
	tidy_tmp_files="no"
fi
if [ -z "$network_analyses" ]; then
	network_analyses="no"
fi
if [ -z "$clusterProfiler_universe" ]; then
	clusterProfiler_universe="detected"
fi
if [ -z "$clusterProfiler_minGSSize" ]; then
	clusterProfiler_minGSSize=10
fi
if [ -z "$clusterProfiler_maxGSSize" ]; then
	clusterProfiler_maxGSSize=500
fi
if [ -z "$bed_mode" ]; then
	bed_mode="no"
fi
if [ -z "$pattern_to_remove" ]; then
	pattern_to_remove="none"
fi
if [ -z "$debug_module" ]; then
	debug_module="all"
fi
### Info on the kraken2 decontamination step and databases:
if [ -z "$kraken2_confidence" ]; then
	kraken2_confidence="0"
fi
if [ -z "$compression_level" ]; then
	compression_level=9
fi
if [ -z "$revigo_threshold_similarity" ]; then
	revigo_threshold_similarity=0.7
fi
if [ -z "$fastp_mode" ]; then
	fastp_mode="no"
fi
if [ -z "$fastp_adapter" ]; then
	fastp_adapter="no"
fi
if [ -z "$fastp_trimming" ]; then
	fastp_trimming="none"
fi
if [ -z "$auto_panther_log" ]; then
	auto_panther_log="no"
fi
if [ -z "$cores_reads_to_subsample" ]; then
	cores_reads_to_subsample=10
fi
if [ -z "$bam_mapq_threshold" ]; then
	bam_mapq_threshold=0
fi
if [ -z "$bam_exclude_flags" ]; then
	bam_exclude_flags=""
fi
if [ -z "$bam_require_flags" ]; then
	bam_require_flags=""
fi
if [ -z "$bam_dedup" ]; then
	bam_dedup="no"
fi
if [ -z "$exploreDE_se" ]; then
	exploreDE_se="no"
fi

seqs_location=$output_folder/$name/raw_reads

# Rename a few:
organism_argument=$organism
number_reads=$number_reads_to_subsample
perform_differential_analyses=$full_differential_analyses
if [ -z "$debug_step" ]; then
	debug_step=$debug_module
fi
minstd=$time_course_std
mestimate=$time_course_fuzz
rev_thr=$revigo_threshold_similarity
indexthreads=$cores_index

if [ ! -z "$kraken2_databases" ]; then
	echo -e "\nPLEASE note some local databases for the decontamination step (kraken2-based) are needed"
	echo -e "These databases may be large, so please be aware that the RAM usage may reach hundreds of GBs. Rerun without the -Dk parameter to skip decontamination if not acceptable"
	echo -e "reanalyzerGSE is now going to check or give you instructions so the databases are downloaded and placed in the corresponding folders"
	IFS=',' read -r -a k2_db_array <<< "$kraken2_databases"
	for k2_db in "${k2_db_array[@]}"; do
		if [ "$(find -L $k2_db -name hash.k2d 2>/dev/null | wc -l)" -gt 0 ]; then
			echo -e "\nGood, Kraken2 database found in $k2_db\n"
		else
			echo -e "\nERROR: Kraken2 database not found in $k2_db (missing hash.k2d)"
			echo -e "\nTo download ready-to-use databases: https://benlangmead.github.io/aws-indexes/k2"
			echo -e "\nTo build your own database with k2, example commands:"
			echo -e "  mkdir -p $k2_db && cd $k2_db && cores=64"
			echo -e "  k2 download-taxonomy --threads \$cores --db $k2_db &"
			echo -e "  k2 download-library archaea --threads \$cores --db $k2_db &"
			echo -e "  k2 download-library viral --threads \$cores --db $k2_db &"
			echo -e "  k2 download-library plasmid --threads \$cores --db $k2_db &"
			echo -e "  k2 download-library bacteria --threads \$cores --db $k2_db & # ~150GB download"
			echo -e "  k2 download-library fungi --threads \$cores --db $k2_db &"
			echo -e "  k2 download-library UniVec_Core --threads \$cores --db $k2_db &"
			echo -e "  # k2 download-library nt --threads \$cores --db $k2_db & # Huge (~450GB), requires ~1TB to build"
			echo -e "  # k2 add-to-library custom.fa --threads \$cores --db $k2_db # For custom sequences"
			echo -e "  k2 build --threads \$cores --db $k2_db # ~5h, ~75GB RAM for standard; ~24h, ~500GB for nt"
			echo -e "  k2 inspect --db $k2_db --threads \$cores > ${k2_db}_k2_inspect.txt"
			echo -e "  k2 clean --db $k2_db # Remove downloaded seqs, keep only DB files (saves space)"
			echo -e "\nAlternatively, download the recommended core_nt or gtdb databases from the link above."
			echo -e "Exiting..."
			exit 1
		fi
	done
	echo -e "Confidence scores to be used: $kraken2_confidence"
	echo -e "Note: reanalyzerGSE now uses k2 daemon mode for fast classification. The old -K/Kraken2_fast option is deprecated."
fi
if [ ! -z "$sortmerna_databases" ]; then
	echo -e "\nPLEASE note some local databases for the step of rRNA removal (sortmerna-based) are needed"
	echo -e "reanalyzerGSE is now going to check or give you instructions so the databases are downloaded and placed in the corresponding folders"
	if [[ -f $sortmerna_databases ]]; then
		echo -e "\nGood, apparently the required database for decontamination based on rRNA removal have been found in "$sortmerna_databases"\n"
	else
		echo -e "\The databases for sortmerna have to be downloaded. The software will exit until the following steps are performed. Please note 'wget' may not complete the download and then uncompressing would give errors. You would need to remove any incomplete file and restart download"
		echo -e "The databases have to be downloaded by the user executing: mkdir -p $(dirname $sortmerna_databases) && cd $(dirname $sortmerna_databases) && wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz && tar xzf database.tar.gz && rm database.tar.gz"
		echo -e "Please note that multiple databases are downloaded. Please choose one (e.g. smr_v4.3_sensitive_db.fasta) and remove the rest"
		exit 1
	fi
fi

### Handle TMPDIR, problematic in some systems or with some runs
if [ -z "$TMPDIR_arg" ]; then
	export TMPDIR=$output_folder/$name/miARma_out0/tmpdir; mkdir -p $TMPDIR
	echo -e "\nTMPDIR changed to $TMPDIR\n"
elif [ "$TMPDIR_arg" == "system" ]; then
	echo -e "\nTMPDIR used is $TMPDIR\n"
else
	export TMPDIR=$TMPDIR_arg; mkdir -p $TMPDIR
	echo -e "\nTMPDIR changed to $TMPDIR\n"
fi
