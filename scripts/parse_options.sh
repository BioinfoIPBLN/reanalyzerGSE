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
		-h*) echo "reanalyzerGSE v3.0.0 - usage: reanalyzerGSE.pk.sh [options]
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
	        -C | -covariables # Please input a comma-separated list for the covariable that will be included in the edgeR model for DGE and in limma::removeBatchEffect as covariate (only one covariable allowed, for example a proven batch effect, and adjusted_counts will be overwritten by ComBat-seq if activated)
	        -T | -target # Protopical target file for attempts to differential gene expression analyses (containing filenames and covariates, automatically built if not provided)
	
	        #### Activate alternative modes:
	        -Dm | -debug_module # For debugging, step to remove the content of the corresponding folders and to resume a failed or incomplete run without repeating (one of 'step1', 'step1a', 'step1b', 'step1c', 'step2', 'step3a', 'step3b', 'step4', 'step5', 'step6', 'step7', 'step8', or 'all' to execute everything, by default)
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
	        -K | -Kraken2_fast # Kraken2 fast mode, consisting on copying the Kraken2 database to /dev/shm (RAM) so execution is faster ('yes' or 'no' by default)" && exit 1;;
		-options) options_file=${arguments[index]} ;;
		-i) input=${arguments[index]} ;;
		-n) name=${arguments[index]} ;;
		-o) output_folder=${arguments[index]} ;;
		-p) cores=${arguments[index]} ;;
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
		-b) batch=${arguments[index]} ;;
		-bv) batch_vector=${arguments[index]} ;;
		-bc) batch_biological_covariates=${arguments[index]} ;;
		-B) bed_mode=${arguments[index]} ;;
		-C) covariables=${arguments[index]} ;;
		-S) stop=${arguments[index]} ;;
		-P) number_parallel=${arguments[index]} ;;
		-R) number_reads_to_subsample=${arguments[index]} ;;
		-t) transcripts=${arguments[index]} ;;
		-T) target=${arguments[index]} ;;
		-A) aligner=${arguments[index]} ;;
		-A) aligner_index_cache=${arguments[index]} ;;
		-K) kraken2_fast=${arguments[index]} ;;
		-Dk) kraken2_databases=${arguments[index]} ;;
		-Ds) sortmerna_databases=${arguments[index]} ;;
		-Des) differential_expr_soft=${arguments[index]} ;;
		-Dm) debug_module=${arguments[index]} ;;
		-Dec) differential_expr_comparisons=${arguments[index]} ;;
		-Dc) deconvolution=${arguments[index]} ;;
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
	esac
done

##### From the external_options file...
if [ ! -z "$options_file" ]; then
	source $options_file
fi

##### Deal with defaults or with the user not providing some...
if [ -z "$input" ] || [ -z "$output_folder" ] || [ -z "$cores" ] || [ -z "$reference_genome" ] || [ -z "$annotation" ]; then
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
echo -e "\ncovariables=$covariables\n"
if [ -z "$deconvolution" ]; then
	deconvolution="no"
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
if [ -z "$kraken2_fast" ]; then
	kraken2_fast="no"
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

seqs_location=$output_folder/$name/raw_reads

# Rename a few:
organism_argument=$organism
number_reads=$number_reads_to_subsample
perform_differential_analyses=$full_differential_analyses
debug_step=$debug_module
minstd=$time_course_std
mestimate=$time_course_fuzz
rev_thr=$revigo_threshold_similarity
indexthreads=$cores_index

if [ ! -z "$kraken2_databases" ]; then
	echo -e "\nPLEASE note some local databases for the decontamination step (kraken2-based) are needed"
	echo -e "These databases may be large, so please be aware that the RAM usage may reach hundreds of GBs. Rerun without the -Dk parameter to skip decontamination if not acceptable"
	echo -e "reanalyzerGSE is now going to check or give you instructions so the databases are downloaded and placed in the corresponding folders"
	if [[ -f $kraken2_databases/taxdump/names.dmp ]] && [[ -f $kraken2_databases/taxdump/nodes.dmp ]] && [ "$(find -L $kraken2_databases -name hash.k2d | wc -l)" -gt 0 ]; then
		echo -e "\nGood, apparently all of the required databases for decontamination based on taxonomic classification have been found in "$kraken2_databases"\n"
	else
		echo -e "\The databases for kraken2 have to be downloaded. The software will exit until the following steps are performed. Please note 'wget' may not complete the download and then uncompressing would give errors. You would need to remove any incomplete file and restart download"
		echo -e "The databases names.dmp, nodes.dmp, merged.dmp and delnodes.dmp have to be downloaded by the user executing: mkdir -p $kraken2_databases/taxdump && cd $kraken2_databases/taxdump && wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip && unzip taxdmp.zip && rm taxdmp.zip"
		echo -e "The kraken2 database has to be built by the user and placed within the $kraken2_databases folder. Please execute the following and read these comments if you want to build a kraken2 database containing the standard kraken2 recommended sequences and the EuPathDB v48 preprocessed sequences (http://ccb.jhu.edu/data/eupathDB/):"
		echo -e "cd $kraken2_databases && cores=64"
		echo -e "kraken2-build --download-taxonomy --threads \$cores --db standard_eupathdb_48_kraken2_db & #### The download commands can be sent to the background so they are simultaneously processed"
		echo -e "kraken2-build --download-library archaea --threads \$cores --db standard_eupathdb_48_kraken2_db &"
		echo -e "kraken2-build --download-library viral --threads \$cores --db standard_eupathdb_48_kraken2_db &"
		echo -e "kraken2-build --download-library plasmid --threads \$cores --db standard_eupathdb_48_kraken2_db &"
		echo -e "kraken2-build --download-library bacteria --threads \$cores --db standard_eupathdb_48_kraken2_db & #### This is around a download of 150GB"
		echo -e "kraken2-build --download-library fungi --threads \$cores --db standard_eupathdb_48_kraken2_db &"
		echo -e "kraken2-build --download-library UniVec_Core --threads \$cores --db standard_eupathdb_48_kraken2_db &"
		echo -e "# kraken2-build --download-library protozoa --threads \$cores --db standard_eupathdb_48_kraken2_db & #### Not required here because protozoa from EuPathDB are going to be already added"
		echo -e "# kraken2-build --download-library nt --threads \$cores --db standard_eupathdb_48_kraken2_db & #### With the nt database results would more precise, but it's huge. It works if you are patient, but it would involve a download of ~450GB, and then building would require to allocate at least 1TB. Then the final database would be of ~450GB of size, and require that much RAM to build and run"
		echo -e "# kraken2-build --add-to-library chr1.fa --threads \$cores --db standard_eupathdb_48_kraken2_db #### If you want to add custom sequences or a particular genome you suspect contamination for... check kraken2 manual to handle taxonomy"
		echo -e "# To add EuPathDB:"
		echo -e "mkdir -p standard_eupathdb_48_kraken2_db/library/added && cd standard_eupathdb_48_kraken2_db/library/added"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/AmoebaDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/CryptoDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/FungiDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/GiardiaDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/MicrosporidiaDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/PiroplasmaDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/PlasmoDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/ToxoDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/TrichDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/TriTrypDB48.tgz &"
		echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB48/seqid2taxid.map &"
		echo -e "for f in *.tgz; do tar -xvzf \$f; done && rm *.tgz \$(ls | grep log)"
		echo "awk '{ print "'"TAXID""\t"$1"\t"$2 }'"' seqid2taxid.map > prelim_map_1.txt && rm seqid2taxid.map"
		echo -e "# To build the kraken2 database:"
		echo -e "kraken2-build --build --threads \$cores --db standard_eupathdb_48_kraken2_db #### This would take around 5 hours and require around 75GB of RAM. If you included nt database, around 24 hours and 500GB of RAM. For the larger database including nt, you may need to include the flag --fast-build to avoid stalling or too much time building"
		echo -e "kraken2-inspect --db standard_eupathdb_48_kraken2_db --threads \$cores > standard_eupathdb_48_kraken2_db_k2_inspect.txt"
		echo -e "kraken2-build --clean --threads \$cores --db standard_eupathdb_48_kraken2_db #### This would remove the downloaded sequences and keep only the kraken2 database. This must be done to save space, but first double check the included sequences and the output of kraken2 inspect. If you execute this command and then you want to make any change to the database, you would have to download everything again"

		echo -e "\n\nAlternatively, you can download ready-to-use databases from https://benlangmead.github.io/aws-indexes/k2"
		echo -e "Please note that if you followed installation instructions, the kraken2-build script has been modified so database building is faster (improved masking, https://github.com/DerrickWood/kraken2/pull/39) and a bug in sequences download from ncbi has been corrected (https://github.com/DerrickWood/kraken2/issues/571)"
		echo -e "Please note that if you have allocated enough RAM and the system is compatible, copying the kraken2 database to the faster RAM, something like /dev/shm/, and then pointing to the library there in kraken2 execution with the flag --memory-mapping would greatly improve speed, particularly if multiple runs (https://github.com/DerrickWood/kraken2/issues/451). reanalyzerGSE includes this mode with the flag -K"
		echo -e "The execution of kraken2 classification with the suggested database (archaea + viral + plasmid + bacteria + fungi + UniVec_Core + EuPathDB) will require ~70GB of RAM."
		echo -e "Please note that the kraken2 version within reanalyzerGSE is frozen. Though it should work, kraken2 is under active development to improve database download and build, sequence databases may undergo crucial changes, and external factors such as updates in NCBI structure may cause errors in the code above. For example, 'core_nt' since to be now more recommended than 'nt'. Please double check you are using the most up-to-date/appropriate databases"
		exit 1
	fi
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
