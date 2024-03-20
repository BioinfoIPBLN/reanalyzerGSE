#!/bin/bash
start=`date +%s`
echo -e "\nCurrent time: $(date)\n"
base64 -d <<<"CiAgX19fICBfX19fICBfX19fICAgICAgX19fXyAgX19fXyAgIF9fICAgX18gXyAgIF9fICAgX18gICAgXyAgXyAgX19fXyAgX19fXyAgX19fXyAKIC8gX18pLyBfX18pKCAgX18pICAgICggIF8gXCggIF9fKSAvIF9cICggICggXCAvIF9cICggICkgICggXC8gKShfXyAgKSggIF9fKSggIF8gXAooIChfIFxcX19fIFwgKSBfKSAgICAgICkgICAvICkgXykgLyAgICBcLyAgICAvLyAgICBcLyAoXy9cICkgIC8gIC8gXy8gICkgXykgICkgICAvCiBcX19fLyhfX19fLyhfX19fKSAgICAoX19cXykoX19fXylcXy9cXy9cXylfXylcXy9cXy9cX19fXy8oX18vICAoX19fXykoX19fXykoX19cXykKCmJ5IEJpb2luZm9ybWF0aWNzIFVuaXQJCQkJSVBCTE4tQ1NJQy4gMjAyMwoKYmlvaW5mb3JtYXRpY2FAaXBiLmNzaWMuZXMJCSAgICAgICAgaHR0cHM6Ly9naXRodWIuY29tL0Jpb2luZm9JUEJMTi9yZWFuYWx5emVyR1NFCgo="
echo -e "doi.org/10.1101/2023.07.12.548663\n\n"


###### 0. Define arguments and variables:
### A string with command options and an array with arguments
options=$@
arguments=($options)

### Get arguments by looping through an index
index=0
for argument in $options; do

### Incrementing index
	index=`expr $index + 1`

### Gather the parameters, default values, exit if essential not provided...
	case $argument in
		-h*) echo "reanalyzerGSE v2.3.0 - usage: reanalyzerGSE.pk.sh [options]
		-h | -help # Type this to get help
		-i | -input # GEO_ID (GSEXXXXXX, separated by comma if more than one), or folder containing raw reads (please provide full absolute path, e.g. /path/folder_name/), or almost any accession from ENA/SRA to download .fastq from (any of the ids with the prefixes PRJEB,PRJNA,PRJDB,ERP,DRP,SRP,SAMD,SAME,SAMN,ERS,DRS,SRS,ERX,DRX,SRX,ERR,DRR,SRR, please separated by commas if more than one id as input)
		-n | -name # Name of the project/folder to create and store results
		-o | -output_folder # Destination folder
		-p | -cores # Number of cores
		-pi | -cores_index # Number of cores for genome indexing in aligning step (by default, same than -p)
		-P | -parallel_number # Number of files to be processed in parallel (10 by default)
		-r | -reference_genome # Reference genome to be used (.fasta file or .gz, absolute pathway)
		-ri | -reference_genome_index # If the reference genome to be used already has an index that would like to reuse, please provide full pathway here (by default the provided genome is indexed)
		-q  | -qc_raw_reads # Whether to perform quality control on the raw reads ('yes' by default, or 'no')
		-a | -annotation # Reference annotation to be used (.gtf file, absolute pathway). If hisat2 is used, a gff file (make sure format is '.gff' and not '.gff3') is accepted (some QC steps like 'qualimap rnaseqqc' may be skipped though). You can provide a comma-separated list of the pathways to different annotation, and multiple/independent quantification/outputs from the same alignments will be generated.
		-A | -aligner # Aligner software to use ('hisat2' or 'star', by default)
		-Ac | -aligner_index_cache # Whether to try and keep the genome index on the cache/loaded RAM so concurrent jobs do not have to reload it and can use it more easily ('no', which will empty cache at the end, or 'yes', by default)
		-t | -transcripts # Reference transcripts to be used (.fasta cDNA file, absolute pathway, only used if '-s' argument not provided so salmon prediction of strandness is required)
		-s | -strandness # Strandness of the library ('yes, 'no', 'reverse'). If not provided and '-t' used, this would be predicted by salmon. Please use this parameter if prediction not correct, see explanations in for example in bit.ly/strandness0 and bit.ly/strandness
		-g | -genes # Genes to highlight their expression in plots (one or several, separated by comma and no space)
		-G | -GSM_filter # GSM ids (one or several, separated by comma and no space) within the GSE entry to restrict the analysis to. An alternative to requesting a stop with -S to reorganize the downloaded files manually
		-R | -reads_to_subsample # Number of reads to subsample to the sequences before the analyses (none by default, a comma-separated list with one number per fastq file/pair of files if paired-end must be provided)
		-f | -filter # Threshold of gene counts to use ('bin' to capture the lower expressed genes, or 'standard', by default). Please provide a comma separated list with the filters to use at each quantification if multiple annotation are provided
		-b | -batch # Batch effect present? (no by default, yes if correction through Combat-seq and model is to be performed, and info is going to be required in other arguments or prompts)
		-bv | -batch_vector # Comma-separated list of numbers for use as batch vector with Combat-seq
		-bc | -batch_biological_covariable # Comma-separated list of numbers for use as batch vector of covariables of biological interest with Combat-seq
		-d | -design_custom # Manually specifying the experimental design for GEO download ('no' by default and if 'yes', please expect an interactive prompt after data download from GEO, and please enter the assignment to groups when asked in the terminal, with a comma-separated list of the same length than the number of samples)
		-D | -design_custom_local # Specifying here the experimental design for the local dataset (by default an interactive prompt will ask for a comma-separated list of the same length than the number of samples, if you want to avoid that manual input please provide the list in this argument)
		-O | -organism # Specifying here the scientific name of the organism for the local dataset (by default an interactive prompt will ask for it, if you want to avoid that manual input please provide the full organism name in this argument, please use underline instead of space)
		-C | -covariables # Please input a comma-separated list for the covariable that should be included in the edgeR model for DGE (for now only one covariable allowed, for example a proven batch effect) 
		-T | -target_file # Protopical target file for attempts to differential gene expression analyses (containing filenames and covariates, automatically built if not provided)
		-S | -stop # Manual stop so the automatically downloaded files can be manually modified ('yes' or 'no', by default)
		-K | -Kraken2_fast_mode # Kraken2 fast mode, consisting on copying the Kraken2 database to /dev/shm (RAM) so execution is faster ('yes' or 'no' by default)
		-Dk | -kraken2_databases # Folder (absolute pathway) containing the database that should be used by Kraken2 (any input here, e.g. 'standard_eupathdb_48_kraken2_db', would activate the kraken2-based decontamination step)
		-Ds | -sortmerna_databases # The database (absolute pathway) that should be used by SortMeRNA (any input here, e.g. '/path/to/rRNA_databases/smr_v4.3_sensitive_db.fasta', would activate the sortmerna-based rRNA removal step)
		-De | -differential_expr_software # Software to be used in the differential expression analyses ('edgeR' by default, or 'DESeq2')
		-Df | -databases_function # Manually provide a comma separated list of databases to be used in automatic functional enrichment analyses of DEGs (check out the R package autoGO::choose_database(), but the most popular GO terms are used by default)
		-Dc | -deconvolution # Whether to perform deconvolution of the bulk RNA-seq data by CDSeq ('yes', which may require few hours to complete, or 'no', by default)
		-Dm | -debug_module # For debugging, step to remove the content of the corresponding folders and to resume a failed or incomplete run without repeating (one of 'step1', 'step1a', 'step1b', 'step1c', 'step2', 'step3a', 'step3b', 'step4', 'step5', 'step6', 'step7', 'step8', or 'all' to execute everything, by default)
		-Dec | -differential_expr_comparisons # Whether to restrict the differential expression analyses to only some of the possible comparisons (a comma-separated list of indexes pointing to the comparisons to keep, which you could get from a preliminar previous run, or 'no', by default)
		-Of | -options_featureCounts_feature # The feature type to use to count in featureCounts (default 'exon')
		-Os | -options_featureCounts_seq # The seqid type to use to count in featureCounts (default 'gene_name')
		-iG | -input_GEO_reads # If you want to combine downloading metadata from GEO with reads from GEO or any database already downloaded, maybe from a previous attempt, please provide an absolute path
		-cG | -compression_level # Specify the compression level to gzip the downloaded fastq files from GEO (numeric '0' to '9', default '9')
		-cP | -clusterProfiler # Whether to perform additional functional enrichment analyses using ClusterProfiler (slow if many significant DEGs or multiple number of comparisons, 'no' or 'yes', by default)
		-cPm | -clusterProfiler_method # Method for adjusting p.value in clusterprofiler iterations (one of 'holm','hochberg','hommel','bonferroni','BH','BY,'none', or 'fdr', by default)
		-Pm | -panther_method # Method for adjusting p.value in panther analyses via rbioapi (one of 'NONE','BONFERRONI', or 'FDR', by default)
		-Tc | -time_course_analyses # Whether to perform additional time-course analyses as a last step ('yes' or 'no', by default)
		-Tcsd | -time_course_std # Standard deviation threshold to filter in time course analyses (numeric, 1 by default)
		-Tcf | -time_course_fuzzi # Fuziness value for the soft clustering approach (by default an estimate is automatically computed but manual testing is encouraged)
		-Ti | -tidy_tmp_files # Space-efficient run, with a last step removing raw reads if downloaded, converting bam to cram, removing tmp files... etc ('yes' or 'no', by default)
		-Txls | -tables_in_xlsx # Convert all tables in results from .txt format, without limitation of size to Excel's .xlsx format, with a limitation of 32,767 characters ('yes' or 'no', by default)
		-Tx | -taxon_id # NCBI's taxon id of the organism, please not it is required for network analyses
		-Gt | -revigo_threshold_similarity # Similarity threshold for Revigo summaries of GO terms (0-1, suggested values are 0.9, 0.7, 0.5, 0.4 for large, medium, small, and tiny levels of similarity, respectively, being default 0.7
		-TMP | -TMPDIR # Directory to export the environmental variable TMPDIR (by default or if left empty an internal folder of the output directory is used, or please enter 'system' to use system's default, or an absolute pathway that will be created if it does not exist)
		-M | -memory_max # Max RAM memory to be used by aligner or JAVA in bytes (by default 257698037760, or 240GB, used)" && exit 1;;
		-i) input=${arguments[index]} ;;
		-n) name=${arguments[index]} ;;
		-o) output_folder=${arguments[index]} ;;
		-p) cores=${arguments[index]} ;;
		-pi) indexthreads=${arguments[index]} ;;
		-M) memory_max=${arguments[index]} ;;
		-s) strand=${arguments[index]} ;;
		-r) reference_genome=${arguments[index]} ;;
		-ri) reference_genome_index=${arguments[index]} ;;
		-g) genes=${arguments[index]} ;;
		-G) GSM_filter=${arguments[index]} ;;
		-Gt) rev_thr=${arguments[index]} ;;
		-a) annotation=${arguments[index]} ;;
		-f) filter=${arguments[index]} ;;
		-d) design_custom=${arguments[index]} ;;
		-D) design_custom_local=${arguments[index]} ;;
		-Df) databases_function=${arguments[index]} ;;
		-b) batch=${arguments[index]} ;;
		-bv) batch_vector=${arguments[index]} ;;
		-bc) batch_biological_covariates=${arguments[index]} ;;
		-C) covariables=${arguments[index]} ;;
		-S) stop=${arguments[index]} ;;
		-P) number_parallel=${arguments[index]} ;;
		-R) number_reads=${arguments[index]} ;;
		-t) transcripts=${arguments[index]} ;;
		-T) target=${arguments[index]} ;;
		-A) aligner=${arguments[index]} ;;
		-A) aligner_index_cache=${arguments[index]} ;;
		-K) kraken2_fast=${arguments[index]} ;;
		-Dk) kraken2_databases=${arguments[index]} ;;
		-Ds) sortmerna_databases=${arguments[index]} ;;
		-De) differential_expr_soft=${arguments[index]} ;;
		-Dm) debug_step=${arguments[index]} ;;
		-Dec) differential_expr_comparisons=${arguments[index]} ;;
		-Dc) deconvolution=${arguments[index]} ;;
		-cP) clusterProfiler=${arguments[index]} ;;
		-Of) optionsFeatureCounts_feat=${arguments[index]} ;;
		-O) organism_argument=${arguments[index]} ;;
		-Os) optionsFeatureCounts_seq=${arguments[index]} ;;
		-iG) input_geo_reads=${arguments[index]} ;;
		-cG) compression_level=${arguments[index]} ;;
		-Ti) tidy_tmp_files=${arguments[index]} ;;
		-Txls) convert_tables_excel=${arguments[index]} ;;
		-Tx) taxonid=${arguments[index]} ;;
		-Tc) time_course=${arguments[index]} ;;
		-Tcsd) minstd=${arguments[index]} ;;
		-Tcf) mestimate=${arguments[index]} ;;
		-TMP) TMPDIR_arg=${arguments[index]} ;;
		-q) qc_raw_reads=${arguments[index]} ;;
		-cPm) clusterProfiler_method=${arguments[index]} ;;
  		-Pm) panther_method=${arguments[index]} ;;
	esac
done

CURRENT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export PATH=$CURRENT_DIR/scripts:$PATH; echo -e "\n\nAdding to PATH the scripts folder..."
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
	echo -e "\nsalmon is going to be used to predict strandness. If not correct, please provide it with '-s'\n"
fi
if [ -z "$genes" ]; then
	genes="none"
else
	echo -e "\ngenes=$genes\n"
fi
if [ ! -z "$GSM_filter" ]; then
	echo -e "\nGSM_filter=$GSM_filter\n"
fi
if [ ! -z "$number_reads" ]; then
	echo -e "\nnumber_reads=$number_reads\n"
fi
if [ -z "$filter" ]; then
	filter=standard
fi
echo -e "\nfilter=$filter\n"
if [ -z "$batch" ]; then
	batch=no
fi
echo -e "\nbatch=$batch\n"
if [ -z "$covariables" ]; then
	covariables="none"
fi
echo -e "\ncovariables=$covariables\n"
if [ -z "$deconvolution" ]; then
	deconvolution=no
fi
if [ -z "$design_custom" ]; then
	design_custom=no
fi
echo -e "\ndesign_custom=$design_custom\n"
if [ -z "$stop" ]; then
	stop=no
fi
echo -e "\nstop=$stop\n"
if [ -z "$memory_max" ]; then
	memory_max=257698037760
fi
echo -e "\nmemory_max=$memory_max\n"
if [ -z "$miarma_path" ]; then
	miarma_path=$CURRENT_DIR/external_software/miARma-seq
fi
if [ -d "$output_folder/$name" ]; then
	echo -e "Please note that $output_folder/$name already exists... reanalyzerGSE is going to attempt a new run or resume running, but you may want to remove the folder, change the destination folder with '-o' or '-n', use downloaded raw data from an external software... etc. Sleeping for a while to give you time to exit if you want, and then continuing..."; sleep 60
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
if [ -z "$indexthreads" ]; then
	indexthreads=$cores
fi
if [ -z "$clusterProfiler" ]; then
	clusterProfiler="yes"
fi
if [ -z "$time_course" ]; then
	time_course="no"
fi
if [ -z "$convert_tables_excel" ]; then
	convert_tables_excel="no"
fi
if [ -z "$mestimate" ]; then
	mestimate=0
fi
if [ -z "$minstd" ]; then
	minstd=1
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
if [ -z "$debug_step" ]; then
	debug_step="all"
fi
export debug_step
### Info on the kraken2 decontamination step and databases:
if [ -z "$kraken2_fast" ]; then
	kraken2_fast="no"
fi
if [ -z "$compression_level" ]; then
	compression_level=9
fi
if [ -z "$rev_thr" ]; then
	rev_thr=0.7
fi
echo -e "\nCompression_level raw reads=$compression_level\n"
seqs_location=$output_folder/$name/raw_reads

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

		echo -e "\n\nAlternatively, you can download the already built database: cd $kraken2_databases && wget --no-check-certificate" '"bit.ly/kraken2_dbs" -O standard_eupathdb_48_kraken2_db.tar.xz && tar -xf standard_eupathdb_48_kraken2_db.tar.xz && rm standard_eupathdb_48_kraken2_db.tar.xz'
		echo -e "\n\nAlternatively, other ready-to-use databases can be downloaded from https://benlangmead.github.io/aws-indexes/k2"
		echo -e "Please note that if you followed installation instructions, the kraken2-build script has been modified so database building is faster (improved masking, https://github.com/DerrickWood/kraken2/pull/39) and a bug in sequences download from ncbi has been corrected (https://github.com/DerrickWood/kraken2/issues/571)"
		echo -e "Please note that if you have allocated enough RAM and the system is compatible, copying the kraken2 database to the faster RAM, something like /dev/shm/, and then pointing to the library there in kraken2 execution with the flag --memory-mapping would greatly improve speed, particularly if multiple runs (https://github.com/DerrickWood/kraken2/issues/451). reanalyzerGSE includes this mode with the flag -K"
		echo -e "The execution of kraken2 classification with the suggested database (archaea + viral + plasmid + bacteria + fungi + UniVec_Core + EuPathDB) will require ~70GB of RAM."
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



###### STEP 1. Download info from GEO and organize, or start processing the raw fastq.gz:
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

if [[ $debug_step == "all" || $debug_step == "step1" ]]; then
	echo -e "\n\nSTEP 1: Starting...\nCurrent date/time: $(date)\n\n"
	if [[ $input == G* ]]; then
	### Download info:
		echo -e "\nDownloading info from GEO for $input...\n"
		R_download_GEO_info.R $input $output_folder
		arrIN=(${input//,/ })
		input=$(for a in "${arrIN[@]}"; do echo "$a"; done | sort | tr '\n' '_' | sed 's,_$,,g')
	
	### Get metadata and process the info:
		echo -e "\nProcessing and downloading more data...\n"
		cd $output_folder/$name/GEO_info
		if [ ! -z "$GSM_filter" ]; then
			echo $GSM_filter > gsm_manual_filter.txt
		fi
		series_matrix=$(zcat *_series_matrix.txt.gz | grep SRP | sed 's,.*SRP,SRP,g' | sed 's,",,g')
		if [ ! -z "$series_matrix" ]; then
			for f in $series_matrix; do
				pysradb metadata $f --detailed --assay --desc --expand --saveto sample_info_pysradb_$f.txt
			done
			for i in $(zcat *_series_matrix.txt.gz | egrep 'SRP|Series_geo_accession' | sed 's,.*SRP,SRP,g' | sed 's,",,g' | sed 's,.*\t,,g' | grep SRP); do
				mv $(ls | grep $i) $(ls | grep $i | sed 's,.txt,_,g')$(zcat *_series_matrix.txt.gz | egrep 'SRP|Series_geo_accession' | sed 's,.*SRP,SRP,g' | sed 's,",,g' | sed 's,.*\t,,g' | grep -B1 $i | grep GSE)".txt"
			done
		R_process_pysradb.R $input $output_folder
		fi
		R_download_GEO_info_process.R $input $output_folder
		sed -i -r -e 's/[^[:alnum:]_\t]/_/g' -e 's/_\+/_/g' -e 's/(_[^_]*)\1+/\1/g' sample_names.txt
		if test -f srx_ids.txt; then
			for i in $(cat srx_ids.txt); do esearch -db sra -query $i | esummary | xtract -pattern DocumentSummary -element Run@acc >> srr_ids.txt; done
			rm srx_ids.txt
		fi
		if test -f srr_ids.txt; then
			for i in $(cat srr_ids.txt | shuf | head -2); do esearch -db sra -query $i | esummary | egrep "SINGLE|PAIRED" >> library_layout.txt; done && uniq library_layout.txt | sed 's/          //g;s/<//g;s/>//g;s\/\\g' > library_layout_info.txt && rm library_layout.txt
			paste -d$'\t' srr_ids.txt sample_names.txt $(for f in $(ls | grep full); do echo $f" "$(sort $f | uniq -c | wc -l); done | grep -v " 1" | cut -d" " -f1 | head -1) > samples_info.txt # the third column is a design containing at least more than one element...
			sed -i -r -e 's/[^[:alnum:]_\t]/_/g' -e 's/_\+/_/g' -e 's/(_[^_]*)\1+/\1/g' -e 's/_([0-9]+)/-\1/g' samples_info.txt # Should be redundant but make sure to remove special characters from the sample names and _1/_2... it's crucial for later steps such as fastqc and miarma-seq
		fi
		if [ ! -s srr_ids.txt ]; then
			echo -e "\nI haven't been able to find SRR accession ids to download the sequences and I'm exiting, please double check manually..."; exit 1
		fi
		echo -e "\nAll available info downloaded from GEO, please check it out in $output_folder/$name/GEO_info\n"
		if [ "$design_custom" == "yes" ]; then
			echo "You have requested to manually provide the experimental design instead of the ones shown above. This is the list of samples:"
			cat sample_names.txt
			echo "Please provide a comma-separated list with the conditions for each sample:"
			read -r design_input
			rm $(ls -d $output_folder/$name/GEO_info/* | grep 'design_possible_')
			echo $design_input | sed 's/,/\n/g' > $output_folder/$name/GEO_info/design_possible_full_1.txt
			cat $output_folder/$name/GEO_info/design_possible_full_1.txt | sort | uniq > $output_folder/$name/GEO_info/design_possible_1.txt
		fi
	
	### Stop and continue with other script if it's a single-cell:
		if [ $(zcat $output_folder/$name/GEO_info/*_series_matrix.txt.gz | egrep -e 'single nuclei|single cell|single-cell|snRNA|scRNA' | wc -l) -gt 0 ] || [ $(zcat $(find . -name "*_series_matrix.txt.gz") | egrep -i -e 'single nuclei|single cell|single-cell|snRNA|scRNA' | wc -l) -gt 0 ]; then
		  	echo -e "\n\nDetected this could be a single-cell RNA-seq study... I can try to do stuff automatically (i.e. try and normalize the raw counts or give an estimated bulk expression taking the average), but errors are expected. \nThe script 'R_process_reanalyzer_GSE_single_cell.R' is a template built from the case example GSE118257, and valid to other GEO entries where pheno data and matrix counts are supplementary files clearly named.\nHowever, manual changes are most likely required to work with other studies... These changes should be possible, so please open an issue or go for it if you have the expertise and this one fails!. For example, it's likely that it's just required to point to the directory of the matrix counts, or to manually specify the columns/names of the conditions/cells\n\n"
			# Print the results to review:
			echo -e "This text in the metadata is what made the pipeline to suggest this could be single-cell:"
			zcat $output_folder/$name/GEO_info/*_series_matrix.txt.gz | egrep -e 'single nuclei|single cell|single-cell|snRNA|scRNA'
			zcat $(find . -name "*_series_matrix.txt.gz") | egrep -i -e 'single nuclei|single cell|single-cell|snRNA|scRNA'
			# Choice:
			echo -e "\nWrite 'yes' to continue with single-cell analyses, or 'no' to continue with normal analyses after reviewing the entry and the statements in the metadata pointing to single-cell..."
			read -r single_cell_choice
			if [ "$single_cell_choice" == "yes" ]; then
				R_process_reanalyzer_GSE_single_cell.R $input $output_folder $genes; exit 1
			elif [ "$single_cell_choice" == "no" ]; then
				echo -e "\nContinuing with bulk RNA-seq analyses...\n"
			fi
		fi
	
	### Stop and continue with other script if it's a microarrays:
		if [ $(zcat $output_folder/$name/GEO_info/*_series_matrix.txt.gz | egrep -i -e 'Expression profiling by array|microarray' | wc -l) -gt 0 ] || [ $(zcat $(find . -name "*_series_matrix.txt.gz") | egrep -i -e 'Expression profiling by array|microarray' | wc -l) -gt 0 ]; then
		  	echo -e "\n\nDetected this could be a microarrays study... trying to do analyze automatically, but errors in this log are expected. The script 'R_process_reanalyzer_GSE_microarrays.R' is already supporting the most frequent arrays and platforms, but it could require to be extended in order to work with other studies... These changes should be possible though, so please open an issue or go for it if you have the expertise and this one fails!\n\n"
			# Choice:
			echo -e "\nWrite 'yes' to continue with microarrays analyses, or 'no' to continue with normal analyses after reviewing the entry pointing to microarrays..."
			read -r microarrays_choice
			if [ "$microarrays_choice" == "yes" ]; then
			  	R_process_reanalyzer_GSE_microarrays.R $input $output_folder $genes; exit 1
			elif [ "$microarrays_choice" == "no" ]; then
				echo -e "\nContinuing with bulk RNA-seq analyses...\n"
			fi
		fi
	
	### Stop if SRR not obtained and not single-cell or microarrays
		if [ ! -s srr_ids.txt ]; then
			echo -e "\nI haven't been able to find SRR accession ids to download the sequences and I'm exiting, please double check manually..."; exit 1
		fi
	
	### Get organism:
		organism=$(zcat $output_folder/$name/GEO_info/*_series_matrix.txt.gz | grep "organism" | awk '{$1=""}1' |tr '"' '\n' | sort -u | sed -r '/^\s*$/d')
		echo $organism > $output_folder/$name/GEO_info/organism.txt
		if [ $(zcat $output_folder/$name/GEO_info/*_series_matrix.txt.gz | grep "organism" | awk '{$1=""}1' |tr '"' '\n' | sort -u | sed -r '/^\s*$/d' | wc -l) -gt 1 ]; then
			echo -e "\n Please keep in mind that two different organisms are detected. You are likely requesting an analysis combining multiple GSEXXXXX, please make sure they are from the same organism. Another possibility is there are multiple series_matrix within the same GSEXXXX id, and you may have requested to stop and manually clarify. Continuing with organism: "
			organism=$(zcat *_series_matrix.txt.gz | grep "organism" | awk '{$1=""}1' |tr '"' '\n' | sort -u | sed -r '/^\s*$/d' | head -1)
			echo $organism > $output_folder/$name/GEO_info/organism.txt
			echo -e "$organism\nPlease request on the next run a stop with parameter '-S' and modify manually the file GEO_info/organism.txt if not required...\n"
		fi
		echo -e "\nSTEP 1 DONE. Current time: $(date)\n"
	fi
	export debug_step="all"
	echo -e "\n\nSTEP 1: DONE\nCurrent date/time: $(date)\n\n"
fi


###### STEP 1a. Download fastq files from the GEO ID provided:
if [[ $debug_step == "all" || $debug_step == "step1a" ]]; then
	echo -e "\n\nSTEP 1a: Starting...\nCurrent date/time: $(date)\n\n"
	if [[ $input == G* ]]; then		
		if [ "$stop" == "yes" ]; then
			echo "You have requested a stop to manually provide the SRR ids, or potentially modify other files that may have not been detected properly from GEO, and were not correct, or you just want to adapt some of them. Please double check or manually modify the files GEO_info/srr_ids.txt, samples_info.txt, sample_names.txt, phenodata_extracted.txt, library_layout_info.txt, organism.txt, design_files, etc. The pipeline is stopped. Please press space to continue or Ctrl + C to exit..."
			read -n1 -s -r -p $'Press space to continue...\n' key
			rm $output_folder/$name/possible_designs_all.txt
			for i in $(ls $output_folder/$name/GEO_info| grep "full"); do echo $i >> $output_folder/$name/possible_designs_all.txt && cat $output_folder/$name/GEO_info/$i >> $output_folder/$name/possible_designs_all.txt && echo -e "\n" >> $output_folder/$name/possible_designs_all.txt; done
			cat $output_folder/$name/GEO_info/phenodata_extracted.txt > $output_folder/$name/phenotypic_data_samples.txt
		fi
		if [ ! -d "$seqs_location" ]; then
			mkdir -p $seqs_location
			echo "Downloading the fastq files from SRR..."		
			if [ -z "$input_geo_reads" ]; then
				download_sra_fq.sh $output_folder/$name/GEO_info/srr_ids.txt $seqs_location $(( number_parallel*2 )) $cores $compression_level		
	### Rename the fastq files (max length name 140 characters)	and handle already downloaded datasets if provided:		
				cd $seqs_location
	   			if [[ "$(cat $output_folder/$name/GEO_info/library_layout_info.txt)" == "SINGLE" ]]; then
					for i in $(cat $output_folder/$name/GEO_info/srr_ids.txt); do echo "mv $(ls | egrep ^$i | head -1) $(cat $output_folder/$name/GEO_info/samples_info.txt | grep $i | cut -f 2 | sed -e 's,%,,g;s,(,,g;s,),,g;s/[_]1/1/g;s/[_]2/2/g;s/replicate_/replicate/g' | awk -F '_GSM' '{ gsub(/-/,"",$1); print substr($1, 1, 140) "_GSM" $2 }')""_1.fastq.gz" && mv $(ls | egrep ^$i | head -1) $(cat $output_folder/$name/GEO_info/samples_info.txt | grep $i | cut -f 2 | sed -e 's,%,,g;s,(,,g;s,),,g;s/[_]1/1/g;s/[_]2/2/g;s/replicate_/replicate/g' | awk -F '_GSM' '{ gsub(/-/,"",$1); print substr($1, 1, 140) "_GSM" $2 }')"_1.fastq.gz"; done
				elif [[ "$(cat $output_folder/$name/GEO_info/library_layout_info.txt)" == "PAIRED" ]]; then
					for i in $(cat $output_folder/$name/GEO_info/srr_ids.txt); do echo "mv $(ls | egrep ^$i | head -1) $(cat $output_folder/$name/GEO_info/samples_info.txt | grep $i | cut -f 2 | sed -e 's,%,,g;s,(,,g;s,),,g;s/[_]1/1/g;s/[_]2/2/g;s/replicate_/replicate/g' | awk -F '_GSM' '{ gsub(/-/,"",$1); print substr($1, 1, 140) "_GSM" $2 }')""_1.fastq.gz" && echo "mv $(ls | egrep ^$i | tail -1) $(cat $output_folder/$name/GEO_info/samples_info.txt | grep $i | cut -f 2 | sed -e 's,%,,g;s,(,,g;s,),,g;s/[_]1/1/g;s/[_]2/2/g;s/replicate_/replicate/g' | awk -F '_GSM' '{ gsub(/-/,"",$1); print substr($1, 1, 140) "_GSM" $2 }')""_2.fastq.gz" && mv $(ls | egrep ^$i | head -1) $(cat $output_folder/$name/GEO_info/samples_info.txt | grep $i | cut -f 2 | sed -e 's,%,,g;s,(,,g;s,),,g;s/[_]1/1/g;s/[_]2/2/g;s/replicate_/replicate/g' | awk -F '_GSM' '{ gsub(/-/,"",$1); print substr($1, 1, 140) "_GSM" $2 }')"_1.fastq.gz" && mv $(ls | egrep ^$i | tail -1) $(cat $output_folder/$name/GEO_info/samples_info.txt | grep $i | cut -f 2 | sed -e 's,%,,g;s,(,,g;s,),,g;s/[_]1/1/g;s/[_]2/2/g;s/replicate_/replicate/g' | awk -F '_GSM' '{ gsub(/-/,"",$1); print substr($1, 1, 140) "_GSM" $2 }')"_2.fastq.gz"; done
				fi
			else
				echo -e "\nSoft linking the already downloaded raw reads from the provided directory: $input_geo_reads\n"
				ln -sf $input_geo_reads/* $seqs_location
			fi		
		fi
		echo -e "\nDONE. Current date/time: $(date)"; time1=`date +%s`; echo -e "Elapsed time (secs): $((time1-start))"; echo -e "Elapsed time (hours): $(echo "scale=2; $((time1-start))/3600" | bc -l)\n"
	
	### Process if any download was not successful or subsampling if required:
		cd $seqs_location
		num_gz_files=$(find . -name "*.gz" | wc -l)
		num_samples=$(cat $output_folder/$name/GEO_info/sample_names.txt | wc -l)
		if [ "$num_gz_files" -eq "$(($num_samples * 2))" ] || [ "$num_gz_files" -eq "$num_samples" ]; then
			echo -e "\nPlease double check this order is the same than the rest of lists printed in the log, and that the info, e.g., the correspondence between GSM and SRR, is correct:"
			cat $output_folder/$name/GEO_info/samples_info.txt
		else
			echo -e "\nRaw reads not downloaded fully? Please double check manually the log files and the folder $seqs_location to assess whether there have been errors with downloading. Retrying all downloads... with another approach\n"
			echo -e "\nIn the future this will automatically detect and only resume the downloads that fail...\n"
			seqs_location=$output_folder/$name/raw_reads
			number_ids=$(echo $input | tr ',' '\n' | wc -l)
			if [ $number_ids -le $number_parallel ]; then
				export cores_parallel=$((cores / number_files))
			else
				export cores_parallel=$((cores / number_parallel))
			fi
			cd $seqs_location; rm -rf *
			echo $input | tr ',' '\n' | parallel --joblog $output_folder/$name/fastq_dl_log_parallel.txt -j $number_parallel --max-args 1 'if [ $(echo {} | egrep -c "PRJEB|PRJNA|PRJDB|ERX|DRX|SRX|ERP|DRP|SRP") -eq 1 ]; then fastq-dl --cpus $cores_parallel --accession {}; fi && 
		 																		   if [ $(echo {} | egrep -c "ERS|DRS|SRS|SAMD|SAME|SAMN|ERR|DRR|SRR") -eq 1 ]; then fastq-dl --provider sra --cpus $cores_parallel --accession {}; fi'
		 	if [ "$num_gz_files" -eq "$(($num_samples * 2))" ] || [ "$num_gz_files" -eq "$num_samples" ]; then
		 		echo "It seems the download has been sucessful, but please double check"
		 	else
		 		echo "Download still failed. Please double check manually, exiting..."; exit 1
		 	fi	 																		   	
		fi		
			
		if [ -z "$number_reads" ]; then
			echo -e "\nAll raw data downloaded and info prepared, proceeding with reanalyses...\n"
		else
			echo -e "\nSubsampling...\n"
			# From the input parameter by the user, obtain a random number allowing a +- 10% window:
			IFS=', ' read -r -a arr <<< "$number_reads"
			IFS=', ' read -r -a arr2 <<< "$(ls | egrep .fastq.gz$ | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq | tr '\n' ',')"
			export -f subsample_reads
			subsample_reads() {
								file=$1
								number=$2							
								ten_percent=$(( number * 10 / 100 ))
								random_shift=$((RANDOM % (2 * ten_percent + 1) - ten_percent))
								number_reads_rand=$((number + random_shift))
								echo "$file to $number +- 10%... to $number_reads_rand"
								seqtk sample -s 123 "$file" "$number_reads_rand" | pigz -p $((cores / 4)) -c --best > "${file}_subsamp"
							  }		
			parallel --verbose -j $cores subsample_reads {} ::: "${arr2[@]}" ::: "${arr[@]}"		
			rm $(ls | grep -v subsamp); for file in $(ls); do mv $file $(echo $file | sed 's,_subsamp,,g'); done
			echo -e "\nAll raw data downloaded and info prepared, subsampling (+-10%) completed. Proceeding with reanalyses...\n"
		fi
	fi
	export debug_step="all"
	echo -e "\n\nSTEP 1a: DONE\nCurrent date/time: $(date)\n\n"
fi

 
### STEP 1b. Process if not required to download from NCBI/GEO but raw reads provided locally:
if [[ $debug_step == "all" || $debug_step == "step1b" ]]; then
	echo -e "\n\nSTEP 1b: Starting...\nCurrent date/time: $(date)\n\n"
	if [[ $input == /* ]]; then
		seqs_location=$output_folder/$name/raw_reads
		if [ ! -d "$seqs_location" ]; then
			mkdir -p $seqs_location		
			if [ $(ls -d $input/* | egrep -c "_R1.fastq.gz$|_R1.fq.gz$|_R2.fastq.gz$|_R2.fq.gz$|_1.fastq.gz$|_1.fq.gz$|_2.fastq.gz$|_2.fq.gz$") -eq 0 ]; then
				echo -e "\nPlease make sure that the input files are named _1.fastq.gz, _R1.fastq.gz, _2.fastq.gz, _R2.fastq.gz\n"
				exit 1
			fi
			for f in $(ls -d $input/*); do ln -sf $f $seqs_location/$(basename $f | sed 's,fq,fastq,g;s,_R1.fastq,_1.fastq,g;s,_R2.fastq,_2.fastq,g'); done
			echo -e "\nProcessing the provided fastq files, renaming to _1.fastq and _2.fastq if necessary...\n"
		fi
	 	cd $seqs_location
		if [[ -z `find $output_folder/$name -name library_layout_info.txt` ]]; then
			if [ $(ls $seqs_location | egrep -c "_1.fastq$|_1.fq$|_R1.fastq$|_R1.fq$|_1.fastq.gz$|_1.fq.gz$|_R1.fastq.gz$|_R1.fq.gz$") -gt 0 ]; then
				echo "SINGLE" > $output_folder/$name/library_layout_info.txt
				if [ $(ls $seqs_location | egrep -c "_2.fastq$|_2.fq$|_R2.fastq$|_R2.fq$|_2.fastq.gz$|_2.fq.gz$|_R2.fastq.gz$|_R2.fq.gz$") -gt 0 ]; then
					echo "PAIRED" > $output_folder/$name/library_layout_info.txt
				fi
			else
				echo "SINGLE" > $output_folder/$name/library_layout_info.txt
			fi
		fi
		if [ -z "$number_reads" ]; then
			echo -e "\nAll raw data downloaded and info prepared, proceeding with reanalyses...\n"
		else
			echo -e "\nSubsampling...\n"
			# From the input parameter by the user, obtain a random number allowing a +- 10% window:
			IFS=', ' read -r -a arr <<< "$number_reads"
			IFS=', ' read -r -a arr2 <<< "$(ls | egrep .fastq.gz$ | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq | tr '\n' ',')"
			for index in "${!arr[@]}"; do		
				number=${arr[index]}
				echo "To $number +- 10%..."
				ten_percent=$(( number * 10 / 100 ))
				random_shift=$((RANDOM % (2 * ten_percent + 1) - ten_percent))
				number_reads_rand=$((number + random_shift))
				echo "... to $number_reads_rand"
				if [[ "$(cat $output_folder/$name/library_layout_info.txt)" == "PAIRED" ]]; then
					export SEQKIT_THREADS=$((cores / 2)); ls | egrep .fastq.gz$ | grep ${arr2[index]} | parallel --verbose -j 2 "seqtk sample -s 123 {} $number_reads_rand | pigz -p $((cores / 2)) -c --fast > {}_subsamp"
				else
	  				export SEQKIT_THREADS=$cores; seqtk sample -s 123 $(ls | egrep .fastq.gz$ | grep ${arr2[index]}) $number_reads_rand | pigz -p $cores -c --fast > $(ls | egrep .fastq.gz$ | grep ${arr2[index]})_subsamp
	  			fi
			done
			rm $(ls | grep -v subsamp); for file in $(ls); do mv $file $(echo $file | sed 's,_subsamp,,g'); done
			echo -e "\nAll raw data downloaded and info prepared, subsampling (+-10%) completed. Proceeding with reanalyses...\n"
		fi
	 	echo -e "This is the content of $seqs_location:\n$(ls -l $seqs_location | awk '{ print $9 }' | tail -n +2)\n"
		if [ -z "$design_custom_local" ]; then
			echo -n "From the ordered list above, please input a comma-separated list with the conditions for each sample. Remember to try and avoid complex names, use as few underlines as possible, avoid names starting with numbers or others that would not be sorted appropriately such as containing spaces, and if reads are paired-end, only one name of condition per pair of reads: "
			read -r design_input
		else
			design_input=$design_custom_local
			echo -e "The used conditions are:\n$(echo $design_input | sed 's_,_\n_g')"
		fi
		mkdir -p $output_folder/$name/GEO_info/
		paste <(ls $seqs_location | egrep '.fq|.fastq' | sed "s/_1.fastq.gz//" | sed "s/_2.fastq.gz//" | uniq) <(paste -d'_' <(ls $seqs_location | egrep '.fq|.fastq' | egrep '.fq|.fastq' | sed "s/_1.fastq.gz//" | sed "s/_2.fastq.gz//" | uniq) <(echo $design_input | sed 's/,/\n/g')) <(echo $design_input | sed 's/,/\n/g') > $output_folder/$name/GEO_info/samples_info.txt
		echo $design_input | sed 's/,/\n/g' > $output_folder/$name/GEO_info/design_possible_full_1.txt
		cat $output_folder/$name/GEO_info/design_possible_full_1.txt | sort | uniq > $output_folder/$name/GEO_info/design_possible_1.txt
		echo $name > $output_folder/$name/GEO_info/study_title.txt
		if [ -z "$organism_argument" ]; then
			echo -n "Please input the scientific name of the organism: "
			read -r organism
		else
			organism=$(echo $organism_argument | sed 's,_, ,g')
			echo "Organism used is $organism"
		fi
		echo $organism > $output_folder/$name/GEO_info/organism.txt
	fi
	export debug_step="all"
	echo -e "\n\nSTEP 1b: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP 1c. Process if required to download from manually provided ids from databases
if [[ $debug_step == "all" || $debug_step == "step1c" ]]; then
	echo -e "\n\nSTEP 1c: Starting...\nCurrent date/time: $(date)\n\n"
	if [[ $input == P*  || $input == E* || $input == D* || $input == S* ]]; then
		seqs_location=$output_folder/$name/raw_reads
		number_ids=$(echo $input | tr ',' '\n' | wc -l)
		if [ $number_ids -le $number_parallel ]; then
			export cores_parallel=$((cores / number_files))
		else
			export cores_parallel=$((cores / number_parallel))
		fi
		if [ ! -d "$seqs_location" ]; then
			mkdir -p $seqs_location; cd $seqs_location
			echo -e "\nDownloading from the input accessions that you manually provided...\n"
			echo $input | tr ',' '\n' | parallel -j $number_parallel --max-args 1 'if [ $(echo {} | egrep -c "PRJEB|PRJNA|PRJDB|ERX|DRX|SRX|ERP|DRP|SRP") -eq 1 ]; then fastq-dl --cpus $cores_parallel --accession {}; fi && 
		 																		   if [ $(echo {} | egrep -c "ERS|DRS|SRS|SAMD|SAME|SAMN|ERR|DRR|SRR") -eq 1 ]; then fastq-dl --provider sra --cpus $cores_parallel --accession {}; fi'
		fi
	fi
	export debug_step="all"
	echo -e "\n\nSTEP 1c: DONE\nCurrent date/time: $(date)\n\n"
fi


### Deal with batch correction... The user has to use certain arguments to manually provide a list or do it interactively:
if [[ $debug_step == "all" ]]; then
	if [ "$batch" == "yes" ]; then
		if [ -z "$batch_vector" ]; then
			echo -e "This is the content of $seqs_location:\n$(ls -l $seqs_location | awk '{ print $9 }' | tail -n +2)\n"
			echo -n "Based on the list above, please input a comma-separated list for the vector for batch separation (use only numbers, and if these are paired-end, only once per pair of reads):"
			read -r batch_vector
		fi
		echo $batch_vector > $output_folder/$name/GEO_info/batch_vector.txt
		echo -e "\nThe comma-separated list for the vector for batch separation is $batch_vector\n"
		if [ -z "$batch_biological_covariates" ]; then
			echo -n "Please input a comma-separated list for the biological covariate, and separate by space if multiple biological variables are to be included (use only numbers): "
			read -r batch_biological_covariates		
		fi
		echo $batch_biological_covariates > $output_folder/$name/GEO_info/batch_biological_variables.txt
		echo -e "\nThe comma-separated list for the vector of biological covariable for batch separation is $batch_biological_covariates\n"
	fi

	### Give info of NCBI's current genome:
	Rscript -e "genomes <- rentrez::entrez_summary(db='genome', id=rentrez::entrez_search(db='genome', term='${organism}[orgn]')\$ids);cat(paste(paste0('\n\nNCBI current genome info: ', date()),genomes\$assembly_name,genomes\$assembly_accession,genomes\$create_date,'\n',sep='\n'))"
	organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's/ \+/_/g;s/__*/_/g') # Get again the organism in case it has been manually modified... and without spaces...
	export debug_step="all"
fi


### STEP 2. Decontamination if required:
if [[ $debug_step == "all" || $debug_step == "step2" ]]; then
	echo -e "\n\nSTEP 2: Starting...\nCurrent date/time: $(date)\n\n"
	if [ ! -z "$kraken2_databases" ]; then
		mkdir -p $output_folder/$name/raw_reads_k2; cd $output_folder/$name/raw_reads_k2
		if [[ $kraken2_fast == "yes" ]]; then
			echo -e "\nPreparing database $kraken2_databases for fast access in RAM...\n"
			cp -ru $kraken2_databases /dev/shm/
			db_k2="--db /dev/shm/$(basename $kraken2_databases) --memory-mapping"
		elif [[ $kraken2_fast == "no" ]]; then
			db_k2="--db $kraken2_databases"
		fi
		echo -e "\nExecuting kraken2...\n"
		if [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "SINGLE" ]]; then
			for f in $(ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$"); do \time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" kraken2 $db_k2 --threads $cores --classified-out $f.classification.txt --unclassified-out $f.classification_unknwn.txt --report $f.report.txt --output $f.kraken2_output.txt --use-names $seqs_location/$f 1>> kraken2_log_out.txt 2>> kraken2_log_warnings_errors.txt; done
		elif [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "PAIRED" ]]; then
			for f in $(ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$" | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq); do \time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S .user=%U" kraken2 $db_k2 --threads $cores --classified-out $f.classification.txt --unclassified-out $f.classification_unknwn.txt --report $f.report.txt --output $f.kraken2_output.txt --use-names --paired $seqs_location/${f}_1.fastq.gz $seqs_location/${f}_2.fastq.gz 1>> kraken2_log_out.txt 2>> kraken2_log_warnings_errors.txt; done
		fi
		echo -e "Please check the files report.txt, kraken2_log_out.txt and kraken2_log_out_warnings_errors.txt"
		echo -e "Processing reports and extracting uncontaminated reads..."
		for f in $(ls | egrep "report.txt$"); do echo -e "\nLog of kraken2:"; echo -e "%_reads_covered\t#_reads_covered\t#_reads_directly_assigned\tRank_code\tTaxon_id\tScientific_name" >> $f.final.txt && cat $f >> $f.final.txt && echo -e "\n\nNumber of classified reads at the genus level: $(cat $f | awk '$4 == "G" {print $2"\t"$5}' | awk '{s+=$1}END{print s}')" >> $f.final.txt && echo -e "\nTaxonomy IDs at the genus level assigned to the reads:" >> $f.final.txt && echo -e "#read\tTaxID\n$(awk '$4 == "G" {print $2"\t"$5}' $f)\n" >> $f.final.txt; done
		for f in $(ls | egrep "kraken2_output.txt$"); do rcf -n $kraken2_databases/taxdump -k $f -o $f.recentrifuge_contamination_report.html -e CSV &>> rcf_log_out.txt; done # Add --sequential if problems with multithreading
		
		if [ -z "$taxonid" ]; then
			taxonid=$(echo $organism | sed 's/_\+/ /g' | taxonkit name2taxid --data-dir $kraken2_databases/taxdump | head -1 | cut -f2)
		fi	
		taxon_name=$(taxonkit list --ids $taxonid -n -r --data-dir $kraken2_databases/taxdump | grep $taxonid)
		echo -e "\nOrganism provided: $organism"; echo -e "\nOrganism provided (taxonid): $taxonid"; echo $taxon_name
		echo -e "\nIf not correct, please rerun and double check that you have provided it explicitely in the prompt... kraken2 output will be filtered to retain that taxa and below"
		echo -e "\nCheck out the logs in the files rcf_log_out.txt and extract_kraken2_log_out.txt"
		
		mkdir -p $seqs_location\_k2
		if [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "SINGLE" ]]; then
			for f in $(ls | egrep ".kraken2_output.txt$"); do extract_kraken_reads.py -k $f -U $(echo $f | sed 's,.kraken2_output.txt,,g') -o $seqs_location\_k2/$f\_1.fastq.gz -t $taxonid -r $(echo $f | sed 's,.kraken2_output.txt,,g').report.txt --include-children &>> extract_kraken2_log_out.txt; done
		elif [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "PAIRED" ]]; then
			for f in $(ls | egrep ".kraken2_output.txt$"); do extract_kraken_reads.py -k $f -s1 $(echo $f | sed 's,.kraken2_output.txt,,g')\_1.fastq.gz -s2 $(echo $f | sed 's,.kraken2_output.txt,,g')\_2.fastq.gz -o $seqs_location\_k2/$f\_1.fastq.gz -o2 $seqs_location\_k2/$f\_2.fastq.gz -t $taxonid -r $(echo $f | sed 's,.kraken2_output.txt,,g').report.txt --include-children &>> extract_kraken2_log_out.txt; done
		fi
		for f in $(ls | grep "k2" | egrep ".fastq.gz$"); do fastqc -t $cores $f; done
	
	fi
	if [ ! -z "$sortmerna_databases" ]; then
		if [ ! -d "$CURRENT_DIR/indexes/$(basename $sortmerna_databases)_sortmerna_index" ]; then
			sortmerna --index 1 --ref $sortmerna_databases --workdir $CURRENT_DIR/indexes/$(basename $sortmerna_databases)_sortmerna_index --threads $cores
		fi
		mkdir -p $seqs_location\_sortmerna; cd $seqs_location\_sortmerna
		# with the argument --paired_out, only the pairs where both reads are coincident (aligning to rRNA or not, are included in the results)
		# I don't include it, so the numbers are exactly the ones in the log, and the properly paired reads can be dealt with later on the mapping
		echo -e "\nExecuting sortmerna and fastqc of the new reads...\n"
		if [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "PAIRED" ]]; then
			for f in $(ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$" | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq); do sortmerna --idx-dir $CURRENT_DIR/indexes/$(basename $sortmerna_databases)_sortmerna_index/idx --ref $sortmerna_databases --reads $seqs_location/${f}_1.fastq.gz --reads $seqs_location/${f}_2.fastq.gz --workdir ${f}_sortmerna_out --fastx --threads $cores --out2 --aligned ${f}_rRNA --other ${f}_no_rRNA -v; done
			rm -rf $(ls | egrep "_sortmerna_out$")
		elif [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "SINGLE" ]]; then
			for f in $(ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$"); do sortmerna --idx-dir $CURRENT_DIR/indexes/$(basename $sortmerna_databases)_sortmerna_index/idx --ref $sortmerna_databases --reads $seqs_location/$f --workdir ${f}_sortmerna_out --fastx --threads $cores --aligned ${f}_rRNA --other ${f}_no_rRNA -v; done
			rm -rf $(ls | egrep "_sortmerna_out$")
		fi
		for f in $(ls | grep "rRNA" | egrep ".fastq|.fq"); do fastqc -t $cores $f; done
	fi	
	export debug_step="all"
	echo -e "\n\nSTEP 2: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP3a. Prepare the data and info for running miARma-seq:
if [[ $debug_step == "all" || $debug_step == "step3a" ]]; then
	echo -e "\n\nSTEP 3a: Starting...\nCurrent date/time: $(date)\n\n"
	### Prepare the salmon index from the trancripts sequences if required and strandness prediction... (if the miarma0.ini does not exist yet, pointing to a previous miarma run)
	if [[ ! -e "$output_folder/$name/miarma0.ini" ]]; then
		if [ -z "$strand" ]; then
			echo -e "\nLooking for indexes or indexing transcripts in $transcripts...\n"
			mkdir -p $output_folder/$name/indexes
			salmon_idx=$CURRENT_DIR/indexes/${organism}_salmon_idx
			if [ ! -d "$salmon_idx" ]; then
				salmon index -p $cores -t $transcripts -i $output_folder/$name/indexes/${organism}_salmon_idx --tmpdir $TMPDIR &> $output_folder/$name/indexes/${organism}_salmon_idx.log
				salmon_idx=$output_folder/$name/indexes/${organism}_salmon_idx
			fi
			kall_idx=$CURRENT_DIR/indexes/${organism}_kallisto_idx
			if [ ! -s "$kall_idx" ]; then
				kallisto index -i $output_folder/$name/indexes/${organism}_kallisto_idx $transcripts &> $output_folder/$name/indexes/${organism}_kallisto_idx.log
				kall_idx=$output_folder/$name/indexes/${organism}_kallisto_idx
			fi
			echo -e "\nPredicting strandness on two random sample...\n"		
			mkdir -p $output_folder/$name/strand_prediction/salmon_out; mkdir -p $output_folder/$name/strand_prediction/how_are_we_stranded_here_out	
			if [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "SINGLE" ]]; then
				salmon quant -i $salmon_idx -l A -r $seqs_location/$(ls $seqs_location | shuf | head -1) -p $cores -o $output_folder/$name/strand_prediction/salmon_out/ --skipQuant &> $output_folder/$name/strand_prediction/salmon_out/salmon_out.log
				cd $output_folder/$name/strand_prediction/how_are_we_stranded_here_out; check_strandedness --gtf $(echo $annotation | sed 'sa,.*aag') --transcripts $transcripts --reads_1 $seqs_location/$(ls $seqs_location | shuf | head -1) --kallisto_index $kall_idx --print_commands &> check_strandedness_out.log
			elif [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "PAIRED" ]]; then
				rand_sample_root=$(ls $seqs_location | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq | shuf | head -1)
				salmon quant -i $salmon_idx -l A -1 $seqs_location/${rand_sample_root}_1.fastq.gz -2 $seqs_location/${rand_sample_root}_2.fastq.gz -p $cores -o $output_folder/$name/strand_prediction/salmon_out/ --skipQuant &> $output_folder/$name/strand_prediction/salmon_out/salmon_out.log
				rand_sample_root=$(ls $seqs_location | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq | shuf | head -1)
				cd $output_folder/$name/strand_prediction/how_are_we_stranded_here_out; check_strandedness --gtf $(echo $annotation | sed 'sa,.*aag') --transcripts $transcripts --reads_1 $seqs_location/${rand_sample_root}_1.fastq.gz --reads_2 $seqs_location/${rand_sample_root}_2.fastq.gz --kallisto_index $kall_idx --print_commands &> check_strandedness_out.log
			fi
			cd $output_folder/$name/strand_prediction/salmon_out/
			salmon_strand=$(grep -r "Automatically detected most likely library type as " | sed 's,.*Automatically detected most likely library type as ,,g' | sort | uniq)
			if [[ "$salmon_strand" == "SR" || "$salmon_strand" == "ISR" ]]; then
				strand=reverse
			elif [[ "$salmon_strand" == "SF" || "$salmon_strand" == "ISF" ]]; then
				strand=yes
			elif [[ "$salmon_strand" == "U" || "$salmon_strand" == "IU" ]]; then
				strand=no
			fi
			strand_second_opinion=$(cat $output_folder/$name/strand_prediction/how_are_we_stranded_here_out/check_strandedness_out.log | grep "Data is likely " | sed 's,*Data is likely ,,g')
			cd $output_folder/$name	
			echo "Salmon prediction 1: $salmon_strand" > $output_folder/$name/strand_info.txt
			echo "Salmon prediction 2: $strand" >> $output_folder/$name/strand_info.txt
			echo -e "how_are_we_stranded_here prediction: $strand_second_opinion" >> $output_folder/$name/strand_info.txt		
			if [ $(grep -c "Data is likely" $output_folder/$name/strand_info.txt) -gt 0 ]; then		
	  			echo "Please double check carefully, based on the kit used in the library preparation, the paper, the GEO entry... because this is crucial for quantification. Please rerun with the argument '-s' in the unlikely case that the prediction by salmon is not correct, or if the second opinion by how_are_we_stranded_here is different (if transcripts from GENCODE or any particular format are used, the latter option may fail to identify the data type though)"
	    			cat $output_folder/$name/strand_info.txt
				rm $(find $output_folder/$name/strand_prediction/how_are_we_stranded_here_out -type f -name "*.bam")
			else
	  			echo "Salmon or how_are_we_stranded_here seem to have failed. This is not acceptable, plese double check. Exiting..."
				exit 1
			fi
		fi
	
	### Prepare other info required by the updated version of miARma...
		echo -e "\nPreparing miARma-seq execution...\n"
		number_files=$(ls $seqs_location | sed 's,_1.fastq.gz*,,g' | sed 's,_2.fastq.gz*,,g' | sort | uniq | wc -l)
		if [ $number_files -le $number_parallel ]; then
			cores_parallel=$((cores / number_files))
		else
			cores_parallel=$((cores / number_parallel))	
		fi
		if [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "SINGLE" ]]; then
			library_layout=Single
		elif [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "PAIRED" ]]; then
			library_layout=Paired
		fi
		read_length_for_miarma=$(zcat $seqs_location/$(ls $seqs_location | shuf | head -1) | head -2 | sed -n '2p' | awk '{print length -1}')
	
	### Prepare the ini file:
		IFS=', ' read -r -a array <<< "$annotation"
		IFS=', ' read -r -a array2 <<< "$optionsFeatureCounts_seq"
		IFS=', ' read -r -a array3 <<< "$optionsFeatureCounts_feat"
		for index in "${!array[@]}"; do
			cd $output_folder/$name
			cp $CURRENT_DIR/external_software/miARma-seq/bakk_miARma1.7.ini miarma$index.ini
			gff=${array[index]}
			sed -i "s,read_length=,read_length=$read_length_for_miarma,g" miarma$index.ini
			sed -i "s,read_dir=,read_dir=$seqs_location,g" miarma$index.ini
			sed -i "s,^threads=,threads=$cores_parallel,g" miarma$index.ini
			sed -i "s,label=,label=$name,g" miarma$index.ini
			sed -i "s,miARmaPath=,miARmaPath=$miarma_path,g" miarma$index.ini
			sed -i "s,output_dir=,output_dir=$output_folder/$name/miARma_out$index,g" miarma$index.ini
			sed -i "s,stats_file=miARma_stat.log,stats_file=$output_folder/$name/miARma_out$index/miARma_stat$index.log,g" miarma$index.ini
			sed -i "s,logfile=miARma_logfile.log,logfile=$output_folder/$name/miARma_out$index/miARma_logfile$index.log,g" miarma$index.ini
			sed -i "s,strand=yes,strand=$strand,g" miarma$index.ini
			sed -i "s,fasta=,fasta=$reference_genome,g" miarma$index.ini
			sed -i "s,gtf=,gtf=$gff,g" miarma$index.ini
			sed -i "s,database=,database=$gff,g" miarma$index.ini
			sed -i "s,seqtype=Paired,seqtype=$library_layout,g" miarma$index.ini
			sed -i "s,organism=mouse,organism=$organism,g" miarma$index.ini
			sed -i "s,indexthreads=,indexthreads=$indexthreads,g" miarma$index.ini
			sed -i "s,parallelnumber=,parallelnumber=$number_parallel,g" miarma$index.ini
			sed -i "s,memorylimit=,memorylimit=$memory_max,g" miarma$index.ini
			if [[ "$aligner" == "star" ]]; then
				if [ -z "$reference_genome_index" ]; then
					sed -i "s,indexname=,indexname=${organism}_$(basename ${reference_genome%.*})_$(basename ${gff%.*})_star_idx,g" miarma$index.ini
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" miarma$index.ini
				else
					sed -i "s,starindex=,starindex=$reference_genome_index,g" miarma$index.ini
					sed -i "s,indexname=,indexname=${organism}_star_idx,g" miarma$index.ini
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" miarma$index.ini
				fi
			elif [[ "$aligner" == "hisat2" ]]; then
				sed -i "s,aligner=star,aligner=hisat2,g" miarma$index.ini
				if [ -z "$reference_genome_index" ]; then
					sed -i "s,indexname=,indexname=${organism}_$(basename ${reference_genome%.*})_$(basename ${gff%.*})_hisat2_idx,g" miarma$index.ini
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" miarma$index.ini
				else
					sed -i "s,hisat2index=,hisat2index=$reference_genome_index,g" miarma$index.ini
					sed -i "s,indexname=,indexname=${organism}_hisat2_idx,g" miarma$index.ini
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" miarma$index.ini
				fi
			fi
			if [ ! -z "$optionsFeatureCounts_seq" ]; then
				sed -i "s,seqid=gene_name,seqid=${array2[index]},g" miarma$index.ini
			fi
			if [ ! -z "$optionsFeatureCounts_feat" ]; then
				sed -i "s,featuretype=exon,featuretype=${array3[index]},g" miarma$index.ini
			fi
			# Final renaming of fastq raw files if SRR present in the filename:
			if [ $(ls $seqs_location | grep -c SRR) -gt 0 ]; then
				for i in $(ls $seqs_location/*); do mv $i $(echo $i | sed 's,_SRR.*_,_,g'); done
			fi
		done
	fi
	echo -e "\nDONE. Current date/time: $(date)"; time1=`date +%s`; echo -e "Elapsed time (secs): $((time1-start))"; echo -e "Elapsed time (hours): $(echo "scale=2; $((time1-start))/3600" | bc -l)\n"
	export debug_step="all"
	echo -e "\n\nSTEP 3a: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP3b. Running miARma-seq:
# 2024: I've modified miARma RNA-seq mode to leverage GNU's parallel and increase speed, introduce limit RAM in aligners and multithreading index, replace the shebang with #!/usr/bin/env perl so it uses the PATH's/environment's one, etc...
# Eventually, WIP nicludes to also improve and integrate the rest of modules of miARma, such as adapter cutting, stats, miRNAs...
if [[ $debug_step == "all" || $debug_step == "step3b" ]]; then
	echo -e "\n\nSTEP 3b: Starting...\nCurrent date/time: $(date)\n\n"
	echo "Please double check all the parameters above, in particular the stranded or the reference genome files and annotation used. Proceeding with miARma execution in..."
	secs=$((1 * 30))
	dir=$output_folder/$name/miARma_out0
	while [ $secs -gt 0 ]; do
		echo -ne "$secs\033[0K\r"
		sleep 1
		: $((secs--))
	done
	for index in "${!array[@]}"; do
		if [ -d "$dir" ] && [ "$(ls -A $dir)" ] && [ "$index" -gt 0 ]; then
			dir2=$(echo $dir | sed "s,out0,out$index,g")
			mkdir -p $dir2; cd $dir2 
		fi
		cd $output_folder/$name
		if [ "$qc_raw_reads" == "no" ]; then
			mkdir -p $output_folder/$name/miARma_out$index/Pre_fastqc_results/_skip_
		fi
		$miarma_path/miARma miarma$index.ini
	done
		
	### Reformat the logs by parallel...
	for f in $(find $output_folder/$name -name "*_log_parallel.txt"); do awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' $f > tmp && mv tmp $f; done
	
 	echo -e "\nmiARma-seq and STEP 4 DONE. Current date/time: $(date)"; time1=`date +%s`; echo -e "Elapsed time (secs): $((time1-start))"; echo -e "Elapsed time (hours): $(echo "scale=2; $((time1-start))/3600" | bc -l)\n"
	export debug_step="all"
	echo -e "\n\nSTEP 3b: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP 4. Process output of miARma. Get figures, final counts, standard DGE, violin plots...
if [[ $debug_step == "all" || $debug_step == "step4" ]]; then
	echo -e "\n\nSTEP 4: Starting...\nCurrent date/time: $(date)\n\n"
	IFS=', ' read -r -a array2 <<< "$filter"
	for index in "${!array[@]}"; do
		annotation_file=${array[index]}
		R_process_reanalyzer_GSE.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index $genes ${array2[index]} $organism $target $differential_expr_soft $covariables $deconvolution $differential_expr_comparisons
		cd $output_folder/$name/final_results_reanalysis$index/DGE/
		tar -cf - $(ls | egrep ".RData$") | pigz -p $cores > allRData.tar.gz; rm -rf $(ls | egrep ".RData$")
	done
	export debug_step="all"
	echo -e "\n\nSTEP 4: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP 5. Time course analyses if required
if [[ $debug_step == "all" || $debug_step == "step5" ]]; then
	echo -e "\n\nSTEP 5: Starting...\nCurrent date/time: $(date)\n\n"
	for index in "${!array[@]}"; do
		if [[ "$time_course" == "yes" ]]; then 
			echo -e "\nPerforming time course analyses."
			R_process_time_course.R $output_folder/$name/final_results_reanalysis$index/ DGE_analysis_comp1.RData edgeR_object_norm $minstd $mestimate
		fi
	done
	export debug_step="all"
	echo -e "\n\nSTEP 5: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP 6. Functional enrichment analyses: clusterProfiler, autoGO, Panther, network analyses...
if [[ $debug_step == "all" || $debug_step == "step6" ]]; then
	echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
	for index in "${!array[@]}"; do
		if [[ "$organism" == "Mus_musculus" || "$organism" == "Homo_sapiens" || "$organism" == "Mus musculus" || "$organism" == "Homo sapiens" ]]; then
			if [ ! -z "$taxonid" ]; then			
				R_network_analyses.R $output_folder/$name/final_results_reanalysis$index/DGE/ $output_folder/$name/final_results_reanalysis$index/RPKM_counts_genes.txt "^DGE_analysis_comp[0-9]+.txt$" $taxonid &> network_analyses_funct_enrichment.log
			fi
			if [[ "$clusterProfiler" == "no" ]]; then
				echo -e "\nSkipping final clusterProfiler execution\n"
			else
				echo -e "\nPerforming clusterProfiler execution for DEGs. The results up to this point are ready to use (including DEGs and expression, even if not annotated), this and the may take long if many significant DEGs or comparisons, but check out the final steps of annotating and tyding and you may not need to wait...\n"
				cd $output_folder/$name/final_results_reanalysis$index/DGE/
				R_clusterProfiler_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index/DGE/ $organism $cores $clusterProfiler_method "^DGE_analysis_comp[0-9]+.txt$" &> clusterProfiler_funct_enrichment.log
				if [[ "$time_course" == "yes" ]]; then 
					cd $output_folder/$name/final_results_reanalysis$index/time_course_analyses
					R_clusterProfiler_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index/time_course_analyses $organism $cores $clusterProfiler_method "^DGE_limma_timecourse.*.txt$" &> clusterProfiler_funct_enrichment.log
				fi
				cd $output_folder/$name/final_results_reanalysis$index
				if [ $(find . -name clusterProfiler_funct_enrichment.log | xargs cat | grep -c "Ensembl site unresponsive") -gt 0 ] || [ $(find . -name clusterProfiler_funct_enrichment.log | xargs cat | grep -c "Error in curl") -gt 0 ] || [ $(find . -name funct_enrichment_analyses.tar.gz -exec tar -tzvf {} \; | grep -c _clusterProfiler/) -eq 0 ]; then
					echo "Apparently at least part of the parallel execution of clusterProfiler failed. Restarting and attempting serial execution. Please note it may be very slow..."
					cd $output_folder/$name/final_results_reanalysis$index/DGE; rm -rf $(ls | grep _clusterProfiler)
					R_clusterProfiler_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index/DGE/ $organism "1" $clusterProfiler_method "^DGE_analysis_comp[0-9]+.txt$" &> clusterProfiler_funct_enrichment_serial.log
					if [[ "$time_course" == "yes" ]]; then 
						cd $output_folder/$name/final_results_reanalysis$index/time_course_analyses; rm -rf $(ls | grep _clusterProfiler)
						R_clusterProfiler_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index/time_course_analyses $organism "1" $clusterProfiler_method "^DGE_limma_timecourse.*.txt$" &> clusterProfiler_funct_enrichment_serial.log
					fi
				fi
			fi
			echo -e "\nPerforming autoGO and Panther execution for the rest of relevant datasets. The rest of the results are ready, this may take long if many genes or comparisons...\n"		
			R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism $cores $databases_function "^DGE_analysis_comp.*\\.txt$|^DGE_limma_timecourse.*.txt$" $panther_method &> autoGO_panther_funct_enrichment.log
		else
			echo "Organism is $organism... Functional analyses apart from human/mouse is not fully supported yet"
			if [ $(egrep -c "GO:|Ontology|tology_term|tology term" $annotation_file) -gt 0 ]; then
				cd $output_folder/$name/final_results_reanalysis$index/DGE/
				echo "However, an automatic approach based on clusterProfiler's enrichr function and automatically extracted GO terms from the annotation can be attempted for DEGs..."
				paste <(egrep "GO:|,GO:|Ontology|tology_term|tology term" $annotation_file | sed 's,.*ID=,,g;s,.*Parent=,,g;s,;.*,,g') <(egrep "GO:|,GO:|Ontology|tology_term|tology term" $annotation_file | sed 's,.*tology_term=,,g') | sort -t $'\t' -k1,1 -k2,2 | awk -F'\t' '!a[$1,$2]++' | awk -F'\t' '{ a[$1] = (a[$1] ? a[$1]","$2 : $2); } END { for (i in a) print i"\t"a[i]; }' | awk -F '\t' '{n=split($2,a,","); for (i=1; i<=n; i++) print $1,a[i]}' | uniq > $output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annotation_file).automatically_extracted_GO_terms.txt
				annotation_go=$output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annotation_file).automatically_extracted_GO_terms.txt
				sed -i '1s/^/source_id Computed_GO_Process_IDs\n/' $annotation_go
				if [ -s "$annotation_go" ]; then
					R_clusterProfiler_enrichr.R $annotation_go $output_folder/$name/final_results_reanalysis$index/RPKM_counts_genes.txt $output_folder/$name/final_results_reanalysis$index/DGE "^DGE_analysis_comp[0-9]+.txt$" &> clusterProfiler_enrichr_funct_enrichment.log
				fi
			else
				echo "For $organism and the annotation $annotation_file, it does not seem there's GO or functional information available..."
			fi
		fi
		
		# Add to the tables of functional enrichment the number of genes up/down:
		cd $output_folder/$name/final_results_reanalysis$index/
		echo "Processing results of functional enrichment analyses..."
		find . \( -name "*.txt" -o -name "*.tsv" -o -name "*.csv" \) | grep funct | parallel -j $cores "file={}; R_enrich_format.R \"\$file\" \$(echo \"\$file\" | sed 's,DGE/.*,DGE/,g')\$(echo \"\$file\" | sed 's,.*DGE_analysis_comp,DGE_analysis_comp,g;s,_pval.*,,g;s,_fdr.*,,g;s,_funct.*,,g;s,_cluster.*,,g' | sed 's,.txt,,g').txt $organism $rev_thr" &> $PWD/enrichment_format.log
	done
	export debug_step="all"
	echo -e "\n\nSTEP 6: DONE\nCurrent date/time: $(date)\n\n"
fi

 
### STEP 7. Annotation: Tables of DEGs, lists of genes, etc
if [[ $debug_step == "all" || $debug_step == "step7" ]]; then
	echo -e "\n\nSTEP 7: Starting...\nCurrent date/time: $(date)\n\n"
	for index in "${!array[@]}"; do
		# All the tables that contain list of genes, annotate them:
		R_annotate_genes.R $output_folder/$name/final_results_reanalysis$index/ "^DGE_analysis_comp\\d+\\.txt$|^DGE_limma_timecourse_T\\d+_vs_T\\d+\\.txt$|mfuzz_elements_clusters|counts|WGCNA_all_modules_|STRINGdb_all_modules_" $organism
	
		# All the tables of DEGs, provide bed files for direct upload in genome browser
		cd $output_folder/$name/final_results_reanalysis$index/
		for file in $(find . -name "DGE_analysis_comp*" | egrep "_fdr_05.txt$|_pval_05.txt$"); do
			cut -f1 "$file" | parallel -j $((cores*2)) "gene={}; foldchange=\$(grep -i \"\$gene\" \"$file\" | cut -f3 | sed -n 's/\(.*[.,][0-9]\{2\}\).*/\1/p'); \
																   grep -i \"=\$gene\" \"$annotation_file\" | head -1 | awk -v id=\"\$gene\" -v fc=\"\$foldchange\" '{ print \$1\"\\t\"\$4\"\\t\"\$5\"\\t\"id\"_\"fc\"\\t.\t\"\$7 }' >> \"$file.bed\""
		done
	done
	export debug_step="all"
	echo -e "\n\nSTEP 7: DONE\nCurrent date/time: $(date)\n\n"
fi


###### STEP 8. Tidy up, prepare for storage if final results have been created and the number of aligned files is equal to the numbers of samples, rename folders, convert tables to xlsx if required... etc
# Compress the folders
if [[ $debug_step == "all" || $debug_step == "step8" ]]; then
	echo -e "\n\nSTEP 8: Starting...\nCurrent date/time: $(date)\n\n"
	for index in "${!array[@]}"; do
	 	cd $output_folder/$name/final_results_reanalysis$index/DGE/
		folders_funct=$(find . -type d \( -name "*_autoGO" -o -name "*_clusterProfiler" -o -name "*_panther" \))
		if [ -n "$folders_funct" ]; then
			tar -cf - $folders_funct | pigz --best -p $cores > funct_enrichment_analyses.tar.gz; rm -rf $folders_funct
		fi
		if [[ "$time_course" == "yes" ]]; then 
			cd $output_folder/$name/final_results_reanalysis$index/time_course_analyses
			folders_funct=$(find . -type d \( -name "*_autoGO" -o -name "*_clusterProfiler" -o -name "*_panther" \))
			if [ -n "$folders_funct" ]; then
				tar -cf - $folders_funct | pigz --best -p $cores > funct_enrichment_analyses.tar.gz; rm -rf $folders_funct
			fi
		fi
	done

	for f in $(find $output_folder -type d -name "final_results_reanaly*"); do
		mv $f $(echo $f"_"$(basename $output_folder))
	done
	
	for f in $(find $output_folder -name "*_QC*"); do
		mv $f $(echo $(dirname $f)"/"$(basename $output_folder)$(echo $f | sed 's,.*_QC,_QC,g'))
	done

	if [ "$tidy_tmp_files" == "yes" ]; then
		num_raw_files=$(cat $output_folder/$name/miARma_out0/Pre_fastqc_results/list_of_files.txt | grep -c "fastq.gz")
		for index in "${!array[@]}"; do
			if [ -d "$output_folder/$name/final_results_reanalysis$index" ] && [[ $(ls $output_folder/$name/miARma_out$index/$aligner\_results | egrep -c ".bam$") -eq $num_raw_files || $(ls $output_folder/$name/miARma_out$index/$aligner\_results | egrep -c ".bam$") -eq $((num_raw_files / 2)) ]]; then
				echo -e "\nTidying up...\n"
				cd $seqs_location
				echo "After execution, raw reads have been removed for the sake of efficient storage. These were... " > readme
				ls -lh >> readme
				rm $(ls | grep -v readme)
	
				cd $output_folder/$name/miARma_out$index/$aligner\_results
				echo "For the sake of efficiente storage: samtools view -@ cores -T ref_genome -C -o xxx.bam.cram xxx.bam && rm xx.bam" >> conversion_bam_to_cram.txt
				find . -type f -name "*.bam" | parallel --verbose -j $number_parallel --max-args 1 samtools view -T $reference_genome -C -@ $((cores / number_parallel)) -o {}.cram {}
				rm -rf $(ls | egrep ".bam$") $TMPDIR
			fi
		done
	fi

	if [ "$aligner" == "star" ] && [ "$aligner_index_cache" == "no" ]; then
		STAR --runThreadN $cores --genomeDir $(find $output_folder/$name/ -name "star_log_parallel.txt" | xargs cat | grep "genomeDir" | sed 's,.*genomeDir ,,g;s, .*,,g' | sort | uniq) --genomeLoad Remove --outFileNamePrefix genomeloading.tmp && rm genomeloading.tmp
	fi

	if [ "$convert_tables_excel" == "yes" ]; then
		R_convert_tables.R $output_folder/$name/ $cores "log_parallel|jquery|bamqc|rnaseqqc|samtools|strand" > /dev/null 2>&1
	fi
	export debug_step="all"
	echo -e "\n\nSTEP 8: DONE\nCurrent date/time: $(date)\n\n"
fi

echo -e "\n\n\nALL STEPS DONE! Best wishes\n\n\n"
