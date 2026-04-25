#!/bin/bash
start=`date +%s`
echo -e "\nCurrent time: $(date)\n"
base64 -d <<<"CiAgX19fICBfX19fICBfX19fICAgICAgX19fXyAgX19fXyAgIF9fICAgX18gXyAgIF9fICAgX18gICAgXyAgXyAgX19fXyAgX19fXyAgX19fXyAKIC8gX18pLyBfX18pKCAgX18pICAgICggIF8gXCggIF9fKSAvIF9cICggICggXCAvIF9cICggICkgICggXC8gKShfXyAgKSggIF9fKSggIF8gXAooIChfIFxcX19fIFwgKSBfKSAgICAgICkgICAvICkgXykgLyAgICBcLyAgICAvLyAgICBcLyAoXy9cICkgIC8gIC8gXy8gICkgXykgICkgICAvCiBcX19fLyhfX19fLyhfX19fKSAgICAoX19cXykoX19fXylcXy9cXy9cXylfXylcXy9cXy9cX19fXy8oX18vICAoX19fXykoX19fXykoX19cXykKCmJ5IEJpb2luZm9ybWF0aWNzIFVuaXQJCQkJSVBCTE4tQ1NJQy4gMjAyMwoKYmlvaW5mb3JtYXRpY2FAaXBiLmNzaWMuZXMJCSAgICAgICAgaHR0cHM6Ly9naXRodWIuY29tL0Jpb2luZm9JUEJMTi9yZWFuYWx5emVyR1NFCgo="
echo "doi.org/10.1101/2023.07.12.548663v2"

###### 0. Define arguments and variables:
### Export a string with command options and an array with arguments and deal with them in parse_options.sh
export options=$@
export arguments=($options)

CURRENT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $CURRENT_DIR/scripts/parse_options.sh
if [ $? -ne 0 ]; then
	echo "Exiting..."; exit 1
fi
mkdir -p $TMPDIR

# Initialize step timing log for Gantt chart (only on full runs, preserve on resume)
STEP_TIMES_FILE="$output_folder/$name/step_times.tsv"
if [[ $debug_step == "all" ]]; then
	mkdir -p "$output_folder/$name"
	echo -e "step\tepoch\tevent" > "$STEP_TIMES_FILE"
fi
_log_step() { mkdir -p "$(dirname "$STEP_TIMES_FILE")" 2>/dev/null; echo -e "$1\t$(date +%s)\t$2" >> "$STEP_TIMES_FILE" 2>/dev/null; }

###### STEP 1. Download info from GEO and organize metadata and so:
if [[ $debug_step == "all" || $debug_step == "step1" ]]; then
	rm -rf $output_folder/*
	if [[ $input == G* ]]; then
		echo -e "\n\nSTEP 1: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_1_Download" "start"
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
			echo "Please provide a comma-separated list with the conditions for each sample, and if more than one separate the comma-separated lists with '/', no spaces:"
			read -r design_input
			rm $(ls -d $output_folder/$name/GEO_info/* | grep 'design_possible_')
			IFS='/' read -ra ADDR <<< "$design_input"
			for i in "${!ADDR[@]}"; do
				# For each comma-separated list, split by ',' and echo to file
				IFS=',' read -ra ITEMS <<< "${ADDR[$i]}"
				for item in "${ITEMS[@]}"; do
					echo "$item"
				done > "$output_folder/$name/GEO_info/design_possible_full_$(($i + 1)).txt"
				cat "$output_folder/$name/GEO_info/design_possible_full_$(($i + 1)).txt" | sort | uniq > "$output_folder/$name/GEO_info/design_possible_$(($i + 1)).txt"
			done
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
_log_step "Step_1_Download" "end"
		echo -e "\nSTEP 1 DONE. Current time: $(date)\n"
	fi
_log_step "Step_1_Download" "end"
  	echo -e "\n\nSTEP 1: DONE\nCurrent date/time: $(date)\n\n"
	export debug_step="all"
fi


###### STEP 1. Download and process fastq files from the GEO ID provided:
if [[ $debug_step == "all" || $debug_step == "step1a" ]]; then
	rm -rf $seqs_location
	mkdir -p $TMPDIR
	if [[ $input == G* ]]; then
 		echo -e "\n\nSTEP 1: Downloading from the $input id provided...\nCurrent date/time: $(date)\n\n"
		if [ "$stop" == "yes" ]; then
			echo "You have requested a stop to manually provide the SRR ids, or potentially modify other files that may have not been detected properly from GEO, and were not correct, or you just want to adapt some of them. Please double check or manually modify the files GEO_info/srr_ids.txt, samples_info.txt, sample_names.txt, phenodata_extracted.txt, library_layout_info.txt, organism.txt, design_files, etc. The pipeline is stopped. Please press space to continue or Ctrl + C to exit..."
			read -n1 -s -r -p $'Press space to continue...\n' key
			rm $output_folder/$name/possible_designs_all.txt
			for i in $(ls $output_folder/$name/GEO_info | grep "full"); do echo $i >> $output_folder/$name/possible_designs_all.txt && cat $output_folder/$name/GEO_info/$i >> $output_folder/$name/possible_designs_all.txt && echo -e "\n" >> $output_folder/$name/possible_designs_all.txt; done
			cat $output_folder/$name/GEO_info/phenodata_extracted.txt > $output_folder/$name/phenotypic_data_samples.txt
		fi
		if [ ! -d "$seqs_location" ]; then # I'm now removing the seqs_location at the beginning of this section, in the context of the new system of resuming by -Dm stepx, so this should always be done
			mkdir -p $seqs_location
			echo "Downloading the fastq files from SRR..."
			if [ -z "$input_geo_reads" ]; then
				download_sra_fq.sh $output_folder/$name/GEO_info/srr_ids.txt $seqs_location $(( number_parallel*2 )) $cores
	### Rename the fastq files (max length name 140 characters) or handle already downloaded datasets if provided:
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
			num_files=$(ls | wc -l); num_samples=$(cat $output_folder/$name/GEO_info/srr_ids.txt | wc -l)
			if [ "$num_files" -lt "$num_samples" ]; then
				echo -e "\nPlease double check manually, is there some issue with the downloaded raw data? Exiting the script...\n"
				exit 1
			fi
		fi
		echo -e "\nDONE. Current date/time: $(date)"; time1=`date +%s`; echo -e "Elapsed time (secs): $((time1-start))"; echo -e "Elapsed time (hours): $(echo "scale=2; $((time1-start))/3600" | bc -l)\n"

	### Process if any download was not successful or subsampling was required:
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

		if [ ! -z "$number_reads" ]; then
			echo -e "\nSubsampling...\n"
			# From the input parameter by the user, obtain a random number allowing a +- 10% window:
			IFS=', ' read -r -a arr <<< "$number_reads"
			IFS=', ' read -r -a arr2 <<< "$(ls | egrep .fastq.gz$ | sed 's,[12].fastq.gz,,g' | sort | uniq | tr '\n' ',')"
			desired_number=${arr[1]}
			apply_random_shift() {
				Rscript -e '
				  modify_number <- function(number) {
				    percentage <- runif(1, 0, 10)  # Random % between 0 and 10
				    change <- ifelse(runif(1) < 0.5, -1, 1)  # Randomly add or subtract
				    change_amount <- number * (percentage / 100) * change
				    return(round(number + change_amount))
				  }
				  cat(modify_number('"$1"'), "\n")
				'
			}
   			# Name-based lookup: match each sample prefix against reads_numbers.txt
			declare -a arr3=()
			for sample_prefix in "${arr2[@]}"; do
				sample_name="${sample_prefix%_}"
				matched_line=$(awk -F'\t' -v name="$sample_name" 'NR>1 && $1 == name' "${arr[0]}")
				if [ -z "$matched_line" ]; then
					echo "WARNING: Sample '$sample_name' not found in ${arr[0]}. Keeping all reads."
					arr3+=("0")
				else
					col2=$(echo "$matched_line" | cut -f2)
					col3=$(echo "$matched_line" | cut -f3)
					desired_number_rand=$(apply_random_shift $desired_number)
					if (( col3 < desired_number_rand )); then
						arr3+=("$col2")
					else
						arr3+=("$((col2 * desired_number_rand / col3))")
					fi
				fi
			done
   			subsample_reads() {
				files=$(ls | grep $1)
				number=$2
				if [ "$number" -eq 0 ]; then
					echo "Skipping subsampling for $1 (not found in reads_numbers.txt, keeping all reads)"
					for file in $files; do cp "$file" "${file}_subsamp"; done
					return
				fi
				for file in $files; do seqtk sample -s 123 "$file" "$number" > "${file}_subsamp"; done
			}
			export -f subsample_reads
			parallel --verbose -j $cores_reads_to_subsample subsample_reads {} ::: "${arr2[@]}" :::+ "${arr3[@]}" # 10 only because of RAM
			rm $(ls | grep -v subsamp); for file in $(ls); do mv $file $(echo $file | sed 's,_subsamp,,g;s,.gz,,g'); done
			pigz --best -p $cores * # gz was lost with seqtk sample
			echo -e "\nSubsampling (+-10%) completed...\n"
		fi
_log_step "Step_1_Download" "end"
		echo -e "\n\nSTEP 1: DONE\nCurrent date/time: $(date)\n\n"
 	fi
	export debug_step="all"
fi


### STEP 1. Process if not required to download from NCBI/GEO the metadata and raw reads provided locally:

### STEP 1a_bis. Alignment Removal (Host Filtering)
if [[ $debug_step == "all" || $debug_step == "step1a_bis" ]]; then
    if [ ! -z "$alignment_removal" ]; then
        echo -e "\n\nSTEP 1a_bis: Alignment Removal (Host Filtering)...\nCurrent date/time: $(date)\n\n"
        
        # Define output directory for clean reads
        clean_seqs_location=$output_folder/$name/clean_reads_no_host
        mkdir -p $clean_seqs_location
        mkdir -p $output_folder/$name/indexes

        # Check/Build HISAT2 index for the removal genome
        removal_index_name=$(basename $alignment_removal)_hisat2_idx
        removal_index_path=$output_folder/$name/indexes/$removal_index_name
        
        if [ ! -f "${removal_index_path}.1.ht2" ]; then
            echo "Building HISAT2 index for alignment removal: $alignment_removal"
            hisat2-build -p $cores $alignment_removal $removal_index_path > $output_folder/$name/indexes/hisat2_build_removal.log 2>&1
        fi

        echo "Mapping reads against removal genome and extracting unmapped..."
        cd $seqs_location

        # Detect pairs or single
        if [ $(ls | egrep -c "_1.fastq.gz$") -gt 0 ]; then
             # Paired-end
             for f in $(ls | egrep "_1.fastq.gz$" | sed 's,_1.fastq.gz,,g'); do
                 echo "Processing $f..."
                 hisat2 -p $cores -x $removal_index_path \
                    -1 ${f}_1.fastq.gz -2 ${f}_2.fastq.gz \
                    --un-conc-gz $clean_seqs_location/${f}_%.fastq.gz \
                    --summary-file $clean_seqs_location/${f}_removal_summary.txt \
                    > /dev/null
 
                 # Rename output to match expected format (_1.fastq.gz instead of _1.fastq.gz) - hisat2 output fits usually but let's ensure
                 mv $clean_seqs_location/${f}_1.fastq.gz $clean_seqs_location/${f}_1.fastq.gz 2>/dev/null || true 
                 mv $clean_seqs_location/${f}_2.fastq.gz $clean_seqs_location/${f}_2.fastq.gz 2>/dev/null || true
             done
        else
             # Single-end
             for f in $(ls | egrep ".fastq.gz$"); do
                 echo "Processing $f..."
                 hisat2 -p $cores -x $removal_index_path \
                    -U $f \
                    --un-gz $clean_seqs_location/$f \
                    --summary-file $clean_seqs_location/${f}_removal_summary.txt \
                    > /dev/null
             done
        fi

        # Update seqs_location to point to clean reads
        echo -e "\nAlignment removal completed. Updating sequences location to: $clean_seqs_location\n"
        seqs_location=$clean_seqs_location
        export seqs_location
    fi
fi

if [[ $debug_step == "all" || $debug_step == "step1b" ]]; then
	mkdir -p $TMPDIR
	if [[ $input == /* ]]; then
		echo -e "\n\nSTEP 1b: Preparing the raw reads and metadata provided locally...\nCurrent date/time: $(date)\n\n"
		_log_step "Step_1b_Fastp" "start"
  		seqs_location=$output_folder/$name/raw_reads
		rm -rf $seqs_location # I'm now removing the seqs_location at the beginning of this section, in the context of the new system of resuming by -Dm stepx, so this should always be done
		if [ ! -d "$seqs_location" ]; then
			mkdir -p $seqs_location
			if [ $(ls -d $input/* | egrep -c "_R1.fastq.gz$|_R1.fq.gz$|_R2.fastq.gz$|_R2.fq.gz$|_1.fastq.gz$|_1.fq.gz$|_2.fastq.gz$|_2.fq.gz$") -eq 0 ]; then
				echo -e "\nPlease make sure that the input files are named _1.fastq.gz, _R1.fastq.gz, _2.fastq.gz, _R2.fastq.gz\n"
				exit 1
			fi
			for f in $(ls -d $input/*); do ln -sf $f $seqs_location/$(basename $f | sed 's,fq,fastq,g;s,_R1.fastq,_1.fastq,g;s,_R2.fastq,_2.fastq,g'); done
			if [ ! -z "$input_filter_regex" ]; then
				echo -e "\nFiltering input files with regex: $input_filter_regex\n"
				cd $seqs_location
				ls | egrep -v "$input_filter_regex" | xargs -r rm -f
				if [ $(ls | wc -l) -eq 0 ]; then
					echo "Error: No files matched the regex '$input_filter_regex'. Exiting..."
					exit 1
				fi
				cd - > /dev/null
			fi
			if [ ! -z "$input_filter_regex_exclude" ]; then
				echo -e "\nExcluding input files matching regex: $input_filter_regex_exclude\n"
				cd $seqs_location
				ls | egrep "$input_filter_regex_exclude" | xargs -r rm -f
				if [ $(ls | wc -l) -eq 0 ]; then
					echo "Error: All files were excluded by the regex '$input_filter_regex_exclude'. Exiting..."
					exit 1
				fi
				cd - > /dev/null
			fi
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
		if [ ! -z "$number_reads" ]; then
			echo -e "\nSubsampling...\n"
			# From the input parameter by the user, obtain a random number allowing a +- 10% window:
			IFS=', ' read -r -a arr <<< "$number_reads"
			IFS=', ' read -r -a arr2 <<< "$(ls | egrep .fastq.gz$ | sed 's,[12].fastq.gz,,g' | sort | uniq | tr '\n' ',')"			
			desired_number=${arr[1]}
			apply_random_shift() {
				Rscript -e '
				  modify_number <- function(number) {
				    percentage <- runif(1, 0, 10)  # Random % between 0 and 10
				    change <- ifelse(runif(1) < 0.5, -1, 1)  # Randomly add or subtract
				    change_amount <- number * (percentage / 100) * change
				    return(round(number + change_amount))
				  }
				  cat(modify_number('"$1"'), "\n")
				'
			}
   			# Name-based lookup: match each sample prefix against reads_numbers.txt
			declare -a arr3=()
			for sample_prefix in "${arr2[@]}"; do
				sample_name="${sample_prefix%_}"
				matched_line=$(awk -F'\t' -v name="$sample_name" 'NR>1 && $1 == name' "${arr[0]}")
				if [ -z "$matched_line" ]; then
					echo "WARNING: Sample '$sample_name' not found in ${arr[0]}. Keeping all reads."
					arr3+=("0")
				else
					col2=$(echo "$matched_line" | cut -f2)
					col3=$(echo "$matched_line" | cut -f3)
					desired_number_rand=$(apply_random_shift $desired_number)
					if (( col3 < desired_number_rand )); then
						arr3+=("$col2")
					else
						arr3+=("$((col2 * desired_number_rand / col3))")
					fi
				fi
			done
   			subsample_reads() {
				files=$(ls | grep $1)
				number=$2
				if [ "$number" -eq 0 ]; then
					echo "Skipping subsampling for $1 (not found in reads_numbers.txt, keeping all reads)"
					return
				fi
				for file in $files; do seqtk sample -s 123 "$file" "$number" > "${file}_subsamp"; done
			}
			export -f subsample_reads
			parallel --verbose -j $cores_reads_to_subsample subsample_reads {} ::: "${arr2[@]}" :::+ "${arr3[@]}" # 10 only because of RAM
			rm $(ls | grep -v subsamp); for file in $(ls); do mv $file $(echo $file | sed 's,_subsamp,,g;s,.gz,,g'); done
			pigz --best -p $cores * # gz was lost with seqtk sample
			echo -e "\nSubsampling (+-10%) completed...\n"
		fi
	 	echo -e "This is the content of $seqs_location:\n$(ls -l $seqs_location | awk '{ print $9 }' | tail -n +2)\n"
		if [ -z "$design_custom_local" ]; then
			echo -n "From the ordered list above, please input a comma-separated list with the conditions for each sample. Remember to try and avoid complex names, use as few underlines as possible, avoid names starting with numbers or others that would not be sorted appropriately such as containing spaces, if reads are paired-end, only one name of condition per pair of reads, and if you want to provide more than one design, separate the comma-separated list with a '/', no spaces: "
			read -r design_input
		else
			design_input=$design_custom_local
			echo -e "The used conditions are (the assignment must match sample names):\n"
			paste \
			  <(ls -l "$seqs_location" | awk '{ print $9 }' | tail -n +2 | sed 's,_[12].fastq.gz,,g' | uniq) \
			  <(echo "$design_input" | sed 's_,_\n_g;s,/,\n\n,g')
		fi
		mkdir -p $output_folder/$name/GEO_info/
		paste <(ls $seqs_location | egrep '.fq|.fastq' | sed 's,_[12].fastq.gz,,g' | uniq) <(paste -d'_' <(ls $seqs_location | egrep '.fq|.fastq' | sed 's,_[12].fastq.gz,,g' | uniq) <(echo $design_input | sed 's*/*\t*g'| cut -f1 | sed 's*,*\n*g')) <(echo $design_input | sed 's*/*\t*g'| cut -f1 | sed 's*,*\n*g') > $output_folder/$name/GEO_info/samples_info.txt
		IFS='/' read -ra ADDR <<< "$design_input"
		for i in "${!ADDR[@]}"; do
			# For each comma-separated list, split by ',' and echo to file
			IFS=',' read -ra ITEMS <<< "${ADDR[$i]}"
			for item in "${ITEMS[@]}"; do
				echo "$item"
			done > "$output_folder/$name/GEO_info/design_possible_full_$(($i + 1)).txt"
			cat "$output_folder/$name/GEO_info/design_possible_full_$(($i + 1)).txt" | sort | uniq > "$output_folder/$name/GEO_info/design_possible_$(($i + 1)).txt"
		done
		echo $name > $output_folder/$name/GEO_info/study_title.txt
		if [ -z "$organism_argument" ]; then
			echo -n "Please input the scientific name of the organism: "
			read -r organism
		else
			organism=$(echo $organism_argument | sed 's,_, ,g')
			echo "Organism used is $organism"
		fi
		echo $organism > $output_folder/$name/GEO_info/organism.txt
_log_step "Step_1b_Fastp" "end"
		echo -e "\n\nSTEP 1b: DONE\nCurrent date/time: $(date)\n\n"
 	fi
	export debug_step="all"
fi


### STEP 1. Process if required to download from manually provided ids from databases
if [[ $debug_step == "all" || $debug_step == "step1c" ]]; then
	if [[ $input == P* || $input == E* || $input == D* || $input == S* ]]; then
		echo -e "\n\nSTEP 1: Downloading from the $input id provided...\nCurrent date/time: $(date)\n\n"
  		seqs_location=$output_folder/$name/raw_reads
		rm -rf $seqs_location # I'm now removing the seqs_location at the beginning of this section, in the context of the new system of resuming by -Dm stepx, so this should always be done
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
_log_step "Step_1_Download" "end"
		echo -e "\n\nSTEP 1: DONE\nCurrent date/time: $(date)\n\n"
 	fi
	export debug_step="all"
fi


### STEP 1. Deal with batch correction... The user has to use certain arguments to manually provide a list or do it interactively:
if [[ $debug_step == "all" || $debug_step == "step1d" ]]; then
	if [ "$batch" == "yes" ]; then
		echo -e "\n\nSTEP 1: Preparing batch effect correction...\nCurrent date/time: $(date)\n\n"
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
	if [ ! -z "$covariables" ]; then
 		echo $covariables > $output_folder/$name/GEO_info/covariables.txt
  	fi
fi

### STEP 1. Give info of NCBI's current genome:
Rscript -e "organism <- '${organism}'; tryCatch({ ids <- rentrez::entrez_search(db='assembly', term=paste0(organism, '[orgn]'))\$ids; if(length(ids)==0) stop('No ids'); assemblies <- rentrez::entrez_summary(db='assembly', id=ids[1]); cat(paste(paste0('\n\nNCBI current assembly info: ', date()), assemblies\$assemblyname, assemblies\$assemblyaccession, assemblies\$submissiondate, '\n', sep='\n')) }, error=function(e) cat(paste0('\n\nNo genome information found in NCBI for: ', organism, '\n\n')))"
organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's/ \+/_/g;s/__*/_/g') # Get again the organism in case it has been manually modified... and without spaces...

### STEP 1. Auto-decompress gzipped reference inputs into the indexes subfolder:
mkdir -p $output_folder/$name/indexes
declare -a files_to_decompress=()
declare -a decompressed_outputs=()

if [[ "$reference_genome" == *.gz ]]; then
	decompressed_genome="$output_folder/$name/indexes/$(basename ${reference_genome%.gz})"
	if [ ! -f "$decompressed_genome" ]; then
		files_to_decompress+=("$reference_genome")
		decompressed_outputs+=("$decompressed_genome")
	else
		echo -e "\nDecompressed reference genome already exists: $decompressed_genome\n"
	fi
	reference_genome="$decompressed_genome"
fi

IFS=', ' read -r -a _annot_array <<< "$annotation"
for _ai in "${!_annot_array[@]}"; do
	if [[ "${_annot_array[$_ai]}" == *.gz ]]; then
		_decompressed_annot="$output_folder/$name/indexes/$(basename ${_annot_array[$_ai]%.gz})"
		if [ ! -f "$_decompressed_annot" ]; then
			files_to_decompress+=("${_annot_array[$_ai]}")
			decompressed_outputs+=("$_decompressed_annot")
		else
			echo -e "\nDecompressed annotation already exists: $_decompressed_annot\n"
		fi
		_annot_array[$_ai]="$_decompressed_annot"
	fi
done
annotation=$(IFS=','; echo "${_annot_array[*]}")

if [[ "$transcripts" == *.gz ]]; then
	decompressed_transcripts="$output_folder/$name/indexes/$(basename ${transcripts%.gz})"
	if [ ! -f "$decompressed_transcripts" ]; then
		files_to_decompress+=("$transcripts")
		decompressed_outputs+=("$decompressed_transcripts")
	else
		echo -e "\nDecompressed transcripts already exists: $decompressed_transcripts\n"
	fi
	transcripts="$decompressed_transcripts"
fi

if [ ${#files_to_decompress[@]} -gt 0 ]; then
	parallel --tmpdir $TMPDIR -j $number_parallel 'echo "Using {1} -> {2}"; gunzip -c "{1}" > "{2}"' ::: "${files_to_decompress[@]}" :::+ "${decompressed_outputs[@]}"
fi

### STEP 1. Deal with fastp if required:
mkdir -p $TMPDIR
cores_fastp=$((cores / number_parallel))
if [ $cores_fastp -lt 1 ]; then cores_fastp=1; fi

# Build fastp-specific options based on config
fastp_extra_opts=""
fastp_labels=()

if [ "$fastp_mode" == "yes" ]; then
	fastp_labels+=("quality filtering")
elif [ "$fastp_mode" == "no" ] && [ "$fastp_adapter" != "no" ]; then
	fastp_extra_opts="--disable_quality_filtering"
fi

if [ "$fastp_adapter" == "yes" ]; then
	fastp_extra_opts="$fastp_extra_opts --detect_adapter_for_pe"
	fastp_labels+=("adapter trimming")
elif [[ $fastp_adapter == /* ]]; then
	fastp_extra_opts="$fastp_extra_opts --adapter_fasta $fastp_adapter"
	fastp_labels+=("custom adapter trimming")
elif [ "$fastp_adapter" == "no" ] && [ "$fastp_mode" == "yes" ]; then
	fastp_extra_opts="$fastp_extra_opts --disable_adapter_trimming"
fi

fastp_label=""
if [ ${#fastp_labels[@]} -gt 0 ]; then
	joined_labels=$(IFS=','; echo "${fastp_labels[*]}" | sed 's/,/ and /g')
	fastp_label="Preprocessing with fastp ($joined_labels)"
fi

# Append any user-provided extra fastp arguments
if [ ! -z "$fastp_extra_args" ]; then
	fastp_extra_opts="$fastp_extra_opts $fastp_extra_args"
fi

if [ ! -z "$fastp_label" ]; then
	echo -e "\n\nSTEP 1: $fastp_label...(output files will be renamed and moved to raw_reads internal folder)\nCurrent date/time: $(date)\n\n"
	mkdir -p $output_folder/$name/fastp_out
	cd $output_folder/$name/fastp_out
	layout_fastp=$(find $output_folder/$name -name library_layout_info.txt | xargs cat)

	if [[ "$layout_fastp" == "SINGLE" ]]; then
		ls -d $seqs_location/*.fastq.gz | \
			parallel --tmpdir $TMPDIR --verbose --joblog $output_folder/$name/fastp_out/fastp_log_parallel.txt -j $number_parallel \
			'fastp --in1 {} --out1 {}_fastp.fastq.gz --dont_overwrite --dont_eval_duplication '$fastp_extra_opts' --thread '$cores_fastp' -z '$compression_level' -h '$output_folder/$name'/fastp_out/{/}_report.html -j '$output_folder/$name'/fastp_out/{/}_report.json &>> '$output_folder/$name'/fastp_out/{/}_fastp_out.log'
	elif [[ "$layout_fastp" == "PAIRED" ]]; then
		ls -d $seqs_location/*.fastq.gz | sed 's,_[12].fastq.gz,,g' | sort | uniq | \
			parallel --tmpdir $TMPDIR --verbose --joblog $output_folder/$name/fastp_out/fastp_log_parallel.txt -j $number_parallel \
			'fastp --in1 {}_1.fastq.gz --in2 {}_2.fastq.gz --out1 {}_1.fastq.gz_fastp.fastq.gz --out2 {}_2.fastq.gz_fastp.fastq.gz --dont_overwrite --dont_eval_duplication '$fastp_extra_opts' --thread '$cores_fastp' -z '$compression_level' -h '$output_folder/$name'/fastp_out/{/}_report.html -j '$output_folder/$name'/fastp_out/{/}_report.json &>> '$output_folder/$name'/fastp_out/{/}_fastp_out.log'
	fi
fi

# Trimming step can run in addition to the adapter/mode step above
if [ "$fastp_trimming" != "none" ]; then
	echo -e "\n\nSTEP 1: Preprocessing with fastp and trimming...\nCurrent date/time: $(date)\n\n"
	mkdir -p $output_folder/$name/fastp_out
	cd $output_folder/$name/fastp_out
	IFS=', ' read -r -a arrfastp <<< "$fastp_trimming"
	layout_fastp=$(find $output_folder/$name -name library_layout_info.txt | xargs cat)

	if [[ "$layout_fastp" == "SINGLE" ]]; then
		ls -d $seqs_location/*.fastq.gz | \
			parallel --tmpdir $TMPDIR --verbose --joblog $output_folder/$name/fastp_out/fastp_trim_log_parallel.txt -j $number_parallel \
			'fastp --in1 {} --out1 {}_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --trim_front1 '"${arrfastp[0]}"' --trim_tail1 '"${arrfastp[1]}"' --thread '$cores_fastp' -z '$compression_level' -h '$output_folder/$name'/fastp_out/{/}_trim_report.html -j '$output_folder/$name'/fastp_out/{/}_trim_report.json &>> '$output_folder/$name'/fastp_out/{/}_fastp_trim_out.log'
	elif [[ "$layout_fastp" == "PAIRED" ]]; then
		ls -d $seqs_location/*.fastq.gz | sed 's,_[12].fastq.gz,,g' | sort | uniq | \
			parallel --tmpdir $TMPDIR --verbose --joblog $output_folder/$name/fastp_out/fastp_trim_log_parallel.txt -j $number_parallel \
			'fastp --in1 {}_1.fastq.gz --in2 {}_2.fastq.gz --out1 {}_1.fastq.gz_fastp.fastq.gz --out2 {}_2.fastq.gz_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --trim_front1 '"${arrfastp[0]}"' --trim_tail1 '"${arrfastp[1]}"' --thread '$cores_fastp' -z '$compression_level' -h '$output_folder/$name'/fastp_out/{/}_trim_report.html -j '$output_folder/$name'/fastp_out/{/}_trim_report.json &>> '$output_folder/$name'/fastp_out/{/}_fastp_trim_out.log'
	fi
fi

if [ $(ls -d $seqs_location/* | grep -c _fastp.fastq.gz) -gt 0 ]; then
	for f in $(ls -d $seqs_location/* | grep _fastp.fastq.gz); do mv $f $(echo $f | sed 's,.fastq.gz_fastp.fastq.gz,.fastq.gz,g'); done
	echo "Files in $seqs_location have been successfully processed by fastp!"
fi


### STEP 2. Decontamination if required:
if [[ $debug_step == "all" || $debug_step == "step2" ]]; then
	if [ ! -z "$kraken2_databases" ]; then
  		echo -e "\n\nSTEP 2: Decontamination starting with Kraken2 (k2 daemon mode)...\nCurrent date/time: $(date)\n\n"
_log_step "Step_2_Decontamination" "start"
    		rm -rf $output_folder/$name/raw_reads_k2
		mkdir -p $output_folder/$name/raw_reads_k2
		cd $output_folder/$name/raw_reads_k2

		IFS=',' read -r -a k2_db_array <<< "$kraken2_databases"
		IFS=',' read -r -a k2_conf_array <<< "$kraken2_confidence"
		layout=$(find $output_folder/$name -name library_layout_info.txt | xargs cat)

		for k2_db in "${k2_db_array[@]}"; do
			db_basename=$(basename $k2_db)
			echo -e "\n\nProcessing database: $db_basename ($k2_db)\n"

			for conf in "${k2_conf_array[@]}"; do
				# Build confidence label: e.g. 0 -> "00", 0.20 -> "020", 0.50 -> "050"
				conf_label=$(echo "$conf" | sed 's/^0$/00/;s/^0\.\([0-9]*\)/0\1/')

				echo -e "\n  Running k2 classify with confidence=$conf (label=$conf_label) on database $db_basename...\n"

				if [[ "$layout" == "SINGLE" ]]; then
					for f in $(ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$"); do
						sample_base=$(basename $f | sed 's,\.\(fastq\|fq\)\(\.gz\)\?$,,')
						k2 classify --use-daemon --db $k2_db --threads $cores \
							$([ "$conf" != "0" ] && echo "--confidence $conf") \
							--report-minimizer-data --use-names \
							--output $PWD/${sample_base}.${conf_label}k2_${db_basename} \
							--report $PWD/${sample_base}.${conf_label}k2_${db_basename}_report.txt \
							--log $PWD/${sample_base}.${conf_label}k2_${db_basename}_log.txt \
							$seqs_location/$f \
							1>> kraken2_log_out.txt 2>> kraken2_log_warnings_errors.txt
					done
				elif [[ "$layout" == "PAIRED" ]]; then
					for f in $(ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$" | sed 's,_R\?[12]\.\(fastq\|fq\)\(\.gz\)\?$,,g' | sort | uniq); do
						k2 classify --use-daemon --db $k2_db --threads $cores \
							$([ "$conf" != "0" ] && echo "--confidence $conf") \
							--report-minimizer-data --use-names --paired \
							--output $PWD/${f}.${conf_label}k2_${db_basename} \
							--report $PWD/${f}.${conf_label}k2_${db_basename}_report.txt \
							--log $PWD/${f}.${conf_label}k2_${db_basename}_log.txt \
							${seqs_location}/${f}_*.fastq.gz \
							1>> kraken2_log_out.txt 2>> kraken2_log_warnings_errors.txt
					done
				fi
			done

			# Compress output files for this DB and stop daemon
			pigz --best -p $cores *k2_${db_basename} 2>/dev/null
			k2 clean --stop-daemon 2>/dev/null

			# Process reports for this DB
			for f in $(ls | egrep "k2_${db_basename}_report.txt$" 2>/dev/null); do
				echo -e "\nLog of kraken2 ($f):"
				echo -e "%_reads_covered\t#_reads_covered\t#_reads_directly_assigned\t#_minimizers_total\t#_minimizers_distinct\tRank_code\tTaxon_id\tScientific_name" >> $f.final.txt
				cat $f >> $f.final.txt
				echo -e "\n\nNumber of classified reads at the genus level: $(cat $f | awk '$6 == "G" {print $2"\t"$7}' | awk '{s+=$1}END{print s}')" >> $f.final.txt
				echo -e "\nTaxonomy IDs at the genus level assigned to the reads:" >> $f.final.txt
				echo -e "#read\tTaxID\n$(awk '$6 == "G" {print $2"\t"$7}' $f)\n" >> $f.final.txt
			done
		done

		echo -e "Please check the files *_report.txt, kraken2_log_out.txt and kraken2_log_warnings_errors.txt"
		echo -e "Processing recentrifuge reports and extracting reads..."

		# Recentrifuge for all output files
		first_db=$(basename ${k2_db_array[0]})
		taxdump_dir=""
		if [[ -d "${k2_db_array[0]}/taxdump" ]]; then
			taxdump_dir="${k2_db_array[0]}/taxdump"
		fi
		if [ ! -z "$taxdump_dir" ]; then
			for f in $(ls | egrep "k2_.*\.gz$" 2>/dev/null); do rcf -n $taxdump_dir -k $f -o $f.recentrifuge_contamination_report.html -e CSV &>> rcf_log_out.txt; done
		fi

		# Extract reads for the organism (use first DB for taxdump if available)
		if [ ! -z "$taxdump_dir" ]; then
			if [ -z "$taxonid" ]; then
				taxonid=$(echo $organism | sed 's/_\+/ /g' | taxonkit name2taxid --data-dir $taxdump_dir | head -1 | cut -f2)
			fi
			taxon_name=$(taxonkit list --ids $taxonid -n -r --data-dir $taxdump_dir | grep $taxonid)
			echo -e "\nOrganism provided: $organism"; echo -e "\nOrganism provided (taxonid): $taxonid"; echo $taxon_name
			echo -e "\nIf not correct, please rerun and double check that you have provided it explicitly in the prompt... kraken2 output will be filtered to retain that taxa and below"
			echo -e "\nCheck out the logs in the files rcf_log_out.txt and extract_kraken2_log_out.txt"

			mkdir -p $seqs_location\_k2
			# Use first confidence score + first DB for read extraction
			conf_label_first=$(echo "${k2_conf_array[0]}" | sed 's/^0$/00/;s/^0\.\([0-9]*\)/0\1/')
			if [[ "$layout" == "SINGLE" ]]; then
				for f in $(ls | egrep "\.${conf_label_first}k2_${first_db}\.gz$"); do
					base=$(echo $f | sed "s,\.${conf_label_first}k2_${first_db}\.gz,,g")
					extract_kraken_reads.py -k $f -U $base -o $seqs_location\_k2/${f}_1.fastq.gz -t $taxonid -r ${base}.${conf_label_first}k2_${first_db}_report.txt --include-children &>> extract_kraken2_log_out.txt
				done
			elif [[ "$layout" == "PAIRED" ]]; then
				for f in $(ls | egrep "\.${conf_label_first}k2_${first_db}\.gz$"); do
					base=$(echo $f | sed "s,\.${conf_label_first}k2_${first_db}\.gz,,g")
					extract_kraken_reads.py -k $f -s1 ${base}_1.fastq.gz -s2 ${base}_2.fastq.gz -o $seqs_location\_k2/${f}_1.fastq.gz -o2 $seqs_location\_k2/${f}_2.fastq.gz -t $taxonid -r ${base}.${conf_label_first}k2_${first_db}_report.txt --include-children &>> extract_kraken2_log_out.txt
				done
			fi
		fi

		for f in $(ls | grep "k2" | egrep ".fastq.gz$"); do fastqc -q -t $cores $f; done
_log_step "Step_2_Decontamination" "end"
  		echo -e "\n\nSTEP 2: DONE\nCurrent date/time: $(date)\n\n"
	fi
	if [ ! -z "$sortmerna_databases" ]; then
		echo -e "\n\nSTEP 2: Decontamination starting with sortmerna...\nCurrent date/time: $(date)\n\n"
_log_step "Step_2_Decontamination" "start"
		mkdir -p $seqs_location\_sortmerna $output_folder/$name/indexes/$(basename $sortmerna_databases)_sortmerna_index
		cd $seqs_location\_sortmerna

		### Build index before parallel execution (uses all cores, done once)
		if [ ! -d "$output_folder/$name/indexes/$(basename $sortmerna_databases)_sortmerna_index/idx" ]; then
			echo "Indexing the provided $sortmerna_databases ..."
			sortmerna --index 1 --ref $sortmerna_databases \
				--workdir $output_folder/$name/indexes/$(basename $sortmerna_databases)_sortmerna_index \
				--threads $cores &>> $output_folder/$name/indexes/sortmerna_index.log
		fi

		rm -rf $seqs_location\_sortmerna/*

		### Compute per-job thread count (mirror fastp logic)
		cores_sortmerna=$((cores / number_parallel))
		if [ $cores_sortmerna -lt 1 ]; then cores_sortmerna=1; fi

		layout_sortmerna=$(find $output_folder/$name -name library_layout_info.txt | xargs cat)
		sortmerna_idx=$output_folder/$name/indexes/$(basename $sortmerna_databases)_sortmerna_index/idx
		sortmerna_out=$seqs_location\_sortmerna

		echo -e "\nExecuting sortmerna in parallel and fastqc of the new reads...\n"

		if [[ "$layout_sortmerna" == "PAIRED" ]]; then
			ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$" \
				| sed 's,_[12].fastq.gz,,g' | sort | uniq \
				| parallel --tmpdir $TMPDIR --verbose \
						   --joblog $sortmerna_out/sortmerna_log_parallel.txt \
						   -j $number_parallel --max-args 1 \
				"sortmerna \
					--idx-dir $sortmerna_idx \
					--ref $sortmerna_databases \
					--reads $seqs_location/{}_1.fastq.gz \
					--reads $seqs_location/{}_2.fastq.gz \
					--workdir $sortmerna_out/{}_sortmerna_workdir \
					--fastx --threads $cores_sortmerna --out2 --index 0 \
					--aligned $sortmerna_out/{}_rRNA \
					--other $sortmerna_out/{}_no_rRNA \
					-v &>> $sortmerna_out/{}_out.log"

		elif [[ "$layout_sortmerna" == "SINGLE" ]]; then
			ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$" \
				| parallel --tmpdir $TMPDIR --verbose \
						   --joblog $sortmerna_out/sortmerna_log_parallel.txt \
						   -j $number_parallel --max-args 1 \
				"sortmerna \
					--idx-dir $sortmerna_idx \
					--ref $sortmerna_databases \
					--reads $seqs_location/{} \
					--workdir $sortmerna_out/{}_sortmerna_workdir \
					--fastx --threads $cores_sortmerna --index 0 \
					--aligned $sortmerna_out/{}_rRNA \
					--other $sortmerna_out/{}_no_rRNA \
					-v &>> $sortmerna_out/{}_out.log"
		fi

		### Clean up per-sample workdirs
		rm -rf $sortmerna_out/*_sortmerna_workdir

		fastqc -q -t $cores $sortmerna_out/*.fq.gz
		mkdir $sortmerna_out/out_noRNA; cd $sortmerna_out/out_noRNA
		ln -sf ../*no_rRNA*.fq.gz .
		for f in $(ls); do mv $f $(basename $f | sed 's,.fq.gz,.fastq.gz,g;s,_fwd,_1,g;s,_rev,_2,g;s,_no_rRNA,,g'); done
		export seqs_location=$sortmerna_out/out_noRNA
_log_step "Step_2_Decontamination" "end"
		echo -e "\n\nSTEP 2: DONE\nCurrent date/time: $(date)\n\n"
	fi
	export debug_step="all"
fi


### STEP 2b. Preliminary rRNA QC (Bowtie2 mapping against rRNA references):
if [[ $debug_step == "all" || $debug_step == "step2b" ]]; then
	if [ ! -z "$rrna_qc_databases" ]; then
		echo -e "\n\nSTEP 2b: rRNA QC ...\nCurrent date/time: $(date)\n\n"
_log_step "Step_2b_rRNA_QC" "start"

		RRNA_QC_DIR=$output_folder/$name/preliminar_rrna_qc
		RRNA_INDEX_DIR=$RRNA_QC_DIR
		RRNA_INDEX_NAME="rrna_ref"
		RRNA_MIN_SCORE=$rrna_qc_min_score

		mkdir -p "$RRNA_QC_DIR"

		# Handle comma-separated reference paths
		IFS=',' read -ra RRNA_REF_ARRAY <<< "$rrna_qc_databases"

		if [ ! -f "$RRNA_INDEX_DIR/${RRNA_INDEX_NAME}.1.bt2" ]; then
			echo "Creating combined FASTA from: $rrna_qc_databases"
			RRNA_COMBINED_FASTA="$RRNA_INDEX_DIR/combined_ref.fasta"
			> "$RRNA_COMBINED_FASTA"

			for ref_file in "${RRNA_REF_ARRAY[@]}"; do
				if [[ "$ref_file" == *.gz ]]; then
					zcat "$ref_file" >> "$RRNA_COMBINED_FASTA"
				else
					cat "$ref_file" >> "$RRNA_COMBINED_FASTA"
				fi
			done

			# Check if RNA (U) or DNA (T) and convert if needed
			FIRST_SEQ_LINE=$(grep -v "^>" "$RRNA_COMBINED_FASTA" | head -n 1)

			if echo "${FIRST_SEQ_LINE}" | grep -q -i "U"; then
				RRNA_CONVERTED_FILE="$RRNA_INDEX_DIR/combined_ref_DNA_converted.fasta.gz"
				echo "RNA detected (Uracil 'U' found). Converting RNA to DNA..."
				awk '/^>/ {print; next} {gsub(/U/,"T"); gsub(/u/,"t"); print}' "$RRNA_COMBINED_FASTA" | gzip > "$RRNA_CONVERTED_FILE"
				rm "$RRNA_COMBINED_FASTA"
				echo "Converted file saved as: ${RRNA_CONVERTED_FILE}"
			elif echo "${FIRST_SEQ_LINE}" | grep -q -i "T"; then
				echo "DNA detected (Thymine 'T' found). No conversion needed."
				RRNA_CONVERTED_FILE="$RRNA_COMBINED_FASTA"
			else
				echo "Warning: Could not detect RNA (U) or DNA (T). Assuming no conversion needed."
				RRNA_CONVERTED_FILE="$RRNA_COMBINED_FASTA"
			fi

			bowtie2-build --threads "$cores" "$RRNA_CONVERTED_FILE" "$RRNA_INDEX_DIR/$RRNA_INDEX_NAME" > "$RRNA_INDEX_DIR/${RRNA_INDEX_NAME}_build.log" 2>&1
		else
			echo "Index $RRNA_INDEX_DIR/$RRNA_INDEX_NAME already exists. Skipping build."
		fi

		# Step 2b: Align R1 reads and count
		RRNA_SUMMARY_FILE="$RRNA_QC_DIR/rRNA_mapping_summary_R1.tsv"
		echo -e "Sample\tTotal_Reads\tSense_Count\tAntisense_Count\tSense_Pct\tAntisense_Pct" > "$RRNA_SUMMARY_FILE"

		for fq in "$seqs_location"/*_1.fastq.gz "$seqs_location"/*_R1*.fastq.gz; do
			# Skip if no fastq files exist
			[ -e "$fq" ] || continue
			# Skip duplicates
			[ -f "$fq" ] || continue

			# Extract sample name
			nm=$(basename "$fq" | sed -E 's/(_R1(_[0-9]+)?|_1)\.fastq\.gz//')

			# Get total read count
			TOTAL_READS=$(zcat "$fq" | wc -l | awk '{print $1/4}')

			TMP_COUNTS="$RRNA_QC_DIR/${nm}_tmp_counts.txt"
			ERR_LOG="$RRNA_QC_DIR/${nm}_bowtie2.err"

			# Run bowtie2
			bowtie2 -x "$RRNA_INDEX_DIR/$RRNA_INDEX_NAME" -U "$fq" \
			  -k 10 --trim5 4 --trim3 4 --very-sensitive \
			  -p "$cores" --no-unal --no-hd --mm 2> "$ERR_LOG" | \
			grep "AS:i:" | \
			awk -v OFS='\t' '{
			    score="";
			    for(i=12; i<=NF; i++){
			        if($i ~ /^AS:i:/){
			            score=substr($i, 6);
			            break;
			        }
			    }
			    if(score != "") print $1, $2, $3, score
			}' > "$TMP_COUNTS"

			# Process counts (filtering for best score and min score)
			counts=$(awk -v min_score="$RRNA_MIN_SCORE" -F'\t' '
			{
			    read_id=$1; flag=$2; score=$4;
			    if((score+0) >= (min_score+0)) {
			        if (!(read_id in best_score) || (score+0) > (best_score[read_id]+0)) {
			            best_score[read_id] = score
			            best_flag[read_id] = flag
			        }
			    }
			}
			END {
			    sense_count=0;
			    antisense_count=0;
			    for (id in best_score) {
			        if (best_flag[id] == 0) {
			            sense_count++
			        } else if (best_flag[id] == 16) {
			            antisense_count++
			        }
			    }
			    print sense_count " " antisense_count
			}' "$TMP_COUNTS")

			SENSE_COUNT=$(echo "$counts" | cut -d' ' -f1)
			ANTISENSE_COUNT=$(echo "$counts" | cut -d' ' -f2)

			if [ "$TOTAL_READS" -gt 0 ]; then
			    SENSE_PCT=$(awk -v s="$SENSE_COUNT" -v t="$TOTAL_READS" 'BEGIN { printf "%.2f", (s/t)*100 }')
			    ANTISENSE_PCT=$(awk -v a="$ANTISENSE_COUNT" -v t="$TOTAL_READS" 'BEGIN { printf "%.2f", (a/t)*100 }')
			else
			    SENSE_PCT=0.00
			    ANTISENSE_PCT=0.00
			fi

			echo -e "${nm}\t${TOTAL_READS}\t${SENSE_COUNT}\t${ANTISENSE_COUNT}\t${SENSE_PCT}\t${ANTISENSE_PCT}" >> "$RRNA_SUMMARY_FILE"

			rm -f "$TMP_COUNTS"
		done

		echo -e "\nResults summarized in: $RRNA_SUMMARY_FILE"
		Rscript $CURRENT_DIR/scripts/R_rrna_qc_plot.R "$RRNA_QC_DIR" 2>&1 | tee -a "$RRNA_QC_DIR/R_rrna_qc_plot.log"

_log_step "Step_2b_rRNA_QC" "end"
		echo -e "\n\nSTEP 2b: DONE\nCurrent date/time: $(date)\n\n"
	fi
	export debug_step="all"
fi


### STEP3a. Prepare the data and info for running miARma-seq:
if [[ $debug_step == "all" || $debug_step == "step3a" ]]; then
	echo -e "\n\nSTEP 3a: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_3a_Prepare" "start"
 	echo -e "Preparing miARma-seq execution"
	cd $output_folder/$name/
	rm -rf mi*
	mkdir -p $TMPDIR
	# If the running is resumed in this step, the above has to be done
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi

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

			mkdir -p $output_folder/$name/strand_prediction/salmon_out
			if [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "SINGLE" ]]; then
				rand_sample=$(ls $seqs_location | shuf | head -1)
				echo -e "\nPredicting strandness with salmon on random sample: $rand_sample\n"
				salmon quant -i $salmon_idx -l A -r $seqs_location/$rand_sample -p $cores -o $output_folder/$name/strand_prediction/salmon_out/ --skipQuant &> $output_folder/$name/strand_prediction/salmon_out/salmon_out.log
			elif [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "PAIRED" ]]; then
				rand_sample_root=$(ls $seqs_location | sed 's,_[12].fastq.gz,,g' | uniq | shuf | head -1)
				rand_sample="${rand_sample_root}_1.fastq.gz / ${rand_sample_root}_2.fastq.gz"
				echo -e "\nPredicting strandness with salmon on random sample: $rand_sample\n"
				salmon quant -i $salmon_idx -l A -1 $seqs_location/${rand_sample_root}_1.fastq.gz -2 $seqs_location/${rand_sample_root}_2.fastq.gz -p $cores -o $output_folder/$name/strand_prediction/salmon_out/ --skipQuant &> $output_folder/$name/strand_prediction/salmon_out/salmon_out.log
			fi
			salmon_meta=$output_folder/$name/strand_prediction/salmon_out/aux_info/meta_info.json
			salmon_strand=$(python3 -c "import json; d=json.load(open('$salmon_meta')); print(d['library_types'][0])" 2>/dev/null)
			salmon_json=$output_folder/$name/strand_prediction/salmon_out/lib_format_counts.json
			if [ -f "$salmon_json" ]; then
				salmon_compat_ratio=$(python3 -c "import json; d=json.load(open('$salmon_json')); print(d.get('compatible_fragment_ratio', 'N/A'))" 2>/dev/null)
			else
				salmon_compat_ratio="N/A"
			fi
			if [[ "$salmon_strand" == "SR" || "$salmon_strand" == "ISR" ]]; then
				strand="reverse"
			elif [[ "$salmon_strand" == "SF" || "$salmon_strand" == "ISF" ]]; then
				strand="yes"
			elif [[ "$salmon_strand" == "U" || "$salmon_strand" == "IU" ]]; then
				strand="no"
			fi
			cd $output_folder/$name
			echo "Salmon library type: $salmon_strand" > $output_folder/$name/strand_info.txt
			echo "Strandedness: $strand" >> $output_folder/$name/strand_info.txt
			echo "Compatible fragment ratio: $salmon_compat_ratio" >> $output_folder/$name/strand_info.txt
			echo "Sample used: $rand_sample" >> $output_folder/$name/strand_info.txt
			if [ $(egrep -c "reverse|yes|no" $output_folder/$name/strand_info.txt) -gt 0 ]; then
				echo "Please double check carefully, based on the kit used in the library preparation, the paper, the GEO entry... because this is crucial for quantification. Please rerun with the argument '-s' in the unlikely case that the prediction by salmon is not correct"
				cat $output_folder/$name/strand_info.txt
			else
				echo -e "Salmon to detect strandedness seems to have failed. This is not acceptable, please double check or provide the parameter -s. Exiting...\nThe random sample used was: $rand_sample"
				exit 1
			fi
		else
  			echo $strand > $output_folder/$name/strand_info.txt
  		fi

	### Prepare other info required by the updated version of miARma...
		echo -e "\nPreparing miARma-seq execution...\n"
		number_files=$(ls $seqs_location | sed 's,_[12].fastq.gz.*,,g' | uniq | wc -l)
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
			elif [[ "$aligner" == "kallisto" ]]; then
				sed -i "s,aligner=star,aligner=kallisto,g" miarma$index.ini
				sed -i "s,fasta=$reference_genome,fasta=$transcripts,g" miarma$index.ini
				if [ -z "$reference_genome_index" ]; then
					sed -i "s,indexname=,indexname=${organism}_$(basename ${transcripts%.*})_kallisto_idx,g" miarma$index.ini
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" miarma$index.ini
				else
					sed -i "s,kallistoindex=,kallistoindex=$reference_genome_index,g" miarma$index.ini
					sed -i "s,indexname=,indexname=${organism}_kallisto_idx,g" miarma$index.ini
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" miarma$index.ini
				fi
			fi
			if [ ! -z "$optionsFeatureCounts_seq" ]; then
				sed -i "s,seqid=gene_name,seqid=${array2[index]},g" miarma$index.ini
			fi
			if [ ! -z "$optionsFeatureCounts_feat" ]; then
				sed -i "s,featuretype=exon,featuretype=${array3[index]},g" miarma$index.ini
			fi

			# ── Validate featureCounts parameters against annotation file ──
			mkdir -p $TMPDIR
			fc_feat_val="${array3[index]:-exon}"
			fc_seq_val="${array2[index]:-gene_name}"
			if [ -f "$gff" ]; then
				# Check feature type (-t) exists in column 3
				available_feats=$(zcat -f "$gff" | awk -F'\t' '!/^#/ && NF>=9 {print $3}' | sort -u | tr '\n' ', ' | sed 's/,$//')
				if ! zcat -f "$gff" | awk -F'\t' -v ft="$fc_feat_val" '!/^#/ && NF>=9 && $3==ft {found=1; exit} END {exit !found}'; then
					echo -e "\n\033[1;31mERROR:\033[0m Feature type '$fc_feat_val' (optionsFeatureCounts_feat / -t) was NOT found in column 3 of annotation file:\n  $gff\n\nAvailable feature types: $available_feats\n\nPlease set 'optionsFeatureCounts_feat' in your YAML config to one of the above (e.g. 'exon' for GTF, 'gene' for some GFF3 files).\n" >&2
					exit 1
				fi
				# Check attribute name (-g) exists as a proper key in column 9 (not substring) for the given feature type
				if ! zcat -f "$gff" | awk -F'\t' -v ft="$fc_feat_val" -v attr="$fc_seq_val" '!/^#/ && NF>=9 && $3==ft { n=split($9,pairs,";"); for(i=1;i<=n;i++){ gsub(/^[ \t]+/,"",pairs[i]); split(pairs[i],kv,/[ =]+/); if(kv[1]==attr){found=1; exit} } } END {exit !found}'; then
					# Extract example attributes from the first data line with the selected feature type
					example_attrs=$(zcat -f "$gff" | awk -F'\t' -v ft="$fc_feat_val" '!/^#/ && NF>=9 && $3==ft {print $9; exit}')
					echo -e "\n\033[1;31mERROR:\033[0m Attribute name '$fc_seq_val' (optionsFeatureCounts_seq / -g) was NOT found in column 9 of annotation file (for feature type '$fc_feat_val'):\n  $gff\n\nExample attributes from your file (for '$fc_feat_val'):\n  $example_attrs\n\nFor GTF files, typical values are 'gene_id' or 'gene_name'.\nFor GFF3 files, typical values are 'ID', 'Name', or 'gene_id' (or 'Parent' if feature type is 'exon').\nPlease set 'optionsFeatureCounts_seq' and 'optionsFeatureCounts_feat' in your YAML config accordingly.\n" >&2
					exit 1
				fi
				echo "Annotation validation OK: feature type '$fc_feat_val' and attribute '$fc_seq_val' found in $gff"
			fi
			# ── End featureCounts parameter validation ──
			if [ "$bam_mapq_threshold" -gt 0 ] 2>/dev/null; then
				sed -i "s,quality=10,quality=$bam_mapq_threshold,g" miarma$index.ini
				sed -i "s,bam_mapq_threshold=,bam_mapq_threshold=$bam_mapq_threshold,g" miarma$index.ini
			fi
			if [ ! -z "$bam_require_flags" ]; then
				sed -i "s,bam_require_flags=,bam_require_flags=$bam_require_flags,g" miarma$index.ini
			fi
			if [ ! -z "$bam_exclude_flags" ]; then
				sed -i "s,bam_exclude_flags=,bam_exclude_flags=$bam_exclude_flags,g" miarma$index.ini
			fi
			if [ ! -z "$bam_dedup" ]; then
				sed -i "s,bam_dedup=no,bam_dedup=$bam_dedup,g" miarma$index.ini
			fi
			if [ ! -z "$bam_custom_filter" ]; then
				bam_custom_filter_escaped=$(printf '%s' "$bam_custom_filter" | sed 's/[\\&]/\\&/g')
				sed -i "s,bam_custom_filter=,bam_custom_filter=$bam_custom_filter_escaped,g" miarma$index.ini
			fi
			if [ ! -z "$bam_normalization" ]; then
				sed -i "s,bam_normalization=,bam_normalization=$bam_normalization,g" miarma$index.ini
			fi
			if [ ! -z "$featureCounts_extra_args" ]; then
				fc_extra_escaped=$(printf '%s' "$featureCounts_extra_args" | sed 's/[\\&]/\\&/g')
				sed -i "s,parameters=-M -O -C -B,parameters=$fc_extra_escaped,g" miarma$index.ini
			fi
			if [ ! -z "$aligner_extra_args" ]; then
				ae_escaped=$(printf '%s' "$aligner_extra_args" | sed 's/[\\&]/\\&/g')
				if [[ "$aligner" == "star" ]]; then
					sed -i "s,starparameters=,starparameters=$ae_escaped,g" miarma$index.ini
				elif [[ "$aligner" == "hisat2" ]]; then
					sed -i "s,hisat2parameters=,hisat2parameters=$ae_escaped,g" miarma$index.ini
				elif [[ "$aligner" == "kallisto" ]]; then
					sed -i "s,kallistoparameters=,kallistoparameters=$ae_escaped,g" miarma$index.ini
				fi
				echo "Aligner extra args for $aligner: $aligner_extra_args"
			fi
			# Final renaming of fastq raw files if SRR present in the filename:
			if [ $(ls $seqs_location | grep -c SRR) -gt 0 ]; then
				for i in $(ls $seqs_location/*); do mv $i $(echo $i | sed 's,_SRR.*_,_,g'); done
			fi
		done
	fi
	export debug_step="all"
_log_step "Step_3a_Prepare" "end"
	echo -e "\n\nSTEP 3a: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP3b. Running miARma-seq:
# miARma RNA-seq mode was modified to leverage GNU's parallel and increase speed, introduce limit RAM in aligners and multithreading index, replace the shebang with #!/usr/bin/env perl so it uses the PATH's/environment's one, etc...
# Eventually, WIP nicludes to also improve and integrate the rest of modules of miARma, such as adapter cutting, stats, miRNAs...
if [[ $debug_step == "all" || $debug_step == "step3b" ]]; then
	echo -e "\n\nSTEP 3b: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_3b_miARma" "start"
	rm -rf $output_folder/$name/miARma_out*
	mkdir -p $TMPDIR
	# If the running is resumed in this step, the above has to be done
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi

	echo -e "miARma configuration .ini:"
        cat miarma$index.ini
	echo -e "\nPlease double check all the parameters above for miARma, in particular the stranded or the reference genome files and annotation used. Proceeding with miARma execution in..."
	secs=$((1 * 15))
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

 	### Clean genome index cache?
  	if [ "$aligner" == "star" ] && [ "$aligner_index_cache" == "no" ]; then
		STAR --runThreadN $cores --genomeDir $(find $output_folder/$name/ -name "star_log_parallel.txt" | xargs cat | grep "genomeDir" | sed 's,.*genomeDir ,,g;s, .*,,g' | sort | uniq) --genomeLoad Remove --outFileNamePrefix genomeloading.tmp
	fi

	echo -e "\nmiARma-seq DONE. Current date/time: $(date)"; time1=`date +%s`; echo -e "Elapsed time (secs): $((time1-start))"; echo -e "Elapsed time (hours): $(echo "scale=2; $((time1-start))/3600" | bc -l)\n"
	export debug_step="all"
_log_step "Step_3b_miARma" "end"
	echo -e "\n\nSTEP 3b: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP 4. Process output of miARma. Get figures, final counts, standard DGE, violin plots...
if [[ $debug_step == "all" || $debug_step == "step4" ]]; then
	# If the running is resumed in this step, this variables has to be created because they would not exist
	rm -rf $output_folder/$name/final_results_* # So it's redone when resuming
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi
	if [ -z "${!array[@]}" ]; then
		IFS=', ' read -r -a array <<< "$annotation"
	fi
	echo -e "\n\nSTEP 4: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_4_R_Process" "start"
 	echo -e "Processing output of miARma-seq, QC figures, plots, DGE if requested..."
	IFS=', ' read -r -a array2 <<< "$filter"
	for index in "${!array[@]}"; do
		annotation_file=${array[index]}
		fc_seq_key=${optionsFeatureCounts_seq:-gene_name}
		fc_feat_type=${optionsFeatureCounts_feat:-exon}
		# Apply counts_custom_gene_filter if provided (filter gene rows from count tables before R processing)
		if [ ! -z "$counts_custom_gene_filter" ]; then
			echo -e "\nApplying counts_custom_gene_filter: $counts_custom_gene_filter\n"
			for count_file in $(find $output_folder/$name/miARma_out$index -name "*_readcount.tab" -o -name "abundance.tsv" 2>/dev/null); do
				cp "$count_file" "${count_file}.bak_before_gene_filter"
				head -1 "$count_file" > "${count_file}.tmp"
				tail -n +2 "$count_file" | eval "$counts_custom_gene_filter" >> "${count_file}.tmp"
				mv "$count_file.tmp" "$count_file"
				echo "  Filtered: $count_file ($(wc -l < "${count_file}.bak_before_gene_filter") -> $(wc -l < "$count_file") lines)"
			done
		fi
  			echo -e "R_process_reanalyzer_GSE.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index $genes ${array2[index]} $organism $target $differential_expr_soft $batch_format $covariables $covariables_format $deconvolution $differential_expr_comparisons $perform_differential_analyses $perform_volcano_venn $pattern_to_remove $annotation_file $fc_seq_key $fc_feat_type $sc_count_matrix $sc_phenotype $bulk_expression_matrix\n\n" > $output_folder/$name/R_process_reanalyzer.log
    		R_process_reanalyzer_GSE.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index $genes ${array2[index]} $organism $target $differential_expr_soft $batch_format $covariables $covariables_format $deconvolution $differential_expr_comparisons $perform_differential_analyses $perform_volcano_venn $pattern_to_remove $annotation_file $fc_seq_key $fc_feat_type $sc_count_matrix $sc_phenotype $bulk_expression_matrix | tee -a $output_folder/$name/R_process_reanalyzer.log
    		echo 'R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter" "edgeR_object" "edgeR_object_norm" $pattern_to_remove $annotation_file $fc_feat_type' > $output_folder/$name/R_qc_figs.log
			R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter" "edgeR_object" "edgeR_object_norm" $pattern_to_remove $annotation_file $fc_feat_type | tee -a $output_folder/$name/R_qc_figs.log
		if [[ -e "$output_folder/$name/final_results_reanalysis$index/counts_adjusted.txt" ]]; then
			echo -e "\n\nRemember that batch effect correction/covariables have been only provided to Combat-Seq/limma for visualization purposes, to include covariables in the DGE model after checking the visualization the argument -C will be used\n\n\nQC_PDF adjusted counts\n\nRemember that you have requested batch effect correction/count adjustment, so you have to mind the figures in this QC_PDF from ComBat-seq/limma counts...\n"
			echo -e 'R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter_adjusted" "edgeR_object_adjusted" "edgeR_object_norm_adjusted" $pattern_to_remove $annotation_file $fc_feat_type' > $output_folder/$name/R_qc_figs_adjusted.log
			R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter_adjusted" "edgeR_object_adjusted" "edgeR_object_norm_adjusted" $pattern_to_remove $annotation_file $fc_feat_type | tee -a $output_folder/$name/R_qc_figs_adjusted.log
		fi
		cd $output_folder/$name/final_results_reanalysis$index/DGE/
		tar -cf - $(ls | egrep ".RData$") | pigz -p $cores > allRData.tar.gz; rm -rf $(ls | egrep ".RData$")
	done
	### Generate SummarizedExperiment for exploreDE app if requested
	if [[ "$exploreDE_se" == "yes" ]]; then
		echo -e "\nGenerating SummarizedExperiment for exploreDE...\n"
		for index in "${!array[@]}"; do
			final_dir=$output_folder/$name/final_results_reanalysis$index
			if [ -f "$final_dir/Raw_counts_genes.txt" ] && [ -f "$final_dir/TPM_counts_genes.txt" ] && [ -f "$final_dir/DGE/list_comp.txt" ]; then
				export ANNOTATION_FILE="${array[index]}"
				Rscript $CURRENT_DIR/scripts/prepare_SE.R \
					"$final_dir/Raw_counts_genes.txt" \
					"$final_dir/TPM_counts_genes.txt" \
					"$output_folder/$name/GEO_info/samples_info.txt" \
					"$final_dir/DGE/list_comp.txt" \
					"$final_dir/DGE" \
					"^DGE_analysis_comp[0-9].txt$" \
					"$name" \
					"$organism" 2>&1 | tee -a "$final_dir/prepare_SE.log"
			else
				echo "Skipping exploreDE SE generation for index $index: required files not found"
			fi
		done
	fi

	export debug_step="all"
_log_step "Step_4_R_Process" "end"
	echo -e "\n\nSTEP 4: DONE\nCurrent date/time: $(date)\n\n"
	if [[ "$perform_differential_analyses" == "no" ]]; then
		echo "Differential analyses not requested, exiting the pipeline..."; exit 1
	fi
fi


### STEP 4b. Splicing analysis if requested (saseR or IsoformSwitchAnalyzeR)
if [[ $debug_step == "all" || $debug_step == "step4b" ]]; then
	_log_step "Step_4b_Splicing" "start"
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi
	if [ -z "${!array[@]}" ]; then
		IFS=', ' read -r -a array <<< "$annotation"
	fi
	if [[ "$splicing_option" == "saser" ]]; then
		echo -e "\n\nSTEP 4b: saseR splicing analysis...\nCurrent date/time: $(date)\n\n"
		for index in "${!array[@]}"; do
			bam_dir=$output_folder/$name/miARma_out$index/${aligner}_results
			saser_out=$output_folder/$name/final_results_reanalysis$index/saseR_splicing
			mkdir -p $saser_out
			library_layout=$(find $output_folder/$name -name library_layout_info.txt | xargs cat)
			samples_info=$output_folder/$name/GEO_info/samples_info.txt
			design_file=$(ls $output_folder/$name/GEO_info/design_possible_full_1.txt 2>/dev/null || echo "none")
			R_saseR_splicing.R \
				"$bam_dir" \
				"${array[index]}" \
				"$saser_out" \
				"$library_layout" \
				"$strand" \
				"$samples_info" \
				"$design_file" \
				"$cores" \
				"$pattern_to_remove" \
				"$differential_expr_comparisons" \
				2>&1 | tee -a $saser_out/saseR_splicing.log
		done
_log_step "Step_4b_Splicing" "end"
		echo -e "\n\nSTEP 4b (saseR): DONE\nCurrent date/time: $(date)\n\n"
	elif [[ "$splicing_option" == "isoformswitchr" ]]; then
		if [[ "$aligner" != "kallisto" ]]; then
			echo -e "\n\033[1;31mERROR:\033[0m IsoformSwitchAnalyzeR (splicing_option=isoformswitchr) requires transcript-level quantification and can only be used with aligner='kallisto'.\n  Current aligner: '$aligner'\n  Please re-run with -A kallisto (or aligner: \"kallisto\" in your YAML config).\n" >&2
			exit 1
		fi
		echo -e "\n\nSTEP 4b: IsoformSwitchAnalyzeR analysis...\nCurrent date/time: $(date)\n\n"
		if [ -z "$transcripts" ]; then
			echo -e "WARNING: No transcript FASTA (-t) provided. Some IsoformSwitchAnalyzeR features may be limited.\n"
			transcripts_arg="none"
		else
			transcripts_arg=$transcripts
		fi
		for index in "${!array[@]}"; do
			# Determine quantification directory: Kallisto or Salmon results
			if [[ "$aligner" == "kallisto" ]]; then
				quant_dir=$output_folder/$name/miARma_out$index/kallisto_results
			else
				# Check if Salmon quantification exists (from strandness prediction or separate run)
				quant_dir=$output_folder/$name/miARma_out$index/salmon_results
				if [ ! -d "$quant_dir" ] || [ $(find $quant_dir -name "quant.sf" 2>/dev/null | wc -l) -eq 0 ]; then
					quant_dir=$output_folder/$name/strand_prediction/salmon_out
				fi
			fi
			isoswitch_out=$output_folder/$name/final_results_reanalysis$index/IsoformSwitchAnalyzeR
			mkdir -p $isoswitch_out
			samples_info=$output_folder/$name/GEO_info/samples_info.txt
			design_file=$(ls $output_folder/$name/GEO_info/design_possible_full_1.txt 2>/dev/null || echo "none")
			R_isoformswitch.R \
				"$quant_dir" \
				"${array[index]}" \
				"$transcripts_arg" \
				"$isoswitch_out" \
				"$samples_info" \
				"$design_file" \
				"$cores" \
				"$pattern_to_remove" \
				"$differential_expr_comparisons" \
				"$aligner" \
				2>&1 | tee -a $isoswitch_out/isoformswitch.log
		done
_log_step "Step_4b_Splicing" "end"
		echo -e "\n\nSTEP 4b (IsoformSwitchAnalyzeR): DONE\nCurrent date/time: $(date)\n\n"
	fi
	export debug_step="all"
fi


### STEP 5. Time course analyses if required
if [[ $debug_step == "all" || $debug_step == "step5" ]]; then
	for index in "${!array[@]}"; do
		if [[ "$time_course" == "yes" ]]; then
			echo -e "\n\nSTEP 5: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_5_QC_Figs" "start"
   			echo -e "\nPerforming time course analyses."
			R_process_time_course.R $output_folder/$name/final_results_reanalysis$index/ DGE_analysis_comp1.RData edgeR_object_norm $minstd $mestimate
_log_step "Step_5_QC_Figs" "end"
   			echo -e "\n\nSTEP 5: DONE\nCurrent date/time: $(date)\n\n"
		fi
	done
	export debug_step="all"	
fi


### STEP 6. Functional enrichment analyses: clusterProfiler, autoGO, Panther, network analyses...
if [[ $debug_step == "all" || $debug_step == "step6" ]]; then
	# If the running is resumed in this step, this variables has to be created because they would not exist
 	# The same may happen in other steps, this is WIP
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi
	if [[ $network_analyses == "yes" ]]; then
		if [ -z "$taxonid" ]; then
			cd $TMPDIR
			mkdir -p taxdump && cd taxdump && rm -rf * && wget -q https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip && unzip -qq taxdmp.zip && rm taxdmp.zip
			taxonid=$(echo $organism | sed 's/_\+/ /g' | taxonkit name2taxid --data-dir $PWD | head -1 | cut -f2)
		fi
	fi
	if [ -z "${!array[@]}" ]; then
		IFS=', ' read -r -a array <<< "$annotation"
	fi

	for index in "${!array[@]}"; do
		cd $output_folder/$name/final_results_reanalysis$index/DGE/
		rm -rf $(find . -type d \( -name "*_autoGO" -o -name "*_clusterProfiler" -o -name "*_panther" -o -name "*funct_enr*" \)) $(find . -type f \( -name "*_autoGO" -o -name "*_clusterProfiler" -o -name "*_panther" -o -name "*funct_enr*" \)) # So it's redone if resuming
		if [ -z "$annotation_file" ]; then
			annotation_file=${array[index]}
		fi

		# Network analyses (Only implemented for Human/Mouse)
		if [[ "$organism" == "Mus_musculus" || "$organism" == "Homo_sapiens" || "$organism" == "Mus musculus" || "$organism" == "Homo sapiens" ]]; then
			if [[ $network_analyses == "yes" ]]; then
				mkdir -p network_analyses && rm -rf network_analyses/* && cd network_analyses
				echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_6_Enrichment" "start"
    				echo -e "\nPerforming network analyses...\n"
				R_network_analyses.R $output_folder/$name/final_results_reanalysis$index/DGE/ $output_folder/$name/final_results_reanalysis$index/RM_counts_genes.txt "^DGE_analysis_comp[0-9]+.txt$" $taxonid &> network_analyses_funct_enrichment.log
			fi
		fi

		# Functional Enrichment Analyses
		if [[ "$functional_enrichment_analyses" == "no" ]]; then
			echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_6_Enrichment" "start"
    			echo -e "\nSkipping functional enrichment analyses\n"
		else
			if [[ "$organism" == "Mus_musculus" || "$organism" == "Homo_sapiens" || "$organism" == "Mus musculus" || "$organism" == "Homo sapiens" ]]; then
				echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_6_Enrichment" "start"
    				echo -e "\nPerforming functional enrichment analyses for DEGs. The results up to this point are ready to use (including DEGs and expression table including gene_ids). This step of funtional enrichment analyses may take long if many significant DEGs, comparisons, or analyses...\n"
				export ANNOTATION_FILE="${array[index]}"
				cd $output_folder/$name/final_results_reanalysis$index/DGE/
				ls | egrep "^DGE_analysis_comp[0-9]+.txt$" | parallel --joblog R_clusterProfiler_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_clusterProfiler_analyses_parallel.R $PWD $organism "1" $clusterProfiler_method $clusterProfiler_full $aPEAR_execution '^{}$' $clusterProfiler_universe $clusterProfiler_minGSSize $clusterProfiler_maxGSSize &> clusterProfiler_{}_funct_enrichment.log"
				echo -e "\nPerforming autoGO and Panther execution... this may take long if many genes or comparisons...\n"
				ls | egrep "^DGE_analysis_comp[0-9]+.txt$" | parallel --joblog R_autoGO_panther_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism "1" $databases_function {} $panther_method $auto_panther_log &> autoGO_panther_{}_funct_enrichment.log"
				if [[ "$time_course" == "yes" ]]; then
					cd $output_folder/$name/final_results_reanalysis$index/time_course_analyses
					ls | egrep "^DGE_limma_timecourse.*.txt$" | parallel --joblog R_clusterProfiler_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_clusterProfiler_analyses_parallel.R $PWD $organism "1" $clusterProfiler_method $clusterProfiler_full $aPEAR_execution '^{}$' $clusterProfiler_universe $clusterProfiler_minGSSize $clusterProfiler_maxGSSize &> clusterProfiler_{}_funct_enrichment.log"
					ls | egrep "^DGE_limma_timecourse.*.txt$" | parallel --joblog R_autoGO_panther_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism "1" $databases_function {} $panther_method $auto_panther_log &> autoGO_panther_{}_funct_enrichment.log"
				fi
			else
				echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_6_Enrichment" "start"
   				echo "Organism is $organism... Functional analyses apart from human/mouse is not fully supported yet"
				# Determine which annotation file to use for functional enrichment
				annot_enrichm=""
				if [ ! -z "$non_reference_funct_enrichm" ]; then
					echo "Using provided non-reference functional enrichment file: $non_reference_funct_enrichm"
					annot_enrichm="$non_reference_funct_enrichm"
				elif [ $(egrep -c "GO:|Ontology|tology_term|tology term" $annotation_file) -gt 0 ]; then
					echo "Using main annotation file for functional enrichment: $annotation_file"
					annot_enrichm="$annotation_file"
				fi

				enrichment_results_found="no"
				if [ ! -z "$annot_enrichm" ]; then
					cd $output_folder/$name/final_results_reanalysis$index/DGE/
					echo "Applying clusterProfiler enrichr with GO terms extracted from the provided annotation..."

					# Check if input is GAF (Gene Association File)
					if [[ "$annot_enrichm" == *.gaf ]] || [[ "$annot_enrichm" == *.gaf.gz ]]; then
						# Detect number of columns to handle standard (17-col) vs simplified (2-col) GAF
						gaf_ncols=$(zcat -f "$annot_enrichm" | grep -v "^!" | head -1 | awk -F'\t' '{print NF}')
						if [ "$gaf_ncols" -le 2 ]; then
							echo "Detected simplified ${gaf_ncols}-column GAF file. Extracting GeneID and GO term columns..."
							# Determine column order: one col has GO:xxxx pattern, the other is the gene ID
							first_col=$(zcat -f "$annot_enrichm" | grep -v "^!" | head -1 | cut -f1)
							if [[ "$first_col" == GO:* ]]; then
								# Format: GO_ID<tab>GeneID -> swap to GeneID<tab>GO_ID
								zcat -f "$annot_enrichm" | grep -v "^!" | awk -F'\t' '{print $2"\t"$1}' | sort -u > $output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annot_enrichm).automatically_extracted_GO_terms.txt
							else
								# Format: GeneID<tab>GO_ID (already correct)
								zcat -f "$annot_enrichm" | grep -v "^!" | cut -f 1,2 | sort -u > $output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annot_enrichm).automatically_extracted_GO_terms.txt
							fi
						else
							echo "Detected standard ${gaf_ncols}-column GAF format. Extracting Gene IDs (col2) and GO terms (col5)..."
							# Standard GAF 2.x: Column 2 = DB Object ID (Gene ID), Column 5 = GO ID
							zcat -f "$annot_enrichm" | grep -v "^!" | cut -f 2,5 | sort -u > $output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annot_enrichm).automatically_extracted_GO_terms.txt
						fi

					elif [[ "$annot_enrichm" == *.gmt ]] || [[ "$annot_enrichm" == *.gmt.gz ]]; then
						echo "Detected GMT format. Transforming to GeneID-TermID format..."
						# GMT format: TermID <tab> Description <tab> Gene1 <tab> Gene2 ...
						# We need: GeneID <tab> TermID
						# Process GMT:
						# 1. Read line
						# 2. Extract Term (col 1)
						# 3. Iterate from col 3 to end (Gene IDs)
						# 4. Print "GeneID \t TermID"
						zcat -f "$annot_enrichm" | awk -F'\t' '{term=$1; for(i=3;i<=NF;i++) print $i"\t"term}' | sort -u > $output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annot_enrichm).automatically_extracted_GO_terms.txt

					else
						# Assume GFF/GTF/GFF3 — extract Gene ID and GO terms
						echo "Detected GFF/GTF format. Extracting Gene IDs and GO terms..."
						zcat -f "$annot_enrichm" | awk -F'\t' '/GO:/ && !/^#/ {
						attrs = $9; gid = ""; go = ""
						if (attrs ~ /gene_id "/) { tmp = attrs; sub(/.*gene_id "/, "", tmp); sub(/".*/, "", tmp); gid = tmp }
						if (gid == "" && attrs ~ /ID=/) { tmp = attrs; sub(/.*ID=/, "", tmp); sub(/[;].*/, "", tmp); gid = tmp }
						if (gid == "" && attrs ~ /Parent=/) { tmp = attrs; sub(/.*Parent=/, "", tmp); sub(/[;].*/, "", tmp); gid = tmp }
						if (attrs ~ /[Oo]ntology_term/) { tmp = attrs; sub(/.*[Oo]ntology_term[= ]*"?/, "", tmp); sub(/"?[;].*/, "", tmp); sub(/"$/, "", tmp); go = tmp }
						if (gid != "" && go != "") { n = split(go, a, ","); for (i = 1; i <= n; i++) if (a[i] ~ /^GO:/) print gid "\t" a[i] }
					}' | sort -u > $output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annot_enrichm).automatically_extracted_GO_terms.txt
					fi

					annotation_go=$output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annot_enrichm).automatically_extracted_GO_terms.txt
					go_data_lines=$(grep "GO:" "$annotation_go" 2>/dev/null | wc -l)
					if [ "$go_data_lines" -lt 2 ]; then
						echo "WARNING: GO term extraction produced $go_data_lines valid entries. The input file may have an unexpected format. Skipping enrichment."
					else
						echo "Extracted $go_data_lines gene-GO associations. Running enrichr..."
						sed -i '1s/^/source_id\tComputed_GO_Process_IDs\n/' $annotation_go
						R_clusterProfiler_enrichr.R $annotation_go $output_folder/$name/final_results_reanalysis$index/RPKM_counts_genes.txt $output_folder/$name/final_results_reanalysis$index/DGE "^DGE_analysis_comp[0-9]+.txt$" &> clusterProfiler_enrichr_funct_enrichment.log
						echo "enrichr execution completed. Please double check the results and the log: clusterProfiler_enrichr_funct_enrichment.log"
					fi
				else
					echo "For $organism and the annotation $annotation_file, no GO or functional information found. Consider providing a GAF, GMT, or GO-annotated GFF/GTF via 'non_reference_funct_enrichm'"
				fi
			fi

			cd $output_folder/$name/final_results_reanalysis$index/DGE/
			error_files=$(grep -l "Err" *_funct_enrichment.log 2>/dev/null | sed 's/.txt_funct_enrichment.log//g' | sort | uniq)
			if [ -n "$error_files" ]; then
			    echo -e "\nFunctional enrichment analyses done!\nYou may want to check out the following logs, which seem to contain some errors:\n"
			    echo "$error_files"

			    # Automatic retry: re-run failed enrichment scripts once
			    echo -e "\nRetrying failed enrichment scripts once..."
			    cd $output_folder/$name/final_results_reanalysis$index/DGE/

			    # Retry clusterProfiler failures
			    failed_cp=$(grep -l "Err" clusterProfiler_*_funct_enrichment.log 2>/dev/null | sed 's/clusterProfiler_//g;s/_funct_enrichment.log//g' | sort | uniq)
			    for ff in $failed_cp; do
			        if [ -f "$ff" ]; then
			            echo "  Retrying clusterProfiler for $ff ..."
			            rm -rf $(echo $ff | sed 's/.txt$//')_funct_enrich_clusterProfiler
			            R_clusterProfiler_analyses_parallel.R $PWD $organism "1" $clusterProfiler_method $clusterProfiler_full $aPEAR_execution "^${ff}$" $clusterProfiler_universe $clusterProfiler_minGSSize $clusterProfiler_maxGSSize &> clusterProfiler_${ff}_funct_enrichment_retry.log
			        fi
			    done

			    # Retry autoGO+Panther failures
			    failed_ago=$(grep -l "Err" autoGO_panther_*_funct_enrichment.log 2>/dev/null | sed 's/autoGO_panther_//g;s/_funct_enrichment.log//g' | sort | uniq)
			    for ff in $failed_ago; do
			        if [ -f "$ff" ]; then
			            echo "  Retrying autoGO+Panther for $ff ..."
			            R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism "1" $databases_function $ff $panther_method $auto_panther_log &> autoGO_panther_${ff}_funct_enrichment_retry.log
			        fi
			    done

			    # Check again after retry
			    retry_errors=$(grep -l "Err" *_funct_enrichment_retry.log 2>/dev/null | sed 's/_funct_enrichment_retry.log//g' | sort | uniq)
			    if [ -n "$retry_errors" ]; then
			        echo -e "\nAfter retry, the following still have errors:"
			        echo "$retry_errors"
			    else
			        echo -e "\nRetry completed successfully — no more errors detected."
			    fi
			else
			    echo -e "\nFunctional enrichment analyses done! No errors detected in logs."
			fi

			# Add to the tables of functional enrichment the number of genes up/down:
			cd $output_folder/$name/final_results_reanalysis$index/
			files_to_process=$(find . \( -name "*.txt" -o -name "*.tsv" -o -name "*.csv" \) | grep funct | grep -v _err.txt)
			if [ -n "$files_to_process" ]; then
				enrichment_results_found="yes"
				cd $output_folder/$name/final_results_reanalysis$index/DGE/
				echo "Formatting $(echo $files_to_process | wc -w) functional enrichment result file(s)..."
				echo $files_to_process | parallel --joblog R_enrich_format_analyses_parallel_log_parallel.txt -j $cores "file={}; R_enrich_format.R \"\$file\" \$(echo \"\$file\" | sed 's,DGE/.*,DGE/,g')\$(echo \"\$file\" | sed 's,.*DGE_analysis_comp,DGE_analysis_comp,g;s,_pval.*,,g;s,_fdr.*,,g;s,_funct.*,,g;s,_cluster.*,,g' | sed 's,.txt,,g').txt $organism $rev_thr" &> $PWD/enrichment_format.log
			else
				echo "No functional enrichment results found. Report will not be rendered."
			fi
		fi

		# Render functional enrichment HTML report (self-contained), only if results exist
		if [[ "$functional_enrichment_analyses" != "no" ]] && [[ "$enrichment_results_found" == "yes" ]] && [ -d "$output_folder/$name/final_results_reanalysis$index/DGE" ]; then
			echo -e "Rendering functional enrichment HTML report..."
			Rscript $CURRENT_DIR/scripts/render_enrichment_report.R \
				"$output_folder/$name/final_results_reanalysis$index/DGE" \
				"$name" \
				"$organism" &> "$output_folder/$name/enrichment_report_render.log"
			if [ $? -eq 0 ] && [ -f "$output_folder/$name/final_results_reanalysis$index/DGE/functional_enrichment_report.html" ]; then
				echo "Done! Report: $output_folder/$name/final_results_reanalysis$index/DGE/functional_enrichment_report.html"
			else
				echo "WARNING: Functional enrichment report rendering failed. Check enrichment_report_render.log"
			fi
		fi
	done
	export debug_step="all"
_log_step "Step_6_Enrichment" "end"
	echo -e "\n\nSTEP 6: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP 7. Annotation: Tables of DEGs, lists of genes, etc
if [[ $debug_step == "all" || $debug_step == "step7" ]]; then
	echo -e "\n\nSTEP 7: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_7_Annotation" "start"
	echo -e "\n\nAnnotating list of genes...\n\n"
	for index in "${!array[@]}"; do
		# Export the annotation file path so R scripts can use it for ENSEMBL->Symbol mapping
		export ANNOTATION_FILE="${array[index]}"
		# All the tables that contain list of genes, annotate them:
		R_annotate_genes.R $output_folder/$name/final_results_reanalysis$index/ "^DGE_analysis_comp\\d+\\.txt$|^DGE_limma_timecourse_T\\d+_vs_T\\d+\\.txt$|mfuzz_elements_clusters|counts|WGCNA_all_modules_|STRINGdb_all_modules_" $organism			

		if [[ "$bed_mode" == "yes" ]]; then
			# All the tables of DEGs, provide bed files for direct upload in genome browser
			cd $output_folder/$name/final_results_reanalysis$index/
			for file in $(find . -name "DGE_analysis_comp*" | egrep "_fdr_05.txt$|_pval_05.txt$"); do
				cut -f1 "$file" | parallel -j $((cores*3)) "gene={}; foldchange=\$(grep -i \"\$gene\" \"$file\" | cut -f3 | sed -n 's/\(.*[.,][0-9]\{2\}\).*/\1/p'); \
															grep -i \"=\$gene\" \"$annotation_file\" | head -1 | awk -v id=\"\$gene\" -v fc=\"\$foldchange\" '{ print \$1\"\\t\"\$4\"\\t\"\$5\"\\t\"id\"_\"fc\"\\t.\t\"\$7 }' >> \"$file.bed\""
			done
		fi
	done
	export debug_step="all"
_log_step "Step_7_Annotation" "end"
	echo -e "\n\nSTEP 7: DONE\nCurrent date/time: $(date)\n\n"
fi


###### STEP 8. Sum up results in a sphinx report
if [[ $debug_step == "all" || $debug_step == "step8" ]]; then
	_log_step "Step_8_Report" "start"
	sphinx_report.sh $output_folder/$name $name
_log_step "Step_8_Report" "end"
 	echo -e "\n\nSTEP 8: Final report DONE\nCurrent date/time: $(date)\n\n"
fi


###### STEP 9. Tidy up, prepare for storage if final results have been created and the number of aligned files is equal to the numbers of samples, rename folders, convert tables to xlsx if required... etc
# Compress the folders
if [[ $debug_step == "all" || $debug_step == "step9" ]]; then
	echo -e "\n\nSTEP 9: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_9_Cleanup" "start"
	echo -e "\n\nTidying up, removing empty folders, temp files, compressing...\n\n"

	# Remove decompressed reference files from the indexes subfolder
	if [ -d "$output_folder/$name/indexes" ]; then
		for _decomp_file in $(find $output_folder/$name/indexes -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.gtf' -o -name '*.gff' \) 2>/dev/null); do
			rm -f "$_decomp_file"
		done
	fi


	cd $output_folder/$name/ && find . -type f \( -name "*_fdr_05.txt" -o -name "*_logneg.txt" -o -name "*_logpos.txt" \) -exec rm -f {} +
	if [ "$convert_tables_excel" == "yes" ]; then
		R_convert_tables.R $output_folder/$name/ $cores "log_parallel|jquery|bamqc|rnaseqqc|samtools|strand" > R_convert_tables.log 2>&1
	fi

	for index in "${!array[@]}"; do
	 	cd $output_folder/$name/final_results_reanalysis$index/DGE/
		find . -type d -empty -delete
		folders_funct=$(find . -type d \( -name "*_autoGO" -o -name "*_clusterProfiler" -o -name "*_panther" \) -o -type f \( -name "*_funct_enrichment.log" -o -name "funct_enrich_*" \))
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
	if [ "$tidy_tmp_files" == "yes" ]; then
		num_raw_files=$(cat $output_folder/$name/miARma_out0/Pre_fastqc_results/list_of_files.txt | grep -c "fastq.gz")
		for index in "${!array[@]}"; do
			if [ -d "$output_folder/$name/final_results_reanalysis$index" ] && [[ $(ls $output_folder/$name/miARma_out$index/$aligner\_results | egrep -c ".bam$") -eq $num_raw_files || $(ls $output_folder/$name/miARma_out$index/$aligner\_results | egrep -c ".bam$") -eq $((num_raw_files / 2)) ]]; then
				echo -e "\nTidying up...\n"
				cd $seqs_location
				echo "After execution, raw reads have been removed for the sake of efficient storage. These were... " > readme
				ls -lh >> readme
				ls | grep -v readme | xargs -r rm -rf

				cd $output_folder/$name/miARma_out$index/$aligner\_results
				echo "For the sake of efficiente storage: samtools view -@ cores -T ref_genome -C -o xxx.bam.cram xxx.bam && rm xx.bam" >> conversion_bam_to_cram.txt
				find . -type f -name "*.bam" | parallel --verbose -j $number_parallel --max-args 1 samtools view -T $reference_genome -C -@ $((cores / number_parallel)) -o {}.cram {}
				rm -rf $(ls | egrep ".bam$") $TMPDIR
			fi
		done
	fi

	### Deploy igvShinyApp.R to final results folders
	if [ -f "$CURRENT_DIR/scripts/igvShinyApp.R" ]; then
		IFS=', ' read -r -a array_annot <<< "$annotation"
		for index in "${!array_annot[@]}"; do
			final_dir_igv=$(find $output_folder/$name -maxdepth 1 -type d -name "final_results_reanalysis${index}_*" | head -1)
			if [ -z "$final_dir_igv" ]; then
				final_dir_igv="$output_folder/$name/final_results_reanalysis${index}_$(basename $output_folder)"
			fi
			if [ -d "$final_dir_igv" ]; then
				bw_dir="$output_folder/$name/miARma_out${index}/${aligner}_results"
				gtf_for_igv="${array_annot[index]}"
				cp "$CURRENT_DIR/scripts/igvShinyApp.R" "$final_dir_igv/igvShinyApp.R"
				sed -i "s|/path/to/reference_genome.fa|${reference_genome}|g" "$final_dir_igv/igvShinyApp.R"
				sed -i "s|/path/to/annotation.gtf|${gtf_for_igv}|g"          "$final_dir_igv/igvShinyApp.R"
				sed -i "s|/path/to/bigwig_folder/|${bw_dir}/|g"              "$final_dir_igv/igvShinyApp.R"
				sed -i "s|GENOME_NAME_PLACEHOLDER|${organism}|g"             "$final_dir_igv/igvShinyApp.R"
			fi
		done
	fi

	export debug_step="all"
_log_step "Step_9_Cleanup" "end"
	echo -e "\n\nSTEP 9: DONE\nCurrent date/time: $(date)\n\n"
fi

# Final Gantt chart with all steps (rendered after everything completes)
if [ -f "$STEP_TIMES_FILE" ]; then
	for qc_dir in $(find "$output_folder/$name" -maxdepth 2 -type d -name "QC_and_others" 2>/dev/null); do
		gantt_out="$qc_dir/pipeline_gantt_chart.pdf"
		echo -e "\nRendering final pipeline Gantt chart..."
		Rscript $CURRENT_DIR/scripts/R_gantt_chart.R "$STEP_TIMES_FILE" "$gantt_out" 2>&1 || \
			echo "WARNING: Gantt chart rendering failed."
	done
fi

echo -e "\n\n\nALL STEPS DONE! Best wishes\n\n\n"
