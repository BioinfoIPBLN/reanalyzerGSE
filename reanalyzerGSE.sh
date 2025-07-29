#!/bin/bash
start=`date +%s`
echo -e "\nCurrent time: $(date)\n"
base64 -d <<<"CiAgX19fICBfX19fICBfX19fICAgICAgX19fXyAgX19fXyAgIF9fICAgX18gXyAgIF9fICAgX18gICAgXyAgXyAgX19fXyAgX19fXyAgX19fXyAKIC8gX18pLyBfX18pKCAgX18pICAgICggIF8gXCggIF9fKSAvIF9cICggICggXCAvIF9cICggICkgICggXC8gKShfXyAgKSggIF9fKSggIF8gXAooIChfIFxcX19fIFwgKSBfKSAgICAgICkgICAvICkgXykgLyAgICBcLyAgICAvLyAgICBcLyAoXy9cICkgIC8gIC8gXy8gICkgXykgICkgICAvCiBcX19fLyhfX19fLyhfX19fKSAgICAoX19cXykoX19fXylcXy9cXy9cXylfXylcXy9cXy9cX19fXy8oX18vICAoX19fXykoX19fXykoX19cXykKCmJ5IEJpb2luZm9ybWF0aWNzIFVuaXQJCQkJSVBCTE4tQ1NJQy4gMjAyMwoKYmlvaW5mb3JtYXRpY2FAaXBiLmNzaWMuZXMJCSAgICAgICAgaHR0cHM6Ly9naXRodWIuY29tL0Jpb2luZm9JUEJMTi9yZWFuYWx5emVyR1NFCgo="
echo "doi.org/10.1101/2023.07.12.548663v2"
echo -e "reanalyzerGSE V3.1.1\n\n"

###### 0. Define arguments and variables:
### Export a string with command options and an array with arguments and deal with them in parse_options.sh
export options=$@
export arguments=($options)

CURRENT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source $CURRENT_DIR/scripts/parse_options.sh
if [ $? -ne 0 ]; then
	echo "Exiting..."; exit 1
fi

###### STEP 1. Download info from GEO and organize metadata and so:
if [[ $debug_step == "all" || $debug_step == "step1" ]]; then
	rm -rf $output_folder/*
	if [[ $input == G* ]]; then
		echo -e "\n\nSTEP 1: Starting...\nCurrent date/time: $(date)\n\n"
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
		echo -e "\nSTEP 1 DONE. Current time: $(date)\n"
	
 	fi
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
			IFS=', ' read -r -a arr2 <<< "$(ls | egrep .fastq.gz$ | sed 's,1.fastq.gz,,g;s,2.fastq.gz,,g' | sort | uniq | tr '\n' ',')"
			desired_number=${arr[1]}
			apply_random_shift() {
				Rscript -e '
				  modify_number <- function(number) {
				    percentage <- runif(1, 0, 10)  # Random % between 0 and 10
				    change <- ifelse(runif(1) < 0.5, -1, 1)  # Randomly add or subtract
				    change_amount <- number * (percentage / 100) * change
				    return(number + change_amount)
				  }
				  cat(modify_number('"$1"'), "\n")
				'
			}
   			desired_numbers=$(while IFS=$'\t' read -r col1 col2 col3; do
					    desired_number_rand=$(apply_random_shift $desired_number)
					    if (( col3 < desired_number_rand )); then
					        result=$col2
					    else
					        result=$((col2 * desired_number_rand / col3))
					    fi
					    echo $result
					 done < <(sed '1d' ${arr[0]}))
			IFS=', ' read -r -a arr3 <<< "$(echo $desired_numbers | tr ' ' ',')"
   			subsample_reads() {
				files=$(ls | grep $1)
				number=$2
				for file in $files; do seqtk sample -s 123 "$file" "$number" > "${file}_subsamp"; done
			}
			export -f subsample_reads
			parallel --verbose -j $cores_reads_to_subsample subsample_reads {} ::: "${arr2[@]}" :::+ "${arr3[@]}" # 10 only because of RAM
			rm $(ls | grep -v subsamp); for file in $(ls); do mv $file $(echo $file | sed 's,_subsamp,,g;s,.gz,,g'); done
			pigz --best -p $cores * # gz was lost with seqtk sample
			echo -e "\nSubsampling (+-10%) completed...\n"
		fi
		echo -e "\n\nSTEP 1: DONE\nCurrent date/time: $(date)\n\n"
 	fi
	export debug_step="all"
fi


### STEP 1. Process if not required to download from NCBI/GEO the metadata and raw reads provided locally:
if [[ $debug_step == "all" || $debug_step == "step1b" ]]; then
	mkdir -p $TMPDIR
	if [[ $input == /* ]]; then
		echo -e "\n\nSTEP 1b: Preparing the raw reads and metadata provided locally...\nCurrent date/time: $(date)\n\n"
  		seqs_location=$output_folder/$name/raw_reads
		rm -rf $seqs_location # I'm now removing the seqs_location at the beginning of this section, in the context of the new system of resuming by -Dm stepx, so this should always be done
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
		if [ ! -z "$number_reads" ]; then
			echo -e "\nSubsampling...\n"
			# From the input parameter by the user, obtain a random number allowing a +- 10% window:
			IFS=', ' read -r -a arr <<< "$number_reads"
			IFS=', ' read -r -a arr2 <<< "$(ls | egrep .fastq.gz$ | sed 's,1.fastq.gz,,g;s,2.fastq.gz,,g' | sort | uniq | tr '\n' ',')"			
			desired_number=${arr[1]}
			apply_random_shift() {
				Rscript -e '
				  modify_number <- function(number) {
				    percentage <- runif(1, 0, 10)  # Random % between 0 and 10
				    change <- ifelse(runif(1) < 0.5, -1, 1)  # Randomly add or subtract
				    change_amount <- number * (percentage / 100) * change
				    return(number + change_amount)
				  }
				  cat(modify_number('"$1"'), "\n")
				'
			}
   			desired_numbers=$(while IFS=$'\t' read -r col1 col2 col3; do
					    desired_number_rand=$(apply_random_shift $desired_number)
					    if (( col3 < desired_number_rand )); then
					        result=$col2
					    else
					        result=$((col2 * desired_number_rand / col3))
					    fi
					    #echo "desired_rand is $desired_number_rand"
					    echo $result
					 done < <(sed '1d' ${arr[0]}))
			IFS=', ' read -r -a arr3 <<< "$(echo $desired_numbers | tr ' ' ',')"
   			subsample_reads() {
				files=$(ls | grep $1)
				number=$2
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
			echo -e "The used conditions are (must match the list above):\n$(echo $design_input | sed 's_,_\n_g;s,/,\n\n,g')"
		fi
		mkdir -p $output_folder/$name/GEO_info/
		paste <(ls $seqs_location | egrep '.fq|.fastq' | sed "s/_1.fastq.gz//" | sed "s/_2.fastq.gz//" | uniq) <(paste -d'_' <(ls $seqs_location | egrep '.fq|.fastq' | egrep '.fq|.fastq' | sed "s/_1.fastq.gz//" | sed "s/_2.fastq.gz//" | uniq) <(echo $design_input | sed 's*/*\t*g'| cut -f1 | sed 's*,*\n*g')) <(echo $design_input | sed 's*/*\t*g'| cut -f1 | sed 's*,*\n*g') > $output_folder/$name/GEO_info/samples_info.txt
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
Rscript -e "organism <- '${organism}'; assemblies <- rentrez::entrez_summary(db='assembly', id=rentrez::entrez_search(db='assembly', term=paste0(organism, '[orgn]'))\$ids[1]); cat(paste(paste0('\n\nNCBI current assembly info: ', date()), assemblies\$assemblyname, assemblies\$assemblyaccession, assemblies\$submissiondate, '\n', sep='\n'))"
organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's/ \+/_/g;s/__*/_/g') # Get again the organism in case it has been manually modified... and without spaces...

### STEP 1. Deal with fastp if required:
if [ "$fastp_mode" == "yes" ]; then
	echo -e "\n\nSTEP 1: Preprocessing with fastp...\nCurrent date/time: $(date)\n\n"
 	mkdir -p $output_folder/$name/fastp_out
	cd $output_folder/$name/fastp_out
	if [[ "$(find $output_folder/$name -name library_layout_info.txt | xargs cat)" == "SINGLE" ]]; then
		for f in $(ls -d $seqs_location/*); do echo "Processing $f"; fastp --in1 $f --out1 $f\_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --disable_adapter_trimming --thread $cores -h $(basename $f)\_report.html -j $(basename $f)\_report.json &>> $(basename $f)\_fastp_out.log; done
	elif [[ "$(find $output_folder/$name -name library_layout_info.txt | xargs cat)" == "PAIRED" ]]; then
		for f in $(ls -d $seqs_location/* | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq); do echo "Processing $f"; fastp --in1 $f\_1.fastq.gz --in2 $f\_2.fastq.gz --out1 $f\_1.fastq.gz_fastp.fastq.gz --out2 $f\_2.fastq.gz_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --disable_adapter_trimming --thread $cores -h $(basename $f)\_report.html -j $(basename $f)\_report.json &>> $(basename $f)\_fastp_out.log; done
	fi
fi

if [ "$fastp_adapter" == "yes" ]; then
	echo -e "\n\nSTEP 1: Preprocessing with fastp to remove adapters...\nCurrent date/time: $(date)\n\n"
 	mkdir -p $output_folder/$name/fastp_out
	cd $output_folder/$name/fastp_out
	if [[ "$(find $output_folder/$name -name library_layout_info.txt | xargs cat)" == "SINGLE" ]]; then
		for f in $(ls -d $seqs_location/*); do echo "Processing $f"; fastp --in1 $f --out1 $f\_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --thread $cores -h $(basename $f)\_report.html -j $(basename $f)\_report.json &>> $(basename $f)\_fastp_out.log; done
	elif [[ "$(find $output_folder/$name -name library_layout_info.txt | xargs cat)" == "PAIRED" ]]; then
		for f in $(ls -d $seqs_location/* | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq); do echo "Processing $f"; fastp --in1 $f\_1.fastq.gz --in2 $f\_2.fastq.gz --out1 $f\_1.fastq.gz_fastp.fastq.gz --out2 $f\_2.fastq.gz_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --detect_adapter_for_pe --thread $cores -h $(basename $f)\_report.html -j $(basename $f)\_report.json &>> $(basename $f)\_fastp_out.log; done
	fi
elif [[ $fastp_adapter == /* ]]; then
	echo -e "\n\nSTEP 1: Preprocessing with fastp to remove adapters...\nCurrent date/time: $(date)\n\n"
 	mkdir -p $output_folder/$name/fastp_out
	cd $output_folder/$name/fastp_out
	if [[ "$(find $output_folder/$name -name library_layout_info.txt | xargs cat)" == "SINGLE" ]]; then
		for f in $(ls -d $seqs_location/*); do echo "Processing $f"; fastp --in1 $f --out1 $f\_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --adapter_fasta $fastp_adapter --thread $cores -h $(basename $f)\_report.html -j $(basename $f)\_report.json &>> $(basename $f)\_fastp_out.log; done
	elif [[ "$(find $output_folder/$name -name library_layout_info.txt | xargs cat)" == "PAIRED" ]]; then
		for f in $(ls -d $seqs_location/* | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq); do echo "Processing $f"; fastp --in1 $f\_1.fastq.gz --in2 $f\_2.fastq.gz --out1 $f\_1.fastq.gz_fastp.fastq.gz --out2 $f\_2.fastq.gz_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --adapter_fasta $fastp_adapter --thread $cores -h $(basename $f)\_report.html -j $(basename $f)\_report.json &>> $(basename $f)\_fastp_out.log; done
	fi
fi

if [ "$fastp_trimming" != "none" ]; then
	echo -e "\n\nSTEP 1: Preprocessing with fastp and trimming...\nCurrent date/time: $(date)\n\n"
 	mkdir -p $output_folder/$name/fastp_out
	cd $output_folder/$name/fastp_out
	IFS=', ' read -r -a arrfastp <<< "$fastp_trimming"
	if [[ "$(find $output_folder/$name -name library_layout_info.txt | xargs cat)" == "SINGLE" ]]; then
		for f in $(ls -d $seqs_location/*); do echo "Processing $f"; fastp --in1 $f --out1 $f\_fastp.fastq.gz --dont_overwrite --dont_eval_duplication  --trim_front1 "${arrfastp[0]}" --trim_tail1 "${arrfastp[1]}" --thread $cores -z $compression_level -h $(basename $f)\_report.html -j $(basename $f)\_report.json &>> $(basename $f)\_fastp_out.log; done
	elif [[ "$(find $output_folder/$name -name library_layout_info.txt | xargs cat)" == "PAIRED" ]]; then
		for f in $(ls -d $seqs_location/* | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq); do echo "Processing $f"; fastp --in1 $f\_1.fastq.gz --in2 $f\_2.fastq.gz --out1 $f\_1.fastq.gz_fastp.fastq.gz --out2 $f\_2.fastq.gz_fastp.fastq.gz --dont_overwrite --dont_eval_duplication  --trim_front1 "${arrfastp[0]}" --trim_tail1 "${arrfastp[1]}" --thread $cores -z $compression_level -h $(basename $f)\_report.html -j $(basename $f)\_report.json &>> $(basename $f)\_fastp_out.log; done
	fi
fi

if [ $(ls -d $seqs_location/* | grep -c _fastp.fastq.gz) -gt 0 ]; then
	for f in $(ls -d $seqs_location/* | grep _fastp.fastq.gz); do mv $f $(echo $f | sed 's,.fastq.gz_fastp.fastq.gz,.fastq.gz,g'); done
 	echo "Files in $seqs_location have been successfully processed by fastp!"
fi


### STEP 2. Decontamination if required:
if [[ $debug_step == "all" || $debug_step == "step2" ]]; then
	if [ ! -z "$kraken2_databases" ]; then
  		echo -e "\n\nSTEP 2: Decontamination starting with Kraken2...\nCurrent date/time: $(date)\n\n"
    		rm -rf $output_folder/$name/raw_reads_k2 # I'm now removing the seqs_location at the beginning of this section, in the context of the new system of resuming by -Dm stepx, so this should always be done
		mkdir -p $output_folder/$name/raw_reads_k2
		cd $output_folder/$name/raw_reads_k2

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
  		echo -e "\n\nSTEP 2: DONE\nCurrent date/time: $(date)\n\n"
	fi
	if [ ! -z "$sortmerna_databases" ]; then
		echo -e "\n\nSTEP 2: Decontamination starting with sortmerna...\nCurrent date/time: $(date)\n\n"
  		if [ ! -d "$CURRENT_DIR/indexes/$(basename $sortmerna_databases)_sortmerna_index" ]; then
			echo "Indexing the provided $sortmerna_databases ..."
   			sortmerna --index 1 --ref $sortmerna_databases --workdir $CURRENT_DIR/indexes/$(basename $sortmerna_databases)_sortmerna_index --threads $cores &>> $output_folder/$name/sortmerna_out/index.log
		fi
		mkdir -p $output_folder/$name/sortmerna_out; mkdir -p $seqs_location\_sortmerna; cd $seqs_location\_sortmerna
  		rm -rf $seqs_location\_sortmerna/* $output_folder/$name/sortmerna_out/* # I'm now removing also at the beginning of this section, in the context of the new system of resuming by -Dm stepx, so this should always be done
    		
		# with the argument --paired_out, only the pairs where both reads are coincident (aligning to rRNA or not, are included in the results)
		# I don't include it, so the numbers are exactly the ones in the log, and the properly paired reads can be dealt with later on the mapping
		echo -e "\nExecuting sortmerna and fastqc of the new reads...\n"
		if [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "PAIRED" ]]; then
			for f in $(ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$" | sed 's,_1.fastq.gz,,g;s,_2.fastq.gz,,g' | sort | uniq); do echo "Processing $f"; sortmerna --idx-dir $CURRENT_DIR/indexes/$(basename $sortmerna_databases)_sortmerna_index/idx --ref $sortmerna_databases --reads $seqs_location/${f}_1.fastq.gz --reads $seqs_location/${f}_2.fastq.gz --workdir ${f}_sortmerna_out --fastx --threads $cores --out2 --aligned ${f}_rRNA --other ${f}_no_rRNA -v &>> $output_folder/$name/sortmerna_out/${f}_out.log; done
			rm -rf $(ls | egrep "_sortmerna_out$")
		elif [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "SINGLE" ]]; then
			for f in $(ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$"); do echo "Processing $f"; sortmerna --idx-dir $CURRENT_DIR/indexes/$(basename $sortmerna_databases)_sortmerna_index/idx --ref $sortmerna_databases --reads $seqs_location/$f --workdir ${f}_sortmerna_out --fastx --threads $cores --aligned ${f}_rRNA --other ${f}_no_rRNA -v &>> $output_folder/$name/sortmerna_out/${f}_out.log; done
			rm -rf $(ls | egrep "_sortmerna_out$")
		fi
		fastqc -t $cores *.fq.gz
  		mkdir out_noRNA; cd out_noRNA; for f in (ls ../* | grep no_RNA*.fq.gz); do ln -sf $f $(echo $f | sed 's,.fq.gz,.fastq.gz,g'); done; export $seqs_location=$seqs_location\_sortmerna/out_noRNA
  		echo -e "\n\nSTEP 2: DONE\nCurrent date/time: $(date)\n\n"
 	fi
	export debug_step="all"
fi


### STEP3a. Prepare the data and info for running miARma-seq:
if [[ $debug_step == "all" || $debug_step == "step3a" ]]; then
	echo -e "\n\nSTEP 3a: Starting...\nCurrent date/time: $(date)\n\n"
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
			kall_idx=$CURRENT_DIR/indexes/${organism}_kallisto_idx
			if [ ! -s "$kall_idx" ]; then
				kallisto index -i $output_folder/$name/indexes/${organism}_kallisto_idx $transcripts &> $output_folder/$name/indexes/${organism}_kallisto_idx.log
				kall_idx=$output_folder/$name/indexes/${organism}_kallisto_idx
			fi
			echo -e "\nPredicting strandness on two random sample...\n"
			mkdir -p $output_folder/$name/strand_prediction/salmon_out; mkdir -p $output_folder/$name/strand_prediction/how_are_we_stranded_here_out
			if [[ $(find $output_folder/$name -name library_layout_info.txt | xargs cat) == "SINGLE" ]]; then
				salmon quant -i $salmon_idx -l A -r $seqs_location/$(ls $seqs_location | l | head -1) -p $cores -o $output_folder/$name/strand_prediction/salmon_out/ --skipQuant &> $output_folder/$name/strand_prediction/salmon_out/salmon_out.log
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
				strand="reverse"
			elif [[ "$salmon_strand" == "SF" || "$salmon_strand" == "ISF" ]]; then
				strand="yes"
			elif [[ "$salmon_strand" == "U" || "$salmon_strand" == "IU" ]]; then
				strand="no"
			fi
			strand_second_opinion=$(cat $output_folder/$name/strand_prediction/how_are_we_stranded_here_out/check_strandedness_out.log | egrep "Data is likely |Data does not fall into" | sed 's,*Data is likely ,,g;s/([^)]*)//g')
			cd $output_folder/$name	
			echo "Salmon prediction 1: $salmon_strand" > $output_folder/$name/strand_info.txt
			echo "Salmon prediction 2: $strand" >> $output_folder/$name/strand_info.txt
			echo -e "how_are_we_stranded_here prediction: $strand_second_opinion" >> $output_folder/$name/strand_info.txt
			if [ $(egrep -c "reverse|yes|no" $output_folder/$name/strand_info.txt) -gt 0 ]; then
				echo "Please double check carefully, based on the kit used in the library preparation, the paper, the GEO entry... because this is crucial for quantification. Please rerun with the argument '-s' in the unlikely case that the prediction by salmon is not correct, or if the second opinion by how_are_we_stranded_here is different (if transcripts from GENCODE or any particular format are used, the latter option may fail to identify the data type though)"
				cat $output_folder/$name/strand_info.txt
				rm $(find $output_folder/$name/strand_prediction/how_are_we_stranded_here_out -type f -name "*.bam")
			else
				echo "Salmon to detect strandedness seems to have failed. This is not acceptable, plese double check or provide the parameter -s. Exiting..."
				exit 1
			fi
		else
  			echo $strand > $output_folder/$name/strand_info.txt
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
	export debug_step="all"
	echo -e "\n\nSTEP 3a: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP3b. Running miARma-seq:
# 2024: I've modified miARma RNA-seq mode to leverage GNU's parallel and increase speed, introduce limit RAM in aligners and multithreading index, replace the shebang with #!/usr/bin/env perl so it uses the PATH's/environment's one, etc...
# Eventually, WIP nicludes to also improve and integrate the rest of modules of miARma, such as adapter cutting, stats, miRNAs...
if [[ $debug_step == "all" || $debug_step == "step3b" ]]; then
	echo -e "\n\nSTEP 3b: Starting...\nCurrent date/time: $(date)\n\n"
	rm -rf $output_folder/$name/miARma_out*
	mkdir -p $TMPDIR
	# If the running is resumed in this step, the above has to be done
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/GEO_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi

	echo "Please double check all the parameters above, in particular the stranded or the reference genome files and annotation used. Proceeding with miARma execution in..."
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
 	echo -e "Processing output of miARma-seq, QC figures, plots, DGE if requested..."
	IFS=', ' read -r -a array2 <<< "$filter"
	for index in "${!array[@]}"; do
		annotation_file=${array[index]}
  		echo -e "R_process_reanalyzer_GSE.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index $genes ${array2[index]} $organism $target $differential_expr_soft $batch_format $covariables $covariables_format $deconvolution $differential_expr_comparisons $perform_differential_analyses $perform_volcano_venn $pattern_to_remove\n\n" >> $output_folder/$name/final_results_reanalysis$index/R_process_reanalyzer.log
    		R_process_reanalyzer_GSE.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index $genes ${array2[index]} $organism $target $differential_expr_soft $batch_format $covariables $covariables_format $deconvolution $differential_expr_comparisons $perform_differential_analyses $perform_volcano_venn $pattern_to_remove | tee -a $output_folder/$name/final_results_reanalysis$index/R_process_reanalyzer.log
    		R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter" "edgeR_object" "edgeR_object_norm" $pattern_to_remove
		if [[ -e "$output_folder/$name/final_results_reanalysis$index/counts_adjusted.txt" ]]; then
			echo -e "\n\nRemember that batch effect correction/covariables have been only provided to Combat-Seq/limma for visualization purposes, to include covariables in the DGE model after checking the visualization the argument -C will be used\n\n\nQC_PDF adjusted counts\n\nRemember that you have requested batch effect correction/count adjustment, so you have to mind the figures in this QC_PDF from ComBat-seq/limma counts...\n"
			R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter_adjusted" "edgeR_object_adjusted" "edgeR_object_norm_adjusted" $pattern_to_remove
		fi
		cd $output_folder/$name/final_results_reanalysis$index/DGE/
		tar -cf - $(ls | egrep ".RData$") | pigz -p $cores > allRData.tar.gz; rm -rf $(ls | egrep ".RData$")
	done
	export debug_step="all"
	echo -e "\n\nSTEP 4: DONE\nCurrent date/time: $(date)\n\n"
	if [[ "$perform_differential_analyses" == "no" ]]; then
		echo "Differential analyses not requested, exiting the pipeline..."; exit 1
	fi
fi


### STEP 5. Time course analyses if required
if [[ $debug_step == "all" || $debug_step == "step5" ]]; then	
	for index in "${!array[@]}"; do
		if [[ "$time_course" == "yes" ]]; then
			echo -e "\n\nSTEP 5: Starting...\nCurrent date/time: $(date)\n\n"
   			echo -e "\nPerforming time course analyses."
			R_process_time_course.R $output_folder/$name/final_results_reanalysis$index/ DGE_analysis_comp1.RData edgeR_object_norm $minstd $mestimate
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
		if [[ "$organism" == "Mus_musculus" || "$organism" == "Homo_sapiens" || "$organism" == "Mus musculus" || "$organism" == "Homo sapiens" ]]; then
			if [[ $network_analyses == "yes" ]]; then
				mkdir -p network_analyses && rm -rf network_analyses/* && cd network_analyses
				echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
    				echo -e "\nPerforming network analyses...\n"
				R_network_analyses.R $output_folder/$name/final_results_reanalysis$index/DGE/ $output_folder/$name/final_results_reanalysis$index/RM_counts_genes.txt "^DGE_analysis_comp[0-9]+.txt$" $taxonid &> network_analyses_funct_enrichment.log
			fi
			if [[ "$functional_enrichment_analyses" == "no" ]]; then
				echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
    				echo -e "\nSkipping functional enrichment analyses\n"
			else
				echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
    				echo -e "\nPerforming functional enrichment analyses for DEGs. The results up to this point are ready to use (including DEGs and expression, that are only lacking annotation). This step of funtional enrichment analyses may take long if many significant DEGs, comparisons, or analyses...\n"
				cd $output_folder/$name/final_results_reanalysis$index/DGE/
				ls | egrep "^DGE_analysis_comp[0-9]+.txt$" | parallel --joblog R_clusterProfiler_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_clusterProfiler_analyses_parallel.R $PWD $organism "1" $clusterProfiler_method $clusterProfiler_full $aPEAR_execution '^{}$' $clusterProfiler_universe $clusterProfiler_minGSSize $clusterProfiler_maxGSSize &> clusterProfiler_{}_funct_enrichment.log"
				echo -e "\nPerforming autoGO and Panther execution... this may take long if many genes or comparisons...\n"
				ls | egrep "^DGE_analysis_comp[0-9]+.txt$" | parallel --joblog R_autoGO_panther_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism "1" $databases_function {} $panther_method $auto_panther_log &> autoGO_panther_{}_funct_enrichment.log"
				if [[ "$time_course" == "yes" ]]; then
					cd $output_folder/$name/final_results_reanalysis$index/time_course_analyses
					ls | egrep "^DGE_limma_timecourse.*.txt$" | parallel --joblog R_clusterProfiler_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_clusterProfiler_analyses_parallel.R $PWD $organism "1" $clusterProfiler_method $clusterProfiler_full $aPEAR_execution '^{}$' $clusterProfiler_universe $clusterProfiler_minGSSize $clusterProfiler_maxGSSize &> clusterProfiler_{}_funct_enrichment.log"
					ls | egrep "^DGE_limma_timecourse.*.txt$" | parallel --joblog R_autoGO_panther_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism "1" $databases_function {} $panther_method $auto_panther_log &> autoGO_panther_{}_funct_enrichment.log"
				fi
			fi
		else
			echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
   			echo "Organism is $organism... Functional analyses apart from human/mouse is not fully supported yet"
			if [ $(egrep -c "GO:|Ontology|tology_term|tology term" $annotation_file) -gt 0 ]; then
				cd $output_folder/$name/final_results_reanalysis$index/DGE/
				echo "However, an automatic approach based on clusterProfiler's enrichr function and automatically extracted GO terms from the annotation can be applied for DEGs..."
				paste <(egrep "GO:|,GO:|Ontology|tology_term|tology term" $annotation_file | sed 's,.*ID=,,g;s,.*Parent=,,g;s,;.*,,g') <(egrep "GO:|,GO:|Ontology|tology_term|tology term" $annotation_file | sed 's,.*tology_term=,,g') | sort -t $'\t' -k1,1 -k2,2 | awk -F'\t' '!a[$1,$2]++' | awk -F'\t' '{ a[$1] = (a[$1] ? a[$1]","$2 : $2); } END { for (i in a) print i"\t"a[i]; }' | awk -F '\t' '{n=split($2,a,","); for (i=1; i<=n; i++) print $1,a[i]}' | uniq > $output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annotation_file).automatically_extracted_GO_terms.txt
				annotation_go=$output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annotation_file).automatically_extracted_GO_terms.txt
				sed -i '1s/^/source_id Computed_GO_Process_IDs\n/' $annotation_go
				if [ -s "$annotation_go" ]; then
					R_clusterProfiler_enrichr.R $annotation_go $output_folder/$name/final_results_reanalysis$index/RM_counts_genes.txt $output_folder/$name/final_results_reanalysis$index/DGE "^DGE_analysis_comp[0-9]+.txt$" &> clusterProfiler_enrichr_funct_enrichment.log
					echo "DONE. Please double check the attempt of executing enrichr with automatically detected GO terms from the annotation"
				fi
			else
				echo "For $organism and the annotation $annotation_file, it does not seem there's GO or functional information available..."
			fi
		fi
		cd $output_folder/$name/final_results_reanalysis$index/DGE/
		echo -e "\nFunctional enrichment analyses done!\nYou may want to check out the following logs, which seem to contain some errors:\n"
		grep Err* $(ls | egrep _funct_enrichment.log) | cut -d":" -f1 | sed 's,.txt_funct_enrichment.log,,g' | sort | uniq

		# Add to the tables of functional enrichment the number of genes up/down:
		cd $output_folder/$name/final_results_reanalysis$index/
		echo "Processing results of functional enrichment analyses if any, for example executing Revigo..."
		files_to_process=$(find . \( -name "*.txt" -o -name "*.tsv" -o -name "*.csv" \) | grep funct | grep -v _err.txt)
		if [ -n "$files_to_process" ]; then
			cd $output_folder/$name/final_results_reanalysis$index/DGE/
			echo $files_to_process | parallel --joblog R_enrich_format_analyses_parallel_log_parallel.txt -j $cores "file={}; R_enrich_format.R \"\$file\" \$(echo \"\$file\" | sed 's,DGE/.*,DGE/,g')\$(echo \"\$file\" | sed 's,.*DGE_analysis_comp,DGE_analysis_comp,g;s,_pval.*,,g;s,_fdr.*,,g;s,_funct.*,,g;s,_cluster.*,,g' | sed 's,.txt,,g').txt $organism $rev_thr" &> $PWD/enrichment_format.log
		else
			echo "No functional enrichment results found, exiting the pipeline..."; exit 1
		fi
	done
	export debug_step="all"
	echo -e "\n\nSTEP 6: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP 7. Annotation: Tables of DEGs, lists of genes, etc
if [[ $debug_step == "all" || $debug_step == "step7" ]]; then
	echo -e "\n\nSTEP 7: Starting...\nCurrent date/time: $(date)\n\n"
	echo -e "\n\nAnnotating list of genes...\n\n"
	for index in "${!array[@]}"; do
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
	echo -e "\n\nSTEP 7: DONE\nCurrent date/time: $(date)\n\n"
fi


###### STEP 8. Tidy up, prepare for storage if final results have been created and the number of aligned files is equal to the numbers of samples, rename folders, convert tables to xlsx if required... etc
# Compress the folders
if [[ $debug_step == "all" || $debug_step == "step8" ]]; then
	echo -e "\n\nSTEP 8: Starting...\nCurrent date/time: $(date)\n\n"
	echo -e "\n\nTidying up, removing empty folders, temp files, compressing...\n\n"
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
				rm $(ls | grep -v readme)

				cd $output_folder/$name/miARma_out$index/$aligner\_results
				echo "For the sake of efficiente storage: samtools view -@ cores -T ref_genome -C -o xxx.bam.cram xxx.bam && rm xx.bam" >> conversion_bam_to_cram.txt
				find . -type f -name "*.bam" | parallel --verbose -j $number_parallel --max-args 1 samtools view -T $reference_genome -C -@ $((cores / number_parallel)) -o {}.cram {}
				rm -rf $(ls | egrep ".bam$") $TMPDIR
			fi
		done
	fi	

	cd $output_folder/$name/ && rm $(find . -name "*_fdr_05.txt" -o -name "*_logneg.txt" -o -name "*_logpos.txt")
	if [ "$convert_tables_excel" == "yes" ]; then
		R_convert_tables.R $output_folder/$name/ $cores "log_parallel|jquery|bamqc|rnaseqqc|samtools|strand" > /dev/null 2>&1
	fi
	export debug_step="all"
	echo -e "\n\nSTEP 8: DONE\nCurrent date/time: $(date)\n\n"
fi


###### STEP 9. Sum up results in a sphinx report
if [[ $debug_step == "all" || $debug_step == "step9" ]]; then
	sphinx_report.sh $output_folder/$name $name
 	echo -e "\n\nSTEP 9: Final report DONE\nCurrent date/time: $(date)\n\n"
fi

echo -e "\n\n\nALL STEPS DONE! Best wishes\n\n\n"
