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
_log_step() { mkdir -p "$(dirname "$STEP_TIMES_FILE")" 2>/dev/null; echo -e "$1\t$(date +%s)\t$2" >> "$STEP_TIMES_FILE" 2>/dev/null; }

###### Step gating (resume with -Dm/debug_step, force-end with -Es/end_step):
# Ordered list of pipeline steps. MUST stay in sync with the block order below,
# including the optional/internal 'step1a_bis', so the after-comparison is correct.
STEP_ORDER=(step1 step1a step1a_bis step1b step1c step1d step2 step2b step3a step3b step4 step4b step5 step6 step7 step8 step9)
# Echo the 0-based position of a step in STEP_ORDER, or -1 if unknown.
_step_pos() { local t=$1 i=0; for s in "${STEP_ORDER[@]}"; do [ "$s" = "$t" ] && { echo "$i"; return; }; i=$((i+1)); done; echo -1; }
# Whether the block for step $1 should run: the resume gate (debug_step, which
# flips to "all" once the resume point is reached) AND the force-end gate (skip
# any step strictly after end_step, so the script falls through to the summary).
run_step() {
	local s=$1
	[[ "$debug_step" == "all" || "$debug_step" == "$s" ]] || return 1
	if [[ -n "$end_step" && "$end_step" != "none" && "$end_step" != "all" ]]; then
		[ "$(_step_pos "$s")" -gt "$(_step_pos "$end_step")" ] && return 1
	fi
	return 0
}
build_alignment_units() {
	unit_ini=(); unit_out=(); unit_fasta=(); unit_gtf=(); unit_reads=(); unit_ri=(); unit_label=()
	local _i; local _annot_arr=()
	IFS=', ' read -r -a _annot_arr <<< "$annotation"
	# In kallisto mode annotation is optional; ensure there is at least one unit
	if [[ "$aligner" == "kallisto" && ${#_annot_arr[@]} -eq 0 ]]; then _annot_arr=(""); fi
	if [ ! -z "$reference_genome_groups" ]; then
		for _i in "${!genome_group_labels[@]}"; do
			unit_ini+=("miarma_genome$_i.ini")
			unit_out+=("$output_folder/$name/miARma_out_genome$_i")
			unit_fasta+=("${genome_group_fastas[$_i]}")
			unit_gtf+=("${_annot_arr[0]}")
			unit_reads+=("$output_folder/$name/reads_group$_i")
			unit_ri+=("")
			unit_label+=("${genome_group_labels[$_i]}")
		done
	else
		for _i in "${!_annot_arr[@]}"; do
			unit_ini+=("miarma$_i.ini")
			unit_out+=("$output_folder/$name/miARma_out$_i")
			unit_fasta+=("$reference_genome")
			unit_gtf+=("${_annot_arr[$_i]}")
			unit_reads+=("$seqs_location")
			unit_ri+=("$reference_genome_index")
			unit_label+=("")
		done
	fi
}
# Validate -Es/end_step and its ordering against -Dm/debug_step (resume point).
if [[ -n "$end_step" && "$end_step" != "none" && "$end_step" != "all" ]]; then
	if [ "$(_step_pos "$end_step")" -lt 0 ]; then
		echo "Error: -Es/-end_step '$end_step' is not a valid step. Choose one of: ${STEP_ORDER[*]} (or 'none')."; exit 1
	fi
	if [[ "$debug_step" != "all" && "$debug_step" != "none" && "$(_step_pos "$debug_step")" -ge 0 \
	      && "$(_step_pos "$end_step")" -lt "$(_step_pos "$debug_step")" ]]; then
		echo "Error: -Es/-end_step '$end_step' comes before -Dm/-debug_module '$debug_step'. The end step must be at or after the resume step."; exit 1
	fi
	echo -e "\nThe pipeline will force-end after '$end_step' (-Es/-end_step); later steps will be skipped.\n"
fi

###### STEP 1. Download info from GEO and organize metadata and so:
if run_step step1; then
	rm -rf $output_folder/*
	# Re-create step_times.tsv after rm -rf (which deletes it)
	if [[ $debug_step == "all" ]]; then
		mkdir -p "$output_folder/$name"
	fi
	if [ -n "$options_file" ] && [ -f "$options_file" ]; then
		mkdir -p "$output_folder/$name"
		cp -f "$options_file" "$output_folder/$name/" 2>/dev/null || true
	fi
	_log_step "Step_1_Download" "start"
	if [[ $input == G* ]]; then
		echo -e "\n\nSTEP 1: Starting...\nCurrent date/time: $(date)\n\n"
 	### Download info:
		echo -e "\nDownloading info from GEO for $input...\n"
		R_download_GEO_info.R $input $output_folder
		arrIN=(${input//,/ })
		input=$(for a in "${arrIN[@]}"; do echo "$a"; done | sort | tr '\n' '_' | sed 's,_$,,g')

	### Get metadata and process the info:
		echo -e "\nProcessing and downloading more data...\n"
		cd $output_folder/$name/reads_study_info
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
			for i in $(cat srx_ids.txt); do
				srr=""
				for attempt in 1 2 3 4; do
					delay=$(( (attempt - 1) * 10 ))
					[ $delay -gt 0 ] && sleep $delay
					srr=$(esearch -db sra -query "$i" 2>/dev/null | esummary 2>/dev/null | xtract -pattern DocumentSummary -element Run@acc 2>/dev/null)
					[ -n "$srr" ] && break
				done
				[ -n "$srr" ] && echo "$srr" >> srr_ids.txt
			done
			rm srx_ids.txt
		fi
		if test -f srr_ids.txt; then
			paste -d$'\t' srr_ids.txt sample_names.txt $(for f in $(ls | grep full); do echo $f" "$(sort $f | uniq -c | wc -l); done | grep -v " 1" | cut -d" " -f1 | head -1) > samples_info.txt # the third column is a design containing at least more than one element...
			sed -i -r -e 's/[^[:alnum:]_\t]/_/g' -e 's/_\+/_/g' -e 's/(_[^_]*)\1+/\1/g' -e 's/_([0-9]+)/-\1/g' samples_info.txt # Should be redundant but make sure to remove special characters from the sample names and _1/_2... it's crucial for later steps such as fastqc and miarma-seq
			
			# Filter metadata files if GSM_filter is provided so fastq-dl and AI metadata streamlining only process requested samples
			if [ -n "$GSM_filter" ] && [ -f gsm_manual_filter.txt ]; then
				echo -e "\nFiltering metadata and SRR accession download list for specified GSM(s): $GSM_filter..."
				python3 -c "
import os, sys
gsm_list = [g.strip() for g in open('gsm_manual_filter.txt').read().strip().split(',') if g.strip()]
if not gsm_list:
    sys.exit(0)

# Filter samples_info.txt
if os.path.exists('samples_info.txt'):
    lines = [l for l in open('samples_info.txt').read().splitlines() if l.strip()]
    kept = [l for l in lines if any(gsm in l for gsm in gsm_list)]
    if kept:
        open('samples_info.txt', 'w').write('\n'.join(kept) + '\n')
        # Update srr_ids.txt and sample_names.txt from filtered samples_info.txt
        srrs = [l.split('\t')[0] for l in kept]
        snames = [l.split('\t')[1] for l in kept]
        open('srr_ids.txt', 'w').write('\n'.join(srrs) + '\n')
        open('sample_names.txt', 'w').write('\n'.join(snames) + '\n')

        # Filter design files
        for f in os.listdir('.'):
            if f.startswith('design_possible_full_'):
                dlines = [l for l in open(f).read().splitlines() if l.strip()]
                # Keep matching indices from original lines
                kept_d = [dlines[i] for i, l in enumerate(lines) if any(gsm in l for gsm in gsm_list) and i < len(dlines)]
                if kept_d:
                    open(f, 'w').write('\n'.join(kept_d) + '\n')
            elif f.startswith('design_possible_') and not f.startswith('design_possible_full_'):
                dlines = [l for l in open(f).read().splitlines() if l.strip()]
                kept_d = list(dict.fromkeys([lines[i].split('\t')[2] if len(lines[i].split('\t'))>2 else dlines[i] for i, l in enumerate(lines) if any(gsm in l for gsm in gsm_list) and i < len(lines)]))
                if kept_d:
                    open(f, 'w').write('\n'.join(kept_d) + '\n')

# Filter phenodata_extracted.txt if present
if os.path.exists('phenodata_extracted.txt'):
    plines = [l for l in open('phenodata_extracted.txt').read().splitlines() if l.strip()]
    if plines:
        pkept = [l for i, l in enumerate(plines) if i == 0 or any(gsm in l for gsm in gsm_list)]
        if pkept:
            open('phenodata_extracted.txt', 'w').write('\n'.join(pkept) + '\n')
"
			fi

			# AI-assisted metadata streamlining (if LLM endpoint enabled and metadata AI requested)
			if [ -n "$llm_endpoint" ] && [ "${ai_do_metadata:-1}" = 1 ] && command -v llm_simplify_metadata.py >/dev/null 2>&1; then
				echo -e "\nRunning AI-assisted sample and condition name streamlining..."
				llm_simplify_metadata.py --info-dir "$output_folder/$name/reads_study_info"
			fi
		fi
		if [ ! -s srr_ids.txt ]; then
			echo -e "\nI haven't been able to find SRR accession ids to download the sequences and I'm exiting, please double check manually..."; exit 1
		fi
		echo -e "\nAll available info downloaded from GEO, please check it out in $output_folder/$name/reads_study_info\n"
		if [ "$design_custom" == "yes" ]; then
			echo "You have requested to manually provide the experimental design instead of the ones shown above. This is the list of samples:"
			cat sample_names.txt
			echo "Please provide a comma-separated list with the conditions for each sample, and if more than one separate the comma-separated lists with '/', no spaces:"
			read -r design_input
			rm $(ls -d $output_folder/$name/reads_study_info/* | grep 'design_possible_')
			IFS='/' read -ra ADDR <<< "$design_input"
			for i in "${!ADDR[@]}"; do
				# For each comma-separated list, split by ',' and echo to file
				IFS=',' read -ra ITEMS <<< "${ADDR[$i]}"
				for item in "${ITEMS[@]}"; do
					echo "$item"
				done > "$output_folder/$name/reads_study_info/design_possible_full_$(($i + 1)).txt"
				cat "$output_folder/$name/reads_study_info/design_possible_full_$(($i + 1)).txt" | sort | uniq > "$output_folder/$name/reads_study_info/design_possible_$(($i + 1)).txt"
			done
		fi

	### Stop and continue with other script if it's a single-cell:
		if [ $(zcat $output_folder/$name/reads_study_info/*_series_matrix.txt.gz | egrep -e 'single nuclei|single cell|single-cell|snRNA|scRNA' | wc -l) -gt 0 ] || [ $(zcat $(find . -name "*_series_matrix.txt.gz") | egrep -i -e 'single nuclei|single cell|single-cell|snRNA|scRNA' | wc -l) -gt 0 ]; then
			echo -e "\n\nDetected this could be a single-cell RNA-seq study... I can try to do stuff automatically (i.e. try and normalize the raw counts or give an estimated bulk expression taking the average), but errors are expected. \nThe script 'R_process_reanalyzer_GSE_single_cell.R' is a template built from the case example GSE118257, and valid to other GEO entries where pheno data and matrix counts are supplementary files clearly named.\nHowever, manual changes are most likely required to work with other studies... These changes should be possible, so please open an issue or go for it if you have the expertise and this one fails!. For example, it's likely that it's just required to point to the directory of the matrix counts, or to manually specify the columns/names of the conditions/cells\n\n"
			# Print the results to review:
			echo -e "This text in the metadata is what made the pipeline to suggest this could be single-cell:"
			zcat $output_folder/$name/reads_study_info/*_series_matrix.txt.gz | egrep -e 'single nuclei|single cell|single-cell|snRNA|scRNA'
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
		if [ $(zcat $output_folder/$name/reads_study_info/*_series_matrix.txt.gz | egrep -i -e 'Expression profiling by array|microarray' | wc -l) -gt 0 ] || [ $(zcat $(find . -name "*_series_matrix.txt.gz") | egrep -i -e 'Expression profiling by array|microarray' | wc -l) -gt 0 ]; then
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
		organism=$(zcat $output_folder/$name/reads_study_info/*_series_matrix.txt.gz | grep "organism" | awk '{$1=""}1' |tr '"' '\n' | sort -u | sed -r '/^\s*$/d')
		echo $organism > $output_folder/$name/reads_study_info/organism.txt
		if [ $(zcat $output_folder/$name/reads_study_info/*_series_matrix.txt.gz | grep "organism" | awk '{$1=""}1' |tr '"' '\n' | sort -u | sed -r '/^\s*$/d' | wc -l) -gt 1 ]; then
			echo -e "\n Please keep in mind that two different organisms are detected. You are likely requesting an analysis combining multiple GSEXXXXX, please make sure they are from the same organism. Another possibility is there are multiple series_matrix within the same GSEXXXX id, and you may have requested to stop and manually clarify. Continuing with organism: "
			organism=$(zcat *_series_matrix.txt.gz | grep "organism" | awk '{$1=""}1' |tr '"' '\n' | sort -u | sed -r '/^\s*$/d' | head -1)
			echo $organism > $output_folder/$name/reads_study_info/organism.txt
			echo -e "$organism\nPlease request on the next run a stop with parameter '-S' and modify manually the file reads_study_info/organism.txt if not required...\n"
		fi
		# An organism provided manually with '-O' always takes precedence over the one detected in
		# the GEO metadata. It is written to organism.txt because every later step re-reads that
		# file, so this is what makes the manual choice effective throughout the whole pipeline.
		if [ ! -z "$organism_argument" ]; then
			organism_from_geo="$organism"
			organism=$(echo $organism_argument | sed 's,_, ,g')
			echo -e "\nOrganism detected in GEO: '$organism_from_geo'. Overriding it with the one provided in '-O': '$organism' (this manual value will be used in all the following steps).\n"
			echo $organism > $output_folder/$name/reads_study_info/organism.txt
		fi
_log_step "Step_1_Download" "end"
		echo -e "\nSTEP 1 DONE. Current time: $(date)\n"
	fi
_log_step "Step_1_Download" "end"
  	echo -e "\n\nSTEP 1: DONE\nCurrent date/time: $(date)\n\n"
	export debug_step="all"
fi


###### STEP 1. Download and process fastq files from the GEO ID provided:
if run_step step1a; then
	rm -rf $seqs_location
	mkdir -p $TMPDIR
	if [[ $input == G* ]]; then
 		echo -e "\n\nSTEP 1: Downloading from the $input id provided...\nCurrent date/time: $(date)\n\n"
		if [ "$stop" == "yes" ]; then
			echo "You have requested a stop to manually provide the SRR ids, or potentially modify other files that may have not been detected properly from GEO, and were not correct, or you just want to adapt some of them. Please double check or manually modify the files reads_study_info/srr_ids.txt, samples_info.txt, sample_names.txt, phenodata_extracted.txt, library_layout_info.txt, organism.txt, design_files, etc. The pipeline is stopped. Please press space to continue or Ctrl + C to exit..."
			read -n1 -s -r -p $'Press space to continue...\n' key
			rm $output_folder/$name/possible_designs_all.txt
			for i in $(ls $output_folder/$name/reads_study_info | grep "full"); do echo $i >> $output_folder/$name/possible_designs_all.txt && cat $output_folder/$name/reads_study_info/$i >> $output_folder/$name/possible_designs_all.txt && echo -e "\n" >> $output_folder/$name/possible_designs_all.txt; done
			cat $output_folder/$name/reads_study_info/phenodata_extracted.txt > $output_folder/$name/phenotypic_data_samples.txt
		fi
		if [ ! -d "$seqs_location" ]; then # I'm now removing the seqs_location at the beginning of this section, in the context of the new system of resuming by -Dm stepx, so this should always be done
			mkdir -p $seqs_location
			echo "Downloading the fastq files from SRR..."
			if [ -z "$input_geo_reads" ]; then
				download_sra_fq.sh $output_folder/$name/reads_study_info/srr_ids.txt $seqs_location $number_parallel $cores $compression_level
				if [ -f "$output_folder/$name/reads_study_info/library_layout_info.txt" ]; then
					cp -f "$output_folder/$name/reads_study_info/library_layout_info.txt" "$output_folder/$name/library_layout_info.txt" 2>/dev/null || true
				fi
	### Rename the fastq files (max length name 140 characters) or handle already downloaded datasets if provided:
				cd $seqs_location
				if [[ "$(cat $output_folder/$name/reads_study_info/library_layout_info.txt)" == "SINGLE" ]]; then
					for i in $(cat $output_folder/$name/reads_study_info/srr_ids.txt); do
						sname=$(awk -F'\t' -v id="$i" '$1 == id { print $2 }' $output_folder/$name/reads_study_info/samples_info.txt)
						if [ -n "$sname" ]; then
							old_file=$(ls | egrep "^${i}[._]" 2>/dev/null | head -1)
							[ -z "$old_file" ] && old_file=$(ls | egrep "^${i}" 2>/dev/null | head -1)
							[ -f "$old_file" ] && mv "$old_file" "${sname}_1.fastq.gz"
						fi
					done
				elif [[ "$(cat $output_folder/$name/reads_study_info/library_layout_info.txt)" == "PAIRED" ]]; then
					for i in $(cat $output_folder/$name/reads_study_info/srr_ids.txt); do
						sname=$(awk -F'\t' -v id="$i" '$1 == id { print $2 }' $output_folder/$name/reads_study_info/samples_info.txt)
						if [ -n "$sname" ]; then
							f1=$(ls | egrep "^${i}_1\." 2>/dev/null | head -1)
							f2=$(ls | egrep "^${i}_2\." 2>/dev/null | head -1)
							[ -z "$f1" ] && f1=$(ls | egrep "^${i}" 2>/dev/null | head -1)
							[ -z "$f2" ] && f2=$(ls | egrep "^${i}" 2>/dev/null | tail -1)
							[ -f "$f1" ] && mv "$f1" "${sname}_1.fastq.gz"
							[ -f "$f2" ] && mv "$f2" "${sname}_2.fastq.gz"
						fi
					done
				fi
				# If layout is SINGLE and files still lack _1.fastq.gz or _2.fastq.gz suffix, add it
				if [[ "$(cat $output_folder/$name/reads_study_info/library_layout_info.txt)" == "SINGLE" ]]; then
					cd $seqs_location
					if [ $(ls *.fastq.gz 2>/dev/null | egrep -c "_[12]\.fastq\.gz$") -eq 0 ] && [ $(ls *.fastq.gz 2>/dev/null | wc -l) -gt 0 ]; then
						echo -e "\nSingle-end reads detected without _1.fastq.gz suffix, adding it...\n"
						for f in *.fastq.gz; do
							mv "$f" "${f%.fastq.gz}_1.fastq.gz"
						done
					fi
				fi
			else
				echo -e "\nSoft linking the already downloaded raw reads from the provided directory: $input_geo_reads\n"
				ln -sf $input_geo_reads/* $seqs_location
			fi
			num_files=$(ls | wc -l); num_samples=$(cat $output_folder/$name/reads_study_info/srr_ids.txt | wc -l)
			if [ "$num_files" -lt "$num_samples" ]; then
				echo -e "\nPlease double check manually, is there some issue with the downloaded raw data? Exiting the script...\n"
				exit 1
			fi
		fi
		echo -e "\nDONE. Current date/time: $(date)"; time1=`date +%s`; echo -e "Elapsed time (secs): $((time1-start))"; echo -e "Elapsed time (hours): $(echo "scale=2; $((time1-start))/3600" | bc -l)\n"

	### Process if any download was not successful or subsampling was required:
		num_gz_files=$(ls -1 $output_folder/$name/raw_reads/*.fastq.gz 2>/dev/null | wc -l)
		num_samples=$(cat $output_folder/$name/reads_study_info/srr_ids.txt | wc -l)

		if [ "$num_gz_files" -eq "$(($num_samples * 2))" ] || [ "$num_gz_files" -eq "$num_samples" ]; then
			echo "Number of fastq.gz files matched expected sample counts."
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
			echo $input | tr ',' '\n' | parallel --halt-on-error 2 --joblog $output_folder/$name/fastq_dl_log_parallel.txt -j $number_parallel --max-args 1 "if [ \$(echo {} | egrep -c 'PRJEB|PRJNA|PRJDB|ERX|DRX|SRX|ERP|DRP|SRP') -eq 1 ]; then fastq-dl --cpus $cores_parallel --silent --force --max-attempts 10 --gzip-level $compression_level --accession {}; fi && 
																																		   if [ \$(echo {} | egrep -c 'ERS|DRS|SRS|SAMD|SAME|SAMN|ERR|DRR|SRR') -eq 1 ]; then fastq-dl --cpus $cores_parallel --silent --force --max-attempts 10 --gzip-level $compression_level --accession {}; fi"
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
			parallel --halt-on-error 2 --verbose -j $cores_reads_to_subsample subsample_reads {} ::: "${arr2[@]}" :::+ "${arr3[@]}" # 10 only because of RAM
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
if run_step step1a_bis; then
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

if run_step step1b; then
	mkdir -p $TMPDIR
	if [[ $input == /* ]]; then
		echo -e "\n\nSTEP 1b: Preparing the raw reads and metadata provided locally...\nCurrent date/time: $(date)\n\n"
		_log_step "Step_1b_Fastp" "start"
  		seqs_location=$output_folder/$name/raw_reads
		rm -rf $seqs_location # I'm now removing the seqs_location at the beginning of this section, in the context of the new system of resuming by -Dm stepx, so this should always be done
		if [ ! -d "$seqs_location" ]; then
			mkdir -p $seqs_location
			if [ $(ls -d $input/* | egrep -c "_R1.fastq.gz$|_R1.fq.gz$|_R2.fastq.gz$|_R2.fq.gz$|_1.fastq.gz$|_1.fq.gz$|_2.fastq.gz$|_2.fq.gz$|_R1_[0-9]+.fastq.gz$|_R1_[0-9]+.fq.gz$|_R2_[0-9]+.fastq.gz$|_R2_[0-9]+.fq.gz$") -eq 0 ]; then
				echo -e "\nPlease make sure that the input files are named _1.fastq.gz, _R1.fastq.gz, _R1_001.fastq.gz, _2.fastq.gz, _R2.fastq.gz, _R2_001.fastq.gz\n"
				exit 1
			fi
			for f in $(ls -d $input/*); do ln -sf $f $seqs_location/$(basename $f | sed 's,fq,fastq,g;s,_R1_[0-9]*.fastq,_1.fastq,g;s,_R2_[0-9]*.fastq,_2.fastq,g;s,_R1.fastq,_1.fastq,g;s,_R2.fastq,_2.fastq,g'); done
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
			parallel --halt-on-error 2 --verbose -j $cores_reads_to_subsample subsample_reads {} ::: "${arr2[@]}" :::+ "${arr3[@]}" # 10 only because of RAM
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
		mkdir -p $output_folder/$name/reads_study_info/
		paste <(ls $seqs_location | egrep '.fq|.fastq' | sed 's,_[12].fastq.gz,,g' | uniq) <(paste -d'_' <(ls $seqs_location | egrep '.fq|.fastq' | sed 's,_[12].fastq.gz,,g' | uniq) <(echo $design_input | sed 's*/*\t*g'| cut -f1 | sed 's*,*\n*g')) <(echo $design_input | sed 's*/*\t*g'| cut -f1 | sed 's*,*\n*g') > $output_folder/$name/reads_study_info/samples_info.txt
		IFS='/' read -ra ADDR <<< "$design_input"
		for i in "${!ADDR[@]}"; do
			# For each comma-separated list, split by ',' and echo to file
			IFS=',' read -ra ITEMS <<< "${ADDR[$i]}"
			for item in "${ITEMS[@]}"; do
				echo "$item"
			done > "$output_folder/$name/reads_study_info/design_possible_full_$(($i + 1)).txt"
			cat "$output_folder/$name/reads_study_info/design_possible_full_$(($i + 1)).txt" | sort | uniq > "$output_folder/$name/reads_study_info/design_possible_$(($i + 1)).txt"
		done
		echo $name > $output_folder/$name/reads_study_info/study_title.txt
		if [ -z "$organism_argument" ]; then
			echo -n "Please input the scientific name of the organism: "
			read -r organism
		else
			organism=$(echo $organism_argument | sed 's,_, ,g')
			echo "Organism used is $organism"
		fi
		echo $organism > $output_folder/$name/reads_study_info/organism.txt
_log_step "Step_1b_Fastp" "end"
		echo -e "\n\nSTEP 1b: DONE\nCurrent date/time: $(date)\n\n"
 	fi
	export debug_step="all"
fi


### STEP 1. Process if required to download from manually provided ids from databases
if run_step step1c; then
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
			echo $input | tr ',' '\n' | parallel --halt-on-error 2 -j $number_parallel --max-args 1 "if [ \$(echo {} | egrep -c 'PRJEB|PRJNA|PRJDB|ERX|DRX|SRX|ERP|DRP|SRP') -eq 1 ]; then fastq-dl --cpus $cores_parallel --silent --force --max-attempts 10 --gzip-level $compression_level --accession {}; fi && 
		 																		   if [ \$(echo {} | egrep -c 'ERS|DRS|SRS|SAMD|SAME|SAMN|ERR|DRR|SRR') -eq 1 ]; then fastq-dl --cpus $cores_parallel --silent --force --max-attempts 10 --gzip-level $compression_level --accession {}; fi"
		fi
_log_step "Step_1_Download" "end"
		echo -e "\n\nSTEP 1: DONE\nCurrent date/time: $(date)\n\n"
 	fi
	export debug_step="all"
fi


### STEP 1. Deal with batch correction... The user has to use certain arguments to manually provide a list or do it interactively:
if run_step step1d; then
	if [ "$batch" == "yes" ]; then
		echo -e "\n\nSTEP 1: Preparing batch effect correction...\nCurrent date/time: $(date)\n\n"
  		if [ -z "$batch_vector" ]; then
			echo -e "This is the content of $seqs_location:\n$(ls -l $seqs_location | awk '{ print $9 }' | tail -n +2)\n"
			echo -n "Based on the list above, please input a comma-separated list for the vector for batch separation (use only numbers, and if these are paired-end, only once per pair of reads):"
			read -r batch_vector
		fi
		echo $batch_vector > $output_folder/$name/reads_study_info/batch_vector.txt
		echo -e "\nThe comma-separated list for the vector for batch separation is $batch_vector\n"
		if [ -z "$batch_biological_covariates" ]; then
			echo -n "Please input a comma-separated list for the biological covariate, and separate by space if multiple biological variables are to be included (use only numbers): "
			read -r batch_biological_covariates
		fi
		echo $batch_biological_covariates > $output_folder/$name/reads_study_info/batch_biological_variables.txt
		echo -e "\nThe comma-separated list for the vector of biological covariable for batch separation is $batch_biological_covariates\n"
	fi
	if [ ! -z "$covariables" ]; then
 		echo $covariables > $output_folder/$name/reads_study_info/covariables.txt
  	fi
fi

### STEP 1. Give info of NCBI's current genome:
Rscript -e "organism <- '${organism}'; tryCatch({ ids <- rentrez::entrez_search(db='assembly', term=paste0(organism, '[orgn]'))\$ids; if(length(ids)==0) stop('No ids'); assemblies <- rentrez::entrez_summary(db='assembly', id=ids[1]); cat(paste(paste0('\n\nNCBI current assembly info: ', date()), assemblies\$assemblyname, assemblies\$assemblyaccession, assemblies\$submissiondate, '\n', sep='\n')) }, error=function(e) cat(paste0('\n\nNo genome information found in NCBI for: ', organism, '\n\n')))"
organism=$(cat $output_folder/$name/reads_study_info/organism.txt | sed 's/ \+/_/g;s/__*/_/g') # Get again the organism in case it has been manually modified... and without spaces...
echo -e "\nYou are using $reference_genome\n"

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

for _gi in "${!genome_group_fastas[@]}"; do
	if [[ "${genome_group_fastas[$_gi]}" == *.gz ]]; then
		_decompressed_gg="$output_folder/$name/indexes/$(basename ${genome_group_fastas[$_gi]%.gz})"
		if [ ! -f "$_decompressed_gg" ]; then
			files_to_decompress+=("${genome_group_fastas[$_gi]}")
			decompressed_outputs+=("$_decompressed_gg")
		else
			echo -e "\nDecompressed genome of group '${genome_group_labels[$_gi]}' already exists: $_decompressed_gg\n"
		fi
		genome_group_fastas[$_gi]="$_decompressed_gg"
	fi
done

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
	parallel --halt-on-error 2 --tmpdir $TMPDIR -j $number_parallel 'echo "Using {1} -> {2}"; gunzip -c "{1}" > "{2}"' ::: "${files_to_decompress[@]}" :::+ "${decompressed_outputs[@]}"
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
	layout_fastp=$(find $output_folder/$name -name library_layout_info.txt 2>/dev/null | head -1 | xargs cat 2>/dev/null | head -n1 | tr -d ' \r\n')

	if [[ "$layout_fastp" == "SINGLE" ]]; then
		ls -d $seqs_location/*.fastq.gz | \
			parallel --halt-on-error 2 --tmpdir $TMPDIR --verbose --joblog $output_folder/$name/fastp_out/fastp_log_parallel.txt -j $number_parallel \
			'fastp --in1 {} --out1 {}_fastp.fastq.gz --dont_overwrite --dont_eval_duplication '$fastp_extra_opts' --thread '$cores_fastp' -z '$compression_level' -h '$output_folder/$name'/fastp_out/{/}_report.html -j '$output_folder/$name'/fastp_out/{/}_report.json &>> '$output_folder/$name'/fastp_out/{/}_fastp_out.log'
	elif [[ "$layout_fastp" == "PAIRED" ]]; then
		ls -d $seqs_location/*.fastq.gz | sed 's,_[12].fastq.gz,,g' | sort | uniq | \
			parallel --halt-on-error 2 --tmpdir $TMPDIR --verbose --joblog $output_folder/$name/fastp_out/fastp_log_parallel.txt -j $number_parallel \
			'fastp --in1 {}_1.fastq.gz --in2 {}_2.fastq.gz --out1 {}_1.fastq.gz_fastp.fastq.gz --out2 {}_2.fastq.gz_fastp.fastq.gz --dont_overwrite --dont_eval_duplication '$fastp_extra_opts' --thread '$cores_fastp' -z '$compression_level' -h '$output_folder/$name'/fastp_out/{/}_report.html -j '$output_folder/$name'/fastp_out/{/}_report.json &>> '$output_folder/$name'/fastp_out/{/}_fastp_out.log'
	fi
fi

# Trimming step can run in addition to the adapter/mode step above
if [ "$fastp_trimming" != "none" ]; then
	echo -e "\n\nSTEP 1: Preprocessing with fastp and trimming...\nCurrent date/time: $(date)\n\n"
	mkdir -p $output_folder/$name/fastp_out
	cd $output_folder/$name/fastp_out
	IFS=', ' read -r -a arrfastp <<< "$fastp_trimming"
	layout_fastp=$(find $output_folder/$name -name library_layout_info.txt 2>/dev/null | head -1 | xargs cat 2>/dev/null | head -n1 | tr -d ' \r\n')

	if [[ "$layout_fastp" == "SINGLE" ]]; then
		ls -d $seqs_location/*.fastq.gz | \
			parallel --halt-on-error 2 --tmpdir $TMPDIR --verbose --joblog $output_folder/$name/fastp_out/fastp_trim_log_parallel.txt -j $number_parallel \
			'fastp --in1 {} --out1 {}_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --trim_front1 '"${arrfastp[0]}"' --trim_tail1 '"${arrfastp[1]}"' --thread '$cores_fastp' -z '$compression_level' -h '$output_folder/$name'/fastp_out/{/}_trim_report.html -j '$output_folder/$name'/fastp_out/{/}_trim_report.json &>> '$output_folder/$name'/fastp_out/{/}_fastp_trim_out.log'
	elif [[ "$layout_fastp" == "PAIRED" ]]; then
		ls -d $seqs_location/*.fastq.gz | sed 's,_[12].fastq.gz,,g' | sort | uniq | \
			parallel --halt-on-error 2 --tmpdir $TMPDIR --verbose --joblog $output_folder/$name/fastp_out/fastp_trim_log_parallel.txt -j $number_parallel \
			'fastp --in1 {}_1.fastq.gz --in2 {}_2.fastq.gz --out1 {}_1.fastq.gz_fastp.fastq.gz --out2 {}_2.fastq.gz_fastp.fastq.gz --dont_overwrite --dont_eval_duplication --trim_front1 '"${arrfastp[0]}"' --trim_tail1 '"${arrfastp[1]}"' --thread '$cores_fastp' -z '$compression_level' -h '$output_folder/$name'/fastp_out/{/}_trim_report.html -j '$output_folder/$name'/fastp_out/{/}_trim_report.json &>> '$output_folder/$name'/fastp_out/{/}_fastp_trim_out.log'
	fi
fi

if [ $(ls -d $seqs_location/* | grep -c _fastp.fastq.gz) -gt 0 ]; then
	for f in $(ls -d $seqs_location/* | grep _fastp.fastq.gz); do mv $f $(echo $f | sed 's,.fastq.gz_fastp.fastq.gz,.fastq.gz,g'); done
	echo "Files in $seqs_location have been successfully processed by fastp!"
fi


### STEP 2. Decontamination if required:
if run_step step2; then
	if [ ! -z "$kraken2_databases" ]; then
  		echo -e "\n\nSTEP 2: Decontamination starting with Kraken2 (k2 daemon mode)...\nCurrent date/time: $(date)\n\n"
_log_step "Step_2_Decontamination" "start"
    		rm -rf $output_folder/$name/raw_reads_k2
		mkdir -p $output_folder/$name/raw_reads_k2
		cd $output_folder/$name/raw_reads_k2

		IFS=',' read -r -a k2_db_array <<< "$kraken2_databases"
		IFS=',' read -r -a k2_conf_array <<< "$kraken2_confidence"
		layout=$(find $output_folder/$name -name library_layout_info.txt 2>/dev/null | head -1 | xargs cat 2>/dev/null | head -n1 | tr -d ' \r\n')

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

		layout_sortmerna=$(find $output_folder/$name -name library_layout_info.txt 2>/dev/null | head -1 | xargs cat 2>/dev/null | head -n1 | tr -d ' \r\n')
		sortmerna_idx=$output_folder/$name/indexes/$(basename $sortmerna_databases)_sortmerna_index/idx
		sortmerna_out=$seqs_location\_sortmerna

		echo -e "\nExecuting sortmerna in parallel and fastqc of the new reads...\n"

		if [[ "$layout_sortmerna" == "PAIRED" ]]; then
			ls $seqs_location | egrep ".fastq$|.fq$|.fastq.gz$|.fq.gz$" \
				| sed 's,_[12].fastq.gz,,g' | sort | uniq \
				| parallel --halt-on-error 2 --tmpdir $TMPDIR --verbose \
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
				| parallel --halt-on-error 2 --tmpdir $TMPDIR --verbose \
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
if run_step step2b; then
	if [ ! -z "$rrna_qc_databases" ]; then
		echo -e "\n\nSTEP 2b: rRNA QC ...\nCurrent date/time: $(date)\n\n"
_log_step "Step_2b_rRNA_QC" "start"

		RRNA_QC_DIR=$output_folder/$name/preliminar_rrna_qc
		RRNA_INDEX_DIR=$output_folder/$name/indexes/rrna_bowtie2_idx
		RRNA_INDEX_NAME="rrna_ref"
		RRNA_MIN_SCORE=$rrna_qc_min_score

		mkdir -p "$RRNA_QC_DIR" "$RRNA_INDEX_DIR"

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
if run_step step3a; then
	echo -e "\n\nSTEP 3a: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_3a_Prepare" "start"
	cd $output_folder/$name/
	rm -rf mi*
	mkdir -p $TMPDIR
	# If the running is resumed in this step, the above has to be done
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/reads_study_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi

	### Prepare the salmon index from the trancripts sequences if required and strandness prediction... (if the first miarma ini does not exist yet, pointing to a previous miarma run)
	if [[ ! -e "$output_folder/$name/miarma0.ini" && ! -e "$output_folder/$name/miarma_genome0.ini" ]]; then
		if [ -z "$strand" ]; then
			echo -e "\nLooking for indexes or indexing transcripts in $transcripts...\n"
			mkdir -p $output_folder/$name/indexes
			salmon_idx=$CURRENT_DIR/indexes/${organism}_salmon_idx
			if [ ! -d "$salmon_idx" ]; then
				if [ ! -d "$output_folder/$name/indexes/${organism}_salmon_idx" ]; then
					salmon index -p $cores -t $transcripts -i $output_folder/$name/indexes/${organism}_salmon_idx --tmpdir $TMPDIR &> $output_folder/$name/indexes/${organism}_salmon_idx.log
				fi
				salmon_idx=$output_folder/$name/indexes/${organism}_salmon_idx
			fi

			mkdir -p $output_folder/$name/strand_prediction/salmon_out
			layout=$(find $output_folder/$name -name library_layout_info.txt 2>/dev/null | xargs cat 2>/dev/null | head -n 1 | tr -d ' \r\n')
			if [[ "$layout" == "SINGLE" ]]; then
				rand_sample=$(ls $seqs_location | shuf | head -1)
				echo -e "\nPredicting strandness with salmon on random sample: $rand_sample\n"
				salmon quant -i $salmon_idx -l A -r $seqs_location/$rand_sample -p $cores -o $output_folder/$name/strand_prediction/salmon_out/ --skipQuant &> $output_folder/$name/strand_prediction/salmon_out/salmon_out.log
			elif [[ "$layout" == "PAIRED" ]]; then
				rand_sample_root=$(ls $seqs_location | sed -E 's/_[12]\.fastq\.gz$//' | grep -v '^\s*$' | sort -u | shuf | head -1)
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
		layout_detected=$(find $output_folder/$name -name library_layout_info.txt 2>/dev/null | head -1 | xargs cat 2>/dev/null | head -n1 | tr -d ' \r\n')
		if [[ "$layout_detected" == "SINGLE" ]]; then
			library_layout=Single
		elif [[ "$layout_detected" == "PAIRED" ]]; then
			library_layout=Paired
		fi
		read_length_for_miarma=$(zcat $seqs_location/$(ls $seqs_location | shuf | head -1) | head -2 | sed -n '2p' | awk '{print length -1}')

		# Final renaming of fastq raw files if SRR present in the filename:
		if [ $(ls $seqs_location | grep -c SRR) -gt 0 ]; then
			for i in $(ls $seqs_location/*); do mv $i $(echo $i | sed 's,_SRR.*_,_,g'); done
		fi

		declare -A gg_index_of_sample=()
		if [ ! -z "$reference_genome_groups" ]; then
			echo -e "\nResolving genome groups (-rG)...\n"
			gg_assignment_file=$output_folder/$name/reads_study_info/genome_group_assignment.tsv
			mkdir -p $output_folder/$name/reads_study_info
			echo -e "sample\tgenome_group\treference_genome" > $gg_assignment_file
			mapfile -t _gg_samples < <(ls $seqs_location | egrep "\.fastq\.gz$|\.fq\.gz$" | sed -E 's,_(R)?[12]\.(fastq|fq)\.gz$,,;s,\.(fastq|fq)\.gz$,,' | sort -u)
			if [ ${#_gg_samples[@]} -eq 0 ]; then
				echo -e "\n\033[1;31mERROR:\033[0m No fastq files found in $seqs_location, so genome groups (-rG) cannot be resolved.\n" >&2; exit 1
			fi
			_gg_unassigned=(); _gg_ambiguous=()
			for _s in "${_gg_samples[@]}"; do
				_hits=()
				for _i in "${!genome_group_labels[@]}"; do
					if echo "$_s" | egrep -q "${genome_group_regexes[$_i]}"; then _hits+=("$_i"); fi
				done
				if [ ${#_hits[@]} -eq 0 ]; then
					_gg_unassigned+=("$_s")
				elif [ ${#_hits[@]} -gt 1 ]; then
					_lab=""; for _h in "${_hits[@]}"; do _lab="$_lab, ${genome_group_labels[$_h]}"; done
					_gg_ambiguous+=("$_s (matches: ${_lab#, })")
				else
					gg_index_of_sample[$_s]=${_hits[0]}
					echo -e "$_s\t${genome_group_labels[${_hits[0]}]}\t${genome_group_fastas[${_hits[0]}]}" >> $gg_assignment_file
				fi
			done
			if [ ${#_gg_unassigned[@]} -gt 0 ] || [ ${#_gg_ambiguous[@]} -gt 0 ]; then
				echo -e "\n\033[1;31mERROR:\033[0m every sample must match exactly one genome group in -rG.\n" >&2
				if [ ${#_gg_unassigned[@]} -gt 0 ]; then
					echo -e "Samples matching NO group (${#_gg_unassigned[@]}):" >&2
					printf '  %s\n' "${_gg_unassigned[@]}" >&2
				fi
				if [ ${#_gg_ambiguous[@]} -gt 0 ]; then
					echo -e "Samples matching MORE THAN ONE group (${#_gg_ambiguous[@]}):" >&2
					printf '  %s\n' "${_gg_ambiguous[@]}" >&2
				fi
				echo -e "\nGroups provided:" >&2
				for _i in "${!genome_group_labels[@]}"; do echo "  ${genome_group_labels[$_i]}: '${genome_group_regexes[$_i]}'" >&2; done
				echo -e "\nPlease adjust the regexes in -rG so the groups partition the samples exactly.\n" >&2
				exit 1
			fi
			echo "Sample to reference genome assignment (also saved in $gg_assignment_file):"
			column -t -s$'\t' $gg_assignment_file 2>/dev/null || cat $gg_assignment_file

			_gg_shared_annot=$(echo "$annotation" | cut -d',' -f1)
			_gg_gtf_extents=$TMPDIR/gg_gtf_extents.txt
			mkdir -p $TMPDIR
			zcat -f "$_gg_shared_annot" | awk -F'\t' '!/^#/ && NF>=9 { if ($5+0 > m[$1]) m[$1]=$5+0 } END { for (k in m) print k"\t"m[k] }' > $_gg_gtf_extents
			if [ ! -s "$_gg_gtf_extents" ]; then
				echo -e "\n\033[1;31mERROR:\033[0m could not read any feature from the annotation $_gg_shared_annot to validate it against the genome groups.\n" >&2; exit 1
			fi
			for _i in "${!genome_group_labels[@]}"; do
				_gg_fa=${genome_group_fastas[$_i]}
				_gg_fai=$output_folder/$name/indexes/$(basename $_gg_fa).gg.fai
				if [ -s "${_gg_fa}.fai" ]; then
					_gg_fai="${_gg_fa}.fai"
				elif [ ! -s "$_gg_fai" ]; then
					if ! { samtools faidx "$_gg_fa" --fai-idx "$_gg_fai" &>/dev/null && [ -s "$_gg_fai" ]; }; then
						if samtools faidx "$_gg_fa" &>/dev/null && [ -s "${_gg_fa}.fai" ]; then
							_gg_fai="${_gg_fa}.fai"
						else
							echo "  Could not index $_gg_fa with samtools, reading contig lengths directly (slower)..."
							awk '/^>/ { if (n != "") print n"\t"l; n=substr($1,2); l=0; next } { l+=length($0) } END { if (n != "") print n"\t"l }' "$_gg_fa" > "$_gg_fai"
						fi
					fi
				fi
				if [ ! -s "$_gg_fai" ]; then
					echo -e "\n\033[1;31mERROR:\033[0m could not determine the contig names/lengths of ${genome_group_labels[$_i]} ($_gg_fa).\n" >&2; exit 1
				fi
				_gg_bad=$(awk -F'\t' 'NR==FNR { len[$1]=$2+0; next }
					{ if (!($1 in len)) { print "  missing contig: "$1; bad++ }
					  else if ($2+0 > len[$1]) { print "  "$1": annotation reaches "$2" but the contig is only "len[$1]" bp"; bad++ }
					} END { if (bad > 10) print "  ... and "bad-10" more" }' "$_gg_fai" "$_gg_gtf_extents" | head -11)
				if [ ! -z "$_gg_bad" ]; then
					echo -e "\n\033[1;31mERROR:\033[0m the annotation $_gg_shared_annot is not valid on the genome of group '${genome_group_labels[$_i]}' ($_gg_fa):\n$_gg_bad\n" >&2
					echo -e "All the genome groups must share one feature space, so the same annotation has to be usable on every genome (same contig names and compatible coordinates). This is the case for variant-personalized genomes built by substituting SNPs, or for a common assembly plus extra sample-specific contigs, but NOT for assemblies with shifted coordinates (indels, different builds). Please provide per-genome annotations lifted to a common gene-ID space, or use one reference for all samples.\n" >&2
					exit 1
				fi
				echo "  Annotation validated against ${genome_group_labels[$_i]} ($(basename $_gg_fa))"
			done

			_gg_sinfo=$output_folder/$name/reads_study_info/samples_info.txt
			if [ -s "$_gg_sinfo" ]; then
				_gg_pairs=$TMPDIR/gg_design_pairs.txt; : > $_gg_pairs
				for _s in "${!gg_index_of_sample[@]}"; do
					_gg_des=$(awk -F'\t' -v s="$_s" '{ for (i=1;i<=NF;i++) if ($i==s) { print $NF; exit } }' $_gg_sinfo)
					[ ! -z "$_gg_des" ] && echo -e "${genome_group_labels[${gg_index_of_sample[$_s]}]}\t$_gg_des" >> $_gg_pairs
				done
				_gg_ngroups=$(cut -f1 $_gg_pairs | sort -u | wc -l)
				_gg_nlevels=$(cut -f2 $_gg_pairs | sort -u | wc -l)
				_gg_shared=$(sort -u $_gg_pairs | cut -f2 | sort | uniq -c | awk '$1>1' | wc -l)
				if [ "$_gg_ngroups" -gt 1 ] && [ "$_gg_nlevels" -gt 1 ] && [ "$_gg_shared" -eq 0 ]; then
					echo -e "\n\033[1;31mWARNING:\033[0m the genome group is perfectly confounded with the experimental design: no condition is present in more than one genome group."
					echo -e "Differential expression computed on the merged matrix would not be able to separate the biological effect from the reference/mapping effect. Consider aligning all the samples to a single common reference for the contrasts you care about, or restrict the interpretation to within-group comparisons."
					echo -e "Genome group vs design:"; sort -u $_gg_pairs | column -t -s$'\t' 2>/dev/null || sort -u $_gg_pairs
					_gg_secs=15
					while [ $_gg_secs -gt 0 ]; do echo -ne "Continuing in $_gg_secs\033[0K\r"; sleep 1; : $((_gg_secs--)); done
					echo ""
				fi
			fi

			### Stage one reads folder per group (symlinks, so no data is duplicated)
			for _i in "${!genome_group_labels[@]}"; do
				rm -rf $output_folder/$name/reads_group$_i; mkdir -p $output_folder/$name/reads_group$_i
			done
			for _f in $(ls $seqs_location | egrep "\.fastq\.gz$|\.fq\.gz$"); do
				_s=$(echo "$_f" | sed -E 's,_(R)?[12]\.(fastq|fq)\.gz$,,;s,\.(fastq|fq)\.gz$,,')
				_i=${gg_index_of_sample[$_s]}
				ln -sf $seqs_location/$_f $output_folder/$name/reads_group$_i/$_f
			done
			for _i in "${!genome_group_labels[@]}"; do
				echo "  Group ${genome_group_labels[$_i]}: $(ls $output_folder/$name/reads_group$_i | wc -l) fastq file(s) staged in reads_group$_i"
			done
		fi

	### Prepare the ini file:
		IFS=', ' read -r -a array <<< "$annotation"
		IFS=', ' read -r -a array2 <<< "$optionsFeatureCounts_seq"
		IFS=', ' read -r -a array3 <<< "$optionsFeatureCounts_feat"
		# In kallisto mode annotation is optional; ensure loop runs at least once
		if [[ "$aligner" == "kallisto" && ${#array[@]} -eq 0 ]]; then
			array=("")
		fi
		build_alignment_units
		for index in "${!unit_ini[@]}"; do
			cd $output_folder/$name
			cp $CURRENT_DIR/external_software/miARma-seq/bakk_miARma1.7.ini ${unit_ini[index]}
			gff=${unit_gtf[index]}
			if [ ! -z "$reference_genome_groups" ]; then fc_opt_index=0; else fc_opt_index=$index; fi
			number_files=$(ls ${unit_reads[index]} | sed 's,_[12].fastq.gz.*,,g' | uniq | wc -l)
			if [ "$number_files" -lt 1 ]; then number_files=1; fi
			if [ $number_files -le $number_parallel ]; then
				cores_parallel=$((cores / number_files))
			else
				cores_parallel=$((cores / number_parallel))
			fi
			if [ $cores_parallel -lt 1 ]; then cores_parallel=1; fi
			sed -i "s,read_length=,read_length=$read_length_for_miarma,g" ${unit_ini[index]}
			sed -i "s,read_dir=,read_dir=${unit_reads[index]},g" ${unit_ini[index]}
			sed -i "s,^threads=,threads=$cores_parallel,g" ${unit_ini[index]}
			sed -i "s,label=,label=$name,g" ${unit_ini[index]}
			sed -i "s,miARmaPath=,miARmaPath=$miarma_path,g" ${unit_ini[index]}
			sed -i "s,output_dir=,output_dir=${unit_out[index]},g" ${unit_ini[index]}
			sed -i "s,stats_file=miARma_stat.log,stats_file=${unit_out[index]}/miARma_stat$index.log,g" ${unit_ini[index]}
			sed -i "s,logfile=miARma_logfile.log,logfile=${unit_out[index]}/miARma_logfile$index.log,g" ${unit_ini[index]}
			sed -i "s,strand=yes,strand=$strand,g" ${unit_ini[index]}
			sed -i "s,fasta=,fasta=${unit_fasta[index]},g" ${unit_ini[index]}
			sed -i "s,gtf=,gtf=$gff,g" ${unit_ini[index]}
			sed -i "s,database=,database=$gff,g" ${unit_ini[index]}
			sed -i "s,seqtype=Paired,seqtype=$library_layout,g" ${unit_ini[index]}
			sed -i "s,organism=mouse,organism=$organism,g" ${unit_ini[index]}
			sed -i "s,indexthreads=,indexthreads=$indexthreads,g" ${unit_ini[index]}
			sed -i "s,parallelnumber=,parallelnumber=$number_parallel,g" ${unit_ini[index]}
			sed -i "s,memorylimit=,memorylimit=$memory_max,g" ${unit_ini[index]}
			if [[ "$aligner" == "star" ]]; then
				if [ -z "${unit_ri[index]}" ]; then
					sed -i "s,indexname=,indexname=${organism}_$(basename ${unit_fasta[index]%.*})_$(basename ${gff%.*})_star_idx,g" ${unit_ini[index]}
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" ${unit_ini[index]}
				else
					sed -i "s,starindex=,starindex=${unit_ri[index]},g" ${unit_ini[index]}
					sed -i "s,indexname=,indexname=${organism}_star_idx,g" ${unit_ini[index]}
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" ${unit_ini[index]}
				fi
			elif [[ "$aligner" == "hisat2" ]]; then
				sed -i "s,aligner=star,aligner=hisat2,g" ${unit_ini[index]}
				if [ -z "${unit_ri[index]}" ]; then
					sed -i "s,indexname=,indexname=${organism}_$(basename ${unit_fasta[index]%.*})_$(basename ${gff%.*})_hisat2_idx,g" ${unit_ini[index]}
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" ${unit_ini[index]}
				else
					sed -i "s,hisat2index=,hisat2index=${unit_ri[index]},g" ${unit_ini[index]}
					sed -i "s,indexname=,indexname=${organism}_hisat2_idx,g" ${unit_ini[index]}
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" ${unit_ini[index]}
				fi
			elif [[ "$aligner" == "kallisto" ]]; then
				sed -i "s,aligner=star,aligner=kallisto,g" ${unit_ini[index]}
				sed -i "s,fasta=${unit_fasta[index]},fasta=$transcripts,g" ${unit_ini[index]}
				if [ -z "${unit_ri[index]}" ]; then
					sed -i "s,indexname=,indexname=${organism}_$(basename ${transcripts%.*})_kallisto_idx,g" ${unit_ini[index]}
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" ${unit_ini[index]}
				else
					sed -i "s,kallistoindex=,kallistoindex=${unit_ri[index]},g" ${unit_ini[index]}
					sed -i "s,indexname=,indexname=${organism}_kallisto_idx,g" ${unit_ini[index]}
					sed -i "s,indexdir=,indexdir=$output_folder/$name/indexes/,g" ${unit_ini[index]}
				fi
			fi
			if [ ! -z "$optionsFeatureCounts_seq" ]; then
				sed -i "s,seqid=gene_name,seqid=${array2[fc_opt_index]},g" ${unit_ini[index]}
			fi
			if [ ! -z "$optionsFeatureCounts_feat" ]; then
				sed -i "s,featuretype=exon,featuretype=${array3[fc_opt_index]},g" ${unit_ini[index]}
			fi

			# ── Validate featureCounts parameters against annotation file ──
			mkdir -p $TMPDIR
			fc_feat_val="${array3[fc_opt_index]:-exon}"
			fc_seq_val="${array2[fc_opt_index]:-gene_name}"
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
				sed -i "s,quality=10,quality=$bam_mapq_threshold,g" ${unit_ini[index]}
				sed -i "s,bam_mapq_threshold=,bam_mapq_threshold=$bam_mapq_threshold,g" ${unit_ini[index]}
			fi
			if [ ! -z "$bam_require_flags" ]; then
				sed -i "s,bam_require_flags=,bam_require_flags=$bam_require_flags,g" ${unit_ini[index]}
			fi
			if [ ! -z "$bam_exclude_flags" ]; then
				sed -i "s,bam_exclude_flags=,bam_exclude_flags=$bam_exclude_flags,g" ${unit_ini[index]}
			fi
			if [ ! -z "$bam_dedup" ]; then
				sed -i "s,bam_dedup=no,bam_dedup=$bam_dedup,g" ${unit_ini[index]}
			fi
			if [ ! -z "$bam_custom_filter" ]; then
				bam_custom_filter_escaped=$(printf '%s' "$bam_custom_filter" | sed 's/[\\&]/\\&/g')
				sed -i "s,bam_custom_filter=,bam_custom_filter=$bam_custom_filter_escaped,g" ${unit_ini[index]}
			fi
			if [ ! -z "$bam_normalization" ]; then
				sed -i "s,bam_normalization=,bam_normalization=$bam_normalization,g" ${unit_ini[index]}
			fi
			if [ "$save_unaligned" == "yes" ]; then
				sed -i "s,save_unaligned=no,save_unaligned=yes,g" ${unit_ini[index]}
			fi
			if [ ! -z "$featureCounts_extra_args" ]; then
				fc_extra_escaped=$(printf '%s' "$featureCounts_extra_args" | sed 's/[\\&]/\\&/g')
				sed -i "s,parameters=-M -O -C -B,parameters=$fc_extra_escaped,g" ${unit_ini[index]}
			fi
			# Extra aligner arguments: each aligner has its own option, so a preset that is only
			# valid for one of them (e.g. hisat2's '--very-sensitive') can never reach another.
			if [[ "$aligner" == "star" ]]; then
				aligner_extra_args_used="$star_extra_args"
			elif [[ "$aligner" == "hisat2" ]]; then
				aligner_extra_args_used="$hisat2_extra_args"
			elif [[ "$aligner" == "kallisto" ]]; then
				aligner_extra_args_used="$kallisto_extra_args"
			else
				aligner_extra_args_used=""
			fi
			if [ ! -z "$aligner_extra_args_used" ]; then
				ae_escaped=$(printf '%s' "$aligner_extra_args_used" | sed 's/[\\&]/\\&/g')
				sed -i "s,${aligner}parameters=,${aligner}parameters=$ae_escaped,g" ${unit_ini[index]}
				echo "Extra args for $aligner: $aligner_extra_args_used"
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
if run_step step3b; then
	echo -e "\n\nSTEP 3b: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_3b_miARma" "start"
	rm -rf $output_folder/$name/miARma_out*
	mkdir -p $TMPDIR
	# If the running is resumed in this step, the above has to be done
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/reads_study_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi
	build_alignment_units

	echo -e "miARma configuration .ini:"
        cat ${unit_ini[0]}
	echo -e "\nPlease double check all the parameters above for miARma, in particular the stranded or the reference genome files and annotation used. Proceeding with miARma execution in..."
	secs=$((1 * 15))
	dir=${unit_out[0]}
	if [ ! -z "$reference_genome_groups" ]; then
		mkdir -p $output_folder/$name/multiqc_out && touch $output_folder/$name/multiqc_out/.rgse_multiqc_deferred
	fi
	while [ $secs -gt 0 ]; do
		echo -ne "$secs\033[0K\r"
		sleep 1
		: $((secs--))
	done
	for index in "${!unit_ini[@]}"; do
		if [ -d "$dir" ] && [ "$(ls -A $dir)" ] && [ "$index" -gt 0 ]; then
			mkdir -p ${unit_out[index]}; cd ${unit_out[index]}
		fi
		cd $output_folder/$name
		if [ ! -z "${unit_label[index]}" ]; then
			echo -e "\n\nAligning and counting genome group '${unit_label[index]}' ($(ls ${unit_reads[index]} | wc -l) fastq file(s)) against $(basename ${unit_fasta[index]})...\n"
		fi
		if [ "$qc_raw_reads" == "no" ]; then
			mkdir -p ${unit_out[index]}/Pre_fastqc_results/_skip_
		fi
		$miarma_path/miARma ${unit_ini[index]}
		miarma_exit=$?
		if [ $miarma_exit -ne 0 ]; then
			echo -e "\n\033[1;31mERROR:\033[0m miARma-seq exited with code $miarma_exit for ${unit_ini[index]}. Please check the log above for details.\n" >&2
			exit $miarma_exit
		fi
	done

	if [ ! -z "$reference_genome_groups" ]; then
		echo -e "\n\nMerging the genome groups into a single quantification...\n"
		merged_dir=$output_folder/$name/miARma_out0
		rc_suffix=$(basename $(find ${unit_out[0]} -maxdepth 1 -type d -name "*_readcount_results" 2>/dev/null | head -1))
		if [ -z "$rc_suffix" ]; then
			echo -e "\n\033[1;31mERROR:\033[0m no *_readcount_results folder found in ${unit_out[0]}. The alignment/counting of the first genome group did not produce counts.\n" >&2; exit 1
		fi
		aligner_results_dir=${aligner}_results
		gg_assignment_file=$output_folder/$name/reads_study_info/genome_group_assignment.tsv
		rm -rf $merged_dir
		mkdir -p $merged_dir/$rc_suffix $merged_dir/$aligner_results_dir $merged_dir/Pre_fastqc_results
		mkdir -p $TMPDIR

		ref_ids=$TMPDIR/gg_feature_ids_ref.txt; cur_ids=$TMPDIR/gg_feature_ids_cur.txt
		common_ids=$TMPDIR/gg_feature_ids_common.txt
		for index in "${!unit_out[@]}"; do
			_one_tab=$(find ${unit_out[index]}/$rc_suffix -maxdepth 1 -name "*.tab" 2>/dev/null | head -1)
			if [ -z "$_one_tab" ]; then
				echo -e "\n\033[1;31mERROR:\033[0m genome group '${unit_label[index]}' produced no count tables in ${unit_out[index]}/$rc_suffix.\n" >&2; exit 1
			fi
			grep -v "^#" "$_one_tab" | tail -n +2 | cut -f1 | sort -u > $cur_ids
			if [ ! -s "$cur_ids" ]; then
				echo -e "\n\033[1;31mERROR:\033[0m could not read any feature id from $_one_tab (or could not write to \$TMPDIR=$TMPDIR). This is a problem reading the counts, not a difference between the genome groups.\n" >&2; exit 1
			fi
			if [ "$index" -eq 0 ]; then
				cp $cur_ids $ref_ids; cp $cur_ids $common_ids
			else
				if ! cmp -s $ref_ids $cur_ids; then
					echo -e "\n\033[1;33mWARNING:\033[0m genome group '${unit_label[index]}' was quantified over a different set of features than '${unit_label[0]}' ($(wc -l < $cur_ids) vs $(wc -l < $ref_ids) features)."
					echo "The merged count matrix will only keep the features present in every group. If this is unexpected, check that the annotation is equally usable on all the group genomes and that featureCounts completed for every sample (see miARma_out_genome*/)."
				fi
				comm -12 $common_ids $cur_ids > $common_ids.tmp && mv $common_ids.tmp $common_ids
			fi
		done
		n_common=$(wc -l < $common_ids)
		if [ "$n_common" -lt 1 ]; then
			echo -e "\n\033[1;31mERROR:\033[0m the genome groups share no feature at all, so their counts cannot be merged into one matrix. Please check that the annotation given in '-a' is the one used for every group. The per-group results are kept in miARma_out_genome*.\n" >&2
			exit 1
		fi
		echo "  Feature space: $n_common feature(s) shared by all the genome groups."

		# Counts are copied (not linked): STEP 4 rewrites them in place when -Fgene is used
		for index in "${!unit_out[@]}"; do
			cp -f ${unit_out[index]}/$rc_suffix/*.tab ${unit_out[index]}/$rc_suffix/*.tab.summary $merged_dir/$rc_suffix/ 2>/dev/null
			# Alignments, coverage tracks and per-sample QC are symlinked (they can be large)
			for _f in $(find ${unit_out[index]}/$aligner_results_dir -maxdepth 1 -mindepth 1 2>/dev/null); do
				_b=$(basename $_f)
				[ "$_b" == "list_multi.txt" ] && continue
				if [ -d "$_f" ]; then
					# bamqc/rnaseqqc folders hold one subfolder per sample, and samples are disjoint
					mkdir -p $merged_dir/$aligner_results_dir/$_b
					for _s in $(find $_f -maxdepth 1 -mindepth 1 2>/dev/null); do ln -sfn $_s $merged_dir/$aligner_results_dir/$_b/$(basename $_s); done
				else
					ln -sf $_f $merged_dir/$aligner_results_dir/$_b
				fi
			done
			for _f in $(find ${unit_out[index]}/Pre_fastqc_results -maxdepth 1 -type f 2>/dev/null); do
				[ "$(basename $_f)" == "list_of_files.txt" ] && continue
				ln -sf $_f $merged_dir/Pre_fastqc_results/$(basename $_f)
			done
			cat ${unit_out[index]}/Pre_fastqc_results/list_of_files.txt >> $merged_dir/Pre_fastqc_results/list_of_files.txt 2>/dev/null
		done
		if [ "$qc_raw_reads" == "no" ]; then mkdir -p $merged_dir/Pre_fastqc_results/_skip_; fi

		n_tabs=$(ls $merged_dir/$rc_suffix/*.tab 2>/dev/null | wc -l)
		n_samples=$(cut -f1 $gg_assignment_file 2>/dev/null | tail -n +2 | sort -u | wc -l)
		if [ "$n_tabs" -ne "$n_samples" ]; then
			echo -e "\n\033[1;31mERROR:\033[0m merged $n_tabs count table(s) but $n_samples sample(s) were assigned to genome groups. Some samples were not quantified; please check the per-group logs in miARma_out_genome*.\n" >&2
			exit 1
		fi
		echo "  Merged $n_tabs per-sample count table(s) from ${#unit_out[@]} genome group(s) into $merged_dir/$rc_suffix"

		# MultiQC was deferred above so that it covers every group; build it once now
		rm -rf $output_folder/$name/multiqc_out; mkdir -p $output_folder/$name/multiqc_out
		echo "  Building MultiQC over all the genome groups..."
		if [ -n "$RGSE_DO_MULTIQC" ] && command -v multiqc_ai.py >/dev/null 2>&1; then
			multiqc_ai.py --analysis-dir $output_folder/$name --out-dir $output_folder/$name/multiqc_out ${MULTIQC_AI_ARGS} > $output_folder/$name/multiqc_out/multiqc.log 2>&1 || {
				echo "AI MultiQC failed - using standard MultiQC" >> $output_folder/$name/multiqc_out/multiqc.log
				multiqc -f $output_folder/$name --ignore 'miARma_stat*' -n multiqc_report -o $output_folder/$name/multiqc_out -p -q >> $output_folder/$name/multiqc_out/multiqc.log 2>&1; }
		else
			multiqc -f $output_folder/$name --ignore 'miARma_stat*' -n multiqc_report -o $output_folder/$name/multiqc_out -p -q > $output_folder/$name/multiqc_out/multiqc.log 2>&1
		fi
	fi

	### Reformat the logs by parallel...
	for f in $(find $output_folder/$name -name "*_log_parallel.txt"); do awk -F"\t" 'NR==1; NR > 1{OFS="\t"; $3=strftime("%Y-%m-%d %H:%M:%S", $3); print $0}' $f > tmp && mv tmp $f; done

 	### Clean genome index cache?
  	if [ "$aligner" == "star" ] && [ "$aligner_index_cache" == "no" ]; then
		for genome_dir_loaded in $(find $output_folder/$name/ -name "star_log_parallel.txt" | xargs cat 2>/dev/null | grep "genomeDir" | sed 's,.*genomeDir ,,g;s, .*,,g' | sort | uniq); do
			STAR --runThreadN $cores --genomeDir $genome_dir_loaded --genomeLoad Remove --outFileNamePrefix genomeloading.tmp
		done
	fi

	### Display alignment percentages per sample on screen and save summary table
	if [ -f "$CURRENT_DIR/scripts/show_alignment_stats.py" ]; then
		python3 $CURRENT_DIR/scripts/show_alignment_stats.py \
			--dir "$output_folder/$name" \
			--aligner "$aligner" \
			--save "$output_folder/$name/reads_study_info/alignment_summary.tsv" 2>/dev/null || true
	fi

	echo -e "\nmiARma-seq DONE. Current date/time: $(date)"; time1=`date +%s`; echo -e "Elapsed time (secs): $((time1-start))"; echo -e "Elapsed time (hours): $(echo "scale=2; $((time1-start))/3600" | bc -l)\n"
	export debug_step="all"
_log_step "Step_3b_miARma" "end"
	echo -e "\n\nSTEP 3b: DONE\nCurrent date/time: $(date)\n\n"
fi


### STEP 4. Process output of miARma. Get figures, final counts, standard DGE, violin plots...
if run_step step4; then
	# If the running is resumed in this step, this variables has to be created because they would not exist
	rm -rf $output_folder/$name/final_results_* # So it's redone when resuming
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/reads_study_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi
	if [ -z "${!array[@]}" ]; then
		IFS=', ' read -r -a array <<< "$annotation"
		# In kallisto mode annotation is optional; ensure loop runs at least once
		if [[ "$aligner" == "kallisto" && ${#array[@]} -eq 0 ]]; then
			array=("")
		fi
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
		# Pre-extract GO terms into the DGE folder so that R_process_reanalyzer_GSE.R
		# can merge them into the final expression tables (otherwise they'd only be
		# extracted in Step 6, too late for the merged tables).
		_pre_go_annot=""
		if [ ! -z "$non_reference_funct_enrichm" ]; then
			_pre_go_annot="$non_reference_funct_enrichm"
		elif [ ! -z "$annotation_file" ] && [ -f "$annotation_file" ] && [ $(zcat -f "$annotation_file" 2>/dev/null | head -500 | egrep -c '(^|[^[:alnum:]])GO:[0-9]{7}|[Oo]ntology_term|go_(process|function|component)') -gt 0 ]; then
			_pre_go_annot="$annotation_file"
		fi
		if [ ! -z "$_pre_go_annot" ]; then
			_pre_go_outdir=$output_folder/$name/final_results_reanalysis$index/DGE
			mkdir -p "$_pre_go_outdir"
			_pre_go_outfile="$_pre_go_outdir/$(basename $_pre_go_annot).automatically_extracted_GO_terms.txt"
			if [ ! -f "$_pre_go_outfile" ]; then
				echo "Pre-extracting GO and functional terms from $_pre_go_annot for merged tables..."
				python3 $CURRENT_DIR/scripts/parse_functional_annotation.py \
					-i "$_pre_go_annot" \
					${annotation_file:+-r "$annotation_file"} \
					-o "$_pre_go_outfile" \
					-k "$_pre_kegg_outfile"
				go_pre_lines=$(grep "GO:" "$_pre_go_outfile" 2>/dev/null | wc -l)
				echo "  Pre-extracted $go_pre_lines gene-GO associations for merged tables."
				kegg_pre_lines=$(grep -v "^source_id" "$_pre_kegg_outfile" 2>/dev/null | wc -l || echo 0)
				if [ "$kegg_pre_lines" -gt 0 ]; then
					echo "  Pre-extracted $kegg_pre_lines gene-KEGG associations for merged tables."
				else
					echo "  No KEGG terms found in annotation."
					rm -f "$_pre_kegg_outfile"
				fi
			fi
		fi
  			echo -e "R_process_reanalyzer_GSE.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index $genes ${array2[index]} $organism $target $differential_expr_soft $batch_format $covariables $covariables_format $deconvolution $differential_expr_comparisons $perform_differential_analyses $perform_volcano_venn $pattern_to_remove $annotation_file $fc_seq_key $fc_feat_type $sc_count_matrix $sc_phenotype $bulk_expression_matrix\n\n" > $output_folder/$name/R_process_reanalyzer.log
    		R_process_reanalyzer_GSE.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index $genes ${array2[index]} $organism $target $differential_expr_soft $batch_format $covariables $covariables_format $deconvolution $differential_expr_comparisons $perform_differential_analyses $perform_volcano_venn $pattern_to_remove $annotation_file $fc_seq_key $fc_feat_type $sc_count_matrix $sc_phenotype $bulk_expression_matrix | tee -a $output_folder/$name/R_process_reanalyzer.log
		if [ ! -d "$output_folder/$name/final_results_reanalysis$index" ]; then
			echo -e "\nWARNING: R_process_reanalyzer_GSE.R did not create $output_folder/$name/final_results_reanalysis$index."
			echo "Check $output_folder/$name/R_process_reanalyzer.log for details. Skipping remaining post-processing for this index."
			continue
		fi
    		echo 'R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter" "edgeR_object" "edgeR_object_norm" $pattern_to_remove $annotation_file $fc_feat_type $split_sections' > $output_folder/$name/R_qc_figs.log
			# Determine if per-section PDF splitting is needed for AI interleaving
			_ai_val="${ai_insights:-all}"
			if [ -n "$LLM_ENDPOINT" ] && [ "$_ai_val" != "no" ] && [ "$_ai_val" != "enrichment" ] && [ "$_ai_val" != "dge" ] && [ "$_ai_val" != "counts" ] && command -v qc_pdf_ai.py >/dev/null 2>&1; then
				split_sections="yes"
			else
				split_sections="no"
			fi
			R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter" "edgeR_object" "edgeR_object_norm" $pattern_to_remove $annotation_file $fc_feat_type $split_sections | tee -a $output_folder/$name/R_qc_figs.log
		if [[ -e "$output_folder/$name/final_results_reanalysis$index/counts_adjusted.txt" ]]; then
			echo -e "\n\nRemember that batch effect correction/covariables have been only provided to Combat-Seq/limma for visualization purposes, to include covariables in the DGE model after checking the visualization the argument -C will be used\n\n\nQC_PDF adjusted counts\n\nRemember that you have requested batch effect correction/count adjustment, so you have to mind the figures in this QC_PDF from ComBat-seq/limma counts...\n"
			echo -e 'R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter_adjusted" "edgeR_object_adjusted" "edgeR_object_norm_adjusted" $pattern_to_remove $annotation_file $fc_feat_type $split_sections' > $output_folder/$name/R_qc_figs_adjusted.log
			R_qc_figs.R $output_folder/$name $output_folder/$name/miARma_out$index $output_folder/$name/final_results_reanalysis$index "edgeR_object_prefilter_adjusted" "edgeR_object_adjusted" "edgeR_object_norm_adjusted" $pattern_to_remove $annotation_file $fc_feat_type $split_sections | tee -a $output_folder/$name/R_qc_figs_adjusted.log
		fi
	done
	### Generate SummarizedExperiment for the exploreLocalDE app if requested
	if [[ "$exploreDE_se" == "yes" ]]; then
		echo -e "\nGenerating SummarizedExperiment for exploreLocalDE...\n"
		for index in "${!array[@]}"; do
			final_dir=$output_folder/$name/final_results_reanalysis$index
			if [ -f "$final_dir/Raw_counts_genes.txt" ] && [ -f "$final_dir/TPM_counts_genes.txt" ] && [ -f "$final_dir/DGE/list_comp.txt" ]; then
				export ANNOTATION_FILE="${array[index]}"
				Rscript $CURRENT_DIR/scripts/prepare_SE.R \
					"$final_dir/Raw_counts_genes.txt" \
					"$final_dir/TPM_counts_genes.txt" \
					"$output_folder/$name/reads_study_info/samples_info.txt" \
					"$final_dir/DGE/list_comp.txt" \
					"$final_dir/DGE" \
					"^DGE_analysis_comp[0-9].txt$" \
					"$name" \
					"$organism" \
					"${_pre_go_outfile:-$non_reference_funct_enrichm}" 2>&1 | tee -a "$final_dir/DGE/deResults_prepare_SE.log"
			else
				echo "Skipping exploreLocalDE SE generation for index $index: required files not found"
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
if run_step step4b; then
	_log_step "Step_4b_Splicing" "start"
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/reads_study_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
	fi
	if [ -z "${!array[@]}" ]; then
		IFS=', ' read -r -a array <<< "$annotation"
		# In kallisto mode annotation is optional; ensure loop runs at least once
		if [[ "$aligner" == "kallisto" && ${#array[@]} -eq 0 ]]; then
			array=("")
		fi
	fi
	if [[ "$splicing_option" == "saser" ]]; then
		echo -e "\n\nSTEP 4b: saseR splicing analysis...\nCurrent date/time: $(date)\n\n"
		for index in "${!array[@]}"; do
			bam_dir=$output_folder/$name/miARma_out$index/${aligner}_results
			saser_out=$output_folder/$name/final_results_reanalysis$index/saseR_splicing
			mkdir -p $saser_out
			library_layout=$(find $output_folder/$name -name library_layout_info.txt 2>/dev/null | head -1 | xargs cat 2>/dev/null | head -n1 | tr -d ' \r\n')
			samples_info=$output_folder/$name/reads_study_info/samples_info.txt
			design_file=$(ls $output_folder/$name/reads_study_info/design_possible_full_1.txt 2>/dev/null || echo "none")
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
			samples_info=$output_folder/$name/reads_study_info/samples_info.txt
			design_file=$(ls $output_folder/$name/reads_study_info/design_possible_full_1.txt 2>/dev/null || echo "none")
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
		echo -e "\n\nSTEP 4b (IsoformSwitchAnalyzeR): DONE\nCurrent date/time: $(date)\n\n"
	fi
	_log_step "Step_4b_Splicing" "end"
	export debug_step="all"
fi


### STEP 5. Time course analyses if required
if run_step step5; then
	for index in "${!array[@]}"; do
		if [[ "$time_course" == "yes" ]]; then
			echo -e "\n\nSTEP 5: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_5_QC_Figs" "start"
   			echo -e "\nPerforming time course analyses."
			R_process_time_course.R $output_folder/$name/final_results_reanalysis$index/ DGE_analysis_comp1.qs2 edgeR_object_norm $minstd $mestimate
_log_step "Step_5_QC_Figs" "end"
   			echo -e "\n\nSTEP 5: DONE\nCurrent date/time: $(date)\n\n"
		fi
	done
	export debug_step="all"	
fi


### STEP 6. Functional enrichment analyses: clusterProfiler, autoGO, Panther, network analyses...
if run_step step6; then
	# If the running is resumed in this step, this variables has to be created because they would not exist
 	# The same may happen in other steps, this is WIP
	if [ -z "$organism" ]; then
		organism=$(cat $output_folder/$name/reads_study_info/organism.txt | sed 's, ,_,g;s,_+,_,g')
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
		# In kallisto mode annotation is optional; ensure loop runs at least once
		if [[ "$aligner" == "kallisto" && ${#array[@]} -eq 0 ]]; then
			array=("")
		fi
	fi
	for index in "${!array[@]}"; do
		if [ ! -d "$output_folder/$name/final_results_reanalysis$index/DGE/" ]; then
			echo -e "\nWARNING: DGE directory not found at $output_folder/$name/final_results_reanalysis$index/DGE/. Skipping enrichment for index $index."
			continue
		fi
		cd $output_folder/$name/final_results_reanalysis$index/DGE/
		rm -rf $(find . -type d \( -name "*_autoGO" -o -name "*_clusterProfiler" -o -name "*_panther" -o -name "*funct_enr*" \)) $(find . -type f \( -name "*_autoGO" -o -name "*_clusterProfiler" -o -name "*_panther" -o -name "*funct_enr*" \)) # So it's redone if resuming
		if [ -z "$annotation_file" ]; then
			annotation_file=${array[index]}
		fi

		# Check if any DGE comparison file contains significant DEGs (FDR < 0.05)
		deg_count=0
		dge_files=(DGE_analysis_comp*.txt)
		if [ -e "${dge_files[0]}" ]; then
			deg_count=$(awk -F'\t' '
				FNR==1 {
					fdr_col=0;
					for(i=1;i<=NF;i++) {
						if($i=="FDR" || $i=="adj.P.Val" || $i=="padj") fdr_col=i
					}
					next
				}
				fdr_col > 0 && $fdr_col != "" && $fdr_col+0 < 0.05 { count++ }
				END { print count+0 }
			' DGE_analysis_comp*.txt 2>/dev/null || echo 0)
		fi

		if [ "$deg_count" -eq 0 ]; then
			echo -e "\n=========================================================================="
			echo "NOTICE: 0 significant DEGs (FDR < 0.05) found across all comparisons."
			echo "Skipping Step 6 functional enrichment & network analyses for index $index."
			echo "==========================================================================\n"
			echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
			_log_step "Step_6_Enrichment" "start"
			_log_step "Step_6_Enrichment" "end"
			continue
		fi

		# Network analyses (WGCNA is organism-agnostic; STRINGdb supports any organism with a valid taxon ID)
		if [[ $network_analyses == "yes" ]]; then
			mkdir -p network_analyses && rm -rf network_analyses/* && cd network_analyses
			echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_6_Enrichment" "start"
    			echo -e "\nPerforming network analyses (WGCNA mode: $wgcna_mode)...\n"
			R_functional_network_analyses.R $output_folder/$name/final_results_reanalysis$index/DGE/ $output_folder/$name/final_results_reanalysis$index/Raw_counts_genes.txt "^DGE_analysis_comp[0-9]+.txt$" $taxonid $wgcna_mode $organism $output_folder/$name/reads_study_info/samples_info.txt &> network_analyses_funct_enrichment.log
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
				ls | egrep "^DGE_analysis_comp[0-9]+.txt$" | parallel --halt-on-error 2 --joblog R_clusterProfiler_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_clusterProfiler_analyses_parallel.R $PWD $organism "1" $clusterProfiler_method $clusterProfiler_full $aPEAR_execution '^{}$' $clusterProfiler_universe $clusterProfiler_minGSSize $clusterProfiler_maxGSSize &> clusterProfiler_{}_funct_enrichment.log"
				echo -e "\nPerforming autoGO and Panther execution... this may take long if many genes or comparisons...\n"
				ls | egrep "^DGE_analysis_comp[0-9]+.txt$" | parallel --halt-on-error 2 --joblog R_autoGO_panther_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism "1" $databases_function {} $panther_method $auto_panther_log &> autoGO_panther_{}_funct_enrichment.log"
				if [[ "$time_course" == "yes" ]]; then
					cd $output_folder/$name/final_results_reanalysis$index/time_course_analyses
					ls | egrep "^DGE_limma_timecourse.*.txt$" | parallel --halt-on-error 2 --joblog R_clusterProfiler_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_clusterProfiler_analyses_parallel.R $PWD $organism "1" $clusterProfiler_method $clusterProfiler_full $aPEAR_execution '^{}$' $clusterProfiler_universe $clusterProfiler_minGSSize $clusterProfiler_maxGSSize &> clusterProfiler_{}_funct_enrichment.log"
					ls | egrep "^DGE_limma_timecourse.*.txt$" | parallel --halt-on-error 2 --joblog R_autoGO_panther_analyses_parallel_log_parallel.txt -j $cores --max-args 1 "R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism "1" $databases_function {} $panther_method $auto_panther_log &> autoGO_panther_{}_funct_enrichment.log"
				fi
			else
				echo -e "\n\nSTEP 6: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_6_Enrichment" "start"
   				echo "Organism '$organism' is not human/mouse, so functional enrichment support is limited: reanalyzerGSE will attempt GO/KEGG over-representation using terms extracted from the provided annotation (dedicated OrgDb-based analyses are human/mouse only)."
				# Determine which annotation file to use for functional enrichment
				annot_enrichm=""
				if [ ! -z "$non_reference_funct_enrichm" ]; then
					echo "Using provided non-reference functional enrichment file: $non_reference_funct_enrichm"
					annot_enrichm="$non_reference_funct_enrichm"
				elif [ $(zcat -f "$annotation_file" 2>/dev/null | head -500 | egrep -c "GO:|Ontology|tology_term|go_(process|function|component)") -gt 0 ]; then
					echo "Using main annotation file for functional enrichment: $annotation_file"
					annot_enrichm="$annotation_file"
				fi

				enrichment_results_found="no"
				if [ ! -z "$annot_enrichm" ]; then
					cd $output_folder/$name/final_results_reanalysis$index/DGE/
					echo "Applying clusterProfiler enrichr with GO terms extracted from the provided annotation..."

					# GO terms were already pre-extracted in Step 4 for merged tables.
					# Just reference the existing file.
					annotation_go=$output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annot_enrichm).automatically_extracted_GO_terms.txt
					go_data_lines=$(grep "GO:" "$annotation_go" 2>/dev/null | wc -l)

					# KEGG terms were already pre-extracted in Step 4 for merged tables.
					# Just reference the existing file.
					annotation_kegg=$output_folder/$name/final_results_reanalysis$index/DGE/$(basename $annot_enrichm).automatically_extracted_KEGG_terms.txt
					if [ ! -f "$annotation_kegg" ]; then
						echo "No KEGG terms found in the annotation."
						annotation_kegg=""
					else
						kegg_data_lines=$(wc -l < "$annotation_kegg" 2>/dev/null || echo 0)
						echo "Found $kegg_data_lines gene-KEGG associations (pre-extracted)."
					fi

					if [ "$go_data_lines" -lt 2 ] && [ -z "$annotation_kegg" ]; then
						echo "WARNING: Neither GO nor KEGG term extraction produced valid entries. The input file may have an unexpected format. Skipping enrichment."
					else
						if [ "$go_data_lines" -ge 2 ]; then
							echo "Extracted $go_data_lines gene-GO associations. Running enrichr..."
							# Add header if not already present (GAF extraction adds it directly)
							if ! head -1 "$annotation_go" | grep -q "^source_id"; then
								sed -i '1s/^/source_id\tComputed_GO_Process_IDs\n/' "$annotation_go"
							fi
						else
							echo "No GO terms found, skipping GO enrichment."
							annotation_go=""
						fi
						R_clusterProfiler_enrichr.R "$annotation_go" $output_folder/$name/final_results_reanalysis$index/RPKM_counts_genes.txt $output_folder/$name/final_results_reanalysis$index/DGE "^DGE_analysis_comp[0-9]+.txt$" "$annotation_kegg" &> clusterProfiler_enrichr_funct_enrichment.log
						echo "enrichr finished (log: clusterProfiler_enrichr_funct_enrichment.log); status is summarised below."
					fi
				else
					echo "For $organism and the annotation $annotation_file, no GO or functional information found. Consider providing a GAF, GMT, or GO-annotated GFF/GTF via 'non_reference_funct_enrichm'"
				fi
			fi

			cd $output_folder/$name/final_results_reanalysis$index/DGE/
			# A log counts as failed only on a real R fatal-error marker (the old
			# broad "Err" match flagged benign lines). Show the actual error lines so
			# the cause (e.g. a missing package) is visible instead of buried.
			enrich_err_re="Execution halted|Error in |Error: |no package called"
			error_logs=$(grep -lE "$enrich_err_re" *_funct_enrichment.log 2>/dev/null | sort -u)
			if [ -n "$error_logs" ]; then
			    echo -e "\nSome enrichment scripts reported errors:"
			    for lg in $error_logs; do
			        echo "  * $lg"
			        grep -hE "$enrich_err_re" "$lg" 2>/dev/null | sed 's/^/      /' | head -4
			    done

			    # Automatic retry, once. Only the per-comparison parallel scripts
			    # (clusterProfiler_<comp>.txt..., autoGO_panther_<comp>.txt...) can be
			    # re-run this way; single-shot scripts (e.g. the enrichr path) cannot.
			    echo -e "\nRetrying failed per-comparison enrichment scripts once..."
			    failed_cp=$(grep -lE "$enrich_err_re" clusterProfiler_*_funct_enrichment.log 2>/dev/null | sed 's/clusterProfiler_//g;s/_funct_enrichment.log//g' | sort | uniq)
			    for ff in $failed_cp; do
			        if [ -f "$ff" ]; then
			            echo "  Retrying clusterProfiler for $ff ..."
			            rm -rf $(echo $ff | sed 's/.txt$//')_funct_enrich_clusterProfiler
			            R_clusterProfiler_analyses_parallel.R $PWD $organism "1" $clusterProfiler_method $clusterProfiler_full $aPEAR_execution "^${ff}$" $clusterProfiler_universe $clusterProfiler_minGSSize $clusterProfiler_maxGSSize &> clusterProfiler_${ff}_funct_enrichment_retry.log
			        fi
			    done
			    failed_ago=$(grep -lE "$enrich_err_re" autoGO_panther_*_funct_enrichment.log 2>/dev/null | sed 's/autoGO_panther_//g;s/_funct_enrichment.log//g' | sort | uniq)
			    for ff in $failed_ago; do
			        if [ -f "$ff" ]; then
			            echo "  Retrying autoGO+Panther for $ff ..."
			            R_autoGO_panther_analyses_parallel.R $output_folder/$name/final_results_reanalysis$index $organism "1" $databases_function $ff $panther_method $auto_panther_log &> autoGO_panther_${ff}_funct_enrichment_retry.log
			        fi
			    done

			    if ls *_funct_enrichment_retry.log >/dev/null 2>&1; then
			        retry_errors=$(grep -lE "$enrich_err_re" *_funct_enrichment_retry.log 2>/dev/null | sort -u)
			        if [ -n "$retry_errors" ]; then
			            echo "  Still failing after retry: $(echo $retry_errors | tr '\n' ' ')"
			        else
			            echo "  Retry succeeded — no errors on re-run."
			        fi
			    else
			        echo "  Nothing could be retried automatically (the failing step is not a per-comparison script); see the errors above."
			    fi
			fi

			# Success is decided solely by whether result files were produced, not by
			# the presence/absence of log errors. Add the up/down gene counts to them.
			cd $output_folder/$name/final_results_reanalysis$index/
			files_to_process=$(find . \( -name "*.txt" -o -name "*.tsv" -o -name "*.csv" \) | grep funct | grep -v '_err.txt\|_aPEAR\|_similarity')
			if [ -n "$files_to_process" ]; then
				enrichment_results_found="yes"
				cd $output_folder/$name/final_results_reanalysis$index/DGE/
				echo -e "\nFunctional enrichment completed: $(echo $files_to_process | wc -w) result file(s) produced. Formatting..."
				echo $files_to_process | parallel --halt-on-error 2 --joblog R_enrich_format_analyses_parallel_log_parallel.txt -j $cores "file={}; R_enrich_format.R \"\$file\" \$(echo \"\$file\" | sed 's,DGE/.*,DGE/,g')\$(echo \"\$file\" | sed 's,.*DGE_analysis_comp,DGE_analysis_comp,g;s,_pval.*,,g;s,_fdr.*,,g;s,_funct.*,,g;s,_cluster.*,,g' | sed 's,.txt,,g').txt $organism $rev_thr" &> $PWD/enrichment_format.log
			else
				enrichment_results_found="no"
				echo -e "\nFunctional enrichment produced no result files, so the HTML report will not be rendered."
				if [ -n "$error_logs" ]; then
				    echo "Reason: the enrichment step(s) above failed — see the error lines. Fix the cause and re-run this step with '-Dm step6'."
				fi
			fi
		fi

		# --- Optional AI insight boxes for the reports (opt-in: needs -llm_endpoint) ---
		# Reuses the LLM_* env vars exported for MultiQC AI. llm_insight.py no-ops
		# without an endpoint anyway; we also gate here to skip the work when off.
		# Insight files (*.ai_insight.md) are picked up by the enrichment .Rmd and by
		# the Sphinx report; generation failures are non-fatal (logged, boxes omitted).
		# Each ai_insights OPTION (design / counts / dge / enrichment) is handled as
		# its own pass. llm_insight.py exits 42 on a >5-min LLM timeout; on that we
		# stop the pass and fill EVERY box for that option with a 'timeout' placeholder
		# (item: on timeout the whole option is skipped and its boxes say so).
		: "${ai_do_design:=0}" "${ai_do_counts:=0}" "${ai_do_dge:=0}" "${ai_do_enrichment:=0}" "${ai_do_qualimap:=0}" "${ai_do_qc_pdf:=0}"
		if [ -n "$LLM_ENDPOINT" ] && [ $((ai_do_design + ai_do_counts + ai_do_dge + ai_do_enrichment + ai_do_qualimap + ai_do_qc_pdf)) -gt 0 ]; then
			ai_fdir="$output_folder/$name/final_results_reanalysis$index"
			ai_rsi="$output_folder/$name/reads_study_info"
			ai_dge_dir="$ai_fdir/DGE"
			if [ -d "$ai_dge_dir" ]; then
				ai_log="$ai_dge_dir/ai_insights.log"; : > "$ai_log"
				echo -e "\nGenerating AI report insights (ai_insights=$ai_insights, model=${llm_model:-<unset>})..."
				# Write the 'timeout, please try again' placeholder box for one output.
				ai_write_timeout() {   # $1 = output .md path
					printf '%s\n\n%s\n' "**AI summary** — the LLM did not respond in time." \
						"ai_insights timeout. Please try again" > "$1"
				}
				# Base DE comparison tables (comp1, comp2, ...) present in the DGE dir.
				ai_comps=()
				for ai_dgef in "$ai_dge_dir"/DGE_analysis_comp*.txt; do
					[ -e "$ai_dgef" ] || continue
					ai_b=$(basename "$ai_dgef" .txt)
					[[ "$ai_b" =~ ^DGE_analysis_comp[0-9]+$ ]] && ai_comps+=("$ai_b")
				done

				# --- design: organism + comparisons + covariables + layout/strand +
				#     sample metadata (the info shown in report sections 1-3). ---
				ai_sinfo="$ai_rsi/samples_info.txt"
				if [ "$ai_do_design" = 1 ] && [ -f "$ai_sinfo" ]; then
					ai_tmp=$(mktemp)
					{ echo "# organism";     cat "$ai_rsi/organism.txt" 2>/dev/null
					  echo; echo "# comparisons"; cat "$ai_dge_dir/list_comp.txt" 2>/dev/null
					  echo; echo "# covariables / potential batch effects"
					  cat "$ai_rsi/covariables.txt" "$ai_rsi/batch_vector.txt" "$ai_rsi/batch_biological_variable.txt" 2>/dev/null
					  echo; echo "# library layout";  cat "$output_folder/$name/library_layout_info.txt" 2>/dev/null
					  echo; echo "# strandedness";    cat "$output_folder/$name/strand_info.txt" 2>/dev/null
					  echo; echo "# samples (name, condition, batch)"; cat "$ai_sinfo"; } > "$ai_tmp"
					llm_insight.py --input "$ai_tmp" --task design --title "$name" \
						--out "$ai_dge_dir/study_design.ai_insight.md" >> "$ai_log" 2>&1
					[ $? -eq 42 ] && ai_write_timeout "$ai_dge_dir/study_design.ai_insight.md"
					rm -f "$ai_tmp"
				fi

				# --- counts: per-sample expression landscape. The full log2(TPM+1)
				#     table is NEVER sent; we send a compact per-sample summary (value
				#     ranges/quartiles + Low/Medium/High category counts) plus the
				#     sample->condition metadata, so the LLM can comment on it. ---
				ai_counts_tbl="$ai_fdir/TPM_counts_genes_log2_1_categ.txt"
				if [ "$ai_do_counts" = 1 ] && [ -f "$ai_counts_tbl" ]; then
					ai_tmp=$(mktemp)
					python3 - "$ai_counts_tbl" "$ai_sinfo" > "$ai_tmp" 2>/dev/null <<'PYEOF'
import sys
tbl = sys.argv[1]; sinfo = sys.argv[2] if len(sys.argv) > 2 else None
with open(tbl, encoding="utf-8", errors="replace") as fh:
    header = fh.readline().rstrip("\n").split("\t"); ncol = len(header)
    num_cols = [i for i, h in enumerate(header)
                if i > 0 and h and not h.endswith("_categ") and not h.endswith("_categ_2")]
    cat_of = {header[i]: i for i, h in enumerate(header) if h.endswith("_categ") and not h.endswith("_categ_2")}
    vals = {i: [] for i in num_cols}
    cats = {i: {} for i in cat_of.values()}
    ngenes = 0
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < ncol:
            continue
        ngenes += 1
        for i in num_cols:
            try: vals[i].append(float(p[i]))
            except ValueError: pass
        for i in cats:
            v = p[i]; cats[i][v] = cats[i].get(v, 0) + 1
def q(s, pr):
    if not s: return float("nan")
    k = (len(s) - 1) * pr; f = int(k); c = min(f + 1, len(s) - 1)
    return s[f] + (s[c] - s[f]) * (k - f)
print("Per-sample expression summary (values are log2(TPM+1)); genes: %d" % ngenes)
for i in num_cols:
    s = sorted(vals[i]); name = header[i]
    ci = cat_of.get(name + "_categ")
    cc = cats.get(ci, {}) if ci is not None else {}
    cs = ", ".join("%s=%d" % (k, cc[k]) for k in ("Low", "Medium", "High") if k in cc)
    if s:
        print("- %s: min=%.2f q25=%.2f median=%.2f q75=%.2f max=%.2f%s"
              % (name, s[0], q(s, .25), q(s, .5), q(s, .75), s[-1],
                 ("; categories: " + cs) if cs else ""))
if sinfo:
    try:
        with open(sinfo, encoding="utf-8", errors="replace") as sh:
            print("\n# Sample metadata (name, condition, ...):\n" + sh.read().strip())
    except Exception:
        pass
PYEOF
					if [ -s "$ai_tmp" ]; then
						llm_insight.py --input "$ai_tmp" --task counts --title "$name" \
							--out "$ai_dge_dir/counts.ai_insight.md" >> "$ai_log" 2>&1
						[ $? -eq 42 ] && ai_write_timeout "$ai_dge_dir/counts.ai_insight.md"
					fi
					rm -f "$ai_tmp"
				fi

				# --- dge: one DEG summary per comparison (from the base DE table). ---
				if [ "$ai_do_dge" = 1 ]; then
					ai_to=0
					for ai_c in "${ai_comps[@]}"; do
						ai_out="$ai_dge_dir/${ai_c}.ai_insight.md"
						if [ "$ai_to" = 1 ]; then ai_write_timeout "$ai_out"; continue; fi
						llm_insight.py --input "$ai_dge_dir/${ai_c}.txt" --task dge --title "$ai_c" \
							--out "$ai_out" >> "$ai_log" 2>&1
						if [ $? -eq 42 ]; then ai_to=1; ai_write_timeout "$ai_out"; fi
					done
				fi

				# --- enrichment: one narrative per comparison, map-reduced (one
				#     sequential LLM call per table, then a synthesis call) over the
				#     PRIMARY over-representation result tables for that comparison.
				#     The tables live in the per-comparison SUBDIRS (clusterProfiler /
				#     panther / autoGO), not at the DGE top level, so recurse into them
				#     (the old top-level glob found nothing -> no insight was made).
				#     We deliberately send only the headline results and skip the noisy
				#     extras: groupGO/level breakdowns (no p-values), fdr_01 duplicates,
				#     panther expr-background / binom / GO_SLIM variants. Kept:
				#       - autoGO GO_*.tsv (what the enrichment report shows)
				#       - clusterProfiler GO/KEGG/REACT over-representation at fdr_05
				#       - panther whole-background Fisher tables (GO BP/MF/CC, PC, PATHWAY, REACTOME) ---
				if [ "$ai_do_enrichment" = 1 ]; then
					ai_to=0
					for ai_c in "${ai_comps[@]}"; do
						ai_out="$ai_dge_dir/${ai_c}.enrichment_ai_insight.md"
						if [ "$ai_to" = 1 ]; then ai_write_timeout "$ai_out"; continue; fi
						ai_rfs=()
						while IFS= read -r ai_rf; do
							[ -n "$ai_rf" ] && ai_rfs+=("$ai_rf")
						done < <(find "$ai_dge_dir" -type f \( -name '*.tsv' -o -name '*.txt' \) \
							\( -path "*${ai_c}_funct_enrich_clusterProfiler/*" -o -path "*${ai_c}_*funct_enrichment_panther/*" \) 2>/dev/null \
							| grep -E '/_autoGO/GO_.*\.tsv$|/GO_overrepresentation_test.*_fdr_fdr_05\.txt$|/KEGG_enrich_fdr_fdr_05\.txt$|/REACT_fdr_fdr_05\.txt$|/enriched_fisher_whole_back_.*_panther\.txt$' \
							| grep -vE 'GO_SLIM|_err\.txt$' | sort)
						if [ ${#ai_rfs[@]} -gt 0 ]; then
							llm_insight.py --input "${ai_rfs[@]}" --task enrichment --title "$ai_c" --max-rows 0 \
								--out "$ai_out" \
								--rel-dir "$ai_dge_dir" \
								--pertable-out "$ai_dge_dir/${ai_c}.enrichment_ai_pertable.tsv" >> "$ai_log" 2>&1
							if [ $? -eq 42 ]; then ai_to=1; ai_write_timeout "$ai_out"; fi
						fi
					done
				fi
				echo "AI insights written under $ai_dge_dir (*.ai_insight.md); log: $ai_log"
			fi
		fi

		# Render functional enrichment HTML report (self-contained), only if results exist
		if [[ "$functional_enrichment_analyses" != "no" ]] && [[ "$enrichment_results_found" == "yes" ]] && [ -d "$output_folder/$name/final_results_reanalysis$index/DGE" ]; then
			echo -e "Rendering functional enrichment HTML report..."
			Rscript $CURRENT_DIR/scripts/render_enrichment_report.R \
				"$output_folder/$name/final_results_reanalysis$index/DGE" \
				"$name" \
				"$organism" &> "$output_folder/$name/final_results_reanalysis$index/DGE/enrichment_report_render.log"
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
if run_step step7; then
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
				cut -f1 "$file" | parallel --halt-on-error 2 -j $((cores*3)) "gene={}; foldchange=\$(grep -i \"\$gene\" \"$file\" | cut -f3 | sed -n 's/\(.*[.,][0-9]\{2\}\).*/\1/p'); \
															grep -i \"=\$gene\" \"$annotation_file\" | head -1 | awk -v id=\"\$gene\" -v fc=\"\$foldchange\" '{ print \$1\"\\t\"\$4\"\\t\"\$5\"\\t\"id\"_\"fc\"\\t.\t\"\$7 }' >> \"$file.bed\""
			done
		fi

		### Bibliographic context for the top DEGs:
		lit_dge_dir="$output_folder/$name/final_results_reanalysis$index/DGE"
		lit_out="$lit_dge_dir/literature"
		lit_log="$lit_dge_dir/lit_gather.log"
		: "${ai_do_literature:=0}"
		if [ "$literature" != "no" ] && [ -d "$lit_dge_dir" ] && command -v lit_gather.py >/dev/null 2>&1; then
			echo -e "\nGathering bibliographic context for the top DEGs (organism: $organism)..."
			lit_gather.py --dge-dir "$lit_dge_dir" --organism "$organism" \
				--out-dir "$lit_out" \
				--cache-dir "$output_folder/$name/reads_study_info/literature_cache" \
				$literature_args > "$lit_log" 2>&1
			lit_rc=$?
			case $lit_rc in
				0) echo "Bibliographic context saved under $lit_out" ;;
				3) echo "No usable bibliographic context: no significant genes, none resolved at NCBI, or no literature linked. See $lit_log" ;;
				4) echo "NCBI E-utilities unreachable; bibliographic context skipped. See $lit_log" ;;
				*) echo "Literature gathering failed (exit $lit_rc); see $lit_log" ;;
			esac

			if [ "$lit_rc" = 0 ] && [ -n "$LLM_ENDPOINT" ] && [ "$ai_do_literature" = 1 ]; then
				lit_write_timeout() {
					printf '%s\n\n%s\n' "**AI summary** — the LLM did not respond in time." \
						"ai_insights timeout. Please try again" > "$1"
				}
				lit_comps=$(python3 -c "import json,sys; print(' '.join(json.load(open(sys.argv[1]))['usable_comparisons']))" \
					"$lit_out/status.json" 2>/dev/null)
				if [ -n "$lit_comps" ]; then
					echo "Generating AI literature synthesis (model=${llm_model:-<unset>})..."
					lit_to=0
					for lit_c in $lit_comps; do
						lit_box="$lit_dge_dir/${lit_c}.literature_ai_insight.md"
						if [ "$lit_to" = 1 ]; then lit_write_timeout "$lit_box"; continue; fi
						llm_insight.py --input "$lit_out/$lit_c/for_llm.txt" --task literature \
							--title "$lit_c" --max-rows 0 \
							--pmid-whitelist "$lit_out/$lit_c/pmids.txt" \
							--out "$lit_box" >> "$lit_log" 2>&1
						if [ $? -eq 42 ]; then lit_to=1; lit_write_timeout "$lit_box"; fi
					done
				fi
			fi
		elif [ "$literature" = "no" ] && [ "$ai_do_literature" = 1 ] && [ -n "$LLM_ENDPOINT" ] && [ "$index" = 0 ]; then
			echo -e "\nNOTE: 'literature' is among the requested ai_insights, but literature gathering is off ('-lit no'), so there is nothing to summarise."
		fi
	done
	export debug_step="all"
_log_step "Step_7_Annotation" "end"
	echo -e "\n\nSTEP 7: DONE\nCurrent date/time: $(date)\n\n"
fi


###### STEP 8. Sum up results in a sphinx report
if run_step step8; then
	_log_step "Step_8_Report" "start"

	# AI Insights for Qualimap and QC PDFs (only when LLM endpoint is provided)
	ai_insights_val="${ai_insights:-all}"
	if [ -n "$LLM_ENDPOINT" ] && [ "$ai_insights_val" != "no" ] && command -v llm_insight.py >/dev/null 2>&1; then
		ai_do_qualimap=1
		ai_do_qc_pdf=1
		case "$ai_insights_val" in
			enrichment|dge|counts) ai_do_qualimap=0; ai_do_qc_pdf=0 ;;
		esac

		for index in "${!array[@]}"; do
			ai_fdir="$output_folder/$name/final_results_reanalysis$index"
			ai_dge_dir="$ai_fdir/DGE"
			ai_log="$ai_dge_dir/ai_insights.log"
			mkdir -p "$ai_dge_dir"

			if [ "$ai_do_qualimap" = 1 ]; then
				echo "Annotating Qualimap HTML reports with AI insights..."
				qualimap_ai.py --analysis-dir "$output_folder/$name" >> "$ai_log" 2>&1 || true
			fi
			if [ "$ai_do_qc_pdf" = 1 ]; then
				echo "Annotating QC figures PDF with AI insights..."
				_label=$(basename "$output_folder/$name")
				# Process each edgeR variant's sections directory independently
				for edger_suffix in "norm" "norm_adjusted"; do
					sections_dir="$ai_fdir/QC_and_others/sections_${edger_suffix}"
					_qc_pdf="$ai_fdir/QC_and_others/${_label}_${edger_suffix}_QC.pdf"
					if [ -d "$sections_dir" ]; then
						# Per-section mode: R created per-section PDFs, Python assembles + interleaves
						qc_pdf_ai.py --tables-dir "$ai_fdir/QC_and_others/tables" --pdf "$_qc_pdf" --sections-dir "$sections_dir" >> "$ai_log" 2>&1 || true
					elif [ -f "$_qc_pdf" ]; then
						# Fallback: monolithic PDF exists, append AI slides at end
						qc_pdf_ai.py --tables-dir "$ai_fdir/QC_and_others/tables" --pdf "$_qc_pdf" >> "$ai_log" 2>&1 || true
					fi
				done
			fi
		done
	fi

	# Render a preliminary Gantt chart with completed steps (1-7) so it exists when sphinx_report.sh builds
	if [ -f "$STEP_TIMES_FILE" ]; then
		for qc_dir in $(find "$output_folder/$name" -maxdepth 2 -type d -name "QC_and_others" 2>/dev/null); do
			gantt_out="$qc_dir/pipeline_gantt_chart.pdf"
			echo -e "\nRendering preliminary pipeline Gantt chart for the Sphinx report..."
			Rscript $CURRENT_DIR/scripts/R_gantt_chart.R "$STEP_TIMES_FILE" "$gantt_out" 2>&1 || \
				echo "WARNING: Gantt chart rendering failed."
		done
	fi

	# Convert tables to xlsx before sphinx so download links work in the report
	if [ "$convert_tables_excel" == "yes" ]; then
		echo -e "Converting tables to xlsx..."
		R_convert_tables.R $output_folder/$name/ $cores "log_parallel|jquery|bamqc|rnaseqqc|samtools|strand" > $output_folder/$name/R_convert_tables.log 2>&1
	fi
	sphinx_report.sh $output_folder/$name $name
_log_step "Step_8_Report" "end"
 	echo -e "\n\nSTEP 8: Final report DONE\nCurrent date/time: $(date)\n\n"
fi



###### STEP 9. Tidy up, prepare for storage if final results have been created and the number of aligned files is equal to the numbers of samples, rename folders, convert tables to xlsx if required... etc
# Compress the folders
if run_step step9; then
	echo -e "\n\nSTEP 9: Starting...\nCurrent date/time: $(date)\n\n"
_log_step "Step_9_Cleanup" "start"
	echo -e "\n\nTidying up, removing empty folders, temp files, compressing...\n\n"

	# Remove decompressed reference files from the indexes subfolder
	if [ -d "$output_folder/$name/indexes" ]; then
		for _decomp_file in $(find $output_folder/$name/indexes -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.gtf' -o -name '*.gff' \) 2>/dev/null); do
			rm -f "$_decomp_file"
		done
	fi


	# Remove intermediate DEG-list text files (DGE_analysis_comp*_fdr_05.txt, logpos/logneg),
	# but spare everything inside the enrichment result folders: those tables are named after
	# the gene set they came from (e.g. GO_overrepresentation_test_BP_fdr_fdr_05_logpos.txt)
	# and match the same globs. Without the exclusion the broad "*_fdr_05.txt" glob deleted the
	# fdr_05 enrichment tables while keeping fdr_01, so the archived report could no longer be
	# re-rendered with both thresholds.
	cd $output_folder/$name/ && find . -type f \( -name "*_fdr_05.txt" -o -name "*_logneg.txt" -o -name "*_logpos.txt" \) \
		! -path "*_funct_enrich_clusterProfiler/*" ! -path "*_funct_enrichment_panther/*" ! -path "*/enrichr_clusterProfiler/*" \
		-exec rm -f {} +
	# Note: xlsx conversion now happens in STEP 8 (before sphinx report), not here

	for index in "${!array[@]}"; do
	 	if [ ! -d "$output_folder/$name/final_results_reanalysis$index/DGE/" ]; then
			continue
		fi
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
		build_alignment_units
		any_tidied=0
		for index in "${!unit_out[@]}"; do
			num_raw_files=$(cat ${unit_out[index]}/Pre_fastqc_results/list_of_files.txt 2>/dev/null | grep -c "fastq.gz")
			[ "$num_raw_files" -lt 1 ] && continue
			if [ ! -z "$reference_genome_groups" ]; then _fr_glob="final_results_reanalysis0*"; else _fr_glob="final_results_reanalysis$index*"; fi
			if [ -z "$(find $output_folder/$name -maxdepth 1 -type d -name "$_fr_glob" 2>/dev/null | head -1)" ]; then
				echo "Skipping the clean up of ${unit_out[index]}: no results folder matching '$_fr_glob' (the analysis did not complete)"
				continue
			fi
			n_bams=$(ls ${unit_out[index]}/$aligner\_results 2>/dev/null | egrep -c ".bam$")
			if [[ "$n_bams" -eq "$num_raw_files" || "$n_bams" -eq $((num_raw_files / 2)) ]]; then
				echo -e "\nTidying up${unit_label[index]:+ (genome group ${unit_label[index]})}...\n"
				cd ${unit_out[index]}/$aligner\_results
				echo "For the sake of efficiente storage: samtools view -@ cores -T ref_genome -C -o xxx.bam.cram xxx.bam && rm xx.bam" >> conversion_bam_to_cram.txt
				echo "Reference genome used for this CRAM conversion: ${unit_fasta[index]}" >> conversion_bam_to_cram.txt
				find . -type f -name "*.bam" | parallel --halt-on-error 2 --verbose -j $number_parallel --max-args 1 samtools view -T ${unit_fasta[index]} -C -@ $((cores / number_parallel)) -o {}.cram {}
				rm -rf $(ls | egrep ".bam$")
				any_tidied=1
			fi
		done
		if [ "$any_tidied" -eq 1 ] && [ -d "$seqs_location" ]; then
			cd $seqs_location
			echo "After execution, raw reads have been removed for the sake of efficient storage. These were... " > readme
			ls -lh | grep -v "^total" | grep -v " readme$" >> readme
			ls | grep -v readme | xargs -r rm -rf
			rm -rf $TMPDIR
		fi
		if [ ! -z "$reference_genome_groups" ] && [ -d "$output_folder/$name/miARma_out0/${aligner}_results" ]; then
			find $output_folder/$name/miARma_out0/${aligner}_results -xtype l -delete 2>/dev/null
			for index in "${!unit_out[@]}"; do
				for _f in $(find ${unit_out[index]}/${aligner}_results -maxdepth 1 -type f -name "*.cram" 2>/dev/null); do
					ln -sf $_f $output_folder/$name/miARma_out0/${aligner}_results/$(basename $_f)
				done
			done
		fi
		# Remove the indexes folder since data has been converted to CRAM
		if [ -d "$output_folder/$name/indexes" ]; then
			echo "Removing indexes folder: $output_folder/$name/indexes"
			rm -rf "$output_folder/$name/indexes"
		fi
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
				if [ ! -z "$reference_genome_groups" ]; then
					{
						echo "LIMITATION: this analysis used SEVERAL reference genomes (-rG), but igvShinyApp.R can only be configured with one."
						echo "The deployed app points at: ${reference_genome} (the -r genome), with coverage tracks from ${bw_dir}/"
						echo "Only the samples aligned to that same genome will be displayed correctly. To browse the others, edit igvShinyApp.R"
						echo "and replace the reference genome and the bigwig folder with the matching pair from the table below."
						echo ""
						echo "group<TAB>reference_genome<TAB>bigwig_folder" | sed 's,<TAB>,\t,g'
						for _gi in "${!genome_group_labels[@]}"; do
							echo -e "${genome_group_labels[$_gi]}\t${genome_group_fastas[$_gi]}\t$output_folder/$name/miARma_out_genome$_gi/${aligner}_results"
						done
						echo ""
						echo "The per-sample assignment is in reads_study_info/genome_group_assignment.tsv"
					} > "$final_dir_igv/igvShinyApp_README.txt"
					echo "Genome groups in use: wrote $final_dir_igv/igvShinyApp_README.txt documenting the single-reference limitation of the IGV app"
				fi
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
	# Also update the final Gantt chart in the built Sphinx HTML report directory if it exists
	if [ -d "$output_folder/$name/sphinx_report/html" ]; then
		echo "Updating final pipeline Gantt chart in the Sphinx report..."
		Rscript $CURRENT_DIR/scripts/R_gantt_chart.R "$STEP_TIMES_FILE" "$output_folder/$name/sphinx_report/html/pipeline_gantt_chart.pdf" 2>&1 || \
			echo "WARNING: Sphinx report Gantt chart update failed."
	fi
fi

# Shrink the built report tree, once nothing else writes into it
if [ "${tidy_report_files:-yes}" != "no" ] && [ -d "$output_folder/$name/sphinx_report/html" ]; then
	$CURRENT_DIR/scripts/prune_sphinx_html.py "$output_folder/$name/sphinx_report/html" \
		--extra-html "$output_folder/$name/final_report.html" || \
		echo "WARNING: report tidying failed. The report itself is unaffected."
fi

if [[ -n "$end_step" && "$end_step" != "none" && "$end_step" != "all" ]]; then
	echo -e "\n\n\nPIPELINE FORCE-ENDED after '$end_step' as requested (-Es/-end_step); any later steps were skipped. Best wishes\n\n\n"
else
	echo -e "\n\n\nALL STEPS DONE! Best wishes\n\n\n"
fi
