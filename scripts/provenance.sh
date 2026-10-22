#!/bin/bash

rgse_option_keys() {
	sed -n 's/^[[:space:]]*-[^)]*)[[:space:]]*\([a-zA-Z_][a-zA-Z0-9_]*\)=\${arguments\[index\]}.*/\1/p' "$CURRENT_DIR_SCRIPTS/parse_options.sh" | sort -u
	printf 'debug_step\norganism\n'
}

rgse_pipeline_commit() {
	local root c
	root=$(dirname "$CURRENT_DIR_SCRIPTS")
	if [ -s "$root/COMMIT" ]; then cat "$root/COMMIT"; return 0; fi
	if command -v git >/dev/null 2>&1 && c=$(git -C "$root" rev-parse HEAD 2>/dev/null); then
		[ -n "$(git -C "$root" status --porcelain --untracked-files=no 2>/dev/null)" ] && c="$c (with uncommitted changes)"
		echo "$c"; return 0
	fi
	echo "unknown"
}

rgse_container() { echo "${APPTAINER_CONTAINER:-${SINGULARITY_CONTAINER:-none}}"; }

record_run_parameters() {
	local dir="$output_folder/$name/provenance" stamp key val
	mkdir -p "$dir" 2>/dev/null || return 0
	stamp=$(date +%Y%m%d_%H%M%S)
	{
		echo "# reanalyzerGSE v${REANALYZER_VERSION}, commit $(rgse_pipeline_commit)"
		echo "# Options in effect for the run started on $(date '+%Y-%m-%d %H:%M:%S'): command line, YAML and defaults merged."
		echo "# Rerun with 'reanalyzerGSE.sh -options <this file>' (set debug_step to 'all' for a full run). llm_endpoint and llm_api_key are never recorded."
		for key in $(rgse_option_keys); do
			case "$key" in
				llm_endpoint|llm_api_key) echo "$key: \"\""; continue ;;
			esac
			val="${!key}"
			val=${val//\\/\\\\}
			val=${val//\"/\\\"}
			echo "$key: \"$val\""
		done
	} > "$dir/run_parameters_$stamp.yaml"
	cp "$dir/run_parameters_$stamp.yaml" "$dir/run_parameters_latest.yaml"
	[ -s "$dir/run_history.tsv" ] || printf 'Started\tVersion\tCommit\tFrom_step\tEnd_step\tHost\tContainer\tParameters\n' > "$dir/run_history.tsv"
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$REANALYZER_VERSION" "$(rgse_pipeline_commit)" \
		"$debug_step" "${end_step:-none}" "$(hostname)" "$(rgse_container)" "run_parameters_$stamp.yaml" >> "$dir/run_history.tsv"
}

write_provenance() {
	local dir="$output_folder/$name/provenance" entry label cmd v f list
	mkdir -p "$dir" 2>/dev/null || return 0
	{
		echo "reanalyzerGSE v${REANALYZER_VERSION}"
		echo "Commit: $(rgse_pipeline_commit)"
		[ -s "$(dirname "$CURRENT_DIR_SCRIPTS")/BUILD_REF" ] && echo "Image built from: $(cat "$(dirname "$CURRENT_DIR_SCRIPTS")/BUILD_REF")"
		echo "Container: $(rgse_container)"
		echo "Host: $(hostname), $(uname -sr)"
		echo "Report built: $(date '+%Y-%m-%d %H:%M:%S')"
		echo ""
		for entry in "hisat2|hisat2 --version" "STAR|STAR --version" "kallisto|kallisto version" "salmon|salmon --version" \
			"samtools|samtools --version" "featureCounts|featureCounts -v" "fastp|fastp --version" "FastQC|fastqc --version" \
			"MultiQC|multiqc --version" "Bowtie2|bowtie2 --version" "Kraken2|kraken2 --version" "SortMeRNA|sortmerna --version" \
			"deepTools|bamCoverage --version" "fastq-dl|fastq-dl --version" "R|Rscript --version" "Python|python3 --version" \
			"Perl|perl -v" "Sphinx|sphinx-build --version" "GNU parallel|parallel --version"; do
			label=${entry%%|*}; cmd=${entry#*|}
			command -v "${cmd%% *}" >/dev/null 2>&1 || continue
			v=$(timeout 60 bash -c "$cmd" 2>&1 | grep -v '^[[:space:]]*$' | head -1)
			printf '%-14s %s\n' "$label" "$v"
		done
	} > "$dir/tool_versions.txt"
	Rscript -e 'p <- c("edgeR","limma","DESeq2","ashr","sva","clusterProfiler","enrichplot","DOSE","ReactomePA","msigdbr","pathview","autoGO","rbioapi","WGCNA","STRINGdb","saseR","IsoformSwitchAnalyzeR","BisqueRNA","CDSeq","Mfuzz","GO.db","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","qs2"); cat(R.version.string, "\n"); if (requireNamespace("BiocManager", quietly = TRUE)) cat("Bioconductor", as.character(BiocManager::version()), "\n"); for (x in p) cat(sprintf("%-22s %s\n", x, tryCatch(as.character(packageVersion(x)), error = function(e) "not installed")))' > "$dir/R_packages.txt" 2>/dev/null
	python3 -m pip list --format=freeze > "$dir/python_packages.txt" 2>/dev/null
	{
		printf 'md5\tbytes\tfile\n'
		list="$reference_genome,$annotation,$transcripts,$non_reference_funct_enrichm,$sortmerna_databases,$rrna_qc_databases,$proteome_faa,$ortho_human_proteome"
		for f in "${genome_group_fastas[@]}"; do list="$list,$f"; done
		for f in $(echo "$list" | tr ',' '\n' | sed '/^[[:space:]]*$/d' | sort -u); do
			[ -f "$f" ] || continue
			printf '%s\t%s\t%s\n' "$(md5sum "$f" | cut -d' ' -f1)" "$(stat -L -c %s "$f")" "$f"
		done
		for f in $(echo "$kraken2_databases" | tr ',' '\n' | sed '/^[[:space:]]*$/d'); do
			[ -f "$f/hash.k2d" ] && printf 'not computed (Kraken2 database, size and date only)\t%s\t%s (%s)\n' "$(stat -L -c %s "$f/hash.k2d")" "$f/hash.k2d" "$(date -r "$f/hash.k2d" '+%Y-%m-%d')"
		done
	} > "$dir/reference_checksums.tsv"
}
