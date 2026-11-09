#!/bin/bash
dry_errors=0; dry_warnings=0
dry_err() { echo -e "  \033[1;31mERROR:\033[0m $1"; dry_errors=$((dry_errors + 1)); }
dry_warn() { echo -e "  \033[1;33mWARNING:\033[0m $1"; dry_warnings=$((dry_warnings + 1)); }
dry_ok() { echo "  ok: $1"; }
dry_file() {
	local label=$1 f=$2
	[ -z "$f" ] && return 0
	if [ -e "$f" ]; then dry_ok "$label: $f"; else dry_err "$label does not exist: $f"; fi
}
dry_tool() {
	local t
	for t in "$@"; do
		if command -v "$t" >/dev/null 2>&1; then dry_ok "$t ($(command -v "$t"))"; else dry_err "$t is not in PATH"; fi
	done
}

echo -e "\n\n========================================================================"
echo "DRY RUN: the configuration is checked and nothing is written to $output_folder/$name"
echo -e "========================================================================\n"

echo "Input ($input):"
dry_mode=""
if [[ $input == G* ]]; then
	dry_mode="GEO"
	dry_ok "GEO series $input: the metadata, the reads and the design come from GEO/SRA at run time"
	[ -n "$organism_argument" ] && dry_ok "organism forced to '$organism_argument'"
elif [[ $input == /* ]]; then
	dry_mode="local"
	if [ ! -d "$input" ]; then
		dry_err "the input folder does not exist: $input"
	else
		mapfile -t dry_fastqs < <(ls "$input" | egrep "_R1.fastq.gz$|_R1.fq.gz$|_R2.fastq.gz$|_R2.fq.gz$|_1.fastq.gz$|_1.fq.gz$|_2.fastq.gz$|_2.fq.gz$|_R1_[0-9]+.fastq.gz$|_R1_[0-9]+.fq.gz$|_R2_[0-9]+.fastq.gz$|_R2_[0-9]+.fq.gz$")
		[ -n "$input_filter_regex" ] && mapfile -t dry_fastqs < <(printf '%s\n' "${dry_fastqs[@]}" | egrep "$input_filter_regex")
		[ -n "$input_filter_regex_exclude" ] && mapfile -t dry_fastqs < <(printf '%s\n' "${dry_fastqs[@]}" | egrep -v "$input_filter_regex_exclude")
		dry_samples=$(printf '%s\n' "${dry_fastqs[@]}" | sed -E 's/\.(fastq|fq)\.gz$//; s/_R?[12](_[0-9]+)?$//' | sort -u | grep .)
		n_dry=$(echo "$dry_samples" | grep -c .)
		if [ "$n_dry" -eq 0 ]; then
			dry_err "no fastq.gz files named *_1/_2, *_R1/_R2 or *_R1_001/_R2_001 in $input"
		else
			dry_ok "$n_dry sample(s) in ${#dry_fastqs[@]} fastq file(s):"
			echo "$dry_samples" | sed 's/^/      /'
			if [ -z "$design_custom_local" ]; then
				dry_warn "no design_custom_local: the run will stop and ask for the design interactively"
			elif [ "$design_custom_local" == "auto" ]; then
				if echo "$dry_samples" | grep -Evq '_[Rr][Ee][Pp][0-9]'; then
					dry_err "design 'auto' needs a _Rep<N> suffix in every sample name"
				else
					dry_ok "design 'auto': $(echo "$dry_samples" | sed -E 's,_[Rr][Ee][Pp][0-9].*$,,' | paste -sd,)"
				fi
			elif [ "$design_custom_local" == "llm" ]; then
				if [ -n "$LLM_ENDPOINT" ] && [ -n "$LLM_MODEL" ]; then dry_ok "design 'llm' (model $LLM_MODEL)"; else dry_err "design 'llm' needs llm_endpoint and llm_model"; fi
			else
				IFS='/' read -ra dry_designs <<< "$design_custom_local"
				for k in "${!dry_designs[@]}"; do
					n_cond=$(echo "${dry_designs[$k]}" | awk -F',' '{print NF}')
					if [ "$n_cond" -ne "$n_dry" ]; then
						dry_err "design $((k + 1)) has $n_cond condition(s) for $n_dry sample(s)"
					else
						dry_ok "design $((k + 1)) (sample, condition):"
						paste <(echo "$dry_samples") <(echo "${dry_designs[$k]}" | tr ',' '\n') | sed 's/^/      /'
					fi
				done
			fi
		fi
	fi
	if [ -z "$organism_argument" ]; then dry_warn "no organism: the run will ask for it interactively"; else dry_ok "organism '$organism_argument'"; fi
else
	dry_mode="accessions"
	echo "  querying ENA for the runs of $input..."
	if ! python3 $CURRENT_DIR/scripts/ena_metadata.py --accessions "$input" --print ${organism_argument:+--organism "${organism_argument//_/ }"} 2>&1 | sed 's/^/      /'; then
		dry_err "ENA returned no usable runs for $input (see above)"
	fi
	if [ -z "$design_custom_local" ]; then
		dry_warn "no design_custom_local: the run will stop and ask for the design interactively, in the sample order printed above"
	else
		dry_ok "design_custom_local: $design_custom_local (one condition per sample, in the sample order printed above; 'auto' drops the replicate numbers of the sample titles)"
	fi
fi

echo "Reference files:"
for f in $(echo "$reference_genome" | tr ',' '\n'); do dry_file "reference genome" "$f"; done
for f in $(echo "$annotation" | tr ',' '\n'); do
	dry_file "annotation" "$f"
	if [ -f "$f" ]; then
		_feat=${optionsFeatureCounts_feat:-exon}; _attr=${optionsFeatureCounts_seq:-gene_name}
		if ! zcat -f "$f" | awk -F'\t' -v ft="$_feat" '!/^#/ && NF>=9 && $3==ft {found=1; exit} END {exit !found}'; then
			dry_err "feature type '$_feat' (optionsFeatureCounts_feat) is not in column 3 of $(basename $f); available: $(zcat -f "$f" | awk -F'\t' '!/^#/ && NF>=9 {print $3}' | sort -u | paste -sd,)"
		elif ! zcat -f "$f" | awk -F'\t' -v ft="$_feat" -v attr="$_attr" '!/^#/ && NF>=9 && $3==ft { n=split($9,pairs,";"); for(i=1;i<=n;i++){ gsub(/^[ \t]+/,"",pairs[i]); split(pairs[i],kv,/[ =]+/); if(kv[1]==attr){found=1; exit} } } END {exit !found}'; then
			dry_err "attribute '$_attr' (optionsFeatureCounts_seq) is not in column 9 of the '$_feat' features of $(basename $f)"
		else
			_nmiss=$(zcat -f "$f" | awk -F'\t' -v ft="$_feat" -v attr="$_attr" '!/^#/ && NF>=9 && $3==ft { ok=0; n=split($9,pairs,";"); for(i=1;i<=n;i++){ gsub(/^[ \t]+/,"",pairs[i]); split(pairs[i],kv,/[ =]+/); if(kv[1]==attr){ok=1; break} } if(!ok) m++ } END {print m+0}')
			if [ "$_nmiss" -gt 0 ]; then
				dry_err "$_nmiss '$_feat' feature(s) of $(basename $f) have no '$_attr' attribute and featureCounts would stop on the first one; set optionsFeatureCounts_seq to 'gene_id' (Ensembl GTFs lack gene_name on unnamed genes) or use a GENCODE GTF"
			else
				dry_ok "featureCounts -t $_feat -g $_attr found in every '$_feat' feature of $(basename $f)"
			fi
		fi
	fi
done
dry_file "transcripts" "$transcripts"
if [ -n "$reference_genome_index" ]; then
	if [ "$aligner" == "bowtie1" ]; then
		ls "${reference_genome_index}".1.ebwt* >/dev/null 2>&1 && dry_ok "bowtie1 index $reference_genome_index" || dry_err "no bowtie1 index at $reference_genome_index (expected ${reference_genome_index}.1.ebwt)"
	else
		dry_file "reference genome index" "$reference_genome_index"
	fi
fi
dry_file "host genome for alignment removal" "$alignment_removal"
dry_file "proteome" "$proteome_faa"
dry_file "human proteome" "$ortho_human_proteome"
dry_file "functional annotation" "$non_reference_funct_enrichm"
dry_file "SortMeRNA database" "$sortmerna_databases"
dry_file "rRNA QC database" "$rrna_qc_databases"
dry_file "miRBase mature FASTA" "$mirna_mature_fa"
dry_file "miRBase hairpin FASTA" "$mirna_hairpin_fa"
[[ $fastp_adapter == /* ]] && dry_file "fastp adapters" "$fastp_adapter"
for d in $(echo "$kraken2_databases" | tr ',' '\n'); do
	if [ -f "$d/hash.k2d" ]; then dry_ok "Kraken2 database $d"; else dry_err "Kraken2 database without hash.k2d: $d"; fi
done
for _i in "${!genome_group_fastas[@]}"; do dry_file "genome of group ${genome_group_labels[$_i]}" "${genome_group_fastas[$_i]}"; done

echo "Output:"
if [ -d "$output_folder/$name" ]; then
	dry_warn "$output_folder/$name exists: step 1 wipes it (resume with -Dm/debug_module to keep it)"
fi
_probe=$output_folder
while [ ! -d "$_probe" ] && [ "$_probe" != "/" ]; do _probe=$(dirname "$_probe"); done
if [ -w "$_probe" ]; then dry_ok "$output_folder is writable ($_probe)"; else dry_err "$output_folder cannot be created or written ($_probe is not writable)"; fi
if [ -n "$TMPDIR_arg" ] && [ "$TMPDIR_arg" != "system" ]; then
	_probe=$TMPDIR_arg
	while [ ! -d "$_probe" ] && [ "$_probe" != "/" ]; do _probe=$(dirname "$_probe"); done
	[ -w "$_probe" ] && dry_ok "TMPDIR $TMPDIR_arg" || dry_err "TMPDIR $TMPDIR_arg cannot be created"
fi

echo "Tools:"
dry_tool pigz parallel fastp samtools featureCounts multiqc Rscript python3 yq perl
[ "$qc_raw_reads" != "no" ] && dry_tool fastqc
case "$aligner" in
	star) dry_tool STAR ;;
	hisat2) dry_tool hisat2 hisat2-build ;;
	kallisto) dry_tool kallisto ;;
	bowtie1) dry_tool bowtie bowtie-build ;;
esac
[ -z "$strand" ] && dry_tool salmon
[ "$dry_mode" != "local" ] && dry_tool fastq-dl
[ "$dry_mode" == "GEO" ] && dry_tool pysradb esearch efetch
[ -n "$kraken2_databases" ] && dry_tool kraken2
[ -n "$sortmerna_databases" ] && dry_tool sortmerna
[ -n "$rrna_qc_databases" ] && dry_tool bowtie2
[ "$bam_normalization" != "None" ] && dry_tool bamCoverage
[ "$repair_pairs" != "no" ] && dry_tool repair.sh
[ "$mirna_novel" == "yes" ] && dry_tool mapper.pl miRDeep2.pl RNAfold
[ "$ortho_detection" != "no" ] && { case "$ortho_detection" in DIAMOND*) dry_tool diamond ;; *) dry_tool blastp makeblastdb ;; esac; }
[ "$network_analyses" == "yes" ] && [ -z "$taxonid" ] && dry_tool taxonkit
_rpk="edgeR,limma,qs2,data.table"
[ "$differential_expr_soft" == "DESeq2" ] && _rpk="$_rpk,DESeq2"
[ "$functional_enrichment_analyses" != "no" ] && _rpk="$_rpk,clusterProfiler,enrichplot,rmarkdown,DT"
[ "$small_rna" == "yes" ] && [ "$functional_enrichment_analyses" != "no" ] && _rpk="$_rpk,multiMiR"
[ "$network_analyses" == "yes" ] && _rpk="$_rpk,WGCNA,STRINGdb"
[ "$exploreDE" != "no" ] && _rpk="$_rpk,anndataR,rhdf5"
_missing_r=$(Rscript -e 'p <- strsplit(commandArgs(TRUE)[1], ",")[[1]]; cat(p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)], sep = " ")' "$_rpk" 2>/dev/null)
if [ -z "$_missing_r" ]; then dry_ok "R packages: $_rpk"; else dry_err "R packages not installed: $_missing_r"; fi

echo "Plan:"
echo "  reads: $dry_mode input, layout $( [ -n "$strand" ] && echo "strand $strand" || echo "strand predicted by salmon" ), $( [ "$small_rna" == "yes" ] && echo "small RNA mode (fastp adapter $small_rna_adapter, $small_rna_length nt, bowtie1)" || echo "aligner $aligner" )"
if [ "$fastp_mode" == "yes" ] || [ "$fastp_adapter" != "no" ] || [ "$fastp_trimming" != "none" ] || [ -n "$fastp_extra_args" ]; then echo "  fastp: quality filtering $fastp_mode, adapter $fastp_adapter, trimming $fastp_trimming, extra '$fastp_extra_args'"; fi
[ -n "$kraken2_databases" ] && echo "  decontamination: Kraken2 ($kraken2_databases, confidence $kraken2_confidence)$( [ -n "$sortmerna_databases" ] && echo " + SortMeRNA" )"
[ -n "$alignment_removal" ] && echo "  host removal: $alignment_removal"
echo "  quantification: featureCounts -t ${optionsFeatureCounts_feat:-exon} -g ${optionsFeatureCounts_seq:-gene_name} ${featureCounts_extra_args}"
echo "  differential expression: $differential_expr_soft (edger_test $edger_test, lfc_threshold $lfc_threshold, filter $filter)$( [ "$covariables" != "none" ] && echo ", covariables $covariables" )$( [ "$batch" == "yes" ] && echo ", batch correction" )$( [ "$differential_expr_comparisons" != "no" ] && echo ", comparisons $differential_expr_comparisons" )"
[ -n "$design_formula" ] && echo "  custom model: $design_formula with contrasts $design_contrasts"
[ "$leave_out_samples" != "no" ] && echo "  leave-out runs: $leave_out_samples"
if [ "$functional_enrichment_analyses" == "no" ]; then
	echo "  functional enrichment: off"
elif [ "$small_rna" == "yes" ]; then
	echo "  functional enrichment: ORA of the $mirna_target_tables multiMiR targets of the DE miRNAs (human, mouse, rat)"
else
	echo "  functional enrichment: ORA (GO, KEGG, Reactome), GSEA $clusterProfiler_gsea, extended analyses $clusterProfiler_full, autoGO/Panther for human and mouse; report images: $report_images"
fi
echo "  network analyses: $network_analyses; time course: $time_course; splicing: ${splicing_option:-no}; exploreDE file: $exploreDE; novel miRNAs: $mirna_novel"
echo "  resources: $cores cores, $number_parallel samples in parallel, $memory_max bytes of RAM, TMPDIR ${TMPDIR_arg:-inside the output folder}"
[ "$debug_step" != "all" ] && echo "  resume from: $debug_step"
[ -n "$end_step" ] && [ "$end_step" != "none" ] && echo "  stop after: $end_step"

echo ""
if [ "$dry_errors" -gt 0 ]; then
	echo -e "\033[1;31mDry run: $dry_errors error(s), $dry_warnings warning(s). Please fix the errors above before running.\033[0m\n"
	exit 1
fi
echo -e "Dry run: no errors, $dry_warnings warning(s). Remove dry_run (or set it to 'no') to start the analysis.\n"
exit 0
