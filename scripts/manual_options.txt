#### Input/output: 
input=
input_geo_reads=
name=project_name
output_folder=
genes=none
TMPDIR_arg=

#### Reference and databases:
reference_genome=
reference_genome_index=
annotation=
transcripts=
kraken2_databases=
sortmerna_databases=
databases_function=GO_Biological_Process_2023,GO_Molecular_Function_2023,GO_Cellular_Component_2023

#### Metadata and sample info:
design_custom_local=
design_custom=no
organism=
taxonid=
number_reads_to_subsample=
batch_vector=
batch_biological_covariates=
covariables=none
target=

#### Activate alternative modes:
debug_step=all
qc_raw_reads=yes
full_differential_analyses=yes
functional_enrichment_analyses=yes
clusterProfiler_full=no
batch=no
bed_mode=no
deconvolution=no
perform_volcano_venn=yes
aPEAR_execution=no
tidy_tmp_files=no
convert_tables_excel=no
time_course=no
network_analyses=no
auto_panther_log=no

#### Processing parameters:
strand=
filter=standard
optionsFeatureCounts_feat=exon
optionsFeatureCounts_seq=gene_name
aligner=star
differential_expr_soft=edgeR
fastp_mode=no
fastp_adapter=no
fastp_trimmnig=none
time_course_std=1
time_course_fuzz=0

#### Filtering out samples/comparisons:
GSM_filter=
stop=no
pattern_to_remove=none
differential_expr_comparisons=no

#### Functional enrichment/networking analyses
clusterProfiler_method=fdr
clusterProfiler_universe=detected
clusterProfiler_minGSSize=10
clusterProfiler_maxGSSize=500
panther_method=FDR
revigo_threshold_similarity=0.7

#### Performance:
cores=30
cores_reads_to_subsample=10
cores_index=30
memory_max=257698037760
number_parallel=10
compression_level=9
aligner_index_cache=yes
kraken2_fast=no
