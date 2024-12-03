#### Input/output: 
input=
input_geo_reads=
name=project_name
output_folder=
number_reads=
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
organism_argument=
taxonid=
batch_vector=
batch_biological_covariates=
covariables=none
differential_expr_comparisons=no
target=

#### Activate alternative modes:
qc_raw_reads=yes
perform_differential_analyses=yes
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
debug_step=all

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
minstd=1
mestimate=0

#### Filtering out samples:
GSM_filter=
stop=no
pattern_to_remove=none

#### Functional enrichment/networking analyses
clusterProfiler_method=fdr
clusterProfiler_universe=detected
clusterProfiler_minGSSize=10
clusterProfiler_maxGSSize=500
panther_method=FDR
rev_thr=0.7

#### Performance:
cores=30
cores_reads_to_subsample=10
indexthreads=30
memory_max=257698037760
number_parallel=10
compression_level=9
aligner_index_cache=yes
kraken2_fast=no
