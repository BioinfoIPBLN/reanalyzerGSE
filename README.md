# reanalyzerGSE
reanalyzerGSE is a pipeline to help with and streamline the transcriptomic analyses of various datasets (i.e. microarrays, RNA-seq, scRNA-seq) by automatically reanalyzing raw data submitted to Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) or provided locally by the user. The pipeline is based on several steps (i.e. quality control, alignment to reference genome, quantification and differential gene expression analyses) implementing standard tools and novel scripts.

## Installation
We suggest two alternatives for installation. Please choose one of:

1) The fastest option is to use the folder 'external_software', which contain some of the required software (i.e. miARma-seq), a script to install dependencies (mainly through miniconda, 'external_software/installation/install.sh'), and a suggestion of the PATH to be set ('source_path'). To perform a conda-based installation and setup everything required to run reanalyzerGSE, please execute:

```
git clone https://github.com/BioinfoIPBLN/reanalyzerGSE
cd reanalyzerGSE/external_software/installation/
bash install.sh 2>&1 | tee install.sh.log # Check log to ensure successful installation of dependencies
```

This should work if you already have miniconda3 installed, and also install miniconda3 if not available. Plese keep in mind that in the install.sh script most of the versions of the tools installed by conda are frozen (by means of multiple '.yml' files), so please open an issue or you may need to manually change installed versions if conda fails with new problems dependency-related.


2) The less recommended option is to manually install the required software.
If you want to manually install the software, check out in the files '.yml' within the folder external_software/installation the list of required tools, which must be in the PATH when running. Please be aware that many scripts (bash, perl...) within the 'scripts' folder are also used, so you may need to manually change the interpreter in the corresponding shebang statements (first line #!) to ensure that everything works in your system. You may also need to make scripts executable ('chmod' command) and to make source or export so that the PATH variable and others are available for all scripts.


## Quick start / Minimal examples
```
source reanalyzerGSE/external_software/source_path # To set up the PATH if you have followed option 1 for installation above
cd reanalyzerGSE/
### Case examples of mouse transcriptomic datasets, analyzed in a machine with 30 cores and 200GB RAM available.
cores=30
pigz -p $cores -dkc $PWD/test_data/GRCm39.primary_assembly.genome.fa.gz # To uncompress the reference genome and annotation
pigz -p $cores -dkc $PWD/test_data/gencode.vM28.annotation.gtf.gz 
pigz -p $cores -dkc $PWD/test_data/Mus_musculus.GRCm39.cdna.all.fa.gz
# A) Local raw data (subset, 500k reads) simultaneously processing 4 samples and highligthing the genes B4galt3 and Chd1
reanalyzerGSE.pk.sh -i $PWD/test_data -r $PWD/test_data/GRCm39.primary_assembly.genome.fa -a $PWD/test_data/gencode.vM28.annotation.gtf -s reverse -o $PWD/test_data_out/test_data_1/ -p $cores -P 4 -g B4galt3,Chd1 -M 214748364800 2>&1 | tee -a $PWD/test_data_out/test_data_1.log

# B) The GEO entry GSE118451, simultaneously processing 6 samples, highligthing the genes B4galt3 and Chd1, and predicting strandness from the transcript sequences
reanalyzerGSE.pk.sh -i GSE118451 -r $PWD/test_data/GRCm39.primary_assembly.genome.fa -a $PWD/test_data/gencode.vM28.annotation.gtf -t $PWD/test_data/Mus_musculus.GRCm39.cdna.all.fa -o $PWD/test_data_out/test_data_2/ -p $cores -P 6 -g B4galt3,Chd1 -M 214748364800 2>&1 | tee -a $PWD/test_data_out/test_data_2.log
```
These test runs will take ~X minutes and ~X minutes, respectively. The generation of the required indexes (i.e. for Salmon and STAR) will be done automatically if required (if not present beforehand in the subfolder 'indexes'), and take around ~X minutes.

Please go through the log files 'test_data_X.log' to get the details on the pipeline processing steps and final output.


## reanalyzerGSE arguments
Please refer to the help page for futher details:
```
ILRA.sh -h

usage: ILRA.sh [options]
		-h | -help # Type this to get help
		-a | -assembly # Name of the long reads assembly to correct (FASTA format, can be gzipped)
		-f | -filter_contig_size # Size threshold to filter the contigs (bp)
		-F | -Do you have even Illumina reads coverage? If yes contigs in the long reads assembly not covered by Illumina short reads by a certain threshold will be filtered out ('no' /'yes' by default)
		-I | -Illumina_reads # Root name of the paired-end Illumina reads (FASTQ format, must be gzipped and name root_1.fastq.gz and root_2.fastq.gz prior execution)
		-C | -Correction_Illumina_reads # Whether illumina reads are provided and all steps of correction should be performed ('no' /'yes' by default)
		-R | -Range_insert_size # Insert size range of the Illumina reads to use in mapping by SMALT (bp)
		-i | -iterations_iCORN2 # Number of iterations to perform in iCORN2
		-r | -reference # Reference file (full pathway, FASTA format)
		-g | -gff_file # Reference annotation file (full pathway, GFF format)
		-c | -corrected_reads # Corrected long reads for circulatization of the contigs containing the strings by -s and -S (FASTQ format, can be gzipped)
		-s | -seq_circularize_1 # Regex pattern to find in the contig names and circularize
		-S | -Seq_circularize_2 # Regex pattern to find in the contig names and circularize
		-L | -Long_reads_technology # Technology of long reads sequencing (pb/ont)
		-T | -Taxon_ID # NCBI taxon ID or a scientific name
		-e | -ending_telomere_seq_1 # Telomere-associated sequence to search in the first strand
		-E | -Ending_telomere_seq_2 # Telomere-associated sequence to search in the second strand
		-o | -output # Output folder (absolute pathway)
		-n | -name # Base name of the output file
		-t | -threads # Number of cores to use in multithreaded steps
		-d | -debug_step # For debug, step to remove the content of the corresponding folder and resume a failed run ('step1', 'step2a', 'step2b', 'step3', 'step4', 'step4i', 'step5', 'step6', or 'step7')
		-D | -databases # Folder for storing the databases for the decontamination step (by default, 'databases' under ILRA main folder)
		-K | -Kraken2_fast_mode # Kraken2 fast mode, consisting on copying the Kraken2 database to /dev/shm (RAM) so execution is faster ('yes' /'no' by default)
		-k | -Kraken2_databases # Folder within the folder databases (-D) containing the database used by Kraken2 (by default, 'standard_eupathdb_48_kraken2_db')
		-b | -block_size # Block size for parallel processing (by default, 10)
		-p | -pilon # Whether to use pilon instead of iCORN2 ('yes'/'no' by default)
		-P | -parts_icorn2_split # Number of parts to split the input sequences of iCORN2 before processing them (0 by default, which means no splitting)
		-A | -abacas2_split # Number of parts to split and process in parallel in ABACAS2 (by default the argument -b, block_size, but may be necessary to decrease due to memory issues)
		-B | -abacas2_blast # Whether to do blast within ABACAS2 to compare with the reference and display in ACT (1 by default, which means blasting, or 0)
		-q | -quality_assesment # Whether to execute a final step for assessing the quality of the corrected assembly, gathering sequences, analyzing telomeres... etc ('no'/'yes' by default)
		-M | -java_memory # Max Java memory (heap space) to be used ('XXg', by default 200g=200GB used)
		-l | -low_memory # Activate low memory mode for iCORN2 ('yes'/'no' by default)
		-m | -mode # Add 'taxon' to execute decontamination based on taxonomic classification by kraken2, add 'blast' to execute decontamination based on BLAST against databases as requested by the DDBJ/ENA/Genbank submission, add 'both' to execute both approaches, and add 'light' to execute ILRA in light mode and skip these steps (default)
```

Parameters are not positional. If you did not provide a required parameter, the pipeline may exit or use default values if possible (check the help page above, the log after execution, or the 'Arguments / Variables' section in the ILRA main script 'ILRA.sh'). For example, if the argument '-f' is not provided, contigs shorter than 5,000 bp by default will be filtered out. 

In general, from an assembly as input (argument '-a'), ILRA is going to provide a polished assembly as output with the argument '-n' as name (the file 'name.ILRA.fasta', in the folder provided by the argument '-o').

Please do provide or not the arguments '-C', '-c', '-R' and '-I' to indicate whether to use short reads to perform error correction (iCORN2 / Pilon iteratively, with argument '-i' as number of iterations), and to find and filter out overlapping contigs (if coverage by Illumina short reads is even, argument 'F'). Please do provide the long reads sequencing technology used with the argument '-L'.

Depending on whether you provided a reference genome (argument '-r'), reordering and renaming of the contigs (ABACAS2 and the argument '-B' to perform blasting) is going to be skipped, and assessment by QUAST would be run without the reference. Similarly, the availability of a reference annotation (argument '-g') would determine the mode to run QUAST, or the presence of certain names in the contigs marking the sequences to circularize (arguments '-s' and '-S') would mean that Circlator is executed or not. The debug mode (argument '-d' makes possible to resumen the execution of ILRA from a particular step). The argument '-p' determine whether Pilon should be used for short reads correction instead of iCORN2 (default 'no'). The argument '-q' determines whether a final extra step for assessing the quality and completeness of the corrected assembly (i.e., QUAST, BUSCO, gathering sequences, looking in the telomeres for the sequenes provided by the arguments '-e' and '-E'...) is included (default 'yes'). The argument 'M' is required to control the maximum RAM memory used, and the argument '-l' activates a 'low memory' mode at the expense of more processing time. The arguments '-b', '-P' and '-A' provide the number of parts to split the sequences and to simultaneously process in ILRA, iCORN2/Pilon and ABACAS2, respectively. The argument '-t' provides the number of cores to be used when multithreading is possible.

Finally, ILRA can be run in alternative modes (argument '-m'): 
* '-m taxon': To perform decontamination based on taxonomic classification, which would be more computationally expensive.
* '-m blast': To perform decontamination and formatting for online submission based on blasting against databases, which would be less computationally expensive.
* '-m both': To perform both. 
* '-m light': To skip decontamination and expedite the process (default if argument not provided).
The location of the downloaded databases to be used in the decontamination step should be provided with the argument '-D' and the location of the kraken2 database should be provided with the argument '-k'. The taxon id of the organism of interest that should be kept when filtering must be provided with the argument '-T'.


## Comments
We used reanalyzerGSE to automatically reanalyze several transcriptomic datasets and interrogate the expression levels of genes of interest.  
Please cite this reference when using reanalyzerGSE for your publications: (PENDING TO UPDATE ONCE WE HAVE THE BIORXIV)

> From contigs to chromosomes: automatic Improvement of Long Read Assemblies (ILRA)
> 
> José L Ruiz, Susanne Reimering, Mandy Sanders, Juan David Escobar-Prieto, Nicolas M. B. Brancucci, Diego F. Echeverry, Abdirahman I. Abdi, Matthias Marti, Elena Gomez-Diaz, Thomas D. Otto
> 
> bioRxiv 2021.07.30.454413; doi: https://doi.org/10.1101/2021.07.30.454413
```
@article {Ruiz2021.07.30.454413,
	author = {Ruiz, Jos{\'e} L and Reimering, Susanne and Sanders, Mandy and Escobar-Prieto, Juan David and Brancucci, Nicolas M. B. and Echeverry, Diego F. and Abdi, Abdirahman I. and Marti, Matthias and Gomez-Diaz, Elena and Otto, Thomas D.},
	title = {From contigs to chromosomes: automatic Improvement of Long Read Assemblies (ILRA)},
	elocation-id = {2021.07.30.454413},
	year = {2021},
	doi = {10.1101/2021.07.30.454413},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Recent advances in long read technologies not only enable large consortia to aim to sequence all eukaryotes on Earth, but they also allow many laboratories to sequence their species of interest. Although there is a promise to obtain {\textquoteright}perfect genomes{\textquoteright} with long read technologies, the number of contigs often exceeds the number of chromosomes significantly, containing many insertion and deletion errors around homopolymer tracks. To overcome these issues, we implemented ILRA to correct long reads-based assemblies, a pipeline that orders, names, merges, and circularizes contigs, filters erroneous small contigs and contamination, and corrects homopolymer errors with Illumina reads. We successfully tested our approach to assemble the genomes of four novel Plasmodium falciparum samples, and on existing assemblies of Trypanosoma brucei and Leptosphaeria spp. We found that correcting homopolymer tracks reduced the number of genes incorrectly annotated as pseudogenes, but an iterative correction seems to be needed to reduce high numbers of homopolymer errors. In summary, we described and compared the performance of a new tool, which improves the quality of long read assemblies. It can be used to correct genomes of a size of up to 300 Mb.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2021/08/01/2021.07.30.454413},
	eprint = {https://www.biorxiv.org/content/early/2021/08/01/2021.07.30.454413.full.pdf},
	journal = {bioRxiv}
}
```


## Support
Please report any [issue](https://github.com/BioinfoIPBLN/reanalyzerGSE/issues) or [contact us](mailto:bioinformatica@ipb.csic.es?subject=[GitHub]%20Source%20reanalyzerGSE%20Support) to request support.
