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

reanalyzerGSE.pk.sh -h

usage: [options]
		-h | -help # Type this to get help
		-i | -GEO_ID # GEO_ID (GSEXXXXXX, separated by comma if more than one) or folder containing raw reads (please provide full absolute path, e.g. /path/folder_name/)
		-n | -name # Name of the project/folder to create and store results
		-o | -output_folder # Destination folder
		-p | -cores # Number of cores
		-P | -number_parallel # Number of files to be processed in parallel (10 by default)		
		-r | -reference_genome # Reference genome to be used (.fasta file or .gz, absolute pathway)
		-a | -annotation_file # Reference annotation to be used (.gtf file, absolute pathway)
		-t | -transcripts # Referece transcripts to be used (.fasta cDNA file, absolute pathway, only used if '-s' argument not provided so salmon prediction of strandness is required)
		-s | -strandness # Strandness of the library ('yes, 'no', 'reverse'). If not provided and '-t' used, this would be predicted by salmon. Please use this parameter if prediction not correct, see explanations in for example in bit.ly/strandness0 and bit.ly/strandness
		-g | -genes # Genes to highlight their expression in plots (one or several, separated by comma and no space)
		-G | -GSM_filter # GSM ids (one or several, separated by comma and no space) within the GSE entry to restrict the analysis to. An alternative to requesting a stop with -S to reorganize the downloaded files manually
		-f | -filter # Threshold of gene counts to use (bin to capture the lower expressed genes, e.g. Cort, or standard, by default)
		-b | -batch # Batch effect present? (no by default, yes if correction through Combat-seq and model is to be performed)
		-d | -design # Manually specifying the experimental design (no by default, if yes the assignment to groups for each sample must be provided in the prompt when asked with a comma-separated list of the same length than the number of samples)
		-S | -stop # Manual stop so the automatically downloaded files can be manually modified (yes or no, by default)
		-M | -memory # Max RAM memory to be used by STAR ('XXXXXXXXXX bytes', by default 257698037760=240GB used)
		-m | -miARma_seq_path # A default is used if not provided...
```

Parameters are not positional. If you did not provide a required parameter, the pipeline may exit or use default values if possible (check the help page above, the log after execution, or the 'arguments and variables' first section in the main script 'reanalyzerGSE.pk.sh'). For example, if the argument '-s' is not provided, strandness will be predicted using Salmon and requiring transcript sequences, so the pipeline would eixt if not provided with the argument '-t'. 

In general, from a folder containing raw sequences (.fastq.gz) or a GEO entry (GSEXXXXX) as input (argument '-i'), reanalyzerGSE is going to provide a standard transcriptomic analysis named with the argument '-n', in the folder provided by the argument '-o'. The argument '-p' provides the number of cores to be used when multithreading is possible and '-P' the number of files to be processed simultaneously when possible, leveraging GNU's parallel.

An improved version of miARma-seq has been included in reanalyzerGSE and used by default (subfolder 'external_software'). if not detected or available, the path to a working version of miARma-seq must be provided with the argument '-m' | -miARma_seq_path # A default is used if not provided...Please do provide or not the arguments '-C', '-c', '-R' and '-I' to indicate whether to use short reads to perform error correction (iCORN2 / Pilon iteratively, with argument '-i' as number of iterations), and to find and filter out overlapping contigs (if coverage by Illumina short reads is even, argument 'F'). Please do provide the long reads sequencing technology used with the argument '-L'.

Depending on whether you provided a reference genome (argument '-r'), reordering and renaming of the contigs (ABACAS2 and the argument '-B' to perform blasting) is going to be skipped, and assessment by QUAST would be run without the reference. Similarly, the availability of a reference annotation (argument '-g') would determine the mode to run QUAST, or the presence of certain names in the contigs marking the sequences to circularize (arguments '-s' and '-S') would mean that Circlator is executed or not. The debug mode (argument '-d' makes possible to resumen the execution of ILRA from a particular step). The argument '-p' determine whether Pilon should be used for short reads correction instead of iCORN2 (default 'no'). The argument '-q' determines whether a final extra step for assessing the quality and completeness of the corrected assembly (i.e., QUAST, BUSCO, gathering sequences, looking in the telomeres for the sequenes provided by the arguments '-e' and '-E'...) is included (default 'yes'). The argument 'M' is required to control the maximum RAM memory used, and the argument '-l' activates a 'low memory' mode at the expense of more processing time. The arguments '-b', '-P' and '-A' provide the number of parts to split the sequences and to simultaneously process in ILRA, iCORN2/Pilon and ABACAS2, respectively. 

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
