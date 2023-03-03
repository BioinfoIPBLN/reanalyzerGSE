# reanalyzerGSE
reanalyzerGSE is a pipeline to help with and streamline the transcriptomic analyses of various datasets (i.e. microarrays, RNA-seq, scRNA-seq) by automatically reanalyzing raw data submitted to Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) or provided locally by the user. The pipeline is based on several steps (i.e. quality control, alignment to reference genome, quantification and differential gene expression analyses) implementing standard tools and novel scripts.

## Installation
We suggest two alternatives for installation. Please choose one of:

1) The fastest option is to use the folder 'external_software', which contain some of the required software (i.e. miARma-seq), a script to install dependencies (mainly through miniconda, 'external_software/installation/install.sh'), and a suggestion of the PATH to be set ('source_path'). To perform a conda-based installation and setup everything required to run reanalyzerGSE, please execute:

```
git clone https://github.com/BioinfoIPBLN/reanalyzerGSE
bash reanalyzerGSE/external_software/installation/install.sh 2>&1 | tee reanalyzerGSE/external_software/installation/install.sh.log # Check log to ensure successful installation of dependencies
```

This should work if you already have miniconda3 installed, and also install miniconda3 if not available. Plese keep in mind that in the install.sh script most of the versions of the tools installed by conda are frozen (by means of multiple '.yml' files), so please open an issue or you may need to manually change installed versions if conda fails with new problems dependency-related.


2) The less recommended option is to manually install the required software.
If you want to manually install the software, check out in the files '.yml' within the folder external_software/installation the list of required tools, which must be in the PATH when running. Please be aware that many scripts (bash, perl...) within the 'scripts' folder are also used, so you may need to manually change the interpreter in the corresponding shebang statements (first line #!) to ensure that everything works in your system. You may also need to make scripts executable ('chmod' command) and to make source or export so that the PATH variable and others are available for all scripts.


## Quick start / Minimal examples
```
bash reanalyzerGSE/external_software/source_path.sh # To get a suggestion of the PATH to export if you have followed option 1 for installation above
# export the PATH following the printed instructions with the previous command
cores=30
cd reanalyzerGSE/test_data
wget -q https://bit.ly/case_examples; unzip -qq case_examples; mv test_data/* .; rm -r test_data case_examples
pigz -p $cores -d *.fa.gz *.gtf.gz; mkdir -p ../references; mv $(ls | egrep '.fa$|gtf$') ../references # To uncompress the reference genome and annotation and move them to a subfolder
cd ../
### Case examples of mouse transcriptomic datasets, analyzed in a machine with 30 cores and 200GB RAM available.
# A) Local raw data (subset, 500k reads) simultaneously processing 4 samples and highligthing the genes B4galt3 and Chd1
reanalyzerGSE.pk.sh -i $PWD/test_data -r $PWD/references/GRCm39.primary_assembly.genome.fa -a $PWD/references/Mus_musculus.GRCm39.109.gtf -s reverse -o $PWD/test_data_out/test_data_1/ -p $cores -P 4 -g B4galt3,Chd1 -M 214748364800 2>&1 | tee -a $PWD/test_data_out/test_data_1.log

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
		-i | -input # GEO_ID (GSEXXXXXX, separated by comma if more than one) or folder containing raw reads (please provide full absolute path, e.g. /path/folder_name/)
		-n | -name # Name of the project/folder to create and store results
		-o | -output_folder # Destination folder
		-p | -cores # Number of cores
		-P | -parallel_number # Number of files to be processed in parallel (10 by default)
		-r | -reference_genome # Reference genome to be used (.fasta file or .gz, absolute pathway)
		-a | -annotation # Reference annotation to be used (.gtf file, absolute pathway)
		-t | -transcripts # Referece transcripts to be used (.fasta cDNA file, absolute pathway, only used if '-s' argument not provided so salmon prediction of strandness is required)
		-s | -strandness # Strandness of the library ('yes, 'no', 'reverse'). If not provided and '-t' used, this would be predicted by salmon. Please use this parameter if prediction not correct, see explanations in for example in bit.ly/strandness0 and bit.ly/strandness
		-g | -genes # Genes to highlight their expression in plots (one or several, separated by comma and no space)
		-G | -GSM_filter # GSM ids (one or several, separated by comma and no space) within the GSE entry to restrict the analysis to. An alternative to requesting a stop with -S to reorganize the downloaded files manually
		-R | -reads_to_subsample # Number of reads to subsample the sequences before the analyses
		-f | -filter # Threshold of gene counts to use ('bin' to capture the lower expressed genes, or 'standard', by default)
		-b | -batch # Batch effect present? (no by default, yes if correction through Combat-seq and model is to be performed, and info is going to be required in prompts)
		-d | -design_custom # Manually specifying the experimental design ('no' by default and if 'yes', please expect an interactive prompt after data download from GEO, and please enter the assignment to groups when asked in the terminal, with a comma-separated list of the same length than the number of samples)
		-S | -stop # Manual stop so the automatically downloaded files can be manually modified ('yes' or no, by default)
		-M | -memory_max # Max RAM memory to be used by STAR in bytes (by default 257698037760, or 240GB, used)
		-m | -miARma_seq_path # By default, the updated version within the main folder is used and reanalyzerGSE may not work with other, but try and provide other path if required
```

Parameters are not positional. If you did not provide a required parameter, the pipeline may exit or use default values if possible (check the help page above, the log after execution, or the 'arguments and variables' first section in the main script 'reanalyzerGSE.pk.sh'). For example, if the argument '-s' is not provided, strandness will be predicted using Salmon and requiring transcript sequences, so the pipeline would exit if not provided with the argument '-t'. Similarly, reference genome and annotation are likely going to be required for the alignment and quantifying steps (arguments '-r' and '-a').

In general, from a folder containing raw sequences (.fastq.gz) or a GEO entry (GSEXXXXX) as input (argument '-i'), reanalyzerGSE is going to provide a standard transcriptomic analysis named with the argument '-n', in the folder provided by the argument '-o'. The argument '-p' provides the number of cores to be used when multithreading is possible, '-P' the number of files to be processed simultaneously when possible, leveraging GNU's parallel, and '-M' the maximum amount of RAM memory available (required mostly for the index generation and alignment step). If the reanalysis of a GEO entry is requested, all the samples in the database will be included by default. The parameter '-G' allows to restrict the analysis to some of them (providing GSMXXXXX ids), while '-S' interrupts the pipeline until the user has made manual changes to the downloaded files, and '-d' allows the user to input a custom experimental design, to be used instead of the automatically-detected one.

An improved version of miARma-seq has been included in reanalyzerGSE and used by default (subfolder 'external_software'). if not detected or available, the path to a working version of miARma-seq must be provided with the argument '-m'. reanalyzer GSE also includes a module to perform batch correction (the design matrix must be provided in a prompt after using the flag '-b'), the possibility to use multiple thresholds when filtering out low expressed genes (argument '-f'), and a module to output plots from the differential gene expression analyses, highlighting the genes provided with the argument '-g'.


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
