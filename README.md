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

This should work if you already have miniconda3 installed, and also install miniconda3 if not available. Plese keep in mind that in the install.sh script most of the versions of the tools installed by conda are frozen (by means of multiple '.yml' files), so please open an issue or you may need to manually change installed versions if conda fails with new dependencies problems.


2) The less recommended option is to manually install the required software.
If you want to manually install the software, check out in the files '.yml' within the folder external_software/installation the list of required tools, which must be in the PATH when running. Please be aware that many scripts (bash, perl...) within the 'scripts' folder are also used, so you may need to manually change the interpreter in the corresponding shebang statements (first line #!) to ensure that everything works in your system. You may also need to make scripts executable ('chmod' command) and to make source or export so that the PATH variable and others are available for all scripts.


## Comments
We used reanalyzerGSE to automatically reanalyze several transcriptomic datasets and interrogate the expression levels of genes of interest.  
Please cite this reference when using reanalyzerGSE for your publications:

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
