# reanalyzerGSE
reanalyzerGSE is a pipeline to help with and streamline the transcriptomic analyses of various datasets (i.e. microarrays, RNA-seq, scRNA-seq) by automatically reanalyzing raw data submitted to Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) or provided locally by the user. The pipeline is based on several steps (i.e. quality control, alignment to reference genome, quantification and differential gene expression analyses) implementing standard tools and novel scripts.

## Installation
We suggest two alternatives for installation. Please choose one of:

1) The fastest option is to use the folder 'external_software/installation', which contain some of the required software, a script to install dependencies (mainly through miniconda, external_software/installation/install.sh), and a suggestion of the PATH to be set. To install and setup everything required to run reanalyzerGSE, please execute:

To perform a conda-based installition and setup everything required to run reanalyzer, please execute:
```
git clone https://github.com/BioinfoIPBLN/reanalyzerGSE
cd reanalyzerGSE/external_software/installation/
bash install.sh 2>&1 | tee install.sh.log # Check log to ensure successful installation of dependencies
```

This should work if you already have miniconda3 installed, and also install miniconda3 if not available. Plese keep in mind that in the install.sh script most of the versions of the tools installed by conda are frozen (by means of multiple '.yml' files), so please open an issue or you may need to manually change installed versions if conda fails with new dependencies problems.


2) The less recommended option is to manually install the required software.
If you want to manually install the software, check out in the files '.yml' within the folder external_software/installation the list of required tools, which must be in the PATH when running. Please be aware that many scripts (bash, perl...) within the 'scripts' folder are also used, so you may need to manually change the interpreter in the corresponding shebang statements (first line #!) to ensure that everything works in your system. You may also need to make scripts executable ('chmod' command) and to make source or export so that the PATH variable and others are available for all scripts.

## Support
Please report any [issue](https://github.com/BioinfoIPBLN/reanalyzerGSE/issues) or [contact us](mailto:bioinformatica@ipb.csic.es?subject=[GitHub]%20Source%20reanalyzerGSE%20Support) to request support.
