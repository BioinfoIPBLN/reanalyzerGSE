# reanalyzerGSE
reanalyzerGSE is a pipeline to assist with and streamline transcriptomic analyses of various datasets (i.e. microarrays, RNA-seq, scRNA-seq) by automatically reanalyzing raw data submitted to public databases like [GEO](https://www.ncbi.nlm.nih.gov/geo/), [ENA](https://www.ebi.ac.uk/ena/browser/home) or [SRA](https://www.ncbi.nlm.nih.gov/sra). Local data can also be provided by the user. The pipeline is based on several steps implementing standard tools and novel scripts. (i.e. data download, quality control, alignment to reference genome, quantification, differential gene expression analyses, functional enrichment analyses...)

## Installation
We suggest alternatives for installation. Please choose one of:

1) An Apptainer/Singularity container (~6 GB) is provided. You can either:

1.1) Use the .def file to create the .sif image by executing:
```
git clone https://github.com/BioinfoIPBLN/reanalyzerGSE
apptainer build reanalyzerGSE.sif reanalyzerGSE/external_software/installation/reanalyzerGSE.def | tee -a reanalyzerGSE.sif.build.log
```
1.2) Download the ready-to-use .sif image:
```
wget -q https://bit.ly/reanalyzer_appt_image -O reanalyzerGSE.sif
```


2) Another option is to use the folder 'external_software', which contain some of the required software (i.e. miARma-seq), and within the 'external_software/installation' folder a wrapper script installs and configures all dependencies (mainly through miniconda and pip, 'external_software/installation/install.sh'). To perform a conda-based installation and setup everything required to run reanalyzerGSE, please execute:
```
git clone https://github.com/BioinfoIPBLN/reanalyzerGSE
bash reanalyzerGSE/external_software/installation/install.sh 2>&1 | tee -a reanalyzerGSE/external_software/installation/install.sh.log # Check log to make sure that installation of all dependencies has been succesful
```

This should work if you already have miniconda3 installed, or install miniconda3 within the reanalyzerGSE folder if not available or if you have kept it out of the PATH. Plese keep in mind that in the 'install.sh' script most of the versions of the tools installed by conda are frozen (by means of multiple '.yml' files corresponding to different environments), so please open an issue or try to install with conda if there are dependency-related problems or any software is not installed.


3) The less recommended option is to manually install the required software.
If you want to manually install the software, check out the list of required tools in the .def file (Apptainer/Singularity) or in the files '.yml' within the folder 'external_software/installation'. Please be aware that many scripts (bash, perl...) within the 'scripts' folder are also used, so you may need to manually change the interpreter in the corresponding statements (first line #!) to ensure that everything works in your system. You may also check out the 'install.sh' script to conform to other needs, such as making scripts executable ('chmod' command).


## Quick start / Minimal examples
Please go to test_data/README and follow instructions.


## reanalyzerGSE arguments
Please be aware that inline parameters can be provided (see 'reanalyzerGSE.sh -h' or the [options script](https://github.com/BioinfoIPBLN/reanalyzerGSE/blob/main/scripts/parse_options.sh)) or a yaml file with all arguments can be provided via '-options' (see [template](https://github.com/BioinfoIPBLN/reanalyzerGSE/blob/main/scripts/config_template.yaml))
```
# Recommended execution:
reanalyzerGSE.sh -options config.yaml | tee -a output.log
```

See also a [quick example](https://github.com/BioinfoIPBLN/reanalyzerGSE/tree/main/test_data)

Parameters are not positional. If you did not provide a required parameter, the pipeline may exit or use default values if possible. For example, if not provided, strandedness will be predicted using Salmon and transcript sequences would be required. Other input such as reference genome and annotation are always required for the alignment and quantifying steps.

Several different pipeline modes can be activated but others are automatic. For example, depending on whether a folder containing raw sequences (.fastq.gz) or a GEO entry (GSEXXXXX) are provided as input, GEO reanalysis or local analysis modes will be triggered.

It is advisable to fine tune parameters such as the maximum RAM available, number of cores to be used in steps allowing multithreading, or number of samples to be processed simultaneously when possible (mostly leveraging GNU's parallel).

An updated version of [miARma-seq](https://github.com/eandresleon/miARma-seq) has been included in reanalyzerGSE [here](https://github.com/BioinfoIPBLN/reanalyzerGSE/tree/main/external_software/miARma-seq).

Please refer to the help ('-h') or contact us for any further clarification.

## Output
Please refer to the [wiki](https://github.com/BioinfoIPBLN/reanalyzerGSE/wiki) for the output of a test run.

## Comments
Please cite this reference when using reanalyzerGSE for your publications:

> Ruiz, J. L., Terrón-Camero, L. C., Castillo-González, J., Fernández-Rengel, I., Delgado, M., Gonzalez-Rey, E., & Andrés-León, E. (2023). reanalyzerGSE: tackling the everlasting lack of reproducibility and reanalyses in transcriptomics. bioRxiv, 2023-07. https://doi.org/10.1101/2023.07.12.548663

```
@article{ruiz2023reanalyzergse,
  title={reanalyzerGSE: tackling the everlasting lack of reproducibility and reanalyses in transcriptomics},
  author={Ruiz, Jose L and Terr{\'o}n-Camero, Laura Carmen and Castillo-Gonz{\'a}lez, Julia and Fern{\'a}ndez-Rengel, Iv{\'a}n and Delgado, Mario and Gonzalez-Rey, Elena and Andr{\'e}s-Le{\'o}n, Eduardo},
  journal={bioRxiv},
  pages={2023--07},
  year={2023},
  publisher={Cold Spring Harbor Laboratory}
}
```


## Support
Please report any [issue](https://github.com/BioinfoIPBLN/reanalyzerGSE/issues) or [contact us](mailto:bioinformatica@ipb.csic.es?subject=[GitHub]%20Source%20reanalyzerGSE%20Support) to request support or any further clarification.
