### This is the previous README file from miARma-Seq. Within reanalyzerGSE, we have updated miARma to v2.0 (mainly to speed it up and boost its performance and memory management), and will be maintained in this repo.
# miARma WebPage #
This is a new repository created in Oct 2021 to store the source code and the guides from miARma-seq. The guides will be available soon in the wiki.

This repository contains the source code published in:
- Andrés-León E, Rojas AM. miARma-Seq, a comprehensive pipeline for the simultaneous study and integration of miRNA and mRNA expression data. Methods. 2019 Jan 1;152:31-40. doi: 10.1016/j.ymeth.2018.09.002. Epub 2018 Sep 22. PMID: 30253202.
- Andrés-León E, Núñez-Torres R, Rojas AM. miARma-Seq: a comprehensive tool for miRNA, mRNA and circRNA analysis. Sci Rep. 2016 May 11;6:25749. doi: 10.1038/srep25749. Erratum in: Sci Rep. 2018 Jan 08;8:46928. PMID: 27167008; PMCID: PMC4863143.

# miARma #

miARma is a fully customizable pipeline for NGS transcriptome analyses. Including gene/transcripts, miRNAs and circRNAs expression measurements.
Created at Computational Biology and Bioinformatics Group (CbBio)
Institute of Biomedicine of Seville. IBIS (Spain)
Modified and Updated at Bioinformatics Unit at IPBLN-CSIC (Institue for Parasitology and Biomedicine Lopez-Neyra, CSIC).
Granada (Spain). 
Copyright (c) 2019 IBIS & IPBLN. All rights reserved.

### miARma 1.7.6 release (07/Oct/19) ###
The latest Bioconductor versions that use R3.5 and R.6 install the packages differently than the other versions. So miARma has been modified to take this into account.
 
### miARma 1.7.5 release (22/Sep/18) ###
New utilities have been included and publisehd in a new [scientific article](https://www.ncbi.nlm.nih.gov/pubmed/30253202). An example to integrate miRNA and mRNA data to infer potential regulation partners is included in a [GitHub repository](https://github.com/eandresleon/miRNA-mRNA_Integration).

 * The possibility of integrating data from miRNAs and mRNAs
 * A statistical correlation (Pearson/Spearman) can be calculated to infer linked miRNA/mRNA expression profiles
 * miRNA/mRNA target prediction based on statstical correlated pairs
 

### miARma 1.7.2 release (18/Dec/17) ###
Minor bugs fixed.
 
 * New Ensembl BiomaRt URL used
 * Order columns from ReadCount section without checking samples names
 * Fixed a bug in the TargetPrediction
 
### miARma 1.7.1 release (21/Aug/17) ###
Minor bugs fixed eg.
 * No aligned reads in hisat2 paired end analysis added.

Added stuff:
* Unaligned files are compressed
* miRDeeparam added to include parameters to miRDeep execution


### miARma 1.7.0 release (09/Aug/17) ###
 * [Hisat v2.1.0](https://ccb.jhu.edu/software/hisat2/index.shtml) has been included as a mRNA aligner.
 * [STAR v020201](https://github.com/alexdobin/STAR/) has been included as a mRNA aligner.

### 1. Included in miARma ###

#### Quality Software ####
* [Fastqc v0.11.5](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#### Trimming Softare ####
* [CutAdapt v1.9.1](https://cutadapt.readthedocs.org/en/stable/)
* [Minion v15-065](ftp://ftp.ebi.ac.uk/pub/contrib/enrightlab/kraken/reaper/src/reaper-latest/doc/minion.html)
* [Reaper v15-065](http://www.ebi.ac.uk/~stijn/reaper/reaper.html)
#### Aligners ####
* [Bowtie v1.1.2](http://bowtie-bio.sourceforge.net/index.shtml)
* [Bowtie v2.2.8](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [TopHat v2.1.1](http://ccb.jhu.edu/software/tophat/index.shtml)
* [BWA v0.7.13](http://bio-bwa.sourceforge.net/)
* [Hisat v2.1.0](https://ccb.jhu.edu/software/hisat2/index.shtml)
* [STAR v020201](https://github.com/alexdobin/STAR/)
#### Entity Quantification####
* [feactureCounts v1.5.0-p1](http://bioinf.wehi.edu.au/featureCounts/)
* [Ciri v1.2](http://sourceforge.net/projects/ciri/files/?source=navbar)
* [Ciri v2.0.1](http://sourceforge.net/projects/ciri/files/?source=navbar)
* [miRDeep v2](https://www.mdc-berlin.de/8551903/en/)
#### Others ####
* [RNAfold v2.2.4](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)
* [Samtools v1.3](http://samtools.sourceforge.net/)

### 2. Pre-requisites ###

miARma-Seq is a tool that provides an easy and common interface to various analysis software. It also intends to reduce to the minimum the number of dependencies. Nevertheless, some basic programs listed below must be correctly installed:

* [Perl v5.6.0 or higher.](http://www.cpan.org/src/5.0/perl-5.6.1.tar.gz)
* [Java JDK v.1.6. or higher.](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
* [R environment v.3.2 or higher.](http://www.r-project.org/)
* [Bioconductor v.1.3 or higher.](https://www.bioconductor.org/install/)

#### Compilers: #####
+ Apple:
    - [Xcode](https://itunes.apple.com/es/app/xcode/id497799835?l=en&mt=12)
+ Linux:
    - [Gcc](https://ftp.gnu.org/gnu/gcc/)
    - [make](https://ftp.gnu.org/gnu/make/)

### How do I get set up? ###

* [miARma can be installed using the following guide]


### Guidelines/How to ###

* [miRNAs guide]
* [mRNAs guide]
* [circRNAs guide]

### Code Documentation ####
* [Perldoc]

### Who do I talk to? ###

* eduardo.andres at csic.es
