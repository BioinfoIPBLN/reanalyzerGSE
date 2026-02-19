## Quick start / Minimal examples

## Analysis of local data:
This is a case example of a mouse transcriptomic dataset, analyzed in a machine with 12 cores and 32GB RAM available. The 4 samples are processed in parlalel and the genes Rpl4 and Krt14 are highlighted in various plots. 
It may depend on the machine, but the analyses should take ~30 min (reference genome indexing will take ~20 min).

Depending on the installation method, you can either:

1) Use the Apptainer image:

```
# wget -q https://bit.ly/reana_apptainer -O reanalyzerGSE.sif # Download the Apptainer image or use the one you created (see Installation instructions in main README)
wget -q https://bit.ly/reana_test -O test_data.tar && tar xf test_data.tar && rm test_data.tar # Download test datasets

wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz # Download reference genome
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.basic.annotation.gtf.gz # Download reference genome annotation
gzip -d *.gz # To uncompress the references

apptainer exec reanalyzerGSE.sif reanalyzerGSE.sh -options options_test.yaml 2>&1 | tee -a out_test.log
```

2) Use the conda-based installation
```
cd reanalyzerGSE # Cloned repo
source external_software/source_path.sh # To set up the PATH

cd test_data
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz # Download reference genome
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.basic.annotation.gtf.gz # Download reference genome annotation
gzip -d *.gz # To uncompress the references

reanalyzerGSE.sh -options options_test.yaml 2>&1 | tee -a out_test.log
```

## Reanalysis of a GEO dataset:
This is a case example of a the mouse transcriptomic dataset GSE118451, analyzed in a machine with 12 cores and 32GB RAM available. The 6 samples are processed in parlalel and the genes Rpl4 and Krt14 are highlighted in various plots. 

Depending on the installation method, you can either:

1) Use the Apptainer image:

```
# wget -q https://bit.ly/reana_apptainer -O reanalyzerGSE.sif # Download the Apptainer image or use the one you created (see Installation instructions in main README)

wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz # Download reference genome
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.basic.annotation.gtf.gz # Download reference genome annotation
gzip -d *.gz # To uncompress the references
wget https://github.com/BioinfoIPBLN/reanalyzerGSE/raw/refs/heads/main/test_data/manual_options_out_test_GEO.txt

apptainer exec reanalyzerGSE.sif reanalyzerGSE.sh -options options_test_GEO.yaml 2>&1 | tee -a out_test_GEO.log
```

2) Use the conda-based installation
```
cd reanalyzerGSE # Cloned repo
source external_software/source_path.sh # To set up the PATH

cd test_data
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz # Download reference genome
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.basic.annotation.gtf.gz # Download reference genome annotation
gzip -d *.gz # To uncompress the references

reanalyzerGSE.sh -options options_test_GEO.yaml 2>&1 | tee -a out_test_GEO.log
```

