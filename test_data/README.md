## Quick start / Minimal examples

This is a case example of a mouse transcriptomic dataset, analyzed in a machine with 12 cores and 32GB RAM available. The 4 samples are processed in parlalel and the genes Rpl4 and Krt14 are highlighted in various plots. 
The analyses will take ~30 min (reference genome indexing will take ~20 min).

Depending on the installation method, you can either:

1) Use the Apptainer image:

```
wget -q https://bit.ly/reana_apptainer_image -O reanalyzerGSE.sif # Download the Apptainer image
wget -q https://bit.ly/reana_test_data -O test_data.tar && tar xf test_data.tar --one-top-level && rm test_data.tar # Download test datasets

cd test_data
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz # Download reference genome
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.basic.annotation.gtf.gz # Download reference genome annotation
gzip -d *.gz # To uncompress the references

apptainer exec ../reanalyzerGSE.sif reanalyzerGSE.pk.sh -options manual_options_out_test.txt 2>&1 | tee -a manual_options_out_test.log
```

2) Use the conda-based installation
```
cd reanalyzerGSE # Cloned repo
source external_software/source_path.sh # To set up the PATH

cd test_data
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz # Download reference genome
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.basic.annotation.gtf.gz # Download reference genome annotation
gzip -d *.gz # To uncompress the references

reanalyzerGSE.pk.sh -options manual_options_out_test.txt 2>&1 | tee -a manual_options_out_test.log
```


