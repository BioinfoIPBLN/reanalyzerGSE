## Quick start / Minimal examples
```
source reanalyzerGSE/external_software/source_path # To set up the PATH if you have followed option 1 for installation above
cores=30
cd reanalyzerGSE/test_data
wget -q https://bit.ly/case_examples; unzip case_examples; mv test_data/* .; rm -r test_data
pigz -p $cores -d *.fa.gz *.gtf.gz # To uncompress the reference genome and annotation
cd ../
cd ../
### Case examples of mouse transcriptomic datasets, analyzed in a machine with 30 cores and 200GB RAM available.
# A) Local raw data (subset, 500k reads) simultaneously processing 4 samples and highligthing the genes B4galt3 and Chd1
reanalyzerGSE.pk.sh -i $PWD/test_data -r $PWD/test_data/GRCm39.primary_assembly.genome.fa -a $PWD/test_data/gencode.vM28.annotation.gtf -s reverse -o $PWD/test_data_out/test_data_1/ -p $cores -P 4 -g B4galt3,Chd1 -M 214748364800 2>&1 | tee -a $PWD/test_data_out/test_data_1.log

# B) The GEO entry GSE118451, simultaneously processing 6 samples, highligthing the genes B4galt3 and Chd1, and predicting strandness from the transcript sequences
reanalyzerGSE.pk.sh -i GSE118451 -r $PWD/test_data/GRCm39.primary_assembly.genome.fa -a $PWD/test_data/gencode.vM28.annotation.gtf -t $PWD/test_data/Mus_musculus.GRCm39.cdna.all.fa -o $PWD/test_data_out/test_data_2/ -p $cores -P 6 -g B4galt3,Chd1 -M 214748364800 2>&1 | tee -a $PWD/test_data_out/test_data_2.log
```
