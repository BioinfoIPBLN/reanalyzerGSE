## Quick start / Minimal examples
```
cd reanalyzerGSE # Installation folder
source external_software/source_path # To set up the PATH

cd test_data
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz # Download reference genome
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.basic.annotation.gtf.gz # Download reference genome annotation
wget -q https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz # Download reference genome transcripts
gzip -d *.gz # To uncompress the references

### Case example of a mouse transcriptomic dataset, analyzed in a machine with 12 cores and 32GB RAM available.
# Local raw data simultaneously processing 4 samples and highligthing the genes Gpatch3 and Tent2
reanalyzerGSE.pk.sh -options manual_options_out_test.txt 2>&1 | tee -a manual_options_out_test.log


```
