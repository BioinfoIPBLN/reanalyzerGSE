#!/bin/bash
file=$1
dest=$2
number_parallel=$3 
cores=$4
#compression_level=$5

total_lines=$(cat $1 | wc -l)
#echo -e "\nProcessing and downloading with fasterq-dump $total_lines files\n"
echo -e "\nProcessing and downloading with fastq-dl $total_lines files\n"

# Make sure cores get distributed:
if [ $number_parallel -le $cores ]; then
	cores_parallel=$((cores / number_parallel))
else
	number_parallel=$cores
	cores_parallel=$((cores / number_parallel))
fi

#echo -e "\nUsing parallel blocks of $number_parallel files with $cores_parallel cores each and $compression_level compression level\n"
echo -e "\nUsing parallel blocks of $number_parallel files with $cores_parallel cores\n"


# Processing:
cd $dest
# https://www.biostars.org/p/359441/
# parallel --verbose --joblog $dest/../download_sra_fq_log_parallel.txt -j $number_parallel download_sra_fq0.sh {} $cores_parallel $compression_level ::: $(cat $file)

parallel --verbose --joblog $dest/../download_sra_fq_log_parallel.txt -j $number_parallel fastq-dl --accession {} --cpus $cores_parallel ::: $(cat $file)

cut -f106 fastq-run-info.tsv | sed '1d' > $(dirname $file)/library_layout_info.txt

rm fastq-run*



