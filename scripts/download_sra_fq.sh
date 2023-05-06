#!/bin/bash
file=$1
dest=$2
number_parallel=$3 
cores=$4

total_lines=$(cat $1 | wc -l)
echo -e "\nProcessing and downloading with fasterq-dump $total_lines files\n"

# Make sure cores get distributed:
if [ $number_parallel -le $cores ]; then
	cores_parallel=$((cores / number_parallel))
else
	number_parallel=$cores
	cores_parallel=$((cores / number_parallel))
fi

echo -e "\nUsing parallel blocks of $number_parallel files with $cores_parallel cores each\n"


# Processing:
cd $dest
# https://www.biostars.org/p/359441/
parallel --verbose -j $number_parallel download_sra_fq0.sh {} $cores_parallel ::: $(cat $file)



