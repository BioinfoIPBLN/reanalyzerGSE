#!/bin/bash
file=$1
cores_parallel=$2
compression_level=$3

fasterq-dump -e $(( cores_parallel *2 )) $file && pigz -$compression_level -p $(( cores_parallel *2 )) $file*
# I added a multiplicator, because nor fasterq-dump neither pigz are very-intensive, but overall faster if pigz ends and it's faster than fasterq-dump
