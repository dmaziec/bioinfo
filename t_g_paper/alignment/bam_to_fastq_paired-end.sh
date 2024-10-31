#!/usr/bin/env bash

# Input arguments
input_file_bam=$1
nthreads=$2
output_file_prefix=$3

# Command line options
collate_opt="-f -r 10000000 -u -O"
# -r
# Fast mode keeps a buffer of alignments in memory so that it can write out most pairs as soon as they are found instead of storing them in temporary files.
# This allows collate to avoid some work and so finish more quickly compared to the standard mode.
# The number of alignments held can be changed using -r, storing more alignments uses more memory but increases the number of pairs that can be written early.
# Setting it to 10,000,000, this will use around 12GB RAM.

fastq_opt="-0 unset_reads.txt -s singletons.txt -n"

# Run
samtools collate -@ $nthreads $collate_opt $input_file_bam | \
samtools fastq -1 ${output_file_prefix}.1.fastq.gz -2 ${output_file_prefix}.2.fastq.gz -@ $nthreads $fastq_opt - || { echo "Cannot convert BAM to FASTQs"; exit 1; }