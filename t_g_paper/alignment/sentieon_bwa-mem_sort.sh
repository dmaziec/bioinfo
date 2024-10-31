#!/usr/bin/env bash

# *******************************************
# Generate an alignment BAM file from paired FASTQ
# files for a single sample.
# The BAM file will be sorted by coordinates.
# *******************************************

## Command line arguments
# Input FASTQ files
fastq_r1=$1
fastq_r2=$2

# Reference data files
reference_fa=$3
reference_bwt=$4

## Other settings
nt=$(nproc) # number of threads to use in computation,
            # set to number of cores in the server

## SymLink to reference files
fasta="reference.fasta"


# BWT reference
ln -s ${reference_bwt}.bwt reference.fasta.bwt
ln -s ${reference_bwt}.ann reference.fasta.ann
ln -s ${reference_bwt}.amb reference.fasta.amb
ln -s ${reference_bwt}.pac reference.fasta.pac
ln -s ${reference_bwt}.sa reference.fasta.sa

if [ -f ${reference_bwt}.alt ]; then
    ln -s ${reference_bwt}.alt reference.fasta.alt
fi

# ******************************************
# 1. Mapping reads with BWA-MEM and
# sort by coordinates.
# The results of this call are dependent on
# the number of threads used.
# To have number of threads independent results,
# add chunk size option -K 10000000.
# ******************************************
( sentieon bwa mem -t $nt -K 10000000 $fasta $fastq_r1 $fastq_r2 || exit 1 ) | samtools sort --no-PG -@ $nt -o sorted.bam - || exit 1

# ******************************************
# 2. Index BAM
# ******************************************
samtools index -@ $nt sorted.bam || exit 1

# ******************************************
# 3. Check BAM integrity.
# ******************************************
py_script="
import sys, os

def check_EOF(filename):
    EOF_hex = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    size = os.path.getsize(filename)
    fb = open(filename, 'rb')
    fb.seek(size - 28)
    EOF = fb.read(28)
    fb.close()
    if EOF != EOF_hex:
        sys.stderr.write('EOF is missing\n')
        sys.exit(1)
    else:
        sys.stderr.write('EOF is present\n')

check_EOF('sorted.bam')
"

python -c "$py_script" || exit 1