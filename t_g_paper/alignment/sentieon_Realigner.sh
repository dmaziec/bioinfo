#!/usr/bin/env bash

# *******************************************
# Realign indels in a sorted BAM file
# for a single sample.
# The sorted BAM file need to be pre-processed
# to add read groups and mark duplicates.
# *******************************************

## Command line arguments
# Input BAM
input_bam=$1

# Reference data files
reference_fa=$2
known_sites_indel=$3

## Other settings
nt=${4:-$(nproc)} # number of threads to use in computation,
                  # set to number of cores in the server if not specified

## SymLink to reference files
fasta="reference.fasta"



# ******************************************
# 1. Indel realignment
# ******************************************
sentieon driver -r $fasta -t $nt -i $input_bam --algo Realigner -k $known_sites_indel realigned.bam || exit 1

# ******************************************
# 2. Check realigned BAM integrity.
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

check_EOF('realigned.bam')
"

python -c "$py_script" || exit 1