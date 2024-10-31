#!/usr/bin/env bash

# *******************************************
# Mark duplicates in a sorted BAM file
# for a single sample.
# The sorted BAM file need to be pre-processed
# to add read groups.
# *******************************************

## Command line arguments
# Input BAM
input_bam=$1

# Other
optical_dup_pix_dist=$2
# The maximum offset between two duplicate clusters in order to consider them optical duplicates.
# This should be set to 100 for (circa 2011+) read names and typical flowcells.
# Structured flow cells (NovaSeq, HiSeq 4000, X) should use ~2500.
# For older conventions, distances could be to some fairly small number (e.g. 5-10 pixels).

## Other settings
nt=$(nproc) # number of threads to use in computation,
            # set to number of cores in the server

# ******************************************
# 1. Mark/remove duplicate reads.
# By ommiting the --rmdup option in Dedup
# we are only marking to match GATK Best Practices.
# ******************************************
sentieon driver -t $nt -i $input_bam --algo LocusCollector --fun score_info score.txt || exit 1
sentieon driver -t $nt -i $input_bam --algo Dedup --optical_dup_pix_dist $optical_dup_pix_dist --score_info score.txt --metrics dedup_metrics.txt deduped.bam || exit 1

# ******************************************
# 2. Check deduped BAM integrity.
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

check_EOF('deduped.bam')
"

python -c "$py_script" || exit 1