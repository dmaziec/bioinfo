#!/usr/bin/env bash

# *******************************************
# Perform base scores and indels recalibration to
# generate an analysis-ready BAM for a single sample.
# The input BAM file need to be sorted and pre-processed
# to add read groups, mark duplicates, and realign indels.
# *******************************************

## Command line arguments
# Input BAM
input_bam=$1

# Reference data files
reference_fa=$2
known_sites_snp=$3
known_sites_indel=$4

## Other settings
nt=$(nproc) # number of threads to use in computation,
            # set to number of cores in the server

## SymLink to reference files
fasta="reference.fasta"
rm reference.fasta* 

# FASTA reference
ln -s ${reference_fa} reference.fasta
ln -s ${reference_fa}.fai reference.fasta.fai
ln -s ${reference_fa}.dict reference.dict

# *****************************************************************************
# 1. Base recalibration - see:
# https://support.sentieon.com/appnotes/arguments/#bqsr-calculate-recalibration
# Not generating RECAL_DATA.TABLE.POST for plotting, just need recal_data.table.
# *****************************************************************************
sentieon driver -r $fasta -t $nt -i $input_bam --algo QualCal -k $known_sites_snp -k $known_sites_indel recal_data.table || exit 1
sentieon driver -r $fasta -t $nt -i $input_bam --read_filter 'QualCalFilter,table=recal_data.table,keep_oq=true' --algo ReadWriter recalibrated.bam || exit 1

# ******************************************
# 2. Check recalibrated BAM integrity.
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

check_EOF('recalibrated.bam')
"

python -c "$py_script" || exit 1