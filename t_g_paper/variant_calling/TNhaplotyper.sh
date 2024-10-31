#!/usr/bin/env bash

# *******************************************
# Script to run TNhaplotyper2 and TNfilter
#   on tumor only data
#
# Arguments:
#   -i input file for tumor in BAM format, index file required
#   -s sample name for tumor as string
#   -r genome reference file in FASTA format, index files required
#   -g population allele frequencies file in compressed (gzip|bgzip) VCF format, index file required
#   [-p panel of normal file in compressed (gzip|bgzip) VCF format, index file required]
#   [-o prefix for the output files, default is input file basename]
#   [-l interval to restrict calculation to (chr:start-end)]
#
# Output:
#   - <output_prefix>[_<interval>].vcf.gz: Output file with RAW variant calls in compressed (bgzip) VCF format
#     <output_prefix>[_<interval>].vcf.gz.tbi: Tabix index file
#   - <output_prefix>[_<interval>]_filtered.vcf.gz: Output file with FILTERED variant calls in compressed (bgzip) VCF format
#     <output_prefix>[_<interval>]_filtered.vcf.gz.tbi: Tabix index file
# *******************************************

export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.01
export PATH=$PATH:/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.01/bin
## Variables
USAGE="Usage: TNhaplotyper2_tumor_only.sh -i <input_file_bam> -s <sample_name> -r <genome_reference_fasta> -g <population_allele_frequencies> [-p panel_of_normal] [-o output_prefix] [-l interval]"

## Functions
check_args()
{
    arg_names=($@)
    for arg_name in ${arg_names[@]}; do
        [ -z ${!arg_name} ] && \
        echo "Missing Argument <${arg_name}>" 1>&2 && \
        echo $USAGE 1>&2 && \
        exit 1
    done
    return 0
}

## Bash command line definition
while getopts 'i:s:r:g:p:o:l:h' opt; do
  case $opt in
    # Required arguments
    i) input_file_bam=${OPTARG} ;;
    s) sample_name=${OPTARG} ;;
    r) genome_reference_fasta=${OPTARG} ;;
    g) population_allele_frequencies=${OPTARG} ;;
    # Optional arguments
    p) panel_of_normal=${OPTARG} ;;
    o) output_prefix=${OPTARG} ;;
    l) interval=${OPTARG} ;;
    ?|h)
      echo $USAGE 1>&2
      exit 1
      ;;
  esac
done
shift $(($OPTIND -1))

## Check arguments
check_args input_file_bam sample_name genome_reference_fasta population_allele_frequencies

## Other settings
nt=$(nproc) # Number of threads to use in computation,
            #   set to number of cores in the server

if ! [ -z ${output_prefix} ]; then
    output=${output_prefix}
  else
    output=$(basename ${input_file_bam} .bam)
fi

if ! [ -z ${interval} ]; then
  output+="_${interval}"
fi

# ******************************************
# 1. Create TNhaplotyper2 command line
# ******************************************
## Basic command
command="sentieon driver -t ${nt} -r ${genome_reference_fasta} -i ${input_file_bam}"

## Add interval if specified
if ! [ -z ${interval} ]; then
  command+=" --interval ${interval}"
fi

## Add TNhaplotyper2 base arguments
command+=" --algo TNhaplotyper2 --tumor_sample ${sample_name} --germline_vcf ${population_allele_frequencies}"

# Add TNhaplotyper2 optional arguments
if ! [ -z ${panel_of_normal} ]; then
  command+=" --pon ${panel_of_normal}"
fi

# Specify output file
command+=" ${output}.vcf.gz"

## Add OrientationBias arguments
command+=" --algo OrientationBias --tumor_sample ${sample_name} ${output}.priors"

## Add ContaminationModel arguments
command+=" --algo ContaminationModel --tumor_sample ${sample_name} -v ${population_allele_frequencies} --tumor_segments ${output}.segments ${output}.contamination"

# ******************************************
# 2. Run TNhaplotyper2
# ******************************************
eval $command || exit 1

# ******************************************
# 2. Run TNfilter
# ******************************************
sentieon driver -t ${nt} -r ${genome_reference_fasta} \
   --algo TNfilter \
   -v ${output}.vcf.gz \
   --tumor_sample ${sample_name} \
   --contamination ${output}.contamination \
   --tumor_segments ${output}.segments \
   --orientation_priors ${output}.priors \
   ${output}_filtered.vcf.gz || exit 1