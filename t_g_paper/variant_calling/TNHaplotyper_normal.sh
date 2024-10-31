#!/usr/bin/env bash

# *******************************************
# Script to run TNhaplotyper2 and TNfilter
#   on tumor/normal paired data
#
# Arguments:
#   -t input file for tumor in BAM format, index file required
#   -n input file for normal in BAM format, index file required
#   -s sample name for tumor as string
#   -r genome reference file in FASTA format, index files required
#   -g population allele frequencies file in compressed (gzip|bgzip) VCF format, index file required
#   [-p panel of normal file in compressed (gzip|bgzip) VCF format, index file required]
#   [-o prefix for the output files, default is tumor input file basename]
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

USAGE="Usage: TNhaplotyper2_tumor_normal.sh -t <tumor_file_bam> -n <normal_file_bam> -s <sample_name> -r <genome_reference_fasta> -g <population_allele_frequencies> [-p panel_of_normal] [-o output_prefix] [-l interval]"

generate_random_string() {
    echo "$(cat /dev/urandom | tr -dc 'a-zA-Z' | fold -w 5 | head -n 1)"
}

process_bam_header() {
    # Variable to store replacement @RG args
    replace_rg_args=""

    # Extract @RG lines from the header of the input BAM
    input_bam="$1"
    rg_lines=$(samtools view -H "$input_bam" | grep "^@RG")

    # Loop through @RG lines and modify ID field
    # Create --replace_rg arguments
    while IFS= read -r rg_line; do
        orig_rg_id=$(echo "$rg_line" | awk -F'\t' '/ID:/ {print $2}' | sed "s/ID://")
        new_rg_id="${orig_rg_id}-$(generate_random_string)"
        new_rg_line=$(echo "$rg_line" | cut -f 2- | sed "s/${orig_rg_id}/${new_rg_id}/")
        formatted_new_rg_line=$(echo -e "$new_rg_line" | sed 's/\t/\\t/g')
        replace_rg_args+=" --replace_rg ${orig_rg_id}='${formatted_new_rg_line}' "
    done <<< "$rg_lines"

    # Return replacement @RG args
    echo "$replace_rg_args"
}

check_args() {
    arg_names=($@)
    for arg_name in ${arg_names[@]}; do
        [ -z ${!arg_name} ] && \
        echo "Missing Argument <${arg_name}>" 1>&2 && \
        echo $USAGE 1>&2 && \
        exit 1
    done
    return 0
}

while getopts 't:n:s:d:r:g:p:o:l:h' opt; do
  case $opt in
    # Required arguments
    t) tumor_file_bam=${OPTARG} ;;
    n) normal_file_bam=${OPTARG} ;;
    s) sample_name=${OPTARG} ;;
    d) normal_name=${OPTARG} ;;
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
check_args tumor_file_bam normal_file_bam sample_name normal_name genome_reference_fasta population_allele_frequencies

## Other settings
nt=$(nproc) # Number of threads to use in computation
if ! [ -z ${output_prefix} ]; then
    output=${output_prefix}
else
    output=$(basename ${tumor_file_bam} .bam)
fi

if ! [ -z ${interval} ]; then
  output+="_${interval}"
fi



input_files=""
replace_args=$(process_bam_header $tumor_file_bam)
input_files+=" $replace_args -i $tumor_file_bam "

replace_args=$(process_bam_header $normal_file_bam)
input_files+=" $replace_args -i $normal_file_bam "

# ******************************************
# 1. Create TNhaplotyper2 command line
# ******************************************
command="sentieon driver -t ${nt} -r ${genome_reference_fasta} ${input_files}"



## Add interval if specified
if ! [ -z ${interval} ]; then
  command+=" --interval ${interval}"
fi

## Add TNhaplotyper2 base arguments
command+=" --algo TNhaplotyper2 --tumor_sample ${sample_name} --normal_sample ${normal_name} --germline_vcf ${population_allele_frequencies}"

# Add optional panel of normal if provided
if ! [ -z ${panel_of_normal} ]; then
  command+=" --pon ${panel_of_normal}"
fi

# Specify output file
command+=" ${output}.vcf.gz"

## Add OrientationBias arguments
command+=" --algo OrientationBias --tumor_sample ${sample_name} ${output}.priors"

## Add ContaminationModel arguments
command+=" --algo ContaminationModel --tumor_sample ${sample_name} --normal_sample ${normal_name} -v ${population_allele_frequencies} --tumor_segments ${output}.segments ${output}.contamination"

# ******************************************
# 2. Run TNhaplotyper2
# ******************************************
eval $command || exit 1

# ******************************************
# 3. Run TNfilter
# ******************************************
sentieon driver -t ${nt} -r ${genome_reference_fasta} \
   --algo TNfilter \
   -v ${output}.vcf.gz \
   --tumor_sample ${sample_name} \
   --normal_sample ${normal_name} \
   --contamination ${output}.contamination \
   --tumor_segments ${output}.segments \
   --orientation_priors ${output}.priors \
   ${output}_filtered.vcf.gz || exit 1
