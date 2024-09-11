#!/bin/bash
#SBATCH -c 16                          # Cores
#SBATCH --mem 64GB
#SBATCH -t 88:00:00                    # Runtime in D-HH:MM format
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -o /n/scratch/users/d/dm334/log/TG_paper_%j.out  # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /n/scratch/users/d/dm334/log/TG_paper_%j.err  # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=END
#------

# Check input arguments
if [ "$#" -ne 6 ]; then
  echo "Usage: $0 <input_bam> <output_folder> <sample> <technology> <force_rerun> <ref_genome>"
  exit 1
fi

input_bam=$1
output_folder=$2
sample=$3
technology=$4
force_rerun=$5
ref_genome=$6 
ERROR_LOG="/n/scratch/users/d/dm334/log/TG_paper_${SLURM_JOB_ID}.err"
STDOUT_LOG="/n/scratch/users/d/dm334/log/TG_paper_${SLURM_JOB_ID}.out"

export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.01
export PATH=$PATH:/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.01/bin

if [ "$ref_genome" = "hg19" ]; then

reference_genome="/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh37d5/human_g1k_v37_decoy.fasta"
known_sites_dbsnp="/n/data1/hms/dbmi/park/dominika/testing/tg_paper/resources/dbsnp_138.hg19.vcf.gz"
known_sites_indel="/n/data1/hms/dbmi/park/dominika/testing/tg_paper/resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
bwt="/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh37d5/human_g1k_v37_decoy.fasta"

fi

if [ "$ref_genome" = "hg38" ]; then
reference_genome="/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38_smaht/hg38_no_alt.fa"
known_sites_dbsnp=""
known_sites_indel="/n/data1/hms/dbmi/park/dominika/testing/smaht/reference_files/known_sites_indel/SMAFIPOL9T5R.vcf.gz"
bwt="/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/GRCh38_smaht/hg38_no_alt.fa"

fi
# Create output directory
#cd /n/scratch/users/d/dm334/tg_paper

mkdir -p "${output_folder}_${technology}"
cd "${output_folder}_${technology}"

rm -rf errorLogs.err
rm -rf stdoutLogs.out 
ln -s $ERROR_LOG errorLogs.err
ln -s $STDOUT_LOG stdoutLogs.out

echo "Sample: " $sample
echo "Technology: " $technology

# Check force_rerun flag for BAM to FASTQ
echo "Running BAM to FASTQ"
if [ "$force_rerun" = "true" ] || ! [ -f "${sample}.1.fastq.gz" ]; then
  /n/data1/hms/dbmi/park/dominika/testing/tg_paper/scripts/bam_to_fastq_paired-end.sh "${input_bam}" 16 "${sample}" || exit 1
fi

# Run Sentieon BWA mem
echo "Running Sentieon BWA mem"
if [ "$force_rerun" = "true" ] || ! [ -f sorted.bam ]; then
  /n/data1/hms/dbmi/park/dominika/testing/tg_paper/scripts/sentieon_bwa-mem_sort.sh "${sample}.1.fastq.gz" "${sample}.2.fastq.gz" $reference_genome $bwt || \
  { echo "Error in Sentieon BWA mem" && rm -f sorted.bam && exit 1; }
fi

# Run AddReadGroups
echo "Running AddReadGroups"
if [ "$force_rerun" = "true" ] || ! [ -f sorted_rg.bam ]; then
  python3 /n/data1/hms/dbmi/park/dominika/testing/tg_paper/scripts/AddReadGroups.py -i sorted.bam -s $sample -x -t 16 || \
  { echo "Error in AddReadGroups" && rm -f sorted_rg.bam && exit 1; }
fi

# Run Sentieon Dedup
echo "Running Sentieon Dedup"
if [ "$force_rerun" = "true" ] || ! [ -f deduped.bam ]; then
  /n/data1/hms/dbmi/park/dominika/testing/tg_paper/scripts/sentieon_Dedup.sh sorted_rg.bam 2500 || \
  { echo "Error in Sentieon Dedup" && rm -f deduped.bam && exit 1; }
fi

echo "END"