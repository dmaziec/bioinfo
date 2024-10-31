#!/bin/bash
#SBATCH -c 16                          # Cores
#SBATCH --mem 64GB
#SBATCH -t 5:00:00                    # Runtime in D-HH:MM format
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -o /n/scratch/users/d/dm334/log/TG_paper_BQSR_%j.out  # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /n/scratch/users/d/dm334/log/TG_paper_BQSR_%j.err  # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=END
#------

# Check input arguments
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <input_bam> <output_folder> <sample> <technology>"
  exit 1
fi

input_bam=$1
output_folder=$2
sample=$3
technology=$4

ERROR_LOG="/n/scratch/users/d/dm334/log/TG_paper_BQSR_${SLURM_JOB_ID}.err"
STDOUT_LOG="/n/scratch/users/d/dm334/log/TG_paper_BQSR_${SLURM_JOB_ID}.out"

export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.01
export PATH=$PATH:/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.01/bin

reference_genome="/n/data1/hms/dbmi/park/DATA/crossplatform/bams/bams_from_cgap/Homo_sapiens_assembly38.fasta"
known_sites_dbsnp="/n/data1/hms/dbmi/park/dominika/testing/cgap/reference_files/known_sites_snps/GAPFI4LJRN98.vcf.gz"
known_sites_indel="/n/data1/hms/dbmi/park/dominika/testing/cgap/reference_files/known_sites_indel/GAPFIAX2PPYB.vcf.gz"
bwt="/n/data1/hms/dbmi/park/DATA/crossplatform/bams/bams_from_cgap/Homo_sapiens_assembly38.fasta"

# Create output directory
cd /n/scratch/users/d/dm334/tg_paper_hg38
mkdir -p "${output_folder}_${technology}"
cd "${output_folder}_${technology}"

pwd

rm -rf errorLogsBQSR.err
rm -rf stdoutLogsBQSR.out 
ln -s $ERROR_LOG errorLogsBQSR.err
ln -s $STDOUT_LOG stdoutLogsBQSR.out

rm recalibrated.bam
echo "Running Sentieon QualCal"
if ! [ -f recalibrated.bam ]; then
  /n/data1/hms/dbmi/park/dominika/testing/tg_paper/scripts/sentieon_QualCal.sh $input_bam $reference_genome $known_sites_dbsnp $known_sites_indel || \
  { echo "Error in Sentieon QualCal" && rm -f recalibrated.bam && exit 1; }
fi