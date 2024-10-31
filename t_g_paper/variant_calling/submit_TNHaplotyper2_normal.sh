#!/bin/bash
#SBATCH -c 16                          # Cores
#SBATCH --mem 32GB
#SBATCH -t 5:00:00                    # Runtime in D-HH:MM format
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -o /n/scratch/users/d/dm334/log/TG_paper_BQSR_%j.out  # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /n/scratch/users/d/dm334/log/TG_paper_BQSR_%j.err  # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=END
#------

module load gcc
module load samtools

# Export
export SENTIEON_LICENSE=license.rc.hms.harvard.edu:8990
export SENTIEON_INSTALL_DIR=/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.01
export PATH=$PATH:/n/data1/hms/dbmi/park/SOFTWARE/Sentieon/sentieon-genomics-202308.01/bin

TUMOR_BAM=$1



NORMAL_BAM=$2
SM_TUMOR=$3
SM_NORMAL=$4
mkdir -p ${SM_TUMOR}_${SM_NORMAL}
cd ${SM_TUMOR}_${SM_NORMAL}

SM_NORMAL=$(samtools view -H $NORMAL_BAM | grep '@RG' | head -n 1 | sed 's/.*SM:\([^ \t]*\).*/\1/')
/n/data1/hms/dbmi/park/dominika/testing/tg_paper/scripts/TNHaplotyper_normal.sh \
-t $TUMOR_BAM -n $NORMAL_BAM -s $SM_TUMOR -d $SM_NORMAL \
-r /n/data1/hms/dbmi/park/DATA/crossplatform/bams/bams_from_cgap/Homo_sapiens_assembly38.fasta \
-g /n/data1/hms/dbmi/park/dominika/testing/smaht/tissues/scripts/ref_files/SMAFIF5Y6TL1.vcf.gz
#/home/dm334/park_home/testing/smaht/tissues/scripts/ref_files/SMAFI3F4TR9Q.vcf.gz
