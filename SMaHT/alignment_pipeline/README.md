## ALIGNMENT.sh 

This is an alignment pipeline for Illumina short reads.
This script implements the following steps:

1. Alignment (bwa-mem)
2. Adding Read Groups
3. Marking duplicates


Input parameters:
- `input_bam` - path to the input BAM
- `output_folder` - output file path where to store the results
- `sample` - sample name
- `technology` - technology (e.g. HiSeq, NovaSeq)
- `force_rerun` - true if you need to overwrite the all the results, false if you need to continue with your analysis that has failed due to timeout/oom 
- `ref_genome` - hg19 or hg38 

## Example run: 
``` sbatch  /n/data1/hms/dbmi/park/dominika/testing/tg_paper/scripts/ALIGNMENT.sh /n/data1/hms/dbmi/park/DATA/crossplatform/bams/bams_from_giab/na12878_novaseq.bam  na12878 na12878 NovaSeq hg19  ```




