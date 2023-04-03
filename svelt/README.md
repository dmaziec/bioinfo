## PCAWG samples 

In order to reproduce the results if the hierarchical clustering on the PCAWG samples available on SVELT,  please run `clustering.py`.  The aim of the tool is to organize the PCAWG samples uploaded to the browser based on the patterns of rearrangement signatures to facilitate the process of structural variants interpretation. 

## CLI 

The help command of the tool: 

```
python3 clustring.py  --help 
Usage: clustring.py [OPTIONS]

  Run hierarchical clustering to organize PCAWG samples on SVELT.

Options:
  --svelt_samples TEXT       Svelt samples in CSV.
  --pubication_samples TEXT  Publication samples in TSV. Source:
                             [https://www.nature.com/articles/s43018-020-0027-5](https://www.nature.com/articles/s43018-020-0027-5)
                             -- Supplementary Table -- Sheet 7
  --help                     Show this message and exit.
``` 





## Input files

|       Parameter|Description                 
|----------------|-------------------------------|
|`svelt_samples` | List of samples used as demo data for SVELT. Cancer type is defined in the `histology_abbreviation` column. Sample UUIDs are defined in the `aliquot_id` column.      |
|`publication_samples`          |  Patters of somatic rearrangements across cancers (Sheet 7. Columns correspond to different patterns of rearrangement signatures. Rows correspond to sample IDs.         |

## Output files and directories

|Output| Description|
|----|----|
|`clusters/`|directory containing lists of ordered samples. Each lists represents one cancer type from the dataset. |
|`clusters/non_clustered_datasets.txt` | file containing cancer types for which SVELT has less than 3 samples. These cancer types were not subject to the clustering algorithm. They are displayed in a default order in the cohort overview. 
|`dataframes/`| directory of TSV files representing matrices used for clustering.|
|`heatmaps/` | directory of generated heat maps in the PDF format for visual inspection on how the algorithm grouped the samples |
|`missing_samples/` | directory containing `missing_samples.tsv` that stores samples demonstrated on SVELT but not included in the analysis due to missing rearrangement types for these samples. 
|`missing_samples_in_publication.txt` | it’s a report containing a compacted list of sample UUIDs excluded from the analysis and their total number for each cancer type | 




