import pandas
import scipy
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os
import click

def create_folder(folder_name):
    if os.path.exists(folder_name) != True:
        os.makedirs(folder_name)
dateTimeObj = datetime.now()
timestampStr = dateTimeObj.strftime("%d_%b_%Y_%H_%M_%S_%f")

output_folder = f"SVELT_clustering/{timestampStr}"
heatmap_folder=f"{output_folder}/heatmaps"
clusters_folder=f"{output_folder}/clusters"
dataframes_folder=f"{output_folder}/dataframes"
missing_samples_folder = f"{output_folder}/missing_samples"
create_folder(output_folder)
create_folder(heatmap_folder)
create_folder(clusters_folder)
create_folder(dataframes_folder)
create_folder(missing_samples_folder)

def create_reports(missing_samples_df, histology):
    report_file = open(f"{output_folder}/missing_samples_in_publication.txt", "a")
    missing_samples = open(f"{missing_samples_folder}/missing_samples.tsv", "a")
    # number of missing samples 
    number_of_samples_missing = len(missing_samples_df)

    # list of missing samples UUIDs
    list_of_missing_samples = ",".join(missing_samples_df["UUID"].tolist())

    
    # save it in the report file
    if number_of_samples_missing > 0:
        report_file.write(f"Histology: {histology}\n")
        report_file.write(f"Number of samples available in SVELT but not available in the publication: {number_of_samples_missing}\n")
        report_file.write(f"UUIDs of the missing samples: {list_of_missing_samples} \n")
        report_file.write("####################### \n")
    
    for missing_uuid in missing_samples_df["UUID"].tolist():
        missing_samples.write(f"{histology}\t{missing_uuid}\n")
    
    report_file.close()
    missing_samples.close()


@click.command()
@click.option("--svelt_samples", help="Svelt samples in CSV.")
@click.option("--pubication_samples", help="Publication samples in TSV. Source: https://www.nature.com/articles/s43018-020-0027-5 -- Supplementary Table -- Sheet 7")
def run_clustering(svelt_samples, pubication_samples):
    """ Run hierarchical clustering to organize PCAWG samples on SVELT."""
    #output folders
    output_folder = f"SVELT_clustering/{timestampStr}"
    heatmap_folder=f"{output_folder}/heatmaps"
    clusters_folder=f"{output_folder}/clusters"
    dataframes_folder=f"{output_folder}/dataframes"
    create_folder(output_folder)
    create_folder(heatmap_folder)
    create_folder(clusters_folder)
    create_folder(dataframes_folder)

    # load svelt samples
    svelt_samples = pandas.read_csv(
        svelt_samples,
        sep=",",
        index_col=False,
    )
    # load publication samples
    publication_data = pandas.read_csv(
        pubication_samples, sep="\t", index_col=False
    )

    # get cancer types 
    histology_abb = svelt_samples.histology_abbreviation.unique().tolist()

    #iterate over cancer types
    for his in histology_abb:
        print(f"processing {his}")

        # select svelt samples with the given cancer type
        svelt_cancer_subset = svelt_samples[svelt_samples["histology_abbreviation"] == his]

        # get samples that are missing in the publication
        svel_publication_missing = svelt_cancer_subset[
            ~svelt_cancer_subset["UUID"].isin(publication_data.UUID.tolist())]

        create_reports(svel_publication_missing, his)
        
        # get corresponding publication samples for the svelt subset samples 
        publication_data_subset = publication_data[
            publication_data["UUID"].isin(svelt_cancer_subset.UUID.tolist())
        ].drop("UUID", axis=1) # drop UUIDs as we want to run clustering on a matrix of numbers

        # run clustering only on a dataset with more than 2 samples 
        if len(publication_data_subset) > 2:

            # hierarchical algorithm and heat map
            cluster_map = sns.clustermap(
                publication_data_subset,
                #annot = True,
                cmap = sns.color_palette("crest", as_cmap=True), #green colors
                yticklabels=publication_data.loc[publication_data_subset.index]["UUID"].tolist() # uuids labels instead of indexes
            )

            # save the dataframe on which we run the algorithm
            publication_data[publication_data["UUID"].isin(svelt_cancer_subset.UUID.tolist())].to_csv(f"{dataframes_folder}/{his}.csv", index=False)

            plt.title(f'Cluster heatmap for {his}', loc='left', fontsize = 20)

            plt.figure(figsize=(6, 6), 
            dpi = 600) 
            cluster_map.savefig(f"{heatmap_folder}/{his}.pdf")
            plt.clf()

            # convert the results to a dendrodram so we can deduce the order of samples
            scipy_dendrogram = scipy.cluster.hierarchy.dendrogram(
                cluster_map.dendrogram_row.linkage,
                labels = publication_data_subset.index)

            # scipy_dendrogram["ivl"] contains the oreder of samples after clustering
            f = open(f"{clusters_folder}/{his}.txt", "a")
            for index in scipy_dendrogram["ivl"]:
                # write the UUIDs
                f.write(f"{publication_data.loc[index]['UUID']}\n")
            f.close()
        
        else:
            f = open(f"{clusters_folder}/not_clustered_datasets.txt", "a")
            f.write(f"{his}\n")
            f.close()

if __name__ == "__main__":
    run_clustering() 