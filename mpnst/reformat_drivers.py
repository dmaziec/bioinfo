#!/usr/bin/env python3

# Authors: Dominika Maziec & Victoria Stevens
# Usage: python3 reformat_drivers.py -m [driver_tsv_filepath] -g [gene_annotation_filepath]

################################################
#   Libraries
################################################
import pandas
import json
import csv
import os
import re
from datetime import datetime
import argparse

################################################
#   Constants
################################################
FRAMESHIFT = "frameshift"
FRAMESHIFT_INSERTION = "frameshift insertion"
FRAMESHIFT_DELETION = "frameshift deletion"
HOMOZYGOUS_DELETION = "homozygous deletion"
DEEP_INTRONIC = "deep intronic"
STOPGAIN = "stopgain"
SPLICING = "splicing"
INDEL = "INDEL"
SV = "SV"
SNV = "SNV"
NONSYNONYMOUS_SNV = "nonsynonymous SNV"

# used for setting the "category" attribute, if applicable
MAPPING = {
    SV: SV,
    SNV: "mutation", # TODO: why not just SNV: "SNV" (same as SNV:SNV)
    NONSYNONYMOUS_SNV: SNV,
    FRAMESHIFT: INDEL,
    FRAMESHIFT_INSERTION: INDEL,
    FRAMESHIFT_DELETION: INDEL,
    STOPGAIN: "SNV or INDEL",
    SPLICING: "SNV or INDEL",
    DEEP_INTRONIC: "SNV or INDEL",  # TODO: NOT SURE
}

# list of all possible columns/attributes within a final, reformatted driver file
# GosCan documentation: https://sehilyi.github.io/goscan/docs/#/datahub?id=drivers-tsv
DRIVER_FILE_ATTRS = [
    "chr", # Required. Can be obtained by looking up gene positions if not explicitly available
    "pos", # Required. Here, we use midpoint, average between start and end positions
    "ref",
    "alt",
    "gene", # Required.
    "category", # options: (germline/somatic) mutation, INDEL, SV, CNV, SNV
    "top_category",
    "transcript_consequence",
    "protein_mutation",
    "allele_fraction",
    "mutation_type", # options: stopgain, nonsynonymous SNV, frameshift deletion, can stay empty for SVs.
                     # If (gene)_somatic_total_cn==0, homozygous deletion
    "biallelic" # Yes, if mutation is concurrent with loss-of-heterozygosity (LOH), when (gene)_somatic_minor_cn==0
]


################################################
#   Helper Functions
################################################
def create_simplified_genes_dictionary(gene_annotation_filepath, mpnst_gene_set):
    """
    Given the filepath of a gene annotation JSON file,
    create a filtered dictionary of dictionaries,
    where each dictionary corresponds to a singular gene and its metadata.
    Here, the file is filtered down to the genes within the
    MPNST consortium set of drivers.

    :param gene_annotation_filepath: filepath of a gene annotation JSON file
    :type gene_annotation_filepath: str
    :param mpnst_gene_set: set of MPNST consortium driver mutation gene names
    :type mpnst_gene_set: set(str)
    :return: dictionary of MPNST genes and metadata
    :rtype: dict of dicts
    """
    genes_file = open(gene_annotation_filepath)
    genes = json.load(genes_file)

    # filtered dict of dicts of possible genes, created from gene annotation file
    genes_dict = {
        gene["gene_symbol"]: {
            "ens_id": gene["ensgid"],
            "chr": gene["chrom"],
            "start": gene["spos"],
            "end": gene["epos"],
            "strand": gene["strand"],
            "mid_point": int(gene["spos"])
            + round((int(gene["epos"]) - int(gene["spos"])) / 2),
        }
        for gene in (genes)
        if gene["gene_symbol"] in mpnst_gene_set
    }
    return genes_dict


def parse_protein_mutation_metadata(proteins):
    """
    Given a string of driver mutations and metadata,
    parse out the transcript consequences and protein mutations,
    where a transcript consequence is of the form "c.{transcript_change}"
    and protein mutation is of the form "p.{protein_mutation}".
    Also handles cases where several mutations are listed in one string,
    comma-separated. Returns the transcript_consquence(s) and
    protein_mutation(s) as strings.

    :param proteins: string of mutation(s) metadata
    :type proteins: str
    :return: transcript consequence(s)
    :rtype: str
    :return: protein mutation(s)
    :rtype: str

    Example:
    parse_protein_mutation_metadata("NF1:NM_000000:exon3:c.2640C<G:p.Asp1255Asn") --> "c.2640C<G", "p.Asp1255Asn"
    """

    transcript_consequences = []
    protein_mutations = []

    def parse_protein_single_string(protein):
        components = protein.split(":")
        for component in components:
            if "c." in component:
                transcript_consequences.append(component)
            elif "p." in component:
                protein_mutations.append(component)

    # if there are several mutations separated by comma
    if "," in proteins:
        single_proteins = proteins.split(",")
        for protein in single_proteins:
            parse_protein_single_string(protein)
    # several proteins separated by semicolon
    elif ";" in proteins:
        single_proteins = proteins.split(";")
        for protein in single_proteins:
            parse_protein_single_string(protein)
    # or just a singular mutation
    else:
        parse_protein_single_string(proteins)

    # get rid of duplicate transcript consequences and protein mutations
    # while still keeping the same order of the lists
    non_dup_transcript_consequences = list(dict.fromkeys(transcript_consequences))
    non_dup_protein_mutations = list(dict.fromkeys(protein_mutations))

    # convert these lists into singular string, comma-separated
    all_transcript_consequences = ", ".join(non_dup_transcript_consequences)
    all_protein_mutations = ", ".join(non_dup_protein_mutations)

    return all_transcript_consequences, all_protein_mutations


def parse_germline_details_NF1(record, gene, category, mutation_type):
    """
    Specific to NF1 germline mutations from MPNST consortium,
    given an unparsed record of the NF1 driver mutation metadata,
    returns a dictionary of the parsed attributes in 
    GosCan-appropriate format.

    :param record: NF1 gene metadata string from MPNST drivers file
    :type record: str
    :param gene: name of gene (in this case, always "NF1")
    :type gene: str
    :param category: "category" attribute of this NF1 germline mutation (germline, etc.)
    :type category: str
    :param mutation_type: "mutation_type" attribute of this NF1 germline mutation
    :type mutation_type: str
    :return: a list of edited results dictionaries of the NF1 driver(s)
    :rtype: list[dict] or None
    """

    def get_chrom_start_end(token):
        """
        Given a token string of NF1 driver mutation metadata,
        parse out the chromosome, start, and end positions,
        and ref and alt allele attributes
        as a string, else returns None.

        :param token: token string of NF1 driver mutation metadata
        :type token: str
        :return: string containing chr, start, end, ref and alt
        :rtype: str or None

        Example:
        get_chrom_start_end("indel chr18 4000000 4000000 T G NF1 NF1:NM_000000...") --> "chr18 4000000 4000000 T G"
        """
        pattern = re.search("chr\d* \d* \d* [ACGT-]* [ACGT-]*", token)
        if pattern:
            return pattern.group()

    def get_gene_protein_str(token, gene):
        """
        Given a token string of (NF1) driver mutation metadata and a gene name (NF1)
        return a re.Match object with the indices of the first occurence
        of the gene name within the token.

        :param token: token string of NF1 driver mutation metadata
        :type token: str
        :param gene: gene name
        :type gene: str
        :return: re.Match object of the match with this gene within token string
        :rtype: re.Match object

        Example:
        get_gene_metadata_str("indel chr18 4000000 4000000 T G NF1 NF1:NM_000000...", "NF1") --> <re.Match object; span=(38, 41), match='NF1'>
        """
        return re.search(f"{gene}", token)

    # final list of dicts that will be returned at the end of fxn
    all_results = []

    # the NF1 germline record may list several drivers, split by ";;"
    records = record.split(";;")

    for r in records:

        # extract chromosome, start/end position, ref and alt allele attributes
        chrom_start_end = get_chrom_start_end(r)
        if chrom_start_end == None:
            return None
        chrom, start, end, ref, alt = chrom_start_end.split(" ")

        results = {}

        results["chr"] = chrom
        results["pos"] = start
        results["ref"] = ref
        results["alt"] = alt
        
        results["gene"] = gene
        results["category"] = category
        results["mutation_type"] = mutation_type
        
        gene_string = get_gene_protein_str(r, gene)
        gene_end = gene_string.end()
        # e.g. if r = "indel chr18 4000000 4000000 T G NF1 NF1:NM_000000.."
        # then proteins = "NF1:NM_000000.."
        # which contains metadata that is then parsed in parse_protein_mutation_metadata()
        proteins = r[gene_end + 1 :].split(" ")[0]
        transcript_consequence, protein_mutation = parse_protein_mutation_metadata(proteins)
        results["transcript_consequence"] = transcript_consequence
        results["protein_mutation"] = protein_mutation
        all_results += [results]

    return all_results


def fill_driver_file_attrs(gene, genes_dictionary, results_dict, category=None, mutation_type=None):
    """
    Given a results dictionary of a gene's driver, fills out its required attributes
    (gene, chr, and pos) as well as category and mutation_type, if applicable.

    :param gene: gene name
    :type gene: str
    :param genes_dictionary: gene annotations dictionary, from CGAP
    :type genes_dictionary: dict of dicts
    :param results_dict: dictionary of given gene's driver, to be filled out
    :type results_dict: dict
    :param category: category attribute for the given gene
    :type category: str or None
    :param mutation_type: mutation_type attribute for the given gene
    :type mutation_type: str or None
    :return: the edited results dictionary of the given gene
    :rtype: dict
    """
    
    # fill out required attrs: gene, chr, pos
    results_dict["gene"] = gene
    results_dict["chr"] = f'chr{genes_dictionary[gene]["chr"]}'
    results_dict["pos"] = genes_dictionary[gene]["mid_point"]
    
    # fill out category and mutation_type, if applicable
    if category is not None:
        results_dict["category"] = category
    if mutation_type is not None:
        results_dict["mutation_type"] = mutation_type
        
    return results_dict


################################################
#
#   reformat_drivers
#
################################################
def reformat_drivers(args):
    """
    Creates sample-wise driver files that are compatible for the GosCan
    visual analytics tool for cancer genomes (https://github.com/sehilyi/goscan),
    given a TSV of sample-wise driver mutations with metadata, and a
    corresponding gene annotation file. Driver file data format that is compatible
    with GosCan can be found at https://sehilyi.github.io/goscan/docs/#/datahub?id=drivers-tsv.
    The newly created driver files will be stored in a timestamped directory
    under a folder "./MPNST_drivers/", created in the same directory in which this script is run.

    Formatting specific to malignant peripheral nerve sheath tumour (MPNST) drivers file
    from the international Genomics of MPNST (GeM) consortium.
    Reference: https://www.biorxiv.org/content/10.1101/2022.05.03.490481v3.full.pdf

    :param mpnst: drivers TSV filepath
    :type mpnst: str
    :param genes: gene annotation JSON filepath (should be decompressed for the script)
    :type genes: str

    Example: 

        $ python3 reformat_drivers.py -m [driver_tsv_filepath] -g [gene_annotation_filepath]

        or

        $ python3 reformat_drivers.py --mpnst [driver_tsv_filepath] --genes [gene_annotation_filepath]

    For CGAP-specific case reformatting MPNST drivers, files in same directory as module:

        $ python3 reformat_drivers.py -m MPNST_drivers.tsv -g gene_inserts_v0.4.6.updated.json
    """

    # timestamp for creating new directory with newly created driver files
    dateTimeObj = datetime.now()
    timestampStr = dateTimeObj.strftime("%d_%b_%Y_%H_%M_%S_%f")

    # create timestamped output folder
    output_folder = f"MPNST_drivers/{timestampStr}"
    if os.path.exists(output_folder) != True:
        os.makedirs(output_folder)

    # extract relevant attributes from raw drivers file
    mpnst_tsv_filtered = pandas.read_csv(
        args["mpnst"], sep="\t", index_col=False
    ).filter(
        regex="(^GeM_ID$)|(^sample$)|(^(germline)$)|(.*_somatic$)|(.*_germline$)|(.*_germline_detailed$)|(.*_somatic_minor_cn$)|(.*_somatic_total_cn$)",
        axis=1,
    )

    # drop columns for following genes: CXorf67, COR, AZF1, ROS, FAM22, ZC3H7, INI1
    # These were checked manually -- aren't available within CGAP-specific gene annotation JSON file
    mpnst_tsv_filtered.drop(
        columns=list(mpnst_tsv_filtered.filter(regex = '(^CXorf67_somatic)|(^COR_somatic)|(^AZF1_somatic)|(^ROS_somatic)|(^FAM22_somatic)|(^ZC3H7_somatic)|(^INI1_somatic)')), 
        inplace = True
    )

    # gene names that are somatic mutations
    mpnst_genes_somatic = set(
        [
            col.split("_")[0]
            for col in list(mpnst_tsv_filtered.columns)
            if (col.endswith("_somatic"))
        ]
    )
    # gene names that are germline mutations
    mpnst_genes_germline = set(
        [
            col.split("_")[0]
            for col in list(mpnst_tsv_filtered.columns)
            if (col.endswith("_germline"))
        ]
    )

    # combine somatic and germline MPNST gene names
    # creates an easy-to-access set of MPNST gene symbols
    mpnst_genes = set(list(mpnst_genes_somatic) + list(mpnst_genes_germline))

    # extract column names referring to minor copy number
    col_mino_total = [ # getting minor copy number to fill out cells without value
        col
        for col in mpnst_tsv_filtered.columns
        if (col.endswith("_somatic_minor_cn") or col.endswith("_somatic_total_cn"))
    ]
    # and fill out NaN cells with negative number
    # negative number because copy number can't be negative
    # but filling with integer rather than NA in order to make 
    # integer comparison later on in this code possible
    mpnst_tsv_filtered[col_mino_total] = mpnst_tsv_filtered[col_mino_total].fillna(
        value=-1000
    )

    # convert pandas dataframe to Python dictionary
    mpnst_tsv_dict = mpnst_tsv_filtered.to_dict(orient="records")

    # use gene annotation JSON file to create dictionary of genes and corresponding metadata on CGAP
    # genes_dict is a dictionary of dictionaries of the form:
    # {
    #   gene_symbol_ 0: {
    #        "ens_id": Ensembl stable ID 
    #        "chr": chromosome
    #        "start": start position (spos)
    #        "end": end position (epos)
    #        "strand": coding or noncoding strand
    #        "mid_point": average of spos and epos
    #   }
    #   gene_symbol_1: {
    #        ...
    #   }   
    #   ... 
    # }
    genes_dict = create_simplified_genes_dictionary(args["genes"], mpnst_genes)


    # for each sample in this MPNST set
    for sample in mpnst_tsv_dict:

        # final list of drivers that will be appended to sample's text file
        # list of dicts, where each dict corresponds to a driver for this sample
        results_final = []

        # checking all columns
        for gene in mpnst_genes:

            # driver-specific dict, which is eventually appended to results_final list
            results = {}
            
            # taking care of somatic drivers
            if gene in mpnst_genes_somatic:

                # flag to add a driver to results_final
                # a driver is only added to results_final if it is one of the specified categories in MAPPING
                # or if it is a homozygous deletion mutation_type
                add = False

                # checking if biallelic based on minor copy number
                if int(sample[f"{gene}_somatic_minor_cn"]) == 0:
                    results["biallelic"] = "yes"
                else:
                    results["biallelic"] = "no"

                # checking if homozygous deletion or another type of somatic mutation
                if int(sample[f"{gene}_somatic_total_cn"]) == 0:
                    results["mutation_type"] = HOMOZYGOUS_DELETION
                    results["category"] = f"(somatic) INDEL"
                    add = True
                elif str(sample[f"{gene}_somatic"]) != "nan":
                    results["mutation_type"] = sample[f"{gene}_somatic"]

                    # setting the "category" column, if applicable
                    # a driver is added to results_final if it is one of the specified categories in MAPPING
                    for key in MAPPING.keys():
                        if str(sample[f"{gene}_somatic"]) in key:
                            results["category"] = f"(somatic) {MAPPING[key]}"
                            add = True 

                # getting the required column values for mutations that are classified as "drivers"
                if add == True:
                    results = fill_driver_file_attrs(gene, genes_dict, results)
                    results_final.append(results)

            # reset results dict -- it is possible for a gene 
            # to acquire both somatic and germline driver mutations (e.g. NF1)
            results = {}

            # now handling germline mutations
            if gene in mpnst_genes_germline:

                mutation_type = sample[f"{gene}_germline"]

                for key in MAPPING.keys():
                    category = f"(germline) {MAPPING[key]}"
                    
                    # if mutation_type is one of the specified categories in MAPPING
                    if str(mutation_type) in key:
                        
                        # handle NF1 -- special parsing
                        if gene == "NF1":
                            details_parsed = parse_germline_details_NF1(
                                sample[f"{gene}_germline_detailed"],
                                gene,
                                category,
                                mutation_type,
                            )
                            if details_parsed is not None:
                                results_final += details_parsed
                            else:
                                results = fill_driver_file_attrs(gene, genes_dict, results, category, mutation_type)
                                results_final.append(results)
                        # handle all other genes
                        else:
                            results = fill_driver_file_attrs(gene, genes_dict, results, category, mutation_type)
                            results_final.append(results)


        if len(results_final) > 0: # catching samples that are only LOH
        
            # fill out attributes that have not been completed
            for elem in results_final:
                for attr in DRIVER_FILE_ATTRS:
                    elem.setdefault(attr, "")

            # and write the drivers to new tab-delimited txt file for this sample
            results_final = [
                dict(s) for s in set(frozenset(d.items()) for d in results_final)
            ]
            with open(f'{output_folder}/{sample["sample"]}.txt', "w") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=DRIVER_FILE_ATTRS, delimiter="\t")
                writer.writeheader()
                writer.writerows(results_final)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat MPNST drivers")
    parser.add_argument("-m", "--mpnst", help="MPNST drivers TSV", required=True)
    parser.add_argument("-g", "--genes", help="CGAP genes JSON", required=True)
    args = vars(parser.parse_args())
    reformat_drivers(args)