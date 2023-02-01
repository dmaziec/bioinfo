import pandas
import json
import pandas
import csv
import os
import re
from datetime import datetime
import argparse

mapping = {
    "SNV": "mutation",
    "frameshift": "INDEL",
    "frameshift insertion": "INDEL",
    "SV": "SV",
    "nonsynonymous SNV": "SNV",
    "stopgain": "SNV or INDEL",
    "splicing": "SNV or INDEL",
    "deep intronic": "SNV or INDEL",  # NOT SURE
}

FRAMESHIFT_INSERTION = "frameshift insertion"
STOPGAIN = "stopgain"
SPLICING = "splicing"
FRAMESHIFT_DELETION = "frameshift deletion"
NONSYNONYMOUS_SNV = "nonsynonymous SNV"

muation_types = [
    FRAMESHIFT_INSERTION,
    STOPGAIN,
    SPLICING,
    FRAMESHIFT_DELETION,
    NONSYNONYMOUS_SNV,
]


def parse_germline_details(record, category, mutation_type):
    def get_mutation_type(token, occurence):
        for mutation in muation_types:
            if re.search(mutation, token):
                return re.findall(mutation, token)[occurence]

    def get_chrom_start_end(token):
        pattern = re.search("chr\d* \d* \d* [ACGT-]* [ACGT-]*", token)
        if pattern:
            return pattern.group()

    def get_gene(token, gene):
        return re.search(f"{gene}", token)

    all_results = []
    records = record.split(";;")

    for r in records:

        chrom_start_end = get_chrom_start_end(r)

        if chrom_start_end == None:
            return None
        results = {}
        chrom, start, end, ref, alt = chrom_start_end.split(" ")
        results["chr"] = chrom
        results["pos"] = start
        results["ref"] = ref
        results["alt"] = alt
        results["category"] = category
        results["mutation_type"] = mutation_type
        mutation_type = get_mutation_type(r, 0)
        results["mutation_type"] = f"germline {mutation_type}"
        gene = get_gene(r, "NF1")
        results["gene"] = gene.group()
        gene_end = gene.end()
        proteins = r[gene_end + 1 :].split(" ")[0]
        results["protein_change"] = proteins
        all_results += [results]

    return all_results


def main(args):
    dateTimeObj = datetime.now()
    timestampStr = dateTimeObj.strftime("%d_%b_%Y_%H_%M_%S_%f")

    output_folder = f"MPNST_drivers/{timestampStr}"
    if os.path.exists(output_folder) != True:
        os.makedirs(output_folder)

    mpnst_tsv_filtered = pandas.read_csv(
        args["mpnst"],
        sep="\t",
        index_col=False).filter(
        regex="(^GeM_ID$)|(^sample$)|(^(germline)$)|(.*_somatic$)|(.*_germline$)|(.*_germline_detailed$)|(.*_somatic_minor_cn$)|(.*_somatic_total_cn$)",
        axis=1)

    mpnst_tsv_filtered.drop(
        inplace=True,
        columns=[
            "CXorf67_somatic",
            "CXorf67_somatic_minor_cn",
            "CXorf67_somatic_total_cn",
            "COR_somatic",
            "COR_somatic_minor_cn",
            "COR_somatic_total_cn",
            "AZF1_somatic",
            "AZF1_somatic_minor_cn",
            "AZF1_somatic_total_cn",
            "ROS_somatic",
            "ROS_somatic_minor_cn",
            "ROS_somatic_total_cn",
            "FAM22_somatic",
            "FAM22_somatic_minor_cn",
            "FAM22_somatic_total_cn",
            "ZC3H7_somatic",
            "ZC3H7_somatic_minor_cn",
            "ZC3H7_somatic_total_cn",
            "INI1_somatic",
            "INI1_somatic_minor_cn",
            "INI1_somatic_total_cn"

        ],
    )  # dropping these genes as I cant find them in cgap genes and they are empty, checked manually

      # just set it to some number insetad of NaN
    mpnst_genes_somatic = set(
        [
            col.split("_")[0]
            for col in list(mpnst_tsv_filtered.columns)
            if (col.endswith("_somatic"))
        ]
    )
    mpnst_genes_germline = set(
        [
            col.split("_")[0]
            for col in list(mpnst_tsv_filtered.columns)
            if (col.endswith("_germline"))
        ]
    )
    mpnst_genes = set(list(mpnst_genes_somatic) + list(mpnst_genes_germline))
    col_mino_total = [col for col in mpnst_tsv_filtered.columns if (col.endswith("_somatic_minor_cn") or col.endswith("_somatic_total_cn"))]
    mpnst_tsv_filtered[col_mino_total] = mpnst_tsv_filtered[col_mino_total].fillna(value=-1000)       # just set it to some number insetad of NaN


    mpnst_tsv_dict = mpnst_tsv_filtered.to_dict(orient="records")

    genes_file = open(args["genes"])
    genes = json.load(genes_file)

    genes_ = {
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
        if gene["gene_symbol"] in mpnst_genes
    }

    for sample in mpnst_tsv_dict:
        results_final = []
        for gene in mpnst_genes:
            results = {}
            # flag if add it to the results, only if we have a specified category or if category is missing we have homozygous deletion
            
            if gene in mpnst_genes_somatic:
                add = False
                if int(sample[f"{gene}_somatic_minor_cn"]) == 0:
                    results["biallelic"] = "YES"
                else:
                    results["biallelic"] = "NO"

                if int(sample[f"{gene}_somatic_total_cn"]) == 0:
                    results["mutation_type"] = "homozygous deletion"
                    add = True
                elif str(sample[f"{gene}_somatic"]) != "nan":
                        results["mutation_type"] = sample[f"{gene}_somatic"]
                        add = True

                for key in mapping.keys():
                    if str(sample[f"{gene}_somatic"]) in key:

                        results[f"category"] = mapping[key]
                        add = True

                if add == True:
                    results["gene"] = gene
                    results["chr"] = f'chr{genes_[gene]["chr"]}'
                    results["pos"] = genes_[gene]["mid_point"]
                    results_final.append(results)
            results = {}
            if gene in mpnst_genes_germline:

                mutation_type = sample[f"{gene}_germline"]

                for key in mapping.keys():
                    category = f"(germline) {mapping[key]}"
                    if str(sample[f"{gene}_germline"]) in key:

                        # check NF1 that has all the details
                        if gene == "NF1":
                            details_parsed = parse_germline_details(
                                sample[f"{gene}_germline_detailed"],
                                category,
                                mutation_type,
                            )
                            if details_parsed != None:
                                results_final += details_parsed
                            else: 
                                results["category"] = category
                                results["mutation_type"] = mutation_type
                                results["gene"] = gene
                                results["chr"] = f'chr{genes_[gene]["chr"]}'
                                results["pos"] = genes_[gene]["mid_point"]
                                results_final.append(results)
                        else:
                            results["category"] = category
                            results["mutation_type"] = mutation_type
                            results["gene"] = gene
                            results["chr"] = f'chr{genes_[gene]["chr"]}'
                            results["pos"] = genes_[gene]["mid_point"]
                            results_final.append(results)

        all_keys = []
        for elem in results_final:
            all_keys += list(elem.keys())

        results_final = [
            dict(s) for s in set(frozenset(d.items()) for d in results_final)
        ]

        field_names = sorted(list(set(all_keys)))
        with open(f'{output_folder}/{sample["sample"]}.tsv', "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=field_names, delimiter="\t")
            writer.writeheader()
            writer.writerows(results_final)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat MPNST drivers")
    parser.add_argument("-m", "--mpnst", help="MPNST drivers TSV", required=True)
    parser.add_argument("-g", "--genes", help="CGAP genes JSON", required=True)
    args = vars(parser.parse_args())
    main(args)
