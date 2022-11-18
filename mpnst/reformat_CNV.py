import pandas
import glob
import os
import shutil
import argparse


def main(args):
    output_folder = "CNV_goscan"
    if os.path.exists(output_folder) != True:
        os.mkdir(output_folder)
    path = args["inputfile"]

    for file_csv in glob.glob(f"{path}/**/*.copynumber.caveman.csv", recursive=True):

        file_name = os.path.basename(file_csv)

        folder_name = file_csv.split("/")[-2]

        ascat_cnv = pandas.read_csv(
            file_csv,
            sep=",",
            names=[
                "ID",
                "chr",
                "startpos",
                "endpos",
                "total_normal",
                "minor_normal",
                "total_tumor",
                "minor_tumor",
            ],
        )

        ascat_cnv["total_tumor"] = ascat_cnv["total_tumor"].astype(int)
        ascat_cnv["minor_tumor"] = ascat_cnv["minor_tumor"].astype(int)

        ascat_cnv["nMajor"] = ascat_cnv["total_tumor"] - ascat_cnv["minor_tumor"]
        goscan_cnv = ascat_cnv[
            ["chr", "startpos", "endpos", "total_tumor", "minor_tumor", "nMajor"]
        ]
        # change the order of columns
        goscan_cnv = goscan_cnv[
            ["chr", "startpos", "endpos", "nMajor", "minor_tumor", "total_tumor"]
        ]
        goscan_cnv.rename(
            columns={"minor_tumor": "nMinor", "total_tumor": "copyNumber"}, inplace=True
        )
        final_folder = f"{output_folder}/{folder_name}"
        if os.path.exists(final_folder) != True:
            os.mkdir(final_folder)
        file_name = file_name.replace("csv", "tsv")
        goscan_cnv.to_csv(
            f"{output_folder}/{folder_name}/{file_name}", sep="\t", index=False
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reformat CNVs from ASCAT into the GoScan CNV format"
    )
    parser.add_argument("-i", "--inputfile", help="input CSV file", required=True)

    args = vars(parser.parse_args())

    main(args)
