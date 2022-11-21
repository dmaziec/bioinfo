import pandas
import glob
import os
import shutil
import argparse
from datetime import datetime


"""
This script is to reformat output CNV data from ASCAT into the format compatible with GosCan 
Under the provided directory (inputdir) files should be named with the suffix: .copynumber.caveman.csv
Input directory will be searched recursively for the files and the files should be stored up to one level from the parent directory e.g. parentdir/ascat/out.copynumber.caveman.csv
the results will stored under the CNV_goscan folder
"""


def main(args):
    dateTimeObj = datetime.now()
    timestampStr = dateTimeObj.strftime("%d_%b_%Y_%H_%M_%S_%f")

    output_folder = f"CNV_goscan/{timestampStr}"
    if os.path.exists(output_folder) != True:
        os.mkdir(output_folder)
    path = args["inputdir"]

    for file_csv in glob.glob(f"{path}/**/*.copynumber.caveman.csv", recursive=True):

        file_name = os.path.basename(file_csv)

        folder_name = file_csv.split("/")[-2]

        ascat_cnv = pandas.read_csv(
            file_csv,
            sep=",",
            names=[
                "ID",
                "chromosome",
                "start",
                "end",
                "total_normal",
                "minor_normal",
                "total_tumor",
                "minor_tumor",
            ],
        )

        ascat_cnv["total_tumor"] = ascat_cnv["total_tumor"].astype(int)
        ascat_cnv["minor_tumor"] = ascat_cnv["minor_tumor"].astype(int)

        ascat_cnv["major_cn"] = ascat_cnv["total_tumor"] - ascat_cnv["minor_tumor"]
        goscan_cnv = ascat_cnv[
            ["chromosome", "start", "end", "total_tumor", "minor_tumor", "major_cn"]
        ]
        # change the order of columns
        goscan_cnv = goscan_cnv[
            ["chromosome", "start", "end", "major_cn", "minor_tumor", "total_tumor"]
        ]
        goscan_cnv.rename(
            columns={"minor_tumor": "minor_cn", "total_tumor": "total_cn"}, inplace=True
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
    parser.add_argument("-i", "--inputdir", help="input CSV directory", required=True)

    args = vars(parser.parse_args())

    main(args)
