from dcicutils import ff_utils
import pandas as pd
import json
import argparse
from pathlib import Path

def generate_command(args):

    with open(f'{Path.home()}/.cgap-keys.json') as json_file:
        keys = json.load(json_file)

    ff_key = keys.get(args['auth_key'])

    ff_key_key = ff_key['key']
    ff_key_secret = ff_key['secret']


    samples = ff_utils.search_metadata(
        f"search/?type=Sample&project.display_title={args['project']}",
        key=ff_key,
    )

    portal_url = ff_key['server']
    files = {}

    for s in samples:

        individual_display_title = s['individual']['display_title']
        tissue_type = s['tissue_type']

        if individual_display_title not in files.keys():
            files[individual_display_title] = {"individual_display_title": individual_display_title}

        for processed_file in s['processed_files']:
            processed_file_uuid = processed_file['uuid']
            processed_file_data = ff_utils.search_metadata(
                f"search/?uuid={processed_file_uuid}",
                key=ff_key)

            extension = processed_file_data[0]['file_format']['display_title']
            output_file_name = processed_file_data[0]["href"].split("/")[-1]
            files[individual_display_title][
                f"{tissue_type}_{extension}"] = f'curl -L {portal_url}{processed_file_data[0]["href"]} --user "{ff_key_key}:{ff_key_secret}" -o {output_file_name}'

            if 'extra_files' in processed_file_data[0].keys():
                output_file_name = processed_file_data[0]["extra_files"][0]["href"].split("/")[-1]
                files[individual_display_title][f"{tissue_type}_{extension}_extra_file"] = f'curl -L {portal_url}{processed_file_data[0]["extra_files"][0]["href"]} --user "{ff_key_key}:{ff_key_secret}" -o {output_file_name}'

    sa_project = ff_utils.search_metadata(
            f"search/?type=SomaticAnalysis&project.display_title={args['project']}",
            key=ff_key)

    for sa in sa_project:
        individual_display_title = sa['individual']['display_title']
        for processed_file in sa['processed_files']:
            processed_files = ff_utils.get_metadata(processed_file['uuid'], key=ff_key)
            for file in processed_files['workflow_run_outputs'][0]['output_files']:
                if processed_file['uuid'] == file['value']['uuid']:

                    output_file_name = processed_files["href"].split("/")[-1]
                    files[individual_display_title][file[
                        'workflow_argument_name']] = f'curl -L {portal_url}{processed_files["href"]} --user "{ff_key_key}:{ff_key_secret}" -o {output_file_name}'
                    if 'extra_files' in processed_files.keys():
                        output_file_name = processed_files["extra_files"][0]["href"].split("/")[-1]
                        files[individual_display_title][
                            f"{file['workflow_argument_name']}_extra_file"] = f'curl -L {portal_url}{processed_files["extra_files"][0]["href"]} --user "{ff_key_key}:{ff_key_secret}" -o {output_file_name}'


    pd.DataFrame.from_dict(files.values()).to_csv(args["output"], index=False, sep ="\t")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Generate a TSV file containing download commands for a somatic project.")
    parser.add_argument("--auth_key", "-a", help="The name of the key in ~/.cgap-keys.json to use")
    parser.add_argument("--project", "-p", help="The name of the project")
    parser.add_argument("--output", "-o", help="Output TSV file")
    args = vars(parser.parse_args())


    generate_command(args)
