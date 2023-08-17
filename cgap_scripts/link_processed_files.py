from dcicutils import ff_utils
import json
import argparse
from pathlib import Path


def get_output_arguments(workflow_runs):

    arguments = {}
    for workflow_run in workflow_runs:
        if 'output' in workflow_run.keys():
            for outp in workflow_run['output']:
                arguments[outp['argument_name']] = outp['file']['uuid']

    return arguments


def link_files(args):

    with open(f'{Path.home()}/.cgap-keys.json') as json_file:
        keys = json.load(json_file)
    ff_key = keys.get(args['auth_key'])


    search = ff_utils.search_metadata(
        f"search/?type=SomaticAnalysis&project.display_title={args['project']}",
        key=ff_key,
    )

    for s in search:
        analysis_uuid = s['uuid']
        processed_files = s['processed_files']
        processed_files_uuid = []
        for file in processed_files:
            processed_files_uuid.append(file['uuid'])

        for metaworkflow_run in s['meta_workflow_runs']:

            if args['metaworkflow_name'] in metaworkflow_run['display_title']:

                search_wfl = ff_utils.search_metadata(
                    f"search/?uuid={metaworkflow_run['uuid']}&type=Item",
                    key=ff_key, )

                output_args_uuids = get_output_arguments(search_wfl[0]['workflow_runs'])
                for argument in args['arguments']:
                    processed_files_uuid.append(output_args_uuids[argument])

        processed_files_uuid = list(set(processed_files_uuid))
        ff_utils.patch_metadata({"processed_files": processed_files_uuid}, obj_id = analysis_uuid, key = ff_key)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Link output files to processed files")
    parser.add_argument("--auth_key", "-a", help="The name of the key in ~/.cgap-keys.json to use")
    parser.add_argument("--project", "-p", help="The name of the project")
    parser.add_argument("--arguments", "-o", help="Argument names", nargs="*")
    parser.add_argument("--metaworkflow_name", "-m", help="Argument names")

    args = vars(parser.parse_args())


    link_files(args)
