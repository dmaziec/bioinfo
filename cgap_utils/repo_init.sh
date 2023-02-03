mkdir -p descriptions
mkdir -p dockerfiles
mkdir -p portal_objects
mkdir -p portal_objects/metaworkflows
mkdir -p portal_objects/workflows
touch portal_objects/file_format.yaml
touch portal_objects/file_reference.yaml
touch portal_objects/software.yaml
basename "$PWD" > PIPELINE
echo "1.0.0" > VERSION