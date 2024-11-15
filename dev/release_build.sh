#!/bin/bash
# Gather files to create a .zip file release

# Prompt for desired version name of the build
read -p "Enter version of new build: " version

# Build docs
cd docs
source ./docs_env/bin/activate
mkdocs build
deactivate
cp -r site ../release/docs
cd ..

# Premake folders
mkdir -p ./release/volumes/postgres-data
mkdir -p ./release/volumes/static-data
mkdir -p ./release/volumes/var-data

tar -czvf ${version}.tar.gz release
zip -r ${version}.zip release

# Cleanup
rm -r ./release/docs
rmdir ./release/volumes/postgres-data
rmdir ./release/volumes/static-data
rmdir ./release/volumes/var-data
rmdir ./release/volumes

echo -e "Success"