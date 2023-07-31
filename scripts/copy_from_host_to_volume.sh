#!/bin/sh 
# For development purposes. DO NOT RUN in deployment container.

# Use docker to copy files from host to volume
echo "Making helper container to facilitate file transfer"
docker run -v barcode-identifier-api_parent-data:/data --name helper busybox true
echo "Performing file transfer."
docker cp ./scripts/. helper:/data/scripts
docker cp ./barcode_blastn/. helper:/data/barcode_blastn
docker cp ./barcode_identifier_api/. helper:/data/barcode_identifier_api
docker cp ./ncbi-blast-* helper:/data
echo "File transfer finished. Destroying container."

docker rm helper
echo "Exiting."