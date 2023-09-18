#!/bin/sh 
# For development purposes. DO NOT RUN in deployment container.

# Use a busybox docker container to interactively view files in volume.
docker run -it --rm -v barcode-identifier-api_parent-data:/vol busybox sh