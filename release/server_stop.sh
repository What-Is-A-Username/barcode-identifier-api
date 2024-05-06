#!/bin/bash

RED_COLOR_ESCAPE='\033[0;31m'
GREEN_COLOR_ESCAPE='\033[0;32m'
NO_COLOR_ESCAPE='\033[0m'
YELLOW_COLOR_ESCAPE='\033[1;33m'

echo -e "${RED_COLOR_ESCAPE}WARNING: This run script will STOP any previously running containers but preserve any existing data" 
read -p "Proceed? (y/n): " answer
if [ "$answer" == "${answer#[Yy]}" ]; then 
    echo "Exiting ..."
    exit
fi
echo -e "${YELLOW_COLOR_ESCAPE}(server_stop.sh) Stopping Docker containers ...${NO_COLOR_ESCAPE}"
echo "-------"
docker compose -f docker-compose.yml -p barrel stop 
echo "-------"
echo -e "${GREEN_COLOR_ESCAPE}(server_stop.sh) End of script.${NO_COLOR_ESCAPE}"
