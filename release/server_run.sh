#!/bin/bash

RED_COLOR_ESCAPE='\033[0;31m'
GREEN_COLOR_ESCAPE='\033[0;32m'
NO_COLOR_ESCAPE='\033[0m'
YELLOW_COLOR_ESCAPE='\033[1;33m'

echo -e "${RED_COLOR_ESCAPE}WARNING: This run script will attempt to start the required containers not yet running. If there"
echo -e "are containers already running, they will be recreated. To get rid of previous data, run ./server_delete.sh"
echo -e "before proceeding with the current script.${NO_COLOR_ESCAPE}" 

read -p "Proceed? (y/n): " answer
if [ "$answer" == "${answer#[Yy]}" ]; then 
    echo "Exiting ..."
    exit
fi

echo -e "${YELLOW_COLOR_ESCAPE}(server_run.sh) Running Docker containers ...${NO_COLOR_ESCAPE}"
echo "-------"
docker compose -f docker-compose.yml -p barrel up barrel barrel_proxy celery_worker -d 
echo "-------"
echo -e "${GREEN_COLOR_ESCAPE}(server_run.sh) End of script.${NO_COLOR_ESCAPE}"