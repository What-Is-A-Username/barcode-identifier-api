#!/bin/bash
# Remove all containers and delete associated data in their volumes

RED_COLOR_ESCAPE='\033[0;31m'
GREEN_COLOR_ESCAPE='\033[0;32m'
NO_COLOR_ESCAPE='\033[0m'
YELLOW_COLOR_ESCAPE='\033[1;33m'

echo -e "${RED_COLOR_ESCAPE}WARNING: This run script will STOP and REMOVE any previously running containers, and DELETE existing data." 
echo -e "To keep data, transfer any files out of the containers and volumes to elsewhere.${NO_COLOR_ESCAPE}"
read -p "Proceed? (y/n): " answer
if [ "$answer" == "${answer#[Yy]}" ]; then 
    echo "Exiting ..."
    exit
fi

echo -e "${YELLOW_COLOR_ESCAPE}(server_delete.sh) Stopping and removing Docker containers ...${NO_COLOR_ESCAPE}"
echo "-------"
docker compose -f docker-compose.yml -p barrel down --volumes
echo "-------"
echo -e "${GREEN_COLOR_ESCAPE}(server_delete.sh) End of script.${NO_COLOR_ESCAPE}"
