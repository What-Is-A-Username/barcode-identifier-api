#!/bin/bash

GREEN_COLOR_ESCAPE='\033[0;32m'
NO_COLOR_ESCAPE='\033[0m'
YELLOW_COLOR_ESCAPE='\033[1;33m'

echo -e "${YELLOW_COLOR_ESCAPE}(server_pull.sh) Pulling and building Docker images from the web using docker-compose.yml ...${NO_COLOR_ESCAPE}"
echo "-------"
docker compose -f docker-compose.yml pull barrel celery_worker
docker compose -f docker-compose.yml -p barrel build barrel_proxy
echo "-------"
echo -e "${GREEN_COLOR_ESCAPE}(server_pull.sh) End of script.${NO_COLOR_ESCAPE}"