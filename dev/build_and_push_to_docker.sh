#!/bin/bash
# Build and push dockers to Docker from dev environment.

RED_COLOR_ESCAPE='\033[0;31m'
NO_COLOR_ESCAPE='\033[0m'

echo -e "${RED_COLOR_ESCAPE}WARNING: This setup script will rebuild images and push to docker.${NO_COLOR_ESCAPE}" 
read -p "Enter version of new build: " version
read -p "Enter Docker Hub username:" username
read -p "Enter Docker Hub repository name: " repo

docker compose -f ./docker-compose-deploy.yml -p barrel build barrel_venv_image
docker compose -f ./docker-compose-deploy.yml -p barrel build celery_worker barrel --no-cache

docker tag barrel-barrel ${username}/${repo}:barrel-${version}
docker push ${username}/${repo}:barrel-${version}

docker tag barrel-celery_worker ${username}/${repo}:celery_worker-${version}
docker push ${username}/${repo}:celery_worker-${version}

echo -e "Success"