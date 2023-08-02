# Docker Containers

Barrel leverages Docker to simplify aspects of the setup process, segregate different components of the workflow, and deliver potential performance improvements.

## Description

### barrel_venv_image

Barrel uses an intermediate image which contains the Python installation used to run both the web application and the worker processes that handle jobs.

This image also includes all the code required to run the app. 

### barrel_rabbitmq

Barrel uses RabbitMQ, a message broker used to queue up user-submitted jobs.

### celery_worker

This image is used to handle the Celery worker process, which processes jobs off the queue provided by barrel_rabbitmq.

### barrel

This is the main image of the app, and houses the actual web application which handle requests and generates responses. 

### barrel_proxy

This is an NGINX image used as a reverse proxy. Essentially, this is the part of the app that first handles incoming web requests. Based on the request URL, it either responds directly with files from the server's filesystem or forwards the requests to the barrel container.

## For server adminstrators

This section describes how to setup a Barrel server. 

Setup involves retrieving (i.e. pulling) Docker images from the web, which will be used to create containers on your server which will run Barrel and its constituent components.

Before starting the below instructions, ensure that you:
- have a Docker account
- have installed Docker on your machine

Ensure that you're logged into Docker. Substitute `<username>` with your docker username
```
docker login -u <username>
```

Download the latest release of Barrel `docker-compose.yml` file from GitHub.

Unzip the file.

Run `setup.sh`, which will prompt you to setup various passwords and secrets that will secure your machine and computer. This automatically creates a `.env` file in the same directory when finished. This information includes server adminstrator credentials, so do not share them publcly! Sharing these passwords and secrets will give parties virtually unlimited access.
```
./setup.sh
```

Retrieve the images from the web by pulling from Docker:
```
```

## For developers

### Building containers in the development environment.

This section is a reference for developers about how to build the Docker containers from the complete source code (i.e. not pulling images from the Docker Hub website). The following assumes that Docker is already installed on the system and you already possess a Docker account.

From the terminal, ensure that you're logged into Docker. Login by running the command, substituting in your username:
```
docker login -u <username>
```

Build and run the images. Below is an example of the command used to build the production image for barrel_venv_image using the docker compose files.
```
docker compose -f ./docker-compose-deploy.yml -p barrel up --build barrel_venv_image
```

Tag the image.
```
docker tag barrel_venv_image username/barrel-dev:barrel_venv_image
```

Push to Docker.
```
docker push username/barrel-dev:barrel_venv_image
``` 

To build the docker containers, build them in the following order:
- barrel_venv_image
- celery_worker
- barrel
- barrel_proxy


