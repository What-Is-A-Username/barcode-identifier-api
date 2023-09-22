# Barrel API

The Barrel API is the backend server API powering the [Barrel app](https://barrel.utsc.utoronto.ca/app).

It is a web service, provided as a Django web API, which allows biological researchers to collaboratively work on DNA barcoding projects through a shared server platform.

The project is still under active development and testing. In the future, downloadable releases will be available at our [main repository](https://github.com/clwillhuang/barrel).

## Quickstart

Below are the steps for running the backend API for Barrel locally, meant for those familiar with development servers such as software developers and bioinformaticians. 
-   If you are instead looking to host your own version of Barrel, [visit and download from official releases](https://github.com/clwillhuang/barrel/releases).
-   If you are instead looking for documentation on how to use Barrel for your own sequence analysis and research, view our documentation 

The present project was developed for Linux Ubuntu. 

Ensure that Docker is installed on the machine. If not, install it over at https://docs.docker.com/engine/install/ubuntu/.

Run the containers using Docker Compose:
```
docker compose -f ./docker-compose-dev.yml -p barrel build barrel_venv_image
docker compose -f ./docker-compose-dev.yml -p barrel up --build celery_worker barrel barrel_proxy --no-cache
```

The primary difference between the development (docker-compose-dev.yml) and the deployment builds (docker-compose-deploy.yml) is that in the development build, the files are in bind mounts so that they are reachable more easily through the local filesystem. For deployment, files are in named volumes managed by Docker.

### Docs (optional)
This repo is also accompanied with documentation built with MKDocs.

First, install the required dependencies:
```
cd docs
python3 -m venv docs_env -r 
source docs_env/bin/activate
pip install -r docs_requirements.txt
```

Then start the MKDocs dev server:
```
mkdocs serve
```

## Usage

DNA barcoding is used widely in biology, especially in biodiversity conservation, to identify specimens originating from unknown species. It is a technique utilizing short DNA sequences or markers from standardized regions of the genetic code to identify species. In animals, researchers typically rely on the cytochrome oxidase 1 gene for this purpose, although there exists many other markers for animals and other taxonomic groups.

The present app, known as "Barrel", provides a server platform for researchers to curate collections of genetic sequences that can be used for species identification using DNA barcoding. Server adminstrators can use the console at `/admin` to manage the service and upload genetic sequences which are organized into separate BLAST databases. End-users can query against the uploaded data by submitting a POST request to the `/blastdbs/<id>/run`, which submits a job to be processed asynchronously. 

The ultimate goal of this project is to create a server solution that can be quickly setup by researchers for their specific purposes, to give them each greater control over the genetic barcoding data they make available to the public.




