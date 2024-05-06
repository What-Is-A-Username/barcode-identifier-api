# For Developers

This page is for developers who want to run the repository in the development environment, in order to add or test changes.

## Setup the documentation generation

Barrel uses [MKDocs](https://www.mkdocs.org/) to generate the documentation for Barrel.

MkDocs runs with Python and the required packages can be downloaded into a new virtual environment using [./docs_requirements.txt](./docs_requirements.txt).

Starting from the main project directory, enter the docs folder
```
cd docs
```

Then create the virtual environment `docs_env`:
```
python3 -m venv docs_env -r 
source docs_env/bin/activate
pip install -r docs_requirements.txt
```

To preview how the docs will look like:
```
mkdocs serve
```

## Build for production

The documentation is made available to site visitors is by building the docs to static site files, and then making the files available to be served by the site.

To do this, build the files:
```
mkdocs build
```

Move the files to the proxy.
```
cp -r ./site ../compose/proxy/docs
```

If the proxy container is already running, ensure that the files are in the container. This is done by rebuilding and restarting the proxy container.

## Building a release

A version release is composed of two major steps
- building the docker images and uploading them to the web (Docker Hub)
- creating the release file, which is a compressed file which users download for installation

### Building and pushing images

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

### Create the release files
The release file includes:
- docker compose file required by users to download and run the docker images
- the configuration files for the server
- files required to show the website documentation
- (optional) files required for the frontend web application

Before building the release, ensure that the virtual environment for building the documentation is already setup. Consult the [docs README](../README.md) for more information.

From the main project directory, simply run the shell script, which will create the `.tar.gz` and `.zip` files.
```
./dev/release_build.sh
```