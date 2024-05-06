# Barrel (Release Version)

This is a release copy of Barrel for distribution.

Barrel is an upcoming collaborative platform, provided as a web API, for simplified barcoding analysis.

## Installation

The instructions below will tell you how to install Barrel on your machine.

Barrel runs within Docker containers, and thus requires an [installation of the Docker Engine](https://docs.docker.com/engine/install/). Install it before moving onto [first-time setup](#first-time-setup) below.

### First-time Setup

First-time setup must be done first before Barrel can be run. Ensure the working directory of your terminal / command prompt is the same folder as this README file.

Run `server_setup.sh` to setup an environment file in this directory containing required passwords and secrets for setting up Barrel. You will be provided two options: a public or private installation. Choose public if you'd like your installation to be used by outside users over the Internet, and private if you'd like to keep it isolated to the same machine.
```
./server_setup.sh
```

Run `server_pull.sh` to start downloading the Docker images from the web, which will be used to create the containers in which Barrel runs. 
```
./server_pull.sh
```

### Run

Given that first-time setup has already been performed, running Barrel at any subsequent time can be done by just running `server_run.sh`. 
```
./server_run.sh
```

Setup time may take a while and is dependent on the available computing resources. It may take upwards of a minute for all aspects to become functional. Once functional, you will be able to access the site at [localhost port 8000](http://localhost:8000). 


## Usage
DNA barcoding is used widely in biology, especially in biodiversity conservation, to identify specimens originating from unknown species. It is a technique utilizing short DNA sequences or markers from standardized regions of the genetic code to identify species. Relying on the natural genetic variation between different species, sequences of unknown origin can be compared against sequences from known species to make an identification.

The present app, known as "Barrel", provides a server platform for researchers to curate collections of genetic sequences that can be used for species identification using DNA barcoding. Server adminstrators can use the console at `/admin` to manage the service and upload genetic sequences which are organized into separate BLAST databases. End-users can query against the uploaded data by submitting a POST request to the `/blastdbs/<id>/run`, which submits a job to be processed asynchronously. 

The ultimate goal of this project is to create a server solution that can be quickly setup by researchers for their specific purposes, to give them each greater control over the genetic barcoding data they make available to the public.

## Source Code
Source code for Barrel is provided across different projects.
- This repository (/barrel) contains the latest official releases for download, as well as the documentation for all aspects of the project.
- The backend API which performs the analysis is found at [/barcode-identifier-api](https://github.com/clwillhuang/barcode-identifier-api)
- The frontend web client which forms the user interface, excluding the administrator console, is over at [/barcode-identifier-app](https://github.com/clwillhuang/barcode-identifier-app)

The source code provided in this repository, and in any other repository related to it (e.g. ) may not be used for publication without the expressed permission from original Barrel developer(s) and contributors. 