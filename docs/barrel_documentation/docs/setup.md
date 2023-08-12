# Setup 

## Step 1: Check your system compatibility

The following setup process has only been tested on Linux (Ubuntu).

If you have a Windows machine, look into installing [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install). 

## Step 2: Install with Docker

Download the release from either GitHub. 

Barrel runs within Docker containers, and thus requires an [installation of the Docker Engine](https://docs.docker.com/engine/install/). Install it before moving onto [first-time setup](#first-time-setup) below.

## Step 3: Download Barrel

Barrel can be downloaded from GitHub from [releases page](#https://github.com/What-Is-A-Username/barrel-releases/releases). 

Download it manually in the browser, or download using `wget` in the terminal.

## Step 4: Run first-time setup scripts

First-time setup must be done first before Barrel can be run. Setup is performed by running a few shell scripts that will setup the required files and installations.

Run `server_setup.sh` to setup an environment file in this directory containing required passwords and secrets for setting up barrel. 
```
./server_setup.sh
```

Run `server_pull.sh` to start downloading the Docker images from the web, which will be used to create the containers in which Barrel runs. 
```
./server_pull.sh
```

## Step 5: Run

**Once first-time installation has been completed, you do not have to do it again**. However, steps 3-4 can be completed again to your discretion if you want certain options to be updated.

### Starting Barrel

Given that first-time setup has already been performed, running Barrel at any subsequent time can be done by just running `server_run.sh`. 
```
./server_run.sh
```

Setup time may take a while and is dependent on the available computing resources. It may take upwards of a minute for all aspects to become functional. Once functional, you will be able to access the site on your localhost. For example, if you set your EXPOSED_PORT value to 8000, go to [http://localhost:8000/app](http://localhost:8000). 

### Stopping Barrel

Everything involved with Barrel can be stopped using the `server_stop.sh` script.
```
./server_stop.sh
```

### Restarting Barrel

If you stopped the application, it can be started back up with:
```
./server_run.sh
```
