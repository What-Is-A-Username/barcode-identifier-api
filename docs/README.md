# Documentation for Barrel

Barrel uses [MKDocs](https://www.mkdocs.org/) to generate the documentation for Barrel.

MkDocs runs with Python and the required packages can be downloaded into a new virtual environment using [./docs_requirements.txt](./docs_requirements.txt).

```
python3 -m venv docs_env -r 
source docs_env/bin/activate
pip install -r docs_requirements.txt
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

