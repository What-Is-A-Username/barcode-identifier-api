# Barrel Documentation

Barrel uses [MKDocs](https://www.mkdocs.org/) to generate the documentation for Barrel.

MkDocs runs with Python and the required packages can be downloaded into a new virtual environment using [./docs_requirements.txt](./docs_requirements.txt).

```
python3 -m venv docs_env -r 
source docs_env/bin/activate
pip install -r docs_requirements.txt
mkdocs serve
```