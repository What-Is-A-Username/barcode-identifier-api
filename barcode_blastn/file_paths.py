import os

def get_ncbi_folder(ncbi_blast_version: str = "ncbi-blast-2.12.0+"):
    return os.path.abspath(f'./{ncbi_blast_version}/bin/')

def get_project_path():
    return os.path.abspath(f".")

def get_data_run_path(run_id: str = ''):
    '''Get absolute path of where to store all run files'''
    return os.path.abspath(f'/var/data/runs/{run_id}/')

def get_data_fishdb_path(db_id: str = ''):
    '''Get absolute path of where to store fishdb'''
    return os.path.abspath(f'/var/data/fishdb/{db_id}/')

def get_static_run_path(run_id: str = ''):
    '''Get absolute path of where to store run files to be served'''
    return os.path.abspath(f'/vol/static/runs/{run_id}/')

def get_static_collection_path():
    '''Get absolute path of where to collect static files'''
    return os.path.abspath(f'/vol/static/static/')

def get_static_media_path():
    '''Get absolute path of where to store media files'''
    return os.path.abspath(f'/vol/static/media/')