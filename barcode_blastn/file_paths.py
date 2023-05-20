import os

from barcode_blastn.models import BlastDb, Library

def get_ncbi_folder(ncbi_blast_version: str = "ncbi-blast-2.12.0+"):
    '''
    Return the absolute path to the /bin/ folder of the local ncbi
    installation in the project directory.
    '''
    return os.path.abspath(f'./{ncbi_blast_version}/bin/')

def get_project_path():
    '''
    Return the absolute path of the project directory.
    '''
    return os.path.abspath(f".")

def get_data_run_path(run_id: str = ''):
    '''
    Get absolute path of the folder containing the main run files.
    '''
    return os.path.abspath(f'/var/data/runs/{run_id}/')

def get_data_library_path(lib: Library):
    return os.path.abspath(f'/var/data/library/{lib.id}/')

def get_data_fishdb_path(db: BlastDb):
    '''
    Get absolute path of where to store BLAST database files made
    when running makeblastdb
    '''
    return os.path.abspath(f'{get_data_library_path(db.library)}/{str(db.id)}/')

def get_static_run_path(run_id: str = ''):
    '''
    Get absolute path of where to store run files to be served
    through the server.
    '''
    return os.path.abspath(f'/vol/static/runs/{run_id}/')

def get_static_collection_path():
    '''
    Get absolute path of where to collect static files
    '''
    return os.path.abspath(f'/vol/static/static/')

def get_static_media_path():
    '''
    Get absolute path of where to store media files to be
    served through the server.
    '''
    return os.path.abspath(f'/vol/static/media/')