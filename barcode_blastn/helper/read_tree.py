import os
from typing import Union
from barcode_blastn.file_paths import get_data_run_path, get_static_run_path
from barcode_blastn.models import BlastRun

def readHitTreeFromFile(run: BlastRun) -> Union[str, None]:
    '''Read the Newick string of the hit tree (tree with query sequences and their hits)
    stored in the .phylotree.ph file.
    '''
    print('Reading hit tree.')
    return readTreeFromFile(run, run.alignment_job_id)

def readDbTreeFromFile(run: BlastRun) -> Union[str, None]:
    '''Read the Newick string of the database tree (tree with query sequences and all database sequences hits)
    stored in the .phylotree.ph file.
    '''
    print('Reading db tree.')
    return readTreeFromFile(run, run.complete_alignment_job_id)

def readTreeFromFile(run: BlastRun, job_id: str) -> Union[str, None]:
    '''Return the tree as a string by reading it from file. Returns None if the tree file could not
    be located
    '''
    run_folder = get_static_run_path(str(run.id))
    run_files = os.listdir(f'{run_folder}')
    try:
        # Identify the tree file based on its unique name and file extension
        phylip_file = [r for r in run_files if r == f"{job_id}.phylotree.ph"][0]
    except IndexError:
        # Set error if file could not be found.
        print("Warning: Tried to find local copy of tree PHYLIP file, but no file was found.")
        run.errors = run.errors + '\nCould not phylotree file to read phylogenetic tree from.'
        run.status = BlastRun.JobStatus.ERRORED
        run.save()
        return ''
    else:
        # Open the file and return its contents
        with open(f'{run_folder}/{phylip_file}') as tree_file:
            tree_text = tree_file.read()
            return tree_text