import os
from barcode_blastn.models import BlastRun

def readHitTreeFromFile(run: BlastRun) -> str:
    '''Read the Newick string of the hit tree (tree with query sequences and their hits)
    stored in the .phylotree.ph file.
    '''
    return readTreeFromFile(str(run.id), run.alignment_job_id)

def readDbTreeFromFile(run: BlastRun) -> str:
    '''Read the Newick string of the database tree (tree with query sequences and all database sequences hits)
    stored in the .phylotree.ph file.
    '''
    return readTreeFromFile(str(run.id), run.complete_alignment_job_id)

def readTreeFromFile(run_id: str, job_id: str) -> str:
    '''Return the tree as a string by reading it from file.
    '''
    run_folder = os.path.abspath(f'./runs/{run_id}')
    run_files = os.listdir(run_folder)
    try:
        phylip_file = [r for r in run_files if r == f"{job_id}.phylotree.ph"][0]
    except IndexError:
        print("Warning: Tried to find local copy of tree PHYLIP file, but no file was found.")
        # TODO: return error if no file is actually found?
        return ''
    else:
        with open(f'{run_folder}/{phylip_file}') as tree_file:
            tree_text = tree_file.read()
        return tree_text