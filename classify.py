from collections import namedtuple
from typing import List, Optional
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceCalculator
from math import log, sqrt
import os

def get_static_run_path(run_id: str = ''):
    '''
    Get absolute path of where to store run files to be served
    through the server.
    '''
    return os.path.abspath(f'/vol/static/runs/{run_id}/')

def calculate_genetic_distance_for_sequence(x_seq: str, y_seq: str) -> float:
    num_transitions = 0
    num_transversions = 0
    num_bases = 0
    for index in range(0, len(x_seq)):
        x = x_seq[index]
        y = y_seq[index]
        if x == '-' or y == '-':
            continue
        change = x+y if x < y else y+x
        if change in ['AG', 'CT']:
            num_transitions = num_transitions + 1
        elif change in ['AC', 'AT', 'CG', 'GT']:
            num_transversions = num_transversions + 1
        num_bases = num_bases + 1

    freq_transitions = num_transitions / num_bases 
    freq_transversions = num_transversions / num_bases

    result = -0.5 * log(1 - 2*freq_transitions - freq_transversions) - 0.25 * log(1 - 2*freq_transversions)
    return result

def classify_genetic_distance(alignment_successful: bool, alignment_job_id, run_id: str):
    alignment_file_path = f'{get_static_run_path(run_id)}/{alignment_job_id}.aln-clustal_num.clustal_num'
    alignment: AlignIO.MultipleSeqAlignment = AlignIO.read(open(alignment_file_path), "clustal")
    print("Alignment length %i" % alignment.get_alignment_length())
    records: List[SeqRecord] = [record for record in alignment]       
    matrix: DistanceMatrix = DistanceMatrix(names=[record.id for record in records])
    for x_index in range(0, len(records) - 1, 1):
        for y_index in range(x_index + 1, len(records)):
            x: SeqRecord = records[x_index]
            y: SeqRecord = records[y_index]
            matrix[x.id, y.id] = calculate_genetic_distance_for_sequence(x.seq, y.seq)
            matrix[y.id, x.id] = matrix[x.id, y.id]
    return matrix

alignment_successful = True 
alignment_job_id = 'clustalo-R20230605-153354-0685-34173264-p1m'
run_id = '3046fd36-6a35-4334-ae94-9cffb8a996e8'
classify_genetic_distance(alignment_successful, alignment_job_id, run_id)