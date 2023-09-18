from typing import List
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix
from math import log

def parse_id_from_msa_id(id: str):
    return id.split('|')

def calculate_genetic_distance_for_sequence(x_seq: str, y_seq: str) -> float:
    '''
    Return the Kimura-2-parameter distance between sequences x_seq and y_seq.
    '''
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

def calculate_genetic_distance(alignment_file_path: str) -> DistanceMatrix:
    '''
    Return a matrix of all pairwise Kimura-2-parameter distances for the sequences contained
    in the alignment file specified by `alignment_file_path`. 
    '''
    # Read sequences and identifiers from MSA file
    alignment: AlignIO.MultipleSeqAlignment = AlignIO.read(open(alignment_file_path), "clustal")
    records: List[SeqRecord] = [record for record in alignment]       
    matrix: DistanceMatrix = DistanceMatrix(names=[record.id for record in records])
    # Loop through all pairs of sequences
    for x_index in range(0, len(records) - 1, 1):
        for y_index in range(x_index + 1, len(records)):
            x: SeqRecord = records[x_index]
            y: SeqRecord = records[y_index]
            matrix[x.id, y.id] = calculate_genetic_distance_for_sequence(x.seq, y.seq)
            matrix[y.id, x.id] = matrix[x.id, y.id]

    return matrix