from collections import namedtuple
from typing import List, Optional
from django.db.models import QuerySet
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from barcode_blastn.controllers.blastdb_controller import compareHits
from barcode_blastn.file_paths import get_data_run_path
from barcode_blastn.models import BlastDb, BlastQuerySequence, BlastRun, HeaderInfo, Hit, NuccoreSequence
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceCalculator
from math import log, sqrt

def parse_id_from_msa_id(id: str):
    return id.split('|')

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

def calculate_genetic_distance(alignment_file_path: str):
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

def annotate_pair_comparison(matrix: DistanceMatrix, run: BlastRun, threshold: float = 0.01):
    '''
    Annotate each query sequence in the matrix using the categorization proposed in
    Janzen et al. 2022.
    
    Names in the matrix each correspond to a reference or query sequence. Reference sequences
    take the form of 'version|species_name' while query sequences take the form of 
    version|species_name|query
    '''
    # 
    names: List[str] = matrix.names
    info_database = [BlastQuerySequence.extract_header_info(name) for name in names]

    # keep a list of all reference species in the db
    refs_in_db: set[str] = set([info.species for info in info_database if not info.is_query])

    # keep a list of all reference AND query ids in the database
    ids_in_db = [info.id for info in info_database]

    seqs: QuerySet[BlastQuerySequence] = run.queries.all()
    seq: BlastQuerySequence

    output_path = get_data_run_path(str(run.id))
    with open(output_path + '/k2p_matrix.phy', 'w') as matrix_handle:
        matrix.format_phylip(matrix_handle)

    debug_handle = open(output_path + '/debug.txt', 'w')
    debug_handle.write(str(names))
    debug_handle.write(str(refs_in_db))
    debug_handle.write(str(ids_in_db))

    classification_path = f'{output_path}/classification.tsv'
    class_handle = open(classification_path, 'w')
    class_handle.write('query_id\ttree_id\tquery_species\treference_species\taccuracy_category\n')
    try:
        for seq in seqs:
            best_hit: Optional[Hit] = None
            hit: Hit
            for hit in seq.hits.all():
                if best_hit is None or compareHits(best_hit, hit):
                    best_hit = hit
            
            query_id = seq.write_tree_identifier()
            query_species = seq.original_species_name if not seq.original_species_name is None else ''
                
            if best_hit is None:
                # If no hits are returned 
                seq.accuracy_category = BlastQuerySequence.QueryClassification.NO_HITS
                class_handle.write(f'{query_id}\t{seq.write_tree_identifier()}\t{query_species}\tNo hits\t{seq.accuracy_category}\n')
                continue
            
            reference: NuccoreSequence = best_hit.db_entry
            ref_id = best_hit.db_entry.write_tree_identifier()
            debug_handle.write(f'{ref_id}\t{query_id}\n')

            divergence = matrix[query_id, ref_id]
            if not isinstance(divergence, float):
                raise ValueError('Distance value is not a float.')

            # sequence information for query sequence
            reference_species = reference.taxon_species.scientific_name if not reference.taxon_species is None else 'Reference_unspecified_species'
            
            result: str = ''
            if divergence < threshold:
                # "Correct ID": Query < 1.0/2.0% divergent from reference, and query name matches reference species name
                if query_species == reference_species:
                    result = BlastQuerySequence.QueryClassification.CORRECT_ID
                # "New ID": Query < 1.0/2.0% divergent from reference, and query not labelled to species, e.g. ‘Gymnotiformes’ or ‘sp.’
                elif 'sp.' in query_species or len(query_species.split(' ')) < 2:
                    result = BlastQuerySequence.QueryClassification.NEW_ID
                # "Incorrect ID": Query < 1.0/2.0% divergent from reference, and query name does not match reference species name
                else:
                    result = BlastQuerySequence.QueryClassification.INCORRECT_ID
            else:
                # "Tentative Correct ID : Query > 1.0/2.0% divergent from reference, and most similar reference name matches query species name
                if query_species == reference_species:
                    result = BlastQuerySequence.QueryClassification.TENTATIVE_CORRECT_ID
                # Unknown ID: Query > 1.0/2.0% divergent from reference, and query not labelled to species, e.g. ‘Gymnotiformes’ or ‘sp.’
                elif 'sp.' in query_species or len(query_species.split(' ')) < 2:
                    result = BlastQuerySequence.QueryClassification.UNKNOWN_ID
                # Tentative Additional Species: Query > 1.0/2.0% divergent from reference, and query labelled as a species not included in reference library
                elif query_species not in refs_in_db:
                    result = BlastQuerySequence.QueryClassification.TENATIVE_ADDITIONAL_SPECIES
                # Incorrect ID without Match": Query > 1.0/2.0% divergent from reference, and most similar reference name does not match reference species name
                else:
                    result = BlastQuerySequence.QueryClassification.INCORRECT_ID_NO_MATCH
            class_handle.write(f'{query_id}\t{seq.write_tree_identifier()}\t{query_species}\t{reference_species}\t{result}\n')
            seq.accuracy_category = result
        BlastQuerySequence.objects.bulk_update(seqs, fields=['accuracy_category'])
    finally:
        class_handle.close()
        debug_handle.close()
        return

            


             


    




 