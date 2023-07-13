from datetime import datetime
from typing import List

from barcode_blastn.file_paths import get_data_run_path
from barcode_blastn.models import (BlastQuerySequence, BlastRun,
                                   Hit, NuccoreSequence)
from Bio.Phylo.TreeConstruction import DistanceMatrix
from django.db.models import QuerySet


def parse_taxonomic_abbreviations(source_organism: str) -> str:
    '''
    Given name of the source organism, remove:
    -   cf. 
    -   aff.

    Examples:
        - 'G cylindricus' => 'G cylindricus'
        - 'G aff. cylindricus' => 'G cylindricus'
        - 'G cf. cylindricus' => 'G cylindricus'
        - 'aff. G cylindricus' => 'G cylindricus'
    '''
    delim = source_organism.split(' ')
    exclude = ['cf.', '.aff']
    new_delim = [d for d in delim if d not in exclude]
    return ' '.join(new_delim)

def annotate_accuracy_category(matrix: DistanceMatrix, run: BlastRun, threshold: float = 0.01) -> bool:
    '''
    Annotate each query sequence in the matrix using the categorization proposed in
    Janzen et al. 2022. Return True if operation was successful, False otherwise.
    
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

    function_result = False

    try:
        for seq in seqs:
            # Find the highest rated hits and store them in best hits
            best_hits: List[Hit] = []
            hit: Hit
            best_hits = seq.best_hits()
            
            # Extract the original species name from the sequence information
            query_id = seq.write_tree_identifier()
            query_species = seq.original_species_name if not seq.original_species_name is None else ''
                
            # If no hits are returned, terminate early and provide a classification
            if len(best_hits) == 0:
                seq.accuracy_category = BlastQuerySequence.QueryClassification.NO_HITS
                class_handle.write(f'{query_id}\t{seq.write_tree_identifier()}\t{query_species}\tNo hits\t{seq.accuracy_category}\n')
                continue
            
            # Create a list, ref_species, which is a list of species names corresponding
            # to the best hits
            references = [b.db_entry for b in best_hits]
            ref_ids = [h.write_tree_identifier() for h in references]
            def extract_ref_species(entry: NuccoreSequence) -> str:
                return entry.taxon_species.scientific_name if not entry.taxon_species is None else 'Reference_unspecified_species'
            ref_species = [extract_ref_species(h) for h in references]
            
            try:
                # If the query species name is one of the best hits, only consider that hit
                index = ref_species.index(query_species)
                reference: NuccoreSequence = references[index]
            except ValueError:
                # If the query species name is missing, just take the very first hit 
                reference: NuccoreSequence = references[0]

            ref_id = reference.write_tree_identifier()
            debug_handle.write(f'{ref_id}\t{query_id}\n')

            # Retrieve the calculated distance between the query and the reference
            divergence = matrix[query_id, ref_id]
            if not isinstance(divergence, float):
                raise ValueError('Distance value is not a float.')

            # sequence information for query sequence
            reference_species = reference.taxon_species.scientific_name if not reference.taxon_species is None else 'Reference_unspecified_species'
            
            query_species = parse_taxonomic_abbreviations(query_species)

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
        function_result = True
    except BaseException as err:
        run.errors = run.errors + '\nErrored while annotating taxonomic assignments.'
        run.status = BlastRun.JobStatus.ERRORED
        run.error_time = datetime.now()
        run.save()
        function_result = False
    finally:
        class_handle.close()
        debug_handle.close()
        return function_result
