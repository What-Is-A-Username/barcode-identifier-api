import json
from decimal import Decimal
import re
from typing import Any, Dict, List, Tuple
from barcode_blastn.helper.compare_hits import compareHits

from barcode_blastn.models import Hit, NuccoreSequence

def parse_results(results_string: str) -> Tuple[List[Dict[str, Any]], Dict[str, List[int]]]:
    '''
    Parse the output text of BLAST into a list of hits that can be used to
    bulk create hits, and return a tuple with first element being a list of all hits,
    and the second element being the corresponding indices of the best hits indexed 
    by the query_accession_version.
    '''
    lines = re.split('\r\n|\n|\r', results_string)
    lines = [l.strip() for l in lines if len(l.strip()) > 0]
    hits = [l for l in lines if not l.startswith('#')]

    # The list of all columns in the BLAST output tables, in the order in which they
    # appear
    model_fields = ['query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']
    # These columns are to be converted to integers
    integer_fields = ['alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end']
    # These columns are to be converted to decimals
    decimal_fields = ['percent_identity', 'evalue', 'bit_score']

    def extract_hit(line):
        '''Convert a list of strings into a dictionary'''
        data = dict(zip(model_fields, line.split('\t')))
        for field in decimal_fields:
            data[field] = Decimal(data[field])
        for field in integer_fields:
            data[field] = int(data[field])
        data['best_hit'] = False
        return data

    last_query_accession_version: str = ''
    current_position: int = 0
    parsed_hits: List[Dict[str, Any]] = []
    best_hits: Dict[str, List[int]] = {}
    for index, hit in enumerate(hits):
        hit_data = extract_hit(hit)

        # If the last subject accession (i.e. query) differs in ID, we know we have reached
        # a new table. So reset the hit position to 1, to indicate that its top of the table
        query_accession_version: str = hit_data.get('query_accession_version', '')
        assert len(query_accession_version) > 0
        if query_accession_version != last_query_accession_version:
            last_query_accession_version = query_accession_version
            current_position = 1
        # Same table as before, so increment the position
        else:
            current_position = current_position + 1
        hit_data['position'] = current_position

        # Compare the present hit with the current best hits
        current_best_hits = best_hits.get(query_accession_version, [])
        if len(current_best_hits) == 0:
            best_hits[query_accession_version] = [index]
        else:
            comparison = compareHits(parsed_hits[current_best_hits[0]], hit_data)
            if comparison == 1:
                best_hits[query_accession_version] = [index]
            elif comparison == 0:
                best_hits[query_accession_version].append(index)
         
        parsed_hits.append(hit_data)

    # Now that the best hits are identified, cache the result 
    for key, value in best_hits.items():
        for index in value:
            parsed_hits[index]['best_hit'] = True


    return (parsed_hits, best_hits)