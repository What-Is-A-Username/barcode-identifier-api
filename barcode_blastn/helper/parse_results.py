import json
from decimal import Decimal
import re

from barcode_blastn.models import NuccoreSequence

def parse_results(results_string: str):
    lines = re.split('\r\n|\n|\r', results_string)
    lines = [l.strip() for l in lines if len(l.strip()) > 0]
    hits = [l for l in lines if not l.startswith('#')]

    model_fields = ['query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']
    integer_fields = ['alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end']

    def extract_hit(line):
        data = dict(zip(model_fields, line.split('\t')))
        data['percent_identity'] = Decimal(data['percent_identity'])
        data['evalue'] = Decimal(data['evalue'])
        data['bit_score'] = Decimal(data['bit_score'])
        for name in integer_fields:
            data[name] = int(data[name])
        return data

    parsed_hits = [extract_hit(h) for h in hits]

    # establish relation between hit and sequence
    
       
    # For outfmt = 15
    # Convert all hits to Hit entries and save them
    # hits = results["BlastOutput2"][0]["report"]["results"]["search"]["hits"]
    # query_id = results["BlastOutput2"][0]["report"]["results"]["search"]["query_id"]

    # hits = [{
    #     'query_accession_version': query_id,
    #     'subject_accession_version': h["description"][0]["title"],
    #     'percent_identity': Decimal(h["hsps"][0]["identity"]),
    #     'alignment_length': h["hsps"][0]["align_len"],
    #     'mismatches': 0,
    #     'gap_opens': h["hsps"][0]["gaps"],
    #     'query_start': h["hsps"][0]["query_from"],
    #     'query_end': h["hsps"][0]["query_to"],
    #     'sequence_start': h["hsps"][0]["hit_from"],
    #     'sequence_end': h["hsps"][0]["hit_to"],
    #     'evalue': Decimal(h["hsps"][0]["evalue"]),
    #     'bit_score': Decimal(h["hsps"][0]["bit_score"])
    # } for h in hits]

    print(parsed_hits)

    return parsed_hits