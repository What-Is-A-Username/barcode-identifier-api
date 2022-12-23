from typing import Dict, Generator, List
from Bio import Entrez 
from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature
from urllib.error import HTTPError
from ratelimit import limits
from ratelimit.decorators import sleep_and_retry
from datetime import datetime
from urllib.error import HTTPError

# NCBI Rate limit without an API key is 3 requests per second
PERIOD = 1 # Time between calls, in number of seconds
MAX_ACCESSIONS_PER_REQUEST = 300 # Limit each request to only return a max number of accessions

class InvalidAccessionNumberError(Exception):
    """Raised when accession number cannot be found in the retrieved GenBank XML"""
    pass

@sleep_and_retry
@limits(calls = 1, period = PERIOD)
def retrieve_gb(accession_numbers: List[str]) -> List[Dict]: 
    """
        Retrieve from GenBank the nucleotide entries for the given accession numbers.
        
    Raises:
        ValueError: If no accession numbers are provided, or the number of accessions is more than the maximum specified by MAX_ACCESSIONS_PER_REQUEST
    """
    if accession_numbers is None:
        raise ValueError(f'List of accession numbers to query is None.')
    if len(accession_numbers) == 0:
        raise ValueError(f'List of accession numbers to query is empty.')
    if len(accession_numbers) > MAX_ACCESSIONS_PER_REQUEST:
        raise ValueError(f'Cannot query the {len(accession_numbers)} found in a single request. Limit per request is {MAX_ACCESSIONS_PER_REQUEST}')

    accessions = ','.join(accession_numbers)

    Entrez.email = "william.huang1212@gmail.com"
    Entrez.max_tries = 1
    Entrez.tool = "barcode_identifier"

    request_time_str = datetime.now().strftime("%H:%M:%S.%f")
    print(f'{request_time_str} | Fetching from GenBank for accessions {accessions}')

    try: 
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accessions)
    except HTTPError:
        raise 
    else:
        response_time_str = datetime.now().strftime("%H:%M:%S.%f")
        print(f'{response_time_str} | Response received from GenBank for accessions {accessions}')
    seq_records : Generator = SeqIO.parse(handle, "genbank")

    all_data = []

    seq_record : SeqRecord
    for seq_record in seq_records:
        current_data : Dict = {
            'accession_number': seq_record.name,
            'definition': seq_record.description,
            'dna_sequence': str(seq_record.seq)
        }
        features : List[ SeqFeature ] = seq_record.features
        qualifiers_to_extract = ['organism', 'organelle', 'mol_type', 'isolate', 'country', 'specimen_voucher', 'lat_lon']
        try: 
            source_feature : SeqFeature = [feature for feature in features if feature.type == 'source'][0]
        except IndexError:
            for qualifier_name in qualifiers_to_extract:
                current_data[qualifier_name] = ''
        else:
            for qualifier_name in qualifiers_to_extract:
                if qualifier_name in source_feature.qualifiers:
                    try:
                        current_data[qualifier_name] = source_feature.qualifiers[qualifier_name][0]
                    except IndexError:
                        current_data[qualifier_name] = ''
                else:
                    current_data[qualifier_name] = ''

        all_data.append(current_data)

    handle.close()

    return all_data
