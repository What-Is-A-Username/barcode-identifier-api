from typing import Dict, Generator, List
from Bio import Entrez 
from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature
from ratelimit import limits
from ratelimit.decorators import sleep_and_retry
from datetime import datetime

# NCBI Rate limit without an API key is 3 requests per second
PERIOD = 1 # Time between calls, in number of seconds
MAX_ACCESSIONS_PER_REQUEST = 300 # Limit each request to only return a max number of accessions

class InvalidAccessionNumberError(Exception):
    """Raised when accession number cannot be found in the retrieved GenBank XML"""
    pass

@sleep_and_retry
@limits(calls = 1, period = PERIOD)
def retrieve_gb(accession_numbers: List[str]) -> List[Dict[str, str]]: 
    """
        Retrieve from GenBank the nucleotide entries for the given accession numbers.
        
    Raises:
        ValueError: If no accession numbers are provided, or the number of accessions is more than the maximum specified by MAX_ACCESSIONS_PER_REQUEST
        HTTPError: If there is a network error when retrieving data from GenBank using Bio.Entrez
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

    # Raises HTTPError if there is a network error
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accessions)
    response_time_str = datetime.now().strftime("%H:%M:%S.%f")
    print(f'{response_time_str} | Response received from GenBank for accessions {accessions}')
    seq_records : Generator = SeqIO.parse(handle, "genbank")

    all_data: List[Dict[str, str]]= []

    seq_record : SeqRecord
    for seq_record in seq_records:
        current_data : Dict = {
            'accession_number': seq_record.name,
            'definition': seq_record.description,
            'dna_sequence': str(seq_record.seq)
        }
        features : List[ SeqFeature ] = seq_record.features
        qualifiers_to_extract = ['organism', 'organelle', 'isolate', 'country', 'specimen_voucher', 'type_material', 'lat_lon']
        try: 
            source_feature : SeqFeature = [feature for feature in features if feature.type == 'source'][0]
        except IndexError:
            # set fields to N/A if no source features found
            for qualifier_name in qualifiers_to_extract:
                current_data[qualifier_name] = 'N/A'
        else:
            for qualifier_name in qualifiers_to_extract:
                if qualifier_name in source_feature.qualifiers:
                    try:
                        # take only the first line/element of the feature
                        current_data[qualifier_name] = source_feature.qualifiers[qualifier_name][0]
                    except IndexError:
                        # set it to 'error' if above code errored
                        current_data[qualifier_name] = 'error'
                else:
                    # set it to empty string if qualifier was not found in the data
                    current_data[qualifier_name] = ''

            # use type material specified in 'note' if no type_material was found 
            if (current_data['type_material'] == '' and 'note' in source_feature.qualifiers):
                notes : str = ''
                try:
                    # include all notes 
                    notes = "\n".join(source_feature.qualifiers['note'])
                except BaseException:
                    notes = 'error'
                else:
                    notes_lower = notes.lower()
                    if 'paratype' in notes_lower or 'holotype' in notes_lower:
                        if notes_lower.startswith('type: ') and len(notes_lower) > 6:
                            # remove "type: " from the beginning if it is present
                            notes = notes[6:]
                        else:
                            # print a warning to the console if the "type: " string was not found at start
                            print(f'Inferred type material from notes since it contained "paratype" or "holotype". It did not start with "type: ". Consider checking /type_material and/or /note info for {seq_record.name} in GenBank.')
                finally:
                    current_data['type_material'] = notes
        finally:
            all_data.append(current_data)

    handle.close()

    return all_data
