from typing import Dict
from urllib.error import HTTPError
import urllib.request
import xml.etree.ElementTree as ET
from ratelimit import limits
from ratelimit.decorators import sleep_and_retry
from rest_framework import status
from datetime import datetime

'''
Given a url, return the XML response as a string
'''

ONE_SECOND = 1 # NCBI Rate limit without an API key is 3 requests per second

@sleep_and_retry
@limits(calls = 1, period = ONE_SECOND)
def retrieve_gb(accession_number: str) -> str: 
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession_number}&rettype=gb&retmode=xml'
    
    request_time = datetime.now()
    request_time_str = request_time.strftime("%H:%M:%S.%f")
    print(f'{request_time_str} | Fetching from GenBank at url "{url}"')

    try:
        response = urllib.request.urlopen(url)
    except HTTPError as e:
        # e.msg = f'Failed to fetch entry from GenBank.'
        raise
    
    response_time = datetime.now()
    response_time_str = response_time.strftime("%H:%M:%S.%f")
    print(f'{response_time_str} | Response received from "{url}"')

    xml_string = response.read().decode('UTF8')
    return xml_string

class InvalidAccessionNumberError(Exception):
    """Raised when accession number cannot be found in the retrieved GenBank XML"""
    pass

def parse_gbx_xml(xml_string: str, accession_number: str) -> Dict:
    qualifiers_to_extract = ['organism', 'organelle', 'mol_type', 'isolate', 'country', 'specimen_voucher']
    
    gbset = ET.fromstring(xml_string)
    gbseq = gbset.find(f"./GBSeq[GBSeq_primary-accession='{accession_number}']")
    if gbseq is None:
        raise InvalidAccessionNumberError(f'Could not locate the data for accession number {accession_number} in the GenBank file.')
    gb_source_quals = gbseq.find('./GBSeq_feature-table/GBFeature[GBFeature_key="source"]/GBFeature_quals')

    current_data = {
        'accession_number': accession_number,
        'definition': gbseq.find('./GBSeq_definition').text,
        'dna_sequence': gbseq.find('./GBSeq_sequence').text
    }

    qualifiers = gb_source_quals.findall('./GBQualifier')
    for q in qualifiers:
        qualifier_name = q.find('./GBQualifier_name').text
        if qualifier_name in qualifiers_to_extract:
            current_data[qualifier_name] = q.find('./GBQualifier_value').text
        
    return current_data

    


