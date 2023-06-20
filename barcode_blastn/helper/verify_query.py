import re

MIN_SEQUENCE_LENGTH = 12
MIN_HEADER_LENGTH = 1

def verify_header(header: str, sequence: str):
    if header is None:
        raise ValueError('Definition missing.')
    elif len(header) < MIN_HEADER_LENGTH:
        raise ValueError(f'Header must be at least {MIN_HEADER_LENGTH} characters in length')  
    elif not bool(re.match('^[A-Za-z0-9.:#*_-]+$', header)):
        raise ValueError(f"The sequence identifier in the header '{header}' has invalid characters. Only letters, digits, hyphens, underscores, periods, colons, asterisks and number signs are allowed.")
    return True

def verify_dna(header: str, sequence: str):
    if sequence is None:
        raise ValueError('Sequence missing.')
    elif len(sequence) < MIN_SEQUENCE_LENGTH:
        raise ValueError(f'Sequence must be at least {MIN_SEQUENCE_LENGTH} nucleotides in length')
    elif not bool(re.match('^[RYSWKMBDHVNCAGTryswkmbdhvncagt.-]+$', sequence)):
        raise ValueError(f"The sequence provided for '>{header}' has non-DNA nucleotide characters.")
    return True