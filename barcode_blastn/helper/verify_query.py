'''
Given a fasta definition and a sequence, verify if the input is valid
'''
def verify_query(definition: str, sequence: str):
    if definition is None or sequence is None:
        return False
    # TODO: Add constant indicating the minimum sequence length processible
    elif len(definition) < 1 or len(sequence) < 6:
        return False 
    # TODO: Check for invalid characters in definition
    # TODO: Check for invalid definition format (e.g. extra >)
    elif definition[0].startswith('>'):
        return False
    else:
        return True