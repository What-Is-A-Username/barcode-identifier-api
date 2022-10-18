import os

# Make an entrez query usable by the NCBI Blast API (https://blast.ncbi.nlm.nih.gov/Blast.cgi), based on the accession numbers present in default.txt
default_db = os.path.abspath('./default.txt')
with open(default_db, 'r') as default:
    p = default.readlines()
p = [f'{an.strip()}[Accession]' for an in p]
full_query = ' OR '.join(p)
out_query = os.path.abspath('./entrez_query.txt')
with open(out_query, 'w') as outf:
    outf.write(full_query)
    outf.close()

