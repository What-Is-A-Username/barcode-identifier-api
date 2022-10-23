# Generate the default database, using the accession numbers found in default.txt

import requests 

db_id = '4f33c746-e566-4cfb-a79d-1d4bcb8cae6d'
url = f'http://127.0.0.1:8000/blastdbs/{db_id}/add/'

with open('./default.txt', 'r') as default_numbers: 
    lines = default_numbers.readlines()
    lines = [l.strip() for l in lines]

lines = lines[:20]

for l in lines:
    post_data = {
        'accession_number': l
    }
    requests.post(url = url, data = post_data)

print(lines)