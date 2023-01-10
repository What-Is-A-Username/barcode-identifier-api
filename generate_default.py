# THIS SCRIPT IS OBSELETE:
#       Adding to a database requires admin privileges, which are not given by this script

# # This python script sends a single request to add a list of accession numbers found in default.txt to the specified database. 

# import requests 

# # Set to the ID of the destination database where you want all the sequences to be added to
# db_id = '4f33c746-e566-4cfb-a79d-1d4bcb8cae6d'
# # Url to the server endpoint which adds sequences to the database
# url = f'http://127.0.0.1:8000/blastdbs/{db_id}/add/'

# with open('./default.txt', 'r') as default_numbers: 
#     lines = default_numbers.readlines()
#     lines = [l.strip() for l in lines]

# lines = lines[:20]

# for l in lines:
#     post_data = {
#         'accession_number': l # send a list of accession numbers in the POST request body
#     }
#     requests.post(url = url, data = post_data)

# # print the list of accession numbers sent; only for debug logging
# print(lines)

# import json 

# with open("current_default.json", "r") as my_file:
#     existing = "".join([line.strip() for line in my_file.readlines()])
# json_data = json.loads(existing)
# existing_sequences = [r["accession_number"] for r in json_data["sequences"]]




# print single comma-separated string, no quotes
with open("default.txt", "r") as desired_seqs:
    lines = desired_seqs.readlines()
    desired = [f'{line.strip()}' for line in lines]
print(",".join(desired))

# print single comma-separated string, with quotes
with open("default.txt", "r") as desired_seqs:
    lines = desired_seqs.readlines()
    desired = [f'"{line.strip()}"' for line in lines]
print(",".join(desired))