import time
from typing import Any, Dict, List, Union
from django.test import TestCase
from rest_framework import status
from dateutil.parser import parse
from datetime import datetime
from django.contrib.auth.models import User
from barcode_blastn.models import BlastDb, BlastRun, Library
from copy import deepcopy
import json

COMMON_DATE = datetime(2023, 1, 1, 0, 0, 0, 0)
USERNAME = 'John Doe'
PASSWORD = '1234FakePassword!!!!'
EMAIL = 'johndoe@email.com'

def remove_id(obj: Union[Dict[str, Any], List]):
    if isinstance(obj, dict):
        if 'id' in obj:
            obj['id'] = ''
        for key, value in obj.items():
            if isinstance(value, dict):
                obj[key] = remove_id(value)
            elif isinstance(value, list):
                obj[key] = [remove_id(o) for o in value]
            elif isinstance(value, str):
                try:
                    parse(obj[key])
                except ValueError:
                    pass
                else:
                    # Set all dates and times to some common time
                    obj[key] = str(COMMON_DATE)
    elif isinstance(obj, list):
        obj = [remove_id(o) for o in obj]
    return obj

class AuthTestCase(TestCase):
    '''Base class for all tests requiring a user with superuser and staff permissions.'''

    owner_serializer = {
        "username": USERNAME,
        "email": EMAIL
    }

    def setUp(self) -> None:
        self.user = User.objects.create_user(USERNAME, password=PASSWORD, email=EMAIL, is_superuser=True, is_staff=True)
        self.user.save()
        self.maxDiff = None
        self.client.login(username=USERNAME, password=PASSWORD)

    def assertEqualIgnoreId(self, response, expected):
        print('ACTUAL:\n')
        print(response)
        print('EXPECTED:\n')
        print(expected)
        cleaned_response = remove_id(response)
        cleaned_expected = remove_id(expected)
        self.assertEqual(cleaned_response, cleaned_expected)

    def tearDown(self) -> None:
        super().tearDown()
        self.client.logout()

class LibraryListTest(AuthTestCase):
    '''
    Setup: Two libraries A and B. A has two databases.
    '''

    # Data to add a reference library thru a POST request. This library will also have databases added
    dummy_library_1 = {
        'custom_name': 'Library of South American Fishes',
        'description': 'A collection of voucher specimens from the literature.',
        'public': True
    }

    dummy_library_2 = {
        'custom_name': 'Library of European Fishes',
        'description': 'A collection of European voucher specimens from the literature.',
        'public': True
    }

    dummy_db_1 = {
        'custom_name': 'January Release',
        'description': 'Sequences gathered as of January 2022',
        'locked': True,
        'accession_numbers': ['ON303390.1']
    }

    dummy_db_2 = {
        'custom_name': 'February Release',
        'description': 'Sequences gathered as of February 2022',
        'locked': True,
        'accession_numbers': ['ON303390.1', 'ON303391.1']
    }

    dummy_db_3 = {
        'custom_name': 'March Release',
        'description': 'Sequences gathered as of March 2022',
        'locked': False,
        'accession_numbers': ['ON303390.1', 'ON303391.1', 'ON303392.1']
    }

    def setUp(self) -> None:
        super().setUp()
        # Create libraries
        library_1 = self.client.post('/libraries/', self.dummy_library_1, content_type='application/json')
        library_2 = self.client.post('/libraries/', self.dummy_library_2, content_type='application/json')
        self.assertEqual(library_1.status_code, status.HTTP_201_CREATED)
        self.assertEqual(library_2.status_code, status.HTTP_201_CREATED)
        self.library = json.loads(library_1.content)
        self.library_id = self.library['id']
        # Create databases
        time.sleep(1)
        db_1 = self.client.post(f'/libraries/{self.library_id}/versions/', self.dummy_db_1, content_type='application/json')
        self.assertEqual(BlastDb.objects.all().count(), 1)
        db_2 = self.client.post(f'/libraries/{self.library_id}/versions/', self.dummy_db_2, content_type='application/json')
        self.assertEqual(BlastDb.objects.all().count(), 2)
        time.sleep(1)
        db_3 = self.client.post(f'/libraries/{self.library_id}/versions/', self.dummy_db_3, content_type='application/json')
        self.assertEqual(BlastDb.objects.all().count(), 3)
       
        self.assertEqual(db_1.status_code, status.HTTP_201_CREATED)
        self.assertEqual(db_2.status_code, status.HTTP_201_CREATED)
        self.assertEqual(db_3.status_code, status.HTTP_201_CREATED)
        self.db = json.loads(db_2.content)
        self.db_id = self.db['id']
        self.db_3 = json.loads(db_3.content)

    post_libraries_201_request = {
        'custom_name': 'Newly Sequenced Species Reference Library',
        'description': 'A collection of new sequences from several species of interest.',
        'public': True
    }

    post_libraries_201_response = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "owner": AuthTestCase.owner_serializer,
            "custom_name": "Newly Sequenced Species Reference Library",
            "description": "A collection of new sequences from several species of interest.",
            "public": True,
            "latest": None
        }

    def test_post_libraries(self):
        response = self.client.post('/libraries/', self.post_libraries_201_request, content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqualIgnoreId(json.loads(response.content), self.post_libraries_201_response)
        self.assertEqual(Library.objects.count(), 3)

    get_libraries_200_response = [
        {
            "owner": AuthTestCase.owner_serializer,
            "id": '66855f2c-f360-4ad9-8c98-998ecb81fa2v',
            'custom_name': 'Library of European Fishes',
            'description': 'A collection of European voucher specimens from the literature.',
            'public': True,
            "latest": None
        },
        {
            "owner": AuthTestCase.owner_serializer,
            "id": '66855f2c-f360-4ad9-8c98-998ecb812222',
            'custom_name': 'Library of South American Fishes',
            'description': 'A collection of voucher specimens from the literature.',
            'public': True,
            "latest": {
                "id": "d83fc783-0a61-4b2f-94a8-d04d094e2a6c",
                "description": "Sequences gathered as of February 2022",
                "version_number": "2.1.1",
                "custom_name": "February Release",
                "sequence_count": 2,
                "created": "2023-06-26T20:04:43.113904Z",
                "locked": True
            }
        },
    ]

    def test_get_libraries(self):
        response = self.client.get('/libraries/', format='application/json')
        self.assertEqual(len(json.loads(response.content)), 2)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqualIgnoreId(json.loads(response.content), self.get_libraries_200_response)

    get_libraries_id_200_response = [r for r in get_libraries_200_response if r['custom_name'] == 'Library of South American Fishes'][0]

    def test_get_libraries_id(self):
        response = self.client.get(f'/libraries/{self.library_id}', content_type='application/json')

        expected = deepcopy(self.get_libraries_id_200_response)
        expected['id'] = self.library_id

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqualIgnoreId(json.loads(response.content), expected)

    patch_libraries_id_204_request = {
        "custom_name": "Updated title for library of South American Fishes",
        "description": "Updated description describing the library.",
        "public": False
    }

    patch_libraries_id_204_response = {
        "custom_name": "Updated title for library of South American Fishes",
        "description": "Updated description describing the library.",
        "public": False
    }

    def test_patch_libraries_id(self):
        response = self.client.patch(f'/libraries/{self.library_id}', self.patch_libraries_id_204_request, content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        response = self.client.get(f'/libraries/{self.library_id}', content_type='application/json')
        data = json.loads(response.content) 
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqualIgnoreId(data, self.patch_libraries_id_204_response)

    def test_delete_libraries_id(self):
        response = self.client.delete(f'/libraries/{self.library_id}', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_204_NO_CONTENT)

        response = self.client.delete(f'/libraries/{self.library_id}', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)
        
        response = self.client.get(f'/libraries/{self.library_id}', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)

    get_libraries_id_versions_200_response = [
        {
            "id": "d83fc783-0a61-4b2f-94a8-d04d094e2a6c",
            "version_number": "2.1.1",
            "custom_name": "February Release",
            "sequence_count": 2,
            "description": "Sequences gathered as of February 2022",
            "locked": True
        },
        {
            "id": "5f377bd9-6a98-4071-ab13-3197077815db",
            "version_number": "1.1.1",
            "custom_name": "January Release",
            "sequence_count": 1,
            "description": "Sequences gathered as of January 2022",
            "locked": True
        },
        {
            "id": "71f5d9f3-6df2-4b49-bad7-1a388e928e95",
            "version_number": "0.0.0",
            "custom_name": "March Release",
            "sequence_count": 3,
            "description": "Sequences gathered as of March 2022",
            "locked": False
        }
    ]

    def test_get_libraries_id_versions(self):
        '''Test retrieving list of two existing versions'''
        response = self.client.get(f'/libraries/{self.library_id}/versions/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = json.loads(response.content)
        self.assertEqual(len(data), 3)
        self.assertEqualIgnoreId(data, self.get_libraries_id_versions_200_response)

    get_blastdbs_id_200_response = {
        "id": "d83fc783-0a61-4b2f-94a8-d04d094e2a6c",
        "library": {
            "id": "ff5da812-1c0e-4682-9423-f08dbbc8bc5f",
            "custom_name": "Library of South American Fishes",
            "description": "A collection of voucher specimens from the literature.",
            "public": True,
            "owner": AuthTestCase.owner_serializer
        },
        "custom_name": "February Release",
        "version_number": "2.1.1",
        "description": "Sequences gathered as of February 2022",
        "locked": True,
        "sequences": [
            {
            "id": "232e3b26-f347-4b09-88aa-2dcb52080d80",
            "accession_number": "ON303390",
            "version": "ON303390.1",
            "organism": "Gymnotus cylindricus",
            "organelle": "mitochondrion",
            "isolate": "2094",
            "country": "Costa Rica",
            "specimen_voucher": "ROM:84772",
            "dna_sequence": "ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTCAACCCGGAGCCCTCCTCGGGGACGACCAAATTTATAATGTAATTGTTACTGCCCACGCTTTCGTAATAATTTTTTTTATAGTAATACCCATTATAATTGGAGGCTTCGGAAATTGATTAACTCCACTAATAATTGGAGCCCCAGACATAGCATTTCCCCGAATAAATAATATAAGCTTTTGACTCCTTCCTCCTTCTTTTTTACTTCTCCTTGCATCATCCGGAGTTGAAGCAGGGGCCGGAACAGGCTGAACAGTATACCCCCCTCTTGCAGGCAATCTTGCCCATGCAGGAGCCTCAGTAGATCTAACTATTTTCTCTCTTCATCTAGCCGGAGTTTCTTCAATTCTAGGATCCATTAACTTTATTACCACAATTATTAATATAAAACCTCCAGCCATTTCTCAATATCAAACCCCACTATTTATTTGATCACTTCTAGTAACCACTGTCCTTTTACTCCTCTCTCTTCCAGTACTAGCTGCTGGTATCACCATACTACTAACAGATCGAAACTTAAACACAACATTCTTTGACCCGGCGGGCGGAGGAGATCCTATTTTATATCAACATTTA",
            "lat_lon": "10.54 N 83.50 W",
            "type_material": "",
            "created": "2023-06-26T20:04:45.016196Z",
            "updated": "2023-06-26T20:04:45.016200Z",
            "genbank_modification_date": "2022-07-04",
            "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
            "taxon_superkingdom": {
                "id": 2759,
                "scientific_name": "Eukaryota"
            },
            "taxon_kingdom": {
                "id": 33208,
                "scientific_name": "Metazoa"
            },
            "taxon_phylum": {
                "id": 7711,
                "scientific_name": "Chordata"
            },
            "taxon_class": {
                "id": 186623,
                "scientific_name": "Actinopteri"
            },
            "taxon_order": {
                "id": 8002,
                "scientific_name": "Gymnotiformes"
            },
            "taxon_family": {
                "id": 30771,
                "scientific_name": "Gymnotidae"
            },
            "taxon_genus": {
                "id": 36670,
                "scientific_name": "Gymnotus"
            },
            "taxon_species": {
                "id": 699532,
                "scientific_name": "Gymnotus cylindricus"
            },
            "annotations": []
            },
            {
            "id": "221e2c35-52f8-43bf-a21f-30c15ccf36c3",
            "accession_number": "ON303391",
            "version": "ON303391.1",
            "organism": "Gymnotus esmeraldas",
            "organelle": "mitochondrion",
            "isolate": "10865",
            "country": "Ecuador",
            "specimen_voucher": "ZOO.A.V.Pe0310",
            "dna_sequence": "ATAGTATTTGGTGCCTGAGCTGGAATAGTTGGCACAGCCTTGAGCCTACTGATCCGAGCAGAACTAAGCCAACCCGGAACCCTCCTAGGCGATGACCAAATTTATAATGTAATCGTTACTGCCCACGCCTTCGTAATGATTTTCTTTATAGTTATACCTATTATGATTGGAGGCTTCGGGAACTGATTAATTCCACTAATAATTGGGGCCCCAGACATAGCATTCCCCCGAATAAATAATATAAGCTTTTGACTTCTCCCCCCTTCTTTTTTACTTCTCCTTGCTTCATCCGGAGTCGAGGCCGGAGCCGGAACAGGCTGAACCGTATATCCGCCTCTTGCAAGCAATCTTGCTCACGCAGGAGCTTCAGTAGATTTGGCTATTTTTTCACTACACCTCGCCGGAATCTCCTCAATCCTCGGGTCTATTAACTTTATTACCACAATTATTAACATAAAACCCCCAGCCATTACCCAATACCAAACCCCCCTATTTATCTGAGCACTTCTAGTGACTACCGTCCTCCTACTCCTCTCTCTTCCGGTACTGGCTGCCGGCATTACTATGCTATTAACAGACCGAAACCTAAACACAACTTTCTTTGACCCTGCAGGCGGGGGAGACCCTATCTTATACCAACACTTA",
            "lat_lon": "3.06 S 79.74 W",
            "type_material": "",
            "created": "2023-06-26T20:04:45.016134Z",
            "updated": "2023-06-26T20:04:45.016146Z",
            "genbank_modification_date": "2022-07-04",
            "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
            "taxon_superkingdom": {
                "id": 2759,
                "scientific_name": "Eukaryota"
            },
            "taxon_kingdom": {
                "id": 33208,
                "scientific_name": "Metazoa"
            },
            "taxon_phylum": {
                "id": 7711,
                "scientific_name": "Chordata"
            },
            "taxon_class": {
                "id": 186623,
                "scientific_name": "Actinopteri"
            },
            "taxon_order": {
                "id": 8002,
                "scientific_name": "Gymnotiformes"
            },
            "taxon_family": {
                "id": 30771,
                "scientific_name": "Gymnotidae"
            },
            "taxon_genus": {
                "id": 36670,
                "scientific_name": "Gymnotus"
            },
            "taxon_species": {
                "id": 2594964,
                "scientific_name": "Gymnotus esmeraldas"
            },
            "annotations": []
            }
        ]
    }

    def test_get_blastdbs_id(self):
        '''Retrieve a BLAST database by ID'''
        response = self.client.get(f'/blastdbs/{self.db_id}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = json.loads(response.content) 
        self.assertEqualIgnoreId(data, self.get_blastdbs_id_200_response)

    patch_blastdbs_id_204_request = {
        "custom_name": "New Updated Title",
        "description": "This is a new description replacing the old one.",
    }

    patch_blastdbs_id_204_response = {
        "custom_name": "New Updated Title",
        "description": "This is a new description replacing the old one.",
        "locked": False
    }

    def test_patch_blastdbs_id(self):
        '''Patch a BLAST database by ID'''
        response = self.client.patch(f'/blastdbs/{self.db_3["id"]}/', self.patch_blastdbs_id_204_request, content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_204_NO_CONTENT)
        response = self.client.get(f'/blastdbs/{self.db_3["id"]}/', content_type='application/json')
        data = json.loads(response.content) 
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqualIgnoreId(data, self.patch_blastdbs_id_204_response)

    def test_delete_blastdbs_id(self):
        '''Delete an existing blastdb'''
        response = self.client.delete(f'/blastdbs/{self.db_id}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_204_NO_CONTENT)

        response = self.client.delete(f'/blastdbs/{self.db_id}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)

        response = self.client.get(f'/blastdbs/{self.db_id}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)

    post_libraries_id_versions_201_request = {
        'custom_name': 'Interim March Release',
        'description': 'February release augmented with two additional sequences.',
        'locked': False,
        'accession_numbers': ['ON303380.1', 'ON303381.1'],
        'base': 'd83fc783-0a61-4b2f-94a8-d04d094e2a6c' 
    }

    post_libraries_id_versions_201_response = {
        "id": "b6b92be2-d219-4f09-b14b-ed904086f5f0",
        "library": {
            "id": "ff5da812-1c0e-4682-9423-f08dbbc8bc5f",
            "custom_name": "Library of South American Fishes",
            "description": "A collection of voucher specimens from the literature.",
            "public": True,
            "owner": AuthTestCase.owner_serializer
        },
        "custom_name": "Interim March Release",
        "version_number": "0.0.0",
        "description": "February release augmented with two additional sequences.",
        "locked": False,
        "sequences": [
            {
                "id": "8a1e0ff9-2305-4069-83a0-a28337b759da",
                "accession_number": "ON303380",
                "version": "ON303380.1",
                "organism": "Gymnotus ardilai",
                "organelle": "mitochondrion",
                "isolate": "8186",
                "country": "Colombia",
                "specimen_voucher": "IAvHP 11510",
                "dna_sequence": "ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAATTAAGTCAACCCGGGGCCCTCCTAGGTGATGATCAAATTTATAATGTAATCGTTACTGCCCACGCTTTCGTAATGATTTTCTTTATGGTAATACCTATCATGATTGGAGGTTTCGGAAACTGATTAATTCCCCTAATAATTGGAGCCCCAGACATAGCATTTCCCCGAATAAATAACATAAGCTTCTGACTTCTCCCCCCTTCTTTTTTACTTCTCCTTGCATCATCCGGAGTTGAAGCTGGAGCTGGAACAGGCTGAACAGTATACCCCCCTCTCGCAGGTAATCTTGCCCATGCAGGTGCTTCAGTAGACTTAACTATCTTCTCCCTTCATCTAGCCGGAGTTTCCTCTATTCTAGGGTCTATTAATTTTATTACCACAATTATTAATATGAAACCTCCAGCTATTTCTCAATATCAAACCCCATTATTTATTTGAGCACTTCTAATTACCACTGTTCTTTTACTCCTTTCCCTTCCCGTGCTAGCCGCTGGAATTACTATGTTATTAACAGATCGAAACTTAAACACAACCTTCTTTGATCCGGCAGGCGGGGGAGATCCCATTCTCTATCAACACTTA",
                "lat_lon": "7.02 N 73.16 W",
                "type_material": "",
                "created": "2023-06-26T20:11:41.619669Z",
                "updated": "2023-06-26T20:11:41.619673Z",
                "genbank_modification_date": "2022-07-04",
                "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
                "taxon_superkingdom": {
                    "id": 2759,
                    "scientific_name": "Eukaryota"
                },
                "taxon_kingdom": {
                    "id": 33208,
                    "scientific_name": "Metazoa"
                },
                "taxon_phylum": {
                    "id": 7711,
                    "scientific_name": "Chordata"
                },
                "taxon_class": {
                    "id": 186623,
                    "scientific_name": "Actinopteri"
                },
                "taxon_order": {
                    "id": 8002,
                    "scientific_name": "Gymnotiformes"
                },
                "taxon_family": {
                    "id": 30771,
                    "scientific_name": "Gymnotidae"
                },
                "taxon_genus": {
                    "id": 36670,
                    "scientific_name": "Gymnotus"
                },
                "taxon_species": {
                    "id": 2138147,
                    "scientific_name": "Gymnotus ardilai"
                },
                "annotations": []
            },
            {
                "id": "7379bb73-9ea2-4532-824b-022b6ae9b2ac",
                "accession_number": "ON303381",
                "version": "ON303381.1",
                "organism": "Gymnotus bahianus",
                "organelle": "mitochondrion",
                "isolate": "7244",
                "country": "Brazil",
                "specimen_voucher": "MZUSP:102898",
                "dna_sequence": "ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGCCAACCCGGAGCCCTCCTAGGCGATGATCAAATTTATAATGTAATTGTTACTGCCCACGCTTTCGTAATGATTTTCTTTATAGTAATGCCTATCATGATTGGAGGCTTCGGAAACTGATTAATTCCACTAATAATTGGAGCTCCAGACATAGCATTTCCCCGAATAAATAACATAAGCTTCTGACTTCTTCCCCCTTCTTTTTTACTTCTCCTTGCATCATCCGGAGTTGAAGCCGGAGCTGGAACAGGCTGAACAGTATACCCCCCTCTCGCAGGCAATCTTGCCCATGCAGGTGCTTCAGTAGACCTAACTATCTTCTCCCTTCACCTAGCCGGAGTTTCCTCTATTCTAGGGTCTATTAATTTTATTACCACAATTATTAATATGAAACCTCCAGCCATTTCTCAGTATCAAACCCCATTATTTATTTGAGCCCTTCTAATCACTACTGTCCTTTTACTCCTCTCCCTTCCCGTATTAGCCGCTGGAATTACTATGTTATTAACAGATCGAAACTTAAACACAACCTTCTTTGACCCAGCAGGCGGTGGAGACCCCATTCTCTACCAACACTTA",
                "lat_lon": "14.76 S 39.09 W",
                "type_material": "",
                "created": "2023-06-26T20:11:41.619506Z",
                "updated": "2023-06-26T20:11:41.619517Z",
                "genbank_modification_date": "2022-07-04",
                "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
                "taxon_superkingdom": {
                    "id": 2759,
                    "scientific_name": "Eukaryota"
                },
                "taxon_kingdom": {
                    "id": 33208,
                    "scientific_name": "Metazoa"
                },
                "taxon_phylum": {
                    "id": 7711,
                    "scientific_name": "Chordata"
                },
                "taxon_class": {
                    "id": 186623,
                    "scientific_name": "Actinopteri"
                },
                "taxon_order": {
                    "id": 8002,
                    "scientific_name": "Gymnotiformes"
                },
                "taxon_family": {
                    "id": 30771,
                    "scientific_name": "Gymnotidae"
                },
                "taxon_genus": {
                    "id": 36670,
                    "scientific_name": "Gymnotus"
                },
                "taxon_species": {
                    "id": 2138148,
                    "scientific_name": "Gymnotus bahianus"
                },
                "annotations": []
            },
            {
                "id": "2036bde1-0ce0-4ea8-a58d-bf6e6a896e6a",
                "accession_number": "ON303390",
                "version": "ON303390.1",
                "organism": "Gymnotus cylindricus",
                "organelle": "mitochondrion",
                "isolate": "2094",
                "country": "Costa Rica",
                "specimen_voucher": "ROM:84772",
                "dna_sequence": "ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTCAACCCGGAGCCCTCCTCGGGGACGACCAAATTTATAATGTAATTGTTACTGCCCACGCTTTCGTAATAATTTTTTTTATAGTAATACCCATTATAATTGGAGGCTTCGGAAATTGATTAACTCCACTAATAATTGGAGCCCCAGACATAGCATTTCCCCGAATAAATAATATAAGCTTTTGACTCCTTCCTCCTTCTTTTTTACTTCTCCTTGCATCATCCGGAGTTGAAGCAGGGGCCGGAACAGGCTGAACAGTATACCCCCCTCTTGCAGGCAATCTTGCCCATGCAGGAGCCTCAGTAGATCTAACTATTTTCTCTCTTCATCTAGCCGGAGTTTCTTCAATTCTAGGATCCATTAACTTTATTACCACAATTATTAATATAAAACCTCCAGCCATTTCTCAATATCAAACCCCACTATTTATTTGATCACTTCTAGTAACCACTGTCCTTTTACTCCTCTCTCTTCCAGTACTAGCTGCTGGTATCACCATACTACTAACAGATCGAAACTTAAACACAACATTCTTTGACCCGGCGGGCGGAGGAGATCCTATTTTATATCAACATTTA",
                "lat_lon": "10.54 N 83.50 W",
                "type_material": "",
                "created": "2023-06-26T20:11:41.619619Z",
                "updated": "2023-06-26T20:11:41.619623Z",
                "genbank_modification_date": "2022-07-04",
                "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
                "taxon_superkingdom": {
                    "id": 2759,
                    "scientific_name": "Eukaryota"
                },
                "taxon_kingdom": {
                    "id": 33208,
                    "scientific_name": "Metazoa"
                },
                "taxon_phylum": {
                    "id": 7711,
                    "scientific_name": "Chordata"
                },
                "taxon_class": {
                    "id": 186623,
                    "scientific_name": "Actinopteri"
                },
                "taxon_order": {
                    "id": 8002,
                    "scientific_name": "Gymnotiformes"
                },
                "taxon_family": {
                    "id": 30771,
                    "scientific_name": "Gymnotidae"
                },
                "taxon_genus": {
                    "id": 36670,
                    "scientific_name": "Gymnotus"
                },
                "taxon_species": {
                    "id": 699532,
                    "scientific_name": "Gymnotus cylindricus"
                },
                "annotations": []
            },
            {
                "id": "5164ddb1-d955-4e7f-94af-df2a87b7af93",
                "accession_number": "ON303391",
                "version": "ON303391.1",
                "organism": "Gymnotus esmeraldas",
                "organelle": "mitochondrion",
                "isolate": "10865",
                "country": "Ecuador",
                "specimen_voucher": "ZOO.A.V.Pe0310",
                "dna_sequence": "ATAGTATTTGGTGCCTGAGCTGGAATAGTTGGCACAGCCTTGAGCCTACTGATCCGAGCAGAACTAAGCCAACCCGGAACCCTCCTAGGCGATGACCAAATTTATAATGTAATCGTTACTGCCCACGCCTTCGTAATGATTTTCTTTATAGTTATACCTATTATGATTGGAGGCTTCGGGAACTGATTAATTCCACTAATAATTGGGGCCCCAGACATAGCATTCCCCCGAATAAATAATATAAGCTTTTGACTTCTCCCCCCTTCTTTTTTACTTCTCCTTGCTTCATCCGGAGTCGAGGCCGGAGCCGGAACAGGCTGAACCGTATATCCGCCTCTTGCAAGCAATCTTGCTCACGCAGGAGCTTCAGTAGATTTGGCTATTTTTTCACTACACCTCGCCGGAATCTCCTCAATCCTCGGGTCTATTAACTTTATTACCACAATTATTAACATAAAACCCCCAGCCATTACCCAATACCAAACCCCCCTATTTATCTGAGCACTTCTAGTGACTACCGTCCTCCTACTCCTCTCTCTTCCGGTACTGGCTGCCGGCATTACTATGCTATTAACAGACCGAAACCTAAACACAACTTTCTTTGACCCTGCAGGCGGGGGAGACCCTATCTTATACCAACACTTA",
                "lat_lon": "3.06 S 79.74 W",
                "type_material": "",
                "created": "2023-06-26T20:11:41.619567Z",
                "updated": "2023-06-26T20:11:41.619572Z",
                "genbank_modification_date": "2022-07-04",
                "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
                "taxon_superkingdom": {
                    "id": 2759,
                    "scientific_name": "Eukaryota"
                },
                "taxon_kingdom": {
                    "id": 33208,
                    "scientific_name": "Metazoa"
                },
                "taxon_phylum": {
                    "id": 7711,
                    "scientific_name": "Chordata"
                },
                "taxon_class": {
                    "id": 186623,
                    "scientific_name": "Actinopteri"
                },
                "taxon_order": {
                    "id": 8002,
                    "scientific_name": "Gymnotiformes"
                },
                "taxon_family": {
                    "id": 30771,
                    "scientific_name": "Gymnotidae"
                },
                "taxon_genus": {
                    "id": 36670,
                    "scientific_name": "Gymnotus"
                },
                "taxon_species": {
                    "id": 2594964,
                    "scientific_name": "Gymnotus esmeraldas"
                },
                "annotations": []
            }
        ]
    }

    def test_post_library_id_versions(self):
        '''Add BLAST database to one of the libraries'''
        req = self.post_libraries_id_versions_201_request
        req['base'] = self.db_id
        response = self.client.post(f'/libraries/{self.library_id}/versions/', self.post_libraries_id_versions_201_request)

        expected = deepcopy(self.post_libraries_id_versions_201_response)
        data = json.loads(response.content)
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(self.library_id, data['library']['id'])
        self.assertEqualIgnoreId(data, expected)

class SequenceTester(AuthTestCase):

    # Data to add a reference library thru a POST request. This library will also have databases added
    dummy_library_1 = {
        'custom_name': 'Library of South American Fishes',
        'description': 'A collection of voucher specimens from the literature.',
        'public': True
    }

    dummy_db_1 = {
        'custom_name': 'January Release',
        'description': 'Sequences gathered as of January 2022',
        'locked': True,
        'accession_numbers': ['ON303390.1']
    }

    dummy_db_2 = {
        'custom_name': 'February Release',
        'description': 'A database currently being curated.',
        'public': False,
        'accession_numbers': []
    }

    def setUp(self) -> None:
        super().setUp()

        library_1 = self.client.post('/libraries/', self.dummy_library_1, content_type='application/json')
        self.assertEqual(library_1.status_code, status.HTTP_201_CREATED)
        self.library = json.loads(library_1.content)
        self.library_id = self.library['id']

        locked_db = self.client.post(f'/libraries/{self.library_id}/versions/', self.dummy_db_1, content_type='application/json')
        unlocked_db = self.client.post(f'/libraries/{self.library_id}/versions/', self.dummy_db_2, content_type='application/json')
        self.assertEqual(locked_db.status_code, status.HTTP_201_CREATED)
        self.assertEqual(unlocked_db.status_code, status.HTTP_201_CREATED)

        self.locked_db = json.loads(locked_db.content)
        self.unlocked_db = json.loads(unlocked_db.content)

        self.nuccore_id = self.locked_db['sequences'][0]['id']

    get_nuccores_id_200_response = {
        "id": "232e3b26-f347-4b09-88aa-2dcb52080d80",
        "annotations": [],
        "owner_database": {
            "id": "d83fc783-0a61-4b2f-94a8-d04d094e2a6c",
            "description": "Sequences gathered as of January 2022",
            "library": {
                "id": "ff5da812-1c0e-4682-9423-f08dbbc8bc5f",
                "custom_name": "Library of South American Fishes",
                "description": "A collection of voucher specimens from the literature.",
                "public": True,
                "owner": AuthTestCase.owner_serializer
            },
            "custom_name": "January Release",
            "version_number": "2.1.1"
        },
        "accession_number": "ON303390",
        "version": "ON303390.1",
        "definition": "Gymnotus cylindricus isolate 2094 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
        "organism": "Gymnotus cylindricus",
        "organelle": "mitochondrion",
        "isolate": "2094",
        "country": "Costa Rica",
        "specimen_voucher": "ROM:84772",
        "lat_lon": "10.54 N 83.50 W",
        "dna_sequence": "ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTCAACCCGGAGCCCTCCTCGGGGACGACCAAATTTATAATGTAATTGTTACTGCCCACGCTTTCGTAATAATTTTTTTTATAGTAATACCCATTATAATTGGAGGCTTCGGAAATTGATTAACTCCACTAATAATTGGAGCCCCAGACATAGCATTTCCCCGAATAAATAATATAAGCTTTTGACTCCTTCCTCCTTCTTTTTTACTTCTCCTTGCATCATCCGGAGTTGAAGCAGGGGCCGGAACAGGCTGAACAGTATACCCCCCTCTTGCAGGCAATCTTGCCCATGCAGGAGCCTCAGTAGATCTAACTATTTTCTCTCTTCATCTAGCCGGAGTTTCTTCAATTCTAGGATCCATTAACTTTATTACCACAATTATTAATATAAAACCTCCAGCCATTTCTCAATATCAAACCCCACTATTTATTTGATCACTTCTAGTAACCACTGTCCTTTTACTCCTCTCTCTTCCAGTACTAGCTGCTGGTATCACCATACTACTAACAGATCGAAACTTAAACACAACATTCTTTGACCCGGCGGGCGGAGGAGATCCTATTTTATATCAACATTTA",
        "translation": "",
        "type_material": "",
        "created": "2023-06-26T20:04:45.016196Z",
        "genbank_modification_date": "2022-07-04",
        "taxid": 699532,
        "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
        "title": "A new taxonomist-curated reference library of DNA barcodes for Neotropical electric fishes (Teleostei: Gymnotiformes)",
        "journal": "Zool J Linn Soc (2022) In press",
        "authors": "Janzen,F.H., Crampton,W.G. and Lovejoy,N.R.",
        "taxon_superkingdom": {
            "id": 2759,
            "scientific_name": "Eukaryota"
        },
        "taxon_kingdom": {
            "id": 33208,
            "scientific_name": "Metazoa"
        },
        "taxon_phylum": {
            "id": 7711,
            "scientific_name": "Chordata"
        },
        "taxon_class": {
            "id": 186623,
            "scientific_name": "Actinopteri"
        },
        "taxon_order": {
            "id": 8002,
            "scientific_name": "Gymnotiformes"
        },
        "taxon_family": {
            "id": 30771,
            "scientific_name": "Gymnotidae"
        },
        "taxon_genus": {
            "id": 36670,
            "scientific_name": "Gymnotus"
        },
        "taxon_species": {
            "id": 699532,
            "scientific_name": "Gymnotus cylindricus"
        }
    }

    def test_get_patch_delete_blastdbs_id_sequences(self):
        '''Test that a sequence can be retrieved, but cannot be patched or deleted from a locked database'''

        response = self.client.get(f'/nuccores/{self.nuccore_id}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        data = json.loads(response.content)
        self.assertEqualIgnoreId(data, self.get_nuccores_id_200_response)

        response = self.client.patch(f'/blastdbs/{self.locked_db["id"]}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)

        response = self.client.delete(f'/blastdbs/{self.locked_db["id"]}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)

    post_blastdbs_id_sequences_201_request = {
        "accession_numbers": ["ON303382.1"],
        "search_term": "ON303383.1[Accession]"
    }

    post_blastdbs_id_sequences_201_response = [
        {
            "id": "ba2e25aa-b781-4a0b-aadb-0ab3827ab302",
            "accession_number": "ON303383",
            "version": "ON303383.1",
            "organism": "Gymnotus cataniapo",
            "organelle": "mitochondrion",
            "isolate": "2063",
            "country": "Venezuela",
            "specimen_voucher": "UF:174332",
            "dna_sequence": "ATAGTATTTGGCGCCTGAGCCGGAATAGTTGGCACAGCCTTAAGCCTCCTCATCCGGGCAGAACTCAGTCAACCCGGGGCCCTCCTTGGCGACGACCAAATTTATAATGTAATTGTTACTGCCCACGCCTTCGTAATAATCTTCTTTATGGTGATACCCATCATGATTGGAGGCTTTGGAAACTGACTAATCCCACTAATAATCGGAGCCCCAGATATAGCATTCCCACGAATAAACAACATGAGCTTTTGACTTCTCCCGCCCTCTTTCCTGCTTCTCCTTGCCTCCTCAGGGGTTGAAGCTGGGGCTGGGACAGGCTGAACCGTATATCCCCCCCTTGCAGGCAACCTTGCCCACGCAGGAGCCTCAGTAGACCTGACTATCTTCTCCCTCCACCTTGCCGGGGTTTCTTCAATTCTTGGGTCTATTAACTTTATTACCACAATTATTAACATGAAACCCCCAGCTATTTCCCAATATCAAACCCCATTGTTTATTTGGGCCTTATTAGTAACCACCGTCCTCTTGCTTCTCTCCCTTCCAGTCTTAGCTGCTGGTATTACTATACTACTAACGGACCGAAATTTAAACACAACTTTCTTCGACCCCGCAGGCGGAGGGGACCCCATCCTATACCAACACCTA",
            "lat_lon": "5.56 N 67.47 W",
            "type_material": "",
            "created": "2023-06-26T21:42:53.548893Z",
            "updated": "2023-06-26T21:42:53.548907Z",
            "genbank_modification_date": "2022-07-04",
            "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
            "taxon_superkingdom": {
                "id": 2759,
                "scientific_name": "Eukaryota"
            },
            "taxon_kingdom": {
                "id": 33208,
                "scientific_name": "Metazoa"
            },
            "taxon_phylum": {
                "id": 7711,
                "scientific_name": "Chordata"
            },
            "taxon_class": {
                "id": 186623,
                "scientific_name": "Actinopteri"
            },
            "taxon_order": {
                "id": 8002,
                "scientific_name": "Gymnotiformes"
            },
            "taxon_family": {
                "id": 30771,
                "scientific_name": "Gymnotidae"
            },
            "taxon_genus": {
                "id": 36670,
                "scientific_name": "Gymnotus"
            },
            "taxon_species": {
                "id": 699527,
                "scientific_name": "Gymnotus cataniapo"
            },
            "annotations": []
        },
        {
            "id": "825d67b0-14ae-4f4a-b9b0-4c178b765f7d",
            "accession_number": "ON303382",
            "version": "ON303382.1",
            "organism": "Gymnotus carapo",
            "organelle": "mitochondrion",
            "isolate": "7005",
            "country": "Suriname",
            "specimen_voucher": "UF:180169",
            "dna_sequence": "ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTCAACCCGGAGCCCTCCTAGGCGATGATCAAATTTATAATGTAATTGTTACTGCCCACGCTTTCGTAATGATTTTCTTTATAGTAATGCCTATCATGATTGGAGGCTTCGGAAACTGATTAATTCCACTAATAATTGGAGCTCCAGACATAGCATTTCCCCGAATAAATAACATAAGCTTCTGACTTCTTCCCCCTTCTTTTTTACTTCTCCTTGCATCATCCGGAGTTGAAGCCGGAGCTGGAACAGGCTGAACAGTATACCCCCCTCTCGCAGGCAATCTTGCCCATGCAGGTGCTTCAGTAGACTTAACTATCTTCTCCCTTCACCTAGCCGGAGTTTCCTCTATTCTAGGGTCTATTAATTTTATTACCACAATTATTAATATGAAACCTCCAGCCATTTCTCAGTATCAAACCCCATTATTTATTTGAGCACTTCTAATTACTACTGTCCTTTTACTCCTCTCCCTTCCCGTACTAGCCGCTGGAATTACTATGTTATTAACAGATCGAAACTTAAACACAACCTTCTTTGACCCAGCAGGCGGTGGAGATCCCATTCTCTATCAACACTTA",
            "lat_lon": "5.25 N 55.10 W",
            "type_material": "",
            "created": "2023-06-26T21:42:53.548958Z",
            "updated": "2023-06-26T21:42:53.548963Z",
            "genbank_modification_date": "2022-07-04",
            "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Gymnotoidei,Gymnotidae,Gymnotus",
            "taxon_superkingdom": {
                "id": 2759,
                "scientific_name": "Eukaryota"
            },
            "taxon_kingdom": {
                "id": 33208,
                "scientific_name": "Metazoa"
            },
            "taxon_phylum": {
                "id": 7711,
                "scientific_name": "Chordata"
            },
            "taxon_class": {
                "id": 186623,
                "scientific_name": "Actinopteri"
            },
            "taxon_order": {
                "id": 8002,
                "scientific_name": "Gymnotiformes"
            },
            "taxon_family": {
                "id": 30771,
                "scientific_name": "Gymnotidae"
            },
            "taxon_genus": {
                "id": 36670,
                "scientific_name": "Gymnotus"
            },
            "taxon_species": {
                "id": 94172,
                "scientific_name": "Gymnotus carapo"
            },
            "annotations": []
        }
    ]

    def test_post_delete_blastdbs_id_sequences(self):
        '''Test POST and DELETE for a sequence to an unlocked database.'''
        response = self.client.post(f'/blastdbs/{self.unlocked_db["id"]}/sequences/', self.post_blastdbs_id_sequences_201_request, content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        data = json.loads(response.content)
        id_1 = data[0]["id"]
        self.assertEqualIgnoreId(data, self.post_blastdbs_id_sequences_201_response)

        response = self.client.delete(f'/nuccores/{id_1}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_204_NO_CONTENT)

        response = self.client.delete(f'/nuccores/{id_1}/', content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)

class RunTester(AuthTestCase):
    dummy_library_1 = {
        'custom_name': 'Library of South American Fishes',
        'description': 'A collection of voucher specimens from the literature.',
        'public': True
    }

    dummy_db_1 = {
        'custom_name': 'January Release',
        'description': 'Sequences gathered as of January 2022',
        'locked': True,
        'accession_numbers': ['ON303419.1', 'ON303462.1', 'ON303402.1', 'ON303378.1']
    }

    def setUp(self) -> None:
        library = self.client.post('/libraries/', self.dummy_library_1, content_type='application/json')
        self.assertEqual(library.status_code, status.HTTP_201_CREATED)
        self.library = json.loads(library.content)

        database = self.client.post(f'/libraries/{self.library["id"]}/versions/', self.dummy_db_1, content_type='application/json')
        self.assertEqual(database.status_code, status.HTTP_201_CREATED)
        self.database = json.loads(database.content)

    post_blastdbs_id_blast_blast_only_201_request = {
        'query_sequence': 'GGCACACTTTACATAGTGTTTGGCGCCTGGGCGGGTATAATTGGTACTGCTCTAAGCCTTCTAATCCGGGCCGAGCTTAACCAACCAGGCACCCTCCTAGAAGACGACCAAATTTATAATGTAGCCGTTACCGCCCATGCCTTCGTAATAATTTTCTTTATAGTTATGCCAATCATAATTGGAGGCTTTGGCAATTGGCTTATTCCCCTAATAATTGCCGCACCAGATATGGCATTCCCACGAATAAACAATATAAGCTTTTGGCTACTTCCTCCGTCATTCTTCCTCCTCCTAGCCTCTGCTGGCTTAGAGGCCGGGGTTGGAACAGGCTGGACCCTATACCCCCCTCTGGCCGGTAATGCAGCACACGCCGGAGCTTCCGTAGACTTAACCATTTTCTCCCTTCACTTGGCCGGTGTTTCATCTATCCTCGGCTCTATTAACTTTATCACTACAATTATTAATATAAAACCCC',
        'job_name': 'Unknown seq',
        'databaseSelect': '4cd539fe-5046-4283-800f-ffd06b48d469',
        'librarySelect': '14a53f56-486e-4005-913b-f73af4eb91a2',
        'create_hit_tree': False,
        'create_db_tree': False,
    }

    def test_post_blastdb_id_blast_blast_only(self):
        '''Test a run using only BLAST (no alignments)'''
        self.post_blastdbs_id_blast_blast_only_201_request['databaseSelect'] = self.database['id']
        self.post_blastdbs_id_blast_blast_only_201_request['librarySelect'] = self.library['id']
        raise NotImplementedError('Not yet implemented.')

    post_blastdbs_id_blast_alignment_201_request = {
        'query_sequence': 'GGCACACTTTACATAGTGTTTGGCGCCTGGGCGGGTATAATTGGTACTGCTCTAAGCCTTCTAATCCGGGCCGAGCTTAACCAACCAGGCACCCTCCTAGAAGACGACCAAATTTATAATGTAGCCGTTACCGCCCATGCCTTCGTAATAATTTTCTTTATAGTTATGCCAATCATAATTGGAGGCTTTGGCAATTGGCTTATTCCCCTAATAATTGCCGCACCAGATATGGCATTCCCACGAATAAACAATATAAGCTTTTGGCTACTTCCTCCGTCATTCTTCCTCCTCCTAGCCTCTGCTGGCTTAGAGGCCGGGGTTGGAACAGGCTGGACCCTATACCCCCCTCTGGCCGGTAATGCAGCACACGCCGGAGCTTCCGTAGACTTAACCATTTTCTCCCTTCACTTGGCCGGTGTTTCATCTATCCTCGGCTCTATTAACTTTATCACTACAATTATTAATATAAAACCCC',
        'job_name': 'Unknown seq',
        'databaseSelect': '4cd539fe-5046-4283-800f-ffd06b48d469',
        'librarySelect': '14a53f56-486e-4005-913b-f73af4eb91a2',
        'create_hit_tree': True,
        'create_db_tree': False,
    }

    post_blastdbs_id_blast_alignment_201_response = {
        "id": "a4baadd8-95fe-495f-8867-e2a456afcfcd",
        "job_name": "Unknown seq",
        "queries": [
            {
                "definition": "query_sequence",
                "query_sequence": "GGCACACTTTACATAGTGTTTGGCGCCTGGGCGGGTATAATTGGTACTGCTCTAAGCCTTCTAATCCGGGCCGAGCTTAACCAACCAGGCACCCTCCTAGAAGACGACCAAATTTATAATGTAGCCGTTACCGCCCATGCCTTCGTAATAATTTTCTTTATAGTTATGCCAATCATAATTGGAGGCTTTGGCAATTGGCTTATTCCCCTAATAATTGCCGCACCAGATATGGCATTCCCACGAATAAACAATATAAGCTTTTGGCTACTTCCTCCGTCATTCTTCCTCCTCCTAGCCTCTGCTGGCTTAGAGGCCGGGGTTGGAACAGGCTGGACCCTATACCCCCCTCTGGCCGGTAATGCAGCACACGCCGGAGCTTCCGTAGACTTAACCATTTTCTCCCTTCACTTGGCCGGTGTTTCATCTATCCTCGGCTCTATTAACTTTATCACTACAATTATTAATATAAAACCCC",
                "hits": [],
                "results_species_name": None,
                "accuracy_category": None,
                "original_species_name": "",
                "write_tree_identifier": "query_sequence||query"
            }
        ],
        "db_used": {
            "id": "4cd539fe-5046-4283-800f-ffd06b48d469",
            "description": "Sequences gathered as of January 2022",
            "library": {
                "id": "14a53f56-486e-4005-913b-f73af4eb91a2",
                "custom_name": "Library of South American Fishes",
                "description": "A collection of voucher specimens from the literature.",
                "public": True,
                "owner": AuthTestCase.owner_serializer
            },
            "custom_name": "January Release",
            "version_number": "1.1.1"
        },
        "start_time": None,
        "status": "QUE",
        "received_time": "2023-06-22T02:10:04.255023Z",
        "end_time": None,
        "error_time": None,
        "create_hit_tree": True,
        "hit_tree": "",
        "alignment_job_id": "",
        "create_db_tree": False,
        "db_tree": "",
        "complete_alignment_job_id": "",
    }

    get_runs_id_alignment_200_response = {
        "id": "cdc7ae5d-5756-4d49-8a80-028221a49846",
        "job_name": "Unknown seq",
        "queries": [
            {
            "definition": "query_sequence",
            "query_sequence": "GGCACACTTTACATAGTGTTTGGCGCCTGGGCGGGTATAATTGGTACTGCTCTAAGCCTTCTAATCCGGGCCGAGCTTAACCAACCAGGCACCCTCCTAGAAGACGACCAAATTTATAATGTAGCCGTTACCGCCCATGCCTTCGTAATAATTTTCTTTATAGTTATGCCAATCATAATTGGAGGCTTTGGCAATTGGCTTATTCCCCTAATAATTGCCGCACCAGATATGGCATTCCCACGAATAAACAATATAAGCTTTTGGCTACTTCCTCCGTCATTCTTCCTCCTCCTAGCCTCTGCTGGCTTAGAGGCCGGGGTTGGAACAGGCTGGACCCTATACCCCCCTCTGGCCGGTAATGCAGCACACGCCGGAGCTTCCGTAGACTTAACCATTTTCTCCCTTCACTTGGCCGGTGTTTCATCTATCCTCGGCTCTATTAACTTTATCACTACAATTATTAATATAAAACCCC",
            "hits": [
                {
                "db_entry": {
                    "id": "0a0ea541-a32f-4cb0-9daf-972e5c26db00",
                    "accession_number": "ON303419",
                    "version": "ON303419.1",
                    "definition": "Microsternarchus bilineatus isolate 2138 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
                    "organism": "Microsternarchus bilineatus",
                    "country": "Brazil",
                    "specimen_voucher": "MCP 45480",
                    "type_material": "",
                    "lat_lon": "3.42 S 64.70 W",
                    "annotations": []
                },
                "query_accession_version": "query_sequence",
                "subject_accession_version": "ON303419.1",
                "percent_identity": "78.233",
                "alignment_length": 464,
                "mismatches": 99,
                "gap_opens": 2,
                "query_start": 13,
                "query_end": 475,
                "sequence_start": 1,
                "sequence_end": 463,
                "evalue": "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000061200000000000000",
                "bit_score": "296.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
                },
                {
                "db_entry": {
                    "id": "64451ee2-69bd-4839-8904-2caa59fab839",
                    "accession_number": "ON303462",
                    "version": "ON303462.1",
                    "definition": "Sternopygus branco isolate 2108 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
                    "organism": "Sternopygus branco",
                    "country": "Brazil",
                    "specimen_voucher": "MCP 32246",
                    "type_material": "paratype of Sternopygus branco",
                    "lat_lon": "3.12 S 64.78 W",
                    "annotations": []
                },
                "query_accession_version": "query_sequence",
                "subject_accession_version": "ON303462.1",
                "percent_identity": "79.570",
                "alignment_length": 465,
                "mismatches": 91,
                "gap_opens": 4,
                "query_start": 13,
                "query_end": 475,
                "sequence_start": 1,
                "sequence_end": 463,
                "evalue": "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000006030000",
                "bit_score": "329.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
                },
                {
                "db_entry": {
                    "id": "16ad29dd-3766-414b-97ea-77cb9b8e42e7",
                    "accession_number": "ON303402",
                    "version": "ON303402.1",
                    "definition": "Gymnotus pedanopterus isolate 2059 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
                    "organism": "Gymnotus pedanopterus",
                    "country": "Venezuela",
                    "specimen_voucher": "UF:174328",
                    "type_material": "",
                    "lat_lon": "3.93 N 67.61 W",
                    "annotations": []
                },
                "query_accession_version": "query_sequence",
                "subject_accession_version": "ON303402.1",
                "percent_identity": "80.130",
                "alignment_length": 463,
                "mismatches": 92,
                "gap_opens": 0,
                "query_start": 13,
                "query_end": 475,
                "sequence_start": 1,
                "sequence_end": 463,
                "evalue": "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000060",
                "bit_score": "346.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
                },
                {
                "db_entry": {
                    "id": "b60e83b7-9920-4047-9454-006947a4abc5",
                    "accession_number": "ON303378",
                    "version": "ON303378.1",
                    "definition": "Gymnotus anguillaris isolate 10545 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
                    "organism": "Gymnotus anguillaris",
                    "country": "Suriname",
                    "specimen_voucher": "ROM:100941",
                    "type_material": "",
                    "lat_lon": "5.53 N 54.42 W",
                    "annotations": []
                },
                "query_accession_version": "query_sequence",
                "subject_accession_version": "ON303378.1",
                "percent_identity": "80.562",
                "alignment_length": 463,
                "mismatches": 90,
                "gap_opens": 0,
                "query_start": 13,
                "query_end": 475,
                "sequence_start": 1,
                "sequence_end": 463,
                "evalue": "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
                "bit_score": "357.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
                }
            ],
            "results_species_name": "Gymnotus anguillaris",
            "accuracy_category": "Unknown ID",
            "original_species_name": "",
            "write_tree_identifier": "query_sequence||query"
            }
        ],
        "db_used": {
            "id": "4cd539fe-5046-4283-800f-ffd06b48d469",
            "description": "Sequences gathered as of January 2022",
            "library": {
            "id": "14a53f56-486e-4005-913b-f73af4eb91a2",
            "custom_name": "Library of South American Fishes",
            "description": "A collection of voucher specimens from the literature.",
            "public": True,
            "owner": AuthTestCase.owner_serializer
            },
            "custom_name": "January Release",
            "version_number": "1.1.1"
        },
        "start_time": "2023-06-27T00:46:29.015260Z",
        "status": "FIN",
        "received_time": "2023-06-27T00:46:28.064601Z",
        "end_time": "2023-06-27T00:47:30.660222Z",
        "error_time": None,
        "create_hit_tree": True,
        "hit_tree": "(\n(\nquery_sequence||query:0.11672,\nON303462.1|Sternopygus_branco:0.09279)\n:0.00295,\nON303419.1|Microsternarchus_bilineatus:0.09935,\n(\nON303402.1|Gymnotus_pedanopterus:0.05705,\nON303378.1|Gymnotus_anguillaris:0.05303)\n:0.02313);\n",
        "alignment_job_id": "clustalo-R20230627-014648-0996-2185536-p1m",
        "create_db_tree": False,
        "db_tree": "",
        "complete_alignment_job_id": ""
    }

    get_runs_id_alignment_status_200_response = {
        "id": "cdc7ae5d-5756-4d49-8a80-028221a49846",
        "job_name": "Unknown seq",
        "received_time": "2023-06-24T16:12:11.328957Z",
        "status": "FIN",
        "start_time": "2023-06-24T16:13:12.892249Z",
        "end_time": "2023-06-24T16:13:14.777021Z",
        "error_time": None
    }


    def test_post_blastdb_id_blast_db_alignment(self):
        '''Test BLAST and hit tree alignment with 4 sequences and 1 query.'''
        # Specify library and database
        self.post_blastdbs_id_blast_alignment_201_request['databaseSelect'] = self.database['id']
        self.post_blastdbs_id_blast_alignment_201_request['librarySelect'] = self.library['id']

        raise NotImplementedError('Request and response data not complete for alignment')

        response = self.client.post(f'/blastdbs/{self.database["id"]}/blast', self.post_blastdbs_id_blast_alignment_201_request, content_type='application/json')
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        data = json.loads(response.content)
        self.assertEqualIgnoreId(data, self.post_blastdbs_id_blast_alignment_201_response)

        self.assertEqual(data['db_used']['id'], self.database['id'])
        self.assertEqual(data['db_used']['library']['id'], self.library['id'])
        run_id = data['id']
        job_name = self.post_blastdbs_id_blast_alignment_201_request["job_name"]

        max_tries = 30
        poll_interval = 1
        tries = 1

        # Validate status 
        bad_statuses = [BlastRun.JobStatus.ERRORED, BlastRun.JobStatus.DENIED, BlastRun.JobStatus.UNKNOWN]
        
        while tries <= max_tries:
            response = self.client.get(f'runs/{run_id}/status', content_type='application/json')
            self.assertEqual(response.status_code, status.HTTP_200_OK)
            status_data = json.loads(response.content)
            self.assertEqual(status_data["job_name"], job_name)
            self.assertEqual(status_data["id"], run_id)
            self.assertNotIn(status_data['status'], bad_statuses)
            if status_data['status'] == BlastRun.JobStatus.FINISHED:
                break
            else:
                self.assertIsNone(status_data['end_time'])
                self.assertIsNone(status_data['error_time'])
            tries = tries + 1
            time.sleep(poll_interval)

        # Validate results
        response = self.client.get(f'runs/{run_id}/')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        results = json.loads(response.content)
        # ignore differences in clustal job ids
        results['alignment_job_id'] = ''
        results['complete_alignment_job_id'] = ''
        self.get_runs_id_alignment_200_response['alignment_job_id'] = ''
        self.get_runs_id_alignment_200_response['complete_alignment_job_id'] = ''
        self.assertEqualIgnoreId(results, self.get_runs_id_alignment_200_response) 
