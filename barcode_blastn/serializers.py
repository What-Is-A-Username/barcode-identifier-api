import json
import os
from typing import Union
from django.contrib.auth.models import User
from rest_framework import serializers
from barcode_blastn.models import BlastQuerySequence, BlastRun, DatabaseShare, Hit, Library, NuccoreSequence, BlastDb

library_title = 'Reference Library'
blast_db_title = 'BLAST Database'
nuccore_title = 'GenBank Accession'
hit_title = 'BLASTN hit'
query_title = 'Query Sequence'
run_title = 'Run'
share_title = 'Share Details'

class LibraryOwnerSerializer(serializers.ModelSerializer):
    f'''
    Show basic information about the owner of a database
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in ['username', 'email']:
            self.fields[name].read_only = True

    class Meta:
        model = User
        ref_name = f'{library_title} owner'
        fields = ['username', 'email']
        example = {
            'username': 'John Doe',
            'email': 'johndoe@example.com'
        }

class LibraryShortSerializer(serializers.ModelSerializer):
    f'''
    Used when returning a list of {library_title}
    '''
    owner = LibraryOwnerSerializer(many=False, read_only=True)

    def __init__(self, instance=None, data=..., **kwargs):
        super().__init__(instance, data, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = False
            self.fields[name].read_only = False
        for f in ['id', 'custom_name', 'description', 'public', 'owner']:
            self.fields[f].read_only = True 
            self.fields[f].required = False

    class Meta:
        model = Library
        ref_name = blast_db_title + ' summary'
        fields = ['id', 'custom_name', 'description', 'public', 'owner']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "Neotropical electric knifefish",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "public": True,
            "owner": LibraryOwnerSerializer.Meta.example
        }
        tags = ['DBS']

class BlastDbTinySerializer(serializers.ModelSerializer):
    f"""
    Information about a {blast_db_title}, when displayed as a list under its parent reference library."""

    class Meta:
        model = BlastDb
        fields = ['id', 'description', 'version_number', 'sequence_count', 'created', 'locked']

class BlastDbShortSerializer(serializers.ModelSerializer):
    f"""
    Information about a {blast_db_title}, when displayed under a {nuccore_title} it includes.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = False

    library = LibraryShortSerializer(many=False, read_only=True)

    class Meta:
        model = BlastDb
        ref_name = blast_db_title + ' summary'
        fields = ['id', 'description', 'library']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "library": LibraryShortSerializer.Meta.example
        }
        tags = ['DBS']

class NuccoreSequenceSerializer(serializers.ModelSerializer):
    f"""
    Complete information about a {nuccore_title} found within a {blast_db_title}. 
    """
    owner_database = BlastDbShortSerializer(many=False, read_only=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = False

    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title + ' summary'
        fields = ['id', 'owner_database', 'accession_number', 'version', 'definition', 'organism', 'organelle', 'isolate', 'country', 'specimen_voucher', 'lat_lon', 'dna_sequence', 'translation', 'type_material', 'created']
        example = {
            "id": "5100cbd8-2fda-4b42-8aa1-10ede078448b",
            "owner_database": BlastDbShortSerializer.Meta.example,
            "accession_number": "ON303341",
            "version": "ON303341.1",
            "definition": "Brachyhypopomus arrayae isolate 12586 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
            "organism": "Brachyhypopomus arrayae",
            "organelle": "mitochondrion",
            "isolate": "12586", 
            "country": "Bolivia",
            "specimen_voucher": "ANSP:197574",
            "lat_lon": "11.03 S 66.09 W",
            "dna_sequence": "ATAGTATTTGGTGCATGAGCTGGGATAGTAGGCACAGCCTTAAGCCTCTTAATCCGAGCAGAACTAAGCCAGCCAGGAGCTCTTATGGGCGACGACCAAATTTACAATGTGATTGTTACTGCGCACGCTTTCGTAATAATTTTCTTCATGGTTATGCCCATTATAATCGGCGGGTTCGGCAACTGATTAATTCCCCTAATACTCGGTGCCCCTGACATGGCATTCCCACGAATAAACAACATAAGCTTCTGACTTCTGCCCCCATCATTCCTTCTACTCCTTGCATCCTCTGGGGTCGAAGCGGGAGCCGGAACCGGCTGAACTGTTTACCCCCCTCTCGCTAGCAACCTCGCCCACGCAGGGGCCTCCGTTGATCTAACTATCTTCTCCCTTCACCTTGCTGGGGTTTCTTCCATCCTTGGCTCTATCAACTTCATTACTACCATTATTAACATGAAACCCCCAGCCATATCTCAGTATCAAACCCCTCTATTTATTTGAGCGCTCCTAATTACCACAGTTCTCCTACTGTTATCCCTTCCCGTACTGGCCGCTGGTATCACCATGCTGCTAACAGACCGAAACCTAAATACAACCTTCTTCGACCCCGCAGGAGGAGGGGACCCCGTCCTTTATCAGCACTTA",
            "translation": "",
            "type_material": "paratype of Brachyhypopomus arrayae",
            "created": "2023-02-20T01:00:31.240961Z"
        }

class NuccoreSequenceAddSerializer(serializers.ModelSerializer):
    '''
    Required fields for adding a sequence accession to the database
    '''
    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title + ' addition'
        fields = ['accession_number']
        example = {'accession_number': 'ON303341'}

class NuccoreSequenceBulkAddSerializer(serializers.Serializer):
    '''
    Serialize data when receiving a request to bulk add sequences to a database
    '''
    accession_numbers = serializers.ListField(
        child=serializers.CharField()
    )

    class Meta:
        ref_name = nuccore_title
        example = {
            "accession_numbers": [
                "GU701771",
                "KF533332"
            ]
        }

class BlastDbSequenceEntryShortSerializer(serializers.ModelSerializer):
    f'''Information about a {nuccore_title} when included in a list of {blast_db_title} contents.'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = True 
            self.fields[name].read_only = True 

    '''
    Show a condensed summary of a sequence, in order to display with a list of all blastdbs
    '''
    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title + ' item'
        fields = ['accession_number', 'version', 'organism', 'specimen_voucher', 'country'] 
        example = {
            "accession_number": "GU701771",
            "version": "ON303341.1",
            "organism": "Gymnotus pantherinus",
            'specimen_voucher': 'LBP-24532',
            "country": "Brazil: Sao Paulo, Upper Parana Basin"
        }

class LibraryCreateSerializer(serializers.ModelSerializer):
    f'''
    Information required for creating a {library_title}'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in [ "custom_name", "description", "id"]:
            self.fields[name].required = True
        self.fields['id'].read_only = True
        self.fields['public'].required = False

    '''
    Create a new blastdb
    '''
    class Meta:
        model = Library
        ref_name = library_title + ' creation'
        fields = ['id', 'custom_name', 'description', 'public']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "Newly Sequenced Species Reference Library",
            "description": "A collection of new sequences from several species of interest.",
            "public": True
        }

class BlastDbSequenceEntrySerializer(serializers.ModelSerializer):
    f'''
    Show a summary of a {nuccore_title}, to be shown when its corresponding {blast_db_title} or {hit_title} is shown.
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = True
            self.fields[name].read_only = True

    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title
        fields = ['id', 'accession_number', 'version', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher', 'lat_lon',  'type_material', 'dna_sequence']
        example = {
            "id": "5100cbd8-2fda-4b42-8aa1-10ede078448b",
            "accession_number": "ON303341",
            "version": "ON303341.1",
            "definition": "Brachyhypopomus arrayae isolate 12586 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
            "organism": "Brachyhypopomus arrayae",
            "isolate": "12586",
            "country": "Bolivia",
            "specimen_voucher": "ANSP:197574",
            "lat_lon": "11.03 S 66.09 W",
            "type_material": "paratype of Brachyhypopomus arrayae",
            "dna_sequence": "ATAGTATTTGGTGCATGAGCTGGGATAGTAGGCACAGCCTTAAGCCTCTTAATCCGAGCAGAACTAAGCCAGCCAGGAGCTCTTATGGGCGACGACCAAATTTACAATGTGATTGTTACTGCGCACGCTTTCGTAATAATTTTCTTCATGGTTATGCCCATTATAATCGGCGGGTTCGGCAACTGATTAATTCCCCTAATACTCGGTGCCCCTGACATGGCATTCCCACGAATAAACAACATAAGCTTCTGACTTCTGCCCCCATCATTCCTTCTACTCCTTGCATCCTCTGGGGTCGAAGCGGGAGCCGGAACCGGCTGAACTGTTTACCCCCCTCTCGCTAGCAACCTCGCCCACGCAGGGGCCTCCGTTGATCTAACTATCTTCTCCCTTCACCTTGCTGGGGTTTCTTCCATCCTTGGCTCTATCAACTTCATTACTACCATTATTAACATGAAACCCCCAGCCATATCTCAGTATCAAACCCCTCTATTTATTTGAGCGCTCCTAATTACCACAGTTCTCCTACTGTTATCCCTTCCCGTACTGGCCGCTGGTATCACCATGCTGCTAACAGACCGAAACCTAAATACAACCTTCTTCGACCCCGCAGGAGGAGGGGACCCCGTCCTTTATCAGCACTTA"
        }

class BlastDbCreateSerializer(serializers.ModelSerializer):
    f"""Information required when creating a {blast_db_title}"""

    accession_numbers = serializers.ListField(child=serializers.CharField(), required=False)
    base = serializers.UUIDField()
    library = LibraryShortSerializer(many=False, read_only=True, required=False)
    sequences = BlastDbSequenceEntrySerializer(many=True, read_only=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in ['id', 'library', 'version_number', 'sequences']:
            self.fields[name].read_only = True
        for name in ['id', 'version_number', 'description', 'locked', 'sequences']:
            self.fields[name].required = True
        for name in ['accession_numbers', 'base']:
            self.fields[name].required = False

    '''
    Create a new blastdb
    '''
    class Meta:
        model = BlastDb
        ref_name = blast_db_title + ' creation'
        fields = ['id', 'library', 'version_number', 'description', 'locked', 'sequences', 'base', 'accession_numbers']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "library": LibraryShortSerializer.Meta.example,
            "description": "A collection of new sequences from several species of interest.",
            "version_number": '254.2.1',
        }

class BlastDbEditSerializer(serializers.ModelSerializer):
    f'''
    Determine the editable fields for a specific {blast_db_title}.
    These fields are editable on the admin page and through PATCH requests.
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in ['description', 'locked']:
            self.fields[name].required = True

    class Meta:
        model = BlastDb
        ref_name = 'Editable' + blast_db_title + ' Fields'
        fields = ['description', 'locked']
        example = {
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "locked": True,
        }

class BlastDbSerializer(serializers.ModelSerializer):
    f'''
    Show detailed information about a specific {blast_db_title}
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        read_only_fields = ['id', 'library', 'sequences', 'version_number']
        for name in read_only_fields:
            self.fields[name].read_only = True
        for name in ['id', 'library', 'description', 'locked', 'sequences', 'version_number']:
            self.fields[name].required = True

    library = LibraryShortSerializer(many=False, read_only=True)
    sequences = BlastDbSequenceEntrySerializer(many=True, read_only=True)

    class Meta:
        model = BlastDb
        ref_name = blast_db_title
        fields = ['id', 'library', 'version_number', 'description', 'locked', 'sequences']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "owner": LibraryOwnerSerializer.Meta.example,
            "locked": True,
            "sequences": BlastDbSequenceEntrySerializer.Meta.example
        }

class LibrarySerializer(serializers.ModelSerializer):
    f'''
    Show detailed information about a specific {library_title}
    '''
    owner = LibraryOwnerSerializer(many=False, read_only=True)
    blastdb_set = BlastDbTinySerializer(many=True, read_only=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        read_only_fields = ['id', 'owner', 'blastdb_set']
        for name in read_only_fields:
            self.fields[name].read_only = True
        for name in ['id', 'custom_name', 'description', 'owner', 'public', 'blastdb_set']:
            self.fields[name].required = True

    class Meta:
        model = Library
        ref_name = library_title
        fields = ['id', 'custom_name', 'description', 'owner', 'public', 'blastdb_set']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "Neotropical electric knifefish",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "owner": LibraryOwnerSerializer.Meta.example,
            "blastdb_set": BlastDbSerializer.Meta.example,
            "public": True,
        }

class LibraryEditSerializer(serializers.ModelSerializer):
    f'''
    Determine the editable fields for a specific {library_title}.
    These fields are editable on the admin page and through PATCH requests.
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in ['custom_name', 'description', 'public']:
            self.fields[name].required = True

    class Meta:
        model = Library
        ref_name = library_title + ' edits'
        fields = ['custom_name', 'description', 'public']
        example = {
            "custom_name": "Neotropical electric knifefish",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "public": True,
        }

class LibraryListSerializer(serializers.ModelSerializer):
    owner = LibraryOwnerSerializer(many=False, read_only=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = False
            self.fields[name].read_only = True

    class Meta:
        model = Library
        ref_name = library_title 
        fields = ['id', 'custom_name', 'description', 'owner', 'public']
        example = [ LibrarySerializer.Meta.example ]

class BlastDbListSerializer(serializers.ModelSerializer):
    '''
    Show a condensed summary of each blastdb, to list out all blastdbs
    '''
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = False
            self.fields[name].read_only = True

    class Meta:
        model = BlastDb
        ref_name = blast_db_title 
        fields = ['id', 'version_number', 'description', 'locked']
        example = [{
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "version_number": '254.2.1',
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "locked": True,
        }]

class HitSerializer(serializers.ModelSerializer):
    '''
    Serialize a hit to be returned with a list of all hits
    '''
    db_entry = BlastDbSequenceEntrySerializer(many=False, read_only=True)

    class Meta:
        model = Hit
        ref_name = hit_title + ' item'
        fields = ['db_entry', 'owner_run', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']
        example = {
            "db_entry": BlastDbSequenceEntrySerializer.Meta.example,
            "query_accession_version": "MG653404.1",
            "subject_accession_version": "MG653404",
            "percent_identity": "100.000",
            "alignment_length": 490,
            "mismatches": 0,
            "gap_opens": 0,
            "query_start": 1,
            "query_end": 490,
            "sequence_start": 1,
            "sequence_end": 490,
            "evalue": "0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
            "bit_score": "905.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        }
        
class NuccoreSequenceHitSerializer(serializers.ModelSerializer):
    '''
    Show a summary of a sequence, when a query registers a hit on it.
    '''
    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title + ' summary'
        fields = ['accession_number', 'version', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher', 'type_material', 'lat_lon']
        example = {
            "db_entry": {
                "accession_number": "ON303423",
                "version": "ON303423.1",
                "definition": "Porotergus duende isolate 2916 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
                "organism": "Porotergus duende",
                "isolate": "2916",
                "country": "Brazil",
                "specimen_voucher": "MCP 37359",
                "type_material": "paratype of Porotergus duende",
                "lat_lon": "3.22 S 54.38 W"
            }
        }

class HitEntrySerializer(serializers.ModelSerializer):
    '''
    Serialize a hit to be returned with a list of all hits
    '''
    db_entry = NuccoreSequenceHitSerializer(many=False, read_only=True)

    class Meta:
        model = Hit
        ref_name = hit_title + ' item'
        fields = ['db_entry', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']
        example = {
            "db_entry": NuccoreSequenceHitSerializer.Meta.example,
            "query_accession_version": "MG653404.1",
            "subject_accession_version": "ON303423",
            "percent_identity": "88.247",
            "alignment_length": 485,
            "mismatches": 57,
            "gap_opens": 0,
            "query_start": 2,
            "query_end": 486,
            "sequence_start": 70,
            "sequence_end": 554,
            "evalue": "0.0000000000000000000000000000000000000000000000000000000000054500000000000000000000000000000000000000",
            "bit_score": "580.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
        }

class BlastQuerySequenceSerializer(serializers.ModelSerializer):
    '''
    Show a sequence submitted in a run
    '''
    class Meta:
        model = BlastQuerySequence
        ref_name = query_title
        fields = ['definition', 'query_sequence']
        example = {
            "definition": "Steatogenys_elegans isolate 8807 cytochrome c oxidase subunit I (COI) gene, partial cds; mitochondrial",
            "query_sequence": "GGCACCCTTTATATAGTGTTTGGTGCCTGAGCCGGAATGGTTGGCACGGCCTTAAGCCTCCTTATTCGAGCCGAGCTAAGCCAACCCGGGGCCCTAATGGGTGATGACCAGATTTACAATGTTA"
        }

class BlastRunRunSerializer(serializers.ModelSerializer):
    '''
    Required fields for submitting a blast run
    '''
    query_sequence = serializers.CharField()
    # query_file = serializers.FileField()

    class Meta:
        model = BlastRun
        ref_name = run_title + ' submission'
        fields = ['id', 'job_name', 'query_sequence', 'create_hit_tree', 'create_db_tree']

def load_run_example():
    '''
    Load example of run result response from .json file
    '''
    f = open(os.path.abspath('./barcode_blastn/run_example.json'))
    data = json.load(f)
    f.close()
    return data

class BlastRunSerializer(serializers.ModelSerializer):
    '''
    Show results of a blast run
    '''
    db_used = BlastDbShortSerializer(many=False, read_only=True)
    hits = HitEntrySerializer(many=True, read_only=True)
    queries = BlastQuerySequenceSerializer(many=True, read_only=True)    

    class Meta:
        model = BlastRun    
        ref_name = run_title
        fields = ['id', 'job_name', 'queries', 'db_used', 'runtime', 'job_status', 'job_start_time', 'job_end_time', 'job_error_time', 'hits', 'create_hit_tree', 'hit_tree', 'alignment_job_id', 'create_db_tree', 'db_tree', 'complete_alignment_job_id']
        example = load_run_example()

class BlastRunStatusSerializer(serializers.ModelSerializer):    
    '''
    Show status of a blast run
    '''
    class Meta:
        model = BlastRun    
        ref_name = run_title + ' status'
        fields = ['id', 'job_name', 'runtime', 'job_status', 'job_start_time', 'job_end_time', 'job_error_time']
        example = {
            "id": "2e5898a3-14da-4e7f-9599-ba1ef35f1e7a",
            "job_name": "two sequences",
            "runtime": "2023-02-21T03:27:11.323521Z",
            "job_status": "FIN",
            "job_start_time": "2023-02-21T03:27:12.427193Z",
            "job_end_time": "2023-02-21T03:28:02.102346Z",
            "job_error_time": None,
        }

class UserSerializer(serializers.ModelSerializer):
    class Meta:
        model = User 
        fields = ['id', 'username', 'email', 'is_staff', 'is_superuser']
        example = {
            "id": 2316,
            "username": "admin",
            "username": "JohnSmith",
            "is_staff": False,
            "is_superuser": False
        }

class DatabaseShareSerializer(serializers.ModelSerializer):
    '''
    Show the action of an admin sharing a database with another user.
    '''

    database = BlastDbShortSerializer(many=False, read_only=True)
    grantee = UserSerializer(many=False, read_only=True)

    class Meta:
        model = DatabaseShare
        ref_name = share_title 
        fields = ['database', 'grantee']

