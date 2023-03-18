import json
import os
from rest_framework import serializers
from barcode_blastn.models import BlastQuerySequence, BlastRun, Hit, NuccoreSequence, BlastDb

blast_db_title = 'Blast Database'
nuccore_title = 'GenBank Accession'
hit_title = 'BLASTN hit'
query_title = 'Query Sequence'
run_title = 'Run'

class BlastDbShortSerializer(serializers.ModelSerializer):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = False

    f"""
    Information about a {blast_db_title}, when displayed under a {nuccore_title} it includes.
    """
    class Meta:
        model = BlastDb
        ref_name = blast_db_title + ' summary'
        fields = ['id', 'custom_name', 'description']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "Neotropical electric knifefish",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database."
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
        fields = ['id', 'owner_database', 'accession_number', 'definition', 'organism', 'organelle', 'isolate', 'country', 'specimen_voucher', 'lat_lon', 'dna_sequence', 'translation', 'type_material', 'created']
        example = {
            "id": "5100cbd8-2fda-4b42-8aa1-10ede078448b",
            "owner_database": BlastDbShortSerializer.Meta.example,
            "accession_number": "ON303341",
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
        fields = ['accession_number', 'organism', 'country'] 
        example = {
            "accession_number": "GU701771",
            "organism": "Gymnotus pantherinus",
            "country": "Brazil: Sao Paulo, Upper Parana Basin"
        }

class BlastDbCreateSerializer(serializers.ModelSerializer):
    f"""Information required when creating a {blast_db_title}"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in [ "custom_name", "description" ]:
            self.fields[name].required = True
        self.fields['id'].required = True
        self.fields['id'].read_only = True

    '''
    Create a new blastdb
    '''
    class Meta:
        model = BlastDb
        ref_name = blast_db_title + ' creation'
        fields = ['id', 'custom_name', 'description']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "Newly Sequenced Species",
            "description": "A collection of new sequences from several species of interest."
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
        fields = ['id', 'accession_number', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher', 'lat_lon',  'type_material', 'dna_sequence']
        example = {
            "id": "5100cbd8-2fda-4b42-8aa1-10ede078448b",
            "accession_number": "ON303341",
            "definition": "Brachyhypopomus arrayae isolate 12586 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
            "organism": "Brachyhypopomus arrayae",
            "isolate": "12586",
            "country": "Bolivia",
            "specimen_voucher": "ANSP:197574",
            "lat_lon": "11.03 S 66.09 W",
            "type_material": "paratype of Brachyhypopomus arrayae",
            "dna_sequence": "ATAGTATTTGGTGCATGAGCTGGGATAGTAGGCACAGCCTTAAGCCTCTTAATCCGAGCAGAACTAAGCCAGCCAGGAGCTCTTATGGGCGACGACCAAATTTACAATGTGATTGTTACTGCGCACGCTTTCGTAATAATTTTCTTCATGGTTATGCCCATTATAATCGGCGGGTTCGGCAACTGATTAATTCCCCTAATACTCGGTGCCCCTGACATGGCATTCCCACGAATAAACAACATAAGCTTCTGACTTCTGCCCCCATCATTCCTTCTACTCCTTGCATCCTCTGGGGTCGAAGCGGGAGCCGGAACCGGCTGAACTGTTTACCCCCCTCTCGCTAGCAACCTCGCCCACGCAGGGGCCTCCGTTGATCTAACTATCTTCTCCCTTCACCTTGCTGGGGTTTCTTCCATCCTTGGCTCTATCAACTTCATTACTACCATTATTAACATGAAACCCCCAGCCATATCTCAGTATCAAACCCCTCTATTTATTTGAGCGCTCCTAATTACCACAGTTCTCCTACTGTTATCCCTTCCCGTACTGGCCGCTGGTATCACCATGCTGCTAACAGACCGAAACCTAAATACAACCTTCTTCGACCCCGCAGGAGGAGGGGACCCCGTCCTTTATCAGCACTTA"
        }

class BlastDbSerializer(serializers.ModelSerializer):
    f'''
    Show detailed information about a specific {blast_db_title}
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        read_only_fields = ['id', 'sequences']
        for name in read_only_fields:
            self.fields[name].read_only = True
        for name in ['id', 'custom_name', 'description', 'locked', 'sequences']:
            self.fields[name].required = True

    sequences = BlastDbSequenceEntrySerializer(many=True, read_only=True)
    class Meta:
        model = BlastDb
        ref_name = blast_db_title
        fields = ['id', 'custom_name', 'description', 'locked', 'sequences']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "Neotropical electric knifefish",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "locked": True,
            "sequences": [
                {
                    "accession_number": "GU701771",
                    "organism": "Gymnotus pantherinus",
                    "country": "Brazil: Sao Paulo, Upper Parana Basin"
                },
                {
                    "accession_number": "KF533332",
                    "organism": "Hypopomus artedi",
                    "country": "Guyana: Mazaruni"
                }
            ]
        }

class BlastDbListSerializer(serializers.ModelSerializer):
    '''
    Show a condensed summary of each blastdb, to list out all blastdbs
    '''
    
    sequences = BlastDbSequenceEntryShortSerializer(many=True, read_only=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = False
            self.fields[name].read_only = True

    class Meta:
        model = BlastDb
        ref_name = blast_db_title 
        fields = ['id', 'custom_name', 'description', 'locked', 'sequences']
        example = [ BlastDbSerializer.Meta.example ]

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
        fields = ['accession_number', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher', 'type_material', 'lat_lon']
        example = {
            "db_entry": {
                "accession_number": "ON303423",
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
    query_file = serializers.FileField()

    class Meta:
        model = BlastRun
        ref_name = run_title + ' submission'
        fields = ['id', 'job_name', 'query_sequence', 'query_file', 'create_hit_tree', 'create_db_tree']

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


