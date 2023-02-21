from rest_framework import serializers
from barcode_blastn.models import BlastQuerySequence, BlastRun, Hit, NuccoreSequence, BlastDb

class BlastDbShortSerializer(serializers.ModelSerializer):
    class Meta:
        model = BlastDb
        fields = ['id', 'custom_name', 'description']

class NuccoreSequenceSerializer(serializers.ModelSerializer):
    owner_database = BlastDbShortSerializer(many=False, read_only=True)

    class Meta:
        model = NuccoreSequence
        fields = ['id', 'owner_database', 'accession_number', 'definition', 'organism', 'organelle', 'mol_type', 'isolate', 'country', 'specimen_voucher', 'lat_lon', 'dna_sequence', 'translation', 'type_material', 'created']

'''
Show a condensed summary of a sequence, in order to display with a list of all blastdbs
'''
class BlastDbSequenceEntryShortSerializer(serializers.ModelSerializer):
    class Meta:
        model = NuccoreSequence
        fields = ['accession_number', 'organism', 'country'] 

'''
Show a condensed summary of each blastdb, to list out all blastdbs
'''
class BlastDbListSerializer(serializers.ModelSerializer):
    sequences = BlastDbSequenceEntryShortSerializer(many=True, read_only=True)

    class Meta:
        model = BlastDb
        fields = ['id', 'custom_name', 'description', 'locked', 'sequences']

'''
Show a summary of a sequence, to be shown when only its blastdb or blastrun is shown.
'''
class BlastDbSequenceEntrySerializer(serializers.ModelSerializer):
    class Meta:
        model = NuccoreSequence
        fields = ['id', 'accession_number', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher', 'lat_lon',  'type_material', 'dna_sequence']

'''
Show detailed information about a specific blastdb
'''
class BlastDbSerializer(serializers.ModelSerializer):
    sequences = BlastDbSequenceEntrySerializer(many=True, read_only=True)

    class Meta:
        model = BlastDb
        fields = ['id', 'custom_name', 'description', 'locked', 'sequences']

'''
Required fields for adding a sequence accession to the database
'''
class NuccoreSequenceAddSerializer(serializers.ModelSerializer):
    class Meta:
        model = NuccoreSequence
        fields = ['accession_number']

'''
Serialize data when receiving a request to bulk add sequences to a database
'''
class NuccoreSequenceBulkAddSerializer(serializers.Serializer):
    accession_numbers = serializers.ListField(
        child=NuccoreSequenceAddSerializer()
    )

'''
Serialize a hit to be returned with a list of all hits
'''
class HitSerializer(serializers.ModelSerializer):
    db_entry = BlastDbSequenceEntrySerializer(many=False, read_only=True)

    class Meta:
        model = Hit
        fields = ['db_entry', 'owner_run', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']

'''
Show a summary of a sequence, when a query registers a hit on it.
'''
class NuccoreSequenceHitSerializer(serializers.ModelSerializer):
    class Meta:
        model = NuccoreSequence
        fields = ['accession_number', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher', 'type_material', 'lat_lon']


'''
Serialize a hit to be returned with a list of all hits
'''
class HitEntrySerializer(serializers.ModelSerializer):
    db_entry = NuccoreSequenceHitSerializer(many=False, read_only=True)

    class Meta:
        model = Hit
        fields = ['db_entry', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']

'''
Show a sequence submitted in a run
'''
class BlastQuerySequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = BlastQuerySequence
        fields = ['definition', 'query_sequence']

'''
Required fields for submitting a blast run
'''
class BlastRunRunSerializer(serializers.ModelSerializer):
    queries = BlastQuerySequenceSerializer(many=True, read_only=True)

    class Meta:
        model = BlastRun
        fields = ['id', 'job_name', 'queries', 'create_hit_tree', 'create_db_tree']

'''
Show results of a blast run
'''
class BlastRunSerializer(serializers.ModelSerializer):
    db_used = BlastDbShortSerializer(many=False, read_only=True)
    hits = HitEntrySerializer(many=True, read_only=True)
    queries = BlastQuerySequenceSerializer(many=True, read_only=True)    

    class Meta:
        model = BlastRun    
        fields = ['id', 'job_name', 'queries', 'db_used', 'runtime', 'job_status', 'job_start_time', 'job_end_time', 'job_error_time', 'hits', 'create_hit_tree', 'hit_tree', 'alignment_job_id', 'create_db_tree', 'db_tree', 'complete_alignment_job_id']

'''
Show status of a blast run
'''
class BlastRunStatusSerializer(serializers.ModelSerializer):    
    class Meta:
        model = BlastRun    
        fields = ['id', 'job_name', 'runtime', 'job_status', 'job_start_time', 'job_end_time', 'job_error_time']


