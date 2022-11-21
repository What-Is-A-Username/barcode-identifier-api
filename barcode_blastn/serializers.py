from rest_framework import serializers
from barcode_blastn.models import BlastRun, Hit, NuccoreSequence, BlastDb

class BlastDbShortSerializer(serializers.ModelSerializer):
    class Meta:
        model = BlastDb
        fields = ['id', 'custom_name', 'description']

class NuccoreSequenceSerializer(serializers.ModelSerializer):
    owner_database = BlastDbShortSerializer(many=False, read_only=True)

    class Meta:
        model = NuccoreSequence
        fields = ['owner_database', 'accession_number', 'definition', 'organism', 'organelle', 'mol_type', 'isolate', 'country', 'specimen_voucher', 'lat_lon', 'dna_sequence', 'translation', 'created']

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
        fields = ['id', 'accession_number', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher', 'lat_lon', 'dna_sequence']

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
Serialize a hit to be returned with a list of all hits
'''
class HitSerializer(serializers.ModelSerializer):
    db_entry = BlastDbSequenceEntrySerializer(many=False, read_only=True)

    class Meta:
        model = Hit
        fields = ['db_entry', 'owner_run', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']

'''
Serialize a hit to be returned with a list of all hits
'''
class HitEntrySerializer(serializers.ModelSerializer):
    db_entry = BlastDbSequenceEntrySerializer(many=False, read_only=True)

    class Meta:
        model = Hit
        fields = ['db_entry', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']

'''
Required fields for submitting a blast run
'''
class BlastRunRunSerializer(serializers.ModelSerializer):
    class Meta:
        model = BlastRun
        fields = ['id', 'job_name', 'query_sequence']

'''
Show results of a blast run
'''
class BlastRunSerializer(serializers.ModelSerializer):

    db_used = BlastDbShortSerializer(many=False, read_only=True)
    hits = HitEntrySerializer(many=True, read_only=True)
    
    class Meta:
        model = BlastRun    
        fields = ['id', 'job_name', 'query_sequence', 'db_used', 'runtime', 'job_status', 'job_start_time', 'job_end_time', 'hits']

'''
Show status of a blast run
'''
class BlastRunStatusSerializer(serializers.ModelSerializer):    
    class Meta:
        model = BlastRun    
        fields = ['id', 'job_name', 'runtime', 'job_status', 'job_start_time', 'job_end_time']


