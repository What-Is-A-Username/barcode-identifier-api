from rest_framework import serializers
from barcode_blastn.models import BlastRun, Hit, NuccoreSequence, BlastDb

class BlastDbShortSerializer(serializers.ModelSerializer):
    class Meta:
        model = BlastDb
        fields = ['id', 'custom_name']

class NuccoreSequenceSerializer(serializers.ModelSerializer):
    owner_database = BlastDbShortSerializer(many=False, read_only=True)

    class Meta:
        model = NuccoreSequence
        fields = ['owner_database', 'accession_number', 'definition', 'organism', 'organelle', 'mol_type', 'isolate', 'country', 'specimen_voucher', 'dna_sequence', 'translation', 'created']

class BlastDbSequenceEntrySerializer(serializers.ModelSerializer):
    class Meta:
        model = NuccoreSequence
        fields = ['id', 'accession_number', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher']

class BlastDbSerializer(serializers.ModelSerializer):
    sequences = BlastDbSequenceEntrySerializer(many=True, read_only=True)

    class Meta:
        model = BlastDb
        fields = ['id', 'custom_name', 'locked', 'sequences']
    
class NuccoreSequenceAddSerializer(serializers.ModelSerializer):
    class Meta:
        model = NuccoreSequence
        fields = ['accession_number']

class HitSerializer(serializers.ModelSerializer):
    class Meta:
        model = Hit
        fields = ['owner_run', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue', 'bit_score']

class BlastRunRunSerializer(serializers.ModelSerializer):
    class Meta:
        model = BlastRun
        fields = ['id', 'job_name', 'query_sequence']

class BlastRunSerializer(serializers.ModelSerializer):

    db_used = BlastDbShortSerializer(many=False, read_only=True)
    hits = HitSerializer(many=True, read_only=True)
    
    class Meta:
        model = BlastRun
        fields = ['id', 'job_name', 'query_sequence', 'db_used', 'hits']

