from rest_framework import serializers
from barcode_blastn.models import BlastRun, NuccoreSequence, BlastDb

class NuccoreSequenceSerializer(serializers.ModelSerializer):
    owner_database = serializers.PrimaryKeyRelatedField(many=False, queryset = BlastDb.objects.all())

    class Meta:
        model = NuccoreSequence
        fields = ['id', 'owner_database', 'accession_number', 'definition', 'organism', 'organelle', 'mol_type', 'isolate', 'country', 'specimen_voucher', 'dna_sequence', 'translation', 'created']

class BlastDbSerializer(serializers.ModelSerializer):
    sequences = serializers.PrimaryKeyRelatedField(many=True, read_only=True)

    class Meta:
        model = BlastDb
        fields = ['id', 'custom_name', 'sequences']

class BlastRunSerializer(serializers.ModelSerializer):

    db_used = serializers.PrimaryKeyRelatedField(many=False, queryset=BlastRun.objects.all())
    
    class Meta:
        model = BlastRun
        fields = ['id', 'job_name', 'query_sequence', 'db_used']
