from rest_framework import serializers
from barcode_tree.models import ResultTree

class ResultTreeCreatorSerializer(serializers.ModelSerializer):
    dummy_field = serializers.CharField()

    class Meta:
        model = ResultTree
        fields = ['dummy_field']

class ResultTreeDetailSerializer(serializers.ModelSerializer):
    class Meta:
        model = ResultTree
        fields = ['owner_run', 'alignment_job_id', 'internal_status']