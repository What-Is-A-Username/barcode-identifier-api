from datetime import datetime
import json
import os
from typing import Union
from django.contrib.auth.models import User
from rest_framework import serializers
from barcode_blastn.models import Annotation, BlastQuerySequence, BlastRun, DatabaseShare, Hit, Library, NuccoreSequence, BlastDb, TaxonomyNode
from barcode_blastn.validators import QueryFileValidator

library_title = 'Reference Library'
blast_db_title = 'BLAST Database Version'
nuccore_title = 'GenBank Accession'
hit_title = 'BLASTN hit'
query_title = 'Query Sequence'
run_title = 'Run'
share_title = 'Share Details'
annotation_title = 'User Annotation'

class TaxonomyNodeSerializer(serializers.ModelSerializer):
    class Meta:
        model = TaxonomyNode
        fields = ['id', 'scientific_name']

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
        fields = ['id', 'description', 'version_number', 'custom_name', 'sequence_count', 'created', 'locked']
        example = []
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "March 2022 Update",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "library": LibraryShortSerializer.Meta.example,
            "version_number": "1.1.1",
            "sequence_count": 2,
            "created": datetime.now(),
            "locked": True
        }

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
        fields = ['id', 'description', 'library', 'custom_name', 'version_number']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "custom_name": "Newly Sequenced Species Reference Library",
            "library": LibraryShortSerializer.Meta.example,
            "version_number": "2.1.1"
        }
        tags = ['DBS']

class AnnotationPosterSerializer(serializers.ModelSerializer):
    f'''
    Show basic information about who posted an annotation
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in ['username']:
            self.fields[name].read_only = True

    class Meta:
        model = User
        ref_name = f'{annotation_title} poster'
        fields = ['username']
        example = {
            'username': 'John Doe',
        }

class AnnotationSerializer(serializers.ModelSerializer):
    '''
    Show annotations when a list of annotations is returned
    '''
    # replies = RecursiveField(many=True)
    poster = AnnotationPosterSerializer(read_only=True)

    class Meta:
        model = Annotation
        ref_name = annotation_title
        # fields = ['replies', 'sequence', 'poster', 'timestamp', 'annotation_type', 'comment']
        fields = ['sequence', 'poster', 'timestamp', 'annotation_type', 'comment']

class NuccoreSequenceSerializer(serializers.ModelSerializer):
    f"""
    Complete information about a {nuccore_title} including the {blast_db_title} it belongs to. 
    """
    owner_database = BlastDbShortSerializer(many=False, read_only=True)
    taxon_genus = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_family = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_species = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_order = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_class = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_phylum = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_kingdom = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_superkingdom = TaxonomyNodeSerializer(many=False, read_only=True)
    annotations = AnnotationSerializer(many=True, read_only=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = False

    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title + ' summary'
        fields = ['id', 'annotations', 'owner_database', 'accession_number', 'version', 'definition', 'organism', 'organelle', 'isolate', 'country', 'specimen_voucher', 'lat_lon', 'dna_sequence', 'translation', 'type_material', 'created', 'genbank_modification_date', 'taxid', 'taxonomy', 'title', 'journal', 'authors', 'taxonomy', 'taxon_superkingdom', 'taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus', 'taxon_species']
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

class NuccoreSequenceBulkAddSerializer(serializers.Serializer):
    '''
    Serialize data when receiving a request to bulk add sequences to a database
    '''
    accession_numbers = serializers.ListField(
        child=serializers.CharField(),
        required=False
    )
    # Allow bulk addition of sequences by uploading file. Limit file size to 512KB
    accession_file = serializers.FileField(allow_empty_file=True, required=False,validators=[QueryFileValidator(max_size=524288)])

    class Meta:
        ref_name = nuccore_title
        example = {
            "accession_numbers": [
                "GU701771",
                "KF533332"
            ]
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

    taxon_genus = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_family = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_species = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_order = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_class = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_phylum = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_kingdom = TaxonomyNodeSerializer(many=False, read_only=True)
    taxon_superkingdom = TaxonomyNodeSerializer(many=False, read_only=True)
    annotations = AnnotationSerializer(many=True, read_only=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name in self.Meta.fields:
            self.fields[name].required = True
            self.fields[name].read_only = True

    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title
        fields = ['id', 'accession_number', 'version', 'organism', 'organelle', 'isolate', 'country', 'specimen_voucher', 'dna_sequence', 'lat_lon', 'type_material', 'created', 'updated', 'genbank_modification_date', 'taxonomy', 'taxon_superkingdom', 'taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus', 'taxon_species', 'annotations'] 
        example = {
            "id": "5100cbd8-2fda-4b42-8aa1-10ede078448b",
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
            "created": "2023-02-20T01:00:31.240961Z",
            "updated": "2023-02-20T01:00:31.240961Z"
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

    '''
    Create a new blastdb
    '''
    class Meta:
        model = BlastDb
        ref_name = blast_db_title + ' creation'
        fields = ['id', 'library', 'version_number', 'custom_name', 'description', 'locked', 'sequences', 'base', 'accession_numbers']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "2022 Official Release",
            "library": LibraryShortSerializer.Meta.example,
            "description": "A collection of new sequences from several species of interest.",
            "version_number": '254.2.1',
        }

class BlastDbEditSerializer(serializers.ModelSerializer):
    f'''
    Determine the editable fields for a specific {blast_db_title}.
    These fields are editable on the admin page and through PATCH requests.
    '''

    class Meta:
        model = BlastDb
        ref_name = 'Editable' + blast_db_title + ' Fields'
        fields = ['description', 'locked', 'custom_name']
        example = {
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "custom_name": "2022 Official Release",
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
        fields = ['id', 'library', 'custom_name', 'version_number', 'description', 'locked', 'sequences']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "owner": LibraryOwnerSerializer.Meta.example,
            "locked": True,
            "sequences": BlastDbSequenceEntrySerializer.Meta.example
        }

class BlastDbSequenceExportSerializer(BlastDbSequenceEntrySerializer):
    '''
    Gather detailed information for exporting the sequence to file when
    its database is exported.
    It includes all fields returned by the superclass (BlastDbSequenceEntrySerializer), but adds publication information
    to the exported file.
    '''
    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title + ' export'
        fields = ['accession_number', 'version', 'organism', 'organelle', 'isolate', 'country', 'specimen_voucher', 'dna_sequence', 'lat_lon', 'type_material', 'created', 'updated', 'genbank_modification_date', 'taxonomy', 'taxon_superkingdom', 'taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus', 'taxon_species', 'authors', 'title', 'journal']

class BlastDbExportSerializer(BlastDbSerializer):
    '''
    Gather detailed information for exporting the database to file.
    It includes all data otherwise added by BlastDbSerializer, but
    includes more information for each sequence.
    '''
    sequences = BlastDbSequenceExportSerializer(many=True, read_only=True)

class LibrarySerializer(serializers.ModelSerializer):
    f'''
    Show detailed information about a specific {library_title}
    '''
    owner = LibraryOwnerSerializer(many=False, read_only=True)
    latest = BlastDbTinySerializer(many=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        read_only_fields = ['id', 'owner']
        for name in read_only_fields:
            self.fields[name].read_only = True
        for name in ['id', 'custom_name', 'description', 'owner', 'public']:
            self.fields[name].required = True

    class Meta:
        model = Library
        ref_name = library_title
        fields = ['id', 'custom_name', 'description', 'owner', 'public', 'latest']
        example = {
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "custom_name": "Neotropical electric knifefish",
            "description": "This BLAST database is a collection of barcodes from 167 species of Neotropical electric knifefish (Teleostei: Gymnotiformes) which was presented by Janzen et al. 2022. All sequences and related feature data are updated daily at midnight (UTC) from NCBI's Genbank database.",
            "owner": LibraryOwnerSerializer.Meta.example,
            "latest": BlastDbTinySerializer.Meta.example,
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
        fields = ['id', 'version_number', 'custom_name', 'sequence_count', 'description', 'locked']
        example = [{
            "id": "66855f2c-f360-4ad9-8c98-998ecb815ff5",
            "version_number": '254.2.1',
            "sequence_count": 167,
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

    annotations = AnnotationSerializer(many=True, read_only=True)

    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title + ' summary'
        fields = ['id', 'accession_number', 'version', 'definition', 'organism', 'country', 'specimen_voucher', 'type_material', 'lat_lon', 'annotations']
        example = {
            "db_entry": {
                "accession_number": "ON303423",
                "version": "ON303423.1",
                "definition": "Porotergus duende isolate 2916 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial",
                "organism": "Porotergus duende",
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
    hits = HitEntrySerializer(many=True, read_only=True)

    class Meta:
        model = BlastQuerySequence
        ref_name = query_title
        fields = ['definition', 'query_sequence', 'hits', 'results_species_name', 'accuracy_category', 'original_species_name', 'write_tree_identifier']
        example = {
            "definition": "Steatogenys_elegans isolate 8807 cytochrome c oxidase subunit I (COI) gene, partial cds; mitochondrial",
            "query_sequence": "GGCACCCTTTATATAGTGTTTGGTGCCTGAGCCGGAATGGTTGGCACGGCCTTAAGCCTCCTTATTCGAGCCGAGCTAAGCCAACCCGGGGCCCTAATGGGTGATGACCAGATTTACAATGTTA"
        }

class BlastRunRunSerializer(serializers.ModelSerializer):
    '''
    Required fields for submitting a blast run
    '''
    query_sequence = serializers.CharField(allow_blank=False, min_length=10)
    query_file = serializers.FileField(max_length=2621440, validators=[QueryFileValidator()])

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['query_file'].required = False
        self.fields['query_sequence'].required = False
        self.fields['id'].read_only = True 

    class Meta:
        model = BlastRun
        ref_name = run_title + ' submission'
        fields = ['id', 'job_name', 'query_sequence', 'create_hit_tree', 'create_db_tree', 'query_file']

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
    queries = BlastQuerySequenceSerializer(many=True, read_only=True)    

    class Meta:
        model = BlastRun    
        ref_name = run_title
        fields = ['id', 'job_name', 'queries', 'db_used', 'start_time', 'status', 'received_time', 'start_time', 'end_time', 'error_time', 'create_hit_tree', 'hit_tree', 'alignment_job_id', 'create_db_tree', 'db_tree', 'complete_alignment_job_id']
        example = load_run_example()

class BlastQuerySequenceTaxonomySerializer(serializers.ModelSerializer):
    '''
    Retrieve query sequences, without also retrieving hits, to gather run data
    pertaining to taxonomic assignment made on the query sequence, as required by 
    BlastRunTaxonomyCSVRenderer as similar.
    '''
    class Meta:
        model = BlastQuerySequence
        fields = ['definition', 'query_sequence', 'results_species_name', 'accuracy_category', 'original_species_name', 'write_tree_identifier']

class BlastRunTaxonomySerializer(serializers.ModelSerializer):
    '''
    Load run data pertaining to taxonomic assignments made, as required by
    BlastRunTaxonomyCSVRenderer and similar.
    '''
    queries = BlastQuerySequenceTaxonomySerializer(many=True, read_only=True)

    class Meta:
        model = BlastRun
        ref_name = f'{run_title} taxonomy'
        fields = ['id', 'job_name', 'queries']

class BlastRunStatusSerializer(serializers.ModelSerializer):    
    '''
    Show status of a blast run
    '''
    class Meta:
        model = BlastRun    
        ref_name = run_title + ' status'
        fields = ['id', 'job_name', 'received_time', 'status', 'start_time', 'end_time', 'error_time']
        example = {
            "id": "2e5898a3-14da-4e7f-9599-ba1ef35f1e7a",
            "job_name": "two sequences",
            "received_time": "2023-02-21T03:26:32.323521Z",
            "start_time": "2023-02-21T03:27:11.323521Z",
            "status": "FIN",
            "start_time": "2023-02-21T03:27:12.427193Z",
            "end_time": "2023-02-21T03:28:02.102346Z",
            "error_time": None,
        }

class UserSerializer(serializers.ModelSerializer):
    '''
    Show basic user information.
    '''
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

class NuccoreAnnotationSerializer(serializers.ModelSerializer):
    '''
    Show a very brief summary of the sequence which is the subject of
    a user annotation.
    '''
    class Meta:
        model = NuccoreSequence
        ref_name = nuccore_title + ' entry'
        fields = ['version', 'organism']

class RecursiveField(serializers.Serializer):
    def to_representation(self, value):
        serializer = self.parent.parent.__class__(value, context=self.context)
        return serializer.data

class AnnotationPoster(serializers.ModelSerializer):
    '''
    Required data to post/create an annotation
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        required_fields = ['comment', 'annotation_type']
        readonly_fields = ['id', 'timestamp', 'poster']
        for field in readonly_fields:
            self.fields[field].read_only = True
            self.fields[field].required = True
        for field in required_fields:
            self.fields[field].required = True

    class Meta:
        model = Annotation
        ref_name = annotation_title
        fields = ['id', 'poster', 'timestamp', 'annotation_type', 'comment']
        