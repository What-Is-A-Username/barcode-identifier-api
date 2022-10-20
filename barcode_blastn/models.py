from django.db import models

class BlastDb(models.Model):
    # The creation datetime of this database
    created = models.DateTimeField(auto_now_add=True)
    # The user-customized title
    custom_name = models.CharField(max_length=255)
    # Locked
    locked = models.BooleanField(default=False)

    class Meta:
        ordering = ['custom_name']

class NuccoreSequence(models.Model):
    owner_database = models.ForeignKey(BlastDb, related_name='sequences',
        on_delete=models.CASCADE)

    created = models.DateTimeField(auto_now_add=True)
    accession_number = models.CharField(max_length=255)
    uid = models.CharField(max_length=2048, blank=True, default='')
    definition = models.CharField(max_length=255, blank=True, default='')
    organism = models.CharField(max_length=255, blank=True, default='')
    organelle = models.CharField(max_length=255, blank=True, default='')
    mol_type = models.CharField(max_length=255, blank=True, default='')
    isolate = models.CharField(max_length=255, blank=True, default='')
    country = models.CharField(max_length=255, blank=True, default='')
    specimen_voucher = models.CharField(max_length=150, blank=True, default='')
    dna_sequence = models.TextField(max_length=10000, blank=True, default='')
    translation = models.TextField(max_length=10000, blank=True, default='')

    class Meta:
        ordering = ['accession_number']

class BlastRun(BlastDb):
    # Reference to the database used
    db_used = models.ForeignKey(BlastDb, related_name='usages', on_delete=models.CASCADE)
    # When was the blastn run?
    runtime = models.DateTimeField(auto_now_add=True)
    # Job name
    job_name = models.CharField(max_length=255, blank=True, default='')
    # Query sequence
    query_sequence = models.TextField(max_length=10000, blank=True, default='')

    # RESULTS
    # Blast version
    blast_version = models.TextField(max_length=100, blank=True, default='')
    # Database name 

    # Error
    errors = models.TextField(max_length=10000, blank=True, default='')

class Hit(models.Model):
    owner_run = models.ForeignKey(BlastRun, related_name='hits', on_delete=models.CASCADE)

    query_accession_version = models.CharField(max_length=128)
    subject_accession_version = models.CharField(max_length=128)
    percent_identity = models.DecimalField(max_digits=6, decimal_places=3)
    alignment_length = models.IntegerField()
    mismatches = models.IntegerField()
    gap_opens = models.IntegerField()
    query_start = models.IntegerField()
    query_end = models.IntegerField()
    sequence_start = models.IntegerField()
    sequence_end = models.IntegerField()
    evalue = models.DecimalField(max_digits=110, decimal_places=100)
    bit_score = models.DecimalField(max_digits=110, decimal_places=100)