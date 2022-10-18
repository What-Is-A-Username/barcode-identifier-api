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
