from unittest.util import _MAX_LENGTH
from django.db import models
from django.utils.translation import gettext_lazy as _
import uuid

class BlastDb(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # The creation datetime of this database
    created = models.DateTimeField(auto_now_add=True)
    # The user-customized title
    custom_name = models.CharField(max_length=255)
    # Locked
    locked = models.BooleanField(default=False)
    # Short description of the database
    description = models.CharField(max_length=1024, blank=True, default='')

    class Meta:
        ordering = ['custom_name']

class NuccoreSequence(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

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
    lat_lon = models.CharField(max_length=64, blank=True, default='')

    class Meta:
        ordering = ['accession_number']



class BlastRun(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # Reference to the database used
    db_used = models.ForeignKey(BlastDb, related_name='usages', on_delete=models.CASCADE)
    # When was the blastn run?
    runtime = models.DateTimeField(auto_now_add=True)
    # Job name
    job_name = models.CharField(max_length=255, blank=True, default='')

    # Query sequence
    # TODO: remove query_sequence reference
    query_sequence = models.TextField(max_length=10000, blank=True, default='')

    class JobStatus(models.TextChoices):
        UNKNOWN = 'UNK', _('UNKNOWN')
        DENIED = 'DEN', _('DENIED')
        QUEUED = 'QUE', _('QUEUED')
        STARTED = 'STA', _('RUNNING')
        ERRORED = 'ERR', _('ERRORED')
        FINISHED = 'FIN', _('FINISHED')

    job_status = models.CharField(max_length=3,choices=JobStatus.choices, default=JobStatus.UNKNOWN)

    job_start_time = models.DateTimeField(blank=True, null=True)
    job_end_time = models.DateTimeField(blank=True, null=True)

    # RESULTS
    # Blast version
    blast_version = models.TextField(max_length=100, blank=True, default='')
    # Database name 

    # Error
    errors = models.TextField(max_length=10000, blank=True, default='')

class BlastQuerySequence(models.Model):
    owner_run = models.ForeignKey(BlastRun, related_name='queries', on_delete=models.CASCADE)
    definition = models.CharField(max_length=255)
    query_sequence = models.CharField(max_length=10000)

class Hit(models.Model):
    owner_run = models.ForeignKey(BlastRun, related_name='hits', on_delete=models.CASCADE)
    db_entry = models.ForeignKey(NuccoreSequence, on_delete=models.CASCADE)

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