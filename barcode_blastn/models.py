from datetime import datetime
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
    type_material = models.CharField(max_length=255, blank=True, default='')

    class Meta:
        ordering = ['accession_number']

class BlastRun(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # Reference to the database used
    db_used = models.ForeignKey(BlastDb, related_name='usages', on_delete=models.CASCADE)
    # When was the run request first received (i.e. added to queue)
    runtime = models.DateTimeField(auto_now_add=True)
    # Job name
    job_name = models.CharField(max_length=255, blank=True, default='')

    # Perform alignment and construct NJ tree of query sequences + hits?
    create_hit_tree = models.BooleanField(default=True)
    # Job ID of alignment using query + hit sequences
    alignment_job_id = models.CharField(max_length=100, blank=True, default='')
    
    # Perform alignment and construct NJ tree of query sequences + all DB sequences?
    create_db_tree = models.BooleanField(default=True)
    # Job ID for alignment using query + all database sequences
    complete_alignment_job_id = models.CharField(max_length=100, blank=True, default='')
    # Newick string of tree with query sequences + hits 
    hit_tree = models.TextField(blank=True, default='')
    # Newick string of tree with query sequences + all database sequences
    db_tree = models.TextField(blank=True, default='')

    class JobStatus(models.TextChoices):
        UNKNOWN = 'UNK', _('UNKNOWN')
        DENIED = 'DEN', _('DENIED')
        QUEUED = 'QUE', _('QUEUED')
        STARTED = 'STA', _('RUNNING')
        ERRORED = 'ERR', _('ERRORED')
        FINISHED = 'FIN', _('FINISHED')

    def throw_error(self, debug_error_message: str = ''):
        '''Designate the current run to error and add debug_error_message string to errors.
        '''
        self.job_error_time = datetime.now()
        self.errors = ('\n' + debug_error_message) if len(self.errors) > 0 else debug_error_message
        self.job_status = self.JobStatus.ERRORED
        self.save()

    # What is the current status of the job?
    job_status = models.CharField(max_length=3,choices=JobStatus.choices, default=JobStatus.UNKNOWN)
    # Time that the job started running
    job_start_time = models.DateTimeField(blank=True, null=True)
    # Time that job successfully finished
    job_end_time = models.DateTimeField(blank=True, null=True)
    # Time that job errored
    job_error_time = models.DateTimeField(blank=True, null=True)

    # Blast version
    blast_version = models.TextField(max_length=100, blank=True, default='')

    # Error for internal debugging
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