from datetime import datetime
from unittest.util import _MAX_LENGTH
from django.db import models
from django.utils.translation import gettext_lazy as _
import uuid

class BlastDb(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, help_text='Unique identifier of this BLAST database')

    # The creation datetime of this database
    created = models.DateTimeField(auto_now_add=True, help_text='Date and time at which database was created')
    # The user-customized title
    custom_name = models.CharField(max_length=255, help_text='Name of BLAST database')
    # Locked
    locked = models.BooleanField(default=False, help_text='Is editing of entry set (adding/removing) in the database locked?')
    # Short description of the database
    description = models.CharField(max_length=1024, blank=True, default='', help_text='Description of BLAST database')

    def __str__(self) -> str:
        return f'"{self.custom_name}" Database ({str(self.id)})'

    class Meta:
        ordering = ['custom_name']
        verbose_name = 'BLAST Database'
        verbose_name_plural = 'BLAST Databases'

class NuccoreSequence(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, help_text='Unique identifier of this sequence entry')

    owner_database = models.ForeignKey(BlastDb, related_name='sequences',
        on_delete=models.CASCADE, help_text='The curated database to which this sequence was added')

    created = models.DateTimeField(auto_now_add=True, help_text='Date and time at which record was last updated from GenBank')
    accession_number = models.CharField(max_length=255, help_text='Accession number on GenBank')
    uid = models.CharField(max_length=2048, blank=True, default='', help_text='Obselete UUID')
    definition = models.CharField(max_length=255, blank=True, default='', help_text='The definition line')
    organism = models.CharField(max_length=255, blank=True, default='', help_text='Scientific name of source organism')
    organelle = models.CharField(max_length=255, blank=True, default='', help_text='Organelle of the source')
    isolate = models.CharField(max_length=255, blank=True, default='', help_text='Isolate of the source specimen')
    country = models.CharField(max_length=255, blank=True, default='', help_text='Origin country of the source specimen')
    specimen_voucher = models.CharField(max_length=150, blank=True, default='', help_text = 'Specimen voucher of the source specimen')
    dna_sequence = models.TextField(max_length=10000, blank=True, default='', help_text='Sequence data')
    translation = models.TextField(max_length=10000, blank=True, default='', help_text='Amino acid translation corresponding to the coding sequence')
    lat_lon = models.CharField(max_length=64, blank=True, default='', help_text='Latitude and longitude from which specimen originated')
    type_material = models.CharField(max_length=255, blank=True, default='', help_text='Specimen type of the source')

    def __str__(self) -> str:
        return f'{self.accession_number}, {str(self.organism)} ({str(self.id)})'

    class Meta:
        ordering = ['accession_number']
        verbose_name = 'GenBank Accession'
        verbose_name_plural = 'GenBank Accessions'

class BlastRun(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, help_text='Unique identifier of the run')

    # Reference to the database used
    db_used = models.ForeignKey(BlastDb, related_name='usages', on_delete=models.CASCADE, help_text='The curated BLAST database against which the query BLAST was run.')
    # When was the run request first received (i.e. added to queue)
    runtime = models.DateTimeField(auto_now_add=True, help_text='Date and time when run first received by server')
    # Job name
    job_name = models.CharField(max_length=255, blank=True, default='', help_text='Job name given by run submission')

    # Perform alignment and construct NJ tree of query sequences + hits?
    create_hit_tree = models.BooleanField(default=True, help_text='Perform alignment and construct "hit tree" of query sequences and hits?')
    # Job ID of alignment using query + hit sequences
    alignment_job_id = models.CharField(max_length=100, blank=True, default='', help_text='External job ID used to construct hit tree')
    
    # Perform alignment and construct NJ tree of query sequences + all DB sequences?
    create_db_tree = models.BooleanField(default=True, help_text='Perform alignment and construct "database tree" of query sequences and all database sequences?')
    # Job ID for alignment using query + all database sequences
    complete_alignment_job_id = models.CharField(max_length=100, blank=True, default='', help_text='External job ID used to construct database tree')

    # Newick string of tree with query sequences + hits 
    hit_tree = models.TextField(blank=True, default='', help_text='Newick/phylip tree string of hit tree.')
    # Newick string of tree with query sequences + all database sequences
    db_tree = models.TextField(blank=True, default='', help_text='Newick/phylip tree string of database tree.')

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
    job_status = models.CharField(max_length=3,choices=JobStatus.choices, default=JobStatus.UNKNOWN, help_text='Current status of the job')
    # Time that the job started running
    job_start_time = models.DateTimeField(blank=True, null=True, help_text='Date and time when job first started running')
    # Time that job successfully finished
    job_end_time = models.DateTimeField(blank=True, null=True, help_text='Date and time when job successfully finished running')
    # Time that job errored
    job_error_time = models.DateTimeField(blank=True, null=True, help_text='Date and time when job encountered an error')

    # Blast version
    blast_version = models.TextField(max_length=100, blank=True, default='', help_text='Version of BLASTn used')

    # Error for internal debugging
    errors = models.TextField(max_length=10000, blank=True, default='', help_text='Error message text')

    def __str__(self) -> str:
        return f'Run {self.id}'

    class Meta:
        ordering = ['runtime']
        verbose_name = 'BLASTN Run'
        verbose_name_plural = 'BLASTN Runs'

class BlastQuerySequence(models.Model):
    owner_run = models.ForeignKey(BlastRun, related_name='queries', on_delete=models.CASCADE, help_text='Job/run in which this query entry appeared')
    definition = models.CharField(max_length=255, help_text='Definition line')
    query_sequence = models.CharField(max_length=10000, help_text='Sequence text')

    class Meta:
        verbose_name = 'BLASTN Query Sequence'
        verbose_name_plural = 'BLASTN Query Sequences'

class Hit(models.Model):
    owner_run = models.ForeignKey(BlastRun, related_name='hits', on_delete=models.CASCADE, help_text='Run in which this hit appeared')
    db_entry = models.ForeignKey(NuccoreSequence, on_delete=models.CASCADE, help_text='BLAST database used in the run')

    query_accession_version = models.CharField(max_length=128, help_text='Sequence identifier of query sequence')
    subject_accession_version = models.CharField(max_length=128, help_text='Sequence identifier of sequence in database')
    percent_identity = models.DecimalField(max_digits=6, decimal_places=3, help_text='Percent identity')
    alignment_length = models.IntegerField(help_text='Alignment length')
    mismatches = models.IntegerField(help_text='Number of mismatches')
    gap_opens = models.IntegerField(help_text='Number of Gap openings')
    query_start = models.IntegerField(help_text='Start of alignment in query')
    query_end = models.IntegerField(help_text='End of alignment in query')
    sequence_start = models.IntegerField(help_text='Start of alignment in subject')
    sequence_end = models.IntegerField(help_text='End of alignment in subject')
    evalue = models.DecimalField(max_digits=110, decimal_places=100, help_text='Expect value')
    bit_score = models.DecimalField(max_digits=110, decimal_places=100, help_text='Bit score')

    class Meta:
        ordering = ['percent_identity']
        verbose_name = 'BLASTN Run Hit'
        verbose_name_plural = 'BLASTN Run Hits'