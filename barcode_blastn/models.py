import uuid
from collections import namedtuple
from datetime import datetime
from typing import Any, List, Union

from django.contrib.auth.models import AbstractBaseUser, AnonymousUser, User
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from django.utils.translation import gettext_lazy as _
from simple_history.models import HistoricalRecords

from barcode_blastn.database_permissions import DatabasePermissions

class BlastDbHistoricalModel(models.Model):
    '''
    Abstract model for tracking the accession numbers or search
    terms used to add or remove sequences from the database.
    '''
    added = models.CharField(max_length=10000, blank=True, default='None')
    deleted = models.CharField(max_length=10000, blank=True, default='None')
    search_terms = models.CharField(max_length=10000, blank=True, default='None')
    filter_options = models.CharField(max_length=10000, blank=True, default='')

    class Meta:
        abstract = True

class TaxonomyNode(models.Model):
    class TaxonomyRank(models.TextChoices):
        SUPERKINGDOM = 'sk'
        KINGDOM = 'k'
        PHYLUM = 'p'
        CLASS = 'c'
        ORDER = 'o'
        FAMILY = 'f'
        GENUS = 'g'
        SPECIES = 's'        

    scientific_name = models.CharField(help_text='Scientific name.', max_length=512)
    rank = models.CharField(help_text='Taxonomic rank', choices=TaxonomyRank.choices, max_length=128)
    last_updated = models.DateTimeField(auto_now=True)

    def __str__(self) -> str:
        return f'{self.scientific_name}: {self.id} ({self.rank})'

    class Meta:
        ordering = ['scientific_name']
        verbose_name = 'NCBI Taxonomy Node'
        verbose_name_plural = 'NCBI Taxonomy Nodes'

class LibraryManager(models.Manager):
    '''
    Model manager for the BlastDb class.
    '''
    def editable(self, user: Union[AbstractBaseUser, AnonymousUser]):
        '''
        Return a queryset of Reference Libraries that are editable by the given user.
        Editability is given if:
        - user is a superuser 
        - user is owner
        - user given explicit edit permission

        Returns an empty queryset if user is not authenticated, the full queryset
        if user is a superuser, and a queryset of all BlastDbs editable by the
        user.
        '''
        if not user.is_authenticated:
            # public users can never edit a database
            return super().none()
        elif isinstance(user, User) and user.is_superuser:
            # give access to all databases if superuser
            return super().get_queryset().all()
        else:
            return super().get_queryset().filter(
                # include databases owned by user
                models.Q(owner=user) |
                # include those with access explicitly given to 
                models.Q(shares=user, databaseshare__permission_level__in=[DatabasePermissions.CAN_EDIT_DB])
            )
        
    def viewable(self, user: Union[AbstractBaseUser, AnonymousUser]):
        '''
        Retrieve a set of databases that the user can view. A user can view if:
        - the reference library is public, and either not explicitly denied or unauthenticated
        - user is superuser
        - user is owner
        - user is given explicit view, run or edit permission
        '''
        if not user.is_authenticated:
            # only return public databases for public users
            return super().get_queryset().filter(public=True)
        elif isinstance(user, User) and user.is_superuser:
            # give access to all databases if superuser
            return super().get_queryset().all()
        else:
            return super().get_queryset().exclude(
                # exclude databases that are public but explicitly denied to user
                models.Q(public=True) &
                models.Q(shares=user, 
                         databaseshare__permission_level=DatabasePermissions.DENY_ACCESS)
            ).filter(
                # include databases owned by user
                models.Q(owner=user) |
                # include those that are public
                models.Q(public=True) |
                # include those with access explicitly given to 
                models.Q(shares=user,
                         databaseshare__permission_level__in=[
                            DatabasePermissions.CAN_VIEW_DB, DatabasePermissions.CAN_EDIT_DB, DatabasePermissions.CAN_RUN_DB
                         ])
            )

    def runnable(self, user: Union[AbstractBaseUser, AnonymousUser]):
        '''
        Retrieve a set of reference libraries that the user can run on. Permission is given if:
        - the user can view
        '''
        # Currently, if a user can view the database and the database is locked
        # they can run it
        return self.viewable(user) 

    def deletable(self, user: Union[AbstractBaseUser, AnonymousUser]):
        '''
        Retrieve a set of reference libraries that is deletable by the user. 
        A reference library is deletable if:
        - user is superuser
        '''
        if not user.is_authenticated:
            # public users cannot delete any databases
            return super().none()
        elif isinstance(user, User) and user.is_superuser:
            # A superuser can delete any database
            return super().get_queryset().all()
        else:
            return super().none()

class Library(models.Model):

    objects = LibraryManager()

    history = HistoricalRecords()

    class MarkerGenes(models.TextChoices):
        MARKER_CO1 = 'CO1', _('CO1')
        MARKER_18S = '18S', _('18S')
        MARKER_16S = '16S', _('16S')
        MARKER_12S = '12S', _('12S')
        MARKER_CYTB = 'CytB', _('Cytb')
        MARKER_ITS = 'ITS', _('ITS')
    
    marker_gene = models.CharField(choices=MarkerGenes.choices, max_length=8)

    # Original creator of the database
    owner = models.ForeignKey(User, on_delete=models.CASCADE)

    # Is the database accessible to the public for viewing and running?
    public = models.BooleanField(default=False, help_text='Is this reference library accessible to the public?')
    
    # Access permissions granted to users by the owner
    shares = models.ManyToManyField(User, related_name='permission', related_query_name='permissions', through='DatabaseShare')

    # Short description of the database
    description = models.CharField(max_length=1024, blank=True, default='', help_text='Description of contents and usage')

    # Unique identifier for the BLAST database
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, help_text='Unique identifier of this Reference Library')

    # The user-customized title
    custom_name = models.CharField(max_length=255, help_text='Name of Reference Library', unique=True)

    # Time of first creation
    created = models.DateTimeField(auto_now_add=True, help_text='Date and time at which reference library was first created')

    def __str__(self) -> str:
        return f'"{self.custom_name}" Library ({str(self.id)})'

    class Meta:
        ordering = ['custom_name']
        verbose_name = 'Reference Library'
        verbose_name_plural = 'Reference Libraries'

    def latest(self):
        '''Retrieve the latest blastdb that is published (locked=True) from a given library'''
        return BlastDb.objects.latest(self)

    def latest_version(self):
        '''Retrieve a string of the latest database that is published (locked=True) for this library'''
        db = self.latest()
        if db is None:
            return ""
        else:
            return f'{db.version_number()} ({db.custom_name})'

class BlastDbManager(models.Manager):
    def latest(self, library: Library):
        '''Retrieve the latest blastdb that is published (locked=True) from a given library'''
        return BlastDb.objects.filter(library=library, locked=True).last()

    def editable(self, user: User):
        '''Return a QuerySet of BLAST databases that are editable by the given user. BLAST databases are editable if their reference library are editable by the same user.'''
        libs: models.QuerySet[Library] = Library.objects.editable(user)
        return BlastDb.objects.filter(library__in=libs)

    def viewable(self, user: Union[AbstractBaseUser, AnonymousUser]):
        '''Return a QuerySet of BLAST databases that are viewable by the given user. BLAST databases are editable if their reference library are viewable by the same user.'''
        libs: models.QuerySet[Library] = Library.objects.viewable(user)
        return BlastDb.objects.filter(library__in=libs)

    def runnable(self, user: Union[AbstractBaseUser, AnonymousUser]):
        '''Return a QuerySet of BLAST databases that are runnable by the given user. BLAST databases are editable if their reference library are runnable by the same user.'''
        libs: models.QuerySet[Library] = Library.objects.runnable(user)
        return BlastDb.objects.filter(library__in=libs, locked=True)

    def deletable(self, user: User):
        '''Return a QuerySet of BLAST databases that are deletable by the given user. BLAST databases are editable if their reference library are deletable by the same user.'''
        libs: models.QuerySet[Library] = Library.objects.deletable(user)
        return BlastDb.objects.filter(library__in=libs)

class BlastDb(models.Model):
    objects = BlastDbManager()

    history = HistoricalRecords(bases=[BlastDbHistoricalModel])

    def save_without_historical_record(self, *args, **kwargs):
        '''Save database instance without creating a historical record'''
        # https://django-simple-history.readthedocs.io/en/latest/querying_history.html#save-without-a-historical-record
        self.skip_history_when_saving = True
        try:
            ret = self.save(*args, **kwargs)
        finally:
            del self.skip_history_when_saving
        return ret

    custom_name = models.CharField(max_length=255, help_text='Name of the database version')
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, help_text='Unique identifier of this database version of the library')

    library = models.ForeignKey(Library, on_delete=models.CASCADE, help_text='The reference library which this database is a version of.')

    # TODO: Reorganize version numbers, deciding whether to integrate GenBank version
    genbank_version = models.SmallIntegerField(default=0, validators=[MaxValueValidator(32767), MinValueValidator(0)], help_text='Version number reflective of changes in sequence data, accession.versions, and the set of accessions included.')
    major_version = models.SmallIntegerField(default=0, validators=[MaxValueValidator(32767), MinValueValidator(0)], help_text='Version number reflective of changes in important metadata such as source information, location, specimen identifiers.')
    minor_version = models.SmallIntegerField(default=0, validators=[MaxValueValidator(32767), MinValueValidator(0)], help_text='Version number reflective of minor changes such as database description, names, references.')

    def version_number(self) -> str:
        return f'{self.genbank_version}.{self.major_version}.{self.minor_version}'

    def sequence_count(self) -> int:
        return NuccoreSequence.objects.filter(owner_database=self).count()

    # The creation datetime of this version
    created = models.DateTimeField(auto_now_add=True, help_text='Date and time at which database was created.')
    # Last modified
    modified = models.DateTimeField(auto_now=True, help_text='Date and time at which database was last modified.')
    # Short description of the version
    description = models.CharField(max_length=1024, blank=True, default='', help_text='Description of this version')
    # Locked
    locked = models.BooleanField(default=False, help_text='Is editing of entry set (adding/removing) in the database locked?')

    blacklisted_accessions = models.JSONField(default=list, help_text='List of accession numbers and versions blacklisted, which will not be added to the database.')

    def __str__(self) -> str:
        if self.locked:
            return f'"{self.library.custom_name}", Version {self.version_number()}, ({self.custom_name}) ({self.id})'
        else:
            return f'"{self.library.custom_name}", Unpublished ({self.custom_name}), ({self.id})'

    class Meta:
        ordering = ['genbank_version', 'major_version', 'minor_version']
        verbose_name = 'BLAST Database Version'
        verbose_name_plural = 'BLAST Database Versions'

class DatabaseShare(models.Model):
    # The database these permissions apply to
    database = models.ForeignKey(Library, on_delete=models.CASCADE, help_text='The database these permissions apply to.')
    # User that the permissions apply to
    grantee = models.ForeignKey(User, on_delete=models.CASCADE, help_text='User that the permissions apply to.')
    # Permission level given to this relationship 
    permission_level = models.CharField(max_length=16, choices=DatabasePermissions.choices, default=DatabasePermissions.DENY_ACCESS, help_text='Access permissions')

    def __str__(self) -> str:
        return f'"{self.permission_level}" permission for "{self.grantee.username}" on library "{self.database.custom_name}"'

    class Meta:
        unique_together = ['database', 'grantee']
        verbose_name = 'BLAST Database Access Permission'
        verbose_name_plural = 'BLAST Database Access Permissions'
        ordering = ['database__custom_name', 'grantee__username']

class NuccoreSequenceManager(models.Manager):
    def viewable(self, user: User):
        '''
        Return a queryset of GenBank accessions that are viewable
        by the given user (i.e. accessions located in databases
        that the user can view)
        '''
        db_qs = BlastDb.objects.viewable(user)
        return self.get_queryset().filter(owner_database__in=db_qs)

class NuccoreSequence(models.Model):

    objects = NuccoreSequenceManager()

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, help_text='Unique identifier of this sequence entry')

    owner_database = models.ForeignKey(BlastDb, related_name='sequences',
        on_delete=models.CASCADE, help_text='The curated database to which this sequence was added')

    accession_number = models.CharField(max_length=255, help_text='Accession number on GenBank')
    version = models.CharField(max_length=63, help_text='The accession.version on GenBank')
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

    taxid = models.IntegerField(help_text='Taxonomic identifier of source organism.', blank=True, default=-2)

    taxon_species = models.ForeignKey(TaxonomyNode, related_name='species_member', related_query_name='species_members', help_text='Species taxa of the source organism', on_delete=models.CASCADE, null=True, blank=True) 
    taxon_genus = models.ForeignKey(TaxonomyNode, related_name='genus_member', related_query_name='genus_members', help_text='Genus of the source organism', on_delete=models.CASCADE, null=True, blank=True) 
    taxon_family = models.ForeignKey(TaxonomyNode, related_name='family_member', related_query_name='family_members', help_text='Family of the source organism', on_delete=models.CASCADE, null=True, blank=True) 
    taxon_order = models.ForeignKey(TaxonomyNode, related_name='order_member', related_query_name='order_members', help_text='Order of the source organism', on_delete=models.CASCADE, null=True, blank=True) 
    taxon_class = models.ForeignKey(TaxonomyNode, related_name='class_member', related_query_name='class_members', help_text='Class of the source organism', on_delete=models.CASCADE, null=True, blank=True) 
    taxon_phylum = models.ForeignKey(TaxonomyNode, related_name='phylum_member', related_query_name='phylum_members', help_text='Phylum of the source organism', on_delete=models.CASCADE, null=True, blank=True) 
    taxon_kingdom = models.ForeignKey(TaxonomyNode, related_name='kingdom_member', related_query_name='kingdom_members', help_text='Kingdom of the source organism', on_delete=models.CASCADE, null=True, blank=True) 
    taxon_superkingdom = models.ForeignKey(TaxonomyNode, related_name='superkingdom_member', related_query_name='superkingdom_members', help_text='Superkingdom of the source organism', on_delete=models.CASCADE, null=True, blank=True) 
    
    keywords = models.CharField(max_length=256, help_text='Comma-separated list of keywords on the sequence record.', blank=True, default='')
    title = models.CharField(max_length=512, help_text='Title of a publication by the authors of the sequence record.', blank=True, default='')
    journal = models.CharField(max_length=512, help_text='The journal of a publication, expressed with the MEDLINE abbreviation.', blank=True, default='')
    authors = models.CharField(max_length=512, help_text='Author list from the publication.', blank=True, default='')
    taxonomy = models.CharField(max_length=512, help_text='Comma-separated list of taxonomic lineage from GenBank', blank=True, default='')
    
    genbank_modification_date = models.DateField(help_text='Modification date of the GenBank entry, shown in the locus field of the flat file.', null=True, blank=True, default=None)

    collection_date = models.CharField(max_length=128, help_text='The date/time of or date/time range in which the specimen was collected.', blank=True, default='')
    collected_by = models.CharField(max_length=128, help_text='Persons or institute that collected the specimen.', blank=True, default='')
    identified_by = models.CharField(max_length=128, help_text='Name of expert who made the taxonomic identification of the species.', blank=True, default='')

    # Time when instance first created
    created = models.DateTimeField(auto_now_add=True, help_text='Date and time at which record was first created')
    # Date and time when data was last updated by this app
    updated = models.DateTimeField(auto_now=True, help_text='Date and time when sequence data was last updated from GenBank')

    # Only keep if CustomSequence queryset doesn't work
    class SequenceSource(models.TextChoices):
        GENBANK = 'GB', _('GenBank')
        IMPORT = 'IM', _('Import')

    # Source of this sequence
    data_source = models.CharField(max_length=2, choices=SequenceSource.choices, help_text='The source of this sequence.', blank=False, null=False, default=SequenceSource.GENBANK)

    def __str__(self) -> str:
        return f'{self.accession_number}, {str(self.organism)} ({str(self.id)})'

    def write_tree_identifier(self) -> str:
        '''Create the pipe-delimited identifier string used to represent this sequence in multiple sequence alignment or
        phylogenetic tree construction jobs.
        
        Raises ValueError if the components used to create the string already have pipes (i.e.: `|`)
        '''
        info = [self.version, 
        self.taxon_species.scientific_name if not self.taxon_species is None else '']
        for i in info:
            if '|' in i:
                raise ValueError('The SeqId or species name specified have existing pipe delimiters.')
        return '|'.join(info).replace(' ', '_')

    @staticmethod
    def extract_header_info(id):
        '''
        Parse the information stored within the identifier string `id`, based on the format
        used in `BlastQuerySequence.write_tree_identifier()`.
        
        HeaderInfo is a named tuple of three elements. The last element indicates whether
        the sequence was a query sequence. 
        '''
        header = id.split('|')
        id = header[0]
        species = header[1].replace('_', ' ')
        # Query sequences are differentiated from other seqs by a third value in the header
        is_query = len(header) > 2 and header[2] == 'query'
        return HeaderInfo(id, species, is_query)

    def sequence_length(self):
        '''
        Return the number of characters in the dna sequence.
        '''
        return len(self.dna_sequence)

    def reference(self):
        '''
        Return a string representation of the reference paper / source on GenBank, including
        the journal, author list, and paper title.
        '''
        return f'{self.journal} ({self.authors})'

    class Meta:
        ordering = ['accession_number']
        verbose_name = 'GenBank Accession'
        verbose_name_plural = 'GenBank Accessions'

class CustomSequence(NuccoreSequence):
    '''
    Used to represent a custom sequence that was added manually, instead of
    GenBank. Proxy model inheriting from NuccoreSequence
    '''
    class Meta:
        proxy = True
        ordering = ['accession_number']

class BlastRunManager(models.Manager):
    def listable(self, user: User):
        '''
        Return a queryset of BlastRun objects that are explicitly 
        viewable in a list for the given user. Primarily,
        this refers to BLASTN runs performed on databases the user has
        edit permission in.

        So, although blast runs are public, this will not necessary 
        return all blast runs within the database.
        '''
        libraries = Library.objects.editable(user)
        databases = BlastDb.objects.filter(library__in=libraries)
        return self.get_queryset().filter(db_used__in=databases)

class BlastRun(models.Model):
    objects = BlastRunManager()

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False, help_text='Unique identifier of the run')

    # Reference to the database used
    db_used = models.ForeignKey(BlastDb, related_name='usages', on_delete=models.CASCADE, help_text='The curated BLAST database against which the query BLAST was run.')
    # Job name
    job_name = models.CharField(max_length=255, blank=True, default='', help_text='A customizable name given by the submitter to refer to and describe this run.')

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

    def query_sequence_count(self):
        '''
        Return the number of query sequences returned by the result
        '''
        return self.queries.count()

    def __str__(self) -> str:
        return f'Run {self.id}'

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
        self.error_time = datetime.now()
        self.errors = ('\n' + debug_error_message) if len(self.errors) > 0 else debug_error_message
        self.status = self.JobStatus.ERRORED
        self.save()

    # What is the current status of the job?
    status = models.CharField(max_length=3,choices=JobStatus.choices, default=JobStatus.UNKNOWN, help_text='Current status of the job')
    # When was the run request first received (i.e. added to queue)
    received_time = models.DateTimeField(auto_now_add=True, help_text='Date and time when run first received by server')
    # Time that the job started running
    start_time = models.DateTimeField(blank=True, null=True, help_text='Date and time when job first started running')
    # Time that job successfully finished
    end_time = models.DateTimeField(blank=True, null=True, help_text='Date and time when job successfully finished running')
    # Time that job errored
    error_time = models.DateTimeField(blank=True, null=True, help_text='Date and time when job encountered an error')

    # Blast version
    blast_version = models.TextField(max_length=100, blank=True, default='', help_text='Version of BLASTn used')

    # Error for internal debugging
    errors = models.TextField(max_length=10000, blank=True, default='', help_text='Error message text')

    class Meta:
        ordering = ['-received_time']
        verbose_name = 'BLASTN Run'
        verbose_name_plural = 'BLASTN Runs'

HeaderInfo = namedtuple('HeaderInfo', 'id species is_query')

class BlastQuerySequence(models.Model):
    owner_run = models.ForeignKey(BlastRun, related_name='queries', on_delete=models.CASCADE, help_text='Job/run in which this query entry appeared')
    definition = models.CharField(max_length=255, help_text='Definition line')
    
    def query_seq_identifier(self) -> str:
        '''Return the sequence identifier from the definition line, which is used by the command-line BLAST+
        when identify sequences and report hits.'''
        return self.definition.split(' ')[0]

    def write_tree_identifier(self) -> str:
        '''Create the pipe-delimited identifier string used to represent this sequence in multiple sequence alignment or
        phylogenetic tree construction jobs.
        
        Raises ValueError if the components used to create the string already have pipes (i.e.: `|`)
        '''
        info = [self.query_seq_identifier(), self.original_species_name, 'query']
        for i in info:
            if '|' in i:
                raise ValueError('The SeqId or species name specified have existing pipe delimiters.')
        return '|'.join(info).replace(' ', '_')

    @staticmethod
    def extract_header_info(id):
        '''
        Parse the information stored within the identifier string `id`, based on the format
        used in `BlastQuerySequence.write_tree_identifier()`.
        
        HeaderInfo is a named tuple of three elements. The last element indicates whether
        the sequence was a query sequence. 
        '''
        header = id.split('|')
        id = header[0]
        species = header[1]
        # Query sequences are differentiated from other seqs by a third value in the header
        is_query = len(header) > 2 and header[2] == 'query'
        return HeaderInfo(id, species, is_query)

    query_sequence = models.CharField(max_length=50000, help_text='Sequence text')
    # The binomial name of the classification made by nucleotide blast. None if assignment not yet done, empty string if no assignment could be made.
    results_species_name = models.CharField(max_length=512, help_text='Binomial species identity assigned from the best BLAST hits. If there are multiple best hits, this is a comma-delimited list.', null=True, default=None)

    # Classification
    class QueryClassification(models.TextChoices):
        CORRECT_ID = 'Correct ID'
        INCORRECT_ID = 'Incorrect ID'
        NEW_ID = 'New ID'
        TENTATIVE_CORRECT_ID = 'Tentative Correct ID'
        INCORRECT_ID_NO_MATCH = 'Incorrect ID without match'
        UNKNOWN_ID = 'Unknown ID'
        TENATIVE_ADDITIONAL_SPECIES = 'Tentative additional species'
        NO_HITS = 'No hits'
    accuracy_category = models.CharField(choices=QueryClassification.choices, max_length=32, help_text='Category assigned after comparing sequence against reference libraries, using categories in Janzen et al. 2022.', null=True, default=None)
    # Original label
    original_species_name = models.CharField(max_length=255, help_text='Original binomial species identity reported from the user query', null=True, default=None)

    def best_hits(self) -> List[Any]:
        '''
        Return a list of best hits 
        '''
        return [h for h in self.hits.filter(best_hit=True)]
        
    def highest_percent_identity(self):
        '''
        Return the percent identity of the highest-performing hit
        '''
        best_hits = self.best_hits()
        percent_identities = [hit.percent_identity for hit in best_hits]
        return max(percent_identities) if len(percent_identities) > 0 else 0

    def evalue(self):
        '''
        Return the evalue of the best-performing hit
        '''
        best_hits = self.best_hits()
        return best_hits[0].evalue if len(best_hits) > 0 else 0

    class Meta:
        verbose_name = 'BLASTN Query Sequence'
        verbose_name_plural = 'BLASTN Query Sequences'
        # Prevent two query sequences in the same run from having the same definition
        unique_together = ('owner_run', 'definition')
        # Display all instances in order of definition
        ordering = ['definition'] 

    def __str__(self) -> str:
        return self.definition

class Hit(models.Model):
    db_entry = models.ForeignKey(NuccoreSequence, on_delete=models.CASCADE, help_text='BLAST database used in the run')
    query_sequence = models.ForeignKey(BlastQuerySequence, related_name='hits', on_delete=models.CASCADE, help_text='Query sequence to which this hit was registered.')

    def owner_run(self):
        return self.query_sequence.owner_run

    query_accession_version = models.CharField(max_length=128, help_text='Sequence identifier of query sequence')
    subject_accession_version = models.CharField(max_length=128, help_text='Sequence identifier of sequence in database')
    percent_identity = models.DecimalField(max_digits=6, decimal_places=3, help_text='Percent identity')
    alignment_length = models.PositiveIntegerField(help_text='Alignment length')
    mismatches = models.PositiveIntegerField(help_text='Number of mismatches')
    gap_opens = models.PositiveIntegerField(help_text='Number of Gap openings')
    query_start = models.PositiveIntegerField(help_text='Start of alignment in query')
    query_end = models.PositiveIntegerField(help_text='End of alignment in query')
    sequence_start = models.PositiveIntegerField(help_text='Start of alignment in subject')
    sequence_end = models.PositiveIntegerField(help_text='End of alignment in subject')
    evalue = models.DecimalField(max_digits=110, decimal_places=100, help_text='Expect value')
    bit_score = models.DecimalField(max_digits=110, decimal_places=100, help_text='Bit score')
    position = models.PositiveSmallIntegerField(help_text='Position of the hit within the results table. 1 = 1st (top) of the table.', validators=[MinValueValidator(1)])
    best_hit = models.BooleanField(help_text='Was this the best hit used for taxonomic assignment?', default=False)

    def first_hit(self):
        return self.position == 1

    class Meta:
        ordering = ['-percent_identity', 'db_entry__version']
        verbose_name = 'BLASTN Run Hit'
        verbose_name_plural = 'BLASTN Run Hits'
        unique_together = ['query_sequence', 'position']

class Annotation(models.Model):
    sequence = models.ForeignKey(NuccoreSequence, related_name='annotations', on_delete=models.CASCADE, help_text='Sequence to which this annotation was made.')
    poster = models.ForeignKey(User, on_delete=models.SET_NULL, null=True, help_text='User which created of the annotation.')
    timestamp = models.DateTimeField(help_text='Creation time of annotation', auto_now_add=True)
    class AnnotationType(models.TextChoices):
        MISIDENTIFICATION = 'Misidentification'
        TENTATIVE_NEW_SPECIES = 'Tentative New Species'
        UNRESOLVED_TAXONOMY = 'Unresolved Taxonomy'
        SEQUENCE_QUALITY = 'Sequence Issue'
        CONTAMINATION = 'Contamination'
        METADATA_QUALITY = 'Metadata Issue'
        VALIDATED_TAXONOMY = 'Validated Taxonomy'
        FLAG = 'Flag'
        KEYWORD = 'Keyword'
        QUESTIONS = 'Questions'
        OTHER = 'Other'
    annotation_type = models.CharField(max_length=32, choices=AnnotationType.choices, help_text='Categorization of content within the annotation.')
    comment = models.CharField(max_length=512, help_text='User-added comment text for annotation', blank=True, default='')

    def __str__(self) -> str:
        user_poster_display = '<deleted user>' \
                            if self.poster is None else self.poster.first_name
        return f'{Annotation.AnnotationType(self.annotation_type).label}, \
            {user_poster_display} ({self.timestamp})'
        

