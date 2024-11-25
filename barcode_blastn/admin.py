from django.db.models.functions import Length
from typing import Any, Dict, List, Optional, Tuple, Union
from django.contrib import admin
from django.db import models
from django import forms
from django.db.models.query import QuerySet
from django.db.models import Q
from django.template.loader import render_to_string
from simple_history.admin import SimpleHistoryAdmin
from django.utils.safestring import mark_safe

from django.http.request import HttpRequest
from django.urls import reverse

from django.utils.html import format_html
from django.contrib.auth.models import User
from barcode_blastn.admin_list_filters import BlastRunLibraryFilter, NuccoreSequencePublicationFilter
from barcode_blastn.controllers.blastdb_controller import add_sequences_to_database, create_blastdb, delete_blastdb, delete_library, filter_sequences_in_database, log_deleted_sequences, save_blastdb, save_custom_sequence, save_sequence, filter_sequences_in_database
from barcode_blastn.controllers.custom_sequence_controller import parse_file_upload_to_custom_sequence
from barcode_blastn.helper.parse_gb import GenBankConnectionError, InsufficientAccessionData, retrieve_gb
from barcode_blastn.models import Annotation, BlastDb, BlastQuerySequence, CustomSequence, DatabaseShare, Library, NuccoreSequence, BlastRun, Hit, TaxonomyNode
from django.forms import BaseInlineFormSet, Form, ModelForm, ValidationError
from barcode_blastn.permissions import DatabaseSharePermissions, HitSharePermission, LibrarySharePermissions, NuccoreSharePermission, RunSharePermissions
from barcode_blastn.renderers import get_letter

from barcode_blastn.serializers import LibraryEditSerializer, library_title, blast_db_title, run_title, nuccore_title, hit_title, custom_nuccore_title

class AdminAuthenticatedHttpRequest(HttpRequest):
    user: User

class NuccoreSequenceInline(admin.TabularInline):
    model = NuccoreSequence     
    show_change_link = False    # Show link to page to edit the sequence
    classes = ['collapse']      # Allow the entries to be collapsed
    extra = 0                   # Show one extra row by default
    fields = ['accession_number', 'version', 'id', 'organism', 'specimen_voucher', 'seq_len', 'annot_count', 'updated', 'link']
    ordering = ['version']

    @admin.display(description='Annotation count')
    def annot_count(self, obj):
        return obj.annot_count
 
    @admin.display(description='Length (bp)')
    def seq_len(self, obj):
        return obj.seq_len

    @admin.display(description='Link')
    def link(self, obj: Union[NuccoreSequence, CustomSequence]):
        if obj.data_source == NuccoreSequence.SequenceSource.GENBANK:
            link = 'admin:barcode_blastn_nuccoresequence_change'
        elif obj.data_source == NuccoreSequence.SequenceSource.IMPORT:
            link = 'admin:barcode_blastn_customsequence_change'
        else:
            raise NotImplementedError(f'Missing data page for {obj.data_source}')
        
        return mark_safe('<a href="%s">View</a>' % \
                            reverse(link,
                            args=(obj.id,)))

    def get_queryset(self, request: HttpRequest) -> QuerySet[NuccoreSequence]:
        # Get all sequences to display, and compute the length and
        # number of annotations per sequence.
        qs = super().get_queryset(request).annotate(seq_len=Length('dna_sequence'), annot_count=models.Count('annotations'))
        return qs

    def get_readonly_fields(self, request: HttpRequest, obj: Union[BlastDb, None] = None):
        base = list(set(
            [field.name for field in self.opts.local_fields] +  
            [field.name for field in self.opts.local_many_to_many]
        ))
        base.remove('accession_number')
        base.append('annot_count')
        base.append('seq_len')
        base.append('link')
        return base

    def has_module_permission(self, request: HttpRequest) -> bool:
        return NuccoreSharePermission.has_module_permission(request.user)

    def has_view_permission(self, request: HttpRequest, obj: Union[BlastDb, None]=None) -> bool:
        return DatabaseSharePermissions.has_view_permission(request.user, obj)
    
    def has_change_permission(self, request: HttpRequest, obj: Union[BlastDb, None] = None) -> bool:
        # Sequences cannot be modified directly; if an accession should be altered,
        # the old should be deleted and a new entry should be added
        return False

    def has_delete_permission(self, request, obj: Union[BlastDb, None] = None) -> bool:
        if obj is None:
            return DatabaseSharePermissions.has_change_permission(request.user, obj)
        else:
            return not obj.locked and DatabaseSharePermissions.has_change_permission(request.user, obj)

    def has_add_permission(self, request, obj: Union[BlastDb, None] = None) -> bool:
        return False

class UserPermissionsInline(admin.TabularInline):
    model = DatabaseShare
    extra = 0

    def get_fields(self, request, obj: Library):
        return ['grantee', 'permission_level', 'id']

    def get_readonly_fields(self, request, obj: Union[Library, None] = None):
        return ['id']

    def has_view_permission(self, request: HttpRequest, obj: Union[Library, None] = None) -> bool:
        return LibrarySharePermissions.has_view_permission(request.user, obj) 

    def has_add_permission(self, request, obj: Union[Library, None] = None) -> bool:
        '''
        Allow only the creator to add permissions.
        '''
        return LibrarySharePermissions.has_delete_permission(request.user, obj)

    def has_delete_permission(self, request, obj: Union[Library, None] = None) -> bool:
        '''
        Allow only the creator to delete permissions.
        '''
        return LibrarySharePermissions.has_delete_permission(request.user, obj)

    def has_change_permission(self, request, obj: Union[Library, None] = None) -> bool:
        '''
        Allow only the creator to change permissions.
        '''
        return LibrarySharePermissions.has_delete_permission(request.user, obj)
        
class VersionInline(admin.TabularInline):
    '''
    Show a version of a BlastDb inline in the Library admin page
    '''
    model = BlastDb
    fields = ['description', 'id', 'version_number', 'sequence_count', 'locked', 'created']
    readonly_fields = ['id', 'version_number', 'sequence_count', 'locked', 'created']
    show_change_link = True

    def sequence_count(self, obj: BlastDb) -> int:
        num_seqs: int = NuccoreSequence.objects.filter(owner_database=obj).count()
        return num_seqs

    def has_add_permission(self, request, obj: Union[Library, None] = None) -> bool:
        return False # Allow only the creator to add permissions.

    def has_delete_permission(self, request: HttpRequest, obj: Union[Library, None] = None) -> bool:
        return False # Disable deletions

    def has_change_permission(self, request: HttpRequest, obj: Union[Library, None] = None) -> bool:
        return False # Disable edits through the inline editor

# TODO: Make model form for creating and updating libraries
# https://stackoverflow.com/questions/24047308/django-rest-framework-serializers-and-django-forms

@admin.register(Library)
class LibraryAdmin(SimpleHistoryAdmin):
    '''
    Admin page for Library instances
    '''
    inlines = [UserPermissionsInline, VersionInline]
    list_display = ('id', 'marker_gene', 'custom_name', 'latest_version', 'owner', 'public', 'sequence_count')
    search_fields = ['description', 'id', 'custom_name', 'marker_gene']
    list_filter = ['public', 'owner__username', 'marker_gene']

    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {library_title} to view or change'}
        return super(LibraryAdmin, self).changelist_view(request, extra_context=extra_context)

    def get_fields(self, request: HttpRequest, obj: Optional[Library] = None):
        fields = LibraryEditSerializer.Meta.fields.copy()
        fields.extend(['owner', 'created'])
        return fields

    def has_module_permission(self, request: HttpRequest) -> bool:
        return LibrarySharePermissions.has_module_permission(request.user)
    
    def has_add_permission(self, request: HttpRequest) -> bool:
        return LibrarySharePermissions.has_add_permission(request.user, None)

    def has_change_permission(self, request: HttpRequest, obj: Optional[Library] = None) -> bool:
        return LibrarySharePermissions.has_change_permission(request.user, obj)

    def has_delete_permission(self, request: HttpRequest, obj: Optional[Library] = None) -> bool:
        return LibrarySharePermissions.has_delete_permission(request.user, obj)

    def has_view_permission(self, request: HttpRequest, obj: Optional[Library] = None) -> bool:
        return LibrarySharePermissions.has_view_permission(request.user, obj)
         
    def get_readonly_fields(self, request: HttpRequest, obj: Optional[Library]) :
        return ['owner', 'created']
    
    def sequence_count(self, obj: Library):
        blastdb: Union[BlastDb, None] = BlastDb.objects.latest(obj)
        if blastdb is None:
            return 0
        return NuccoreSequence.objects.filter(owner_database=blastdb).count()

    def delete_model(self, request: HttpRequest, obj: Library) -> None:
        delete_library(obj)

    def save_model(self, request: Any, obj: Any, form: Any, change: Any) -> None:
        obj.owner = request.user
        return super().save_model(request, obj, form, change)

    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        if isinstance(request.user, User):
            return Library.objects.viewable(request.user)
        else:
            return Library.objects.none()

class BlastDbForm(ModelForm):
    
    min_length = forms.IntegerField(min_value=-1, max_value=10000, initial=-1, help_text='Filter out sequences with character count less than the minimum sequence length. Set to -1 to ignore this option.', required=False)
    max_length = forms.IntegerField(min_value=-1, max_value=10000, initial=-1, help_text='Filter out sequences with character count greater than the minimum sequence length. Set to -1 to ignore this option.', required=False)
    max_ambiguous_bases = forms.IntegerField(min_value=-1, max_value=10000, initial=-1, help_text='Filter out sequences with a greater number of ambiguous characters (Ns) than the maximum. Set to -1 to ignore this option.', required=False)

    # Allow filtering by blacklisting accessions
    blacklist = forms.CharField(
        widget=forms.Textarea,
        min_length=0,
        max_length=10000,
        help_text='Filter out the given accession numbers or versions, by preventing these from being added and removing records already added. Leave empty to ignore this option.',
        required=False
    )

    # If True, filter if taxonomy missing
    require_taxonomy = forms.BooleanField(help_text='Filter out records without taxonomic information between the superkingdom to species levels (inclusive). Uncheck to ignore this option.', required=False, initial=False)
   
   
    # Methods for adding assessions from GenBank
    accession_list_upload = forms.FileField(help_text='Add large numbers of sequences at once by uploading a .txt file, with each accession number on a separate line.', required=False)
    accession_list_text = forms.CharField(widget=forms.Textarea, help_text='Add multiple sequences by pasting in a list of accession numbers, one per line.', required=False)
    search_term = forms.CharField(widget=forms.Textarea, help_text='Add sequences by GenBank search terms.', required=False)
    base_database = forms.ModelChoiceField(queryset=BlastDb.objects.none(), help_text='Inherit all accession numbers from an existing BLAST database in order to populate the current.', required=False)

    # Methods for importing sequences
    dna_fasta_upload = forms.FileField(help_text='Import custom DNA sequences from FASTA file.', required=False)


    def __init__(self, *args, **kwargs) -> None:
        super(BlastDbForm, self).__init__(*args, **kwargs)
        user = getattr(self, 'user', None)
        if user:
            self.fields['base_database'].queryset = BlastDb.objects.viewable(user)
        else:
            raise ValueError('No user')

    def clean(self) -> Dict[str, Any]:
        cleaned_data = super().clean()
        existing = NuccoreSequence.objects.filter(owner_database=self.instance)

        # Check that the accession numbers or accession.versions in the list_text are not duplicates
        accession_file_upload = cleaned_data.get('accession_list_upload', None)
        if not accession_file_upload is None:
            accession_file_upload.seek(0)
            file_numbers = [line.decode('utf-8') for line in accession_file_upload.readlines()]
            # remove whitespace and newline characters
            file_numbers = [t.strip() for t in file_numbers]
            file_numbers = [t for t in file_numbers if len(t) > 0]
            # query the database for any entries with matching accession.versions or accession numbers
            exists_duplicate = existing.filter(Q(accession_number__in=file_numbers) | Q(version__in=file_numbers))
            if exists_duplicate.exists():
                raise ValidationError({'accession_list_upload': 'The file given for "accession list upload" contains accession numbers or accession.versions already in the database.'})

        # Check that the accession numbers or accession.versions in the list_upload are not duplicates
        accession_list_text: str = cleaned_data.get('accession_list_text', '')
        if len(accession_list_text) > 0:
            # remove whitespace and newline characters
            list_text = accession_list_text.replace('\r\n', '\n').split('\n')
            list_text = [t.strip() for t in list_text]
            list_text = [t for t in list_text if len(t) > 0]
            # query the database for any entries with matching accession.versions or accession numbers
            exists_duplicate = existing.filter(Q(accession_number__in=list_text) | Q(version__in=list_text))
            if exists_duplicate.exists():
                raise ValidationError({'accession_list_text': 'The raw text list given for "accession list text" contains accession numbers or accession.versions already in the database.'})
        
        return cleaned_data

@admin.register(BlastDb)
class BlastDbAdmin(SimpleHistoryAdmin):
    '''
    Admin page for BlastDb instances.
    '''
    inlines = [NuccoreSequenceInline]
    form = BlastDbForm
    list_display = ('custom_name', 'library', 'version_number', 'sequences_admin_count', 'id', 'locked')
    search_fields = ['custom_name', 'id', 'library__custom_name', 'library__owner__username', 'library__description', 'description']
    list_filter = ['library__custom_name', 'genbank_version', 'locked']
    history_list_display = ['added', 'filter_options', 'deleted', 'search_terms', 'locked', 'changed_fields']

    def added(self, obj: BlastDb):
        '''Return a column with text wrapping to show added sequence identifiers in history'''
        return format_html("<div style='width: 100px; word-wrap: break-word'>{text}</div>", text=getattr(obj, 'added', ''))

    def deleted(self, obj: BlastDb):
        '''Return a column with text wrapping to show deleted sequence identifiers in history'''
        return format_html("<div style='width: 100px; word-wrap: break-word'>{text}</div>", text=getattr(obj, 'deleted', ''))

    def search_terms(self, obj: BlastDb):
        '''Return a column with text wrapping to show search terms in history'''
        return format_html("<div style='width: 100px; word-wrap: break-word'>{text}</div>", text=getattr(obj, 'search_terms', ''))

    def changed_fields(self, obj):
        if obj.prev_record:
            delta = obj.diff_against(obj.prev_record)
            return delta.changed_fields
        return None

    def library_is_public(self, obj: BlastDb):
        return not obj.library is None and obj.library.public

    def library_owner(self, obj: BlastDb):
        if not obj.library is None:
            return obj.library.owner
        else:
            return None

    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {blast_db_title} to view or change'}
        return super(BlastDbAdmin, self).changelist_view(request, extra_context=extra_context)

    def delete_model(self, request: HttpRequest, obj: BlastDb) -> None:
        '''Delete instance of database'''
        delete_blastdb(obj)

    def get_form(self, request: Any, obj: Optional[BlastDb] = ..., change: bool = ..., **kwargs: Any) -> Any:
        form = super().get_form(request, obj, change, **kwargs)
        form.user = request.user
        return form

    def get_fieldsets(self, request: AdminAuthenticatedHttpRequest, obj: Optional[BlastDb] = None) -> List[Tuple[str, Dict[str, Any]]]:
        f: List[Tuple[str, Dict[str, Any]]]
        f = [('Details', { 'fields': ['library', 'library_owner', 'custom_name', 'description', 'version_number']})]
        # Only allow filter options in creation mode, or in edit mode if database is not locked
        if obj is None or not obj.locked:
            f.append(('Filter', { 'fields': ['min_length', 'max_length', 'max_ambiguous_bases', 'blacklist', 'require_taxonomy']}))
        if obj is None:
            f.extend([
                ('Visibility', { 'fields': ['locked', 'library_is_public']}),
                ('Download GenBank accessions', { 'fields': ('base_database', 'accession_list_upload', 'accession_list_text', 'search_term')}),
                ('Import custom sequences', { 'fields': ('dna_fasta_upload',)})
            ])
        elif not obj.locked and isinstance(request.user, User) and DatabaseSharePermissions.has_change_permission(request.user, obj=obj):
            f.extend([
                ('Download additional GenBank accessions', { 'fields': ('accession_list_upload', 'accession_list_text', 'search_term')}),
                ('Import additional custom sequences', { 'fields': ('dna_fasta_upload',)}),
            ])
        return f

    def save_formset(self, request: AdminAuthenticatedHttpRequest, form: Any, formset: BaseInlineFormSet, change: Any) -> None:
        super_results = formset.save()

        # Save right now, instead of in save_model()
        obj = form.save(commit=False)

        # Has the user indicated for the database to be locked?
        will_lock = form.cleaned_data.get('locked', False)

        # If sequences are deleted, add it to change history
        if len(formset.deleted_objects) > 0:
            log_deleted_sequences(formset.deleted_objects, obj)
       
        if will_lock:
            # if user wants to lock the database, perform locking after saving the formset
            obj.locked = True 
        if will_lock or len(formset.deleted_objects) > 0:
            # Save to lock database and/or log sequence deletions
            save_blastdb(obj, request.user, perform_lock=will_lock)
        
        return super_results

    def save_model(self, request: AdminAuthenticatedHttpRequest, obj: BlastDb, form: Any, change: Any) -> None:
        # Since the model is saved first before any changes in sequence, don't lock the
        # database yet. The database will be locked if needed by save_formset
        obj.locked = False

        accessions = []

        base = form.cleaned_data.get('base_database', None)

        # Add accessions from the file upload
        accession_file_upload = form.cleaned_data.get('accession_list_upload')
        if accession_file_upload:
            accession_file_upload.seek(0)
            file_numbers = [line.decode('utf-8') for line in accession_file_upload.readlines()]
            file_numbers = [t.strip() for t in file_numbers]
            file_numbers = [t for t in file_numbers if len(t) > 0]
            accessions.extend(file_numbers)

        # Add other accessions specified manually
        accession_list_text = form.cleaned_data.get('accession_list_text')
        if len(accession_list_text) > 0:
            list_text = accession_list_text.replace('\r\n', '\n').split('\n')
            list_text = [t.strip() for t in list_text]
            list_text = [t for t in list_text if len(t) > 0]
            accessions.extend(list_text)
        
        blastdb_fields = {
            'custom_name': form.cleaned_data.get('custom_name'),
            'description': form.cleaned_data.get('description'),
        }
        filter_fields = ['min_length', 'max_length', 'max_ambiguous_bases', 'blacklist', 'require_taxonomy']
        filter_args: dict[str, Any] = {
            field: form.cleaned_data.get(field) for field in filter_fields if field in form.cleaned_data
        }
        raw_blacklist = filter_args.get('blacklist', '').strip()
        filter_args['blacklist'] = raw_blacklist.split('\n') if len(raw_blacklist) > 0 else []
        
        search_term = form.cleaned_data.get('search_term', None)
        print('Submitting save request from admin with ', change, search_term, accessions)
        if not change:
            create_blastdb(additional_accessions=accessions, user=request.user, base=base, database=obj, search_term=search_term, **blastdb_fields, library=obj.library, **filter_args)
        else:
            obj._change_reason = 'Save changes'
            if len(accessions) > 0 or (not search_term is None and len(search_term) > 0):
                add_sequences_to_database(obj, user=request.user, desired_numbers=accessions, search_term=search_term, **filter_args)     
            else:
                filter_sequences_in_database(obj, user=request.user, **filter_args)  

        dna_fasta_upload = form.cleaned_data.get('dna_fasta_upload')
        if dna_fasta_upload:
            seqs = parse_file_upload_to_custom_sequence(dna_fasta_upload, obj)
            save_custom_sequence(seqs, request.user, True)

    def get_readonly_fields(self, request, obj: Union[BlastDb, None] = None):
        if obj is None: # is adding
            return ['version_number', 'library_is_public', 'library_owner']
        else:
            return ['library', 'version_number', 'library_is_public', 'library_owner']

    @admin.display(
        ordering="sequences__count",
        description="Sequences",
    )
    def sequence_count(self, obj: BlastDb):
        return obj.sequence_count()

    def formfield_for_foreignkey(self, db_field: models.ForeignKey[BlastDb], request: HttpRequest, **kwargs):
        if request is None:
            return super().formfield_for_foreignkey(db_field, request, **kwargs)
        if db_field.name == 'library':
            if isinstance(request.user, User):
                # restrict choices of library to those that the user can edit
                kwargs['queryset'] = Library.objects.editable(request.user)
            else:
                kwargs['queryset'] = Library.objects.none()
        return super().formfield_for_foreignkey(db_field, request, **kwargs)
        
    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        qs : QuerySet[BlastDb]
        if not isinstance(request.user, User) or not request.user.is_authenticated:
            qs = BlastDb.objects.none()
        else:
            qs = BlastDb.objects.viewable(request.user)
            # Introduce "seqences__count" by counting the number of related sequences
        qs = qs.annotate(sequences__count=models.Count('sequences'))
        return qs

    @admin.display(ordering="sequences__count", description='Sequences')
    def sequences_admin_count(self, obj: BlastDb):
        # Return value computed in self.get_queryset()
        return obj.sequences__count

    def has_module_permission(self, request: HttpRequest) -> bool:
        return DatabaseSharePermissions.has_module_permission(request.user)

    def has_view_permission(self, request: HttpRequest, obj: Union[BlastDb, None] = None) -> bool:
        return DatabaseSharePermissions.has_view_permission(request.user, obj)

    def has_add_permission(self, request: HttpRequest) -> bool:
        return DatabaseSharePermissions.has_add_permission(request.user, None)

    def has_delete_permission(self, request: HttpRequest, obj: Union[BlastDb, None] = None) -> bool:
        return DatabaseSharePermissions.has_delete_permission(request.user, obj)

    def has_change_permission(self, request: HttpRequest, obj: Union[BlastDb, None] = None) -> bool:
        return DatabaseSharePermissions.has_change_permission(request.user, obj)        

    def changeform_view(self, request: HttpRequest, object_id: Union[str, None], form_url: str, extra_context: Optional[Dict[str, bool]]) -> Any:
        extra_context = extra_context or {}
        extra_context['show_save_and_continue'] = True # show Save and Continue button
        extra_context['show_save_and_add_another'] = True # hide Save and Add Another button
        extra_context['show_save'] = True # hide save button
        extra_context['show_delete'] = True # show delete button
        return super().changeform_view(request, object_id, form_url, extra_context)

class NuccoreAdminModifyForm(ModelForm):
    '''
    Form for changing and editing accessions.
    '''

    def clean(self):
        '''
        Validate new data in the form before saving it.

        If accession number does not correspond to entry, raise a ValidationError.

        Return nothing if no validation errors.
        '''

        cleaned_data = super().clean()
        instance: NuccoreSequence = self.instance

        # Only check data if we are adding a new entry
        if instance._state.adding:
            accession_number: Any = cleaned_data.get('accession_number', None)
            # check if the form data contains an accession_number
            if accession_number is None:
                raise ValidationError(f'Could not locate an accession number in the submitted form data.')

            owner_database: Union[BlastDb, None] = cleaned_data.get('owner_database', None)

            # check that the sequence will be assigned to a database
            if owner_database is None:
                raise ValidationError(f'No owner database specified. Sequence must be associated with a database when being created. If no databases are selectable, contact a superuser to give you edit permission to a database.')
            else:
                if owner_database.locked:
                    raise ValidationError(f'The database {owner_database} ({owner_database.id}) is locked, so changes cannot be made to its sequences and new sequences cannot be added to it.')

            # check if the accession already exists in the database
            duplicate_exists: bool = False
            if self.instance.pk is None: # if this is an addition
                duplicate_exists = NuccoreSequence.objects.filter(accession_number = accession_number, owner_database=owner_database).exists()
            
            if duplicate_exists:
                raise ValidationError(f'Error: Sequence entry for accession number {accession_number} already exists in the same database.')
        
            try:
                retrieve_gb(accession_numbers=[accession_number], raise_if_missing=True)
            except GenBankConnectionError as err:
                raise ValidationError('Encountered connection error while requesting data from GenBank. Please try again.') from err  
            except InsufficientAccessionData as err:
                raise ValidationError(f'Some number of following the accession numbers and search terms do not match any record. Accessions: {", ".join(err.missing_accessions)}. Terms: {err.term}.')
            
        return None

class AnnotationInline(admin.TabularInline):
    '''
    Inline for showing annotation under a sequence
    '''
    model = Annotation
    fields = ['user', 'timestamp', 'annotation_type', 'comment']
    readonly_fields = ['user', 'timestamp']
    extra = 0

    def user(self, obj: Optional[Annotation]):
        if obj is None:
            return '<unspecified_user>'
        elif obj.poster is None:
            return '<deleted_user>'
        else:
            return obj.poster.username

    def has_module_permission(self, request: HttpRequest) -> bool:
        return NuccoreSharePermission.has_module_permission(request.user)
    
    def has_view_permission(self, request: HttpRequest, obj: Union[NuccoreSequence, None]=None) -> bool:
        return NuccoreSharePermission.has_view_permission(request.user, obj)
    
    def has_change_permission(self, request: HttpRequest, obj: Union[NuccoreSequence, None] = None) -> bool:
        # Annotations cannot be modified after initial creation
        return False

    def has_delete_permission(self, request, obj: Union[NuccoreSequence, None] = None) -> bool:
        return NuccoreSharePermission.has_delete_permission(request.user, obj)

    def has_add_permission(self, request, obj: Union[NuccoreSequence, None] = None) -> bool:
        return request.user.is_authenticated

@admin.register(CustomSequence)
class CustomSequenceAdmin(admin.ModelAdmin):
    # Queryset will only include custom sequence
    source: NuccoreSequence.SequenceSource = NuccoreSequence.SequenceSource.IMPORT

    list_display = ('accession_number', 'version', 'organism', 'specimen_voucher', 'id', 'owner_database_link', 'seq_length')
    fields_excluded = ['uuid', 'genbank_modification_date']
    inlines = [AnnotationInline]
    list_filter = ['owner_database', NuccoreSequencePublicationFilter]

    def get_search_fields(self, request: HttpRequest) -> List[str]:
        search_fields = ['version', 'organism', 'specimen_voucher', 'id', 'taxonomy', 'owner_database__custom_name']
        # Allow search of taxon scientific names
        for taxon_level in ['species', 'genus', 'family', 'order', 'class', 'kingdom', 'superkingdom']:
            search_fields.extend([f'taxon_{taxon_level}__scientific_name', f'taxon_{taxon_level}__id'])
        return search_fields

    @admin.display(ordering='dna_sequence__len', description='Length (bp)')
    def seq_length(self, obj: NuccoreSequence):
        # Return value computed in self.get_queryset()
        return obj.dna_sequence__len

    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {custom_nuccore_title} to view or change'}
        return super(CustomSequenceAdmin, self).changelist_view(request, extra_context=extra_context)

    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        qs: QuerySet[NuccoreSequence]
        # Restrict the sequences visible based on which databases they belong to, and which libraries the user can see
        if not isinstance(request.user, User) or not request.user.is_authenticated:
            qs = CustomSequence.objects.none()
        else:
            qs = CustomSequence.objects.viewable(request.user).filter(data_source=self.source)
        # Compute the lengths of all sequences
        qs = qs.annotate(dna_sequence__len=Length('dna_sequence'))
        return qs

    def formfield_for_foreignkey(self, db_field, request: Optional[HttpRequest], **kwargs):
        if request is None:
            return super().formfield_for_foreignkey(db_field, request, **kwargs)
        if db_field.name == 'owner_database':
            if not request.user is None:
                if isinstance(request.user, User):
                    # Set the available choices of database the accession can belong to
                    # based on the user's permission
                    kwargs['queryset'] = BlastDb.objects.editable(request.user).filter(locked=False)
        elif db_field.name in ['taxon_species', 
            'taxon_genus', 'taxon_family', 'taxon_order', 'taxon_class', 
            'taxon_phylum', 'taxon_kingdom', 'taxon_superkingdom']:
            kwargs['queryset'] = TaxonomyNode.objects.filter(rank=get_letter(db_field.name))
        return super().formfield_for_foreignkey(db_field, request, **kwargs) 

    def get_fields(self, request, obj=None):
        fields = [
            'owner_database', 'organism', 'version', 'definition', 
            'organelle', 'accession_number', 'specimen_voucher',
            'id', 'isolate', 'country', 'dna_sequence', 'collected_by', 
            'collection_date',  'identified_by', 'lat_lon', 'type_material',
            'created', 'updated', 'taxonomy', 'taxon_species', 
            'taxon_genus', 'taxon_family', 'taxon_order', 'taxon_class', 
            'taxon_phylum', 'taxon_kingdom', 'taxon_superkingdom', 'title', 
            'authors', 'journal'
        ]
        return fields

    def get_readonly_fields(self, request, obj: Optional[NuccoreSequence]=None):
        fields = [ 'id', 'created', 'updated', 'data_source']
        if not obj is None:
            fields.extend(['owner_database', 'accession_number'])
        return fields

    def get_fieldsets(self, request: HttpRequest, obj: Optional[NuccoreSequence] = None) -> List[Tuple[Optional[str], Dict[str, Any]]]:
        summary_fields = ['id', 'accession_number', 'version', 'definition', 'taxonomy', 'owner_database']
        # Omit data source from add form, but keep it readonly in change
        if obj is not None:
            summary_fields.append('data_source')
        return [
            ('Summary', { 'fields': summary_fields }), 
            ('Source Information', {'fields': ['specimen_voucher', 'type_material', 'organelle', 'isolate', 'country', 'collected_by', 'collection_date', 'lat_lon', 'identified_by']}),
            ('History', { 'fields': ['created', 'updated']}),
            ('Taxonomy', { 'fields': ['taxonomy', 'taxon_species', 'taxon_genus', 'taxon_family', 'taxon_order', 'taxon_class', 'taxon_phylum', 'taxon_kingdom', 'taxon_superkingdom']}),
            ('Reference', { 'fields': ['title', 'authors', 'journal']}),
            ('Sequence', { 'fields': ['dna_sequence'] })
        ]
     
    @admin.display(description='Database')
    def owner_database_link(self, obj):
        url = reverse("admin:%s_%s_change" % ('barcode_blastn', 'blastdb'), args=(obj.owner_database.id,))
        return format_html("<a href='{url}'>{obj}</a>", url=url, obj=obj.owner_database)

    def has_module_permission(self, request: HttpRequest) -> bool:
        return NuccoreSharePermission.has_module_permission(request.user)

    def has_view_permission(self, request: HttpRequest, obj: Union[NuccoreSequence, None]=None) -> bool:
        return NuccoreSharePermission.has_view_permission(request.user, obj)

    def has_change_permission(self, request, obj: Union[NuccoreSequence, None]=None):
        return NuccoreSharePermission.has_change_permission(request.user, obj)

    def has_delete_permission(self, request, obj: Union[NuccoreSequence, None]=None):
        return NuccoreSharePermission.has_delete_permission(request.user, obj)

    def has_add_permission(self, request, obj: Union[NuccoreSequence, None]=None):
        return NuccoreSharePermission.has_add_permission(request.user, obj)
    

    def delete_model(self, request: HttpRequest, obj: NuccoreSequence) -> None:
        super().delete_model(request, obj)
        log_deleted_sequences([obj], obj.owner_database)
        obj.owner_database.save()

    def save_formset(self, request: Any, form: Any, formset: Any, change: Any) -> None:
        instances = formset.save(commit=False)
        for instance in instances:
            # Do something with `instance`
            if instance._state.adding:
                instance.poster = request.user
            instance.save() 
        formset.save_m2m()
        return super().save_formset(request, form, formset, change)

    def save_model(self, request: Any, obj: Any, form: Any, change: bool) -> None:
        save_custom_sequence([obj], user=request.user)

    def changeform_view(self, request: HttpRequest, object_id: Union[str, None], form_url: str, extra_context: Optional[Dict[str, bool]]) -> Any:
        extra_context = extra_context or {}

        extra_context['show_save_and_continue'] = True # show Save and Continue button
        extra_context['show_save_and_add_another'] = object_id is None # show Save and Add Another button
        extra_context['show_save'] = object_id is None # hide save button
        extra_context['show_delete'] = object_id is not None # show delete button
        return super().changeform_view(request, object_id, form_url, extra_context)

@admin.register(NuccoreSequence)
class NuccoreAdmin(CustomSequenceAdmin):
    '''
    Admin page for showing accession instances.
    '''
    source = NuccoreSequence.SequenceSource.GENBANK

    form = NuccoreAdminModifyForm

    fields_excluded = ['uuid']


    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {nuccore_title} to view or change'}
        return super(NuccoreAdmin, self).changelist_view(request, extra_context=extra_context)

    def get_readonly_fields(self, request, obj: Optional[NuccoreSequence]=None):
        fields = [
            'organism', 'version', 'definition', 'organelle',
            'specimen_voucher', 'data_source', 'id', 'isolate', 'country',
            'dna_sequence', 'collected_by', 'collection_date', 
            'identified_by', 'lat_lon', 'type_material',
            'created', 'updated', 'taxonomy', 'taxon_species', 
            'taxon_genus', 'taxon_family', 'taxon_order', 'taxon_class', 
            'taxon_phylum', 'taxon_kingdom', 'taxon_superkingdom','title', 
            'authors', 'journal', 
        ]
        if not obj is None:
            fields.extend(['owner_database', 'accession_number'])
        fields.extend(['taxon'])
        return fields

    def get_fieldsets(self, request: HttpRequest, obj: Optional[NuccoreSequence] = None) -> List[Tuple[Optional[str], Dict[str, Any]]]:
        is_adding = obj is None
        if is_adding:
            return [
                ('Summary', { 'fields': ['id', 'accession_number', 'owner_database']})
            ]
        else:
            return [
                ('Summary', { 'fields': ['id', 'accession_number', 'version', 'definition', 'taxonomy', 'owner_database', 'data_source']}), 
                ('Source Information', {'fields': ['specimen_voucher', 'type_material', 'organelle', 'isolate', 'country', 'collected_by', 'collection_date', 'lat_lon', 'identified_by']}),
                ('History', { 'fields': ['created', 'updated']}),
                ('Taxonomy', { 'fields': ['taxonomy', 'taxon_species', 'taxon_genus', 'taxon_family', 'taxon_order', 'taxon_class', 'taxon_phylum', 'taxon_kingdom', 'taxon_superkingdom']}),
                ('Reference', { 'fields': ['title', 'authors', 'journal']}),
                ('Sequence', { 'fields': ['dna_sequence'] })
            ]

    def save_model(self, request: Any, obj: Any, form: Any, change: bool) -> None:
        if not change:
            # Because the fields of sequences cannot be modified,
            # only call save if this is a newly added sequence
            save_sequence(obj, user=request.user, commit=True)
     
class BlastQuerySequenceInline(admin.StackedInline):
    '''
    Row to show query sequences under a BLAST run.
    '''
    model = BlastQuerySequence
    extra = 0
    show_change_link = True

    def has_change_permission(self, request, obj=None):
        return False
    def has_add_permission(self, request, obj=None) -> bool:
        return False
    def has_delete_permission(self, request, obj=None) -> bool:
        return False

class HitInline(admin.TabularInline):
    '''
    Row to show a database hit under a BLAST run.
    '''
    model = Hit
    fields = ['db_entry', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue_value', 'bit_score_value']
    readonly_fields = ['evalue_value', 'bit_score_value']
    extra = 0
    show_change_link = True

    def evalue_value(self, obj):
        return '%.4E' % obj.evalue
    def bit_score_value(self, obj):
        return '%.4E' % obj.bit_score
    def has_change_permission(self, request, obj=None):
        return False
    def has_add_permission(self, request, obj=None) -> bool:
        return False
    def has_delete_permission(self, request, obj=None) -> bool:
        return False

@admin.register(BlastQuerySequence)
class BlastQuerySequenceAdmin(admin.ModelAdmin):
    inlines = [HitInline]
    fields = ['definition', 'query_sequence', 'original_species_name', 'results_species_name', 'accuracy_category']
    list_display = ['definition', 'length', 'number_of_hits', 'owner_run']
    list_display_links = ['definition']
    search_fields = ['definition']

    def number_of_hits(self, obj: Optional[BlastQuerySequence]):
        if not obj is None:
            return Hit.objects.filter(query_sequence=obj).count()
        else:
            return 0

    def length(self, obj: Optional[BlastQuerySequence]):
        if not obj is None:
            return len(obj.query_sequence)
        else:
            return 0

    def has_module_permission(self, request: HttpRequest) -> bool:
        return RunSharePermissions.has_module_permission(request.user)

    def has_view_permission(self, request: HttpRequest, obj: Union[BlastQuerySequence, None]=None) -> bool:
        return RunSharePermissions.has_view_permission(request.user, obj.owner_run if not obj is None else None)

    def has_change_permission(self, request, obj: Union[BlastQuerySequence, None]=None):
        return RunSharePermissions.has_change_permission(request.user, obj.owner_run if not obj is None else None)

    def has_delete_permission(self, request, obj: Union[BlastQuerySequence, None]=None):
        return RunSharePermissions.has_delete_permission(request.user, obj.owner_run if not obj is None else None)

    def has_add_permission(self, request, obj: Union[BlastQuerySequence, None]=None):
        return RunSharePermissions.has_add_permission(request.user, obj.owner_run if not obj is None else None)

@admin.register(BlastRun)
class BlastRunAdmin(admin.ModelAdmin):
    '''
    Admin page for showing run information.
    '''
    inlines = [BlastQuerySequenceInline]
    show_change_link = True
    list_display = ('id', 'query_sequences', 'job_name', 'received_time', 'status', 'db_used_link', 'create_hit_tree', 'create_db_tree')
    search_fields = ['id', 'job_name', 'db_used__custom_name', 'alignment_job_id', 'complete_alignment_job_id', 'errors']
    list_filter = ['status', 'blast_version', 'create_hit_tree', 'create_db_tree', BlastRunLibraryFilter]

    @admin.display(ordering='queries__count', description='Queries')
    def query_sequences(self, obj: BlastDb):
        # Return value computed in self.get_queryset()
        return obj.queries__count

    @admin.display(ordering='db_used__custom_name', description='Database used')
    def db_used_link(self, obj: Optional[BlastRun]):
        if obj is None:
            return ''
        else:
            url = reverse("admin:%s_%s_change" % ('barcode_blastn', 'blastdb'), args=(obj.db_used.id,))
            return format_html("<a href='{url}'>{obj}</a>", url=url, obj=obj.db_used)

    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {run_title} to view or change'}
        return super(BlastRunAdmin, self).changelist_view(request, extra_context=extra_context)
    
    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        qs : QuerySet[BlastDb]
        if not request.user.is_authenticated:
            qs = BlastRun.objects.none()
        elif isinstance(request.user, User):
            qs = BlastRun.objects.listable(request.user)
        else:
            qs = BlastRun.objects.none()    
        qs = qs.annotate(queries__count=models.Count('queries'))
        return qs
                
    def has_module_permission(self, request: HttpRequest) -> bool:
        # allow module access if they can access databases
        return RunSharePermissions.has_module_permission(request.user)

    def has_view_permission(self, request: HttpRequest, obj: Union[BlastRun, None]=None) -> bool:
        return RunSharePermissions.has_view_permission(request.user, obj)

    def has_change_permission(self, request, obj: Union[BlastRun, None]=None):
        return RunSharePermissions.has_change_permission(request.user, obj) 

    def has_delete_permission(self, request, obj: Union[BlastRun, None]=None):
        return RunSharePermissions.has_delete_permission(request.user, obj) 

    def has_add_permission(self, request):
        return RunSharePermissions.has_add_permission(request.user, None) 

    def get_readonly_fields(self, request, obj=None):
        return list(set(
            [field.name for field in self.opts.local_fields] +
            [field.name for field in self.opts.local_many_to_many]
        ))

@admin.register(Hit)
class HitAdmin(admin.ModelAdmin):
    '''
    Admin page for showing Hit instances.
    '''
    inlines = []
    show_change_link = True
    list_display = ['id', 'query_link', 'db_entry_link', 'owner_run_link']
    search_fields = ('id', 'db_entry__version', 'db_entry__organism', 'query_sequence__definition')

    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {hit_title} to view or change'}
        return super(HitAdmin, self).changelist_view(request, extra_context=extra_context)

    @admin.display(description='Query')
    def query_link(self, obj):
        url = reverse("admin:%s_%s_change" % ('barcode_blastn', 'blastquerysequence'), args=(obj.query_sequence.id,))
        return format_html("<a href='{url}'>{obj}</a>", url=url, obj=obj.query_sequence) 

    @admin.display(description='Reference')
    def db_entry_link(self, obj):
        url = reverse("admin:%s_%s_change" % ('barcode_blastn', 'nuccoresequence'), args=(obj.db_entry.id,))
        source = obj.db_entry.taxon_species.scientific_name if not obj.db_entry.taxon_species is None else obj.db_entry.organism
        return format_html("<a href='{url}'>{obj}</a>", url=url, obj=f'{obj.db_entry.version}, {source}') 

    @admin.display(description='Run')
    def owner_run_link(self, obj):
        url = reverse("admin:%s_%s_change" % ('barcode_blastn', 'blastrun'), args=(obj.owner_run().id,))
        return format_html("<a href='{url}'>{obj}</a>", url=url, obj=obj.owner_run())

    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        if not request.user.is_authenticated:
            return Hit.objects.none()
        elif isinstance(request.user, User):
            return Hit.objects.filter(query_sequence__owner_run__in=BlastRun.objects.listable(request.user))
        else:
            return Hit.objects.none()

    def has_module_permission(self, request: HttpRequest) -> bool:
        # allow module access if they can access databases
        return HitSharePermission.has_module_permission(request.user)

    def has_view_permission(self, request: HttpRequest, obj: Union[Hit, None]=None) -> bool:
        return HitSharePermission.has_view_permission(request.user, obj)

    def has_change_permission(self, request, obj: Union[Hit, None]=None):
        # prohibit runs from being modified through admin panel
        return HitSharePermission.has_change_permission(request.user, obj) 

    def has_delete_permission(self, request, obj: Union[Hit, None]=None):
        # runs can only be deleted if the user can delete the run
        return HitSharePermission.has_delete_permission(request.user, obj) 

    def has_add_permission(self, request, obj: Union[Hit, None]=None):
        # prohibit runs from being added through admin panel
        return HitSharePermission.has_add_permission(request.user, obj) 

    def get_readonly_fields(self, request, obj=None):
        return list(set(
            [field.name for field in self.opts.local_fields] +
            [field.name for field in self.opts.local_many_to_many]
        ))

@admin.register(TaxonomyNode)
class TaxonomyNodeAdmin(admin.ModelAdmin):
    list_display = ('scientific_name', 'rank', 'id')
    readonly_fields = ('scientific_name', 'rank', 'id')
    search_fields = ('scientific_name', 'id')
    list_filter = ('rank',)

    def has_module_permission(self, request: HttpRequest) -> bool:
        return True 
    
    def has_view_permission(self, request: HttpRequest, obj: Optional[TaxonomyNode] = None) -> bool:
        return isinstance(request.user, User) and request.user.is_superuser

    def has_add_permission(self, request: HttpRequest, obj: Optional[TaxonomyNode] = None) -> bool:
        return False
    
    def has_delete_permission(self, request: HttpRequest, obj: Optional[TaxonomyNode] = None) -> bool:
        return isinstance(request.user, User) and request.user.is_superuser

    def has_change_permission(self, request: HttpRequest, obj: Optional[TaxonomyNode] = None) -> bool:
        return isinstance(request.user, User) and request.user.is_superuser

    def changeform_view(self, request: HttpRequest, object_id: Union[str, None], form_url: str, extra_context: Optional[Dict[str, bool]]) -> Any:
        extra_context = extra_context or {}

        extra_context['show_save_and_continue'] = False # show Save and Continue button
        extra_context['show_save_and_add_another'] = False # hide Save and Add Another button
        extra_context['show_save'] = False # hide save button
        extra_context['show_delete'] = True # show delete button
        return super().changeform_view(request, object_id, form_url, extra_context)