import io
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union
from django.contrib import admin
from django.db import models
from django import forms
from django.db.models.query import QuerySet
from django.db.models import Q
from django.template.response import TemplateResponse

from django.http.request import HttpRequest
from django.urls import reverse

from django.utils.html import format_html
from django.contrib.auth.models import User
from barcode_blastn.controllers.blastdb_controller import add_sequences_to_database, create_blastdb, delete_blastdb, delete_library, save_blastdb, save_sequence
from barcode_blastn.helper.parse_gb import GenBankConnectionError, InsufficientAccessionData, retrieve_gb
from barcode_blastn.models import Annotation, BlastDb, BlastQuerySequence, DatabaseShare, Library, NuccoreSequence, BlastRun, Hit, TaxonomyNode
from django.forms import BaseInlineFormSet, ModelForm, ValidationError
from barcode_blastn.permissions import DatabaseSharePermissions, HitSharePermission, LibrarySharePermissions, NuccoreSharePermission, RunSharePermissions

from barcode_blastn.serializers import LibraryEditSerializer, library_title, blast_db_title, run_title, nuccore_title, hit_title

@admin.register(TaxonomyNode)
class TaxonomyNodeAdmin(admin.ModelAdmin):
    list_display = ('scientific_name', 'rank', 'id')
    readonly_fields = ('scientific_name', 'rank', 'id')

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

class SequenceFormset(BaseInlineFormSet):

    def save(self, commit: bool = False) -> Any:
        return super().save(commit)

    def save_new(self, form, commit=True):
        '''
        Callback when a new sequence is added through a save button on the admin form.
        '''
        obj = super(SequenceFormset, self).save_new(form, commit=False)
        # TODO: Catch errors and/or enable validation
        return save_sequence(obj, change=False, commit=commit, raise_if_missing=True) 

    def clean(self) -> None:
        # TODO: validate in bulk
        super().clean()
        valid_forms = [
            form
            for form in self.forms
            if form.is_valid() and form not in self.deleted_forms
        ]
        accession_numbers = [form.cleaned_data.get('accession_number', '') for form in valid_forms]
        if any([len(acc) == 0 for acc in accession_numbers]):
            raise ValidationError(f'One or more accessions to be added had a missing accession number')
        try:
            retrieve_gb(accession_numbers=accession_numbers, raise_if_missing=True)
        except InsufficientAccessionData as exc:
            raise ValidationError(f'One or more accession numbers do not match a GenBank record: {", ".join(exc.missing_accessions)}')
        except ValueError:
            pass
        except BaseException as exc:
            raise ValidationError('Error validating accession numbers with GenBank.')
        else:
            return None

class NuccoreSequenceInline(admin.TabularInline):
    model = NuccoreSequence     
    formset = SequenceFormset   # Specify the form to handle edits and additions
    show_change_link = True     # Show link to page to edit the sequence
    classes = ['collapse']      # Allow the entries to be collapsed
    extra = 0                   # Show one extra row by default
    fields = ['accession_number', 'version', 'id', 'updated', 'organism', 'specimen_voucher', 'annotation_count']

    def annotation_count(self, obj):
        return Annotation.objects.filter(sequence=obj).count()

    def get_readonly_fields(self, request: HttpRequest, obj: Union[BlastDb, None] = None):
        base = list(set(
            [field.name for field in self.opts.local_fields] +  
            [field.name for field in self.opts.local_many_to_many]
        ))
        base.remove('accession_number')
        base.append('annotation_count')
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
        if obj is None:
            return DatabaseSharePermissions.has_change_permission(request.user, obj)
        else:
            return not obj.locked and DatabaseSharePermissions.has_change_permission(request.user, obj)

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
class LibraryAdmin(admin.ModelAdmin):
    '''
    Admin page for Library instances
    '''
    inlines = [UserPermissionsInline, VersionInline]
    list_display = ('custom_name', 'owner', 'public', 'latest_version', 'sequence_count', 'id')

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
    accession_list_upload = forms.FileField(help_text='Add large numbers of sequences at once by uploading a .txt file, with each accession number on a separate line.', required=False)
    accession_list_text = forms.CharField(widget=forms.Textarea, help_text='Add multiple sequences by pasting in a list of accession numbers, one per line.', required=False)
    base_library = forms.ModelChoiceField(queryset=BlastDb.objects.none(), help_text='Inherit all accession numbers from an existing BLAST database in order to populate the current.', required=False)

    def __init__(self, *args, **kwargs) -> None:
        super(BlastDbForm, self).__init__(*args, **kwargs)
        user = getattr(self, 'user', None)
        if user:
            self.fields['base_library'].queryset = BlastDb.objects.runnable(user)
        else:
            raise ValueError('No user')

    def clean(self) -> Dict[str, Any]:
        cleaned_data = super().clean()
        existing = NuccoreSequence.objects.filter(owner_database=self.instance)

        accessions_to_add = []

        # Check that the accesson numbers or accession.versions in the list_text are not duplicates
        accession_file_upload = cleaned_data.get('accession_list_upload')
        if accession_file_upload:
            accession_file_upload.seek(0)
            file_numbers = [line.decode('utf-8') for line in accession_file_upload.readlines()]
            # remove whitespace and newline characters
            file_numbers = [t.strip() for t in file_numbers]
            file_numbers = [t for t in file_numbers if len(t) > 0]
            # query the database for any entries with matching accession.versions or accession numbers
            exists_duplicate = existing.filter(Q(accession_number__in=file_numbers) | Q(version__in=file_numbers))
            if exists_duplicate.exists():
                raise ValidationError({'accession_list_upload': 'The file given for "accession list upload" contains accession numbers or accession.versions already in the database.'})
            accessions_to_add.extend(file_numbers)
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
            accessions_to_add.extend(list_text)

        # Check that accessions actually correspond with GenBank records
        try:
            retrieve_gb(accession_numbers=accessions_to_add, raise_if_missing=True)
        except InsufficientAccessionData as exc:
            raise ValidationError(f'One or more accession numbers or accession.versions do not match a GenBank record: {", ".join(exc.missing_accessions)}')
        except ValueError:
            pass
        except BaseException as exc:
            raise ValidationError('None of the accession numbers or accession.versions match a GenBank record.')

        return cleaned_data

@admin.register(BlastDb)
class BlastDbAdmin(admin.ModelAdmin):
    '''
    Admin page for BlastDb instances.
    '''
    inlines = [NuccoreSequenceInline]
    form = BlastDbForm
    list_display = ('custom_name', 'library', 'version_number', 'sequence_count', 'id', 'locked')
    fieldsets = [('Details', { 'fields': ['library', 'library_owner', 'custom_name', 'description', 'version_number']}), ('Visibility', { 'fields': ['locked', 'library_is_public']})]


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

    def get_fieldsets(self, request: HttpRequest, obj: Optional[BlastDb] = None) -> List[Tuple[Optional[str], Dict[str, Any]]]:
        f = super().get_fieldsets(request, obj).copy()
        if obj is None:
            f.append(('Accessions Included', { 'fields': ('base_library', 'accession_list_upload', 'accession_list_text')}))
        else:
            if not obj.locked and isinstance(request.user, User) and DatabaseSharePermissions.has_change_permission(request.user, obj=obj):
                f.append(('Add additional accessions', { 'fields': ('accession_list_upload', 'accession_list_text')}))
        return f

    def save_formset(self, request: Any, form: Any, formset: Any, change: Any) -> None:
        super_results = super().save_formset(request, form, formset, change)
        will_lock = form.cleaned_data['locked'] if 'locked' in form.cleaned_data else False
        if will_lock:
            # if user wants to lock the database, perform locking after saving the formset
            obj = form.save(commit=False)
            obj.locked = True 
            # directly call save_model on superclass to avoid locked setting to False
            save_blastdb(obj, perform_lock=True)
        return super_results

    def save_model(self, request: Any, obj: BlastDb, form: Any, change: Any) -> None:
        # Since the model is saved first before any changes in sequence, don't lock the
        # database yet. The database will be locked if needed by save_formset
        obj.locked = False

        accessions = []

        base = form.cleaned_data.get('base_library', None)

        accession_file_upload = form.cleaned_data.get('accession_list_upload')
        if accession_file_upload:
            accession_file_upload.seek(0)
            file_numbers = [line.decode('utf-8') for line in accession_file_upload.readlines()]
            file_numbers = [t.strip() for t in file_numbers]
            file_numbers = [t for t in file_numbers if len(t) > 0]
            accessions.extend(file_numbers)

        accession_list_text = form.cleaned_data.get('accession_list_text')
        if len(accession_list_text) > 0:
            list_text = accession_list_text.replace('\r\n', '\n').split('\n')
            list_text = [t.strip() for t in list_text]
            list_text = [t for t in list_text if len(t) > 0]
            accessions.extend(list_text)
        
        blastdb_fields = {
            'custom_name': form.cleaned_data.get('custom_name'),
            'description': form.cleaned_data.get('description')
        }
        if not change:
            create_blastdb(additional_accessions=accessions, base=base, database=obj, **blastdb_fields, library=obj.library)
        else:
            if len(accessions) > 0:
                add_sequences_to_database(obj, desired_numbers=accessions)
            else:
                obj.save()

    def get_readonly_fields(self, request, obj: Union[BlastDb, None] = None):
        if obj is None: # is adding
            return ['version_number', 'library_is_public', 'library_owner']
        else:
            return ['library', 'version_number', 'library_is_public', 'library_owner']

    def sequence_count(self, obj):
        num_seqs: int = NuccoreSequence.objects.filter(owner_database=obj).count()
        return num_seqs

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
        if not isinstance(request.user, User) or not request.user.is_authenticated:
            return BlastDb.objects.none()
        else:
            return BlastDb.objects.viewable(request.user)

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

    def __init__(self, *args, **kwargs):
        super(NuccoreAdminModifyForm, self).__init__(*args, **kwargs)
        instance = getattr(self, 'instance', None)
        # if instance and instance.pk:
        #     self.fields['accession_number'].widget.attrs['readonly'] = True
        #     self.fields['owner_database'].widget.attrs['readonly'] = True

    def clean(self):
        '''
        Validate new data in the form before saving it.

        If accession number does not correspond to entry, raise a ValidationError.

        Return nothing if no validation errors.
        '''

        cleaned_data = super().clean()
        instance: NuccoreSequence = self.instance
        # isAdding = instance._state.adding if not instance is None else 'instance none'
        # pk = instance.pk is None
        # raise ValidationError(str(instance is None) + ' ' + str(isAdding) + ' ' + str(pk) + ' ' + str(self.instance.accession_number is None))

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
                raise ValidationError(f'Could not retrieve a record for {accession_number}, likely due to accession number not matching with an existing and unique record. Check the accession number provided.')  
            
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

@admin.register(NuccoreSequence)
class NuccoreAdmin(admin.ModelAdmin):
    '''
    Admin page for showing accession instances.
    '''
    list_display = (
        'accession_number', 'version', 'organism', 'specimen_voucher', 'id', 'owner_database_link'
    )
    fields_excluded = ['uuid']
    form = NuccoreAdminModifyForm
    inlines = [AnnotationInline]

    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {nuccore_title} to view or change'}
        return super(NuccoreAdmin, self).changelist_view(request, extra_context=extra_context)

    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        # Restrict the sequences visible based on which databases they belong to, and
        # which libraries the user can see
        if not isinstance(request.user, User) or not request.user.is_authenticated:
            return NuccoreSequence.objects.none()
        else:
            return NuccoreSequence.objects.viewable(request.user)

    def formfield_for_foreignkey(self, db_field, request: Optional[HttpRequest], **kwargs):
        if request is None:
            return super().formfield_for_foreignkey(db_field, request, **kwargs)
        if db_field.name == 'owner_database':
            if not request.user is None:
                if isinstance(request.user, User):
                    # Set the available choices of database the accession can belong to
                    # based on the user's permission
                    kwargs['queryset'] = BlastDb.objects.editable(request.user)
        return super().formfield_for_foreignkey(db_field, request, **kwargs) 

    def get_fields(self, request, obj=None):
        fields = [
            'owner_database',
            'organism',
            'accession_number',
            'version',
            'definition',
            'organelle',
            'specimen_voucher',
            'id',
            'isolate',
            'country',
            'dna_sequence',
            'lat_lon',
            'type_material',
            'translation',
            'created',
            'updated'
        ]
        return fields

    def get_readonly_fields(self, request, obj: Optional[NuccoreSequence]=None):
        fields = [
            'organism', 'version', 'definition', 'organelle',
            'specimen_voucher', 'id', 'isolate', 'country',
            'dna_sequence', 'lat_lon', 'type_material',
            'created', 'updated', 'taxonomy', 'taxon_species', 
            'taxon_genus', 'taxon_family', 'taxon_order', 'taxon_class', 
            'taxon_phylum', 'taxon_kingdom', 'taxon_superkingdom','title', 
            'authors', 'journal'
        ]
        if not obj is None:
            fields.extend(['owner_database', 'accession_number'])
        fields.extend(['taxon'])
        return fields

    def get_fieldsets(self, request: HttpRequest, obj: Optional[NuccoreSequence] = None) -> List[Tuple[Optional[str], Dict[str, Any]]]:
        return [('Summary', { 'fields': ['id', 'accession_number', 'version', 'definition', 'taxonomy', 'owner_database']}), 
        ('Source Information', {'fields': ['specimen_voucher', 'type_material', 'organelle', 'isolate', 'country', 'lat_lon']}),
        ('History', { 'fields': ['created', 'updated']}),
        ('Taxonomy', { 'fields': ['taxonomy', 'taxon_species', 'taxon_genus', 'taxon_family', 'taxon_order', 'taxon_class', 'taxon_phylum', 'taxon_kingdom', 'taxon_superkingdom']}),
        ('Reference', { 'fields': ['title', 'authors', 'journal']}),
        ('Sequence', { 'fields': ['dna_sequence'] })
        ]
     
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
    
    def save_model(self, request: Any, obj: Any, form: Any, change: bool) -> None:
        if not change:
            save_sequence(obj, change=change, commit=True)

    def save_formset(self, request: Any, form: Any, formset: Any, change: Any) -> None:
        instances = formset.save(commit=False)
        for instance in instances:
            # Do something with `instance`
            if instance._state.adding:
                instance.poster = request.user
            instance.save() 
        formset.save_m2m()
        return super().save_formset(request, form, formset, change)

    def changeform_view(self, request: HttpRequest, object_id: Union[str, None], form_url: str, extra_context: Optional[Dict[str, bool]]) -> Any:
        extra_context = extra_context or {}

        extra_context['show_save_and_continue'] = True # show Save and Continue button
        extra_context['show_save_and_add_another'] = object_id is None # show Save and Add Another button
        extra_context['show_save'] = object_id is None # hide save button
        extra_context['show_delete'] = object_id is not None # show delete button
        return super().changeform_view(request, object_id, form_url, extra_context)

class BlastQuerySequenceInline(admin.StackedInline):
    '''
    Row to show query sequences under a BLAST run.
    '''
    model = BlastQuerySequence
    extra = 0

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
    fields = ['db_entry', 'query_sequence', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue_value', 'bit_score_value']
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

@admin.register(BlastRun)
class BlastRunAdmin(admin.ModelAdmin):
    '''
    Admin page for showing run information.
    '''
    inlines = [BlastQuerySequenceInline]
    show_change_link = True
    list_display = (
        'id', 'job_name', 'received_time', 'status'
    )

    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {run_title} to view or change'}
        return super(BlastRunAdmin, self).changelist_view(request, extra_context=extra_context)
    
    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        if not request.user.is_authenticated:
            return BlastRun.objects.none()
        elif isinstance(request.user, User):
            return BlastRun.objects.listable(request.user)
        else:
            return BlastRun.objects.none()        

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
    list_display = (
        'db_entry', 'id', 'owner_run_link'
    )

    def changelist_view(self, request, extra_context: Optional[Dict[str, str]] = None):
        # Customize the title at the top of the change list
        extra_context = {'title': f'Select a {hit_title} to view or change'}
        return super(HitAdmin, self).changelist_view(request, extra_context=extra_context)

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