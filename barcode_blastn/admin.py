from typing import Any, Dict, List, Optional, Union
from django.contrib import admin
from django.db import models
from django.db.models.query import QuerySet

from django.http.request import HttpRequest
from django.urls import reverse

from django.utils.html import format_html
from django.contrib.auth.models import User
from barcode_blastn.controllers.blastdb_controller import delete_blastdb, delete_library, save_sequence
from barcode_blastn.helper.parse_gb import GenBankConnectionError, InsufficientAccessionData, retrieve_gb
from barcode_blastn.models import BlastDb, BlastQuerySequence, DatabaseShare, Library, NuccoreSequence, BlastRun, Hit
from django.forms import BaseInlineFormSet, ModelForm, ValidationError
from barcode_blastn.permissions import DatabaseSharePermissions, HitSharePermission, LibrarySharePermissions, NuccoreSharePermission, RunSharePermissions

from barcode_blastn.serializers import BlastDbEditSerializer, LibraryEditSerializer, NuccoreSequenceSerializer

class SequenceFormset(BaseInlineFormSet):
    def save_new(self, form, commit=True):
        '''
        Callback when a new sequence is added through a save button on the admin form.
        '''
        obj = super(SequenceFormset, self).save_new(form, commit=False)
        return save_sequence(obj, change=False, commit=commit, raise_if_missing=True) 

    def save_existing(self, form, instance, commit=True):
        '''
        Callback when a sequence is edited through a save button in the admin form.
        '''
        obj = super(SequenceFormset, self).save_existing(form, instance, commit=False)
        return save_sequence(obj, change=True, commit=commit, raise_if_missing=True)

class NuccoreSequenceInline(admin.TabularInline):
    model = NuccoreSequence     
    formset = SequenceFormset   # Specify the form to handle edits and additions
    show_change_link = True     # Show link to page to edit the sequence
    classes = ['collapse']      # Allow the entries to be collapsed
    extra = 0                   # Show one extra row by default

    def get_fields(self, request: HttpRequest, obj: Union[BlastDb, None]):
        return [('accession_number'), ('id', 'created'), ('organism', 'definition', 'organelle'), ('specimen_voucher')]

    def get_readonly_fields(self, request: HttpRequest, obj: Union[BlastDb, None] = None):
        base = list(set(
            [field.name for field in self.opts.local_fields] +  
            [field.name for field in self.opts.local_many_to_many]
        ))
        if not obj or obj and not obj.locked: # if db not locked, exclude accession number
            base = [b for b in base if b != 'accession_number']
        return base
    
    def has_view_permission(self, request: HttpRequest, obj: Union[BlastDb, None]=None) -> bool:
        return DatabaseSharePermissions.has_view_permission(request.user, obj)
    
    def has_change_permission(self, request: HttpRequest, obj: Union[BlastDb, None] = None) -> bool:
        return DatabaseSharePermissions.has_change_permission(request.user, obj)

    def has_delete_permission(self, request, obj: Union[BlastDb, None] = None) -> bool:
        return self.has_change_permission(request, obj)

    def has_add_permission(self, request, obj: Union[BlastDb, None] = None) -> bool:
        return self.has_change_permission(request, obj)

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

@admin.register(Library)
class LibraryAdmin(admin.ModelAdmin):
    '''
    Admin page for Library instances
    '''
    inlines = [UserPermissionsInline, VersionInline]
    list_display = ('custom_name', 'owner', 'public', 'sequence_count', 'id')
    fields = LibraryEditSerializer.Meta.fields

    def sequence_count(self, obj: Library):
        blastdb: Union[BlastDb, None] = BlastDb.objects.latest(obj)
        if blastdb is None:
            return 0
        return NuccoreSequence.objects.filter(owner_database=blastdb).count()

    def delete_model(self, request: HttpRequest, obj: Library) -> None:
        delete_library(obj)

@admin.register(BlastDb)
class BlastDbAdmin(admin.ModelAdmin):
    '''
    Admin page for BlastDb instances.
    '''
    inlines = [NuccoreSequenceInline]
    list_display = ('library', 'version_number', 'sequence_count', 'id')

    def delete_model(self, request: HttpRequest, obj: BlastDb) -> None:
        '''Delete instance of database'''
        delete_blastdb(obj)
    
    def get_fields(self, request: HttpRequest, obj: Optional[BlastDb]):
        f: List[str] = BlastDbEditSerializer.Meta.fields.copy()
        additional_fields = ['library']
        if obj and obj.locked:
            f.extend(['genbank_version', 'major_version', 'minor_version'])
        for field in additional_fields:
            if not field in f:
                f.append(field)
        return f 

    def save_model(self, request: Any, obj: BlastDb, form: Any, change: Any) -> None:
        # TODO: Use logic in blastdb_controller.py
        raise NotImplementedError()

    def get_readonly_fields(self, request, obj: Union[BlastDb, None] = None):
        base = list(set(
            [field.name for field in self.opts.local_fields] +
            [field.name for field in self.opts.local_many_to_many]
        ))
        # specify which fields are editable
        editable_fields = BlastDbEditSerializer.Meta.fields.copy()
        base = [b for b in base if b not in editable_fields]
        # these fields are only readonly when user is modifying an existing page
        readonly_only_when_changing = ['library']
        if obj and obj.id:
            for key in readonly_only_when_changing:
                if not key in base:
                    base.append(key)
        else:
            for key in readonly_only_when_changing:
                if key in base:
                    base.remove(key)
        return base

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
        return super().formfield_for_foreignkey(db_field, request, **kwargs)
        
    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        if not isinstance(request.user, User) or not request.user.is_authenticated:
            return BlastDb.objects.none()
        else:
            return BlastDb.objects.viewable(request.user)

    def has_module_permission(self, request: HttpRequest) -> bool:
        super().has_module_permission
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
        accession_number: Any = cleaned_data.get("accession_number")
        # check if the form data contains an accession_number
        if accession_number is None:
            raise ValidationError(f'Could not locate an accession number in the submitted form data.')

        owner_database: Union[BlastDb, None] = cleaned_data.get('owner_database')

        # check that the sequence will be assigned to a database
        if owner_database is None:
            raise ValidationError(f'No owner database specified. Sequence must be associated with a database when being created. If no databases are selectable, contact a superuser to give you edit permission to a database.')
        else:
            if owner_database.locked:
                raise ValidationError(f'The database {owner_database} ({owner_database.id}) is locked, so changes cannot be made to its sequences and new sequences cannot be added to it.')

        # check if the accession already exists in the database
        duplicate_exists: bool = False
        if self.instance is None: # if this is an addition
            duplicate_exists = NuccoreSequence.objects.filter(accession_number = accession_number, owner_database_id = owner_database.id).exists()
        else: # if this is a modification of an existing entry
            duplicate_exists = NuccoreSequence.objects.filter(accession_number = accession_number, owner_database_id = owner_database.id).exclude(pk=self.instance.id).exists()
        
        if duplicate_exists:
            raise ValidationError(f'Error: Sequence entry for accession number {accession_number} already exists in the same database.')
        
        try:
            retrieve_gb(accession_numbers=[accession_number], raise_if_missing=True)
        except GenBankConnectionError as err:
            raise ValidationError('Encountered connection error while requesting data from GenBank. Please try again.') from err  
        except InsufficientAccessionData as err:
            raise ValidationError(f'Could not retrieve a record for {accession_number}, likely due to accession number not matching with an existing and unique record. Check the accession number provided.')  
        
        return None

@admin.register(NuccoreSequence)
class NuccoreAdmin(admin.ModelAdmin):
    '''
    Admin page for showing accession instances.
    '''
    list_display = (
        'accession_number', 'organism', 'specimen_voucher', 'id', 'owner_database_link'
    )
    fields_excluded = ['uuid']
    form = NuccoreAdminModifyForm

    # TODO: Make the owner_database field read_only if database is locked 

    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
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
                    # set the available choices of database the accession can belong to
                    kwargs['queryset'] = BlastDb.objects.editable(request.user)
        return super().formfield_for_foreignkey(db_field, request, **kwargs) 

    def get_fields(self, request, obj=None):
        fields = super(NuccoreAdmin, self).get_fields(request, obj)
        if 'uid' in fields:
            fields.remove('uid')
        if 'translation' in fields:
            fields.remove('translation')
        return fields

    def get_readonly_fields(self, request, obj: Optional[NuccoreSequence]=None):
        fields = [
            'organism',
            'definition',
            'organelle',
            'specimen_voucher',
            'id',
            'isolate',
            'country',
            'dna_sequence',
            'lat_lon',
            'type_material',
            'translation'
        ]
        return fields
     
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
        save_sequence(obj, change=change, commit=True)

    def changeform_view(self, request: HttpRequest, object_id: Union[str, None], form_url: str, extra_context: Optional[Dict[str, bool]]) -> Any:
        extra_context = extra_context or {}

        extra_context['show_save_and_continue'] = True # show Save and Continue button
        extra_context['show_save_and_add_another'] = True # show Save and Add Another button
        extra_context['show_save'] = True # hide save button
        extra_context['show_delete'] = True # show delete button
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
    fields = ['db_entry', 'owner_run', 'query_accession_version', 'subject_accession_version', 'percent_identity', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'sequence_start', 'sequence_end', 'evalue_value', 'bit_score_value']
    readonly_fields = ['evalue_value', 'bit_score_value']

    show_change_link = True

    def evalue_value(self, obj):
        return '%.4E' % obj.evalue

    def bit_score_value(self, obj):
        return '%.4E' % obj.bit_score

    extra = 0

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
    inlines = [HitInline, BlastQuerySequenceInline]
    show_change_link = True
    list_display = (
        'id', 'job_name', 'runtime'
    )
    
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
class HitAdmin(BlastRunAdmin):
    '''
    Admin page for showing Hit instances.
    '''
    inlines = []
    show_change_link = True
    list_display = (
        'db_entry', 'id', 'owner_run_link'
    )

    def owner_run_link(self, obj):
        url = reverse("admin:%s_%s_change" % ('barcode_blastn', 'blastrun'), args=(obj.owner_run.id,))
        return format_html("<a href='{url}'>{obj}</a>", url=url, obj=obj.owner_run)

    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        if not request.user.is_authenticated:
            return Hit.objects.none()
        elif isinstance(request.user, User):
            return Hit.objects.filter(owner_run__in=BlastRun.objects.listable(request.user))
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