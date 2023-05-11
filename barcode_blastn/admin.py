from typing import Any, Dict, List, Union
from django.contrib import admin
from django.db import models
from django.db.models.query import QuerySet

from django.http.request import HttpRequest
from django.urls import reverse

from django.utils.html import format_html
from django.contrib.auth.models import User
from barcode_blastn.helper.parse_gb import retrieve_gb
from barcode_blastn.models import BlastDb, BlastQuerySequence, DatabaseShare, NuccoreSequence, BlastRun, Hit
from django.forms import BaseInlineFormSet, ModelForm, ValidationError
from barcode_blastn.database_permissions import DatabasePermissions
from barcode_blastn.permissions import DatabaseSharePermissions, HitSharePermission, NuccoreSharePermission, RunSharePermissions

from barcode_blastn.serializers import NuccoreSequenceSerializer

def fetch_data(accession_number: str) -> Dict[str, str]: 
    '''
    Fetch data from GenBank using the accession_number in obj, relying on function retrieve_gb.

    Raises an ValueError if there was an issue retrieving the data

    Returns a single dictionary of key-value pairs to populate a NuccoreSequence with.
    '''
    if accession_number is None:
        raise ValueError('Error no accession number found')
    currentData = {}
    genbank_data: Union[List[Dict], None] = []
    try:
        genbank_data = retrieve_gb(accession_numbers = [accession_number])
    except BaseException:
        raise ValueError('Error retrieving accession')
    
    try:
        assert(not genbank_data is None)
        assert(len(genbank_data) <= 1)
    except AssertionError:
        raise ValueError('More than one accession was found for this accession number. Differentiating between these accessions is currently not supported')
    
    try:
        assert(len(genbank_data) > 0)
    except AssertionError:
        raise ValueError('No data was found for this accession number')

    return genbank_data[0]

class SequenceFormset(BaseInlineFormSet):
    def save_new(self, form, commit=True):
        '''
        Callback when a new sequence is added through a save button on the admin form.
        '''
        obj = super(SequenceFormset, self).save_new(form, commit=False)
        accession_number = obj.accession_number

        # check if there is a duplicate
        try:
            NuccoreSequence.objects.get(accession_number = accession_number, owner_database_id = obj.owner_database.id)
        except NuccoreSequence.DoesNotExist:
            pass
        else:
            obj.accession_number = f'Error: duplicate for {accession_number}'
            return obj

        # fetch GenBank data
        currentData = fetch_data(accession_number=accession_number)

        # check that the GenBank data is valid
        try:
            assert not (currentData is None)
            serializer = NuccoreSequenceSerializer(data = currentData)
            if serializer.is_valid():
                saved = serializer.save(owner_database = obj.owner_database)
                return saved
                
        except AssertionError:
            obj.accession_number = f'Error: no data for {accession_number}'
        if commit:
            obj.save()
        return obj

    def save_existing(self, form, instance, commit=True):
        '''
        Callback when a sequence is edited through a save button in the admin form.
        '''
        obj = super(SequenceFormset, self).save_new(form, commit=False)
        accession_number = obj.accession_number
        try:
            duplicate = NuccoreSequence.objects.get(accession_number = accession_number, owner_database_id = obj.owner_database.id)
        except NuccoreSequence.DoesNotExist:
            pass
        else:
            obj.accession_number = f'Error: duplicate for {accession_number}'
            return obj
        
        currentData = fetch_data(accession_number=accession_number)

        try:
            assert not (currentData is None) # assert that there must be data retrieved
            for key, value in currentData.items():
                setattr(obj, key, str(value))
        except AssertionError:
            obj.accession_number = f'Error: no data for {accession_number}'
        if commit:
            obj.save()

        return obj

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

    def get_fields(self, request, obj: BlastDb):
        return ['grantee', 'perms', 'id']

    def get_readonly_fields(self, request, obj: Union[BlastDb, None] = None):
        return ['id']

    def has_view_permission(self, request: HttpRequest, obj: Union[BlastDb, None] = None) -> bool:
        return DatabaseSharePermissions.has_view_permission(request.user, obj) 

    def has_add_permission(self, request, obj: Union[BlastDb, None] = None) -> bool:
        '''
        Allow only the creator to add permissions.
        '''
        return DatabaseSharePermissions.has_delete_permission(request.user, obj)

    def has_delete_permission(self, request, obj: Union[BlastDb, None] = None) -> bool:
        '''
        Allow only the creator to delete permissions.
        '''
        return DatabaseSharePermissions.has_delete_permission(request.user, obj)

    def has_change_permission(self, request, obj: Union[BlastDb, None] = None) -> bool:
        '''
        Allow only the creator to change permissions.
        '''
        return DatabaseSharePermissions.has_delete_permission(request.user, obj)
        
@admin.register(BlastDb)
class BlastDbAdmin(admin.ModelAdmin):
    '''
    Admin page for BlastDb instances.
    '''
    # form = BlastDbForm

    inlines = [UserPermissionsInline, NuccoreSequenceInline]
    list_display = ('custom_name', 'owner', 'public', 'sequence_count', 'id')
    fields = ['custom_name', 'description', 'public', 'locked', 'owner']

    def sequence_count(self, obj):
        num_seqs: int = NuccoreSequence.objects.filter(owner_database=obj).count()
        return num_seqs
        
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

    def save_model(self, request, obj, form, change) -> None:
        if not(obj and obj.id):
            # set the owner to the current if the obj is being first created
            # (i.e. when no id exists)
            obj.owner = request.user
        return super().save_model(request, obj, form, change)

    def get_readonly_fields(self, request, obj: Union[BlastDb, None] = None):
        base = list(set(
            [field.name for field in self.opts.local_fields] +
            [field.name for field in self.opts.local_many_to_many]
        ))
        # specify which fields are editable
        excluded_fields = ['custom_name', 'description', 'public']
        if not obj or obj and not obj.locked:
            excluded_fields.append('sequences')
        base = [b for b in base if b not in excluded_fields and b != 'locked']
    
        return base

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

        # check if the accession already exists in the database
        duplicate_exists: bool = False
        if self.instance is None: # if this is an addition
            duplicate_exists = NuccoreSequence.objects.filter(accession_number = accession_number, owner_database_id = owner_database.id).exists()
        else: # if this is a modification of an existing entry
            duplicate_exists = NuccoreSequence.objects.filter(accession_number = accession_number, owner_database_id = owner_database.id).exclude(pk=self.instance.id).exists()
        
        if duplicate_exists:
            raise ValidationError(f'Error: Sequence entry for accession number {accession_number} already exists in the same database.')
        
        try:
            currentData = fetch_data(accession_number=accession_number)
        except ValueError as err:
            raise ValidationError('Could not retrieve data for this accession number, either due to server error or accession number not matching with a unique accession. Try again.') from err

        return        

@admin.register(NuccoreSequence)
class NuccoreAdmin(admin.ModelAdmin):
    '''
    Admin page for showing accession instances.
    '''
    list_display = (
        'accession_number', 'organism', 'specimen_voucher', 'id', 
        'owner_database_link'
    )
    fields_excluded = ['uuid']
    form = NuccoreAdminModifyForm

    # TODO: Make the owner_database field read_only if database is locked 

    def get_queryset(self, request: HttpRequest) -> QuerySet[Any]:
        if not isinstance(request.user, User) or not request.user.is_authenticated:
            return NuccoreSequence.objects.none()
        else:
            return NuccoreSequence.objects.viewable(request.user)

    def formfield_for_foreignkey(self, db_field, request: HttpRequest, **kwargs):
        if db_field.name == 'owner_database':
            if not request.user is None:
                if isinstance(request.user, User):
                    user: User = request.user
                    # set the available choices of database the accession can belong to
                    kwargs['queryset'] = BlastDb.objects.filter(
                        ~models.Q(shares=user, databaseshare__perms__startswith=DatabasePermissions.DENY_ACCESS) # omit databases the user is denied from
                        & 
                        (
                            models.Q(owner=request.user) |
                            models.Q(shares=user, databaseshare__perms__in=[DatabasePermissions.CAN_EDIT_DB])
                        ) # add databases the user can edit
                    )
        return super().formfield_for_foreignkey(db_field, request, **kwargs) 

    def get_fields(self, request, obj=None):
        fields = super(NuccoreAdmin, self).get_fields(request, obj)
        if 'uid' in fields:
            fields.remove('uid')
        if 'translation' in fields:
            fields.remove('translation')
        return fields

    def get_readonly_fields(self, request, obj=None):
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
        # allow module access if they can access databases
        return NuccoreSharePermission.has_module_permission(request.user)

    def has_view_permission(self, request: HttpRequest, obj: Union[NuccoreSequence, None]=None) -> bool:
        return NuccoreSharePermission.has_view_permission(request.user, obj)

    def has_change_permission(self, request, obj: Union[NuccoreSequence, None]=None):
        return NuccoreSharePermission.has_change_permission(request.user, obj)

    def has_delete_permission(self, request, obj: Union[NuccoreSequence, None]=None):
        return NuccoreSharePermission.has_delete_permission(request.user, obj)

    def has_add_permission(self, request, obj: Union[NuccoreSequence, None]=None):
        # sequences can only be added if the user can edit the database
        return NuccoreSharePermission.has_add_permission(request.user, obj)
    
    def save_model(self, request: Any, obj: Any, form: Any, change: Any) -> None:
        accession_number = obj.accession_number
        
        currentData = fetch_data(accession_number=accession_number)
        for key, value in currentData.items():
            setattr(obj, key, str(value))

        super(NuccoreAdmin, self).save_model(request, obj, form, change)

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