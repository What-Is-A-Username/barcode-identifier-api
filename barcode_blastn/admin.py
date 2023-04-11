from django.contrib import admin
from django.urls import reverse
from django.utils.html import format_html
from barcode_blastn.helper.parse_gb import retrieve_gb
from barcode_blastn.models import BlastDb, BlastQuerySequence, NuccoreSequence, BlastRun, Hit
from django.forms import BaseInlineFormSet

from barcode_blastn.serializers import NuccoreSequenceSerializer

@admin.register(Hit)
class HitAdmin(admin.ModelAdmin):
    list_display = (
        'db_entry', 'id', 'owner_run_link'
    )

    def owner_run_link(self, obj):
        url = reverse("admin:%s_%s_change" % ('barcode_blastn', 'blastrun'), args=(obj.owner_run.id,))
        return format_html("<a href='{url}'>{obj}</a>", url=url, obj=obj.owner_run)

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return True 

    def has_add_permission(self, request, obj=None):
        return False 

@admin.register(NuccoreSequence)
class NuccoreAdmin(admin.ModelAdmin):
    list_display = (
        'accession_number', 'organism', 'specimen_voucher', 'id', 'owner_database_link'
    )
    fields_excluded = ['uuid']

    def get_fields(self, request, obj=None):
        fields = super(NuccoreAdmin, self).get_fields(request, obj)
        if 'uid' in fields:
            fields.remove('uid')
        if 'translation' in fields:
            fields.remove('translation')
        return fields

    def owner_database_link(self, obj):
        url = reverse("admin:%s_%s_change" % ('barcode_blastn', 'blastdb'), args=(obj.owner_database.id,))
        return format_html("<a href='{url}'>{obj}</a>", url=url, obj=obj.owner_database)

    def has_change_permission(self, request, obj=None):
        return False

    def has_delete_permission(self, request, obj=None):
        return True 

    def has_add_permission(self, request, obj=None):
        return False 

class BlastQuerySequenceInline(admin.StackedInline):
    model = BlastQuerySequence
    extra = 0

    def has_change_permission(self, request, obj=None):
        return False
    def has_add_permission(self, request, obj=None) -> bool:
        return False
    def has_delete_permission(self, request, obj=None) -> bool:
        return True

class HitInline(admin.TabularInline):
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
        return True

@admin.register(BlastRun)
class BlastRunAdmin(admin.ModelAdmin):
    inlines = [HitInline, BlastQuerySequenceInline]
    show_change_link = True
    list_display = (
        'id', 'job_name', 'runtime'
    )
    def has_add_permission(self, request) -> bool:
        return False 

    def get_readonly_fields(self, request, obj=None):
        return list(set(
            [field.name for field in self.opts.local_fields] +
            [field.name for field in self.opts.local_many_to_many]
        ))

class SequenceFormset(BaseInlineFormSet):
    def fetch_data(self, obj: NuccoreSequence): 
        accession_number = obj.accession_number
        currentData = {}
        genbank_data = []
        try:
            genbank_data = retrieve_gb(accession_numbers = [accession_number])
            assert len(genbank_data) > 0  
        except BaseException:
            obj.organism = f'Error: could not find data for "{accession_number}".'

        try:
            print(genbank_data)
            assert len(genbank_data) == 1
        except AssertionError:
            return {'organism': f'Error: could not find data for "{accession_number}".'}
        else:
            return genbank_data[0]

    def save_new(self, form, commit=True):
        obj = super(SequenceFormset, self).save_new(form, commit=False)
        accession_number = obj.accession_number
        try:
            duplicate = NuccoreSequence.objects.get(accession_number = accession_number, owner_database_id = obj.owner_database.id)
        except NuccoreSequence.DoesNotExist:
            pass
        else:
            obj.organism = 'Error: duplicate accession number in database.'
            return obj

        currentData = self.fetch_data(obj)

        try:
            assert not (currentData is None)
            serializer = NuccoreSequenceSerializer(data = currentData)
            if serializer.is_valid():
                saved = serializer.save(owner_database = obj.owner_database)
                return saved
        except AssertionError:
            obj.organism = 'Error: could not retrieve data for new entry.'
        if commit:
            obj.save()
        return obj

    def save_existing(self, form, instance, commit=True):
        obj = super(SequenceFormset, self).save_new(form, commit=False)
        accession_number = obj.accession_number
        try:
            duplicate = NuccoreSequence.objects.get(accession_number = accession_number, owner_database_id = obj.owner_database.id)
        except NuccoreSequence.DoesNotExist:
            pass
        else:
            obj.organism = 'Error: duplicate accession number in database.'
            return obj
        
        currentData = self.fetch_data(obj)

        try:
            assert not (currentData is None)
            for key, value in currentData.items():
                setattr(obj, key, str(value))
        except AssertionError:
            pass
        if commit:
            obj.save()
        return obj

class NuccoreSequenceInline(admin.TabularInline):
    model = NuccoreSequence
    extra = 1
    formset = SequenceFormset
    fields = [('accession_number'), ('id', 'created'), ('organism', 'definition', 'organelle'), ('specimen_voucher')]
    show_change_link = True 

    def get_readonly_fields(self, request, obj=None):
        base = list(set(
            [field.name for field in self.opts.local_fields] +
            [field.name for field in self.opts.local_many_to_many]
        ))
        if not obj or obj and not obj.locked: # if db not locked, exclude accession number
            base = [b for b in base if b != 'accession_number']
        return base

    def has_delete_permission(self, request, obj=None):
        return obj and not obj.locked

    def has_add_permission(self, request, obj=None):
        return not obj or obj and not obj.locked
        
@admin.register(BlastDb)
class BlastDbAdmin(admin.ModelAdmin):

    inlines = [NuccoreSequenceInline]
    list_display = (
        'custom_name', 'id'
    )
    def get_readonly_fields(self, request, obj=None):
        base = list(set(
            [field.name for field in self.opts.local_fields] +
            [field.name for field in self.opts.local_many_to_many]
        ))
        excluded_fields = ['custom_name', 'description']
        if not obj or obj and not obj.locked:
            excluded_fields.append('sequences')
        base = [b for b in base if b not in excluded_fields and b != 'locked']
    
        return base