from django.contrib import admin

from barcode_blastn.models import BlastDb, NuccoreSequence, BlastRun, Hit

# Register your models here.

admin.site.register(BlastDb)
admin.site.register(NuccoreSequence)
admin.site.register(BlastRun)
admin.site.register(Hit)