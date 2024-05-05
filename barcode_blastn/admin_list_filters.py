from typing import List, Optional, Tuple
from django.contrib import admin
from django.db.models import QuerySet

from barcode_blastn.models import Library, NuccoreSequence

class NuccoreSequencePublicationFilter(admin.SimpleListFilter):
    '''
    Allow sequences to be filtered by the reference publication
    they belong to.
    '''
    title = 'Reference Publication'
    parameter_name = 'reference'
    def lookups(self, request, model_admin) -> List[Tuple[str, str]]:
        '''
        Return a list of valid values to filter by.
        '''
        qs: QuerySet[NuccoreSequence] = model_admin.get_queryset(request)
        distinct = qs.order_by().values('title', 'journal', 'authors').distinct()
        options = [('all', 'All')] + [(str(d['title']), f"{d['journal']}, {d['title']}") for d in distinct]
        return options

    def queryset(self, request, queryset: QuerySet[NuccoreSequence]) -> Optional[QuerySet[NuccoreSequence]]:
        '''
        Return a filtered queryset based on the value specified for
        reference filter.
        '''
        param = self.value()
        if param is None:
            return queryset
        else:
            return queryset.filter(title=param)

class BlastRunLibraryFilter(admin.SimpleListFilter):
    '''
    Allow BlastRuns to be filtered by the reference library they were run on.
    '''
    title = 'Reference Library'
    parameter_name = 'library'

    def lookups(self, request, model_admin) -> List[Tuple[str, str]]:
        '''
        '''
        qs: QuerySet[NuccoreSequence] = model_admin.get_queryset(request)
        # Return a list of all libraries utilized by the presently-shown runs
        libraries = Library.objects.distinct().filter(blastdb__usages__in=qs)
        options = [(str(library.id), library.custom_name,) for library in libraries]
        return options

    def queryset(self, request, queryset: QuerySet[NuccoreSequence]) -> Optional[QuerySet[NuccoreSequence]]:
        param = self.value()
        if param is None:
            return queryset
        else:
            return queryset.filter(db_used__library__id=param)