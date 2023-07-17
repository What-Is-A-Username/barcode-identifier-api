from rest_framework.pagination import PageNumberPagination

class BlastRunQueryPagination(PageNumberPagination):
    '''
    Pagination for lists of BlastQuerySequences.
    '''
    page_size = 100
    page_size_query_param = 'page_size'
    max_page_size = 1000

class BlastRunHitPagination(PageNumberPagination):
    '''
    Pagination for lists of Hits.
    '''
    page_size = 50
    page_size_query_param = 'page_size'
    max_page_size = 100

class BlastDbSequencePagination(PageNumberPagination):
    '''
    Pagination for lists of NuccoreSequences under \
        a BLAST database
    '''
    page_size = 100
    page_size_query_param = 'page_size'
    max_page_size = 1000