from rest_framework import status, generics, mixins 
from urllib import request

# Create your views here.

class ResultTreeDetail(mixins.CreateModelMixin, mixins.DestroyModelMixin, generics.GenericAPIView):
    '''
        View and create phylogenic trees for blast output sequences.
    '''

    def get(self, request, *args, **kwargs):
        '''
            Retrieve tree data, if it exists.
        '''

        # Check if a ResultTree object has been made 

        # If there is no object, return nothing

        # If there is an object, then based on status,
            # check if tree result already exists in the file system
                # return it if it exists
            # check if multiple sequence alignment exists in system (ALTERNATIVE: use celery worker)
                # if job not yet submitted, submit to simple phylogeny
                # if job submitted, poll for result
                    # if response has result, write result to file system and return the response
                    # else return status
            # throw error 

        pass

    def post(self, request, *args, **kwargs):
        '''
            Submit a request to perform tree construction
        '''

        # If a Result tree object already exists for this run, don't do anything and return a response

        # Gather sequences (query sequences + hit sequences) into a single FASTA

        # Submit the fasta for multiple sequence alignment https://www.ebi.ac.uk/Tools/msa/ with ClustalOmega

        # return success

        pass

