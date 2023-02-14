from io import StringIO
from rest_framework import status, generics, mixins 
from urllib import request
from barcode_blastn.models import BlastRun, Hit
from barcode_tree.models import ResultTree
from rest_framework.response import Response
from Bio import SeqIO
import shutil
from Bio.Seq import Seq
import os
import json
from barcode_tree.modified_clustalo import getMultipleAlignmentResult, serviceGetStatus, submitMultipleAlignmentAsync
from barcode_tree.modified_simple_phylogeny import checkSimplePhylogenyStatus, getSimplePhylogenyOutput, readTreeFromFile, submitSimplePhylogenyAsync

from barcode_tree.serializers import ResultTreeCreatorSerializer, ResultTreeDetailSerializer 

# Create your views here.

class ResultTreeDetail(mixins.DestroyModelMixin, generics.GenericAPIView):
    queryset = ResultTree.objects.all()
    serializer_class = ResultTreeCreatorSerializer
    '''
        View and create phylogenic trees for blast output sequences.
    '''

    def get(self, request, *args, **kwargs):
        '''
            Retrieve tree data, if it exists.
        '''

        run_id = kwargs['run']

        # Check if a ResultTree object has been made 
        try:
            tree : ResultTree = ResultTree.objects.get(owner_run=run_id)
        except ResultTree.DoesNotExist:
            # If there is no object, return nothing
            return Response(status=status.HTTP_404_NOT_FOUND);       
    
        if tree.internal_status == ResultTree.TreeStatus.ALIGNING:
            # If the last status indicated that multiple sequence alignment was occurring

            # Check if the job finished

            operationFinished = serviceGetStatus(tree.alignment_job_id)
            # if operation is finished, change status and prepare data for tree construction
            if operationFinished:
                tree.internal_status = ResultTree.TreeStatus.PROCESSING
                tree.save()
                # save clustal files locally
                try:
                    getMultipleAlignmentResult(job_id=tree.alignment_job_id, run_id=run_id)
                except BaseException:
                    tree.internal_status = ResultTree.TreeStatus.ERRORED
                    tree.save()
                    return Response(status=status.HTTP_500_INTERNAL_SERVER_ERROR)
                else:
                    # TODO: Confirm operation was successful by checking for network error, exceptions
                    
                    # move files to be downloaded
                    destination_dir = f'/var/www/runs/{run_id}'
                    parent_folder = os.path.abspath(f'./runs/{run_id}')
                    files = os.listdir(parent_folder)
                    files_to_transfer = ['.aln-clustal_num.clustal_num', '.phylotree.ph', '.pim.pim', '.sequence.txt']
                    for file in files:
                        if any(file.endswith(extension) for extension in files_to_transfer):
                            shutil.copy(f'{parent_folder}/{file}', destination_dir)
                            
                    
                    tree.internal_status = ResultTree.TreeStatus.FINISHED
                    tree.save()

        serializer = ResultTreeDetailSerializer(tree)
        data = serializer.data
        # read the tree string from file, return an empty string otherwise
        data['tree'] = readTreeFromFile(run_id) if tree.internal_status == ResultTree.TreeStatus.FINISHED else ''
        response = Response(data, status = status.HTTP_200_OK)

        return response

    '''
        View and create phylogenic trees for blast output sequences.
    '''
    def post(self, request, *args, **kwargs):
        '''
            Submit a request to perform tree construction
        '''

        run_id = kwargs['run']

        # Check that the run exists
        try:
            run : BlastRun = BlastRun.objects.get(id=run_id)
        except BlastRun.DoesNotExist:
            return Response(
                {'message': 'Could not find a run result with the given id.'},
                status=status.HTTP_400_BAD_REQUEST
            )

        # If a Result tree object already exists for this run, don't do anything and return a response
        try:
            tree : ResultTree = ResultTree.objects.get(owner_run=run_id)
        except ResultTree.DoesNotExist:
            pass
        else:
            return Response(
                {'message': 'Already started performing tree construction. Cannot start tree construction more than once.' }, 
                status=status.HTTP_400_BAD_REQUEST)

        # Gather sequences (query sequences + hit sequences) into a single FASTA
        run_folder = os.path.abspath(f'./runs/{run_id}')

        # Gather query sequences as a list
        s = list(SeqIO.parse(f'{run_folder}/query.fasta', 'fasta'))
        
        hit : Hit
        # Query over all hits
        hits = list(run.hits.all())  # type: ignore
        for hit in hits:
            s.append(SeqIO.SeqRecord(seq=Seq(hit.db_entry.dna_sequence), id=hit.db_entry.accession_number, description=hit.db_entry.definition))

        sequence_string = StringIO()
        SeqIO.write(s, sequence_string, 'fasta')

        # The following code submits the fasta file to ClustalO at EMBL-EBI for multiple sequence alignment, and is adapted from clustalo.py at https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation
        # ClustalO:
        # TODO: Add citation to web interface
        #   Sievers F., Wilm A., Dineen D., Gibson T.J., Karplus K., Li W., Lopez R., McWilliam H., Remmert M., Söding J., Thompson J.D. and Higgins D.G. (2011)
        # Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. 
        # Mol. Syst. Biol. 7:539
        # PMID:  21988835 
        # DOI:  10.1038/msb.2011.75
        # TODO: Find citation for EMBL-EBI

        # This runs the tool, mimicking the following command line calls from the original Python wrapper
        # python clustalo.py --asyncjob --email email@domain.com ./aggregate.fasta
        # This checks if its done
        # python clustalo.py --status --jobid clustalo-R20230207-011805-0893-96632125-p2m
        # This outputs the data 
        # python clustalo.py --polljob --jobid clustalo-R20230207-011805-0893-96632125-p2m
        # job_id: str = submitClustalOAsyncJob(f'runs/{run_id}/aggregate.fasta', run_id)
        job_id: str = submitMultipleAlignmentAsync(sequence=sequence_string.getvalue(), run_id=run_id)

        tree = ResultTree(owner_run=run, internal_status=ResultTree.TreeStatus.ALIGNING, alignment_job_id=job_id)
        tree.save()

        serializer = ResultTreeDetailSerializer(tree)

        # return success
        return Response(serializer.data, status=status.HTTP_201_CREATED); 

