# TODO: Install Python 3.10.6 to match EC2

import io
from celery import signature
from barcode_identifier_api.celery import app
from typing import Dict, List
import uuid
import os
import shutil
from Bio import SeqIO
from barcode_blastn.helper.parse_gb import retrieve_gb
from barcode_blastn.helper.verify_query import verify_dna
from barcode_blastn.models import BlastQuerySequence, BlastRun, Hit, NuccoreSequence, BlastDb
from barcode_blastn.permissions import IsAdminOrReadOnly
from barcode_blastn.renderers import BlastRunCSVRenderer, BlastRunTxtRenderer, BlastDbCSVRenderer, BlastDbFastaRenderer, BlastRunFastaRenderer
from barcode_blastn.serializers import BlastDbCreateSerializer, BlastRunRunSerializer, BlastRunSerializer, BlastRunStatusSerializer, HitSerializer, NuccoreSequenceAddSerializer, NuccoreSequenceBulkAddSerializer, NuccoreSequenceSerializer, BlastDbSerializer, BlastDbListSerializer
from rest_framework import status, generics, mixins
from rest_framework.parsers import JSONParser, FormParser, MultiPartParser
from rest_framework.permissions import IsAdminUser, AllowAny
from rest_framework.renderers import JSONRenderer, BrowsableAPIRenderer, TemplateHTMLRenderer
from rest_framework.response import Response
from rest_framework.serializers import Serializer
from urllib.error import HTTPError
from barcode_blastn.tasks import run_blast_command

from drf_yasg.utils import swagger_auto_schema
from drf_yasg import openapi

tag_runs = 'Runs/Jobs'
tag_blastdbs = 'BLAST Databases'
tag_sequences = 'GenBank Accessions'
tag_admin = 'Admin Tools'

class NuccoreSequenceList(mixins.ListModelMixin, generics.GenericAPIView):
    '''
    List all the sequences in the server, irrespective of database
    '''
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer
    permission_classes = [IsAdminUser]

    @swagger_auto_schema(
        operation_summary='Global view of all sequence database entries.',
        operation_description='Return a list of all accession numbers across all databases.',
        tags = [tag_admin],
        responses={
            '200': openapi.Response(
                description="All sequence entries",
                # TODO: Schema and example should be an array
                schema=NuccoreSequenceSerializer(many=True),
                examples={
                    'application/json': [NuccoreSequenceSerializer.Meta.example]
                }
            )
        }
    )
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

class NuccoreSequenceBulkAdd(generics.CreateAPIView):   
    '''
    Bulk add sequences to a specified database
    '''
    serializer_class = NuccoreSequenceBulkAddSerializer
    permission_classes = [ IsAdminUser ]
    queryset = NuccoreSequence.objects.all()
    
    @swagger_auto_schema(
        operation_summary='Bulk-add accessions.',
        operation_description='From a list of accession numbers, bulk add them to an existing database. List must contain between 1-100 accession numbers inclusive.',
        tags = [tag_blastdbs, tag_sequences],
        # TODO: Add request schema and example
        request_body=openapi.Schema(
            type=openapi.TYPE_OBJECT,
            properties={
                'accession_numbers': openapi.Schema(
                    type=openapi.TYPE_ARRAY, 
                    description='List of GenBank accession numbers to add',
                    items=openapi.Schema(
                        type=openapi.TYPE_STRING,
                        example='GU701771',
                        description='GenBank accession number'
                    )
                ),
            }
        ),
        responses={
            '200': openapi.Response(
                description='No new sequences added because all already exist in the database.',
                schema=NuccoreSequenceSerializer(many=True),
                examples={
                    'application/json': [ NuccoreSequenceSerializer.Meta.example ]
                }
            ),
            '201': 'Sequences successfully added to the database.',
            '400': 'Bad parameters in request. Example reasons: list may be empty or too long (>100); accession numbers may be invalid or duplicate.',
            '404': 'BLAST database matching the specified ID was not found.',
            '500': 'Unexpected error.',
            '502': 'Encountered error connecting to GenBank.'
        }
    )
    def post(self, request, *args, **kwargs):
        '''
        Create a new accession number and add it to an existing database
        '''
        desired_numbers = request.data['accession_numbers']
        if len(desired_numbers) == 0:
            return Response({'message': 'The list of numbers to add cannot be empty.', 'accession_numbers': str(desired_numbers)}, status.HTTP_400_BAD_REQUEST)

        pk = kwargs['pk']
        # Query for the Blast database
        try:
            db = BlastDb.objects.get(id=pk)
        except BlastDb.DoesNotExist:
            return Response({'message': 'Database does not exist', 'requested': pk}, status.HTTP_404_NOT_FOUND)
        
        # Check what accession numbers are duplicate (i.e. already existing in the datab)
        existing = NuccoreSequence.objects.filter(
            owner_database_id = pk
        ).filter(
            accession_number__in = desired_numbers
        )
            
        # return an OK status if all the numbers to be added already exist
        if existing.count() == len(desired_numbers):
            return Response({'message': "Every accession number specified to be added already exists in the database.", 'accession_numbers': desired_numbers}, status=status.HTTP_200_OK)

        existing_numbers = [dup.accession_number for dup in existing]
        missing_numbers = [number for number in desired_numbers if number not in existing_numbers ]

        # retrieve data
        genbank_data : List[Dict]
        try:
            genbank_data = retrieve_gb(accession_numbers = missing_numbers)
        except ValueError:
            return Response({
                    'message': f"Bulk addition of sequences requires a list of accession numbers and the list length must have 1 to 100 elements (inclusive).", 
                    "accession_numbers": desired_numbers 
                }, 
                status=status.HTTP_400_BAD_REQUEST)
        except HTTPError:
            return Response({
                    'message': f"The GenBank database could not be connected to.", 
                    "accession_numbers": desired_numbers, 
                    "queried_numbers": missing_numbers
                }, 
                status=status.HTTP_502_BAD_GATEWAY)
        except BaseException:
            return Response({'message': f"Failed to parse GenBank database response when requested for a list of accession numbers", 
                    "accession_numbers": desired_numbers, 
                    "queried_numbers": missing_numbers
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR)
                
        # check if we are missing any data
        failed_queries : List[str]
        failed_queries = [entry["accession_number"] for entry in genbank_data]
        # stop execution if data not complete
        if len(failed_queries) < len(missing_numbers):
            return Response({
                    'message': f"The GenBank database failed to return some accession numbers so no numbers were added at all.", 
                    "accession_numbers": desired_numbers,
                    "queried_numbers": failed_queries
                }, 
                status=status.HTTP_400_BAD_REQUEST)
        
        # save the new sequences using a bulk operation
        all_new_sequences = [NuccoreSequence(owner_database=db, **data) for data in genbank_data]
        created_sequences = NuccoreSequence.objects.bulk_create(all_new_sequences)

        # serialize the list of created objects so they can be sent back
        serialized_data = NuccoreSequenceSerializer(created_sequences, many = True).data

        return Response(serialized_data, status = status.HTTP_201_CREATED)

'''
Add a sequence to a specified database
'''
class NuccoreSequenceAdd(generics.CreateAPIView):
    serializer_class = NuccoreSequenceAddSerializer
    permission_classes = [IsAdminUser]
    queryset = NuccoreSequence.objects.all()
    
    @swagger_auto_schema(
        operation_summary='Add a sequence entry.',
        operation_description='Add a sequence to a BLAST database, using GenBank accession data corresponding to the accession number given.',
        request_body=NuccoreSequenceAddSerializer,
        tags = [tag_blastdbs, tag_sequences],
        responses={
            '201': openapi.Response(
                description='Sequences successfully added to the database.',
                schema=NuccoreSequenceSerializer(),
                examples={
                    'application/json': NuccoreSequenceSerializer.Meta.example
                }
            ),
            '400': 'Bad parameters in request. Example reasons: accession number was not found on GenBank; an entry with the same accession number already exists in the BLAST database.',
            '404': 'BLAST database matching the specified ID was not found.',
            '500': 'Unexpected error.',
            '502': 'Encountered error connecting to GenBank.'
        }
    )
    def post(self, request, *args, **kwargs):
        '''
        Create a new accession number and add it to an existing database
        '''
        accession_number = request.data['accession_number']

        pk = kwargs['pk']
        # Query for the Blast database
        try:
            db = BlastDb.objects.get(id=pk)
        except BlastDb.DoesNotExist:
            return Response({'message': 'Database does not exist', 'requested': pk}, status.HTTP_404_NOT_FOUND)
        
        # Check if the database already has this accession number
        try:
            duplicate = NuccoreSequence.objects.get(accession_number = accession_number, owner_database_id = pk)
        except NuccoreSequence.DoesNotExist:
            pass
        else:
            return Response({'message': "Entry already exists for the accession number specified.", 'accession_number': accession_number}, status.HTTP_400_BAD_REQUEST)

        currentData = {}
        try:
            genbank_data = retrieve_gb(accession_numbers = [accession_number])
            assert len(genbank_data) > 0  
        except AssertionError:
            return Response({'message': f"The GenBank database did not return an entry for the accession number {accession_number}"}, status=status.HTTP_502_BAD_GATEWAY)
        except HTTPError:
            return Response({'message': f"The GenBank database does not contain an entry for the accession number {accession_number}"}, status=status.HTTP_400_BAD_REQUEST)
        except BaseException:
            return Response({'message': f"Failed to parse retrieved GenBank database entry for the accession number {accession_number}"}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        
        currentData = genbank_data[0]
        currentData['owner_database'] = db.id
        all_data = request.POST.copy()
        all_data.update(currentData)

        serializer = NuccoreSequenceSerializer(data = all_data)
        if serializer.is_valid():
            serializer.save(owner_database = db)
            return Response(serializer.data, status = status.HTTP_201_CREATED)

        return Response(serializer.errors, status = status.HTTP_500_INTERNAL_SERVER_ERROR)

class NuccoreSequenceDetail(mixins.DestroyModelMixin, generics.RetrieveAPIView):
    '''
    Get or delete a sequence specified by ID
    '''
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer
    permission_classes = [IsAdminOrReadOnly]

    @swagger_auto_schema(
        operation_summary='Get information about a sequence entry.',
        operation_description='Get information about a sequence entry.',
        tags = [tag_sequences],
        responses={
            '200': openapi.Response(
                description='Successfully returned sequence information',
                schema=NuccoreSequenceSerializer(read_only=True),
                examples={
                    'application/json': NuccoreSequenceSerializer.Meta.example
                }
            ),
            '404': 'Sequence does not exist',
            '500': 'Unexpected error.',
        }
    )
    def get(self, request, *args, **kwargs):
        db_primary_key = kwargs['pk']

        try:
            seq = NuccoreSequence.objects.get(id = db_primary_key)
        except NuccoreSequence.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        
        serializer = NuccoreSequenceSerializer(seq)
        return Response(serializer.data, status = status.HTTP_200_OK)

    @swagger_auto_schema(
        operation_summary='Delete a sequence entry.',
        operation_description='Delete a sequence entry from a BLAST database.',
        tags = [tag_blastdbs, tag_sequences],
        responses={
            '204': 'Deletion successful.',
            '400': 'Deletion failed because entry belongs to a database that is locked for editing.',
            '404': 'Sequence does not exist',
            '500': 'Unexpected error.',
        }
    )
    def delete(self, request, *args, **kwargs):
        pk = kwargs['pk']
        try:
            seq : NuccoreSequence = NuccoreSequence.objects.get(id=pk)
        except NuccoreSequence.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        db : BlastDb = seq.owner_db   # type: ignore
        if db.locked:
            return Response({'message': 'Cannot remove sequences from a locked database.'}, status=status.HTTP_400_BAD_REQUEST)
        return self.destroy(request, *args, **kwargs)

class BlastDbList(mixins.ListModelMixin, generics.CreateAPIView):
    '''
    Return a list of all blast databases
    '''
    queryset = BlastDb.objects.all()
    serializer_class = BlastDbCreateSerializer
    permission_classes = [IsAdminOrReadOnly]

    def get_serializer_class(self):
        if self.request is None:
            return BlastDbListSerializer
        elif self.request.method == 'POST':
            return BlastDbCreateSerializer
        else:
            return BlastDbListSerializer

    @swagger_auto_schema(
        operation_summary='Get all databases.',
        operation_description='Return a list of all BLAST databases publicly available for queries.',
        tags = [tag_blastdbs],
        responses={
            '200': openapi.Response(
                description='A list of all accession numbers saved to all databases.',
                schema=BlastDbListSerializer(many=True, read_only=True),
                examples={
                    'application/json': BlastDbListSerializer.Meta.example,
                }
            )
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        List all blast databases
        '''
        return self.list(request, *args, **kwargs)

    @swagger_auto_schema(
        operation_summary='Create a database.',
        operation_description='Create a new publicly accessible BLAST database available for queries.',
        tags = [tag_admin, tag_blastdbs],
        request_body=BlastDbCreateSerializer(),
        responses={
            '200': openapi.Response(
                description='Creation was successful.',
                schema=BlastDbCreateSerializer(),
                examples={
                    'application/json': BlastDbCreateSerializer.Meta.example
                }
            )
        }
    )
    def post(self, request, *args, **kwargs):
        '''
        Create a new blast database.
        '''
        serializer = self.get_serializer(data=request.data)
        if serializer.is_valid():
            custom_name = serializer.validated_data['custom_name']  # type: ignore
            description = serializer.validated_data['description']  # type: ignore
            blast_db = BlastDb(custom_name=custom_name, description=description)
            blast_db.save()
            return Response(BlastDbCreateSerializer(blast_db).data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class BlastDbDetail(mixins.RetrieveModelMixin, mixins.UpdateModelMixin, mixins.DestroyModelMixin, generics.GenericAPIView):
    '''
    Retrieve the details of a given blast database
    '''
    queryset = BlastDb.objects.all()
    serializer_class = BlastDbSerializer
    permission_classes = [IsAdminOrReadOnly]
    renderer_classes = [JSONRenderer, BrowsableAPIRenderer, TemplateHTMLRenderer, BlastDbCSVRenderer, BlastDbFastaRenderer]

    def get_serializer_class(self):
        return super().get_serializer_class()

    @swagger_auto_schema(
        operation_summary='Get information about a BLAST database.',
        operation_description='Return all available information of the BLAST database specified by the ID given.',
        tags = [tag_blastdbs],
        responses={
            '200': openapi.Response(
                description='Information on the BLAST database matching the given ID.',
                schema=BlastDbSerializer(),
                examples={
                    'application/json': BlastDbSerializer.Meta.example
                }
            ),
            '404': 'Database with the ID does not exist',
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        Retrieve data in the blastdb.
        '''

        db_primary_key = kwargs['pk']

        try:
            db = BlastDb.objects.get(id = db_primary_key)
        except BlastDb.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        
        response = Response(self.get_serializer_class()(db).data, status=status.HTTP_200_OK, template_name='blastdb.html')

        # based on the media type file returned, specify attachment and file name
        file_name_root = db.custom_name
        if request.accepted_media_type.startswith('text/csv'):
            response['Content-Disposition'] = f'attachment; filename={file_name_root}.csv";'
        elif request.accepted_media_type.startswith('text/x-fasta'):
            response['Content-Disposition'] = f'attachment; filename="{file_name_root}.fasta";'

        return response

    @swagger_auto_schema(
        operation_summary='Update information of a BLAST database.',
        operation_description='Edit information of a BLAST database. This method does not allow adding/removing/editing of the sequence entries within.',
        tags = [tag_blastdbs],
        responses={
            '200': openapi.Response(
                description='Database updated successfully.',
                schema=BlastDbSerializer,
                examples={
                    'application/json': BlastDbSerializer.Meta.example
                }
            ),
            '400': 'Bad request parameters.',
            '404': 'Database with the ID does not exist',
        }
    )
    def patch(self, request, *args, **kwargs):
        '''
        Partially update this blastdb.
        '''
        try:
            db = BlastDb.objects.get(id = kwargs['pk'])
        except BlastDb.DoesNotExist:
            return Response("Resource does not exist", status = status.HTTP_404_NOT_FOUND)
        
        if db.locked and ('locked' in request.data and request.data['locked'] != db.locked):
            return Response("This entry is locked and cannot be patched.", status = status.HTTP_400_BAD_REQUEST)

        return self.partial_update(request, *args, **kwargs)

    @swagger_auto_schema(
        operation_summary='Delete a BLAST database.',
        operation_description='Delete the BLAST database specified by the given ID.',
        tags = [tag_blastdbs],
        responses={
            '204': 'Deletion successful.',
            '404': 'BLAST database with the given ID does not exist',
            '500': 'Unexpected error.',
        }
    )
    def delete(self, request, *args, **kwargs):
        '''
        Delete the blastdb
        '''
        try:
            database_id = str(kwargs['pk'])
            db = BlastDb.objects.get(id = database_id)
        except BlastDb.DoesNotExist:
            return Response(f"Resource does not exist", status = status.HTTP_404_NOT_FOUND)
    
        try:
            local_db_folder = os.path.abspath(f'./fishdb/{database_id}/')
            if len(database_id) > 0 and os.path.exists(local_db_folder):
                shutil.rmtree(local_db_folder, ignore_errors=True)
        except BaseException:
            return Response(status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        return self.destroy(request, *args, **kwargs)

class BlastRunList(mixins.ListModelMixin, mixins.CreateModelMixin, generics.GenericAPIView):
    '''
    List all blast runs
    '''

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [IsAdminUser]

    @swagger_auto_schema(
        operation_summary='List all runs.',
        operation_description='Return a list of runs/jobs/queries saved by the API, including queued, running, and completed jobs.',
        tags = [tag_admin],
        responses={
            '200': openapi.Response(
                description='A list of all jobs saved.',
                schema = BlastRunSerializer(many=True),
                examples={
                    'application/json': [ BlastRunSerializer.Meta.example ]
                }
            )
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        List all blast databases
        '''
        return self.list(request, *args, **kwargs)

class BlastRunRun(generics.CreateAPIView):
    '''
    Run a blast query
    '''
    serializer_class = BlastRunRunSerializer
    permission_classes = [AllowAny]
    queryset = BlastRun.objects.all()
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    
    @swagger_auto_schema(
        operation_summary='Submit a run.',
        operation_description='Submit query sequence(s) and associated data to begin a run. Sequences can be supplied as either raw text in "query_sequence" or in an uploaded file named "query_file". Unless the submission is erroneous, the response will contain a unique run ID used to keep track of the status (the run is run asynchronously) and look up results when complete.\n\nA run will perform a BLASTN query of each query sequence against the sequences found within the BLAST database. Based on the values of the "create_hit_tree" and "create_db_tree" parameters, then up to 2 multiple alignment jobs will be performed using the query sequences and sequences from the hits or entire database, respectively.',
        # TODO: Add request / parameter schema

        request_body=openapi.Schema(
            type=openapi.TYPE_OBJECT,
            title='Parameters',
            properties={
                'job_name': openapi.Schema(
                    type=openapi.TYPE_STRING,
                    example='two sequences'
                ),
                'query_sequence': openapi.Schema(
                    type=openapi.TYPE_STRING,
                    example='>MG653404.1 Compsaraia iara isolate SA217 cytochrome c oxidase subunit 1 (COI) gene, partial cds; mitochondrial\nCCAACCAGGCGCCCTCCTGGGAGACGACCAAATTTACAATGTGGTCGTTACCGCCCATGCCTTCGTAATAATTTTCTTTATAGTAATGCCAATTATAATCGGAGGCTTTGGCAATTGACTTATCCCCTTAATAATTGCCGCGCCCGATATGGCATTCCCACGAATAAATAATATAAGCTTCTGACTGCTTCCCCCATCATTCTTCCTCCTACTTGCCTCTGCCGGGTTAGAGGCCGGAGTCGGGACAGGCTGAACGCTTTACCCCCCTCTTGCCGGTAATGCAGCACACGCTGGAGCCTCTGTAGACCTAACCATTTTCTCCCTTCACTTGGCCGGTGTCTCATCTATCCTCGGATCTATTAACTTTATCACTACAATTATTAATATGAAACCCCCAACAATATCCCAATACCAGCTTCCATTATTTATTTGATCCTTACTAGTAACCACAGTACTTCTACTACTCTCCCTTCCTGTTCTAGCTGCTGGA\n>MK401918.1 Porotergus gymnotus isolate ANSP189647 cytochrome oxidase subunit I (COI) gene, partial cds; mitochondrial\nGGCACACTATACATGGTGTTNGGCGCCTGGGCGGGTATAATTGGTACTGCTCTTANGCTTCTAATCCGGGCCGAGCTAAATCAACCGGGCACCCTCCTAGAAGACGACCAAATTTACAATGTGGCCGTCACTGCCCATGCCTTTGTAATAATTTTCTTTATAGTTATGCCAATCATAATTGGAGGCTTTGGCAATTGGCTTATCCCCCTAATAATTGCCGCGCCAGACATAGCATTCCCACGAATAAATAACATAAGCTTCTGGCTACTCCCCCCATCATTCTTCCTGCTCCTCGCCTCTGCTGGCTTAGAGGCCGGAGTTGGAACAGGTTGGACCCTATACCCCCCTCTTGCCGCTAATGCAGCACACGCCGGAGCTTCCGTAGACCTAACTATCTTCTCCCTTCACTTGGCGGGTGTTTCATCCATCCTCGGCTCCATTAACTTTATTACTACAATTATTAATATAAAACCTCCAACAATATCCCAATACCAACTCCCACTCTTTATCTGGTCCCTGCTGGTTACTACCGTGCTTCTACTACTCTCCCTTCCTGTCCTAGCTGCTGGTATTACCATGCTACTCACAGACCGAAATCTAAACACAGCATTCTTCGACCCAACGGGCGGCGGTGACCCAATTCTGTACCAACACTTGTTCTGGTT'
                ),
                'query_file': openapi.Schema(
                    type=openapi.TYPE_FILE
                ),
                'create_hit_tree': openapi.Schema(
                    type=openapi.TYPE_BOOLEAN,
                    example=True
                ),
                'create_db_tree': openapi.Schema(
                    type=openapi.TYPE_BOOLEAN,
                    example=False
                )
            }),
        tags = [tag_runs],
        responses={
            '400': openapi.Response(
                description='Bad request parameters. An accompanying message may specify the error with the request.',
                schema=openapi.Schema(
                    type=openapi.TYPE_OBJECT,
                    properties={
                        'message': openapi.Schema(
                            type=openapi.TYPE_STRING,
                            description='Helper message clarifying the cause of the error.'
                        )
                    }
                ),
                examples={
                    'application/json': {
                        'message': 'Request did not have raw sequence text or file upload specified to run with.'
                    }
                }
                ),
            '404': 'The BLAST database specified by the ID does not exist.',
            '201': openapi.Response(
                description='Run was successfully added to queue and a new unique run identifier is returned. Users should now use the given ID to continually check the status of the run to check whether it has completed.',
                schema=BlastRunSerializer(),
                examples={'application/json': BlastRunSerializer.Meta.example}
            ),
            '500': 'Unexpected error. No new run was created.'
        }
    )
    def post(self, request, *args, **kwargs):      
        '''
        Run blast query
        '''

        db_used = kwargs['pk']

        # Bad request if no query_sequence is given
        if not 'query_sequence' in request.data and not 'query_file' in request.FILES:
            return Response({
                'message': 'Request did not have raw sequence text or file upload specified to run with.'
            }, status = status.HTTP_400_BAD_REQUEST)

        # Bad request if database could not be found
        try:
            odb = BlastDb.objects.get(id = db_used)
        except BlastDb.DoesNotExist:
            return Response({
                'message': f"The database specified (id = {db_used}) does not exist."
            }, status = status.HTTP_400_BAD_REQUEST)

        # Bad request if the database does not have the minimum number of sequences
        MINIMUM_NUMBER_OF_SEQUENCES = 1
        sequences = NuccoreSequence.objects.filter(owner_database = odb.id)
        num_sequences_found = len(sequences)
        if num_sequences_found < MINIMUM_NUMBER_OF_SEQUENCES:
            return Response({
                'message': f"Cannot begin a blastn query on a database of less than {MINIMUM_NUMBER_OF_SEQUENCES} sequences. This database only contains {num_sequences_found} sequences."
            }, status = status.HTTP_400_BAD_REQUEST)

        # Bad request if create_hit_tree and create_db_tree are not boolean values
        create_hit_tree = False
        create_db_tree = False 
        if 'create_hit_tree' in request.data:
            val = request.data['create_hit_tree']
            if val is bool:
                create_hit_tree = val
            elif isinstance(val, str):
                create_hit_tree = (val == 'true')
            else:
                return Response({'message': 'Could not parse parameters for create_hit_tree'},
                status=status.HTTP_400_BAD_REQUEST)
        if 'create_db_tree' in request.data:
            val = request.data['create_db_tree']
            if val is bool:
                create_db_tree = val 
            elif isinstance(val, str):
                create_db_tree = (val == 'true')
            else:
                return Response({'message': 'Could not parse parameters for create_db_tree'},
                status=status.HTTP_400_BAD_REQUEST)

        query_sequences : List = []
        # obtain the query sequence, either from raw text or from the file upload
        # take precedence for the raw text
        if 'query_sequence' in request.data and len(request.data['query_sequence']) > 0:
            full_query : str = request.data['query_sequence']
            # try parsing with fasta first
            with io.StringIO(full_query) as query_string_io:
                query_records = SeqIO.parse(query_string_io, 'fasta')
                query_record : SeqIO.SeqRecord
                
                for query_record in query_records:
                    seq = str(query_record.seq).strip()
                    if not verify_dna(seq):
                        query_string_io.close()
                        return Response({
                            'message': f"The sequence provided for '>{query_record.description}' has non-DNA nucleotide characters."
                        }, status = status.HTTP_400_BAD_REQUEST) 
                    query_sequences.append({
                        'definition': query_record.description,
                        'query_sequence': seq
                    })

            # if we found no fasta entries, then parse the whole thing as a sequence
            if len(query_sequences) == 0:
                # replace all whitespace
                full_query = full_query.strip().replace('\r', '').replace('\n', '').replace(' ', '')
                if not verify_dna(full_query):
                    query_string_io.close()
                    return Response({
                        'message': f"Tried to parse input as single entry by removing whitespace, but the sequence provided has non-nucleotide characters."
                    }, status = status.HTTP_400_BAD_REQUEST) 
                query_sequences.append({
                    'definition': 'query_sequence', 
                    'query_sequence': full_query
                })
        else:
            query_file = request.FILES.get('query_file', False)
            # return an error if no file of name 'query_file' was found in the request
            if not query_file:
                return Response({'message': 'No submitted sequences were found. Either upload a .fasta file or paste raw text for a single sequence.'}, status = status.HTTP_400_BAD_REQUEST)

            # TODO: Check file type uploaded

            # parse the file with Biopython
            try:
                query_file_io = io.StringIO(query_file.file.read().decode('UTF-8'))
                query_file_seqs = SeqIO.parse(query_file_io, 'fasta')
            except BaseException:
                return Response({'message': 'The fasta file received was unable to be parsed with Biopython SeqIO using the "fasta" format. Double check the fasta contents.'}, status = status.HTTP_400_BAD_REQUEST)

            parsed_entry: SeqIO.SeqRecord
            for parsed_entry in query_file_seqs:
                seq = str(parsed_entry.seq)
                if not verify_dna(seq):
                    return Response({
                        'message': f"The sequence provided for '>{parsed_entry.description}' has non-DNA nucleotide characters."
                    }, status = status.HTTP_400_BAD_REQUEST) 
                query_sequences.append({
                    'definition': parsed_entry.description,
                    'query_sequence': seq
                })

        if len(query_sequences) == 0:
            return Response({'message': 'Found zero sequences in the request to query with. Sequences can be sent as text in a form or a file.'}, status = status.HTTP_400_BAD_REQUEST)

        print(query_sequences)

        job_name = request.data['job_name'] if 'job_name' in request.data else ''

        # path to parent folder of blast tools and django app
        project_root = os.path.abspath(os.path.dirname('./'))
        # path to ncbi.../bin/
        ncbi_blast_version = 'ncbi-blast-2.12.0+'
        blast_root =  f'{project_root}/{ncbi_blast_version}/bin'

        # path to output the database
        fishdb_root = project_root + '/fishdb'
        fishdb_path = fishdb_root + '/' + str(odb.id)
        if not os.path.exists(fishdb_path):
            os.mkdir(fishdb_path)

        # path to output the results
        results_uuid = uuid.uuid4()
        results_root = project_root + '/runs'
        results_path = results_root + '/' + str(results_uuid)
        if not os.path.exists(results_path):
            os.mkdir(results_path)
            # ensure that celery worker has access to file to write results.txt later
            shutil.chown(results_path, group='celery') 
            os.chmod(results_path, 0o774)

        # static path to download results
        static_path = '/var/www/runs/' + str(results_uuid)
        if not os.path.exists(static_path):
            os.mkdir(static_path)
            shutil.chown(static_path, group='celery')
            os.chmod(static_path, 0o775)

        print("Run id: " + str(results_uuid))

        if not odb.locked:
            print('Database was not locked. Creating database locally ...')
            
            # Make the database from scratch
            if os.path.exists(fishdb_path):
                try:
                    # delete the old folder
                    shutil.rmtree(fishdb_path, ignore_errors = False)
                except BaseException as base_exception:
                    return Response({
                        'message': "Server errored making the database.",
                        'error_type': type(base_exception)
                    }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
            
            print('Gathering sequences into fasta ...')
            os.mkdir(fishdb_path)
            fasta_file = fishdb_path + f'/database.fasta'
            with open(fasta_file, 'w') as my_file:
                for x in sequences:
                    identifier = x.accession_number
                    dna_sequence = x.dna_sequence
                    my_file.write('>' + identifier + '\n' + dna_sequence + '\n')
            
            my_file.close()

            print('Creating db now ...')
            command = '{} -in {} -dbtype nucl -out {} -title {} -parse_seqids'.format(blast_root + '/makeblastdb', fasta_file, fishdb_path + '/database', 'database')
        
            # Lock the database
            odb.locked = True
            odb.save()

            os.system(command)
        else:
            print('Database was locked. Will use existing local database.')

        # Perform blast search
        print('Generating query fasta file ...')
        query_file = results_path + '/query.fasta'
        with open(query_file, 'w') as tmp:
            for query_entry in query_sequences:
                tmp.write(f'>{query_entry["definition"]}\n')
                tmp.write(f'{query_entry["query_sequence"]}\n')

        tmp.close()

        print('Running BLAST search ...')

        # TODO: remove query_sequence reference
        run_details = BlastRun(id = results_uuid, db_used = odb, job_name = job_name, blast_version = ncbi_blast_version, errors = '', job_status=BlastRun.JobStatus.QUEUED, job_start_time = None, job_end_time = None, job_error_time = None, create_hit_tree = create_hit_tree, create_db_tree = create_db_tree)
        run_details.save()

        # make query sequence objects
        all_query_sequences = [BlastQuerySequence(**query_sequence, owner_run = run_details) for query_sequence in query_sequences]
        BlastQuerySequence.objects.bulk_create(all_query_sequences)

        run_details_id = run_details.id

        message_properties = {
            'MessageGroupId': str(run_details_id),
            'MessageDeduplicationId': str(run_details_id),
        }

        run_id_str = str(run_details_id)

        run_blast_command.apply_async(
            queue='BarcodeQueue.fifo',
            **message_properties, 
            kwargs={
                'ncbi_blast_version': ncbi_blast_version,
                'fishdb_id': str(odb.id),
                'run_id': run_id_str
            }, 
            link=signature('barcode_blastn.tasks.performAlignment', 
                kwargs={
                    'run_id': run_id_str
                }
            )  # type: ignore
        )      # type: ignore

        # create response 
        serializer = BlastRunSerializer(run_details)
        return Response(serializer.data, status=status.HTTP_201_CREATED)

class BlastRunDetail(mixins.DestroyModelMixin, generics.GenericAPIView):
    '''
    View the results of any given run
    '''

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [IsAdminOrReadOnly]
    
    @swagger_auto_schema(
        operation_summary='Get run results.',
        operation_description='Get run information and results of a BLAST run. The job_status indicates the status of the run. The results returned are complete once the job_status is "FIN" (i.e. run is finished). The run information includes the BLASTN hits and multiple alignment trees.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Run information successfully retrieved.',
                schema=BlastRunSerializer(),
                examples={
                    'application/json': BlastRunSerializer.Meta.example
                }
            ),
            '404': 'Run with the ID does not exist',
        }
    )
    def get(self, request, *args, **kwargs):
        db_primary_key = kwargs['pk']

        try:
            run = BlastRun.objects.get(id = db_primary_key)
        except BlastRun.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        
        serializer = BlastRunSerializer(run)
        response = Response(serializer.data, status=status.HTTP_200_OK)

        return response

class BlastRunInputDownload(generics.GenericAPIView):
    '''
    Return the sequences submitted to a run as a .fasta file
    '''
    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [IsAdminOrReadOnly]
    renderer_classes = [BlastRunFastaRenderer]
    
    @swagger_auto_schema(
        operation_summary='Get query sequence fasta file.',
        operation_description='Returns the original query sequences from the run submission, formatted as a FASTA file attachment.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Query sequences retrieved successfully and a fasta file attachment is returned.',
                schema=openapi.Schema(
                    type=openapi.TYPE_STRING,
                    format='binary'
                ),
            ),
            '404': 'Run data corresponding to the specified ID does not exist'
        }
    )
    def get(self, request, *args, **kwargs):
        db_primary_key = kwargs['pk']

        try:
            run = BlastRun.objects.get(id = db_primary_key)
        except BlastRun.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        
        file_name = 'query.fasta'
        query_path = os.path.abspath(f'./runs/{run.id}/{file_name}')
        if os.path.exists(query_path) and os.path.isfile(query_path):
            # return file in response
            file_handler = open(query_path, 'r')
            response = Response(file_handler, content_type="text/x-fasta", status=status.HTTP_200_OK)
            response['Content-Disposition'] = f'attachment; filename="{file_name}";'
            return response
        else:
            # return 404 error if the query.fasta file does not exist
            return Response(status=status.HTTP_404_NOT_FOUND)

class BlastRunDetailDownload(generics.GenericAPIView):
    '''
    Download the run results in either .txt or .csv format
    '''
    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [IsAdminOrReadOnly]
    renderer_classes = [BlastRunCSVRenderer, BlastRunTxtRenderer]
    
    @swagger_auto_schema(
        operation_summary='Get BLASTN results as a file',
        operation_description='Returns BLASTN results as a file attachment. The file format is .txt if Accept header specifies text/plain, and .csv if it is text/csv',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Results retrieved successfully and a file attachment is returned.',
                schema = openapi.Schema(
                    type= openapi.TYPE_STRING,
                    format='binary'
                ),
            ),
            '404': 'Run data corresponding to the specified ID does not exist'
        },
    )
    def get(self, request, *args, **kwargs):
        db_primary_key = kwargs['pk']

        try:
            run = BlastRun.objects.get(id = db_primary_key)
        except BlastRun.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        
        serializer = BlastRunSerializer(run)
       
        response = Response(serializer.data, status=status.HTTP_200_OK)
        file_name_root = serializer.data['job_name']
        if file_name_root is None or len(file_name_root) == 0:
            file_name_root = 'results'       
        
        if request.accepted_media_type.startswith('text/csv'):
            response['Content-Disposition'] = f'attachment; filename="{file_name_root}.csv";'
        elif request.accepted_media_type.startswith('text/plain'):
            response['Content-Disposition'] = f'attachment; filename="{file_name_root}.txt";'

        return response

class BlastRunStatus(generics.RetrieveAPIView):
    '''
    Check the status of a run
    '''
    queryset = BlastRun.objects.all()
    serializer_class = BlastRunStatusSerializer
    permission_classes = [IsAdminOrReadOnly]

    @swagger_auto_schema(
        operation_summary='Get status of run',
        operation_description='Returns a minimal set of information useful for polling/checking the status of run.\n\nHow to interpret job_status:\n"QUE" (Queued): The run is currently waiting in the queue for its turn to be processed.\n"STA" (Started): The run is currently being processed.\n"FIN" (Finished): The run successfully completed and complete results are now visible.\n"ERR" (Errored): The run encountered an unexpected error and terminated.\n"UNK" (Unknown): The status is unknown, likely because there was an unexpected database or server error. Please submit another run.\n"DEN" (Denied): The run submission was received by the server, but it was denied from being processed.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Results retrieved successfully and a file attachment is returned.',
                schema=BlastRunStatusSerializer(),
                examples={
                    'application/json': BlastRunStatusSerializer.Meta.example
                }
            ),
            '404': 'Run data corresponding to the specified ID does not exist'
        },
    )
    def get(self, request, *args, **kwargs):
        return super().get(request, *args, **kwargs)
