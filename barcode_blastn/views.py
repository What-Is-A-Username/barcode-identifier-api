# TODO: Install Python 3.10.6 to match EC2

import io
from celery import chain, signature
from barcode_identifier_api.celery import app
from typing import Dict, List
from urllib import request
import uuid
import os
import shutil
from Bio import SeqIO
from barcode_blastn.helper.parse_gb import retrieve_gb
from barcode_blastn.helper.verify_query import verify_query, verify_dna
from barcode_blastn.helper.read_tree import readDbTreeFromFile, readHitTreeFromFile
from barcode_blastn.models import BlastQuerySequence, BlastRun, Hit, NuccoreSequence, BlastDb
from barcode_blastn.permissions import IsAdminOrReadOnly
from barcode_blastn.renderers import BlastRunCSVRenderer, BlastRunTxtRenderer, BlastDbCSVRenderer, BlastDbFastaRenderer, BlastRunFastaRenderer
from barcode_blastn.serializers import BlastRunRunSerializer, BlastRunSerializer, BlastRunStatusSerializer, HitSerializer, NuccoreSequenceAddSerializer, NuccoreSequenceBulkAddSerializer, NuccoreSequenceSerializer, BlastDbSerializer, BlastDbListSerializer
from rest_framework import status, generics, mixins
from rest_framework.parsers import JSONParser, FormParser, MultiPartParser
from rest_framework.permissions import IsAdminUser, AllowAny
from rest_framework.renderers import JSONRenderer, BrowsableAPIRenderer, TemplateHTMLRenderer
from rest_framework.response import Response
from urllib.error import HTTPError
from django.db.models import QuerySet

from barcode_blastn.tasks import run_blast_command, performAlignment
# from barcode_tree.tasks import performAlignment

'''
List all the sequences in the server, irrespective of database
'''
class NuccoreSequenceList(mixins.ListModelMixin, generics.GenericAPIView):
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer
    permission_classes = [AllowAny]

    '''
    List all accession numbers saved to all databases
    '''
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

'''
Bulk add sequences to a specified database
'''
class NuccoreSequenceBulkAdd(generics.CreateAPIView):
    serializer_class = NuccoreSequenceBulkAddSerializer
    permission_classes = [ IsAdminUser ]
    queryset = NuccoreSequence.objects.all()
    
    '''
    Create a new accession number and add it to an existing database
    '''
    def post(self, request, *args, **kwargs):
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
                status=status.HTTP_400_BAD_REQUEST)
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
    
    '''
    Create a new accession number and add it to an existing database
    '''
    def post(self, request, *args, **kwargs):
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
            return Response({'message': f"The GenBank database did not return an entry for the accession number {accession_number}"}, status=status.HTTP_400_BAD_REQUEST)
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

        return Response(serializer.errors, status = status.HTTP_400_BAD_REQUEST)

'''
Get or delete a sequence specified by ID
'''
class NuccoreSequenceDetail(mixins.DestroyModelMixin, generics.RetrieveAPIView):
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer
    permission_classes = [IsAdminOrReadOnly]

    '''
    Delete an accession number from the database
    '''
    def delete(self, request, *args, **kwargs):
        return self.destroy(request, *args, **kwargs)

'''
Get a sequence by accession number
'''
class NuccoreSequenceDetailByAccession(mixins.RetrieveModelMixin, generics.GenericAPIView):
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer
    permission_classes = [IsAdminOrReadOnly]

    '''
    Retrieve all sequences with an accession number specified in the list
    '''
    def get(self, request, *args, **kwargs):
        accession_number = kwargs['an'].split(',')

        try:
            seq = NuccoreSequence.objects.filter(accession_number__in=accession_number)
        except NuccoreSequence.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)

        serializer = NuccoreSequenceSerializer(seq, many=True)
        response = Response(serializer.data, status=status.HTTP_200_OK)
        return response

'''
Return a list of all blast databases
'''
class BlastDbList(mixins.ListModelMixin, mixins.CreateModelMixin, generics.GenericAPIView):
    queryset = BlastDb.objects.all()
    serializer_class = BlastDbListSerializer
    permission_classes = [IsAdminOrReadOnly]

    '''
    List all blast databases
    '''
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

    '''
    Create a new blast database.
    '''
    def post(self, request, *args, **kwargs):
        return self.create(request, *args, **kwargs)

'''
Retrieve the details of a given blast database
'''
class BlastDbDetail(mixins.RetrieveModelMixin, mixins.UpdateModelMixin, mixins.DestroyModelMixin, generics.GenericAPIView):
    queryset = BlastDb.objects.all()
    serializer_class = BlastDbSerializer
    permission_classes = [IsAdminOrReadOnly]
    renderer_classes = [JSONRenderer, BrowsableAPIRenderer, TemplateHTMLRenderer, BlastDbCSVRenderer, BlastDbFastaRenderer]

    '''
    Retrieve data in the blastdb.
    '''
    def get(self, request, *args, **kwargs):

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

    '''
    Partially update this blastdb.
    '''
    def patch(self, request, *args, **kwargs):
        try:
            db = BlastDb.objects.get(id = kwargs['pk'])
        except BlastDb.DoesNotExist:
            return Response("Resource does not exist", status = status.HTTP_404_NOT_FOUND)
        
        if db.locked and (not 'locked' in request.data or request.data['locked'] == db.locked):
            return Response("This entry is locked and cannot be patched.", status = status.HTTP_400_BAD_REQUEST)

        return self.partial_update(request, *args, **kwargs)

    '''
    Delete the blastdb
    '''
    def delete(self, request, *args, **kwargs):
        try:
            database_id = str(kwargs['pk'])
            db = BlastDb.objects.get(id = database_id)
        except BlastDb.DoesNotExist:
            return Response(f"Resource does not exist", status = status.HTTP_404_NOT_FOUND)

        local_db_folder = os.path.abspath(f'./fishdb/{database_id}/')
        if len(database_id) > 0 and os.path.exists(local_db_folder):
            shutil.rmtree(local_db_folder, ignore_errors=True)

        return self.destroy(request, *args, **kwargs)

'''
List all blast runs
'''
class BlastRunList(mixins.ListModelMixin, mixins.CreateModelMixin, generics.GenericAPIView):

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [IsAdminUser]

    '''
    List all blast databases
    '''
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

'''
Run a blast query
'''
class BlastRunRun(mixins.CreateModelMixin, generics.GenericAPIView):
    serializer_class = BlastRunRunSerializer
    permission_classes = [AllowAny]
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    
    '''
    Run blast query
    '''
    def post(self, request, *args, **kwargs):      

        db_used = kwargs['pk']

        # Bad request if no query_sequence is given
        if not 'query_sequence' in request.data and not 'query_file' in request.FILES:
            return Response({
                'message': 'Request did not have raw sequence text or file upload specified to blast with.'
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
                    })
            
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

        # django_rq.enqueue(run_blast_command, blast_root=blast_root, fishdb_path=fishdb_path, query_file=query_file, run_details=run_details, results_path=results_path)

        # create response 
        serializer = BlastRunSerializer(run_details)
        return Response(serializer.data, status=status.HTTP_201_CREATED)

'''
View the results of any given run
'''
class BlastRunDetail(mixins.DestroyModelMixin, generics.GenericAPIView):

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [IsAdminOrReadOnly]
    renderer_classes = [JSONRenderer, BrowsableAPIRenderer, BlastRunTxtRenderer]
    
    def get(self, request, *args, **kwargs):
        db_primary_key = kwargs['pk']

        try:
            run = BlastRun.objects.get(id = db_primary_key)
        except BlastRun.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        
        serializer = BlastRunSerializer(run)
        response = Response(serializer.data, status=status.HTTP_200_OK)

        return response

    '''
    Delete an accession number from the database
    '''
    def delete(self, request, *args, **kwargs):
        run_id = str(kwargs['pk'])
        local_run_folder = os.path.abspath(f'./runs/{run_id}/')
        if len(run_id) > 0 and os.path.exists(local_run_folder):
            shutil.rmtree(local_run_folder, ignore_errors=True)

        return self.destroy(request, *args, **kwargs)

'''
Return the sequences submitted to a run as a .fasta file
'''
class BlastRunInputDownload(generics.GenericAPIView):
    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [IsAdminOrReadOnly]
    renderer_classes = [BlastRunFastaRenderer]
    
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

'''
Download the run results in either .txt or .csv format
'''
class BlastRunDetailDownload(generics.GenericAPIView):
    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [IsAdminOrReadOnly]
    renderer_classes = [BlastRunCSVRenderer, BlastRunTxtRenderer]
    
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

'''
Access a particular hit
'''
class HitDetail(mixins.ListModelMixin, generics.GenericAPIView):
    queryset = Hit.objects.all()
    serializer_class = HitSerializer
    permission_classes = [IsAdminOrReadOnly]

    '''
    List all accession numbers saved to all databases
    '''
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

'''
Check the status of a run
'''
class BlastRunStatus(generics.RetrieveAPIView):
    queryset = BlastRun.objects.all()
    serializer_class = BlastRunStatusSerializer
    permission_classes = [IsAdminOrReadOnly]
