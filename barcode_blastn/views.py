import uuid
import os
import re
import json
import shutil
import django_rq
from django.http.response import FileResponse, HttpResponse
from barcode_blastn.helper.parse_gb import InvalidAccessionNumberError, parse_gbx_xml, retrieve_gb
from rest_framework.response import Response
from barcode_blastn.models import BlastRun, Hit, NuccoreSequence, BlastDb
from barcode_blastn.serializers import BlastRunRunSerializer, BlastRunSerializer, HitSerializer, NuccoreSequenceAddSerializer, NuccoreSequenceSerializer, BlastDbSerializer
from rest_framework import status, generics, mixins
from urllib.error import HTTPError
from rest_framework.parsers import MultiPartParser
from django.template import loader

from barcode_blastn.helper.run_blast import run_blast_command

class NuccoreSequenceList(mixins.ListModelMixin, generics.GenericAPIView):
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer

    '''
    List all accession numbers saved to all databases
    '''
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)
        
class NuccoreSequenceAdd(generics.CreateAPIView):
    serializer_class = NuccoreSequenceAddSerializer
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
            xml_string = retrieve_gb(accession_number = accession_number)
            currentData = parse_gbx_xml(xml_string = xml_string, accession_number = accession_number)
        except HTTPError:
            return Response({'message': f"The GenBank database does not contain an entry for this accession number.", 'accession_number': accession_number}, status=status.HTTP_400_BAD_REQUEST)
        except BaseException:
            return Response({'message': f"Failed to parse retrieved GenBank database entry.", 'accession_number': accession_number}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        
        currentData['owner_database'] = db.id
        all_data = request.POST.copy()
        all_data.update(currentData)

        serializer = NuccoreSequenceSerializer(data = all_data)
        if serializer.is_valid():
            serializer.save(owner_database = db)
            return Response(serializer.data, status = status.HTTP_201_CREATED)

        return Response(serializer.errors, status = status.HTTP_400_BAD_REQUEST)

class NuccoreSequenceDetail(mixins.DestroyModelMixin, generics.RetrieveAPIView):
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer

    '''
    Delete an accession number from the database
    '''
    def delete(self, request, *args, **kwargs):
        return self.destroy(request, *args, **kwargs)

# TODO: Remove uploading feature; one POST should not be able to create this many resources at once. Front end should handle upload of multiple accession numbers
class NuccoreSequenceListUpload(mixins.CreateModelMixin, generics.GenericAPIView):
    parser_classes = [MultiPartParser]

    def post(self, request, filename, format = None):

        uploaded_file = request.FILES['file'].file

        try:
            content = uploaded_file.read().decode('UTF-8')
        except BaseException as e:
            return Response({'message': 'Unexpected error reading the uploaded file.', 
            'error_type': type(e)}, status=status.HTTP_400_BAD_REQUEST)

        lines = re.split('\r\n|\n|\r', content)
        lines = [l.strip() for l in lines if len(l.strip()) > 0]

        return Response({'accession_numbers': json.dumps(lines)}, 
        status=status.HTTP_201_CREATED)

class BlastDbList(mixins.ListModelMixin, mixins.CreateModelMixin, generics.GenericAPIView):
    queryset = BlastDb.objects.all()
    serializer_class = BlastDbSerializer

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

class BlastDbDetail(mixins.RetrieveModelMixin, mixins.UpdateModelMixin, mixins.DestroyModelMixin, generics.GenericAPIView):
    queryset = BlastDb.objects.all()
    serializer_class = BlastDbSerializer

    '''
    Retrieve data in the blastdb.
    If get_format = 'text', return a .txt file of every accession number, line by line
    '''
    def get(self, request, *args, **kwargs):

        db_primary_key = kwargs['pk']

        try:
            db = BlastDb.objects.get(id = db_primary_key)
        except BlastDb.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)

        sequences = NuccoreSequence.objects.filter(owner_database = db_primary_key)
        accession_numbers = '\n'.join([s.accession_number for s in sequences])

        if 'get_format' in request.query_params:
            gformat = request.query_params['get_format']

            # Return downloadable text file
            if gformat == 'file':
                response = FileResponse(accession_numbers, content_type='text/plain')
                response['Content-Disposition'] = 'attachment; filename=db.txt;'
                return response
            # Return plain text html
            elif gformat == 'text':
                template = loader.get_template('blastdb.html')
                page = template.render(context = { 'list': accession_numbers})

                return HttpResponse(page, content_type='text/html', status=status.HTTP_200_OK)
        
        serializer = BlastDbSerializer(db)
        
        return Response(serializer.data, status=status.HTTP_200_OK)
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
        # TODO: Test this method
        try:
            database_id = str(kwargs['pk'])
            db = BlastDb.objects.get(id = database_id)
        except BlastDb.DoesNotExist:
            return Response("Resource does not exist", status = status.HTTP_404_NOT_FOUND)

        local_db_folder = os.path.abspath(f'./fishdb/{database_id}/')
        if len(database_id) > 0 and os.path.exists(local_db_folder):
            shutil.rmtree(local_db_folder, ignore_errors=True)

        return self.destroy(request, *args, **kwargs)

class BlastRunList(mixins.ListModelMixin, mixins.CreateModelMixin, generics.GenericAPIView):

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer

    '''
    List all blast databases
    '''
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

class BlastRunRun(mixins.CreateModelMixin, generics.GenericAPIView):
    serializer_class = BlastRunRunSerializer

    '''
    Run blast query
    '''
    def post(self, request, *args, **kwargs):      

        db_used = kwargs['pk']

        # Bad request if no query_sequence is given
        if not 'query_sequence' in request.data:
            return Response({
                'message': 'No nucleotide query sequence given to blast with.'
            }, status = status.HTTP_400_BAD_REQUEST)
        query_sequence = request.data['query_sequence']

        job_name = request.data['job_name'] if 'job_name' in request.data else ''

        # Bad request if database could not be found
        try:
            odb = BlastDb.objects.get(id = db_used)
        except BlastDb.DoesNotExist:
            return Response({
                'message': "The database specified does not exist.",
                'database': db_used
            }, status = status.HTTP_400_BAD_REQUEST)

        # Bad request if the database does not have the minimum number of sequences
        MINIMUM_NUMBER_OF_SEQUENCES = 1
        sequences = NuccoreSequence.objects.filter(owner_database = odb.id)
        if len(sequences) < MINIMUM_NUMBER_OF_SEQUENCES:
            return Response({
                'message': f"Cannot begin a blastn query on a database of less than {MINIMUM_NUMBER_OF_SEQUENCES} sequences.",
                'current_size': len(sequences),
                'min_size': MINIMUM_NUMBER_OF_SEQUENCES
            }, status = status.HTTP_400_BAD_REQUEST)

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

        if not odb.locked:
            print('Database was not locked. Creating database locally ...')
            # Make the database
            
            if os.path.exists(fishdb_path):
                try:
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
                    print(x.organism)
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
        # else
            print('Database was locked. Will use existing local database.')

        # Perform blast search
        print('Generating query fasta file ...')
        query_file = results_path + '/query.fasta'
        with open(query_file, 'w') as tmp:
            tmp.write('>query_sequence\n')
            tmp.write(query_sequence)
        tmp.close()

        print('Running BLAST search ...')

        run_details = BlastRun(id = results_uuid, db_used = odb, job_name = job_name, blast_version = ncbi_blast_version, errors = '', query_sequence = query_sequence, job_status=BlastRun.JobStatus.QUEUED, job_start_time = None, job_end_time = None)
        run_details.save()

        django_rq.enqueue(run_blast_command, blast_root=blast_root, fishdb_path=fishdb_path, query_file=query_file, run_details=run_details, results_path=results_path)

        # create response 
        # details = BlastRun.objects.get(pk=run_details.pk)
        serializer = BlastRunSerializer(run_details)
        return Response(serializer.data, status=status.HTTP_201_CREATED)


class BlastRunDetail(mixins.DestroyModelMixin, generics.GenericAPIView):

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    
    # TODO: Return blast results as JSON, HTML text, and .txt
    def get(self, request, *args, **kwargs):
        # self.get(request, *args, **kwargs)
        db_primary_key = kwargs['pk']

        try:
            run = BlastRun.objects.get(id = db_primary_key)
        except BlastRun.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        
        serializer = BlastRunSerializer(run)

        return(Response(serializer.data, status=status.HTTP_200_OK))

    '''
    Delete an accession number from the database
    '''
    def delete(self, request, *args, **kwargs):
        run_id = str(kwargs['pk'])
        print(run_id)
        local_run_folder = os.path.abspath(f'./runs/{run_id}/')
        if len(run_id) > 0 and os.path.exists(local_run_folder):
            shutil.rmtree(local_run_folder, ignore_errors=True)

        return self.destroy(request, *args, **kwargs)

class HitDetail(mixins.ListModelMixin, generics.GenericAPIView):
    queryset = Hit.objects.all()
    serializer_class = HitSerializer

    '''
    List all accession numbers saved to all databases
    '''
    def get(self, request, *args, **kwargs):
        pk = kwargs['pk']
        return self.list(request, *args, **kwargs)
