import os
import re
import json
import shutil
import subprocess
from django.http.response import FileResponse, HttpResponse
from barcode_blastn.helper.parse_gb import InvalidAccessionNumberError, parse_gbx_xml, retrieve_gb
from rest_framework.response import Response
from barcode_blastn.models import BlastRun, NuccoreSequence, BlastDb
from barcode_blastn.serializers import BlastRunSerializer, NuccoreSequenceSerializer, BlastDbSerializer
from rest_framework import status, generics, mixins
from urllib.error import HTTPError
from rest_framework.parsers import MultiPartParser
from django.template import loader

# TODO: Move GET and POST to /blastdbs/<pk>/add
class NuccoreSequenceList(mixins.ListModelMixin, mixins.CreateModelMixin, generics.GenericAPIView):
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer

    '''
    List all accession numbers saved to all databases
    '''
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

    '''
    Create a new accession number and add it to an existing database
    '''
    def post(self, request, *args, **kwargs):
        currentData = {}

        accession_number = request.data['accession_number']

        # TODO: Check for duplicates in the SAME db 
        try:
            duplicate = NuccoreSequence.objects.get(accession_number = accession_number)
        except NuccoreSequence.DoesNotExist:
            pass
        else:
            return Response({'message': "Entry already exists for the accession number specified.", 'accession_number': accession_number}, status.HTTP_400_BAD_REQUEST)
        
        try:
            xml_string = retrieve_gb(accession_number = accession_number)
        except HTTPError:
            return Response({'message': f"The GenBank database does not contain an entry for this accession number.", 'accession_number': accession_number}, status=status.HTTP_400_BAD_REQUEST)
        except BaseException:
            return Response({'message': f"Failed to parse retrieved GenBank database entry.", 'accession_number': accession_number}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        try:
            currentData = parse_gbx_xml(xml_string = xml_string, accession_number = accession_number)
        except BaseException:
            return Response({'message': f"The GenBank entry for accession number could not be parsed properly.", 'accession_number': accession_number}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        all_data = request.POST.copy()
        all_data.update(currentData)

        serializer = NuccoreSequenceSerializer(data = all_data)
        if serializer.is_valid():
            serializer.save()
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
        
        if db.locked:
            return Response("This entry is locked and cannot be patched.", status = status.HTTP_400_BAD_REQUEST)

        return self.partial_update(request, *args, **kwargs)

    '''
    Delete the blastdb
    '''
    def delete(self, request, *args, **kwargs):
        # TODO: If a local blastdb file has been made, delete it
        return self.destroy(request, *args, **kwargs)

class BlastRunList(mixins.ListModelMixin, mixins.CreateModelMixin, generics.GenericAPIView):

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer

    '''
    List all blast databases
    '''
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

    '''
    Run blast query
    '''
    def post(self, request, *args, **kwargs):      
        # Bad request if no query_sequence is given
        if not 'query_sequence' in request.data:
            return Response({
                'message': 'No nucleotide query sequence given to blast with.'
            }, status = status.HTTP_400_BAD_REQUEST)
        query_sequence = request.data['query_sequence']

        job_name = request.data['job_name'] if 'job_name' in request.data else ''

        # Bad request if no database is given
        if not 'database' in request.data:
            return Response({
                'message': "No blast database specified."
            }, status = status.HTTP_400_BAD_REQUEST)

        # Bad request if database could not be found
        try:
            odb = BlastDb.objects.get(id = request.data['database'])
        except BlastDb.DoesNotExist:
            return Response({
                'message': "The database specified does not exist.",
                'database': request.data['database']
            }, status = status.HTTP_400_BAD_REQUEST)

        # Bad request if the database does not have the minimum number of sequences
        MINIMUM_NUMBER_OF_SEQUENCES = 1
        sequences = NuccoreSequence.objects.filter(owner_database = odb.id)
        if len(sequences) < MINIMUM_NUMBER_OF_SEQUENCES:
            return Response({
                'message': "Cannot begin a blastn query on a database of less than " + MINIMUM_NUMBER_OF_SEQUENCES + " sequences.",
                'current_size': len(sequences),
                'min_size': MINIMUM_NUMBER_OF_SEQUENCES
            }, status = status.HTTP_400_BAD_REQUEST)

        # path to parent folder of blast tools and django app
        project_root = os.path.abspath(os.path.dirname('./'))
        # path to ncbi.../bin/
        blast_root = project_root + '/ncbi-blast-2.12.0+/bin'

        fishdb_path = project_root + '/fishdb'

        # TODO: Generate unique name for each database
        db_name = 'fish'

        if not odb.locked:
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
            fasta_file = fishdb_path + '/{}.fasta'.format(db_name)
            with open(fasta_file, 'w') as my_file:
                for x in sequences:
                    print(x.organism)
                    identifier = '_'.join(x.organism.split(' '))
                    dna_sequence = x.dna_sequence
                    my_file.write('>' + identifier + '\n' + dna_sequence)
            my_file.close()

            print('Creating db now ...')
            command = '{} -in {} -dbtype nucl -out {} -title {}'.format(blast_root + '/makeblastdb', fasta_file, fishdb_path + '/' + db_name, db_name)

            # TODO: Lock the blastdb once the db has been created

            os.system(command)

        # Perform blast search
        print('Generating query fasta file ...')
        with open(fishdb_path + '/query.fasta', 'w') as tmp:
            tmp.write(query_sequence)
        tmp.close()

        print('Running BLAST search ...')
        blast_command = '{} -db {} -outfmt 7 -query {}'.format(blast_root + '/blastn', fishdb_path + '/' + db_name, fishdb_path + '/query.fasta')

        process = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        print('BLAST search completed.')

        out, err = process.stdout.read().decode(), process.stderr.read().decode()

        print('Writing results to file ...')

        with open(fishdb_path + '/results.txt', "w") as results_file:
            results_file.write(out + '\n')
            if err:
                results_file.write('Errors occurred when processing the query:\n')
                results_file.write(err)
        results_file.close()

        outlines = out.split('\n')
        blast_version = outlines[0].replace('# ', '')
        errors = err 

        run_details = BlastRun(blastdb = odb, job_name = job_name, blast_version = blast_version, errors = errors)

        run_details.save()

        return Response("Successfully ran blast query.", status=status.HTTP_201_CREATED)

class BlastRunDetail(generics.RetrieveAPIView):

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer

    # TODO: Return blast results as JSON, HTML text, and .txt
    def get(self, request, *args, **kwargs):
        return self.get(request, *args, **kwargs)
