# TODO: Install Python 3.10.6 to match EC2
# TODO: Add endpoint to give create, edit and remove DatabaseShare permissions
# TODO: Guard endpoints to only allow signed in users

from datetime import datetime
import io
import os
import shutil
import uuid
from typing import Any, Dict, List
from urllib.error import HTTPError

from barcode_identifier_api.celery import app
from Bio import SeqIO
from celery import signature
from django.contrib.auth import login
from django.contrib.auth.models import User
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from knox.auth import TokenAuthentication
from knox.views import LoginView as KnoxLoginView
from knox.views import LogoutAllView as KnoxLogoutAllView
from knox.views import LogoutView as KnoxLogoutView
from rest_framework import generics, mixins, status
from rest_framework.authtoken.serializers import AuthTokenSerializer
from rest_framework.exceptions import PermissionDenied
from rest_framework.parsers import FormParser, JSONParser, MultiPartParser
from rest_framework.permissions import AllowAny, IsAuthenticated
from rest_framework.renderers import (BrowsableAPIRenderer, JSONRenderer,
                                      TemplateHTMLRenderer)
from rest_framework.response import Response

from barcode_blastn.controllers.blastdb_controller import (
    AccessionsAlreadyExist, AccessionsNotFound, DatabaseLocked, InsufficientAccessionData,
    add_sequences_to_database, create_blastdb, delete_blastdb, delete_library, delete_sequences_in_database,
    save_blastdb, update_sequences_in_database)
from barcode_blastn.file_paths import (get_data_fishdb_path, get_data_run_path,
                                       get_ncbi_folder, get_static_run_path)
from barcode_blastn.helper.parse_gb import (MAX_ACCESSIONS_PER_REQUEST,
                                            AccessionLimitExceeded,
                                            GenBankConnectionError)
from barcode_blastn.helper.verify_query import verify_dna
from barcode_blastn.models import (BlastDb, BlastQuerySequence, BlastRun,
                                   Library, NuccoreSequence)
from barcode_blastn.permissions import (BlastDbEndpointPermission, BlastRunEndpointPermission,
                                        DatabaseSharePermissions,
                                        LibraryEndpointPermission,
                                        LibrarySharePermissions,
                                        NuccoreSequenceEndpointPermission,
                                        NuccoreSharePermission)
from barcode_blastn.renderers import (BlastDbCSVRenderer, BlastDbCompatibleRenderer, BlastDbFastaRenderer,
                                      BlastRunCSVRenderer,
                                      BlastRunFastaRenderer,
                                      BlastRunTxtRenderer)
from barcode_blastn.serializers import (BlastDbCreateSerializer,
                                        BlastDbEditSerializer,
                                        BlastDbListSerializer,
                                        BlastDbSequenceEntrySerializer,
                                        BlastDbSerializer,
                                        BlastRunRunSerializer,
                                        BlastRunSerializer,
                                        BlastRunStatusSerializer,
                                        LibraryCreateSerializer,
                                        LibraryEditSerializer,
                                        LibraryListSerializer,
                                        LibrarySerializer,
                                        NuccoreSequenceBulkAddSerializer,
                                        NuccoreSequenceSerializer,
                                        UserSerializer)

tag_libraries = 'Libraries'
tag_runs = 'Runs/Jobs'
tag_blastdbs = 'BLAST Databases'
tag_sequences = 'GenBank Accessions'
tag_admin = 'Admin Tools'
tag_users = 'User Authentication'

class UserDetailView(generics.GenericAPIView):
    '''
    Return the user details associated with the authenticated
    user
    '''
    authentication_classes = (TokenAuthentication,)
    permission_classes = (IsAuthenticated,)

    @swagger_auto_schema(
        operation_summary='Get user details.',
        operation_description='Retrieve information about the current user signed-in through token authentication.',
        tags = [tag_users],
        responses={
            '200': openapi.Response(
                description='Successfully returned sequence information',
                schema=UserSerializer(read_only=True),
                examples={
                    'application/json': UserSerializer.Meta.example
                }
            ),
            '403': 'User could not be authenticated.'
        }
    )
    def get(self, request, format=None):
        user: User = request.user
        return Response(UserSerializer(request.user).data)

class LogoutAllView(KnoxLogoutAllView):
    @swagger_auto_schema(
        security=[{'Basic': []}, {'Bearer': []}],
        operation_summary='User logoff all',
        operation_description='Sign out of the user account by signing out of all tokens corresponding to the user specified by the provided authentication token. The token is provided in the request header with the format `Authorization: Bearer <token>`, where <token> is the full token.',
        tags = [tag_users],
    )
    def post(self, request, format=None):
        return super().post(request, format)

class LogoutView(KnoxLogoutView):
    '''
    Sign out of the signed in user
    '''
    @swagger_auto_schema(
        security=[{'Basic': []}, {'Bearer': []}],
        operation_summary='User logoff',
        operation_description='Sign out of the user account by signing out the provided authentication token. The token is provided in the request header with the format `Authorization: Bearer <token>`, where <token> is the full token.',
        tags = [tag_users],
    )
    def post(self, request, format=None):
        return super().post(request, format)

class LoginView(KnoxLoginView):
    '''
    Sign in the user and return the authentication token.
    '''
    permission_classes = (AllowAny,)

    def get_user_serializer_class(self):
        return UserSerializer

    @swagger_auto_schema(
        operation_summary='User login',
        operation_description='Sign in with the specified credentials and return the authentication token.',
        tags = [tag_users],
        request_body=openapi.Schema(
            type=openapi.TYPE_OBJECT,
            properties={
                'username': openapi.Schema(type=openapi.TYPE_STRING),
                'password': openapi.Schema(type=openapi.TYPE_STRING),
            }
        )
    )
    def post(self, request, format=None):
        serializer = AuthTokenSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        user = serializer.validated_data['user']  # type: ignore
        login(request, user)
        response = super(LoginView, self).post(request, format=None)

        return response

class NuccoreSequenceAdd(mixins.UpdateModelMixin, mixins.DestroyModelMixin, generics.CreateAPIView):   
    '''
    Bulk add sequences to a specified database
    '''
    serializer_class = NuccoreSequenceBulkAddSerializer
    permission_classes = [NuccoreSequenceEndpointPermission]
    queryset = NuccoreSequence.objects.all()
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    
    @swagger_auto_schema(
        operation_summary='Add accession numbers to database.',
        operation_description=f'From a list of accession numbers, add them to an existing database. List must contain between 1-{MAX_ACCESSIONS_PER_REQUEST} accession numbers inclusive.',
        tags = [tag_blastdbs, tag_sequences],
        manual_parameters=[openapi.Parameter(
            name='id',
            description='Id of BLAST database to add accession numbers to',
            in_=openapi.IN_PATH, 
            type=openapi.TYPE_STRING,
            format=openapi.FORMAT_UUID
        )],
        request_body=NuccoreSequenceBulkAddSerializer,
        responses={
            '201': openapi.Response(
                description='All accession numbers added to the database.',
                schema=BlastDbSequenceEntrySerializer(many=True),
                examples={
                    'application/json': [ BlastDbSequenceEntrySerializer.Meta.example ]
                }
            ),
            '400': 'Bad parameters in request. Example reasons: list may be empty or too long (>100); accession numbers may be invalid or duplicate; database may be locked',
            '403': 'Insufficient permissions.',
            '404': 'BLAST database matching the specified ID was not found.',
            '500': 'Unexpected error.',
            '502': 'Encountered error connecting to GenBank.'
        }
    )
    def post(self, request, *args, **kwargs):
        '''
        Create a new accession number and add it to an existing database
        '''
        pk = kwargs['pk']
        db: BlastDb
        # Query for the Blast database
        try:
            db = BlastDb.objects.get(id=pk)
        except BlastDb.DoesNotExist:
            return Response({'message': 'Database does not exist', 'requested': pk}, status.HTTP_404_NOT_FOUND)
        
        # Desired numbers collects the accession numbers to add to the database
        desired_numbers = []
        # Check if accession number file provided
        accession_file = request.FILES.get('accession_file', None)
        # return an error if no file of name 'query_file' was found in the request
        if not accession_file is None:
            query_file_io = io.StringIO(accession_file.file.read().decode('UTF-8'))
            file_accessions: List[str] = query_file_io.readlines()
            desired_numbers.extend(file_accessions)

        # Check if accession numbers provided
        serialized_data = NuccoreSequenceBulkAddSerializer(data=request.data)
        if serialized_data.is_valid():
            desired_numbers.extend(request.data.get('accession_numbers', []))

        # Respond with error if accessions empty
        if len(desired_numbers) == 0:
            return Response({'message': 'The list of numbers to add cannot be empty.', 'accession_numbers': str(desired_numbers)}, status.HTTP_400_BAD_REQUEST)

        # Deny user if user has insufficient permissions
        if not NuccoreSharePermission.has_add_permission(user=request.user, obj=None):
            return Response(status=status.HTTP_403_FORBIDDEN)
        
        # Respond with an error if database is locked
        if db.locked:
            return Response({'message': 'The database is locked and its accession numbers cannot be added, edited or removed.'}, status=status.HTTP_400_BAD_REQUEST)

        try:
            created_sequences = add_sequences_to_database(db, desired_numbers)
        except DatabaseLocked as exc:
            return Response({'message': f"Cannot modify sequences in a locked database."}, status=status.HTTP_400_BAD_REQUEST)
        except AccessionsAlreadyExist as exc:
            return Response({
                    'message': "Every accession number specified to be added already exists in the database.", 
                    'accession_numbers': exc.accession_numbers}, 
                    status=status.HTTP_400_BAD_REQUEST)
        except AccessionLimitExceeded as exc:
            return Response({
                    'message': f"Bulk addition of sequences requires a list of accession numbers and the list length must have 1 to {exc.max_accessions} elements (inclusive).", 
                    "accession_numbers": desired_numbers 
                }, 
                status=status.HTTP_400_BAD_REQUEST)
        except GenBankConnectionError as exc:
            return Response({
                    'message': f"The GenBank database could not be connected to.", 
                    "accession_numbers": desired_numbers, 
                    "queried_numbers": exc.queried_accessions
                }, 
                status=status.HTTP_502_BAD_GATEWAY)
        except InsufficientAccessionData as exc:
            return Response({
                    'message': f"The GenBank database failed to return some accession numbers so no numbers were added at all.", 
                    "accession_numbers": desired_numbers,
                    "missing_accessions": exc.missing_accessions
                }, 
                status=status.HTTP_400_BAD_REQUEST)
        else:
            # serialize the list of created objects so they can be sent back
            serialized_data = BlastDbSequenceEntrySerializer(created_sequences, many = True).data
            return Response(serialized_data, status = status.HTTP_201_CREATED)

    @swagger_auto_schema(
        operation_summary='Update accessions from GenBank data.',
        operation_description=f'From a list of accession numbers, update their accessions by refetching data from GenBank. List must contain between 1-{MAX_ACCESSIONS_PER_REQUEST} accession numbers inclusive.',
        tags = [tag_blastdbs, tag_sequences],
        manual_parameters=[openapi.Parameter(
            name='id',
            description='Id of BLAST database to update accession numbers in',
            in_=openapi.IN_PATH, 
            type=openapi.TYPE_STRING,
            format=openapi.FORMAT_UUID
        )],
        request_body=NuccoreSequenceBulkAddSerializer,
        responses={
            '200': openapi.Response(
                description='All accession numbers updated.',
                schema=BlastDbSequenceEntrySerializer(many=True),
                examples={
                    'application/json': [ BlastDbSequenceEntrySerializer.Meta.example ]
                }
            ),
            '400': 'Bad parameters in request. Example reasons: list may be empty or too long (>100); accession numbers may not exist in database yet; database may be locked',
            '403': 'Insufficient permissions.',
            '404': 'BLAST database matching the specified ID was not found.',
            '500': 'Unexpected error.',
            '502': 'Encountered error connecting to GenBank.'
        }
    )
    def patch(self, request, *args, **kwargs):
        '''
        Update certain accession numbers.
        '''
        pk = kwargs['pk']
        db: BlastDb
        # Query for the Blast database
        try:
            db = BlastDb.objects.get(id=pk)
        except BlastDb.DoesNotExist:
            return Response({'message': 'Database does not exist', 'requested': pk}, status.HTTP_404_NOT_FOUND)
        
        # Check if accession numbers provided
        if not 'accession_numbers' in request.data:
            return Response({'message': 'Missing required list of accession numbers to update.', 'accession_numbers': []}, status.HTTP_400_BAD_REQUEST)
        desired_numbers = request.data['accession_numbers']

        # Deny user if user has insufficient permissions
        if not NuccoreSharePermission.has_change_permission(user=request.user, obj=None):
            return Response(status=status.HTTP_403_FORBIDDEN)
        
        # Respond with an error if database is locked
        if db.locked:
            return Response({'message': 'The database is locked and its accession numbers cannot be added, edited or removed.'}, status=status.HTTP_400_BAD_REQUEST)

        try:
            updated_sequences = update_sequences_in_database(db, desired_numbers)
        except DatabaseLocked as exc:
            return Response({'message': f"Cannot modify sequences in a locked database."}, status=status.HTTP_400_BAD_REQUEST)
        except AccessionsNotFound as exc:
            return Response({'message': "Every accession number specified to be added already exists in the database.", 'accession_numbers': desired_numbers}, status=status.HTTP_200_OK)
        except AccessionLimitExceeded as exc:
            return Response({
                    'message': f"Bulk addition of sequences requires a list of accession numbers and the list length must have 1 to {exc.max_accessions} elements (inclusive).", 
                    "accession_numbers": desired_numbers 
                }, 
                status=status.HTTP_400_BAD_REQUEST)
        except GenBankConnectionError as exc:
            return Response({
                    'message': f"The GenBank database could not be connected to.", 
                    "accession_numbers": desired_numbers, 
                    "queried_numbers": exc.queried_accessions
                }, 
                status=status.HTTP_502_BAD_GATEWAY)
        except InsufficientAccessionData as exc:
            return Response({
                    'message': f"The GenBank database failed to return some accession numbers so no numbers were added at all.", 
                    "accession_numbers": desired_numbers,
                    "missing_accessions": exc.missing_accessions
                }, 
                status=status.HTTP_400_BAD_REQUEST)
        else:
            # serialize the list of updated objects so they can be sent back
            serialized_data = BlastDbSequenceEntrySerializer(updated_sequences, many = True).data
            return Response(serialized_data, status = status.HTTP_202_ACCEPTED)

    @swagger_auto_schema(
        operation_summary='Delete accessions from database.',
        operation_description=f'Delete all accession numbers from a given BLAST database.',
        tags = [tag_blastdbs, tag_sequences],
        manual_parameters=[openapi.Parameter(
            name='id',
            description='Id of BLAST database to update accession numbers in',
            in_=openapi.IN_PATH, 
            type=openapi.TYPE_STRING,
            format=openapi.FORMAT_UUID
        )],
        responses={
            '202': openapi.Response(
                description='Sequences deleted.',
                schema=openapi.Schema(
                        type=openapi.TYPE_OBJECT,
                        properties={
                            'deleted': openapi.Schema(
                                type=openapi.TYPE_STRING,
                                description='The number of objects deleted',
                                example=1
                            ),
                        }
                    ),
                examples={
                    'application/json': {
                        'deleted': 2
                    }
                }
            ),
            '400': 'Bad parameters in request. Database may be locked',
            '403': 'Insufficient permissions.',
            '404': 'BLAST database matching the specified ID was not found.',
            '500': 'Unexpected error.',
            '502': 'Encountered error connecting to GenBank.'
        }
    )
    def delete(self, request, *args, **kwargs):
        '''
        Delete all existing accessions
        '''
        pk = kwargs['pk']
        db: BlastDb
        # Query for the Blast database
        try:
            db = BlastDb.objects.get(id=pk)
        except BlastDb.DoesNotExist:
            return Response({'message': 'Database does not exist', 'requested': pk}, status.HTTP_404_NOT_FOUND)
        
        # Check if accession numbers provided
        if not 'accession_numbers' in request.data:
            return Response({'message': 'Missing required list of accession numbers to add.', 'accession_numbers': []}, status.HTTP_400_BAD_REQUEST)
        desired_numbers = request.data['accession_numbers']

        # Deny user if user has insufficient permissions
        if not NuccoreSharePermission.has_delete_permission(user=request.user, obj=None):
            return Response(status=status.HTTP_403_FORBIDDEN)
        
        # Respond with an error if database is locked
        if db.locked:
            return Response({'message': 'The database is locked and its accession numbers cannot be added, edited or removed.'}, status=status.HTTP_400_BAD_REQUEST)

        try:
            deleted = delete_sequences_in_database(db, desired_nums=desired_numbers)
        except DatabaseLocked as exc:
            return Response({'message': f"Cannot modify sequences in a locked database."}, status=status.HTTP_400_BAD_REQUEST)
        else:
            return Response({'deleted': deleted}, status = status.HTTP_202_ACCEPTED)

class NuccoreSequenceDetail(mixins.DestroyModelMixin, generics.RetrieveAPIView):
    '''
    Get or delete a sequence specified by ID
    '''
    queryset = NuccoreSequence.objects.all()
    serializer_class = NuccoreSequenceSerializer
    permission_classes = (NuccoreSequenceEndpointPermission,)

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
            '403': 'Insufficient permissions.',            
            '404': 'Sequence does not exist.',
            '500': 'Unexpected error.',
        }
    )
    def get(self, request, *args, **kwargs):
        db_primary_key = kwargs['pk']

        try:
            seq = NuccoreSequence.objects.get(id = db_primary_key)
        except NuccoreSequence.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        d: BlastDb = BlastDb.objects.get()

        try:
            self.check_object_permissions(request, obj=seq)
        except PermissionDenied:
            # Check if user has permissions to access object
            return Response(status=status.HTTP_403_FORBIDDEN)
        except:
            serializer = NuccoreSequenceSerializer(seq)
            return Response(serializer.data, status = status.HTTP_200_OK)

    @swagger_auto_schema(
        operation_summary='Delete a sequence entry.',
        operation_description='Delete a sequence entry from a BLAST database.',
        tags = [tag_blastdbs, tag_sequences],
        responses={
            '204': 'Deletion successful.',
            '400': 'Deletion failed because entry belongs to a database that is locked for editing.',
            '403': 'Insufficient permissions.',
            '404': 'Sequence does not exist',
            '500': 'Unexpected error.',
        }
    )
    def delete(self, request, *args, **kwargs):
        pk = kwargs['pk']
        try:
            seq : NuccoreSequence = NuccoreSequence.objects.get(id=pk)
            self.check_object_permissions(request, seq)
        except NuccoreSequence.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        except PermissionDenied:
            return Response(status=status.HTTP_403_FORBIDDEN)
        else:
            return self.destroy(request, *args, **kwargs)

class LibraryListView(mixins.ListModelMixin, generics.CreateAPIView):
    '''
    Return a list of all reference libraries.
    '''
    queryset = Library.objects.all()

    def get_serializer_class(self):
        if self.request is None:
            return LibrarySerializer
        elif self.request.method == 'POST':
            return LibraryCreateSerializer
        else:
            return LibrarySerializer

    @swagger_auto_schema(
        operation_summary='Get all reference libraries.',
        operation_description='Return a list of all reference libraries available to the user.',
        tags = [tag_libraries],
        responses={
            '200': openapi.Response(
                description='A list of all reference libraries viewable and runnable.',
                schema=LibraryListSerializer(many=True, read_only=True),
                examples={
                    'application/json': LibraryListSerializer.Meta.example,
                }
            )
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        List all reference libraries
        '''

        queryset = Library.objects.viewable(request.user)

        page = self.paginate_queryset(queryset)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        serializer = self.get_serializer(queryset, many=True)
        return Response(serializer.data)

    @swagger_auto_schema(
        operation_summary='Create a reference library.',
        operation_description='Create a reference library.',
        tags = [tag_libraries],
        request_body=LibraryCreateSerializer,
        responses={
            '200': openapi.Response(
                description='Creation was successful.',
                schema=LibraryCreateSerializer(),
                examples={
                    'application/json': LibraryCreateSerializer.Meta.example
                }
            ),
            '403': 'Insufficient permissions.',
        }
    )
    def post(self, request, *args, **kwargs):
        '''
        Create a reference library
        '''
        # Check that user can create a reference library
        if not LibrarySharePermissions.has_add_permission(request.user, None):
            return Response(status=status.HTTP_403_FORBIDDEN)

        serializer = self.get_serializer(data=request.data)
        if serializer.is_valid():
            serializer.validated_data['owner'] = request.user
            library_instance: Library = serializer.save()
            blast_db: BlastDb = BlastDb(library=library_instance)
            blast_db.save()
            library_data = LibrarySerializer(library_instance)
            return Response(library_data.data, status=status.HTTP_201_CREATED) 
            
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class LibraryDetailView(mixins.RetrieveModelMixin, mixins.UpdateModelMixin, mixins.DestroyModelMixin, generics.GenericAPIView):
    '''
    Retrieve the details of a given blast database
    '''
    queryset = Library.objects.all()
    permission_classes = [LibraryEndpointPermission]
    lookup_url_kwarg = 'library' # name of path parameter for primary key lookup

    def get_serializer_class(self):
        if not self.request or self.request.method == 'GET' or self.request.method == 'DELETE':
            return LibrarySerializer
        elif self.request.method == 'PATCH':
            return LibraryEditSerializer
        else:
            return LibrarySerializer

    @swagger_auto_schema(
        operation_summary='Get information about a reference library.',
        operation_description='Return data about a reference library.',
        tags = [tag_blastdbs],
        responses={
            '200': openapi.Response(
                description='Reference library updated successfully.',
                schema=LibrarySerializer(),
                examples={
                    'application/json': LibrarySerializer.Meta.example
                }
            ),
            '400': 'Bad request parameters.',
            '403': 'Insufficient permissions.',
            '404': 'Reference library with the given ID does not exist',
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        Retrieve data in the reference library.
        '''

        db_primary_key = kwargs['library']

        try:
            db = Library.objects.get(id = db_primary_key)
        except Library.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)

        # Check if user has access to this database
        try:
            self.check_object_permissions(request, db)
        except PermissionDenied:
            return Response(status=status.HTTP_403_FORBIDDEN)

        return Response(self.get_serializer_class()(db).data, status=status.HTTP_200_OK)

    @swagger_auto_schema(
        operation_summary='Update information of an existing reference library.',
        operation_description='Update an existing reference library with the content in the request.',
        tags = [tag_libraries],
        request_body=LibraryEditSerializer,
        responses={
            '200': openapi.Response(
                description='Reference library updated successfully.',
                schema=LibraryEditSerializer(),
                examples={
                    'application/json': LibraryEditSerializer.Meta.example
                }
            ),
            '400': 'Bad request parameters.',
            '403': 'Insufficient permissions.',
            '404': 'Reference library with the given ID does not exist.',
        }
    )
    def patch(self, request, *args, **kwargs):
        '''
        Partially update this reference library.
        '''
        try:
            db = Library.objects.get(id = kwargs['library'])
        except Library.DoesNotExist:
            return Response("Resource does not exist", status = status.HTTP_404_NOT_FOUND)
        
        # Check if user has access to this library
        try:
            self.check_object_permissions(request, db)
        except PermissionDenied:
            return Response(status=status.HTTP_403_FORBIDDEN)

        return self.partial_update(request, *args, **kwargs)

    @swagger_auto_schema(
        operation_summary='Delete a reference library.',
        operation_description='Delete the reference library matching the id specified by "library".',
        tags = [tag_libraries],
        responses={
            '204': 'Deletion successful.',
            '403': 'Insufficient permissions.',
            '404': 'Reference library with the given ID does not exist.',
            '500': 'Unexpected error.',
        }
    )
    def delete(self, request, *args, **kwargs):
        '''
        Delete the reference library
        '''
        # check if database exists
        try:
            libraryId = str(kwargs['library'])
            library: Library = Library.objects.get(id=libraryId)
        except Library.DoesNotExist:
            return Response(status = status.HTTP_404_NOT_FOUND)
        
        # check if user has delete permissions
        try:
            self.check_object_permissions(request, library)
        except PermissionDenied:
            return Response(status=status.HTTP_403_FORBIDDEN)

        try: 
            delete_library(library)
        except BaseException:
            return Response(status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)

class LibraryBlastDbList(mixins.ListModelMixin, generics.CreateAPIView):
    '''
    Return a list of all blast databases under a library
    '''
    queryset = BlastDb.objects.all()
    serializer_class = BlastDbCreateSerializer

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

        library_pk: str = kwargs['library']     
        try:
            library: Library = Library.objects.get(id=library_pk)
        except Library.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)

        if not LibrarySharePermissions.has_view_permission(request.user, library):
            return Response(status=status.HTTP_403_FORBIDDEN) 

        queryset = BlastDb.objects.viewable(request.user).filter(library=library_pk).reverse()

        page = self.paginate_queryset(queryset)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        serializer = self.get_serializer(queryset, many=True)
        return Response(serializer.data)

    @swagger_auto_schema(
        operation_summary='Create a database for this reference library.',
        operation_description="""
        Create a new publicly accessible BLAST database available for queries.
        
        Accession numbers can be immediately added to the database by specifying the numbers with `accession_numbers` and/or using the `base` parameter to specify the ID of a BLAST database from which to copy accession numbers from. Thus, resulting database with contain accession numbers found in the BLAST database specified by `base` and those found in `accession_numbers`.
        """,
        tags = [tag_admin, tag_blastdbs],
        # Request body must be manually specified because drf_yasg issues
        request_body=openapi.Schema(
            type=openapi.TYPE_OBJECT,
            properties={
                'description': openapi.Schema(
                    type=openapi.TYPE_STRING,
                    description='Description of the reference library version',
                    example='Sequences gathered as of January 2022'
                ),
                'locked': openapi.Schema(
                    type=openapi.TYPE_BOOLEAN,
                    description='Whether the version will be published and thus locked for future edits',
                    example=False
                ),
                'base': openapi.Schema(
                    type=openapi.TYPE_STRING,
                    format=openapi.FORMAT_UUID,
                    description='The ID of a BLAST database from which to copy accession numbers from into the new database.',
                    example='123e4567-e89b-12d3-a456-426614174000'
                ),
                'accession_numbers': openapi.Schema(
                    type=openapi.TYPE_ARRAY,
                    items=openapi.Schema(type=openapi.TYPE_STRING),
                    description='List of accession numbers to add to the database',
                    example=['ON303390', 'ON303391']
                )
            }
        ),
        responses={
            '201': openapi.Response(
                description='Creation was successful.',
                schema=BlastDbCreateSerializer(),
                examples={
                    'application/json': BlastDbCreateSerializer.Meta.example
                }
            ),
            '400': 'Bad request.',
            '403': 'Insufficient permissions.',
            '500': 'Unspecified error.',
            '502': 'Error connecting to GenBank.',
        }
    )
    def post(self, request, *args, **kwargs):
        '''
        Create a new blast database.
        '''
        library_pk: str = kwargs['library']        
        try:
            library: Library = Library.objects.get(id=library_pk)
        except Library.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)

        # Check that user can create a database in the reference library
        if not LibrarySharePermissions.has_change_permission(request.user, library):
            return Response(status=status.HTTP_403_FORBIDDEN)

        # TODO: Move logic to blastdb_controller.py
        updated_data = request.data.copy()
        updated_data['library'] = library.id
        serializer = self.get_serializer_class()(data=request.data)
        if serializer.is_valid():
            try:
                assert isinstance(serializer.validated_data, dict)
            except AssertionError:
                return Response(status=status.HTTP_500_INTERNAL_SERVER_ERROR)
            # Ensure that accession_numbers and base are parsed in order to be
            # passed to create_blastdb
            serializer_data = serializer.validated_data
            additional_accessions = serializer_data.pop('accession_numbers', [])
            base = serializer_data.pop('base', None)
            if not base is None:
                try:
                    base = BlastDb.objects.get(id=base)
                except BlastDb.DoesNotExist:
                    return Response({'message': f'Database with id {base} does not exist'}, status=status.HTTP_400_BAD_REQUEST)

            try:
                new_database: BlastDb = create_blastdb(additional_accessions=additional_accessions, base=base, **serializer_data, library=library)
            except AccessionLimitExceeded as exc:
                return Response({
                        'message': f"Bulk addition of sequences requires a list of accession numbers and the list length must have 1 to {exc.max_accessions} elements (inclusive)." 
                    }, 
                    status=status.HTTP_400_BAD_REQUEST)
            except GenBankConnectionError as exc:
                return Response({
                        'message': f"The GenBank database could not be connected to.", 
                        "queried_numbers": exc.queried_accessions
                    }, 
                    status=status.HTTP_502_BAD_GATEWAY)
            except InsufficientAccessionData as exc:
                return Response({
                        'message': f"The GenBank database failed to return some accession numbers so no numbers were added at all.", 
                        "missing_accessions": exc.missing_accessions
                    }, 
                    status=status.HTTP_400_BAD_REQUEST)
            except BaseException as exc:
                raise exc
            else:
                # serialize the list of created objects so they can be sent back
                return Response(BlastDbSerializer(BlastDb.objects.get(id=new_database.id)).data, status = status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class BlastDbDetail(mixins.RetrieveModelMixin, mixins.UpdateModelMixin, mixins.DestroyModelMixin, generics.GenericAPIView):
    '''
    Retrieve the details of a given blast database
    '''
    queryset = BlastDb.objects.all()
    renderer_classes = [JSONRenderer, BrowsableAPIRenderer, TemplateHTMLRenderer, BlastDbCSVRenderer, BlastDbFastaRenderer]
    permission_classes = [BlastDbEndpointPermission]

    def get_serializer_class(self):
        if not self.request or self.request.method == 'GET':
            return BlastDbSerializer
        elif self.request.method == 'PATCH':
            return BlastDbEditSerializer
        elif self.request.method == 'DELETE':
            return BlastDbEditSerializer
        else:
            return BlastDbSerializer

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
            '403': 'Insufficient permissions.',
            '404': 'Database with the ID does not exist',
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        Retrieve data in the blastdb.
        '''

        db_primary_key = kwargs['pk']

        try:
            db: BlastDb = BlastDb.objects.get(id = db_primary_key)
            # Check if user has access permissions to this database
            self.check_object_permissions(request, db)
        except BlastDb.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        except PermissionDenied:
            return Response(status=status.HTTP_403_FORBIDDEN)

        response = Response(self.get_serializer_class()(db).data, status=status.HTTP_200_OK, template_name='blastdb.html')

        # based on the media type file returned, specify attachment and file name
        file_name_root = get_data_fishdb_path(db)
        if request.accepted_media_type.startswith('text/csv'):
            response['Content-Disposition'] = f'attachment; filename={file_name_root}.csv";'
        elif request.accepted_media_type.startswith('text/x-fasta'):
            response['Content-Disposition'] = f'attachment; filename="{file_name_root}.fasta";'

        return response

    @swagger_auto_schema(
        operation_summary='Update information of an existing BLAST database.',
        operation_description='Edit information of a BLAST database. This method does not allow adding/removing/editing of the sequence entries within.',
        tags = [tag_blastdbs],
        request_body=BlastDbEditSerializer(),
        responses={
            '200': openapi.Response(
                description='Database updated successfully.',
                schema=BlastDbSerializer(),
                examples={
                    'application/json': BlastDbSerializer.Meta.example
                }
            ),
            '400': 'Bad request parameters.',
            '403': 'Insufficient permissions.',
            '404': 'Database with the ID does not exist',
        }
    )
    def patch(self, request, *args, **kwargs):
        '''
        Partially update this blastdb.
        '''        
        try:
            db: BlastDb = BlastDb.objects.get(id = kwargs['pk'])
            # Check if user has access to this database
            self.check_object_permissions(request, db)
        except BlastDb.DoesNotExist:
            return Response("Resource does not exist", status = status.HTTP_404_NOT_FOUND)
        except PermissionDenied:
            return Response(status=status.HTTP_403_FORBIDDEN)

        updated_data: Dict[str, Any] = request.data.copy()
        new_lock: bool = False
        if not db.locked and 'locked' in updated_data:
            # We want to set 'locked' to false first
            # so we can lock it later on if needed
            locked = updated_data.pop('locked', 'False')[0]
            new_lock = (locked.lower() == 'true')
            updated_data['locked'] = False
        
        serializer = self.get_serializer_class()(db, data=updated_data, partial=True) 
        if serializer.is_valid():
            updated = serializer.save()
        else:
            return Response(status=status.HTTP_400_BAD_REQUEST)

        if new_lock:
            # If the request locks the database, we need to save the version
            updated = save_blastdb(updated, perform_lock=True)
        return Response(status=status.HTTP_204_NO_CONTENT)

    @swagger_auto_schema(
        operation_summary='Delete a BLAST database.',
        operation_description='Delete the BLAST database specified by the given ID.',
        tags = [tag_blastdbs],
        responses={
            '204': 'Deletion successful.',
            '403': 'Insufficient permissions.',
            '404': 'BLAST database with the given ID does not exist',
            '500': 'Unexpected error.',
        }
    )
    def delete(self, request, *args, **kwargs):
        '''
        Delete the BLAST database
        '''
        # check if database exists
        try:
            database_id = str(kwargs['pk'])
            db = BlastDb.objects.get(id = database_id)
        except BlastDb.DoesNotExist:
            return Response(status = status.HTTP_404_NOT_FOUND)
        
        # check if user has delete permissions
        try:
            self.check_object_permissions(request, db)
        except PermissionDenied as exc:
            return Response({'message': exc.detail}, status=status.HTTP_403_FORBIDDEN)
        
        try:
            delete_blastdb(db)
        except BaseException:
            return Response(status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        else:
            return Response(status=status.HTTP_204_NO_CONTENT)

class BlastDbExport(generics.GenericAPIView):
    serializer_class = BlastDbSerializer
    renderer_classes = [JSONRenderer, BrowsableAPIRenderer, BlastDbFastaRenderer, BlastDbCompatibleRenderer]

    def get_queryset(self):
        return BlastDb.objects.viewable(self.request.user)

    def get_renderer_context(self):
        context = super().get_renderer_context()
        context['fasta_format'] = self.request.GET.get('fasta_format', '')
        return context
    
    def get(self, request, *args, **kwargs):
        '''
        Export blastdb to compatible formats for taxonomic assignment.
        '''

        db_primary_key = kwargs['pk']

        try:
            db: BlastDb = BlastDb.objects.get(id = db_primary_key)
            # Check if user has access permissions to this database
            self.check_object_permissions(request, db)
        except BlastDb.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)
        except PermissionDenied:
            return Response(status=status.HTTP_403_FORBIDDEN)

        response = Response(self.get_serializer_class()(db).data, status=status.HTTP_200_OK, template_name='blastdb.html')

        # based on the media type file returned, specify attachment and file name
        if request.accepted_media_type.startswith('text/x-fasta'):
            ff = self.request.GET.get('fasta_format', 'fasta')
            response['Content-Disposition'] = f'attachment; filename="{db.custom_name}.fasta";'
        elif request.accepted_media_type.startswith('application/zip'):
            ff = self.request.GET.get('fasta_format', 'fasta')
            response['Content-Disposition'] = f'attachment; filename="{db.custom_name}.zip";'
        else:
            # If the accepted media type cannot be chosen (ie content negotiation)
            # fails, then return a 406 response
            return Response(status=status.HTTP_406_NOT_ACCEPTABLE)
        return response

class BlastRunList(mixins.ListModelMixin, mixins.CreateModelMixin, generics.GenericAPIView):
    '''
    List all BLAST runs. Used only for database administration purposes.
    '''

    serializer_class = BlastRunSerializer

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
        List all blast runs viewable
        '''
        return self.list(request, *args, **kwargs)
    
    def get_queryset(self):
        if isinstance(self.request.user, User):
            return BlastRun.objects.listable(self.request.user)
        else:
            return BlastRun.objects.none()

class BlastRunRun(generics.CreateAPIView):
    '''
    Run a blast query
    '''
    serializer_class = BlastRunRunSerializer
    permission_classes = [BlastRunEndpointPermission]
    queryset = BlastRun.objects.all()
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    
    @swagger_auto_schema(
        operation_summary='Submit a run.',
        operation_description='Submit query sequence(s) and associated data to begin a run. Sequences can be supplied as either raw text in "query_sequence" or in an uploaded file named "query_file". Unless the submission is erroneous, the response will contain a unique run ID used to keep track of the status (the run is run asynchronously) and look up results when complete.\n\nA run will perform a BLASTN query of each query sequence against the sequences found within the BLAST database. Based on the values of the "create_hit_tree" and "create_db_tree" parameters, then up to 2 multiple alignment jobs will be performed using the query sequences and sequences from the hits or entire database, respectively.',
        # TODO: Add request / parameter schema

        # request_body=openapi.Schema(
        #     type=openapi.TYPE_OBJECT,
        #     title='Parameters',
        #     properties={
        #         'job_name': openapi.Schema(
        #             type=openapi.TYPE_STRING,
        #             example='two sequences'
        #         ),
        #         'query_sequence': openapi.Schema(
        #             type=openapi.TYPE_STRING,
        #             example='>MG653404.1 Compsaraia iara isolate SA217 cytochrome c oxidase subunit 1 (COI) gene, partial cds; mitochondrial\nCCAACCAGGCGCCCTCCTGGGAGACGACCAAATTTACAATGTGGTCGTTACCGCCCATGCCTTCGTAATAATTTTCTTTATAGTAATGCCAATTATAATCGGAGGCTTTGGCAATTGACTTATCCCCTTAATAATTGCCGCGCCCGATATGGCATTCCCACGAATAAATAATATAAGCTTCTGACTGCTTCCCCCATCATTCTTCCTCCTACTTGCCTCTGCCGGGTTAGAGGCCGGAGTCGGGACAGGCTGAACGCTTTACCCCCCTCTTGCCGGTAATGCAGCACACGCTGGAGCCTCTGTAGACCTAACCATTTTCTCCCTTCACTTGGCCGGTGTCTCATCTATCCTCGGATCTATTAACTTTATCACTACAATTATTAATATGAAACCCCCAACAATATCCCAATACCAGCTTCCATTATTTATTTGATCCTTACTAGTAACCACAGTACTTCTACTACTCTCCCTTCCTGTTCTAGCTGCTGGA\n>MK401918.1 Porotergus gymnotus isolate ANSP189647 cytochrome oxidase subunit I (COI) gene, partial cds; mitochondrial\nGGCACACTATACATGGTGTTNGGCGCCTGGGCGGGTATAATTGGTACTGCTCTTANGCTTCTAATCCGGGCCGAGCTAAATCAACCGGGCACCCTCCTAGAAGACGACCAAATTTACAATGTGGCCGTCACTGCCCATGCCTTTGTAATAATTTTCTTTATAGTTATGCCAATCATAATTGGAGGCTTTGGCAATTGGCTTATCCCCCTAATAATTGCCGCGCCAGACATAGCATTCCCACGAATAAATAACATAAGCTTCTGGCTACTCCCCCCATCATTCTTCCTGCTCCTCGCCTCTGCTGGCTTAGAGGCCGGAGTTGGAACAGGTTGGACCCTATACCCCCCTCTTGCCGCTAATGCAGCACACGCCGGAGCTTCCGTAGACCTAACTATCTTCTCCCTTCACTTGGCGGGTGTTTCATCCATCCTCGGCTCCATTAACTTTATTACTACAATTATTAATATAAAACCTCCAACAATATCCCAATACCAACTCCCACTCTTTATCTGGTCCCTGCTGGTTACTACCGTGCTTCTACTACTCTCCCTTCCTGTCCTAGCTGCTGGTATTACCATGCTACTCACAGACCGAAATCTAAACACAGCATTCTTCGACCCAACGGGCGGCGGTGACCCAATTCTGTACCAACACTTGTTCTGGTT'
        #         ),
        #         'query_file': openapi.Schema(
        #             type=openapi.TYPE_FILE
        #         ),
        #         'create_hit_tree': openapi.Schema(
        #             type=openapi.TYPE_BOOLEAN,
        #             example=True
        #         ),
        #         'create_db_tree': openapi.Schema(
        #             type=openapi.TYPE_BOOLEAN,
        #             example=False
        #         )
        #     }),
        tags = [tag_runs]
        # responses={
        #     '400': openapi.Response(
        #         description='Bad request parameters. An accompanying message may specify the error with the request.',
        #         schema=openapi.Schema(
        #             type=openapi.TYPE_OBJECT,
        #             properties={
        #                 'message': openapi.Schema(
        #                     type=openapi.TYPE_STRING,
        #                     description='Helper message clarifying the cause of the error.'
        #                 )
        #             }
        #         ),
        #         examples={
        #             'application/json': {
        #                 'message': 'Request did not have raw sequence text or file upload specified to run with.'
        #             }
        #         }
        #         ),
        #     '404': 'The BLAST database specified by the ID does not exist.',
        #     '201': openapi.Response(
        #         description='Run was successfully added to queue and a new unique run identifier is returned. Users should now use the given ID to continually check the status of the run to check whether it has completed.',
        #         schema=BlastRunSerializer(),
        #         examples={'application/json': BlastRunSerializer.Meta.example}
        #     ),
        #     '500': 'Unexpected error. No new run was created.'
        # }
    )
    def post(self, request, *args, **kwargs):      
        '''
        Run blast query
        '''

        db_used = kwargs['pk']

        print("Beginning run on database", str(db_used))

        # Bad request if database could not be found
        try:
            odb = BlastDb.objects.get(id = db_used)
        except BlastDb.DoesNotExist:
            return Response({
                'message': f"The database specified (id = {db_used}) does not exist."
            }, status = status.HTTP_400_BAD_REQUEST)
        
        # Check whether user has permissions to run on this database
        if not DatabaseSharePermissions.has_run_permission(request.user, obj=odb):
            return Response({
                'message': 'Insufficient permissions to run BLAST on this database.'
            }, status=status.HTTP_403_FORBIDDEN)

        # Bad request if no query_sequence is given
        if not 'query_sequence' in request.data and not 'query_file' in request.FILES:
            return Response({
                'message': 'Request did not have raw sequence text or file upload specified to run with.'
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

        print("Beginning to parse sequences")

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

        job_name = request.data['job_name'] if 'job_name' in request.data else ''
        
        # path to output the results
        results_uuid = uuid.uuid4()
        results_path = get_data_run_path(str(results_uuid))
        if not os.path.exists(results_path):
            os.mkdir(results_path)
            # ensure that celery worker has access to file to write results.txt later
            # shutil.chown(results_path, group='appgroup', user='appuser') 
            os.chmod(results_path, 0o774)

        # static path to download results
        static_path = get_static_run_path(run_id=str(results_uuid))
        if not os.path.exists(static_path):
            os.mkdir(static_path)
            # shutil.chown(static_path, group='appgroup', user='appuser')
            os.chmod(static_path, 0o775)

        print("Run id: " + str(results_uuid))

        if not odb.locked:
            return Response({
                'message': 'The specified database is not published yet for queries',
            }, status=status.HTTP_400_BAD_REQUEST)
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
        ncbi_blast_version = 'ncbi-blast-2.12.0+'
        run_details = BlastRun(id = results_uuid, db_used = odb, job_name = job_name, blast_version = ncbi_blast_version, errors = '', status=BlastRun.JobStatus.QUEUED, start_time = None, end_time = None, error_time = None, create_hit_tree = create_hit_tree, create_db_tree = create_db_tree)
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

        # How to chain for Docker-hosted RabbitMQ 
        app.send_task('barcode_blastn.tasks.run_blast_command', 
            queue='BarcodeQueue.fifo',
            **message_properties,
            kwargs={
                    'ncbi_blast_version': ncbi_blast_version,
                    'fishdb_id': str(odb.id),
                    'run_id': run_id_str
        },
        chain=[signature('barcode_blastn.tasks.performAlignment', 
            queue='BarcodeQueue.fifo',
                kwargs={
                    'run_id': run_id_str
                }
            )])

        # How to chain for Amazon SQS
        # run_blast_command.apply_async(
        #     queue='BarcodeQueue.fifo',
        #     **message_properties, 
        #     kwargs={
        #         'ncbi_blast_version': ncbi_blast_version,
        #         'fishdb_id': str(odb.id),
        #         'run_id': run_id_str
        #     }, 
        #     link=signature('barcode_blastn.tasks.performAlignment', 
        #         kwargs={
        #             'run_id': run_id_str
        #         }
        #     )  # type: ignore
        # )      # type: ignore

        # create response 
        serializer = BlastRunSerializer(run_details)
        return Response(serializer.data, status=status.HTTP_201_CREATED)

class BlastRunDetail(mixins.DestroyModelMixin, generics.GenericAPIView):
    '''
    View the results of any given run
    '''

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    permission_classes = [AllowAny,]
    
    @swagger_auto_schema(
        operation_summary='Get run results.',
        operation_description='Get run information and results of a BLAST run. The status indicates the status of the run. The results returned are complete once the status is "FIN" (i.e. run is finished). The run information includes the BLASTN hits and multiple alignment trees.',
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
    permission_classes = (AllowAny,)
    renderer_classes = (BlastRunFastaRenderer,)
    
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
        query_path = get_data_run_path(str(run.id)) + '/' + file_name
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
    permission_classes = (AllowAny,)
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
    permission_classes = (AllowAny,)

    @swagger_auto_schema(
        operation_summary='Get status of run',
        operation_description='Returns a minimal set of information useful for polling/checking the status of run.\n\nHow to interpret status:\n"QUE" (Queued): The run is currently waiting in the queue for its turn to be processed.\n"STA" (Started): The run is currently being processed.\n"FIN" (Finished): The run successfully completed and complete results are now visible.\n"ERR" (Errored): The run encountered an unexpected error and terminated.\n"UNK" (Unknown): The status is unknown, likely because there was an unexpected database or server error. Please submit another run.\n"DEN" (Denied): The run submission was received by the server, but it was denied from being processed.',
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
