# TODO: Install Python 3.10.6 to match EC2

import io
import os
import re
import uuid
from typing import Any, Dict, List
from barcode_blastn.ordering import CustomHitOrderingFilter, CustomSequenceOrderingFilter
from barcode_blastn.pagination import BlastDbSequencePagination, BlastRunHitPagination, BlastRunQueryPagination
from barcode_blastn.tests import LibraryListTest, SequenceTester, LibraryCreateTest

from barcode_identifier_api.celery import app
from Bio import SeqIO
from celery import signature, chain
from barcode_blastn.tasks import run_blast_command
from django.contrib.auth import login
from django.contrib.auth.models import User
from django.db.models import QuerySet, Count, When, Case, F, Value, PositiveSmallIntegerField, CharField
from django.db.models.functions import Substr, StrIndex, Length, Trim
from django.http import Http404
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
from rest_framework.renderers import JSONRenderer
from rest_framework.response import Response
from rest_framework.filters import OrderingFilter

from barcode_blastn.controllers.blastdb_controller import (
    AccessionsAlreadyExist, AccessionsNotFound, DatabaseLocked, InsufficientAccessionData,
    add_sequences_to_database, create_blastdb, delete_blastdb, delete_library, delete_sequences_in_database, filter_sequences_in_database, log_deleted_sequences,
    save_blastdb, update_sequences_in_database)
from barcode_blastn.file_paths import (get_data_fishdb_path, get_data_run_path,
                                       get_ncbi_folder, get_static_run_path)
from barcode_blastn.helper.parse_gb import (ACCESSIONS_PER_REQUEST,
                                            AccessionLimitExceeded,
                                            GenBankConnectionError, retrieve_gb)
from barcode_blastn.helper.verify_query import verify_dna, verify_header
from barcode_blastn.models import (BlastDb, BlastQuerySequence, BlastRun, Hit,
                                   Library, NuccoreSequence)
from barcode_blastn.permissions import (BlastDbEndpointPermission, BlastRunEndpointPermission,
                                        DatabaseSharePermissions,
                                        LibraryEndpointPermission,
                                        LibrarySharePermissions,
                                        NuccoreSequenceEndpointPermission,
                                        NuccoreSharePermission)
from barcode_blastn.renderers import (BlastDbCSVRenderer, BlastDbCompatibleRenderer, BlastDbFastaRenderer, BlastDbTSVRenderer, BlastDbXMLRenderer,
                                      BlastRunHitsCSVRenderer,
                                      BlastRunHitsTSVRenderer,
                                      BlastRunFastaRenderer,
                                      BlastRunHitsTxtRenderer, BlastRunTaxonomyCSVRenderer, BlastRunTaxonomyTSVRenderer)
from barcode_blastn.serializers import (BlastDbCreateSerializer,
                                        BlastDbEditSerializer, BlastDbExportSerializer, BlastDbHistoricalSerializer,
                                        BlastDbListSerializer,
                                        BlastDbSequenceEntrySerializer,
                                        BlastDbSerializer, BlastQuerySequenceSerializer, BlastQuerySequenceShortSerializer, BlastRunResultsShortSerializer,
                                        BlastRunRunSerializer,
                                        BlastRunSerializer,
                                        BlastRunStatusSerializer, BlastRunTaxonomySerializer, HitSerializer, HitTaxonomySerializer,
                                        LibraryCreateSerializer,
                                        LibraryEditSerializer,
                                        LibraryListSerializer,
                                        LibrarySerializer,
                                        NuccoreSequenceBulkAddSerializer,
                                        NuccoreSequenceSerializer,
                                        UserSerializer)

tag_libraries = 'Libraries'
tag_runs = 'Runs'
tag_blastdbs = 'BLAST Databases'
tag_sequences = 'GenBank Accessions'
tag_admin = 'Admin Tools'
tag_users = 'User Authentication'

class UserDetailView(generics.RetrieveAPIView):
    '''
    Return the user details associated with the authenticated
    user
    '''
    authentication_classes = (TokenAuthentication,)
    permission_classes = (IsAuthenticated,)
    serializer_class = (UserSerializer,)

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
    def get(self, request, *args, **kwargs):
        return super(UserDetailView, self).get(request, *args, **kwargs)

    def get_object(self):
        return self.request.user

class LogoutAllView(KnoxLogoutAllView):
    @swagger_auto_schema(
        security=[{'Basic': []}, {'Bearer': []}],
        operation_summary='User logoff all',
        operation_description='Sign out of the user account by signing out of all tokens corresponding to the user specified by the provided authentication token. The token is provided in the request header with the format `Authorization: Bearer <token>`, where <token> is the full token.',
        tags = [tag_users],
    )
    def post(self, request, format=None):
        return super(LogoutAllView, self).post(request, format)

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
        return super(LogoutView, self).post(request, format)

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

class BlastDbSequenceList(mixins.UpdateModelMixin, mixins.DestroyModelMixin, mixins.ListModelMixin, generics.CreateAPIView):   
    '''
    Retrieve, filter and bulk add sequences under a specified database.
    '''
    serializer_class = BlastDbSequenceEntrySerializer
    permission_classes = [NuccoreSequenceEndpointPermission]
    queryset = NuccoreSequence.objects.all()
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    pagination_class = BlastDbSequencePagination
    ordering_fields = '__all__'
    ordering = ['version']
    filter_backends = [CustomSequenceOrderingFilter]

    def get_queryset(self):
        pk = self.kwargs.get('pk')
        try:
            database = BlastDb.objects.get(id=pk)
        except BlastDb.DoesNotExist:
            raise Http404
        if DatabaseSharePermissions.has_view_permission(self.request.user, database):
            return database.sequences.all()
        else:
            raise PermissionDenied

    @swagger_auto_schema(
        operation_summary='Retrieve all sequences in the database.',
        operation_description=f'Get an object featuring a paginated list of all sequences in the reference database.',
        tags = [tag_blastdbs, tag_sequences],
        responses={
            '200': openapi.Response(
                description='Sequences successfully retrieved.',
                schema=NuccoreSequenceBulkAddSerializer(many=True),
            ),
            '404': 'BLAST database matching the specified ID was not found.',
        }
    )
    def get(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)
        
    @swagger_auto_schema(
        operation_summary='Add accession numbers to database.',
        operation_description=f'From a list of accession numbers, add them to an existing database. List must contain between 1-{ACCESSIONS_PER_REQUEST} accession numbers inclusive.',
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
                    'application/json': SequenceTester.post_blastdbs_id_sequences_201_response
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
        serializer = NuccoreSequenceBulkAddSerializer(data=request.data)
        if not serializer.is_valid() or not isinstance(serializer.validated_data, dict):
            return Response(status=status.HTTP_400_BAD_REQUEST)
        # Subset kwargs 
        min_length = serializer.validated_data.get('min_length', -1)
        max_length = serializer.validated_data.get('min_length', -1)
        max_ambiguous_bases = serializer.validated_data.get('max_ambiguous_bases', -1)
        blacklist = serializer.validated_data.get('blacklist', [])
        require_taxonomy = serializer.validated_data.get('require_taxonomy', False)
        desired_numbers.extend(serializer.validated_data.get('accession_numbers', []))
        search_term = serializer.validated_data.get('search_term', None)

        # Deny user if user has insufficient permissions
        if not NuccoreSharePermission.has_add_permission(user=request.user, obj=None):
            return Response(status=status.HTTP_403_FORBIDDEN)
        
        # Respond with an error if database is locked
        if db.locked:
            return Response({'message': 'The database is locked and its accession numbers cannot be added, edited or removed.'}, status=status.HTTP_400_BAD_REQUEST)

        try:
            created_sequences = add_sequences_to_database(db, request.user, desired_numbers=desired_numbers, search_term=search_term, \
                min_length=min_length, max_length=max_length, max_ambiguous_bases=max_ambiguous_bases, blacklist=blacklist, \
                require_taxonomy=require_taxonomy)
        except DatabaseLocked as exc:
            return Response({'message': f"Cannot modify sequences in a locked database."}, status=status.HTTP_400_BAD_REQUEST)
        except AccessionsAlreadyExist as exc:
            return Response({
                    'message': "Every accession number specified to be added already exists in the database.", 
                    'accession_numbers': exc.accession_numbers}, 
                    status=status.HTTP_400_BAD_REQUEST)
        except AccessionLimitExceeded as exc:
            return Response({
                    'message': f"Bulk addition of sequences limited to 1 to {exc.max_accessions} records per operation (inclusive). Current number: {exc.curr_accessions}",
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
                    'message': f"Could not find a matching record for at least one accession number, or no records could be returned from the query.", 
                    "accession_numbers": desired_numbers,
                    "missing_accessions": exc.missing_accessions,
                    "term": exc.term,
                }, 
                status=status.HTTP_400_BAD_REQUEST)
        else:
            # serialize the list of created objects so they can be sent back
            serializer = BlastDbSequenceEntrySerializer(created_sequences, many = True).data
            return Response(serializer, status = status.HTTP_201_CREATED)

    @swagger_auto_schema(
        operation_summary='Update accessions from GenBank data.',
        operation_description=f'From a list of accession numbers, update their accessions by refetching data from GenBank. List must contain between 1-{ACCESSIONS_PER_REQUEST} accession numbers inclusive.',
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
                    'application/json': SequenceTester.post_blastdbs_id_sequences_201_response
                }
            ),
            '400': 'Bad parameters in request. Example reasons: list is too long (>100); accession numbers may not exist in database yet; database may be locked',
            '403': 'Insufficient permissions.',
            '404': 'BLAST database matching the specified ID was not found.',
            '500': 'Unexpected error.',
            '502': 'Encountered error connecting to GenBank.'
        }
    )
    def patch(self, request, *args, **kwargs):
        '''
        Filter accession numbers
        '''
        pk = kwargs['pk']
        db: BlastDb
        # Query for the Blast database
        try:
            db = BlastDb.objects.get(id=pk)
        except BlastDb.DoesNotExist:
            return Response({'message': 'Database does not exist', 'requested': pk}, status.HTTP_404_NOT_FOUND)
        
        # Deny user if user has insufficient permissions
        if not NuccoreSharePermission.has_change_permission(user=request.user, obj=None):
            return Response(status=status.HTTP_403_FORBIDDEN)
        
        # Respond with an error if database is locked
        if db.locked:
            return Response({'message': 'The database is locked and its accession numbers cannot be added, edited or removed.'}, status=status.HTTP_400_BAD_REQUEST)

        serializer = NuccoreSequenceBulkAddSerializer(data=request.data)
        if serializer.is_valid() and isinstance(serializer.validated_data, dict):
            # Subset kwargs 
            min_length = serializer.validated_data.get('min_length', -1)
            max_length = serializer.validated_data.get('min_length', -1)
            max_ambiguous_bases = serializer.validated_data.get('max_ambiguous_bases', -1)
            blacklist = serializer.validated_data.get('blacklist', [])
            require_taxonomy = serializer.validated_data.get('require_taxonomy', False)
            updated = filter_sequences_in_database(db, user=request.user, min_length=min_length, max_length=max_length, \
                max_ambiguous_bases=max_ambiguous_bases, blacklist=blacklist, \
                require_taxonomy=require_taxonomy)
            return Response(updated, status=status.HTTP_200_OK)
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

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
                                type=openapi.TYPE_NUMBER,
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
            deleted = delete_sequences_in_database(db, user=request.user, desired_nums=desired_numbers)
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
    lookup_url_kwarg = 'pk'

    @swagger_auto_schema(
        operation_summary='Get information about a sequence entry.',
        operation_description='Get information about a sequence entry.',
        tags = [tag_sequences],
        responses={
            '200': openapi.Response(
                description='Successfully returned sequence information',
                schema=NuccoreSequenceSerializer(read_only=True),
                examples={
                    'application/json': SequenceTester.get_nuccores_id_200_response
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

        try:
            self.check_object_permissions(request, obj=seq)
        except PermissionDenied:
            # Check if user has permissions to access object
            return Response(status=status.HTTP_403_FORBIDDEN)
        else:
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
            database = seq.owner_database
            log_deleted_sequences([seq], database)
            save_blastdb(database, request.user, perform_lock=False)
            return self.destroy(request, *args, **kwargs)

class LibrariesList(mixins.ListModelMixin, generics.CreateAPIView):
    '''
    Return a list of all reference libraries.
    '''
    queryset = Library.objects.all()
    filter_backends = [OrderingFilter]
    ordering_fields = '__all__'
    ordering = ['custom_name']

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
                    'application/json': LibraryListTest.get_libraries_200_response
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
            '201': openapi.Response(
                description='Creation was successful.',
                schema=LibraryCreateSerializer,
                examples={
                    'application/json': LibraryCreateTest.post_libraries_201_request
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
                    'application/json': LibraryListTest.get_libraries_id_200_response
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
                    'application/json': LibraryListTest.patch_libraries_id_204_response
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
    Return a list of all runnable blast databases under a library
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

    def get_queryset(self):
        library_pk: str = self.kwargs.get('library')
        try:
            library: Library = Library.objects.get(id=library_pk)
        except Library.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)

        if not LibrarySharePermissions.has_view_permission(self.request.user, library):
            return Response(status=status.HTTP_403_FORBIDDEN) 
        
        return BlastDb.objects.runnable(self.request.user).filter(library=library_pk).order_by('genbank_version', 'major_version', 'minor_version') 

    @swagger_auto_schema(
        operation_summary='Get all databases.',
        operation_description='Return a list of all BLAST databases publicly available for queries.',
        tags = [tag_blastdbs],
        responses={
            '200': openapi.Response(
                description='List all database versions under a reference library.',
                schema=BlastDbListSerializer(many=True, read_only=True),
                examples={
                    'application/json': LibraryListTest.get_libraries_id_versions_200_response,
                }
            )
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        List all blast databases under a specific library
        '''
        return super().list(request, *args, **kwargs)

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
                    # example=LibraryListTest.post_libraries_id_versions_request['description']
                ),
                'locked': openapi.Schema(
                    type=openapi.TYPE_BOOLEAN,
                    description='Whether the version will be published and thus locked for future edits',
                    # example=LibraryListTest.post_libraries_id_versions_request['locked']
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
                    # example=LibraryListTest.post_libraries_id_versions_request['accession_numbers']
                )
            }
        ),
        responses={
            '201': openapi.Response(
                description='Creation was successful.',
                schema=BlastDbCreateSerializer(),
                examples={
                    'application/json': LibraryListTest.post_libraries_id_versions_201_response
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
                new_database: BlastDb = create_blastdb(additional_accessions=additional_accessions, user=request.user, base=base, **serializer_data, library=library)
            except AccessionLimitExceeded as exc:
                return Response({
                        'message': f"Bulk addition of sequences limited to 1 to {exc.max_accessions} accessions (inclusive). Current number: {exc.curr_accessions}", 
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
                        'message': f"Could not find a matching record for at least one accession number, or no records could be returned from the query.", 
                        "missing_accessions": exc.missing_accessions,
                        "term": exc.term,
                    }, 
                    status=status.HTTP_400_BAD_REQUEST)
            except BaseException as exc:
                raise exc
            else:
                # serialize the list of created objects so they can be sent back
                return Response(BlastDbSerializer(BlastDb.objects.get(id=new_database.id)).data, status = status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class BlastDbDetail(mixins.UpdateModelMixin, mixins.DestroyModelMixin, generics.RetrieveAPIView):
    '''
    Retrieve the details of a given blast database
    '''
    queryset = BlastDb.objects.all()
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
                    'application/json': LibraryListTest.get_blastdbs_id_200_response
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
        return super().get(request, *args, **kwargs)

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
                    'application/json': LibraryListTest.patch_blastdbs_id_204_response
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
        was_locked: bool = db.locked
        if not db.locked and 'locked' in updated_data:
            # We want to set 'locked' to false first
            # so we can lock it later on if needed
            locked = updated_data.pop('locked', False)
            updated_data['locked'] = locked
        
        serializer = self.get_serializer_class()(db, data=updated_data, partial=True) 
        if serializer.is_valid():
            updated = serializer.save()
        else:
            return Response(status=status.HTTP_400_BAD_REQUEST)

        if was_locked and serializer.data.get('locked', False):
            # If the request locks the database, we need to save the version
            updated = save_blastdb(updated, request.user, perform_lock=True)
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

class BlastDbHistory(mixins.RetrieveModelMixin, generics.GenericAPIView):
    '''
    Retrieve a list representing complete history changelog of a database.
    '''
    serializer_class = BlastDbHistoricalSerializer
    queryset = BlastDb.objects.all()
    permission_classes = [BlastDbEndpointPermission]

    @swagger_auto_schema(
        operation_summary='Show edit history of the database.',
        operation_description='Show a time log of the edit actions of the database including the removal and addition of sequences, as well as the editing of database metadata.',
        tags = [tag_blastdbs],
        responses={
            '200': 'Success.',
            '403': 'Insufficient permissions.',
            '404': 'BLAST database with the given ID does not exist',
            '500': 'Unexpected error.',
        }
    )
    def get(self, request, *args, **kwargs):
        db_id: str = kwargs.pop('pk', None)
        if db_id is None:
            return Response(status=status.HTTP_404_NOT_FOUND)
        try:
            database: BlastDb = BlastDb.objects.get(id=db_id)
            self.check_object_permissions(request, database)
        except BlastDb.DoesNotExist:
            return Response(status=status.HTTP_404_NOT_FOUND)

        return Response(self.get_serializer_class()(database.history.all(), many=True).data, status=status.HTTP_200_OK)

class BlastDbExport(generics.RetrieveAPIView):
    '''
    Export the reference library to JSON and different file formats.
    '''
    serializer_class = BlastDbExportSerializer
    renderer_classes = [JSONRenderer, BlastDbCSVRenderer, BlastDbTSVRenderer, BlastDbXMLRenderer, BlastDbFastaRenderer, BlastDbCompatibleRenderer]
    queryset = BlastDb.objects.all()
    permission_classes = [BlastDbEndpointPermission]
    lookup_url_kwarg = 'pk'

    def get_renderer_context(self):
        context = super().get_renderer_context()
        context['export_format'] = self.request.GET.get('export_format', '')
        return context
    
    @swagger_auto_schema(
        operation_summary='Export BLAST database to file.',
        operation_description='Export the sequences and/or metadata of the database to different file types.\n \
            The file type exported will depend on the `format` parameter in the request.\n\
            Exports can also be made to be compatible with popular bioinformatics tools such as SINTAX, DADA2 and QIIME2, by using the `export_format` parameter. Consult the documentation for all supported combinations of parameters.',
        tags = [tag_blastdbs],
        manual_parameters=[
            openapi.Parameter(
                name='format',
                description='Specify the file format. For .zip and .fasta files, additional settings can be set with the `export_format` parameter.',
                required=True,
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                enum=['csv', 'tsv', 'fasta', 'json', 'xml', 'zip']
            ),
            openapi.Parameter(
                name='export_format',
                description='Only applicable for `.zip` or `.fasta` file exports. This specifies the format of the file export to be compatible with a third-party bioinformatics tool. Leave blank for the default minimalistic fasta format.',
                required=False,
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                enum=['', 'qiime2', 'dada2tax', 'dada2sp', 'sintax', 'rdp', 'mothur']
            ),
        ],
        responses={
            '200': 'Export successful. Response should contain a file for download.',
            '403': 'Insufficient permissions.',
            '404': 'BLAST database with the given ID does not exist',
            '500': 'Unexpected error. Database with the given ID may also not exist',
            '502': 'Unexpected error.',
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        Export blastdb to compatible formats for taxonomic assignment.
        '''
        response = super(BlastDbExport, self).get(request, *args, **kwargs)
        db: BlastDb = self.get_object()

        # based on the media type file to be returned, specify attachment and file name
        if request.accepted_media_type.startswith('text/csv'):
            response['Content-Disposition'] = f'filename="barrel_database_{db.custom_name}_{db.id}.csv";'
        elif request.accepted_media_type.startswith('text/plain'):
            response['Content-Disposition'] = f'filename="barrel_database_{db.custom_name}_{db.id}.txt";'
        elif request.accepted_media_type.startswith('text/tsv'):
            response['Content-Disposition'] = f'filename="barrel_database_{db.custom_name}_{db.id}.tsv";'
        elif request.accepted_media_type.startswith('application/xml'):
            response['Content-Disposition'] = f'filename="barrel_database_{db.custom_name}_{db.id}.xml";'
        elif request.accepted_media_type.startswith('application/x-fasta'):
            response['Content-Disposition'] = f'filename="barrel_database_{db.custom_name}_{db.id}.fasta";'
        elif request.accepted_media_type.startswith('application/zip'):
            response['Content-Disposition'] = f'attachment; filename="barrel_database_{db.custom_name}_{db.id}.zip";'
        elif request.accepted_media_type.startswith('application/json'):
            return response
        else:
            # If the accepted media type cannot be chosen (ie content negotiation)
            # fails, then return a 406 response
            return Response(status=status.HTTP_406_NOT_ACCEPTABLE)
        return response

class BlastDbSummary(generics.GenericAPIView):
    lookup_url_kwarg = 'pk'
    queryset = BlastDb.objects.all()
    permission_classes = [BlastDbEndpointPermission]

    @swagger_auto_schema(
        operation_summary='Show some general statistics about the database\'s sequences.',
        operation_description='Show aggregate summaries of the database\'s sequences, including the number of sequences from each taxon, country of origin, and publication.',
        tags = [tag_blastdbs],
        
        responses={
            '200': openapi.Response(description='Success.', schema=openapi.Schema(
                title='Summary data',
                description='Number of sequences per each taxon, country and publication.',
                type=openapi.TYPE_OBJECT,
                properties={
                    **{
                        f'taxon_{i}__scientific_name': openapi.Schema(
                            title=f'Summary by {i}',
                            type=openapi.TYPE_ARRAY,
                            items=openapi.Items(
                                type=openapi.TYPE_OBJECT,
                                properties={
                                    f'taxon_{i}__scientific_name': openapi.Schema(title=f'Name of {i}', type=openapi.TYPE_STRING),
                                    f'taxon_{i}__id': openapi.Schema(title='Taxon id',type=openapi.TYPE_INTEGER),
                                    'count': openapi.Schema(title='Number of entries',type=openapi.TYPE_INTEGER),
                                }
                            )
                        )
                        for i in ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
                    },
                    'country': openapi.Schema(
                        title='Summary by country',
                        type=openapi.TYPE_ARRAY,
                        items=openapi.Items(
                            type=openapi.TYPE_OBJECT,
                            properties={
                                f'country': openapi.Schema(title=f'Full description of collection locale', type=openapi.TYPE_STRING),
                                **{
                                    f'{i}_name': openapi.Schema(title=f'Name of {i}, as extracted from country value.', type=openapi.TYPE_STRING)
                                    for i in ['country', 'region', 'locality']
                                },
                                'count': openapi.Schema(title='Number of entries', type=openapi.TYPE_INTEGER),
                            }
                        )
                    ),
                    'title': openapi.Schema(
                        title='Summary by title of reference',
                        type=openapi.TYPE_ARRAY,
                        items=openapi.Items(
                            type=openapi.TYPE_OBJECT,
                            properties={
                                f'title': openapi.Schema(title=f'Title of the reference resource, such as an article', type=openapi.TYPE_STRING),
                                f'journal': openapi.Schema(title=f'Journal name of the reference resource.', type=openapi.TYPE_STRING),
                                'count': openapi.Schema(title='Number of entries', type=openapi.TYPE_INTEGER),
                            }
                        )
                    ),
                }
            )),
            '403': 'Insufficient permissions.',
            '404': 'BLAST database with the given ID does not exist',
            '500': 'Unexpected error.',
        }
    )
    def get(self, request, *args, **kwargs):
        db = self.get_object()
        sequences : QuerySet[NuccoreSequence] = db.sequences.all()
        taxon_levels = ['taxon_superkingdom', 'taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus', 'taxon_species']
        fields = [f'{taxa}__scientific_name' for taxa in taxon_levels]
        fields.extend([f'annotations__annotation_type', 'country', 'title'])
        data = {}
        for field in fields:
            verbose_name: str = field
            # verbose_name = verbose_names.get(field, field)
            if field == 'country':
                data[verbose_name] = sequences.annotate(
                    colon_ind=StrIndex(field, Value(':'), output_field=PositiveSmallIntegerField()), 
                    comma_ind=StrIndex(field, Value(','), output_field=PositiveSmallIntegerField())
                ).annotate(
                    country_name=Case(
                        When(colon_ind__gt=1,then=Trim(Substr(field, Value(1), F('colon_ind') - Value(1)))),
                        # default=Value('Unknown', output_field=CharField()),
                        output_field=CharField()
                    ),
                    region_name=Case(
                        # has ; and ,
                        When(colon_ind__gt=1, comma_ind__gt=1, then=Trim(Substr(field, F('colon_ind') + 1, F('comma_ind') - Value(1) - F('colon_ind')))),
                        # has ; only
                        When(colon_ind__gt=1, then=Trim(Substr(field, F('colon_ind') + 1, Length(field)))),
                        default=Value('Unknown'),
                        output_field=CharField()
                    ),
                    locality_name=Case(
                        When(colon_ind__gt=1, comma_ind__gt=1, then=Trim(Substr(field, F('comma_ind') + 1, Length(field) - Value(1)))),
                        default=Value('Unknown'),
                        output_field=CharField()
                    )
                ).values('country', 'country_name', 'region_name', 'locality_name').annotate(count=Count('country'))
            elif field.startswith('taxon_'):
                data[verbose_name] = sequences.values(field, field.replace('__scientific_name', '__id', 1)).annotate(count=Count(field))
            elif field == 'title':
                data[verbose_name] = sequences.values('title', 'journal').annotate(count=Count('title'))
            else:
                data[verbose_name] = sequences.values(field).annotate(count=Count(field))
        return Response(data, status=status.HTTP_200_OK)

class BlastRunList(mixins.CreateModelMixin, generics.ListAPIView):
    '''
    List all BLAST runs. Used only for database administration purposes.
    '''

    serializer_class = BlastRunResultsShortSerializer

    def get_queryset(self):
        if isinstance(self.request.user, User):
            return BlastRun.objects.listable(self.request.user)
        else:
            return BlastRun.objects.none()

    @swagger_auto_schema(
        operation_summary='List all runs.',
        operation_description='Return a list of runs/jobs/queries saved by the API, including queued, running, and completed jobs.',
        tags = [tag_admin],
        responses={
            '200': openapi.Response(
                description='A list of all jobs saved.',
                schema = BlastRunResultsShortSerializer(many=True),
            )
        }
    )
    def get(self, request, *args, **kwargs):
        '''
        List all blast runs viewable
        '''
        return super().get(request, *args, **kwargs)

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
        request_body=BlastRunRunSerializer,
        tags = [tag_runs],
        responses={
            '400': openapi.Response(
                description='Bad request parameters. An accompanying message may specify the error with the request.',
                schema=BlastRunSerializer,
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

        # Validate the request
        serializer = BlastRunRunSerializer(data=request.data)
        if serializer.is_valid():
            query_string = serializer.data.get('query_sequence', '').strip()
            query_file = request.FILES.get('query_file', None)
            query_identifiers = serializer.data.get('query_identifiers', '').strip()
            query_identifiers_file = request.FILES.get('query_identifiers_file', None)
            job_name = serializer.data.get('job_name', '')
            create_db_tree = serializer.data.get('create_db_tree', False) 
            create_hit_tree = serializer.data.get('create_hit_tree', False)
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        
        # Bad request if no data is given
        if len(query_string) == 0 and query_file is None and \
            len(query_identifiers) == 0 and query_identifiers_file is None:
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

        print("Beginning to parse sequences")

        query_sequences : List = []
        unique_headers: set[str] = set([])
        # Parse any raw sequences given
        if len(query_string) > 0 or not query_file is None:
            if len(query_string) > 0:
                # try parsing with fasta first
                with io.StringIO(query_string) as query_string_io:
                    query_no = 1
                    query_records = SeqIO.parse(query_string_io, 'fasta')
                    query_record : SeqIO.SeqRecord
                    
                    for query_record in query_records:
                        seq = str(query_record.seq).strip()
                        if len(query_record.description) == 0:
                            query_record.description = f'query_sequence_{query_no}'
                            query_no = query_no + 1
                        query_sequences.append({
                            'definition': query_record.description,
                            'query_sequence': seq
                        })

                # if we found no fasta entries, then parse the whole thing as a sequence
                if len(query_sequences) == 0:
                    # replace all whitespace
                    query_string = query_string.strip().replace('\r', '').replace('\n', '').replace(' ', '')
                    query_sequences.append({
                        'definition': 'query_sequence', 
                        'query_sequence': query_string
                    })            
            if not query_file is None:
                # parse the file with Biopython
                try:
                    query_file_io = io.StringIO(query_file.file.read().decode('UTF-8'))
                    query_file_seqs = SeqIO.parse(query_file_io, 'fasta')
                except BaseException:
                    return Response({'message': 'The fasta file received was unable to be parsed. \
                                     Double check the text format and ensure it has correctly \
                                     formatted headers for each record.'}, 
                                     status = status.HTTP_400_BAD_REQUEST)

                parsed_entry: SeqIO.SeqRecord
                for parsed_entry in query_file_seqs:
                    query_sequences.append({
                        'definition': parsed_entry.description,
                        'query_sequence': str(parsed_entry.seq)
                    })

            # For user-provided sequences, clean the definition and identify any species labels in the definition
            for query_entry in query_sequences:
                # only keep the seqid 
                new_definition = query_entry['definition'].split(' ')[0]
                # parse the seqid to see if organism is specified 
                metadata = new_definition.split('\t')
                
                if len(metadata) >= 2:
                    # set the species name by replacing all non alpha-numeric characters and non-periods with spaces
                    query_entry['original_species_name'] = re.sub('[^0-9a-zA-Z.]+', ' ', metadata[1])
                else:
                    query_entry['original_species_name'] = ''
                new_definition = metadata[0]
                query_entry['definition'] = new_definition

                if new_definition in unique_headers:
                    return Response({'message': f"The sequence identifier in the header '{new_definition}'\
                            is not unique. All submitted sequences must have unique sequence identifiers \
                            (all text before the first whitespace)."}, status = status.HTTP_400_BAD_REQUEST) 
                unique_headers.add(new_definition)
        
        # Parse any identifiers given
        if len(query_identifiers) > 0 or not query_identifiers_file is None:
            ids = []
            # Add identifiers from text
            if len(query_identifiers) > 0:
                ids = query_identifiers.replace('\r\n', '\n').strip().split('\n')
                ids.extend([id.strip() for id in ids])

            # Add identifiers from file upload
            if not query_identifiers_file is None:
                query_file_io = io.StringIO(query_identifiers_file.file.read().decode('UTF-8'))
                file_ids = query_file_io.readlines()
                ids.extend([id.strip() for id in file_ids])

            # Filter out any empty ids
            ids = [id for id in ids if len(id) > 0]
                            
            try:
                # Retrieve sequences from GenBank based on the identifiers provided
                identifier_sequences = retrieve_gb(accession_numbers=ids, raise_if_missing=False)
            except ValueError:
                return Response({'message': f'Query identifiers could not be parsed from the query identifier text or \
                    file.'}, 
                    status=status.HTTP_502_BAD_GATEWAY)
            except InsufficientAccessionData as exc:
                return Response({'message': f'One or more of the following identifiers could not be retrieved from \
                    GenBank: {exc.missing_accessions}'}, status=status.HTTP_400_BAD_REQUEST)
            except GenBankConnectionError as exc:
                return Response({'message': f'Encountered unexpected error connecting to GenBank'}, 
                    status=status.HTTP_502_BAD_GATEWAY)
            except AccessionLimitExceeded as exc:
                return Response({'message': f'Too many accessions submitted in a single operation. \
                                 Consider splitting this query into different jobs. Absolute maximum \
                                 is {exc.max_accessions}.'}, status=status.HTTP_400_BAD_REQUEST)
            else:
                if identifier_sequences is None:
                    return Response({'message': f'No identifiers could not be retrieved from \
                    GenBank.'}, status=status.HTTP_400_BAD_REQUEST)
                
                # Find existing database entries if they match the query sequences in accession.version
                existing: QuerySet[NuccoreSequence] = NuccoreSequence.objects.filter(owner_database=odb)

                for identifier_sequence in identifier_sequences:
                    seqid = identifier_sequence.get('version')

                    # Do not allow query sequences with accession.versions already in the database.
                    # TODO: Fix bug where query can be identifiers in the reference database
                    if existing.filter(version=seqid).exists():
                        continue 

                    if seqid in unique_headers:
                        return Response({'message': f"The sequence identifier in the header '{seqid}'\
                            is not unique. All submitted sequences must have unique sequence identifiers \
                            (all text before the first space)."}, status = status.HTTP_400_BAD_REQUEST)
                    
                    query_sequences.append({
                        'definition': seqid,
                        'query_sequence': identifier_sequence.get('dna_sequence'),
                        'original_species_name': identifier_sequence.get('organism')
                    })

        # Raise an error if nothing could be retrieved based on sequences or identifiers given
        if len(query_sequences) == 0:
            return Response({'message': 'Could not successfully parse any query sequences or identifiers in the request to query with. Check the formatting and identifiers provided. If identifiers were given corresponding to records already in the database, they were ignored.'}, status = status.HTTP_400_BAD_REQUEST)
        
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

        # Generate query fasta file to blast with
        print('Generating query fasta file ...')
        query_file = results_path + '/query.fasta'
        blast_tmp = open(query_file, 'w')
        
        for query_entry in query_sequences:
            defn = query_entry.get('definition')
            seq = query_entry.get('query_sequence')
            try:
                verify_header(defn, seq)
                verify_dna(defn, seq)
            except ValueError as value_error:
                return Response({'message': value_error.args}, status = status.HTTP_400_BAD_REQUEST) 
            # Write the header and sequence to fasta 
            blast_tmp.write(f'>{defn}\n{seq}\n')

        blast_tmp.close()

        # Perform blast search
        print('Running BLAST search ...')

        ncbi_blast_version = 'ncbi-blast-2.12.0+'
        run_details = BlastRun(id = results_uuid, db_used = odb, job_name = job_name, blast_version = ncbi_blast_version, errors = '', status=BlastRun.JobStatus.QUEUED, start_time = None, end_time = None, error_time = None, create_hit_tree = create_hit_tree, create_db_tree = create_db_tree)
        run_details.save()

        # make query sequence objects
        all_query_sequences = [BlastQuerySequence(**query_sequence, owner_run = run_details) for query_sequence in query_sequences]
        BlastQuerySequence.objects.bulk_create(all_query_sequences)

        # If we are performing alignment, create the fasta file
        if create_db_tree or create_hit_tree:
            aln_tmp = open(results_path + '/alignment_query.fasta', 'w')
            for query_seq in all_query_sequences:
                aln_tmp.write(f'>{query_seq.write_tree_identifier()}\n')
                aln_tmp.write(f'{query_seq.query_sequence}\n')
            aln_tmp.close()

        run_details_id = run_details.id

        run_id_str = str(run_details_id)

        # How to chain for Docker-hosted RabbitMQ 
        res = chain(run_blast_command.s(ncbi_blast_version, str(odb.id), run_id_str)).apply_async(queue='BarcodeQueue.fifo')

        # How to chain for Amazon SQS
        # message_properties = {
        #     'MessageGroupId': str(run_details_id),
        #     'MessageDeduplicationId': str(run_details_id),
        # }
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

class BlastRunDetail(mixins.DestroyModelMixin, generics.RetrieveAPIView):
    '''
    View the results of any given run
    '''

    queryset = BlastRun.objects.all()
    serializer_class = BlastRunResultsShortSerializer
    permission_classes = [AllowAny,]
    
    def get_object(self):
        db_primary_key = self.kwargs.get('pk')
        try:
            run = BlastRun.objects.get(id = db_primary_key)
        except BlastRun.DoesNotExist:
            raise Http404()
        else:
            return run

    @swagger_auto_schema(
        operation_summary='Get run results.',
        operation_description='Get a condensed summary of the run information. The results returned are complete once the status is "FIN" (i.e. run is finished), and may be incomplete, empty or missing otherwise. The run information includes the submitted sequences and Newick strings for any phylogenetic trees. To browse the hits, separately access each query sequence individually using its identifier.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Run information successfully retrieved.',
                schema=BlastRunResultsShortSerializer(),
            ),
            '404': 'Run with the ID does not exist',
        }
    )
    def get(self, request, *args, **kwargs):
        return super().get(request, *args, **kwargs)

class BlastRunQueryList(generics.ListAPIView):
    '''
    Return a paginated list of query sequences for a given run, \
        minus with data of all the resulting hits.
    '''
    pagination_class = BlastRunQueryPagination
    serializer_class = BlastQuerySequenceShortSerializer
    filter_backends = [OrderingFilter]
    ordering_fields = '__all__'
    ordering = ['definition']

    def get_queryset(self):
        pk: str = self.kwargs.get('pk')
        try:
            run = BlastRun.objects.get(id=pk)
        except BlastRun.DoesNotExist:
            raise Http404()
        return run.queries.all()

    @swagger_auto_schema(
        operation_summary='Get query sequences submitted',
        operation_description='Get a list of the query sequences submitted in a run, including the submitted sequences and tentative species identity for each sequence. To read the hit information, retrieve each individual query sequence separately using its id.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Information successfully retrieved.',
                schema=BlastQuerySequenceShortSerializer,
            ),
            '404': 'Run with the ID does not exist',
        }
    )
    def get(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)

class BlastQuerySequenceDetail(generics.RetrieveAPIView):
    '''
    Return a query sequence with hits.
    '''
    serializer_class = BlastQuerySequenceSerializer

    def get_queryset(self):
        return BlastQuerySequence.objects.all()

    def get_object(self):
        pk: str = self.kwargs.get('pk')
        try:
            run = BlastRun.objects.get(id=pk)
        except BlastQuerySequence.DoesNotExist:
            raise Http404()

        query: str = self.kwargs.get('query')
        try:
            obj = BlastQuerySequence.objects.get(owner_run=run, id=query)
        except BlastQuerySequence.DoesNotExist:
            raise Http404()
        return obj

    @swagger_auto_schema(
        operation_summary='Get BLAST query sequence information',
        operation_description='Get information about a query sequence submitted for a BLAST run, including hits. \
            All hits are returned in a single list within the returned object, so the volume of the response may be large. \
                To retrieve the hits, it is recommended that the `runs/{id}/query/{query}/hits` endpoint is used, since it will paginate \
                    results across different requests.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Information successfully retrieved.',
                schema=BlastQuerySequenceSerializer(many=True),
            ),
            '404': 'Run with the ID does not exist',
        }
    )
    def get(self, request, *args, **kwargs):
        return super().get(request, *args, **kwargs)

class BlastQueryHitList(generics.ListAPIView):
    '''
    Return a paginated list of all hits for a given query sequence
    '''
    pagination_class = BlastRunHitPagination
    filter_backends = [CustomHitOrderingFilter]
    ordering_fields = '__all__'
    ordering = ['query_accession_version']
    serializer_class = HitSerializer

    def get_queryset(self):
        pk: str = self.kwargs.get('pk')
        try:
            run = BlastRun.objects.get(id=pk)
        except BlastRun.DoesNotExist:
            raise Http404()

        query: str = self.kwargs.get('query')
        try:
            obj = BlastQuerySequence.objects.get(owner_run=run, id=query)
        except BlastQuerySequence.DoesNotExist:
            raise Http404()
        return obj.hits.all() 

    @swagger_auto_schema(
        operation_summary='Get BLAST hit information',
        operation_description='Get a list of BLAST hits for the corresponding query sequence, with each hit including all calculated statistics such as alignment length, percent identity, bit score, etc.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Information successfully retrieved.',
                schema=HitSerializer(many=True),
            ),
            '404': 'Run with the ID does not exist',
        }
    )
    def get(self, request, *args, **kwargs):
        return super().get(request, *args, **kwargs)

class BlastRunInputDownload(generics.GenericAPIView):
    '''
    Return the sequences submitted to a run as a .fasta file
    '''
    serializer_class = BlastRunSerializer
    permission_classes = (AllowAny,)
    renderer_classes = [BlastRunFastaRenderer]
    queryset = BlastRun.objects.all()
    
    @swagger_auto_schema(
        operation_summary='Get query sequence fasta file.',
        operation_description='Returns the original query sequences from the run submission, formatted as a FASTA file attachment.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Query sequences retrieved successfully and a fasta file attachment is returned.',
                schema=openapi.Schema(
                    type=openapi.TYPE_FILE,
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
        
        response = Response(self.get_serializer_class()(run).data, status=status.HTTP_200_OK)
        if request.accepted_media_type.startswith('application/x-fasta'):
            response['Content-Disposition'] = f'attachment; filename="barrel_{run.id}.input.fasta";'
        else:
            return Response(status=status.HTTP_406_NOT_ACCEPTABLE)
        return response

class BlastRunDetailDownload(generics.GenericAPIView):
    '''
    Download the run results in either .txt or .csv format
    '''
    queryset = BlastRun.objects.all()
    serializer_class = BlastRunSerializer
    renderer_classes = [BlastRunHitsCSVRenderer, BlastRunHitsTSVRenderer, BlastRunHitsTxtRenderer]
    
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
        if request.accepted_media_type.startswith('text/csv'):
            response['Content-Disposition'] = f'attachment; filename="barrel_{run.id}.hits.csv";'
        elif request.accepted_media_type.startswith('text/tsv'):
            response['Content-Disposition'] = f'attachment; filename="barrel_{run.id}.hits.tsv";'
        elif request.accepted_media_type.startswith('text/plain'):
            response['Content-Disposition'] = f'attachment; filename="barrel_{run.id}.hits.txt";'

        return response

class BlastRunTaxonomyDownload(generics.GenericAPIView):
    '''
    Download the taxonomic classification results.
    '''
    serializer_class = BlastRunTaxonomySerializer
    renderer_classes = [JSONRenderer, BlastRunTaxonomyCSVRenderer, BlastRunTaxonomyTSVRenderer]

    def get_queryset(self):
        pk = self.kwargs.get('pk', '')
        # Get all query sequences under the run specified by pk
        try:
            run: BlastRun = BlastRun.objects.get(id=pk)
        except BlastRun.DoesNotExist:
            raise Http404()
        return run
    
    @swagger_auto_schema(
        operation_summary='Get taxonomic classification results',
        operation_description='Returns taxonomic assignments/classifications made by the application according to BLAST hit results. Whether taxonomic assignments are made for a run is determined by the run parameters.',
        tags = [tag_runs],
        responses={
            '200': openapi.Response(
                description='Returns a file.',
                schema = openapi.Schema(
                    type= openapi.TYPE_STRING,
                    format='binary'
                ),
            ),
            '404': 'Run information or taxonomic assignment data corresponding to the specified ID does not exist'
        },
    )
    def get(self, request, *args, **kwargs):
        best_hits = self.get_queryset()
        pk = self.kwargs.get('pk', '')
        serializer = self.get_serializer(best_hits)
        response = Response(serializer.data, status=status.HTTP_200_OK)
        if request.accepted_media_type.startswith('text/csv'):
            response['Content-Disposition'] = f'attachment; filename="barrel_{pk}.taxonomy.csv";'
        elif request.accepted_media_type.startswith('text/tsv'):
            response['Content-Disposition'] = f'attachment; filename="barrel_{pk}.taxonomy.tsv";'

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
