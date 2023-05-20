from abc import ABC, abstractmethod
from typing import Callable, Generic, List, Type, TypeVar, Union, get_args
from django.db import models
from django.contrib.auth.models import (AbstractBaseUser, AbstractUser,
                                        AnonymousUser, User)
from django.http.request import HttpRequest
from rest_framework import permissions

from barcode_blastn.models import (BlastDb, BlastRun,
                                   DatabaseShare, Hit, Library, NuccoreSequence)

from barcode_blastn.database_permissions import DatabasePermissions                                   

class IsAdminOrReadOnly(permissions.BasePermission):
    def has_permission(self, request, view):
        # allow GET, HEAD, OPTIONS to be used by anyone
        if request.method in permissions.SAFE_METHODS:
            return True
        
        # only allow other requests if the user is admin
        return request.user.is_staff

T = TypeVar('T', bound=models.Model)
class CustomPermissions(Generic[T], ABC):
    '''
    Determine user permissions for a user-defined object
    '''
    @staticmethod
    @abstractmethod
    def has_module_permission(user: Union[AbstractBaseUser, AnonymousUser]) -> bool:
        """
        Return True if the given request has any permission in the given
        app label.

        Can be overridden by the user in subclasses. In such case it should
        return True if the given request has permission to view the module on
        the admin index page and access the module's index page. Overriding it
        does not restrict access to the add, change or delete views. Use
        `ModelAdmin.has_(add|change|delete)_permission` for that.
        """
        pass

    @staticmethod
    @abstractmethod
    def has_add_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[T, None]) -> bool:
        """
        Return True if the given request has permission to add an object.
        Can be overridden by the user in subclasses.

        It should return True if the given user has permission to add the given 
        obj instance. If `obj` is None, it should return True if the user
        has permission to add any object of the given type.
        """
        pass

    @staticmethod
    @abstractmethod
    def has_view_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[T, None]) -> bool:
        """
        Return True if the given request has permission to view the given
        Django model instance. The default implementation doesn't examine the
        `obj` parameter.

        If overridden by the user in subclasses, it should return True if the
        given request has permission to view the `obj` model instance. If `obj`
        is None, it should return True if the request has permission to view
        any object of the given type.
        """
        pass

    @staticmethod
    @abstractmethod
    def has_change_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[T, None]) -> bool:
        """
        Return True if the given request has permission to change the given
        Django model instance, the default implementation doesn't examine the
        `obj` parameter.

        Can be overridden by the user in subclasses. In such case it should
        return True if the given request has permission to change the `obj`
        model instance. If `obj` is None, this should return True if the given
        request has permission to change *any* object of the given type.
        """
        pass

    @staticmethod
    @abstractmethod
    def has_delete_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[T, None]) -> bool:
        """
        Return True if the given request has permission to change the given
        Django model instance, the default implementation doesn't examine the
        `obj` parameter.

        Can be overridden by the user in subclasses. In such case it should
        return True if the given request has permission to delete the `obj`
        model instance. If `obj` is None, this should return True if the given
        request has permission to delete *any* object of the given type.
        """
        pass

class LibrarySharePermissions(CustomPermissions[Library]):
    '''
    Given a user, an optionally an object, determine whether an action is permissible.
    '''

    @staticmethod 
    def has_explicit_permission_given(user: AbstractUser, database: Library, acceptable_permissions: List[str], prohibited_permissions: List[str] = [DatabasePermissions.DENY_ACCESS]):
        '''
        Check if a user has the required level of permissions.

        Returns True in any of the following cases: 
        - user's permissions matches one of the options in necessary_permissions and does NOT match any of the options specified by prohibited_permissions

        Importantly, it returns False if:
        - user is not an instance that can be stored in a DatabaseShare permission
        - user does not have any permissions assigned for the database.
        '''
        storedperms : DatabaseShare
        if not isinstance(user, User):
            return False
        try:
            storedperms = DatabaseShare.objects.get(database=database, grantee=user)
        except DatabaseShare.DoesNotExist:
            return False    # deny if no explicit permission given 
        return storedperms.permission_level in acceptable_permissions and not storedperms.permission_level in prohibited_permissions # return true if the current permission matches with any of those given

    @staticmethod 
    def has_module_permission(user: Union[AbstractBaseUser, AnonymousUser]) -> bool:
        return isinstance(user, AbstractUser) and (user.is_staff or user.is_superuser)

    @staticmethod
    def has_add_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[Library, None]) -> bool:
        return isinstance(user, AbstractUser) and user.is_superuser

    @staticmethod
    def has_run_permission(user: Union[AbstractBaseUser, AnonymousUser], library: Union[Library, None]) -> bool:
        '''
        Check if the user has permission to run BLAST on the database.
        '''
        if library is None: # for generic models, give access by default
            return True
        elif not user.is_authenticated:
            # allow all public databases to be run
            return library.public
        else:
            if user == library.owner:
                return True 
            elif not isinstance(user, AbstractUser):
                return library.public
            elif user.is_superuser:
                return True 
            elif library.public:
                # signed in users can query public databases so long as they aren't explicitly denied
                return not LibrarySharePermissions.has_explicit_permission_given(
                    user,
                    library,
                    acceptable_permissions=[DatabasePermissions.DENY_ACCESS],
                    prohibited_permissions=[]
                ) 
            return LibrarySharePermissions.has_explicit_permission_given(
                user, 
                library, 
                acceptable_permissions=[DatabasePermissions.CAN_RUN_DB, DatabasePermissions.CAN_EDIT_DB],
                prohibited_permissions=[DatabasePermissions.DENY_ACCESS]) # check if user has these permissions

    @staticmethod
    def has_view_permission(user: Union[AbstractBaseUser, AnonymousUser], database: Union[Library, None]) -> bool:
        if database is None: # for generic models, give access by default
            return True
        elif not user or not user.is_authenticated:
            return database.public
        else:
            if user == database.owner:
                return True
            elif not isinstance(user, AbstractUser):
                return False
            elif database.public:
                return not LibrarySharePermissions.has_explicit_permission_given(
                    user,
                    database,
                    acceptable_permissions=[DatabasePermissions.DENY_ACCESS],
                    prohibited_permissions=[]
                )
            elif (user.is_superuser):
                return True 
            else:
                return LibrarySharePermissions.has_explicit_permission_given(
                    user, 
                    database, 
                    acceptable_permissions=[DatabasePermissions.CAN_VIEW_DB, DatabasePermissions.CAN_EDIT_DB, DatabasePermissions.CAN_RUN_DB],
                    prohibited_permissions=[DatabasePermissions.DENY_ACCESS]) # check if user has these permissions

    @staticmethod
    def has_change_permission(user: Union[AbstractBaseUser, AnonymousUser], database: Union[Library, None]) -> bool:
        if not user.is_authenticated: # user account needed for editing
            return False 
        elif database is None: # for generic models, give access so long as there is a deletable database
            return isinstance(user, User) and Library.objects.editable(user).exists()
        else:
            if not isinstance(user, User):
                return False
            elif database.owner == user or user.is_superuser:
                return True 
            return LibrarySharePermissions.has_explicit_permission_given(user, 
                database, 
                acceptable_permissions=[DatabasePermissions.CAN_EDIT_DB],
                prohibited_permissions=[DatabasePermissions.DENY_ACCESS]) # check if user has these permissions

    @staticmethod
    def has_delete_permission(user: Union[AbstractBaseUser, AnonymousUser], database: Union[Library, None]) -> bool:
        '''
        Return whether the user can delete a reference library
        '''
        return isinstance(user, User) and user.is_superuser

class DatabaseSharePermissions(CustomPermissions[BlastDb]):
    @staticmethod 
    def has_explicit_permission_given(user: AbstractUser, database: BlastDb, acceptable_permissions: List[str], prohibited_permissions: List[str] = [DatabasePermissions.DENY_ACCESS]):
        pass

    @staticmethod
    def has_module_permission(user: Union[AbstractBaseUser, AnonymousUser]) -> bool:
        return LibrarySharePermissions.has_module_permission(user)

    @staticmethod
    def has_add_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastDb, None]) -> bool:
        if obj is None:
            return isinstance(user, AbstractUser) and user.is_authenticated
        else:
            return LibrarySharePermissions.has_add_permission(user, obj.library)

    @staticmethod
    def has_view_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastDb, None]) -> bool:
        if obj is None:
            return True
        else:
            return LibrarySharePermissions.has_view_permission(user, obj.library)
 
    @staticmethod
    def has_run_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastDb, None]) -> bool:
        if obj is None:
            return isinstance(user, AbstractUser) and user.is_authenticated
        else:
            return LibrarySharePermissions.has_run_permission(user, obj.library)

    @staticmethod
    def has_change_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastDb, None]) -> bool:
        if obj is None:
            return isinstance(user, AbstractUser) and user.is_authenticated
        else:
            return LibrarySharePermissions.has_change_permission(user, obj.library)

    @staticmethod
    def has_delete_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastDb, None]) -> bool:
        if obj is None:
            return isinstance(user, AbstractUser) and user.is_authenticated
        else:
            return LibrarySharePermissions.has_delete_permission(user, obj.library)

class RunSharePermissions(CustomPermissions[BlastRun]):
    '''
    Determine user permissions for a BlastRun
    '''

    @staticmethod
    def defer_to_database( 
            permission_function: Callable[[Union[AbstractBaseUser, AnonymousUser], Union[BlastDb, None]], bool], 
            user: Union[AbstractBaseUser, AnonymousUser], 
            obj: Union[BlastRun, None]
        ):
        '''
        Return True if request has permission to the database associated with the given blastrun 
        '''
        if obj is None:
            return permission_function(user, None)
        else:
            return permission_function(user, obj.db_used)

    @staticmethod
    def has_module_permission(user: Union[AbstractBaseUser, AnonymousUser]) -> bool:
        return DatabaseSharePermissions.has_module_permission(user)

    @staticmethod
    def has_add_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastRun, None]) -> bool:
        if obj is None:
            # runs cannot be added in general
            return False 
        else:
            # runs can only be added if there is run access to that database
            return RunSharePermissions.defer_to_database(DatabaseSharePermissions.has_run_permission, user, obj)

    @staticmethod
    def has_view_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastRun, None]) -> bool:
        # runs can only be viewed if the user can edit the database
        return RunSharePermissions.defer_to_database(DatabaseSharePermissions.has_change_permission, user, obj)

    @staticmethod
    def has_change_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastRun, None]) -> bool:
        # prohibit runs from being modified after creation
        return False

    @staticmethod
    def has_delete_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[BlastRun, None]) -> bool:
        # runs can only be deleted if the user can delete the database
        return RunSharePermissions.defer_to_database(DatabaseSharePermissions.has_delete_permission, user, obj)

class NuccoreSharePermission(CustomPermissions[NuccoreSequence]):
    @staticmethod
    def defer_to_database( 
            permission_function: Callable[[Union[AbstractBaseUser, AnonymousUser], Union[BlastDb, None]], bool], 
            user: Union[AbstractBaseUser, AnonymousUser], 
            obj: Union[NuccoreSequence, None]
        ):
        '''
        Return True if request has permission to the database associated with the given blastrun 
        '''
        if obj is None:
            return permission_function(user, None)
        else:
            return permission_function(user, obj.owner_database)

    @staticmethod
    def has_module_permission(user: Union[AbstractBaseUser, AnonymousUser]) -> bool:
        return DatabaseSharePermissions.has_module_permission(user)

    @staticmethod
    def has_add_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[NuccoreSequence, None]) -> bool:
        '''
        Check whether the given sequence `obj` can be added to its owner_database.

        If obj is None, return True if any sequence can be added at all to any database. Otherwise, return True if user has permissions to add to database, database is not locked, and database doesn't already have the accession_number.
        '''
        # sequences can only be added if the user can edit the database, database not locked and no existing entry matches the accession number
        if obj is None:
            if isinstance(user, User):
                return BlastDb.objects.editable(user).filter(locked=False).exists()
            else:
                return False
        else:
            return NuccoreSharePermission.defer_to_database(
                DatabaseSharePermissions.has_change_permission, user, obj) \
                and not obj.owner_database.locked \
                and not obj.owner_database.sequences.filter(accession_number=obj.accession_number).exists()

    @staticmethod
    def has_view_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[NuccoreSequence, None]) -> bool:
        # runs can only be viewed if the user can view the database
        return NuccoreSharePermission.defer_to_database(DatabaseSharePermissions.has_view_permission, user, obj)

    @staticmethod
    def has_change_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[NuccoreSequence, None]) -> bool:
        # sequences can only be edited if the user can edit the database and database not locked
        if obj is None:
            return False
        else:
            return NuccoreSharePermission.defer_to_database(
                DatabaseSharePermissions.has_change_permission, user, obj) \
                and not obj.owner_database.locked

    @staticmethod
    def has_delete_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[NuccoreSequence, None]) -> bool:
        # sequences can only be deleted if the user can edit the database and database not locked
        if obj is None:
            return NuccoreSharePermission.defer_to_database(
                DatabaseSharePermissions.has_change_permission, user, obj)
        elif isinstance(user, User) and user.is_superuser:
            return NuccoreSharePermission.defer_to_database(
                DatabaseSharePermissions.has_delete_permission, user, obj) 
        else:
            return NuccoreSharePermission.defer_to_database(
                DatabaseSharePermissions.has_change_permission, user, obj) \
                and not obj.owner_database.locked

class HitSharePermission(CustomPermissions[Hit]):
    
    @staticmethod
    def defer_to_run( 
            permission_function: Callable[[Union[AbstractBaseUser, AnonymousUser], Union[BlastRun, None]], bool], 
            user: Union[AbstractBaseUser, AnonymousUser], 
            obj: Union[Hit, None]
        ):
        '''
        Return True if request has permission to the database associated with the given blastrun 
        '''
        if obj is None:
            return permission_function(user, None)
        else:
            return permission_function(user, obj.owner_run)

    @staticmethod
    def has_module_permission(user: Union[AbstractBaseUser, AnonymousUser]) -> bool:
        return DatabaseSharePermissions.has_module_permission(user)

    @staticmethod
    def has_add_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[Hit, None]) -> bool:
        return HitSharePermission.defer_to_run(RunSharePermissions.has_add_permission, user, obj)

    @staticmethod
    def has_view_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[Hit, None]) -> bool:
        return HitSharePermission.defer_to_run(RunSharePermissions.has_view_permission, user, obj)

    @staticmethod
    def has_change_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[Hit, None]) -> bool:
        return HitSharePermission.defer_to_run(RunSharePermissions.has_change_permission, user, obj)

    @staticmethod
    def has_delete_permission(user: Union[AbstractBaseUser, AnonymousUser], obj: Union[Hit, None]) -> bool:
        return HitSharePermission.defer_to_run(RunSharePermissions.has_delete_permission, user, obj)

# Define a generic type that must be a subclass of CustomPermissions
SharePermission = TypeVar('SharePermission', bound=CustomPermissions) 
EndpointModel = TypeVar('EndpointModel', bound=models.Model)   
class BaseEndpointPermission(Generic[SharePermission, EndpointModel], permissions.BasePermission):
    '''
    Permission class used for endpoints. Subclasses 
    '''

    @classmethod
    def get_share_permission_type(cls) -> Type[SharePermission]:
        '''
        Retrieve the type of SharePermission used
        '''
        # Requires Python >= 3.8
        return get_args(cls.__orig_bases__[0])[0]

    def determine_obj_permission(self, share_permission_type: Type[SharePermission], request, view, obj: Union[Type[EndpointModel], None]):
        '''
        Restrict view and edit access based on permissions determined by CustomPermissions instance
        '''       
        if request.method in permissions.SAFE_METHODS:
            return share_permission_type.has_view_permission(request.user, obj)
            # return share_permission_type.has_view_permission(request.user, obj)
        elif request.method in ['PUT', 'PATCH']:
            return share_permission_type.has_change_permission(request.user, obj)
        elif request.method in ['POST']:
            return share_permission_type.has_add_permission(request.user, obj)
        elif request.method == 'DELETE':
            return share_permission_type.has_delete_permission(request.user, obj) 
        else:
            return share_permission_type.has_view_permission(request.user, obj)
    
    def has_permission(self, request, view):
        return super().has_permission(request, view)

    def has_object_permission(self, request, view, obj):
        return self.determine_obj_permission(self.get_share_permission_type(), request, view, obj)
     
class NuccoreSequenceEndpointPermission(BaseEndpointPermission[NuccoreSharePermission, NuccoreSequence]):
    '''
    Permission class used for NuccoreSequence endpoints
    '''
    message = 'Insufficient permissions for this action'

class BlastDbEndpointPermission(BaseEndpointPermission[DatabaseSharePermissions, BlastDb]):
    '''
    Permission class used for BlastDb endpoints
    '''
    message = 'Insufficient permissions for this action.'

class BlastRunEndpointPermission(BaseEndpointPermission[RunSharePermissions, BlastRun]):
    '''
    Permission class used for BlastRun endpoints
    '''
    message = 'Insufficient permissions for this action.'

class LibraryEndpointPermission(BaseEndpointPermission[LibrarySharePermissions, Library]):
    '''
    Permission class used for Library endpoints
    '''
    message = 'Insufficient permissions for this action.'