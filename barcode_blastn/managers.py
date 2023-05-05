from typing import Callable, Generic, List, TypeVar, Union
from django.db import models
from django.contrib.auth.models import (AbstractBaseUser, AbstractUser,
                                        AnonymousUser, User)
from django.http.request import HttpRequest
from rest_framework import permissions

from barcode_blastn.database_permissions import DatabasePermissions

class BlastDbManager(models.Manager):
    '''
    Model manager for the BlastDb class.
    '''
    def editable(self, user: User):
        '''
        Return a queryset of BlastDb objects that are editable by the given user.

        Returns an empty queryset if user is not authenticated, the full queryset
        if user is a superuser, and a queryset of all BlastDbs editable by the
        user
        '''
        if not user.is_authenticated:
            # public users can never edit a database
            return super().none()
        elif user.is_superuser:
            # give access to all databases if superuser
            return super().get_queryset().all()
        else:
            return super().get_queryset().filter(
                # include databases owned by user
                models.Q(owner=user) |
                # include those with access explicitly given to 
                models.Q(shares=user, databaseshare__perms__in=[DatabasePermissions.CAN_EDIT_DB])
            )
        
    def viewable(self, user: User):
        '''
        Retrieve a set of databases that the user can view
        '''
        if not user.is_authenticated:
            # only return public databases for public users
            return super().get_queryset().filter(public=True)
        elif user.is_superuser:
            # give access to all databases if superuser
            return super().get_queryset().all()
        else:
            return super().get_queryset().exclude(
                # exclude databases that are public but explicitly denied to user
                models.Q(public=True) &
                models.Q(shares=user, 
                         databaseshare__perms=DatabasePermissions.DENY_ACCESS)
            ).filter(
                # include databases owned by user
                models.Q(owner=user) |
                # include those that are public
                models.Q(public=True) |
                # include those with access explicitly given to 
                models.Q(shares=user,
                         databaseshare__perms__in=[
                            DatabasePermissions.CAN_VIEW_DB, DatabasePermissions.CAN_EDIT_DB, DatabasePermissions.CAN_RUN_DB
                         ])
            )

    def runnable(self, user: User):
        '''
        Retrieve a set of databases that the user can run on
        '''
        # Currently, if a user can view the database, they can run it
        return self.viewable(user)

    def deletable(self, user: User):
        '''
        Retrieve a database set that is deletable by the user.
        '''
        if not user.is_authenticated:
            # public users cannot delete any databases
            return super().none()
        elif user.is_superuser:
            # A superuser can delete any database
            return super().get_queryset().all()
        else:
            return super().get_queryset().filter(owner=user)
        