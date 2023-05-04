from typing import Callable, Generic, List, TypeVar, Union
from django.db import models
from django.contrib.auth.models import (AbstractBaseUser, AbstractUser,
                                        AnonymousUser, User)
from django.http.request import HttpRequest
from rest_framework import permissions

from barcode_blastn.database_permissions import DatabasePermissions

class BlastDbManager(models.Manager):
    def editable(self, user: User):
        '''
        Return a queryset of BlastDb objects that are editable by the given user
        '''
        if not user.is_authenticated:
            return super().none()
        elif user.is_superuser:
            return super().get_queryset()
        else:
            return super().get_queryset().filter(
                models.Q(owner=user) |
                models.Q(shares=user, databaseshare__perms__in=[DatabasePermissions.CAN_EDIT_DB])
            )