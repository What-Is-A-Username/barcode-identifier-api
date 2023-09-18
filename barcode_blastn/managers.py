from typing import Callable, Generic, List, TypeVar, Union
from django.db import models
from django.contrib.auth.models import (AbstractBaseUser, AbstractUser,
                                        AnonymousUser, User)
from django.http.request import HttpRequest
from rest_framework import permissions

from barcode_blastn.database_permissions import DatabasePermissions

