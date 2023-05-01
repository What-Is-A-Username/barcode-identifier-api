from rest_framework import permissions
from django.contrib.auth.models import User
from barcode_blastn.models import BlastDb, DatabaseShare

class IsAdminOrReadOnly(permissions.BasePermission):
    def has_permission(self, request, view):
        # allow GET, HEAD, OPTIONS to be used by anyone
        if request.method in permissions.SAFE_METHODS:
            return True
        
        # only allow other requests if the user is admin
        return request.user.is_staff

class BlastDbAccessPermission(permissions.BasePermission):
    message = 'Insufficient permissions for this action.'

    def has_permission(self, request, view):
        return super().has_permission(request, view)

    def has_object_permission(self, request, view, obj: BlastDb):
        '''
        Restrict view and edit access based on permissions for each BlastDb instance
        '''       
        if not obj:
            user: User = request.user
            if not user.is_authenticated:
                return False 
            else:
                return user.is_staff or user.is_superuser
        else:
            if request.method in permissions.SAFE_METHODS:
                return DatabaseShare.has_view_and_run_permission(request.user, obj)
            elif request.method in ['PUT', 'PATCH']:
                return DatabaseShare.has_edit_permission(request.user, obj)
            else:
                return DatabaseShare.has_delete_permission(request.user, obj)