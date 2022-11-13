from rest_framework import permissions

class IsAdminOrReadOnly(permissions.BasePermission):
    def has_permission(self, request, view):
        # allow GET, HEAD, OPTIONS to be used by anyone
        if request.method in permissions.SAFE_METHODS:
            return True
        
        # only allow other requests if the user is admin
        return request.user.is_staff
