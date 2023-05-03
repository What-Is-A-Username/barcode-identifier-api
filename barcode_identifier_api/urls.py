"""barcode_identifier_api URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.urls import path, include
from django.contrib import admin

admin.site.site_title = 'Admin'
admin.site.site_header = 'Barcode Identifier API Adminstration'
admin.site.index_title = 'Database Adminstration Dashboard'

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('barcode_blastn.urls')),
    path('', include('barcode_tree.urls'))
]

