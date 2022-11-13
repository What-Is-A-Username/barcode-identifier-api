from django.urls import path, include
from django.urls.conf import re_path
from barcode_blastn import views


urlpatterns = [
    path('blastdbs/', views.BlastDbList.as_view()),
    path('blastdbs/<uuid:pk>/', views.BlastDbDetail.as_view()),
    path('blastdbs/<uuid:pk>/add/', views.NuccoreSequenceAdd.as_view()),
    path('blastdbs/<uuid:pk>/run/', views.BlastRunRun.as_view()),
    path('nuccores/', views.NuccoreSequenceList.as_view()),
    path('nuccores/<uuid:pk>/', views.NuccoreSequenceDetail.as_view()),
    path('runs/', views.BlastRunList.as_view()),
    path('runs/<uuid:pk>/', views.BlastRunDetail.as_view()),
    path('django-rq/', include('django_rq.urls'))
]