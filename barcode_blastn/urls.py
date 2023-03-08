from django.urls import path, re_path
from barcode_blastn import views
from drf_yasg.views import get_schema_view
from drf_yasg import openapi
from rest_framework.permissions import AllowAny

schema_view = get_schema_view(
    openapi.Info(
        title="Barcode Identifier API",
        default_version="v0.0.1",
        description="This web-based API compares DNA sequence input against curated sequence libraries using a streamlined multi-step workflow with nucleotide BLAST (BLASTN), multiple sequence alignment and tree construction.",
        # TODO: Add terms of service?
        # terms_of_service="https://www.google.com/policies/terms/",
        # TODO: Contact required?
        # contact=openapi.Contact(email="contact@email.com"),
        # TODO: Is licence required?
        # license=openapi.License(name='MIT License")
    ),
    public=True,
    permission_classes=[AllowAny],
)

urlpatterns = [
    re_path(r'^swagger(?P<format>\.json|\.yaml)$', schema_view.without_ui(cache_timeout=0), name='schema-json'),
    re_path(r'^swagger/$', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
    re_path(r'^redoc/$', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc'),
    path('blastdbs/', views.BlastDbList.as_view()),
    path('blastdbs/<uuid:pk>/', views.BlastDbDetail.as_view()),
    path('blastdbs/<uuid:pk>/add/', views.NuccoreSequenceAdd.as_view()),
    path('blastdbs/<uuid:pk>/bulk/', views.NuccoreSequenceBulkAdd.as_view()),
    path('blastdbs/<uuid:pk>/run/', views.BlastRunRun.as_view()),
    path('nuccores/', views.NuccoreSequenceList.as_view()),
    path('nuccores/<uuid:pk>/', views.NuccoreSequenceDetail.as_view()),
    path('runs/', views.BlastRunList.as_view()),
    path('runs/<uuid:pk>/', views.BlastRunDetail.as_view()),
    path('runs/<uuid:pk>/status/', views.BlastRunStatus.as_view()),
    path('runs/<uuid:pk>/download/', views.BlastRunDetailDownload.as_view()),
    path('runs/<uuid:pk>/input-download/', views.BlastRunInputDownload.as_view())
]