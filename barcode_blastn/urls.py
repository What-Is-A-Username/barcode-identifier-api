from django.urls import path, re_path, include
from barcode_blastn import views
from drf_yasg.views import get_schema_view
from knox import views as knox_views
from drf_yasg import openapi
from rest_framework.permissions import AllowAny

from django.conf.urls.static import static
from django.conf import settings

schema_view = get_schema_view(
    openapi.Info(
        title="Barrel API",
        default_version="v0.0.2",
        description="This web-based API compares DNA sequence input against curated sequence libraries using a streamlined multi-step workflow with nucleotide BLAST (BLASTN), multiple sequence alignment and tree construction.",
        # TODO: Add terms of service?
        # terms_of_service="https://www.google.com/policies/terms/",
        contact=openapi.Contact(
            name="Lab Inquiries",
            url="http://www.utsc.utoronto.ca/~lovejoy/",
            email="nathan.lovejoy@utoronto.ca"),
        # TODO: Is licence required?
        # license=openapi.License(name='MIT License")
    ),
    public=True,
    permission_classes=[AllowAny],
)

urlpatterns = [
        # Swagger Docs 
        re_path(r'^swagger(?P<format>\.json|\.yaml)$', schema_view.without_ui(cache_timeout=0), name='schema-json'),
        re_path(r'^swagger/$', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
        re_path(r'^redoc/$', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc'),

        # Django-Rest-Knox Authentication
        path(r'login/', views.LoginView.as_view(), name='knox_login'),
        path(r'logout/', views.LogoutView.as_view(), name='knox_logout'),
        path(r'logoutall/', views.LogoutAllView.as_view(), name='knox_logoutall'),

        # Token Authentication Testing 
        path('users/', views.UserDetailView.as_view()),

        # Reference Libraries
        path('libraries/', views.LibrariesList.as_view()),
        path('libraries/<uuid:library>/versions', views.LibraryBlastDbList.as_view()),
        path('libraries/<uuid:library>', views.LibraryDetailView.as_view()),

        # Library Versions
        path('blastdbs/<uuid:pk>', views.BlastDbDetail.as_view()),
        path('blastdbs/<uuid:pk>/history', views.BlastDbHistory.as_view()),
        path('blastdbs/<uuid:pk>/sequences', views.BlastDbSequenceList.as_view()),
        path('blastdbs/<uuid:pk>/run', views.BlastRunRun.as_view()),
        path('blastdbs/<uuid:pk>/export', views.BlastDbExport.as_view()),
        path('blastdbs/<uuid:pk>/summary', views.BlastDbSummary.as_view()),

        # nuccore sequence
        path('nuccores/<uuid:pk>', views.NuccoreSequenceDetail.as_view()),

        # runs
        path('runs/', views.BlastRunList.as_view()),
        path('runs/<uuid:pk>', views.BlastRunDetail.as_view()),
        path('runs/<uuid:pk>/queries', views.BlastRunQueryList.as_view()),
        path('runs/<uuid:pk>/queries/<str:query>', views.BlastQuerySequenceDetail.as_view()),
        path('runs/<uuid:pk>/queries/<str:query>/hits', views.BlastQueryHitList.as_view()),
        path('runs/<uuid:pk>/status', views.BlastRunStatus.as_view()),
        path('runs/<uuid:pk>/download', views.BlastRunDetailDownload.as_view()),
        path('runs/<uuid:pk>/download/input', views.BlastRunInputDownload.as_view()),
        path('runs/<uuid:pk>/download/taxonomy', views.BlastRunTaxonomyDownload.as_view()),
    ]

# give access to media files if in debug mode (i.e. when nginx not serving them)
if settings.DEBUG:
    urlpatterns += static(
        settings.MEDIA_URL,
        document_root=settings.MEDIA_ROOT
    )