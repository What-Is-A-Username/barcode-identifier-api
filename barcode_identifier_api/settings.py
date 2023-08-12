"""
Django settings for barcode_identifier_api project.

Generated by 'django-admin startproject' using Django 4.1.1.

For more information on this file, see
https://docs.djangoproject.com/en/4.1/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/4.1/ref/settings/
"""
import os
from pathlib import Path
from kombu.utils.url import safequote

# Fall back for local development environments where env variables not set from docker compose

if os.environ.get("DB_SECRET_KEY", "defaultkey") == "defaultkey":
    from dotenv import load_dotenv
    print("Loaded local environment variables from .env")
    load_dotenv('.env')

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/4.1/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = os.environ.get('SECRET_KEY')

# SECURITY WARNING: don't run with debug turned on in production!
# DEBUG = bool(int(os.environ.get('DEBUG', 0)))
DEBUG = 1

# TODO: Change when deploying to production
ALLOWED_HOSTS = [
'*'
]

ALLOWED_HOSTS.extend(filter(None,os.environ.get('ALLOWED_HOSTS', '').split(',')))

# Application definition
INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'rest_framework',
    'knox',
    'simple_history',
    'corsheaders',
    'drf_yasg',
    # Django extensions can be installed and used to generate 
    # graph models: https://django-extensions.readthedocs.io/en/latest/graph_models.html
    # 'django_extensions', 
    'barcode_blastn.apps.BarcodeBlastnConfig'
]

AUTHENTICATION_BACKENDS = (
    ('django.contrib.auth.backends.ModelBackend'),
)

REST_FRAMEWORK = {
    'DEFAULT_AUTHENTICATION_CLASSES': (
        'knox.auth.TokenAuthentication',
        'rest_framework.authentication.SessionAuthentication'
    ),
    'DEFAULT_RENDERER_CLASSES': (
        'rest_framework.renderers.JSONRenderer',
        # Only enable browsable API renderer in development,
        # to prevent scraping by bots and search engines.
        # 'rest_framework.renderers.BrowsableAPIRenderer'
    )
}

REST_KNOX = {
    'AUTH_HEADER_PREFIX': 'Bearer'
}

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'simple_history.middleware.HistoryRequestMiddleware' # Used to track which user made changes to libraries and databases
]

ROOT_URLCONF = 'barcode_identifier_api.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'barcode_identifier_api.wsgi.application'


# Database
# https://docs.djangoproject.com/en/4.1/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': os.environ.get("DB_NAME"),
        'USER': os.environ.get("DB_USER"),
        'PASSWORD': os.environ.get("DB_PASSWORD"),
        'HOST': os.environ.get("DB_HOST"), # set to name of container if using docker, or "localhost" otherwise 
        'PORT': int(os.environ.get('DB_PORT', 5432)),
    }
}

# TODO: Refine CORS_HEADERS settings once the domain is known.
# Docs: https://github.com/adamchainz/django-cors-headers#configuration
CORS_ALLOWED_ORIGINS = [
    # Allow local React app
    'http://localhost:3000',
    'https://localhost:3000',
    # Allow production website
    'http://localhost:8000',
    'https://localhost:8000',
]

CORS_ALLOW_CREDENTIALS = True


# Password validation
# https://docs.djangoproject.com/en/4.1/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/4.1/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/4.1/howto/static-files/

# Old URL
# STATIC_URL = 'static/'

STATIC_URL = '/static/static/'
MEDIA_URL = '/static/media/'

MEDIA_ROOT = '/vol/static/media'
STATIC_ROOT = '/vol/static/static'

# STATIC_ROOT = os.path.join(BASE_DIR, "static/")

# Default primary key field type
# https://docs.djangoproject.com/en/4.1/ref/settings/#default-auto-field

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

# CELERY
AWS_ACCESS_KEY_ID = safequote(os.environ.get('AWS_ACCESS_KEY_ID', ''))

AWS_SECRET_ACCESS_KEY = safequote(os.environ.get('AWS_SECRET_ACCESS_KEY', ''))

CELERY_BROKER_URL = os.environ.get('CELERY_BROKER_URL')

# CELERY_BROKER_TRANSPORT_OPTIONS = {
#     'region': 'us-east-2',
#     'visibility_timeout': 1800,
#     'polling_interval': 1,
#     'wait_time_seconds': 10,
#     'predefined_queues': {
#         'BarcodeQueue.fifo': {
#             'url': os.environ.get('CELERY_BROKER_TRANSPORT_BARCODE_URL'),
#             'access_key_id': os.environ.get('CELERY_BROKER_TRANSPORT_ACCESS_KEY_ID'),
#             'secret_access_key': os.environ.get('CELERY_BROKER_TRANSPORT_SECRET_ACCESS_KEY'),
#             'backoff_policy': {1: 10, 2: 15, 3: 20},
#             'backoff_tasks': ['barcode_blastn.tasks.run_blast_command'],
#         }
#     }
# }

CELERY_USER = "appuser" 
CELERY_GROUP = "appgroup"

# CELERY_ROUTES = {
#     'barcode_blastn.tasks.run_blast_command': 'BarcodeQueue.fifo',
# }

# CELERY_TASK_DEFAULT_QUEUE = 'BarcodeQueue.fifo'

# SWAGGER API DOC SETTINGS
# https://drf-yasg.readthedocs.io/en/stable/settings.html
SWAGGER_SETTINGS = {
    'USE_SESSION_AUTH': True, # hide authentication
    'DEFAULT_MODEL_DEPTH': -1, # hide models widget
    'DEFAULT_MODEL_RENDERING': 'example', # show examples by default
    'SECURITY_DEFINITIONS': {
        'Bearer': {
            'type': 'apiKey',
            'name': 'Authorization',
            'in': 'header'
        }
    },
}


# FILE UPLOAD
FILE_UPLOAD_MAX_MEMORY_SIZE = 2621440 # = 2.5MB
DATA_UPLOAD_MAX_MEMORY_SIZE = 2621440 # = 2.5MB
DATA_UPLOAD_MAX_NUMBER_FILES = 1 # Only allow one file upload via POST in a multipart/form-data request
FILE_UPLOAD_MAX_NUMBER_FILES = 0 # Ensure that we disable file uploads to the file system

# DJANGO-SIMPLE-HISTORY
SIMPLE_HISTORY_REVERT_DISABLED = True # Prevent users from reverting to earlier versions

GRAPH_MODELS = {
    'all_applications': True,
    'group_models': True
}