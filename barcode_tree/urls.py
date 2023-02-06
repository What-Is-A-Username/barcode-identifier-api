from django.urls import path, include
from django.urls.conf import re_path
from barcode_tree import views
urlpatterns = [
    path('trees/<uuid:pk>/', views.ResultTreeDetail.as_view()),
]