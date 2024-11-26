from typing import Any, List
from django.dispatch import receiver
from django.db.models.signals import post_save
from simple_history.signals import pre_create_historical_record

from barcode_blastn.models import Annotation, CustomSequence, NuccoreSequence

@receiver(pre_create_historical_record)
def add_history_blastdb_sequences(sender, instance, **kwargs):
    '''Signal handler/callback for also logging the original GenBank search terms
    and accessions'''

    # Import object storing HistoricalData created by django-simple-history
    from barcode_blastn.models import HistoricalBlastDb, HistoricalLibrary

    if sender == HistoricalBlastDb:
        history_instance = kwargs['history_instance']
        # Move data from instance to the record
        for attr in ['added', 'deleted', 'search_terms', 'filter_options']:
            setattr(history_instance, attr, getattr(instance, attr, ''))
    elif sender == HistoricalLibrary:
        history_instance = kwargs['history_instance']
    else:
        raise NotImplementedError(f'Signal for {str(sender)}')

@receiver(post_save, sender=NuccoreSequence)
def save_annotations(sender, instance, created, **kwargs):
    '''
    Save annotations based on `create_annotations` on the instance.
    
    If no poster is specified for the annotations of an instance, then the value
    specified by `annotation_user` is used as the default.
    '''
    annotations = getattr(instance, 'create_annotations', [])
    default_user = getattr(instance, 'annotation_user')
    new_annotations: List[Annotation] = []
    annotation: dict[str, Any]
    for annotation in annotations:
        annotation.setdefault('poster', default_user)
        new_annotations.append(Annotation(sequence=instance, **annotation))
    Annotation.objects.bulk_create(new_annotations)

@receiver(post_save, sender=CustomSequence)
def save_annotations(sender, instance, created, **kwargs):
    '''
    Save annotations based on `create_annotations` on the instance.
    
    If no poster is specified for the annotations of an instance, then the value
    specified by `annotation_user` is used as the default.
    '''
    annotations = getattr(instance, 'create_annotations', [])
    default_user = getattr(instance, 'annotation_user')
    new_annotations: List[Annotation] = []
    annotation: dict[str, Any]
    for annotation in annotations:
        annotation.setdefault('poster', default_user)
        new_annotations.append(Annotation(sequence=instance, **annotation))
    Annotation.objects.bulk_create(new_annotations)