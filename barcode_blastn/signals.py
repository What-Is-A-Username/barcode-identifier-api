from django.dispatch import receiver
from simple_history.signals import pre_create_historical_record

@receiver(pre_create_historical_record)
def add_history_blastdb_sequences(sender, instance, **kwargs):
    '''Signal handler/callback for also logging the original GenBank search terms
    and accessions'''

    # Import object storing HistoricalData created by django-simple-history
    from barcode_blastn.models import HistoricalBlastDb, HistoricalLibrary

    if sender == HistoricalBlastDb:
        history_instance = kwargs['history_instance']
        # Move data from instance to the record
        for attr in ['added', 'deleted', 'search_terms']:
            setattr(history_instance, attr, getattr(instance, attr, ''))
    elif sender == HistoricalLibrary:
        history_instance = kwargs['history_instance']
    else:
        raise NotImplementedError(f'Signal for {str(sender)}')