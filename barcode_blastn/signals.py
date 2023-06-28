from django.dispatch import receiver
from simple_history.signals import pre_create_historical_record

@receiver(pre_create_historical_record)
def add_history_blastdb_sequences(sender, instance, **kwargs):
    '''Signal handler/callback for also logging the original GenBank search terms
    and accessions'''

    # Import object storing HistoricalData created by django-simple-history
    from barcode_blastn.models import HistoricalBlastDb

    if sender == HistoricalBlastDb:
        history_instance = kwargs['history_instance']
        # Move data from instance to the record
        history_instance.ids = getattr(instance, 'ids', '')
        history_instance.search_terms = getattr(instance, 'search_terms', '')
    else:
        raise NotImplementedError(f'Signal for {str(sender)}')