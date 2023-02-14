from ctypes import alignment
from django.db import models
from barcode_blastn.models import BlastRun
from django.utils.translation import gettext_lazy as _

class ResultTree(models.Model):
    '''
        Stores a tree generated from the sequence input
    '''
    # The run from which the results were gathered
    owner_run = models.ForeignKey(BlastRun, related_name='tree', on_delete=models.CASCADE)
    
    # Job ID of the asynchronous request to Multiple Sequence Alignment with ClustalOmega
    alignment_job_id = models.CharField(max_length=100, blank=True, default='')

    class TreeStatus(models.TextChoices):
        # Waiting for web tool to return an aligned sequence
        ALIGNING = 'ALN', _('ALIGNING')
        # Waiting for alignment output to be saved
        PROCESSING = 'PRO', _('PROCESSING')
        # Waiting for Simple Phylogeny to construct tree (Obselete)
        CONSTRUCTING = 'CON', _('CONSTRUCTING')
        # Finalizing and cleaning up tree data (Obselete)
        CLEANING = 'CLE', _('CLEANING')
        # Finished constructing tree and saving it locally
        FINISHED = 'FIN', _('FINISHED')
        # Unexpected error 
        ERRORED = 'ERR', _('ERRORED')

    internal_status = models.CharField(max_length=3,choices=TreeStatus.choices, default=TreeStatus.ALIGNING)
        
        
        
