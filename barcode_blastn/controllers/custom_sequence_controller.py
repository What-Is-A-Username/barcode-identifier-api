from io import StringIO
from typing import List
from django.core.files.uploadedfile import UploadedFile

from barcode_blastn.models import BlastDb, CustomSequence, NuccoreSequence
from Bio import SeqIO


def parse_file_upload_to_custom_sequence(upload: UploadedFile, database: BlastDb) -> List[CustomSequence]:
    """
    Parse the uploaded file into custom sequence objects without saving them.
    """
    if upload is None:
        return []
    upload.seek(0)
    seqs: List[CustomSequence] = []
    record : SeqIO.SeqRecord
    data = [line.decode('utf-8') for line in upload.readlines()]
    for record in SeqIO.parse(StringIO('\n'.join(data)), "fasta"):
        accession, definition = record.description.split(' ', 1)
        version = accession
        if '.' in accession:
            accession = accession.rsplit('.', 1)[0]
        seqs.append(
            CustomSequence(
                dna_sequence=str(record.seq),
                accession_number=accession,
                version=version,
                owner_database=database,
                definition=definition,
                data_source=NuccoreSequence.SequenceSource.IMPORT
            )
        )
    return seqs