# from barcode_blastn.models import BlastDb, NuccoreSequence
# import uuid
# from barcode_blastn.serializers import NuccoreSequenceSerializer, NuccoreSequenceSaveSerializer

# db = BlastDb.objects.all()[0]
# # ns = NuccoreSequence(accession_number = '123', owner_database = db)
# # nss = NuccoreSequence(accession_number = '345', owner_database_id = db.id)

# data = {'accession_number': '123', 'owner_database': { 'id': 'd8d1f06f-4bf3-4f8e-9182-42c16bdaecdc' }}

# data = {'accession_number': '123', 'owner_database': 'd8d1f06f-4bf3-4f8e-9182-42c16bdaecdc'}

# serializer = NuccoreSequenceSaveSerializer(data = data)
# if serializer.is_valid():
#     serializer.save(owner_database = db)

# ref = {'accession_number': '123', 'owner_database': db }
# serializer = NuccoreSequenceSaveSerializer(db, ref = data)
# if serializer.is_valid():
#     serializer.save()
