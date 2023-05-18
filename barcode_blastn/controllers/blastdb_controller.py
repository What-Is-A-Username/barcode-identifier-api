from datetime import datetime, timezone
import os
import shutil
from typing import Any, Dict, List, Optional, Tuple
from barcode_blastn.file_paths import get_data_fishdb_path
from barcode_blastn.helper.parse_gb import AccessionLimitExceeded, GenBankConnectionError, InsufficientAccessionData, retrieve_gb
from barcode_blastn.models import BlastDb, Library, NuccoreSequence
from barcode_blastn.serializers import NuccoreSequenceSerializer 
from django.db.models import QuerySet

class SequenceUpdateSummary:
    no_change = []
    accession_version_changed = []
    metadata_changed = []
    deleted = []
    added = []

def delete_blastdb(blast_db: Optional[BlastDb]) -> None:
    '''
    Delete a BLAST database. Also remove its files from the filesystem.

    Raises ValueError if value is None

    Raises OSError if error encountered while deleting files.
    '''
    if blast_db is None:
        raise ValueError('Library cannnot be None')

    # delete the database from the file system
    database_id: str = str(blast_db.id)
    try:
        local_db_folder = get_data_fishdb_path(database_id)
        if len(database_id) > 0 and os.path.exists(local_db_folder):
            shutil.rmtree(local_db_folder, ignore_errors=True)
    except BaseException as exc:
        raise OSError(exc)
    else:
        blast_db.delete()

def delete_library(library: Optional[Library]):
    '''
    Shared logic for deleting a reference library

    Raises ValueError if value is None

    Raises OSError if error encountered while deleting files.
    '''
    if library is None:
        raise ValueError('Library cannnot be None')

    dbs = BlastDb.objects.filter(library=library)
    db: BlastDb
    for db in dbs:
        delete_blastdb(db)
    library.delete()     

def save_blastdb(obj: BlastDb, perform_lock: bool = False) -> BlastDb:
    '''
    Save the BlastDb. The BlastDb should have accession numbers **added and saved**, reference library indicated.

    The save process performed includes:
    - locking the database, if perform_lock is True
    - assign a version number if it was previously locked, based on the latest library version
    - calling `.save()` on the instance

    Raises:
    
    '''
    if perform_lock:
        lastPublished: Optional[BlastDb] = BlastDb.objects.latest(obj.library)
        version_nums = (1,1,1) # if this is the first published version of the library, assign version 1.1.1 
        if lastPublished is not None:
            sequence_summary: SequenceUpdateSummary = calculate_update_summary(last=lastPublished, current=obj)
            if len(sequence_summary.deleted) > 0 or len(sequence_summary.added) > 0 or len(sequence_summary.accession_version_changed) > 0:
                version_nums = (lastPublished.genbank_version + 1, 1, 1)
            elif len(sequence_summary.metadata_changed) > 0:
                version_nums = (lastPublished.genbank_version, lastPublished.major_version + 1, 1)
            else:
                version_nums = (lastPublished.genbank_version, lastPublished.major_version, lastPublished.minor_version + 1)
        
        obj.genbank_version = version_nums[0]
        obj.major_version = version_nums[1]
        obj.minor_version = version_nums[2]
        obj.locked = True

    obj.save()

    return obj

def create_blastdb(additional_accessions: List[str], base: Optional[BlastDb] = None, **kwargs) -> BlastDb:
    '''
    Create an entirely new blastdb with the accession numbers in additional_accessions. Also add the data from base.
    Data for every accession, from both additional_accessions and base, will be refetched from GenBank.
    The new blastdb and accessions will be saved. The saved instance will be returned.

    Returns the new database.

    Raises various exceptions (ValueError, GenBankConnectionError, AccessionLimitExceeded, InsufficientAccessionData) if errors encountered while fetching data from GenBank.
    '''
    genbank_version = kwargs.pop('genbank_version', 0)
    major_version = kwargs.pop('major_version', 0)
    minor_version = kwargs.pop('minor_version', 0)
    created = kwargs.pop('created', datetime.now())
    locked = kwargs.pop('locked', False) # ignore value of locked for now
    new_database: BlastDb = BlastDb(genbank_version=genbank_version, major_version=major_version, minor_version=minor_version, created=created, locked=False, **kwargs)
    new_database.save()
    
    accessions_to_add = additional_accessions
    if base is not None:
        seqs = NuccoreSequence.objects.filter(owner_database=base)
        accessions_to_add.extend([s.accession_number for s in seqs])
    
    if len(accessions_to_add) > 0:
        add_sequences_to_database(new_database, desired_numbers=accessions_to_add)

    if locked:
        new_database = save_blastdb(new_database, perform_lock=locked)
    return new_database

def calculate_update_summary(last: BlastDb, current: BlastDb):
    fields_to_check = ['definition', 'dna_sequence', 'organism', 'organelle', 'isolate', 'country', 'specimen_voucher', 'type_material', 'lat_lon']

    last_sequences: QuerySet[NuccoreSequence] = NuccoreSequence.objects.filter(owner_database=last)
    current_sequences: QuerySet[NuccoreSequence] = NuccoreSequence.objects.filter(owner_database=current)

    summary: SequenceUpdateSummary = SequenceUpdateSummary()

    an_to_last: Dict[str, NuccoreSequence] = {}
    an_to_current: Dict[str, NuccoreSequence] = {}
    seq: NuccoreSequence
    for seq in current_sequences:
        an_to_current[seq.accession_number] = seq
    for seq in last_sequences:
        an_to_last[seq.accession_number] = seq
        if seq.accession_number not in an_to_current:
            summary.deleted.append(seq.accession_number)
    
    for current_seq in current_sequences:
        if current_seq.accession_number not in an_to_last:
            summary.added.append(current_seq.accession_number)
        else:
            instance = an_to_last[current_seq.accession_number]
            if current_seq.version != instance.version or current_seq.dna_sequence != instance.dna_sequence:
                summary.accession_version_changed.append(current_seq.accession_number)
            else:
                has_any_field_changed: bool = False
                for field in fields_to_check:
                    if not has_any_field_changed and getattr(current_seq, field) != getattr(instance, field):
                        has_any_field_changed = True 
                
                if has_any_field_changed:
                    summary.metadata_changed.append(instance.accession_number)
                else:
                    summary.no_change.append(instance.accession_number)
    return summary

def update_from_genbank(sequence_instances: QuerySet[NuccoreSequence]) -> SequenceUpdateSummary:
    '''
    Update a query set of existing and saved NuccoreSequence objects by retrieving new GenBank data and saving the data back to the database.
    If new GenBank data cannot be found for an existing instance accession number, then the instance is deleted.

    Returns a summary of the sequences updated.
    '''
    all_numbers: List[str] = [seq.accession_number for seq in sequence_instances]
    new_data = retrieve_gb(all_numbers)
    
    seq : NuccoreSequence
    # make a dictionary mapping genbank accession -> genbank entry
    an_to_gb : Dict[str, Dict[str, str]] = {}
    update_summary: SequenceUpdateSummary = SequenceUpdateSummary()
    for genbank_data in new_data:
        an_to_gb[genbank_data['accession_number']] = genbank_data

    new_created_time = datetime.now(timezone.utc)
    genbank_dict: Dict[str, Any]
    instance: NuccoreSequence

    # using the fetched data as a base, add updated values for 'created' and ensure that the id is present
    fields_to_update = list(new_data[0].keys())
    fields_to_update.append('created')
    to_update: List[NuccoreSequence] = []
    to_delete: List[NuccoreSequence] = []
    for instance in sequence_instances:
        if instance.accession_number not in an_to_gb: # genbank data missing accession
            update_summary.deleted.append(instance.accession_number)
            to_delete.append(instance)
        else:
            genbank_dict = an_to_gb[instance.accession_number]
            instance.created = new_created_time
            has_any_field_changed: bool = False
            has_version_changed: bool = instance.version != genbank_dict['version']

            for field in fields_to_update:
                if not has_any_field_changed and genbank_dict[field] != getattr(instance, field):
                    has_any_field_changed = True 
                setattr(instance, field, genbank_dict[field])
            
            if has_version_changed:
                update_summary.accession_version_changed.append(instance.accession_number)
            elif has_any_field_changed:
                update_summary.metadata_changed.append(instance.accession_number)
            else:
                update_summary.no_change.append(instance.accession_number)
            to_update.append(instance)

    NuccoreSequence.objects.bulk_update(to_update, fields=fields_to_update, batch_size=100) 
    for delete_seq in to_delete:
        delete_seq.delete() 
    return update_summary     

class AccessionsAlreadyExist(BaseException): 
    '''
    A set of accession numbers already exist in a database.
    '''
    accession_numbers: List[str]
    def __init__(self, accession_numbers: List[str]) -> None:
        self.accession_numbers = accession_numbers

class AccessionsNotFound(BaseException):
    '''
    A set of accession numbers could not be located in the database.
    '''
    accession_numbers: List[str]
    def __init__(self, accession_numbers: List[str]) -> None:
        self.accession_numbers = accession_numbers

class DatabaseLocked(BaseException): ...

def delete_sequences_in_database(database: BlastDb, desired_nums: List[str]) -> int:
    '''
    Remove all sequences.

    Returns the number of sequence objects deleted.

    Raises
        DatabaseLocked: If database is locked for editing.
    '''
    if database.locked:
        raise DatabaseLocked()
    desired_nums = list(set(desired_nums))
    to_delete: QuerySet[NuccoreSequence] = NuccoreSequence.objects.filter(owner_database=database, accession_number__in=desired_nums)
    result = to_delete.delete()

    return result[1]['barcode_blastn.NuccoreSequence']

def update_sequences_in_database(database: BlastDb, desired_numbers: List[str]) -> List[NuccoreSequence]:
    ''''
    Bulk update the given accession numbers in the database and return the resulting list of NuccoreSequences. If a number doesn't exists, then an AccessionsNotFound Error is raised.

    Returns a tuple of two lists. The first list contains the saved NuccoreSequence objects, as returned by `.bulk_create()`. The second list contains the updated objects, as returned by `.bulk_update()`.

    Raises:
        DatabaseLocked: If database is locked for editing.
        ValueError: If no accession numbers are specified
        AccessionsNotFound: If a given accession number is not present in the database.
        AccessionLimitExceeded: If the number of accessions to add exceeds the maximum allowed
        GenbankConnectionError: Could not connect to GenBank or the request sent was bad
        InsufficientAccessionData: If all accession numbers could not be identified by GenBank 
    '''
    if database.locked:
        raise DatabaseLocked()
    if len(desired_numbers) == 0:
        raise ValueError('No accessions given')
    desired_numbers = list(set(desired_numbers))
    # Check what accession numbers are duplicate (i.e. already existing in the datab)
    existing: QuerySet[NuccoreSequence] = NuccoreSequence.objects.filter(owner_database=database, accession_number__in=desired_numbers) # type: ignore        
    # raise error if not all numbers exist
    if existing.count() < len(desired_numbers):
        e: NuccoreSequence
        existing_nums: List[str] = [e.accession_number for e in existing]
        raise AccessionsNotFound([d for d in desired_numbers if d not in existing_nums])

    # retrieve data
    genbank_data = retrieve_gb(accession_numbers=desired_numbers)
    keys = [d['accession_number'] for d in genbank_data]
    # map accession numbers -> retrieved data dictionary
    acc_to_data: Dict[str, Dict[str, str]] = dict(zip(keys, genbank_data))

    # save the new sequences using a bulk operation
    to_update: List[NuccoreSequence] = []
    record: NuccoreSequence
    for record in existing:
        for key, value in acc_to_data[record.accession_number].items():
            setattr(record, key, str(value))
        record.updated = datetime.now()
        to_update.append(record)

    fields_to_update = list(genbank_data[0].keys())
    fields_to_update.append('updated')
    NuccoreSequence.objects.bulk_update(existing, fields=fields_to_update)
    return to_update


def add_sequences_to_database(database: BlastDb, desired_numbers: List[str]) -> List[NuccoreSequence]:
    ''''
    Add a list of accession numbers to a database by bulk creating and return the resulting list of NuccoreSequences. If a number already exists, then an AccessionsAlreadyExist Error is raised.

    Returns a tuple of two lists. The first list contains the saved NuccoreSequence objects, as returned by `.bulk_create()`. The second list contains the updated objects, as returned by `.bulk_update()`.

    Raises:
        DatabaseLocked: If database is locked for editing.
        ValueError: If no accession numbers are specified
        AccessionsAlreadyExist: If an accession already exists in the database
        AccessionLimitExceeded: If the number of accessions to add exceeds the maximum allowed
        GenbankConnectionError: Could not connect to GenBank or the request sent was bad
        InsufficientAccessionData: If all accession numbers could not be identified by GenBank 
    '''
    if database.locked:
        raise DatabaseLocked()
    if len(desired_numbers) == 0:
        raise ValueError('No accessions given')
    desired_numbers = list(set(desired_numbers))
    # Check what accession numbers are duplicate (i.e. already existing in the datab)
    existing: QuerySet[NuccoreSequence] = NuccoreSequence.objects.filter(owner_database=database, accession_number__in=desired_numbers) # type: ignore        
    # raise error if any the numbers to be added already exist
    if existing.count() > 0:
        e: NuccoreSequence
        existing_nums: List[str] = [e.accession_number for e in existing]
        raise AccessionsAlreadyExist(existing_nums)

    # retrieve data
    genbank_data = retrieve_gb(accession_numbers=desired_numbers)

    # save the new sequences using a bulk operation
    to_create: List[NuccoreSequence] = []
    for data in genbank_data:
        an: str = data['accession_number']
        try:
            old_version: NuccoreSequence = existing.get(accession_number=an)
        except NuccoreSequence.DoesNotExist:
            to_create.append(NuccoreSequence(owner_database=database, **data))
        else:
            raise AccessionsAlreadyExist([an])

    created_sequences = []
    created_sequences = NuccoreSequence.objects.bulk_create(to_create)
    return created_sequences

def save_sequence(obj: NuccoreSequence, change: bool=True, commit: bool = False, raise_if_missing: bool = False):
    '''
    Save an accession by fetching data from GenBank again and populating the instance. Used for both saving and updating an accession.
    Obj should already have id, accession_number and owner_database values.
    If change is to False, raise AccessionsAlreadyExist if a query returns a entry in the same database with the same accession_number
    If commit is set to True, also save the instance. If set to False, you will need to call save manually after this function.

    Raises:
        ValueError does not have accession_number and owner_database values.

        AccessionsAlreadyExist: If change is False, and an existing accession with the same accession number already exists in the same database.

        AssertionError: If sequence claims to be editing an existing obj but that obj does not exist yet in the database

        GenBankConnectionError: Error connecting to GenBank.

        InsufficientAccessionData: If data from GenBank is insufficient for populating the fields.
    '''
    #TODO: Decide whether to keep 'change' parameter
    if not obj.accession_number or len(obj.accession_number) == 0:
        raise ValueError()
    if not obj.owner_database:
        raise ValueError()
    if obj.owner_database.locked:
        raise ValueError('Database locked')

    accession_number = obj.accession_number

    # check if there is a duplicate
    try:
        existing = NuccoreSequence.objects.filter(owner_database=obj.owner_database).exclude(id=obj.id).get(accession_number=accession_number)
    except NuccoreSequence.DoesNotExist:
        pass
    else:
        # raise error if duplicate accession already exists
        raise AccessionsAlreadyExist([e.accession_number for e in existing]) 

    # fetch GenBank data
    try:
        currentData = retrieve_gb(accession_numbers=[accession_number],
                                  raise_if_missing=raise_if_missing)[0]
    except GenBankConnectionError as exc:
        raise exc
    except InsufficientAccessionData as exc:
        raise exc

    # check that the GenBank data is valid
    serializer = NuccoreSequenceSerializer(obj, data=currentData, partial=True)
    if serializer.is_valid():
        if commit:
            # only save if commit specified
            return serializer.save()
        else:
            # only update instance key values and don't save
            for key, value in currentData.items():
                setattr(obj, key, str(value))
            return obj
    else:
        # raise an error if the genbank data format is not valid
        raise AssertionError(serializer.errors)

        
        
    
