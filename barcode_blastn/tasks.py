import os
from ratelimit import limits
from shlex import quote
from datetime import datetime, timezone
import subprocess
from typing import Dict, List
from barcode_blastn.helper.parse_gb import retrieve_gb
from barcode_blastn.helper.parse_results import parse_results

from barcode_blastn.models import BlastDb, BlastRun, Hit, NuccoreSequence 
from django.db.models import QuerySet
from celery import shared_task

# update the blast run with the errors and print to console
def raise_error(run: BlastRun, err: str) -> None:
    print(err)
    run.errors = run.errors + "\n" + err
    run.job_status = BlastRun.JobStatus.ERRORED
    run.save()
    return

@shared_task(time_limit=30)  # hard time limit of 30 seconds 
def run_blast_command(ncbi_blast_version: str, fishdb_id: str, run_id: str) -> None:   
    '''
        Performs a BLASTn search using the specified ncbi_blast_version, database with id fishdb_id, for the run of id run_id
    '''
    print('Beginning queued BLAST search ...')

    project_root = os.path.abspath(os.path.dirname('./'))
    print(f'Working directory is {project_root}')
    print(f'Using blast version {ncbi_blast_version}')
    print(f'Using database id {fishdb_id} and run id {run_id}')
    
    # TODO: throw error if blastrun DNE
    try:
        run_details : BlastRun = BlastRun.objects.get(id = run_id)
    except BlastRun.DoesNotExist:
        print(f"Critical error: BlastRun of id {run_id} does not exist.")
        return

    run_details.job_status = BlastRun.JobStatus.STARTED
    run_details.job_start_time = datetime.now()
    run_details.save()

    blast_root =  f'{project_root}/{ncbi_blast_version}/bin'
    command_app = blast_root + '/blastn'
    # check if blast executable are present
    if not os.path.exists(command_app) or not os.path.isfile(command_app):
        raise_error(run_details, f"Critical error: Failed to find BLASTn executable at {command_app}")
        return

    db_path = project_root + '/fishdb/' + fishdb_id + '/database' 
    # check if the database file exists
    if not os.path.exists(db_path + '.fasta') or not os.path.isfile(db_path + '.fasta'):
        raise_error(run_details, f"Critical error: Failed to find BLAST database at {db_path}")
        return

    run_folder = f'{project_root}/runs/{run_id}'
    results_file = f'{run_folder}/results.txt'
    errors_file = f'{run_folder}/errors.txt'
    query_file = f'{run_folder}/query.fasta'

    # check if run status folder exists
    if not os.path.exists(run_folder) or not os.path.isdir(run_folder):
        raise_error(run_details, f"Critical error: Failed to find run folder at {run_folder}")
        return
    # check if query file exists
    if not os.path.exists(query_file) or not os.path.isfile(query_file):
        raise_error(run_details, f"Critical error: Failed to find query sequence file at {query_file}")
        return

    outfmt = 7

    # build the command to run in the shell
    blast_command = f'{quote(command_app)} -db {quote(db_path)} -outfmt {outfmt} -query {quote(query_file)}'

    process = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    print('BLAST search completed.')

    out, err  = 'Unable to find out', 'Unable to find err'
    try:
        if process.stdout is not None:
            out = process.stdout.read().decode()
        if process.stderr is not None:
            err = process.stderr.read().decode()
    except BaseException:
        raise_error(run_details, f"Errored while retrieving shell output.")

    print('Writing results to file ...')

    try:
        with open(results_file, "w") as results_out:
            results_out.write(out + '\n')
        if err:
            with open(errors_file, "w") as errors_out:
                errors_out.write(err)
    except BaseException:
        raise_error(run_details, f"Errored while handling results output.")
        return
       
    # bulk create hits
    try:
        parsed_data = parse_results(out)
        accessions = [hit['subject_accession_version'] for hit in parsed_data]
        accession_entries = NuccoreSequence.objects.filter(accession_number__in=accessions)
        all_hits = [Hit(**hit, owner_run = run_details, db_entry = accession_entries.get(accession_number=hit['subject_accession_version'])) for hit in parsed_data]
        Hit.objects.bulk_create(all_hits)
    except BaseException:
        raise_error(run_details, f"Errored while bulk creating hit results.")
        return

    # update run information
    try:
        run_details.errors = err
        run_details.job_status = BlastRun.JobStatus.FINISHED
        run_details.job_end_time = datetime.now()
        run_details.save()
    except BaseException:
        raise_error(run_details, f"Errored while updating run status.")
        return

    print('Queued BLAST search completed.')
    return

PERIOD = 3500 # Time between calls, in number of seconds
@limits(calls = 1, period = PERIOD)
@shared_task(time_limit=30)  # hard time limit of 30 seconds 
def update_database() -> None:
    print("update_database ran at " + datetime.now().strftime("%H:%M:%S.%f"))
    all_dbs : QuerySet = BlastDb.objects.all() 
    db : BlastDb
    for db in all_dbs:
        all_seqs: QuerySet = NuccoreSequence.objects.filter(owner_database=db.id)
        all_numbers: List[str] = [seq.accession_number for seq in all_seqs]
        print("Updating database " + str(db.id) + " ...")

        new_data = retrieve_gb(all_numbers)
        fields_to_update = list(new_data[0].keys())
        fields_to_update.append('created')
        seq : NuccoreSequence

        # make a dictionary mapping accession -> uid
        an_to_uid : Dict = {}
        for seq in all_seqs:
            an_to_uid[seq.accession_number] = str(seq.id)

        new_created_time = datetime.now(timezone.utc)

        # using the fetched data as a base, add updated values for 'created' and ensure that the id is present
        def make_updated_dict(old_dict: Dict):
            old_dict['created'] = new_created_time
            old_dict['id'] = an_to_uid[old_dict['accession_number']]
            return old_dict

        updated_data = [make_updated_dict(entry) for entry in new_data if entry['accession_number'] in an_to_uid]
        
        print(f"Submitting updates to database ...")
        
        NuccoreSequence.objects.bulk_update([NuccoreSequence(**data) for data in updated_data], fields=fields_to_update, batch_size=100)

        print(f"Finished updating database {str(db.id)}")
    print("Finished updating all databases.")
        
