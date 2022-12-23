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

@shared_task(time_limit=30)  # hard time limit of 30 seconds 
def run_blast_command(ncbi_blast_version, fishdb_id, run_id):   
    print('Beginning queued BLAST search ...')

    project_root = os.path.abspath(os.path.dirname('./'))
    print(f'Working directory is {project_root}')
    print(f'Using blast version {ncbi_blast_version}')
    print(f'Using database id {fishdb_id} and run id {run_id}')
    
    run_details : BlastRun = BlastRun.objects.get(id = run_id)

    run_details.job_status = BlastRun.JobStatus.STARTED
    run_details.job_start_time = datetime.now()
    run_details.save()

    blast_root =  f'{project_root}/{ncbi_blast_version}/bin'
    db_path = project_root + '/fishdb/' + fishdb_id + '/database' 

    run_folder = f'{project_root}/runs/{run_id}'
    results_file = f'{run_folder}/results.txt'
    query_file = f'{run_folder}/query.fasta'

    if not os.path.exists(run_folder) or not os.path.exists(query_file):
        run_details.job_status = BlastRun.JobStatus.ERRORED

    command_app = blast_root + '/blastn'
    outfmt = 7
    blast_command = f'{quote(command_app)} -db {quote(db_path)} -outfmt {outfmt} -query {quote(query_file)}'

    process = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    print('BLAST search completed.')

    out, err = 'Unable to find out', 'Unable to find err'
    if process.stdout is not None:
        out = process.stdout.read().decode()
    if process.stderr is not None:
        err = process.stderr.read().decode()

    print('Writing results to file ...')

    try:
        with open(results_file, "w") as results_out:
            results_out.write(out + '\n')
            if err:
                results_out.write('\nErrors occurred when processing the query:\n')
                results_out.write(err)
    except BaseException:
        run_details.job_status = BlastRun.JobStatus.ERRORED
    else:
        results_out.close()

    # bulk create hits
    try:
        parsed_data = parse_results(out)
        accessions = [hit['subject_accession_version'] for hit in parsed_data]
        accession_entries = NuccoreSequence.objects.filter(accession_number__in=accessions)
        all_hits = [Hit(**hit, owner_run = run_details, db_entry = accession_entries.get(accession_number=hit['subject_accession_version'])) for hit in parsed_data]
        Hit.objects.bulk_create(all_hits)

        run_details.errors = err
        run_details.job_status = BlastRun.JobStatus.FINISHED
        run_details.job_end_time = datetime.now()
        run_details.save()
    except BaseException:
        run_details.job_status = BlastRun.JobStatus.ERRORED

    # TODO: Add a soft time limit to the function, and update the job status as ERRORED if reached. Docs: https://docs.celeryq.dev/en/stable/userguide/workers.html#time-limits

    print('Queued BLAST search completed.')

PERIOD = 3500 # Time between calls, in number of seconds
@limits(calls = 1, period = PERIOD)
@shared_task(time_limit=30)  # hard time limit of 30 seconds 
def update_database():
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
        
