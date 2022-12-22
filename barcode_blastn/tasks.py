import os
from shlex import quote
from datetime import datetime
import subprocess
from barcode_blastn.helper.parse_results import parse_results

from barcode_blastn.models import BlastRun, Hit, NuccoreSequence 

from celery import shared_task

@shared_task(time_limit=30)  # hard time limit of 30 seconds 
def run_blast_command(ncbi_blast_version, fishdb_id, run_id):   
    print('Beginning queued BLAST search ...')

    project_root = os.path.abspath(os.path.dirname('./'))
    print(f'Working directory is {project_root}')
    print(f'Using blast version {ncbi_blast_version}')
    print(f'Using database id {fishdb_id} and run id {run_id}')
    
    run_details = BlastRun.objects.get(id = run_id)

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
