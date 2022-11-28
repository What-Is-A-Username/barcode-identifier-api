from datetime import datetime
import subprocess
from time import sleep
from barcode_blastn.helper.parse_results import parse_results

from barcode_blastn.models import BlastRun, Hit, NuccoreSequence 

from celery import shared_task

@shared_task
def run_blast_command(blast_root, fishdb_path, query_file, run_details_id, results_path):   
    print('Beginning queued BLAST search ...')
    
    run_details = BlastRun.objects.get(id = run_details_id)

    run_details.job_status = BlastRun.JobStatus.STARTED
    run_details.job_start_time = datetime.now()
    run_details.save()

    command_app = blast_root + '/blastn'
    db_path = fishdb_path + '/database'
    outfmt = 7
    blast_command = f'{command_app} -db {db_path} -outfmt {outfmt} -query {query_file}'

    # TODO: avoid shell=True or add escaping https://docs.python.org/3/library/shlex.html

    process = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    sleep(10)

    print('BLAST search completed.')

    out, err = 'Unable to find out', 'Unable to find err'
    if process.stdout is not None:
        out = process.stdout.read().decode()
    if process.stderr is not None:
        err = process.stderr.read().decode()

    print('Writing results to file ...')

    with open(results_path + '/results.txt', "w") as results_file:
        results_file.write(out + '\n')
        if err:
            results_file.write('Errors occurred when processing the query:\n')
            results_file.write(err)
    results_file.close()


    # bulk create hits
    parsed_data = parse_results(out)
    accessions = [hit['subject_accession_version'] for hit in parsed_data]
    accession_entries = NuccoreSequence.objects.filter(accession_number__in=accessions)
    all_hits = [Hit(**hit, owner_run = run_details, db_entry = accession_entries.get(accession_number=hit['subject_accession_version'])) for hit in parsed_data]
    Hit.objects.bulk_create(all_hits)

    run_details.errors = err
    run_details.job_status = BlastRun.JobStatus.FINISHED
    run_details.job_end_time = datetime.now()
    run_details.save()

    print('Queued BLAST search completed.')

    # TODO: Result is currently kept for 500 seconds from the rqworker. See if we can decrease this value to decrease demand on resources
