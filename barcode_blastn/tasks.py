import copy
from io import StringIO
import os
import shutil
from typing import Tuple
from ratelimit import limits
from shlex import quote
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime, timezone
import subprocess
from typing import Dict, List
from barcode_blastn.file_paths import get_data_fishdb_path, get_data_run_path, get_ncbi_folder, get_static_run_path
from barcode_blastn.helper.modified_clustalo import getMultipleAlignmentResult, submitMultipleAlignmentAsync
from barcode_blastn.helper.parse_gb import retrieve_gb
from barcode_blastn.helper.parse_results import parse_results
from rest_framework import status
from barcode_blastn.helper.read_tree import readDbTreeFromFile, readHitTreeFromFile
from barcode_blastn.models import BlastDb, BlastRun, Hit, NuccoreSequence 
from django.db.models import QuerySet
from barcode_identifier_api.celery import app

# Max time
HARD_TIME_LIMIT_IN_SECONDS = 600

# update the blast run with the errors and print to console
def raise_error(run: BlastRun, err: str) -> None:
    print(err)
    run.errors = run.errors + "\n" + err
    run.status = BlastRun.JobStatus.ERRORED
    run.save()
    return

@app.task(time_limit=HARD_TIME_LIMIT_IN_SECONDS)  # hard time limit of 30 seconds 
def run_blast_command(ncbi_blast_version: str, fishdb_id: str, run_id: str) -> bool:   
    '''Performs a BLASTn search using the specified ncbi_blast_version, database with id fishdb_id, for the run of id run_id

    Return true if operation was successful.
    '''
    print('Beginning queued BLAST search ...')

    project_root = os.path.abspath(os.path.dirname('./'))
    print(f'Working directory is {project_root}')
    print(f'Using blast version {ncbi_blast_version}')
    print(f'Using database id {fishdb_id} and run id {run_id}')
    
    # TODO: throw error if blastrun DNE
    try:
        run_details : BlastRun = BlastRun.objects.get(id = run_id)
    except BlastRun.DoesNotExist as exc:
        print(f"Critical error: BlastRun of id {run_id} does not exist.")
        raise RuntimeError('BlastRun of specified id does not exist')

    run_details.status = BlastRun.JobStatus.STARTED
    run_details.start_time = datetime.now()
    run_details.save()

    blast_root = get_ncbi_folder(ncbi_blast_version=ncbi_blast_version)
    command_app = blast_root + '/blastn'
    # check if blast executable are present
    if not os.path.exists(command_app) or not os.path.isfile(command_app):
        raise_error(run_details, f"Critical error: Failed to find blastn executable at {command_app}")
        raise FileNotFoundError('Failed to find blastn executable')

    db_path = get_data_fishdb_path(fishdb_id) + '/database'
    # check if the database file exists
    if not os.path.exists(db_path + '.fasta') or not os.path.isfile(db_path + '.fasta'):
        raise_error(run_details, f"Critical error: Failed to find BLAST database at {db_path}")
        raise FileNotFoundError('Failed to find BLAST database')

    run_folder = get_data_run_path(run_id=run_id)
    results_file = f'{run_folder}/results.txt'
    errors_file = f'{run_folder}/errors.txt'
    query_file = f'{run_folder}/query.fasta'

    # check if run status folder exists
    if not os.path.exists(run_folder) or not os.path.isdir(run_folder):
        raise_error(run_details, f"Critical error: Failed to find run folder at {run_folder}")
        raise FileNotFoundError('Failed to find run folder')
    # check if query file exists
    if not os.path.exists(query_file) or not os.path.isfile(query_file):
        raise_error(run_details, f"Critical error: Failed to find query sequence file at {query_file}")
        raise FileNotFoundError('Failed to find query sequence file')

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
    except BaseException as exc:
        raise_error(run_details, f"Errored while retrieving shell output.")
        raise RuntimeError('Errored while retrieving shell output.') from exc

    print(out)

    print('Writing results to file ...')

    try:
        with open(results_file, "w") as results_out:
            results_out.write(out + '\n')
        if err:
            with open(errors_file, "w") as errors_out:
                errors_out.write(err)
    except BaseException as exc:
        raise_error(run_details, f"Errored while handling results output.")
        raise RuntimeError('Errored while handling results output.') from exc
       
    # bulk create hits
    try:
        parsed_data = parse_results(out)
        accessions = [hit['subject_accession_version'] for hit in parsed_data]
        accession_entries = NuccoreSequence.objects.filter(accession_number__in=accessions, owner_database=run_details.db_used)
        all_hits = [Hit(**hit, owner_run = run_details, db_entry = accession_entries.get(accession_number=hit['subject_accession_version'])) for hit in parsed_data]
        Hit.objects.bulk_create(all_hits)
    except BaseException as exc:
        raise_error(run_details, f"Errored while bulk creating hit results.")
        raise RuntimeError('Errored while bulk creating hit results.') from exc

    # update run information
    try:
        run_details.errors = err
        run_details.save()
    except BaseException as exc:
        raise_error(run_details, f"Errored while updating run errors.")
        raise RuntimeError('Errored while updating run errors.') from exc

    print('Queued BLAST search completed.')
    return True

PERIOD = 3500 # Time between calls, in number of seconds
@limits(calls = 1, period = PERIOD)
@app.task(time_limit=HARD_TIME_LIMIT_IN_SECONDS)  # hard time limit of 30 seconds 
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
        

PERIOD = 10 # Time between calls, in number of seconds
@limits(calls = 1, period = PERIOD)
@app.task(time_limit=HARD_TIME_LIMIT_IN_SECONDS)
def performAlignment(blast_successful: bool, run_id: str) -> Tuple[str, int]:
    '''Perform tree construction for the given run specified by run_id.
        Return a three item tuple, representing the status message and status number.
        The job was successful if status number is 201 or 200.
    '''

    if not blast_successful:
        return ('BLAST returned error code', status.HTTP_500_INTERNAL_SERVER_ERROR)
    else:
        print(f"Starting alignment submission for run {run_id}")

    # Check that the run exists
    try:
        run : BlastRun = BlastRun.objects.get(id=run_id)
    except BlastRun.DoesNotExist:
        return ('Could not find a run result with the given id.', status.HTTP_400_BAD_REQUEST)
        
    # Return if no construction needs to occur
    if not run.create_db_tree and not run.create_hit_tree:
        run.status = BlastRun.JobStatus.FINISHED
        run.end_time = datetime.now()
        run.save()
        return ('', status.HTTP_200_OK)

    # Gather sequences (query sequences + hit sequences) into a single list
    run_folder = get_data_run_path(run_id=run_id)
    query_sequences: List[SeqIO.SeqRecord] = list(SeqIO.parse(f'{run_folder}/query.fasta', 'fasta'))
    for i in range(len(query_sequences)):
        query_sequences[i].id += '|query'

    def gatherSequences(sequences: List[NuccoreSequence], existingSequences: List[SeqIO.SeqRecord] = []) -> List[SeqIO.SeqRecord]:
        '''Creates a new List of SeqRecords created by making a new deep copy of existingSequences and SeqRecords
        made from the NuccoreSequences in sequences
        '''
        concatSequences = copy.deepcopy(existingSequences)
        sequences_added = [q.id for q in concatSequences]
        for seq in sequences:
            seq_id = f'{seq.accession_number}|{seq.organism}'.replace(' ', '_')
            if seq_id not in sequences_added:
                concatSequences.append(SeqIO.SeqRecord(
                    seq=Seq(seq.dna_sequence), 
                    id= seq_id,
                    description=seq.definition
                ))
                sequences_added.append(seq_id)
        return concatSequences
    
    if run.create_hit_tree:
        # Construct a list of all hit and query sequences
        queryAndHitSequences = query_sequences[:]
        # Query over all hits
        hits = [h.db_entry for h in list(run.hits.all())] # type: ignore
        queryAndHitSequences = gatherSequences(hits, query_sequences)

        sequence_string = StringIO()
        SeqIO.write(queryAndHitSequences, sequence_string, 'fasta')

        hit_job_result = completeAlignment(sequence_string.getvalue(), run_id)
        run.alignment_job_id = hit_job_result[0]
        if hit_job_result[-1] != status.HTTP_200_OK:
            run.throw_error(hit_job_result[1])
            return (hit_job_result[1], hit_job_result[2]) 
    
    if run.create_db_tree:
    # Construct a list of query sequences and all database sequences
        allSequences = query_sequences[:]
        blast_db = run.db_used
        db_sequences = list(blast_db.sequences.all()) # type:ignore
        allSequences = gatherSequences(db_sequences, query_sequences)
        
        sequence_string = StringIO()
        SeqIO.write(allSequences, sequence_string, 'fasta')

        all_job_result = completeAlignment(sequence_string.getvalue(), run_id)
        run.complete_alignment_job_id = all_job_result[0]
        if all_job_result[-1] != status.HTTP_200_OK:
            run.throw_error(all_job_result[1])
            return (all_job_result[1], all_job_result[2]) 

    try:
        run.status = BlastRun.JobStatus.FINISHED
        run.db_tree = readDbTreeFromFile(run) if run.create_db_tree else ''
        run.hit_tree = readHitTreeFromFile(run) if run.create_hit_tree else ''
        run.end_time = datetime.now()
        run.save()
    except BaseException as exc:
        raise_error(run, f"Errored while updating run status.")
        raise RuntimeError('Errored while updating run status.') from exc

    return ('', status.HTTP_201_CREATED)




def completeAlignment(sequence_string: str, run_id: str) -> Tuple[str, str, int]:
    '''
        Submit a .fasta string for submission as a synchronous job to Clustal Omega.
        Return a three-item tuple, containing the job ID, status message, and status number.
        The job was successful if status message is 200 (OK).
    '''

    # The following code submits the fasta file to ClustalO at EMBL-EBI for multiple sequence alignment, and is adapted from clustalo.py at https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation
    # ClustalO:
    # TODO: Add citation to web interface
    #   Sievers F., Wilm A., Dineen D., Gibson T.J., Karplus K., Li W., Lopez R., McWilliam H., Remmert M., Söding J., Thompson J.D. and Higgins D.G. (2011)
    # Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. 
    # Mol. Syst. Biol. 7:539
    # PMID:  21988835 
    # DOI:  10.1038/msb.2011.75
    # TODO: Find citation for EMBL-EBI

    # This runs the tool, mimicking the following command line calls from the original Python wrapper
    # python clustalo.py --asyncjob --email email@domain.com ./aggregate.fasta
    # This checks if its done
    # python clustalo.py --status --jobid clustalo-R20230207-011805-0893-96632125-p2m
    # This outputs the data 
    # python clustalo.py --polljob --jobid clustalo-R20230207-011805-0893-96632125-p2m
    # job_id: str = submitClustalOAsyncJob(f'runs/{run_id}/aggregate.fasta', run_id)

    job_id: str
    try:
        job_id = submitMultipleAlignmentAsync(sequence=sequence_string, run_id=run_id)
    except BaseException:
        return ('', 'Encountered unexpected error submitting job', status.HTTP_500_INTERNAL_SERVER_ERROR)
        # tree.internal_status = ResultTree.TreeStatus.ERRORED
        # tree.save()

    # tree.internal_status = ResultTree.TreeStatus.PROCESSING
    # tree.save()
    # save clustal files locally
    try:
        getMultipleAlignmentResult(job_id=job_id, run_id=run_id)
    except BaseException:
        return (job_id, 'Encountered unexpected error retrieving results', status.HTTP_500_INTERNAL_SERVER_ERROR)
    else:
        # TODO: Confirm operation was successful by checking for network error, exceptions
        
        # move files to be downloaded
        try:
            destination_dir = get_static_run_path(run_id)
            parent_folder = get_data_run_path(run_id)
            files = os.listdir(parent_folder)
            files_to_transfer = ['.aln-clustal_num.clustal_num', '.phylotree.ph', '.pim.pim', '.sequence.txt']
            for file in files:
                if any([file.endswith(extension) for extension in files_to_transfer]):
                    print(f"Moving {file}")
                    shutil.copy(f'{parent_folder}/{file}', destination_dir)
            # tree.internal_status = ResultTree.TreeStatus.FINISHED
        except BaseException as exc:
            return (job_id, 'Failed to process output files', status.HTTP_500_INTERNAL_SERVER_ERROR)
        # except BaseException:
            # tree.internal_status = ResultTree.TreeStatus.ERRORED
        # tree.save()

    # return success
    return (job_id, 'Successfully ran alignment job', status.HTTP_200_OK)

