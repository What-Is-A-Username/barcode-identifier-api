import copy
from io import StringIO
import os
import shutil
from typing import Any, Tuple
from ratelimit import limits
from shlex import quote
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime, timezone
import subprocess
from typing import Dict, List
from barcode_blastn.file_paths import get_data_fishdb_path, get_data_run_path, get_ncbi_folder, get_static_run_path
from barcode_blastn.helper.assign_accuracy import annotate_accuracy_category
from barcode_blastn.helper.calculate_distance import calculate_genetic_distance
from barcode_blastn.helper.modified_clustalo import getMultipleAlignmentResult, submitMultipleAlignmentAsync
from barcode_blastn.helper.parse_gb import retrieve_gb, save_taxonomy
from barcode_blastn.helper.parse_results import parse_results
from rest_framework import status
from barcode_blastn.helper.read_tree import readDbTreeFromFile, readHitTreeFromFile
from barcode_blastn.models import BlastDb, BlastQuerySequence, BlastRun, Hit, NuccoreSequence 
from django.db.models import QuerySet
from barcode_identifier_api.celery import app

# Max time for any task to be running for
HARD_TIME_LIMIT_IN_SECONDS = 900

def raise_error(run: BlastRun, err: str) -> None:
    '''
    Update the blast run with the error message and status
    specified by `err` and print it to console.
    '''
    print(err)
    run.errors = run.errors + "\n" + err
    run.status = BlastRun.JobStatus.ERRORED
    run.save()
    return

def mark_run_complete(run: BlastRun) -> None:
    '''Mark the run status as complete'''
    run.status = BlastRun.JobStatus.FINISHED
    run.end_time = datetime.now()
    run.save()
    return 

@app.task(time_limit=HARD_TIME_LIMIT_IN_SECONDS)  
def run_blast_command(ncbi_blast_version: str, fishdb_id: str, run_id: str) -> bool:   
    '''Performs a BLASTn search using the specified ncbi_blast_version, database with id fishdb_id, for the run of id run_id

    Return true if operation was successful.
    '''
    print('Beginning queued BLAST search ...')

    project_root = os.path.abspath(os.path.dirname('./'))
    print(f'Working directory is {project_root}')
    print(f'Using blast version {ncbi_blast_version}')
    print(f'Using database id {fishdb_id} and run id {run_id}')
    
    try:
        run_details : BlastRun = BlastRun.objects.get(id = run_id)
    except BlastRun.DoesNotExist as exc:
        raise RuntimeError(f'BlastRun of specified id {run_id} does not exist')

    run_details.status = BlastRun.JobStatus.STARTED
    run_details.start_time = datetime.now()
    run_details.save()

    blast_root = get_ncbi_folder(ncbi_blast_version=ncbi_blast_version)
    command_app = blast_root + '/blastn'
    # check if blast executable are present
    if not os.path.exists(command_app) or not os.path.isfile(command_app):
        raise_error(run_details, f"Critical error: Failed to find blastn executable at {command_app}")
        raise FileNotFoundError('Failed to find blastn executable')

    try:
        db: BlastDb = BlastDb.objects.get(id=fishdb_id)
    except BlastDb.DoesNotExist as exc:
        raise RuntimeError('BlastDb of specified id does not exist')
    db_path = get_data_fishdb_path(db) + '/database'
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
       
    print('Results written. Now saving hits ...')
    
    # bulk create hits
    try:
        best_hits: Dict[str, List[int]]
        parsed_data: List[Dict[str, Any]]
        parsed_data, best_hits = parse_results(out)
        accessions = [hit['subject_accession_version'] for hit in parsed_data]
        accession_entries = NuccoreSequence.objects.filter(version__in=accessions, owner_database=run_details.db_used)

        # Map query id -> query object
        query_entries: List[BlastQuerySequence] = list(run_details.queries.all())
        query_entry_ids = [q.query_seq_identifier() for q in query_entries]
        query_entries_dict: dict[str, BlastQuerySequence] = dict(zip(query_entry_ids, query_entries))

        for i in range(len(parsed_data)):
            parsed_data[i].update({
                'query_sequence': query_entries_dict[parsed_data[i]['query_accession_version']],
                'db_entry': accession_entries.get(version=parsed_data[i]['subject_accession_version'])
            })
        all_hits = [Hit(**hit) for hit in parsed_data]
        hit_objects = Hit.objects.bulk_create(all_hits)
    except BaseException as exc:
        raise_error(run_details, f"Errored while bulk creating hit results.")
        raise RuntimeError('Errored while bulk creating hit results.') from exc
    
    print('Hits saved. Performing taxonomy assignments based on best hit ...')
        
    # assign identities to query sequences
    for query in query_entries:
        best_hit: List[int] = best_hits.get(query.query_seq_identifier(), [])
        if len(best_hit) > 0:
            taxa = []
            for hit_index in best_hit:
                hit = parsed_data[hit_index]
                assert hit.get('best_hit') == True
                seq: NuccoreSequence = hit.get('db_entry', None)
                assert not seq is None
                if seq.taxon_species is None:
                    taxa.append(f'{seq.version}_unspecified_species')
                else:
                    taxa.append(seq.taxon_species.scientific_name)
            taxa = list(set(taxa))
            query.results_species_name = ', '.join(taxa)
        else:
            # If there was no hits on the reference sequences,
            # leave the result blank
            query.results_species_name = ''
    BlastQuerySequence.objects.bulk_update(query_entries, fields=['results_species_name'])

    print('Updated and saved assignments. Cleaning up')

    # update run information
    try:
        run_details.errors = err
        run_details.save()
    except BaseException as exc:
        raise_error(run_details, f"Errored while updating run errors.")
        raise RuntimeError('Errored while updating run errors.') from exc

    print('Queued BLAST search and assignment completed.')
    print('Next parameters: ', run_details.create_hit_tree, run_details.create_db_tree)
    
    if run_details.create_hit_tree or run_details.create_db_tree:
        performAlignment(True, run=run_details)
    mark_run_complete(run_details)
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
        new_data = save_taxonomy(taxonomy_info=new_data, user=None)
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
        
def performAlignment(blast_successful: bool, run: BlastRun) -> bool:
    '''Perform tree construction for the given run specified by run_id.
        Return True if the operation was successful, False otherwise.
    '''
    run_id = str(run.id)
    print(f'Performing alignment for {run_id}')

    if not blast_successful:
        return False
    else:
        print(f"Starting alignment submission for run {run_id}")
        
    # Return if no construction needs to occur
    if not run.create_db_tree and not run.create_hit_tree:
        return True

    # Gather sequences (query sequences + hit sequences) into a single list
    run_folder = get_data_run_path(run_id=run_id)
    query_sequences: List[SeqIO.SeqRecord] = list(SeqIO.parse(f'{run_folder}/alignment_query.fasta', 'fasta'))

    def gatherSequences(reference_sequences: List[NuccoreSequence], query_sequences: List[SeqIO.SeqRecord] = []) -> List[SeqIO.SeqRecord]:
        '''Creates a new List of SeqRecords created by making a new deep copy of existingSequences and SeqRecords
        made from the NuccoreSequences in sequences
        '''
        concatSequences = copy.deepcopy(query_sequences)
        sequences_added = [q.description for q in concatSequences]
        for seq in reference_sequences:
            seq_id = seq.write_tree_identifier()
            if seq_id not in sequences_added:
                concatSequences.append(SeqIO.SeqRecord(
                    seq=Seq(seq.dna_sequence), 
                    id=seq_id,
                    description=seq.definition
                ))
                sequences_added.append(seq_id)
        return concatSequences
    
    if run.create_hit_tree:
        # Construct a list of all hit and query sequences
        queryAndHitSequences = query_sequences[:]
        # Query over all hits
        hits = []
        hit_query = Hit.objects.filter(query_sequence__owner_run_id=run_id)
        hits = [h.db_entry for h in hit_query]
        queryAndHitSequences = gatherSequences(hits, query_sequences)

        sequence_string = StringIO()
        SeqIO.write(queryAndHitSequences, sequence_string, 'fasta')

        hit_job_result = completeAlignment(sequence_string.getvalue(), run_id)
        run.alignment_job_id = hit_job_result[0]
        if hit_job_result[-1] != status.HTTP_200_OK:
            run.throw_error(hit_job_result[1])
            raise_error(run, 'Encountered error creating hit alignment/tree')
            return False
    
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
            raise_error(run, 'Encountered error creating db alignment/tree')
            return False

    try:
        db_tree = readDbTreeFromFile(run) if run.create_db_tree else ''
        if db_tree is None:
            raise_error(run, 'Could not locate result file to read db tree string from.')
            return False
        else:
            run.db_tree = db_tree
            run.save()
        hit_tree = readHitTreeFromFile(run) if run.create_hit_tree else ''
        if hit_tree is None:
            raise_error(run, 'Could not locate result file to read hit tree string from.')
            return False
        else:
            run.hit_tree = hit_tree
            run.save()
    except BaseException as exc:
        raise_error(run, f"Errored while reading alignment trees from file.")
        return False
    else:
        print('Successfully finished tree construction')
        if run.create_hit_tree:
            classify_genetic_distance(alignment_successful=True, run=run)
        return True

def completeAlignment(sequence_string: str, run_id: str) -> Tuple[str, str, int]:
    '''
        Submit a .fasta string for submission as a synchronous job to Clustal Omega.
        Return a three-item tuple, containing the job ID, status message, and status number.
        The job was successful if status message is 200 (OK).
    '''
    print("Considering tree alignment for " + run_id)
    # The following code submits the fasta file to ClustalO at EMBL-EBI for multiple sequence alignment, and is adapted from clustalo.py at https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation
    # ClustalO:
    #   Sievers F., Wilm A., Dineen D., Gibson T.J., Karplus K., Li W., Lopez R., McWilliam H., Remmert M., Söding J., Thompson J.D. and Higgins D.G. (2011)
    # Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. 
    # Mol. Syst. Biol. 7:539
    # PMID:  21988835 
    # DOI:  10.1038/msb.2011.75

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

    # save clustal files locally
    try:
        getMultipleAlignmentResult(job_id=job_id, run_id=run_id)
    except BaseException:
        return (job_id, 'Encountered unexpected error retrieving results', status.HTTP_500_INTERNAL_SERVER_ERROR)
    else:       
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
        except BaseException as exc:
            return (job_id, 'Failed to process output files', status.HTTP_500_INTERNAL_SERVER_ERROR)
        
    # return success
    return (job_id, 'Successfully ran alignment job', status.HTTP_200_OK)

# @limits(calls = 1, period = 10)
# @app.task(time_limit=HARD_TIME_LIMIT_IN_SECONDS)
def classify_genetic_distance(alignment_successful: bool, run: BlastRun) -> bool:
    run_id = str(run.id)
    if alignment_successful:
        print(f'Running db accuracy classification for {run_id}')
    else:
        print(f'Skipping db classification for {run_id} because alignment was unsuccessful')
        return True

    # we can only calculate distances with a complete alignment of hits with query sequences
    if not run.create_hit_tree:
        print('Will not classify distances because no alignment performed.')
        return True

    alignment_file_path = f'{get_static_run_path(run_id)}/{run.alignment_job_id}.aln-clustal_num.clustal_num'
    if os.path.exists(alignment_file_path) and os.path.isfile(alignment_file_path):
        try:
            annotation_result = annotate_accuracy_category(run, alignment_file_path)
            if annotation_result:
                print('Successfully finished annotating accuracy')
            return annotation_result
        except BaseException:
            raise_error(run=run, err='Errored while annotating accuracy category.')
            return False
            
    else:
        raise_error(run, f'Could not find multiple alignment file at {alignment_file_path}\n' + run.errors)
        return False