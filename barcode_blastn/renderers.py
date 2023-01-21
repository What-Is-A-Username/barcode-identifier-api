import csv
import io
import os
from django.http.response import HttpResponse
from rest_framework.renderers import BaseRenderer
from barcode_blastn.models import BlastRun

from barcode_blastn.serializers import BlastDbSequenceEntrySerializer,  BlastRunSerializer, HitSerializer

'''
Return the entries of blastdb in FASTA format (accession number + sequence)
'''
class BlastDbFastaRenderer(BaseRenderer):
    media_type = 'text/x-fasta'
    format = 'fasta'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        fasta_file = []
        for sequence in data['sequences']:
            fasta_file.append(f'>{sequence["accession_number"]}\n{sequence["dna_sequence"]}\n')
        return ''.join(fasta_file).encode(self.charset)

'''
Return input file of run in FASTA format 
'''
class BlastRunFastaRenderer(BaseRenderer):
    media_type = 'text/x-fasta'
    format = 'fasta'
    charset = 'utf-8'
    
    def render(self, data, accepted_media_type=None, renderer_context=None):
        return data

'''
Return the blastdb in CSV format
'''
class BlastDbCSVRenderer(BaseRenderer):
    media_type = 'text/csv'
    format = 'csv'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        response = io.StringIO()
       
        # fields displayed for each sequence entry
        sequence_fields = BlastDbSequenceEntrySerializer.Meta.fields[:]
        sequence_fields = [field for field in sequence_fields if not field in ['id']]

        comment_writer = csv.writer(response)
        comment_writer.writerow([f'# Barcode Identifier API'])
        comment_writer.writerow([f"# Database: {data['custom_name']}"])
        comment_writer.writerow([f"# Description: {data['description']}"])
        
        # fields displayed for each sequence
        writer = csv.DictWriter(response, fieldnames=sequence_fields, extrasaction='ignore', dialect='unix')

        writer.writeheader()
        writer.writerows(data['sequences'])

        return response.getvalue().encode(self.charset)

'''
Return the blast run results in txt format
'''
class BlastRunTxtRenderer(BaseRenderer):
    media_type = 'text/plain'
    format = 'txt'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        out_lines = ['# Barcode Identifier API\n']
        run_id = data['id']
        
        # get a list of fields to display as comments at the start of the file
        comment_fields = [ field for field in BlastRunSerializer.Meta.fields if field != 'hits' and field != 'owner_run']

        # write the comments
        for comment in comment_fields:
            if comment == 'db_used':
                out_lines.append(f'# db_used-id={data[comment]["id"]}\n')
                out_lines.append(f'# db_used-custom_name={data[comment]["custom_name"]}\n')
                out_lines.append(f'# db_used-description={data[comment]["description"]}\n')
            else:
                val : str = data[comment]
                if val is None:
                    val = 'NULL'
                out_lines.append(f'# {comment}={val}\n')

        out_lines.append('\n')

        run_folder = os.path.abspath(f'./runs/{run_id}/')
        # only print results if run is finished and file exists
        if data['job_status'] == BlastRun.JobStatus.FINISHED:
            results_file = f'{run_folder}/results.txt'
            if os.path.exists(results_file) and os.path.isfile(results_file):
                with open(results_file, 'r') as results_txt_file:
                    out_lines.extend(results_txt_file.readlines())

                # omit "# Database" line to prevent leaking of file paths
                out_lines = [line for line in out_lines if not line.startswith('# Database: ')]
            else:
                out_lines.append('Error: Could not find results file generated from this run')

        else: 
            out_lines.append('Error: Job has not yet completed so hits cannot be compiled at this time.')

        return ''.join(out_lines).encode(self.charset)
        

'''
Return the information of a blast run in CSV format
'''
class BlastRunCSVRenderer(BaseRenderer):
    media_type = 'text/csv'
    format = 'csv'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        response = io.StringIO()

        # get a list of fields to display as comments at the start of the file
        comment_fields = [ field for field in BlastRunSerializer.Meta.fields if field != 'hits']

        # write the comments
        comment_writer = csv.writer(response)
        comment_writer.writerow([f'# Barcode Identifier API'])
        for comment in comment_fields:
            if comment == 'db_used':
                comment_writer.writerow([f'# db_used-id={data[comment]["id"]}'])
                comment_writer.writerow([f'# db_used-custom_name={data[comment]["custom_name"]}'])
                comment_writer.writerow([f'# db_used-description={data[comment]["description"]}'])
            else:
                val : str = data[comment]
                if val is None:
                    val = 'NULL'
                comment_writer.writerow([f'# {comment}={val}'])
        
        # get a list of displayed fields to display for each hit
        hit_fields = HitSerializer.Meta.fields[:]
        hit_fields = [field for field in hit_fields if field != 'owner_run' and field != 'db_entry']
        db_entry_fields = ['accession_number', 'definition', 'organism', 'isolate', 'country', 'specimen_voucher']
        hit_fields.extend(db_entry_fields)

        # add the fields under db_entry to the parent dictionary for output
        for hit in data['hits']:
            db_entry_values = [hit['db_entry'][field] for field in db_entry_fields]
            hit.update(zip(db_entry_fields, db_entry_values))

        # make a DictWriter to write data
        writer = csv.DictWriter(response, fieldnames=hit_fields, extrasaction='ignore', dialect='unix')

        # write csv header row + all hits
        writer.writeheader()
        writer.writerows(data['hits'])

        return response.getvalue().encode(self.charset)

class BlastRunHTMLRenderer(BlastRunTxtRenderer):
    media_type = 'text/plain'
    format = 'html'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        return super().render(data, accepted_media_type, renderer_context)