import csv
import io
import os
import zipfile
import xmltodict
from typing import Any, Dict, List, Optional, Set, Tuple
from rest_framework.renderers import BaseRenderer
from barcode_blastn.file_paths import get_data_run_path
from barcode_blastn.models import BlastRun

from barcode_blastn.serializers import BlastDbSequenceExportSerializer,  BlastRunSerializer, HitSerializer

def get_letter(rank: str) -> str:
    '''
    Return the letter abbreviation of the taxonomic rank.
    E.g. `rank = "taxon_family" => "f"`
    '''
    assert rank.startswith('taxon_')
    if rank == 'taxon_superkingdom':
        return 'sk'
    else:
        return rank[6]

def renderDada2Taxonomy(data) -> List[str]:
    taxa_ranks = ['taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus']   
    sequence_file = []
    for sequence in data['sequences']:
        taxa = [sequence[rank] for rank in taxa_ranks]
        taxa_to_print = []
        for t in taxa:
                # stop printing as soon as rank unspecified
            if t is None:
                break
            taxa_to_print.append(t['scientific_name'].replace(' ', '_'))
        lineage = ';'.join(taxa_to_print)
        sequence_file.append(f'>{lineage};\n{sequence["dna_sequence"]}\n')
    return sequence_file

def renderDada2Species(data) -> List[str]:
    sequence_file = []
    for sequence in data['sequences']:
        id: str = sequence["version"].replace('.', '_')
        species = sequence['taxon_species']
        if species is None:
            genus = sequence['taxon_species']
            if not genus is None:
                species = f'{genus["scientific_name"]} {id}'
            else:
                species = f'{id} unidentified_species'
        else:
            species = species['scientific_name']
        sequence_file.append(f'>{id} {species}\n{sequence["dna_sequence"]}\n')
    return sequence_file

def renderSintax(data) -> List[str]:
    taxa_ranks = ['taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus']   
    sequence_file = []
    for sequence in data['sequences']:
        id: str = sequence["version"].replace('.', '_')
        taxa_to_print: List[str] = []
        for rank in taxa_ranks:
            if not sequence[rank] is None:
                name = sequence[rank]["scientific_name"].replace(' ', '_')
                taxa_to_print.append(f'{get_letter(rank)}:{name}')
        lineage = ','.join(taxa_to_print)
        sequence_file.append(f'>{id};tax={lineage};\n{sequence["dna_sequence"]}\n')
    return sequence_file

def renderDefaultFasta(data) -> List[str]:
    sequence_file = []
    for sequence in data['sequences']:
        sequence_file.append(f'>{sequence["accession_number"]}\n{sequence["dna_sequence"]}\n')
    return sequence_file

class BlastDbXMLRenderer(BaseRenderer):
    '''
    Return the entries of the blastdb in XML format.
    '''
    media_type = 'application/xml'
    format = 'xml'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        export_data = {'Database': data}
        return xmltodict.unparse(export_data).encode(self.charset)

class BlastDbFastaRenderer(BaseRenderer):
    '''
    Return the entries of blastdb in FASTA format.
    '''
    media_type = 'application/x-fasta'
    format = 'fasta'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context: Optional[Dict[Any, Any]]=None):
        fasta_format = renderer_context.get('fasta_format', '') if not renderer_context is None else ''
        
        sequence_file: List[str]
        if fasta_format == 'dada2tax':
            sequence_file = renderDada2Taxonomy(data)
        elif fasta_format == 'dada2sp':
            sequence_file = renderDada2Species(data)
        elif fasta_format == 'sintax':
            sequence_file = renderSintax(data)
        else:
            sequence_file = renderDefaultFasta(data)
            
        return ''.join(sequence_file).encode(self.charset)

class BlastDbCompatibleRenderer(BaseRenderer):
    '''
    Return the entries of blastdb in ZIP file format
    '''
    
    media_type = 'application/zip'
    format = 'zip'
    charset = 'utf-8'

    def get_rdp_rank_name(self, rank: str) -> str:
        assert rank.startswith('taxon_')
        return rank[6:]

    def render(self, data, accepted_media_type=None, renderer_context=None):
        fasta_format = renderer_context.get('fasta_format', '') if not renderer_context is None else ''

        zip_buffer = io.BytesIO()
        zip_files: List[Tuple[str,str]] = [] # a list of tuples, (filename, content)

        if fasta_format == 'qiime2':
            sequence_file = []
            taxonomy_file = []
            taxa_ranks = ['taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus']
            for sequence in data['sequences']:
                id: str = sequence["version"].replace('.', '_')
                found_empty: bool = False
                taxa_to_print = []
                sequence_file.append(f'>{id}\n{sequence["dna_sequence"]}\n')
                for rank in taxa_ranks:
                    if found_empty:
                        taxa_to_print.append(f'{get_letter(rank)}__')
                    elif sequence[rank] is None:
                        taxa_to_print.append(f'{get_letter(rank)}__')
                        found_empty = True
                    else:
                        taxa_to_print.append(f'{get_letter(rank)}__{sequence[rank]["scientific_name"].replace(" ", "_")}\n')
                taxonomy_file.append(f'{id}\t{"; ".join(taxa_to_print)}\n')
            zip_files.append((f'{data["custom_name"]}__q2.fasta', ''.join(sequence_file)))
            zip_files.append((f'{data["custom_name"]}__q2taxonomy.txt', ''.join(taxonomy_file)))
                
        elif fasta_format == 'dada2tax':
            sequence_file = renderDada2Taxonomy(data)
            zip_files.append((f'{data["custom_name"]}__dada2taxonomy.fasta', ''.join(sequence_file)))

        elif fasta_format == 'dada2sp':
            sequence_file = renderDada2Species(data)
            zip_files.append((f'{data["custom_name"]}__dada2species.fasta', ''.join(sequence_file)))

        elif fasta_format == 'sintax':
            sequence_file = renderSintax(data)
            zip_files.append((f'{data["custom_name"]}__sintax.fasta', ''.join(sequence_file)))

        elif fasta_format == 'rdp':
            rdp_taxa_ranks = ['taxon_superkingdom', 'taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus', 'taxon_species']
            sequence_file = []
            taxonomy_file = []
            taxonomy_file.append('1*Root*-1*0*rootrank\n')
            # set of taxomic names added already. key = taxonomic_name, id = taxid
            rdpTaxaAdded: Dict[str, int] = {}
            # next unused taxid
            nextRdpTaxId = -2
            for sequence in data['sequences']:
                id: str = sequence["version"].replace('.', '_')
                taxa_to_print: List[str] = ['Root']
                last_rank_id = 1
                for index, rank in enumerate(rdp_taxa_ranks):
                    taxon_name = taxa_to_print[-1] if sequence[rank] is None else sequence[rank]["scientific_name"].replace(' ', '_')
                    taxon_name = f'{get_letter(rank)}__{taxon_name}'
                    taxa_to_print.append(taxon_name)
                    if taxon_name not in rdpTaxaAdded:
                        # For taxon names that haven't been encountered before,
                        # write them to the taxonomy file
                        if not sequence[rank] is None:
                            taxid = sequence[rank]["id"]
                        else:
                            taxid = nextRdpTaxId
                            nextRdpTaxId = nextRdpTaxId - 1
                        taxonomy_file.append(f'{taxid}*{taxon_name}*{last_rank_id}*{index+1}*{self.get_rdp_rank_name(rank)}\n')
                        rdpTaxaAdded[taxon_name] = taxid
                        
                        last_rank_id = taxid 
                    else:
                        last_rank_id = rdpTaxaAdded[taxon_name]

                lineage = ';'.join(taxa_to_print)
                sequence_file.append(f'>{id}\t{lineage}\n{sequence["dna_sequence"]}\n')
            
            db_name = data["custom_name"].replace(' ', '_')
            zip_files.append((f'{db_name}__rdp.fasta', ''.join(sequence_file)))
            zip_files.append((f'{db_name}__rdp.txt', ''.join(taxonomy_file)))
            
        elif fasta_format == 'mothur':
            mothur_taxa_ranks = ['taxon_superkingdom', 'taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus', 'taxon_species']
            sequence_file = []
            taxonomy_file = []
            for sequence in data['sequences']: 
                id: str = sequence["version"].replace('.', '_')
                taxa_to_print = ['Root']
                for rank in mothur_taxa_ranks:
                    # stop printing as soon as rank unspecified
                    if sequence[rank] is None: 
                        break
                    name = sequence[rank]["scientific_name"].replace(' ', '_')
                    taxa_to_print.append(f"{get_letter(rank)}__{name}")
                lineage = ';'.join(taxa_to_print)
                sequence_file.append(f'>{id}\t{lineage};\n{sequence["dna_sequence"]}\n')
                taxonomy_file.append(f'{id}\t{lineage};\n')
            zip_files.append((f'{data["custom_name"]}__mothur.pds.fasta', ''.join(sequence_file)))
            zip_files.append((f'{data["custom_name"]}__mothur.pds.tax', ''.join(taxonomy_file)))

        else:
            sequence_file = renderDefaultFasta(data)
            zip_files.append((f'{data["custom_name"]}.fasta', ''.join(sequence_file)))
        
        with zipfile.ZipFile(zip_buffer, 'w') as zip_file:
            for file in zip_files:
                zip_file.writestr(file[0], file[1])
        
        return zip_buffer.getvalue()

def BlastDbDictWriter(data, delimiter: str = ','):
    response = io.StringIO()
    # fields displayed for each sequence entry
    sequence_fields = BlastDbSequenceExportSerializer.Meta.fields[:]
    sequence_fields = [field for field in sequence_fields if not field in ['id']]

    # For the following keys, flatten the keys so that there are no nested objects
    taxon_keys = [field for field in sequence_fields if field.startswith('taxon_')]
    def mutate_taxonomy(sequence):
        '''Remove nexted taxonomy objects from a sequence dictionary, and
        add information back as keys.'''
        for taxon_key in taxon_keys:
            # remove taxonomy object
            taxon_data = sequence.pop(taxon_key, None)
            # add back information as keys
            if taxon_data is None:
                sequence[f'{taxon_key}_id'] = None
                sequence[f'{taxon_key}_scientific_name'] = None
            else:
                sequence[f'{taxon_key}_id'] = taxon_data['id']
                sequence[f'{taxon_key}_scientific_name'] = taxon_data['scientific_name']
        return sequence
    sequence_data = [mutate_taxonomy(sequence) for sequence in data['sequences']]
    # Ensure that the newly added keys are added
    for taxon_key in taxon_keys:
        sequence_fields.extend([f'{taxon_key}_id', f'{taxon_key}_scientific_name'])
        sequence_fields.remove(taxon_key)
    # Write to file
    writer = csv.DictWriter(response, fieldnames=sequence_fields, extrasaction='ignore', dialect='unix', delimiter=delimiter)
    writer.writeheader()
    writer.writerows(sequence_data)
    return response

class BlastDbTSVRenderer(BaseRenderer):
    '''
    Return the entries of blastdb in TSV format.
    '''
    media_type = 'text/tsv'
    format = 'tsv'
    charset = 'utf-8'
    def render(self, data, accepted_media_type=None, renderer_context=None):
        response = BlastDbDictWriter(data=data, delimiter='\t')
        return response.getvalue().encode(self.charset)

class BlastRunFastaRenderer(BaseRenderer):
    '''
    Return input file of run in FASTA format 
    '''
    media_type = 'application/x-fasta'
    format = 'fasta'
    charset = 'utf-8'
    
    def render(self, data, accepted_media_type=None, renderer_context: Optional[Dict[Any, Any]] = None):
        sequence_file = []
        for sequence in data['queries']:
            if not sequence["original_species_name"] is None:
                sequence_file.append(f'>{sequence["definition"]}\t{sequence["original_species_name"].replace(" ", "_")}\n{sequence["query_sequence"]}\n')
            else:
                sequence_file.append(f'>{sequence["definition"]}\n{sequence["query_sequence"]}\n')
        return ''.join(sequence_file).encode(self.charset)

class BlastDbCSVRenderer(BaseRenderer):
    '''
    Return the blastdb in CSV format
    '''
    media_type = 'text/csv'
    format = 'csv'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        response = BlastDbDictWriter(data=data, delimiter=',')
        return response.getvalue().encode(self.charset)

class BlastRunHitsTxtRenderer(BaseRenderer):
    '''
    Return the blast run hits in txt format
    '''
    media_type = 'text/plain'
    format = 'txt'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        out_lines = []
        run_id = data['id']
        run_folder = get_data_run_path(run_id)

        # only print results if run is finished and file exists
        if data['status'] == BlastRun.JobStatus.FINISHED:
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
        
def blast_run_hit_dict_writer(data, delimiter=',') -> io.StringIO:
    '''
    Given data for all BLAST hits from a blast run, write the data
    to the response with columns delimited by `delimiter` character using
    `csv.DictWriter()`. Return the response.
    '''
    response = io.StringIO()
        
    # get a list of displayed fields to display for each hit
    hit_fields = HitSerializer.Meta.fields[:]
    hit_fields = [field for field in hit_fields if field != 'owner_run' and field != 'db_entry']
    db_entry_fields = ['accession_number', 'organism', 'version', 'country', 'specimen_voucher']
    hit_fields.extend(db_entry_fields)

    # make a DictWriter to write data
    writer = csv.DictWriter(response, delimiter=delimiter, fieldnames=hit_fields, extrasaction='ignore', dialect='unix')
    writer.writeheader()
    # add the fields under db_entry to the parent dictionary for output
    for query in data['queries']:
        for hit in query['hits']:
            hit.update(hit['db_entry'])
                
        writer.writerows(query['hits'])
    return response

class BlastRunHitsCSVRenderer(BaseRenderer):
    '''
    Return the information of a blast run hits in CSV format
    '''
    media_type = 'text/csv'
    format = 'csv'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        return blast_run_hit_dict_writer(data, delimiter=',').getvalue().encode(self.charset)

class BlastRunHitsTSVRenderer(BaseRenderer):
    '''
    Return the information of a blast run hits in TSV format
    '''
    media_type = 'text/tsv'
    format = 'tsv'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        return blast_run_hit_dict_writer(data, delimiter='\t').getvalue().encode(self.charset)

def blast_run_taxonomy_dict_writer(data, delimiter=','):
    '''
    Given data on taxonomic assignments from a blast run, write the data
    to the response with columns delimited by `delimiter` character using
    `csv.DictWriter()`. Return the response.
    '''
    response = io.StringIO()
    # get a list of displayed fields to display for each hit
    query_fields = ['definition', 'tree_query_id', 'original_species_name', 'results_species_name', 'accuracy_category', 'highest_percent_identity', 'evalue']

    def populate_fields(old_dict: dict):
        query_id = old_dict.pop('write_tree_identifier')
        old_dict['tree_query_id'] = query_id
        return old_dict

    finalData = [populate_fields(q) for q in data['queries']]
    # make a DictWriter to write data
    writer = csv.DictWriter(response, delimiter=delimiter, fieldnames=query_fields, extrasaction='ignore', dialect='unix')
    writer.writeheader()
    # add the fields under db_entry to the parent dictionary for output
    writer.writerows(finalData)
    return response

class BlastRunTaxonomyCSVRenderer(BaseRenderer):
    '''
    Return information on taxonomic assignments from run in CSV format.
    '''
    media_type = 'text/csv'
    format = 'csv'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        return blast_run_taxonomy_dict_writer(data, delimiter=',').getvalue().encode(self.charset)

class BlastRunTaxonomyTSVRenderer(BaseRenderer):
    '''
    Return information on taxonomic assignments from run in TSV format.
    '''
    media_type = 'text/tsv'
    format = 'tsv'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        return blast_run_taxonomy_dict_writer(data, delimiter='\t').getvalue().encode(self.charset)

class BlastRunHTMLRenderer(BlastRunHitsTxtRenderer):
    '''
    Return BLAST text results as an HTML page
    '''
    media_type = 'text/plain'
    format = 'html'
    charset = 'utf-8'

    def render(self, data, accepted_media_type=None, renderer_context=None):
        return super().render(data, accepted_media_type, renderer_context)