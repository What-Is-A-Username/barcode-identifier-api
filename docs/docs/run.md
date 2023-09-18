# Running queries with Barrel

Queries can be submitted to Barrel via the web app, on the `/blast` page.

<img src='../images/blast_run.png' width='80%' alt='Drawing'>
</img>

## Query sequences

Query sequences can be uploaded in four ways:

**Method 1**. Upload of sequences in a FASTA format 

```
>query1 Gymnotus cylindricus
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
AACCCGGAGCCCTCCTCGGGGACGACCAAATTTATAATGTAATTGTTACTGCCCACGC
>query2 Gymnotus esmeraldas
ATAGTATTTGGTGCCTGAGCTGGAATAGTTGGCACAGCCTTGAGCCTACTGATCCGAGCAGAACTAAGCC
AACCCGGAACCCTCCTAGGCGATGACCAAATTTATAATGTAATCGTTACTGCCCACGC
```

**Method 2**. Pasting raw text of a FASTA file, or a single string of nucleotides to be interpreted as a single query sequence

```
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
AACCCGGAGCCCTCCTCGGGGACGACCAAATTTATAATGTAATTGTTACTGCCCACGC
```
   
*Note: this sequence will be identified as "query_sequence" by Barrel.*

**Method 3**. Upload of GenBank accession numbers or accession versions in a file, with one identifier per line
```
ON303390
ON303391
```
   

**Method 4**. Pasting a list of GenBank accession numbers of accession versions in a file, one per line.
```
ON303390
ON303391
``` 
   

## Nucleotide BLAST Parameters

Each run queries a single BLAST database specified in this section. BLAST is run with the command line arguments `-outfmt 7`, `-query` and `-db`, which supply input and determine output. 

The web app allows you to narrow down to the desired database using the marker gene and reference library it belongs to. Recall that each reference library and each BLAST database is uniquely identified by its universal unique identifier (UUID). Example: `ae47cd0a-b7b7-4dfe-9223-8e026d2980c5`.

All jobs will run BLAST. The BLAST top hit for each sequence, as determined by percent identity, is used for taxonomic assignment.

## Multiple Alignment and Tree Parameters

### Create hit tree
If enabled (checked or set to "true"), a multiple sequence alignment job will be run with ClustalOmega (hosted at EMBL-EBI) with default parameters. The submitted sequences for alignment will consist of all query sequences, along with all hits they collectively returned in the database.

Additionally, Barrel will calculate the Kimura-2-parameter distance between all query sequences and their best hits. If the query sequences were supplied with tentative species identities, the distances are used to roughly estimate an ["accuracy category"](./results.md#accuracy-category).

### Create database tree
If enabled (checked or set to "true"), a multiple sequence alignment job will be run with ClustalOmega (hosted at EMBL-EBI) with default parameters. The submitted set of sequences is made from all the query sequences and all the sequences in the database.

## Job Name

This is an optional text field that can be added to describe the job for future reference. It will appear on results webpages and some of the result exports. This can be left blank.

