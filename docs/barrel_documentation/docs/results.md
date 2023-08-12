## Reading Barrel Results

### Hits

The BLASTN run will return a set of hits for each query sequence supplied. 

<img src='../images/hits.png' width='80%' alt='Drawing'/>

This is displayed in a table with the following columns:

-   **Percent Identity, Alignment Length, Evalue, Bit Score**: As computed by BLAST
-   **Query ID**: The Seq ID line corresponding to which of the query sequences the hit is for.
-   **Organism, Country, Specimen Voucher, Type, Latitude/Longitude**: As mined from GenBank, corresponding to source organism, `/country`, `/specimen_voucher`, `/type_material` and `/lat_lon`.
-   **Info**: Click to open a comprehensive summary of the sequence, including its modification dates and publication information.

### Multiple Alignments and Trees

The number of alignments performed will depend the multiple [alignment options selected](#multiple-alignment-and-tree-parameters).

<img src='../images/tree.png' width='80%' alt='Interactive phylogenetic tree of hits and query sequences'/>

For *each* tree there are three downloadable files, each directly obtained from the job ran on ClustalOmega at EMBL-EBI and unaltered:

-   Multiple sequence alignment (MSA): A `.clustal_num` file containing  an alignment with nucleotides numbered.
-   Phylogenetic tree: `.ph` file containing the phylogenetic tree in Newick (Phylip) format. Note that the web app also allows you to copy the raw Newick string.
-   Sequences: A `.txt` file containing the sequences submitted to perform the multiple sequence alignment, in FASTA format.

In all three files, the sequence identifiers of query sequences will be altered by appending `|query` at the end to denote these as query sequences. All non-query sequences (i.e. sequences from the database, such as hits) will have sequence identifiers in the format of `<accession_number>|<organism_name>`, where organism name is generally `<genus>_<species>`.

More information and additional examples of these files are available in the [documentation for ClustalOmega at EMBL-EBI](https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation).

### Taxonomy Assignment

Barrel will attempt to assign taxonomic identity based on the top BLAST hit as determined by percent identity, and this assignment is considered the "Classification Result". If multiple hits have identical values, then this result is a comma-separated list of all unique organisms with the highest percent-identity.

If a species identity was initially given with the query sequences, then the parsed identity is visible under "Original Classification". 

#### Accuracy Category

- "Correct ID": Query < 1.0/2.0% divergent from reference, and query name matches reference species name
- "New ID": Query < 1.0/2.0% divergent from reference, and query not labelled to species, e.g. "Gymnotiformes" or "sp."
- "Incorrect ID": Query < 1.0/2.0% divergent from reference, and query name does not match reference species name
- "Tentative Correct ID : Query > 1.0/2.0% divergent from reference, and most similar reference name matches query species name
- "Unknown ID": Query > 1.0/2.0% divergent from reference, and query not labelled to species, e.g. "Gymnotiformes" or "sp."
- "Tentative Additional Species": Query > 1.0/2.0% divergent from reference, and query labelled as a species not included in reference library
- "Incorrect ID without Match": otherwise

<img src='../images/taxonomy.png' width='80%' alt='Interactive phylogenetic tree of hits and query sequences'/>
