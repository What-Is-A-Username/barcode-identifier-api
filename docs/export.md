# Exporting databases

Reference library databases can be exported to formats compatible with popular taxonomic assignment tools and other applications.

**Table of Contents**:
- [Exporting databases](#exporting-databases)
- [How to export a reference library](#how-to-export-a-reference-library)
- [Supported Export Formats](#supported-export-formats)
  - [Basic Exports](#basic-exports)
    - [.csv](#csv)
    - [.tsv](#tsv)
    - [.fasta](#fasta)
    - [.xml](#xml)
      - [Example:](#example)
    - [.json](#json)
      - [Example:](#example-1)
  - [QIIME 2](#qiime-2)
    - [QIIME2 `.fasta` file](#qiime2-fasta-file)
    - [QIIME2 `.txt` file](#qiime2-txt-file)
  - [DADA2](#dada2)
    - [DADA2 assignTaxonomy(...)](#dada2-assigntaxonomy)
    - [DADA2 assignSpecies(...)](#dada2-assignspecies)
  - [sintax](#sintax)
  - [RDP](#rdp)
    - [RDP training set sequence file](#rdp-training-set-sequence-file)
    - [RDP taxonomy file](#rdp-taxonomy-file)
  - [Mothur](#mothur)
    - [Mothur `.pds.fasta` file](#mothur-pdsfasta-file)
    - [Mothur `.pds.tax` file](#mothur-pdstax-file)

# How to export a reference library

The most direct way to export a reference library is to visit the API endpoint link corresponding to the database and export format desired. The link takes the form of `https://<domain>.com/blastdbs/<database-id>/export?export_format=<export_format>&format=<format>`, and it will allow you to download the requested files.

- `export_format` is used to indicate the format requested
- `format` is used to indicate the file type requested (can be `zip` or `fasta`, depending on the export format requested)

*Example*: To export a database with ID `66855f2c-f360-4ad9-8c98-998ecb815ff5` from `https://barcodebarrel.com` in order to use it for the classifier offered in Mothur, you would visit:
```
https://barcodebarrel.com/blastdbs/66855f2c-f360-4ad9-8c98-998ecb815ff5/export?export_format=mothur&format=zip
```

| Supported Export Format                                | Summary                                                            | `export_format` | `format`           |
| ------------------------------------------------------ | ------------------------------------------------------------------ | ------------- | ---------------- |
| [Basic Export](#)                                      | File containing sequence and/or database metadata. | Blank | `csv`, `tsv`, `fasta`, `json` or `xml` |
| [QIIME2](#qiime-2)                                     | Zip file containing sequence and taxonomy files for use by QIIME2. | `qiime2`      | `zip`            |
| [DADA2 (`assignTaxonomy(...)`)](#dada2-assigntaxonomy) | Zip or .fasta file containing sequence data for use by DADA2.      | `dada2tax`    | `zip` or `fasta` |
| [DADA2 (`assignSpecies(...)`)](#dada2-assignspecies)   | Zip or .fasta file containing sequence data for use by DADA2.      | `dada2sp`     | `zip` or `fasta` |
| [SINTAX](#sintax)                                      | Zip or .fasta file containing sequence data for use by SINTAX.     | `sintax`      | `fasta` or `zip` |
| [RDP](#rdp)                                            | Zip file containing sequence and taxonomy files for use by RDP.    | `rdp`         | `zip`            |
| [Mothur](#mothur)                                      | Zip file containing sequence and taxonomy files for use by Mothur. | `mothur`      | `zip`            |

# Supported Export Formats

## Basic Exports

The database can be output to several quintessential file formats: `.csv`, `.tsv`, `.fasta`, `.xml` and `.json`.

### .csv

The database can be exported to a comma-separated values file (.csv) which contains only the sequence data. Other database data, like user data, names and descriptions, is not exported. Info for the following columns is provided for each sequence:
- `accession_number`
- `version` (accession.version)
- `organism`
- `organelle`
- `isolate` 
- `country`
- `specimen_voucher`
- `dna_sequence`
- `lat_lon`
- `type_material`
- `created`
- `updated`
- `genbank_modification_date`
- `taxonomy`
- `taxon_superkingdom_id`
- `taxon_superkingdom_scientific_name`
- `taxon_kingdom_id`
- `taxon_kingdom_scientific_name`
- `taxon_phylum_id`
- `taxon_phylum_scientific_name`
- `taxon_class_id`
- `taxon_class_scientific_name`
- `taxon_order_id`
- `taxon_order_scientific_name`
- `taxon_family_id`
- `taxon_family_scientific_name`
- `taxon_genus_id`
- `taxon_genus_scientific_name`
- `taxon_species_id`
- `taxon_species_scientific_name`
- `authors`
- `title`
- `journal`

### .tsv

The tab-separated values file (.tsv) exported contains sequence data using the same columns as the [.csv export](#csv). Other database data, like user data, names and descriptions, is not included.

### .fasta

The .fasta file exported contains only the sequence and identifier data. Each sequence's header is just the ID.
```
>ON303390_1
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

### .xml

The .xml file exported contains information about each sequence as well as metadata about the database and reference library themselves. Thus, it is the most comprehensive export format alongside [.json](#json)

#### Example:
```
<?xml version="1.0" encoding="UTF-8"?>
<Database>
	<id>e3e61e46-3aaa-42d2-bcb6-89c7fd7eddb9</id>
	<library>
		<id>1f739dde-ca5e-4744-a6cb-8d79cc1db907</id>
		<custom_name>c</custom_name>
		<description>d</description>
		<public>false</public>
		<owner>
			<username>admin</username>
			<email>admin@example.com</email>
		</owner>
	</library>
	<custom_name>with superkingdom</custom_name>
	<version_number>1.1.1</version_number>
	<description/>
	<locked>true</locked>
	<sequences>
		<accession_number>ON303324</accession_number>
		<version>ON303324.1</version>
		<organism>Adontosternarchus balaenops</organism>
		<organelle>mitochondrion</organelle>
		<isolate>2612</isolate>
		<country>Peru</country>
		<specimen_voucher>UF:116559</specimen_voucher>
		<dna_sequence>ATAGTATTTGGCGCCTGAGCCGGTATAATTGGAACTGCTCTCAGCCTATTAATTCGAGCTGAGCTCAACCAACCCGGCACCCTCTTAGAAGATGACCAAATTTACAACGTGGCAGTCACCGCTCATGCCTTCGTAATAATTTTCTTTATAGTTATGCCCATCATAATTGGAGGCTTTGGCAATTGACTTATTCCTTTAATAATTGCCGCACCAGACATGGCATTCCCCCGAATAAATAACATAAGCTTCTGACTACTCCCCCCATCATTCTTCCTTCTCCTTGCCTCCGCCGGCTTAGAAGCTGGGGTTGGAACAGGCTGAACCCTATACCCCCCTCTTGCTGGCAACGCCGCACACGCCGGAGCTTCTGTTGACTTAACCATTTTCTCCCTTCACCTTGCCGGTGTCTCCTCCATTCTCGGCTCCATCAACTTTATTACCACAATTATTAATATAAAACCCCCCACAATAACTCAATACCAACTTCCATTATTCATTTGATCCCTACTAGTAACTACCGTACTCCTACTACTTTCTCTCCCTGTCCTAGCTGCTGGCATTACCATACTTCTTACAGACCGAAATTTAAACACAGCATTCTTTGATCCTACAGGAGGAGGAGACCCCATCCTGTACCAACACCTA</dna_sequence>
		<lat_lon>3.74 S 73.25 W</lat_lon>
		<type_material/>
		<created>2023-06-02T17:35:09.262135Z</created>
		<updated>2023-06-02T17:35:09.262139Z</updated>
		<genbank_modification_date>2022-07-04</genbank_modification_date>
		<taxonomy>Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Sternopygoidei,Apteronotidae,Adontosternarchus</taxonomy>
		<taxon_superkingdom>
			<id>2759</id>
			<scientific_name>Eukaryota</scientific_name>
		</taxon_superkingdom>
		<taxon_kingdom>
			<id>33208</id>
			<scientific_name>Metazoa</scientific_name>
		</taxon_kingdom>
		<taxon_phylum>
			<id>7711</id>
			<scientific_name>Chordata</scientific_name>
		</taxon_phylum>
		<taxon_class>
			<id>186623</id>
			<scientific_name>Actinopteri</scientific_name>
		</taxon_class>
		<taxon_order>
			<id>8002</id>
			<scientific_name>Gymnotiformes</scientific_name>
		</taxon_order>
		<taxon_family>
			<id>30766</id>
			<scientific_name>Apteronotidae</scientific_name>
		</taxon_family>
		<taxon_genus>
			<id>36683</id>
			<scientific_name>Adontosternarchus</scientific_name>
		</taxon_genus>
		<taxon_species>
			<id>1740075</id>
			<scientific_name>Adontosternarchus balaenops</scientific_name>
		</taxon_species>
		<authors>Janzen,F.H., Crampton,W.G. and Lovejoy,N.R.</authors>
		<title>A new taxonomist-curated reference library of DNA barcodes for Neotropical electric fishes (Teleostei: Gymnotiformes)</title>
		<journal>Zool J Linn Soc (2022) In press</journal>
	</sequences>
	<sequences>
		<accession_number>ON303325</accession_number>
		<version>ON303325.1</version>
		<organism>Adontosternarchus clarkae</organism>
		<organelle>mitochondrion</organelle>
		<isolate>2906</isolate>
		<country>Brazil</country>
		<specimen_voucher>MCP 39341</specimen_voucher>
		<dna_sequence>ATAGTGTTTGGTGCCTGAGCCGGCATAATTGGAACTGCTCTCAGCCTGTTGATTCGGGCAGAGCTCAACCAACCTGGCACTCTCTTAGAAGACGACCAAATTTATAACGTAGCCGTTACCGCTCATGCCTTCGTAATAATTTTCTTTATAGTTATGCCAATCATAATTGGAGGCTTTGGCAATTGACTTATCCCCCTAATAATTGCCGCCCCAGACATGGCATTCCCACGAATAAATAACATAAGCTTTTGACTACTTCCCCCATCATTTTTCCTCCTTCTCGCCTCTGCTGGCTTAGAAGCTGGAGTAGGAACAGGCTGAACCTTATACCCCCCTCTTGCTGGCAACGCCGCACACGCCGGAGCTTCCGTAGACCTAACCATTTTCTCTCTCCACCTTGCCGGTGTCTCCTCCATCCTTGGCTCTATCAACTTTATTACTACCATCATTAATATGAAACCCCCCACAATAACCCAATATCAACTCCCACTATTTATCTGATCCCTGCTAGTAACCACCGTACTCCTTCTCCTTTCTCTCCCCGTTCTAGCTGCCGGCATTACCATACTTCTTACGGACCGAAATTTAAACACAGCATTCTTTGACCCCACAGGAGGAGGAGACCCCATCCTGTACCAACACTTA</dna_sequence>
		<lat_lon>3.12 S 64.78 W</lat_lon>
		<type_material/>
		<created>2023-06-02T17:35:09.260837Z</created>
		<updated>2023-06-02T17:35:09.260841Z</updated>
		<genbank_modification_date>2022-07-04</genbank_modification_date>
		<taxonomy>Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Sternopygoidei,Apteronotidae,Adontosternarchus</taxonomy>
		<taxon_superkingdom>
			<id>2759</id>
			<scientific_name>Eukaryota</scientific_name>
		</taxon_superkingdom>
		<taxon_kingdom>
			<id>33208</id>
			<scientific_name>Metazoa</scientific_name>
		</taxon_kingdom>
		<taxon_phylum>
			<id>7711</id>
			<scientific_name>Chordata</scientific_name>
		</taxon_phylum>
		<taxon_class>
			<id>186623</id>
			<scientific_name>Actinopteri</scientific_name>
		</taxon_class>
		<taxon_order>
			<id>8002</id>
			<scientific_name>Gymnotiformes</scientific_name>
		</taxon_order>
		<taxon_family>
			<id>30766</id>
			<scientific_name>Apteronotidae</scientific_name>
		</taxon_family>
		<taxon_genus>
			<id>36683</id>
			<scientific_name>Adontosternarchus</scientific_name>
		</taxon_genus>
		<taxon_species>
			<id>597343</id>
			<scientific_name>Adontosternarchus clarkae</scientific_name>
		</taxon_species>
		<authors>Janzen,F.H., Crampton,W.G. and Lovejoy,N.R.</authors>
		<title>A new taxonomist-curated reference library of DNA barcodes for Neotropical electric fishes (Teleostei: Gymnotiformes)</title>
		<journal>Zool J Linn Soc (2022) In press</journal>
	</sequences>
</Database>
```

### .json

Contains the same information as the [.xml](#xml) export, but in Javascript Object Notation (JSON) format.

#### Example:
```
{
  "id": "e3e61e46-3aaa-42d2-bcb6-89c7fd7eddb9",
  "library": {
    "id": "1f739dde-ca5e-4744-a6cb-8d79cc1db907",
    "custom_name": "c",
    "description": "d",
    "public": false,
    "owner": {
      "username": "admin",
      "email": "admin@example.com"
    }
  },
  "custom_name": "with superkingdom",
  "version_number": "1.1.1",
  "description": "",
  "locked": true,
  "sequences": [
    {
      "accession_number": "ON303324",
      "version": "ON303324.1",
      "organism": "Adontosternarchus balaenops",
      "organelle": "mitochondrion",
      "isolate": "2612",
      "country": "Peru",
      "specimen_voucher": "UF:116559",
      "dna_sequence": "ATAGTATTTGGCGCCTGAGCCGGTATAATTGGAACTGCTCTCAGCCTATTAATTCGAGCTGAGCTCAACCAACCCGGCACCCTCTTAGAAGATGACCAAATTTACAACGTGGCAGTCACCGCTCATGCCTTCGTAATAATTTTCTTTATAGTTATGCCCATCATAATTGGAGGCTTTGGCAATTGACTTATTCCTTTAATAATTGCCGCACCAGACATGGCATTCCCCCGAATAAATAACATAAGCTTCTGACTACTCCCCCCATCATTCTTCCTTCTCCTTGCCTCCGCCGGCTTAGAAGCTGGGGTTGGAACAGGCTGAACCCTATACCCCCCTCTTGCTGGCAACGCCGCACACGCCGGAGCTTCTGTTGACTTAACCATTTTCTCCCTTCACCTTGCCGGTGTCTCCTCCATTCTCGGCTCCATCAACTTTATTACCACAATTATTAATATAAAACCCCCCACAATAACTCAATACCAACTTCCATTATTCATTTGATCCCTACTAGTAACTACCGTACTCCTACTACTTTCTCTCCCTGTCCTAGCTGCTGGCATTACCATACTTCTTACAGACCGAAATTTAAACACAGCATTCTTTGATCCTACAGGAGGAGGAGACCCCATCCTGTACCAACACCTA",
      "lat_lon": "3.74 S 73.25 W",
      "type_material": "",
      "created": "2023-06-02T17:35:09.262135Z",
      "updated": "2023-06-02T17:35:09.262139Z",
      "genbank_modification_date": "2022-07-04",
      "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Sternopygoidei,Apteronotidae,Adontosternarchus",
      "taxon_superkingdom": {
        "id": 2759,
        "scientific_name": "Eukaryota"
      },
      "taxon_kingdom": {
        "id": 33208,
        "scientific_name": "Metazoa"
      },
      "taxon_phylum": {
        "id": 7711,
        "scientific_name": "Chordata"
      },
      "taxon_class": {
        "id": 186623,
        "scientific_name": "Actinopteri"
      },
      "taxon_order": {
        "id": 8002,
        "scientific_name": "Gymnotiformes"
      },
      "taxon_family": {
        "id": 30766,
        "scientific_name": "Apteronotidae"
      },
      "taxon_genus": {
        "id": 36683,
        "scientific_name": "Adontosternarchus"
      },
      "taxon_species": {
        "id": 1740075,
        "scientific_name": "Adontosternarchus balaenops"
      },
      "authors": "Janzen,F.H., Crampton,W.G. and Lovejoy,N.R.",
      "title": "A new taxonomist-curated reference library of DNA barcodes for Neotropical electric fishes (Teleostei: Gymnotiformes)",
      "journal": "Zool J Linn Soc (2022) In press"
    },
    {
      "accession_number": "ON303325",
      "version": "ON303325.1",
      "organism": "Adontosternarchus clarkae",
      "organelle": "mitochondrion",
      "isolate": "2906",
      "country": "Brazil",
      "specimen_voucher": "MCP 39341",
      "dna_sequence": "ATAGTGTTTGGTGCCTGAGCCGGCATAATTGGAACTGCTCTCAGCCTGTTGATTCGGGCAGAGCTCAACCAACCTGGCACTCTCTTAGAAGACGACCAAATTTATAACGTAGCCGTTACCGCTCATGCCTTCGTAATAATTTTCTTTATAGTTATGCCAATCATAATTGGAGGCTTTGGCAATTGACTTATCCCCCTAATAATTGCCGCCCCAGACATGGCATTCCCACGAATAAATAACATAAGCTTTTGACTACTTCCCCCATCATTTTTCCTCCTTCTCGCCTCTGCTGGCTTAGAAGCTGGAGTAGGAACAGGCTGAACCTTATACCCCCCTCTTGCTGGCAACGCCGCACACGCCGGAGCTTCCGTAGACCTAACCATTTTCTCTCTCCACCTTGCCGGTGTCTCCTCCATCCTTGGCTCTATCAACTTTATTACTACCATCATTAATATGAAACCCCCCACAATAACCCAATATCAACTCCCACTATTTATCTGATCCCTGCTAGTAACCACCGTACTCCTTCTCCTTTCTCTCCCCGTTCTAGCTGCCGGCATTACCATACTTCTTACGGACCGAAATTTAAACACAGCATTCTTTGACCCCACAGGAGGAGGAGACCCCATCCTGTACCAACACTTA",
      "lat_lon": "3.12 S 64.78 W",
      "type_material": "",
      "created": "2023-06-02T17:35:09.260837Z",
      "updated": "2023-06-02T17:35:09.260841Z",
      "genbank_modification_date": "2022-07-04",
      "taxonomy": "Eukaryota,Metazoa,Chordata,Craniata,Vertebrata,Euteleostomi,Actinopterygii,Neopterygii,Teleostei,Ostariophysi,Gymnotiformes,Sternopygoidei,Apteronotidae,Adontosternarchus",
      "taxon_superkingdom": {
        "id": 2759,
        "scientific_name": "Eukaryota"
      },
      "taxon_kingdom": {
        "id": 33208,
        "scientific_name": "Metazoa"
      },
      "taxon_phylum": {
        "id": 7711,
        "scientific_name": "Chordata"
      },
      "taxon_class": {
        "id": 186623,
        "scientific_name": "Actinopteri"
      },
      "taxon_order": {
        "id": 8002,
        "scientific_name": "Gymnotiformes"
      },
      "taxon_family": {
        "id": 30766,
        "scientific_name": "Apteronotidae"
      },
      "taxon_genus": {
        "id": 36683,
        "scientific_name": "Adontosternarchus"
      },
      "taxon_species": {
        "id": 597343,
        "scientific_name": "Adontosternarchus clarkae"
      },
      "authors": "Janzen,F.H., Crampton,W.G. and Lovejoy,N.R.",
      "title": "A new taxonomist-curated reference library of DNA barcodes for Neotropical electric fishes (Teleostei: Gymnotiformes)",
      "journal": "Zool J Linn Soc (2022) In press"
    },
  ]
}
```


## QIIME 2

The reference library can be exported to QIIME 2 for use by the [q2-feature-classifier](https://library.qiime2.org/plugins/q2-feature-classifier/3/) plugin used for taxonomic classification. 

The exported file is a zip file containing two files, a `.fasta` file containing the sequences and a `.txt` file containing the taxonomic information.

There is a [tutorial in the QIIME2 documentation](https://docs.qiime2.org/2023.5/tutorials/feature-classifier/) that shows how to use q2-feature-classifier.

### QIIME2 `.fasta` file 
Each sequence is exported as a single line, with the header containing only the accession.
```
>ON303390_1
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

### QIIME2 `.txt` file
The taxonomy file contains one line per sequence. Each line consists of two tab-separated columns, the first column containing the ID and the second column containing the lineage in the form of `k__Kingdom; p__Phylum; c__Class; o__Order; f__Family; g__Genus; s__species`.
```
ON303390_1 k__Metazoa; p__Chordata; c__Actinopterygii; o__Gymnotiformes; f__Gymnotidae; g__Gymnotus; s__cylindricus
```

If a taxonomic rank is missing, that rank and all the ranks below it are printed without a specified name:
```
ON303390_1 k__Metazoa; p__Chordata; c__Actinopterygii; o__Gymnotiformes; f__; g__; s__
```

## DADA2

A reference library can be exported to be used as a training dataset for [dada2](https://benjjneb.github.io/dada2/index.html), used for `assignTaxonomy(...)` or `assignSpecies(...)`.

### DADA2 assignTaxonomy(...) 

This format is for using the `assignTaxonomy(...)` function within DADA2. The header of each sequence is a semicolon-separated list of Kingdom, Phylum, Class, Order, Family, Genus:
```
>Metazoa;Chordata;Actinopterygii;Gymnotiformes;Gymnotidae;Gymnotus;
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```
For sequences with incomplete information, the list stops at the first missing rank:
```
>Metazoa;Chordata;Actinopterygii;Gymnotiformes;
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

The resulting .fasta file exported can then be used when running DADA2 in R:
```
library(dada2)
# load our query sequences in from the fasta file
fs = dada2:::derepFasta('sequence.fasta')
seqTable = dada2::makeSequenceTable(fs)

set.seed(100)
# run classification procedure by passing in query table and exported fasta file
taxonomic_assignments <- assignTaxonomy(seqTable, "my_database.fasta", multithread=TRUE)
```

### DADA2 assignSpecies(...)

This format is for using the `assignSpecies(...)` function within DADA2.

The header of each sequence takes the form of `ID Genus Species`:
```
>ON303390_1 Gymnotus cylindricus
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

The resulting .fasta file exported can then be used when running DADA2 in R:
```
library(dada2)
# load our query sequences in from the fasta file
fs = dada2:::derepFasta('sequence.fasta')
seqTable = dada2::makeSequenceTable(fs)

set.seed(100)
# run classification procedure by passing in query table and exported fasta file
species_assignments <- dada2::assignSpecies(seqTable, "dada2species.fasta")
```

## sintax

The reference database is exported as a single .fasta file that is SINTAX compatible SINTAX algorithm, so it can be run with software such as [USEARCH](https://www.drive5.com/usearch/manual/cmd_sintax.html).

Example with USEARCH:
```
# Reference database was exported as my_database.fasta.
# Make the database file in USEARCH so that it can be loaded quickly:
./usearch11.exe -makeudb_usearch my_database.fasta -output my_database.udb
# Run the SINTAX classification on sequences within query.fasta:
./usearch11.exe -sintax query.fasta -db my_database.udb -tabbedout reads.sintax -strand both -sintax_cutoff 0.8
```

The header of each sequence takes the form of `ID;tax=k:Kingdom,p:Phylum,c:Class,o:Order,f:Family,g:Genus`.
```
>ON303390_1;tax=k:Metazoa,p:Chordata,c:Actinopterygii,o:Gymnotiformes,f:Gymnotidae,g:Gymnotus;
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

If information for taxonomic ranks are missing for a sequence, only those ranks are omitted from that particular header.
```
>ON303390_1;tax=k:Metazoa,p:Chordata,o:Gymnotiformes,g:Gymnotus;
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

More information on using SINTAX within USEARCH is [here](https://www.drive5.com/usearch/manual/tax_annot.html).

## RDP

This export format will produce a zip file which contains both the sequence file and the taxonomy file required for a new training set in RDP based on the reference library.

The format is consistent with that produced by [mkCOInr](https://mkcoinr.readthedocs.io/en/latest/content/io.html#rdp-classifier-taxid-file).

Example usage:
```
# Suppose we extracted the zip file contents to /my_database/
# First create the training files
java -jar classifier.jar train -o output_training_files -s /my_database.fasta -t my_database/my_database.txt
# Run the classifier on our query sequences in query.fasta, 
# writing results to classifications.txt
java -jar classifier.jar classify -c 0.8 -f allrank -t output_training_files/rRNAClassifier.properties -o classifications.txt query.fasta
```

### RDP training set sequence file

This is a fasta file, with the header of each sequence of the form `ID   sk__Superkingdom;k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus`. 

```
>ON303390_1 sk__Eukaryota__;k__Metazoa;p__Chordata;c__Actinopterygii;o__Gymnotiformes;f__Gymnotidae;g__Gymnotus;
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

If information for taxonomic ranks are missing for a sequence, all ranks under it are renamed according to the last rank specified. The species rank will be made unique using the identifier.
```
>ON303390_1 sk__Eukaryota__;k__Metazoa;p__Chordata;c__Actinopterygii;o__Gymnotiformes;f__o__Gymnotiformes;g__f__o__Gymnotiformes;g__f__o__Gymnotiformes__ON303390_1
``` 

Note that the ID is always separated from the taxonomy data by a tab, and not spaces.

### RDP taxonomy file

This is a file describing the taxonomy with each line containing a sequence of values separated by asterisks (*) in the form of `TaxID*Taxon_name_TaxID*Parent_TaxID*Taxonomic_rank_index*Taxonomic_rank`.

The `Taxonomic_rank` is one of root, superkingdom, kingdom, phylum, class, order, family, genus, species. The `Taxonomic_rank_index` is the number corresponding to each `Taxonomic_rank`, with root = 0, superkingdom = 1 ... species = 8.

```
699532*s__Gymnotus cylindricus*36670*8*species
36670*g__Gymnotus*30771*7*genus
30771*f__Gymnotidae*1489620*6*family
...
33208*k__Metazoa*2759*3*kingdom
0*Root*-1*0*rootrank
```

## Mothur

The reference library can be exported for use in the [`classify.seqs(...)` function](https://mothur.org/wiki/classify.seqs/). The files are exported as a zip file which contains a `.pds.fasta` file and a `.pds.tax` file. Simply unzip and specify these files in the `classify.seqs(...)` call.

Example:
```
mothur > classify.seqs(fasta=sequence.fasta,reference=my_database.pds.fasta,taxonomy=my_database.pds.tax)
```

### Mothur `.pds.fasta` file

This is a fasta file, with the header of each sequence of the form `ID   Root;sk__Superkingdom;k__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species`. 

```
>ON303390_1 Root;sk__Eukaryota;k__Metazoa;p__Chordata;c__Actinopterygii;o__Gymnotiformes;f__Gymnotidae;g__Gymnotus;s__Gymnotus_cylindricus;
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

If information for taxonomic ranks are missing for a sequence, that rank and all ranks below it are omitted.
```
>ON303390_1 Root;sk__Eukaryota;k__Metazoa;p__Chordata;c__Actinopterygii;o__Gymnotiformes;
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
``` 

Note that the ID is always separated from the taxonomy data by a tab character, and not spaces.

### Mothur `.pds.tax` file

This is a taxonomy file comprised of two tab-separated columns, first for ID and second for taxonomic information. 

The taxonomic information is written as `Kingdom;Phylum;Class;Order;Family;Genus`.
```
ON303390_1  k__Metazoa;p__Chordata;c__Actinopterygii;o__Gymnotiformes;f__Gymnotidae;g__Gymnotus;
```

If there is missing information for a taxonomic rank, that rank and all ranks beneath it will be renamed systematically.
```
ON303390_1  k__Metazoa,p__Chordata;c__Actinopterygii;o__Gymnotiformes;f__o__Gymnotiformes;g__f__o__Gymnotiformes;
```




