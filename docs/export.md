# Exporting databases

Reference library databases can be exported to formats compatible with popular taxonomic assignment tools:
- [Exporting databases](#exporting-databases)
- [.fasta](#fasta)
- [.csv](#csv)
- [.tsv](#tsv)
- [DADA2](#dada2)
  - [assignTaxonomy(...)](#assigntaxonomy)
  - [assignSpecies(...)](#assignspecies)
- [sintax](#sintax)
- [RDP](#rdp)
  - [RDP training set sequence file](#rdp-training-set-sequence-file)
    - [RDP taxonomy file](#rdp-taxonomy-file)
- [Mothur](#mothur)
  - [Mothur `.pds.fasta` file](#mothur-pdsfasta-file)
  - [Mothur `.pds.tax` file](#mothur-pdstax-file)

# .fasta

With the most basic .fasta export, each sequence's header is just the ID.
```
>ON303390_1
ATAGTATTTGGTGCCTGAGCCGGAATAGTTGGCACAGCTTTAAGCCTCCTTATCCGAGCAGAACTAAGTC
```

# .csv

# .tsv

# DADA2

A reference library can be exported to be used as a training dataset for [dada2](https://benjjneb.github.io/dada2/index.html).

## assignTaxonomy(...) 

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

## assignSpecies(...)

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

# sintax

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

# RDP

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

## RDP training set sequence file

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

# Mothur

The reference library can be exported for use in the [`classify.seqs(...)` function](https://mothur.org/wiki/classify.seqs/). The files are exported as a zip file which contains a `.pds.fasta` file and a `.pds.tax` file. Simply unzip and specify these files in the `classify.seqs(...)` call.

Example:
```
mothur > classify.seqs(fasta=sequence.fasta,reference=my_database.pds.fasta,taxonomy=my_database.pds.tax)
```

## Mothur `.pds.fasta` file

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

## Mothur `.pds.tax` file

This is a taxonomy file comprised of two tab-separated columns, first for ID and second for taxonomic information. 

The taxonomic information is written as `Kingdom;Phylum;Class;Order;Family;Genus`.
```
ON303390_1  k__Metazoa;p__Chordata;c__Actinopterygii;o__Gymnotiformes;f__Gymnotidae;g__Gymnotus;
```

If there is missing information for a taxonomic rank, that rank and all ranks beneath it will be renamed systematically.
```
ON303390_1  k__Metazoa,p__Chordata;c__Actinopterygii;o__Gymnotiformes;f__o__Gymnotiformes;g__f__o__Gymnotiformes;
```




