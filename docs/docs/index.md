# Welcome to Barrel 

Barrel is an upcoming web app designed to enable researchers to create, share and maintain barcode reference libraries collaboratively over the web. It consists of a backend application programming interface (API) as well as a modern and responsive frontend web interface.

**Barrel is still under active development, so documentation pages and features may not be fully complete. We are currently aiming for a release at the end of 2023.**

At a glance, Barrel allows users to:

- Mine sequence data from GenBank as reference barcode sequences.
- Organize barcode reference libraries by version and library. 
- Collaborate with other users to curate libraries with both public and private visibility.
- Filter mined sequences by GenBank search query, sequence quality, sequence length, taxonomy
- Run BLAST asynchronously on uploaded data and files to assign taxonomy to query sequences.
- Run multiple sequence alignment with database and query sequences included.
- Export databases to popular taxonomic assignment software.

Barrel thus provides an intuitive and user-friendly interface for users to collaborate to create specialized DNA barcoding reference libraries based on custom search queries and criteria. It also provides a platform for users to quickly share reference libraries and quickly make queries to determine taxonomic identities and reference library utility.

**Barrel can be configured to run locally, or as a web server.** When run as a web server, fellow researchers can see your publicly published reference libraries and run their sequences.

For more detailed API documentation, view our [interactive API documentation](/app/api-docs)