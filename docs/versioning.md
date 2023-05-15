# Versioning

Continuous maintenance and updates are key to ensuring high-quality barcode reference libraries. Newly publish data and quality control findings frequently necessitate the need for corrections or updates to species identity and specimen metadata, and addition and removal of reference barcodes. 

To fulfill these needs, our application tracks database versions each with a meaningful version number.

## Version Number

Versions are assigned numbers in a meaningful format to make them human readable and useful: `<genbank_version>.<major version>.<minor version>`. 

**Genbank version** changes whenever an [accession.version](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/#VersionB) in the database changes, or a barcode is added or removed from the library. As described by GenBank, the accession.version is incremented whenever sequence data changes.

**Major version** changes when there is a change to the source information. Examples of database changes that necessitate a new major version include:

- changes in sequence information not indicated by the change in accession.version
- changes to source information of sample:
  - `organism`
  - `organelle`
  - `marker_gene`
  - `definition`
  - `country`
  - latitude/longitude (`lat_lon`)
  - specimen identifiers (`specimen_voucher`)
  - type information (`type_material`)
  - voucher information (`specimen_voucher`)
- type information changes 

**Minor version** changes when there is a change that will unlikely have the ability to change the results of future runs. These changes can include:
- description or name of the database is changed
- changes to specimen metadata:
  - `reference` (authorship, cited paper)
- specimen voucher or isolate is changed
- publication is changed

Versions do not change if:
- databases are made public or non-public by editors
- user permissions on databases are changed 

## Making updated versions

To construct a database, users with edit permissions in a reference library can create a new database under it. This database can inherit accession numbers from an existing database and can also have sequences added or deleted so long as the database is not 'locked'. The process of creating a new database will refetch data from GenBank.

After changes are made by the database maintainers, it must be published (locked) before it can be queried by users. Upon publishing, the application will:
- run `makeblastdb` with the sequence data 
- prevent users from making further changes
- assigns a version number based on the differences of this version from the most recent version

Once a database is locked, it should never be unlocked. This is since run results on that database relate back to the sequence data of the database, so its contents must remain present and unaltered for run results to persist.

Since the development copy is based on the latest version, versions can only follow one after another in a linear fashion. Maintainers who want to make a new version by branching or forking from an older version are recommended to create a separate reference library to do so.

## Using different versions

Unless specified otherwise, the API will retrieve the latest version of a database. Users can submit queries on older versions specifying the version number.