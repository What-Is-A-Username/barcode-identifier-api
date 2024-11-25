# Barrel Changelog

## Version 0.0.2

### Feature Changes

-   Summary graphs now display properly when browsing a specific BLAST database.
-   Fixed incorrect static file serving for run results files (.clustal_num, .ph., .txt).
-   Fixed download links for downloading tree results.
-   Fixed spacing of 'Run A Query' buttons to fix unresponsive clicks.
-   Added version number to website footer.

## Version 0.0.3

### Feature Changes
-   BLAST run query submissions will ignore records already in the reference database; i.e. query sequences will never include reference sequences. 
-   Revamped server setup and installation process to allow user to specify between installations for private (e.g. only personal usage) and public (i.e. over the Internet) usage 
    -   Private usage will only make the app accessible on the localhost at user-specified ports (EXPOSED_PORT). This use case is meant for users who wish to install the tools for sole use on the same machine.
    -   Public usage will ensure that the app is only accessible over encrypted HTTPS (port 443) to allow user credentials and information to be securely communicated. This requires a SSL/TLS certificate, which must be obtained independently by the server admin and provided in the setup process.
        -   The app will listen on ports 80 and 443, so both must be free on the machine.
-   Removed example links from homepage, since link will differ between builds.

### Bug fixes
-   Fixed bug where database page crashed if there was a species with no taxonomy
-   Fixed bug where error message was not visible after submitting a BLAST run 
-   Fixed bug where pressing 'Run a Query' to the run submission page did not properly select the library and database.
-   Fixed bug where API could never successfully retrieve query identifiers file submitted in a run.

### Advanced changes (for developers and adminstrators)
-   Migrated website bundler from create-react-app to vite for the build pack for increased performance for both development and end-user.
-   Removed several unused dependencies.

## Known Issues
-   Library and database browse links on run submission page are broken.
-   Interactive API documentation tags for some /blastdb/ endpoints are erroneous.

## Version 0.0.4

### Bug Fixes
-   Cleaned unnecessary files from repo.

### Changes for Developers and Admins
-   Added helper scripts to push updates to docker.
-   Added build script in anticipation of future full release
-   Added additional docs for developers
-   Blast is now longer installed in the main barrel container
-   Docs now served as static files

## Version 0.0.5

### Feature Changes
-   Admin console now allows the uploading of custom sequences
    -   For this purpose, a new page is added for "Custom Sequences" which is separate from "GenBank Accessions". Custom Sequences refer to any reference library entries that were imported from file or manually added via the admin console
    -   Most fields of custom sequences can be edited. However, no checks for history or locked databases yet are implemented, so changes are not tracked
    -   The file used for uploading custom sequences must be a FASTA file.
        - The header can include the accession number, version, and definition in the format of `<accession>.<version> <definition>` (ie. accession and version cannot contain spaces and are period-separated, and the definition is the rest).
-   Sequences in the database now have a `data_source` attribute which indicate where the sequence originated (from GenBank or custom import/creation)
-   When viewing a blast database summary in the web app, number of sequences is now also tabulated by origin of sequence (in addition to previously available options such as journal title, taxa and country)

### Bug Fixes
-   Fixed links on Blast run page to correctly direct user to browse database and reference libraries when clicked.
-   

### Changes for Developers and Admins
-   Added new environment variable ENTREZ_EMAIL, which will be attached to every request to NCBI. This previously used hard-code or used the db admin as the fallback.
