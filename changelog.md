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

### Feature changes
-   Added documentation pages to website, using Read The Docs styles. Raw markdown content is found in the `docs` folder of this repository
-   Added accuracy classification between supposed species identity and BLAST result species, using Kimura-2-parameter genetic distance between sequences.

### Bug Fixes
-   Cleaned unnecessary files from repo.

### Changes for Developers and Admins
-   Added helper scripts to push updates to docker.
-   Added scripts for server deployment, but these scripts will be moved to another repo for full release.