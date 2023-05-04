# Permissions

The BLAST database allows for different levels of permissions to be assigned to site users, which will be summarized here.

## Users and Accounts

A user can be either **authenticated** or a **public**. An authenticated user refers to any user that is signed in into an account. Otherwise, they are public user. 

Authenticated users can additionally be designated as **staff** and/or **superusers**. These can only be designated by superusers.

## Database-level settings
Every BLAST database has an **owner**, which refers to an authenticated user that first created the BLAST database instance.

### Permissions Levels
The database owner can assign different levels of **database permissions** to authenticated users:
-   `Can Edit Db`
-   `Can Run Db`
-   `Can View Db`
-   `Deny Access`

## Database Actions
The user details and assigned database permission level will be used to determine what actions a given user can perform.
-   [viewing an existing BLAST database](#viewing-an-existing-blast-database)
-   [running BLAST on an existing BLAST database](#running-blast-on-an-existing-blast-database)
-   [editing an existing BLAST database](#editing-an-existing-blast-database)
-   [deleting an existing BLAST database](#deleting-an-existing-blast-database)
-   [creating a new BLAST database](#creating-a-new-blast-database)

### Viewing an existing BLAST database
Requires any one of the following conditions to be met:
-   [user can run the database](#running-blast-on-an-existing-blast-database)
-   user is an authenticated user and has been given any of the following database permission levels:
    - `Can Edit Db`
    - `Can Run Db`
    - `Can View Db`

### Running BLAST on an existing BLAST database
Requires any ONE of the following conditions to be met:
-   [the user can edit the database](#editing-an-existing-blast-database)
-   user is public and the database is public
-   user is authenticated and possesses any ONE of following database permission levels:
    - `Can Edit Db`
    - `Can Run Db`

### Editing an existing BLAST database
Requires any ONE of the following conditions to be met:
-   [the user can delete the database](#deleting-an-existing-blast-database)
-   the user is authenticated and possesses any ONE of the following database permission levels:
    - `Can Edit Db`

### Deleting an existing BLAST database
Requires that the user is a superuser.

### Creating a new BLAST database
Requires that the user is a superuser.
