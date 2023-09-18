# Permissions

The BLAST database allows for different levels of permissions to be assigned to site users, which will be summarized here.

The [Permission Table](#permission-table) summarizes what actions a user can perform based on user and database settings. For explanations of the terminology used, follow the links or directly read from the [Definitions](#definitions) section.

## Permission Table

<table>
	<tbody text-align=center align=center>
		<tr>
			<td rowspan=3><a href="#public-vs-private">BLAST Database Visibility</a></td>
			<td rowspan=3><a href="#database-actions">Action</a></td>
			<td rowspan=3><a href="#users-and-accounts">Public user</a></td>
			<td colspan=6><a href="#users-and-accounts">Authenticated user</a></td>
		</tr>
		<tr>
			<td rowspan=2>No additional permissions given</td>
			<td colspan=3><a href="#permission-levels">User Permission Given</a></td>
			<td rowspan=2><a href="#ownership">Owner</a></td>
			<td rowspan=2><a href="#users-and-accounts">Superuser</a></td>
		</tr>
		<tr>
			<td>Can View Db</td>
			<td>Can Run Db</td>
			<td>Can Edit Db</td>
		</tr>
		<tr>
			<td>Public</td>
			<td>View a BLAST database</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
		</tr>
		<tr>
			<td>Public</td>
			<td>Running BLAST on an existing BLAST database</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
		</tr>
		<tr>
			<td>Private or Public</td>
			<td>Editing existing BLAST database</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
		</tr>
		<tr>
			<td>Private or Public</td>
			<td>Granting user permissions to an existing BLAST database</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>✔️</td>
			<td>✔️</td>
		</tr>
		<tr>
			<td>Private or Public</td>
			<td>Deleting existing BLAST database</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>✔️</td>
		</tr>
		<tr>
			<td>Private or Public</td>
			<td>Create a new BLAST database</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>❌</td>
			<td>✔️</td>
		</tr>
		<tr>
			<td>Private</td>
			<td>View a BLAST database</td>
			<td>❌</td>
			<td>❌</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
		</tr>
		<tr>
			<td>Private</td>
			<td>Running BLAST on an existing BLAST database</td>
			<td>❌</td>
			<td>❌</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
			<td>✔️</td>
		</tr>
	</tbody>
</table>

## Definitions

### Users and Accounts

A user can be either **authenticated** or a **public**. An authenticated user refers to any user that is signed in into an account. Otherwise, they are public user. 

Authenticated users can additionally be designated as **staff** and/or **superusers**. These can only be designated by superusers.

### BLAST database-level settings

#### Ownership
Every BLAST database has an **owner**, which refers to an authenticated user that first created the BLAST database instance.

#### Public vs. private
People with edit permissions of a BLAST database can designate it as **public**, allowing public users to view and run BLAST with the database. Otherwise, the BLAST database is considered **private**.

#### Permission Levels
The database owner and any superusers can assign different levels of **database permissions** to authenticated users:

-   `Can Edit Db`
-   `Can Run Db`
-   `Can View Db`
-   `Deny Access`

### Database Actions

The user details and assigned database permission level will be used to determine what actions a given user can perform.

#### Viewing an existing BLAST database

**This action allows a user to**

- view associated database information (sequences, database name, database description, version, etc.)

**This action requires the user to meet any ONE of the following conditions**

-   [user can run the database](#running-blast-on-an-existing-blast-database)
-   user is an authenticated user and has been given any of the following database permission levels:
	- `Can Edit Db`
	- `Can Run Db`
	- `Can View Db`

#### Running BLAST on an existing BLAST database

**This action allows a user to**

- submit a BLAST query on the database

**This action requires the user to meet any ONE of the following conditions**

-   [the user can edit the database](#editing-an-existing-blast-database)
-   user is public and the database is public
-   user is authenticated and possesses any ONE of following database permission levels:
	- `Can Edit Db`
	- `Can Run Db`

#### Editing an existing BLAST database

**This action allows a user to**

- edit associated database information (database name, database description, version, etc.)
- add, edit and remove accession numbers from the database

**This action requires the user to meet any ONE of the following conditions**

-   [the user can delete the database](#deleting-an-existing-blast-database)
-   the user is authenticated and is the database owner
-   the user is authenticated and possesses any ONE of the following database permission levels:
	- `Can Edit Db`

#### Granting user permissions to a BLAST database

**This action allows a user to**

- add, edit, and delete a [permission level](#permission-levels) for another user

**This action requires the user to meet any ONE of the following conditions**

-   [the user can delete the database](#deleting-an-existing-blast-database)
-   the user is authenticated and is the database owner

#### Deleting an existing BLAST database

**This action allows a user to**
- delete the BLAST database (also deletes associated sequences, runs)

**This action requires the user to be a superuser**

#### Creating a new BLAST database
**This action allows a user to**
- create a BLAST database instance (also deletes associated sequences, runs) and be designated the BLAST database owner

*Tip: Although a superuser is needed to create the initial BLAST database instance, they can grant a non-superuser user edit permissions so that this other user can add/edit/remove sequences to maintain the BLAST database.*

**This action requires the user to be a superuser**

