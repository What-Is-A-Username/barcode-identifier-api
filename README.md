# Barcode Identifier API

The Barcode Identifier API is a work-in-progress web API that allows users to BLAST against sequences in customized databases.

Built in Python with the Django Rest Framework.

Documentation will be made available when the project nears completion.

## First time installation

Clone from the repository:
```
git clone https://github.com/<username>/barcode-identifier-api.git
cd barcode-identifier-api/
```

Or, if the repo is private, follow the following format:
```
git clone https://username@github.com/<username>/barcode-identifier-api.git
cd barcode-identifier-api/
```

Install the required Linux packages through the terminal:
```
sudo apt install postgresql postgresql-contrib
sudo apt install libpq-dev
```

This package was built using Python 3.8.10. I recommend using the pyenv package to maintain multiple python versions installed on my machine. WARNING: use of pyenv will change how `python` is resolved in the terminal. You can check where `python` currently points to in the terminal with the command `which python`. The following uses pyenv to install and make a virtual environment for Python 3.8.10:
```
# check version
python --version

# recommend installing pyenv 
curl https://pyenv.run | bash
exec $SHELL

# install 3.8.10 using pyenv
pyenv install 3.8.10
pyenv global 3.8.10

# make the virtual env
python -m venv env_barcode
source env_barcode/bin/activate

# install all python packages required
pip install -r requirements.txt
``` 

Set up the PostGre SQL database:
```
sudo service postgresql start
sudo su - postgres
psql
CREATE DATABASE barcode_identifier_api
CREATE USER admin with PASSWORD '<password>'
GRANT ALL PRIVILEGES ON DATABASE barcode_identifier_api TO admin
\q
exit
```

Optionally, if you want the database to be populated with data from a local `db.sql` file:
```
psql -d barcode_identifier_db -f db.sql -U admin -h 127.0.0.1
```

Set-up the redis-server, and the default rqworker queue.
Ubuntu:
```
sudo apt-get install rabbitmq-server
sudo service rabbiqmq-server start
```

### Running the API
Ensure that you have already activated the virtual environment.

Open a terminal and run the worker queue used to perform the BLAST searches:
```
celery -A barcode_identifier_api worker --loglevel=INFO
```

Open another terminal and run the API server
```
python manage.py runserver
```

## Downloading database data
In the event that the database should be dumped/downloaded to a file, run the following in the terminal to create a `db.sql` file which can be transferred.
```
pg_dump -U admin -h 127.0.0.1 barcode_identifier_db > db.sql
```

## Troubleshooting

### I get the error `psql: error: FATAL:  Peer authentication failed for user "admin"` in the terminal when I try to run `psql`
If you are accessing `psql` through the terminal, try specifying the port number and the database to connect to. Example:
```
psql -d barcode_identifier_db -U admin -h 127.0.0.1
```

In the event that Postgre SQL still fails to authenticate you (e.g.: `psql: error: FATAL:  Peer authentication failed for user "admin"` or `psql: error: FATAL:  role "linux_username" does not exist`), one remedy is to change the authentication settings of Postgre SQL (requires root access) so that the `postgres` username is usable. Details on how to do so are [online](https://stackoverflow.com/a/26735105). After, running `psql` in the terminal requires you to specify the username explicitly with `-U postgres`. You can check whether you can log into the database by running `psql -U <username> <database>`.

### I get the error `django.db.utils.ProgrammingError: permission denied for table django_migrations` when running `python manage.py makemigrations`
This arises when the local PostGre SQL database is not owned by the the `admin` user which the API access the database from. We have to regrant permissions to the `admin` user, from some user second_user with higher permissions to the database (this second user is likely `postgres`).
```
# log into the database with a user with greater permissions
psql -U <second_user> -h 127.0.0.1 barcode_identifier_db

# in postgres:
GRANT ALL ON ALL TABLES IN SCHEMA public to admin;
GRANT ALL ON ALL SEQUENCES IN SCHEMA public to admin;
GRANT ALL ON ALL FUNCTIONS IN SCHEMA public to admin;
```
Now the error should be resolved, and commands like `python manage.py runserver` should work.

### I get errors when syncing the local database with the `db.sql` file with `psql -f db.sql`
One method of resolving it, albeit very blunt and destructive, is to remove all the data in the table and let the database be constructed from scratch using the `db.sql` file.
```
python manage.py flush
python manage.py migrate barcode_blastn zero
psql -d barcode_identifier_db -f db.sql -U admin -h 127.0.0.1 
```

Alternatively, you can drop the database and recreate it:
```
sudo su postgres
psql
drop database your_database_name;
create database your_database_name with owner user_you_use_in_django;
\q
exit
python manage.py makemigrations
python manage.py migrate
psql -h 127.0.0.1 -d barcode_identifier_db -U admin -f db.sql
```





