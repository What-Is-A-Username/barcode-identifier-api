import time
from typing import Any, Optional 

from psycopg2 import OperationalError as Psycopg2OperationalError

from django.db.utils import OperationalError
from django.core.management.base import BaseCommand

class Command(BaseCommand):
    """Django command to wait for the database to start."""

    def handle(self, *args: Any, **options: Any) -> Optional[str]:
        self.stdout.write('Waiting for database...')
        db_conn = None
        while not db_conn:
            try:
                self.check(databases=['default'])  # type: ignore
                db_conn = True
            except (Psycopg2OperationalError, OperationalError):
                self.stdout.write('Database unavailable, waiting 1 second...')
                time.sleep(1)

        self.stdout.write(self.style.SUCCESS('Database available!'))