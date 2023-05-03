from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand

class Command(BaseCommand):
	help = 'Create a superuser (admin-user) if none currently exists for the given username'

	def add_arguments(self, parser):
		parser.add_argument('--username', help='Username of superuser')
		parser.add_argument('--email', help='Email of superuser')
		parser.add_argument('--password', help='Password of superuser')
		return

	def handle(self, *args, **options):
		User = get_user_model()
		if not User.objects.filter(username=options['username']).exists():
			User.objects.create_superuser(username=options['username'], email=options['email'], password=options['password'])
		return

