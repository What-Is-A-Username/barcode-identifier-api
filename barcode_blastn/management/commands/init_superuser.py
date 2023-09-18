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
		self.stdout.write(f'Verifying that superuser {options["username"]} exists with password {options["password"]}')
		User = get_user_model()
		if not User.objects.filter(username=options['username']).exists():
			self.stdout.write(f'Creating superuser {options["username"]} because no existing user found')
			User.objects.create_superuser(username=options['username'], email=options['email'], password=options['password'])
		else:
			ad = User.objects.get(username=options['username'])
			if not ad.is_superuser:
				self.stdout.write(f'\tWarning! {options["username"]} exists but is not a superuser')
			else:
				self.stdout.write(f'\tSuccess! {options["username"]} exists and is superuser')
			if not ad.is_staff:
				self.stdout.write(f'\tWarning! {options["username"]} exists but is not staff')
			else:
				self.stdout.write(f'\tSuccess! {options["username"]} exists and is staff')
		self.stdout.write(f'Superuser verified')
		return

