from django.apps import AppConfig

class BarcodeBlastnConfig(AppConfig):
    name = 'barcode_blastn'

    def ready(self) -> None:
        print('Config for `barcode_blastn` is ready.')
        import barcode_blastn.signals