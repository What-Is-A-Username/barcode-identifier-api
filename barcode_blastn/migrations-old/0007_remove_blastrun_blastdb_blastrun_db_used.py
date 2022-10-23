# Generated by Django 4.1.1 on 2022-10-11 02:32

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('barcode_blastn', '0006_blastrun_blast_version_blastrun_errors_and_more'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='blastrun',
            name='blastdb',
        ),
        migrations.AddField(
            model_name='blastrun',
            name='db_used',
            field=models.ForeignKey(default=3, on_delete=django.db.models.deletion.CASCADE, related_name='usages', to='barcode_blastn.blastdb'),
            preserve_default=False,
        ),
    ]