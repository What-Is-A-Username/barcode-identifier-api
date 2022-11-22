# Generated by Django 4.1.2 on 2022-11-21 00:02

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('barcode_blastn', '0006_alter_hit_db_entry'),
    ]

    operations = [
        migrations.AddField(
            model_name='blastdb',
            name='description',
            field=models.CharField(blank=True, default='', max_length=1024),
        ),
        migrations.AddField(
            model_name='nuccoresequence',
            name='lat_lon',
            field=models.CharField(blank=True, default='', max_length=64),
        ),
    ]