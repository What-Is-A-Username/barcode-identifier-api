# Generated by Django 4.1.3 on 2023-02-18 21:11

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('barcode_blastn', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='blastrun',
            name='alignment_job_id',
            field=models.CharField(blank=True, default='', max_length=100),
        ),
        migrations.AddField(
            model_name='blastrun',
            name='complete_alignment_job_id',
            field=models.CharField(blank=True, default='', max_length=100),
        ),
        migrations.AddField(
            model_name='blastrun',
            name='create_db_tree',
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name='blastrun',
            name='create_hit_tree',
            field=models.BooleanField(default=True),
        ),
    ]
