# Generated by Django 4.1.3 on 2023-02-18 21:11

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('barcode_tree', '0003_remove_resulttree_tree_job_id'),
    ]

    operations = [
        migrations.AddField(
            model_name='resulttree',
            name='complete_alignment_job_id',
            field=models.CharField(blank=True, default='', max_length=100),
        ),
    ]