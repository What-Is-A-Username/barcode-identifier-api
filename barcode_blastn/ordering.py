from rest_framework.filters import OrderingFilter

from barcode_blastn.serializers import BlastDbSequenceEntrySerializer

class CustomSequenceOrderingFilter(OrderingFilter):
    def get_ordering(self, request, queryset, view):
        ordering = super().get_ordering(request, queryset, view)
        # Specify nested data to be organized by the scientific name, instead
        # of other fields such as the ID.
        field_map = {
            'taxon_species': 'taxon_species__scientific_name',
            'taxon_genus': 'taxon_genus__scientific_name',
            'taxon_family': 'taxon_family__scientific_name',
            'taxon_order': 'taxon_order__scientific_name',
            'taxon_phylum': 'taxon_phylum__scientific_name',
            'taxon_kingdom': 'taxon_kingdom__scientific_name',
            'taxon_superkingdom': 'taxon_superkingdom__scientific_name',
        }
        if not ordering is None:
            return [field_map.get(o, o) for o in ordering]
        return field_map.items()

class CustomHitOrderingFilter(OrderingFilter):
    def get_valid_fields(self, queryset, view, context={}):
        valid_fields = super().get_valid_fields(queryset, view, context=context)
        exclusions = ['taxon_superkingdom', 'taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus', 'taxon_species', 'annotations']
        extra_fields = [
            (f'db_entry__{field}', f'reference {field}') for field in BlastDbSequenceEntrySerializer.Meta.fields if field not in exclusions
        ]
        valid_fields.extend(extra_fields)
        return valid_fields

    def get_ordering(self, request, queryset, view):
        ordering = super().get_ordering(request, queryset, view)
        exclusions = ['taxon_superkingdom', 'taxon_kingdom', 'taxon_phylum', 'taxon_class', 'taxon_order', 'taxon_family', 'taxon_genus', 'taxon_species', 'annotations']
        extra_fields = [(f'db_entry__{field}', f'db_entry__{field}') for field in BlastDbSequenceEntrySerializer.Meta.fields if field not in exclusions]
        field_map = dict((k, v) for k, v in extra_fields)
        if not ordering is None:
            return [field_map.get(o, o) for o in ordering]
        return field_map.items()
