from rest_framework.filters import OrderingFilter

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
        return [field_map.get(o, o) for o in ordering]