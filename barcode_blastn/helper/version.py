from typing import Union
from barcode_blastn.models import Library, BlastDb

def get_new_version(lib: Library, new: BlastDb) -> tuple:
    last: Union[BlastDb, None] = BlastDb.objects.latest(library=lib)
    # TODO: Compute whether new major and minor versions are needed
    if last is None:
        return (254, 1, 1)
    else:
        return (254, last.major_version + 1, last.minor_version)

def assign_new_version(lib: Library, new: BlastDb):
    version_nos = get_new_version(lib, new)
    new.genbank_version = version_nos[0]
    new.major_version = version_nos[1]
    new.minor_version = version_nos[2]
    return
