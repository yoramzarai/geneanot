
from .Gene_annotation import Gene_cls, ensembl_gff3_df
from .annotation_file_utils import (
    suggested_annotation_file_name,
    get_annotation_files_info_in_folder, 
    get_ensembl_release, 
    update_local_release_to_latest
)
from .genomic_sequences_utils import fetch_seq
from importlib.metadata import version

__all__ = [
    "Gene_cls", 
    "ensembl_gff3_df",
    "suggested_annotation_file_name",
    "get_annotation_files_info_in_folder",
    "get_ensembl_release",
    "update_local_release_to_latest",
    "fetch_seq",
]
# from geneanot.Gene_annotation import *
# from geneanot.annotation_file_utils import *
# from .Gene_annotation import *
# from .annotation_file_utils import *

# an example of an entry point. See [project.scripts] in pyproject.toml
def main() -> None:
    print(f"Welcome to the vertebrates gene annotation package geneanot version {version('geneanot')}.\n\
Please visit https://github.com/yoramzarai/geneanot for documentation.")