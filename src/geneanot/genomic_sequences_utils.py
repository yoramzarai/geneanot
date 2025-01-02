# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments,too-many-lines,too-many-instance-attributes,too-many-locals
# type: ignore   # for Pylance
""" 
Utils related to extracting genomic sequences.
"""
from pathlib import Path
import numpy as np
from .Fasta_segment import Fasta_segment
from .ensembl_rest_utils import REST_API
import geneanot.translation as tran


def fetch_seq(source_info: Path | str, start_p: int, end_p: int, rev: bool = False, species: str = 'homo_sapiens', rest_assembly: str = 'GRCh38') -> str:
    """Fetch sequence either from a local Fasta file or from the Ensembl REST API.

    The fetch algorithm is:

                if isinstance(source_info, Path):
                    # fetch from the chromosome Fasta file source_info
                else:  # isinstance(source_info, str)
                    # fetch from Ensembl REST API from chromosome number source_info

    Args:
        source_info (Path | str): a chromosome fasta file (of type Path) or the chromosome number of type str (e.g., '1' or 'Y').
        start_p (int): 1-based start coordinate
        end_p (int): 1-based end coordinate
        rev (bool, optional): reverse-complement the sequence. Defaults to False.
        species (str, optional): species (relevant for isinstance(source_info, str)). Defaults to 'homo_sapiens'.
        rest_assembly (str, optional): REST assembly (relevant for isinstance(source_info, str). Defaults to 'GRCh38'.

    Returns:
        str: the fetched sequence.
    """
    if isinstance(source_info, Path):
        # load from a local chromosome Fasta file 
        return extract_fasta_seq(str(source_info), start_p, end_p, rev=rev)
    elif isinstance(source_info, str):
        # fetch from Ensembl REST API
        return extract_chromosome_seq(source_info, start_p, end_p, rev=rev, species=species, rest_assembly=rest_assembly)
    return ''

def extract_fasta_seq(
        fasta_file: str,
        start_p: int,   # 1-based
        end_p: int,  # 1-based
        rev: bool = False) -> str:
    """
    Extracts a Fasta sequence.

    This function can be used ONLY with Fasta files where ALL rows
    (excluding the headers, and possibly the last row) have the same number of characters!

    Args:
        fasta_file (str): Fasta file name.
        start_p (int): 1-based start position
        end_p (int): 1-based end position
        rev (bool, optional): True to get the reveresed-complement sequence (i.e. the sequence in the negative strand). Defaults to False.

    Returns:
        str: The extracted sequence.
    """
    # read_segment gets 0-based position, and size of the sequence
    seq = Fasta_segment().read_segment(fasta_file, start_p-1, end_p - start_p + 1)
    return tran.reverse_complement(seq) if rev else seq

def extract_chromosome_seq(
        chrm: str,  # without 'chr' prefix, e.g., '3', or 'X'
        start_p: int,
        end_p: int,
        rev: bool = False,
        species: str = 'homo_sapiens',
        rest_assembly: str = 'GRCh38'
        ) -> str:
    """Extracts a chromosome sequence using Ensembl REST API.

    Args:
        chrm (str): chromosome number. E.g., '1' or 'X'.
        start_p (int): 1-based start position.
        end_p (int): 1-based end position.
        rev (bool, optional): True to extract the reveresed-complement sequence 
                             (i.e. the sequence in the negative strand). Defaults to False.
        species (str, optional): Defaults to 'homo_sapiens'.
        assembly (str, optional): Determines the REST URL. Defaults to 'GRCh38'.

    Returns:
        str: The extracted sequence.
    """
    if 'chr' in chrm:
        chrm = chrm[3:]
    strand = -1 if rev else 1
    return REST_API(assembly=rest_assembly).sequence_region_endpoint_base(chrm, start_p, end_p, strand=strand, species=species, content_type='text/plain')

def pos2seg_info(pos: int, seg_start: list, seg_end: list, rev: bool) -> tuple[int, int, int] | None:
    """
    Given a position pos, and segments defined by start and end lists, the function returns
    the tuple (x, y, z), where
    x - the 0-based segment number containing pos
    y - the 0-based position within the segment x corresponding to pos
    z - the 0-based offset from the first segment, when concatenating all segments, corresponding to pos

    Note that seg_start[i] > seg_end[i] for rev==True, otherwise seg_start[i] < seg_end[i].
    """
    # converting so that start < end
    s_start, s_end = (seg_end, seg_start) if rev else (seg_start, seg_end)

    try:
        # 0-based
        seg_indx = np.argwhere((np.array(s_start) <= pos) & (pos <= np.array(s_end)))[0][0]
    except IndexError:
        return None

    seg_sizes = [e - s + 1 for s, e in zip(s_start, s_end)]
    dist_from_seg = (seg_start[seg_indx] - pos) if rev else (pos - seg_start[seg_indx])
    dist_from_first_seg = sum(seg_sizes[:seg_indx]) + dist_from_seg

    return int(seg_indx), int(dist_from_seg), int(dist_from_first_seg)