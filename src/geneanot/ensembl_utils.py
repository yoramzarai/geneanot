# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments,too-many-lines
# type: ignore  # for Pylance
"""
Ensembl related utils.
"""
import re

# pattern for Ensembl stable ID
Ensb_pattern: str = r"ENS([A-Z]*)([A-Z])(\d{11})"
# Pattern for Ensembl stable ID with version
# Ensb_pattern_ver: str = r"ENS([A-Z]*)([A-Z])(\d{11}).(\d*)"
# Pattern for Ensembl stable ID with or without a version
Ensb_pattern_ver: str = r"ENS([A-Z]*)([A-Z])(\d{11})(\.\d+)?"

# See http://www.ensembl.org/info/genome/stable_ids/prefixes.html
Ensb_feature_mapping: dict[str,str] = {
    "E": "Exon",
    "G": "Gene",
    "T": "Transcript",
    "P": "Protein",
    "R": "Regulatory Feature",
}

def parse_ensembl_id(ensb_id: str) -> dict:
    """Parses Ensembl stable ID.
    The function does NOT support Ensembl protein family (FM) and gene tree (GT) feature prefixes.

    Args:
        ensb_id (str): Ensembl stable ID (with or without a version).

    Returns:
        dict: contains the feature, species, ID and version fields.
    """
    if (res := re.match(Ensb_pattern_ver, ensb_id)) is None:
        return {}
    feature_prefix, ver = res.group(2), res.group(4)

    return {
        "feature_prefix": feature_prefix,
        "species_prefix": res.group(1) or None,  # None implies Homo sapiens
        "ID": res.group(3),
        "feature_name": Ensb_feature_mapping.get(feature_prefix, '?'),
        "version": ver[1:] if ver else None
    }

def is_id(ensb_id: str) -> bool:
    """Checks if a valid Ensembl ID.

    Args:
        ensb_id (str): An ID.

    Returns:
        bool: True if an Ensembl ID
    """
    return bool(parse_ensembl_id(ensb_id))

def _is_feature_id(ensb_id: str, feature_prefix: str) -> bool:
    """Checks if ensembl ID correponds to a given feature prefix.

    Args:
        ensb_id (str): Ensembl stable ID.
        feature_prefix (str): one of the feature prefixes.

    Returns:
        bool: True if ensb_id corresponds to a feature feature_prefix.
    """
    if not (p := parse_ensembl_id(ensb_id)):
        return False
    return p['feature_prefix'] == feature_prefix

def is_id_gene_id(ensb_id: str) -> bool:
    """Checks if ensb_id is a gene ID."""
    return _is_feature_id(ensb_id, 'G')

def is_id_transcript_id(ensb_id: str) -> bool:
    """Checks if ensb_id is a transcript ID."""
    return _is_feature_id(ensb_id, 'T')

def is_id_protein_id(ensb_id: str) -> bool:
    """Checks if ensb_id is a protein ID."""
    return _is_feature_id(ensb_id, 'P')

def is_id_exon_id(ensb_id: str) -> bool:
    """Checks if ensb_id is an exon ID."""
    return _is_feature_id(ensb_id, 'E')

def is_id_regulatory_id(ensb_id: str) -> bool:
    """Checks if ensb_id is a regulatory ID."""
    return _is_feature_id(ensb_id, 'R')
