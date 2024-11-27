# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments,too-many-lines,too-many-instance-attributes
# type: ignore   # for Pylance
"""
Utils related to the GFF3 annotation file.
"""
import re
import pathlib
import subprocess

import geneanot.ensembl_rest_utils as erut

def suggested_annotation_file_name(species: str = 'homo_sapiens',
                                   file_type: str = 'gff3', 
                                   rest_assembly: str = 'GRCh38') -> tuple[str,str,str] | None:
    """Returns the annotation file name based on the latest Ensembl release, and the provided species and file_type."""
    rapi = erut.REST_API(assembly=rest_assembly)
    try:
        assembly_name = rapi.get_assembly_info(species=species)["default_coord_system_version"]
        release_number = str(rapi.get_release_info()['releases'][0])
    except KeyError as ke:
        print(f"Error in retreiving assembly and release informtion: {ke}")
        return None
    return f"{species.capitalize()}.{assembly_name}.{release_number}.{file_type}.gz", assembly_name, release_number

def list_files_in_folder(folder: pathlib.Path,
                         file_pattern: str) -> list[pathlib.Path]:
    """Returns a list of files (with a given signature) found in a folder.

    Args:
        folder (pathlib.Path): search folder.
        file_pattern (str, optional): file regex signature. E.g. 'Homo_sapiens.GRCh38.(\\d+).gff3.gz'.

    Returns:
        list[pathlib.Path]: list of annotation files found.
    """
    return [x for x in folder.glob("*.*") if re.match(file_pattern, x.name)]

def get_ensembl_release() -> int:
    """Returns Ensembl current release."""
    return int(erut.REST_API().get_release_info()['releases'][0])

def get_annotation_files_info_in_folder(folder: pathlib.Path,
                                        gff3_pattern: str = 'Homo_sapiens.GRCh38.XXX.gff3.gz') -> dict[int,str]:
    """Get annotation files info in local folder.

    Args:
        folder (pathlib.Path): local folder.
        gff3_pattern (str, optional): annotation file pattern. Defaults to 'Homo_sapiens.GRCh38.XXX.gff3.gz'. XXX is a place holder for the version.

    Returns:
        dict[int,str]: key is version number, value is correponding file name.
    """
    # create a regex pattern
    #gff3_match_pattern = gff3_pattern.replace('XXX', '(\d+)').replace('.', '\.')
    # gff3_match_pattern = gff3_pattern.replace('XXX', '(\d+)')
    gff3_match_pattern = gff3_pattern.replace('XXX', r'(\d+)')

    if (files := list_files_in_folder(folder, file_pattern=gff3_match_pattern)):
        return {int(re.search(f"{gff3_match_pattern}", file.name).group(1)): file.name for file in files}
    return {}

def report_annotation_files_in_folder_and_in_ensembl(local_gff3_folder: pathlib.Path,
                                                     gff3_pattern: str = 'Homo_sapiens.GRCh38.XXX.gff3.gz') -> tuple[dict, int]:
    """Reports all local annotation files, and Ensembl remote version.

    Args:
        local_gff3_folder (pathlib.Path): local folder.
        gff3_pattern (str, optional): annotation file pattern. Defaults to 'Homo_sapiens.GRCh38.XXX.gff3.gz'.

    Returns:
        tuple[dict int]: [0] - local annotation information. key is version number, value is correponding file name.
                         [1] - latest Ensembl version.
    """
    # key is version number, value is correponding file name.
    local_version_info = get_annotation_files_info_in_folder(local_gff3_folder, gff3_pattern=gff3_pattern)

    # latest annotattion version in Ensembl
    ensm_rel = get_ensembl_release()

    return local_version_info, ensm_rel

def update_local_release_to_latest(local_gff3_folder: pathlib.Path,
                                   enable_download: bool,
                                   gff3_pattern: str = 'Homo_sapiens.GRCh38.XXX.gff3.gz',
                                   species: str = 'homo_sapiens',
                                   ) -> tuple[bool, str, str]:
    """Checks local GFF3 annotation file version vs Ensembl latest version, and downloads Ensembl's file if needed (and download is enabled).

    1. Chceks for annotation files in local folder and report.
    2. Checks Ensembl current release and report.
    3. If no local annotation files found - download Ensembl file. Use the 2nd returned output as the annotation file.
    4. If local annotation files found, but latest version older than Ensemblversion: if enable_download == True, download Ensembl version, (use the 2nd returned output as the annotation file). If enable_download == False, download is not done (use 3rd returned output as the annotation file).
    5. If local annotation files found, and latest version is identical to Ensembl version, report that no update is needed (use the 3rd returned output as the annotation file).

    Args:
        local_gff3_folder (pathlib.Path): local folder holding (or to hold) GFF3 annotation files.
        gff3_pattern (str, optional): GFF3 file pattern. Defaults to 'Homo_sapiens.GRCh38.XXX.gff3.gz', where XXX is a place holder for Ensembl version.
        species (str, optional): Defaults to 'homo_sapiens'.

    Returns:
        tuple[bool, str, str]: [0] - True if downloaded the latest Ensembl GFF3 file to local folder, otherwise False.
                               [1] - latest Ensembl GFF3 file.
                               [2] - latest GFF3 file found in local folder (PRIOR to the download, if [0] == True). 
    """
    local_version_info, ensm_rel = report_annotation_files_in_folder_and_in_ensembl(local_gff3_folder, gff3_pattern=gff3_pattern)

    # local files
    if local_version_info:
        max_local_version = max(local_version_info.keys())
        local_gff3_max_ver_file = local_version_info[max_local_version]
        print(f"Local annotation releases found: {', '.join(sorted(map(str, local_version_info.keys())))}. Latest is {max_local_version} ({local_gff3_max_ver_file}).")
    else:
        max_local_version, local_gff3_max_ver_file = None, ''
        print("No local annotation files found.")

    # Ensembl file
    print(f"Ensembl latest release is {ensm_rel}.")
    ensembl_gff3_latest_file = gff3_pattern.replace('XXX', str(ensm_rel))

    if ensm_rel == max_local_version:
        print("Ensembl and local releases match.")
        return False, ensembl_gff3_latest_file, local_gff3_max_ver_file

    ensembl_url: str = f'rsync://ftp.ebi.ac.uk/ensemblorg/pub/current_gff3/{species}'
    # cases to consider now are no local version, or local version < ensembl version
    if (max_local_version is None) or enable_download:
        # no local version, or (local version < ensembl version and enable_download == True). Thus need to download Ensembl file
        if max_local_version is not None:
            print("Local latest release is older than latest Ensembl release.")
        print(f"Downloading Ensembl latest annotation file ({ensembl_gff3_latest_file})...")

        ret = subprocess.run(f"rsync -av {ensembl_url}/{ensembl_gff3_latest_file} {str(local_gff3_folder)}", shell=True, capture_output=True, check=False)
        if (ret_success := (ret.returncode == 0)):
            print(f"Ensembl {ensembl_gff3_latest_file} file downloaded to {local_gff3_folder} successfully.")
        else:
            print(f"Error in downloading Ensembl {ensembl_gff3_latest_file} file from {ensembl_url}:")
            print(f"STDOUT:\n{ret.stdout.decode()}")
            print(f"STDERR:\n{ret.stderr.decode()}")
        return ret_success, ensembl_gff3_latest_file, local_gff3_max_ver_file

    # local version < ensembl version and enable_download == False
    print(f"Download did not occur (since {enable_download=}). Use local {local_gff3_max_ver_file} annotation file.")
    return False, ensembl_gff3_latest_file, local_gff3_max_ver_file
