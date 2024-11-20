# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments,too-many-lines
# type: ignore   # for Pylance
"""
Utils for GFF3 gene annotation files from Ensembl.

See the GFF3 file in the Gene annotation section in
https://useast.ensembl.org/Homo_sapiens/Info/Index

See also the GFF3 readme for a description of the different GFF3 fields.
for example: https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/README.

Update on 6/5/24 - see below (search for 6/5/24).

Update on 11/20/24 - using ensembl_utils.
"""
import pathlib
import pandas as pd

import geneanot.ensembl_utils as eu

"""
GFF3 Type attributes
All values in the "Type" column (third column) of the GFF3 file:

['gene' 'mRNA' 'exon' 'five_prime_UTR' 'CDS' 'three_prime_UTR' 'lnc_RNA'
 'transcript' 'biological_region' 'ncRNA_gene' 'snRNA' 'pseudogene'
 'pseudogenic_transcript' 'unconfirmed_transcript' 'miRNA' 'snoRNA' 'rRNA'
 'chromosome' 'scRNA' 'C_gene_segment' 'J_gene_segment' 'V_gene_segment'
 'D_gene_segment' 'processed_transcript' 'scaffold' 'tRNA']
"""
# all Type values in GFF3 file
Gene_type_values: list = ['gene', 'ncRNA_gene', 'pseudogene']

# gene's transcript Types to process
# added 'transcript', 'unconfirmed_transcript', 'processed_transcript', on 6/5/24
Types_transcript_processed: list = ['mRNA', 'lnc_RNA', 'pseudogenic_transcript',
                                    'transcript', 'unconfirmed_transcript', 'processed_transcript',
                                    'ncRNA', 'miRNA', 'snoRNA', 'snRNA', 'scRNA', 'rRNA', 'tRNA']

# sub-transcript types to process
Types_in_transcript: list = ['CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR']


"""
Functions
"""
def load_ensembl_human_gff3_annotation_file(file: pathlib.Path, gene_type_values: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Loads the GFF3 file into a dataframe, and returns it and a subset of it containing only rows with column Type values in Gene_type_values."""
    try:
        df = pd.read_csv(file, sep='\t', comment='#', low_memory=False,
                      names=['Chrm', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'],
                      dtype={
                        'Chrm': str,
                         'Source': str,
                         'Type': str,
                         'Start': int,  #  1-based
                         'End': int,    #  1-based.
                         'Score': str,
                         'Strand': str,  # '+' for positive strand, '-' for negative strand.
                         'Phase': str,  # some of the rows have non-integer type
                         'Attributes': str
                     })
    except FileNotFoundError:
        print(f"Ensembl annotation file {file} not found !!")
        raise

    return df, get_gff3_df_gene_type_values(df, gene_type_values, type_col_name='Type', attr_col_name='Attributes')


def get_gff3_df_gene_type_values(df: pd.DataFrame, gene_type_values: list, type_col_name: str = 'Type', attr_col_name: str = 'Attributes') -> pd.DataFrame:
    """
    Given the GFF3 dataframe, the function returns a subset of it containing only rows with column Type values in Gene_type_values.
    Note that 'gene_id' and 'Name' strings are taken from the GFF3 file.
    """
    def extract_geneid_name_from_attributes(x: pd.Series, attributes_col_name: str) -> tuple[str, str]:
        """Extracts X and Y from 'gene_id=X' and 'Name=Y' fields in the Attributes column."""
        flds = dict([v.split('=') for v in x[attributes_col_name].split(';')])
        return flds.get('gene_id', ''), flds.get('Name', '')
    # ------------------------------------------------------------------------
    df_type = df.query(f"{type_col_name} in {gene_type_values}")  # dataframe containing only rows with Type value equal to one of the gene_type_values items

    # extract gene_id and Name from the Attributes column to two separate columns
    gid_nm = df_type.apply(extract_geneid_name_from_attributes, args=(attr_col_name, ), axis=1)
    # convert to dataframe
    d = pd.DataFrame(gid_nm.to_list(), columns=['gene_id', 'Name'], index=gid_nm.index)
    return df_type.merge(d, left_index=True, right_index=True)

def display_gene_annotation_structure(gene_s: dict, transcripts_label: str = 'transcripts'):
    """Displays structure"""
    for k, v in gene_s.items():
        if k == transcripts_label:
            for i, y in enumerate(v, start=1):
                print(f"Transcript #{i}:")
                for k1, v1 in y.items():
                    print(f"\t{k1}: {v1}")
        else:
            print(f"{k}: {v}")

def gather_transcript_attributes(
        gene_rev: bool,
        mrna_name: str, transcript: str, transcript_v: str, ccdsid: str, mrna_start: int, mrna_end: int,
        mrna_ver: int, mrna_biotype: str, mrna_type: str,
        c_exon_start: list, c_exon_end: list, c_5UTR_start: list, c_5UTR_end: list,
        c_3UTR_start: list, c_3UTR_end: list, c_CDS_start: list, c_CDS_end: list,
        c_CDS_phase: list, c_CDS_protein_id: list, exon_id: list, exon_phase: list,
        exon_end_phase: list, exon_ver: list) -> dict:
    """
    If gene encoded on the negative strand (i.e., gene_rev == True), the order of all the start and end from the GFF3 file
    is also from low to high (as in the case of gene encoded on the positive strand). However, the actual order
    should be reversed.

    For example, given the following exon_start and exon_end of a gene encoded on the negative strand
    (this is what is given in the GFF3 file):

        exon_start = [100, 250, 360]
        exon_end   = [140, 310, 420],

    then relative to the reverse strand, the first exon is from 420 to 360, the second exon is from 310 to 250, and the third (last)
    exon is from 140 to 100.
    Thus:

        a_exon_start = exon_end[::-1] = [420, 310, 140]  # the actual exon start list
        a_exon_end = exon_start[::-1] = [360, 250, 100]  # the actual exon end list

    This is also the order and values of the transcripts in the Ensembl web
    (see, e.g., https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000138413;r=2:208236229-208266074;t=ENST00000345146
    for the IDH1 gene, which is encoded on the negtive strand.)
    """
    # all protein IDs are the same for the same CDS.
    p_id = '' if c_CDS_protein_id == [] else c_CDS_protein_id[0]

    if gene_rev:
        a_mrna_start = mrna_end
        a_mrna_end = mrna_start

        a_c_exon_start = c_exon_end[::-1]
        a_c_exon_end = c_exon_start[::-1]
        a_exon_id = exon_id[::-1]

        # these two should not be switched but only reversed
        a_exon_phase = exon_phase[::-1]  # exon_end_phase[::-1]
        a_exon_end_phase = exon_end_phase[::-1]  # exon_phase[::-1]

        a_exon_ver = exon_ver[::-1]

        a_c_CDS_start = c_CDS_end[::-1]
        a_c_CDS_end = c_CDS_start[::-1]
        a_c_CDS_phase = c_CDS_phase[::-1]

        a_c_5UTR_start = c_5UTR_end[::-1]
        a_c_5UTR_end = c_5UTR_start[::-1]

        a_c_3UTR_start = c_3UTR_end[::-1]
        a_c_3UTR_end = c_3UTR_start[::-1]
    else:
        a_mrna_start = mrna_start
        a_mrna_end = mrna_end

        a_c_exon_start = c_exon_start
        a_c_exon_end = c_exon_end
        a_exon_id = exon_id

        a_exon_phase = exon_phase
        a_exon_end_phase = exon_end_phase

        a_exon_ver = exon_ver

        a_c_CDS_start = c_CDS_start
        a_c_CDS_end = c_CDS_end
        a_c_CDS_phase = c_CDS_phase

        a_c_5UTR_start = c_5UTR_start
        a_c_5UTR_end = c_5UTR_end

        a_c_3UTR_start = c_3UTR_start
        a_c_3UTR_end = c_3UTR_end

    return {
        'transcript_name': mrna_name,
        'transcript_id': transcript,
        'transcript_v_id': transcript_v,
        'transcript_ccdsid': ccdsid,
        'transcript_start': a_mrna_start,
        'transcript_end': a_mrna_end,
        'transcript_ver': mrna_ver,
        'transcript_type': mrna_type,
        'transcript_biotype': mrna_biotype,

        'exon_start': a_c_exon_start,
        'exon_end': a_c_exon_end,
        'exon_id': a_exon_id,
        'exon_phase': a_exon_phase,
        'exon_end_phase': a_exon_end_phase,
        'exon_ver': a_exon_ver,

        'CDS_start': a_c_CDS_start,
        'CDS_end': a_c_CDS_end,
        'CDS_phase': a_c_CDS_phase,
        'protein_id': p_id,

        '5UTR_start': a_c_5UTR_start,
        '5UTR_end': a_c_5UTR_end,

        '3UTR_start': a_c_3UTR_start,
        '3UTR_end': a_c_3UTR_end
    }

def parse_ensembl_gene_anot_gene(df: pd.DataFrame, types_transcript_processed: list, types_in_transcript: list, verbose: bool = True) -> dict:
    """
    This function parses a dataframe (with columns from the GFF3 file), whose first column has
    a Type 'gene' or 'ncRNA_gene' and generates the annotation dictionary of the corresponding gene.

    df - a dataframe containing the GFF3 data.
    types_transcript_processed - types of transcripts (from the 'type' column) to process
                                 (e.g., mRNA, lnc_RNA, ncRNA, miRNA)
    types_in_transcript - types of sub-transcript (from the 'type' column)
                          to process (e.g., exon, CDS, etc.).
    """
    gn_info = df.iloc[0]

    if (gn_info['Type'] != 'gene') and (gn_info['Type'] != 'ncRNA_gene'):
        if verbose:
            print("First row in df must have Type of gene or ncRNA_gene !!")
        return {}

    # convert Attributes string to dict
    d = dict([tuple(x.split('=')) for x in gn_info['Attributes'].split(';')])

    gene_chrm, gene_start, gene_end, gene_rev = gn_info['Chrm'], int(gn_info['Start']), int(gn_info['End']), (
        False if gn_info['Strand'] == '+' else True)
    try:
        gene_name = d['Name']
    except KeyError:
        # non-coding genes do not have "Name=<name>" string in the Attributes column.
        # In this case use gene_id as gene name
        # changed on 6/5/24 - using None for gene name in these cases
        #gene_name = d['gene_id']
        gene_name = None
    gene_ID = d['gene_id']
    gene_type = d['biotype']
    try:
        gene_desc = d['description']
    except KeyError:
        gene_desc = ''
    try:
        gene_ver = int(d['version'])
    except KeyError:
        gene_ver = -1

    # parse subsequent rows and get transcript information
    num_mrnas, gene_transcript_strct, in_mrna = 0, [], False

    # this is not really necessary, as these are initialized at a transcript level, before
    # being used in a sub-transcript level.
    c_exon_start, c_exon_end, c_5UTR_start, c_5UTR_end, c_3UTR_start, c_3UTR_end, c_CDS_start, c_CDS_end = [], [], [], [], [], [], [], []
    c_CDS_phase, c_CDS_protein_id, exon_id, exon_phase, exon_end_phase, exon_ver = [], [], [], [], [], []
    mrna_name, transcript, transcript_v, ccdsid = '', '', '', ''
    mrna_start, mrna_end, mrna_ver, mrna_biotype, mrna_type = '', '', '', '',  ''

    # recall that the first row is the gene attributes
    for _, r in df.iloc[1:].iterrows():
        # The assumption here is that the beginning of a gene in the GFF3 file is defined
        # by Type='gene' (protein coding genes), Type='ncRNA_gene' (non-protein coding genes),
        # or Type='pseudogene'.
        if isinstance(r['Type'], str) and 'gene' in r['Type']:
            # a new gene, need to exit
            break

        # any processed Type values
        if r['Type'] in (types_in_transcript + types_transcript_processed):

            # Type is a transcript type (mRNA, ncMRNA, tRNA, etc), but not a sub-mRNA type (e.g. CDS, exon, ..)
            if r['Type'] in types_transcript_processed:
                in_mrna = True
                # starting a new transcript

                if num_mrnas > 0:
                    # gather previous transcript attributes
                    gene_transcript_strct.append(gather_transcript_attributes(
                        gene_rev,
                        mrna_name, transcript, transcript_v, ccdsid, mrna_start, mrna_end, mrna_ver, mrna_biotype,
                        mrna_type,
                        c_exon_start, c_exon_end, c_5UTR_start, c_5UTR_end,
                        c_3UTR_start, c_3UTR_end, c_CDS_start, c_CDS_end,
                        c_CDS_phase, c_CDS_protein_id, exon_id, exon_phase, exon_end_phase, exon_ver
                    ))

                num_mrnas += 1

                # parse the new transcript
                d1 = dict([tuple(x.split('=')) for x in r['Attributes'].split(';')])

                # verify that the new transcript corresponds to the gene
                try:
                    if d1['Parent'] != f"gene:{gene_ID}":
                        if verbose:
                            print(f"Transcript (type={r['Type']}), ({d1}) does not seems to belongs to {gene_ID} !! Exiting ...")
                        return {}
                except KeyError:
                    pass

                # the mrna_ names below should actually be called transcript_
                transcript = d1['transcript_id']
                transcript_v = f"{transcript}.{d1['version']}"
                mrna_type = r['Type']
                try:
                    mrna_name = d1['Name']
                except KeyError:
                    mrna_name = d1['transcript_id']
                mrna_biotype = d1['biotype']
                mrna_start = int(r['Start'])
                mrna_end = int(r['End'])
                mrna_ver = int(d1['version'])
                try:
                    ccdsid = d1['ccdsid']
                except KeyError:
                    ccdsid = None

                # print(f"\t{transcript=}, {transcript_v=}, {mrna_name=}, {ccdsid}, {mrna_start=}, {mrna_end=}")

                # for the next rows in the dataframe, which consist of the transcript sub-types (exon, CDS, etc)
                c_exon_start, c_exon_end, c_5UTR_start, c_5UTR_end, c_3UTR_start, c_3UTR_end, c_CDS_start, c_CDS_end = [], [], [], [], [], [], [], []
                c_CDS_phase, c_CDS_protein_id, exon_id, exon_phase, exon_end_phase, exon_ver = [], [], [], [], [], []
            else:
                # Type is in types_in_transcript. Updating current transcript
                # don't parse if the transcript is not in types_transcript_processed
                if not in_mrna:
                    continue

                # verify that this belongs to the current mRNA
                d2 = dict([tuple(x.split('=')) for x in r['Attributes'].split(';')])
                try:
                    if d2['Parent'] != f"transcript:{transcript}":
                        if verbose:
                            print(f"{r['Type']} ({d2}) does not seems to belong to {transcript} !! Exiting ...")
                        return {}
                except KeyError:
                    pass

                match r['Type']:
                    case 'five_prime_UTR':
                        c_5UTR_start.append(int(r['Start']))
                        c_5UTR_end.append(int(r['End']))
                    case 'exon':
                        c_exon_start.append(int(r['Start']))
                        c_exon_end.append(int(r['End']))

                        exon_id.append(d2['exon_id'])
                        exon_end_phase.append(int(d2['ensembl_end_phase']))
                        exon_phase.append(int(d2['ensembl_phase']))
                        exon_ver.append(int(d2['version']))
                    case 'CDS':
                        c_CDS_start.append(int(r['Start']))
                        c_CDS_end.append(int(r['End']))
                        c_CDS_phase.append(int(r['Phase']))

                        c_CDS_protein_id.append(d2['protein_id'])  # this should be the same for all CDSs of a mRNA
                    case 'three_prime_UTR':
                        c_3UTR_start.append(int(r['Start']))
                        c_3UTR_end.append(int(r['End']))
                    case _:
                        print(f"Type {r['Type']} in {gene_name=}, ({mrna_name}) is not supported !!")
                        continue
        else:
            # print(f"Type={r['Type']} currently not supported.")
            # don't process any (next) rows with Type in types_in_transcript, until encountering
            # a row with Type in types_transcript_processed
            in_mrna = False

    # gather last transcript attributes
    if num_mrnas > 0:
        gene_transcript_strct.append(gather_transcript_attributes(
            gene_rev,
            mrna_name, transcript, transcript_v, ccdsid, mrna_start, mrna_end, mrna_ver, mrna_biotype, mrna_type,
            c_exon_start, c_exon_end, c_5UTR_start, c_5UTR_end,
            c_3UTR_start, c_3UTR_end, c_CDS_start, c_CDS_end,
            c_CDS_phase, c_CDS_protein_id, exon_id, exon_phase, exon_end_phase, exon_ver
        ))

    # construct final gene structure
    return {
        # changed on 6/5/24
        #'gene': gene_name,
        'gene_name': gene_name,
        'chrm': gene_chrm,
        'gene_start': gene_start,
        'gene_end': gene_end,
        'rev': gene_rev,
        'gene_ID': gene_ID,
        'gene_type': gene_type,
        'gene_desc': gene_desc,
        'gene_ver': gene_ver,
        'num_transcripts': len(gene_transcript_strct),
        'transcripts': gene_transcript_strct
    }


def get_df_start_with_gene(df: pd.DataFrame, df_genes: pd.DataFrame, gene: str, verbose: bool = True) -> pd.DataFrame | None:
    """
    The function returns the dataframe where the first row contains Type gene
    of the corresponding gene.

    df - GFF3 file loaded to a dataframe
    df_genes - a subset of df containing only rows with Type values that is in Gene_type_values (see the function load_ensembl_human_gff3_annotation_file)
    """
    fld = 'gene_id' if eu.is_id(gene) else 'Name'
    if (df_gene := df_genes.query(f"{fld} == '{gene}'")).empty:
        if verbose:
            print(f"{gene=} not found in annotation dataframe !!")
        return None

    # this should work as well, instead of the 3 lines below (since the indexes in df_genes are a subset of the indexes in df)
    #return df.loc[df_gene.index[0]:]

    # the label index of the row that contains Type of 'gene', corresponding to gene
    idx = df_gene.index

    # the corresponding index in df
    i_idx = df.index.get_indexer(idx)

    # first row is the row with Type value in gene_types_to_include of the corresponding gene
    return df.iloc[i_idx[0]:]


def extract_gene_dict(df: pd.DataFrame, df_genes: pd.DataFrame, gene: str, verbose: bool = True) -> dict:
    """Finds and extracts gene data from GFF3 dataframe."""
    if (df_gene_all := get_df_start_with_gene(df, df_genes, gene, verbose=verbose)) is None:
        if verbose:
            print(f"No entry for {gene=} found in annotation dataframe!!")
        return {}

    # verify
    d = dict([tuple(x.split('=')) for x in df_gene_all.iloc[0]['Attributes'].split(';')])
    try:
        # first_row_gene = d['Name']
        # update on 6/4/24 - to enable 'gene' to be either a gene name (Hugo symbol) or a gene ID (i.e., ENS<species prefix>G, where <species prefix> is empty for Homo sapiens)
        first_row_gene = d['gene_id'] if eu.is_id(gene) else d['Name']

    except KeyError:
        # non-coding genes have no "Name" field. In these cases we use the gene ID as the name.
        first_row_gene = d['gene_id']

    if gene != first_row_gene:
        if verbose:
            print(f"Error in extracting gene information (first row in extracted dataframe (with {first_row_gene}) does not match {gene=} !!")
        return {}

    # generate the annotation dictionary
    return parse_ensembl_gene_anot_gene(df_gene_all, Types_transcript_processed, Types_in_transcript, verbose=verbose)
