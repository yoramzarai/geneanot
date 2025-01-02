# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments,too-many-lines,too-many-instance-attributes,too-many-locals
# type: ignore   # for Pylance
"""
Vertebrates gene Ensembl annotation.

Based on Ensembl GFF3 annotation file.
See the GFF3 readme for a description of the different GFF3 fields,
for example: https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/README
"""
from pathlib import Path
from dataclasses import dataclass, field
import numpy as np
import pandas as pd

import geneanot.translation as tran
#from . import translation as tran
import geneanot.ensembl_gene_annotations_utils as egna
from geneanot.excel_utils import dfs_to_excel_file, create_excel_description_sheet
#from .excel_utils import dfs_to_excel_file, create_excel_description_sheet
from geneanot.genomic_sequences_utils import fetch_seq, pos2seg_info

def ensembl_gff3_df(file: Path, gene_type_values: list = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse the annotation file into two dataframes, which can be used to instantiate the annotation class.

    file: GFF3 file path.
    """
    if gene_type_values is None:
        gene_type_values = egna.Gene_type_values 
    return egna.load_ensembl_human_gff3_annotation_file(file, gene_type_values)

@dataclass
class Transcript_gff3_cls:
    """
    Transcript class derived from a gene name and the GFF3 Ensembl annotation file.

    This is basically the base class for the gene annotation class.

    The class supports two chromosome data access modes:
    1. local - user provides the corresponding chromosome Fasta file.
    2. remote - the class uses the Ensembl REST API to extract sequence from the chromosome.

    In case of a local access mode, the provided Fasta file MUST contain equal bps per rows in all sequence
    rows (other than possibly the last row).

    Instantiating this class requires the following inputs:
    1. gene name or gene ID (ENS<species prefix>G, where <species prefix> is empty for Homo sapiens).
    2. Either the GFF3 file (Path), or a tuple of the GFF3 dataframe and its subset dataframe
       (containing only rows with Type value defined by the user (deafult is egna.Gene_type_values)).
       Use the function ensembl_gff3_df to generate this tuple if instantiating using the tuple.
    3. Optional - verbose flag.
    4. Optional - species
    5. Optional - chromosome Fasta file.

    All other class members are set by the __post__init__ method.

    The pharse "transcript" below (in the context of a sequence or cof a boundary) refers to
    the primary transcript (pre-RNA sequence).
    """
    gene: str  # gene name or gene ID
    gff3_source: Path | tuple  # Ensembl GFF3 file name, or a tuple of the GFF3 dataframe and its gene type subset
    species: str = 'homo_sapiens'  # needed when using remote chromosome sequence extraction
    chrm_fasta_file: Path | None = None  # if None using remote sequence extraction, otherwise using this chromosome fasta file
    verbose: bool = True  # set to False at instantiation to suppress prints

    # user should not set this. I need to change this when a new REST API URL is available.
    _rest_assembly: str = field(init=False, default='GRCh38')

    # these are set by the __post_init__ method.
    gff3_df: pd.DataFrame = field(init=False)
    gff3_df_gene_type: pd.DataFrame = field(init=False)

    gene_ID: str = field(init=False)
    gene_name: str = field(init=False)
    chrm: str = field(init=False)
    rev: bool = field(init=False)
    gene_start: int = field(init=False)
    gene_end: int = field(init=False)
    gene_type: str = field(init=False)
    gene_desc: str = field(init=False)
    gene_ver: int = field(init=False)
    transcripts: list = field(init=False)
    source: str = field(init=False)

    transcripts_info: dict = field(init=False)  # keys are transcripts, values are transcript details
    exon_intron_maps: dict = field(init=False)  # keys are transcripts, values are maps in dataframe format

    __protein_coding_labels_in_biotype: list | None = None

    def __post_init__(self):
        if isinstance(self.gff3_source, Path):  # instantiating by GFF3 file name
            self.source = str(self.gff3_source)
            if self.verbose:
                print(f"Loading {self.gff3_source} to a dataframe ... ", end='')
            self.gff3_df, self.gff3_df_gene_type = egna.load_ensembl_human_gff3_annotation_file(self.gff3_source, egna.Gene_type_values)
            if self.verbose:
                print("Done.")
        else:  # instantiating by a tuple containing the GFF3 dataframe and its subset
            self.gff3_df, self.gff3_df_gene_type = self.gff3_source
            self.source = f'Input DataFrame (ID={id(self.gff3_df)}, size={self.gff3_df.shape}). (Gene type DataFrame (ID={id(self.gff3_df_gene_type)}, size={self.gff3_df_gene_type.shape}))'

        # extract the gene gene from the GFF3 dataframe and generate the description dictionary
        if not (d := egna.extract_gene_dict(self.gff3_df, self.gff3_df_gene_type, self.gene, verbose=self.verbose)):
            raise ValueError(f"Error in extracting the gene {self.gene} from GFF3 dataframe !!")

        self.chrm = d['chrm']  # the chromosome number in str form (e.g. '3')
        self.gene_start = d['gene_start']
        self.gene_end = d['gene_end']
        self.rev = d['rev']
        self.gene_ID = d['gene_ID']
        self.gene_name = d['gene_name']
        self.gene_type = d['gene_type']
        self.gene_desc = d['gene_desc']
        self.gene_ver = d['gene_ver']
        self.transcripts = d['transcripts']

        if self.__protein_coding_labels_in_biotype is None:
            # see https://www.ensembl.org/info/genome/genebuild/biotypes.html
            self.__protein_coding_labels_in_biotype = [
                "nonsense_mediated_decay",
                "protein_coding",
            ]

        # creates transcript information for each transcript. This populates self.transcripts_info
        self.__gen_transcripts_info()

        # creates exon-intron map for each transcript. This populates self.exon_intron_maps
        self.__create_exon_intron_map()

    def access_mode(self) -> str:
        """Chromosome data access mode."""
        return 'remote' if self.chrm_fasta_file is None else 'local'

    # def _extract_sequence(self, start_p: int, end_p: int, rev: bool = False) -> str:
    #     """Extracts sequence from teh chromosome.

    #     Args:
    #     start_p (int): 1-based start position.
    #     end_p (int): 1-based end position.
    #     rev (bool, optional): True to extract the reveresed-complement sequence 
    #                          (i.e. the sequence in the negative strand). Defaults to False.
    # Returns:
    #     str: The extracted sequence.
    #     """
    #     # if self.chrm_fasta_file is None:
    #     if self.access_mode() == 'remote':
    #         # use remote chromosome file via Ensembl REST API
    #         return extract_chromosome_seq(self.chrm, start_p, end_p, rev=rev, species=self.species, rest_assembly=self._rest_assembly)
    #     # use local fasta file
    #     return extract_fasta_seq(self.chrm_fasta_file, start_p, end_p, rev=rev)


    def _extract_sequence(self, start_p: int, end_p: int, rev: bool = False) -> str:
        """Extracts sequence from the chromosome.

        Args:
        start_p (int): 1-based start position.
        end_p (int): 1-based end position.
        rev (bool, optional): True to extract the reveresed-complement sequence 
                             (i.e. the sequence in the negative strand). Defaults to False.
    Returns:
        str: The extracted sequence.
        """
        source_info: Path | str = self.chrm if self.access_mode() == 'remote' else self.chrm_fasta_file
        return fetch_seq(source_info, start_p, end_p, rev=rev, species=self.species, rest_assembly=self._rest_assembly)

    def __len__(self) -> int:
        "Number of transcripts."
        return len(self.transcripts)

    def __check_transcript_exists(self, transcript_id: str) -> bool:
        if transcript_id in self.transcripts_info:
            return True
        if self.verbose:
            print(f"{transcript_id=} not in {self.gene} transcript list !!.")
        return False

    def is_protein_coding_transcript(self, transcript_id: str) -> bool:
        """Returns True [False] is transcript_id is [is not] a protein coding transcript."""
        if not self.__check_transcript_exists(transcript_id):
            return False
        try:
            return (
                self.transcripts_info[transcript_id]["transcript_biotype"]
                in self.__protein_coding_labels_in_biotype
            )
        except KeyError:
            print(f"No biotype information for {transcript_id}. Assuming non-coding transcript !!")
            return False

    def __gen_transcripts_info(self):
        """
        Reads the data from self.transcripts and creates a dictionary with transcript ID
        as keys and a dictionary as value, where the keys and values are taken from the self.transcripts.

        Note that for genes encoded on the negative start, all the chromosome position
        lists (e.g., exon_start, exon_end, CDS_start, CDS_end, CDS_phase, ..) contain
        decreasing values, otherwise (genes encoded on the positive strand) these
        lists contain increasing values.
        """
        key_names = [
            "transcript_name",
            "transcript_v_id",
            "transcript_ccdsid",
            "transcript_start",
            "transcript_end",
            "transcript_ver",
            "transcript_type",
            "transcript_biotype",
            "exon_start",
            "exon_end",
            "exon_id",
            "exon_phase",
            "exon_end_phase",
            "exon_ver",
            "CDS_start",
            "CDS_end",
            "CDS_phase",
            "protein_id",
            "5UTR_start",
            "5UTR_end",
            "3UTR_start",
            "3UTR_end",
        ]

        m = -1 if self.rev else 1
        self.transcripts_info = {}
        for x in self.transcripts:
            # compute ORF start and end from the beginning of the transcript
            # (i.e., from the beginning of the first exon)
            if x["transcript_biotype"] in self.__protein_coding_labels_in_biotype:
                start_codon_first_bp_pos = x["CDS_start"][0] + m * x["CDS_phase"][0]
                end_codon_last_bp_pos = x["CDS_end"][-1]

                if (
                    tmp := pos2seg_info(
                        start_codon_first_bp_pos,
                        x["exon_start"],
                        x["exon_end"],
                        self.rev,
                    )
                ) is None:
                    orf_start = -1
                else:
                    # the 0-based offset of the start_codon_first_bp_pos from the beginning of the mRNA (i.e.,
                    # from the beginning of a sequence containing all exons)
                    orf_start = tmp[2]

                if (
                    tmp := pos2seg_info(
                        end_codon_last_bp_pos, x["exon_start"], x["exon_end"], self.rev
                    )
                ) is None:
                    orf_end = -1
                else:
                    # the 0-based offset of the end_codon_last_bp_pos from the beginning of the mRNA (i.e.,
                    # from the beginning of a sequence containing all exons)
                    orf_end = tmp[2]

                orf_d = {"ORF_start_offset": orf_start, "ORF_end_offset": orf_end}
                # (1-based) chromosome positions
                orf_chrm_d = {
                    "ORF_start_chrm_pos": start_codon_first_bp_pos,
                    "ORF_end_chrm_pos": end_codon_last_bp_pos,
                }
            else:
                orf_d = {"ORF_start_offset": None, "ORF_end_offset": None}
                orf_chrm_d = {"ORF_start_chrm_pos": None, "ORF_end_chrm_pos": None}

            self.transcripts_info[x["transcript_id"]] = (
                {k: x[k] for k in key_names} | orf_d | orf_chrm_d
            )

    def __create_exon_intron_map(self):
        self.exon_intron_maps = {
            transcript: self.__create_transcript_exon_intron_map(transcript)
            for transcript in self.transcripts_info.keys()
        }

    def __create_transcript_exon_intron_map(self, transcript_id: str) -> pd.DataFrame | None:
        """
        Created an exon-intron map of a transcript
        """
        if not self.__check_transcript_exists(transcript_id):
            return None

        t_info = self.transcripts_info[transcript_id]

        exons = list(zip(t_info["exon_start"], t_info["exon_end"]))

        delta = -1 if self.rev else 1
        intr_start_list = [x + delta for x in t_info["exon_end"][:-1]]
        intr_end_list = [x - delta for x in t_info["exon_start"][1:]]
        introns = list(zip(intr_start_list, intr_end_list))

        ex_lbl = [f"Exon {x}" for x in range(1, len(exons) + 1)]
        in_lbl = [f"Intron {x}" for x in range(1, len(introns) + 1)]

        name = [y for x in zip(ex_lbl, in_lbl) for y in x] + [ex_lbl[-1]]
        region = [y for x in zip(exons, introns) for y in x] + [exons[-1]]

        return pd.DataFrame(
            {
                "name": name,
                "region": region,
                "region_size": [np.abs(r[1] - r[0]) + 1 for r in region],
            }
        )

    def info(self):
        """Prints gene information."""
        print(
            f"Gene={self.gene}, Gene name={self.gene_name}, Gene ID={self.gene_ID}\nSource={self.source}\nPositive strand={not self.rev}."
        )
        print(f"Type={self.gene_type}, version={self.gene_ver}.")
        print(f"Description={self.gene_desc}")
        print(f"Gene region=chr{self.chrm}:{self.gene_start:,}-{self.gene_end:,}.")
        print(f"{len(self.transcripts_info)} transcripts.")

    def transcript_info(self, transcript_id: str, verbose: bool = False):
        """Prints transcript information."""
        if self.__check_transcript_exists(transcript_id):
            t = self.transcripts_info[transcript_id]
            print(f"{transcript_id} information:")
            print(
                f"Name={t['transcript_name']}, ID={t['transcript_v_id']}, version={t['transcript_ver']}"
            )
            print(f"Type={t['transcript_type']}, biotype={t['transcript_biotype']}")
            print(f"Start={t['transcript_start']:,}, end={t['transcript_end']:,}")
            print(f"protein ID = {t['protein_id']}")
            if t["transcript_ccdsid"] is not None:
                print(f"ccdsid={t['transcript_ccdsid']}")
            if self.is_protein_coding_transcript(transcript_id):
                print(
                    f"ORF start offset = {t['ORF_start_offset']:,}, ORF end offset = {t['ORF_end_offset']:,} (0-based offsets from the beginning of the RNA)"
                )
                print(
                    f"ORF start chromosome position = {t['ORF_start_chrm_pos']:,}, ORF end chromosome position = {t['ORF_end_chrm_pos']:,}"
                )
            else:
                print("No ORF.")

            if verbose:
                print("Exon regions:")
                for i, (start, end, eid, ep, eep) in enumerate(
                    zip(
                        t["exon_start"],
                        t["exon_end"],
                        t["exon_id"],
                        t["exon_phase"],
                        t["exon_end_phase"],
                    ),
                    start=1,
                ):
                    print(
                        f"\t{i:2}. {start:,} - {end:,}, ID={eid}, phase={ep}, end phase={eep}"
                    )

                print("CDS regions:")
                for i, (start, end, p) in enumerate(
                    zip(t["CDS_start"], t["CDS_end"], t["CDS_phase"]), start=1
                ):
                    print(f"\t{i:2}. {start:,} - {end:,}, phase={p}")

                print("5UTR regions:")
                for i, (start, end) in enumerate(
                    zip(t["5UTR_start"], t["5UTR_end"]), start=1
                ):
                    print(f"\t{i:2}. {start:,} - {end:,}")

                print("3UTR regions:")
                for i, (start, end) in enumerate(
                    zip(t["3UTR_start"], t["3UTR_end"]), start=1
                ):
                    print(f"\t{i:2}. {start:,} - {end:,}")
            else:
                print(f"{len(t['exon_start'])} exons.")

    def show_all_transcript_IDs(self):
        """Prints transcript IDs."""
        print(f"{self.gene} transcript IDs:")
        for i, t in enumerate(self.transcripts_info.keys(), start=1):
            print(f"{i:2}. {t}")

    def exon_intron_map(self, transcript_id: str) -> pd.DataFrame | None:
        """Returns exon and intron map. This is basically the primary transcript map."""
        if not self.__check_transcript_exists(transcript_id):
            return None

        return self.exon_intron_maps[transcript_id]

    def exon_intron_map_to_excel(self, transcript_id: str, excel_file: str, usr_desc: dict = None, desc_sheet_pos: int = None):
        """Writes exon map to excel file"""
        if self.__check_transcript_exists(transcript_id):
            dfs_to_excel_file([self.exon_intron_maps[transcript_id]], excel_file, ["Exon-Intron Map"])

            # adding description
            desc = {
                "Description": f"Exon-Intron map of {self.gene} gene, {transcript_id} transcript.",
                "Gene info": f"chr{self.chrm}:{self.gene_start:,}-{self.gene_end:,}, {'negative' if self.rev else 'positive'} strand.",
                "Source": self.source
            }
            if usr_desc is not None:
                desc |= usr_desc
            create_excel_description_sheet(excel_file, desc, sheet_name="Description", sheet_pos=desc_sheet_pos)

    def exon_intron_map_to_html(self, transcript_id: str, html_file: str):
        """Writes exon map to HTML file"""
        if self.__check_transcript_exists(transcript_id):
            self.exon_intron_maps[transcript_id].to_html(html_file, index=False)

    def exon_intron_map_to_csv(self, transcript_id: str, csv_file: str, sep=','):
        """Writes exon map to CSV file"""
        if self.__check_transcript_exists(transcript_id):
            self.exon_intron_maps[transcript_id].to_csv(csv_file, sep=sep, index=False)

    def transcript_boundaries(self, transcript_id: str) -> tuple[int, int] | None:
        """
        Returns (transcript_start_site, transcript_end_site) of a (primary) transcript.
        Note that for a gene encoded on the negative strand,
        transcript_start_site>transcript_end_site.
        """
        if not self.__check_transcript_exists(transcript_id):
            return None

        t_info = self.transcripts_info[transcript_id]
        return t_info["transcript_start"], t_info["transcript_end"]

    def seq(self, transcript_id: str) -> str | None:
        """
        Generates the (primary) transcript sequence (i.e. pre-rna sequence).
        """
        if not self.__check_transcript_exists(transcript_id):
            return None

        ts, te = self.transcript_boundaries(transcript_id)

        # for rev=True, transc_start [transc_end] is the end [start] of the transcript on the negative strand
        (transc_start, transc_end) = (te, ts) if self.rev else (ts, te)

        #seq = extract_fasta_seq(chrm_path, transc_start, transc_end, rev=False)
        seq = self._extract_sequence(transc_start, transc_end, rev=False)
        return tran.reverse_complement(seq) if self.rev else seq

    def exon_intron_seq(self, name: str, number: int, transcript_id: str) -> tuple[str, str] | None:
        """
        Returns the sequence corresponding to exon (if name=='Exon') or intron (if name=='Intron') number <number>,
        (1 for first exon/intron), and the lable which contains <name>_<number>:<chromosome>:<start 1-based position>:<end 1-based position>.

        Note: for genes encoded on the negative DNA strand, the returned Intron/Exon sequence
        is already the DNA reveresed-complemented sequence (i.e., the sequence as appears in the RNA).
        """
        if not self.__check_transcript_exists(transcript_id):
            return None

        if name not in ['Intron', 'Exon']:
            print("name input must be either 'Intron' or 'Exon' !!")
            return None

        if (df_q := self.exon_intron_maps[transcript_id].query(f"name == '{name} {number}'")).empty:
            print(f"No {name} #{number} found in exon/intron table of {transcript_id=} !!")
            return None

        region = df_q.iloc[0]['region']
        #seq = extract_fasta_seq(chrm_path, np.min(region), np.max(region), rev=False)
        seq = self._extract_sequence(np.min(region), np.max(region), rev=False)
        if self.rev:
            seq = tran.reverse_complement(seq)

        return seq, f"{name}_{number}:chr{self.chrm}:{region[0]}:{region[1]}"

    def modified_transcript(self, exon_list: list, intron_list: list, transcript_id: str) -> str | None:
        """
        Generates a modified (primary) transcript.

        exon_list - a list of (1-based) exon numbers indicating exons to retain in the emodified transcript. The list order is irrelevant.
        intron_list - a list of (1-based) intron numbers indicating introns to retain in the modified transcript. The list order is irrelevant.

        A transcript contains the following exon/intron arrangment: E_1 I_1 E_2 I_2 ... I_(n-1), E_n,
        where E_i [I_i] is the i'th Exon [Intron].
        In the returned modified transcript, E_j exon and/or I_k intron are missing if j \notin exon_list
        and/or k \notin intron_list.

        For example, given a transcript that contains 3 exons: E_1 I_1 E_2 I_2 E_3,
        then for exon_list=[1,2,3], intron_list = [2], the modified transcript is: E_1 E_2 I_2 E_3.
        As another example, exon_list=[1,2,..n] and intron_list=[] yields the mRNA.
        """
        if not self.__check_transcript_exists(transcript_id):
            return None

        num_exons = len(self.transcripts_info[transcript_id]['exon_start'])
        # verify:
        if exon_list and (min(exon_list) < 1 or max(exon_list) > num_exons):
            print(f"exon list out of range (values should be between 1 and {num_exons}) !!")
            return None

        if intron_list and (min(intron_list) < 1 or max(intron_list) > (num_exons - 1)):
            print(f"Intron list out of range (values should be between 1 and {num_exons -1}) !!")
            return None

        list_info = {'Exon': exon_list, 'Intron': intron_list}
        return ''.join(
            [
                self.exon_intron_seq(name, index, transcript_id)[0] if index in list_info[name] else ''
                for index in range(1, num_exons+1)
                for name in ['Exon', 'Intron']
            ]
        )

    def chrm_pos_info(self, transcript_id: str, chrm_pos: int) -> dict | None:
        """
        Given a 1-based chromosome position, the function returns
        the dictionary {'region': <the name of the region containing chrm_pos position>,
        'region_pos': <the 1-based position of chrm_pos in the region>,
        'bp': <the basepair value at the strand in which th gene is encoded>}
        """
        if not self.__check_transcript_exists(transcript_id):
            return None

        df = self.exon_intron_maps[transcript_id]
        (low_indx, high_indx) = (1, 0) if self.rev else (0, 1)
        df_pos = df.loc[
            df.apply(
                lambda x: x["region"][low_indx] <= chrm_pos <= x["region"][high_indx],
                axis=1,
            )
        ]

        if not df_pos.empty:
            return {
                "region": "_".join(df_pos.iloc[0]["name"].split()),
                "region_pos": abs(chrm_pos - df_pos.iloc[0]["region"][0]) + 1,
                # added on 12/9/24 - distance from boundary. 0 implies right at the region boundary
                "dist_from_region_boundary": (abs(chrm_pos - df_pos.iloc[0]["region"][0]), abs(chrm_pos - df_pos.iloc[0]["region"][1])),
            } | {"bp": self._extract_sequence(chrm_pos, chrm_pos, rev=self.rev)}
        print(f"Chromosome position {chrm_pos:,} out of {transcript_id} bound.")
        return None

    def chrm_pos2rna_pos(self, transcript_id: str, chrm_pos: int) -> tuple[int, str] | None:
        """
        Given a 1-based chromosome position, then if the position is within the exons, the function returns
        the corresponding 1-based position within the RNA sequence (i.e. the corresponding 1-based position
        within a sequence that contains all exons), and the corresponding
        bp value. If the chromosome position is outside the exons, the function returns None.
        """
        if not self.__check_transcript_exists(transcript_id):
            return None

        t_info = self.transcripts_info[transcript_id]
        if (tmp := pos2seg_info(chrm_pos, t_info["exon_start"], t_info["exon_end"], self.rev)) is None:
            return None
        return tmp[-1] + 1, self._extract_sequence(chrm_pos, chrm_pos, rev=self.rev)

    def rna_pos2chrm_pos(self, transcript_id: str, pos: int) -> tuple[int, str] | None:
        """
        Given a 1-base RNA position pos, the function returns its corresponding chromosome coordinate.
        If pos is outsize the RNA size, the returned valiue is -1.
        """
        if not self.__check_transcript_exists(transcript_id):
            return None

        t_info = self.transcripts_info[transcript_id]
        exon_start, exon_end = t_info["exon_start"], t_info["exon_end"]

        exon_sizes = [abs(x - y) + 1 for x, y in zip(exon_start, exon_end)]
        csum = np.insert(np.cumsum(exon_sizes), 0, 0)  # [0, len(exon1), len(exon1)+len(exon2), ..]

        if pos < 1 or pos > csum[-1]:
            return -1, ''  # out of range

        # index of exon containing RNA position pos
        #exon_index =  np.array([csum[i] < pos <= csum[i+1] for i in range(csum.shape[0]-1)]).nonzero()[0][0]  # original (result in a lint error)
        exon_index = np.asarray([csum[i] < pos <= csum[i+1] for i in range(csum.shape[0]-1)]).nonzero()[0][0]  # YZ, 11/14/24
        # 0-based offset of pos from the beginning of the corresponding exon
        # (recall that csum[exon_index] is the cumulative sum up to the previous exon, as csum[0]=0, thus (csum[exon_index] + 1) is the beginning of the corresponding exon)
        offset = pos - (csum[exon_index] + 1)
        chrm_pos = exon_start[exon_index] - offset if self.rev else exon_start[exon_index] + offset

        # return (int(chrm_pos), extract_fasta_seq(chrm_path, chrm_pos, chrm_pos, rev=self.rev)) if chrm_path is not None else int(chrm_pos)
        return (int(chrm_pos), self._extract_sequence(chrm_pos, chrm_pos, rev=self.rev))

    def __repr__(self):
        return f"Transcript_gff3_cls(gene={self.gene}, gff3_source=gff3_file|(gff3_df, gff3_df_gene_type))"

    def __str__(self):
        return "Transcript data of all transcripts of a gene based on Ensembl GFF3 annotation file.\n"


@dataclass
class Gene_gff3_cls(Transcript_gff3_cls):
    """
    RNA class derived from Transcript_gff3_cls class.

    This class adds few other query methods, the RNA (mRNA), and in case of a protein coding transcript also the
    ORF, AA, 3'UTR, and 5'UTR sequences.
    """
    exon_maps: dict = None  # keys are transcripts, values are exon maps in df format

    def __post_init__(self):
        super().__post_init__()

        # create exon map for each transcript. This populates self.exon_maps
        self.__create_exon_map()

    def __check_transcript_exists(self, transcript_id: str) -> bool:
        return self._Transcript_gff3_cls__check_transcript_exists(transcript_id)

    def get_sequence_from_start_end_segments(self, start: list[int], end: list[int], rev: bool, offset: int = 0) -> str:
        """
        Get the chromosome sequence, at the corresponding strand (based on rev input) from
        the concatenation of segments defined by start and end list.

        rev == True  implies start[i]>end[i], start[j]<start[i], end[j]<end[i], for j>i.
        rev == False implies start[i]<end[i], start[j]>start[i], end[j]>end[i], for j>1.
        """
        # convert to increasing lists
        a_start, a_end = (end[::-1], start[::-1]) if rev else (start, end)
        seq = "".join(
            [
                #extract_fasta_seq(chrm_path, a_s, a_e, rev=False)
                self._extract_sequence(a_s, a_e, rev=False)
                for a_s, a_e in zip(a_start, a_end)
            ]
        )
        return tran.reverse_complement(seq)[offset:] if rev else seq[offset:]

    def show_exon_map(self, transcript_id: str):
        """Prints exon map."""
        if self.__check_transcript_exists(transcript_id):
            print(self.exon_maps[transcript_id])

    def exon_map(self, transcript_id: str) -> pd.DataFrame | None:
        """Returns exon map."""
        if self.__check_transcript_exists(transcript_id):
            return self.exon_maps[transcript_id]
        return None

    def exon_map_to_excel(self, transcript_id, excel_file: str, usr_desc: dict | None = None, desc_sheet_pos: int | None = None):
        """Writes exon map to excel file"""
        if self.__check_transcript_exists(transcript_id):
            dfs_to_excel_file([self.exon_maps[transcript_id]], excel_file, ["Exon Map"])

            # adding description
            desc = {
                "Description": f"Exon map of {self.gene} gene, {transcript_id} transcript.",
                "Gene info": f"chr{self.chrm}:{self.gene_start:,}-{self.gene_end:,}, {'negative' if self.rev else 'positive'} strand.",
                "Source": self.source
            }
            if usr_desc is not None:
                desc |= usr_desc
            create_excel_description_sheet(excel_file, desc, sheet_name="Description", sheet_pos=desc_sheet_pos)

    def exon_map_to_html(self, transcript_id, html_file: str):
        """Writes exon map to HTML file"""
        if self.__check_transcript_exists(transcript_id):
            self.exon_maps[transcript_id].to_html(html_file, index=False)

    def exon_map_to_csv(self, transcript_id, csv_file: str, sep: str = ','):
        """Writes exon map to CSV file"""
        if self.__check_transcript_exists(transcript_id):
            self.exon_maps[transcript_id].to_csv(csv_file, sep=sep, index=False)

    def __create_exon_map(self):
        self.exon_maps = {
            transcript: self.__create_transcript_exon_map(transcript)
            for transcript in self.transcripts_info.keys()
        }

    def __create_transcript_exon_map(self, transcript_id: str) -> pd.DataFrame | None:
        """
        Creates an exon map (mRNA map) of a protein-coding transcript.

        Given a transcript ID of a gene, we generate a map between
        exons and the corresponding exon's NT regions, and AA and ORF's NT regions.

        Let NT_ORF be the ORF's NT sequence, and let AA_ORF be the corresponding amino-acid sequence
        (i.e. AA_ORF=translate(NT_OTF)).
        For example, these can be obtained by executing:
                        AA_ORF = self.AA(transcript_id, chrm_path),
                        NT_ORF = self.ORF(transcript_id, chrm_path).

        We return a dataframe with the following columns:
        1. 'exon_number' - 1-based exon number.
        2. 'exon_region' - (exon_start, exon_end)
        3. 'exon_size' - number of NTs in the exon
        4. 'mRNA_NT_region' - [mRNA_NT_index_in_the_beginning_of_the_exon, mRNA_NT_index_in_the_end_of_the_exon]
        5. 'ORF_AA_region' - [AA_index_in_the_beginning_of_the_exon, AA_index_in_the_end_of_the_exon],
            where AA_index is the index of the AA in AA_ORF (1 indicates the first AA in AA_ORF).
        6. 'ORF_NT_region' - [ORF_NT_index_in_the_beginning_of_the_exon, ORF_NT_index_in_the_end_of_the_exon],
            where ORF_NT_index is the index of the NT in NT_ORF (1 indicates the first NT in NT_ORF).
            This includes the stop codon.
        7. 'next_exon_frame_alignment' - <a value of 0, 1 or 2>. This is the number of NTs required in the beginning
            of the next exon to complete the last AA (last codon) in the current exon.
            For example, a value of 1 means that the first NT of the next exon is the last NT
            of the last codon in the current exon, and the frame alignment in the next exon starts from the second NT.
        8. 'exon_ID' - the Ensembl exon ID.

        Note that if the gene is encoded on the negative DNA strand, then the first exon of the transcript is
        the exon with the largest coordinates, and in this case exon_start>exon_end.
        """
        if not self.is_protein_coding_transcript(transcript_id):
            # print(f"INFO: {transcript_id} is not a protein coding transcript.")
            return None

        t_info = self.transcripts_info[transcript_id]
        exon_info = [
            {
                "exon": (ex_start, ex_end),
                "exon_size": ex_start - ex_end + 1
                if self.rev
                else ex_end - ex_start + 1,
            }
            for ex_start, ex_end in zip(t_info["exon_start"], t_info["exon_end"])
        ]

        """
        The value of exon_info[i]['exon'] is a tuple with the exon start and exon end.
        Note that if self.rev==True, then exon_start>exon_end

        Also, note that the first dictionary in the exon_info list corresponds to the first exon (which,
        in the case of rev=True, corresponds to the exon with the largest coordinates), and the last
        dictionary in the exon_info list corresponds to the last exon (which, in case of rev=True,
        corresponds to the exon with the smallest coordinates).

        Gather AA and NT regions per exon.

        Note that in general, the ORF is generated as follows:

        For self.rev==False: we collect all exons (starting from the smallest coordinates) to one string of nucleotides,
        referred to as seq, and then we perform:
        1. orf = seq[cds_start_nt:cds_end_nt + 1]

        For self.rev==True: we collect all exons (again starting from the smallest coordinates) to one
        string of nucleotides, referred to as tmp_seq, and then we perform:
        1. seq = reverse_complement(tmp_seq)
        2. orf = seq[cds_start:cds_end+1]

        Thus, in both cases, the seq variable corresponds to the order of exons in exon_info.

        The assumption here is that CDS_start and CDS_end completely overlap with a corresponding
        exon_start_and exon_end, OTHER than (possibly) the first and last pairs of
        CDS_start and CDS_end.
        """
        m = -1 if self.rev else 1

        # the chromosome position of the first bp of the start codon
        start_codon_first_bp_pos = t_info["CDS_start"][0] + m * t_info["CDS_phase"][0]
        # finding in which exon start_codon_first_bp_pos is located, and its offset from the beginning of that exon
        if (
            tmp := pos2seg_info(
                start_codon_first_bp_pos,
                t_info["exon_start"],
                t_info["exon_end"],
                self.rev,
            )
        ) is None:
            if self.verbose:
                print(f"{start_codon_first_bp_pos=} can not be found in {transcript_id=} exons !!")
            return None
        exon_containing_start_codon = tmp[0]

        # the chromosome position of the third bp of the last codon
        end_codon_last_bp_pos = t_info["CDS_end"][-1]
        if (
            tmp := pos2seg_info(
                end_codon_last_bp_pos,
                t_info["exon_start"],
                t_info["exon_end"],
                self.rev,
            )
        ) is None:
            if self.verbose:
                print(f"{end_codon_last_bp_pos=} can not be found in {transcript_id=} exons !!")
            return None
        exon_containing_last_codon = tmp[0]

        exon_aa_info, aa_region, nt_region = [], [], []
        cds_indx = 0
        for ex_indx, ex in enumerate(exon_info):
            if ex_indx < exon_containing_start_codon:
                # this exon is upstream of the start codon
                aa_region = [0, 0]
                nt_region = [0, 0]
                last_aa_nt_required = (
                    0  # number of NTs required in the next exon to complete the last AA
                )
            elif ex_indx > exon_containing_last_codon:
                # this exon is downstream of the stop codon
                aa_region = [0, 0]
                nt_region = [0, 0]
                last_aa_nt_required = (
                    0  # number of NTs required in the next exon to complete the last AA
                )
            else:
                # this exon contains CDS (ORF) bps
                if cds_indx == 0:
                    # this is the first exon containing CDS data (ORF)
                    aa_nts_in_exon = (
                        m * (t_info["CDS_end"][cds_indx] - start_codon_first_bp_pos) + 1
                    )
                    aa_aas_in_exon = (
                        aa_nts_in_exon - 1
                    ) // 3 + 1  # number of aa in exon
                    nt_region = [1, aa_nts_in_exon]
                    aa_region = [1, aa_aas_in_exon]
                else:
                    # this is the x'th (x>1) exon containing CDS data
                    aa_nts_in_exon = (
                        m
                        * (t_info["CDS_end"][cds_indx] - t_info["CDS_start"][cds_indx])
                        + 1
                    )
                    nt_region = [nt_region[1] + 1, nt_region[1] + aa_nts_in_exon]

                    # the CDS phase contains the number of bps in this CDS segment that are part
                    # of the last codon in the previous segment.
                    aa_aas_in_exon = (
                        (aa_nts_in_exon - t_info["CDS_phase"][cds_indx]) - 1
                    ) // 3 + 1  # number of aa
                    offset = int(
                        t_info["CDS_phase"][cds_indx] == 0
                    )  # 1 if CDS_phase==0, otherwise 0
                    aa_region = [aa_region[1] + offset, aa_region[1] + aa_aas_in_exon]

                # set last_aa_nt_required value
                if cds_indx == (len(t_info["CDS_start"]) - 1):
                    # cds_indx is the last CDS segment index
                    last_aa_nt_required = 0  # recall that this is for the next exon
                else:
                    # cds_indx is not the last CDS segment index
                    # use the phase of the next CDS segment
                    last_aa_nt_required = t_info["CDS_phase"][cds_indx + 1]

                cds_indx += 1

            exon_aa_info.append(
                {
                    "exon_region": ex["exon"],
                    "exon_size": ex["exon_size"],
                    "ORF_NT_region": nt_region,
                    "ORF_AA_region": aa_region,
                    "next_exon_frame_alignment": last_aa_nt_required,
                }
            )

        df = (
            pd.DataFrame(exon_aa_info, index=np.arange(1, len(exon_aa_info) + 1))
            .reset_index()
            .rename(columns={"index": "exon_number"})
        )

        # add mRNA NT region per exon
        end = df["exon_size"].cumsum()
        start = [1] + [x + 1 for x in end[:-1]]
        df.insert(3, "mRNA_NT_region", list(map(list, zip(start, end))))

        # add exon ID
        df.insert(df.shape[1], 'exon_ID', t_info['exon_id'])

        '''
        The maximal value in ORF_AA_region[1] is y:=Z+1, where Z is the number of amino-acids in the ORF. The reason
        for the +1 is due to the stop codon. We change it now to be Z instead of Z+1.
        These are the possible cases of ORF_AA_region where ORF_AA_region[1]==y:
        1. ORF_AA_region = [x, y], x < y. This needs to change to [x, y-1].
        2. ORF_AA_region = [y, y]. This needs to change to [0, 0], since this exon contains only the stop codon from the ORF
        (either 1, 2 or 3 NTs of the stop codon, but none of the codons preceding the stop codon),
        and the NTs that are downstream of the ORF.
        '''
        #y = df.apply(lambda x: x['ORF_AA_region'][1], axis=1).max().max()
        y = df.apply(lambda x: x['ORF_AA_region'][1], axis=1).max()  # changed from above on 6/26/23.
        # change ORF_AA_region in rows containing ORF_AA_region[1] == y
        for r_lbl, r in df.loc[df.apply(lambda x: x['ORF_AA_region'][1] == y, axis=1)].iterrows():
            df.at[r_lbl, 'ORF_AA_region'] = [0, 0] if r['ORF_AA_region'][0] == y else [r['ORF_AA_region'][0], y - 1]

        return df

    def start_and_stop_codons_pos(self, transcript_id: str) -> tuple[tuple[int, int], tuple[int, int]] | None:
        """
        Returns the tuple (a,b), where:
        1. a = (start_codon_first_bp_chromosome_position, start_codon_first_bp_rna_position)
        2. b = (stop_codon_first_bp_chromosome_position, stop_codon_first_bp_rna_position)
        All positions are 1-based.
        """
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        t_info = self.transcripts_info[transcript_id]
        m = -1 if self.rev else 1

        # ORF_start_offset and ORF_end_offset are 0-based offsets from the beginning of the RNA.
        # ORF_start_chrm_pos and ORF_end_chrm_pos are 1-based chromosome positions.
        start_codon_poss = (
            t_info["ORF_start_chrm_pos"],
            t_info["ORF_start_offset"] + 1,  # output is 1-based
        )
        stop_codon_poss = (
            t_info["ORF_end_chrm_pos"] - m * 2,  # t_info['ORF_end_chrm_pos'] is the position of the last bp of the stop codon
            (t_info["ORF_end_offset"] + 1) - 2,  # ORF_end_offset is the last bp offset from the beginning of the RNA
        )
        return start_codon_poss, stop_codon_poss

    def orf_start_and_end_exon_info(self, transcript_id: str) -> dict | None:
        """
        Computes the exon number and offset of the first ORF NT and the last ORF NT.
        Returns the dictionary:
            {'start': start_info, 'end': end_info},
        where *_info is the dictionary:
        {'exon_number': <the 1-based exon number containing the first/last ORF NT>,
         'exon_offset': <the 0-based position of the first/last ORF NT in the corresponding exon>}
        """
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        t_info = self.transcripts_info[transcript_id]
        if (tmp_start := pos2seg_info(t_info['ORF_start_chrm_pos'], t_info['exon_start'], t_info['exon_end'], self.rev)) is None:
            print(f"Can not find ORF_start_chrm_pos={t_info['ORF_start_chrm_pos']:,} in exons !!")
            return None
        if (tmp_end := pos2seg_info(t_info['ORF_end_chrm_pos'], t_info['exon_start'], t_info['exon_end'], self.rev)) is None:
            print(f"Can not find ORF_end_chrm_pos={t_info['ORF_end_chrm_pos']:,}, in exons !!")
            return None

        return {
            'start': {'exon_number': tmp_start[0] + 1, 'exon_offset': tmp_start[1]},  # exon_number is 1-based. exon_offset is 0-based.
            'end' : {'exon_number': tmp_end[0] + 1, 'exon_offset': tmp_end[1]},
        }

    def rna(self, transcript_id: str) -> str | None:
        """Generates the RNA nucleotide sequence of a transcript."""
        if not self.__check_transcript_exists(transcript_id):
            return None

        t_info = self.transcripts_info[transcript_id]
        # return get_sequence_from_start_end_segments(t_info["exon_start"], t_info["exon_end"], self.rev, chrm_path)
        return self.get_sequence_from_start_end_segments(t_info["exon_start"], t_info["exon_end"], self.rev)
        # return tran.DNA2RNA(
        #    get_sequence_from_start_end_segments(t_info['exon_start'], t_info['exon_end'], self.rev, chrm_path)
        # )

    def ORF(self, transcript_id: str) -> str | None:
        """Generates the ORF nucleotide sequence of a transcript."""
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        t_info = self.transcripts_info[transcript_id]
        orf_offset = t_info["CDS_phase"][0]
        # return get_sequence_from_start_end_segments(
        #     t_info["CDS_start"],
        #     t_info["CDS_end"],
        #     self.rev,
        #     chrm_path,
        #     offset=orf_offset,
        # )
        return self.get_sequence_from_start_end_segments(
            t_info["CDS_start"],
            t_info["CDS_end"],
            self.rev,
            offset=orf_offset,
        )
        # return tran.DNA2RNA(
        #    get_sequence_from_start_end_segments(t_info['CDS_start'], t_info['CDS_end'], self.rev, chrm_path,
        #                                         offset=orf_offset)
        # )

    def AA(self, transcript_id: str) -> str | None:
        """Generates the amino-acid sequence of a transcript."""
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        orf = self.ORF(transcript_id)[:-1]  # remove the stop codon
        return tran.convert_aa_32IUPAC(tran.translate(orf))

    def UTR5(self, transcript_id: str) -> str | None:
        """Generates the 5'UTR nucleotide sequence of a transcript."""
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        t_info = self.transcripts_info[transcript_id]
        # return get_sequence_from_start_end_segments(t_info["5UTR_start"], t_info["5UTR_end"], self.rev, chrm_path)
        return self.get_sequence_from_start_end_segments(t_info["5UTR_start"], t_info["5UTR_end"], self.rev)
        # return tran.DNA2RNA(
        #    get_sequence_from_start_end_segments(t_info['5UTR_start'], t_info['5UTR_end'], self.rev, chrm_path)
        # )

    def UTR3(self, transcript_id: str) -> str | None:
        """Generates the 3'UTR nucleotide sequence of a transcript."""
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        t_info = self.transcripts_info[transcript_id]
        # return get_sequence_from_start_end_segments(t_info["3UTR_start"], t_info["3UTR_end"], self.rev, chrm_path)
        return self.get_sequence_from_start_end_segments(t_info["3UTR_start"], t_info["3UTR_end"], self.rev)
        # return tran.DNA2RNA(
        #    get_sequence_from_start_end_segments(t_info['3UTR_start'], t_info['3UTR_end'], self.rev, chrm_path)
        # )

    def exon_nt_segment(self, transcript_id: str, exon_number: int, nt_number: int) -> tuple[str, int] | None:
        """
        Given an exon number and a nucleotide number within the exon (all 1-based), the function returns
        the segment type that the nucleotide belongs to (i.e. '5UTR', 'ORF', or '3UTR'),
        and the 1-based position of the nucleotide within that segment.
        """
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        df = self.exon_map(transcript_id)
        if exon_number < 0 or exon_number > df.shape[0]:
            print(
                f"{exon_number=} out of boundary ({transcript_id} contains {df.shape[0]} exons) !!"
            )
            return None

        exon_info = df.query(f"exon_number=={exon_number}").iloc[0]

        if nt_number < 0 or nt_number > exon_info["exon_size"]:
            print(
                f"{nt_number=} out of boundary (exon #{exon_number} contains {exon_info['exon_size']} bps) !!"
            )
            return None

        # the 1-based position of nt_number in the mRNA
        nt_pos_in_mrna = exon_info["mRNA_NT_region"][0] + nt_number - 1

        t_info = self.transcripts_info[transcript_id]
        # the 1-based position of ORF start and ORF end in the mRNA
        # (ORF_start_offset and ORF_end_offset are 0-based offset from the begiinning of the mRNA)
        orf_start_in_mrna, orf_end_in_mrna = (
            t_info["ORF_start_offset"] + 1,
            t_info["ORF_end_offset"] + 1,
        )

        if nt_pos_in_mrna < orf_start_in_mrna:
            segment = "5UTR"
            pos_in_segment = nt_pos_in_mrna
        elif orf_start_in_mrna <= nt_pos_in_mrna <= orf_end_in_mrna:
            segment = "ORF"
            pos_in_segment = nt_pos_in_mrna - orf_start_in_mrna + 1
        else:
            segment = "3UTR"
            pos_in_segment = nt_pos_in_mrna - orf_end_in_mrna

        return segment, pos_in_segment

    def exon_nt_info(self, transcript_id: str, exon_number: int, nt_number: int) -> dict | None:
        """
        Given an exon number and a nucleotide number within the exon (all 1-based), the function returns
        the following information about the nucleotide:
        1. The segment it belongs to (i.e. 5UTR, ORF or 3UTR).
        2. The 1-based position of the nucleotide within the segment.
        3. The nucleotide.

        In case the segment is ORF, the function also returns:
        4. The 1-based codon number within the ORF that this nucleotide belongs to.
        5. The corresponding 1-based nucleotide number within the codon.
        6. The codon.
        7. The amino acid.
        """
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        if (s_tup := self.exon_nt_segment(transcript_id, exon_number, nt_number)) is None:
            return None
        segment, pos_in_segment = s_tup

        info = {"segment": segment, "pos_in_segment": pos_in_segment}

        match segment:
            case "5UTR" | "3UTR":
                utr5 = self.UTR5(transcript_id)
                utr3 = self.UTR3(transcript_id)
                nt = (
                    utr5[pos_in_segment - 1]
                    if segment == "5UTR"
                    else utr3[pos_in_segment - 1]
                )

                info |= {"NT": nt}
            case "ORF":
                orf = self.ORF(transcript_id)
                nt = orf[pos_in_segment - 1]

                # codon number and nt position in codon (all 0-based)
                codon_num, nt_in_codon = (pos_in_segment - 1) // 3, (
                    pos_in_segment - 1
                ) % 3
                codon = orf[codon_num * 3 : (codon_num + 1) * 3]
                aa = tran.convert_aa_32IUPAC(tran.translate_RNA_codon(codon))

                info |= {
                    "NT": nt,
                    "codon_number": codon_num + 1,  # 1-based
                    "nt_in_codon": nt_in_codon + 1,  # 1-based
                    "codon": codon,
                    "aa": aa,
                }
            case _:
                print(f"Unknown {segment=} !!")
                return None

        return info

    def chrm_pos_info(self, transcript_id: str, chrm_pos: int) -> dict | None:
        """
        If position in intron, provides chrm_pos_info() of the Transcript_gff3_cls method.
        If position in exon, provides the details in exon_nt_info() method in addition to
        chrm_pos_info() of the Transcript_gff3_cls method.

        chrm_pos - 1-based position in the chromosome.
        """
        if (base_info := super().chrm_pos_info(transcript_id, chrm_pos)) is None:
            return None

        reg_type, reg_num = base_info["region"].split("_")
        if reg_type == "Exon":
            base_info |= self.exon_nt_info(transcript_id, int(reg_num), base_info["region_pos"])
            # base info contain 'bp' (from exon_nt_info) and 'NT' (from chrm_pos_info), which are the same
            base_info.pop('bp', None)

        return base_info

    def rna_pos2chrm_info(self, transcript_id: str, rna_pos: int) -> dict | None:
        """
        Given a 1-based RNA position rna_pos, the function returns the details 
        in exon_nt_info() method in addition to chromosome position from the 
        menthod rna_pos2chrm_pos() of the Transcript_gff3_cls method.
        """
        if (a := super().rna_pos2chrm_pos(transcript_id, rna_pos)) == -1:
            return None
        return {'chrm_pos': a[0]} | self.chrm_pos_info(transcript_id, a[0])

    def aa_exon_info(self, transcript_id: str, aa_number: int) -> dict | None:
        """
        Given an AA (1-based) position within the chain of AA (ORF), the function returns
        the following information:
        1. The corresponding codon
        2. The corresponding AA
        3. The position of each NT of the codon in the exons, in a form of a list of
           three elements, each of the form 'exon<exon_number>:<1-based position in the exon>'.
        4. The 1-based position of each NT of the codon in the chromosome
        5. The 1-based position of each NT of the codon in the mRNA
        """
        # added YZ, 11/14/24
        def check_pos(x: pd.Series, p: int):
            return x["mRNA_NT_region"][0] <= p <= x["mRNA_NT_region"][1]
        
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        aa = self.AA(transcript_id)

        if not 1 <= aa_number <= len(aa):
            print(f"AA position {aa_number} out of bound (transcript contains {len(aa)} AAs) !!")
            return None

        orf = self.ORF(transcript_id)
        codon = orf[(aa_number - 1) * 3 : aa_number * 3]

        info = {"codon": codon, "AA": aa[aa_number - 1]}
        # determine exon and chromosome position of each NT in codon
        t_info = self.transcripts_info[transcript_id]
        df = self.exon_map(transcript_id)
        m = -1 if self.rev else 1
        # the 1-based mRNA offset of the first NT of the codon (i.e. offset relative to the beginning of exon 1)
        # note that for rev==True, the beginning of an exon is larger than the end of the exon
        # recall that ORF_start_offset is a 0-based value.
        mrna_start_pos = t_info["ORF_start_offset"] + (aa_number - 1) * 3 + 1
        # pos here is actually an offset, thus +1, and +2 hold regardless of self.rev
        mrna_poss = [mrna_start_pos, mrna_start_pos + 1, mrna_start_pos + 2]
        pos_info, chrm_info = [], []

        for pos in mrna_poss:
            # original:
            # df_pos = df.loc[
            #     df.apply(
            #         lambda x: x["mRNA_NT_region"][0] <= pos <= x["mRNA_NT_region"][1],
            #         axis=1,
            #     )
            # ].iloc[0]

            # added YZ, 11/14/24
            df_pos = df.loc[
                df.apply(
                    check_pos, args=(pos,),
                    axis=1,
                )
            ].iloc[0]
            exon_pos = pos - df_pos["mRNA_NT_region"][0] + 1

            pos_info.append(f"Exon_{df_pos['exon_number']}:{exon_pos}")
            # the chromosome position is a function of self.rev
            chrm_info.append(df_pos["exon_region"][0] + m * (exon_pos - 1))

        # Note: for genes encoded on the negative strand (i.e. self.rev==True), exon_start>exon_end and so
        # exon_pos in pos_info corresponds to exon_start-exon_pos+1 chromosome coordinate
        return info | {"codon_exon_pos": pos_info, "codon_chromosome_pos": chrm_info, 'mrna_pos': mrna_poss}

    def AA_mut_to_DNA_SNP_mut(self, aa_mut: str, transcript_id: str) -> dict | None:
        """
        Given a mutation in an amino-acid (aa_mut), in the format:
          <reference AA><AA 1-based position in ORF><mutated AA> (e.g. 'C231S'),

        the function converts aa_mut to (possibly many) corresponding DNA mutations.

        The output is a dictionary, where the key is a mutated codon (that encodes <mutated AA>),
        and the value is a dictionary containing the DNA mutation information:
        1. key='start_pos', value=the 1-based position in the (positive strand) chromosome
        2. key='reference_allele', value=the positive strand reference DNA allele
        3. key='alternative_allele', value=the positive strand mutated DNA allele.

        The DNA mutation can be either SNP, DNP, or TNP.
        """
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        aa_r, aa_pos, aa_m = aa_mut[0], int(aa_mut[1:-1]), aa_mut[-1]
        # aa_pos is the 1-based position of aa in the chain (i.e. the codon 1-based position in the ORF)

        if aa_r == aa_m:
            return None  # no mutation

        if (aa_info := self.aa_exon_info(transcript_id, aa_pos)) is None:
            return None  # aa_pos not valid

        if aa_r != aa_info["AA"]:
            print(f"Position {aa_pos} in the aa chain is {aa_info['AA']}, but the input aa mutation ({aa_mut}) implies that it is {aa_r} !!")
            return None

        ref_codon = aa_info["codon"]

        # codons that encodes aa_m
        possible_mut_codons = [
            tran.RNA2DNA(x)
            for x in tran.RNA_aa_synon_table[tran.convert_aa_IUPAC23(aa_m)]
        ]

        dna_mutations = {}
        for mut_codon in possible_mut_codons:
            # find (0-based) positions of change in the codon
            d = [
                (i, x, y)
                for i, (x, y) in enumerate(zip(ref_codon, mut_codon))
                if x != y
            ]

            mut_i = {}
            match len(d):
                case 1:  # SNP
                    # codon_chromosome_pos is a list of the 3 chromosome positions corresponding to the codon.
                    # It is in descending order for genes encoded on the negative strand, otherwise ascending order
                    # Example:
                    #   for negative strand: codon='AGT', codon_chromosome_pos[12, 11, 10], implying that
                    #   'A' is encoded on the negative strand at chromosome position 12, ...
                    #   for positive strand: codon='AGT', codon_chromosome_pos[71, 72, 72] implying that
                    #   'A' is encoded on the positive strand at chromosome position 71, ...
                    mut_i["start_pos"] = aa_info["codon_chromosome_pos"][d[0][0]]
                    mut_i["reference_allele"] = (
                        tran.reverse_complement(d[0][1].upper())
                        if self.rev
                        else d[0][1]
                    )
                    mut_i["alternative_allele"] = (
                        tran.reverse_complement(d[0][2].upper())
                        if self.rev
                        else d[0][2]
                    )
                case 2:  # DNP or TNP
                    # if both positions are consecutive, then it is DNP, otherwise TNP (as we cant
                    # represent two non-consecutive changes with a SNP or a DNP mutation).
                    if d[1][0] - d[0][0] == 1:  # DNP
                        if d[0][0] == 0:  # changes in the first and second bps in codon
                            # if the gene is encoded on the negative strand, then second codon position,
                            # which corresponds to codon_chromosome_pos[1], is smaller than the first codon
                            # position (codon_chromosome_pos[0]), and so codon_chromosome_pos[1] is the mutation
                            # position on the positive strand.
                            indx_in_chromosome_pos = 1 if self.rev else 0
                            # pick the first two codon's bps
                            codon_start_pos, codon_end_pos = (
                                0,
                                2,
                            )  # end_pos is exclusive
                        else:  # changes in the second and third bps in codon
                            # if the gene is encoded on the negative strand, then third codon position,
                            # which corresponds to codon_chromosome_pos[2], is smaller than the second codon
                            # position (codon_chromosome_pos[1]), and so codon_chromosome_pos[2] is the mutation
                            # position on the positive strand.
                            indx_in_chromosome_pos = 2 if self.rev else 1
                            # pick the second and third codon's bps
                            codon_start_pos, codon_end_pos = (
                                1,
                                3,
                            )  # end_pos is exclusive

                        mut_i["start_pos"] = aa_info["codon_chromosome_pos"][
                            indx_in_chromosome_pos
                        ]
                        r_bps = ref_codon[codon_start_pos:codon_end_pos].upper()
                        m_bps = mut_codon[codon_start_pos:codon_end_pos].upper()
                        mut_i["reference_allele"] = (
                            tran.reverse_complement(r_bps) if self.rev else r_bps
                        )
                        mut_i["alternative_allele"] = (
                            tran.reverse_complement(m_bps) if self.rev else m_bps
                        )
                    else:  # TNP (first and third codon's bps changed, thus we report TNP mutation (where the second bp
                        # is the same in the reference and alternative alleles)
                        # recall that mut_i refers to the positive strand
                        mut_i["start_pos"] = min(aa_info["codon_chromosome_pos"])
                        mut_i["reference_allele"] = (
                            tran.reverse_complement(ref_codon.upper())
                            if self.rev
                            else ref_codon.upper()
                        )
                        mut_i["alternative_allele"] = (
                            tran.reverse_complement(mut_codon.upper())
                            if self.rev
                            else mut_codon.upper()
                        )
                case 3:  # TNP (all three bps in the codon change)
                    # recall that mut_i refers to the positive strand
                    mut_i["start_pos"] = min(aa_info["codon_chromosome_pos"])
                    mut_i["reference_allele"] = (
                        tran.reverse_complement(ref_codon.upper())
                        if self.rev
                        else ref_codon.upper()
                    )
                    mut_i["alternative_allele"] = (
                        tran.reverse_complement(mut_codon.upper())
                        if self.rev
                        else mut_codon.upper()
                    )
                case _:
                    # we verified above that aa_r != aa_m, thus d can not be empty, nor have more than 3 elements.
                    pass

            dna_mutations[mut_codon] = mut_i

        return dna_mutations

    def DNA_SNP_mut_to_AA_mut(
        self,
        dna_ref_allele: str,
        dna_mut_allele: str,
        dna_chrm_pos: int,
        transcript_id: str,
    ) -> str | None:
        """
        Given a DNA mutation within the ORF, defined by:
        1. The positive strand reference allele bp: dna_ref_allele
        2. The positive strand mutated (alternative) allele bp: dna_mut_allele
        3. The 1-based position in the chromosome,

        the function returns the corresponding AA mutation, if the mutation falls within the ORF,
        else None. The AA mutation format is: <reference AA><AA 1-based position in ORF><mutated AA> (e.g. 'C231S').
        If the DNA mutation is outside the ORF, the function returns None.

        INS and DEL DNA mutations are not supported (see the following two methods for INS and DEL).

        SNP, and SNP-like mutations (i.e., DNP, ONP, TNP, ...) are supported. However, we only
        take into consideration the reference and mutated allele bps that are within a codon
        boundary, i.e., we ignore the bps that spill to the next codon. This means that
        we only report the first AA mutation in case of multiple consecutive AA mutations due
        to alleles that span multiple codons.
        """
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        # INS and DEL not supported
        if (dna_ref_allele == "-") or (dna_mut_allele == "-"):
            return None

        # For DNP, TNP, ..., in case of gene encoded on the negative strand, the 'start_position' in the
        # negative strand is dna_chrm_pos + len(dna_mut_allele) - 1, as dna_chrm_pos is the first bp of
        # the DNP, TNP, ... on the positive strand.
        strand_chrm_pos = (
            (dna_chrm_pos + len(dna_mut_allele) - 1) if self.rev else dna_chrm_pos
        )

        if (seg_info := self.chrm_pos_info(transcript_id, strand_chrm_pos)) is None:
            return None  # dna_chrm_pos or transcript_id is invalid
        if "segment" not in seg_info or seg_info["segment"] != "ORF":
            return None  # DNA mutation is not in the ORF.

        mut_aa_pos = seg_info["codon_number"]
        mut_aa_ref = seg_info["aa"]

        # convert alleles to the correct strand
        r = (
            tran.reverse_complement(dna_ref_allele.upper())
            if self.rev
            else dna_ref_allele
        )
        m = (
            tran.reverse_complement(dna_mut_allele.upper())
            if self.rev
            else dna_mut_allele
        )

        # the information in seg_info is in RNA alphabet
        # r = tran.DNA2RNA(r)
        # m = tran.DNA2RNA(m)

        codon = seg_info["codon"]
        codon_pos = seg_info["nt_in_codon"] - 1  # nt_in_codon is 1-based

        # verify that r corresponds to the right bp in the codon.
        # restrict bps to codon boundary
        r_check = r[: 3 - codon_pos]
        to_pos = min(3, codon_pos + len(r))

        if codon[codon_pos:to_pos] != r_check:
            print(
                f"Warning: Position number {codon_pos + 1} in {codon} does not match (positive strand) {dna_ref_allele} !!"
            )

        # generate the mutated codon
        m_in_codon = m[: 3 - codon_pos]  # restrict to codon boundary
        codon_bps = list(codon)
        codon_bps[codon_pos:to_pos] = m_in_codon
        mut_codon = "".join(codon_bps)
        mut_aa_mut = tran.convert_aa_32IUPAC(tran.translate_RNA_codon(mut_codon))

        return f"{mut_aa_ref}{mut_aa_pos}{mut_aa_mut}"

    def DNA_INS_mut_to_affected_AA(self, dna_mut_allele: str, dna_chrm_pos: int, transcript_id: str) -> dict | None:
        """
        Given an insertion mutation, defined by its (exclusive) chromosome (1-based) position and the
        positive strand insertion bps, the function returns a dictionary with key, value:
        1. key=aa_mut, value=<reference AA><AA 1-based position><mutated AA>. This is the first affected AA.
        2. key=codon_mut, value=the mutated (first) codon
        3. key=num_DS_affected_aa, value=number of AAs, downstream from <AA 1-based position>, that are
           directly affected by the insertion sequence. We do not consider here the effect
           of frame shift due to the insertion on all other AA in the ORF.

        Recall that in INS, the bp at dna_chrm_pos is not altered and dna_mut_allele sequence
        is inserted in the positive strand following the bp at dna_chrm_pos.
        """

        """
        In INS, the bp at dna_chrm_pos in the positive strand is not altered.

        Example:
        Given the reference positive strand sequence:

                seq (positive strand) =  -----AG+++++,

        where A is at chromosome position X. Assume an insertion of dna_mut_allele=QQQQ at
        position dna_chrm_pos=X, i.e., the mutated sequence (in the positive strand) is:

              mut_seq (positive strand) =  ----AQQQQG++++,

        then, for gene encoded on the negative strand, the sequence up to C=reverese_complement(G)
        are not altered (these are the ++++ bps downstream from G on the positive strand, up to G),
        and the bps following C (on the negative strand) are reverese_complement(----AQQQQ), so
        in the negative strand, the mutated sequence is ++++CqqqqT----, where
        qqqqT----- = reverese_complement(----AQQQQ).

        Thus, for gene encoded on the negative strand, the chromosome position at which the insertion
        occurs on the negative strand is X+1 (containing C in the example above). Position X+1 in the
        negative strand is not altered, but the following positions (X, X-1, ..) are.

        For gene encoded on the positive strand, the chromosome position at which the insertion
        occurs on the positive strand is X (containing A in the example above). Position X in the
        positive strand is not altered, but the following positions (X+1, X+2, ..) are.
        """
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        strand_chrm_pos = (dna_chrm_pos + 1) if self.rev else dna_chrm_pos

        if (seg_info := self.chrm_pos_info(transcript_id, strand_chrm_pos)) is None:
            return None  # dna_chrm_pos or transcript_id is invalid
        if "segment" not in seg_info or seg_info["segment"] != "ORF":
            return None  # DNA mutation is not in the ORF.

        m = (
            tran.reverse_complement(dna_mut_allele.upper())
            if self.rev
            else dna_mut_allele
        )
        codon = seg_info["codon"]
        codon_pos = seg_info["nt_in_codon"] - 1  # nt_in_codon is 1-based

        # generate the mutated codon (in INS, dna_chrm_pos is not altered)
        mut_codon = (codon[: codon_pos + 1] + m + codon[len(m) :])[:3]
        mut_aa = tran.convert_aa_32IUPAC(tran.translate(mut_codon))

        # number of bps from m used in the current codon
        bps_used_by_mutated_current_codon = 2 - codon_pos
        # number of downstream AA affected from the insertion
        num_ds_affected_aa = int(
            np.ceil((len(m) - bps_used_by_mutated_current_codon) / 3)
        )

        return {
            "aa_mut": f"{seg_info['aa']}{seg_info['codon_number']}{mut_aa}",
            "codon_mut": mut_codon,
            "num_DS_affected_aa": num_ds_affected_aa,
        }

    def DNA_DEL_mut_to_affected_AA(self, dna_ref_allele: str, dna_chrm_pos: int, transcript_id: str) -> dict | None:
        """
        Given a deletion mutation, defined by its (inclusive) chromosome (1-based) position on the
        chromosome (dna_chrm_pos), and the deleted sequence from the positive strand
        (dna_ref_allele), the function returns a dictionary with keys, values:
        1. key=aa_pos, value=the AA position of the first affected AA.
        2. key=aa, value=the (reference) AA at position aa_pos.
        3. key=num_DS_affected_aa, value=the number of AAs, downstream from aa_pos, that are
           directly affected by the deleted sequence. We do not consider here the effect
           of frame shift due to the deletion on all other AA in the ORF.

        'dna_chrm_pos' is the chromosome position of dna_ref_alllele[0] (the first deleted bp).
        'dna_ref_allele' the positive strand sequence deleted.

        In DEL, the first deleted bp on the positive strand is at position dna_chrm_pos.
        Thus, all bps from dna_chrm_pos up to dna_chrm_pos+len(dna_ref_allele)-1 are deleted.

        Example:

              seq (positive strand) =     ----ACGTCAG++++,

        thus

              seq (reversed complement) = ++++CTGACGT----

        where the first A (in the positive strand) is in chromosome position X-1 (so C is in position X).

        Assume a deletion of dna_ref_allele='CGT' at position dna_chrm_pos=X. Then

            mut_seq (positive strand) =     ----ACAG++++

            mut_seq (reverse complement) =  ++++CTGT----

        Thus, for gene encoded on the negative strand, the first deleted bp is A=reverese_complement(T)
        at chromosome position X+len(CGT)-1=X+2, and the following CG=reverese_complement(CG)
        are also deleted.
        """
        if not self.is_protein_coding_transcript(transcript_id):
            print(f"{transcript_id} is not a protein coding transcript !!")
            return None

        strand_chrm_pos = (
            (dna_chrm_pos + len(dna_ref_allele) - 1) if self.rev else dna_chrm_pos
        )

        # note: if strand_pos_pos is outside of the ORF, but some of the bps deleted are within,
        # then the following will result in None, even though there is an effect on the ORF.
        if (seg_info := self.chrm_pos_info(transcript_id, strand_chrm_pos)) is None:
            return None  # dna_chrm_pos or transcript_id is invalid
        if "segment" not in seg_info or seg_info["segment"] != "ORF":
            return None  # DNA mutation is not in the ORF.

        codon_pos = seg_info["nt_in_codon"] - 1  # nt_in_codon is 1-based

        # number of bps from dna_ref_allele that are within the current codon
        # recall that the first deleted bp is at dna_chrm_pos
        # which corresponds to (0-based) codon_pos in the codon.
        bps_deleted_from_current_codon = 3 - codon_pos

        # following codons/AA affected from the deletion
        num_ds_affected_aa = int(
            np.ceil((len(dna_ref_allele) - bps_deleted_from_current_codon) / 3)
        )

        return {
            "aa_pos": seg_info["codon_number"],
            "aa": seg_info["aa"],
            "num_DS_affected_aa": num_ds_affected_aa,
        }

    def __repr__(self):
        return f"Gene_gff3_cls(gene={self.gene}, gff3_file|(gff3_df, gff3_df_gene_type))"

    def __str__(self):
        return "Gene annotation class based on Ensembl GFF3 annotation file.\n"


class Gene_cls(Gene_gff3_cls):
    """
    Gene annotation class based on gene name.

    User should use this class.

    Instantiating this class requires the following inputs:
    1. gene name or gene ID (ENS<species prefix>G, where <species prefix> is empty for Homo sapiens).
    2. Either the GFF3 file (Path), or a tuple of the GFF3 dataframe and its subset dataframe
       (containing only rows with Type value defined by the user (deafult is egna.Gene_type_values)).
       Use the function ensembl_gff3_df to generate this tuple if instantiating using the tuple.
    3. Optional - verbose flag.
    4. Optional - species
    5. Optional - chromosome Fasta file.
    """
    def __init__(self, gene: str, gff3_source: Path | tuple, 
                 species: str = 'homo_sapiens',
                 chrm_fasta_file: Path | None = None,
                 verbose: bool = True) -> None:
        super().__init__(gene, gff3_source, species=species, chrm_fasta_file=chrm_fasta_file, verbose=verbose)

    def __repr__(self) -> str:
        return f"Gene_cls('{self.gene}, gff3_file|(gff3_df, gff3_df_gene_type)')"

    def __str__(self) -> str:
        return"Gene annotation class based on gene name and annotation data (annotation file or dataframes).\n"
