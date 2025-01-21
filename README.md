# Vertebrates Gene Annotation

A Python package that annotates genes and transcripts in vertebrates species based on Ensembl.

# Python

Requires python>=3.10

# Installing
`geneanot` is available on PyPI:
```console
pip install geneanot
```

Running it provides a documentation link:

```console
$ geneanot

Welcome to the vertebrates gene annotation package geneanot.
Please visit https://github.com/yoramzarai/geneanot for documentation.
```

# Requirements
A designated local folder is required to hold the Ensembl annotation file.

Please consult the [usage notebook](https://github.com/yoramzarai/geneanot/blob/main/Scripts/usage_examples.ipynb) for more information.

## Chromosome Data
Some of the methods require the chromosome sequence data. Two chromosome data **access modes** are supported:
- `local` - user provides the corresponding chromosome Fasta file
- `remote` - `genenot` uses the Ensembl REST API to extract sequence from the chromosome

`local` access implies faster sequence retrival from the chromosome, whereas `remote` access requires network connection and implies a slower sequence retrival.

The provided chromosome Fasta file, in case of a `local` access mode, MUST contain equal number of bps per row in all sequence
rows (other than possibly the last row). The chromosome Fasta file can be either gzip compressed (with a `.gz` suffix) or uncompressed, and can be downloaded from Ensembl (Ensembl chromosome Fasta files contain equal number of bps per sequence row).

Please consult the [usage notebook](https://github.com/yoramzarai/geneanot/blob/main/Scripts/usage_examples.ipynb) for more information.

# Usage
See the [usage notebook](https://github.com/yoramzarai/geneanot/blob/main/Scripts/usage_examples.ipynb) for a detailed usage description.

Here are few basic usage examples (assuming the code is executed from the `.Scripts/` folder). We start with annotating a Homo sapiens gene (Homo sapiens is the default assumed species in `geneanot`).

## Instantiating an annotation object
```python
from pathlib import Path
import geneanot as u

# folder to hold Ensembl annotation file.
Annotation_folder: Path = Path('./../AnnotationDB')

# get/update annotation file
download_done, ensembl_file, local_file = u.update_local_release_to_latest(Annotation_folder, enable_download=True)
annotation_full_file = Annotation_folder / (ensembl_file if download_done else local_file)

# instantiate an annotation object (can use gene name or gene ID) in remote access mode
gA = u.Gene_cls('EGFR', annotation_full_file, verbose=True)

# or - instantiate (can use gene name or gene ID) in local access mode
chrm_fasta_file: str = 'Homo_sapiens.GRCh38.dna_sm.chromosome.7.fa'  # change to your path
gA = u.Gene_cls('EGFR', annotation_full_file, chrm_fasta_file=chrm_fasta_file, verbose=True)
```

## Basic gene annotation
```python
# show gene basic information
gA.info()

# show all transcripts
print(f"\n{gene} contains {len(gA)} transcripts.")
gA.show_all_transcript_IDs()

# gene start and end coordinates
print(f"\nGene start = {gA.gene_start:,}, gene end = {gA.gene_end:,}")

print(f"\n{gene=} encoded in chromosome {gA.chrm} on the {'negative' if gA.rev else 'positive'} DNA strand.")

# the various gene attributes printed from the above commands are also available via gA
print(gA.gene_name, gA.gene_ID, gA.gene_desc, gA.gene_type, gA.gene_ver, gA.rev, gA.chrm, gA.gene_start, gA.gene_end, len(gA.transcripts_info), sep='\n')
```

## Transcript annotation
```python
transcript_id: str = 'ENST00000275493'

# show transcript basic information
gA.transcript_info(transcript_id, verbose=True)

# the transcript features are accessible via gA.transcripts_info[transcript_id], for example
t_info = gA.transcripts_info[transcript_id]
print(t_info['transcript_name'], t_info['transcript_v_id'], t_info['transcript_ver'])
```

## Transcript and RNA tables

The transcript table lists the exons and introns in the transcript. Following is a <ins> **partial** </ins> transcript table of the transcript `ENST00000275493`.

<img src="https://github.com/yoramzarai/geneanot/blob/main/metadata/figs/partial_transcript_table.png" alt="" width="380" height="800">

The RNA table lists the exons in the RNA, and in case of protein-coding transcripts, maps the ORF and the AAs to the exons. Following is the RNA table of the transcript `ENST00000275493`.

![](https://github.com/yoramzarai/geneanot/blob/main/metadata/figs/mRNA_table.png)


Here `mRNA_NT_region` gives the mRNA bp count, `ORF_NT_region` gives the ORF bp count, `ORF_AA_region` gives the protein AA count, and `next_exon_frame_alignment` gives the codon phase (i.e., number of bps required in the next exon to complete the last codon in the current exon).

```python
# Transcript table
df = gA.exon_intron_map(transcript_id)
display(df)

# the transcript table can be written to excel, CSV or HTML file
gA.exon_intron_map_to_csv(transcript_id, './../Reports/transcript_table.csv')
gA.exon_intron_map_to_excel(transcript_id, './../Reports/transcript_table.xlsx', usr_desc={"Description": "Transcript table", "Transcript": transcript_id})
gA.exon_intron_map_to_html(transcript_id, './../Reports/transcript_table.html')

# RNA table
df = gA.exon_map(transcript_id)
display(df)

# this table can also be written to excel, CSV or HTML file
gA.exon_map_to_csv(transcript_id, './../Reports/mRNA_table.csv')
gA.exon_map_to_excel(transcript_id, './../Reports/mRNA_table.xlsx', usr_desc={"Description": "mRNA table", "Transcript": transcript_id})
gA.exon_map_to_html(transcript_id, './../Reports/mRNA_table.html')
```

## Sequences
```python
pre_RNA_seq = gA.seq(transcript_id).upper()
print(f"pre-RNA contains {len(pre_RNA_seq):,} bps.")

rna_seq = gA.rna(transcript_id).upper()
print(f"\nrna:\n{rna_seq}\n{len(rna_seq):,} bps.")

orf_seq = gA.ORF(transcript_id).upper()
print(f"\norf:\n{orf_seq}\n{len(orf_seq):,} bps.")

aa_seq = gA.AA(transcript_id)
print(f"\nprotein:\n{aa_seq}\n{len(aa_seq):,} AAs.")

utr5_seq = gA.UTR5(transcript_id).upper()
print(f"\n5'UTR:\n{utr5_seq}\n{len(utr5_seq):,} bps.")

utr3_seq = gA.UTR3(transcript_id).upper()
print(f"\n3'UTR=\n{utr3_seq}\n{len(utr3_seq):,} bps.")

exon_seq, seq_info = gA.exon_intron_seq('Exon', 5, transcript_id)
print(f"\nExon =\n{exon_seq}\n{seq_info}")

intron_seq, seq_info = gA.exon_intron_seq('Intron', 5, transcript_id)
print(f"\nIntron =\n{intron_seq}\n{seq_info.upper()}")

# A modified transcript
use_exon_list: list[int] = [1, 3]  # exon numbers to use
use_intron_list: list[int] = [2]  # intron numbers to use
m_seq = gA.modified_transcript(use_exon_list, use_intron_list, transcript_id)
print(f"\nA modified transcript containing exons {', '.join(map(str,use_exon_list))} and introns {', '.join(map(str,use_intron_list))} contains {len(m_seq.upper()):,} bps.")
```

## Queries
```python
# query the start and stop codon positions in the chromosome and in the RNA
# --------------------------------------------------------------------------
start_codon_info, stop_codon_info = gA.start_and_stop_codons_pos(transcript_id)

# query the start and end ORF positions in the exons
# --------------------------------------------------
orf_start_end_info = gA.orf_start_and_end_exon_info(transcript_id)

# query a chromosome position
chrm_pos: int = 55_157_663
# ---------------------------
chrm_info = gA.chrm_pos_info(transcript_id, chrm_pos)
print(chrm_info)

# query a RNA position
rna_pos: int = 685
# ----------------
chrm_info = gA.rna_pos2chrm_info(transcript_id, rna_pos)
print(chrm_info)

# map a chromosome position to the RNA position
chrm_pos: int = 55_211_628 
# ------------------------
rna_pos = gA.chrm_pos2rna_pos(transcript_id, chrm_pos)

# map a RNA position to the chromosome position
rna_pos: int = 685
# ----------------
chrm_p = gA.rna_pos2chrm_pos(transcript_id, rna_pos)
print(f"{rna_pos=} --> {chrm_p=}")

# query an exon position
exon_number: int = 7
bp_index_in_exon: int = 47
# -------------------
info = gA.exon_nt_info(transcript_id, exon_number, bp_index_in_exon)
print(info)

# query the RNA segment (UTR, ORF) of an exon position
exon_number: int = 28
bp_index_in_exon: int = 400
# --------------------------
exon_segment_info := gA.exon_nt_segment(transcript_id, exon_number, bp_index_in_exon)

# query an amino-acid position
aa_number: int = 163
# ------------------
aa_info = gA.aa_exon_info(transcript_id, aa_number)
print(aa_info)

# query an AA variant given a DNA variant
ref_allele, var_allele, chromosome_pos = 'G', 'C', 55_152_609
# ------------------------------------------------------------
aa_var = gA.DNA_SNP_mut_to_AA_mut(ref_allele, var_allele, chromosome_pos, transcript_id)
print(f"chr{gA.chrm}:{chromosome_pos}:{ref_allele}>{var_allele} --> {aa_var=}")

# query (all) DNA variants given an AA variant
aa_var: str = 'C231S'
# -------------------
dna_all_muts = gA.AA_mut_to_DNA_SNP_mut(aa_var, transcript_id)
print(f"{aa_var=} corresponds to the following DNA variants:")
for i, (codon, var_info) in enumerate(dna_all_muts.items(), start=1):
    print(f"{i}. {codon=}, {var_info}")
```


## Annotating multiple genes (faster instantiation)
```python
# parse annotation file into annotation dataframes
gff3_dfs = u.ensembl_gff3_df(annotation_full_file)  

# now instantiate with the annotation dataframes (faster)
# can instantiate in local mode by providing chrm_fasta_file (see above)
g_a1 = u.Gene_cls('EGFR', gff3_dfs)
g_a1.info()
g_a2 = u.Gene_cls('BRCA1', gff3_dfs)
g_a2.info()
g_a3 = u.Gene_cls('IDH1', gff3_dfs)
g_a3.info()
...
```

## Annotating other vertebrates species
```python
# Example: Mus_musculus

species: str = 'mus_musculus'  # required for non Homo sapiens species.

# The suggest annotation file name (used to extract annotation signature, which is required next)
suggest_annotation_file_name, _, release_n = u.suggested_annotation_file_name(species=species)

# update/download annotation file
download_done, ensembl_file, local_file = u.update_local_release_to_latest(
    Annotation_folder, 
    enable_download=True, 
    gff3_pattern=suggest_annotation_file_name.replace(release_n, 'XXX'),
    species=species)

annotation_full_file = Annotation_folder / (ensembl_file if download_done else local_file)

# instantiate annotation class (can use gene name or gene ID) in remote mode (see above for a local mode example)
gA = u.Gene_cls('ENSMUSG00000017167', annotation_full_file, species=species, verbose=True)
gA.info()

print(f"\n{gene} contains {len(gA)} transcripts.")
gA.show_all_transcript_IDs()

transcript_id: str = 'ENSMUST00000103109'
df = gA.exon_map(transcript_id)
print(df.to_string())

# protein sequence
aa_seq = gA.AA(transcript_id)
print(f"\nprotein:\n{aa_seq}\n{len(aa_seq):,} AAs.")
```

## Extra
`geneanot` supports also an arbitrary sequence retrievel from a chromosome, using either a chromosome Fasta file, or remotely via the Ensembl REST API. 

```python
# fetching using a chromosome Fasta file
chrm_fasta_file: Path = Path('./../Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa')
start_p: int = 122_989_200  # 1-based start coordinate
end_p: int = 122_989_229  # 1-based end coordinate
rev: bool = False  # True to fetch the reveresed complement sequence

seq = u.fetch_seq(chrm_fasta_file, start_p, end_p, rev=rev)
print(seq)


# fetching remotely (this requires a network connection)
chrm: str = '5'  # chromosome number, e.g., '1' or 'Y'
start_p: int = 52_120_100  # 1-based start coordinate
end_p: int = 52_120_120  # 1-based end coordinate
species: str = "Danio_rerio"  # species
rev: bool = True  # Troe to fetch the reveresed complement sequence

seq = u.fetch_seq(chrm, start_p, end_p, rev=rev, species=species)
print(seq)
```

# References
- [Ensembl](http://www.ensembl.org/index.html)
- [Ensembl GFF3 File Format](https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/README)
- [GFF3 Format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
- [Ensembl REST](https://rest.ensembl.org)
