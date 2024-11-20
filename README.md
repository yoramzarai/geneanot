# Eukaryotes Gene Annotation

A Python package that annotates eukaryotes genes and transcripts based on Ensembl.

# Python Version
Python>=3.10

# Installing

GeneAnot is available on PyPI:
```console
pip install geneanot
```

# Requirements
- A designated local folder to hold the Ensembl annotation file
- Some of the methods (e.g., sequence retrieval methods) require the chromosome Fasta file

Please consult the [usage notebook](/Scripts/usage_examples.ipynb) for more information.

# Usage
See the [usage notebook](/Scripts/usage_examples.ipynb) for a detailed usage description.

Here are few basic usage examples (assuming the code is executed from the `.Scripts/` folder):

Instantiating an annotation object:
```python
from pathlib import Path
import geneanot as u

# folder to hold Ensembl annotation file.
Annotation_folder: Path = Path('./../AnnotationDB')

# get/update annotation file
download_done, ensembl_file, local_file = u.update_local_release_to_latest(Annotation_folder, enable_download=True)
annotation_full_file = Annotation_folder / (ensembl_file if download_done else local_file)

# instantiate annotation class (can use gene name or gene ID)
gA = u.Gene_cls('EGFR', annotation_full_file, verbose=True)
```

Basic gene annotation:
```python
# show gene basic information
gA.info()

# show all transcripts
print(f"\n{gene} contains {len(gA)} transcripts.")
gA.show_all_transcript_IDs()

# gene start and end coordinates
print(f"\nGene start = {gA.gene_start:,}, gene end = {gA.gene_end:,}")

print(f"{gene=} encoded on the {'negative' if gA.rev else 'positive'} DNA strand.")

# the various gene attributes above are available using gA
print(gA.gene_name, gA.gene_ID, gA.gene_desc, gA.gene_type, gA.gene_ver, gA.rev, gA.chrm, gA.gene_start, gA.gene_end, len(gA.transcripts_info), sep='\n')
```

Annotation of a transcript:
```python
# show transcript basic information
gA.transcript_info('ENST00000275493', verbose=True)

# the transcript features are accessible via gA.transcripts_info[transcript_id], for example
t_info = gA.transcripts_info[transcript_id]
print(t_info['transcript_name'], t_info['transcript_v_id'], t_info['transcript_ver'])
```

Transcript and mRNA tables
```python
# Transcript table
df = gA.exon_intron_map(transcript_id)
display(df)

# the transcript table can be written to excel, CSV or HTML file
gA.exon_intron_map_to_csv(transcript_id, './../Reports/transcript_table.csv')
gA.exon_intron_map_to_excel(transcript_id, './../Reports/transcript_table.xlsx', usr_desc={"Description": "Transcript table", "Transcript": transcript_id})
gA.exon_intron_map_to_html(transcript_id, './../Reports/transcript_table.html')

# mRNA table
df = gA.exon_map(transcript_id)
display(df)

# similarly to exon_intron_map() above, this table can also we written to excel, CSV or HTML file
gA.exon_map_to_csv(transcript_id, './../Reports/mRNA_table.csv')
gA.exon_map_to_excel(transcript_id, './../Reports/mRNA_table.xlsx', usr_desc={"Description": "mRNA table", "Transcript": transcript_id})
gA.exon_map_to_html(transcript_id, './../Reports/mRNA_table.html')
```

Sequences
```python
chrm_fasta_file: str = f'./../Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.{gA.chrm}.fa'  # replace with your fasta file

pre_mRNA_seq = gA.seq(transcript_id, chrm_fasta_file).upper()
print(f"pre-mRNA contains {len(pre_mRNA_seq):,} bps.")

rna_seq = gA.rna(transcript_id, chrm_fasta_file).upper()
print(f"\nrna=\n{rna_seq}")

orf_seq = gA.ORF(transcript_id, chrm_fasta_file).upper()
print(f"\norf=\n{orf_seq}")

aa_seq = gA.AA(transcript_id, chrm_fasta_file)
print(f"\nAA=\n{aa_seq}")

utr5_seq = gA.UTR5(transcript_id, chrm_fasta_file).upper()
print(f"\nutr5=\n{utr5_seq}")

utr3_seq = gA.UTR3(transcript_id, chrm_fasta_file).upper()
print(f"\nutr3=\n{utr3_seq}")

exon_seq, seq_info = gA.exon_intron_seq('Exon', 5, transcript_id, chrm_fasta_file)
print(f"\nExon =\n{exon_seq}\n{seq_info}")

intron_seq, seq_info = gA.exon_intron_seq('Intron', 5, transcript_id, chrm_fasta_file)
print(f"\nIntron =\n{intron_seq}\n{seq_info.upper()}")
```

Queries
```python
# query a chromosome position
chrm_pos: int = 55_157_663
# ---------------------------
chrm_info = gA.chrm_pos_info(transcript_id, chrm_pos, chrm_fasta_file)
print(chrm_info)

# query RNA position
rna_pos: int = 685
# ----------------
chrm_info = gA.rna_pos2chrm_info(transcript_id, rna_pos, chrm_fasta_file)
print(chrm_info)

# query AA variant based on DNA variant
ref_allele, var_allele, chromosome_pos = 'G', 'C', 55_152_609
# ------------------------------------------------------------
aa_var = gA.DNA_SNP_mut_to_AA_mut(ref_allele, var_allele, chromosome_pos, transcript_id, chrm_fasta_file)
print(f"chr{gA.chrm}:{chromosome_pos}:{ref_allele}>{var_allele} --> {aa_var=}")

# query (all) DNA variants based on an AA variant
aa_var: str = 'C231S'
# -------------------
dna_all_muts = gA.AA_mut_to_DNA_SNP_mut(aa_var, transcript_id, chrm_fasta_file)
print(f"{aa_var=} corresponds to the following DNA variant:")
for i, (codon, var_info) in enumerate(dna_all_muts.items(), start=1):
    print(f"{i}. {codon=}, {var_info}")
```


Annotating multiple genes (faster instantiation)
```python
# parse annotation file into annotation dataframes
gff3_dfs = u.ensembl_gff3_df(annotation_full_file)  

# now instantiate with the annotation dataframes (faster)
g_a1 = u.Gene_cls('EGFR', gff3_dfs)
g_a1.info()
g_a2 = u.Gene_cls('BRCA1', gff3_dfs)
g_a2.info()
g_a3 = u.Gene_cls('IDH1', gff3_dfs)
g_a3.info()
```

Annotating other Eukaryotes Species
```python
# Example: Mus_musculus

# Annotation file signature
annotation_file_signature: str = 'Mus_musculus.GRCm39.XXX.gff3.gz'
# Ensembl FTP URL
ensembl_url: str = 'rsync://ftp.ebi.ac.uk/ensemblorg/pub/current_gff3/mus_musculus'

# update/download annotation file
download_done, ensembl_file, local_file = u.update_local_release_to_latest(Annotation_folder, 
                                                                           enable_download=True, 
                                                                           gff3_pattern=annotation_file_signature,
                                                                           ensembl_url=ensembl_url)
annotation_full_file = Annotation_folder / (ensembl_file if download_done else local_file)

# instantiate annotation class (can use gene name or gene ID)
gA = u.Gene_cls('ENSMUSG00000017167', annotation_full_file, verbose=True)
gA.info()

print(f"\n{gene} contains {len(gA)} transcripts.")
gA.show_all_transcript_IDs()

transcript_id: str = 'ENSMUST00000103109'
df = gA.exon_map(transcript_id)
print(df.to_string())
```
