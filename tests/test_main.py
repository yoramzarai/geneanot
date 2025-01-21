""" 
Activate environment and run "pytest".

The "expected" in the tests below are based on release 113 (Homo_sapiens.GRCh38.113.gff3.gz).
"""
from pathlib import Path

import pytest

import geneanot as u

human_gff3_dfs = u.ensembl_gff3_df(Path('AnnotationDB/Homo_sapiens.GRCh38.113.gff3.gz'))

def get_homo_sapiens_object(gene: str) -> u.Gene_cls:
    g = u.Gene_cls(gene, human_gff3_dfs, verbose=False)
    # if the chromosome file exists, use it (local mode - faster retrival), otherwise use remote access.
    chrm_fasta_file = Path(f'Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.{g.chrm}.fa')
    g.chrm_fasta_file = chrm_fasta_file if Path(chrm_fasta_file).is_file() else None
    return g

def get_mus_musculus_object(gene: str) -> u.Gene_cls:
    annotation_full_file: Path = Path("AnnotationDB/Mus_musculus.GRCm39.113.gff3.gz")
    g = u.Gene_cls(gene, annotation_full_file, species="mus_musculus", verbose=False)
    chrm_fasta_file = Path(f'Chromosome/Mus_musculus.GRCm39.dna_sm.chromosome.{g.chrm}.fa')
    g.chrm_fasta_file = chrm_fasta_file if chrm_fasta_file.is_file() else None
    return g


# Homo sapiens 
# ============
@pytest.mark.parametrize(
    "gene, expected",
    [
        # homo_sapiens
        ("EGFR", (55_019_017, 55_211_628, '7', False)),
        ("IDH1", (208_236_229, 208_266_074, '2', True)),
    ],
)
def test_gene_start_end_chrm_rev(gene: str, expected: tuple[int, int, str, bool]) -> None:
    g = get_homo_sapiens_object(gene)
    assert (g.gene_start, g.gene_end, g.chrm, g.rev) == expected


@pytest.mark.parametrize(
    "gene, expected",
    [
        # homo_sapiens
        ("EGFR", 13),
        ("IDH1", 9),
    ],
)
def test_gene_num_transcripts(gene: str, expected: int) -> None:
    g = get_homo_sapiens_object(gene)
    assert len(g) == expected


@pytest.mark.parametrize(
    "gene, transcript, expected",
    [
        # homo_sapiens
        ("EGFR", "ENST00000275493", 192_612),
        ("IDH1",  "ENST00000345146", 18_843),
    ],
)
def test_gene_transcript_premrna_size(gene: str, transcript: str, expected: int) -> None:
    g = get_homo_sapiens_object(gene)
    assert len(g.seq(transcript)) == expected


@pytest.mark.parametrize(
    "gene, transcript, expected_partial_protein_seq",
    # sequences here are not complete (first AAs)
    [
        # these two are encoded on the positive strand
        ("EGFR", "ENST00000275493", "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKL"),
        #("AEBP1", "ENST00000223357", "MAAVRGAPLLSCLLALLALCPGGRPQTVLTDDEIEEFLEGFLSELEPEPREDDVEAPPPPEPTPRVRKAQAGGKPGKRPGTAAE"),
        # IDH1 encoded on the negative strand
        ("IDH1", "ENST00000345146", "MSKKISGGSVVEMQGDEMTRIIWELIKEKLIFPYVELDLHSYDLGIENRDATNDQVTKDAAEAIKKHNVGVKCATITPDEKRVEEFKLKQMWKSPNGTIRNILGGTVFRE")
    ],
)
def test_gene_transcript_protein_seq(gene: str, transcript: str, expected_partial_protein_seq: str) -> None:
    g = get_homo_sapiens_object(gene)
    #g.chrm_fasta_file = (f"Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.{g.chrm}.fa")
    if (protein_seq := g.AA(transcript)) is None:
        raise ValueError
    assert protein_seq[:len(expected_partial_protein_seq)].upper() == expected_partial_protein_seq.upper()


@pytest.mark.parametrize(
    "gene, transcript, expected_partial_utr3_seq",
    # sequences here are not complete (first bps)
    [
        ("EGFR", "ENST00000275493", "CCACGGAGGATAGTATGAGCCCTAAAAATCCAGACTCTTTCGATACCCAGGACCAAGCCACAGCAGGTCCTCCATCCCAACAGCCATGCCCGCATTAGCTCTTAGACCCACAGACTGGTTTTGCAACGTTTA"),
        ("IDH1", "ENST00000345146", "GTTCATACCTGAGCTAAGAAGGATAATTGTCTTTTGGTAACTAGGTCTACAGGTTTACATTTTTCTGTGTTACACTCAAGGATAAAGGCAAAATCAATTTTGTAATTTGTTTAGAAGCCAGAGTTTATCTTTTCTATAAGTTTACAGCCT")
    ],
)
def test_gene_transcript_UTR3_seq(gene: str, transcript: str, expected_partial_utr3_seq: str) -> None:
    g = get_homo_sapiens_object(gene)
    if (seq := g.UTR3(transcript)) is None:
        raise ValueError
    assert seq[:len(expected_partial_utr3_seq)].upper() == expected_partial_utr3_seq.upper()


@pytest.mark.parametrize(
    "gene, transcript, expected_partial_utr5_seq",
    # sequences here are not complete (first bps)
    [
        ("EGFR", "ENST00000275493", "AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCGGCCGAGGCGGCCGGAGTCCCGAGCTAGCCCCGGCGGCCGC"),
        ("IDH1", "ENST00000345146", "GGGTTTCTGCAGAGTCTACTTCAGAAGCGGAGGCACTGGGAGTCCGGTTTGGGATTGCCAGGCTGTGGTTGTGAGTCTGAGCTTGTGAGCGGCTGTGGCGCCCCAACTCTTCGCCAGCATATCATCCCGG")
    ],
)
def test_gene_transcript_UTR5_seq(gene: str, transcript: str, expected_partial_utr5_seq: str) -> None:
    g = get_homo_sapiens_object(gene)
    if (seq := g.UTR5(transcript)) is None:
        raise ValueError
    assert seq[:len(expected_partial_utr5_seq)].upper() == expected_partial_utr5_seq.upper()

@pytest.mark.parametrize(
    "gene, transcript, expected",
    [
        ("EGFR", "ENST00000275493", (55_019_017, 55_211_628)),   # gene encoded in the positive strand
        ("IDH1", "ENST00000345146", (208_255_071, 208_236_229)), # gene encoed on the negative strand
    ],
)
def test_transcript_start_and_end(gene: str, transcript: str, expected: tuple[int,int]) -> None:
    g = get_homo_sapiens_object(gene)
    t_info = g.transcripts_info[transcript]
    assert t_info['transcript_start'] == expected[0] and t_info['transcript_end'] == expected[1]

@pytest.mark.parametrize(
    "gene, transcript, expected",
    [
        ("EGFR", "ENST00000275493", (55_019_278, 55_205_615)),
        ("IDH1", "ENST00000345146",(208_251_551, 208_237_081))
    ],
)
def test_transcript_start_end_chrm_codon(gene: str, transcript: str, expected: tuple[int,int]) -> None:
    g = get_homo_sapiens_object(gene)
    if (start_end_stop := g.start_and_stop_codons_pos(transcript)) is None:
        raise ValueError(f"{transcript=} not recognized !!")
    start_codon_info, stop_codon_info = start_end_stop
    assert start_codon_info[0] == expected[0] and stop_codon_info[0] == expected[1]

@pytest.mark.parametrize(
    "gene, transcript, chrm_pos, expected",
    [
        ("EGFR", "ENST00000275493", 55_157_663, {
            'region': 'Exon_11', 'region_pos': 1, 'dist_from_region_boundary': (0, 90), 'segment': 'ORF', 'pos_in_segment': 1208, 'NT': 'G', 'codon_number': 403, 'nt_in_codon': 2, 'codon': 'GGG', 'aa': 'G'}),
        ("IDH1", "ENST00000345146", 208_243_594, {'region': 'Exon_6', 'region_pos': 11, 'dist_from_region_boundary': (10, 167), 'segment': 'ORF', 'pos_in_segment': 531, 'NT': 'T', 'codon_number': 177, 'nt_in_codon': 3, 'codon': 'GGT', 'aa': 'G'}),
    ],
)
def test_query_chrm_pos(gene: str, transcript: str, chrm_pos: int, expected: dict) -> None:
    g = get_homo_sapiens_object(gene)
    assert g.chrm_pos_info(transcript, chrm_pos) == expected

@pytest.mark.parametrize(
    "gene, transcript, chrm_pos, expected",
    [
        ("EGFR", "ENST00000275493", 55_211_628, (9905, 'A')),
        ("IDH1", "ENST00000345146", 208_236_229, (2_318, 'A')),
    ],
)
def test_map_chrm_pos(gene: str, transcript: str, chrm_pos: int, expected: tuple[int,str]) -> None:
    g = get_homo_sapiens_object(gene)
    if (rna_pos := g.chrm_pos2rna_pos(transcript, chrm_pos)) is None:
        raise ValueError(f"{chrm_pos=:,} does not overlap with RNA.")
    assert rna_pos == expected

@pytest.mark.parametrize(
    "gene, transcript, rna_pos, expected",
    [
        ("EGFR", "ENST00000275493", 685, (55143488, 'G')),
        ("IDH1", "ENST00000345146", 685, (208245377, 'A')),
    ],
)
def test_map_rna_pos(gene: str, transcript: str, rna_pos: int, expected: tuple[int,str]) -> None:
    g = get_homo_sapiens_object(gene)
    assert g.rna_pos2chrm_pos(transcript, rna_pos) == expected

@pytest.mark.parametrize(
    "gene, transcript, rna_pos, expected",
    [
        ("EGFR", "ENST00000275493", 685, {'chrm_pos': 55143488, 'region': 'Exon_3', 'region_pos': 184, 'dist_from_region_boundary': (183, 0), 'segment': 'ORF', 'pos_in_segment': 424, 'NT': 'G', 'codon_number': 142, 'nt_in_codon': 1, 'codon': 'GAA', 'aa': 'E'}),
        ("IDH1", "ENST00000345146", 685, {'chrm_pos': 208245377, 'region': 'Exon_5', 'region_pos': 48, 'dist_from_region_boundary': (47, 58), 'segment': 'ORF', 'pos_in_segment': 462, 'NT': 'A', 'codon_number': 154, 'nt_in_codon': 3, 'codon': 'ATA', 'aa': 'I'}),
    ],
)
def test_query_rna_pos(gene: str, transcript: str, rna_pos: int, expected: dict) -> None:
    g = get_homo_sapiens_object(gene)
    assert g.rna_pos2chrm_info(transcript, rna_pos) == expected

@pytest.mark.parametrize(
    "gene, transcript, exon_number, nt_number, expected",
    [
        ("EGFR", "ENST00000275493", 7, 47, {'segment': 'ORF', 'pos_in_segment': 794, 'NT': 'C', 'codon_number': 265, 'nt_in_codon': 2, 'codon': 'CCC', 'aa': 'P'}),
        ("IDH1", "ENST00000345146", 7, 47, {'segment': 'ORF', 'pos_in_segment': 745, 'NT': 'A', 'codon_number': 249, 'nt_in_codon': 1, 'codon': 'AGG', 'aa': 'R'}),
    ],
)
def test_query_exon_pos(gene: str, transcript: str, exon_number: int, nt_number: int, expected: dict) -> None:
    g = get_homo_sapiens_object(gene)
    assert g.exon_nt_info(transcript, exon_number, nt_number) == expected

@pytest.mark.parametrize(
    "gene, transcript, exon_number, nt_number, expected",
    [
        ("EGFR", "ENST00000275493", 28, 400, ('3UTR', 38)),
        ("IDH1", "ENST00000345146", 5, 74, ('ORF', 488)),
    ],
)
def test_query_exon_segment(gene: str, transcript: str, exon_number: int, nt_number: int, expected: tuple[str,int]) -> None:
    g = get_homo_sapiens_object(gene)
    assert g.exon_nt_segment(transcript, exon_number, nt_number) == expected

@pytest.mark.parametrize(
    "gene, transcript, aa_number, expected",
    [
        ("EGFR", "ENST00000275493", 163, {'codon': 'CAG', 'AA': 'Q', 'codon_exon_pos': ['Exon_4:63', 'Exon_4:64', 'Exon_4:65'], 'codon_chromosome_pos': [55146668, 55146669, 55146670], 'mrna_pos': [748, 749, 750]}),
        ("IDH1", "ENST00000345146", 49, {'codon': 'CGT', 'AA': 'R', 'codon_exon_pos': ['Exon_4:23', 'Exon_4:24', 'Exon_4:25'], 'codon_chromosome_pos': [208248638, 208248637, 208248636], 'mrna_pos': [368, 369, 370]}),
    ],
)
def test_query_aa_pos(gene: str, transcript: str, aa_number: int, expected: dict) -> None:
    g = get_homo_sapiens_object(gene)
    assert g.aa_exon_info(transcript, aa_number) == expected

@pytest.mark.parametrize(
    "gene, transcript, ref_allele, var_allele, pos, expected_aa_var",
    [
        ("EGFR", "ENST00000275493", "G", "C", 55_152_609, "C231S"),
        ("IDH1", "ENST00000345146", "C", "T", 208_248_455, "E110K"),
    ],
)
def test_DNA2AA_var(gene: str, transcript: str, ref_allele: str, var_allele: str, pos: int, expected_aa_var: str) -> None:
    g = get_homo_sapiens_object(gene)
    assert g.DNA_SNP_mut_to_AA_mut(ref_allele, var_allele, pos, transcript) == expected_aa_var

@pytest.mark.parametrize(
    "gene, transcript, aa_var, expected_dna_var",
    [
        ("EGFR", "ENST00000275493", "C231S", {'TCT': {'start_pos': 55152609, 'reference_allele': 'GC', 'alternative_allele': 'CT'}, 'TCC': {'start_pos': 55152609, 'reference_allele': 'G', 'alternative_allele': 'C'}, 'TCA': {'start_pos': 55152609, 'reference_allele': 'GC', 'alternative_allele': 'CA'}, 'TCG': {'start_pos': 55152609, 'reference_allele': 'GC', 'alternative_allele': 'CG'}, 'AGT': {'start_pos': 55152608, 'reference_allele': 'TGC', 'alternative_allele': 'AGT'}, 'AGC': {'start_pos': 55152608, 'reference_allele': 'T', 'alternative_allele': 'A'}}),
        ("IDH1", "ENST00000345146", "E110K", {'AAA': {'start_pos': 208248455, 'reference_allele': 'C', 'alternative_allele': 'T'}, 'AAG': {'start_pos': 208248453, 'reference_allele': 'TTC', 'alternative_allele': 'CTT'}}),
    ],
)
def test_AA2DNA_var(gene: str, transcript: str, aa_var: str, expected_dna_var: dict) -> None:
    g = get_homo_sapiens_object(gene)
    if (dna_all_muts := g.AA_mut_to_DNA_SNP_mut(aa_var, transcript)) is None:
        raise ValueError(f"Error with {aa_var=} in {transcript=} !!")
    assert dna_all_muts == expected_dna_var

# Mus musculus 
# ============
@pytest.mark.parametrize(
    "gene, expected",
    [
        # mus musculus
        ("ENSMUSG00000017167", 2),
        ("Drg1", 7),
        ("Pik3c2a", 7),
    ],
)
def test_mouse_gene_num_transcripts(gene: str, expected: int) -> None:
    g = get_mus_musculus_object(gene)
    assert len(g) == expected

@pytest.mark.parametrize(
    "gene, expected",
    [
        ("ENSMUSG00000017167", (101_061_349, 101_081_550, '11', False)),
        ("Drg1", (3_137_360, 3_216_415, '11', True)),
        ("Pik3c2a", (115_936_500, 116_042_684, '7', True)),
    ],
)
def test_mouse_gene_start_end_chrm_rev(gene: str, expected: tuple[int, int, str, bool]) -> None:
    g = get_mus_musculus_object(gene)
    assert (g.gene_start, g.gene_end, g.chrm, g.rev) == expected

@pytest.mark.parametrize(
    "gene, transcript, expected_partial_protein_seq",
    [
        # encoded on the positive strand
        ("Cntnap1", "ENSMUST00000103109", "MMSLRLFSILLATVVSGAWGWGYYGCNEELVGPLYARSLGASSYYGLFTTARFARLHGISGWSPRIGDPNPWLQIDLMKKHRIRAVATQGAFNSWDWVTRYMLLYGDRVDSWTPFYQKGHN"),
        # encoded on the negative strand
        ("ENSMUSG00000020457", "ENSMUST00000020741", "MSGTLAKIAEIEAEMARTQKNKATAHHLGLLKARLAKLRRELITPKGGGGGGPGEGFDVAKTGDARIGFVGFPSVGKSTLLSNLAGVYSEVAAYEFTTLTTVPGVIRYKGAKIQLLDLPGIIEGAKDGKGRGRQVIAVARTCNLILIVLDVLKPLGHKKIIENELEGFGIRLNSKPPNIGFKKKDKGGINLTATCPQSELDAETVKSILAEYKIHNADVTLRSDATADDLIDVVEGNRVYIPCIYVLNKIDQISIEELDIIYKVPHCVPISAHHRWNFDDLLEKIWDYLKLVRIYTKPKG"),
        ("Pik3c2a", "ENSMUST00000206219",
        "MAQISNNSEFKQCSSSHPEPIRTKDVNKAEALQMEAEALAKLQKDRQMTDSPRGFELSSSTRQRTQGFNKQDYDLMVFPELDSQKRAVDIDVEKLTQAELEKILLDDNFETRKPPALPVTPVLSPSFSTQLYLRPSGQRGQWPPGLCGPSTYTLPSTYPSAYSKQATFQNGFSPRMPTFPSTESVYLRLPGQSPYFSYPLTPATPFHPQGSLPVYRPLVSPDMAKLFEKIASTSEFLKNGKARTDLEIANSKASVCNLQISPKSEDINKFDW")
    ],
)
def test_mouse_gene_transcript_partial_protein_seq(gene: str, transcript: str, expected_partial_protein_seq: str) -> None:
    g = get_mus_musculus_object(gene)
    if (protein_seq := g.AA(transcript)) is None:
        raise ValueError
    assert protein_seq[:len(expected_partial_protein_seq)].upper() == expected_partial_protein_seq.upper()


# Chromosome sequence fetch
# =========================
@pytest.mark.parametrize(
    "chrm_info, start_p, end_p, rev, species, expected_seq",
    [
        (Path('Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa'), 122_989_200, 122_989_200 + 29, False, 'homo_sapiens', "gcacccactgcgatctcagagcaacctggg"),
        (Path('Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa'), 22_345_678, 22_345_790, True, 'homo_sapiens', "TTTGTTTATTGAGTTGTTTATTTGATTTGTTTATCCTTTCAAATATTCCAGCCTTCAATTTAGAGAATTTGTTTGAAATATATTATAAAAACAGAAGATAGGAGAAAGGAAAA"),
        (Path('Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa'), 60_000_000, 60_000_123, False, 'homo_sapiens', "ACTATCCTAAATCAATTAAGTCAGATCTTTGGGGTATAGAACCCTGGAATTTCTCTTTCAGGGGAAAAAAATGGTTGACTCAGGTGATTCTGATACTGAGCCAGGATGGAGAAGCAAACAAAAT"),
        ('11', 16_341_123, 16_341_143, False, 'mus_musculus', "TCAACATCCTTAATCATCAGG"),
        ('5', 52_120_100, 52_120_120, True, 'Danio_rerio', "TTTTAGAGGTTGAACAGCCAC"),
        ('14', 80_237_279, 80_237_310, True, 'Pan troglodytes', 'AGTAAATCTCCAATAAAGTTATCGTCTGTTCA'),
        ('Y', 1_100_800, 1_100_850, True, 'homo_sapiens', "ACTTACCTCTCTTAAAACCTCATTGGAATAAAACTTAAAAGGAATAAAAAA")
    ]
)
def test_fetch_seq(chrm_info: Path | str, start_p: int, end_p: int, rev: bool, species: str, expected_seq: str) -> None:
    assert u.fetch_seq(chrm_info, start_p, end_p, rev=rev, species=species).upper() == expected_seq.upper()