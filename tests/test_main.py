from pathlib import Path
import pytest
import geneanot.Gene_annotation as ga

@pytest.mark.parametrize(
    "gene, expected",
    [
        # homo_sapiens
        ("EGFR", 13),
        ("IDH1", 9),
    ],
)
def test_gene_num_transcripts(gene: str, expected: int) -> None:
    annotation_full_file: Path = Path('AnnotationDB/Homo_sapiens.GRCh38.113.gff3.gz')
    g = ga.Gene_cls(gene, annotation_full_file, species='homo_spaiens', verbose=False)
    assert len(g) == expected


@pytest.mark.parametrize(
    "gene, transcript, expected_partial_protein_seq",
    [
        ("EGFR", "ENST00000275493", "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKL"),
        ("AEBP1", "ENST00000223357", "MAAVRGAPLLSCLLALLALCPGGRPQTVLTDDEIEEFLEGFLSELEPEPREDDVEAPPPPEPTPRVRKAQAGGKPGKRPGTAAE"),
    ],
)
def test_gene_transcript_protein_seq(gene: str, transcript: str, expected_partial_protein_seq: str) -> None:
    annotation_full_file: Path = Path("AnnotationDB/Homo_sapiens.GRCh38.113.gff3.gz")
    g = ga.Gene_cls(gene, annotation_full_file, species="homo_sapiens", verbose=False)
    #g.chrm_fasta_file = (f"Chromosome/Homo_sapiens.GRCh38.dna_sm.chromosome.{g.chrm}.fa")
    if (protein_seq := g.AA(transcript)) is None:
        raise ValueError
    assert protein_seq[:len(expected_partial_protein_seq)].upper() == expected_partial_protein_seq.upper()
