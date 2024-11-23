from pathlib import Path
import pytest
import geneanot.Gene_annotation as ga

@pytest.mark.parametrize(
        "gene, expected",
        [
            # homo_sapiens
            ('EGFR', 13),
            ('IDH1', 9)

        ],
)

def test_gene_num_transcripts(gene: str, expected: int) -> None:
    annotation_full_file: Path = Path('AnnotationDB/Homo_sapiens.GRCh38.113.gff3.gz')
    g = ga.Gene_cls(gene, annotation_full_file, species='homo_spaiens', verbose=False)
    assert len(g) == expected