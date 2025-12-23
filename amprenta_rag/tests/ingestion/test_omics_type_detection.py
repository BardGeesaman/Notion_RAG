from __future__ import annotations

from amprenta_rag.ingestion import omics_type_detection as otd

def test_detect_omics_type_from_filename():
    assert otd.detect_omics_type("test_lipidomics.csv") == "lipidomics"
    assert otd.detect_omics_type("my_metabolites.csv") == "metabolomics"
    assert otd.detect_omics_type("protein_data.tsv") == "proteomics"
    assert otd.detect_omics_type("rnaseq_genes.csv") == "transcriptomics"

def test_detect_omics_type_from_header(tmp_path):
    # Lipidomics
    p1 = tmp_path / "unknown_1.csv"
    p1.write_text("Lipid Name,Value\nC16:0,1.2")
    assert otd.detect_omics_type(str(p1)) == "lipidomics"

    # Metabolomics
    p2 = tmp_path / "unknown_2.csv"
    p2.write_text("Compound ID,Intensity\nC001,100")
    assert otd.detect_omics_type(str(p2)) == "metabolomics"

    # Proteomics
    p3 = tmp_path / "unknown_3.csv"
    p3.write_text("UniProt ID,Abundance\nP12345,50")
    assert otd.detect_omics_type(str(p3)) == "proteomics"

    # Transcriptomics
    p4 = tmp_path / "unknown_4.csv"
    p4.write_text("Gene Symbol,Log2FoldChange\nTP53,2.5")
    assert otd.detect_omics_type(str(p4)) == "transcriptomics"

def test_detect_omics_type_failure(tmp_path):
    # Empty file
    p = tmp_path / "empty.csv"
    p.touch()
    assert otd.detect_omics_type(str(p)) is None

    # No matches
    p2 = tmp_path / "random.csv"
    p2.write_text("Col1,Col2\n1,2")
    assert otd.detect_omics_type(str(p2)) is None

def test_detect_omics_type_read_error(monkeypatch):
    # Passing a directory path raises IsADirectoryError on open
    assert otd.detect_omics_type(".") is None
