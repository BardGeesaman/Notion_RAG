"""
Integration tests for batch omics ingestion and type detection.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List

import pytest
from unittest.mock import patch

from amprenta_rag.ingestion.omics_type_detection import detect_omics_type as simple_detect_omics_type
from scripts import batch_ingest_omics


# ---------------------------------------------------------------------------
# 1. Omics type detection (amprenta_rag.ingestion.omics_type_detection)
# ---------------------------------------------------------------------------


def test_detect_omics_type_filename_patterns(tmp_path):
    """Filename-only heuristics should detect omics type without reading headers."""
    files = {
        "sample_lipid_data.csv": "lipidomics",
        "study_metabol_results.csv": "metabolomics",
        "prote_quant.tsv": "proteomics",
        "transcript_counts.txt": "transcriptomics",
    }

    for filename, expected in files.items():
        file_path = tmp_path / filename
        file_path.write_text("dummy\n", encoding="utf-8")
        detected = simple_detect_omics_type(str(file_path))
        assert detected == expected


@pytest.mark.parametrize(
    "headers,expected",
    [
        (["Lipid Species", "Intensity"], "lipidomics"),
        (["Metabolite", "Compound", "Area"], "metabolomics"),
        (["Protein", "UniProt", "LFQ"], "proteomics"),
        (["Gene", "log2FoldChange", "padj"], "transcriptomics"),
    ],
)
def test_detect_omics_type_header_inspection(tmp_path, headers, expected):
    """
    Header-based heuristics should detect type when filename is ambiguous.
    """
    file_path = tmp_path / "ambiguous.csv"
    file_path.write_text(",".join(headers) + "\n", encoding="utf-8")

    detected = simple_detect_omics_type(str(file_path))
    assert detected == expected


def test_detect_omics_type_returns_none_for_ambiguous(tmp_path):
    """Ambiguous files without clear signals should return None."""
    file_path = tmp_path / "data.csv"
    file_path.write_text("id,value\n1,2\n", encoding="utf-8")

    detected = simple_detect_omics_type(str(file_path))
    assert detected is None


def test_detect_omics_type_with_real_sample_headers():
    """
    Sanity check against real sample files in the repository.

    Uses existing test data files where available; if files are missing,
    the test is skipped rather than failing.
    """
    root = Path(__file__).parent.parent.parent.parent

    candidates = [
        (root / "sample_lipidomics_data.csv", "lipidomics"),
        (root / "test_data" / "test_lipidomics_features.csv", "lipidomics"),
        (root / "test_data" / "test_metabolomics_features.csv", "metabolomics"),
        (root / "test_data" / "test_proteomics_features.csv", "proteomics"),
        (root / "test_data" / "test_transcriptomics_features.csv", "transcriptomics"),
    ]

    existing: List[Dict[str, str]] = [
        {"path": str(path), "expected": expected}
        for path, expected in candidates
        if path.exists()
    ]

    if not existing:
        pytest.skip("No real sample header files found in repository")

    for item in existing:
        detected = simple_detect_omics_type(item["path"])
        assert detected is not None


# ---------------------------------------------------------------------------
# 2. Batch ingestion orchestrator (scripts.batch_ingest_omics)
# ---------------------------------------------------------------------------


def _write_csv(path: Path, headers: List[str]):
    path.write_text(",".join(headers) + "\n", encoding="utf-8")


def test_batch_ingest_mixed_omics_parallel(tmp_path, capsys):
    """
    Mixed directory with multiple omics types should:
    - Detect correct type for each file via batch_ingest_omics.detect_omics_type
    - Process files in parallel without errors
    - Produce a summary table in stdout
    """
    # Create mixed omics files
    files_and_types = [
        ("lipid_1.csv", "lipidomics", ["Lipid Species", "Intensity"]),
        ("lipid_2.csv", "lipidomics", ["Lipid", "Abundance"]),
        ("metab_1.csv", "metabolomics", ["Metabolite", "Area"]),
        ("metab_2.csv", "metabolomics", ["Compound", "Peak"]),
        ("prote_1.csv", "proteomics", ["Protein", "UniProt"]),
        ("transcript_1.csv", "transcriptomics", ["Gene", "log2FoldChange"]),
    ]

    for name, _, headers in files_and_types:
        _write_csv(tmp_path / name, headers)

    # Stub ingestion functions to avoid real external work
    def _fake_ingest(file_path: str, **kwargs):
        return f"PAGE-{Path(file_path).stem}"

    argv = [
        "batch_ingest_omics",
        "--directory",
        str(tmp_path),
        "--parallel",
        "--max-workers",
        "4",
    ]

    with patch.object(sys, "argv", argv), patch(
        "scripts.batch_ingest_omics.ingest_lipidomics_file", side_effect=_fake_ingest
    ) as mock_lipid, patch(
        "scripts.batch_ingest_omics.ingest_metabolomics_file", side_effect=_fake_ingest
    ) as mock_metab, patch(
        "scripts.batch_ingest_omics.ingest_proteomics_file", side_effect=_fake_ingest
    ) as mock_prot, patch(
        "scripts.batch_ingest_omics.ingest_transcriptomics_file", side_effect=_fake_ingest
    ) as mock_tx, patch(
        "scripts.batch_ingest_omics.sys.exit"
    ) as mock_exit:
        batch_ingest_omics.main()
        captured = capsys.readouterr()

    # No non-zero exit expected
    mock_exit.assert_not_called()

    # Each type-specific ingestion function should have been called at least once
    assert mock_lipid.call_count == 2
    assert mock_metab.call_count == 2
    assert mock_prot.call_count == 1
    assert mock_tx.call_count == 1

    # Summary header and per-type lines should appear
    stdout = captured.out
    assert "BATCH INGESTION SUMMARY" in stdout
    assert "LIPIDOMICS" in stdout.upper()
    assert "METABOLOMICS" in stdout.upper()
    assert "PROTEOMICS" in stdout.upper()
    assert "TRANSCRIPTOMICS" in stdout.upper()


def test_batch_ingest_error_handling_continues_and_summarizes(tmp_path, capsys):
    """
    When one file raises an ingestion error:
    - Batch continues processing other files
    - Error is captured and summarized
    - Exit code is non-zero
    """
    ok_file = tmp_path / "lipid_ok.csv"
    error_file = tmp_path / "lipid_fail.csv"
    _write_csv(ok_file, ["Lipid Species", "Intensity"])
    _write_csv(error_file, ["Lipid Species", "Intensity"])

    def _ok_ingest(file_path: str, **kwargs):
        return f"PAGE-{Path(file_path).stem}"

    def _fail_ingest(file_path: str, **kwargs):
        raise RuntimeError("Simulated ingestion error")

    argv = [
        "batch_ingest_omics",
        "--directory",
        str(tmp_path),
        "--parallel",
        "--max-workers",
        "2",
    ]

    exit_codes: List[int] = []

    def _fake_exit(code: int):
        exit_codes.append(code)
        # Do not actually exit during tests
        return None

    def _side_effect(file_path: str, **kwargs):
        # Route calls based on filename to simulate one success and one failure
        if Path(file_path).name == ok_file.name:
            return _ok_ingest(file_path, **kwargs)
        return _fail_ingest(file_path, **kwargs)

    with patch.object(sys, "argv", argv), patch(
        "scripts.batch_ingest_omics.ingest_lipidomics_file",
        side_effect=_side_effect,
    ), patch(
        "scripts.batch_ingest_omics.sys.exit", side_effect=_fake_exit
    ):
        batch_ingest_omics.main()
        captured = capsys.readouterr()

    # Should have attempted to exit with non-zero code due to failures
    assert exit_codes, "Expected sys.exit to be called"
    assert exit_codes[-1] == 1

    stdout = captured.out
    # Summary should include both successful and failed sections
    assert "Successfully Ingested" in stdout
    assert "Failed Files" in stdout
    assert "Simulated ingestion error" in stdout


def test_manual_type_override_uses_forced_type(tmp_path, capsys):
    """
    --omics-type flag should override automatic detection.

    Verify that detect_omics_type is not called and that the specified type is used.
    """
    file_path = tmp_path / "ambiguous_data.csv"
    _write_csv(file_path, ["id", "value"])

    argv = [
        "batch_ingest_omics",
        "--file",
        str(file_path),
        "--omics-type",
        "metabolomics",
    ]

    def _fake_ingest_metabolomics(file_path: str, **kwargs):
        return f"PAGE-METAB-{Path(file_path).stem}"

    with patch.object(sys, "argv", argv), patch(
        "scripts.batch_ingest_omics.ingest_metabolomics_file",
        side_effect=_fake_ingest_metabolomics,
    ) as mock_metab, patch(
        "scripts.batch_ingest_omics.detect_omics_type"
    ) as mock_detect, patch(
        "scripts.batch_ingest_omics.sys.exit"
    ) as mock_exit:
        batch_ingest_omics.main()
        captured = capsys.readouterr()

    # Auto-detection should not be invoked when --omics-type is provided
    mock_detect.assert_not_called()

    # Metabolomics ingestion should have been used
    assert mock_metab.call_count == 1

    # No failure exit expected
    mock_exit.assert_not_called()
    assert "METABOLOMICS" in captured.out.upper()


