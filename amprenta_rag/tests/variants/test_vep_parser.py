from __future__ import annotations

from pathlib import Path

import pytest

from amprenta_rag.variants.vep_parser import parse_vep_tsv


def test_parse_vep_tsv_extracts_location_consequence_and_predictions(tmp_path: Path):
    p = tmp_path / "vep.tsv"
    p.write_text(
        "\n".join(
            [
                "Location\tAllele\tGene\tSYMBOL\tConsequence\tIMPACT\tSIFT\tPolyPhen\tgnomAD_AF\tCADD_PHRED",
                "1:12345\tG\tENSG000001\tTP53\tmissense_variant\tMODERATE\tdeleterious(0.02)\tpossibly_damaging(0.8)\t0.001\t20.5",
                "",
            ]
        ),
        encoding="utf-8",
    )
    out = parse_vep_tsv(str(p))
    assert len(out) == 1
    r = out[0]
    assert r["chromosome"] == "1"
    assert r["position"] == 12345
    assert r["alt_allele"] == "G"
    assert r["consequence"] == "missense_variant"
    assert r["impact"] == "MODERATE"
    assert r["sift_prediction"] == "deleterious"
    assert r["polyphen_prediction"] == "possibly_damaging"
    assert r["gnomad_af"] == pytest.approx(0.001)
    assert r["cadd_phred"] == pytest.approx(20.5)
    assert r["gene_symbol"] == "TP53"


