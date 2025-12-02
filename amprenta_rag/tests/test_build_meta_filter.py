# tests/test_build_meta_filter.py

from amprenta_rag.query.rag_query_engine import build_meta_filter


def test_build_meta_filter_single_fields():
    filt = build_meta_filter(
        disease="ALS",
        target="SPTLC1",
        lipid=None,
        signature="ALS-CSF-Core-6Ceramides",
    )

    assert filt["diseases"] == {"$in": ["ALS"]}
    assert filt["targets"] == {"$in": ["SPTLC1"]}
    assert filt["lipid_signatures"] == {"$in": ["ALS-CSF-Core-6Ceramides"]}


def test_build_meta_filter_lipid_or_raw():
    filt = build_meta_filter(
        disease=None,
        target=None,
        lipid="Cer(d18:1/16:0)",
        signature=None,
    )

    assert "$or" in filt
    ors = filt["$or"]
    assert {"lipids": {"$in": ["Cer(d18:1/16:0)"]}} in ors
    assert {"lipids_raw": {"$in": ["Cer(d18:1/16:0)"]}} in ors