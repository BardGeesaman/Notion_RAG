# tests/test_metadata_classifier.py

from amprenta_rag.metadata.classify_literature import (
    _build_notion_updates_from_classification, _classify_literature_doc)


def test_classifier_separates_targets_and_lipids():
    # Synthetic doc that clearly mixes gene targets and lipid terms
    title = "SPTLC1 variants and ceramide dysregulation in ALS"
    abstract = (
        "SPTLC1 and DEGS1 are key enzymes in sphingolipid metabolism. "
        "Ceramide and sphingomyelin levels are altered in ALS patient CSF. "
        "We explore SPTLC1 gain-of-function variants and their impact on ceramide levels."
    )

    cls = _classify_literature_doc(title, abstract)
    updates = _build_notion_updates_from_classification(cls)

    targets_prop = updates["Targets"]["multi_select"]
    lipids_prop = updates["Lipid Species (raw)"]["multi_select"]

    targets = {t["name"] for t in targets_prop}
    lipids = {l["name"] for l in lipids_prop}

    # Targets: only molecular targets (genes/proteins/etc.)
    # e.g., SPTLC1, DEGS1 – *not* "ceramide" or "sphingomyelin"
    assert not any("ceramide" == t.lower() for t in targets)
    assert not any("sphingomyelin" == t.lower() for t in targets)

    # Lipids: contain lipid terms
    lipid_str = " ".join(l.lower() for l in lipids)
    assert ("ceramide" in lipid_str) or ("sphingomyelin" in lipid_str)
