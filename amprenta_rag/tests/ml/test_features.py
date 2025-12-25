from __future__ import annotations

import pytest

from amprenta_rag.ml.admet.features import DESCRIPTOR_NAMES, FEATURE_NAMES, get_feature_name


def test_feature_names_length() -> None:
    assert len(FEATURE_NAMES) == 2054


def test_morgan_bit_names() -> None:
    assert FEATURE_NAMES[0] == "MorganBit_0000"
    assert FEATURE_NAMES[2047] == "MorganBit_2047"


def test_descriptor_names() -> None:
    assert FEATURE_NAMES[-6:] == [
        "MolWt",
        "MolLogP",
        "TPSA",
        "NumHDonors",
        "NumHAcceptors",
        "NumRotatableBonds",
    ]
    assert DESCRIPTOR_NAMES == FEATURE_NAMES[-6:]


def test_get_feature_name_bounds() -> None:
    assert get_feature_name(0) == "MorganBit_0000"
    assert get_feature_name(2053) == "NumRotatableBonds"


def test_get_feature_name_out_of_range() -> None:
    with pytest.raises(IndexError):
        get_feature_name(2054)


