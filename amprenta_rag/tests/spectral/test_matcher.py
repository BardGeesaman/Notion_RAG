from __future__ import annotations

from amprenta_rag.spectral.matcher import match_spectrum


def test_match_spectrum_precursor_filter_and_ranking():
    # Query peaks overlap strongly with ref1
    query_mz = 760.585
    query_peaks = [(100.0, 10.0), (150.0, 20.0), (200.0, 5.0)]

    refs = [
        {
            "id": "11111111-1111-1111-1111-111111111111",
            "lipid_name": "PC 34:1",
            "precursor_mz": 760.585,
            "spectrum": {"mz": [100.0, 150.0, 300.0], "intensity": [10.0, 20.0, 1.0]},
        },
        {
            "id": "22222222-2222-2222-2222-222222222222",
            "lipid_name": "PE 36:2",
            "precursor_mz": 760.800,  # outside 10ppm tolerance
            "spectrum": {"mz": [100.0, 150.0, 200.0], "intensity": [1.0, 1.0, 1.0]},
        },
        {
            "id": "33333333-3333-3333-3333-333333333333",
            "lipid_name": "PC 34:2",
            "precursor_mz": 760.590,  # within tolerance, but weaker overlap
            "spectrum": {"mz": [100.0, 250.0], "intensity": [10.0, 10.0]},
        },
    ]

    out = match_spectrum(query_mz, query_peaks, refs, precursor_tolerance_ppm=10, fragment_tolerance_da=0.01, top_k=10)
    assert len(out) >= 1
    # PE 36:2 should be filtered out by precursor tolerance
    assert all(m.lipid_name != "PE 36:2" for m in out)
    # Best hit should be PC 34:1 due to stronger peak overlap
    assert out[0].lipid_name == "PC 34:1"
    assert out[0].score >= out[-1].score


