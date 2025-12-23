from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from amprenta_rag.multi_omics.result_parser import parse_mofa_output


def test_parse_mofa_output_minimal_hdf5(tmp_path: Path):
    h5py = pytest.importorskip("h5py")

    p = tmp_path / "mofa_model.hdf5"
    with h5py.File(p, "w") as h5:
        exp = h5.create_group("expectations")

        # Z: samples x factors
        exp.create_dataset("Z", data=np.array([[0.1, 0.2], [0.3, 0.4]]))

        # W per view: features x factors
        wgrp = exp.create_group("W")
        wgrp.create_dataset("transcriptomics", data=np.array([[1.0, -1.0], [0.5, 0.2]]))

        # Optional names
        h5.create_dataset("samples", data=np.array([b"S1", b"S2"]))
        h5.create_dataset("transcriptomics/features", data=np.array([b"G1", b"G2"]))

        # Optional variance
        h5.create_dataset("variance_explained", data=np.array([[0.1, 0.2]]))

    out = parse_mofa_output(str(p))
    assert out["factors"] == [0, 1]
    assert "scores" in out and out["scores"].shape == (2, 2)
    assert "loadings" in out and "transcriptomics" in out["loadings"]
    assert out["loadings"]["transcriptomics"].shape == (2, 2)


