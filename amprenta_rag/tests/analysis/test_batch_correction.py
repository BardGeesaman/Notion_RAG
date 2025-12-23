from __future__ import annotations

import pandas as pd
import pytest


def test_correct_batch_effects_df_validation():
    from amprenta_rag.analysis.batch_correction import correct_batch_effects_df

    df = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]}, index=["f1", "f2"])
    with pytest.raises(ValueError, match="batch_labels length"):
        correct_batch_effects_df(df, batch_labels=["x"], method="combat")


def test_correct_batch_effects_df_calls_combat(monkeypatch):
    from amprenta_rag.analysis import batch_correction as bc

    df = pd.DataFrame(
        {"s1": [1.0, 2.0], "s2": [10.0, 20.0], "s3": [1.2, 2.2], "s4": [9.5, 19.0]},
        index=["f1", "f2"],
    )
    batch = ["A", "B", "A", "B"]

    def fake_combat(m, batch):  # noqa: ANN001
        # simple centering by batch mean (not real ComBat), for test determinism
        out = m.copy()
        for b in sorted(set(batch)):
            cols = [c for c, lab in zip(m.columns, batch) if lab == b]
            out[cols] = out[cols] - out[cols].mean(axis=1).to_numpy()[:, None]
        return out

    monkeypatch.setattr(bc, "_combat_correct", fake_combat)
    corrected, stats = bc.correct_batch_effects_df(df, batch_labels=batch, method="combat")
    assert corrected.shape == df.shape
    assert stats["method"] == "combat"
    assert stats["batches"] == {"A": 2, "B": 2}


