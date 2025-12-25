from __future__ import annotations

import numpy as np
import pytest

xgb = pytest.importorskip("xgboost")
sklearn = pytest.importorskip("sklearn")
_ = (xgb, sklearn)


def test_in_domain_check():
    from amprenta_rag.ml.admet.applicability import ApplicabilityChecker

    X_train = np.ones((50, 32), dtype=np.float32)
    chk = ApplicabilityChecker(threshold=0.3).fit(X_train)

    X = np.ones((5, 32), dtype=np.float32)
    in_domain, sim = chk.check(X)
    assert in_domain.dtype == bool
    assert sim.shape == (5,)
    assert in_domain.all()
    assert (sim >= 0.99).all()


def test_out_of_domain_check():
    from amprenta_rag.ml.admet.applicability import ApplicabilityChecker

    X_train = np.ones((50, 32), dtype=np.float32)
    chk = ApplicabilityChecker(threshold=0.3).fit(X_train)

    X = np.zeros((5, 32), dtype=np.float32)
    in_domain, sim = chk.check(X)
    assert (~in_domain).all()
    assert (sim <= 1e-6).all()


def test_widen_ci_ood_only():
    from amprenta_rag.ml.admet.applicability import ApplicabilityChecker

    chk = ApplicabilityChecker(threshold=0.3)

    ci_low = np.array([0.0, 0.0, 0.0], dtype=float)
    ci_high = np.array([1.0, 1.0, 1.0], dtype=float)
    in_domain = np.array([True, False, True])

    new_low, new_high = chk.widen_ci(ci_low, ci_high, in_domain)
    # only index 1 widens by 2x width: [0,1] -> [-0.5,1.5]
    assert new_low[0] == 0.0 and new_high[0] == 1.0
    assert new_low[2] == 0.0 and new_high[2] == 1.0
    assert new_low[1] == -0.5 and new_high[1] == 1.5


