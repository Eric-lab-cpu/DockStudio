"""Unit tests for enrichment ROC/AUC/EF (v2.0 feature #16)."""

import numpy as np

from dockstudio.core import enrich


def test_roc_auc_perfect_separation():
    """Perfectly separated actives/decoys must give AUC = 1."""
    # more negative affinity = better; actives are much better than decoys
    active = [-9.5, -9.1, -8.8, -8.2]
    decoy = [-6.0, -5.5, -5.0, -4.4]
    res = enrich.roc_auc_ef(active, decoy)
    assert res["auc"] == pytest_approx(1.0)
    assert res["n_active"] == 4 and res["n_decoy"] == 4


def test_roc_auc_reverse_is_zero():
    active = [-4.4, -5.0, -5.5, -6.0]
    decoy = [-8.2, -8.8, -9.1, -9.5]
    res = enrich.roc_auc_ef(active, decoy)
    assert res["auc"] == pytest_approx(0.0)


def test_roc_auc_half_random():
    rng = np.random.default_rng(0)
    active = list(rng.normal(-7.0, 1.0, 60))
    decoy = list(rng.normal(-7.0, 1.0, 60))
    res = enrich.roc_auc_ef(active, decoy)
    assert 0.25 < res["auc"] < 0.75  # near-random but deterministic
    assert len(res["roc_points"]) >= 2


def test_roc_points_monotone():
    active = [-9.0, -8.0, -7.0]
    decoy = [-6.0, -5.0]
    res = enrich.roc_auc_ef(active, decoy)
    # ROC points start at (0,0)
    assert res["roc_points"][0] == (0.0, 0.0)
    # last point should reach (1,1)
    assert res["roc_points"][-1][0] == pytest_approx(1.0)
    assert res["roc_points"][-1][1] == pytest_approx(1.0)


def pytest_approx(v, rel=1e-6):
    # small local helper keeps the test file dependency-light
    import pytest
    return pytest.approx(v, rel=rel)
