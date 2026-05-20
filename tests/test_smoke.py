"""Smoke tests for the limma core linear-model workflow."""
from __future__ import annotations

import numpy as np
import pytest

import pylimma


def _make_data(seed=0, n_genes=400, n_per_group=5):
    """Synthetic log-expression matrix with planted DE genes."""
    rng = np.random.default_rng(seed)
    n = 2 * n_per_group
    base = rng.normal(8.0, 2.0, n_genes)
    is_de = rng.uniform(size=n_genes) < 0.15
    delta = np.where(is_de, rng.choice([-1, 1], n_genes) * 2.0, 0.0)
    expr = np.zeros((n_genes, n))
    for g in (0, 1):
        mu = base + (delta if g else 0.0)
        expr[:, g * n_per_group:(g + 1) * n_per_group] = rng.normal(
            mu[:, None], 0.6, (n_genes, n_per_group))
    design = np.zeros((n, 2))
    design[:, 0] = 1.0
    design[n_per_group:, 1] = 1.0
    de_genes = {f"gene_{i}" for i in np.flatnonzero(is_de)}
    return expr, design, de_genes


def test_lmfit_ebayes_toptable():
    expr, design, de_genes = _make_data()
    fit = pylimma.lmFit(expr, design)
    fit = pylimma.eBayes(fit)
    res = pylimma.topTable(fit, coef=1, number=np.inf)
    assert len(res) == expr.shape[0]
    # planted DE genes should dominate the top of the ranking
    top = set(res["gene"].iloc[:len(de_genes)])
    assert len(top & de_genes) / max(len(de_genes), 1) > 0.6


def test_voom_pipeline_runs():
    rng = np.random.default_rng(1)
    counts = rng.poisson(200, (300, 8)).astype(float)
    design = np.zeros((8, 2)); design[:, 0] = 1.0; design[4:, 1] = 1.0
    v = pylimma.voom(counts, design)
    # voom returns an EList — feed its expression matrix + weights to lmFit.
    fit = pylimma.eBayes(pylimma.lmFit(v.E, design, weights=v.weights))
    res = pylimma.topTable(fit, coef=1, number=np.inf)
    assert len(res) == 300


def test_contrasts_fit():
    expr, design, _ = _make_data(seed=2)
    fit = pylimma.lmFit(expr, design)
    contrast = np.array([[0.0], [1.0]])
    fit2 = pylimma.eBayes(pylimma.contrasts_fit(fit, contrast))
    res = pylimma.topTable(fit2, coef=0, number=5)
    assert len(res) == 5


def test_rstyle_aliases_exist():
    for nm in ("lmFit", "eBayes", "topTable", "decideTests",
               "duplicateCorrelation", "removeBatchEffect"):
        assert hasattr(pylimma, nm), nm
