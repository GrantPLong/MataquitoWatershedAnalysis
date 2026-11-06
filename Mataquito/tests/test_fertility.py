"""Tests for mataquito.fertility — quartz and zircon fertility and flux calculations."""

import numpy as np
import pytest
from mataquito.fertility import (
    THRESHOLD,
    ct3_star_area,
    load_mixing_results,
    quartz_fertility,
    ratio_stats,
    zircon_fertility,
)

MIX_XLSX = "results/mixing_coefficients/dPy-mix_ct-5_ct-6_to_ct-10_results_vmax.xlsx"


def _erosion_draws(n=1000, seed=42):
    rng = np.random.default_rng(seed)
    E_a = rng.normal(400, 50, n)
    E_b = rng.normal(30, 5, n)
    E_c = rng.normal(250, 30, n)
    return E_a, E_b, E_c


class TestLoadMixingResults:
    def test_counts_and_pairing(self):
        """Filtered coefficient arrays are paired and no longer than the ensemble."""
        r = load_mixing_results(MIX_XLSX, "CT-10")
        assert r["parents"] == ["CT-5", "CT-6"]
        assert 0 < r["n_kept"] <= r["n_total"]
        assert len(r["wa"]) == len(r["wb"]) == r["n_kept"]

    def test_coefficients_in_range(self):
        """All kept coefficients lie within (THRESHOLD, 1 - THRESHOLD)."""
        r = load_mixing_results(MIX_XLSX, "CT-10")
        for w in (r["wa"], r["wb"]):
            assert np.all(w > THRESHOLD)
            assert np.all(w < 1 - THRESHOLD)

    def test_best_fit_sums_to_one(self):
        r = load_mixing_results(MIX_XLSX, "CT-10")
        assert r["best_a"] + r["best_b"] == pytest.approx(1.0, abs=1e-3)


class TestCT3StarArea:
    def test_residual_area(self):
        A = {"CT-9": 6000.0, "CT-8": 4800.0}
        assert ct3_star_area(A, "CT-8", "CT-9") == pytest.approx(1200.0)

    def test_parent_larger_than_child_raises(self):
        with pytest.raises(ValueError):
            ct3_star_area({"CT-9": 100.0, "CT-8": 200.0}, "CT-8", "CT-9")


class TestQuartzFertility:
    def test_positive_output_only(self):
        """Only physical (positive) fertility and flux ratios are returned."""
        Q, F, mask = quartz_fertility(*_erosion_draws(), 1500, 2500, 14.6, 5.9, 11.5)
        assert len(Q) > 0
        assert np.all(Q > 0)
        assert np.all(F > 0)

    def test_mask_matches_output(self):
        Q, F, mask = quartz_fertility(*_erosion_draws(), 1500, 2500, 14.6, 5.9, 11.5)
        assert len(Q) == len(F) == mask.sum()

    def test_flux_ratio_independent_of_area(self):
        """The quartz flux ratio does not depend on the source areas."""
        draws = _erosion_draws()
        _, F1, _ = quartz_fertility(*draws, 1500, 2500, 14.6, 5.9, 11.5)
        _, F2, _ = quartz_fertility(*draws, 300, 9000, 14.6, 5.9, 11.5)
        np.testing.assert_allclose(F1, F2)


class TestZirconFertility:
    def _run(self):
        rng = np.random.default_rng(42)
        wa = rng.uniform(0.3, 0.7, 500)
        wb = 1.0 - wa
        ER_A = rng.normal(400, 50, 500)
        ER_B = rng.normal(30, 5, 500)
        return zircon_fertility(wa, wb, ER_A, ER_B, 1500, 2500, num_samples=200, seed=17)

    def test_output_structure(self):
        result = self._run()
        for key in ("ZA_ZB", "ZB_ZA", "FZA_FZB", "FZB_FZA"):
            assert key in result
            assert "median" in result[key]

    def test_reproducible(self):
        assert self._run()["ZA_ZB"]["median"] == self._run()["ZA_ZB"]["median"]


class TestRatioStats:
    def test_log_median_matches_linear(self):
        x = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
        st = ratio_stats(x)
        assert st["median"] == pytest.approx(2.0)
        assert 10 ** st["log10_median"] == pytest.approx(st["median"])
