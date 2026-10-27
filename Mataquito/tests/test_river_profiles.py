"""Tests for mataquito.river_profiles — k_sn regression and K solution."""

import numpy as np
import pytest
from mataquito.river_profiles import get_ksn, compute_r2, solve_K


class TestGetKsn:
    def test_recovers_slope(self):
        """A noiseless linear chi-elevation profile returns its slope as k_sn."""
        chi = np.linspace(0, 100, 50)
        elev = 2.5 * chi + 10
        ksn, se = get_ksn(chi, elev, 'synthetic', verbose=False)
        assert ksn == pytest.approx(2.5)
        assert se == pytest.approx(0.0, abs=1e-10)

    def test_too_few_pixels(self):
        ksn, se = get_ksn(np.arange(3.0), np.arange(3.0), 'short', verbose=False)
        assert ksn is None and se is None

    def test_zero_chi_spread(self):
        ksn, se = get_ksn(np.ones(10), np.arange(10.0), 'flat', verbose=False)
        assert ksn is None and se is None


class TestComputeR2:
    def test_perfect_fit(self):
        chi = np.linspace(0, 10, 20)
        assert compute_r2(chi, 3 * chi) == pytest.approx(1.0)

    def test_degenerate_is_nan(self):
        assert np.isnan(compute_r2(np.ones(10), np.arange(10.0)))


class TestSolveK:
    def test_k_value(self):
        """K = U / (k_sn * Ao^m), with U converted from m/Myr to m/yr."""
        res = solve_K(100.0, 10.0, 2.0, 0.0, 'test', ao_m=1e6 ** 0.4, verbose=False)
        expected = 100e-6 / (2.0 * 1e6 ** 0.4)
        assert res['K'] == pytest.approx(expected)
        # k_sn SE is zero, so relative K error equals relative U error
        assert res['sigma_K'] / res['K'] == pytest.approx(0.1)

    def test_quadrature_uncertainty(self):
        res = solve_K(100.0, 10.0, 2.0, 0.2, 'test', ao_m=1.0, verbose=False)
        assert res['sigma_K'] / res['K'] == pytest.approx(np.sqrt(0.1 ** 2 + 0.1 ** 2))

    def test_missing_ksn_returns_nan(self):
        res = solve_K(100.0, 10.0, None, None, 'skip', ao_m=1.0, verbose=False)
        assert np.isnan(res['K']) and np.isnan(res['sigma_K'])
