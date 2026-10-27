"""
river_profiles.py — chi-elevation, k_sn and erodibility functions for notebook 08.

Ported from the Rapel Basin ``rapel_utils.py`` (river-profile subset only; the
tau / knickpoint-segmentation helpers are not included).

Functions only. **Analysis parameters live in the notebook**, not here — theta,
A0 and erosion rates are defined in the notebook setup cell. Nothing in this
module reads a notebook global; everything it needs arrives as an argument.

Contents
--------
Grids and loading
    mean_pixel_dimension_grid, extract_watershed
Network traversal
    fix_indices_recursive, upstream_list
Chi-elevation analysis
    get_ksn, compute_r2, solve_K
Plotting primitives
    plot_upstream_profile, fit_and_draw
"""

import numpy as np
from scipy.stats import linregress

__all__ = [
    'mean_pixel_dimension_grid', 'extract_watershed',
    'fix_indices_recursive', 'upstream_list',
    'get_ksn', 'compute_r2', 'solve_K',
    'plot_upstream_profile', 'fit_and_draw',
]

# Project-wide constants that are not analysis choices.
AREA_THRESHOLD_DEFAULT = 1e6    # m², 1 km² — minimum channel-forming drainage area
MYR_TO_YR = 1e-6                # unit conversion, not a parameter


# ═════════════════════════════════════════════════════════════════════════════
# Grids and loading
# ═════════════════════════════════════════════════════════════════════════════

def mean_pixel_dimension_grid(area):
    """Per-pixel along-step length in metres, as a grid TopoAnalysis can consume.

    Required as the ``de`` argument to ``map_values_to_recursive_list``.
    """
    from TopoAnalysis import dem as d

    mean_px = d.BaseSpatialGrid()
    mean_px._copy_info_from_grid(area, True)
    mean_px._griddata = area._mean_pixel_dimension()
    return mean_px


def extract_watershed(outlet_coords, fd, area, elev_grid, area_grid,
                      theta, ao, area_threshold=AREA_THRESHOLD_DEFAULT):
    """Chi-elevation-area arrays for the watershed upstream of ``outlet_coords``.

    Chi is solved from this outlet alone, so chi = 0 at the outlet and each
    catchment sits on its own axis. Membership is chi > 0 combined with the
    area threshold.

    Returns
    -------
    chi_vals, elev_vals, area_vals : 1-D arrays, NaN pixels removed
    """
    from TopoAnalysis import dem as d

    chi_sample = d.Chi(flow_direction=fd, area=area, theta=theta, Ao=ao,
                       outlets=[outlet_coords])
    mask = (chi_sample._griddata > 0) & (area_grid >= area_threshold)

    chi_vals  = chi_sample._griddata[mask]
    elev_vals = elev_grid[mask]
    area_vals = area_grid[mask]

    valid = ~np.isnan(chi_vals) & ~np.isnan(elev_vals)
    return chi_vals[valid], elev_vals[valid], area_vals[valid]


# ═════════════════════════════════════════════════════════════════════════════
# Network traversal
# ═════════════════════════════════════════════════════════════════════════════

def fix_indices_recursive(ld_item):
    """Unwrap single-element tuple indices to (row, col) form, in place.

    ``map_values_to_recursive_list`` returns some indices as 1-tuples; downstream
    code assumes (row, col).
    """
    if isinstance(ld_item['index'], tuple) and len(ld_item['index']) == 1:
        ld_item['index'] = ld_item['index'][0]
    for next_item in (ld_item.get('next') or []):
        fix_indices_recursive(next_item)


def _as_rc(idx):
    """Normalise a node index to (row, col)."""
    return idx if isinstance(idx, tuple) else (idx, 0)


def upstream_list(outlet_coords, fd, area, elevation, mean_px=None):
    """Recursive upstream network from an outlet, with indices normalised.

    ``mean_px`` is built on demand if not supplied, but building it once in the
    notebook and passing it in avoids repeating that work per outlet.
    """
    if mean_px is None:
        mean_px = mean_pixel_dimension_grid(area)
    ld = fd.map_values_to_recursive_list(
        outlet_coords, elevation=elevation, area=area, de=mean_px)
    fix_indices_recursive(ld)
    return ld


# ═════════════════════════════════════════════════════════════════════════════
# Chi-elevation analysis
# ═════════════════════════════════════════════════════════════════════════════

def get_ksn(chi_vals, elev_vals, label, min_pixels=5, verbose=True):
    """k_sn and its standard error from a chi-elevation regression.

    Guarded on both pixel count and chi spread: ``linregress`` raises if every x
    is identical. Such a reach returns (None, None) rather than crashing, so
    callers must tolerate None.
    """
    if len(chi_vals) < min_pixels or np.ptp(chi_vals) == 0:
        print(f'  {label}: insufficient pixels or zero chi spread')
        return None, None

    slope, intercept, r, _, se = linregress(chi_vals, elev_vals)
    if verbose:
        print(f'  {label}:')
        print(f'    n       = {len(chi_vals):,}')
        print(f'    k_sn    = {slope:.4f}')
        print(f'    R²      = {r ** 2:.4f}')
        print(f'    SE      = {se:.6f}')
    return slope, se


def compute_r2(chi_vals, elev_vals, min_pixels=5):
    """R² of the chi-elevation linear fit; NaN if degenerate.

    Same guard as ``get_ksn`` — too few points or zero chi variance gives NaN.
    """
    if len(chi_vals) < min_pixels or np.ptp(chi_vals) == 0:
        return np.nan
    _, _, r, _, _ = linregress(chi_vals, elev_vals)
    return r ** 2


def solve_K(U_mMyr, sigma_mMyr, ksn_val, ksn_se_val, label, ao_m, verbose=True):
    """Erodibility K from steady-state stream power with n = 1.

        K = U / (k_sn * Ao^m)

    ``ao_m`` is Ao**theta and must be passed explicitly — it differs between
    concavities.

    Uncertainty is propagated in quadrature from the erosion rate and the k_sn
    standard error. A reach skipped by ``get_ksn`` arrives with ksn=None and
    returns NaN rather than raising.
    """
    if ksn_val is None or ksn_se_val is None:
        print(f'  {label}: skipped — no k_sn')
        return {'K': np.nan, 'sigma_K': np.nan,
                'U': U_mMyr * MYR_TO_YR, 'ksn': np.nan}

    U   = U_mMyr * MYR_TO_YR
    sig = sigma_mMyr * MYR_TO_YR
    K   = U / (ksn_val * ao_m)

    rel_err_U   = sig / U
    rel_err_ksn = ksn_se_val / ksn_val
    sigma_K = K * np.sqrt(rel_err_U ** 2 + rel_err_ksn ** 2)

    if verbose:
        print(f'  {label}:')
        print(f'    U       = {U_mMyr} m/Myr')
        print(f'    k_sn    = {ksn_val:.4f}')
        print(f'    Ao^m    = {ao_m:.4e}')
        print(f'    K       = {K:.4e} m^(1-2m)/yr')
        print(f'    sigma_K = {sigma_K:.4e}')

    return {'K': K, 'sigma_K': sigma_K, 'U': U, 'ksn': ksn_val}


# ═════════════════════════════════════════════════════════════════════════════
# Plotting primitives
# ═════════════════════════════════════════════════════════════════════════════

def plot_upstream_profile(ax, ld_list, minimum_area=AREA_THRESHOLD_DEFAULT,
                          color='k', starting_distance=0.0, linewidth=0.8):
    """Draw a distance-elevation profile of the whole upstream network.

    Walks every branch above ``minimum_area``, accumulating real-space distance
    (D8 step length, x1.414 on diagonals).
    """
    def plot_recursive(current_distance, ld_item):
        for next_item in (ld_item.get('next') or []):
            if next_item['area'] < minimum_area:
                continue
            r0, c0 = _as_rc(ld_item['index'])
            r1, c1 = _as_rc(next_item['index'])
            diag = (r0 != r1) and (c0 != c1)
            next_distance = current_distance + ld_item['de'] * (1.414 if diag else 1.0)
            ax.plot([current_distance / 1000, next_distance / 1000],
                    [ld_item['elevation'], next_item['elevation']],
                    color=color, linewidth=linewidth)
            plot_recursive(next_distance, next_item)

    plot_recursive(starting_distance, ld_list)


def fit_and_draw(ax, chi_v, elev_v, color, label, linewidth=1):
    """Fit elevation on chi, draw the regression line, return (slope, r2).

    Returns (None, None) if there are too few points to fit.
    """
    if len(chi_v) < 2:
        return None, None
    slope, intercept, r, _, _ = linregress(chi_v, elev_v)
    x = np.array([chi_v.min(), chi_v.max()])
    ax.plot(x, slope * x + intercept, color=color, linewidth=linewidth,
            label=f'{label}  k_sn={slope:.3f}  R²={r ** 2:.3f}')
    return slope, r ** 2
