"""Monte Carlo erosion rate sampling and mass-balance (flux) ordering.

REVISED VERSION -- replaces erosion.py.

What changed relative to erosion.py, and why
--------------------------------------------
1. RETURNS BOTH THE ORDERED AND THE UNORDERED ARRAYS, so the notebook can run
   the whole fertility analysis twice and report the sensitivity of every result
   to the mass-balance constraint.

   Why this matters. Of the five flux-ordering constraints, only ONE actually
   binds on the measured data:

       constraint            upstream flux   downstream flux   retained
       -------------------   -------------   ---------------   --------
       CT-1  -> CT-6           3.17e4           7.67e4          ~100%
       CT-5+CT-6 -> CT-10      6.95e5           1.20e6          ~100%
       CT-5+CT-6 -> CT-11      6.95e5           1.41e6          ~100%
       CT-5+CT-6 -> CT-8       6.95e5           1.30e6          ~100%
       CT-7  -> CT-5           6.43e5           5.79e5          ~21%   <-- binds

   The measured CT-7 and CT-5 rates are mutually inconsistent with downstream-
   increasing flux (E_7*A_7 > E_5*A_5), so the filter discards ~79% of draws and
   shifts CT-5's median from 387 to 413 m/Myr, a +6.8% change. Because CT-5 is
   the fast-eroding parent in every upstream scenario, that 6.8% propagates into
   a ~40-50% change in the quartz fertility and flux ratios. The constraint is
   physically necessary and is retained as the primary analysis, but the size of
   its effect must be reportable, hence the unordered arrays.

   CT-5 is the ONLY array that moves. CT-6 shifts by <0.1%; CT-8, CT-10 and CT-11
   are 100% retained and are numerically identical ordered or not; CT-3 and CT-4
   have no upstream constraint at all. The sensitivity therefore only bites the
   four CT-5/CT-6-parent scenarios.

2. PAIRED TRIPLES. Downstream code needs (E_parent_a, E_parent_b, E_child) as
   arrays of equal length whose index i is the SAME Monte Carlo iteration. The
   joint CT-5+CT-6 -> child filter selects a common index set, so the surviving
   CT-5, CT-6 and child draws stay index-aligned. This module therefore returns
   the aligned triples directly rather than leaving the notebook to slice them.

3. Local `np.random.default_rng(seed)` instead of the global `np.random.seed`.
   Deterministic as before, but no longer mutates global random state. Numbers
   will differ from the legacy implementation (PCG64 vs MT19937 stream).

4. The unused CT-4 -> CT-8 ordering block was removed. Nothing consumed it, and
   it encoded a topology that contradicts the manuscript's stated constraint set.


What "flux" means here -- read this before citing E*A as a mass flux
--------------------------------------------------------------------
A ^10Be-derived catchment-averaged erosion rate is a QUARTZ-FLUX-WEIGHTED mean,
not an area-weighted mean: the nuclide concentration of the sediment is the
quartz-flux-weighted average of its sources. Consequently E*A is a quartz-
weighted flux PROXY, not a bulk sediment mass flux, and the ordering constraint
is a physically motivated PRIOR rather than a statement of mass conservation.

This is why the apparent fluxes at CT-10/CT-11/CT-8.2 exceed the sum of the
parent fluxes by 50-70%: the Teno (CT-5) is both faster-eroding and (per this
study's result) more quartz-fertile, which pulls the apparent downstream erosion
rate upward. That excess is the fertility signal, not a mass-balance violation.
"""

import numpy as np


def generate_mc_samples(df, N=100_000, seed=17):
    """Draw N erosion-rate samples per site from a normal distribution.

    Mean = measured erosion rate; standard deviation = EXTERNAL uncertainty.
    External uncertainty includes production-rate scaling, spallation and
    muogenic contributions (Balco et al., 2008) and is the conservative choice.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'Sample_ID', 'Erosion_rate',
        'Erosion_rate_uncertainty_external'.
    N : int
        Draws per site.
    seed : int
        Seed for the local Generator.

    Returns
    -------
    dict
        {Sample_ID: np.ndarray of shape (N,)}
    """
    rng = np.random.default_rng(seed)
    ER = {}
    for _, row in df.iterrows():
        ER[row['Sample_ID']] = rng.normal(
            loc=row['Erosion_rate'],
            scale=row['Erosion_rate_uncertainty_external'],
            size=N,
        )
    return ER


def flux_order_samples(ER, A, N=100_000, seed=17):
    """Apply the downstream-increasing-flux constraint through the flow network.

    Constraint set (as stated in the manuscript Methods):
        CT-7        -> CT-5
        CT-1        -> CT-6
        CT-5 + CT-6 -> CT-10, CT-11, CT-8

    CT-3 (hydrologically independent Coastal Cordillera) and CT-4 (excluded as a
    poorly mixed sample) carry no upstream constraint. CT-2 and CT-9 are not
    ordered: the Central Depression lies between CT-8 and the outlet, and the
    ordering statement is not meaningful across a reach where deposition and
    bypass occur.

    Procedure. For each Monte Carlo iteration, the flux proxy E*A is evaluated.
    Iterations in which the downstream flux does not exceed the sum of upstream
    fluxes are discarded. Surviving CT-5 and CT-6 draws are resampled with
    replacement back to N before being used in the CT-5+CT-6 -> child step, so
    that the child filter operates on a full-size, correctly-weighted parent
    distribution.

    Parameters
    ----------
    ER : dict
        {Sample_ID: np.ndarray}, from generate_mc_samples().
    A : dict
        {Sample_ID: area_km2}.
    N : int
        Target array length after resampling.
    seed : int
        Seed for the local Generator.

    Returns
    -------
    dict with keys:
        'CT5_ord', 'CT7_ord'   : surviving CT-7 -> CT-5 pairs (index-aligned)
        'CT6_ord', 'CT1_ord'   : surviving CT-1 -> CT-6 pairs (index-aligned)
        'CT5_res', 'CT6_res'   : resampled back to N
        'CT3_res', 'CT4_res'   : resampled (unconstrained; for length parity)
        'CT10_res','CT11_res','CT8_res' : ordered children, resampled to N
        'triples'              : {child: (E_CT5, E_CT6, E_child)} index-aligned
        'retained'             : {constraint_name: fraction retained}
    """
    rng = np.random.default_rng(seed)
    out = {}
    retained = {}

    # -- CT-7 -> CT-5 --------------------------------------------------------
    # This is the constraint that actually binds (~21% retained). The measured
    # values violate it: E_7*A_7 = 6.43e5 > E_5*A_5 = 5.79e5.
    m75 = ER['CT-5'] * A['CT-5'] > ER['CT-7'] * A['CT-7']
    out['CT5_ord'] = ER['CT-5'][m75]
    out['CT7_ord'] = ER['CT-7'][m75]
    retained['CT-7 -> CT-5'] = float(m75.mean())

    # -- CT-1 -> CT-6 --------------------------------------------------------
    m16 = ER['CT-6'] * A['CT-6'] > ER['CT-1'] * A['CT-1']
    out['CT6_ord'] = ER['CT-6'][m16]
    out['CT1_ord'] = ER['CT-1'][m16]
    retained['CT-1 -> CT-6'] = float(m16.mean())

    # -- Resample the ordered parents back up to N ---------------------------
    # CT-5 and CT-6 drain independent catchments and were filtered against
    # different upstream samples, so they are resampled independently.
    CT5_res = rng.choice(out['CT5_ord'], size=N, replace=True)
    CT6_res = rng.choice(out['CT6_ord'], size=N, replace=True)
    out['CT5_res'] = CT5_res
    out['CT6_res'] = CT6_res

    # -- CT-5 + CT-6 -> each child ------------------------------------------
    # The joint filter selects a single index set, so the surviving CT-5, CT-6
    # and child draws remain index-aligned. That alignment is what makes the
    # downstream elementwise fertility calculation a valid Monte Carlo.
    upstream_flux = CT5_res * A['CT-5'] + CT6_res * A['CT-6']
    triples = {}
    for child in ['CT-10', 'CT-11', 'CT-8']:
        m = ER[child] * A[child] >= upstream_flux
        triples[child] = (CT5_res[m], CT6_res[m], ER[child][m])
        retained[f'CT-5 + CT-6 -> {child}'] = float(m.mean())
        # Also expose the ordered child on its own, resampled to N, for use as a
        # PARENT in the downstream CT-x/CT-3 scenarios. Key names drop the hyphen:
        # 'CT-10' -> 'CT10_res', 'CT-8' -> 'CT8_res'.
        key = child.replace('-', '') + '_res'
        out[key] = rng.choice(ER[child][m], size=N, replace=True)

    # CT-4 has no upstream constraint (excluded as poorly mixed), so its "ordered"
    # triple simply pairs the ordered parents against the raw CT-4 draws. This is
    # a change from the legacy notebook, which used UNORDERED CT-5 and CT-6 for
    # CT-4 while every other scenario used ordered ones -- an inconsistency that
    # made CT-4's numbers non-comparable with the rest of the table.
    triples['CT-4'] = (CT5_res, CT6_res, ER['CT-4'])
    retained['CT-4 (unconstrained)'] = 1.0

    out['triples'] = triples
    out['retained'] = retained

    # -- Unconstrained resamples, for length parity in downstream scenarios ---
    out['CT3_res'] = rng.choice(ER['CT-3'], size=N, replace=True)
    out['CT4_res'] = rng.choice(ER['CT-4'], size=N, replace=True)

    return out


def raw_triples(ER):
    """Unordered counterparts of `flux_order_samples()['triples']`.

    Used for the flux-ordering sensitivity test: identical scenario structure,
    but with the raw (unfiltered) Monte Carlo erosion rates. Every array is
    length N and index-aligned by construction, since no filtering occurs.
    """
    return {
        child: (ER['CT-5'], ER['CT-6'], ER[child])
        for child in ['CT-4', 'CT-10', 'CT-11', 'CT-8']
    }
