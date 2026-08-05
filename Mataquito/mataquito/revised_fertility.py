"""Quartz and zircon mineral fertility + mineral flux calculations.

REVISED VERSION -- replaces fertility.py.

What changed relative to fertility.py, and why
---------------------------------------------
1. MIXING COEFFICIENTS ARE NOW READ FROM THE detritalPy-mix RESULTS WORKBOOK
   (`dPy-mix_*_results_{metric}.xlsx`) rather than from `mix_coeffs_all_*.csv`.

   Reason: the CSVs are ambiguous. The batch script writes
   `mix_coeffs_all_ct-5_ct-6_to_ct-8.2_{metric}.csv` from BOTH the combined
   four-child scenario AND the single-child scenario. Those are different
   bootstrap ensembles (different md5-derived seeds), and whichever scenario ran
   last silently overwrote the other. Nothing on disk records which one won.

   The results workbook is unambiguous (one file per scenario per metric) and it
   additionally contains, at full precision:
       - sheet '<CHILD>'   : the full bootstrap ensemble (one row per iteration)
       - sheet 'Best-fit'  : the deterministic best-fit mixing coefficients
       - sheet 'Model_fit' : best-fit objective value, crit, and two fit
                             diagnostics (see load_mixing_results docstring)
   Reading it therefore removes BOTH hard-coded blocks that used to live in the
   notebook (the best-fit lists on the violin plots, and the `fit_results` dict
   used by the goodness-of-fit gate). Those hard-coded values had already gone
   stale: the notebook carried vmax crit=0.081 for CT-4 while the current run
   reports crit=0.110.

2. NEW: MINERAL FLUX RATIOS, reported separately from FERTILITY RATIOS.

   These are different physical quantities and the manuscript conflated them.

       fertility ratio  Q_A/Q_B  = concentration of quartz per unit rock in A
                                   relative to B.  Dimensionless.  Says nothing
                                   about how much quartz each tributary delivers.

       flux ratio    F_QA/F_QB   = (Q_A * E_A * A_A) / (Q_B * E_B * A_B)
                                 = the ratio of quartz MASS DELIVERED per unit
                                   time.  THIS is "how much more quartz comes
                                   from A than from B".

   The zircon analogues are Z_A/Z_B (fertility) and F_ZA/F_ZB (flux).

3. `quartz_fertility` no longer takes A_c. It never used it -- the child area
   cancels out of the mass balance -- and carrying it in the signature implied
   otherwise.

4. Random number generation uses a local `np.random.default_rng(seed)` instead
   of `np.random.seed(seed)`. The old code reseeded the GLOBAL numpy RNG on every
   single call to `zircon_fertility` (nine times per notebook run), silently
   resetting the random state for anything downstream. Results are still fully
   deterministic, but the numbers differ from the legacy implementation because
   the PCG64 generator produces a different stream than the legacy MT19937.

5. NEW: `ct3_star_area()` centralises the CT-3* area convention so it cannot
   drift between scenarios (it previously did -- three scenarios used the raw
   CT-3 watershed area, a ~7x error relative to the intended convention).


The equations
-------------
Quartz (manuscript Eq. 5). From the flux-weighted nuclide mass balance, assuming
steady state, negligible radioactive decay, and constant attenuation length and
rock density across catchments:

    Q_A / Q_B = [ A_B * (P_B*E_C - P_C*E_B) ] / [ A_A * (P_C*E_A - P_A*E_C) ]

Quartz flux ratio. Substituting the above into F_QA/F_QB = Q_A E_A A_A / (Q_B E_B A_B),
the AREAS CANCEL EXACTLY:

    F_QA / F_QB = [ E_A * (P_B*E_C - P_C*E_B) ] / [ E_B * (P_C*E_A - P_A*E_C) ]

This is worth noting: the quartz flux ratio is INDEPENDENT of the source areas,
and is therefore immune to the "child basin area != A_A + A_B" approximation.

Zircon (manuscript Eq. 11):

    Z_A / Z_B = (w_z1 / w_z2) * (E_B * A_B) / (E_A * A_A)

Zircon flux ratio. By the definition of the mixing coefficient (the fraction of
the child's zircon supplied by each parent):

    F_ZA / F_ZB = Z_A E_A A_A / (Z_B E_B A_B) = w_z1 / w_z2

i.e. the zircon flux ratio IS the mixing coefficient ratio, with no erosion rates
and no areas involved at all. Like the quartz flux ratio, it is immune to the
area approximation.


Sign flips (important)
----------------------
Write the quartz ratios in terms of two shared kernels:

    num_k = P_B*E_C - P_C*E_B          (vanishes when E_C/E_B -> P_C/P_B)
    den_k = P_C*E_A - P_A*E_C          (vanishes when E_A/E_C -> P_A/P_C)

    Q_A/Q_B     = (A_B * num_k) / (A_A * den_k)
    F_QA/F_QB   = (E_A * num_k) / (E_B * den_k)

EITHER kernel can straddle zero within the Monte Carlo cloud, and which one does
differs by scenario. A draw for which either kernel changes sign yields a NEGATIVE
(unphysical) ratio and is masked out. Because Q and F share both kernels and
differ only by the strictly positive factor (E_A*A_A)/(E_B*A_B), they always share
a sign -- one mask serves both.

Measured distance of each kernel from zero, in units of its own Monte Carlo
standard deviation (parents CT-5 / CT-6):

    child     num_k         den_k         draws retained
    -------   -----------   -----------   --------------
    CT-4      0.14 sigma    10.9 sigma    ~56%   <-- NUMERATOR is the problem
    CT-10     10.1 sigma     1.9 sigma    ~97%
    CT-8.2    10.2 sigma     1.3 sigma    ~91%
    CT-11     10.4 sigma     0.6 sigma    ~74%

For CT-4 it is the NUMERATOR that nearly vanishes, because E_CT-4 ~= E_CT-6 AND
P_CT-4 ~= P_CT-6: the child is nearly indistinguishable from parent B in both
erosion rate and production rate, so the two-endmember system is very nearly
degenerate and the ^10Be data cannot constrain even the SIGN of the quartz
fertility contrast. That is precisely why CT-4 is excluded from the reported
quartz results (manuscript Table 2, "CT-4 excluded"). The scenario is still
COMPUTED so the failure is visible rather than hidden, but the retention fraction
must always be reported alongside it.

For CT-10/CT-11/CT-8.2 it is the DENOMINATOR that comes closest to zero, which is
why those distributions have long right tails. Their retention is high enough for
the results to stand, but it is not 100% and must be reported.
"""

import numpy as np
import pandas as pd

# Mixing coefficients outside (THRESHOLD, 1 - THRESHOLD) are discarded: a parent
# contributing <1% of the child's zircon is not physically meaningful and drives
# the w_a/w_b ratio to extreme values. Note that because w_a + w_b = 1 exactly,
# the four conditions (a > T, a < 1-T, b > T, b < 1-T) collapse to the single
# condition T < a < 1-T. Both are checked anyway for explicitness.
THRESHOLD = 0.01


# ---------------------------------------------------------------------------
# CT-3* source area convention
# ---------------------------------------------------------------------------

def ct3_star_area(A, parent_a, child):
    """Return the CT-3* source area (km^2) for a given parent/child pair.

    CT-3* is NOT the CT-3 watershed. CT-3 is a small (189.2 km^2) Coastal
    Cordillera tributary that is used as a LITHOLOGIC PROXY for the entire
    portion of the child basin that is not covered by the Principal Cordillera
    parent. The area that proxy stands for is the residual:

        A_CT3* = A(child) - A(parent_a)

    This is the convention already used in the manuscript for the
    CT-8.2 & CT-3 -> CT-9 scenario, and it is applied here uniformly to every
    scenario in which CT-3 appears as a parent.

    Applying it consistently gives:

        parent_a   child    A_CT3*                        km^2
        --------   -----    --------------------------    --------
        CT-8.2     CT-9     A(CT-9)  - A(CT-8)            1239.321
        CT-8.2     CT-2     A(CT-2)  - A(CT-8)             809.733
        CT-4       CT-9     A(CT-9)  - A(CT-4)            1482.999
        CT-10      CT-9     A(CT-9)  - A(CT-10)           1324.774
        CT-11      CT-9     A(CT-9)  - A(CT-11)           1276.382

    Two corrections relative to the legacy notebook:

    (a) The CT-8.2 -> CT-2 case previously used
            A(CT-2) - A(CT-8) + A(CT-3) = 998.938 km^2
        The trailing `+ A(CT-3)` double-counted: it made the tiled area
        A(CT-8) + A_CT3* = A(CT-2) + A(CT-3) = 5949.3 km^2, which OVER-tiles the
        actual CT-2 basin (5760.1 km^2) by exactly A(CT-3). CT-3's drainage lies
        outside the CT-2 basin polygon, so it cannot be added to it. Under the
        proxy interpretation the correct area is simply the residual.

    (b) The CT-4/CT-10/CT-11 -> CT-9 scenarios previously used the RAW CT-3
        watershed area (189.2 km^2), a factor of ~7 too small, and inconsistent
        with the CT-8.2 case sitting in the same results table.

    Parameters
    ----------
    A : dict
        {Sample_ID: area_km2}, from sample_data.get_areas().
    parent_a : str
        Sample_ID of the Principal Cordillera parent ('CT-8', 'CT-4', 'CT-10',
        'CT-11'). Note CT-8.2 is stored as 'CT-8' in the sample spreadsheet.
    child : str
        Sample_ID of the child ('CT-9' or 'CT-2').

    Returns
    -------
    float
        CT-3* source area in km^2.
    """
    area = A[child] - A[parent_a]
    if area <= 0:
        raise ValueError(
            f"CT-3* area for parent {parent_a} -> child {child} is non-positive "
            f"({area:.3f} km^2). The parent basin is not contained within the "
            f"child basin; the residual convention does not apply."
        )
    return area


# ---------------------------------------------------------------------------
# Mixing model results loader
# ---------------------------------------------------------------------------

def load_mixing_results(xlsx_path, child, threshold=THRESHOLD):
    """Load one child's mixing results from a detritalPy-mix results workbook.

    Reads three sheets from a single `dPy-mix_*_results_{metric}.xlsx`:

        '<child>'    -- bootstrap ensemble. Columns are the two parent names
                        followed by the objective metric ('vmax' or 'r2-kde').
                        One row per bootstrap iteration (10,000 rows).
        'Best-fit'   -- deterministic best-fit mixing coefficients, one row per
                        child. These come from a grid search on the UNPERTURBED
                        data and so are seed-independent.
        'Model_fit'  -- one row per child, up to four diagnostic columns:
                          col 1: best-fit objective value
                          col 2: crit -- 95th percentile of comparisons between
                                 the child and BOOTSTRAPPED VERSIONS OF ITSELF.
                                 The noise floor: the best fit should beat it.
                          col 3: fraction of BOOTSTRAPPED BEST-FIT MIXTURES worse
                                 than the best-fit value. This is "worse than
                                 BEST-FIT", NOT "worse than crit" -- the legacy
                                 notebook mislabelled it.
                          col 4: fraction of RESAMPLED OBSERVATIONS that yielded a
                                 better objective value than the mixture model.
                                 detritalPy-mix's own note: a sample that is
                                 INDISTINGUISHABLE from the mixture should sit at
                                 ~0.5. Near 1.0 means almost any resample of the
                                 child fits it better than the mixture does -- the
                                 mixture is a poor description. This is the most
                                 directly interpretable fit diagnostic in the
                                 workbook, and the legacy notebook ignored it.
                                 (Some older workbooks omit it; NaN if absent.)

    Parameters
    ----------
    xlsx_path : str or Path
        Path to the results workbook.
    child : str
        Child sample name, which is also the sheet name (e.g. 'CT-8.2').
    threshold : float
        Mixing coefficients outside (threshold, 1 - threshold) are discarded.

    Returns
    -------
    dict with keys:
        'parents'          : [name_a, name_b]
        'metric'           : 'vmax' or 'r2-kde'
        'wa', 'wb'         : filtered bootstrap coefficient arrays (paired)
        'n_total'          : bootstrap iterations before filtering
        'n_kept'           : bootstrap iterations after filtering
        'best_a', 'best_b' : deterministic best-fit coefficients (full precision)
        'best_fit'         : best-fit objective value
        'crit'             : critical objective value (noise floor)
        'pct_boot_worse_than_best_fit'      : Model_fit col 3
        'pct_resampled_better_than_mixture' : Model_fit col 4 (~0.5 = ideal;
                                              NaN if the workbook omits it)
        'gof_pass'         : bool. For vmax (lower is better): best_fit < crit.
                             For r2-kde (higher is better): best_fit > crit.
    """
    # --- bootstrap ensemble -------------------------------------------------
    boot = pd.read_excel(xlsx_path, sheet_name=child)
    # Columns: [parent_a, parent_b, metric_name]. The metric name is the third
    # column header, so we read the metric off the file rather than trusting the
    # filename -- this makes a mislabelled file impossible to use silently.
    parent_a, parent_b = boot.columns[0], boot.columns[1]
    metric = str(boot.columns[2])

    a = boot[parent_a].to_numpy(dtype=float)
    b = boot[parent_b].to_numpy(dtype=float)
    n_total = len(a)

    keep = (a > threshold) & (a < 1 - threshold) & (b > threshold) & (b < 1 - threshold)
    wa, wb = a[keep], b[keep]

    # --- deterministic best fit --------------------------------------------
    bf = pd.read_excel(xlsx_path, sheet_name='Best-fit').set_index('Sample')
    best_a = float(bf.loc[child, parent_a])
    best_b = float(bf.loc[child, parent_b])

    # --- model fit / goodness of fit ---------------------------------------
    # Column headers in this sheet embed the metric name and are verbose, so
    # index positionally:
    #   [Sample, best-fit, crit, %boot-worse-than-best-fit, %resampled-better]
    mf = pd.read_excel(xlsx_path, sheet_name='Model_fit')
    mf = mf[mf.iloc[:, 0] == child]
    best_fit = float(mf.iloc[0, 1])
    crit = float(mf.iloc[0, 2])
    pct_boot_worse = float(mf.iloc[0, 3])
    # Column 4 is the resampled-observations diagnostic (~0.5 = indistinguishable
    # from the mixture). Tolerate its absence in older workbooks.
    if mf.shape[1] > 4:
        pct_resampled_better = float(mf.iloc[0, 4])
    else:
        pct_resampled_better = float('nan')

    # vmax / Kuiper: LOWER is better, so a good fit sits BELOW the noise floor.
    # r2-kde:        HIGHER is better, so a good fit sits ABOVE the noise floor.
    if metric.lower().startswith('vmax'):
        gof_pass = best_fit < crit
    else:
        gof_pass = best_fit > crit

    return {
        'parents': [parent_a, parent_b],
        'metric': metric,
        'wa': wa,
        'wb': wb,
        'n_total': n_total,
        'n_kept': int(keep.sum()),
        'best_a': best_a,
        'best_b': best_b,
        'best_fit': best_fit,
        'crit': crit,
        'pct_boot_worse_than_best_fit': pct_boot_worse,
        'pct_resampled_better_than_mixture': pct_resampled_better,
        'gof_pass': bool(gof_pass),
    }


# ---------------------------------------------------------------------------
# Quartz
# ---------------------------------------------------------------------------

def quartz_fertility(E_a, E_b, E_c, A_a, A_b, P_a, P_b, P_c):
    """Quartz fertility ratio Q_A/Q_B and quartz flux ratio F_QA/F_QB.

    Both are computed elementwise over paired Monte Carlo erosion-rate draws, so
    E_a, E_b, E_c must be arrays of equal length whose index i corresponds to the
    same Monte Carlo iteration.

    Fertility (manuscript Eq. 5):

        Q_A/Q_B = [A_b * (P_b*E_c - P_c*E_b)] / [A_a * (P_c*E_a - P_a*E_c)]

    Flux ratio (areas cancel; see module docstring):

        F_QA/F_QB = [E_a * (P_b*E_c - P_c*E_b)] / [E_b * (P_c*E_a - P_a*E_c)]

    Both share the two kernels (P_b*E_c - P_c*E_b) and (P_c*E_a - P_a*E_c), EITHER
    of which can straddle zero within the Monte Carlo cloud. Draws for which the
    resulting ratio is negative are unphysical and are masked out; `mask` records
    which draws survived so the retention fraction can be reported. See the module
    docstring for which kernel is the binding one in each scenario -- it is the
    NUMERATOR for CT-4 and the DENOMINATOR for CT-10/CT-11/CT-8.2.

    NOTE: A_c does not appear. The child area cancels out of the mass balance.

    Parameters
    ----------
    E_a, E_b, E_c : np.ndarray
        Paired erosion rate draws (m/Myr) for parent A, parent B, child C.
    A_a, A_b : float
        Source areas (km^2) for parent A and parent B.
    P_a, P_b, P_c : float
        Surface production rates (atoms/g/yr).

    Returns
    -------
    Qa_Qb : np.ndarray
        Quartz fertility ratio, physical draws only.
    Fa_Fb : np.ndarray
        Quartz flux ratio, same draws (identical mask).
    mask : np.ndarray of bool
        True where the draw was retained.
    """
    E_a = np.asarray(E_a, dtype=float)
    E_b = np.asarray(E_b, dtype=float)
    E_c = np.asarray(E_c, dtype=float)
    if not (len(E_a) == len(E_b) == len(E_c)):
        raise ValueError(
            f"E_a, E_b, E_c must be paired and equal length; got "
            f"{len(E_a)}, {len(E_b)}, {len(E_c)}"
        )

    # Shared kernels. Splitting them out makes the two possible sign flips
    # explicit: EITHER can straddle zero, and which one does differs by scenario.
    num_k = P_b * E_c - P_c * E_b   # vanishes when E_c/E_b -> P_c/P_b (CT-4's problem)
    den_k = P_c * E_a - P_a * E_c   # vanishes when E_a/E_c -> P_a/P_c (CT-10/11/8.2)

    with np.errstate(divide='ignore', invalid='ignore'):
        Qa_Qb = (A_b * num_k) / (A_a * den_k)
        Fa_Fb = (E_a * num_k) / (E_b * den_k)

    # A negative fertility is unphysical. Because F = Q * (E_a A_a)/(E_b A_b) and
    # erosion rates and areas are positive, Q and F always share a sign -- one
    # mask serves both.
    mask = np.isfinite(Qa_Qb) & (Qa_Qb > 0)

    return Qa_Qb[mask], Fa_Fb[mask], mask


# ---------------------------------------------------------------------------
# Zircon
# ---------------------------------------------------------------------------

def zircon_fertility(wa, wb, ER_A, ER_B, A_A, A_B, num_samples=10000, seed=17):
    """Zircon fertility ratio Z_A/Z_B and zircon flux ratio F_ZA/F_ZB.

    Fertility (manuscript Eq. 11):

        Z_A/Z_B = (w_a / w_b) * (E_B * A_B) / (E_A * A_A)

    Flux ratio:

        F_ZA/F_ZB = w_a / w_b

    The zircon flux ratio requires no erosion rates and no areas -- the mixing
    coefficient IS the fraction of the child's zircon supplied by each parent, so
    their ratio IS the zircon flux ratio, by construction. It is therefore
    completely independent of the ^10Be data and of the source-area assumptions.
    It is reported explicitly because it, not Z_A/Z_B, is the quantity that
    answers "how much zircon comes from each tributary".

    Two independent random draws are combined per Monte Carlo sample:
      - one bootstrap iteration of the mixing coefficients (from detritalPy-mix)
      - one Monte Carlo draw of each parent's erosion rate
    These are independent sources of uncertainty, so pairing them at random is
    correct. Both are drawn with replacement using a LOCAL Generator, so this
    function no longer perturbs the global numpy random state.

    Parameters
    ----------
    wa, wb : np.ndarray
        Paired, filtered bootstrap mixing coefficient arrays.
    ER_A, ER_B : np.ndarray
        Erosion rate draws (m/Myr) for parents A and B. Need not be the same
        length as wa/wb, and need not be the same length as each other.
    A_A, A_B : float
        Source areas (km^2).
    num_samples : int
        Number of Monte Carlo draws to generate.
    seed : int
        Seed for the local Generator.

    Returns
    -------
    dict
        {'ZA_ZB': stats, 'ZB_ZA': stats, 'FZA_FZB': stats, 'FZB_FZA': stats}
        where each `stats` dict carries samples, median, p25, p75 and their
        log10 equivalents.
    """
    rng = np.random.default_rng(seed)

    wa = np.asarray(wa, dtype=float)
    wb = np.asarray(wb, dtype=float)
    ER_A = np.asarray(ER_A, dtype=float)
    ER_B = np.asarray(ER_B, dtype=float)

    # Independent resampling of each uncertainty source, with replacement.
    i_w = rng.integers(0, len(wa), size=num_samples)
    i_a = rng.integers(0, len(ER_A), size=num_samples)
    i_b = rng.integers(0, len(ER_B), size=num_samples)

    w_ratio = wa[i_w] / wb[i_w]            # this IS the zircon flux ratio
    ER_A_s = ER_A[i_a]
    ER_B_s = ER_B[i_b]

    ZA_ZB = w_ratio * (ER_B_s * A_B) / (ER_A_s * A_A)

    def stats(x):
        lx = np.log10(x)
        return {
            'samples': x,
            'median': float(np.median(x)),
            'p25': float(np.percentile(x, 25)),
            'p75': float(np.percentile(x, 75)),
            'log10_samples': lx,
            'log10_median': float(np.median(lx)),
            'log10_p25': float(np.percentile(lx, 25)),
            'log10_p75': float(np.percentile(lx, 75)),
        }

    return {
        'ZA_ZB': stats(ZA_ZB),
        'ZB_ZA': stats(1.0 / ZA_ZB),
        'FZA_FZB': stats(w_ratio),
        'FZB_FZA': stats(1.0 / w_ratio),
    }


# ---------------------------------------------------------------------------
# Shared summary helper
# ---------------------------------------------------------------------------

def ratio_stats(x):
    """Median / IQR of a ratio, in linear and log10 space.

    A note on the back-transform, because the legacy notebook's comments were
    misleading on this point: log10 is strictly increasing, so it preserves rank
    order, so it preserves RANK statistics exactly:

        10 ** median(log10(x))  ==  median(x)          (exactly)
        10 ** percentile(log10(x), q)  ==  percentile(x, q)   (exactly)

    There is therefore no such thing as "back-transforming the median from log
    space" -- it is a no-op, and the linear median reported here is identical to
    the linear median of the raw samples. (This would NOT hold for the mean.)

    log10 is used purely as a DISPLAY transform: it symmetrises a ratio about 1,
    so that "A is 5x B" and "B is 5x A" are equidistant from zero rather than
    being 5 and 0.2.
    """
    lx = np.log10(x)
    return {
        'n': int(len(x)),
        'median': float(np.median(x)),
        'p25': float(np.percentile(x, 25)),
        'p75': float(np.percentile(x, 75)),
        'log10_median': float(np.median(lx)),
        'log10_p25': float(np.percentile(lx, 25)),
        'log10_p75': float(np.percentile(lx, 75)),
    }


def fmt_lin(st):
    """Format a ratio_stats dict as 'median (p25 - p75)' in LINEAR space."""
    return f"{st['median']:.4f} ({st['p25']:.4f} - {st['p75']:.4f})"


def fmt_log(st):
    """Format a ratio_stats dict as 'median (p25 - p75)' in LOG10 space."""
    return f"{st['log10_median']:.4f} ({st['log10_p25']:.4f} - {st['log10_p75']:.4f})"
