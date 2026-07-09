"""
run_mixing_scenarios.py

Batch runner for detritalPy-mix top-down mixture modeling of detrital zircon
U-Pb ages, across every parent/child scenario for the Mataquito dataset.

This replaces manually editing parent_list / child_list in the notebook and
re-running cell-by-cell. Each scenario is run once with objective_metric='vmax'
and once with 'r2-kde', producing the same figures/CSV/xlsx outputs the
notebook produced, saved to MIXING_FIG_DIR / MIXING_RES_DIR.

Expected location: MataquitoWatershedAnalysis/Mataquito/mataquito/run_mixing_scenarios.py
"""

import os
import sys
import io
import time
import argparse
import traceback
import contextlib
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # no interactive windows -- figures are only saved to disk
from matplotlib import pyplot as plt

import detritalpy
from detritalpy import detritalFuncs as dFunc
from detritalpy import detritalMixer as dMix


# ─────────────────────────────────────────────────────────────────────────────
# Directory structure (same layout as the notebook, resolved from this file's
# location instead of Path.cwd(), since a script's cwd isn't guaranteed to be
# where the file lives)
# ─────────────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent   # .../Mataquito
DATA_DIR = os.path.join(BASE_DIR, 'data')
FIG_DIR  = os.path.join(BASE_DIR, 'figures')
RES_DIR  = os.path.join(BASE_DIR, 'results')

ZIRCON_DATA_DIR = os.path.join(DATA_DIR, 'geochronology', 'zircon_upb')
MIXING_FIG_DIR  = os.path.join(FIG_DIR, 'zircon', 'mixing')
MIXING_RES_DIR  = os.path.join(RES_DIR, 'mixing_coefficients')

for folder in [MIXING_FIG_DIR, MIXING_RES_DIR]:
    os.makedirs(folder, exist_ok=True)

LOG_PATH = os.path.join(MIXING_RES_DIR, 'run_log.txt')


# ─────────────────────────────────────────────────────────────────────────────
# Scenario definitions
# Each entry: name (used for filenames/log tags), parent_list, child_list.
# This is the full set of parent/child combinations from the notebook's
# commented-out options, excluding the three-parent mixtures (those were
# tests, not part of the real sweep).
# ─────────────────────────────────────────────────────────────────────────────
SCENARIOS = [
    # --- Parents CT-7, CT-1: each child run individually ---
    {'name': 'ct7ct1_to_CT4',   'parent_list': ['CT-7', 'CT-1'], 'child_list': ['CT-4']},
    {'name': 'ct7ct1_to_CT10',  'parent_list': ['CT-7', 'CT-1'], 'child_list': ['CT-10']},
    {'name': 'ct7ct1_to_CT11',  'parent_list': ['CT-7', 'CT-1'], 'child_list': ['CT-11']},
    {'name': 'ct7ct1_to_CT8_2', 'parent_list': ['CT-7', 'CT-1'], 'child_list': ['CT-8.2']},
    {'name': 'ct7ct1_to_CT2',   'parent_list': ['CT-7', 'CT-1'], 'child_list': ['CT-2']},
    {'name': 'ct7ct1_to_CT9',   'parent_list': ['CT-7', 'CT-1'], 'child_list': ['CT-9']},

    # --- Parents CT-5, CT-6: batched Principal Cordillera children, plus each
    #     individually, plus CT-2 and CT-9 ---
    {'name': 'ct5ct6_batch_CT4_CT10_CT11_CT8_2', 'parent_list': ['CT-5', 'CT-6'],
     'child_list': ['CT-4', 'CT-10', 'CT-11', 'CT-8.2']},
    {'name': 'ct5ct6_to_CT4',   'parent_list': ['CT-5', 'CT-6'], 'child_list': ['CT-4']},
    {'name': 'ct5ct6_to_CT10',  'parent_list': ['CT-5', 'CT-6'], 'child_list': ['CT-10']},
    {'name': 'ct5ct6_to_CT11',  'parent_list': ['CT-5', 'CT-6'], 'child_list': ['CT-11']},
    {'name': 'ct5ct6_to_CT8_2', 'parent_list': ['CT-5', 'CT-6'], 'child_list': ['CT-8.2']},
    {'name': 'ct5ct6_to_CT2',   'parent_list': ['CT-5', 'CT-6'], 'child_list': ['CT-2']},
    {'name': 'ct5ct6_to_CT9',   'parent_list': ['CT-5', 'CT-6'], 'child_list': ['CT-9']},

    # --- Child CT-2, two-parent mixtures with CT-3 ---
    {'name': 'CT2_ct4ct3',   'parent_list': ['CT-4', 'CT-3'],   'child_list': ['CT-2']},
    {'name': 'CT2_ct10ct3',  'parent_list': ['CT-10', 'CT-3'],  'child_list': ['CT-2']},
    {'name': 'CT2_ct11ct3',  'parent_list': ['CT-11', 'CT-3'],  'child_list': ['CT-2']},
    {'name': 'CT2_ct8_2ct3', 'parent_list': ['CT-8.2', 'CT-3'], 'child_list': ['CT-2']},

    # --- Child CT-9, two-parent mixtures with CT-3 ---
    {'name': 'CT9_ct4ct3',   'parent_list': ['CT-4', 'CT-3'],   'child_list': ['CT-9']},
    {'name': 'CT9_ct10ct3',  'parent_list': ['CT-10', 'CT-3'],  'child_list': ['CT-9']},
    {'name': 'CT9_ct11ct3',  'parent_list': ['CT-11', 'CT-3'],  'child_list': ['CT-9']},
    {'name': 'CT9_ct8_2ct3', 'parent_list': ['CT-8.2', 'CT-3'], 'child_list': ['CT-9']},
    {'name': 'CT9_ct2ct3',   'parent_list': ['CT-2', 'CT-3'],   'child_list': ['CT-9']},
]

# Trim this list to run a subset without touching SCENARIOS above.
ACTIVE_SCENARIOS = [s['name'] for s in SCENARIOS]

OBJECTIVE_METRICS = ['vmax', 'r2-kde']


# ─────────────────────────────────────────────────────────────────────────────
# Fixed run parameters (identical across every scenario)
# ─────────────────────────────────────────────────────────────────────────────
X1, X2, XDIF = 0, 4500, 1          # age range / interval (Myr) for distribution calcs
BW, BW_X = [1.5, 5], [25]          # KDE bandwidth for mixing calc; only matters for 'r2-kde'
SIGMA = '2sigma'                   # input data uncertainty level

N_BOOTSTRAP_ITERATIONS_DEFAULT = 10000  # overridable with --iterations
DO_PERTURB_RESAMPLED_AGES = True
N_GRAINS_TO_RESAMPLE = None

# Samples with older/more dispersed age spectra (Coastal Cordillera / outlet)
# need a wider plotting window than the Principal Cordillera-only scenarios.
XAXIS2_TRIGGER_SAMPLES = {'CT-3', 'CT-2', 'CT-9'}


def get_xaxis2(parent_list, child_list):
    """Plotting x-axis max (Ma) for fig1/fig3: 400 if CT-3/CT-2/CT-9 appear
    anywhere in the scenario, 25 otherwise."""
    samples_in_scenario = set(parent_list) | set(child_list)
    if samples_in_scenario & XAXIS2_TRIGGER_SAMPLES:
        return 400
    return 25


# ─────────────────────────────────────────────────────────────────────────────
# Logging: print to console and append to a persistent logfile, so progress
# can be checked even if the terminal isn't being watched during a multi-hour run.
# ─────────────────────────────────────────────────────────────────────────────
def log(msg):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = f'[{timestamp}] {msg}'
    print(line)
    with open(LOG_PATH, 'a') as f:
        f.write(line + '\n')


def write_results_block(results_file, header, lines):
    """Append a results section (header + indented lines) to a scenario's
    summary text file. Also echoed to console so you can watch it live."""
    print(header)
    with open(results_file, 'a') as f:
        f.write(header + '\n')
        for content_line in lines:
            print(f'    {content_line}')
            f.write(f'    {content_line}\n')
        f.write('\n')


def elapsed_min(start_time):
    return (time.time() - start_time) / 60.0


# ─────────────────────────────────────────────────────────────────────────────
# Confirmation: show the full run plan before touching any data
# ─────────────────────────────────────────────────────────────────────────────
def confirm_run_plan(scenarios, active_names, metrics, n_iterations):
    active = [s for s in scenarios if s['name'] in active_names]
    total_runs = len(active) * len(metrics)

    print(f"\n{total_runs} runs queued ({len(active)} scenarios x {len(metrics)} objective metrics: {', '.join(metrics)})")
    print(f"Bootstrap iterations per run: {n_iterations}\n")
    for s in active:
        xaxis2 = get_xaxis2(s['parent_list'], s['child_list'])
        parents_str = ','.join(s['parent_list'])
        children_str = ','.join(s['child_list'])
        print(f"  {s['name']:<28} parents=[{parents_str}]  children=[{children_str}]  xaxis2={xaxis2}")

    answer = input("\nProceed? [y/n]: ").strip().lower()
    return answer.startswith('y')


# ─────────────────────────────────────────────────────────────────────────────
# Core scenario runner -- one full pass (best fit -> bootstrap -> figures ->
# exports) for one (scenario, objective_metric) combination.
# ─────────────────────────────────────────────────────────────────────────────
def run_scenario_metric(scenario, objective_metric, main_byid_df, n_bootstrap_iterations):
    name = scenario['name']
    parent_list = scenario['parent_list']
    child_list = scenario['child_list']
    tag = f"[{name}] {objective_metric}"

    xaxis_1 = 0
    xaxis_2 = get_xaxis2(parent_list, child_list)

    parent_str = "_".join(p.lower() for p in parent_list)
    child_str = "_".join(c.lower() for c in child_list)
    metric_tag = objective_metric.replace('-', '')

    # Separate text file for this scenario+metric's printed results (best-fit
    # coefficients, bootstrap summary, model fit) -- kept apart from run_log.txt,
    # which is just step-by-step progress/timing.
    results_file = os.path.join(MIXING_RES_DIR, f'dPy-mix_{parent_str}_to_{child_str}_summary_{metric_tag}.txt')
    open(results_file, 'w').close()  # start fresh each run

    plt.close('all')  # avoid figure accumulation across runs
    t_scenario = time.time()

    # --- Best-fit mixture ---
    log(f"{tag}: starting best-fit optimization")
    t0 = time.time()
    mix_coeffs_bf, obj_func_val, best_mixed_dist = dMix.find_best_fit_mix(
        parent_list, child_list, main_byid_df, sigma=SIGMA, objective_metric=objective_metric,
        x1=X1, x2=X2, xdif=XDIF, bw=BW, bw_x=BW_X, verbose=False)
    log(f"{tag}: best-fit done ({elapsed_min(t0):.1f} min)")

    # Best-fit mixing coefficients per child -- same content dMix.print_best_fit_mixture
    # shows in the notebook, captured here so it lands in the results file.
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        dMix.print_best_fit_mixture(parent_list, child_list, objective_metric, obj_func_val, mix_coeffs_bf)
    write_results_block(results_file, "Best-fit results", captured.getvalue().strip().split('\n'))

    # --- Bootstrap mixture solutions ---
    log(f"{tag}: starting bootstrap ({n_bootstrap_iterations} iterations)")
    t0 = time.time()
    mix_coeffs_all, obj_vals_all, child_modelled_distributions = dMix.bootstrap_solve_mixture_models(
        parent_list, child_list, main_byid_df,
        sigma=SIGMA, objective_metric=objective_metric,
        nBootstrapIterations=n_bootstrap_iterations,
        doPerturbResampledAges=DO_PERTURB_RESAMPLED_AGES,
        nGrainsToResample=N_GRAINS_TO_RESAMPLE,
        x1=X1, x2=X2, xdif=XDIF, bw=BW, bw_x=BW_X,
        verbose=False, update_freq=100)
    log(f"{tag}: bootstrap done ({elapsed_min(t0):.1f} min)")

    # --- Self-comparison (null distribution) ---
    t0 = time.time()
    selfCompMetrics_bs_set, childDists_bs_set = dMix.bootstrapped_self_comparisons_many_samples(
        main_byid_df, child_list,
        doPerturbResampledAges=DO_PERTURB_RESAMPLED_AGES,
        objective_metric=objective_metric, sigma=SIGMA,
        x1=X1, x2=X2, xdif=XDIF, bw=BW, bw_x=BW_X,
        nBootstrapIterations=n_bootstrap_iterations)
    obj_func_crit, worse_than_crit = dMix.calc_model_fit(
        child_list, obj_func_val, obj_vals_all, selfCompMetrics_bs_set, objective_metric=objective_metric)
    log(f"{tag}: self-comparison done ({elapsed_min(t0):.1f} min)")

    # Bootstrapped mixing coefficient summary: median and 95% CI per parent per child
    summary_lines = []
    for i, child in enumerate(child_list):
        summary_lines.append(f"{child}:")
        coeffs = np.array(mix_coeffs_all[i])
        for j, parent in enumerate(parent_list):
            lo, med, hi = np.percentile(coeffs[:, j], [2.5, 50, 97.5])
            summary_lines.append(f"  {parent}: median={med:.2f}, 95% CI [{lo:.2f}, {hi:.2f}], width={hi-lo:.2f}")
    write_results_block(results_file, "Bootstrapped mixing coefficient summary", summary_lines)

    fit_lines = []
    for i, child in enumerate(child_list):
        fit_lines.append(
            f"{child}: best-fit {objective_metric}={obj_func_val[i]:.3f}, "
            f"crit={obj_func_crit[i]:.3f}, worse_than_crit={worse_than_crit[i]:.3f}")
    write_results_block(results_file, "Model fit summary", fit_lines)

    # --- Export bootstrapped coefficients to CSV (one file per child) ---
    for i, child in enumerate(child_list):
        coeffs = np.array(mix_coeffs_all[i])
        df = pd.DataFrame(coeffs, columns=parent_list)
        filename = f"mix_coeffs_all_{parent_str}_to_{child.lower()}_{metric_tag}.csv"
        df.to_csv(os.path.join(MIXING_RES_DIR, filename), index=False)

    # --- fig: histograms of bootstrapped mixing coefficients ---
    colors = ['green', 'firebrick']
    fig_hist, axs_hist = plt.subplots(len(child_list), 1, figsize=(6, 3 * len(child_list)))
    if len(child_list) == 1:
        axs_hist = [axs_hist]
    for i, child in enumerate(child_list):
        coeffs = np.array(mix_coeffs_all[i])
        ax = axs_hist[i]
        for j, parent in enumerate(parent_list):
            color = colors[j % len(colors)]
            ax.hist(coeffs[:, j], bins=50, alpha=0.8, color=color, edgecolor='black',
                    linewidth=0.5, label=parent)
        ax.set_xlabel('Mixing coefficient value')
        ax.set_ylabel('Frequency')
        ax.set_title(f'Histogram of Mix Coefficients: {child}')
        ax.legend()
        ax.grid(True)
    fig_hist.tight_layout()
    fig_hist.savefig(os.path.join(MIXING_FIG_DIR, f'dPy-mix_{parent_str}_to_{child_str}_coefficient_histograms_{metric_tag}.pdf'))

    # --- fig1: age distribution + mixture (plotMix) ---
    fig1 = dMix.plotMix(
        main_byid_df, parent_list, child_list, plotType='KDE', bw=BW, bw_x=BW_X, x1=X1, x2=X2, xdif=XDIF,
        fillParent=True, parent_colors='Default', child_colors='Default',
        color_by_age=False, agebins=[0, 20, 100, 200], agebinsc=['red', 'orange', 'green'], agebinsc_alpha=0.5,
        xaxis_1=xaxis_1, xaxis_2=xaxis_2, w1=8, w2=4, c=4, plotPie=True, plotMixResults=True,
        best_plot_type='best-fit', plotResultType='violin', violin_width=0.2, sigma=SIGMA,
        mix_coeffs_all=mix_coeffs_all, mix_coeffs_bf=mix_coeffs_bf, obj_func_val=obj_func_val,
        best_mixed_dist=best_mixed_dist, obj_vals_all=obj_vals_all, objective_metric=objective_metric)
    fig1.savefig(os.path.join(MIXING_FIG_DIR, f'dPy-mix_{parent_str}_to_{child_str}_age_distribution_{metric_tag}.pdf'))

    # --- fig2: objective function comparisons (model vs. self-comparison null) ---
    fig2, axs2 = plt.subplots(len(selfCompMetrics_bs_set), 1, figsize=(5.0, 3.0 * len(child_list)), sharex=True)
    dMix.plot_many_bootstrapped_metric_comparisons_model_observations(
        obj_vals_all, selfCompMetrics_bs_set, objective_metric, child_list, main_byid_df,
        obj_func_crit, worse_than_crit, obj_func_val,
        modelled_colors=None, self_compared_colors=None,
        axs=axs2, plotWidth=5.0, subplotHeight=3.0, doAddSummaryTitle=True, alpha=0.5, bins=50)
    fig2.savefig(os.path.join(MIXING_FIG_DIR, f'dPy-mix_{parent_str}_to_{child_str}_objective_function_distributions_{metric_tag}.pdf'))

    # --- fig3: modelled vs. observed distribution comparison ---
    fig3, axs3 = plt.subplots(len(child_list), 1, figsize=(6.0, 3.0 * len(child_list)), sharex=True, sharey='col')
    dMix.plot_child_bootstrappedmodel_distribution_comparison(
        main_byid_df, child_modelled_distributions, child_list,
        xaxis_1=xaxis_1, xaxis_2=xaxis_2, x1=X1, x2=X2, xdif=XDIF, bw=BW, bw_x=BW_X,
        objective_metric=objective_metric, confidence_interval=95.0, fill_alpha=0.25,
        plot_self_comparisons=True, childDists_bs_set=childDists_bs_set,
        axs=axs3, plotWidth=6.0, subplotHeight=3.0, child_colors='Default',
        model_color='black', sigma=SIGMA)
    fig3.savefig(os.path.join(MIXING_FIG_DIR, f'dPy-mix_{parent_str}_to_{child_str}_distribution_comparisons_{metric_tag}.pdf'))

    # --- fig4: mixing coefficients vs. stratigraphy -- skipped for single-child scenarios ---
    if len(child_list) > 1:
        fig4, ax4 = dMix.plot_bootstrapped_mixturecoefficients_stratigraphy(
            parent_list, child_list, mix_coeffs_bf, mix_coeffs_all,
            confidence_interval=95, ax=None, plotWidth=15.0, plotHeight=4.0,
            parent_colors='Default', do_plot_errorbars=True,
            yAxisValues=None, doFlipXY=False, plot_alpha=0.3,
            separate_parents=True, best_plot_type='best-fit')
        fig4.savefig(os.path.join(MIXING_FIG_DIR, f'dPy-mix_{parent_str}_to_{child_str}_spatial_comparison_{metric_tag}.pdf'))

    # --- Export full results workbook ---
    file_name = os.path.join(MIXING_RES_DIR, f'dPy-mix_{parent_str}_to_{child_str}_results_{metric_tag}.xlsx')
    dMix.export_results(
        parent_list, child_list, main_byid_df, objective_metric, XDIF, BW, mix_coeffs_bf,
        obj_func_val, file_name=file_name, verbose=False, version=detritalpy.__version__,
        bootstrap=True, nBootstrapIterations=n_bootstrap_iterations,
        doPerturbResampledAges=DO_PERTURB_RESAMPLED_AGES, nGrainsToResample=N_GRAINS_TO_RESAMPLE,
        mix_coeffs_all=mix_coeffs_all, obj_vals_all=obj_vals_all, obj_func_crit=obj_func_crit,
        worse_than_crit=worse_than_crit, selfCompMetrics_bs_set=selfCompMetrics_bs_set)

    log(f"{tag}: figures + exports saved")
    log(f"{tag}: SCENARIO COMPLETE ({elapsed_min(t_scenario):.1f} min)")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    valid_names = {s['name'] for s in SCENARIOS}

    parser = argparse.ArgumentParser(description="Run detritalPy-mix scenarios.")
    parser.add_argument('scenarios', nargs='*',
                         help="Scenario names to run (default: everything in ACTIVE_SCENARIOS). "
                              "e.g. ct7ct1_to_CT4 ct5ct6_to_CT4")
    parser.add_argument('--iterations', '-n', type=int, default=N_BOOTSTRAP_ITERATIONS_DEFAULT,
                         help=f"Bootstrap iterations per run (default: {N_BOOTSTRAP_ITERATIONS_DEFAULT}). "
                              "Use a small number for a quick test.")
    args = parser.parse_args()

    if args.scenarios:
        unknown = [n for n in args.scenarios if n not in valid_names]
        if unknown:
            print(f"Unknown scenario name(s): {unknown}")
            print(f"Valid names: {sorted(valid_names)}")
            sys.exit(1)
        active_names = args.scenarios
    else:
        active_names = ACTIVE_SCENARIOS

    n_bootstrap_iterations = args.iterations

    if not confirm_run_plan(SCENARIOS, active_names, OBJECTIVE_METRICS, n_bootstrap_iterations):
        print("Aborted -- edit ACTIVE_SCENARIOS (or pass scenario names/--iterations as arguments) and rerun.")
        sys.exit(0)

    log("Loading dataset...")
    dataToLoad = [os.path.join(ZIRCON_DATA_DIR, 'CT-1-11_ZrUPb_datasheets.xlsx')]
    main_df, main_byid_df, samples_df, analyses_df = dFunc.loadDataExcel(dataToLoad, dataSheet='ZrUPb')
    log("Dataset loaded.")

    active = [s for s in SCENARIOS if s['name'] in active_names]
    total_runs = len(active) * len(OBJECTIVE_METRICS)
    run_count = 0

    for scenario in active:
        for metric in OBJECTIVE_METRICS:
            run_count += 1
            log(f"--- Run {run_count}/{total_runs}: {scenario['name']} / {metric} ---")
            try:
                run_scenario_metric(scenario, metric, main_byid_df, n_bootstrap_iterations)
            except Exception:
                log(f"[{scenario['name']}] {metric}: FAILED -- see traceback below")
                log(traceback.format_exc())
                continue

    log(f"All runs finished ({run_count}/{total_runs} attempted).")


if __name__ == '__main__':
    main()
