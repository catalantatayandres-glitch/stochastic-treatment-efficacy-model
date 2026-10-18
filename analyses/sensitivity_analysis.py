"""
sensitivity_analysis.py: One-at-a-time parameter exploration (fixed high-effort treatment)
================================================================================

Purpose
-------
Tests 3-5 alternative values for EACH ecological/disease parameter (I0, N0,
p_C0, r_I, r_N, K_I, K_N, m_I, m_N, alpha_IN, alpha_NI, beta, mu_A,
gamma_N, invasive_carrier_fitness_cost, prob_recruit_carrier,
migration_rate), one parameter at a time, while holding treatment fixed
at a strong "high-effort" setting (efficacy=100%, overall_encounter_rate=80%,
mix_fraction_kill=0.8) and all OTHER ecological parameters at their
baseline value.

This is a classic one-at-a-time (OAT) sensitivity exploration -- distinct
from a global (Sobol) analysis, which varies multiple parameters
SIMULTANEOUSLY to capture interactions. OAT is simpler to read and
communicate (one clean chart per parameter), at the cost of not capturing
interaction effects between parameters -- a Sobol analysis remains the more
rigorous tool for that. Use this script for a clear "what matters, and by
how much" team presentation; use a Sobol analysis if you need statistically
rigorous interaction-aware sensitivity indices.

Baseline parameters
--------------------
Updated to use the literature-grounded values established for this
project (see the project wiki Model page), with two axes
deliberately FIXED rather than tested at their own "baseline" guess:
  - Habitat is fixed at the "Favors Native" scenario (K_I=5000, K_N=6400
    -- see the habitat definitions), since that is the literature-grounded density
    range's lower-invasive/higher-native end.
  - Competition and transmission are fixed at their most MODERATE level
    (alpha_IN=0.3, alpha_NI=0.5, beta=0.5), matching the "Balanced"
    dominance scenario used in the habitat scenarios, since neither has a
    literature value of its own to anchor to.
N0 is derived from the fixed habitat as N0 = 0.2 x K_N = 1280 (same
0.2-fraction convention used throughout), used as N0's baseline
value here -- N0 itself remains an independently TESTED parameter below
(distinct from testing K_N itself), since "how sensitive is the outcome
to the starting population size alone" is a different question from "how
sensitive is it to the habitat's capacity."
I0=1000, p_C0=0.2, m_I=0.08, m_N=0.01, r_I=0.45, r_N=0.12, mu_A=0.8, etc.
Treatment is fixed at efficacy_I=efficacy_N_protect=1.0,
overall_encounter_rate=0.8, mix_fraction_kill=0.8 for every single run in
this script -- only the listed ecological/disease parameters vary.

Each parameter's tested value range is informed by the literature
discussion earlier in this project where available (e.g. p_C0 spanning
the documented 0-80% prevalence range, mu_A spanning the documented
resistant-to-severe range), and is otherwise a reasonable illustrative
spread around the baseline for parameters with no PRECISELY fitted
literature value (e.g. beta, alpha_IN, alpha_NI, gamma_N,
invasive_carrier_fitness_cost, migration_rate -- migration_rate does
have a qualitative literature basis, Bubb et al. 2004/2006, but not a
directly usable "% per month" statistic; see the model engine (engine.py) module docstring).

Update in this version
-----------------------
`plot_parameter_dashboard` now uses a SHARED y-axis range across all 17
panels (instead of each panel auto-scaling to its own data), so the
steepness/impact of different parameters can be compared directly by eye.
By default the shared range runs from 0 to the single highest native
population value found anywhere across all parameters and replicates; it
can also be fixed manually via the new `y_range` argument.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, Any, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. Import the model engine
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))  # lets the script run without `pip install -e .`

from stochastic_treatment_efficacy_model import engine  # noqa: E402

FIG_DIR = REPO_ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# 2. Baseline parameters (exactly as specified) + fixed high-effort treatment
# ---------------------------------------------------------------------------

BASE_PARAMS: Dict[str, Any] = {
    # --- Initial conditions ---
    "I0": 1000,   # = 0.2 * K_I, derived from the fixed "Favors Native" habitat (K_I=5000)
    "N0": 1280,   # = 0.2 * K_N, derived from the fixed "Favors Native" habitat (K_N=6400)
    "p_C0": 0.2,  # literature value (France survey national average) -- was previously 0.5 in error
    "N_I0": 0,

    # --- Demography ---
    "r_I": 0.45,
    "r_N": 0.12,
    "K_I": 5000.0,   # habitat FIXED at "Favors Native" scenario (see the habitat definitions)
    "K_N": 6400.0,   # habitat FIXED at "Favors Native" scenario (see the habitat definitions)
    "m_I": 0.08,
    "m_N": 0.01,

    # --- Competition: FIXED at the most moderate level (no literature value exists) ---
    "alpha_IN": 0.3,
    "alpha_NI": 0.5,

    # --- Disease ---
    "beta": 0.5,   # transmission FIXED at the most moderate level (no literature value exists)
    "transmission_mode": "frequency",
    "epsilon": 1e-6,
    "mu_A": 0.8,
    "gamma_N": 0.0,
    "invasive_carrier_fitness_cost": 0.0,
    "prob_recruit_carrier": 0.2,

    # --- Treatment: FIXED at high effort for this whole script ---
    "efficacy_I": 1.0,
    "efficacy_N_protect": 1.0,
    "overall_encounter_rate": 0.8,
    "mix_fraction_kill": 0.8,
    "migration_rate": 0.0,  # migration off by default -- see the model engine (engine.py) module docstring
    "treatment_months": list(range(1, 13)),

    # --- Simulation settings ---
    "n_years": 10,
    "steps_per_year": 12,
    "seed": 42,
}


# ---------------------------------------------------------------------------
# 3. Parameter test ranges (3-5 values each, baseline included)
# ---------------------------------------------------------------------------

PARAMETER_RANGES: Dict[str, List[float]] = {
    "I0": [300, 1000, 2500, 5000],            # baseline=1000 (=0.2 x K_I of Favors Native habitat)
    "N0": [400, 800, 1280, 2000],             # baseline=1280 (=0.2 x K_N of Favors Native habitat)
    "p_C0": [0.1, 0.2, 0.5, 0.8],           # France average ~0.2, documented range up to 0.8
    "r_I": [0.25, 0.35, 0.45, 0.60],
    "r_N": [0.06, 0.09, 0.12, 0.18],
    "K_I": [5000, 6400, 8000, 11000],       # baseline=5000 (Favors Native habitat); density-proxy uncertainty range
    "K_N": [1400, 3000, 5000, 6400],        # baseline=6400 (Favors Native habitat); refuge-site estimate up to near-parity
    "m_I": [0.04, 0.08, 0.12, 0.16],
    "m_N": [0.005, 0.01, 0.02, 0.03],
    "alpha_IN": [0.1, 0.3, 0.5, 0.7],
    "alpha_NI": [0.3, 0.5, 0.7, 0.9],
    "beta": [0.1, 0.3, 0.5, 1.0, 2.0],       # fully uncalibrated -- wide illustrative range
    "mu_A": [0.3, 0.5, 0.65, 0.8, 0.95],      # documented resistant-to-severe range, baseline=0.8
    "gamma_N": [0.0, 0.05, 0.1, 0.2],        # not well documented -- exploratory
    "invasive_carrier_fitness_cost": [0.0, 0.01, 0.02, 0.05],
    "prob_recruit_carrier": [0.1, 0.2, 0.35, 0.5],
    "migration_rate": [0.0, 0.02, 0.05, 0.1, 0.2],  # baseline=0.0 (off); 0.05 = the model engine (engine.py)'s illustrative literature-informed "on" default (Bubb et al. 2004/2006, no precisely fitted value -- see the model engine (engine.py) module docstring)
}


# ---------------------------------------------------------------------------
# 4. One-at-a-time sweep
# ---------------------------------------------------------------------------

def run_one_at_a_time_sweep(
    base_params: Dict[str, Any] = None,
    parameter_ranges: Dict[str, List[float]] = None,
    n_replicates: int = 30,
    base_seed: int = 50000,
) -> pd.DataFrame:
    """
    For each parameter in `parameter_ranges`, vary ONLY that parameter
    across its listed values (holding every other parameter at
    `base_params`, including the fixed high-effort treatment settings),
    running `n_replicates` stochastic replicates per value.

    Parameters
    ----------
    base_params : dict
        Baseline parameter dictionary (defaults to BASE_PARAMS).
    parameter_ranges : dict
        Mapping from parameter name -> list of test values (defaults to
        PARAMETER_RANGES).
    n_replicates : int
        Stochastic replicates averaged per tested value.
    base_seed : int
        Starting seed; incremented for every simulation run.

    Returns
    -------
    pd.DataFrame
        One row per (parameter, value, replicate), with columns
        'parameter', 'value', 'is_baseline', 'I_final', 'N_final',
        'native_persisted'.
    """
    if base_params is None:
        base_params = BASE_PARAMS
    if parameter_ranges is None:
        parameter_ranges = PARAMETER_RANGES

    records = []
    seed_counter = base_seed

    for param_name, values in parameter_ranges.items():
        baseline_value = base_params[param_name]
        for value in values:
            p = base_params.copy()
            p[param_name] = value
            for _ in range(n_replicates):
                res = engine.simulate_combined_treatment(p, seed=seed_counter)
                seed_counter += 1
                records.append({
                    "parameter": param_name,
                    "value": value,
                    "is_baseline": np.isclose(value, baseline_value),
                    "I_final": res["I"].iloc[-1],
                    "N_final": res["N"].iloc[-1],
                    "native_persisted": res["N"].iloc[-1] > 0,
                })
        print(f"  ...completed parameter '{param_name}' ({len(values)} values x {n_replicates} replicates)")

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# 5. Dashboard: one dose-response panel per parameter (SHARED y-axis)
# ---------------------------------------------------------------------------

def plot_parameter_dashboard(
    sweep_results: pd.DataFrame,
    base_params: Dict[str, Any] = None,
    y_range: tuple | None = None,
) -> None:
    """
    Grid of small panels, one per parameter, each showing median final
    native population (with 25th-75th percentile band) against that
    parameter's tested values. The baseline value is marked with a
    vertical dashed line.

    All panels share the SAME y-axis range, so the steepness/impact of
    different parameters can be compared directly at a glance.

    Parameters
    ----------
    y_range : tuple(float, float), optional
        Manually fix the shared y-axis to (ymin, ymax). If None (default),
        the range runs from 0 to the single highest native population
        value found anywhere across all parameters/replicates (i.e. the
        true max of the raw `N_final` column, not just the percentile
        band), with a small padding margin on top.
    """
    if base_params is None:
        base_params = BASE_PARAMS

    parameters = list(sweep_results["parameter"].unique())
    n_params = len(parameters)
    n_cols = 4
    n_rows = int(np.ceil(n_params / n_cols))

    # --- Pre-compute per-parameter summaries once, so we can reuse them
    # for plotting without recomputing. ---
    summaries = {}
    for param_name in parameters:
        sub = sweep_results[sweep_results["parameter"] == param_name]
        summary = sub.groupby("value").agg(
            N_p25=("N_final", lambda x: np.percentile(x, 25)),
            N_median=("N_final", "median"),
            N_p75=("N_final", lambda x: np.percentile(x, 75)),
        ).reset_index().sort_values("value")
        summaries[param_name] = summary

    # --- Determine shared y-axis range: 0 to the highest population value
    # found anywhere in the raw results (across all parameters and
    # replicates), with a small padding margin on top. ---
    if y_range is None:
        max_pop = sweep_results["N_final"].max()
        pad = 0.05 * max_pop if max_pop > 0 else 1.0
        y_min, y_max = 0, max_pop + pad
    else:
        y_min, y_max = y_range

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.6 * n_cols, 3.8 * n_rows))
    axes = axes.flatten()

    for i, param_name in enumerate(parameters):
        ax = axes[i]
        summary = summaries[param_name]

        ax.fill_between(summary["value"], summary["N_p25"], summary["N_p75"], color="steelblue", alpha=0.2)
        ax.plot(summary["value"], summary["N_median"], marker="o", color="steelblue", linewidth=2, markersize=5)

        baseline_value = base_params[param_name]
        ax.axvline(baseline_value, color="gray", linestyle="--", linewidth=1.2, alpha=0.7)

        ax.set_ylim(y_min, y_max)  # <-- shared range applied to every panel

        ax.set_title(param_name, fontsize=11, fontweight="bold")
        ax.set_xlabel("value tested", fontsize=8.5)
        ax.set_ylabel("median native pop.", fontsize=8.5)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.25)

    for j in range(n_params, len(axes)):
        axes[j].axis("off")

    fig.suptitle(
        "One-at-a-time parameter exploration -- effect on NATIVE population (year 10)\n"
        "Treatment fixed at: efficacy=100%, coverage=80%, mix_fraction_kill=0.8  |  "
        "dashed line = baseline value  |  band = 25th-75th percentile across replicates  |  "
        "y-axis shared across all panels",
        fontsize=12.5,
    )
    fig.subplots_adjust(hspace=0.55, wspace=0.4, top=0.90, bottom=0.05, left=0.05, right=0.98)
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()


# ---------------------------------------------------------------------------
# 6. Tornado chart: overall impact range per parameter
# ---------------------------------------------------------------------------

def plot_tornado_chart(sweep_results: pd.DataFrame, metric: str = "N_final", title_suffix: str = "native population") -> None:
    """
    Horizontal bar chart ranking parameters by how much their tested range
    shifts the median outcome -- a classic "tornado diagram" showing which
    parameters matter most for the given metric.

    Parameters
    ----------
    sweep_results : pd.DataFrame
        Output of `run_one_at_a_time_sweep`.
    metric : str
        Which outcome column to summarize ('N_final' or 'I_final').
    title_suffix : str
        Human-readable label for the chart title/axis.
    """
    summary = sweep_results.groupby(["parameter", "value"]).agg(
        median_outcome=(metric, "median"),
    ).reset_index()

    ranges = summary.groupby("parameter").agg(
        min_outcome=("median_outcome", "min"),
        max_outcome=("median_outcome", "max"),
    ).reset_index()
    ranges["impact_range"] = ranges["max_outcome"] - ranges["min_outcome"]
    ranges = ranges.sort_values("impact_range", ascending=True)

    plt.figure(figsize=(9, 7.5))
    y_pos = np.arange(len(ranges))
    plt.barh(y_pos, ranges["impact_range"], color="darkorange", alpha=0.85)
    plt.yticks(y_pos, ranges["parameter"], fontsize=10)
    plt.xlabel(f"Range in median final {title_suffix}\n(max tested value - min tested value)", fontsize=10.5)
    plt.title(
        f"Parameter impact ranking (tornado chart) -- {title_suffix}\n"
        "Treatment fixed at: efficacy=100%, coverage=80%, mix_fraction_kill=0.8",
        fontsize=12, fontweight="bold",
    )
    plt.grid(alpha=0.3, axis="x")
    for i, (idx, row) in enumerate(ranges.iterrows()):
        plt.text(row["impact_range"], i, f"  {row['impact_range']:,.0f}", va="center", fontsize=9)
    plt.tight_layout()
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()


def plot_wiki_tornado_chart(sweep_results: pd.DataFrame, metric: str = "N_final", save_path: str = None) -> None:
    """
    Publication-quality tornado chart for the wiki's "Other Explorations"
    section: ranks every tested parameter by how much its tested range
    shifts the median final native population, most influential at top.
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 13,
        "axes.linewidth": 1.1,
    })

    summary = sweep_results.groupby(["parameter", "value"]).agg(
        median_outcome=(metric, "median"),
    ).reset_index()
    ranges = summary.groupby("parameter").agg(
        min_outcome=("median_outcome", "min"),
        max_outcome=("median_outcome", "max"),
    ).reset_index()
    ranges["impact_range"] = ranges["max_outcome"] - ranges["min_outcome"]
    ranges = ranges.sort_values("impact_range", ascending=True)

    fig, ax = plt.subplots(figsize=(9.5, 8))
    y_pos = np.arange(len(ranges))
    ax.barh(y_pos, ranges["impact_range"], color="#C1780A", alpha=0.9)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ranges["parameter"], fontsize=12)
    ax.set_xlabel("Range in median final native population\n(highest tested value \u2212 lowest tested value)", fontsize=13, labelpad=8)
    ax.set_title("Which parameters move the outcome most?", fontsize=15.5, fontweight="bold", pad=14)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, axis="x", linewidth=0.7)
    ax.tick_params(labelsize=12)

    for i, (idx, row) in enumerate(ranges.iterrows()):
        ax.text(row["impact_range"], i, f"  {row['impact_range']:,.0f}", va="center", fontsize=11)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  ...saved {save_path}")
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()


# ---------------------------------------------------------------------------
# 7. Sanity checks
# ---------------------------------------------------------------------------

def run_sanity_checks() -> None:
    """Run a small, fast version of the full pipeline to confirm it works end-to-end."""
    import matplotlib
    original_backend = matplotlib.get_backend()
    matplotlib.use("Agg")
    try:
        small_ranges = {
            "p_C0": [0.1, 0.5],
            "mu_A": [0.3, 0.95],
            "K_I": [5000, 11000],
        }
        sweep = run_one_at_a_time_sweep(BASE_PARAMS, small_ranges, n_replicates=5, base_seed=1)
        assert len(sweep) == sum(len(v) for v in small_ranges.values()) * 5, "Unexpected row count."
        assert (sweep["I_final"] >= 0).all() and (sweep["N_final"] >= 0).all(), "Negative populations found."
        assert set(sweep["parameter"].unique()) == set(small_ranges.keys())

        plot_parameter_dashboard(sweep, BASE_PARAMS)
        plot_tornado_chart(sweep, "N_final", "native population")

        plt.close("all")
    finally:
        matplotlib.use(original_backend)

    print("All sensitivity_analysis.py (one-at-a-time parameter exploration) sanity checks passed.")


# ---------------------------------------------------------------------------
# 8. Main execution example
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    run_sanity_checks()

    print("\nRunning full one-at-a-time parameter sweep (this takes a moment)...")
    sweep_results = run_one_at_a_time_sweep(BASE_PARAMS, PARAMETER_RANGES, n_replicates=30)

    print("\nGenerating tornado chart (wiki figure)...")
    plot_wiki_tornado_chart(sweep_results, "N_final",
                            save_path=os.path.join(FIG_DIR, "sensitivity_tornado_wiki.png"))
