"""
validation_bubb2006.py: Validation against the Bubb et al. (2006) displacement -- final consolidated version
================================================================================

Purpose
-------
QUALITATIVE VALIDATION, not quantitative fitting to an external dataset
(competition and transmission have no independent literature value for
this species/pathogen system -- see the project wiki Model page,
so a formal statistical fit isn't possible). This script:

  1. Sweeps a GRID of competition x transmission combinations, under our
     Neutral habitat, no treatment, over the same 2-year window as a
     real, documented displacement, and reports which combination(s)
     come closest to matching that real outcome (FIGURE 1).

  2. Runs the SELECTED scenario (chosen after reviewing the grid --
     moderate competition, beta=0.9, see "Scenario selection" below) as
     a full trajectory, showing BOTH native and invasive population
     curves over the 2-year window (FIGURE 2).

The real-world reference case
--------------------------------------------------------------------------
Bubb, D., Thom, T.J. & Lucas, M.C. (2006) Movement, dispersal and refuge
use of co-occurring introduced and native crayfish. Freshwater Biology,
51, 1359-1368.

Invasive and native crayfish were co-occurring at one river site in 2003.
A resurvey of the SAME site in 2005 (2 years later) recorded ZERO native
crayfish, while invasive density had increased. Target for this
validation: ~100% native decline within 2 years, no treatment, alongside
a growing invasive population.

Scenario selection (why moderate competition, not the closest-fitting grid cell)
--------------------------------------------------------------------------
The single closest match across the full grid was HIGH competition +
beta=0.9 (98% decline). We deliberately did NOT select it: our habitat is
fixed at Neutral throughout this project specifically to test our
strategies under realistic, non-extreme conditions, and pairing a
Neutral habitat with a high competition coefficient would be internally
inconsistent with every other result on this page. MODERATE competition
+ beta=0.9 (95% decline) is the combination we adopted instead -- still
a strong match to the real outcome, and consistent with the baseline
used everywhere else in this project. This is also the beta value used
in the main results section.

Design notes
--------------------------------------------------------------------------
Habitat is FIXED at Neutral (K_I=6400, K_N=5500) throughout -- the SAME
habitat used in the main results section, never swept. N0=0.2xK_N and
I0=0.2xK_I, the same convention used throughout the main results.

Output
------
FIGURE 1: heatmap, competition (rows) x transmission (columns), colored
by % decline in native population over the 2-year window.

FIGURE 2: trajectory over the 2-year window for the SELECTED scenario
(moderate competition, beta=0.9, no treatment), showing BOTH native and
invasive population curves.

Both figures are publication-quality (300 DPI, clean typography, saved
as PNG files next to this script).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, Any, List, Tuple

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
# 2. Validation window -- matches the real 2003-2005 survey gap exactly
# ---------------------------------------------------------------------------

VALIDATION_YEARS = 2
STEPS_PER_YEAR = 12
NEAR_EXTINCT_THRESHOLD = 5


# ---------------------------------------------------------------------------
# 3. Fixed scenario -- Neutral habitat, IDENTICAL to the main results
# ---------------------------------------------------------------------------

K_I_FIXED = 6400.0
K_N_FIXED = 5500.0

N0_FIXED = 0.2 * K_N_FIXED
I0_FIXED = 0.2 * K_I_FIXED

BASE_PARAMS: Dict[str, Any] = {
    "I0": I0_FIXED,
    "N0": N0_FIXED,
    "p_C0": 0.2,
    "N_I0": 0,

    "r_I": 0.45,
    "r_N": 0.12,
    "K_I": K_I_FIXED,
    "K_N": K_N_FIXED,
    "m_I": 0.08,
    "m_N": 0.01,

    "alpha_IN": None,   # PLACEHOLDER: set per competition scenario
    "alpha_NI": None,
    "beta": None,       # PLACEHOLDER: set per transmission value
    "transmission_mode": "frequency",
    "epsilon": 1e-6,
    "mu_A": 0.8,
    "gamma_N": 0.0,
    "invasive_carrier_fitness_cost": 0.0,
    "prob_recruit_carrier": 0.2,
    "migration_rate": 0.0,

    "efficacy_I": 1.0,
    "efficacy_N_protect": 1.0,
    "overall_encounter_rate": 0.0,
    "mix_fraction_kill": 0.0,
    "treatment_months": list(range(1, 13)),

    "n_years": VALIDATION_YEARS,
    "steps_per_year": STEPS_PER_YEAR,
    "seed": 42,
}

COMPETITION_VALUES_TO_TEST: Dict[str, Tuple[float, float]] = {
    "Low":      (0.1, 0.3),
    "Moderate": (0.3, 0.5),
    "High":     (0.7, 0.9),
}
COMPETITION_ORDER: List[str] = ["Low", "Moderate", "High"]
BETA_VALUES_TO_TEST: List[float] = [0.5, 0.6, 0.7, 0.8, 0.9]

# --- The scenario selected for Figure 2, after reviewing the grid --
# see "Scenario selection" in the module docstring. Change these to
# highlight a different combination. ---
SELECTED_COMPETITION = "Moderate"
SELECTED_BETA = 0.9


def build_params(comp_name: str, beta_value: float) -> Dict[str, Any]:
    p = BASE_PARAMS.copy()
    alpha_in, alpha_ni = COMPETITION_VALUES_TO_TEST[comp_name]
    p["alpha_IN"] = alpha_in
    p["alpha_NI"] = alpha_ni
    p["beta"] = beta_value
    return p


# ---------------------------------------------------------------------------
# 4. Publication styling (applied to both figures)
# ---------------------------------------------------------------------------

def _apply_publication_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 13,
        "axes.linewidth": 1.1,
    })


def _clean_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, linewidth=0.7)
    ax.tick_params(labelsize=12)


# ---------------------------------------------------------------------------
# 5. FIGURE 1: validation grid (competition x transmission)
# ---------------------------------------------------------------------------

def run_grid(n_replicates: int = 30, base_seed: int = 210000) -> pd.DataFrame:
    """
    Run the full (competition x transmission) grid, no treatment, over
    the 2-year validation window.

    Returns
    -------
    pd.DataFrame with columns 'competition', 'beta', 'replicate', 'N_final', 'pct_decline'.
    """
    records = []
    seed_counter = base_seed
    total = len(COMPETITION_ORDER) * len(BETA_VALUES_TO_TEST)
    done = 0

    for comp_name in COMPETITION_ORDER:
        for beta_value in BETA_VALUES_TO_TEST:
            p = build_params(comp_name, beta_value)
            for r in range(n_replicates):
                res = engine.simulate_combined_treatment(p, seed=seed_counter)
                seed_counter += 1
                n_final = res["N"].iloc[-1]
                pct_decline = 100 * (N0_FIXED - n_final) / N0_FIXED if N0_FIXED > 0 else np.nan
                records.append({
                    "competition": comp_name, "beta": beta_value, "replicate": r,
                    "N_final": n_final, "pct_decline": pct_decline,
                })
            done += 1
            print(f"  ...completed {comp_name} competition x beta={beta_value} ({done}/{total})")

    return pd.DataFrame(records)


def plot_fit_heatmap(grid_results: pd.DataFrame, save_path: str = None) -> None:
    """FIGURE 1: heatmap, competition (rows) x transmission (columns), colored by median % decline."""
    _apply_publication_style()

    summary = grid_results.groupby(["competition", "beta"]).agg(
        median_pct_decline=("pct_decline", "median"),
    ).reset_index()
    grid = summary.pivot(index="competition", columns="beta", values="median_pct_decline")
    grid = grid.loc[COMPETITION_ORDER, BETA_VALUES_TO_TEST]

    fig, ax = plt.subplots(figsize=(9.5, 6), dpi=100)
    im = ax.imshow(grid.values, aspect="auto", cmap="RdYlGn", vmin=0, vmax=100, origin="upper")

    ax.set_xticks(range(len(BETA_VALUES_TO_TEST)))
    ax.set_xticklabels(BETA_VALUES_TO_TEST, fontsize=12)
    ax.set_yticks(range(len(COMPETITION_ORDER)))
    ax.set_yticklabels(COMPETITION_ORDER, fontsize=12)
    ax.set_xlabel("\u03b2 (transmission)", fontsize=14, labelpad=8)
    ax.set_ylabel("Competition level", fontsize=14, labelpad=8)

    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            val = grid.values[i, j]
            is_selected = (COMPETITION_ORDER[i] == SELECTED_COMPETITION and BETA_VALUES_TO_TEST[j] == SELECTED_BETA)
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center", fontsize=11.5,
                     color="white" if val > 60 or val < 25 else "black",
                     fontweight="bold" if is_selected else "normal")
            if is_selected:
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="#1F3864", linewidth=3)
                ax.add_patch(rect)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Median % decline in native population over 2 years\n(100% = matches real-site outcome)", fontsize=11)
    cbar.ax.tick_params(labelsize=11)

    fig.suptitle(
        "Validation grid: does the model reproduce the documented displacement?",
        fontsize=15.5, fontweight="bold", color="#1F3864", y=0.99,
    )
    fig.text(0.5, 0.905,
        "Reference: Bubb, Thom & Lucas (2006) -- 100% native decline over 2 years, no treatment\n"
        f"Neutral habitat (K_I={K_I_FIXED:.0f}, K_N={K_N_FIXED:.0f}, fixed, same as main results)  |  "
        "outlined cell = scenario selected for Figure 2",
        ha="center", va="top", fontsize=10.5,
    )
    fig.subplots_adjust(top=0.78, bottom=0.12, left=0.13, right=0.98)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  ...saved {save_path}")
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()

    summary["distance_from_target"] = (100 - summary["median_pct_decline"]).abs()
    best_row = summary.loc[summary["distance_from_target"].idxmin()]
    print("\n" + "=" * 70)
    print(f"Closest grid match: {best_row['competition']} competition, beta={best_row['beta']} "
          f"({best_row['median_pct_decline']:.0f}% decline)")
    selected_row = summary[(summary["competition"] == SELECTED_COMPETITION) & (summary["beta"] == SELECTED_BETA)].iloc[0]
    print(f"Selected for Figure 2: {SELECTED_COMPETITION} competition, beta={SELECTED_BETA} "
          f"({selected_row['median_pct_decline']:.0f}% decline) -- see module docstring for why.")
    print("=" * 70 + "\n")


# ---------------------------------------------------------------------------
# 6. FIGURE 2: selected-scenario trajectory, native AND invasive
# ---------------------------------------------------------------------------

def run_selected_scenario(n_replicates: int = 50, base_seed: int = 215000) -> pd.DataFrame:
    """
    Run the SELECTED (competition, beta) scenario, tracking BOTH native
    and invasive population trajectories over the 2-year window.

    Returns
    -------
    pd.DataFrame with columns 'replicate', 'month', 'year', 'N', 'I'.
    """
    p = build_params(SELECTED_COMPETITION, SELECTED_BETA)
    frames = []
    for r in range(n_replicates):
        res = engine.simulate_combined_treatment(p, seed=base_seed + r)
        sub = res[["month", "year", "N", "I"]].copy()
        sub["replicate"] = r
        frames.append(sub)
    print(f"  ...completed {SELECTED_COMPETITION} competition, beta={SELECTED_BETA} ({n_replicates} replicates)")
    return pd.concat(frames, ignore_index=True)


def plot_selected_scenario(trajectories: pd.DataFrame, save_path: str = None) -> None:
    """FIGURE 2: native AND invasive population over the 2-year window, for the selected scenario."""
    _apply_publication_style()

    native_summary = trajectories.groupby("year").agg(
        p10=("N", lambda x: np.percentile(x, 10)),
        median=("N", "median"),
        p90=("N", lambda x: np.percentile(x, 90)),
    ).reset_index()
    invasive_summary = trajectories.groupby("year").agg(
        p10=("I", lambda x: np.percentile(x, 10)),
        median=("I", "median"),
        p90=("I", lambda x: np.percentile(x, 90)),
    ).reset_index()

    fig, ax = plt.subplots(figsize=(10, 6.5), dpi=100)

    NATIVE_COLOR = "#1C7293"
    INVASIVE_COLOR = "#C1292E"

    ax.fill_between(native_summary["year"], native_summary["p10"], native_summary["p90"], color=NATIVE_COLOR, alpha=0.15, linewidth=0)
    ax.plot(native_summary["year"], native_summary["median"], color=NATIVE_COLOR, linewidth=2.6, label="Native population")

    ax.fill_between(invasive_summary["year"], invasive_summary["p10"], invasive_summary["p90"], color=INVASIVE_COLOR, alpha=0.15, linewidth=0)
    ax.plot(invasive_summary["year"], invasive_summary["median"], color=INVASIVE_COLOR, linewidth=2.6, label="Invasive population")

    ax.axhline(NEAR_EXTINCT_THRESHOLD, color="gray", linestyle=":", linewidth=1.2, label=f"Near-extinct threshold ({NEAR_EXTINCT_THRESHOLD})")

    ax.set_xlabel("Year", fontsize=14, labelpad=8)
    ax.set_ylabel("Population", fontsize=14, labelpad=8)
    ax.set_xlim(0, VALIDATION_YEARS)
    ax.set_ylim(bottom=0)

    _clean_axes(ax)
    legend = ax.legend(fontsize=12, loc="center left", frameon=True, framealpha=0.95, edgecolor="#DDDDDD")
    legend.get_frame().set_linewidth(0.8)

    fig.suptitle(
        f"Selected validation scenario: {SELECTED_COMPETITION.lower()} competition, \u03b2 = {SELECTED_BETA}",
        fontsize=15.5, fontweight="bold", color="#1F3864", y=0.99,
    )
    fig.text(0.5, 0.905,
        "Reference: Bubb, Thom & Lucas (2006) -- co-occurring site in 2003, zero natives / increased invasive density by 2005\n"
        "Neutral habitat (fixed, same as main results). No treatment applied.",
        ha="center", va="top", fontsize=10,
    )
    fig.subplots_adjust(top=0.78, bottom=0.11, left=0.10, right=0.97)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  ...saved {save_path}")
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()

    final_native = native_summary["median"].iloc[-1]
    final_invasive = invasive_summary["median"].iloc[-1]
    pct_decline = 100 * (N0_FIXED - final_native) / N0_FIXED
    pct_growth = 100 * (final_invasive - I0_FIXED) / I0_FIXED if I0_FIXED > 0 else float("nan")
    print(f"  Native: {pct_decline:.0f}% decline over {VALIDATION_YEARS} years (from {N0_FIXED:.0f} to {final_native:.0f})")
    print(f"  Invasive: {pct_growth:+.0f}% change over {VALIDATION_YEARS} years (from {I0_FIXED:.0f} to {final_invasive:.0f})")
    print(f"  Real site: natives -> 0, invasive density increased, over the same 2-year window.")


# ---------------------------------------------------------------------------
# 7. Sanity checks
# ---------------------------------------------------------------------------

def run_sanity_checks() -> None:
    """Run a small, fast version of the full pipeline to confirm it works end-to-end."""
    import matplotlib
    original_backend = matplotlib.get_backend()
    matplotlib.use("Agg")
    global COMPETITION_ORDER, BETA_VALUES_TO_TEST
    original_comp, original_beta = COMPETITION_ORDER, BETA_VALUES_TO_TEST
    try:
        COMPETITION_ORDER = ["Low", "Moderate", "High"]
        BETA_VALUES_TO_TEST = [0.5, 0.9]

        small_grid = run_grid(n_replicates=3, base_seed=1)
        expected_rows = len(COMPETITION_ORDER) * len(BETA_VALUES_TO_TEST) * 3
        assert len(small_grid) == expected_rows, "Unexpected row count."
        assert (small_grid["N_final"] >= 0).all(), "Negative native populations found."
        assert (small_grid["pct_decline"] <= 100.01).all(), "% decline should never exceed 100."

        import tempfile
        test_path1 = os.path.join(tempfile.gettempdir(), "_stage21_fig1_check.png")
        plot_fit_heatmap(small_grid, save_path=test_path1)
        assert os.path.exists(test_path1)
        os.remove(test_path1)

        small_selected = run_selected_scenario(n_replicates=3, base_seed=1)
        assert (small_selected["N"] >= 0).all() and (small_selected["I"] >= 0).all(), "Negative populations found."
        assert "I" in small_selected.columns and "N" in small_selected.columns

        test_path2 = os.path.join(tempfile.gettempdir(), "_stage21_fig2_check.png")
        plot_selected_scenario(small_selected, save_path=test_path2)
        assert os.path.exists(test_path2)
        os.remove(test_path2)

        plt.close("all")
    finally:
        COMPETITION_ORDER, BETA_VALUES_TO_TEST = original_comp, original_beta
        matplotlib.use(original_backend)

    print("All validation_bubb2006.py (Bubb et al. 2006 validation, final consolidated) sanity checks passed.")


# ---------------------------------------------------------------------------
# 8. Main execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    run_sanity_checks()


    print(f"\nRunning full validation grid ({len(COMPETITION_ORDER)} competition levels x "
          f"{len(BETA_VALUES_TO_TEST)} transmission values)...")
    grid_results = run_grid(n_replicates=30)
    plot_fit_heatmap(grid_results, save_path=os.path.join(FIG_DIR, "validation_figure1_grid.png"))

    print(f"\nRunning selected scenario ({SELECTED_COMPETITION} competition, beta={SELECTED_BETA})...")
    selected_traj = run_selected_scenario(n_replicates=50)
    plot_selected_scenario(selected_traj, save_path=os.path.join(FIG_DIR, "validation_figure2_selected_scenario.png"))
