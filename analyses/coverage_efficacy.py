"""
coverage_efficacy.py: Coverage x Efficacy -- absolute final native population, per strategy
================================================================================

Purpose
-------
Extends coverage_efficacy.py (which showed the RNAi ADVANTAGE, i.e. the difference
between strategies) to instead show the ACTUAL final native population
for each strategy directly, across the full coverage x efficacy grid.
Three panels side by side (Defensive (protein), Mixed, Disarming (RNAi)), each colored
by median final native population -- not a comparison between them.

Fixed scenario (same as main results / validation / coverage_efficacy.py)
--------------------------------------------------------------------------
    Habitat      : Neutral   (K_I=6400, K_N=5500)
    Competition  : Moderate  (alpha_IN=0.3, alpha_NI=0.5)
    Transmission : beta = 0.9
    N0=0.2xK_N, I0=0.2xK_I

What is swept
--------------------------------------------------------------------------
    Coverage (overall_encounter_rate) : 0.2, 0.4, 0.6, 0.8, 1.0
    Efficacy (efficacy_I AND efficacy_N_protect, varied together) :
              0.2, 0.4, 0.6, 0.8, 1.0

Output
------
ONE figure, 3 heatmap panels (Defensive (protein) / Mixed / Disarming (RNAi)), each
efficacy (rows) x coverage (columns), colored by median final native
population, sharing the same color scale across all 3 panels so they
are directly comparable to each other.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, Any, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))  # lets the script run without `pip install -e .`

from stochastic_treatment_efficacy_model import engine  # noqa: E402

FIG_DIR = REPO_ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

K_I_FIXED = 6400.0
K_N_FIXED = 5500.0
ALPHA_IN_FIXED = 0.3
ALPHA_NI_FIXED = 0.5
BETA_FIXED = 0.9

N0_FIXED = 0.2 * K_N_FIXED
I0_FIXED = 0.2 * K_I_FIXED

BASE_PARAMS: Dict[str, Any] = {
    "I0": I0_FIXED, "N0": N0_FIXED,
    "p_C0": 0.2, "N_I0": 0,
    "r_I": 0.45, "r_N": 0.12,
    "K_I": K_I_FIXED, "K_N": K_N_FIXED,
    "m_I": 0.08, "m_N": 0.01,
    "alpha_IN": ALPHA_IN_FIXED, "alpha_NI": ALPHA_NI_FIXED,
    "beta": BETA_FIXED,
    "transmission_mode": "frequency",
    "epsilon": 1e-6,
    "mu_A": 0.8,
    "gamma_N": 0.0,
    "invasive_carrier_fitness_cost": 0.0,
    "prob_recruit_carrier": 0.2,
    "migration_rate": 0.0,
    "efficacy_I": None,
    "efficacy_N_protect": None,
    "overall_encounter_rate": None,
    "mix_fraction_kill": None,
    "treatment_months": list(range(1, 13)),
    "n_years": 10,
    "steps_per_year": 12,
    "seed": 42,
}

COVERAGE_VALUES: List[float] = [0.2, 0.4, 0.6, 0.8, 1.0]
EFFICACY_VALUES: List[float] = [0.2, 0.4, 0.6, 0.8, 1.0]

STRATEGIES: List[Dict[str, Any]] = [
    {"label": "Protein Supplementation", "mix": 0.0},
    {"label": "Mixed (mix=0.5)",          "mix": 0.5},
    {"label": "RNAi strategy",            "mix": 1.0},
]


def build_params(coverage: float, efficacy: float, mix: float) -> Dict[str, Any]:
    p = BASE_PARAMS.copy()
    p["overall_encounter_rate"] = coverage
    p["efficacy_I"] = efficacy
    p["efficacy_N_protect"] = efficacy
    p["mix_fraction_kill"] = mix
    return p


def run_sweep(n_replicates: int = 30, base_seed: int = 290000) -> pd.DataFrame:
    """
    Run all 3 strategies across the full coverage x efficacy grid.

    Returns
    -------
    pd.DataFrame with columns 'coverage', 'efficacy', 'strategy', 'replicate', 'N_final'.
    """
    records = []
    seed_counter = base_seed
    total = len(COVERAGE_VALUES) * len(EFFICACY_VALUES) * len(STRATEGIES)
    done = 0

    for coverage in COVERAGE_VALUES:
        for efficacy in EFFICACY_VALUES:
            for strat in STRATEGIES:
                p = build_params(coverage, efficacy, strat["mix"])
                for r in range(n_replicates):
                    res = engine.simulate_combined_treatment(p, seed=seed_counter)
                    seed_counter += 1
                    records.append({
                        "coverage": coverage, "efficacy": efficacy, "strategy": strat["label"],
                        "replicate": r, "N_final": res["N"].iloc[-1],
                    })
                done += 1
                print(f"  ...completed coverage={coverage} x efficacy={efficacy} x {strat['label']} ({done}/{total})")

    return pd.DataFrame(records)


def plot_wiki_figure(results: pd.DataFrame, save_path: str = None) -> None:
    """3 heatmap panels (Defensive (protein) / Mixed / Disarming (RNAi)), absolute median final native population, shared color scale."""
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 13, "axes.linewidth": 0})

    summary = results.groupby(["coverage", "efficacy", "strategy"]).agg(median_N=("N_final", "median")).reset_index()
    vmax = summary["median_N"].max()

    fig, axes = plt.subplots(1, 3, figsize=(17, 6.2))
    cmap = plt.get_cmap("YlGnBu")
    mesh = None

    for i, strat in enumerate(STRATEGIES):
        ax = axes[i]
        sub = summary[summary["strategy"] == strat["label"]]
        grid = sub.pivot(index="efficacy", columns="coverage", values="median_N")
        grid = grid.loc[sorted(EFFICACY_VALUES, reverse=True), COVERAGE_VALUES]
        n_rows, n_cols = grid.shape

        mesh = ax.pcolormesh(grid.values, cmap=cmap, vmin=0, vmax=vmax, edgecolors="white", linewidth=4)

        ax.set_xlim(0, n_cols)
        ax.set_ylim(0, n_rows)
        ax.invert_yaxis()
        ax.set_xticks(np.arange(n_cols) + 0.5)
        ax.set_xticklabels([f"{int(c*100)}%" for c in COVERAGE_VALUES], fontsize=11.5)
        if i == 0:
            ax.set_yticks(np.arange(n_rows) + 0.5)
            ax.set_yticklabels([f"{int(e*100)}%" for e in sorted(EFFICACY_VALUES, reverse=True)], fontsize=11.5)
            ax.set_ylabel("Efficacy", fontsize=14, labelpad=10)
        else:
            ax.set_yticks([])
        ax.set_xlabel("Coverage", fontsize=14, labelpad=8)
        ax.set_title(strat["label"], fontsize=14.5, fontweight="bold", color="#1F3864", pad=12)
        ax.tick_params(length=0)

        for r in range(n_rows):
            for c in range(n_cols):
                val = grid.values[r, c]
                text_color = "white" if val > vmax * 0.55 else "#1a1a1a"
                ax.text(c + 0.5, r + 0.5, f"{val:,.0f}", ha="center", va="center",
                         fontsize=12.5, fontweight="bold", color=text_color)

        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.subplots_adjust(top=0.78, bottom=0.12, left=0.06, right=0.90, wspace=0.08)
    cbar_ax = fig.add_axes([0.925, 0.14, 0.018, 0.60])
    cbar = fig.colorbar(mesh, cax=cbar_ax)
    cbar.set_label("Median final native population", fontsize=12.5, labelpad=12)
    cbar.ax.tick_params(labelsize=11.5)
    cbar.outline.set_visible(False)

    fig.suptitle("Final native population across coverage and efficacy, by strategy",
                 fontsize=17, fontweight="bold", color="#1F3864", y=1.04)
    fig.text(0.44, 0.965,
        "Neutral habitat, moderate competition, \u03b2 = 0.9 (same as main results)",
        ha="center", va="top", fontsize=11, color="#444444")

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  ...saved {save_path}")
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()


def run_sanity_checks() -> None:
    import matplotlib
    original_backend = matplotlib.get_backend()
    matplotlib.use("Agg")
    global COVERAGE_VALUES, EFFICACY_VALUES
    original_cov, original_eff = COVERAGE_VALUES, EFFICACY_VALUES
    try:
        COVERAGE_VALUES = [0.2, 1.0]
        EFFICACY_VALUES = [0.2, 1.0]
        small_results = run_sweep(n_replicates=3, base_seed=1)
        expected_rows = len(COVERAGE_VALUES) * len(EFFICACY_VALUES) * len(STRATEGIES) * 3
        assert len(small_results) == expected_rows, "Unexpected row count."
        assert (small_results["N_final"] >= 0).all(), "Negative native populations found."
        plot_wiki_figure(small_results)
        plt.close("all")
    finally:
        COVERAGE_VALUES, EFFICACY_VALUES = original_cov, original_eff
        matplotlib.use(original_backend)
    print("All coverage_efficacy.py (absolute population, coverage x efficacy) sanity checks passed.")


if __name__ == "__main__":
    run_sanity_checks()
    print(f"\nRunning coverage x efficacy sweep ({len(COVERAGE_VALUES)}x{len(EFFICACY_VALUES)} grid, {len(STRATEGIES)} strategies)...")
    results = run_sweep(n_replicates=30)
    save_path = os.path.join(FIG_DIR, "other_explorations_coverage_efficacy_absolute.png")
    print("\nGenerating figure...")
    plot_wiki_figure(results, save_path=save_path)
