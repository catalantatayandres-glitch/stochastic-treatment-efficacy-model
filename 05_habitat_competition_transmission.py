"""
05_habitat_competition_transmission.py: Habitat x Competition x Transmission -- 27-scenario summary (wiki figure)
================================================================================

Purpose
-------
For the wiki's "Other Explorations" section: does the RNAi-outperforms-
protein result hold across the full range of habitat, competition, and
transmission conditions we established for this model, not just the
single Neutral/moderate/beta=0.9 scenario used in the main Results?

This script tests all 27 combinations (3 habitats x 3 competition levels
x 3 transmission levels) and summarizes them in ONE compact figure: 3
heatmap panels side by side (one per habitat), each a 3x3 grid
(competition x transmission), colored by how much better the RNAi
(disarming (RNAi)) strategy performs than the protein supplementation
(defensive (protein)) strategy, at the final year.

Habitat set (3, matching the main results' OTHER habitats -- Neutral is
already covered in depth in the main Results/Validation sections, so
this panel focuses on the remaining ones)
--------------------------------------------------------------------------
    Favors Native   : K_I=5000,  K_N=6400
    Balanced        : K_I=8000,  K_N=3000
    Favors Invasive : K_I=11000, K_N=1400

Competition and transmission (3 levels each, the standard low/moderate/
high values used throughout this project)
--------------------------------------------------------------------------
    Competition : Low (0.1/0.3), Moderate (0.3/0.5), High (0.7/0.9)
    Transmission: Low (0.2), Moderate (0.5), High (2.0)

Both strategies tested at 100% efficacy, 50% coverage (same convention
used throughout this project). N0=0.2xK_N, I0=0.2xK_I per habitat.

Output
------
ONE figure: 3 heatmap panels (one per habitat), each showing, for every
(competition, transmission) combination, the RNAi advantage: median
final native population under disarming (RNAi) MINUS median final native
population under defensive (protein). Positive (green) = RNAi wins; negative
(red) = protein wins; the actual value shown is the size of that gap.
"""

from __future__ import annotations

import os
import importlib.util
from typing import Dict, Any, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. Dynamically import the model engine (01_model_engine.py)
# ---------------------------------------------------------------------------

def _load_stage6():
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "01_model_engine.py")
    spec = importlib.util.spec_from_file_location("stage6_combined_treatment", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


stage6 = _load_stage6()

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(FIG_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# 2. Scenario definitions
# ---------------------------------------------------------------------------

HABITAT_SCENARIOS: Dict[str, Dict[str, float]] = {
    "Favors Native":   {"K_I": 5000.0,  "K_N": 6400.0},
    "Balanced":        {"K_I": 8000.0,  "K_N": 3000.0},
    "Favors Invasive": {"K_I": 11000.0, "K_N": 1400.0},
}
HABITAT_ORDER = list(HABITAT_SCENARIOS.keys())

COMPETITION_SCENARIOS: Dict[str, Dict[str, float]] = {
    "Low":      {"alpha_IN": 0.1, "alpha_NI": 0.3},
    "Moderate": {"alpha_IN": 0.3, "alpha_NI": 0.5},
    "High":     {"alpha_IN": 0.7, "alpha_NI": 0.9},
}
COMPETITION_ORDER = list(COMPETITION_SCENARIOS.keys())

TRANSMISSION_SCENARIOS: Dict[str, Dict[str, float]] = {
    "Low":      {"beta": 0.2},
    "Moderate": {"beta": 0.5},
    "High":     {"beta": 2.0},
}
TRANSMISSION_ORDER = list(TRANSMISSION_SCENARIOS.keys())

N0_FRACTION_OF_KN: float = 0.2
I0_FRACTION_OF_KI: float = 0.2

BASE_PARAMS: Dict[str, Any] = {
    "I0": None,       # PLACEHOLDER: computed per habitat
    "N0": None,
    "p_C0": 0.2,
    "N_I0": 0,

    "r_I": 0.45,
    "r_N": 0.12,
    "K_I": None,      # PLACEHOLDER: set per habitat
    "K_N": None,
    "m_I": 0.08,
    "m_N": 0.01,

    "alpha_IN": None,  # PLACEHOLDER: set per competition scenario
    "alpha_NI": None,
    "beta": None,      # PLACEHOLDER: set per transmission scenario
    "transmission_mode": "frequency",
    "epsilon": 1e-6,
    "mu_A": 0.8,
    "gamma_N": 0.0,
    "invasive_carrier_fitness_cost": 0.0,
    "prob_recruit_carrier": 0.2,
    "migration_rate": 0.0,

    "efficacy_I": 1.0,
    "efficacy_N_protect": 1.0,
    "overall_encounter_rate": 0.5,
    "mix_fraction_kill": None,  # PLACEHOLDER: 0 (protect) or 1 (kill) per run
    "treatment_months": list(range(1, 13)),

    "n_years": 10,
    "steps_per_year": 12,
    "seed": 42,
}


def build_params(habitat_name: str, comp_name: str, trans_name: str, mix_fraction_kill: float) -> Dict[str, Any]:
    p = BASE_PARAMS.copy()
    p.update(HABITAT_SCENARIOS[habitat_name])
    p.update(COMPETITION_SCENARIOS[comp_name])
    p.update(TRANSMISSION_SCENARIOS[trans_name])
    p["N0"] = N0_FRACTION_OF_KN * p["K_N"]
    p["I0"] = I0_FRACTION_OF_KI * p["K_I"]
    p["mix_fraction_kill"] = mix_fraction_kill
    return p


# ---------------------------------------------------------------------------
# 3. Monte Carlo runner -- disarming (RNAi) and defensive (protein), all 27 scenarios
# ---------------------------------------------------------------------------

def run_all_scenarios(n_replicates: int = 30, base_seed: int = 240000) -> pd.DataFrame:
    """
    Run disarming (RNAi) (mix=1) and defensive (protein) (mix=0) for all 27 (habitat,
    competition, transmission) combinations.

    Returns
    -------
    pd.DataFrame
        One row per (habitat, competition, transmission, strategy,
        replicate), with column 'N_final'.
    """
    records = []
    seed_counter = base_seed
    total = len(HABITAT_ORDER) * len(COMPETITION_ORDER) * len(TRANSMISSION_ORDER)
    done = 0

    for habitat_name in HABITAT_ORDER:
        for comp_name in COMPETITION_ORDER:
            for trans_name in TRANSMISSION_ORDER:
                for strategy_name, mix in [("Defensive (protein)", 0.0), ("Disarming (RNAi)", 1.0)]:
                    p = build_params(habitat_name, comp_name, trans_name, mix)
                    for r in range(n_replicates):
                        res = stage6.simulate_combined_treatment(p, seed=seed_counter)
                        seed_counter += 1
                        records.append({
                            "habitat": habitat_name, "competition": comp_name, "transmission": trans_name,
                            "strategy": strategy_name, "replicate": r, "N_final": res["N"].iloc[-1],
                        })
                done += 1
                print(f"  ...completed {habitat_name} x {comp_name} x {trans_name} ({done}/{total})")

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# 4. Figure: 3 heatmap panels (one per habitat), RNAi advantage
# ---------------------------------------------------------------------------

def plot_wiki_summary_figure(results: pd.DataFrame, save_path: str = None) -> None:
    """
    ONE figure, 3 heatmap panels side by side (one per habitat), each a
    3x3 grid (competition x transmission), colored by the RNAi advantage:
    median final native population under disarming (RNAi) minus under defensive (protein).
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 13,
        "axes.linewidth": 1.1,
    })

    summary = results.groupby(["habitat", "competition", "transmission", "strategy"]).agg(
        median_N=("N_final", "median"),
    ).reset_index()
    pivoted = summary.pivot_table(index=["habitat", "competition", "transmission"], columns="strategy", values="median_N").reset_index()
    pivoted["advantage"] = pivoted["Disarming (RNAi)"] - pivoted["Defensive (protein)"]

    vmax = pivoted["advantage"].abs().max()
    vmax = max(vmax, 1)

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6), sharey=True)

    for i, habitat_name in enumerate(HABITAT_ORDER):
        ax = axes[i]
        sub = pivoted[pivoted["habitat"] == habitat_name]
        grid = sub.pivot(index="competition", columns="transmission", values="advantage")
        grid = grid.loc[COMPETITION_ORDER, TRANSMISSION_ORDER]

        im = ax.imshow(grid.values, aspect="auto", cmap="RdYlGn", vmin=-vmax, vmax=vmax, origin="upper")

        ax.set_xticks(range(len(TRANSMISSION_ORDER)))
        ax.set_xticklabels(TRANSMISSION_ORDER, fontsize=11.5)
        if i == 0:
            ax.set_yticks(range(len(COMPETITION_ORDER)))
            ax.set_yticklabels(COMPETITION_ORDER, fontsize=11.5)
            ax.set_ylabel("Competition", fontsize=13)
        ax.set_xlabel("Transmission", fontsize=13)
        ax.set_title(habitat_name, fontsize=14, fontweight="bold")

        for r in range(grid.shape[0]):
            for c in range(grid.shape[1]):
                val = grid.values[r, c]
                ax.text(c, r, f"{val:+,.0f}", ha="center", va="center", fontsize=10.5,
                         color="white" if abs(val) > vmax * 0.5 else "black")

    fig.suptitle("Does RNAi still outperform protein supplementation everywhere?",
                 fontsize=16, fontweight="bold", color="#1F3864", y=1.02)
    fig.text(0.5, 0.965,
        "RNAi advantage (disarming (RNAi) minus defensive (protein) median final native population), across all 27 tested scenarios\n"
        "Green = RNAi wins by more; red = protein wins by more. Efficacy=100%, coverage=50%.",
        ha="center", va="top", fontsize=10.5)

    cbar = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cbar.set_label("RNAi advantage\n(disarming (RNAi) \u2212 defensive (protein))", fontsize=11)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  ...saved {save_path}")
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()

    print("\n" + "=" * 70)
    n_wins = (pivoted["advantage"] > 0).sum()
    n_ties = (pivoted["advantage"] == 0).sum()
    n_losses = (pivoted["advantage"] < 0).sum()
    n_total = len(pivoted)
    print(f"RNAi (disarming) strictly outperformed protein (defensive (protein)) in {n_wins}/{n_total} scenarios.")
    print(f"Tied in {n_ties}/{n_total} scenarios (both strategies collapsed to zero -- a floor effect, not a protein advantage).")
    print(f"RNAi underperformed protein in {n_losses}/{n_total} scenarios.")
    print("=" * 70)


# ---------------------------------------------------------------------------
# 5. Sanity checks
# ---------------------------------------------------------------------------

def run_sanity_checks() -> None:
    """Run a small, fast version of the full pipeline to confirm it works end-to-end."""
    import matplotlib
    original_backend = matplotlib.get_backend()
    matplotlib.use("Agg")
    try:
        small_results = run_all_scenarios(n_replicates=2, base_seed=1)
        expected_rows = len(HABITAT_ORDER) * len(COMPETITION_ORDER) * len(TRANSMISSION_ORDER) * 2 * 2
        assert len(small_results) == expected_rows, "Unexpected row count."
        assert (small_results["N_final"] >= 0).all(), "Negative native populations found."

        plot_wiki_summary_figure(small_results)
        plt.close("all")
    finally:
        matplotlib.use(original_backend)

    print("All 05_habitat_competition_transmission.py (habitat x competition x transmission summary) sanity checks passed.")


# ---------------------------------------------------------------------------
# 6. Main execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    run_sanity_checks()

    print(f"\nRunning all 27 scenarios x 2 strategies (this takes a while)...")
    results = run_all_scenarios(n_replicates=30)

    here = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(FIG_DIR, "other_explorations_habitat_competition_transmission.png")
    print("\nGenerating summary figure...")
    plot_wiki_summary_figure(results, save_path=save_path)
