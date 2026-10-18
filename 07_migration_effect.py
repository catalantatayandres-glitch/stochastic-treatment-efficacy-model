"""
07_migration_effect.py: Migration effect, under the same validated scenario as main results
================================================================================

Purpose
-------
Final "Other Explorations" panel: does migration change the outcome
under the exact same conditions used throughout the main results and
validation (Neutral habitat, moderate competition, beta=0.9)? Shows
native population over 10 years, migration on vs. off, for the RNAi
(disarming (RNAi)) strategy.

Fixed scenario (identical to main results / validation)
--------------------------------------------------------------------------
    Habitat      : Neutral   (K_I=6400, K_N=5500)
    Competition  : Moderate  (alpha_IN=0.3, alpha_NI=0.5)
    Transmission : beta = 0.9
    Strategy     : RNAi (disarming), coverage=50%, efficacy=100%
                   (same convention as main results)
    N0=0.2xK_N, I0=0.2xK_I

Migration conditions compared
--------------------------------------------------------------------------
    "No migration"   : migration_rate = 0.0 (the model engine (01_model_engine.py) default, off)
    "With migration" : migration_rate = 0.05 (the model engine (01_model_engine.py)'s illustrative
                        literature-informed default -- see the model engine (01_model_engine.py)
                        module docstring, Bubb et al. 2004/2006)

Output
------
One figure: native population over 10 years, 2 curves (No migration /
With migration), same scenario as main results.
"""

from __future__ import annotations

import os
import importlib.util
from typing import Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
    "migration_rate": None,  # PLACEHOLDER: set per condition
    "efficacy_I": 1.0,
    "efficacy_N_protect": 1.0,
    "overall_encounter_rate": 0.5,
    "mix_fraction_kill": 1.0,  # RNAi (disarming) -- the strategy used in main results
    "treatment_months": list(range(1, 13)),
    "n_years": 10,
    "steps_per_year": 12,
    "seed": 42,
}

MIGRATION_CONDITIONS: Dict[str, Dict[str, Any]] = {
    "No migration":   {"migration_rate": 0.0,  "color": "#1C7293"},
    "With migration": {"migration_rate": 0.05, "color": "#C1292E"},
}
MIGRATION_ORDER = list(MIGRATION_CONDITIONS.keys())


def build_params(mig_name: str) -> Dict[str, Any]:
    p = BASE_PARAMS.copy()
    p["migration_rate"] = MIGRATION_CONDITIONS[mig_name]["migration_rate"]
    return p


def run_trajectories(n_replicates: int = 50, base_seed: int = 300000) -> pd.DataFrame:
    frames = []
    seed_counter = base_seed
    for mig_name in MIGRATION_ORDER:
        p = build_params(mig_name)
        for r in range(n_replicates):
            res = stage6.simulate_combined_treatment(p, seed=seed_counter)
            seed_counter += 1
            sub = res[["month", "year", "N"]].copy()
            sub["migration_condition"] = mig_name
            sub["replicate"] = r
            frames.append(sub)
        print(f"  ...completed '{mig_name}' ({n_replicates} replicates)")
    return pd.concat(frames, ignore_index=True)


def plot_figure(trajectories: pd.DataFrame, save_path: str = None) -> None:
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 13, "axes.linewidth": 1.1})
    fig, ax = plt.subplots(figsize=(9.5, 6.5))

    for mig_name in MIGRATION_ORDER:
        color = MIGRATION_CONDITIONS[mig_name]["color"]
        sub = trajectories[trajectories["migration_condition"] == mig_name]
        summary = sub.groupby("year").agg(
            p10=("N", lambda x: np.percentile(x, 10)),
            median=("N", "median"),
            p90=("N", lambda x: np.percentile(x, 90)),
        ).reset_index()
        ax.fill_between(summary["year"], summary["p10"], summary["p90"], color=color, alpha=0.15, linewidth=0)
        ax.plot(summary["year"], summary["median"], color=color, linewidth=2.6, label=mig_name)

    ax.set_xlabel("Year", fontsize=14, labelpad=8)
    ax.set_ylabel("Native population", fontsize=14, labelpad=8)
    ax.set_xlim(0, 10)
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, linewidth=0.7)
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=12, loc="best", frameon=True, framealpha=0.95, edgecolor="#DDDDDD")

    fig.suptitle("Does migration change the outcome?", fontsize=16, fontweight="bold", color="#1F3864", y=0.99)
    fig.text(0.5, 0.905,
        "RNAi strategy, Neutral habitat, moderate competition, \u03b2 = 0.9 (same as main results)",
        ha="center", va="top", fontsize=10.5, color="#444444")
    fig.subplots_adjust(top=0.82, bottom=0.11, left=0.11, right=0.97)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  ...saved {save_path}")
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()


def run_sanity_checks() -> None:
    import matplotlib
    original_backend = matplotlib.get_backend()
    matplotlib.use("Agg")
    try:
        small_traj = run_trajectories(n_replicates=3, base_seed=1)
        expected_rows_per_condition = (BASE_PARAMS["n_years"] * BASE_PARAMS["steps_per_year"] + 1)
        expected_rows = len(MIGRATION_ORDER) * 3 * expected_rows_per_condition
        assert len(small_traj) == expected_rows, "Unexpected row count."
        assert (small_traj["N"] >= 0).all(), "Negative native populations found."
        plot_figure(small_traj)
        plt.close("all")
    finally:
        matplotlib.use(original_backend)
    print("All 07_migration_effect.py sanity checks passed.")


if __name__ == "__main__":
    run_sanity_checks()
    print("\nRunning migration on/off trajectories at the validated scenario...")
    trajectories = run_trajectories(n_replicates=50)
    here = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(FIG_DIR, "other_explorations_migration_effect.png")
    print("\nGenerating figure...")
    plot_figure(trajectories, save_path=save_path)
