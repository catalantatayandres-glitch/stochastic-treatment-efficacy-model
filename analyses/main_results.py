"""
main_results.py: Wiki figure -- native population by treatment strategy, main results
================================================================================

Purpose
-------
Produces ONE standalone, publication-quality figure for direct embedding
in the iGEM wiki -- the main results figure for this project. Shows
native crayfish population over the full 10-year simulation, with 6
curves: a "No treatment" baseline, plus 5 mix_fraction_kill strategies
(0, 0.25, 0.5, 0.75, 1), giving a finer view of the trend across the
full protect-to-kill range than just the 3 endpoints/midpoint.

Fixed scenario
--------------------------------------------------------------------------
    Habitat      : Neutral   (K_I=6400, K_N=5500 -- see the habitat definitions)
    Competition  : Moderate  (alpha_IN=0.3, alpha_NI=0.5)
    Transmission : beta = 0.9
    N0           : 0.2 x K_N (= 1,100)
    I0           : 0.2 x K_I (= 1,280)

Competition (Moderate) and transmission (beta=0.9) are the SAME values
selected and justified in the validation section (validation_bubb2006.py) -- chosen
there because they reproduce the shape of a real, documented native
displacement (Bubb et al. 2006) while staying consistent with testing
under a Neutral, non-extreme habitat. This script uses those same
validated parameters for the main results, rather than an independently
chosen value.

All other parameters are the literature-grounded values established for
this project (p_C0=0.2, mu_A=0.8, m_I=0.08, m_N=0.01, r_I=0.45, r_N=0.12,
gamma_N=0.0, invasive_carrier_fitness_cost=0.0, prob_recruit_carrier=0.2)
-- see the project wiki Model page for the full basis of each.

The 6 curves (identical convention used throughout this project):
efficacy fixed at 100% and coverage fixed at 50% whenever treatment is
active; mix_fraction_kill is the only thing that changes between the 5
active curves.
    No treatment          : overall_encounter_rate = 0.0
    Defensive (protein) (mix=0)   : coverage=0.5, mix_fraction_kill=0.0
    Mixed (mix=0.25)       : coverage=0.5, mix_fraction_kill=0.25
    Mixed (mix=0.5)        : coverage=0.5, mix_fraction_kill=0.5
    Mixed (mix=0.75)       : coverage=0.5, mix_fraction_kill=0.75
    Disarming (RNAi) (mix=1)      : coverage=0.5, mix_fraction_kill=1.0

Output
------
One publication-styled PNG file (300 DPI, ready to embed directly in
the wiki), saved in the SAME FOLDER as this script:
    wiki_figure_transmission_0.9.png
x-axis = year (0-10), y-axis = median native population, 6 lines
with 10th-90th percentile bands, large/clean fonts, minimal-clutter style.
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
# 2. Fixed ecological scenario (change these constants if needed)
# ---------------------------------------------------------------------------

K_I_FIXED = 6400.0    # Neutral habitat
K_N_FIXED = 5500.0    # Neutral habitat
ALPHA_IN_FIXED = 0.3   # Moderate competition
ALPHA_NI_FIXED = 0.5   # Moderate competition

N0_FIXED = 0.2 * K_N_FIXED  # = 1,100
I0_FIXED = 0.2 * K_I_FIXED  # = 1,280

# --- The 2 transmission values compared (one figure each) ---
TRANSMISSION_VALUES: List[float] = [0.9]  # single value decided for the wiki, matches the validation section

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

    "alpha_IN": ALPHA_IN_FIXED,
    "alpha_NI": ALPHA_NI_FIXED,

    "beta": None,  # PLACEHOLDER: set per figure from TRANSMISSION_VALUES
    "transmission_mode": "frequency",
    "epsilon": 1e-6,
    "mu_A": 0.8,
    "gamma_N": 0.0,
    "invasive_carrier_fitness_cost": 0.0,
    "prob_recruit_carrier": 0.2,
    "migration_rate": 0.0,

    "efficacy_I": 1.0,
    "efficacy_N_protect": 1.0,
    "overall_encounter_rate": None,   # PLACEHOLDER: set per treatment curve
    "mix_fraction_kill": None,
    "treatment_months": list(range(1, 13)),

    "n_years": 10,
    "steps_per_year": 12,
    "seed": 42,
}

# --- The 6 curves shown: "No treatment" baseline + 5 mix_fraction_kill
# strategies (0, 0.25, 0.5, 0.75, 1), a finer view than the original 3
# (0, 0.5, 1) so the trend across the full mix range is visible. Colors
# form a clean gradient from protect (blue) to kill (red). ---
TREATMENT_CURVES: List[Dict[str, Any]] = [
    {"label": "No treatment",   "color": "#888888", "linestyle": "--",
     "overall_encounter_rate": 0.0, "mix_fraction_kill": 0.0},
    {"label": "Protein Supplementation (mix=0)",    "color": "#1C7293", "linestyle": "-",
     "overall_encounter_rate": 0.5, "mix_fraction_kill": 0.0},
    {"label": "Mixed (mix=0.25)",  "color": "#4FA37D", "linestyle": "-",
     "overall_encounter_rate": 0.5, "mix_fraction_kill": 0.25},
    {"label": "Mixed (mix=0.5)",   "color": "#F9A03F", "linestyle": "-",
     "overall_encounter_rate": 0.5, "mix_fraction_kill": 0.5},
    {"label": "Mixed (mix=0.75)",  "color": "#E0631E", "linestyle": "-",
     "overall_encounter_rate": 0.5, "mix_fraction_kill": 0.75},
    {"label": "RNAi strategy (mix=1)",       "color": "#C1292E", "linestyle": "-",
     "overall_encounter_rate": 0.5, "mix_fraction_kill": 1.0},
]
CURVE_LABELS = [c["label"] for c in TREATMENT_CURVES]


def build_params(beta_value: float, curve: Dict[str, Any]) -> Dict[str, Any]:
    """Assemble the full parameter dict for one (transmission value, treatment curve) cell."""
    p = BASE_PARAMS.copy()
    p["beta"] = beta_value
    p["overall_encounter_rate"] = curve["overall_encounter_rate"]
    p["mix_fraction_kill"] = curve["mix_fraction_kill"]
    return p


# ---------------------------------------------------------------------------
# 3. Monte Carlo runner
# ---------------------------------------------------------------------------

def run_trajectories_for_beta(
    beta_value: float,
    n_replicates: int = 50,
    base_seed: int = 200000,
) -> pd.DataFrame:
    """
    Run full monthly native-population trajectories for all 4 treatment
    curves, at the given fixed transmission (beta) value.

    Returns
    -------
    pd.DataFrame
        One row per (curve_label, replicate, month), with columns
        'curve_label', 'replicate', 'month', 'year', 'N'.
    """
    frames = []
    seed_counter = base_seed + int(beta_value * 100000)  # keep the two figures' seed streams independent

    for curve in TREATMENT_CURVES:
        p = build_params(beta_value, curve)
        for r in range(n_replicates):
            res = engine.simulate_combined_treatment(p, seed=seed_counter)
            seed_counter += 1
            sub = res[["month", "year", "N"]].copy()
            sub["curve_label"] = curve["label"]
            sub["replicate"] = r
            frames.append(sub)
        print(f"  ...completed '{curve['label']}' (beta={beta_value}, {n_replicates} replicates)")

    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# 4. Publication-quality figure
# ---------------------------------------------------------------------------

def plot_wiki_figure(trajectories: pd.DataFrame, beta_value: float, save_path: str) -> None:
    """
    Single, publication-styled figure: native population over 10 years,
    4 treatment-strategy lines with 10th-90th percentile bands. Saved as
    a 300 DPI PNG ready to embed directly in the wiki.
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 13,
        "axes.linewidth": 1.1,
    })

    fig, ax = plt.subplots(figsize=(10, 6.5))  # normal screen DPI for the interactive window; 300 DPI is applied only at save time below

    curve_lookup = {c["label"]: c for c in TREATMENT_CURVES}
    for label in CURVE_LABELS:
        c = curve_lookup[label]
        sub = trajectories[trajectories["curve_label"] == label]
        summary = sub.groupby("year").agg(
            N_p10=("N", lambda x: np.percentile(x, 10)),
            N_median=("N", "median"),
            N_p90=("N", lambda x: np.percentile(x, 90)),
        ).reset_index()

        ax.fill_between(summary["year"], summary["N_p10"], summary["N_p90"], color=c["color"], alpha=0.15, linewidth=0)
        ax.plot(summary["year"], summary["N_median"], color=c["color"], linewidth=2.6,
                 linestyle=c["linestyle"], label=label)

    ax.set_xlabel("Year", fontsize=14, labelpad=8)
    ax.set_ylabel("Native population", fontsize=14, labelpad=8)
    ax.set_xlim(0, 10)
    ax.set_ylim(bottom=0)
    ax.set_title("Impact of different treatment strategies",
                 fontsize=15.5, fontweight="bold", pad=14)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, linewidth=0.7)
    ax.tick_params(labelsize=12)

    legend = ax.legend(fontsize=12, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                        frameon=True, framealpha=0.95, edgecolor="#DDDDDD", borderaxespad=0)
    legend.get_frame().set_linewidth(0.8)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"  ...saved {save_path}")
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()


# ---------------------------------------------------------------------------
# 5. Sanity checks
# ---------------------------------------------------------------------------

def run_sanity_checks() -> None:
    """Run a small, fast version of the full pipeline to confirm it works end-to-end."""
    import matplotlib
    original_backend = matplotlib.get_backend()
    matplotlib.use("Agg")
    try:
        small_traj = run_trajectories_for_beta(0.5, n_replicates=3, base_seed=1)
        expected_rows_per_condition = (BASE_PARAMS["n_years"] * BASE_PARAMS["steps_per_year"] + 1)
        expected_rows = len(TREATMENT_CURVES) * 3 * expected_rows_per_condition
        assert len(small_traj) == expected_rows, "Unexpected row count."
        assert (small_traj["N"] >= 0).all(), "Negative native populations found."
        assert set(small_traj["curve_label"].unique()) == set(CURVE_LABELS)

        import tempfile
        test_path = os.path.join(tempfile.gettempdir(), "_stage20_sanity_check.png")
        plot_wiki_figure(small_traj, 0.5, test_path)
        assert os.path.exists(test_path), "Figure was not saved to disk."
        os.remove(test_path)

        plt.close("all")
    finally:
        matplotlib.use(original_backend)

    print("All main_results.py (wiki transmission-scenario figures) sanity checks passed.")


# ---------------------------------------------------------------------------
# 6. Main execution
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    run_sanity_checks()

    print(f"\nFixed scenario: Neutral habitat (K_I={K_I_FIXED:.0f}, K_N={K_N_FIXED:.0f}), "
          f"moderate competition, N0={N0_FIXED:.0f}, I0={I0_FIXED:.0f}")

    for beta_value in TRANSMISSION_VALUES:
        print(f"\nRunning trajectories for beta={beta_value} (this takes a moment)...")
        trajectories = run_trajectories_for_beta(beta_value, n_replicates=50)

        save_path = os.path.join(FIG_DIR, f"wiki_figure_transmission_{beta_value}.png")
        print(f"Generating publication figure for beta={beta_value}...")
        plot_wiki_figure(trajectories, beta_value, save_path)
