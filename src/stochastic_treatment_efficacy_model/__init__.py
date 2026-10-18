"""Stochastic, individual-count model comparing a *disarming* (RNAi) and a
*defensive* (protein supplementation) treatment strategy against crayfish
plague (Aphanomyces astaci) in co-occurring invasive and native crayfish.

Quick use::

    from stochastic_treatment_efficacy_model import DEFAULT_PARAMS, simulate_combined_treatment

    params = {**DEFAULT_PARAMS, "mix_fraction_kill": 1.0}   # pure RNAi strategy
    trajectory = simulate_combined_treatment(params, seed=0)  # one row per month
"""

from .engine import (
    DEFAULT_PARAMS,
    run_mix_sweep,
    simulate_combined_treatment,
    validate_params,
)

__all__ = [
    "DEFAULT_PARAMS",
    "run_mix_sweep",
    "simulate_combined_treatment",
    "validate_params",
]

__version__ = "1.0.0"
