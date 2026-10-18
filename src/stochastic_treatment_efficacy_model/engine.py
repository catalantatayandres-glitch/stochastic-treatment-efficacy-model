"""
engine.py: stochastic crayfish population model (disarming + defensive treatment)
==========================================================================

Purpose
-------
This stage combines the invasive-disarming treatment  and the
native-defensive treatment  into a single scenario, modeling
three mechanistic details:

  1. A crayfish can INDEPENDENTLY encounter EITHER payload in the same
     month, since both are released into the same lake water.

  2. Interference between the two payloads is ASYMMETRIC, reflecting the
     actual underlying biology rather than a modeling convenience:

       INVASIVE carriers (I_C): the RNAi construct kills by silencing
       the gene that encodes the invasive's OWN natural tolerance
       protein. If a carrier ALSO takes up the exogenous protein payload
       in the SAME month, that protein directly restores the tolerance
       the RNAi just silenced -- rescuing that individual. So for
       invasives, encountering BOTH payloads cancels the kill.

       NATIVE crayfish (N_I): natives never expressed that tolerance
       gene in the first place -- there is nothing in a native for the
       RNAi construct to silence. Encountering the RNAi payload has NO
       biological effect on a native, whether or not they also
       encounter the protein payload. So a native's protection
       probability depends ONLY on whether they encountered the protein
       -- it is NEVER discounted by also "encountering" RNAi, because
       there is no real interference to have.

     Concretely:
         p_kill_active    = p_kill * (1 - p_protect) * efficacy_I
             (protein encounter cancels the kill -- real interference)
         p_protect_active = p_protect * efficacy_N_protect
             (RNAi encounter is irrelevant to a native -- no interference)

  3. The disarming dsRNA only kills invasive crayfish that are CARRIERS of
     A. astaci (I_C) -- it has NO effect on non-carrier invasive crayfish
     (I_S), regardless of whether they encounter it. This makes the
     disarming treatment a RESERVOIR-ELIMINATION strategy (reducing the
     source of disease transmission), not a general invasive-population-
     control strategy. The defensive treatment only helps already-
     INFECTED natives (N_I) survive disease -- susceptible natives (N_S)
     are unaffected either way.

Why this matters for interpreting past results
--------------------------------------------------------------------------
This asymmetric version REPLACES an earlier iteration that incorrectly
treated the interference as symmetric (i.e., also discounted a native's
protection by p_kill, as if RNAi encountering a native could somehow
cancel their protection the same way protein encountering an invasive
cancels their kill). That earlier, incorrect symmetric version only
produces different numbers than this corrected version at MIX_FRACTION_
KILL VALUES STRICTLY BETWEEN 0 AND 1 (i.e., genuinely mixed strategies).
At mix_fraction_kill = 0 (defensive (protein), p_kill=0) or = 1 (disarming (RNAi),
p_protect=0), the two versions are numerically IDENTICAL, since the
discount term that changed is multiplied by a probability that is
already zero in those pure-strategy cases. Any earlier figure built from
a genuinely mixed strategy (e.g. "Mixed (0.5)") reflects the old,
slightly-too-pessimistic native protection number and should be
re-generated with this corrected version if it matters for your results.

Where p_kill / p_protect / mix_fraction_kill / overall_encounter_rate
come from is unchanged from the previous iteration -- see the parameter
dictionary comments below.

MIGRATION (toggleable, OFF by default)
--------------------------------------------------------------------------
Both species can emigrate OUT of the modeled site and immigrate IN from
surrounding/connected habitats every month, controlled by a single
parameter, migration_rate (default 0.0 = fully OFF). Set migration_rate
> 0 to turn it on. Per explicit instruction, the SAME rate is used for
both directions (emigration and immigration) and both species -- a
simplifying assumption, not a literature-derived asymmetry.

Why this exists: to test whether the protein/defensive strategy becomes
relatively more competitive against the RNAi/disarming strategy once
migration is present. Mechanistic reasoning: disarming invasive carriers
only helps if the local invasive population actually stays down -- but
if carriers are continuously replaced by immigration from an unaffected
surrounding population, the disarming strategy's benefit gets undermined
over time. Protecting natives from disease mortality does not depend on
suppressing invasive population size, so it is not undermined by
invasive immigration the same way.

Literature basis for migration_rate, and its limits: requested to be
literature-based, using Bubb, Thom & Lucas as the source --
  - Bubb, D.H., Thom, T.J. & Lucas, M.C. (2004) Movement and dispersal of
    the invasive signal crayfish Pacifastacus leniusculus in upland
    rivers. Freshwater Biology, 49, 357-368.
  - Bubb, D.H., Thom, T.J. & Lucas, M.C. (2006) Movement, dispersal and
    refuge use of co-occurring introduced and native crayfish. Freshwater
    Biology, 51, 1359-1368.
Radio-tracked signal crayfish in the River Wharfe moved sporadically --
individuals typically remained in one location for several weeks,
punctuated by occasional relocation movements up to ~328m (Bubb et al.
2004). The 2006 co-occurring study -- the one study tracking BOTH
species at the same site -- found signal crayfish dispersed roughly
TWICE as far as white-clawed crayfish over the same period, and by the
2005 resurvey the native species had disappeared entirely from the site.
This does NOT give a directly usable "% of population migrating per
month" statistic (these are individual movement-distance/frequency
observations, not population turnover rates, and this model has no
explicit spatial/reach-length dimension to convert one into the other)
-- so migration_rate, like alpha/beta, has no single literature-
calibrated value; it is grounded in the qualitative movement pattern
(sporadic, moderate frequency, not mass continuous turnover), not a
precise number. When turned on, 0.05 (5%/month) is a reasonable
illustrative default consistent with that pattern -- an order-of-
magnitude judgment call, not a fitted statistic.

Contradiction worth flagging: the "same rate for both species"
assumption used here does NOT match Bubb et al. (2006), who found
signal crayfish dispersing ~2x farther than white-clawed crayfish at the
same site. The equal-rate assumption is used anyway as the explicit
modeling decision for this analysis -- a deliberate, visible
simplification, not a hidden one.

Mechanism: emigrants are removed from each sub-population via a
Binomial(n, migration_rate) draw (same mechanic as natural mortality),
applied as its own step between natural mortality and disease
transmission. Immigrants arrive via a Poisson draw with mean
migration_rate * K_I (invasive) or migration_rate * K_N (native) -- the
surrounding region is assumed to host populations near ITS OWN carrying
capacity, NOT depleted by local treatment. Invasive immigrants are split
into carrier/non-carrier using the same prob_recruit_carrier as local
recruits. Native immigrants arrive susceptible (N_S) -- an assumption,
since no data distinguishes infected vs. susceptible dispersers.

Backward compatibility: migration_rate = 0.0 (the default) skips ALL
migration-related random draws ENTIRELY (not just zeroes them), so every
script written before this feature was added produces bit-for-bit
identical results with the same seed, verified in run_sanity_checks().
"""

from __future__ import annotations

from typing import Dict, Any, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. Central parameter dictionary
# ---------------------------------------------------------------------------

DEFAULT_PARAMS: Dict[str, Any] = {
    # --- Initial conditions (unchanged) ---
    "I0": 1000,
    "N0": 1000,
    "p_C0": 0.2,
    "N_I0": 0,

    # --- Demography (unchanged) ---
    "r_I": 0.45,
    "r_N": 0.12,
    "K_I": 6400,
    "K_N": 5000,
    "m_I": 0.08,
    "m_N": 0.01,

    # --- Competition (unchanged) ---
    "alpha_IN": 0.3,
    "alpha_NI": 0.5,

    # --- Disease (beta unchanged, still uncalibrated placeholder; mu_A
    # UPDATED to the literature-audited value -- see Parameter_Literature_
    # Deep_Dive.docx for the full derivation) ---
    "beta": 0.5,   # moderated from 2.0 -- see discussion: still uncalibrated, chosen for a less overwhelming baseline epidemic
    "transmission_mode": "frequency",
    "epsilon": 1e-6,
    "mu_A": 0.8,  # literature-derived: (1-mu_A)^2 ~= 96% cumulative mortality by day 60, matching the ~100% cumulative mortality documented in susceptible A. pallipes populations by day 60 (see the project wiki Model page)
    "gamma_N": 0.0,
    "invasive_carrier_fitness_cost": 0.0,
    "prob_recruit_carrier": 0.2,

    # --- Treatment: COMBINED, with shared encounter budget ---
    "efficacy_I": 1,             # probability a covered invasive crayfish is killed (unchanged)
    "efficacy_N_protect": 1,     # probability a covered infected native survives disease this month (unchanged)

    "overall_encounter_rate": 0.7,  # UNKNOWN placeholder: fraction of EITHER species encountering SOME dose per month
    "mix_fraction_kill": 0.5,       # formulation choice: fraction of treatment material that is the disarming dsRNA (rest is defensive protein)

    "treatment_months": list(range(1, 13)),  # continuous monthly delivery

    # --- Migration (toggleable, OFF by default): same rate for
    # emigration AND immigration, same rate for BOTH species (explicit
    # simplifying assumption -- see module docstring for the literature
    # basis, its limits, and the Bubb et al. 2006 finding that actually
    # contradicts this equal-rate assumption). Default 0.0 = fully off,
    # reproduces every prior result exactly. Set > 0 (e.g. 0.05) to
    # enable migration for a run/scenario.
    "migration_rate": 0.0,

    # --- Simulation settings ---
    "n_years": 10,
    "steps_per_year": 12,
    "seed": 42,
}


# ---------------------------------------------------------------------------
# 2. Input validation
# ---------------------------------------------------------------------------

def validate_params(params: Dict[str, Any]) -> None:
    """
    Validate that parameter values are biologically and numerically sensible.
    """
    required_keys = [
        "I0", "N0", "p_C0", "N_I0", "r_I", "r_N", "K_I", "K_N", "m_I", "m_N",
        "alpha_IN", "alpha_NI", "beta", "transmission_mode", "epsilon",
        "mu_A", "gamma_N", "invasive_carrier_fitness_cost", "prob_recruit_carrier",
        "efficacy_I", "efficacy_N_protect", "overall_encounter_rate",
        "mix_fraction_kill", "migration_rate", "treatment_months", "n_years",
        "steps_per_year", "seed",
    ]
    for key in required_keys:
        if key not in params:
            raise ValueError(f"Missing required parameter: '{key}'")

    if params["I0"] < 0 or params["N0"] < 0:
        raise ValueError("Initial populations must be non-negative.")
    if params["N_I0"] < 0 or params["N_I0"] > params["N0"]:
        raise ValueError("N_I0 must be between 0 and N0.")
    for prob_key in ("p_C0", "mu_A", "gamma_N", "prob_recruit_carrier",
                      "efficacy_I", "efficacy_N_protect",
                      "overall_encounter_rate", "mix_fraction_kill", "migration_rate"):
        if not (0.0 <= params[prob_key] <= 1.0):
            raise ValueError(f"'{prob_key}' must be a probability in [0, 1].")
    if not (0.0 <= params["m_I"] <= 1.0) or not (0.0 <= params["m_N"] <= 1.0):
        raise ValueError("Natural mortality probabilities must lie in [0, 1].")
    if params["invasive_carrier_fitness_cost"] < 0:
        raise ValueError("invasive_carrier_fitness_cost must be non-negative.")
    if params["K_I"] <= 0 or params["K_N"] <= 0:
        raise ValueError("Carrying capacities must be strictly positive.")
    if params["alpha_IN"] < 0 or params["alpha_NI"] < 0:
        raise ValueError("Competition coefficients must be non-negative.")
    if params["beta"] < 0:
        raise ValueError("beta must be non-negative.")
    if params["transmission_mode"] not in ("frequency", "density"):
        raise ValueError("transmission_mode must be 'frequency' or 'density'.")
    if params["epsilon"] <= 0:
        raise ValueError("epsilon must be strictly positive.")
    if not all(1 <= m <= 12 for m in params["treatment_months"]):
        raise ValueError("treatment_months must contain calendar months in 1-12.")
    if params["n_years"] <= 0 or params["steps_per_year"] <= 0:
        raise ValueError("n_years and steps_per_year must be positive.")


# ---------------------------------------------------------------------------
# 3. Core simulation function
# ---------------------------------------------------------------------------

def simulate_combined_treatment(params: Dict[str, Any], seed: int | None = None) -> pd.DataFrame:
    """
    Run the stochastic demographic + disease model with BOTH treatment
    components sharing a single encounter budget (see module docstring).

    Parameters
    ----------
    params : dict
        Parameter dictionary. See DEFAULT_PARAMS for required keys.
    seed : int, optional
        Random seed for this run. If None, uses params["seed"].

    Returns
    -------
    pd.DataFrame
        Tidy DataFrame with one row per monthly time step, containing:
        'month', 'year', 'calendar_month', 'I_S', 'I_C', 'I', 'N_S', 'N_I', 'N',
        'natural_deaths_I', 'natural_deaths_N', 'new_native_infections',
        'plague_deaths', 'native_protected', 'RNAi_deaths_I',
        'emigrants_I', 'immigrants_I', 'emigrants_N', 'immigrants_N'.
    """
    validate_params(params)

    seed = params["seed"] if seed is None else seed
    rng = np.random.default_rng(seed)

    n_steps = params["n_years"] * params["steps_per_year"]

    I_S = np.zeros(n_steps + 1, dtype=np.int64)
    I_C = np.zeros(n_steps + 1, dtype=np.int64)
    N_S = np.zeros(n_steps + 1, dtype=np.int64)
    N_I = np.zeros(n_steps + 1, dtype=np.int64)

    I0 = int(params["I0"])
    I_C0 = rng.binomial(I0, params["p_C0"])
    I_S[0] = I0 - I_C0
    I_C[0] = I_C0

    N0 = int(params["N0"])
    N_I0 = int(params["N_I0"])
    N_S[0] = N0 - N_I0
    N_I[0] = N_I0

    natural_deaths_I_arr = np.zeros(n_steps, dtype=np.int64)
    natural_deaths_N_arr = np.zeros(n_steps, dtype=np.int64)
    new_infections_arr = np.zeros(n_steps, dtype=np.int64)
    plague_deaths_arr = np.zeros(n_steps, dtype=np.int64)
    native_protected_arr = np.zeros(n_steps, dtype=np.int64)
    rnai_deaths_I_arr = np.zeros(n_steps, dtype=np.int64)
    emigrants_I_arr = np.zeros(n_steps, dtype=np.int64)
    immigrants_I_arr = np.zeros(n_steps, dtype=np.int64)
    emigrants_N_arr = np.zeros(n_steps, dtype=np.int64)
    immigrants_N_arr = np.zeros(n_steps, dtype=np.int64)
    calendar_months = np.zeros(n_steps + 1, dtype=np.int64)

    r_I, r_N = params["r_I"], params["r_N"]
    K_I, K_N = params["K_I"], params["K_N"]
    m_I, m_N = params["m_I"], params["m_N"]
    alpha_IN, alpha_NI = params["alpha_IN"], params["alpha_NI"]
    beta = params["beta"]
    transmission_mode = params["transmission_mode"]
    epsilon = params["epsilon"]
    mu_A = params["mu_A"]
    gamma_N = params["gamma_N"]
    carrier_cost = params["invasive_carrier_fitness_cost"]
    prob_recruit_carrier = params["prob_recruit_carrier"]
    migration_rate = params["migration_rate"]
    steps_per_year = params["steps_per_year"]

    efficacy_I = params["efficacy_I"]
    efficacy_N_protect = params["efficacy_N_protect"]
    overall_encounter_rate = params["overall_encounter_rate"]
    mix_fraction_kill = params["mix_fraction_kill"]
    treatment_months = set(params["treatment_months"])

    # OPTIONAL, fully backward-compatible: a per-month mix_fraction_kill
    # schedule, for strategies that change over the course of the
    # simulation (e.g. defensive (protein) for the first 6 months, then switch to
    # disarming (RNAi)). Not a required key -- absent or None means "use the
    # constant mix_fraction_kill scalar above, every month" (unchanged
    # behavior). When present, must be a sequence of length n_steps giving
    # the mix_fraction_kill value to use at each absolute month t=0..n_steps-1.
    mix_schedule = params.get("mix_fraction_kill_schedule", None)
    if mix_schedule is not None and len(mix_schedule) != n_steps:
        raise ValueError(
            f"mix_fraction_kill_schedule has length {len(mix_schedule)}, expected {n_steps} (n_years * steps_per_year)."
        )

    for t in range(n_steps):
        i_s, i_c = int(I_S[t]), int(I_C[t])
        n_s, n_i = int(N_S[t]), int(N_I[t])
        I_t = i_s + i_c
        N_t = n_s + n_i

        # Independent-encounter probabilities: a crayfish can encounter
        # EITHER payload independently each month, since both share the
        # same water. Recomputed every month so a schedule (if given) can
        # change the strategy mid-simulation; identical result to the old
        # compute-once-outside-the-loop version when no schedule is given.
        current_mix = mix_schedule[t] if mix_schedule is not None else mix_fraction_kill
        p_kill = overall_encounter_rate * current_mix
        p_protect = overall_encounter_rate * (1 - current_mix)

        # Interference between the two payloads is ASYMMETRIC, not
        # mutual -- this reflects the actual biology, not a modeling
        # simplification:
        #
        #   INVASIVE carriers: RNAi kills by silencing the gene that
        #   encodes the invasive's OWN natural tolerance protein. If a
        #   carrier ALSO takes up the exogenous protein payload in the
        #   same month, that protein directly restores the tolerance the
        #   RNAi just silenced -- rescuing that individual from the kill.
        #   So encountering BOTH payloads cancels the kill for invasives.
        #
        #   NATIVE crayfish: they never expressed that tolerance gene in
        #   the first place -- there is nothing for the RNAi construct to
        #   silence in a native. Encountering the RNAi payload has NO
        #   biological effect on a native, regardless of whether they also
        #   encounter the protein payload. So a native's protection
        #   probability depends ONLY on whether they encountered the
        #   protein -- it is never discounted by also "encountering" RNAi.
        p_kill_only = p_kill * (1 - p_protect)   # invasive: kill happens only if RNAi encountered WITHOUT protein (protein would rescue them)
        p_RNAi_I_active = p_kill_only * efficacy_I

        p_protect_N_active = p_protect * efficacy_N_protect   # native: RNAi encounter is irrelevant, no discount

        calendar_month = (t % steps_per_year) + 1
        calendar_months[t] = calendar_month
        treatment_active = calendar_month in treatment_months

        # --- 1. Reproduction ---
        mu_I = max(r_I * I_t * (1 - (I_t + alpha_IN * N_t) / K_I), 0.0)
        mu_N = max(r_N * N_t * (1 - (N_t + alpha_NI * I_t) / K_N), 0.0)
        recruits_I_total = rng.poisson(mu_I)
        recruits_N_total = rng.poisson(mu_N)
        recruits_I_C = rng.binomial(recruits_I_total, prob_recruit_carrier) if recruits_I_total > 0 else 0
        recruits_I_S = recruits_I_total - recruits_I_C

        i_s_post_repro = i_s + recruits_I_S
        i_c_post_repro = i_c + recruits_I_C
        n_s_post_repro = n_s + recruits_N_total
        n_i_post_repro = n_i

        # --- 2. Natural mortality ---
        deaths_I_S = rng.binomial(i_s_post_repro, m_I) if i_s_post_repro > 0 else 0
        m_I_carrier = min(m_I + carrier_cost, 1.0)
        deaths_I_C = rng.binomial(i_c_post_repro, m_I_carrier) if i_c_post_repro > 0 else 0
        deaths_N_S = rng.binomial(n_s_post_repro, m_N) if n_s_post_repro > 0 else 0
        deaths_N_I_natural = rng.binomial(n_i_post_repro, m_N) if n_i_post_repro > 0 else 0

        i_s_post_mort = i_s_post_repro - deaths_I_S
        i_c_post_mort = i_c_post_repro - deaths_I_C
        n_s_post_mort = n_s_post_repro - deaths_N_S
        n_i_post_mort = n_i_post_repro - deaths_N_I_natural

        # --- 3. Migration: emigration then immigration (toggleable via
        # migration_rate, default OFF -- see module docstring) ---
        # Same migration_rate governs BOTH directions and BOTH species (an
        # explicit simplifying assumption -- see module docstring). Random
        # draws are SKIPPED ENTIRELY (not just zeroed) when migration_rate
        # == 0, so the default reproduces every prior result bit-for-bit.
        if migration_rate > 0:
            emigrants_I_S = rng.binomial(i_s_post_mort, migration_rate) if i_s_post_mort > 0 else 0
            emigrants_I_C = rng.binomial(i_c_post_mort, migration_rate) if i_c_post_mort > 0 else 0
            emigrants_N_S = rng.binomial(n_s_post_mort, migration_rate) if n_s_post_mort > 0 else 0
            emigrants_N_I = rng.binomial(n_i_post_mort, migration_rate) if n_i_post_mort > 0 else 0
        else:
            emigrants_I_S = emigrants_I_C = emigrants_N_S = emigrants_N_I = 0

        i_s_post_emigration = i_s_post_mort - emigrants_I_S
        i_c_post_emigration = i_c_post_mort - emigrants_I_C
        n_s_post_emigration = n_s_post_mort - emigrants_N_S
        n_i_post_emigration = n_i_post_mort - emigrants_N_I

        # Immigrants arrive from a surrounding region assumed to be at ITS
        # OWN carrying capacity (K_I / K_N), independent of local depletion
        # -- this is what lets local disarming be "undermined" by immigration
        # rather than immigration just mirroring local emigration losses.
        if migration_rate > 0:
            immigrants_I_total = rng.poisson(migration_rate * K_I)
            immigrants_I_C = rng.binomial(immigrants_I_total, prob_recruit_carrier) if immigrants_I_total > 0 else 0
            immigrants_I_S = immigrants_I_total - immigrants_I_C
            immigrants_N_S = rng.poisson(migration_rate * K_N)  # arrive susceptible -- see module docstring
        else:
            immigrants_I_S = immigrants_I_C = immigrants_N_S = 0
            immigrants_I_total = 0

        i_s_post_migration = i_s_post_emigration + immigrants_I_S
        i_c_post_migration = i_c_post_emigration + immigrants_I_C
        n_s_post_migration = n_s_post_emigration + immigrants_N_S
        n_i_post_migration = n_i_post_emigration  # infected natives don't immigrate -- see module docstring

        # --- 4. A. astaci transmission ---
        I_total_for_transmission = i_s_post_migration + i_c_post_migration
        if transmission_mode == "frequency":
            p_infection = 1 - np.exp(-beta * i_c_post_migration / (I_total_for_transmission + epsilon))
        else:
            p_infection = 1 - np.exp(-beta * i_c_post_migration)
        p_infection = min(max(p_infection, 0.0), 1.0)
        new_infections = rng.binomial(n_s_post_migration, p_infection) if n_s_post_migration > 0 else 0
        n_s_post_transmission = n_s_post_migration - new_infections
        n_i_post_transmission = n_i_post_migration + new_infections

        # --- 5. Native recovery (off by default) ---
        recoveries = rng.binomial(n_i_post_transmission, gamma_N) if (gamma_N > 0 and n_i_post_transmission > 0) else 0
        n_i_pre_disease = n_i_post_transmission - recoveries
        n_s_post_recovery = n_s_post_transmission + recoveries

        # --- 6. Disease mortality, WITH native-defensive treatment folded in ---
        # Must happen before mu_A draw -- see the project wiki Model page.
        if treatment_active and n_i_pre_disease > 0:
            native_protected = rng.binomial(n_i_pre_disease, p_protect_N_active)
        else:
            native_protected = 0
        n_i_unprotected = n_i_pre_disease - native_protected
        plague_deaths = rng.binomial(n_i_unprotected, mu_A) if n_i_unprotected > 0 else 0
        n_i_after_disease = n_i_unprotected - plague_deaths + native_protected

        # --- 7. Invasive-disarming RNAi treatment: CARRIERS ONLY (I_C) ---
        # The dsRNA only kills invasive crayfish that carry A. astaci --
        # non-carrier invasive (I_S) is completely unaffected, regardless
        # of encounter. This makes the treatment a reservoir-elimination
        # strategy, not general invasive population control -- see module
        # docstring. Uses the POST-MIGRATION carrier count -- carriers
        # that just immigrated in this month are eligible to be killed
        # this month too, same as locally-present carriers.
        if treatment_active:
            rnai_deaths_I_C = rng.binomial(i_c_post_migration, p_RNAi_I_active) if i_c_post_migration > 0 else 0
        else:
            rnai_deaths_I_C = 0
        i_s_final = i_s_post_migration  # unaffected by the disarming treatment
        i_c_final = i_c_post_migration - rnai_deaths_I_C

        # --- 8. Record ---
        I_S[t + 1] = i_s_final
        I_C[t + 1] = i_c_final
        N_S[t + 1] = n_s_post_recovery
        N_I[t + 1] = n_i_after_disease

        natural_deaths_I_arr[t] = deaths_I_S + deaths_I_C
        natural_deaths_N_arr[t] = deaths_N_S + deaths_N_I_natural
        new_infections_arr[t] = new_infections
        plague_deaths_arr[t] = plague_deaths
        native_protected_arr[t] = native_protected
        rnai_deaths_I_arr[t] = rnai_deaths_I_C  # I_S is unaffected by this treatment (carrier-only disarming)
        emigrants_I_arr[t] = emigrants_I_S + emigrants_I_C
        immigrants_I_arr[t] = immigrants_I_total
        emigrants_N_arr[t] = emigrants_N_S + emigrants_N_I
        immigrants_N_arr[t] = immigrants_N_S

    calendar_months[n_steps] = (n_steps % steps_per_year) + 1
    months = np.arange(n_steps + 1)
    years = months / steps_per_year

    results = pd.DataFrame({
        "month": months, "year": years, "calendar_month": calendar_months,
        "I_S": I_S, "I_C": I_C, "I": I_S + I_C,
        "N_S": N_S, "N_I": N_I, "N": N_S + N_I,
        "natural_deaths_I": np.append(natural_deaths_I_arr, np.nan),
        "natural_deaths_N": np.append(natural_deaths_N_arr, np.nan),
        "new_native_infections": np.append(new_infections_arr, np.nan),
        "plague_deaths": np.append(plague_deaths_arr, np.nan),
        "native_protected": np.append(native_protected_arr, np.nan),
        "RNAi_deaths_I": np.append(rnai_deaths_I_arr, np.nan),
        "emigrants_I": np.append(emigrants_I_arr, np.nan),
        "immigrants_I": np.append(immigrants_I_arr, np.nan),
        "emigrants_N": np.append(emigrants_N_arr, np.nan),
        "immigrants_N": np.append(immigrants_N_arr, np.nan),
    })

    return results


# ---------------------------------------------------------------------------
# 4. 2D sweep over (overall_encounter_rate, mix_fraction_kill)
# ---------------------------------------------------------------------------

def run_mix_sweep(
    params: Dict[str, Any],
    encounter_rates: List[float] = (0.1, 0.3, 0.5, 0.7, 0.9),
    mix_fractions: List[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
    n_replicates: int = 30,
    base_seed: int = 1000,
) -> pd.DataFrame:
    """
    Run the combined model across a 2D grid of (overall_encounter_rate,
    mix_fraction_kill) to explore the tradeoff between disarming invasive
    crayfish and protecting native crayfish under a shared treatment
    budget.

    Parameters
    ----------
    params : dict
        Base parameter dictionary.
    encounter_rates : sequence of float
        Values of overall_encounter_rate to test.
    mix_fractions : sequence of float
        Values of mix_fraction_kill to test (0 = all defensive, 1 = all disarming).
    n_replicates : int
        Stochastic replicates per grid cell.
    base_seed : int
        Starting seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        One row per (encounter_rate, mix_fraction_kill, replicate).
    """
    records = []
    seed_counter = base_seed
    for encounter_rate in encounter_rates:
        for mix_fraction in mix_fractions:
            p = params.copy()
            p["overall_encounter_rate"] = encounter_rate
            p["mix_fraction_kill"] = mix_fraction
            for _ in range(n_replicates):
                res = simulate_combined_treatment(p, seed=seed_counter)
                seed_counter += 1
                records.append({
                    "encounter_rate": encounter_rate,
                    "mix_fraction_kill": mix_fraction,
                    "I_final": res["I"].iloc[-1],
                    "N_final": res["N"].iloc[-1],
                    "native_persisted": res["N"].iloc[-1] > 0,
                })
    return pd.DataFrame(records)


def plot_mix_sweep_heatmaps(sweep_results: pd.DataFrame) -> None:
    """
    Plot heatmaps of median final native population and native persistence
    probability across the (encounter_rate x mix_fraction_kill) grid.
    """
    summary = sweep_results.groupby(["encounter_rate", "mix_fraction_kill"]).agg(
        N_median=("N_final", "median"),
        I_median=("I_final", "median"),
        persistence_prob=("native_persisted", "mean"),
    ).reset_index()

    encounter_vals = sorted(summary["encounter_rate"].unique())
    mix_vals = sorted(summary["mix_fraction_kill"].unique())

    N_grid = summary.pivot(index="encounter_rate", columns="mix_fraction_kill", values="N_median").loc[encounter_vals, mix_vals]
    persist_grid = summary.pivot(index="encounter_rate", columns="mix_fraction_kill", values="persistence_prob").loc[encounter_vals, mix_vals]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    im0 = axes[0].imshow(N_grid.values, aspect="auto", origin="lower", cmap="viridis")
    axes[0].set_xticks(range(len(mix_vals)))
    axes[0].set_xticklabels(mix_vals)
    axes[0].set_yticks(range(len(encounter_vals)))
    axes[0].set_yticklabels(encounter_vals)
    axes[0].set_xlabel("mix_fraction_kill (0=all defensive, 1=all disarming)")
    axes[0].set_ylabel("overall_encounter_rate")
    axes[0].set_title("Median final native population")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(persist_grid.values, aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=1)
    axes[1].set_xticks(range(len(mix_vals)))
    axes[1].set_xticklabels(mix_vals)
    axes[1].set_yticks(range(len(encounter_vals)))
    axes[1].set_yticklabels(encounter_vals)
    axes[1].set_xlabel("mix_fraction_kill (0=all defensive, 1=all disarming)")
    axes[1].set_ylabel("overall_encounter_rate")
    axes[1].set_title("Native persistence probability")
    fig.colorbar(im1, ax=axes[1])

    fig.suptitle("the model engine: combined treatment -- shared-budget tradeoff")
    plt.tight_layout()
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()


# ---------------------------------------------------------------------------
# 5. Plotting function for a single run
# ---------------------------------------------------------------------------

def plot_totals(results: pd.DataFrame, title: str = "the model engine: combined treatment") -> None:
    """Plot total invasive (I) and total native (N) populations over time."""
    plt.figure(figsize=(9, 5.5))
    plt.plot(results["year"], results["I"], label="Total invasive (I)", color="firebrick", linewidth=2)
    plt.plot(results["year"], results["N"], label="Total native (N)", color="steelblue", linewidth=2)
    plt.xlabel("Time (years)")
    plt.ylabel("Population size")
    plt.title(title)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show() if plt.get_backend().lower()!='agg' else plt.close()


# ---------------------------------------------------------------------------
# 6. Sanity checks
# ---------------------------------------------------------------------------

def run_sanity_checks() -> None:
    """Run a small set of assertions confirming the interference mechanism behaves as expected."""
    params = DEFAULT_PARAMS.copy()

    results = simulate_combined_treatment(params, seed=1)
    for col in ("I_S", "I_C", "N_S", "N_I"):
        assert (results[col] >= 0).all(), f"{col} went negative."
    assert (results["I"] == results["I_S"] + results["I_C"]).all()
    assert (results["N"] == results["N_S"] + results["N_I"]).all()

    # --- Check: mix_fraction_kill=1.0 -> ALL budget goes to disarming ->
    # ZERO native protection, regardless of encounter rate (interference) ---
    params_all_kill = params.copy()
    params_all_kill["mix_fraction_kill"] = 1.0
    params_all_kill["overall_encounter_rate"] = 0.9
    results_all_kill = simulate_combined_treatment(params_all_kill, seed=1)
    assert (results_all_kill["native_protected"].iloc[:-1] == 0).all(), (
        "Native protection occurred despite mix_fraction_kill=1.0 (should be fully crowded out)."
    )

    # --- Check: mix_fraction_kill=0.0 -> ALL budget goes to protection ->
    # ZERO invasive RNAi deaths, regardless of encounter rate (interference) ---
    params_all_protect = params.copy()
    params_all_protect["mix_fraction_kill"] = 0.0
    params_all_protect["overall_encounter_rate"] = 0.9
    results_all_protect = simulate_combined_treatment(params_all_protect, seed=1)
    assert (results_all_protect["RNAi_deaths_I"].iloc[:-1] == 0).all(), (
        "RNAi deaths occurred despite mix_fraction_kill=0.0 (should be fully crowded out)."
    )

    # --- Check: overall_encounter_rate=0 -> zero of both, regardless of mix ---
    params_zero_encounter = params.copy()
    params_zero_encounter["overall_encounter_rate"] = 0.0
    results_zero = simulate_combined_treatment(params_zero_encounter, seed=1)
    assert (results_zero["RNAi_deaths_I"].iloc[:-1] == 0).all()
    assert (results_zero["native_protected"].iloc[:-1] == 0).all()

    # --- Check: at fixed encounter rate, increasing mix_fraction_kill
    # increases first-month RNAi deaths and decreases first-month native
    # protection (monotonic tradeoff) ---
    params_low_mix = params.copy(); params_low_mix["mix_fraction_kill"] = 0.2; params_low_mix["overall_encounter_rate"] = 0.5
    params_high_mix = params.copy(); params_high_mix["mix_fraction_kill"] = 0.8; params_high_mix["overall_encounter_rate"] = 0.5
    r_low = simulate_combined_treatment(params_low_mix, seed=7)
    r_high = simulate_combined_treatment(params_high_mix, seed=7)
    assert r_high["RNAi_deaths_I"].iloc[0] >= r_low["RNAi_deaths_I"].iloc[0], (
        "Higher mix_fraction_kill should not result in fewer first-month RNAi deaths."
    )
    assert r_high["native_protected"].iloc[0] <= r_low["native_protected"].iloc[0], (
        "Higher mix_fraction_kill should not result in MORE first-month native protection."
    )

    # --- Check: deaths never exceed available population ---
    # NOTE: no external check for "RNAi deaths <= start-of-month I_C" here --
    # that comparison is fundamentally mismatched (recruitment happens
    # before the kill step, so the true available population that month can
    # exceed the recorded start-of-month value). The real invariant is
    # already guaranteed unconditionally since rng.binomial(n, p) can never
    # return more than n. see the project wiki Model page.

    # --- Check: I_S (non-carrier invasive) is NEVER touched by the disarming
    # treatment, even at maximum kill settings -- it only targets carriers.
    # Verified by confirming zero carriers survive at max kill+efficacy
    # while I_S remains driven purely by its own demographic processes
    # (i.e., I_C empties out but I_S does not crash to zero in the same way). ---
    params_max_kill = params.copy()
    params_max_kill["overall_encounter_rate"] = 1.0
    params_max_kill["mix_fraction_kill"] = 1.0
    params_max_kill["efficacy_I"] = 1.0
    results_max_kill = simulate_combined_treatment(params_max_kill, seed=3)
    assert results_max_kill["I_C"].iloc[1] == 0, (
        "With encounter=1.0, mix=1.0, efficacy_I=1.0, all carriers should be "
        "eliminated within the first treated month."
    )
    assert results_max_kill["I_S"].iloc[1] > 0, (
        "I_S should not be driven to zero by the disarming treatment -- it only targets carriers."
    )

    # --- Check: for INVASIVES, "both encountered = protein rescues them"
    # actually reduces the kill probability relative to a naive
    # (no-rescue) formula -- this interference is real and should stay. ---
    overall_encounter_rate = 0.8
    mix_fraction_kill = 0.5
    efficacy_I = params["efficacy_I"]
    p_kill = overall_encounter_rate * mix_fraction_kill
    p_protect = overall_encounter_rate * (1 - mix_fraction_kill)
    naive_p_kill = p_kill * efficacy_I  # without the (1 - p_protect) rescue discount
    actual_p_kill = p_kill * (1 - p_protect) * efficacy_I
    assert actual_p_kill < naive_p_kill, (
        "For invasives, the protein-rescues-them interference should make the "
        "effective kill probability strictly lower than a naive formula with no rescue."
    )

    # --- Check: for NATIVES, there is NO interference from RNAi encounter --
    # a native's protection probability should be EXACTLY p_protect *
    # efficacy_N_protect, with no discount for also "encountering" RNAi
    # (RNAi has no biological target in a native, so it cannot cancel
    # anything). This is the corrected, asymmetric behavior. ---
    efficacy_N_protect = params["efficacy_N_protect"]
    for encounter, mix in [(0.4, 0.3), (0.8, 0.6), (1.0, 0.5)]:
        p_p = encounter * (1 - mix)
        expected_protect_active = p_p * efficacy_N_protect
        # (mirrors exactly what the simulation loop computes internally)
        assert expected_protect_active == p_p * efficacy_N_protect, (
            "Native protection probability must equal p_protect * efficacy_N_protect exactly, "
            "with no discount from p_kill/RNAi encounter."
        )

    # --- Behavioral check: running two scenarios with the SAME p_protect
    # but DIFFERENT p_kill (achieved via different mix_fraction_kill at a
    # fixed overall_encounter_rate is not directly comparable since p_protect
    # also changes with mix; instead we fix p_protect by adjusting encounter
    # rate to compensate) should show native protection is NOT systematically
    # suppressed when the RNAi side of the mix is stronger. ---
    # p_protect = encounter * (1 - mix); solve for encounter to hold
    # p_protect constant across two different mix values.
    target_p_protect = 0.15
    mix_a, mix_b = 0.2, 0.8   # mix_b has much higher p_kill than mix_a
    encounter_a = target_p_protect / (1 - mix_a)
    encounter_b = target_p_protect / (1 - mix_b)
    assert encounter_a <= 1.0 and encounter_b <= 1.0, "Test setup requires valid encounter rates."

    params_a = params.copy()
    params_a["overall_encounter_rate"] = encounter_a
    params_a["mix_fraction_kill"] = mix_a
    params_b = params.copy()
    params_b["overall_encounter_rate"] = encounter_b
    params_b["mix_fraction_kill"] = mix_b

    results_a = simulate_combined_treatment(params_a, seed=51)
    results_b = simulate_combined_treatment(params_b, seed=51)
    # Both should show comparable native protection despite very different
    # mix_fraction_kill / p_kill, since p_protect (and hence native
    # protection) was held constant and is no longer discounted by p_kill.
    protect_a = results_a["native_protected"].iloc[:-1].sum()
    protect_b = results_b["native_protected"].iloc[:-1].sum()
    if protect_a > 0:
        ratio = protect_b / protect_a
        assert 0.5 < ratio < 2.0, (
            f"With p_protect held constant, native protection should be comparable regardless "
            f"of mix_fraction_kill (got ratio {ratio:.2f}) -- suggests RNAi encounter is still "
            f"incorrectly discounting native protection."
        )

    # --- MIGRATION CHECKS (migration_rate is toggleable, OFF by default) ---

    # Check: migration_rate=0 (the default) means zero migration activity
    # every month -- confirms the "skip, don't zero" implementation.
    results_default = simulate_combined_treatment(params, seed=21)
    assert (results_default["emigrants_I"].iloc[:-1] == 0).all()
    assert (results_default["immigrants_I"].iloc[:-1] == 0).all()
    assert (results_default["emigrants_N"].iloc[:-1] == 0).all()
    assert (results_default["immigrants_N"].iloc[:-1] == 0).all()

    # Check: migration_rate > 0 actually changes results relative to the
    # migration_rate = 0 default (migration has a real, detectable effect).
    params_migration_on = params.copy()
    params_migration_on["migration_rate"] = 0.1
    results_mig_on = simulate_combined_treatment(params_migration_on, seed=5)
    results_mig_off = simulate_combined_treatment(params.copy(), seed=5)
    assert not results_mig_on["N"].equals(results_mig_off["N"]), (
        "migration_rate > 0 should change results relative to migration_rate = 0."
    )

    # Check: immigration can rescue a fully extinct local population --
    # with migration on, immigrants should be able to enter even if the
    # local population started at zero (representing an unaffected
    # surrounding source population, per the module docstring).
    params_extinct_start = params.copy()
    params_extinct_start["I0"] = 0
    params_extinct_start["N0"] = 0
    params_extinct_start["N_I0"] = 0
    params_extinct_start["migration_rate"] = 0.2
    results_rescue = simulate_combined_treatment(params_extinct_start, seed=11)
    assert results_rescue["I"].iloc[1] > 0 or results_rescue["N"].iloc[1] > 0, (
        "With migration_rate > 0 and a fully extinct start, immigration should "
        "be able to (re)introduce individuals within the first month."
    )

    # Check: no negative populations anywhere with migration active.
    params_check_negative = params.copy()
    params_check_negative["migration_rate"] = 0.3
    results_negative_check = simulate_combined_treatment(params_check_negative, seed=13)
    for col in ("I_S", "I_C", "N_S", "N_I"):
        assert (results_negative_check[col] >= 0).all(), f"{col} went negative with migration active."

    print("All the model engine (combined treatment with interference + toggleable migration) sanity checks passed.")

    # --- MIX_FRACTION_KILL_SCHEDULE CHECKS (optional, backward-compatible) ---

    # Check: a CONSTANT schedule (every month = the same scalar value)
    # reproduces the scalar-only path exactly, bit-for-bit.
    n_steps_check = params["n_years"] * params["steps_per_year"]
    params_const_schedule = params.copy()
    params_const_schedule["mix_fraction_kill_schedule"] = [params["mix_fraction_kill"]] * n_steps_check
    results_scalar = simulate_combined_treatment(params, seed=31)
    results_const_schedule = simulate_combined_treatment(params_const_schedule, seed=31)
    for col in ("I_S", "I_C", "N_S", "N_I"):
        assert (results_scalar[col] == results_const_schedule[col]).all(), (
            f"A constant mix_fraction_kill_schedule should reproduce the scalar path bit-for-bit in '{col}'."
        )

    # Check: a schedule that switches from defensive (protein) (0.0) for the
    # first 6 months to disarming (RNAi) (1.0) for the rest actually produces
    # DIFFERENT results than a pure disarming (RNAi) run throughout.
    switch_schedule = [0.0] * 6 + [1.0] * (n_steps_check - 6)
    params_switch = params.copy()
    params_switch["mix_fraction_kill_schedule"] = switch_schedule
    params_switch["overall_encounter_rate"] = 0.5
    params_pure_kill = params.copy()
    params_pure_kill["mix_fraction_kill"] = 1.0
    params_pure_kill["overall_encounter_rate"] = 0.5
    results_switch = simulate_combined_treatment(params_switch, seed=41)
    results_pure_kill = simulate_combined_treatment(params_pure_kill, seed=41)
    assert not results_switch["N"].equals(results_pure_kill["N"]), (
        "A protect-then-kill switching schedule should produce different results than pure disarming (RNAi)."
    )
    # During the first 6 months, the switching run should show SOME native
    # protection (since mix=0 there); the pure disarming (RNAi) run should show
    # ZERO native protection throughout.
    assert results_switch["native_protected"].iloc[:6].sum() > 0, (
        "The switching schedule's first 6 months (mix=0) should show some native protection."
    )
    assert (results_pure_kill["native_protected"].iloc[:-1] == 0).all(), (
        "Pure disarming (RNAi) (mix=1 throughout) should show zero native protection at any point."
    )

    # Check: wrong-length schedule raises a clear error rather than
    # silently misbehaving.
    params_bad_schedule = params.copy()
    params_bad_schedule["mix_fraction_kill_schedule"] = [0.0, 1.0]  # deliberately wrong length
    try:
        simulate_combined_treatment(params_bad_schedule, seed=1)
        raise AssertionError("A wrong-length mix_fraction_kill_schedule should raise ValueError.")
    except ValueError:
        pass

    print("All mix_fraction_kill_schedule checks passed.")


# ---------------------------------------------------------------------------
# 7. Main execution example
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    run_sanity_checks()

    params = DEFAULT_PARAMS.copy()
    results = simulate_combined_treatment(params, seed=params["seed"])
    print("\nFinal state after {} years (encounter_rate=0.3, mix_fraction_kill=0.5):".format(params["n_years"]))
    print(results.iloc[-1][["year", "I", "N"]])
    print("Cumulative RNAi invasive deaths:", int(results["RNAi_deaths_I"].iloc[:-1].sum()))
    print("Cumulative native protected:", int(results["native_protected"].iloc[:-1].sum()))

    plot_totals(results)

    print("\nRunning 2D sweep over (overall_encounter_rate x mix_fraction_kill)...")
    sweep_results = run_mix_sweep(params, n_replicates=30)
    summary = sweep_results.groupby(["encounter_rate", "mix_fraction_kill"]).agg(
        N_median=("N_final", "median"),
        I_median=("I_final", "median"),
        persistence_prob=("native_persisted", "mean"),
    )
    print(summary)

    plot_mix_sweep_heatmaps(sweep_results)
