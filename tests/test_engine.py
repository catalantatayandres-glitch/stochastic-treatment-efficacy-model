import pytest

from stochastic_treatment_efficacy_model import DEFAULT_PARAMS, engine, simulate_combined_treatment


def test_engine_sanity_checks():
    """The engine's built-in invariants (interference, migration, schedules)."""
    engine.run_sanity_checks()


def test_output_shape_and_conservation():
    res = simulate_combined_treatment(DEFAULT_PARAMS, seed=0)
    n_steps = DEFAULT_PARAMS["n_years"] * DEFAULT_PARAMS["steps_per_year"]
    assert len(res) == n_steps + 1
    assert (res["I"] == res["I_S"] + res["I_C"]).all()
    assert (res["N"] == res["N_S"] + res["N_I"]).all()
    assert (res[["I_S", "I_C", "N_S", "N_I"]] >= 0).all().all()


def test_same_seed_is_reproducible():
    a = simulate_combined_treatment(DEFAULT_PARAMS, seed=123)
    b = simulate_combined_treatment(DEFAULT_PARAMS, seed=123)
    assert a.equals(b)


@pytest.mark.parametrize("key, value", [("mu_A", 1.5), ("K_I", 0), ("transmission_mode", "airborne")])
def test_invalid_parameters_are_rejected(key, value):
    with pytest.raises(ValueError):
        simulate_combined_treatment({**DEFAULT_PARAMS, key: value}, seed=0)
