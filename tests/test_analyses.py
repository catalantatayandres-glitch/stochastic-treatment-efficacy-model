"""Fast end-to-end smoke tests: each analysis runs its own reduced pipeline."""
import importlib

import pytest

ANALYSES = [
    "main_results",
    "validation_bubb2006",
    "sensitivity_analysis",
    "scenario_robustness",
    "coverage_efficacy",
    "migration_effect",
]


@pytest.mark.parametrize("module_name", ANALYSES)
def test_analysis_sanity_checks(module_name):
    module = importlib.import_module(module_name)
    module.run_sanity_checks()
