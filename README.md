# ASTACI – Stochastic crayfish population model (iGEM 2026 UNILausanne)

Monthly stochastic model (120 steps, 10 years) comparing a **disarming strategy** (RNAi, `mix_fraction_kill = 1`) with a **defensive strategy** (protein supplementation, `mix_fraction_kill = 0`).

## Run
```
pip install -r requirements.txt
python 03_main_results.py
```
Figures are written to `figures/`. Scripts 02–07 import `01_model_engine.py`, so keep all files in the same folder.

| Script | Wiki section | Output |
|---|---|---|
| `01_model_engine.py` | How the Model Works / Parameters | engine (`simulate_combined_treatment`) + sanity checks |
| `02_sensitivity_analysis.py` | Other explorations: sensitivity (one-at-a-time) | `sensitivity_tornado_wiki.png` |
| `03_main_results.py` | Results | `wiki_figure_transmission_0.9.png` |
| `04_validation_bubb2006.py` | Validation (Bubb, Thom & Lucas 2006) | `validation_figure1_grid.png`, `validation_figure2_selected_scenario.png` |
| `05_habitat_competition_transmission.py` | Other explorations: 27 scenarios | `other_explorations_habitat_competition_transmission.png` |
| `06_coverage_efficacy.py` | Other explorations: coverage × efficacy | `other_explorations_coverage_efficacy_absolute.png` |
| `07_migration_effect.py` | Other explorations: migration | `other_explorations_migration_effect.png` |

Reference scenario: Neutral habitat (K_I = 6400, K_N = 5500), moderate competition (α_IN = 0.3, α_NI = 0.5), β = 0.9, coverage 50 %, efficacy 100 %, I0 = 0.2·K_I, N0 = 0.2·K_N.
