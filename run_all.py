"""Regenerate every figure in figures/ by running each analysis in turn.

    python run_all.py                  # all analyses
    python run_all.py main_results     # just one (or several) by name
"""
import os
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

# Ordered as on the wiki: Results -> Validation -> Other explorations
ANALYSES = [
    "main_results",
    "validation_bubb2006",
    "sensitivity_analysis",
    "scenario_robustness",
    "coverage_efficacy",
    "migration_effect",
]


def main() -> int:
    selected = sys.argv[1:] or ANALYSES
    unknown = set(selected) - set(ANALYSES)
    if unknown:
        print(f"Unknown analysis: {', '.join(sorted(unknown))}. Choose from: {', '.join(ANALYSES)}")
        return 2

    env = {**os.environ, "MPLBACKEND": "Agg"}  # save figures without opening windows
    for name in selected:
        print(f"\n=== {name} ===", flush=True)
        start = time.time()
        result = subprocess.run([sys.executable, str(REPO_ROOT / "analyses" / f"{name}.py")], env=env)
        if result.returncode != 0:
            print(f"{name} failed (exit code {result.returncode}).")
            return result.returncode
        print(f"=== {name} done in {time.time() - start:.0f}s ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
