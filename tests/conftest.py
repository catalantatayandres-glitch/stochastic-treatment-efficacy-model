import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # never open plot windows during tests

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(REPO_ROOT / "analyses"))
