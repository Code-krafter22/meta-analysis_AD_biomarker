#!/usr/bin/env python3
"""Regenerate grouped validation for both cohorts and build Figure 5."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[4]
validation = root / "scripts" / "04_external_validation" / "05_grouped_validation.py"
figure = Path(__file__).with_name("01_grouped_cv_figure.py")

if not os.getenv("MSBB_EXPR") or not os.getenv("MSBB_META"):
    raise SystemExit(
        "Set MSBB_EXPR and MSBB_META to authorized local MSBB files. "
        "Restricted MSBB source data must not be committed to this repository."
    )

# Keep plotting caches outside the repository and avoid user-home writes.
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="figure5_mpl_"))

for cohort in ("GSE125583", "MSBB"):
    env = os.environ.copy()
    env.update({"VALIDATION_MODE": "grouped_cv", "VALIDATION_COHORT": cohort})
    subprocess.run([sys.executable, str(validation)], cwd=root, env=env, check=True)

subprocess.run([sys.executable, str(figure)], cwd=root, check=True)
print("Figure 5 reproducibility workflow completed.")
