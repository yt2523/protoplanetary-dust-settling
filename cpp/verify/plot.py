#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lightweight plotting: overlay original Python vs C++ profiles.

Reads C++ profile.csv files produced by the acceptance harness and overlays the
corresponding original Python reference profiles.  Saves one PNG under
validation/.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify  # noqa: E402

ROOT = verify.ROOT
CASES = os.path.join(ROOT, "validation", "cases")

# The reference classes are imported lazily; resolve them here.
verify._import_sources(verify.SRC)


def load_cpp(tag):
    prof = np.loadtxt(os.path.join(CASES, tag, "profile.csv"), delimiter=",", skiprows=1)
    return prof[:, 0], prof[:, 1]


fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# laminar representative
tag = "lam_1000000_100_0.5_100_1_4"
z_c, rho_c = load_cpp(tag)
p_diag, p_prof = verify.run_python_laminar(1000000, 100, 0.5, 100, 1.0, 4, 1e3, False)
axes[0].plot(z_c, rho_c, "-", color="tab:blue", label="C++")
axes[0].plot(p_prof[0], p_prof[1], "--", color="black", lw=2, label="Python")
axes[0].set_title("Laminar dust profile (zmax=4)")
axes[0].set_xlabel("z / H")
axes[0].set_ylabel("rho_dust")
axes[0].legend()

# turbulence crossing (stop=True -> early stop at formation ratio>1)
tag = "turb_cross_stop"
z_c, rho_c = load_cpp(tag)
p_diag, p_prof = verify.run_python_turbulence(10000, 400, 0.3, 0.2, 100, 1000.0, 50, 1e-5, 1e3, 1.0, True)
axes[1].plot(z_c, rho_c, "-", color="tab:red", label="C++")
axes[1].plot(p_prof[0], p_prof[1], "--", color="black", lw=2, label="Python")
axes[1].set_title("Turbulence, capped + L clamped (zmax=50)")
axes[1].set_xlabel("z / H")
axes[1].set_ylabel("rho_dust")
axes[1].legend()

fig.tight_layout()
out = os.path.join(ROOT, "validation", "profiles_compare.png")
fig.savefig(out, dpi=130)
print("saved", out)
