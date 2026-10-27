#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lightweight plotting: overlay original Python vs C++ profiles.

Reads the SAVED profiles produced by the acceptance harness — the C++
``profile.csv`` and the Python reference ``profile_python.csv`` — and overlays
them.  This does NOT rerun either solver (review fix #6/#7: lightweight, no
source import, no simulation).  Saves one PNG under validation/.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify  # noqa: E402  (only for ROOT; no solver is imported/run)

ROOT = verify.ROOT
CASES = os.path.join(ROOT, "validation", "cases")


def load_cpp(tag):
    prof = np.loadtxt(os.path.join(CASES, tag, "profile.csv"), delimiter=",", skiprows=1)
    return prof[:, 0], prof[:, 1]


def load_python(tag):
    prof = np.loadtxt(os.path.join(CASES, tag, "profile_python.csv"),
                      delimiter=",", skiprows=1)
    return prof[:, 0], prof[:, 1]


fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# laminar representative
tag = "lam_1000000_100_0.5_100_1_4"
z_c, rho_c = load_cpp(tag)
z_p, rho_p = load_python(tag)
axes[0].plot(z_c, rho_c, "-", color="tab:blue", label="C++")
axes[0].plot(z_p, rho_p, "--", color="black", lw=2, label="Python")
axes[0].set_title("Laminar dust profile (zmax=4)")
axes[0].set_xlabel("z / H")
axes[0].set_ylabel("rho_dust")
axes[0].legend()

# turbulence crossing (stop=True -> early stop at formation ratio>1)
tag = "turb_cross_stop"
z_c, rho_c = load_cpp(tag)
z_p, rho_p = load_python(tag)
axes[1].plot(z_c, rho_c, "-", color="tab:red", label="C++")
axes[1].plot(z_p, rho_p, "--", color="black", lw=2, label="Python")
axes[1].set_title("Turbulence, capped + L clamped (zmax=50)")
axes[1].set_xlabel("z / H")
axes[1].set_ylabel("rho_dust")
axes[1].legend()

fig.tight_layout()
out = os.path.join(ROOT, "validation", "profiles_compare.png")
fig.savefig(out, dpi=130)
print("saved", out)
