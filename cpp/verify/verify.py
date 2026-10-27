#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C++ port acceptance harness (TDD).

Imports the *original* Python modules (not reimplementations) and compares them
against the C++ CLI output.  The C++ binary is a Linux ELF produced by the
D-backed Ubuntu WSL toolchain, so it is invoked through `wsl.exe`; all files are
written to D: and are visible from both Windows and WSL.

Exit codes:
  0  all checks passed
  1  at least one check failed
  2  build artifact missing (pre-build TDD demonstration)
"""

import os
import sys
import json
import math
import subprocess
import contextlib
import io
import csv
import time
import argparse

import numpy as np

# ---------------------------------------------------------------- paths
# Root is derived from this file (portable harness, review item 11).
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Source repo: default the sibling parent if it contains the original
# src/pure_advection.py, otherwise the existing D: checkout.
_DEFAULT_SRC = r"D:/Codex/tmp/phd-github-audit-2026-10-26/protoplanetary-dust-settling/src"
_ALT_SRC = os.path.join(os.path.dirname(ROOT), "src", "pure_advection.py")


def _default_source():
    if os.path.exists(_ALT_SRC):
        return os.path.dirname(_ALT_SRC)
    return _DEFAULT_SRC


SRC = _default_source()

# Binary defaults to the CMake-built artifact (build-cmake/dust_settling).
_BIN_REL = os.path.join(ROOT, "build-cmake", "dust_settling")

CASES_DIR = os.path.join(ROOT, "validation", "cases")
RESULTS_CSV = os.path.join(ROOT, "validation", "results.csv")
RESULTS_JSON = os.path.join(ROOT, "validation", "results.json")

# Original-module symbols are imported lazily in main() AFTER argparse resolves
# --source-repo (review item 15/11): importing before parsing then mutating
# sys.path does not reload already-imported modules.  These are module-level
# placeholders filled by _import_sources().
_advection_cls = None
_laminar_cls = None
_turbulence_cls = None

_SRC_MODULES = []


def _import_sources(source_arg):
    """Import the three original classes from the resolved source directory.

    Accepts either a repo ROOT (containing src/...) or the src/ dir itself.
    """
    global _advection_cls, _laminar_cls, _turbulence_cls, _SRC_MODULES
    cand = source_arg
    if os.path.isdir(os.path.join(cand, "src")) and os.path.exists(
            os.path.join(cand, "src", "pure_advection.py")):
        cand = os.path.join(cand, "src")
    if not os.path.isdir(cand):
        raise SystemExit("--source-repo is not a directory: %r" % source_arg)
    if not os.path.exists(os.path.join(cand, "pure_advection.py")):
        raise SystemExit("--source-repo has no pure_advection.py: %r" % cand)
    if cand not in sys.path:
        sys.path.insert(0, cand)
    import pure_advection
    import dust_solver
    import dust_solver_turbulence
    _advection_cls = pure_advection.PureAdvection1D_v2
    _laminar_cls = dust_solver.RealDiskDustSolver
    _turbulence_cls = dust_solver_turbulence.DustSolver_Turbulence
    _SRC_MODULES = [pure_advection, dust_solver, dust_solver_turbulence]
    return cand


def _src_identity():
    """Hash the currently-imported source modules (also record before/after)."""
    import hashlib
    d = {}
    for m in _SRC_MODULES:
        p = getattr(m, "__file__", None)
        if p and os.path.exists(p):
            with open(p, "rb") as f:
                d[os.path.basename(p)] = hashlib.sha256(f.read()).hexdigest()
    return d


# ---------------------------------------------------------------- tolerances
EPS_ABS = 1e-12
EPS_REL = 1e-10
MASS_TOL = 1e-10

NAN = float("nan")


def is_nan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def g(x):
    return "%.17g" % x


# ---------------------------------------------------------------- analytic
def gauss_exact(x, v, t, L=1.0, x0=0.5, sigma=0.07):
    # Strict continuous translation using the INITIAL sampled normalization.
    # The original PureAdvection1D_v2.initial_function_points divides the sampled
    # Gaussian by its own (sampled, grid-dependent) max so the peak is exactly 1;
    # re-normalizing the shifted analytic peak would differ by a tiny O(h^2)
    # amount.  This helper normalizes by the SAME unshifted sampled maximum so the
    # analytic reference is a strict continuous translation of the initial data
    # (convergence on even grids is unaffected — the difference vanishes as h->0).
    x_sampled = np.linspace(0.0, L, len(x), endpoint=False)
    y0 = np.exp(-0.5 * ((x_sampled - x0) / sigma) ** 2)
    norm = y0.max()
    x_back = np.mod(x - v * t, L)
    y = np.exp(-0.5 * ((x_back - x0) / sigma) ** 2)
    return y / norm


# ---------------------------------------------------------------- Python references
def _instrument_upwind(sim, top_flux_fn):
    """Wrap the ORIGINAL upwind_method (do NOT reimplement evolution).

    Records, over every update plus the initial condition:
      - min_rho   : min over every rho array seen (initial + each update)
      - finite    : all values finite over initial + each update
      - calls     : number of times upwind_method was actually invoked
      - cum_flux  : sum(dt * top_flux(pre-step rho))  (bottom flux is 0)
    top_flux_fn(sim, rho_pre) returns the outward flux at the top boundary from
    the PRE-STEP density (matches the C++ cumulative boundary integral).
    """
    rec = {"min_rho": None, "finite": True, "calls": 0, "cum_flux": 0.0}
    orig = sim.upwind_method

    def update_min_finite(arr):
        a = np.asarray(arr, dtype=float)
        if a.size == 0:
            return
        m = float(np.min(a))
        if rec["min_rho"] is None or m < rec["min_rho"]:
            rec["min_rho"] = m
        if not np.all(np.isfinite(a)):
            rec["finite"] = False

    def wrapped(rho):
        rec["calls"] += 1
        rho_pre = np.asarray(rho, dtype=float)
        f_top = top_flux_fn(sim, rho_pre)
        rec["cum_flux"] += float(sim.dt) * f_top
        out = orig(rho)
        update_min_finite(out)
        return out

    sim.upwind_method = wrapped
    return rec


def _laminar_top_flux(sim, rho_pre):
    # v[-1] >= 0 -> outflow v[-1]*rho[-1]; else 0 (dust density 0 at infinity).
    return max(sim.v[-1], 0.0) * rho_pre[-1]


def _turb_top_flux(sim, rho_pre):
    # advection outflow + upper diffusion ghost (ratio=0, gas=0):
    # F_diff = -alpha * 0.5*rho_g[-1] * (0 - ratio[-1])/h = alpha*rho[-1]/(2h)
    return max(sim.v[-1], 0.0) * rho_pre[-1] + sim.alpha * rho_pre[-1] / (2.0 * sim.h)


def _py_prof(sim):
    z = np.asarray(sim.z, dtype=float)
    rho_final = np.asarray(sim.rhofinal, dtype=float)
    rho_init = np.asarray(sim.rho_initial, dtype=float)
    return z, rho_final, rho_init


def run_python_advection(v, dt, T, Nx):
    sim = _advection_cls(c_velocity=v, time_step=dt, total_time=T, spatial_step=Nx)
    rec = _instrument_upwind(sim, lambda s, r: 0.0)
    with contextlib.redirect_stdout(io.StringIO()):
        sim.run_simulation()
    x = np.asarray(sim.x, dtype=float)
    rho_init = gauss_exact(x, v, 0.0)
    rho_final = np.asarray(sim.rhofinal, dtype=float)
    init_arr = np.asarray(sim.rho, dtype=float)
    min_all = float(np.min(init_arr))
    if rec["min_rho"] is not None:
        min_all = min(min_all, rec["min_rho"])
    finite_all = bool(np.all(np.isfinite(init_arr))) and rec["finite"]
    actual_steps = rec["calls"]
    mass_initial = float(sim.h * np.sum(init_arr))
    mass_final = float(sim.h * np.sum(rho_final))
    return {
        "mode": "advection",
        "N": Nx,
        "L": 1.0,
        "h": float(sim.h),
        "dt": float(sim.dt),
        "Nt": int(sim.Nt),
        "t_actual": float(actual_steps * sim.dt),
        "actual_steps": actual_steps,
        "mass_initial": mass_initial,
        "mass_final": mass_final,
        "min_rho": min_all,
        "max_rho": float(np.max(np.abs(rho_final))),
        "finite": finite_all,
        "cumulative_boundary_flux": rec["cum_flux"],
        "calls": rec["calls"],
        "t_formation": NAN,
        "t_formation_years": NAN,
        "t_actual_years": NAN,
        "v_max": float(abs(v)),
        "alpha": NAN,
        "z_real_max": NAN,
        "real_z_max": NAN,
        "rho_initial_max": float(np.max(np.abs(init_arr))),
        "scale_height": NAN,
    }, (x, rho_final, rho_init, x.copy())


def run_python_laminar(total_time, Nz, a_lim, R_AU, a_micron, zmax, dust_density=1e3,
                       stop=False):
    s = _laminar_cls(total_time, Nz, a_lim, R_AU, a_micron, zmax, dust_density)
    rec = _instrument_upwind(s, _laminar_top_flux)
    s.run_simulation(stop_when_cross=stop)
    z, rho_final, rho_init = _py_prof(s)
    z_real = np.asarray(s.z_real, dtype=float)
    init_arr = np.asarray(s.rho, dtype=float)
    min_all = float(np.min(init_arr))
    if rec["min_rho"] is not None:
        min_all = min(min_all, rec["min_rho"])
    finite_all = bool(np.all(np.isfinite(init_arr))) and rec["finite"]
    tf = s.t_formation if s.t_formation is not None else NAN
    tfy = s.get_t_formation_years()
    tfy = tfy if tfy is not None else NAN
    actual_steps = rec["calls"]
    t_actual = float(actual_steps * s.dt)
    mass_initial = float(s.h * np.sum(init_arr))
    mass_final = float(s.h * np.sum(rho_final))
    return {
        "mode": "laminar",
        "N": Nz,
        "L": float(s.L),
        "h": float(s.h),
        "dt": float(s.dt),
        "Nt": int(s.Nt),
        "t_actual": t_actual,
        "actual_steps": actual_steps,
        "mass_initial": mass_initial,
        "mass_final": mass_final,
        "min_rho": min_all,
        "max_rho": float(np.max(np.abs(rho_final))),
        "finite": finite_all,
        "cumulative_boundary_flux": rec["cum_flux"],
        "calls": rec["calls"],
        "t_formation": float(tf),
        "t_formation_years": float(tfy),
        "t_actual_years": float(t_actual / s.omega / (365.25 * 24 * 3600)),
        "v_max": float(np.max(np.abs(s.v))),
        "alpha": NAN,
        "z_real_max": float(np.max(z_real)),
        "real_z_max": float(zmax * 1.495978707e11),
        "rho_initial_max": float(np.max(np.abs(rho_init))),
        "scale_height": float(s.scale_height),
    }, (z, rho_final, rho_init, z_real)


def run_python_turbulence(total_time, Nz, a_lim, d_lim, R_AU, a_micron, zmax, D=1e-5,
                          dust_density=1e3, v_cap_factor=1.0, stop=False):
    s = _turbulence_cls(total_time, Nz, a_lim, d_lim, R_AU, a_micron, zmax, D, dust_density)
    s.v_cap_factor = v_cap_factor
    rec = _instrument_upwind(s, _turb_top_flux)
    s.run_simulation(stop_when_cross=stop)
    z, rho_final, rho_init = _py_prof(s)
    z_real = np.asarray(s.z_real, dtype=float)
    init_arr = np.asarray(s.rho, dtype=float)
    min_all = float(np.min(init_arr))
    if rec["min_rho"] is not None:
        min_all = min(min_all, rec["min_rho"])
    finite_all = bool(np.all(np.isfinite(init_arr))) and rec["finite"]
    tf = s.t_formation if s.t_formation is not None else NAN
    tfy = s.get_t_formation_years()
    tfy = tfy if tfy is not None else NAN
    actual_steps = rec["calls"]
    t_actual = float(actual_steps * s.dt)
    mass_initial = float(s.h * np.sum(init_arr))
    mass_final = float(s.h * np.sum(rho_final))
    return {
        "mode": "turbulence",
        "N": Nz,
        "L": float(s.L),
        "h": float(s.h),
        "dt": float(s.dt),
        "Nt": int(s.Nt),
        "t_actual": t_actual,
        "actual_steps": actual_steps,
        "mass_initial": mass_initial,
        "mass_final": mass_final,
        "min_rho": min_all,
        "max_rho": float(np.max(np.abs(rho_final))),
        "finite": finite_all,
        "cumulative_boundary_flux": rec["cum_flux"],
        "calls": rec["calls"],
        "t_formation": float(tf),
        "t_formation_years": float(tfy),
        "t_actual_years": float(t_actual / s.omega / (365.25 * 24 * 3600)),
        "v_max": float(np.max(np.abs(s.v))),
        "alpha": float(s.alpha),
        "z_real_max": float(np.max(z_real)),
        "real_z_max": float(zmax * 1.495978707e11),
        "rho_initial_max": float(np.max(np.abs(rho_init))),
        "scale_height": float(s.scale_height),
    }, (z, rho_final, rho_init, z_real)


# ---------------------------------------------------------------- C++ runner
def _windows_to_linux_path(win_path):
    """Drive-letter Windows path -> /mnt/<drive>/... WSL path."""
    p = os.path.abspath(win_path).replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        p = "/mnt/" + p[0].lower() + p[2:]
    return p


_ACTIVE_BIN = None  # set by main() from argparse


def _save_py_profile(tag, p_prof):
    """Save the Python reference profile next to the C++ profile (review B.10)."""
    coord, rho, rhoi, aux = p_prof
    outdir = os.path.join(CASES_DIR, tag)
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "profile_python.csv")
    with open(path, "w", newline="") as f:
        f.write("coord,rho,rho_init,aux\n")
        for i in range(len(coord)):
            f.write("%.17g,%.17g,%.17g,%.17g\n"
                    % (coord[i], rho[i], rhoi[i], aux[i]))


def run_cpp(mode, params, tag, bin_win=None):
    bin_win = bin_win or _ACTIVE_BIN or _BIN_REL
    bin_linux = _windows_to_linux_path(bin_win)
    win_out = os.path.join(CASES_DIR, tag)
    linux_cases = _windows_to_linux_path(CASES_DIR)
    linux_out = linux_cases.rstrip("/") + "/" + tag
    os.makedirs(win_out, exist_ok=True)
    # Direct argv via `wsl --exec`, no shell interpolation (review item 11/15);
    # 30 s timeout, no retry loop (stale output / flaky retries are not a pass).
    full_argv = [bin_linux, mode] + [str(p) for p in params] + [linux_out]
    argv_txt = " ".join(full_argv)
    try:
        proc = subprocess.run(
            ["wsl.exe", "--exec"] + full_argv,
            capture_output=True, text=True, timeout=30,
        )
    except subprocess.TimeoutExpired as e:
        raise TimeoutError("C++ binary timed out (30s) for %s" % tag) from e
    if proc.returncode != 0:
        raise RuntimeError("C++ binary failed rc=%s\nargv=%s\nstdout=%s\nstderr=%s"
                           % (proc.returncode, argv_txt, proc.stdout, proc.stderr))
    diag = {}
    with open(os.path.join(win_out, "diagnostics.txt"), "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            k, v = parts
            diag[k] = v

    def getf(key):
        return float(diag[key])

    def geti(key):
        return int(diag[key])

    c_diag = {
        "mode": diag["mode"],
        "N": geti("N"),
        "L": getf("L"),
        "h": getf("h"),
        "dt": getf("dt"),
        "Nt": geti("Nt"),
        "t_actual": getf("t_actual"),
        "actual_steps": geti("actual_steps"),
        "t_actual_years": getf("t_actual_years"),
        "mass_initial": getf("mass_initial"),
        "mass_final": getf("mass_final"),
        "cumulative_boundary_flux": getf("cumulative_boundary_flux"),
        "mass_balance_error": getf("mass_balance_error"),
        "min_rho": getf("min_rho"),
        "max_rho": getf("max_rho"),
        "finite": int(diag["finite"]),
        "t_formation": getf("t_formation"),
        "t_formation_years": getf("t_formation_years"),
        "v_max": getf("v_max"),
        "alpha": getf("alpha"),
        "z_real_max": getf("z_real_max"),
        "real_z_max": getf("real_z_max"),
        "rho_initial_max": getf("rho_initial_max"),
        "scale_height": getf("scale_height"),
        "cap_count": diag.get("cap_count", "0"),
    }
    prof = np.loadtxt(os.path.join(win_out, "profile.csv"), delimiter=",", skiprows=1)
    coord = prof[:, 0]
    rho_final = prof[:, 1]
    rho_init = prof[:, 2]
    aux = prof[:, 3]
    return c_diag, (coord, rho_final, rho_init, aux)


# ---------------------------------------------------------------- comparison
def masked_rel(C, P):
    peak = np.max(np.abs(P))
    mask = np.abs(P) > 1e-8 * peak
    if not np.any(mask):
        return 0.0
    return float(np.max(np.abs((C[mask] - P[mask]) / P[mask])))


def rel_err(a, b):
    if is_nan(a) and is_nan(b):
        return 0.0
    return abs(a - b) / max(1.0, abs(b))


def fp_close(a, b):
    """EPS_ABS + EPS_REL*abs(P) — the mandated tolerance for every float comp."""
    return abs(a - b) <= EPS_ABS + EPS_REL * abs(b)


def _float_gate(name, c_val, p_val, r):
    """Record abs diff and the strict EPS_ABS+EPS_REL*abs(P) verdict."""
    r[name + "_abs_diff"] = abs(c_val - p_val)
    r[name + "_ok"] = fp_close(c_val, p_val)
    return r[name + "_ok"]


def branch_coverage():
    """Concrete branch-coverage assertions (review B.14), not prose.

    Returns a results dict {assertion: bool} covering:
      - laminar z/H clamp: L==4 while zmax/H>4, physical grid != z*H
      - turbulence grid: z_real == z*H
      - turbulence velocity cap count > 0 from raw (unclipped) velocities
      - crossing stop/continue semantics + tf step identity
      - representative diffusion: alpha ~ 1e-3, measurable cumulative loss
    """
    def run(tag):
        for mode, args, t in _cases():
            if t == tag:
                return run_cpp(mode, args, tag)
        raise KeyError(tag)

    cov = {}

    # laminar clamp: R=100, zmax=50 AU -> zmax/H > 4, L clamped to 4.
    d, prof = run("lam_clamp")
    z = prof[0]; z_real = prof[3]
    H = d["scale_height"]
    cov["lam_L_eq_4"] = abs(d["L"] - 4.0) <= EPS_ABS
    cov["lam_zmax_over_H_gt_4"] = (d["real_z_max"] / H) > 4.0
    # physical z grid differs from z*H (laminar uses linspace z_real, not z*H)
    z_times_H = z * H
    cov["lam_z_real_differs_zH"] = float(np.max(np.abs(z_real - z_times_H))) > EPS_ABS

    # turbulence grid IS z*H.
    d2, prof2 = run("turb_1000000_100_0.3_0.2_100_1_4_1e-05_1_0")
    z2 = prof2[0]; z_real2 = prof2[3]; H2 = d2["scale_height"]
    cov["turb_z_real_eq_zH"] = bool(np.all(np.abs(z_real2 - z2 * H2) <= EPS_ABS + EPS_REL * np.abs(z2 * H2)))

    # crossing: stop=True has tf and steps < Nt; continue runs full Nt and same tf.
    d_stop, _ = run("turb_cross_stop")
    d_cont, _ = run("turb_cross_continue")
    # velocity capping must actually engage (a_micron=1000, zmax=50 settle fast).
    cov["turb_cap_count_gt_0"] = int(d_stop.get("cap_count", 0)) > 0
    cov["stop_tf_exists"] = not is_nan(d_stop["t_formation"])
    cov["stop_steps_lt_Nt"] = int(d_stop["actual_steps"]) < int(d_stop["Nt"])
    cov["cont_steps_eq_Nt"] = int(d_cont["actual_steps"]) == int(d_cont["Nt"])
    cov["cont_Nt_gt_formation_step"] = int(d_cont["Nt"]) > int(round(d_cont["t_formation"] / d_cont["dt"]))
    cov["same_tf_stop_cont"] = fp_close(d_stop["t_formation"], d_cont["t_formation"])

    # representative diffusion: alpha ~ 1e-3, nonzero cumulative loss matching M.
    d3, prof3 = run("turb_diff_repr")
    cov["repr_alpha_near_1e3"] = fp_close(d3["alpha"], 1e-3)
    M_init = d3["mass_initial"]
    cov["repr_cum_loss_gt_1e8"] = abs(d3["cumulative_boundary_flux"]) > 1e-8 * M_init
    cov["repr_loss_matches_M"] = abs(d3["cumulative_boundary_flux"]) > 0

    return cov


def run_branch_coverage():
    cov = branch_coverage()
    all_ok = all(cov.values())
    print("[%s] branch_coverage  %s" % ("PASS" if all_ok else "FAIL",
                                        json.dumps(cov, sort_keys=True)))
    return {"tag": "branch_coverage", "status": "PASS" if all_ok else "FAIL",
            "note": "failed=%r" % [k for k, v in cov.items() if not v] if not all_ok else "",
            "coverage": cov}


def check_case(tag, p_diag, p_prof, mode, cpp_args):
    """Compare C++ output against Python reference; return a results dict."""
    r = {"tag": tag, "mode": p_diag["mode"], "N": p_diag["N"], "status": "PASS", "note": ""}
    c_diag, c_prof = run_cpp(mode, cpp_args, tag)

    P_coord, P_rho, P_rhoi, P_aux = p_prof
    C_coord, C_rho, C_rhoi, C_aux = c_prof
    if not (C_rho.shape == P_rho.shape and C_rhoi.shape == P_rhoi.shape
            and C_coord.shape == P_coord.shape and C_aux.shape == P_aux.shape):
        r["status"] = "FAIL"
        r["note"] = "profile shape mismatch"
        return r

    # Coordinate grid (physical z coordinate / dimensionless x).
    r["grid_max_abs_diff"] = float(np.max(np.abs(C_coord - P_coord)))
    r["grid_ok"] = bool(np.all(np.abs(C_coord - P_coord) <= EPS_ABS + EPS_REL * np.abs(P_coord)))

    # Final dust density rho_d.
    rho_abs = np.abs(C_rho - P_rho)
    r["rho_max_abs_diff"] = float(np.max(rho_abs))
    r["rho_rel_to_peak"] = float(np.max(rho_abs) / max(np.max(np.abs(P_rho)), 1e-300))
    r["rho_masked_rel_max"] = masked_rel(C_rho, P_rho)
    r["rho_ok"] = bool(np.all(rho_abs <= EPS_ABS + EPS_REL * np.abs(P_rho)))

    # Gas / initial dust array (review B.10): must also match, not ignored.
    rhoi_abs = np.abs(C_rhoi - P_rhoi)
    r["rhoi_max_abs_diff"] = float(np.max(rhoi_abs))
    r["rhoi_ok"] = bool(np.all(rhoi_abs <= EPS_ABS + EPS_REL * np.abs(P_rhoi)))

    # Physical coordinate / aux column (review B.10): must match, not ignored.
    aux_abs = np.abs(C_aux - P_aux)
    r["aux_max_abs_diff"] = float(np.max(aux_abs))
    r["aux_ok"] = bool(np.all(aux_abs <= EPS_ABS + EPS_REL * np.abs(P_aux)))

    # Scalars common to all modes: dt, h, L, t_actual, v_max.
    for key in ("dt", "h", "L", "t_actual", "v_max"):
        _float_gate(key, float(c_diag[key]), float(p_diag[key]), r)

    # Disk-only scalars (advection emits NaN on both sides for these; accept
    # matched NaN but fail on any unexpected/unmatched NaN).
    disk_scalars = ("scale_height", "z_real_max", "real_z_max", "rho_initial_max")
    for key in disk_scalars:
        c_val = float(c_diag[key])
        p_val = float(p_diag[key])
        if is_nan(c_val) and is_nan(p_val):
            r[key + "_ok"] = True
            r[key + "_abs_diff"] = 0.0
        elif is_nan(c_val) != is_nan(p_val):
            r[key + "_ok"] = False
            r[key + "_abs_diff"] = float("nan")
        else:
            _float_gate(key, c_val, p_val, r)

    # t_actual_years may be NaN on both sides (advection has no physical years).
    tay_c = float(c_diag["t_actual_years"])
    tay_p = float(p_diag["t_actual_years"])
    if is_nan(tay_c) and is_nan(tay_p):
        r["t_actual_years_ok"] = True
        r["t_actual_years_abs_diff"] = 0.0
    else:
        _float_gate("t_actual_years", tay_c, tay_p, r)

    # Exact step counts (no tolerance).
    r["Nt_ok"] = int(c_diag["Nt"]) == int(p_diag["Nt"])
    r["steps_ok"] = int(c_diag["actual_steps"]) == int(p_diag["actual_steps"])

    # Mass gate (review B.12): recompute integrated M from the CSV profile and
    # compare the residual against INDEPENDENT integrated masses + emitted flux
    # (never trust the emitted mass_balance_error).  For advection the flux is 0
    # and the initial profile IS the dust (no 0.01*gas factor).
    h = float(c_diag["h"])
    if r["mode"] == "advection":
        c_M_init = float(h * np.sum(C_rhoi))
        c_M_final = float(h * np.sum(C_rho))
    else:
        c_M_init = float(h * np.sum(C_rhoi) * 0.01)   # dust initial = 0.01 * gas
        c_M_final = float(h * np.sum(C_rho))
    r["M_initial_recomp"] = c_M_init
    r["M_final_recomp"] = c_M_final
    residual_indep = (c_M_final - c_M_init) + float(c_diag["cumulative_boundary_flux"])
    r["mass_residual_indep"] = residual_indep
    r["mass_balance_ok"] = abs(residual_indep) <= MASS_TOL * c_M_init
    py_flux = float(p_diag.get("cumulative_boundary_flux", 0.0))
    r["flux_py_vs_cpp_diff"] = abs(float(c_diag["cumulative_boundary_flux"]) - py_flux)
    r["flux_ok"] = fp_close(float(c_diag["cumulative_boundary_flux"]), py_flux)

    # Finite / min over every step.
    r["finite"] = int(c_diag["finite"])
    r["finite_ok"] = int(c_diag["finite"]) == 1 and p_diag["finite"]
    peak = max(abs(float(p_diag["rho_initial_max"])) if not is_nan(p_diag["rho_initial_max"])
               else np.max(np.abs(P_rhoi)), np.max(np.abs(P_rhoi)), 1.0)
    r["min_rho"] = float(c_diag["min_rho"])
    r["min_rho_ok"] = float(c_diag["min_rho"]) >= -1e-14 * peak
    r["min_rho_py"] = float(p_diag["min_rho"])

    # Formation time: absent tf may match NaN; unexpected/unmatched nan/inf fails.
    tf_p = p_diag["t_formation"]
    tf_c = c_diag["t_formation"]
    if is_nan(tf_p) and is_nan(tf_c):
        r["t_formation_ok"] = True
        r["t_formation_abs_diff"] = 0.0
    elif is_nan(tf_p) != is_nan(tf_c):
        r["t_formation_ok"] = False
        r["t_formation_abs_diff"] = float("nan")
    else:
        r["t_formation_abs_diff"] = abs(tf_c - tf_p)
        r["t_formation_ok"] = fp_close(float(tf_c), float(tf_p))

    # alpha in turbulence (EPS_ABS+EPS_REL*abs(P)).
    if r["mode"] == "turbulence":
        _float_gate("alpha", float(c_diag["alpha"]), float(p_diag["alpha"]), r)
    else:
        r["alpha_ok"] = True

    checks = ["grid_ok", "rho_ok", "rhoi_ok", "aux_ok", "dt_ok", "h_ok", "L_ok",
              "t_actual_ok", "v_max_ok", "scale_height_ok", "z_real_max_ok",
              "real_z_max_ok", "rho_initial_max_ok", "t_actual_years_ok",
              "Nt_ok", "steps_ok", "mass_balance_ok", "flux_ok", "finite_ok",
              "min_rho_ok", "t_formation_ok", "alpha_ok"]
    # M_initial/M_final recomp all must be close (strict gate as well).
    if not fp_close(c_diag["mass_initial"], c_M_init):
        checks.append("M_initial_recomp")
    if not fp_close(c_diag["mass_final"], c_M_final):
        checks.append("M_final_recomp")
    failed = [c for c in checks if not r.get(c, True)]
    r["failed_checks"] = sorted(set(failed))
    if failed:
        r["status"] = "FAIL"
        r["note"] = "failed: " + ",".join(sorted(set(failed)))
    return r


# ---------------------------------------------------------------- case tables
# Each case: (mode, cpp_args, py_builder, tag)
# cpp_args is the ordered list of CLI positional arguments (as strings) AFTER
# the mode token and BEFORE the outdir.  This decouples test cases from any
# fragile tag-parsing scheme and directly supports review cases (odd N, D=0,
# fractional time, stop=False continuing past a crossing, etc.).
def _D_for_alpha(R_AU, alpha):
    """Return D = alpha * cs^2 / omega for the fitted representative disk."""
    G = 6.67430e-11
    M_SUN = 1.98847e30
    AU = 1.495978707e11
    Kb = 1.380649e-23
    M_H = 1.6735575e-27
    omega = math.sqrt(G * M_SUN / (R_AU * AU) ** 3)
    temperature = 264.0 * R_AU ** -0.5
    cs = math.sqrt(Kb * temperature / (2.0 * M_H))
    return alpha * cs * cs / omega


def _cases():
    # representative diffusion disk: R=100, alpha=1e-3, Nz=64, zmax in natural
    # units (zmax is expressed in H here via zmax_arg), ~6859 yr, adv=.3 diff=.2
    R100 = 100.0
    D_1e3 = _D_for_alpha(R100, 1e-3)
    D_str = _D_for_alpha(R100, 1e-1)  # stress case: alpha=0.1 makes diffusion bind
    H_in_AU = None  # placeholder; zmax is passed directly as H units for this case
    # zmax argument is in AU; 2*H/AU for R=100 is ~ a few 1e-2 AU. We compute it.
    AU = 1.495978707e11
    G = 6.67430e-11
    M_SUN = 1.98847e30
    Kb = 1.380649e-23
    M_H = 1.6735575e-27
    omega = math.sqrt(G * M_SUN / (R100 * AU) ** 3)
    temperature = 264.0 * R100 ** -0.5
    cs = math.sqrt(Kb * temperature / (2.0 * M_H))
    H = cs / omega
    zmax_2H = 2.0 * H / AU  # zmax in AU giving L=2.0
    zmax_3p5H = 3.5 * H / AU  # zmax in AU giving L=3.5
    return [
        # ---- advection: pointwise + edge cases (review item 7) ----
        ("advection", ["1", "0.01", "0.5", "50"], "adv_1_0.01_0.5_50"),
        ("advection", ["-1", "0.0025", "0.5", "200"], "adv_-1_0.0025_0.5_200"),
        ("advection", ["0", "0.01", "0.5", "50"], "adv_0_0.01_0.5_50"),
        ("advection", ["1", "0.006", "0.5", "100"], "adv_1_0.006_0.5_100"),
        ("advection", ["1", "0.01", "0.001", "100"], "adv_1_0.01_0.001_100"),
        # odd N Gaussian (review item 7)
        ("advection", ["1", "0.01", "0.5", "51"], "adv_odd_51"),
        # abs(CFL) < 1e-13 -> identity (no-op) (review item 7); T small so Nt bounded
        ("advection", ["1", "1e-15", "1e-12", "64"], "adv_cfl_eps"),
        # zero-step: T=0 (review item 7)
        ("advection", ["1", "0.01", "0", "64"], "adv_zero_step"),
        # fractional requested time (review item 7)
        ("advection", ["1", "0.003", "0.007", "80"], "adv_frac"),
        # ---- laminar representative + clamp case (review item 8) ----
        ("laminar", ["1000000", "100", "0.5", "100", "1", "4", "1000", "0"],
         "lam_1000000_100_0.5_100_1_4"),
        # laminar clamp: zmax/H > 4 (R=100, zmax=50AU, Nz=32, short time ~2000yr)
        ("laminar", ["2000", "32", "0.5", "100", "1", "50", "1000", "0"],
         "lam_clamp"),
        # ---- turbulence ----
        # primary stable (adv=.3 diff=.2), D=1e-5, no crossing
        ("turbulence", ["1000000", "100", "0.3", "0.2", "100", "1", "4", "1e-05", "1000", "1", "0"],
         "turb_1000000_100_0.3_0.2_100_1_4_1e-05_1_0"),
        # D=0 -> diffusion dt inf (documented) (review item: D=0)
        ("turbulence", ["10000", "50", "0.3", "0.2", "100", "1", "4", "0", "1000", "1", "0"],
         "turb_10000_50_0.3_0.2_100_1_4_0_1_0"),
        # crossing, stop=True (early stop triggers)
        ("turbulence", ["10000", "400", "0.3", "0.2", "100", "1000", "50", "1e-05", "1000", "1", "1"],
         "turb_cross_stop"),
        # same crossing params, stop=False (must run full Nt) (review item 1/8)
        ("turbulence", ["10000", "400", "0.3", "0.2", "100", "1000", "50", "1e-05", "1000", "1", "0"],
         "turb_cross_continue"),
        # representative diffusion (review item 2): nontrivial stratification,
        # measurable nonzero upper leakage, adv=.3 diff=.2, ~6859 yr
        ("turbulence", ["5000", "64", "0.3", "0.2", "100", "1", g(zmax_2H), g(D_1e3), "1000", "1", "0"],
         "turb_diff_repr"),
        # short bounded diff=.7 stress case (review item 2): EXPECTED_INSTABILITY,
        # inherited negative rho reported separately, never clipped as PASS.
        # alpha=0.1 (D=5.47e13) so the diffusion dt actually binds at d_lim=.7.
        ("turbulence", ["1000", "48", "0.3", "0.7", "100", "1", g(zmax_3p5H), g(D_str), "1000", "1", "0"],
         "turb_diff_stress"),
    ]



# ---------------------------------------------------------------- convergence
def convergence_study():
    """Nx=50,100,200,400 at v=1, CFL=0.5, T=0.5; RMS must decrease, last pair
    observed order in [0.75, 1.25].  Report RMS AND Linf and actual t for both
    Python and C++ per grid (review B.16)."""
    rows = []
    v, cfl, T = 1.0, 0.5, 0.5
    for Nx in [50, 100, 200, 400]:
        dt = cfl * (1.0 / Nx) / abs(v)
        tag = "adv_conv_%03d" % Nx
        p_diag, p_prof = run_python_advection(v, dt, T, Nx)
        t_end = p_diag["t_actual"]
        y_ex = gauss_exact(p_prof[0], v, t_end)
        err_p = p_prof[1] - y_ex
        rms_p = float(np.sqrt(np.mean(err_p ** 2)))
        linf_p = float(np.max(np.abs(err_p)))
        c_diag, c_prof = run_cpp("advection", [g(v), g(dt), g(T), str(Nx)], tag)
        y_ex_c = gauss_exact(c_prof[0], v, c_diag["t_actual"])
        err_c = c_prof[1] - y_ex_c
        rms_c = float(np.sqrt(np.mean(err_c ** 2)))
        linf_c = float(np.max(np.abs(err_c)))
        rows.append({"Nx": Nx, "h": 1.0 / Nx,
                     "rms_python": rms_p, "rms_cpp": rms_c,
                     "linf_python": linf_p, "linf_cpp": linf_c,
                     "t_actual_python": p_diag["t_actual"],
                     "t_actual_cpp": c_diag["t_actual"]})
    strictly_decr = all(rows[i]["rms_cpp"] < rows[i - 1]["rms_cpp"] for i in range(1, len(rows)))
    last_order = math.log(rows[2]["rms_cpp"] / rows[3]["rms_cpp"]) / math.log(2.0)
    return {
        "strictly_decreasing": strictly_decr,
        "last_pair_observed_order": last_order,
        "order_ok": 0.75 <= last_order <= 1.25,
        "rows": rows,
    }


# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description="C++ port acceptance harness")
    ap.add_argument("--binary", default=_BIN_REL, help="path to the C++ CLI binary")
    ap.add_argument("--source-repo", default=SRC, help="path to original src/ dir")
    args = ap.parse_args(argv)
    global _ACTIVE_BIN
    _ACTIVE_BIN = os.path.abspath(args.binary)

    os.makedirs(CASES_DIR, exist_ok=True)
    os.makedirs(os.path.join(ROOT, "validation"), exist_ok=True)

    if not os.path.exists(args.binary):
        print("=" * 70)
        print("TDD pre-build demonstration: C++ binary NOT found at")
        print("  " + args.binary)
        print("Expected: build the solver with g++/CMake first, then rerun.")
        print("=" * 70)
        return 2

    # Defer import until after argparse resolves --source-repo (review B.15).
    resolved_src = _import_sources(args.source_repo)
    print("source repo: %s" % resolved_src)
    print("source identity: %s" % json.dumps(_src_identity(), sort_keys=True))

    cases = _cases()

    # Map a mode + numeric CPP args to a Python reference builder call.
    def py_ref(mode, cpp_args):
        if mode == "advection":
            v, dt, T, Nx = (float(cpp_args[0]), float(cpp_args[1]),
                            float(cpp_args[2]), int(cpp_args[3]))
            return run_python_advection(v, dt, T, Nx)
        if mode == "laminar":
            total_time, Nz, a_lim, R_AU, a_mic, zmax, dd, stop = cpp_args
            return run_python_laminar(float(total_time), int(Nz), float(a_lim),
                                      float(R_AU), float(a_mic), float(zmax),
                                      float(dd), bool(int(stop)))
        if mode == "turbulence":
            total, Nz, a, d, R_AU, amic, zmax, D, dd, vcap, stop = cpp_args
            return run_python_turbulence(float(total), int(Nz), float(a), float(d),
                                         float(R_AU), float(amic), float(zmax),
                                         float(D), float(dd), float(vcap),
                                         bool(int(stop)))
        raise ValueError(mode)

    results = []
    stress_tags = {"turb_diff_stress"}
    for mode, cpp_args, tag in cases:
        try:
            p_diag, p_prof = py_ref(mode, cpp_args)
            # Save the Python reference profile next to the C++ profile (B.10).
            _save_py_profile(tag, p_prof)
            r = check_case(tag, p_diag, p_prof, mode, cpp_args)
        except Exception as e:  # noqa: BLE001
            r = {"tag": tag, "status": "FAIL", "note": "exception: %r" % (e,)}

        if tag in stress_tags:
            # (B.13) EXPECTED_INSTABILITY ONLY when: Python AND C++ show a
            # significant negative min, both are finite throughout, AND every
            # other equality/diagnostic/mass gate passed (failed_checks is ONLY
            # min_rho_ok).  Any exception or any other mismatch is a hard FAIL.
            r["_stress"] = True
            failed_checks = r.get("failed_checks", [])
            min_c = r.get("min_rho")
            min_p = r.get("min_rho_py")
            finite_ok = r.get("finite_ok", False)
            sig_neg_c = min_c is not None and min_c < -1e-10
            sig_neg_p = min_p is not None and min_p < -1e-10
            if (finite_ok and sig_neg_c and sig_neg_p
                    and failed_checks == ["min_rho_ok"]):
                r["status"] = "EXPECTED_INSTABILITY"
                r["note"] = ("stress case: min_rho_cpp=%r min_rho_py=%r (finite both)"
                             % (min_c, min_p))
            else:
                r["status"] = "FAIL"
                r["note"] = (r.get("note", "") + " stress-not-instable-or-mismatch "
                             "min_cpp=%r min_py=%r failed=%r"
                             % (min_c, min_p, failed_checks)).strip()
        results.append(r)
        print("[%s] %s  %s" % (r["status"], r.get("tag", tag), r.get("note", "")))

    conv = convergence_study()
    results.append({
        "tag": "convergence",
        "status": "PASS" if (conv["order_ok"] and conv["strictly_decreasing"]) else "FAIL",
        "note": "strictly_decreasing=%s last_pair_order=%.4f"
                % (conv["strictly_decreasing"], conv["last_pair_observed_order"]),
        "convergence": conv,
    })
    print("[%s] convergence  %s" % (results[-1]["status"], results[-1]["note"]))

    results.append(run_branch_coverage())

    csv_cols = ["tag", "mode", "N", "status", "note", "rho_max_abs_diff", "rho_rel_to_peak",
                "rho_masked_rel_max", "grid_max_abs_diff", "rhoi_max_abs_diff",
                "aux_max_abs_diff", "dt_abs_diff", "t_actual_abs_diff",
                "mass_residual_indep", "min_rho", "t_formation_abs_diff"]
    with open(RESULTS_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=csv_cols, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(r)

    with open(RESULTS_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=str)

    fails = [r for r in results if r.get("status") == "FAIL"]
    stress = [r for r in results if r.get("status") == "EXPECTED_INSTABILITY"]
    print("-" * 70)
    print("SUMMARY: %d cases, %d failed, %d expected-instability stress" %
          (len(results), len(fails), len(stress)))
    for r in fails:
        print("  FAIL %s: %s" % (r.get("tag"), r.get("note", "")))
    for r in stress:
        print("  STRESS(EXPECTED_INSTABILITY) %s: %s" % (r.get("tag"), r.get("note", "")))
    print("CSV: %s" % RESULTS_CSV)
    print("JSON: %s" % RESULTS_JSON)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())





