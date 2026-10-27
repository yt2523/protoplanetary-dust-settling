#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Driver: run the predeclared experiments and save CSV/JSON/PNG artifacts.

Audit-corrected timeline: legacy runs FIRST and its ACTUAL natural endpoint is
taken as the target for positive/refinement runs, so graph/control are aligned.
Plots read saved artifacts (profiles/CSV) rather than rerunning unaligned cases.
"""
import os
import sys
import json
import csv
import math
import datetime
import hashlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import settling_experiments as se

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def _stress_params():
    R = 100.0
    om, cs, H = se._disk(R)
    zmax = 3.5 * H / se.AU
    D = se.D_for_alpha(R, 0.1)
    return dict(total_time=1000, Nz=48, a_lim=0.3, d_lim=0.7, R_AU=R,
                a_micron=1, zmax=zmax, D=D, dust_density=1000, v_cap_factor=1.0)


def _stress_params_physical():
    R = 100.0
    om, cs, H = se._disk(R)
    zmax = 2.0 * H / se.AU
    D = se.D_for_alpha(R, 0.001)
    return dict(total_time=1000, Nz=64, a_lim=0.3, d_lim=0.7, R_AU=R,
                a_micron=1, zmax=zmax, D=D, dust_density=1000, v_cap_factor=1.0)


def _param_sha(p):
    key = "%r|%r|%r|%r|%r|%r|%r|%r|%r" % (
        p["total_time"], p["Nz"], p["a_lim"], p["d_lim"], p["R_AU"],
        p["a_micron"], p["zmax"], p["D"], p["dust_density"])
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


def main():
    ident = se.src_identity()
    meta = {
        "timestamp": datetime.datetime.now().isoformat(),
        "source_path": ident["path"],
        "source_sha256": ident["sha256"],
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    }
    rows = []
    json_bundle = {"meta": meta, "runs": [], "summary": {}}

    stress_params = _stress_params()

    # ---- legacy FIRST: fix its ACTUAL natural endpoint -----------------------
    e_leg = se.TurbulenceExperiment(timestep="legacy", top="vacuum_ghost",
                                    safety=0.9, **stress_params)
    e_leg.run()
    t_target = e_leg.t_actual  # ACTUAL natural endpoint of legacy floor run

    # ---- stress cases on the matched endpoint --------------------------------
    for timestep in ("legacy", "positive"):
        for top in ("vacuum_ghost", "closed_total_flux"):
            if timestep == "legacy" and top == "closed_total_flux":
                continue
            e = se.TurbulenceExperiment(timestep=timestep, top=top, safety=0.9,
                                        target_time_natural=t_target,
                                        **stress_params)
            e.run()
            d = e.diagnostics()
            d["tag"] = "stress_%s_%s" % (timestep, top)
            d["param_sha256"] = _param_sha(stress_params)
            rows.append(d)
            json_bundle["runs"].append(d)

    # ---- save stress profile CSV (paired exact times, dimensionless headers) --
    _save_stress_profiles(stress_params, t_target)

    # ---- time refinement on same matched endpoint ----------------------------
    pro = {}
    ref_prof = None
    for s in (0.9, 0.45, 0.225):
        e = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                    safety=s, target_time_natural=t_target,
                                    **stress_params)
        e.run()
        d = e.diagnostics()
        d["tag"] = "refine_s%.3f" % s
        d["param_sha256"] = _param_sha(stress_params)
        rows.append(d)
        json_bundle["runs"].append(d)
        pro[s] = e.profile()[1]
        if s == 0.225:
            ref_prof = e.profile()[1]

    # ---- time-refinement RMS + Linf -----------------------------------------
    rms_09 = se.rms(pro[0.9], ref_prof)
    rms_045 = se.rms(pro[0.45], ref_prof)
    linf_09 = se.linf(pro[0.9], ref_prof)
    linf_045 = se.linf(pro[0.45], ref_prof)
    json_bundle["summary"]["time_refinement"] = {
        "t_target_natural": t_target,
        "rms_s0.9_vs_s0.225": rms_09,
        "rms_s0.45_vs_s0.225": rms_045,
        "linf_s0.9_vs_s0.225": linf_09,
        "linf_s0.45_vs_s0.225": linf_045,
        "rms_decreases": rms_045 < rms_09,
    }

    # ---- equilibrium + observed orders --------------------------------------
    eq_meta = _run_equilibria(json_bundle, rows)

    # ---- controlled boundary pair (alpha=.001, R100, N64, L2) ---------------
    _run_boundary_pair(json_bundle, rows)

    _write_artifacts(rows, json_bundle)
    _plot_from_artifacts(pro, json_bundle, eq_meta)

    print("artifacts written to", OUT)
    print("  csv:", os.path.join(OUT, "results.csv"))
    print("  json:", os.path.join(OUT, "results.json"))
    print("  source sha256:", json_bundle["meta"]["source_sha256"])
    print("  time_refinement rms_decreases:",
          json_bundle["summary"]["time_refinement"]["rms_decreases"])
    print("  boundary dt_equal:",
          json_bundle["summary"]["boundary_pair"]["dt_equal"])
    return 0



def _save_stress_profiles(params, t_target):
    for timestep in ("legacy", "positive"):
        for top in ("vacuum_ghost", "closed_total_flux"):
            if timestep == "legacy" and top == "closed_total_flux":
                continue
            e = se.TurbulenceExperiment(timestep=timestep, top=top, safety=0.9,
                                        target_time_natural=t_target, **params)
            e.run()
            z, rho, _ = e.profile()
            key = "%s_%s" % (timestep, top)
            np.savetxt(os.path.join(OUT, "stress_profile_%s.csv" % key),
                       np.column_stack([z, rho]),
                       header="z_H,rho_dust_natural", delimiter=",",
                       comments="")
    with open(os.path.join(OUT, "stress_paired_times.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "t_actual_natural", "t_actual_years", "dt", "Nt"])
        for timestep in ("legacy", "positive"):
            for top in ("vacuum_ghost", "closed_total_flux"):
                if timestep == "legacy" and top == "closed_total_flux":
                    continue
                e = se.TurbulenceExperiment(timestep=timestep, top=top, safety=0.9,
                                            target_time_natural=t_target, **params)
                e.run()
                d = e.diagnostics()
                w.writerow(["%s_%s" % (timestep, top), d["t_actual"],
                            d["t_actual_years"], d["dt"], d["Nt"]])


def _run_equilibria(json_bundle, rows):
    rms_vals = []
    linf_vals = []
    Ns = [32, 64, 128]
    peak_c = None
    for N in Ns:
        sim = se.make_equilibrium_sim(N, L=2.0, S=0.001, alpha=0.01)
        z = sim.z
        g = sim.rho_initial
        q = se.discrete_closed_equilibrium(sim, 0.001)
        rho_d = q * g
        rho_d = rho_d / (sim.h * np.sum(rho_d))
        rho_c = se.analytic_closed_equilibrium(z, g, 0.001, 0.01)
        rho_c = rho_c / (sim.h * np.sum(rho_c))
        rms_vals.append(se.rms(rho_d, rho_c))
        linf_vals.append(se.linf(rho_d, rho_c))
        if N == max(Ns):
            peak_c = np.max(np.abs(rho_c))
        rows.append({
            "tag": "equilibrium_N%03d" % N, "N": N, "L": float(sim.L),
            "h": float(sim.h), "rms_disc_cont": rms_vals[-1],
            "linf_disc_cont": linf_vals[-1],
            "mass_discrete": float(sim.h * np.sum(rho_d)),
            "mass_continuous": float(sim.h * np.sum(rho_c)),
        })
        np.savetxt(os.path.join(OUT, "eq_profile_N%03d.csv" % N),
                   np.column_stack([z, rho_d, rho_c]),
                   header="z,rho_discrete,rho_continuous", delimiter=",",
                   comments="")
    orders_rms = -np.log2(np.array(rms_vals[1:]) / np.array(rms_vals[:-1]))
    orders_linf = -np.log2(np.array(linf_vals[1:]) / np.array(linf_vals[:-1]))
    eq_meta = {
        "rms": rms_vals, "linf": linf_vals, "orders_rms": orders_rms.tolist(),
        "orders_linf": orders_linf.tolist(),
        "finest_normalized_rms_over_peak": rms_vals[-1] / peak_c,
    }
    json_bundle["summary"]["equilibrium"] = eq_meta
    return eq_meta



def _run_boundary_pair(json_bundle, rows):
    p = _stress_params_physical()
    e_leg = se.TurbulenceExperiment(timestep="legacy", top="vacuum_ghost", **p)
    e_leg.run()
    t_target = e_leg.t_actual
    ev = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                 safety=0.9, target_time_natural=t_target, **p)
    ec = se.TurbulenceExperiment(timestep="positive", top="closed_total_flux",
                                 safety=0.9, target_time_natural=t_target, **p)
    ev._setup(); ev._choose_step()
    ec._setup(); ec._choose_step()
    dt_common = min(ev.bound_dt, ec.bound_dt)
    Nt_common = int(math.ceil(t_target / dt_common))
    dt_common = t_target / Nt_common
    ev2 = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                  safety=0.9, target_time_natural=t_target,
                                  dt=dt_common, **p)
    ec2 = se.TurbulenceExperiment(timestep="positive", top="closed_total_flux",
                                  safety=0.9, target_time_natural=t_target,
                                  dt=dt_common, **p)
    ev2.run(); ec2.run()
    dv = ev2.diagnostics(); dc = ec2.diagnostics()
    dv["tag"] = "boundary_vacuum"; dv["param_sha256"] = _param_sha(p)
    dc["tag"] = "boundary_closed"; dc["param_sha256"] = _param_sha(p)
    rows.append(dv); rows.append(dc)
    json_bundle["runs"].append(dv); json_bundle["runs"].append(dc)
    json_bundle["summary"]["boundary_pair"] = {
        "alpha": float(ev2.sim.alpha), "N": p["Nz"],
        "common_dt": dt_common, "t_target_natural": t_target,
        "dt_equal": abs(dv["dt"] - dc["dt"]),
        "t_actual_equal": abs(dv["t_actual"] - dc["t_actual"]),
        "closed_mass_residual": dc["mass_residual"],
        "vacuum_cum_top_flux": dv["cum_top_flux"],
    }
    zv, rv, gv = ev2.profile()
    zc, rc, gc = ec2.profile()
    np.savetxt(os.path.join(OUT, "boundary_pair_profiles.csv"),
               np.column_stack([zv, rv, rc]),
               header="z_H,rho_vacuum,rho_closed", delimiter=",", comments="")


def _write_artifacts(rows, json_bundle):
    csv_path = os.path.join(OUT, "results.csv")
    if rows:
        keys = sorted(set().union(*(r.keys() for r in rows)))
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    json_path = os.path.join(OUT, "results.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_bundle, f, indent=2, default=str)


def _plot_from_artifacts(pro, json_bundle, eq_meta):
    # stress plot (legacy vs positive, vacuum)
    stress_params = _stress_params()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    leg = np.genfromtxt(os.path.join(OUT, "stress_profile_legacy_vacuum_ghost.csv"),
                        delimiter=",", skip_header=1)
    pos = np.genfromtxt(os.path.join(OUT, "stress_profile_positive_vacuum_ghost.csv"),
                        delimiter=",", skip_header=1)
    ax.plot(leg[:, 0], leg[:, 1], label="legacy (vacuum)", lw=1.5)
    ax.plot(pos[:, 0], pos[:, 1], label="positive (vacuum)", lw=1.5)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("z (H)"); ax.set_ylabel("dust density (natural)")
    ax.set_title("Stress case: legacy vs positive timestep")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stress_legacy_vs_positive.png"), dpi=130)
    plt.close(fig)

    # time refinement plot from saved profiles
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ref = pro[0.225]
    _prof_path = os.path.join(OUT, "stress_profile_positive_vacuum_ghost.csv")
    _z_over_H = np.loadtxt(_prof_path, delimiter=",", skiprows=1)[:, 0]
    assert _z_over_H.size == ref.size
    _pro = {0.9: pro[0.9], 0.45: pro[0.45], 0.225: pro[0.225]}
    np.savetxt(os.path.join(OUT, "time_refinement_profiles.csv"),
               np.column_stack([_z_over_H, _pro[0.9], _pro[0.45], _pro[0.225]]),
               delimiter=",", header="z_over_H,rho_safety_0p9,rho_safety_0p45,rho_safety_0p225", comments="")
    _tab = np.loadtxt(os.path.join(OUT, "time_refinement_profiles.csv"), delimiter=",", skiprows=1)
    z = _tab[:, 0]
    pro = {0.9: _tab[:, 1], 0.45: _tab[:, 2], 0.225: _tab[:, 3]}
    ax.plot(z, pro[0.9], label="safety=0.9", lw=1.5)
    ax.plot(z, pro[0.45], label="safety=0.45", lw=1.5)
    ax.plot(z, pro[0.225], label="safety=0.225 (ref)", lw=1.5)
    ax.set_xlabel("z (H)"); ax.set_ylabel("dust density (natural)")
    ax.set_title("Positive timestep: time refinement")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(OUT, "time_refinement.png"), dpi=130)
    plt.close(fig)

    # equilibrium convergence plot from saved CSV
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for N in (32, 64, 128):
        arr = np.genfromtxt(os.path.join(OUT, "eq_profile_N%03d.csv" % N),
                            delimiter=",", skip_header=1)
        ax.plot(arr[:, 0], arr[:, 1], label="discrete N=%d" % N, lw=1.5)
        ax.plot(arr[:, 0], arr[:, 2], label="continuous N=%d" % N,
                ls="--", lw=1.2)
    ax.set_xlabel("z"); ax.set_ylabel("normalized rho")
    ax.set_title("Discrete vs continuous closed equilibrium")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(os.path.join(OUT, "equilibrium_convergence.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

