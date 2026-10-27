#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the independent settling experiments (literature-driven round 1).

Groups map to research/PLAN.md predeclared metrics; no full-physics scans and
no tolerance loosening.  This module imports the original solver indirectly via
settling_experiments (which calls the original upwind_method), and never starts
any GUI / long relaxation: equilibrium uses the exact discrete recurrence.

Audit fixes incorporated (see REVIEW.md): meaningful stratified matrix fixtures,
matched ACTUAL legacy endpoint for positive/refinement, real serialized validity
diagnostics, source identity regression, controlled boundary pair.
"""
import os
import math
import shutil
import hashlib
import tempfile

import numpy as np
import pytest

import settling_experiments as se


# --------------------------------------------------------------------------- #
# shared builders
# --------------------------------------------------------------------------- #
def _stress_params():
    """Baseline turb_diff_stress parameters (verify.py): N=48, a=.3, d=.7,
    R=100, a_micron=1, zmax=3.5H (AU), D=D_for_alpha(100,0.1), dd=1000, vcap=1.
    """
    R = 100.0
    om, cs, H = se._disk(R)
    zmax = 3.5 * H / se.AU
    D = se.D_for_alpha(R, 0.1)
    return dict(total_time=1000, Nz=48, a_lim=0.3, d_lim=0.7, R_AU=R,
                a_micron=1, zmax=zmax, D=D, dust_density=1000, v_cap_factor=1.0)


def _stress_params_physical():
    """Representative physical input for boundary isolation: alpha=.001, R=100,
    N=64, L=2 (matching PLAN boundary-sensitivity prescription)."""
    R = 100.0
    om, cs, H = se._disk(R)
    zmax = 2.0 * H / se.AU  # 2 scale heights -> L=2 natural units
    D = se.D_for_alpha(R, 0.001)
    return dict(total_time=1000, Nz=64, a_lim=0.3, d_lim=0.7, R_AU=R,
                a_micron=1, zmax=zmax, D=D, dust_density=1000, v_cap_factor=1.0)


def _make_sim(Nz, a_lim=0.3, d_lim=0.7, R=100.0, a_micron=10.0, D=1e-5,
              zmax_hu=3.5, dd=1000, vcap=1.0):
    om, cs, H = se._disk(R)
    zmax = zmax_hu * H / se.AU
    cls = se._import_source()
    sim = cls(1000.0, Nz, a_lim, d_lim, R, a_micron, zmax, D, dd)
    sim.v_cap_factor = vcap
    sim.calculate_natural_unit()
    sim.initial_dust_points()
    return sim


def _make_matrix_sim(N, S=0.05, alpha=0.01):
    """Stratified fixture for matrix tests: L=2, nonflat gas g=exp(-z^2/2),
    nonzero settling v=-S*z/g, alpha=.01.  Diffusion contribution far above the
    1e-12 abs gate (unlike the D=1e-5 -> alpha~1e-20 fixture)."""
    sim = se.make_equilibrium_sim(N, L=2.0, S=S, alpha=alpha)
    sim.alpha = float(alpha)
    sim.v = -S * sim.z / sim.rho_initial
    return sim


def _absrel_ok(got, ref):
    return abs(got - ref) <= se.EPS_ABS + se.EPS_REL * abs(ref)


# --------------------------------------------------------------------------- #
# Audit 1: source identity (first-call SHA valid; switch resolves requested file)
# --------------------------------------------------------------------------- #
def test_source_identity_first_call_valid():
    ident = se.src_identity()
    assert ident["path"] is not None
    assert ident["sha256"] is not None and len(ident["sha256"]) == 64
    assert ident["path"].endswith("dust_solver_turbulence.py")
    assert os.path.exists(ident["path"])


def test_source_identity_switches_to_temp_copy():
    ident0 = se.src_identity()
    src_file = ident0["path"]
    assert src_file and os.path.exists(src_file)
    data = open(src_file, "rb").read()
    tmpdir = tempfile.mkdtemp(prefix="dust_src_audit_")
    try:
        cp_path = os.path.join(tmpdir, "dust_solver_turbulence.py")
        with open(cp_path, "wb") as f:
            f.write(data)
        cls = se._import_source(tmpdir)
        ident1 = se.src_identity()
        assert ident1["path"] == cp_path
        saved = open(cp_path, "rb").read()
        assert ident1["sha256"] == hashlib.sha256(saved).hexdigest()
        assert cls is se._TURBULENCE_CLS
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# --------------------------------------------------------------------------- #
# Group 1: positivity diagonal / adjacency (matrix) checks
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("Nz", [8, 10, 16])
@pytest.mark.parametrize("top", ["vacuum_ghost", "closed_total_flux"])
def test_matrix_diagonal_matches_derived(Nz, top):
    sim = _make_matrix_sim(Nz, S=0.05, alpha=0.01)
    assert (np.max(sim.rho_initial) / np.min(sim.rho_initial)) > 2.0
    sim.choose_time_step()
    if top == "closed_total_flux":
        A, B = se.closed_operator_matrix(sim)
    else:
        A, B = se.operator_matrix(sim)
    lam = se.derived_diagonal_lambda(sim, top)
    diag = np.diag(B)
    assert diag.shape == lam.shape
    for i in range(Nz):
        assert _absrel_ok(diag[i], lam[i]), (i, diag[i], lam[i])
    # diffusion contribution genuinely above tolerance gate
    diffusive = sim.alpha / (sim.h ** 2)
    assert diffusive > 10.0 * se.EPS_ABS


@pytest.mark.parametrize("Nz", [8, 10, 16])
@pytest.mark.parametrize("top", ["vacuum_ghost", "closed_total_flux"])
def test_matrix_offdiag_nonnegative_gain(Nz, top):
    sim = _make_matrix_sim(Nz, S=0.05, alpha=0.01)
    sim.choose_time_step()
    if top == "closed_total_flux":
        A, B = se.closed_operator_matrix(sim)
    else:
        A, B = se.operator_matrix(sim)
    G = se.derived_gain_matrix(sim, top)
    for i in range(Nz):
        for j in range(Nz):
            if i == j:
                continue
            assert _absrel_ok(B[i, j], -G[i, j]), (i, j, B[i, j], G[i, j])
            assert B[i, j] <= 0.0 + 1e-12
    assert np.all(np.diag(B) >= -1e-12)


@pytest.mark.parametrize("top", ["vacuum_ghost", "closed_total_flux"])
def test_matrix_dt_above_bound_gives_negative_diagonal(top):
    sim = _make_matrix_sim(8, S=0.05, alpha=0.01)
    sim.choose_time_step()
    lam = se.derived_diagonal_lambda(sim, top)
    lam_max = np.max(lam)
    sim.dt = 1.01 / lam_max
    if top == "closed_total_flux":
        A = se.closed_update_matrix(sim)
    else:
        A = se.update_matrix(sim)
    assert np.any(np.diag(A) < 0.0)



# --------------------------------------------------------------------------- #
# Group 2: known stress case (legacy reproduces negative, positive stays finite)
# --------------------------------------------------------------------------- #
def _run_stress(timestep, top="vacuum_ghost", safety=0.9, target=None):
    p = _stress_params()
    e = se.TurbulenceExperiment(timestep=timestep, top=top, safety=safety,
                                target_time_natural=target, **p)
    e.run()
    return e


def test_stress_legacy_reproduces_negative():
    e = _run_stress("legacy", "vacuum_ghost")
    d = e.diagnostics()
    assert d["finite_all"] is True
    assert d["min_rho_seen"] < -1e-3
    assert abs(d["mass_residual"]) <= se.MASS_TOL * d["mass_initial"]


def test_stress_positive_finite_nonnegative():
    # run legacy first to get its ACTUAL natural endpoint, then match it
    e_leg = _run_stress("legacy", "vacuum_ghost")
    t_target = e_leg.t_actual
    e = _run_stress("positive", "vacuum_ghost", safety=0.9, target=t_target)
    d = e.diagnostics()
    assert d["finite_all"] is True
    # rho_initial_max is ALREADY the dust peak (0.01); do NOT rescale again
    init_peak = d["rho_initial_max"]
    assert d["min_rho_seen"] >= -1e-14 * init_peak
    assert d["dt"] * d["lambda_max"] <= 0.9 * (1 + 1e-12)
    assert abs(d["mass_residual"]) <= se.MASS_TOL * d["mass_initial"]
    # actual end time matches legacy (matched natural endpoint)
    assert abs(d["t_actual"] - t_target) < 1e-10 * max(t_target, 1.0)


# --------------------------------------------------------------------------- #
# Group 3: stable time refinement (positive safety .9/.45/.225, same end time)
# --------------------------------------------------------------------------- #
def test_time_refinement_rms_decreases():
    e_leg = _run_stress("legacy", "vacuum_ghost")
    t_target = e_leg.t_actual
    pro = {}
    for s in (0.9, 0.45, 0.225):
        e = _run_stress("positive", "vacuum_ghost", safety=s, target=t_target)
        pro[s] = e.profile()[1]
    ref = pro[0.225]
    rms_09 = se.rms(pro[0.9], ref)
    rms_045 = se.rms(pro[0.45], ref)
    assert rms_045 < rms_09
    for s in (0.9, 0.45, 0.225):
        assert np.all(np.isfinite(pro[s]))



# --------------------------------------------------------------------------- #
# Group 4: diffusion isolation (PLAN L2, N32, alpha=.001, v=0, q=.01)
# --------------------------------------------------------------------------- #
def _diffusion_sim(Nz, alpha=0.001):
    sim = se.make_equilibrium_sim(Nz, L=2.0, S=0.0, alpha=alpha)
    sim.alpha = float(alpha)
    sim.v = np.zeros_like(sim.v)      # pure diffusion
    sim.rho = 0.01 * sim.rho_initial.copy()  # q = 0.01
    lam = se.derived_diagonal_lambda(sim, "vacuum_ghost")
    sim.dt = 0.9 / np.max(lam)
    sim.Nt = 1
    return sim


def test_diffusion_closed_preserves_constant_ratio():
    sim = _diffusion_sim(32)
    lam = se.derived_diagonal_lambda(sim, "closed_total_flux")
    dt = 0.9 / np.max(lam)
    sim.dt = dt  # explicit assignment so top-flux compensation uses same dt
    rho = np.asarray(sim.rho, dtype=float).copy()
    rho0 = rho.copy()
    mass0 = sim.h * np.sum(rho0)
    for _ in range(100):
        ft = se.top_flux(sim, rho)
        out = sim.upwind_method(rho)
        out = out.copy()
        out[-1] += (sim.dt / sim.h) * ft  # closed, same sim.dt
        rho = out
    for i in range(rho.size):
        assert _absrel_ok(float(rho[i]), float(rho0[i]))
    assert abs(sim.h * np.sum(rho) - mass0) <= se.MASS_TOL * mass0


def test_diffusion_vacuum_single_step_loss_matches():
    sim = _diffusion_sim(32)
    lam = se.derived_diagonal_lambda(sim, "vacuum_ghost")
    dt = 0.9 / np.max(lam)
    sim.dt = dt
    rho = np.asarray(sim.rho, dtype=float).copy()
    out = sim.upwind_method(rho)  # vacuum (no correction)
    expected_loss = dt * sim.alpha * rho[-1] / (2.0 * sim.h)
    actual_loss = sim.h * np.sum(rho - out)
    assert expected_loss > 0.0
    assert abs(actual_loss - expected_loss) <= se.EPS_ABS + se.EPS_REL * expected_loss
    assert not np.isclose(sim.alpha, 0.0)


def test_diffusion_weighted_not_plain_laplace():
    # constant ratio -> zero weighted diffusion; ordinary Laplacian on rho
    # gives a change far above abs+rel tolerance (wrong update demonstrably
    # FAILS stationarity), so lap>0 alone is insufficient as a "pass".
    sim = _diffusion_sim(32)
    rho = np.asarray(sim.rho, dtype=float).copy()  # = 0.01*g
    out = sim.upwind_method(rho)
    assert np.allclose(out[:-1], rho[:-1], atol=1e-12, rtol=1e-12)
    lap = np.gradient(np.gradient(rho, sim.h), sim.h)
    peak = float(np.max(np.abs(rho)))
    assert np.max(np.abs(lap)) > 1e-8 * peak  # ordinary-Laplace foil change



# --------------------------------------------------------------------------- #
# Group 5: continuous / discrete equilibrium
# --------------------------------------------------------------------------- #
def _eq_sim(N):
    return se.make_equilibrium_sim(N, L=2.0, S=0.001, alpha=0.01)


def test_equilibrium_discrete_stationary():
    sim = _eq_sim(64)
    q = se.discrete_closed_equilibrium(sim, 0.001)
    rho = q * sim.rho_initial
    rho = rho / (sim.h * np.sum(rho))
    ft = se.top_flux(sim, rho)
    out = sim.upwind_method(rho)
    out = out.copy()
    out[-1] += (sim.dt / sim.h) * ft
    assert np.max(np.abs(out - rho)) <= 1e-10
    assert abs(sim.h * np.sum(out) - sim.h * np.sum(rho)) <= 1e-13


def test_equilibrium_converges_to_continuous():
    rms_vals = []
    linf_vals = []
    Ns = [32, 64, 128]
    peak_c = None
    for N in Ns:
        sim = _eq_sim(N)
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
    assert rms_vals[1] < rms_vals[0]
    assert rms_vals[2] < rms_vals[1]
    assert linf_vals[1] < linf_vals[0]
    assert linf_vals[2] < linf_vals[1]
    assert rms_vals[-1] / peak_c < 0.02


def test_equilibrium_evolve_stays_stationary():
    sim = _eq_sim(64)
    q = se.discrete_closed_equilibrium(sim, 0.001)
    rho = q * sim.rho_initial
    rho = rho / (sim.h * np.sum(rho))
    rho0 = rho.copy()
    mass0 = sim.h * np.sum(rho0)
    lam = se.derived_diagonal_lambda(sim, "closed_total_flux")
    dt = 0.9 / np.max(lam)
    for _ in range(200):
        ft = se.top_flux(sim, rho)
        out = sim.upwind_method(rho)
        out = out.copy()
        out[-1] += (dt / sim.h) * ft
        rho = out
    for i in range(rho.size):
        assert _absrel_ok(float(rho[i]), float(rho0[i]))
    assert abs(sim.h * np.sum(rho) - mass0) <= se.MASS_TOL * mass0



# --------------------------------------------------------------------------- #
# Group 6: physical applicability diagnostics (real serialized validity)
# --------------------------------------------------------------------------- #
def test_diagnostics_record_ts_mach_cap():
    p = _stress_params()
    e = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                safety=0.9, **p)
    e.run()
    d = e.diagnostics()
    sim = e.sim
    ts_real = np.sqrt(np.pi / 8.0) * sim.dust_density * sim.real_dust_a \
        / (sim.real_rho * sim.sound_speed)
    tau_s = sim.omega * ts_real
    raw_v_real = -sim.omega ** 2 * sim.z_real * ts_real
    mach_raw = np.abs(raw_v_real) / sim.sound_speed
    vclipped_real = sim.v * sim.omega * sim.scale_height
    mach_used = np.abs(vclipped_real) / sim.sound_speed
    assert d["tau_s_min"] == pytest.approx(np.min(tau_s), rel=1e-9, abs=0.0)
    assert d["tau_s_max"] == pytest.approx(np.max(tau_s), rel=1e-9, abs=0.0)
    assert d["mach_raw_max"] == pytest.approx(np.max(mach_raw), rel=1e-9)
    assert d["mach_used_max"] == pytest.approx(np.max(mach_used), rel=1e-9)
    v_cap = float(sim.v_cap_factor) * float(sim.sound_speed)
    cap_ref = int(np.sum(np.abs(raw_v_real) > v_cap))
    assert d["cap_count"] == cap_ref
    assert d["cap_count"] >= 0
    assert d["source_sha256"] == se.src_identity()["sha256"]
    assert d["source_path"] == se.src_identity()["path"]
    assert len(d["param_sha256"]) == 64
    assert d["finite_all"] is True


def test_diagnostics_static_cap_active_configuration():
    R = 100.0
    om, cs, H = se._disk(R)
    zmax = 3.5 * H / se.AU
    p = dict(total_time=10, Nz=48, a_lim=0.3, d_lim=0.7, R_AU=R,
             a_micron=300.0, zmax=zmax, D=se.D_for_alpha(R, 0.1),
             dust_density=1000, v_cap_factor=1.0)
    e = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                safety=0.9, **p)
    e.run()
    d = e.diagnostics()
    assert d["cap_count"] > 0


def test_diagnostics_valid_nocap_configuration():
    R = 100.0
    om, cs, H = se._disk(R)
    zmax = 3.5 * H / se.AU
    p = dict(total_time=10, Nz=48, a_lim=0.3, d_lim=0.7, R_AU=R,
             a_micron=0.01, zmax=zmax, D=se.D_for_alpha(R, 0.1),
             dust_density=1000, v_cap_factor=1.0)
    e = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                safety=0.9, **p)
    e.run()
    d = e.diagnostics()
    assert d["cap_count"] == 0
    assert d["finite_all"] is True


# --------------------------------------------------------------------------- #
# Audit 6: controlled boundary case (alpha=.001, R100, N64, L2 representative)
# --------------------------------------------------------------------------- #
def test_boundary_isolation_common_time_dt():
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
    ev2.run()
    ec2.run()
    dv = ev2.diagnostics()
    dc = ec2.diagnostics()
    assert dv["N"] == dc["N"]
    assert dv["h"] == dc["h"]
    assert dv["mass_initial"] == pytest.approx(dc["mass_initial"], rel=1e-12)
    assert dv["dt"] == pytest.approx(dc["dt"], rel=1e-12)
    assert dv["t_actual"] == pytest.approx(dc["t_actual"], rel=1e-12)
    assert dv["t_actual"] == pytest.approx(t_target, rel=1e-10)
    assert abs(dc["mass_residual"]) <= se.MASS_TOL * dc["mass_initial"]
    assert dv["cum_top_flux"] > 0.0
    assert pytest.approx(dv["alpha"], rel=1e-9) != 0.1



# --------------------------------------------------------------------------- #
# meaningful unit tests (cheap, no simulation)
# --------------------------------------------------------------------------- #
def test_top_flux_ghost_formula():
    sim = _make_sim(16)
    rho = np.asarray(sim.rho, dtype=float)
    rho[-1] = 2.5
    ft = se.top_flux(sim, rho)
    expected = max(sim.v[-1], 0.0) * rho[-1] + sim.alpha * rho[-1] / (2.0 * sim.h)
    assert abs(ft - expected) < 1e-15


def test_positive_step_time_alignment():
    p = _stress_params()
    e = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                safety=0.9, **p)
    e._setup()
    e._choose_step()
    assert abs(e.dt * e.Nt - e.sim.T) < 1e-12 * e.sim.T
    assert e.bound_dt <= 0.9 / e.lambda_max + 1e-12


def test_legacy_floor_termination_untouched():
    p = _stress_params()
    e = se.TurbulenceExperiment(timestep="legacy", top="vacuum_ghost", **p)
    e._setup()
    e._choose_step()
    assert e.Nt == int(e.sim.T / e.sim.dt)
    assert abs(e.dt - e.sim.dt) < 1e-15


def test_safeguard_positive_mode_validation():
    with pytest.raises(AssertionError):
        se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                safety=0.0, **_stress_params())
    with pytest.raises(AssertionError):
        se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                safety=1.5, **_stress_params())


def test_safeguard_nt_bound():
    p = _stress_params()
    e = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                safety=0.9, **p)
    e._setup()
    e._choose_step()
    assert e.Nt <= 100000


def test_stress_closed_vs_vacuum_mass_conservation():
    p = _stress_params()
    e_leg = se.TurbulenceExperiment(timestep="legacy", top="vacuum_ghost", **p)
    e_leg.run()
    t_target = e_leg.t_actual
    ev = se.TurbulenceExperiment(timestep="positive", top="vacuum_ghost",
                                 safety=0.9, target_time_natural=t_target, **p)
    ev.run()
    ec = se.TurbulenceExperiment(timestep="positive", top="closed_total_flux",
                                 safety=0.9, target_time_natural=t_target, **p)
    ec.run()
    dv = ev.diagnostics()
    dc = ec.diagnostics()
    assert abs(dc["mass_residual"]) <= se.MASS_TOL * dc["mass_initial"]
    assert dv["cum_top_flux"] > 0.0
    assert abs(dv["t_actual"] - dc["t_actual"]) < 1e-10 * max(dv["t_actual"], 1.0)

