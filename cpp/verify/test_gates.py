#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Focused unit tests for the review-fixed acceptance gates (no simulation).

These prove, in isolation, that a deliberately-mutated C++ diagnostic flips the
``verify.check_case`` gates to FAIL — the exact defect fixed in review #1/#2
(the old code read the recomputed mass as a truthy NUMBER, so a mismatch could
never be detected).  ``run_cpp`` is monkey-patched so no PyBullet/CMake binary or
WSL call is needed; the comparison logic is exercised directly.

Run with pytest or plain python.
"""

import os
import sys

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)
import verify  # noqa: E402


def _base_p_diag(mode):
    return {
        "mode": mode, "N": 50, "L": 1.0, "h": 0.02, "dt": 0.01,
        "t_actual": 0.5, "actual_steps": 50, "Nt": 50,
        "t_formation": float("nan"), "t_formation_years": float("nan"),
        "t_actual_years": float("nan"), "mass_initial": 1.0, "mass_final": 1.0,
        "min_rho": 1e-12, "max_rho": 1.0, "finite": True,
        "cumulative_boundary_flux": 0.0, "v_max": 1.0, "alpha": float("nan"),
        "z_real_max": float("nan"), "real_z_max": float("nan"),
        "rho_initial_max": 1.0, "scale_height": float("nan"),
    }


def _make_cpp_diag(p_diag, **overrides):
    d = dict(p_diag)
    d["finite"] = int(d["finite"])
    d.update(overrides)
    return d


def _make_prof(N=50):
    x = np.linspace(0.0, 1.0, N, endpoint=False)
    rhoi = np.exp(-0.5 * ((x - 0.5) / 0.07) ** 2)
    rhoi /= rhoi.max()
    return x, rhoi.copy(), rhoi.copy(), x.copy()


def _run_check(monkeypatch, mode, cpp_overrides, p_overrides=None):
    """Monkeypatch run_cpp to emit a (possibly mutated) C++ result, then check_case."""
    p_diag = _base_p_diag(mode)
    if p_overrides:
        p_diag.update(p_overrides)
    p_prof = _make_prof()
    coord, rho, rhoi, aux = p_prof
    h = 0.02
    # Recompute the masses exactly as check_case does (advection: dust==initial;
    # disks: dust = 0.01 * gas), so a NON-mutated C++ result passes the gate and
    # only an explicit mutation flips it to FAIL.
    if mode == "advection":
        m_init = h * float(np.sum(rhoi))
        m_final = h * float(np.sum(rho))
    else:
        m_init = h * float(np.sum(rhoi)) * 0.01
        m_final = h * float(np.sum(rho))
    p_diag["mass_initial"] = m_init
    p_diag["mass_final"] = m_final
    c_diag = _make_cpp_diag(p_diag, **cpp_overrides)

    def fake_run_cpp(m, params, tag, bin_win=None):
        return c_diag, p_prof

    monkeypatch.setattr(verify, "run_cpp", fake_run_cpp)
    return verify.check_case("mut", p_diag, p_prof, mode, ["1", "0.01", "0.5", "50"])


def test_emitted_mass_mutation_fails(monkeypatch):
    """A deliberately-increased emitted mass_final must FAIL (review fix #1).

    The old code appended the recomputed NUMERIC value into the checks list and
    then r.get(name, True) returned that truthy number, so it could never fail.
    """
    r = _run_check(monkeypatch, "advection",
                   cpp_overrides={"mass_final": 1.0 + 1e-4})
    assert r["emitted_mass_final_ok"] is False
    assert r["status"] == "FAIL"
    assert "emitted_mass_final_ok" in r["failed_checks"]


def test_emitted_mass_match_passes(monkeypatch):
    r = _run_check(monkeypatch, "advection", cpp_overrides={})
    assert r["emitted_mass_initial_ok"] is True
    assert r["emitted_mass_final_ok"] is True


def test_min_rho_mismatch_fails(monkeypatch):
    """A monitored-minimum mismatch between Python and C++ must FAIL (fix #2)."""
    r = _run_check(monkeypatch, "advection",
                   cpp_overrides={"min_rho": -0.01},
                   p_overrides={"min_rho": 0.0})
    assert r["min_rho_match_ok"] is False
    assert r["status"] == "FAIL"
    assert "min_rho_match_ok" in r["failed_checks"]


def test_min_rho_dust_peak_gate(monkeypatch):
    """min_rho gate uses the initial DUST peak (0.01*gas), no max(...,1) floor.

    For a disk (laminar), the gas peak is ~1 so the dust peak is ~0.01 and the
    tolerance is 1e-14 * 0.01 = 1e-16.  A minimum of -5e-16 is BELOW that tighter
    bound and must FAIL, which proves the old loose 1.0-floor (1e-14 tolerance)
    was removed: under the old gate -5e-16 would have silently passed.
    """
    r = _run_check(monkeypatch, "laminar",
                   cpp_overrides={"min_rho": -5e-16},
                   p_overrides={"min_rho": -5e-16})
    assert abs(r["dust_peak"] - 0.01) < 1e-6  # 0.01 * gas_peak(=1.0)
    assert r["min_rho_ok"] is False  # -5e-16 < -1e-16 (tighter than old floor)
    assert r["min_rho_match_ok"] is True  # both sides agree


def test_min_rho_within_dust_peak_passes(monkeypatch):
    """A value within the dust-peak bound must PASS (e.g. -5e-17 > -1e-16)."""
    r = _run_check(monkeypatch, "laminar",
                   cpp_overrides={"min_rho": -5e-17},
                   p_overrides={"min_rho": -5e-17})
    assert r["min_rho_ok"] is True


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))
