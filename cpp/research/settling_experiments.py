#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent Python settling experiments (literature-driven first round).

Inherits / calls the ORIGINAL DustSolver_Turbulence.upwind_method (read-only
source checkout) so the weighted gas-density diffusion stencil is NOT
reimplemented.  Only two explicit experiment switches are added:

  - timestep : "legacy"   -> original choose_time_step() min(adv,diff) + floor
                             termination (baseline, faithful).
               "positive" -> joint positivity-safe step bound from the update
                             matrix diagonal with safety factor, then
                             ceil(T/dt_bound) steps of dt=T/Nt (exact alignment).
  - top      : "vacuum_ghost" -> original upper boundary (ghost ratio=0,gas=0).
               "closed_total_flux" -> explicit zero total flux at the top face,
                             implemented by adding back the original top-face
                             loss to the last cell AFTER the original update.

No clip of negative densities anywhere; no Sc modification; no laminar
grid change.  Reference source:
  D:/Codex/tmp/phd-github-audit-2026-10-26/protoplanetary-dust-settling/src

The update is linear (fixed gas g>0 and fixed velocity v), so the evolution
operator is extracted exactly by applying upwind_method to basis vectors.
"""
import os
import sys
import math
import hashlib

import numpy as np

# --------------------------------------------------------------------------- #
# Source resolution: prefer the ORIGINAL parent repo src/ in the clone layout,
# then explicit root/src via CLI/env, then a phd clone as a last-resort fallback.
# Loading uses the exact file path (not stale import_module cache) by registering
# the module file under a session-unique name and loading it with an open file
# handle.  The resolved path is recorded for provenance and sha256 hashing.
# --------------------------------------------------------------------------- #
_HERE = os.path.dirname(os.path.abspath(__file__))
# research/ lives at <repo>/cpp/research, so <repo>/src is two levels up.
_PARENT_SRC = os.path.abspath(os.path.join(_HERE, "..", "..", "src"))
_DEFAULT_SRC = r"D:/Codex/tmp/phd-github-audit-2026-10-26/protoplanetary-dust-settling/src"

_TURBULENCE_CLS = None
_SRC_MODULE = None
_SRC_FILE = None
_SRC_SHA256 = None


def _resolve_source_dir(source_dir=None):
    """Return the source directory to load, preferring parent repo src."""
    if source_dir:
        cand = os.path.abspath(source_dir)
    elif os.environ.get("DUST_SOLVER_SRC"):
        cand = os.path.abspath(os.environ["DUST_SOLVER_SRC"])
    elif os.path.isdir(_PARENT_SRC) and \
            os.path.exists(os.path.join(_PARENT_SRC, "dust_solver_turbulence.py")):
        cand = _PARENT_SRC
    else:
        cand = _DEFAULT_SRC
    if not os.path.isdir(cand):
        raise OSError("source dir not found: %r" % cand)
    if not os.path.exists(os.path.join(cand, "dust_solver_turbulence.py")):
        raise OSError("no dust_solver_turbulence.py in %r" % cand)
    return cand


def _load_module_from_path(path):
    """Load a .py module by exact path, bypassing stale import_module cache.

    Returns a fresh module object each call (session-unique name), so switching
    to a different source file resolves the requested path deterministically.
    """
    import importlib.util
    path = os.path.abspath(path)
    if not os.path.exists(path):
        raise OSError("source file not found: %r" % path)
    mod_id = "dust_solver_turbulence_%d" % (abs(hash(path)) % (2 ** 32))
    spec = importlib.util.spec_from_file_location(mod_id, path)
    if spec is None or spec.loader is None:
        raise ImportError("cannot create spec for %r" % path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[mod_id] = m
    spec.loader.exec_module(m)
    # normalise __file__ so the recorded path is the requested absolute path
    m.__file__ = path
    return m


def _import_source(source_dir=None):
    global _TURBULENCE_CLS, _SRC_MODULE, _SRC_FILE, _SRC_SHA256
    cand = _resolve_source_dir(source_dir)
    path = os.path.join(cand, "dust_solver_turbulence.py")
    m = _load_module_from_path(path)
    _SRC_MODULE = m
    _SRC_FILE = path
    _SRC_SHA256 = hashlib.sha256(open(path, "rb").read()).hexdigest()
    _TURBULENCE_CLS = m.DustSolver_Turbulence
    return _TURBULENCE_CLS


def src_identity(source_dir=None):
    """Return resolved {path, sha256} of the loaded source (empty sha never).

    On the very first call (nothing loaded yet) this loads the preferred source
    so the returned sha256 is always valid rather than null.
    """
    m = _SRC_MODULE
    if m is None:
        _import_source(source_dir)
        m = _SRC_MODULE
    p = getattr(m, "__file__", None)
    if p and os.path.exists(p):
        data = open(p, "rb").read()
        sha = hashlib.sha256(data).hexdigest()
        return {"path": os.path.abspath(p), "sha256": sha}
    return {"path": None, "sha256": None}


def source_sha256(source_dir=None):
    return src_identity(source_dir)["sha256"]


# --------------------------------------------------------------------------- #
# constants (same physical constants as original / verify.py)
# --------------------------------------------------------------------------- #
G = 6.67430e-11
M_SUN = 1.98847e30
AU = 1.495978707e11
KB = 1.380649e-23
M_H = 1.6735575e-27
YEAR = 365.25 * 24.0 * 3600.0


def _disk(R_AU):
    omega = math.sqrt(G * M_SUN / (R_AU * AU) ** 3)
    temperature = 264.0 * R_AU ** -0.5
    cs = math.sqrt(KB * temperature / (2.0 * M_H))
    H = cs / omega
    return omega, cs, H


def D_for_alpha(R_AU, alpha):
    """D = alpha * cs^2 / omega (matches verify.py _D_for_alpha)."""
    omega, cs, H = _disk(R_AU)
    return alpha * cs * cs / omega


# --------------------------------------------------------------------------- #
# operator extraction: exact linear update matrix A (rho_new = A @ rho)
# --------------------------------------------------------------------------- #
def update_matrix(sim):
    rho0 = np.zeros(sim.v.size, dtype=float)
    N = rho0.size
    cols = []
    for j in range(N):
        e = rho0.copy()
        e[j] = 1.0
        cols.append(np.asarray(sim.upwind_method(e), dtype=float))
    return np.column_stack(cols)


def operator_matrix(sim):
    """A = I - dt*B; B[i,i]>=0 loss rate, B[i,j]<=0 gain for j!=i."""
    A = update_matrix(sim)
    B = (np.eye(A.shape[0]) - A) / sim.dt
    return A, B


def closed_update_matrix(sim):
    """Update matrix for the closed_total_flux model: original update, then the
    top-face flux loss added back to the last cell (zero total top flux).
    top_flux is linear in rho[-1], so this is A_closed = A_vacuum + corr on
    the last row / last column.
    """
    A = update_matrix(sim)
    rho0 = np.zeros(sim.v.size, dtype=float)
    e = rho0.copy()
    e[sim.v.size - 1] = 1.0
    f_top = top_flux(sim, e)  # linear coefficient of rho[-1] in top flux
    corr = (sim.dt / sim.h) * f_top
    A_closed = A.copy()
    A_closed[sim.v.size - 1, sim.v.size - 1] += corr
    return A_closed


def closed_operator_matrix(sim):
    A = closed_update_matrix(sim)
    B = (np.eye(A.shape[0]) - A) / sim.dt
    return A, B


def _verify_matrix_closed(sim):
    """Sanity: closed matrix equals applying closed update to basis vectors."""
    A1 = closed_update_matrix(sim)
    rho0 = np.zeros(sim.v.size, dtype=float)
    cols = []
    for j in range(sim.v.size):
        e = rho0.copy()
        e[j] = 1.0
        out = sim.upwind_method(e)
        out = out.copy()
        out[-1] += (sim.dt / sim.h) * top_flux(sim, e)
        cols.append(np.asarray(out, dtype=float))
    A2 = np.column_stack(cols)
    return np.max(np.abs(A1 - A2))


# --------------------------------------------------------------------------- #
# derived closed-form coefficients (Codex algebra in PLAN.md, re-derived here).
# g = rho_initial (gas), q = rho/g, g_face = arithmetic mean.
# --------------------------------------------------------------------------- #
def _face_left(g, i):
    return 0.5 * (g[i] + g[i - 1])


def derived_diagonal_lambda(sim, top):
    """Vector of diagonal loss rates lambda_i for the chosen boundary.

    top="vacuum_ghost" -> original; top="closed_total_flux" -> zero top flux.
      interior i:  lam_i = [max(v[i+1],0)+max(-v[i],0)]/h
                            + alpha*(g_face_R + g_face_L)/(h^2*g_i)
      i=0 (lower, F_left=0):  lam_0 = max(v[1],0)/h + alpha*g_face_R/(h^2*g_0)
      i=N-1 vacuum:  lam = |v[-1]|/h + alpha*(g_face_L + g[-1]/2)/(h^2*g[-1])
      i=N-1 closed:  lam = max(-v[-1],0)/h + alpha*g_face_L/(h^2*g[-1])
    """
    v = np.asarray(sim.v, dtype=float)
    g = np.asarray(sim.rho_initial, dtype=float)
    h = float(sim.h)
    alpha = float(sim.alpha)
    N = v.size
    lam = np.zeros(N, dtype=float)
    if N >= 3:
        i = np.arange(1, N - 1)
        vr = v[i + 1]
        vl = v[i]
        gr = 0.5 * (g[i + 1] + g[i])
        gl = 0.5 * (g[i] + g[i - 1])
        lam[i] = (np.maximum(vr, 0.0) + np.maximum(-vl, 0.0)) / h \
            + alpha * (gr + gl) / (h * h * g[i])
    # lower boundary
    lam[0] = (max(v[1], 0.0)) / h
    gr0 = 0.5 * (g[1] + g[0])
    lam[0] += alpha * gr0 / (h * h * g[0])
    # top boundary
    gl = _face_left(g, N - 1)
    if top == "closed_total_flux":
        lam[N - 1] = max(-v[N - 1], 0.0) / h + alpha * gl / (h * h * g[N - 1])
    else:
        lam[N - 1] = abs(v[N - 1]) / h \
            + alpha * (gl + 0.5 * g[N - 1]) / (h * h * g[N - 1])
    return lam


def derived_gain_matrix(sim, top):
    """Predicted non-negative off-diagonal gain rates G[i,j] = max(0,-B[i,j]).

    The numerical B must equal -G on off-diagonal entries (gain only), and
    G[i,j]>=0 (pure positivity-preserving adjacency structure).
    """
    v = np.asarray(sim.v, dtype=float)
    g = np.asarray(sim.rho_initial, dtype=float)
    h = float(sim.h)
    alpha = float(sim.alpha)
    N = v.size
    G = np.zeros((N, N), dtype=float)
    if N >= 3:
        for i in range(1, N - 1):
            vr = v[i + 1]
            vl = v[i]
            gr = 0.5 * (g[i + 1] + g[i])
            gl = 0.5 * (g[i] + g[i - 1])
            G[i, i + 1] = max(-vr, 0.0) / h + alpha * gr / (h * h * g[i + 1])
            G[i, i - 1] = max(vl, 0.0) / h + alpha * gl / (h * h * g[i - 1])
    # lower boundary: only right-face gain
    vr0 = v[1]
    gr0 = 0.5 * (g[1] + g[0])
    G[0, 1] = max(-vr0, 0.0) / h + alpha * gr0 / (h * h * g[1])
    # top boundary: only left-face gain (both top models have no top inflow)
    vl = v[N - 1]
    gl = _face_left(g, N - 1)
    G[N - 1, N - 2] = max(vl, 0.0) / h + alpha * gl / (h * h * g[N - 2])
    return G


def top_flux(sim, rho_pre):
    """Outward top-face flux (advection outflow + diffusion ghost loss)."""
    r = np.asarray(rho_pre, dtype=float)
    v_last = float(sim.v[-1])
    return max(v_last, 0.0) * r[-1] + sim.alpha * r[-1] / (2.0 * sim.h)


# --------------------------------------------------------------------------- #
# equilibrium helpers (analytic continuous + exact discrete zero-flux recurrence)
# --------------------------------------------------------------------------- #
def analytic_closed_equilibrium(z, g, S, alpha):
    """rho = C * g * exp[-(S/alpha)*(exp(z^2/2)-1)], unnormalized (C=1).

    Valid only for constant alpha, v=-S*z/g, no velocity clip, zero total flux.
    g = exp(-z^2/2) gas profile (natural units).
    """
    return g * np.exp(-(S / alpha) * (np.exp(0.5 * z ** 2) - 1.0))


def discrete_closed_equilibrium(sim, S):
    """Exact discrete zero-flux recurrence (test-only, not a new integrator).

    q_j = q_{j-1} / [1 - h*v_j*g_j/(alpha*g_face_j)], g_face_j = (g_j+g_{j-1})/2.
    Uses v = -S*z/g (unclipped); requires sim to provide z, g, alpha, h, and a
    negative v consistent with that formula. Returns q (concentration) profile.
    """
    g = np.asarray(sim.rho_initial, dtype=float)
    z = np.asarray(sim.z, dtype=float)
    h = float(sim.h)
    alpha = float(sim.alpha)
    v = -S * z / g
    N = g.size
    q = np.zeros(N, dtype=float)
    q[0] = 1.0
    for j in range(1, N):
        gface = 0.5 * (g[j] + g[j - 1])
        denom = 1.0 - h * v[j] * g[j] / (alpha * gface)
        q[j] = q[j - 1] / denom
    return q


# --------------------------------------------------------------------------- #
# experiment runner
# --------------------------------------------------------------------------- #
class TurbulenceExperiment:
    """Thin wrapper around the original solver with two experiment switches.

    Timeline control: legacy keeps the original floor termination (Nt=int(T/dt));
    positive defaults to the requested natural endpoint. To compare with legacy,
    pass target_time_natural=legacy_actual_time explicitly. It uses ceil(T/dt_bound)
    with global dt=T/Nt (exact end alignment). The original
    requested physical inputs (total_time years, ...) are preserved in metadata.
    """

    def __init__(self, total_time, Nz, a_lim, d_lim, R_AU, a_micron, zmax,
                 D=1e-5, dust_density=1e3, v_cap_factor=1.0,
                 timestep="legacy", top="vacuum_ghost", safety=0.9,
                 source_dir=None, stop_when_cross=False,
                 target_time_natural=None, dt=None, Nt=None,
                 wallclock_budget=30.0):
        cls = _import_source(source_dir)
        self.sim = cls(total_time, Nz, a_lim, d_lim, R_AU, a_micron, zmax,
                       D, dust_density)
        self.sim.v_cap_factor = v_cap_factor
        assert timestep in ("legacy", "positive")
        assert top in ("vacuum_ghost", "closed_total_flux")
        self.timestep = timestep
        self.top = top
        self.safety = float(safety)
        self.stop_when_cross = stop_when_cross
        self.target_time_natural = (None if target_time_natural is None
                                    else float(target_time_natural))
        self._forced_dt = dt
        self._forced_Nt = Nt
        self.wallclock_budget = float(wallclock_budget)
        if self.timestep == "positive":
            _validate_positive_mode(self.safety, self.wallclock_budget)
        self.bound_dt = None
        self.lambda_max = None
        self.dt = None
        self.Nt = None
        self.t_actual = None
        self.actual_steps = None
        self.min_rho_seen = None
        self.finite_all = True
        self.cum_top_flux = 0.0
        self.cap_count = 0
        self.wallclock_elapsed = 0.0

    def _setup(self):
        self.sim.calculate_natural_unit()
        self.sim.initial_dust_points()

    def _choose_step(self):
        if self.timestep == "legacy":
            self.sim.choose_time_step()
            self.dt = self.sim.dt
            self.Nt = self.sim.Nt
            self.bound_dt = None
            self.lambda_max = None
        else:
            lam = derived_diagonal_lambda(self.sim, self.top)
            self.lambda_max = float(np.max(lam))
            if self.lambda_max <= 0.0:
                raise ValueError("max diagonal rate <=0, cannot bound dt")
            self.bound_dt = self.safety / self.lambda_max
            if self._forced_dt is not None:
                # smaller of forced dt and safety bound keeps a shared
                # vacuum/closed dt positivity-safe
                self.dt = min(float(self._forced_dt), self.bound_dt)
            else:
                self.dt = self.bound_dt
            if self.target_time_natural is not None:
                T = self.target_time_natural
            else:
                T = float(self.sim.T)
            self.Nt = int(math.ceil(T / self.dt))
            self.Nt = max(self.Nt, 1)
            self.dt = T / self.Nt  # global exact-end step
            self.sim.dt = self.dt
            self.sim.Nt = self.Nt
        if self._forced_Nt is not None:
            self.Nt = int(self._forced_Nt)

    def _step(self, rho):
        rho = np.asarray(rho, dtype=float)
        f_top = top_flux(self.sim, rho)
        out = self.sim.upwind_method(rho)
        if self.top == "closed_total_flux":
            out = out.copy()
            out[-1] += (self.sim.dt / self.sim.h) * f_top
            # closed: nothing actually leaves the top face (net flux zero)
        else:
            self.cum_top_flux += float(self.sim.dt) * f_top
        m = float(np.min(out))
        if self.min_rho_seen is None or m < self.min_rho_seen:
            self.min_rho_seen = m
        if not np.all(np.isfinite(out)):
            self.finite_all = False
        return out

    def run(self):
        import time as _time
        self._setup()
        self._choose_step()
        if self.Nt > 100000:
            raise ValueError("Nt=%d exceeds safeguard 100000" % self.Nt)
        rho = np.asarray(self.sim.rho, dtype=float).copy()
        m0 = float(np.min(rho))
        self.min_rho_seen = m0
        if not np.all(np.isfinite(rho)):
            self.finite_all = False
        # independent cap count from raw (pre-clip) settling velocity
        self.cap_count = _compute_cap_count(self.sim)
        steps = 0
        t_start = _time.monotonic()
        for n in range(self.Nt):
            rho = self._step(rho)
            steps += 1
            if self.sim.t_formation is None:
                ratio = rho / self.sim.rho_initial
                if np.any(ratio > 1.0):
                    self.sim.t_formation = (n + 1) * self.dt
                    if self.stop_when_cross:
                        break
            if n > 0 and (n % 500 == 0) and \
                    (_time.monotonic() - t_start) > self.wallclock_budget:
                raise TimeoutError("wallclock budget %.1fs exceeded at step %d"
                                   % (self.wallclock_budget, n))
        self.sim.rhofinal = rho.copy()
        self.actual_steps = steps
        self.t_actual = steps * self.dt
        self.wallclock_elapsed = _time.monotonic() - t_start
        return rho

    def diagnostics(self):
        s = self.sim
        init = np.asarray(s.rho, dtype=float)
        fin = np.asarray(s.rhofinal, dtype=float)
        mass_init = s.h * np.sum(init)
        mass_final = s.h * np.sum(fin)
        mass_residual = mass_final - mass_init + self.cum_top_flux
        diag = self._physical_diagnostics()
        ident = src_identity()
        return {
            "mode": "turbulence_experiment",
            "timestep": self.timestep,
            "top": self.top,
            "safety": self.safety,
            "N": int(s.Nz),
            "L": float(s.L),
            "h": float(s.h),
            "dt": float(self.dt),
            "Nt": int(self.Nt),
            "actual_steps": int(self.actual_steps),
            "t_actual": float(self.t_actual),
            "t_actual_years": float(self.t_actual / s.omega / YEAR),
            "target_time_natural": (self.target_time_natural
                                    if self.target_time_natural is not None
                                    else float(s.T)),
            "bound_dt": (float(self.bound_dt) if self.bound_dt is not None else None),
            "lambda_max": (float(self.lambda_max) if self.lambda_max is not None else None),
            "mass_initial": float(mass_init),
            "mass_final": float(mass_final),
            "cum_top_flux": float(self.cum_top_flux),
            "mass_residual": float(mass_residual),
            "min_rho_seen": float(self.min_rho_seen),
            "max_rho_final": float(np.max(fin)),
            "max_abs_rho_final": float(np.max(np.abs(fin))),
            "finite_all": bool(self.finite_all),
            "t_formation": (float(s.t_formation) if s.t_formation is not None else None),
            "alpha": float(s.alpha),
            "v_max": float(np.max(np.abs(s.v))),
            "rho_initial_max": float(np.max(np.abs(init))),
            "cap_count": int(self.cap_count),
            "tau_s_min": float(diag["tau_s_min"]),
            "tau_s_max": float(diag["tau_s_max"]),
            "mach_raw_max": float(diag["mach_raw_max"]),
            "mach_used_max": float(diag["mach_used_max"]),
            "source_path": ident["path"],
            "source_sha256": ident["sha256"],
            "param_sha256": self._param_hash(),
            "wallclock_elapsed": float(self.wallclock_elapsed),
            "scale_height": float(s.scale_height),
            "omega": float(s.omega),
        }

    def _param_hash(self):
        import hashlib as _hl
        key = "%r|%r|%r|%r|%r|%r|%r|%r|%r|%r|%r" % (
            self.timestep, self.top, self.safety,
            self.sim.real_total_time, self.sim.Nz, self.sim.a_lim,
            self.sim.d_lim, self.sim.real_R, self.sim.real_dust_a,
            self.sim.D, self.sim.dust_density)
        return _hl.sha256(key.encode("utf-8")).hexdigest()

    def _physical_diagnostics(self):
        s = self.sim
        ts_real = np.sqrt(np.pi / 8.0) * s.dust_density * s.real_dust_a \
            / (s.real_rho * s.sound_speed)
        tau_s = s.omega * ts_real
        raw_v_real = -s.omega ** 2 * s.z_real * ts_real
        mach_raw = np.abs(raw_v_real) / s.sound_speed
        vclipped_real = s.v * s.omega * s.scale_height
        mach_used = np.abs(vclipped_real) / s.sound_speed
        return {
            "tau_s_min": float(np.min(tau_s)),
            "tau_s_max": float(np.max(tau_s)),
            "mach_raw_max": float(np.max(mach_raw)),
            "mach_used_max": float(np.max(mach_used)),
        }

    def profile(self):
        return (np.asarray(self.sim.z, dtype=float),
                np.asarray(self.sim.rhofinal, dtype=float),
                np.asarray(self.sim.rho_initial, dtype=float))


def _validate_positive_mode(safety, wallclock_budget):
    assert 0.0 < safety <= 1.0, "positive mode safety must be in (0,1]"
    assert wallclock_budget > 0.0, "wallclock budget must be > 0"


def _compute_cap_count(sim):
    """Independent count of velocity-cap activations from the RAW settling
    velocity (pre-clip), computed without mutating solver state."""
    ts_real = np.sqrt(np.pi / 8.0) * sim.dust_density * sim.real_dust_a \
        / (sim.real_rho * sim.sound_speed)
    raw_v_real = -(sim.omega ** 2) * sim.z_real * ts_real
    v_cap = float(sim.v_cap_factor) * float(sim.sound_speed)
    return int(np.sum(np.abs(raw_v_real) > v_cap))


# --------------------------------------------------------------------------- #
# tolerance helpers (shared with tests)
# --------------------------------------------------------------------------- #
EPS_ABS = 1e-12
EPS_REL = 1e-10
MASS_TOL = 1e-10


def err_is_ok(py_val, ref_val, eps_abs=EPS_ABS, eps_rel=EPS_REL):
    return abs(py_val - ref_val) <= eps_abs + eps_rel * abs(ref_val)


def rms(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean((a - b) ** 2)))


def linf(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b)))


# --------------------------------------------------------------------------- #
# synthetic equilibrium test setup: no settling clip, v = -S*z/g, const alpha,
# zero total flux (closed) top.  Returns a sim configured for upwind_method with
# these overrides (g, z, v, alpha, h set directly; D kept for provenance).
# --------------------------------------------------------------------------- #
def make_equilibrium_sim(N, L, S, alpha, source_dir=None):
    cls = _import_source(source_dir)
    # a_micron etc. are irrelevant to the stencil; use defaults so constructor
    # is valid, then override the fields that upwind_method actually uses.
    sim = cls(1.0, N, 0.3, 0.2, 100.0, 1.0, L * (1e-3), D=1e-5, dust_density=1e3)
    sim.calculate_natural_unit()
    sim.initial_dust_points()
    # override grid/lattice to the synthetic equilibrium problem
    sim.L = float(L)
    sim.h = L / N
    sim.z = np.linspace(0.0, L, N, endpoint=False)
    sim.rho_initial = np.exp(-0.5 * sim.z ** 2)
    sim.rho = 0.01 * sim.rho_initial.copy()
    sim.alpha = float(alpha)
    sim.v = -S * sim.z / sim.rho_initial  # negative, unclipped
    sim.v_cap_factor = 1.0
    # choose a dt consistent with a positivity step (not used for equilibrium
    # profiling but kept valid); set dt/ Nt sensible
    lam = derived_diagonal_lambda(sim, "closed_total_flux")
    sim.dt = 0.9 / np.max(lam)
    sim.Nt = 1
    return sim





