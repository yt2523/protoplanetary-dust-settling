from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    HERE = Path(__file__).resolve().parent
    from settling_experiments import _resolve_source_dir
    src = Path(_resolve_source_dir()) / "dust_solver.py"

    spec = importlib.util.spec_from_file_location("dust_solver", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    case = HERE.parent / "validation" / "cases" / "lam_1000000_100_0.5_100_1_4"

    a = np.loadtxt(case / "profile.csv", delimiter=",", skiprows=1)
    z = a[:, 0]
    finaldust = a[:, 1]
    gas = a[:, 2]
    physicalz = a[:, 3]

    diag = {
        k: v
        for line in (case / "diagnostics.txt").read_text().splitlines()
        if line.strip()
        for k, v in [line.split(maxsplit=1)]
    }

    s = mod.RealDiskDustSolver(0, 100, 0.5, 100, 1, 4, 1000)
    s.calculate_natural_unit()
    s.initial_dust_points()
    Nz = s.Nz
    zz = s.z
    z_real = s.z_real
    rho_initial = s.rho_initial
    rho = s.rho
    v = s.v
    h = s.h

    assert a.shape == (100, 4)
    assert np.all(np.isfinite(a))
    assert np.all(finaldust >= 0)
    assert np.all(gas > 0)
    assert np.all(np.diff(gas) <= 0)
    assert np.all(np.isfinite(v))
    assert np.all(v <= 0)

    np.testing.assert_allclose(z, s.z, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(physicalz, s.z_real, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(gas, s.rho_initial, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(s.rho, 0.01 * gas, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(h, float(diag["h"]), rtol=1e-10, atol=1e-12)

    initial = s.rho
    final = a[:, 1]

    initialtail = h * np.cumsum(initial[::-1])[::-1]
    finaltail = h * np.cumsum(final[::-1])[::-1]

    assert np.all(finaltail <= initialtail + 1e-12 + 1e-10 * np.abs(initialtail))

    bound = 0.01 * np.cumsum(gas[::-1])[::-1] / gas
    assert np.all(bound < 1)
    assert np.all(final / gas <= bound + 1e-12 + 1e-10 * np.abs(bound))

    M0 = h * np.sum(initial)
    Mf = h * np.sum(final)
    F = float(diag["cumulative_boundary_flux"])
    assert np.abs(Mf - M0 + F) <= 1e-12 + 1e-10 * np.abs(M0)

    np.testing.assert_allclose(M0, float(diag["mass_initial"]), rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(Mf, float(diag["mass_final"]), rtol=1e-10, atol=1e-12)

    s200 = mod.RealDiskDustSolver(0, 200, 0.5, 100, 1, 4, 1000)
    s200.calculate_natural_unit()
    s200.initial_dust_points()
    g200 = s200.rho_initial
    bound200 = 0.01 * np.cumsum(g200[::-1])[::-1] / g200
    assert np.max(bound200) > 1

    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)

    payload = {
        "sourceSHA": _sha256(src),
        "maxboundN100": float(np.max(bound)),
        "maxboundN200": float(np.max(bound200)),
        "finalmaxratio": float(np.max(final / gas)),
        "massresidual": float(Mf - M0 + F),
        "tails_pass": bool(np.all(finaltail <= initialtail + 1e-12 + 1e-10 * np.abs(initialtail))),
        "note": (
            "N100 strict>1 impossible under one-way laminar, no diffusion/inflow, "
            "nonnegative and monotone gas; N200 capacity permits does not predict crossing."
        ),
    }
    (results / "grid_capacity.json").write_text(json.dumps(payload, indent=2))

    z_over_H = z
    ratio_bound = bound
    final_ratio = final / gas
    csv_rows = np.column_stack((z_over_H, initialtail, finaltail, ratio_bound, final_ratio))
    header = "z_over_H,initial_tail,final_tail,ratio_bound,final_ratio"
    np.savetxt(
        results / "grid_capacity.csv",
        csv_rows,
        delimiter=",",
        header=header,
        comments="",
    )

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
