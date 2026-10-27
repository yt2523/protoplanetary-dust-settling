from pathlib import Path
import numpy as np
import settling_experiments as se

HERE = Path(__file__).resolve().parent


def test_diagnostics():
    e = se.TurbulenceExperiment(
        total_time=1,
        Nz=8,
        a_lim=.5,
        d_lim=.3,
        R_AU=100,
        a_micron=1,
        zmax=2,
        D=1e-5,
    )
    e.run()
    e.sim.rhofinal = np.array([-2., 1.] + [0.] * 6)
    d = e.diagnostics()
    assert d["max_rho_final"] == 1
    assert d["max_abs_rho_final"] == 2


def test_saved_profiles():
    root = HERE / "results"
    time = np.loadtxt(root / "time_refinement_profiles.csv", delimiter=",", skiprows=1)
    stress = np.loadtxt(
        root / "stress_profile_positive_vacuum_ghost.csv", delimiter=",", skiprows=1
    )

    with open(root / "time_refinement_profiles.csv") as f:
        header = f.readline().strip().split(",")
    assert header == [
        "z_over_H",
        "rho_safety_0p9",
        "rho_safety_0p45",
        "rho_safety_0p225",
    ]

    assert time.shape == (48, 4)
    assert np.all(np.isfinite(time))
    assert time[-1, 0] < 3.5

    assert np.allclose(time[:, 0], stress[:, 0], rtol=1e-10, atol=1e-12)
    assert (time[:, 1:] >= -1e-16).all()
