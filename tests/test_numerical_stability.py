import _path_setup  # noqa: F401  (adds ../src to sys.path)
import numpy as np
import matplotlib.pyplot as plt

from dust_solver_turbulence import DustSolver_Turbulence


#index
R_AU_TEST = 100         # Test radius in AU
A_MICRON_TEST  = 1            # Dust grain size in micron
ZMAX_AU_TEST = 4            # Vertical domain upper bound (AU)
TOTAL_TIME_YR  =3.0e6   # Shortened run time for test



ADV_   =1.0  # max advective number allowed
diff=1.0       # max diffusion number allowed



def run_one_case(
    Nz,
    cfl_factor   =0.9,
    v_cap_factor=0.5,
    stop_when_cross=False,
):
   # Initialise solver
    solver = DustSolver_Turbulence( total_time=TOTAL_TIME_YR, z_spatial_points=Nz,     advective_number_limitation=ADV_,
        diffusion_number_limitation=diff,    R_AU=R_AU_TEST,   a_micron=A_MICRON_TEST,
        zmax=ZMAX_AU_TEST,
    )

    solver.v_cap_factor = v_cap_factor
    solver.calculate_natural_unit()
    solver.initial_dust_points()
    solver.choose_time_step()

    solver.dt *= cfl_factor
    solver.Nt = int(solver.T/ solver.dt)
    rho = solver.rho.copy()
    rho_g = solver.rho_initial.copy()

    for n in range(solver.Nt):
        rho = solver.upwind_method( rho)

        if solver.t_formation is None:
            ratio = rho/ rho_g
            if np.any(ratio >1.0):
                solver.t_formation = (n +1) *solver.dt
                if stop_when_cross:
                    break

    solver.rhofinal =rho.copy()
    mid_ratio = (rho / rho_g)[0]
    return solver.z.copy(), rho.copy(), mid_ratio, solver.get_t_formation_years()


#  Resolution test 
def test_resolution():
    Nz_list = [250]
    results = []

    plt.figure()
    for Nz in Nz_list:
        z, rho, mid_ratio, t_cross = run_one_case(
            Nz=Nz,
            cfl_factor=0.9,
            v_cap_factor=0.5,
            stop_when_cross=False,
        )
        results.append((Nz, mid_ratio, t_cross))
        plt.plot(z, rho, label=f"Nz={Nz}")

    plt.xlabel("z / H")
    plt.ylabel("Dust density (arb. units)")
    plt.title("Resolution Test")
    plt.legend()
    plt.tight_layout()
    plt.savefig("resolution_profiles.png", dpi=200)

    print("\n[Resolution Test] Nz vs mid-plane ratio and t_cross:")
    print("Nz mid_ratio_final t_cross_years")
    for Nz, mid_ratio, t_cross in results:
        print(f"{Nz:3d}{mid_ratio:16.6e}  {t_cross if t_cross is not None else'None'}")

#  CFL test (cfl_factor = 0.3, 0.5, 0.7, 0.9)

def test_cfl_():
    Nz_fixed = 200
    cfl_list = [0.3, 0.5, 0.7, 0.9]
    results = []

    print("\n[CFL Test] cfl_factor vs mid-plane ratio & t_cross (Nz=200):")
    print("cfl_factor mid_ratio_final  t_cross_years")

    for cf in cfl_list:
        _, _, mid_ratio, t_cross = run_one_case(
            Nz=Nz_fixed,v_cap_factor=0.5,
            cfl_factor=cf, 
            stop_when_cross=False,
        )
        results.append((cf, mid_ratio, t_cross))
        print(f"{cf:10.2f}  {mid_ratio:16.6e}  {t_cross if t_cross is not None else 'None'}")


#  Velocity cap test
def test_vcap_():
    Nz_fixed = 200
    cap_list = [0.3, 0.5, 0.7]
    results = []

    print("\n[v_cap Test] v_cap_factor vs mid_ratio & t_cross (Nz=200):")
    print("v_cap_factor mid_ratio_final t_cross_years")

    for vf in cap_list:
        _, _, mid_ratio, t_cross = run_one_case( Nz=Nz_fixed, cfl_factor=0.9,
            v_cap_factor= vf,
            stop_when_cross =False, )
        results.append((vf, mid_ratio, t_cross))
        print(f"{vf:12.2f}  {mid_ratio:16.6e}  {t_cross if t_cross is not None else 'None'}")


# run
if __name__ == "__main__":
     test_resolution()
    # test_cfl_()
    # test_vcap_()
