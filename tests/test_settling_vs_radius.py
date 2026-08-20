import _path_setup  # noqa: F401  (adds ../src to sys.path)
import numpy as np
from dust_solver import RealDiskDustSolver

# radius list
R_list = [   10, 30, 50, 70, 100]
t_settle_list =[]

for R in R_list:
    solver = RealDiskDustSolver(
        total_time= 1e7,             
        z_spatial_points=200,
        advective_number_limitation= 0.7,
        R_AU=R,
        a_micron =1.0,
        zmax=R,                     
        dust_density=1e3
    )

    solver.run_simulation(stop_when_cross= True)
    t_years = solver.get_t_formation_years()

    if t_years is None:
        print(f"R = {R:3} AU :NO CROSS (<1e6 yr)")
        t_settle_list.append(np.nan)
    else:
        print(f"R = {R:3} AU :t_cross = {t_years:.6e} yr")
        t_settle_list.append(t_years)

print("\n--- Raw array (for copy-paste) ---")
print(np.array(t_settle_list))

