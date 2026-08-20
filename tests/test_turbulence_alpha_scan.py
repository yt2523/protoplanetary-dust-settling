import _path_setup  # noqa: F401  (adds ../src to sys.path)
from dust_solver_turbulence import DustSolver_Turbulence

import matplotlib.pyplot as plt

R = 100
a_m = 1
total_time_yr = 5e6
Nz = 200
adv_lim = 0.7
diff_lim = 0.7

AU_m = 1.495978707e11

tmp = DustSolver_Turbulence(
    total_time=total_time_yr,
    z_spatial_points=Nz,
    advective_number_limitation=adv_lim,  
    diffusion_number_limitation=diff_lim,
    R_AU=R,
    a_micron=a_m,
    zmax=1.0,   # 占位
    D=0.0
)
tmp.calculate_natural_unit()
H_m   = tmp.scale_height     # [m]
omega = tmp.omega   
cs    = tmp.sound_speed

zmax_target_H = 3.5  #onl run in the range
zmax_AU = zmax_target_H * H_m / AU_m  

alpha_list = [1e-5, 5e-5,1e-4,5e-4,1e-3,   5e-3]

# plot z/H 
plt.figure(figsize=(10,5), dpi=150)

for alpha in alpha_list:
    D = alpha * cs**2 / omega #get d

    solver = DustSolver_Turbulence(
        total_time=total_time_yr,z_spatial_points=Nz,
        advective_number_limitation=adv_lim, diffusion_number_limitation=diff_lim,
        R_AU=R,
        a_micron=a_m,
        zmax=zmax_AU,
        D=D
    )

    solver.calculate_natural_unit()
    solver.initial_dust_points()
    solver.choose_time_step()
    rho_final = solver.run_simulation( stop_when_cross=False)

    ratio_final = rho_final / solver.rho_initial

    print(f"α={alpha:.0e}, D={D:.2e} m²/s, zmax={zmax_AU:.2f} AU, dt={solver.dt:.2e}")
    # 这里的 solver.z 就是 z/H
    plt.plot(solver.z, ratio_final, label=f"α={alpha:.0e}", linewidth=2)

plt.xlabel("z/H", fontsize=15)
plt.ylabel("ρ_d/ ρ_g", fontsize=15)
plt.title("Dust Vertical shape with Turbulence at 100 AU(z/H)", fontsize=15)
plt.legend( fontsize=12)
plt.tight_layout()
plt.show()
