import _path_setup  # noqa: F401  (adds ../src to sys.path)
from dust_solver_turbulence import DustSolver_Turbulence
import matplotlib.pyplot as plt


# index
R_AU = 30                
alpha_target = 1e-4      
total_time_yr     = 1e6     
Nz =250                  
adv_lim =0.7    
diff_lim =0.7           


a_list_micron = [1e4,1e5,1e6]  #list 
AU_m = 1.495978707e11 


tmp_solver = DustSolver_Turbulence(
    total_time=total_time_yr,
    z_spatial_points=Nz,
    advective_number_limitation=adv_lim,
    diffusion_number_limitation=diff_lim,
    R_AU=R_AU,  
    a_micron=1.0, 
    zmax=1.0,                
    D=0.0
)
tmp_solver.calculate_natural_unit()
H_m = tmp_solver.scale_height 
omega = tmp_solver.omega  
cs = tmp_solver.sound_speed
zmax_target_H = 4 
zmax_AU = zmax_target_H* H_m / AU_m

D_target = alpha_target * cs**2 / omega
print(f"R={R_AU} AU, H={H_m:.3e} m,α={alpha_target:.1e} -> D={D_target:.3e} m^2/s, zmax={zmax_AU:.3f} AU")

results = []   # (a_micron, t_form_years, crossed_flag)

   
profiles = {}  # a_micron -> (z, ratio_final)

for a_micron in a_list_micron:
    print(f"Running size a ={a_micron:.3g} micron ...")
    
    solver = DustSolver_Turbulence(
        total_time=total_time_yr,
        z_spatial_points=Nz,
        advective_number_limitation=adv_lim,
        diffusion_number_limitation=diff_lim,
        R_AU=R_AU,
        a_micron=a_micron,
        zmax=zmax_AU,
        D=D_target
    )
  
    solver.calculate_natural_unit()
    solver.initial_dust_points()
    solver.choose_time_step()
   

    solver.run_simulation(stop_when_cross=True)
    t_years = solver.get_t_formation_years()

    if t_years is None:
        print(f"  -> dust/gas never exceeded 1 (within {total_time_yr:.1e} yr)")
        crossed = False
    else:
        print(f"  -> dust/gas>1 achieved at t ≈ {t_years:.3e} yr")
        crossed = True
  
    results.append((a_micron, t_years, crossed))

    ratio_final = solver.rhofinal / solver.rho_initial
    profiles[a_micron] = (solver.z.copy(), ratio_final.copy())


print(f"(R = {R_AU} AU, alpha = {alpha_target:.1e}, t_max = {total_time_yr:.1e} yr)\n")
print("{:>10s}  {:>15s}  {:>10s}".format("a [micron]", "t_form [yr]", "cross?"))
print("-" * 40)
for a_micron, t_years, crossed in results:
    if t_years is None:
        t_str = " >{:.1e}".format(total_time_yr)
    else:
        t_str = "{:.3e}".format(t_years)
    print("{:10.3g}  {:>15s}  {:>10s}".format(a_micron, t_str, "YES" if crossed else "NO"))
