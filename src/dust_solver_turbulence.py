import numpy as np
import matplotlib.pyplot as plt


class DustSolver_Turbulence():
    def __init__(self, total_time, z_spatial_points, advective_number_limitation, diffusion_number_limitation, R_AU, a_micron, zmax, D=0.00001,dust_density=1e3):

        # variable in each time step
        self.Nt =  None  # number of time steps
        self.dt = None  # time step in each iteration

        # parameters_for natural units
        self.a_lim = advective_number_limitation  # max accpetatble advective number
        self.d_lim = diffusion_number_limitation   # max  accpetatble diffusion_number
        self.omega = None #Keplerian angular frequency
        self.temperature = None  # disk temperature at R
        self.sound_speed = None  # sound speed at R
        self.scale_height = None  # disk scale height H at R
        self.rho_g0 = None # midplane gas density at R 
        self.alpha = None  # diffision parameter in natural units

        # parameters in loop (natural units)
        self.T =None  # total time of simulation        t_nu=t*2pi
        self.L =None  # length of the domain          z_nu= z_nu= z/H
        self.Nz =z_spatial_points  # number of spatial points 
        self.h = None  # spatial step size
        # generation of list in time step
        self.z = np.array([])  # discretize x points
        self.rho = np.array([])  # initial function points for dust density
        self.rho_initial = np.array([])  # initial gas density points            for dust/gas ratio
        self.rhofinal = np.array([])  # final function points after simulation
        self.v = np.array([])  #  velocity array at each spatial point

        # parameters for real disk
        self.real_R = R_AU  # radius in AU for a planery real disk
        self.real_dust_a = a_micron*1e-6  # dust size in micron(radius assumed)  micron to m
        self.real_z_max = zmax*1.495978707e11  # planery real disk maximum z in m
        self.real_total_time = total_time  # total time in year
        self.dust_density = dust_density  # inner dust internal density in kg m^-3
        self.real_z = np.array([])  # real z points in m
        self.real_rho = np.array([])  # real dust density points in kg m^-3
        self.D = D  # turbulence diffusion constant in m^2 s^-1


        #record time evolution(when dust/gas>1,planet formation happens)
        self.t_formation = None   # first time when dust/gas > 1 (dimensionless time)
        self.v_cap_factor = 1

 
    def  year_to_second(self, time_year):    
        #Convert year to second 
        time_seconds = time_year * 365.25 * 24 * 3600      
        return time_seconds 

    def  calculate_natural_unit(self):
        #Compute Ω, H, rho_g0 in real units, then trans to natural units
        #index
        G= 6.67430e-11     # m^3 kg^-1 s^-2
        M_SUN= 1.98847e30    # kg, 1 solar mass
        AU= 1.495978707e11      # m
        Kb = 1.380649e-23       # J K^-1
        M_H = 1.6735575e-27       # kg

        # calculate scale number
        self.omega = np.sqrt(G * M_SUN / (self.real_R * AU)**3) #Keplerian angular frequency
        self.temperature = 264 * (self.real_R)**(-0.5)  #disk temperature at R
        self.sound_speed = np.sqrt(Kb * self.temperature / (2* M_H))  # sound speed at R
        self.scale_height = self.sound_speed / self.omega  # disk scale height H at R
        self.rho_g0 = 5e-7 * (self.real_R)**(-3) # midplane gas density at R in kg m^-3
        #self.rho_g0 = 1.4e-9 * (self.real_R)**(-11/4)
        self.alpha = self.D*  self.omega /(self.sound_speed)**2# diffision parameter in natural units

        # to natural units
        real_total_time_ = self.year_to_second(self.real_total_time)  # total time in second
        self.T = real_total_time_ *self.omega  # total time in natural units  t_nu=t*omega

        max_L = 4.0
        self.L = self.real_z_max / self.scale_height # length of the zdomain in natural units z_nu= z/H
    #    if self.L >=4.0:
    #        raise Warning(f"L {self.L:.3f}is too large")
        self.L = min(self.L, max_L) #add a limit to L to avoid too large L( in vacuum)

        self.h = self.L / self.Nz
        if self.h <= 0:
            raise ValueError("Spatial step size is not valid.check initial_function_points")
        
        return None

    def choose_time_step(self):
        # choose time step based on advective number limitation
        v_max= np.max(np.abs(self.v))  # maximum absolute velocity
        if v_max == 0:
            raise ValueError("All velocities are zero, cannot choose dt from CFL.")
        max_dt_adv = self.a_lim * self.h / v_max
        # choose time step based on diffusion number limitation (0.5 for classic)
        density_max= max(self.rho_initial) #gas density maximum
        max_dt_diff = self.d_lim * self.h**2 / (self.alpha * density_max)
        #update dt Nt
        self.dt = min(max_dt_adv, max_dt_diff)
        self.Nt = int(self.T / self.dt)  # update number of time steps
        return None

    def initial_dust_points(self):
        #set up dimensionless grid and initial dust density profile
        #generate initial function points

        #set grid （natural units）
        self.z = np.linspace(0.0, self.L, self.Nz, endpoint=False)
        AU = 1.495978707e11
        self.z_real = self.z_real = self.z * self.scale_height 
        # real_z for plotting in m

        # get rho initial points(in natural units)
        # rho_real(z_real) = rho_g0 * exp(-z_real^2 / (2H^2))
        # rho(z) = exp(- z^2 / 2)    for gas density in natural units
        self.rho = []  # reset initial function points
        self.rho_initial = np.exp(-0.5 * self.z**2)  #for gas density in natural units
        self.real_rho = self.rho_g0 * np.exp(-0.5 * (self.z_real / self.scale_height)**2) #for real gas density in kg m^-3
        self.rho = 0.01 * self.rho_initial.copy()   # rho_d = 0.01 * rho_g


        # settling time in position z 
        ts = (np.sqrt(np.pi / 8.0) *self.dust_density * self.real_dust_a /(self.real_rho * self.sound_speed))
        v_real = -(self.omega**2) * self.z_real * ts

        # speed limite？？？
        v_cap = self.v_cap_factor * self.sound_speed
        v_real = np.clip(v_real, -v_cap, v_cap)
        self.v = v_real / (self.omega * self.scale_height) # v list in natural units


        if len(self.rho) != len(self.z):
            raise ValueError(" check initial_function_points.")

        self.rho = np.array(self.rho)  # store as numpy array
        self.v = np.array(self.v)  # store as numpy array
        self.z= np.array(self.z)  # store as numpy array
        return None
    
    def upwind_method(self, rho):
        # implement upwind method+diffusion
        # use explicit Euler method for the time stepping
        # upwind from right to left for v<0, left to right for v>0
        # v is a list of velocities for each i, v[0]= 0 at boundary,left boundary v= v[i], right boundary v=v[i+1]
        # self.dt is time step
        # self.h is spatial step size
        rho = np.asarray(rho)
        N = rho.size
        if N == 0:
            raise ValueError("rho array is empty.")
        if self.v.size != N:
            raise ValueError("velocity array length not equal to rho array length, check upwind_method.")

        rho_new = np.empty_like(rho)   

        ratio = rho / self.rho_initial
        h = self.h  
        dt_over_h = self.dt / h

        #i = 0 
        if self.v[1] >= 0.0:
            F_right_adv = self.v[1] * rho[0]
        else:
            F_right_adv = self.v[1] * rho[1]
        F_left_adv =0.0 

        right_slope =(ratio[1] -ratio[0]) / h
        rho_g_right = 0.5 * (self.rho_initial[1] +self.rho_initial[0])
        F_right_diff = - self.alpha * rho_g_right *right_slope   
        F_left_diff = 0.0

        F_right_total = F_right_adv + F_right_diff
        F_left_total = F_left_adv +F_left_diff

        rho_new[0]= rho[0] -dt_over_h *(F_right_total - F_left_total)
   
        if N > 2:
            i = np.arange(1, N-1)

            v_right = self.v[i+1] 
            v_left  = self.v[i]  

            #advection
            F_right_adv = np.where(v_right >= 0.0,
                                v_right *rho[i],
                                v_right * rho[i+1])
            F_left_adv = np.where(v_left >= 0.0,
                                v_left *rho[i-1],
                                v_left * rho[i])
   
            # diffusion
            right_slope = (ratio[i+1] - ratio[i]) / h
            left_slope= (ratio[i] - ratio[i-1]) / h
   
            rho_g_right = 0.5 * (self.rho_initial[i+1] + self.rho_initial[i])   
            rho_g_left = 0.5 * (self.rho_initial[i]   + self.rho_initial[i-1])

            F_right_diff = - self.alpha * rho_g_right *right_slope
            F_left_diff = - self.alpha * rho_g_left  *left_slope

            F_right_total = F_right_adv + F_right_diff
            F_left_total = F_left_adv  + F_left_diff

            rho_new[i] = rho[i] - dt_over_h *(F_right_total -F_left_total)   

        #i = N-1
        i_last = N - 1

        if self.v[i_last] >= 0.0:
            F_right_adv = self.v[i_last] * rho[i_last]  
            F_left_adv  = self.v[i_last] * rho[i_last-1]
        else:
            F_right_adv = 0.0
            F_left_adv  = self.v[i_last] * rho[i_last]

        right_slope = (0.0 - ratio[i_last]) / h
        rho_g_right = 0.5 * (0.0 + self.rho_initial[i_last])
        left_slope  = (ratio[i_last] - ratio[i_last-1]) / h
        rho_g_left  = 0.5 * (self.rho_initial[i_last] + self.rho_initial[i_last-1])   

        F_right_diff = - self.alpha * rho_g_right *right_slope
        F_left_diff  = - self.alpha *rho_g_left  * left_slope
  
        F_right_total = F_right_adv + F_right_diff
        F_left_total  = F_left_adv  + F_left_diff
   
        rho_new[i_last] =rho[i_last] - dt_over_h *(F_right_total - F_left_total)

        return rho_new
    
    # def run_simulation(self):

    #     self.calculate_natural_unit()   
     
    #     self.initial_dust_points()
    #     rho = self.rho.copy()

    #     self.choose_time_step()        

    #     #time loop
    #     for _ in range(self.Nt):
    #        rho = self.upwind_method(rho)

    #     self.rhofinal = rho.copy()
    #     return self.rhofinal
    
    def run_simulation(self, stop_when_cross=False):

        self.calculate_natural_unit()
        self.initial_dust_points()
        rho = self.rho.copy()

        self.choose_time_step()

        for n in range(self.Nt):
            rho = self.upwind_method(rho)

            # if self.t_formation is None:
            #     ratio = rho / self.rho_initial   # dust/gas ratio
            #     # ---- 只看中平面附近三格 i = 0,1,2 ----
            #     n_core = min(3, len(ratio))   
            #     core_ratio_mean = np.mean(ratio[:n_core])
            #     if core_ratio_mean > 1.0:
            #         self.t_formation = (n + 1) * self.dt  # dimensionless time Ω t
            #         if stop_when_cross:
            #             break

            if self.t_formation is None:
                ratio = rho / self.rho_initial   # dust/gas ratio
                if np.any(ratio > 1.0):
                    self.t_formation = (n + 1) * self.dt  # Ω t
                    if stop_when_cross:
                        break   

        self.rhofinal = rho.copy()
        return self.rhofinal
    
    def get_t_formation_years(self):
            if self.t_formation is None:
                return None
            
            t_seconds = self.t_formation / self.omega       # s
            t_years = t_seconds / (365.25 * 24 * 3600) # yr
            return t_years