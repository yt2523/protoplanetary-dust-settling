import numpy as np
import matplotlib.pyplot as plt

# constant advection velocity v<0
class PureAdvection1D:
    def __init__(self, c_velocity, time_step, total_time,spatial_step):
        self.x = np.array([]) #discretize x points
        self.rho = np.array([]) # initial function points
        self.rhofinal = np.array([]) # final function points after simulation

        self.h = 0 # spatial step size
        self.v = c_velocity # constant advection velocity (- for the queation)
        self.dt = time_step # time step in each iteration
        self.T = total_time # total time of simulation
        self.Nt = int(self.T / self.dt) # number of time steps
        self.Nx = spatial_step # number of spatial steps

    def initial_function_points(self):
        # generate initial function points

        # variables
        Length=1 # length of the domain
        N= self.Nx    # number of spatial steps
        self.h = Length / N # spatial step size
        self.rho = [] # reset initial function points
        if self.h <= 0 :
            raise ValueError("Spatial step size is not valid.check initial_function_points")

        #discretize x points
        self.x = np.linspace(0, Length, N, endpoint=False)  
        if len(self.x) == 0:
            raise ValueError("x points array is empty. check initial_function_points") 
        
        #calculate initial condition points, generate square wave 
        for i in self.x:
            if 0.3<= i and i<= 0.6: # square wave between 0.3< x<0.6
                add_rho = 1.0
            else:
                add_rho = 0.0  
            self.rho.append(add_rho)

        if len(self.rho) != len(self.x):
            raise ValueError(" check initial_function_points.")
        
        self.rho = np.array(self.rho) #store as numpy array
        return None
    
    def upwind_method(self, rho):
        # implement upwind method for pure advection equation
        # use explicit Euler method for the time stepping
        # upwind from right to left

        rho_new = [] # new array after one time step

        a = -self.v * self.dt / self.h # c, v is negative, so add a minus sign

        if a < 0 or a > 1:
            raise ValueError("invalid advective number, check upwind_method.")

        for i in range(0, len(rho)-1):
            rho_new_add = rho[i]*(1-a) + rho[i+1] *a
            rho_new.append(rho_new_add)

        # periodic boundary condition
        rho_new_add = rho[-1]*(1-a) + rho[0] *a
        rho_new.append(rho_new_add)

        if len(rho_new) != len(rho):
            raise ValueError("new array length not equal to previous one,check upwind_method.")
        
        rho_new = np.array(rho_new) #store as numpy array?

        return rho_new
    
    def run_simulation(self):
        # run this section only

        # initial condition points at t0
        self.initial_function_points()
        rho = self.rho.copy()

        # time stepping 
        stepping = self.Nt

        while stepping > 0:
            rho = self.upwind_method(rho)
            stepping -= 1   #countor decrease for each time step

        self.rhofinal = np.array(rho)

        return self.rhofinal


# various advection velocity v<0, natural units formular
class PureAdvection1D_v2(PureAdvection1D): #inheritance from previous class
    def __init__(self, c_velocity, time_step, total_time,spatial_step):
        super().__init__(c_velocity, time_step, total_time,spatial_step) #inherit init from parent class    
    
    def upwind_method(self, rho):
        # implement upwind method for pure advection equation
        # use explicit Euler method for the time stepping
        # velocity can be positive or negative

        rho_new = [] # new array after one time step

        a = self.v* self.dt / self.h # advective number with sign
        abs_a = np.abs(a) # absolute advective number

        #consideration for abs_a
        if abs_a < 1e-13:    # no change if abs_a near zero
            return rho.copy() 
        
        if abs_a > 1:
            raise ValueError (f"invalid advective number {abs_a:.3f}, upwind_method do not work.")
        elif abs_a >0.75:
            print (f"advective number {abs_a:.3f} too large, may cause big error.")
        else:
            print (f"advective number {abs_a:.3f} acceptable, proceed simulation.")
        
        if a < 0:  # for v<0 ,upwind from right to left
            for i in range(0, len(rho)-1):
                rho_new_add = rho[i]*(1-abs_a) + rho[i+1] *abs_a
                rho_new.append(rho_new_add)
            
            # periodic boundary condition(v<0 add at the end)
            rho_new_add = rho[-1]* (1-abs_a) + rho[0] *abs_a
            rho_new.append(rho_new_add)

        if a > 0:  # for v>0, upwind from left to right
            
            # periodic boundary condition(v>0 add at the start)
            rho_new_add = rho[0]*(1-abs_a) + rho[-1] *abs_a
            rho_new.append(rho_new_add)

            for i in range(1, len(rho)):
                rho_new_add = rho[i]*(1-abs_a) + rho[i-1] *abs_a
                rho_new.append(rho_new_add)

        if len(rho_new) != len(rho):
            raise ValueError("new array length not equal to previous one,check upwind_method.")
        
        rho_new = np.array(rho_new) #store as numpy array?

        return rho_new
    
    def initial_function_points(self):
        # gaussian initial condition instead of square wave

        L  = 1.0          # domain length
        N  = self.Nx
        self.h = L / N    # spatial step
        if self.h <= 0:
            raise ValueError("Spatial step size invalid, check initial_function_points")

        # x grid, periodic domain [0, L)
        self.x = np.linspace(0.0, L, N, endpoint=False)
        if len(self.x) == 0:
            raise ValueError("x array empty, check initial_function_points")

        # gaussian parameters（中心和宽度随便选了一个比较稳的）
        x0    = 0.5       # centre of the bump
        sigma = 0.07      # width, 调大调小都行

        # ρ(x) = exp( -(x-x0)^2 / (2 σ^2) )
        rho_g = np.exp(-0.5 * ((self.x - x0) / sigma)**2)

        # 可以统一归一化下，峰值=1
        rho_g /= rho_g.max()

        self.rho = rho_g.copy()

        return None