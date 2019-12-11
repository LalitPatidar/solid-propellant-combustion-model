import numpy as np
import scipy.integrate              
        
class solidPhase:
    def __init__(self):
        self._T = 298.15
        
    @property    
    def T(self):
        return self._T
        
    @T.setter    
    def T(self,value):
        self._T = value
        
    @property    
    def conductivity(self):        
        # in CGS units 
        return 1.5e-3 - 0.115e-5*self.T
        
    @conductivity.setter
    def conductivity(self,value):
        self._conductivity = value
        
    @property
    def Cp(self):
        return 4.980e-2 + 0.660e-3*self.T
        
    @Cp.setter
    def Cp(self,value):
        self._Cp = value        
    
####################################################################################################################################
class ReactorOde(object):
    def __init__(self, solid,mass_flux):
        self.solid = solid
        self.mass_flux = mass_flux

    def __call__(self, x, y):
        """
        the ODE function, y' = f(t,y) 
        y[0]                = T                 # temperature
        y[1]                = dT/dx             # temperature gradient
        """
        
        mdot = self.mass_flux  

        self.solid.T = y[0]
        
        # governing equations
        dTdx = y[1]
        
        d2Tdx2 = (mdot*self.solid.Cp*y[1])/self.solid.conductivity
                       
        return np.hstack((dTdx,d2Tdx2))

####################################################################################################################################

def main(parameters):
    solid = solidPhase()

    # Initialize with temperature and density 
    # temperature dependent viscosity is implement via temperature setter
    # similarly temperature dependent density could be implemented
       
    mdot = parameters["mdot"]
    Ts = parameters["Ts"]
    Tmelt = parameters["Tmelt"]
    Tbeta2delta = parameters["Tbeta2delta"]
    Tinit = parameters["Tinit"]
      
    solid.T = Tbeta2delta
    
    # dTdx from energy balance of solid phase
    
    # Specific heat of solid HMX in cal/g-K
    # Reference: Beckstead et al. Progress in energy science and combustion 33 (2007) 497-551 
    # Cp_solid = 4.980e-2 + 0.660e-3*T
    # integral_Cp_dt = 4.980e-2*(Tmelt-Tinit) + 0.660e-3*((Tmelt**2)-(Tinit**2))/2.0
    
    # Enthalpy of beta->delta transition of HMX in kcal/mol
    # Reference: Beckstead et al. Progress in energy science and combustion 33 (2007) 497-551
    dHphase = 2.35 
    MW_HMX = 296.155
    
    integral_Cp_dt = 4.980e-2*(Tbeta2delta-Tinit) + 0.660e-3*((Tbeta2delta**2)-(Tinit**2))/2.0 + dHphase*1000/MW_HMX
    
    dTdx = mdot*integral_Cp_dt/solid.conductivity
    
    y0 = np.hstack((solid.T, dTdx))
    
    # Set up objects representing the ODE and the solver
    ode = ReactorOde(solid,mdot)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True, atol=1e-12,rtol=1e-6,nsteps=10000)
    solver.set_initial_value(y0, 0.0)
    
    # Integrate the equations, keeping T(t) and Y(k,t)
    dx = 0.00001

    header = '%15s %15s %15s\n' % ('x(cm)','Temp(K)', 'dT/dx')
    
    with open('delta_layer_solution.txt','w') as File:
        File.writelines(header)
  
    while solver.successful() and solver.y[0] < Tmelt:
        
        solid.T = solver.y[0]
        #Ts_liquid = solver.y[0]
        dTdx = solver.y[1]
        
        line_solid = '%15.6f %15.3f %15.3E\n' %(solver.t, solver.y[0], solver.y[1])
        
        with open('delta_layer_solution.txt','a') as File:
            File.writelines(line_solid)
                      
        solver.integrate(solver.t + dx)
