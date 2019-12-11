import numpy as np
cimport numpy as np
np.import_array()
import scipy.integrate
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Phase.h":
    cdef cppclass Phase:
        Phase(string,string) except +
        double getT()
        void setT(double)
        double getP()
        void setP(double)
        double get_density()
        void set_density(double)
        double get_viscosity()
        void set_viscosity(double)
        double get_conductivity()
        void set_conductivity(double)
                
        int get_n_species()
        int get_n_reactions()
        string species_name(int)
        string reaction_name(int)
        int species_index(string)
        vector[string] get_reactants_list(int)
        vector[string] get_products_list(int)
        
        void setX(vector[double])
        vector[double] getX()
        void setY(vector[double])
        vector[double] getY()
        
        vector[double] get_molecular_weights()
        double get_mean_molecular_weight()
        vector[double] get_concentrations()
        
        vector[double] get_forward_rate_constants()
        vector[double] get_backward_rate_constants()
        
        #vector[double] get_kdiff_f()
        #vector[double] get_kdiff_b()
        
        #vector[vector[double]] get_kdiff()
        
        vector[double] get_net_forward_rate_constants()
        vector[double] get_net_backward_rate_constants()

        vector[double] get_Eaf()
        #vector[double] get_Eab()
        
        vector[double] get_net_rates_of_production()
        
        vector[double] get_Cps()
        vector[double] get_enthalpies()
        double get_Cpbar()              
        
cdef class PyPhase:
    cdef Phase* liquidphase
        
    #cdef np.ndarray _kdiff_f, _kdiff_b

    def __cinit__(self,file1,file2):
        self.liquidphase = new Phase(file1.encode('utf-8'),file2.encode('utf-8'))
        
        
    @property    
    def T(self):
        return self.liquidphase.getT()
        
    @T.setter    
    def T(self,double value):
        self.liquidphase.setT(value)
        self.liquidphase.set_density(self.density)
        self.liquidphase.set_viscosity(self.viscosity)
        self.liquidphase.set_conductivity(self.conductivity)
        
    @property    
    def P(self):
        return self.liquidphase.getP()
        
    @P.setter    
    def P(self,double value):
        self.liquidphase.setP(value)
        
    @property    
    def density(self):
        T = [550.0,600.0,650.0,700.0,750.0,800.0]
        rho = [1.6509,1.6144,1.5869,1.5545,1.5201,1.4882]
        # density in g/cm3
        return np.interp(self.T,T,rho)
        #return self.liquidphase.get_density()
        
    @density.setter    
    def density(self,double value):
        self.liquidphase.set_density(value)
        
    @property    
    def viscosity(self):
        #return self.liquidphase.get_viscosity()
        T = [550.0,600.0,650.0,700.0,750.0,800.0]
        eta = [0.45,0.12,0.04,0.0220,0.01,0.0055]
        # viscosity in g/cm-s
        return 10.0*np.interp(self.T,T,eta)
        #return 0.12
        
    @viscosity.setter
    def viscosity(self,double value):
        self.liquidphase.set_viscosity(value)
        
    @property    
    def conductivity(self):
        #T = [550.0,600.0,650.0,700.0,750.0,800.0]
        #kappa = [9.28e-4,7.51e-4,6.90e-4,6.24e-4,6.30e-4,6.40e-4]
        
        # conductivity in W/m.s obtained from cal/cm.s.K above (Bedrov. et. al for HMX)
        #return 418.4*np.interp(self.T,T,kappa)
        # in CGS units 
        return 1.5e-3 - 0.115e-5*self.T
        
    @conductivity.setter
    def conductivity(self,double value):
        self.liquidphase.set_conductivity(value)
                        
    @property
    def n_species(self):
        return self.liquidphase.get_n_species()
        
    @property
    def n_reactions(self):
        return self.liquidphase.get_n_reactions()
        
    def species_name(self, int k):
        return self.liquidphase.species_name(k).decode('ASCII')
        
    def reaction_name(self, int k):
        return self.liquidphase.reaction_name(k).decode('ASCII')
        
    def species_index(self, name):
        return self.liquidphase.species_index(name.encode('utf-8'))
    
    @property
    def Eaf(self):
        return self.liquidphase.get_Eaf()
    """    
    @property
    def Eab(self):
        return self.liquidphase.get_Eab()
    """
    @property
    def X(self):
        return self.liquidphase.getX()
        
    @X.setter
    def X(self, np.ndarray values):
        self.liquidphase.setX(values)
        
    @property
    def Y(self):
        return self.liquidphase.getY()
        
    @Y.setter
    def Y(self, np.ndarray values):
        self.liquidphase.setY(values)
        
    @property
    def molecular_weights(self):
        return self.liquidphase.get_molecular_weights()
        
    @property
    def mean_molecular_weight(self):
        return self.liquidphase.get_mean_molecular_weight()
        
    @property
    def concentrations(self):
        return self.liquidphase.get_concentrations()    
            
    @property
    def forward_rate_constants(self):
        return self.liquidphase.get_forward_rate_constants()
        
    @property
    def backward_rate_constants(self):
        return self.liquidphase.get_backward_rate_constants()
    """
    @property
    def kdiff_f(self):
        return self.liquidphase.get_kdiff_f()
        
    @property
    def kdiff_b(self):
        return self.liquidphase.get_kdiff_b()
        
    @property
    def kdiff(self):
        return np.asarray(self.liquidphase.get_kdiff())
    """    
    @property
    def net_forward_rate_constants(self):
        return self.liquidphase.get_net_forward_rate_constants()
        
    @property
    def net_backward_rate_constants(self):
        return self.liquidphase.get_net_backward_rate_constants()
        
    def reactants_list(self,int k):
        return [item.decode('utf-8') for item in list(self.liquidphase.get_reactants_list(k))]
        
    def products_list(self,int k):
        return [item.decode('utf-8') for item in list(self.liquidphase.get_products_list(k))]
        
    @property
    def net_rates_of_production(self):
        return self.liquidphase.get_net_rates_of_production()
        
    @property
    def Cp(self):
        return self.liquidphase.get_Cps()
        
    @property
    def H(self):
        return self.liquidphase.get_enthalpies()
        
    @property
    def mean_Cp(self):
        return self.liquidphase.get_Cpbar()
    
####################################################################################################################################
cdef class ReactorOde(object):
    cdef object liquid
    cdef float mass_flux
    def __init__(self, liquid,mass_flux):
        self.liquid = liquid
        self.mass_flux = mass_flux

    def __call__(self, float x, np.ndarray[np.float64_t,ndim=1] y):
        """
        the ODE function, y' = f(t,y) 
        y[0]                = T                 # temperature
        y[1]                = dT/dx             # temperature gradient
        y[2:2+n_species]    = Y                 # Species mass fractions
        """
        
        mdot = self.mass_flux
        cdef int n_species = self.liquid.n_species

        cdef int i    

        self.liquid.T = y[0]

        y = np.asarray([max([y[i],0.0]) for i in range(len(y))])

        y[2:2+n_species] = y[2:2+n_species]/sum(y[2:2+n_species])
        self.liquid.Y = y[2:2+n_species]
               
        cdef np.ndarray[np.float64_t,ndim=1] wdot = np.asarray(self.liquid.net_rates_of_production)
        
        # governing equations
        cdef float dTdx = y[1]
        
        cdef float d2Tdx2 = (mdot*self.liquid.mean_Cp*y[1]/self.liquid.mean_molecular_weight + np.sum(np.multiply(wdot,self.liquid.H)))/self.liquid.conductivity
        
        cdef np.ndarray[np.float64_t,ndim=1] dYdx = np.multiply(wdot,self.liquid.molecular_weights)/mdot
                       
        return np.hstack((dTdx,d2Tdx2, dYdx))

####################################################################################################################################

cpdef main(file1,file2,dict parameters):
    liquid = PyPhase(file1,file2)

    # Initialize with temperature and density 
    # temperature dependent viscosity is implement via temperature setter
    # similarly temperature dependent density could be implemented
       
    mdot = parameters["mdot"]
    Ts = parameters["Ts"]
    Tmelt = parameters["Tmelt"]
    Tinit = parameters["Tinit"]
      
    liquid.T = Tmelt
    cdef int n_species = liquid.n_species
    liquid.X = np.asarray([1.0 if liquid.species_name(i)=='HMX' else 0.0 for i in range(n_species)])
    
    # dTdx from energy balance of solid phase
    
    # Specific heat of solid HMX in cal/g-K
    # Reference: Beckstead et al. Progress in energy science and combustion 33 (2007) 497-551 
    # Cp_solid = 4.980e-2 + 0.660e-3*T
    # integral_Cp_dt = 4.980e-2*(Tmelt-Tinit) + 0.660e-3*((Tmelt**2)-(Tinit**2))/2.0
    
    # Enthalpy of melting of HMX in kcal/mol plus enthalpy of phase change for integral enegery balance in solid phase
    # Reference: Beckstead et al. Progress in energy science and combustion 33 (2007) 497-551
    dHmelting = 16.7 + 2.35     # kcal/mol
    MW_HMX = 296.155            # g/mol
    
    integral_Cp_dt = 4.980e-2*(Tmelt-Tinit) + 0.660e-3*((Tmelt**2)-(Tinit**2))/2.0 + dHmelting*1000/MW_HMX
    
    dTdx = mdot*integral_Cp_dt/liquid.conductivity
    
    y0 = np.hstack((liquid.T, dTdx, liquid.Y))
    
    # Set up objects representing the ODE and the solver
    ode = ReactorOde(liquid,mdot)
    solver = scipy.integrate.ode(ode)
    solver.set_integrator('vode', method='bdf', with_jacobian=True, atol=1e-12,rtol=1e-6,nsteps=10000)
    solver.set_initial_value(y0, 0.0)
    
    # Integrate the equations, keeping T(t) and Y(k,t)
    dx = 0.00001

    mass_fractions_liquid = []

    iHMX = liquid.species_index('HMX')
    header = '%15s %15s %15s' % ('x(cm)','Temp(K)', 'dT/dx')
    
    for i in range(n_species):
        header = header + '%15s' %(liquid.species_name(i))
    header = header + '\n'
    
    with open('mass_fractions_liquid.txt','w') as File:
        File.writelines(header)
  
    while solver.successful() and solver.y[0] < Ts:
        
        liquid.T = solver.y[0]
        Ts_liquid = solver.y[0]
        dTdx = solver.y[1]
        
        Yjc = np.asarray([max([solver.y[i],0.0]) for i in range(len(solver.y))])
        
        liquid.Y = Yjc[2:2+n_species]/sum(Yjc[2:2+n_species])
        
        line_liquid = '%15.6f %15.3f %15.3E' %(solver.t, solver.y[0], solver.y[1])
        
        h = np.divide(liquid.H,liquid.molecular_weights)
        
        for i in range(n_species):
            line_liquid = line_liquid + '%15.3E' %(liquid.Y[i])
            
        line_liquid = line_liquid + '\n'

        #print ('%15.3f %15.3E' %(solver.t,liquid.Y[iHMX]))
        
        with open('mass_fractions_liquid.txt','a') as File:
            File.writelines(line_liquid)
                      
        solver.integrate(solver.t + dx)
        
    #print('No. of species = %d'%(liquid.n_species))
    #print('No. of reactions = %d'%(liquid.n_reactions))
    
    Ysurf = {}
    for i in range(n_species):
        Ysurf[liquid.species_name(i)] = Yjc[i+2]/sum(Yjc[2:2+n_species])
        
    return Ts_liquid, Ysurf, liquid.conductivity, dTdx, sum(np.multiply(Yjc[2:2+n_species]/sum(Yjc[2:2+n_species]),h))
    #"""
