__author__ = "Oyvind Evju <oyvinev@simula.no>"
__date__ = "2013-04-15"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from dolfin import *
from headflow import *

#parameters["reorder_dofs_serial"] = False

# Problem definition
class Aneurysm(NSProblem):
    "Template of a typical implementation of a NSProblem subclass"

    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        '''
        Initiate problem. The following are required:
        - self.mesh
        - self.mu
        - self.rho
        - self.T
        - self.f
        etc.       
        '''
        
        #print parameters["reorder_dofs_serial"]
        
        self.mesh = Mesh(self.params.mesh_file)
        #self.mu = 0.00345
        #self.rho = 0.00106
        
        
        
        self.f = [Constant(0), Constant(0), Constant(0)]
        
        self.params.T = 3*self.params.period
        #self.T = 1.6e-3
        
        #factor = 1000
        profile = [0.4, 1.6, 1.4, 1.0, 0.8, 0.6, 0.55, 0.5, 0.5, 0.45, 0.4]
        #profile = [p*factor for p in profile]
        time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        
        #self.inflow = Womersley(zip(time, profile), self.mesh, 1, self.mu/self.rho, scale_to=1000)
        #self.inflow = Pouseille(zip(time, profile), self.mesh, 1, scale_to=1000)
        self.inflow = Pouseille(zip(time, profile), self.mesh, 1)


    @classmethod
    def default_user_params(cls):
        # Add default parameters (print NSProblem.default_params() for NSProblem parameters)
        params = ParamDict(
            period = 0.8,
            mesh_file="mesh_37k.xml.gz",
            boundary_mesh_file="boundary_mesh_37k.xml.gz",
            mu=0.00345,
            rho=0.00106,
            dt=1e-3,
        )
        return params

    def initial_conditions(self, V, Q):
        #import ipdb; ipdb.set_trace()
        # Return initial conditions as list of scalars (velocity) and scalar (pressure)
        return [Constant(0), Constant(0), Constant(0)], Constant(0)
        

    def boundary_conditions(self, V, Q, t):
        #import ipdb; ipdb.set_trace()
        # Return boundary conditions as lists. 
        for e in self.inflow: e.set_t(t)
        
        bcu = [(self.inflow, 1), ([Constant(0), Constant(0), Constant(0)], 0)]
        bcp = [(Constant(0), 2), (Constant(0), 3)]
        

        return bcu, bcp
        
    def __str__(self):
        # Name of problem
        return "Dog aneurysm"

    
