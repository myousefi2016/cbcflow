#!/usr/bin/env python
__author__ = "Oyvind Evju <oyvinev@simula.no>"
__date__ = "2013-04-15"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

#parameters["reorder_dofs_serial"] = False

c0 = Constant(0)

class DogAneurysm(NSProblem):
    "Template of a typical implementation of a NSProblem subclass"

    def __init__(self, params=None):
        """Initialize problem. The following are required:
        - self.mesh
        etc.
        """
        NSProblem.__init__(self, params)

        #print parameters["reorder_dofs_serial"]

        mesh = Mesh(self.params.mesh_file)
        self.initialize_geometry(mesh)

        self.params.T = self.params.period * self.params.num_periods

        #factor = 1000
        profile = [0.4, 1.6, 1.4, 1.0, 0.8, 0.6, 0.55, 0.5, 0.5, 0.45, 0.4]
        #profile = [p*factor for p in profile]
        time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        #self.inflow = Womersley(zip(time, profile), self.mesh, 1, self.mu/self.rho, scale_to=1000)
        #self.inflow = Pouseille(zip(time, profile), self.mesh, 1, scale_to=1000)
        self.inflow = Pouseille(zip(time, profile), mesh, 1)

    @classmethod
    def default_user_params(cls):
        """Add default parameters for this problem.

        (print NSProblem.default_params()  to see all NSProblem parameters.
        """
        params = ParamDict(
            mesh_file="mesh_37k.xml.gz",
            boundary_mesh_file="boundary_mesh_37k.xml.gz",

            period = 0.8,
            num_periods = 3,
            T = 3*0.8,

            dt=1e-3,

            mu=0.00345,
            rho=0.00106,
        )
        return params

    def initial_conditions(self, V, Q):
        "Return initial conditions as list of scalars (velocity) and scalar (pressure)."
        return [c0, c0, c0], c0

    def boundary_conditions(self, V, Q, t):
        "Return boundary conditions as lists."
        for e in self.inflow: e.set_t(t)

        bcu = [(self.inflow, 1), ([c0, c0, c0], 0)]
        bcp = [(c0, 2), (c0, 3)]

        return bcu, bcp

if __name__ == "__main__":
    p = DogAneurysm()
    show_problem(p)
