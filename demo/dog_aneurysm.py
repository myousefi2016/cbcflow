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

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            dt=1e-3,
            period=1.0,
            num_periods=0.3,
            T=None,
            # Physical parameters
            mu=0.00345,
            rho=0.00106,
            )
        params.update(
            # Spatial discretization parameters
            mesh_file="../data/dog_mesh_37k.xml.gz",
            boundary_mesh_file="../data/dog_boundary_mesh_37k.xml.gz",
            )
        return params

    def initial_conditions(self, spaces, controls):
        "Return initial conditions as list of scalars (velocity) and scalar (pressure)."
        icu = [c0, c0, c0]
        icp = c0
        return (icu, icp)

    def boundary_conditions(self, spaces, u, p, t, controls):
        "Return boundary conditions as lists."

        factor = 1000
        profile = [0.4, 1.6, 1.4, 1.0, 0.8, 0.6, 0.55, 0.5, 0.5, 0.45, 0.4]
        profile = [p*factor for p in profile]
        time = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

        #inflow = Womersley(zip(time, profile), self.mesh, 1, self.mu/self.rho, scale_to=1000)
        #inflow = Pouseille(zip(time, profile), self.mesh, 1, scale_to=1000)

        inflow = Pouseille(zip(time, profile), self.mesh, 1)
        for e in inflow: e.set_t(float(t))

        bcu = [(inflow, 1),
               ([c0, c0, c0], 0)]
        bcp = [(c0, 2),
               (c0, 3)]
        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        inflow = bcu[0][0]
        for e in inflow: e.set_t(float(t))

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(DogAneurysm)
