#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *
from os import path

files = [path.join(path.dirname(path.realpath(__file__)),"../../../data/dog_mesh_37k.xml.gz"),
         path.join(path.dirname(path.realpath(__file__)),"../../../data/dog_mesh_97k.xml.gz"),
        ]

c0 = Constant(0)

class BifurcationAneurysm(NSProblem):
    "Template of a typical implementation of a NSProblem subclass"

    def __init__(self, params=None):
        """Initialize problem. The following are required:
        - self.mesh
        etc.
        """
        NSProblem.__init__(self, params)

        mesh = Mesh(files[self.params.refinement_level])
        self.initialize_geometry(mesh)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            dt=1e-3,
            #period=1.0,
            #num_periods=0.3,
            T=2,
            # Physical parameters
            mu=0.00345,
            rho=0.00106,
            )
        params.update(
            # Spatial discretization parameters
            refinement_level=0,
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
        #inflow = Poiseuille(zip(time, profile), self.mesh, 1, scale_to=1000)

        inflow = Poiseuille(zip(time, profile), self.mesh, 1)
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


def main():
    problem = BifurcationAneurysm()
    scheme = IPCS_Stable()

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    plot_and_save = dict(plot=True, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        ]
    postproc = NSPostProcessor({"casedir": casedir})
    postproc.add_fields(fields)

    solver = NSSolver(problem, scheme, postproc)
    solver.solve()

if __name__ == "__main__":
    main()
