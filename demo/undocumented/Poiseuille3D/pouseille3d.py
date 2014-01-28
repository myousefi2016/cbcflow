#!/usr/bin/env python
__author__ = "Martin Sandve Alnes <martinal@simula.no>"
__date__ = "2013-08-12"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from cbcflow import *
from cbcflow.dol import *

import numpy as np

LENGTH = 10.0
RADIUS = 0.5

class Pouseille3D(NSProblem):
    "3D pipe test problem with known stationary analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Load mesh
        mesh = Mesh(self.params.mesh_filename)

        # We know that the mesh contains markers with these id values
        self.wall_boundary_id = 0
        self.left_boundary_id = 2
        self.right_boundary_id = 1

        # Setup analytical solution constants
        Q = self.params.Q
        nu = self.params.mu / self.params.rho
        self.alpha = Q / (0.5 * pi * RADIUS**4)
        self.beta = 2.0 * nu * self.alpha

        # Toggle to test using Pouseille-shaped bcs with transient flow rate
        if 1:
            print "Using stationary bcs. Analytical solution should hold."
            print "Expected peak velocity:", self.alpha * RADIUS**2
            print "Expected total pressure drop:", self.beta * LENGTH
            self.Q_coeffs = [(0.0, Q), (1.0, Q)]
        else:
            print "Using transient bcs. Analytical solution will not hold."
            T = self.params.T
            P = self.params.period
            tvalues = np.linspace(0.0, P)
            Qvalues = Q * (0.3 + 0.7*np.sin(pi*((P-tvalues)/P)**2)**2)
            self.Q_coeffs = zip(tvalues, Qvalues)

        # Store mesh and markers
        self.initialize_geometry(mesh)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-4,
            period=0.8,
            num_periods=0.1,
            # Physical parameters
            rho=1.0,
            mu=1.0/30.0,
            )
        params.update(
            # Spatial parameters
            mesh_filename="../data/pipe_0.2.xml.gz",
            # Analytical solution parameters
            Q=1.0,
            )
        return params

    def analytical_solution(self, spaces, t):
        ux = Expression("alpha * (radius*radius - x[1]*x[1] - x[2]*x[2])", alpha=1.0, radius=0.5)
        ux.alpha = self.alpha
        ux.radius = RADIUS
        uy = Constant(0.0)
        uz = Constant(0.0)
        u = [ux, uy, uz]

        p = Expression("-beta * x[0]", beta=1.0)
        p.beta = self.beta

        return (u, p)

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip bcs
        d = len(u)
        u0 = [Constant(0.0)] * d
        noslip = (u0, self.wall_boundary_id)

        # Create Pouseille inflow bcs
        uin = make_pouseille_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, None, self.facet_domains)
        #uin = pouseille(self.Q_coeffs, self.mesh, self.facet_domains, self.left_boundary_id) # TODO
        for ucomp in uin:
            ucomp.set_t(t)
        inflow = (uin, self.left_boundary_id)

        # Create outflow bcs for pressure
        pa = Constant(-self.beta*LENGTH)
        outflow = (pa, self.right_boundary_id)

        # Return bcs in two lists
        bcu = [noslip, inflow]
        bcp = [outflow]
        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        uin = bcu[1][0]
        for ucomp in uin:
            ucomp.set_t(t)


def main():
    problem = Pouseille3D()
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
