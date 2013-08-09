#!/usr/bin/env python
__author__ = "Martin Sandve Alnes <kent-and@simula.no>"
__date__ = "2013-07-28"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

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
        self.Q = 1.0
        self.Upeak = self.Q / (0.5 * pi * RADIUS**2)
        self.U = self.Upeak / RADIUS**2
        nu = self.params.mu / self.params.rho
        self.beta = 2.0 * nu * self.U

        print
        print "Expected peak velocity:", self.Upeak
        print "Expected total pressure drop:", self.beta*LENGTH
        print

        # Store mesh and markers
        self.initialize_geometry(mesh)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-3,
            period=0.8,
            num_periods=1.0,
            # Physical parameters
            rho=1.0,
            mu=1.0/30.0,
            )
        params.update(
            # Spatial parameters
            mesh_filename="../data/pipe_0.2.xml.gz",
            )
        return params

    def analytical_solution(self, spaces, t):
        ux = Expression("U*(radius*radius - x[1]*x[1] - x[2]*x[2])", U=1.0, radius=RADIUS)
        ux.U = self.U
        c0 = Constant(0)
        u0 = [ux, c0, c0]

        p0 = Expression("-beta*x[0]", beta=1.0, length=LENGTH)
        p0.beta = self.beta

        return (u0, p0)

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Toggle to compare built-in BC class and direct use of analytical solution:
        if 0:
            print "Using analytical_solution as bcs."
            ua, pa = self.analytical_solution(spaces, t)
        else:
            if 0:
                print "Using stationary bcs."
                coeffs = [(0.0, self.Q), (1.0, self.Q)]
            else:
                print "Using transient bcs."
                T = self.params.T
                P = self.params.period
                tv = np.linspace(0.0, P)
                Q = self.Q * (0.3 + 0.7*np.sin(pi*((P-tv)/P)**2)**2)
                coeffs = zip(tv, Q)
            # Make pouseille objects
            ua = make_pouseille_bcs(coeffs, self.mesh, self.left_boundary_id, None, self.facet_domains)
            for uc in ua:
                uc.set_t(t)
            pa = Constant(-self.beta*LENGTH)

        # Create no-slip and inflow boundary condition for velocity
        c0 = Constant(0)
        bcu = [
            ([c0, c0, c0], self.wall_boundary_id),
            (ua, self.left_boundary_id),
            ]

        # Create outflow boundary conditions for pressure
        bcp = [(pa, self.right_boundary_id)]

        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        ua = bcu[1][0]
        for uc in ua:
            uc.set_t(t)

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(Pouseille3D)
