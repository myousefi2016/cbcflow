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

class Womersley3D(NSProblem):
    "3D pipe test problem with known transient analytical solution."

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
        self.nu = self.params.mu / self.params.rho

        # Setup transient flow rate coefficients
        if 0:
            print "Using stationary bcs."
            self.Q_coeffs = [(0.0, self.Q), (1.0, self.Q)]
        else:
            print "Using transient bcs."
            T = self.params.T
            P = self.params.period
            tv = np.linspace(0.0, P)
            Q = self.Q * (0.3 + 0.7*np.sin(pi*((P-tv)/P)**2)**2)
            self.Q_coeffs = zip(tv, Q)

        # FIXME: This is the pouseille data, update this and analytical_solution
        self.Upeak = self.Q / (0.5 * pi * RADIUS**2)
        self.U = self.Upeak / RADIUS**2
        self.beta = 2.0 * self.nu * self.U

        print
        print "NB! This demo is work in progress."
        #print "Expected peak velocity:", self.Upeak
        #print "Expected total pressure drop:", self.beta*LENGTH
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
        # Create womersley objects
        ua = make_womersley_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains)
        for uc in ua:
            uc.set_t(t)
        pa = Expression("-beta*x[0]", beta=self.beta)
        return (ua, pa)

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        ua, pa = self.analytical_solution(spaces, t)

        # Create no-slip and inflow boundary condition for velocity
        c0 = Constant(0)
        bcu = [
            # Note the ordering, ua should be zero on the outer edges!
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
    demo_main(Womersley3D)
