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
        Q = self.params.Q
        self.nu = self.params.mu / self.params.rho

        # Beta is the Pouseille pressure drop if the flow rate is stationary Q
        self.beta = 4.0 * self.nu * Q / (pi * RADIUS**4)

        # Setup transient flow rate coefficients
        if 0:
            print "Using stationary bcs."
            self.Q_coeffs = [(0.0, self.Q), (1.0, self.Q)]
        else:
            print "Using transient bcs."
            T = self.params.T
            P = self.params.period
            tvalues = np.linspace(0.0, P)
            #Qfloor, Qpeak = 1.0, 0.0
            #Qfloor, Qpeak = 0.3, 0.7
            Qfloor, Qpeak = -0.2, 1.0
            Qvalues = Q * (Qfloor + (Qpeak-Qfloor)*np.sin(pi*((P-tvalues)/P)**2)**2)
            self.Q_coeffs = zip(tvalues, Qvalues)

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
            # Analytical solution parameters
            Q=1.0,
            )
        return params

    def analytical_solution(self, spaces, t):
        # Create womersley objects
        ua = make_womersley_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains)
        #ua = womersley(self.Q_coeffs, self.mesh, self.facet_domains, self.left_boundary_id, self.nu) # TODO
        for uc in ua:
            uc.set_t(t)
        pa = Expression("-beta * x[0]", beta=1.0)
        pa.beta = self.beta # TODO: This is not correct unless stationary...
        return (ua, pa)

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip bcs
        d = len(u)
        u0 = [Constant(0.0)] * d
        noslip = (u0, self.wall_boundary_id)

        # Get other bcs from analytical solution functions
        ua, pa = self.analytical_solution(spaces, t)

        # Create inflow boundary conditions for velocity
        inflow = (ua, self.left_boundary_id)

        # Create outflow boundary conditions for velocity
        u_outflow = (ua, self.right_boundary_id)

        # Create outflow boundary conditions for pressure
        p_outflow = (pa, self.right_boundary_id)

        # Return bcs in two lists
        bcu = [noslip, inflow]
        bcp = []
        if 1: # Switch between pressure or dbc at outlet
            bcp += [p_outflow]
        else:
            bcu += [u_outflow]
        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls):
        bcu, bcp = bcs
        uin = bcu[1][0]
        for ucomp in uin:
            ucomp.set_t(t)

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(Womersley3D)
