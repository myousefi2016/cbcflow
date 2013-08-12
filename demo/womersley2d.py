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

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > LENGTH*(1.0 - DOLFIN_EPS)

class Womersley2D(NSProblem):
    "3D pipe test problem with known transient analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh
        N = self.params.N
        M = int(N*LENGTH/(2*RADIUS) + 0.5)
        mesh = UnitSquareMesh(M, N)
        x = mesh.coordinates()[:,0]
        y = mesh.coordinates()[:,1]
        x = LENGTH*x
        y = RADIUS*2*(y - 0.5)
        mesh.coordinates()[:,0] = x
        mesh.coordinates()[:,1] = y

        # We will apply markers with these id values
        self.wall_boundary_id = 0
        self.left_boundary_id = 1
        self.right_boundary_id = 2
        self.undefined_boundary_id = 3

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(self.undefined_boundary_id)
        DomainBoundary().mark(facet_domains, self.wall_boundary_id)
        Left().mark(facet_domains, self.left_boundary_id)
        Right().mark(facet_domains, self.right_boundary_id)

         # Setup analytical solution constants
        self.Q = 1.0
        self.nu = self.params.mu / self.params.rho

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
        self.initialize_geometry(mesh, facet_domains=facet_domains)

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
            N=32,
            )
        return params

    def analytical_solution(self, spaces, t):
        ux = Expression("U*(radius*radius - x[1]*x[1])", U=1.0, radius=RADIUS) # FIXME: Insert womersley function here
        ux.U = self.U
        uy = Constant(0)
        u0 = [ux, uy]

        p0 = Expression("-beta*x[0]", beta=1.0)
        p0.beta = self.beta

        return (u0, p0)

    def analytical_solution(self, spaces, t):
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
        # Create womersley objects
        ua = make_womersley_bcs(coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains)
        for uc in ua:
            uc.set_t(t)
        pa = Expression("-beta*x[0]", beta=self.beta)
        return (ua, pa)

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Toggle to compare built-in BC class and direct use of analytical solution:
        if 1:
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
            # Create womersley objects
            ua = make_womersley_bcs(coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains)
            for uc in ua:
                uc.set_t(t)
            pa = Constant(-self.beta*LENGTH)

        # Create no-slip and inflow boundary condition for velocity
        c0 = Constant(0)
        bcu = [
            ([c0, c0], self.wall_boundary_id), # TODO: Change ordering, should be the same if ua is actually zero on boundary
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
    demo_main(Womersley2D)
