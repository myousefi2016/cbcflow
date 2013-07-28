#!/usr/bin/env python
__author__ = "Martin Sandve Alnes <kent-and@simula.no>"
__date__ = "2013-07-28"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

from numpy import array

LENGTH = 10.0
RADIUS = 0.5

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > LENGTH*(1.0 - DOLFIN_EPS)

class Pouseille2D(NSProblem):
    "3D pipe test problem with known stationary analytical solution."

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

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(4)
        DomainBoundary().mark(facet_domains, 0)
        Left().mark(facet_domains, 1)
        Right().mark(facet_domains, 2)

        # Setup analytical solution constants
        self.Upeak = 5.0
        self.U = self.Upeak / RADIUS**2
        nu = self.params.mu / self.params.rho
        self.beta = 2.0 * nu * self.U

        print
        print "Expected peak velocity:", self.Upeak
        print "Expected total pressure drop:", self.beta*LENGTH
        print

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=0.5,
            dt=1.0/100,
            # Physical parameters
            rho=1.0,
            mu=1.0/10.0,
            )
        params.update(
            # Spatial parameters
            N=32,
            )
        return params

    def analytical_solution(self, spaces, t):
        ux = Expression("U*(radius*radius - x[1]*x[1])", U=1.0, radius=RADIUS)
        ux.U = self.U
        uy = Constant(0)
        u0 = [ux, uy]

        p0 = Expression("-beta*x[0]", beta=1.0, length=LENGTH)
        p0.beta = self.beta

        return (u0, p0)

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        ua, pa = self.analytical_solution(spaces, t)

        # Create no-slip and inflow boundary condition for velocity
        c0 = Constant(0)
        bcu = [
            ([c0, c0], 0),
            (ua, 1),
            ]

        # Create outflow boundary conditions for pressure
        bcp = [(pa, 2)]

        return (bcu, bcp)

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(Pouseille2D)
