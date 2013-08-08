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

        # Setup analytical solution constants # FIXME: This is the pouseille data
        self.Upeak = 5.0
        self.U = self.Upeak / RADIUS**2
        self.nu = self.params.mu / self.params.rho
        self.beta = 2.0 * self.nu * self.U

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
        ux = Expression("U*(radius*radius - x[1]*x[1])", U=1.0, radius=RADIUS) # FIXME: Insert womersley function here
        ux.U = self.U
        uy = Constant(0)
        u0 = [ux, uy]

        p0 = Expression("-beta*x[0]", beta=1.0)
        p0.beta = self.beta

        return (u0, p0)

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Toggle to compare built-in BC class and direct use of analytical solution:
        if 0:
            ua, pa = self.analytical_solution(spaces, t)
        else:
            # FIXME: FIRST make womersley match the stationary velocity provided here.
            # TODO: THEN test a transient flow rate.
            coeffs = [(0.0, self.Upeak), (1.0, self.Upeak)]
            ua = make_womersley_bcs(coeffs, self.mesh, self.left_boundary_id, self.nu, None, self.facet_domains)
            pa = Constant(-self.beta*LENGTH)

        # Create no-slip and inflow boundary condition for velocity
        c0 = Constant(0)
        bcu = [
            ([c0, c0], self.wall_boundary_id),
            (ua, self.left_boundary_id),
            ]

        # Create outflow boundary conditions for pressure
        bcp = [(pa, self.right_boundary_id)]

        return (bcu, bcp)

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(Womersley2D)
