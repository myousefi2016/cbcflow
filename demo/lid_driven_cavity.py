#!/usr/bin/env python
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

class Lid(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > DOLFIN_EPS and x[0] < 1.0 - DOLFIN_EPS and x[1] > 1.0 - DOLFIN_EPS

c0 = Constant(0.0)
c1 = Constant(1.0)

class LidDrivenCavity(NSProblem):
    "2D lid-driven cavity test problem with known reference value."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        N = self.params.N
        mesh = UnitSquareMesh(N, N)

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(2)
        DomainBoundary().mark(facet_domains, 0)
        Lid().mark(facet_domains, 1)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            dt=0.001,
            T=0.5,
            # Physical parameters
            rho=1.0,
            mu=1.0/1000.0,
            )
        params.update(
            # Spatial parameters
            N=32,
            )
        return params

    def initial_conditions(self, spaces, controls):
        u0 = [c0, c0]
        p0 = c0
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        g = [c1, c0]
        g0 = [c0, c0]
        bcu = [(g0, 0), (g, 1)]
        bcp = []
        return (bcu, bcp)

# Old code:
"""
    def functional(self, t, u, p):
        # Only check final time
        if t < self.T:
            return 0
        else:
            # Compute stream function and report minimum
            psi = StreamFunction(u)
            vals  = psi.vector().array()
            vmin = MPI.min(vals.min())

            headflow_print("Stream function has minimal value"  % vmin)

            return vmin

    def reference(self, t):
        # Only check final time
        if t < self.T:
            return 0.0
        return -0.061076605
"""
# Old code:
"""
def StreamFunction(u):
    "Stream function for a given 2D velocity field."

    # Fetch a scalar function (sub-)space
    try:
        V = u.function_space()
        V = V.sub(0).collapse()
    except AttributeError:
        V = u[0].function_space()

    # Check dimension
    mesh = V.mesh()
    if not mesh.topology().dim() == 2:
        error("Stream-function can only be computed in 2D.")

    # Define variational problem
    q   = TestFunction(V)
    psi = TrialFunction(V)
    a   = dot(grad(q), grad(psi))*dx()
    L   = dot(q, (u[1].dx(0) - u[0].dx(1)))*dx()

    # Define boundary condition
    g  = Constant(0)
    bc = DirichletBC(V, g, DomainBoundary())

    # Compute solution
    psi = Function(V)
    solve(a == L, psi, bc)

    return psi
"""

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(LidDrivenCavity)
