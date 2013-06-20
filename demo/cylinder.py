#!/usr/bin/env python
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 10.0)

class Cylinder(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (sqrt((x[0]-2.0)**2+(x[1]-0.5)**2) < 0.12+DOLFIN_EPS)

class FlowAroundACylinder(NSProblem):
    "Flow around a cylinder in 2D."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh
        r = Rectangle(0,0, 10, 1)
        c = Circle(2.0, 0.5, 0.12)
        mesh = Mesh(r-c, self.params.N)

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(4)
        DomainBoundary().mark(facet_domains, 0)
        Cylinder().mark(facet_domains, 0)
        LeftBoundary().mark(facet_domains, 1)
        RightBoundary().mark(facet_domains, 2)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            # Spatial parameters
            N=64,
            # Time parameters
            T=25.0,
            dt=0.05,
            # Physical parameters
            rho=1.0,
            mu=1.0/1000.0,
            )
        return params

    def initial_conditions(self, spaces, controls):
        c0 = Constant(0)
        u0 = [c0, c0]
        p0 = c0
        return (u0, p0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        c0 = Constant(0)
        c1 = Constant(1)

        # Create no-slip boundary condition for velocity
        bcu1 = ([c1, c0], 1)
        bcu2 = ([c0, c0], 0)

        # Create boundary conditions for pressure
        bcp1 = (c0, 2)

        # Collect and return
        bcu = [bcu1, bcu2]
        bcp = [bcp1]
        return (bcu, bcp)

    '''
    OLD FUNCTIONALITY
    '''
    '''
    def functional(self, t, u, p):
        # Only check final time
        if t < self.T:
            return 0
        else:
            # Compute stream function and report minimum
            psi = StreamFunction(u)
            vals  = psi.vector().array()
            vmin = MPI.min(vals.min())

            headflow_print("Stream function has minimal value %s" % vmin)

            return vmin

    def reference(self, t):
        # Only check final time
        if t < self.T:
            return 0.0
        return -0.061076605
    '''

'''
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

'''

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(FlowAroundACylinder)
