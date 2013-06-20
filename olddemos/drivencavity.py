#!/usr/bin/env python
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

def boundaryvalue(x):
    if x[0] > DOLFIN_EPS and x[0] < 1.0 - DOLFIN_EPS and x[1] > 1.0 - DOLFIN_EPS:
        return [1.0, 0.0]
    else:
        return [0.0, 0.0]

class BoundaryValueVec(Expression):
    def value_shape(self):
        return (2,)
    def eval(self, values, x):
        values[:] = boundaryvalue(x)

class BoundaryValueComp(Expression):
    def __init__(self, component, **kwargs):
        Expression.__init__(self, **kwargs)
        self.component = component
    def eval(self, values, x):
        values[0] = boundaryvalue(x)[self.component]

c0 = Constant(0)

class DrivenCavity(NSProblem):
    "2D lid-driven cavity test problem with known reference value."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)
        N = self.params.N
        mesh = UnitSquareMesh(N, N)

        # Store mesh and markers
        self.initialize_geometry(mesh)

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            N=16,

            dt=0.01,
            T=2.5,

            rho=1.0,
            mu=1.0/1000.0,
            )
        return params

    def initial_conditions(self, V, Q):
        u0 = [c0, c0]
        p0 = c0
        return (u0, p0)

    def boundary_conditions(self, V, Q, t):
        element = FiniteElement("CG", triangle, 1)
        g = [BoundaryValueComp(d, element=element) for d in range(2)]
        bcu = [(g, DomainBoundary())]

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
    p = DrivenCavity()
    show_problem(p)
