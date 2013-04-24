__author__ = "Kent-Andre Mardal <kent-and@simula.no>"
__date__ = "2008-04-03"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.
# Modified by Kristian Valen-Sendstad, 2008-2010.
# Modified by Harish Narayanan, 2009.
# Modified by Mikael Mortensen, 2009.
# Modified by Martin Alnaes, 2013.

from headflow import *
from headflow.dol import *

from numpy import array

class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1 - DOLFIN_EPS

class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

class Channel(NSProblem):
    "2D channel test problem with known analytical solution."

    def __init__(self, params):
        NSProblem.__init__(self, params)

        # Create mesh
        N = self.params.N
        self.mesh = UnitSquareMesh(N, N)

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            # Spatial parameters
            N=16,
            # Time parameters
            T=0.5,
            dt=1.0/80,
            # Physical parameters
            rho=1.0,
            mu=1.0/8.0,
            )
        return params

    def initial_conditions(self, V, Q):
        u0 = [Constant(0), Constant(0)]
        p0 = Expression("1 - x[0]")
        return (u0, p0)

    def boundary_conditions(self, V, Q, t):

        # Create no-slip boundary condition for velocity
        g_noslip = [Constant(0), Constant(0)]
        bcu = [(g_noslip, NoslipBoundary())]

        # Create boundary conditions for pressure
        bcp = [(self.pressure_bc(Q), InflowBoundary()),
               (self.pressure_bc(Q), OutflowBoundary())]

        return (bcu, bcp)

    def pressure_bc(self, Q):
        element = FiniteElement("CG", triangle, 1)
        return Expression("1 - x[0]", element=element)

    # Old code: TODO: Use these to validate
    """
    def functional(self, t, u, p):
        if t < self.T:
            return 0
        else:
            return self.uEval(u, 0, (1.0, 0.5))

    def reference(self, t):
        if t < self.T:
            return 0
        else:
            num_terms = 10000
            u = 1.0
            c = 1.0
            for n in range(1, 2*num_terms, 2):
                a = 32.0 / (DOLFIN_PI**3*n**3)
                b = (1/8.0)*DOLFIN_PI**2*n**2
                c = -c
                u += a*exp(-b*t)*c
            return u

    def tolerance(self, problem):
        return 1e-11
    """
