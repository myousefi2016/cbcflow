__author__ = "Kent-Andre Mardal <kent-and@simula.no>"
__date__ = "2008-04-03"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.
# Modified by Kristian Valen-Sendstad, 2008-2010.
# Modified by Harish Narayanan, 2009.
# Modified by Mikael Mortensen, 2009.

from problembase import *
from numpy import array

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1 - DOLFIN_EPS

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

# Problem definition
class Problem(ProblemBase):
    "2D channel test problem with known analytical solution."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        N = options["N"]
        self.mesh = UnitSquare(N, N)

        # Create right-hand side function with pressure gradient as body force
        self.f = Constant((0, 0))

        # Set viscosity (Re = 8)
        self.nu = 1.0 / 8.0

        # Set end-time
        self.T = 0.5

    def initial_conditions(self, V, Q):

        u0 = Constant((0, 0))
        p0 = Expression("1 - x[0]")

        return u0, p0

    def boundary_conditions(self, V, Q, t):

        # Create no-slip boundary condition for velocity
        bv = DirichletBC(V, Constant((0.0, 0.0)), NoslipBoundary())

        # Create boundary conditions for pressure
        bp0 = DirichletBC(Q, self.pressure_bc(Q), InflowBoundary())
        bp1 = DirichletBC(Q, self.pressure_bc(Q),  OutflowBoundary())

        bcu   = [bv]
        bcp   = [bp0, bp1]

        return bcu, bcp

    def pressure_bc(self, Q):
        element = FiniteElement("CG", triangle, 1)
        return Expression("1 - x[0]", element=element)

    def functional(self, t, u, p):
        if t < self.T:
            return 0
        else:
            x = array((1.0, 0.5))
            values = array((0.0 , 0.0))
            u.eval(values, x)
            return values[0]

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

    def __str__(self):
        return "Channel"
