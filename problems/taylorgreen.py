__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-09-15"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Harish Narayanan, 2009.
# Modified by Anders Logg, 2010.

from math import pi
from problembase import *

class PeriodicBoundaryX(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < (-1.0 + DOLFIN_EPS) and x[0] > (-1.0 - DOLFIN_EPS) and on_boundary

    def map(self, x, y):
        y[0] = x[0] - 2.0
        y[1] = x[1]

class PeriodicBoundaryY(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < (-1.0 + DOLFIN_EPS) and x[1] > (-1.0 - DOLFIN_EPS) and on_boundary

    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1] - 2.0

# Problem definition
class Problem(ProblemBase):
    "2D periodic problem with known analytical solution."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        N = options["N"]
        self.mesh = UnitSquare(N, N)
        self.scale = 2*(self.mesh.coordinates() - 0.5)
        self.mesh.coordinates()[:, :] = self.scale

        # The body force term
        self.f = self.uConstant((0, 0))

        # Set viscosity
#        self.nu = 1.0 / 10.0 # A higher viscosity is used for the hk-refinement test to create a more time dependent solution.
        self.nu = 1.0 / 100.0

        # Set the final time
        self.T = 0.5

        # The time step may be fixed by setting self.dt=0.5 (Used to run hk-refinement)
#        self.dt = 0.5

        # Characteristic velocity in the domain (used to determine timestep if dt is not set)
        self.U = 1.0

        self.analytical_u = ('-(cos(pi*(x[0]))*sin(pi*(x[1]))) * exp(-2.0*nu*pi*pi*t)',
                             ' (cos(pi*(x[1]))*sin(pi*(x[0]))) * exp(-2.0*nu*pi*pi*t)')
        self.analytical_p = '-0.25*(cos(2*pi*(x[0])) + cos(2*pi*(x[1]))) * exp(-4.0*nu*pi*pi*t)'

    def initial_conditions(self, V, Q):

        # Use analytical solutions at t = 0 as initial values
        self.exact_u = self.uExpr(self.analytical_u, nu=self.nu, t=0.0, degree=3)
        self.exact_p = self.uExpr(self.analytical_p, nu=self.nu, t=0.0, degree=3)

        return self.exact_u + self.exact_p

    def boundary_conditions(self, V, Q, t):

        # Periodic boundaries
        px = PeriodicBoundaryX()
        py = PeriodicBoundaryY()

        # Periodic boundary conditions for velocity
        bux = PeriodicBC(V, px)
        buy = PeriodicBC(V, py)
        bcu = [(bux, buy)]
        if self.options['segregated']:
            bcu *= 2

        # Periodic boundary conditions for pressure
        bpx = PeriodicBC(Q, px)
        bpy = PeriodicBC(Q, py)
        bcp = [(bpx, bpy)]

        return bcu + bcp

    def update(self, t, u, p):
        for expr in self.exact_u + self.exact_p:
            expr.t = t

# Logg: One option is to subtract the reference from the functional and return 0 in the reference.

    def functional(self, t, u, p):
        if t < self.T:
            return 0
        else:
            if not self.options['segregated']:
                u = [u]
            val = 0
            for expr in u:
                val += 0.5*sqr(norm(expr, mesh=self.mesh))
            return val

    def reference(self, t):
        if t < self.T:
            return 0
        else:
            val = 0
            for expr in self.exact_u:
                expr.t = t
                val += 0.5*sqr(norm(expr, mesh=self.mesh))
            return val

    def __str__(self):
        return "Taylor-Green"
