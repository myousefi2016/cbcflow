__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2009-10-01"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.

from problembase import *
from numpy import array

# Constants related to the geometry
bmarg = 1.e-3 + DOLFIN_EPS
xmin = 0.0
xmax = 2.2
ymin = 0.0
ymax = 0.41
xcenter = 0.2
ycenter = 0.2
radius = 0.05

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < xmin + bmarg

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx*dx + dy*dy)
        return on_boundary and \
               (x[1] < ymin + bmarg or x[1] > ymax - bmarg or \
                r < radius + bmarg)

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > xmax - bmarg

# Problem definition
class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        refinement_level = options["refinement_level"]
        if refinement_level > 5:
            raise RuntimeError, "No mesh available for refinement level %d" % refinement_level

        self.mesh = Mesh("data/cylinder_%d.xml.gz" % refinement_level)


        # Create right-hand side function
        self.f = self.uConstant((0,0))

        # Set viscosity (Re = 1000)
        self.nu = 1.0 / 1000.0

        # Characteristic velocity in the domain (used to determine timestep)
        self.U = 3.5

        # Set end time
        self.T  = 8.0

    def initial_conditions(self, V, Q):
        u0 = self.uConstant((0,0))
        p0 = self.uConstant(0)

        return u0 + p0

    def boundary_conditions(self, V, Q, t):

        # Create inflow boundary condition
        b0 = InflowBoundary()
        self.g0 = self.uExpr(('4*Um*(x[1]*(ymax-x[1]))*sin(pi*t/8.0)/(ymax*ymax)', '0.0'),
                             Um=1.5, ymax=ymax, t=t)
        bc0 = [DirichletBC(V, g0, b0) for g0 in self.g0]

        # Create no-slip boundary condition
        b1 = NoslipBoundary()
        self.g1 = self.uConstant((0, 0))
        bc1     = [DirichletBC(V, g1, b1) for g1 in self.g1]

        # Create outflow boundary condition for pressure
        b2 = OutflowBoundary()
        self.g2 = Constant(0)
        bc2     = [DirichletBC(Q, self.g2, b2)]

        # Collect boundary conditions
        bcu = zip(bc0, bc1)
        bcp = zip(bc2)

        return bcu + bcp

    def update(self, t, u, p):
        for g0 in self.g0:
            g0.t = t

    def functional(self, t, u, p):
        if t < self.T:
            return 0.0

        x1 = array((xcenter - radius - DOLFIN_EPS, ycenter))
        x2 = array((xcenter + radius + DOLFIN_EPS, ycenter))

        return self.eval(p, x1) - self.eval(p, x2)

    def reference(self, t):
        if t < self.T:
            return 0.0
        return -0.111444953719

    def __str__(self):
        return "Cylinder"
