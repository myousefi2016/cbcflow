__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-11-21"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Harish Narayanan, 2009.
# Modified by Anders Logg, 2010.

from problembase import *
from scipy import *
from numpy import array
from math import pi

# Inflow boundary
class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < DOLFIN_EPS

# Outflow boundary
class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > 0.05 - DOLFIN_EPS

# Define the aneurysm region, everything outside the cylinder
class AneurysmCutoff(Expression):
    def eval(self, values, x):
        r = sqrt(x[1]**2 + x[2]**2)
        # FIXME: Is this well-defined?
        if r < 0.002 + 0.0001:
            values[0] = 0.0
        else:
            values[0] = 1.0

# Problem definition
class Problem(ProblemBase):
    "3D artery with a saccular aneurysm."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        refinement_level = options["refinement_level"]
        if refinement_level > 4:
            raise RuntimeError, "No mesh available for refinement level %d" % refinement_level
        self.mesh = Mesh("data/aneurysm_%d.xml.gz" % refinement_level)

        # The body force term
        self.f = self.uConstant((0, 0, 0))

        # Set viscosity
        self.nu = 3.5 / 1.025e6

        # Characteristic velocity in the domain (used to determine timestep)
        self.U = 2.5

        # Set end-time
        self.T = 0.05

    def initial_conditions(self, V, Q):
        if self.options['segregated']:
            return [Constant(0)] * 4
        else:
            return Constant((0, 0, 0)), Constant(0)

    def boundary_conditions(self, V, Q, t):

        # Mark domains, 0 = noslip, 1 = inflow, 2 = outflow, 3 = rest
        boundary_markers = MeshFunction("uint", V.mesh(), 2)
        boundary_markers.set_all(3)
        DomainBoundary().mark(boundary_markers, 0)
        Inflow().mark(boundary_markers, 1)
        Outflow().mark(boundary_markers, 2)

        # Create no-slip boundary condition for velocity
        self.g_noslip = self.uConstant((0, 0, 0))
        bc_noslip = [DirichletBC(V, g, boundary_markers, 0) for g in self.g_noslip]

        # Create inflow boundary condition for velocity
        inflow_exprs = ('1.0*(sin(30*t))*(1.0-(x[1]*x[1]+x[2]*x[2])/(r*r))',
                        '0.0*(sin(30*t))*(1.0-(x[1]*x[1]+x[2]*x[2])/(r*r))',
                        '0.0*(sin(30*t))*(1.0-(x[1]*x[1]+x[2]*x[2])/(r*r))')
        self.g_inflow = self.uExpr(inflow_exprs, r=0.002, t=t, degree=3)
        bc_inflow = [DirichletBC(V, g, boundary_markers, 1) for g in self.g_inflow]

        # Create outflow boundary condition for pressure
        self.g_outflow = Constant(0)
        bc_outflow = [DirichletBC(Q, self.g_outflow, boundary_markers, 2)]

        # Collect boundary conditions
        bcu = zip(bc_noslip, bc_inflow)
        bcp = zip(bc_outflow)

        return bcu + bcp

    def update(self, t, u, p):
        for g in self.g_inflow:
            g.t = t

    def functional(self, t, u, p):
         if t < self.T:
             return 0.0

         return self.uEval(u, 0, (0.025, -0.006, 0.0))

    def reference(self, t):
        """The reference value was computed using on a fine mesh
        (level 6). Values obtained for different refinement levels
        are listed below for Chorin and IPCS.

              Chorin                 IPCS
        ----------------------------------------
        -0.0325040608617000  -0.0333250879034000
        -0.0470001557641000  -0.0458749339862000
        -0.0370348732066000  -0.0364138324117000
        -0.0359768558469000  -0.0358236703894000
        -0.0356064894317000  -0.0354277722246000
        -0.0355250220872000  -0.0353312047875000
        -0.0356105862451000  -0.0354251625379000

        The reference value is taken as the average of the values
        for Chorin and IPCS on the finest mesh.
        """
        if t < self.T:
            return 0.0

        return -0.0355

    def __str__(self):
        return "Aneurysm"
