__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2009-01-01"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.

from problembase import *
from numpy import array

# Constants related to the geometry
bmarg = 1.e-3 + DOLFIN_EPS
xmin = 0.0;
xmax = 2.1;
ymin = 0.0;
ymax = 1.4;
zmin = 0.0;
zmax = 0.4;
xcenter = 0.5;
ycenter = 0.7;
radius = 0.05;

# Inflow and no-slip boundary
class InflowBoundaryTopBottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool(on_boundary and x[0] < xmax - bmarg and x[0]< bmarg \
                    or x[1] < ymin + bmarg or x[1]  > ymax - bmarg \
                or x[2] < zmin + bmarg or x[2]  > zmax - bmarg)

# Cylinder Boundary
class CylinderBoundary(SubDomain):
    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx*dx + dy*dy)
        return bool(on_boundary and r < radius + bmarg)

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary  and x[0] > (xmax - bmarg)

# No-slip boundary value for velocity
class InflowBoundaryTopBottomBoundaryValue(Expression):
    def value_shape(self):
        return (3,)
    def eval(self, values, x):
        if x[0] < bmarg:
                values[0] = 1.0
                values[1] = 0.0
                values[2] = 0.0
        else:
                values[0] = 0.0
                values[1] = 0.0
                values[2] = 0.0

def InflowBoundaryTopBottomBoundaryValueComponent(component):
    if component > 0:
        return Constant(0)
    else:
        class c(Expression):
            def eval(self, values, x):
                if x[0] < bmarg:
                    values[0] = 1.0
                else:
                    values[0] = 0.0
        return c()

# Problem definition
class Problem(ProblemBase):
    "Unicorn test problem unicorn-0.1.0/ucsolver/icns/bench2D3D."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        self.mesh = Mesh("data/cylinder_3d_bmk.xml.gz")

        # Set viscosity (Re = 1000)
        self.nu = 1.0 / 3900.0

        # Set time step
        self.dt = 0.021650635094

        # Set end time
        self.T = 10*0.021650635094

        # Create right-hand side function
        self.f =  self.uConstant((0, 0, 0))

    def initial_conditions(self, V, Q):

        u0 = self.uConstant((0, 0, 0))
        p0 = [Constant(0)]

        return u0 + p0

    def boundary_conditions(self, V, Q, t):

        # Create inflow boundary condition for velocity
        if self.options['segregated']:
            b0 = [InflowBoundaryTopBottomBoundaryValueComponent(d) for d in range(3)]
            bc0 = [DirichletBC(V, b, InflowBoundaryTopBottomBoundary()) for b in b0]
        else:
            b0 = InflowBoundaryTopBottomBoundaryValue()
            bc0 = [DirichletBC(V, b0, InflowBoundaryTopBottomBoundary())]

        # Create no-slip boundary condition
        b1 = CylinderBoundary()
        g1 = self.uConstant((0, 0, 0))
        bc1 = [DirichletBC(V, g, b1) for g in g1]

        # Create outflow boundary condition for pressure
        b2 = OutflowBoundary()
        g2 = Constant(0)
        bc2 = [DirichletBC(Q, g2, b2)]

        # Collect boundary conditions
        bcu = zip(bc0, bc1)
        bcp = zip(bc2)

        return bcu + bcp

    def functional(self, t, u, p):
        "Return value of functional of interest"
        if t > self.T - DOLFIN_EPS:
            if not self.options['segregated']:
                u = [u]
            val = 0
            for uc in u:
                val += sqr(norm(uc.vector()))
            return sqrt(val) + norm(p.vector())
        return 0.0

    def reference(self, t):
        "Return reference value for functional"
        if t > self.T - DOLFIN_EPS:
            # Reference value for verification against Unicorn (but not for solving NS)
            return 106.583962337 + 38.2586340271
        return 0.0

    def __str__(self):
        return "G2ref"
