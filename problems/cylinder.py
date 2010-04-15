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

# Define the Cylinder region
class CylinderCutoff(Expression):
	def eval(self, values, x):

            dx = x[0] - xcenter
            dy = x[1] - ycenter
            r = sqrt(dx*dx + dy*dy)
            if (r <= radius + bmarg ):
		values[0] = 1.0
            else:
        	values[0] = 0.0

# Calculate symmetric gradient of the velocity field
def epsilon(u):
    return 0.5*(grad(u) + transpose(grad(u)))

# Problem definition
class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        refinement_level = options["refinement_level"]
        if refinement_level > 6:
            raise RuntimeError, "No mesh available for refinement level %d" % refinement_level
        self.mesh = Mesh("data/cylinder_%d.xml.gz" % refinement_level)
        #self.mesh = Mesh("data/cylinder_new_%d.xml.gz" % refinement_level)

        # Create right-hand side function
        self.f =  Constant((0, 0))

        # Set viscosity (Re = 1000)
        self.nu = 1.0 / 1000.0

        # Characteristic velocity in the domain (used to determinde timestep)
        self.U = 2.0

        # Set end time
        self.T  = 8.0

	self.scalar = FunctionSpace(self.mesh, "CG", 1)
	self.cutoff = CylinderCutoff(V=self.scalar)
	self.First = True

        # Option : Compute lift and drag forces on cylinder
        self.comp_forces = True

    def initial_conditions(self, V, Q):

        u0 = Constant((0, 0))
        p0 = Constant(0)

        return u0, p0

    def boundary_conditions(self, V, Q, t):

        # Create inflow boundary condition
        self.g0 = Expression(('4*Um*(x[1]*(ymax-x[1]))*sin(PI*t/8.0)/(ymax*ymax)', '0.0'))
        self.g0.Um   = 1.5
        self.g0.ymax = ymax
        self.g0.PI   = DOLFIN_PI
        self.g0.t = t
        print 'Time in bc is:', t
	self.b0 = InflowBoundary()
        bc0 = DirichletBC(V, self.g0, self.b0)

        # Create no-slip boundary condition
        self.b1 = NoslipBoundary()
        self.g1 = Constant((0, 0))
        bc1     = DirichletBC(V, self.g1, self.b1)

        # Create outflow boundary condition for pressure
        self.b2 = OutflowBoundary()
	self.g2 = Constant(0)
        bc2     = DirichletBC(Q, self.g2, self.b2)

        # Collect boundary conditions
        bcu = [bc0, bc1]
        bcp = [bc2]

        return bcu, bcp

    def update(self, t, u, p):
	self.g0.t = t
	pass

    def functional(self, t, u, p):
         if t < self.T:
             return 0
         else:
             x1 = array((xcenter - radius - DOLFIN_EPS, ycenter))
             x2 = array((xcenter + radius + DOLFIN_EPS, ycenter))

             x1 = array((xcenter - radius - DOLFIN_EPS, ycenter))
             x2 = array((xcenter + radius + DOLFIN_EPS, ycenter))
             values1 = array((0.0))
             values2 = array((0.0))
             p.eval(values1, x1)
             p.eval(values2, x2)
         return values1 - values2

    def reference(self, t):
        return 0.0

    def __str__(self):
        return "Cylinder"
