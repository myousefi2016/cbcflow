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

# Inflow boundary
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

# Define symmetric gradient
def epsilon(u):
    return 0.5*(grad(u) + (grad(u).T))

# Problem definition
class Problem(ProblemBase):
    "3D artery with a saccular aneurysm."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        refinement_level = options["refinement_level"]
        if refinement_level > 6:
            raise RuntimeError, "No mesh available for refinement level %d" % refinement_level
        self.mesh = Mesh("data/aneurysm_%d.xml.gz" % refinement_level)

        # The body force term
        self.f = Constant((0, 0, 0))

	# Set viscosity
        self.nu = 3.5 / 1.025e6
        self.U = 2.0

        # Set end-time
        self.T = DOLFIN_PI/60.0
        self.First = True

    def initial_conditions(self, V, Q):
        u0 = Constant((0, 0, 0))
        p0 = Constant(0)
        return u0, p0

    def boundary_conditions(self, V, Q, t):

        # Mark domains, 0 = noslip, 1 = inflow, 2 = outflow, 3 = rest
        boundary_markers = MeshFunction("uint", V.mesh(), 2)
        boundary_markers.set_all(3)
        DomainBoundary().mark(boundary_markers, 0)
        Inflow().mark(boundary_markers, 1)
        Outflow().mark(boundary_markers, 2)

        # Create no-slip boundary condition for velocity
        self.g0 = Constant((0, 0, 0))
        bc0 = DirichletBC(V, self.g0, boundary_markers, 0)

         # Create inflow boundary condition for velocity
        self.g1 = Expression(('norm0*(sin(30*t))*(1.0 - (x[1]*x[1] + x[2]*x[2]) / (r*r))',
                              'norm1*(sin(30*t))*(1.0 - (x[1]*x[1] + x[2]*x[2]) / (r*r))',
                              'norm2*(sin(30*t))*(1.0 - (x[1]*x[1] + x[2]*x[2]) / (r*r))'),
                             defaults = {'norm0': 1.0, 'norm1': 0.0, 'norm2': 0.0,  'r': 0.002},
                             degree=3)
        self.g1.t = t
        bc1 = DirichletBC(V, self.g1, boundary_markers, 1)

        # Create outflow boundary condition for pressure
	self.g2 = Constant(0)
        bc2 = DirichletBC(Q, self.g2, boundary_markers, 2)

        # Collect boundary conditions
        bcu = [bc0, bc1]
        bcp = [bc2]

        return bcu, bcp

    def pressure_bc(self, Q):
        return 0

    def update(self, t, u, p):
        self.g1.t = t

    def functional(self, t, u, p):

        # Only compute functional at end-time
        if t < self.T:
            return 0

        # Cutoff and scaling functions
        self.DG0 = FunctionSpace(self.mesh, "DG", 0)
        self.cutoff = AneurysmCutoff(V=self.DG0)
        self.scaling = FacetArea(self.mesh)

        # Tangent vector
        n = FacetNormal(self.mesh)
        t = as_vector([-n[1], n[0], 0])

        # Tangential component of shear stress
        sigma = 2.0*self.nu*epsilon(u)
        wss = inner(sigma*n, t)

        # Average tangential shear stress
        self.wss = self.cutoff/self.scaling*wss*ds

        # Function for plotting/debugging
        v = TestFunction(self.DG0)
        self.wss_plot = Function(self.DG0)
        self.wss_plot_form = v*self.cutoff/self.scaling*wss*ds
        self.file = File("wss.pvd")

        # Plot shear stress
        assemble(self.wss_plot_form, tensor=self.wss_plot.vector(), mesh=self.mesh)
        self.file << self.wss_plot

        # Compute shear stress
        wss = assemble(self.wss, mesh=self.mesh)
        print "WSS: %f" % wss

        return wss

    def reference(self, t):
        return 0.0

    def __str__(self):
        return "Aneurysm"
