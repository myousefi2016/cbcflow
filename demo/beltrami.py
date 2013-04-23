__author__ = "Harish Narayanan <harish@simula.no>"
__date__ = "2009-01-20"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian Valen-Sendstad, 2009.
# Modified by Anders Logg, 2010.
# Modified by Martin Alnaes, 2013.

from headflow import *
from headflow.dol import *

from math import pi, e

# Problem definition
class Beltrami(NSProblem):
    "3D test problem with known analytical solution."

    def __init__(self, params):
        NSProblem.__init__(self, params)

        # Create mesh
        # We start with a UnitCubeMesh and modify it to get the mesh we
        # want: (-1, 1) x (-1, 1) x (-1, 1)

        N = self.params.N
        self.mesh = UnitCubeMesh(N, N, N)
        self.scale  = 2*(self.mesh.coordinates() - 0.5)
        self.mesh.coordinates()[:, :] = self.scale

        # The body force term
        self.f = [Constant(0), Constant(0), Constant(0)]

        # Set the viscosity
        self.mu = 1.0
        self.rho = 1.0
        nu = self.mu/self.rho

        # Characteristic velocity in the domain (used to determine timestep)
        self.U = 1.0

        # Set final time
        self.T = 0.5

        # The analytical solution
        # Velocity
        self.analytical_u = \
            ('-((a*(pow(E,a*x[2])*cos(a*x[0] + d*x[1]) + pow(E,a*x[0])*sin(a*x[1] + d*x[2])))/pow(E,pow(d,2)*t*nu))',
             '-((a*(pow(E,a*x[0])*cos(a*x[1] + d*x[2]) + pow(E,a*x[1])*sin(d*x[0] + a*x[2])))/pow(E,pow(d,2)*t*nu))',
             '-((a*(pow(E,a*x[1])*cos(d*x[0] + a*x[2]) + pow(E,a*x[2])*sin(a*x[0] + d*x[1])))/pow(E,pow(d,2)*t*nu))')
        # Pressure
        self.analytical_p = \
            ('-(rho/2.0)*(pow(a,2)*(pow(E,2*a*x[0]) + pow(E,2*a*x[1]) + pow(E,2*a*x[2]) + 2*pow(E,a*(x[1] + x[2]))*cos(d*x[0] + a*x[2])*sin(a*x[0] + d*x[1]) + 2*pow(E,a*(x[0] + x[1]))*cos(a*x[1] + d*x[2])*sin(d*x[0] + a*x[2]) + 2*pow(E,a*(x[0] + x[2]))*cos(a*x[0] + d*x[1])*sin(a*x[1] + d*x[2])))/(pow(E,pow(d,2)*t*nu))')

        # Common parameters pertinent to the functional forms above
        self.u_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e,             'nu': nu, 't': 0.0}
        self.p_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e, 'rho': 1.0, 'nu': nu, 't': 0.0}
        
        #self.exact_u = Expression(self.analytical_u, **self.u_params)
        self.exact_u = [Expression(self.analytical_u[d], **self.u_params) for d in xrange(3)]
        self.exact_p = Expression(self.analytical_p, **self.p_params)
        
        

    @classmethod
    def default_user_params(cls):
        params = ParamDict(N=20)
        return params

    def initial_conditions(self, V, Q):
        # Use analytical solutions at t = 0 as initial values
        return self.exact_u, self.exact_p

    def boundary_conditions(self, V, Q, t):
        # Set exact velocity as boundary conditions

        for e in self.exact_u: e.t = t
        
        self.exact_p.t = t
        
        bcu = [(self.exact_u, DomainBoundary())]
        bcp = []

        return bcu, bcp

    '''
    OLD FUNCTIONALITY
    '''
    '''
    def update(self, t, u, p):
        for expr in self.exact_u + self.exact_p:
            expr.t = t

    def functional(self, t, u, p):
        if t < self.T:
            return 0.0
        else:
            for expr in self.exact_u + self.exact_p:
                expr.t = t
            if not self.params.segregated:
                u = [u]
            error = 0
            for exact_u, calc_u in zip(self.exact_u, u):
                error += sqr(errornorm(exact_u, calc_u) / norm(exact_u, mesh=self.mesh))
            return sqrt(error/len(u))

    def reference(self, t):
        return 0.0

    def __str__(self):
        return "Beltrami"
    '''
