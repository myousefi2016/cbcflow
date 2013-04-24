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
        scale  = 2*(self.mesh.coordinates() - 0.5)
        self.mesh.coordinates()[:, :] = scale

        self.exact_u, self.exact_p = self.analytical_solution(t=0.0)

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            N=20,
            T=0.5,
            dt=0.05,
            rho=1.0,
            mu=1.0,
            )
        return params

    def analytical_solution(self, t):
        # The analytical solution
        # Velocity
        analytical_u = \
            ('-((a*(pow(E,a*x[2])*cos(a*x[0] + d*x[1]) + pow(E,a*x[0])*sin(a*x[1] + d*x[2])))/pow(E,pow(d,2)*t*nu))',
             '-((a*(pow(E,a*x[0])*cos(a*x[1] + d*x[2]) + pow(E,a*x[1])*sin(d*x[0] + a*x[2])))/pow(E,pow(d,2)*t*nu))',
             '-((a*(pow(E,a*x[1])*cos(d*x[0] + a*x[2]) + pow(E,a*x[2])*sin(a*x[0] + d*x[1])))/pow(E,pow(d,2)*t*nu))')

        # Pressure
        analytical_p = \
            ('-(rho/2.0)*(pow(a,2)*(pow(E,2*a*x[0]) + pow(E,2*a*x[1]) + pow(E,2*a*x[2]) + 2*pow(E,a*(x[1] + x[2]))*cos(d*x[0] + a*x[2])*sin(a*x[0] + d*x[1]) + 2*pow(E,a*(x[0] + x[1]))*cos(a*x[1] + d*x[2])*sin(d*x[0] + a*x[2]) + 2*pow(E,a*(x[0] + x[2]))*cos(a*x[0] + d*x[1])*sin(a*x[1] + d*x[2])))/(pow(E,pow(d,2)*t*nu))')

        # Common parameters pertinent to the functional forms above
        u_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e,             'nu': 1.0, 't': 0.0}
        p_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e, 'rho': 1.0, 'nu': 1.0, 't': 0.0}

        # Compile expressions
        exact_u = [Expression(analytical_u[d], **u_params) for d in xrange(3)]
        exact_p = Expression(analytical_p, **p_params)

        # Set configured physical parameters
        nu = self.params.mu / self.params.rho
        exact_u.nu = nu
        exact_p.nu = nu
        exact_p.rho = self.params.rho

        # Set time
        exact_u.t = t
        exact_p.t = t

        return (exact_u, exact_p)

    def initial_conditions(self, V, Q):
        for e in self.exact_u: e.t = 0.0
        self.exact_p.t = 0.0
        return (exact_u, exact_p)

    def boundary_conditions(self, V, Q, t):
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
    '''
