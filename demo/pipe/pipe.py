#!/usr/bin/env python
__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-04-24"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

c0 = Constant(0)

def as_scalar_space(V):
    if V.num_sub_spaces() == 0:
        return V
    else:
        return V.sub(0).collapse()

class Pipe(NSProblem):
    "3D pipe test problem with known analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Get 3D pipe mesh from file
        mesh = Mesh("../../data/pipe_0.2.xml.gz")
        self.initialize_geometry(mesh)

        # Known properties of the mesh
        self.length = 10.0
        self.radius = 0.5

        # Set end time based on period and number of periods NB! Overrides given T!
        self.params.T = self.params.period * self.params.num_periods

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            # Physical parameters
            rho=1.0,
            mu=0.035,

            # Pressure gradient amplitude
            beta=5.0,
            period=0.8,
            num_periods=3,

            # Time parameters
            T=0.8*3,
            dt=1e-3,

            # Control parameters
            alpha=1e-4,
            pdim=1,
            )
        return params

    def controls(self, V, Q):
        V = as_scalar_space(V)
        u0 = [Function(V) for i in range(3)]

        R = FunctionSpace(self.mesh, "Real", 0)
        pdim = self.params.pdim
        #p_out_coeffs = [Constant(0.0) for i in range(pdim)]
        p_out_coeffs = [Function(R) for i in range(pdim)]

        return u0, p_out_coeffs

    def initial_conditions(self, V, Q, controls=None):
        u0, p_out_coeffs = controls
        # Ignoring p_out_coeffs, returning u0 as initial condition

        p0 = Expression("-beta * x[0] * 0.3", beta=1.0)
        p0.beta = self.params.beta

        return (u0, p0)

    def boundary_conditions(self, V, Q, t, controls=None):
        """Return boundary conditions.

        Returns (bcu, bcp) on the format:

          bcu = [([u0, u1, u2], domainid), ...]
          bcp = [(p, domainid), ...]
        """

        # Create no-slip boundary condition for velocity
        g_noslip = [c0, c0, c0]
        bcu = [
            (g_noslip, 0),
            ]

        if 0:
            u0, p_out_coeffs = controls
            # FIXME: Express p1 in terms of p_out_coeffs

        # Create boundary conditions for pressure
        if 0:
            # Works with and without penalty formulation:
            # Expression version, time is kept updated automatically through reference to Constant t:
            p1 = Expression("(-beta * length) * (minflow + (1.0-minflow)*pow(sin(t*period*DOLFIN_PI),2))",
                            beta=1.0, length=10.0, minflow=0.5, period=1.0, t=t)
            p1.length = self.length
            p1.minflow = 0.3
            p1.beta = self.params.beta
            p1.period = self.params.period
        else:
            # Works only with penalty formulation:
            # UFL expression version, time is kept updated automatically through reference to Constant t:
            p1 = (-self.params.beta * self.length) * (0.3 + 0.7*sin(t*self.params.period*pi)**2)

        bcp = [
            (c0, 1),
            (p1, 2),
            ]

        return (bcu, bcp)

    def update_boundary_conditions(self, V, Q, t):
        """Update functions returned by boundary_conditions.

        Called every timestep in scheme before applying the bcs
        returned by boundary_conditions().

        If the bc functions are stationary or reference the
        time constant properly this implementation may be empty.
        """
        pass

    def observation(self, V, t):
        # Quadratic profile
        zx = ("(beta/(4*nu))*(r*r-x[1]*x[1]-x[2]*x[2])", "0.0", "0.0")
        zx = Expression(zx, beta=0.0, nu=1.0, r=0.5)
        zx.beta = self.params.beta
        zx.nu = self.params.mu / self.params.rho
        zx.r = self.radius

        # Transient pulse
        minflow = 0.3
        zt = Expression("minflow + (1.0-minflow)*pow(sin(2*DOLFIN_PI*t/period),2)",
                        minflow=minflow, t=t, period=self.params.T)

        # Set observation to pulsating quadratic profile
        z = zt*zx

        return z

    def J(self, V, Q, t, u, p, controls):
        # Interpret controls argument
        u0, p_out_coeffs = controls
        u0 = as_vector(u0)
        p_out_coeffs = as_vector(p_out_coeffs)

        # Control is on outflow boundary which is 2 in this problem
        dss = self.ds(2)

        # Define distance functional
        z = self.observation(V, t)
        Jdist = (u0 - z)**2*self.dx()*dt

        # Setup priors
        p_out_coeffs_prior = as_vector([0.0 for i in xrange(len(p_out_coeffs))])
        p_out_coeffs_shifted = as_vector([p_out_coeffs[i+1] for i in xrange(len(p_out_coeffs)-1)] + [p_out_coeffs[0]])
        u0_prior = 0.0*u0

        # Define regularization functional
        alpha = Constant(self.params.alpha)
        Jreg = (
            + alpha * (p_out_coeffs-p_out_coeffs_prior)**2
            #+ alpha * (p_out_coeffs_shifted-p_out_coeffs)**2
            + alpha * (u0-u0_prior)**2
            #+ alpha * grad(u0-u0_prior)**2
            ) * dss*dt[START_TIME]

        return Jdist + Jreg

if __name__ == "__main__":
    p = Pipe()
    show_problem(p)
