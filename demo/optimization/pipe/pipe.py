#!/usr/bin/env python
__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-04-24"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

c0 = Constant(0, name="zero")

from headflow.core.utils import as_scalar_space

class Pipe(NSProblem):
    "3D pipe test problem with known analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Get 3D pipe mesh from file
        mesh = Mesh("../../../data/pipe_0.2.xml.gz")
        self.initialize_geometry(mesh)

        # Known properties of the mesh
        self.length = 10.0
        self.radius = 0.5

        # Set end time based on period and number of periods NB! Overrides given T!
        if self.params.num_timesteps:
            self.params.T = self.params.dt * (self.params.num_timesteps+0.5)
        else:
            self.params.T = self.params.period * self.params.num_periods

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            # Physical parameters
            rho=1.0,
            mu=0.035,

            # Pressure gradient amplitude
            beta=0.15,

            # Time parameters
            T=None, # Computed
            dt=1e-3,
            period=0.8,

            # Time to run simulation, periods is used if timesteps are not given
            num_periods=3,
            num_timesteps=0,

            # Control parameters
            alpha=1e-4,
            pdim=5,
            )
        return params

    def controls(self, V, Q):
        # Velocity initial condition control
        V = as_scalar_space(V)
        u0 = [Function(V, name="ui_%d"%i) for i in range(V.cell().d)]
        for u0c in u0:
            u0c.interpolate(Expression("0.0"))

        # Pressure initial condition control
        p0 = Function(Q, name="p0")
        #p0e = Expression("-beta * x[0] * 0.3", beta=1.0)
        #p0e.beta = self.params.beta
        p0e = Expression("0.0")
        p0.interpolate(p0e)

        # Coefficients for pressure bcs
        p_out_coeffs = [Constant(0.0, name="pc%d"%i) for i in range(self.params.pdim)]

        return u0, p0, p_out_coeffs

    def initial_conditions(self, V, Q, controls):
        # Extract initial conditions from controls
        u0, p0, p_out_coeffs = controls
        return (u0, p0)

    def pressure_basis(self, t):
        # TODO: Configurable basis
        if 1: # Fourier basis
            n = (self.params.pdim-1)//2
            assert self.params.pdim == 2*n+1
            return (
                [1.0]
                + [sin(self.params.period*pi*k*t) for k in xrange(1,n+1)]
                + [cos(self.params.period*pi*k*t) for k in xrange(1,n+1)]
                )
        elif 0: # Polynomials?
            pass

    def boundary_conditions(self, V, Q, t, controls):
        """Return boundary conditions.

        Returns (bcu, bcp) on the format:

          bcu = [([u0, u1, u2], domainid), ...]
          bcp = [(p, domainid), ...]
        """

        # Create no-slip boundary condition for velocity
        bcu = [
            ([c0, c0, c0], 0),
            ]

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

        elif 0:
            # Works only with penalty formulation:
            # UFL expression version, time is kept updated automatically through reference to Constant t:
            p1 = (-self.params.beta * self.length) * (0.3 + 0.7*sin(t*self.params.period*pi)**2)

        else:
            # Expressed in terms of p_out_coeffs controls
            u0, p0, p_out_coeffs = controls
            p1 = sum(p_out_coeffs[k] * N for k,N in enumerate(self.pressure_basis(t)))

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
                        minflow=minflow, t=t, period=self.params.period)

        # Set observation to pulsating quadratic profile
        z = zt*zx

        return z

    def J(self, V, Q, t, u, p, controls):
        # Interpret controls argument
        u0, p0, p_out_coeffs = controls
        u0 = as_vector(u0)
        p_out_coeffs = as_vector(p_out_coeffs)

        # Pressure control is on outflow boundary which is 2 in this problem
        control_boundaries = 2
        dsc = self.ds(control_boundaries)
        wall_boundaries = 0
        dsw = self.ds(wall_boundaries)

        # Define distance functional
        z = self.observation(V, t)
        Jdist = (u - z)**2*self.dx()*dt

        # Define cyclic distance functional
        #Jcyc= (u - u0)**2*self.dx()*dt[FINAL_TIME]

        # Setup priors
        u0_prior = 0.0*u0
        p0_prior = 0.0*p0
        p_out_coeffs_prior = as_vector([0.0 for i in xrange(len(p_out_coeffs))])

        # A couple of ways to penalize dp1/dt
        p_out_coeffs_shifted = as_vector([p_out_coeffs[i+1] for i in xrange(len(p_out_coeffs)-1)] + [p_out_coeffs[0]])
        t = variable(t)
        p1 = sum(p_out_coeffs[k] * N for k,N in enumerate(self.pressure_basis(t)))
        p1_t = diff(p1, t)

        # Define regularization functional
        alpha = Constant(self.params.alpha, name="alpha")
        Jreg = (
            # Penalize initial velocity everywhere
            + alpha * (u0-u0_prior)**2
            #+ alpha * div(u0)**2
            + alpha * grad(u0-u0_prior)**2

            # Penalize initial pressure everywhere
            + alpha * (p0-p0_prior)**2
            #+ alpha * grad(p0-p0_prior)**2

            ) * dx*dt[START_TIME] + (

            # Penalize initial velocity hard to be zero on walls
            + u0**2

            ) * dsw*dt[START_TIME] + (

            # Penalize time dependent pressure control
            + alpha * (p_out_coeffs-p_out_coeffs_prior)**2
            #+ alpha * (p_out_coeffs_shifted-p_out_coeffs)**2
            #+ alpha * (p1)**2
            #+ alpha * (p1_t)**2

            ) * dsc*dt[START_TIME]

        return Jdist + Jreg

if __name__ == "__main__":
    p = Pipe()
    show_problem(p)
