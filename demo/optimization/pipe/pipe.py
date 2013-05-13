#!/usr/bin/env python
__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-04-24"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

c0 = Constant(0, name="zero")

from headflow.core.utils import as_scalar_space

class Problem(NSProblem):
    "3D pipe optimization test problem with known analytical solution."

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
        jp = ParamDict(
            # Regularization strength
            alpha=1e-4,

            # Regularization term toggles for initial velocity
            alpha_u_prior = 1,
            alpha_u_div   = 0,
            alpha_u_grad  = 1,
            alpha_u_wall  = 1,

            # Regularization term toggles for pressure bcs
            alpha_p_prior   = 1,
            alpha_p_shifted = 0,
            alpha_p_basis   = 0,
            alpha_p_dt      = 0,

            # Toggle cyclic term in cost functional
            cyclic = 0,
            )
        params = ParamDict(
            # Physical parameters
            rho=1.0,
            mu=0.035,

            # Time parameters
            T=None, # Computed
            dt=1e-3,
            period=0.8,

            # Time to run simulation, periods is used if timesteps are not given
            num_periods=2,
            num_timesteps=0,

            # Control parameters
            pdim=1,

            # Regularization parameters
            J=jp,
            )
        return params

    def controls(self, V, Q):
        # Velocity initial condition control
        V = as_scalar_space(V)
        u0 = [Function(V, name="ui_%d"%i) for i in range(V.cell().d)]
        for u0c in u0:
            u0c.interpolate(Expression("0.0"))

        # Coefficients for pressure bcs
        p_out_coeffs = [Constant(0.0, name="pc%d"%i) for i in range(self.params.pdim)]

        return u0, p_out_coeffs

    def initial_conditions(self, V, Q, controls):
        # Extract initial conditions from controls
        u0, p_out_coeffs = controls

        # Pressure initial condition control # TODO: Does not seem to matter, remove!
        p0 = Function(Q, name="p0")
        p0e = Expression("0.0")
        p0.interpolate(p0e)

        return (u0, p0)

    def pressure_basis(self, t):
        # TODO: Configurable basis
        if 1: # Fourier basis
            n = (self.params.pdim-1)//2
            assert self.params.pdim == 2*n+1
            omega = Constant(self.params.period*pi)
            return (
                [1.0]
                + [sin(omega*k*t) for k in xrange(1,n+1)]
                + [cos(omega*k*t) for k in xrange(1,n+1)]
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

        # Create boundary conditions for pressure expressed in terms of p_out_coeffs controls
        u0, p_out_coeffs = controls
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
        zx = ("upeak*(r*r-x[1]*x[1]-x[2]*x[2])", "0.0", "0.0")
        zx = Expression(zx, upeak=1.0, r=0.5)
        zx.r = self.radius
        zx.upeak = 1.0

        # Transient pulse, NB! this references the actual Constant t passed here!
        zt = Expression("minflow + (1.0-minflow)*pow(sin(2*DOLFIN_PI*t/period),2)",
                        minflow=0.5, t=t, period=1.0)
        zt.period = self.params.period
        zt.minflow = 0.3

        # Set observation to pulsating quadratic profile
        z = zt*zx

        return z

    def J(self, V, Q, t, u, p, controls):
        # Interpret controls argument
        u0, p_out_coeffs = controls
        u0 = as_vector(u0)
        p_out_coeffs = as_vector(p_out_coeffs)

        # Get parameters
        jp = self.params.J
        alpha = Constant(jp.alpha, name="alpha")

        # Pressure control is on outflow boundary which is 2 in this problem
        control_boundaries = 2
        dsc = self.ds(control_boundaries)
        wall_boundaries = 0
        dsw = self.ds(wall_boundaries)

        # Define distance from observation
        z = self.observation(V, t)
        Jdist = (u - z)**2*self.dx()*dt

        # Add cyclic distance functional
        Jdist += jp.cyclic * (u - u0)**2*self.dx()*dt[FINAL_TIME]

        # Setup priors
        u0_prior = 0.0*u0
        p_out_coeffs_prior = as_vector([0.0 for i in xrange(len(p_out_coeffs))])

        # A couple of ways to penalize dp1/dt
        p_out_coeffs_shifted = as_vector([p_out_coeffs[i+1] for i in xrange(len(p_out_coeffs)-1)] + [p_out_coeffs[0]])
        t = variable(t)
        p1 = sum(p_out_coeffs[k] * N for k,N in enumerate(self.pressure_basis(t)))
        p1_t = diff(p1, t)

        # Regularization for initial velocity
        Jreg = (
            # Penalize initial velocity everywhere
            + alpha * jp.alpha_u_prior * (u0-u0_prior)**2
            + alpha * jp.alpha_u_div   * div(u0)**2
            + alpha * jp.alpha_u_grad  * grad(u0-u0_prior)**2
            ) * dx*dt[START_TIME]
        Jreg += (
            # Penalize initial velocity hard to be zero on walls
            + jp.alpha_u_wall * u0**2
            ) * dsw*dt[START_TIME]

        # Regularization for boundary conditions
        Jreg += (
            # Penalize time dependent pressure control
            + alpha * jp.alpha_p_prior   * (p_out_coeffs-p_out_coeffs_prior)**2
            + alpha * jp.alpha_p_shifted * (p_out_coeffs_shifted-p_out_coeffs)**2
            + alpha * jp.alpha_p_basis   * (p1)**2
            + alpha * jp.alpha_p_dt      * (p1_t)**2
            ) * dsc*dt[START_TIME]

        return Jdist + Jreg

if __name__ == "__main__":
    p = Problem()
    show_problem(p)
