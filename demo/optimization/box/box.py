#!/usr/bin/env python
__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-04-24"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

c0 = Constant(0, name="zero")

from headflow.core.utils import as_scalar_space

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1.0-DOLFIN_EPS and on_boundary

class Box(NSProblem):
    "3D box test problem with known analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Build cube mesh
        n = 3
        mesh = UnitCubeMesh(n, n, n)
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(3)
        DomainBoundary().mark(facet_domains, 0)
        Left().mark(facet_domains, 1)
        Right().mark(facet_domains, 2)

        # Canonical initialization of geometry
        self.initialize_geometry(mesh, facet_domains=facet_domains)

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
            pdim=1,
            )
        return params

    def controls(self, V, Q):
        # Velocity initial condition control
        V = as_scalar_space(V)
        u0 = [Function(V, name="ui_%d"%i) for i in range(V.cell().d)]

        # Pressure initial condition control
        p0 = Function(Q, name="p0")

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

        # Create boundary conditions for pressure, expressed in terms of p_out_coeffs controls
        u0, p0, p_out_coeffs = controls
        p1 = sum(p_out_coeffs[k] * N for k,N in enumerate(self.pressure_basis(t)))

        # Collect pressure bcs, fixed to 0 at one side
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
        # Flat profile, unit flow rate
        zx = ("1.0", "0.0", "0.0")
        zx = Expression(zx)

        # Transient pulse
        minflow = 0.3
        zt = Expression("minflow + (1.0-minflow)*pow(sin(2*DOLFIN_PI*t/period),2)",
                        minflow=minflow, t=t, period=self.params.period)

        # Set observation to pulsating flat profile
        z = zt*zx

        return z

    def J(self, V, Q, t, u, p, controls):
        # Interpret controls argument
        u0, p0, p_out_coeffs = controls
        u0 = as_vector(u0)
        p_out_coeffs = as_vector(p_out_coeffs)

        # Control is on outflow boundary which is 2 in this problem
        dss = self.ds(2)

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
            + alpha * (u0-u0_prior)**2
            #+ alpha * div(u0)**2
            #+ alpha * grad(u0-u0_prior)**2

            + alpha * (p0-p0_prior)**2
            #+ alpha * grad(p0-p0_prior)**2

            + alpha * (p_out_coeffs-p_out_coeffs_prior)**2
            #+ alpha * (p_out_coeffs_shifted-p_out_coeffs)**2
            #+ alpha * (p1)**2
            #+ alpha * (p1_t)**2
            ) * dss*dt[START_TIME]

        return Jdist + Jreg

if __name__ == "__main__":
    p = Box()
    show_problem(p)
