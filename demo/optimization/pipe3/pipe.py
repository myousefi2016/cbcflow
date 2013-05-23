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

        self.beta = 4*(self.params.mu/self.params.rho) / self.radius**2
        self.omega = Constant(self.params.period*pi, name="omega")

        # Set end time based on period and number of periods NB! Overrides given T!
        if self.params.num_timesteps:
            self.params.T = self.params.dt * self.params.num_timesteps
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
            alpha_u_curl  = 0,
            alpha_u_grad  = 0,
            alpha_u_wall  = 0,

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
            pdim=3,

            # Regularization parameters
            J=jp,
            )
        return params

    def pressure_basis(self, t): # TODO: Configurable basis
        if 1: # Fourier basis
            n = (self.params.pdim-1)//2
            assert self.params.pdim == 2*n+1
            omega = self.omega
            basis = (
                [1.0]
                + [sin(omega*k*t) for k in xrange(1,n+1)]
                + [cos(omega*k*t) for k in xrange(1,n+1)]
                )
            print "BASIS:", basis
            return basis
        elif 0: # Polynomials?
            pass

    def x_observation_expression(self, t):
        # Quadratic profile times a transient pulse
        ze = "(minflow + (1.0-minflow)*pow(sin(2.0*DOLFIN_PI*t/period),2)) * upeak * (r*r-x[1]*x[1]-x[2]*x[2]) / (r*r)"
        ze = Expression(ze, upeak=1.0, r=0.5, period=1.0, minflow=0.3, t=0.0, name="ze")
        ze.r = self.radius
        ze.period = self.params.period
        ze.minflow = 0.3
        ze.upeak = 1.0
        ze.t = float(t)
        return ze

    def _observation_expression(self, t):
        # STATIONARY SOLUTION FOR DEBUGGING
        ze = "upeak * (r*r-x[1]*x[1]-x[2]*x[2]) / (r*r)"
        ze = Expression(ze, upeak=1.0, r=0.5, name="ze")
        return ze

    def observations(self, spaces, t):
        "Return a list of observation functions that may need updating each timestep."
        U = spaces.U
        d = spaces.d

        z = as_vector([Function(U, name="z_%d" % i) for i in xrange(d)])
        ze = self._observation_expression(t)
        z[0].interpolate(ze)

        # Return observations tuple
        observations = (z,)
        return observations

    def velocity_observation(self, spaces, t, observations):
        z, = observations
        return z

    def update_observations(self, spaces, t, observations):
        "Update functions in list returned by auxilliary_functions() for this timestep."
        # Interpolate ze into x-component of z
        z, = observations
        ze = self._observation_expression(t)
        z[0].interpolate(ze)

    def set_controls(self, controls):
        "Set controls to be returned next time controls() is called."
        self._controls = controls

    def controls(self, spaces):
        U = spaces.U
        d = spaces.d

        # Velocity initial condition control
        u0 = [Function(U, name="ui_%d"%i) for i in xrange(d)]

        # Coefficients for pressure bcs
        p_out_coeffs = [Constant(0.0, name="pc%d"%i) for i in xrange(self.params.pdim)]
        #print '\n'.join([str(pc) for pc in p_out_coeffs])
        #crash

        if hasattr(self, "_controls"):
            # Set given control values
            m_u = self._controls[:d]
            m_p = self._controls[d:]
            for i, u0c in enumerate(u0):
                u0c.interpolate(m_u[i])
            for i, pc in enumerate(p_out_coeffs):
                pc.assign(m_p[i])
        else:
            # Set initial control values
            # STATIONARY SOLUTION FOR DEBUGGING
            ze = "upeak * (r*r-x[1]*x[1]-x[2]*x[2]) / (r*r)"
            ze = Expression(ze, upeak=1.0, r=0.5, name="ze")
            u0[0].interpolate(ze)
            for i in range(1,d):
                u0[i].interpolate(Expression("0.0"))
                p_out_coeffs[0].assign(-self.length*self.beta)

        # Return controls tuple
        controls = (u0, p_out_coeffs)
        return controls

    def initial_conditions(self, spaces, controls):
        # Extract initial conditions from controls
        u0, p_out_coeffs = controls

        # Pressure initial condition control
        p0 = Function(spaces.Q, name="pinit")
        p0e = Expression("-beta*(x[0]+length)", beta=self.beta, length=self.length)
        p0.interpolate(p0e)

        # Return ics tuple
        ics = (u0, p0)
        return ics

    def boundary_conditions(self, spaces, u, p, t, controls):
        """Return boundary conditions.

        Returns (bcu, bcp) on the format:

          bcu = [([u0, u1, u2], domainid), ...]
          bcp = [(p, domainid), ...]
        """
        d = spaces.d

        # Create no-slip boundary condition for velocity
        bcu = [
            ([c0]*d, 0),
            ]

        # Create boundary conditions for pressure expressed in terms of p_out_coeffs controls
        u0, p_out_coeffs = controls
        p1 = sum(p_out_coeffs[k] * N for k,N in enumerate(self.pressure_basis(t)))

        bcp = [
            (c0, 2),
            (p1, 1),
            ]

        # Return bcs tuple
        bcs = (bcu, bcp)
        return bcs

    def update_boundary_conditions(self, spaces, u, p, t, bcs):
        """Update functions returned by boundary_conditions.

        Called every timestep in scheme before applying the bcs
        returned by boundary_conditions().

        If the bc functions are stationary or reference the
        time constant properly this implementation may be empty.
        """
        # TODO: Update bcs?
        #p0, p1 = bcs[1][0][0], bcp[1][1][0]
        pass

    def J(self, spaces, t, u, p, controls, observations):

        # Interpret controls argument
        u0, p_out_coeffs = controls
        u0 = as_vector(u0)
        p_out_coeffs = as_vector(p_out_coeffs)

        # Get parameters
        jp = self.params.J
        alpha = Constant(jp.alpha, name="alpha")

        # Wall is 0 in this problem
        wall_boundaries = 0
        dsw = self.ds(wall_boundaries)

        # Pressure control is on outflow boundary which is 2 in this problem
        control_boundaries = 2
        dsc = self.ds(control_boundaries)

        # Define distance from observation
        z = self.velocity_observation(spaces, t, observations)
        Jdist = (u - z)**2*self.dx()*dt

        # TODO: Use this for scaling of J?
        #z2 = assemble(z**2*self.dx(), mesh=self.mesh, annotate=False)
        #u2 = assemble(u**2*self.dx(), mesh=self.mesh, annotate=False)

        # Add cyclic distance functional
        #u02 = as_vector([Function(u0c) for u0c in u0])
        #Jdist += jp.cyclic * (u - u02)**2*self.dx()*dt[FINISH_TIME]

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
            + alpha * jp.alpha_u_curl  * curl(u0)**2
            + alpha * jp.alpha_u_grad  * grad(u0-u0_prior)**2
            ) * dx*dt[START_TIME]
        Jreg += (
            # Penalize initial velocity hard to be zero on walls
            + alpha * jp.alpha_u_wall * u0**2
            ) * dsw*dt[START_TIME]

        # Regularization for boundary conditions
        Jreg += (
            # Penalize time dependent pressure control
            + alpha * jp.alpha_p_prior   * (p_out_coeffs-p_out_coeffs_prior)**2
            + alpha * jp.alpha_p_shifted * (p_out_coeffs_shifted-p_out_coeffs)**2
            + alpha * jp.alpha_p_basis   * (p1)**2
            + alpha * jp.alpha_p_dt      * (p1_t)**2
            ) * dsc*dt[START_TIME]

        # Add distance and regularization to get total cost functional
        Jtot = Jdist + Jreg

        if 0:
            Q = spaces.Q
            g = Function(Q, name='dummy')
            vol = assemble((g**2 + Constant(1.0))*dx(), annotate=False)
            Jtot = (1.0/vol)*(g**2 + t)*dx()*dt[FINISH_TIME] + -(1.0/vol)*(g**2 + t)*dx()*dt[START_TIME]
            Jtot = (1.0/vol)*(g**2 + t)*dx()*dt

            Jtot = (1.0/vol)*(g**2 + t/0.0105)*dx()*dt

        if 0:
            print
            print
            print float(t)
            print Jtot
            print

        return Jtot

if __name__ == "__main__":
    p = Problem()
    show_problem(p)
