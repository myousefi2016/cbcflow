#!/usr/bin/env python
__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-04-24"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

c0 = Constant(0)

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
            )
        return params

    def initial_conditions(self, V, Q):
        u0 = [c0, c0, c0]
        p0 = Expression("-beta * x[0] * 0.3", beta=1.0)
        p0.beta = self.params.beta
        return (u0, p0)

    def boundary_conditions(self, V, Q, t):
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

if __name__ == "__main__":
    p = Pipe()
    show_problem(p)
