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
        mesh = Mesh("../data/pipe_0.2.xml.gz")
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

    def standard_boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        g_noslip = [c0, c0, c0]
        bcu = [
            (g_noslip, 0),
            ]

        # Create boundary conditions for pressure
        p1 = (-self.params.beta * self.length) * (0.3 + 0.7*sin(t*self.params.period*pi)**2)
        bcp = [
            (c0, 1),
            (Constant(p1), 2),
            ]

        return (bcu, bcp)

    def xdirichlet_boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        g_noslip = [c0, c0, c0]
        bcu = [
            (g_noslip, 0),
            ]

        # Create boundary conditions for pressure
        bcp = [
            (c0, 1),
            ]

        return (bcu, bcp)

    boundary_conditions = standard_boundary_conditions
    #boundary_conditions = dirichlet_boundary_conditions

    def xpenalty_boundary_conditions(self, V, Q, t):
        bcu = []

        # Create boundary conditions for pressure
        p1 = (-self.params.beta * self.length) * (0.3 + 0.7*sin(t*self.params.period*pi)**2)
        bcp = [
            (p1, 2),
            ]

        return (bcu, bcp)

if __name__ == "__main__":
    p = Pipe()
    show_problem(p)
