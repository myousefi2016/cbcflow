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

        # Mesh and known mesh properties
        self.mesh = Mesh("cylinder_0.2.xml.gz")
        self.length = 10.0
        self.radius = 0.5

    @classmethod
    def default_user_params(cls):
        params = ParamDict(
            # Time parameters
            T=0.5,
            dt=0.01,

            # Physical parameters
            rho=1.0,
            mu=0.035,

            beta=5.0,
            )
        return params

    def initial_conditions(self, V, Q):
        u0 = [c0, c0, c0]
        p0 = c0 #Expression("1 - x[0]")
        return (u0, p0)

    def boundary_conditions(self, V, Q, t):
        # Create no-slip boundary condition for velocity
        g_noslip = [c0, c0]
        bcu = [
            (g_noslip, 0),
            ]

        # Create boundary conditions for pressure
        bcp = [
            (c0, 1),
            (Constant(-self.params.beta*self.length*t), 2),
            ]

        return (bcu, bcp)

if __name__ == "__main__":
    p = Pipe()
    show_problem(p)
