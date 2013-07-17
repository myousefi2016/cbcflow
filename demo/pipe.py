#!/usr/bin/env python
__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-04-24"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from headflow import *
from headflow.dol import *

c0 = Constant(0.0)

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

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-3,
            period=0.8,
            num_periods=3,

            # Physical parameters
            rho=1.0,
            mu=0.035,
            )
        params.update(
            # Pressure gradient amplitude
            beta=5.0,
            )
        return params

    def initial_conditions(self, spaces, controls):
        u0 = [c0, c0, c0]
        p0 = Expression("-beta * x[0] * 0.3", beta=1.0)
        p0.beta = self.params.beta
        return (u0, p0)

    def _pressure_drop(self, t):
        e = Expression("-beta * length * (0.3 + 0.7 * sin(t*period*DOLFIN_PI) * sin(t*period*DOLFIN_PI))",
                        beta=1.0, length=10.0, t=t, period=1.0)
        e.beta = self.params.beta
        e.length = self.length
        e.period = self.params.period
        return e

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip boundary condition for velocity
        g_noslip = [c0, c0, c0]
        bcu = [(g_noslip, 0)]

        # Create boundary conditions for pressure
        p1 = self._pressure_drop(t)
        bcp = [(c0, 1), (p1, 2)]

        return (bcu, bcp)

if __name__ == "__main__":
    from demo_main import demo_main
    demo_main(Pipe)
