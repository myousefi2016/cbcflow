#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *

c0 = Constant(0.0)

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

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-3,
            period=0.8,
            num_periods=1.0,

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

def main():
    problem = Pipe()
    scheme = IPCS_Stable()

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    plot_and_save = dict(plot=True, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        ]
    postproc = NSPostProcessor({"casedir": casedir})
    postproc.add_fields(fields)

    solver = NSSolver(problem, scheme, postproc)
    solver.solve()

if __name__ == "__main__":
    main()
