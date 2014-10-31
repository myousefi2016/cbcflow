#!/usr/bin/env python

from cbcflow import *
from cbcflow.dol import *
from cbcpost import ParamDict, PostProcessor

import numpy as np

LENGTH = 10.0
RADIUS = 0.5

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > LENGTH*(1.0 - DOLFIN_EPS)

class Poiseuille2D(NSProblem):
    "2D pipe test problem with known stationary analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh
        refinements = [4,8,16,32,64]
        N = refinements[self.params.refinement_level]
        M = int(N*LENGTH/(2*RADIUS) + 0.5)
        mesh = UnitSquareMesh(M, N)
        x = mesh.coordinates()[:,0]
        y = mesh.coordinates()[:,1]
        x = LENGTH*x
        y = RADIUS*2*(y - 0.5)
        mesh.coordinates()[:,0] = x
        mesh.coordinates()[:,1] = y

        # We will apply markers with these id values
        self.wall_boundary_id = 0
        self.left_boundary_id = 1
        self.right_boundary_id = 2
        self.undefined_boundary_id = 3

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(self.undefined_boundary_id)
        DomainBoundary().mark(facet_domains, self.wall_boundary_id)
        Left().mark(facet_domains, self.left_boundary_id)
        Right().mark(facet_domains, self.right_boundary_id)

        # Setup analytical solution constants
        Q = self.params.Q
        nu = self.params.mu / self.params.rho
        self.alpha = 2.0 * Q / (pi * RADIUS**4)
        self.beta = 2.0 * nu * self.alpha

        # Toggle to test using Poiseuille-shaped bcs with transient flow rate
        if 1:
            #print "Using stationary bcs. Analytical solution should hold."
            #print "Expected peak velocity:", self.alpha * RADIUS**2
            #print "Expected total pressure drop:", self.beta * LENGTH
            self.Q_coeffs = [(0.0, Q), (1.0, Q)]
        else:
            print "Using transient bcs. Analytical solution will not hold."
            T = self.params.T
            P = self.params.period
            tvalues = np.linspace(0.0, P)
            Qvalues = Q * (0.3 + 0.7*np.sin(pi*((P-tvalues)/P)**2)**2)
            self.Q_coeffs = zip(tvalues, Qvalues)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-2,
            period=0.8,
            num_periods=0.1,
            # Physical parameters
            rho = 10.0,
            mu=1.0/30.0,
            )
        params.update(
            # Spatial parameters
            refinement_level=0,
            # Analytical solution parameters
            Q=1.0,
            )
        return params

    def analytical_solution(self, spaces, t):
        A = 2*RADIUS
        Q = self.params.Q
        mu = self.params.mu

        dpdx = 3*Q*mu/(A*RADIUS**2)

        ux = Expression("(radius*radius - x[1]*x[1])/(2*mu)*dpdx", radius=RADIUS, mu=mu, dpdx=dpdx, degree=2)
        uy = Constant(0.0)

        u = [ux, uy]

        p = Expression("dpdx * (length-x[0])", dpdx=dpdx, length=LENGTH)

        return (u, p)

    def observations(self, spaces, t):
        return []

    def controls(self, spaces):
        return []

    def cost_functionals(self, spaces, t, observations, controls):
        return []

    def initial_conditions(self, spaces, controls):
        return self.analytical_solution(spaces, 0.0)

    def boundary_conditions(self, spaces, u, p, t, controls):
        # Create no-slip bcs
        d = len(u)
        u0 = [Constant(0.0) for i in range(d)]
        noslip = (u0, self.wall_boundary_id, "strong")

        # Create Poiseuille inflow bcs
        uin = make_poiseuille_bcs(self.Q_coeffs, self.mesh, self.left_boundary_id, self.params.Q, self.facet_domains)
        for ucomp in uin:
            ucomp.set_t(t)
        inflow = (uin, self.left_boundary_id, "nietche")

        # Create outflow bcs for pressure
        outflow = (Constant(0.0), self.right_boundary_id, "natural")

        # Return bcs in two lists
        bcu = [noslip, inflow]
        bcp = [outflow]
        return (bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls, cost_functionals):
        # TODO: Drop this
        # Update time in boundary condition expressions
        bcu, bcp = bcs
        noslip, inflow = bcu
        uin = inflow[0]
        for ucomp in uin:
            ucomp.set_t(t)

        # TODO: Update observations

        # TODO: Update initial guess for controls at time t
        if controls:
            m0 = controls[0]
            V = m0.function_space()
            m = Function(V)
            #m.interpolate(initial_control_guess)
            controls.append(m)
            m0.assign(m)

        # TODO: Add contribution to cost functionals at time t

def main():
    set_log_level(100)

    problem = Poiseuille2D(
        ParamDict(
            dt=1e-3,
            T=1e-2,
            num_periods=None,
            refinement_level=1,
            )
        )

    scheme = CoupledPicard(
        ParamDict(
            nietche=ParamDict(
                enable=True,
                symmetrize=True,
                stabilize=False,
                gamma=10.0,
                ),
            scale_by_dt=True,
            enable_convection=True,
            )
        )

    casedir = "results_demo_%s_%s" % (problem.shortname(), scheme.shortname())
    plot_and_save = dict(plot=True, save=True)
    fields = [
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        ]
    postproc = PostProcessor({"casedir": casedir})
    postproc.add_fields(fields)

    solver = NSSolver(problem, scheme, postproc)
    solver.solve()


if __name__ == "__main__":
    main()
