#!/usr/bin/env python

import dolfin_adjoint as da

from collections import namedtuple
from cbcflow import *
from cbcflow.dol import *
from cbcpost import ParamDict, PostProcessor

import numpy as np

# TODO: Add types like these to cbcflow to document BC interface?
InitialConditions = namedtuple("InitialConditions", ["icu", "icp"])
VelocityBC = namedtuple("VelocityBC", ["functions", "region", "method"])
PressureBC = namedtuple("PressureBC", ["function", "region"])
BoundaryConditions = namedtuple("BoundaryConditions", ["bcu", "bcp"])

# Create some types to make tuple usage safer and more readable
Observations = namedtuple("Observations", [])
Controls = namedtuple("Controls", ["icu", "uin"])
CostFunctionals = namedtuple("CostFunctionals", [])
VelocityBoundaryConditions = namedtuple("VelocityBoundaryConditions", ["inflow", "noslip"])
PressureBoundaryConditions = namedtuple("PressureBoundaryConditions", ["outflow"])


def create_geometry(M, N, LENGTH, RADIUS):
    # Create mesh
    mesh = UnitSquareMesh(M, N)
    x = mesh.coordinates()[:,0]
    y = mesh.coordinates()[:,1]
    x = LENGTH*x
    y = RADIUS*2*(y - 0.5)
    mesh.coordinates()[:,0] = x
    mesh.coordinates()[:,1] = y

    # Define boundary markers
    boundary_ids = ParamDict(
        wall = 0,
        left = 1,
        right = 2,
        undefined = 3,
        )

    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0] < DOLFIN_EPS

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0] > LENGTH*(1.0 - DOLFIN_EPS)

    # Create boundary markers
    facet_domains = FacetFunction("size_t", mesh)
    facet_domains.set_all(boundary_ids.undefined)
    DomainBoundary().mark(facet_domains, boundary_ids.wall)
    Left().mark(facet_domains, boundary_ids.left)
    Right().mark(facet_domains, boundary_ids.right)

    return mesh, facet_domains, boundary_ids


class Womersley2D(NSProblem):
    "2D pipe test problem with known stationary analytical solution."

    def __init__(self, params=None):
        NSProblem.__init__(self, params)

        # Create mesh
        refinements = [4, 8, 16, 32, 64]
        N = refinements[self.params.refinement_level]
        M = int(N * self.params.length / (2 * self.params.radius) + 0.5)
        mesh, facet_domains, self.boundary_ids = create_geometry(M, N, self.params.length, self.params.radius)

        # Store mesh and markers
        self.initialize_geometry(mesh, facet_domains=facet_domains)

        # Initialize Womersley solution
        self.init_womersley()

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()
        params.replace(
            # Time parameters
            T=None,
            dt=1e-2,
            period=0.8,
            num_periods=1.0,
            # Physical parameters
            rho = 1.0,
            mu=1.0/30.0,
            )
        params.update(
            # Spatial parameters
            refinement_level=2,
            length=10.0,
            radius=0.5,
            # Analytical solution parameters
            Q=1.0,
            Qfloor=0.25,
            Qpeak=1.0,
            transient_Q=True,
            num_womersley_coefficients=25,
            )
        return params

    def init_womersley(self):
        # Get parameters
        nu = self.params.mu / self.params.rho
        P = self.params.period
        Q = self.params.Q
        Qfloor = self.params.Qfloor
        Qpeak = self.params.Qpeak

        if self.params.transient_Q:
            # Setup coefficients for a transient flow rate
            time_values = np.linspace(0.0, P)
            t = np.mod((P - time_values) / P, P)
            time_profile = Qfloor + (Qpeak - Qfloor) * np.sin(pi*t**3)**2
            self.Q_coeffs = zip(time_values, Q * time_profile)
        else:
            # Setup coefficients for a constant flow rate
            self.Q_coeffs = [(0.0, Q), (1.0, Q)]

        # Create womersley objects
        self.womersley = make_womersley_bcs(self.Q_coeffs, self.mesh, self.boundary_ids.left, nu, None, self.facet_domains,
                                            "Q", num_fourier_coefficients=self.params.num_womersley_coefficients)

    def observations(self, spaces, t):
        return []

    def controls(self, spaces):
        d = spaces.d
        icu = [Function(spaces.U) for i in range(d)]
        uin = [] #[[Function(spaces.U) for i in range(d)] for timestep in [0]]
        return Controls(icu, uin)

    def cost_functionals(self, spaces, t, observations, controls):
        return CostFunctionals()

    def initial_conditions(self, spaces, controls):
        d = spaces.d

        # Get/create IC functions
        icu = controls.icu
        icp = Function(spaces.Q)

        # Interpolate womersley solution into entire initial condition (only possible for this simple pipe case!)
        for i in range(d):
            self.womersley[i].set_t(0.0)
            icu[i].interpolate(self.womersley[i])

        return InitialConditions(icu, icp)

    def boundary_conditions(self, spaces, u, p, t, controls):
        d = len(u)

        # Create inflow BC control functions to be returned and used in scheme forms
        uin = [Function(spaces.U) for i in range(d)]
        inflow = VelocityBC(uin, self.boundary_ids.left, "nietche")

        # Create no-slip bcs
        u0 = [Constant(0.0) for i in range(d)]
        noslip = VelocityBC(u0, self.boundary_ids.wall, "strong")

        # Create outflow bcs for pressure
        outflow = PressureBC(Constant(0.0), self.boundary_ids.right)

        # Return bcs in two lists
        bcu = VelocityBoundaryConditions(inflow, noslip)
        bcp = PressureBoundaryConditions(outflow)
        return BoundaryConditions(bcu, bcp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls, cost_functionals):
        d = len(u)

        # Update time in boundary condition expressions
        for component in self.womersley:
            component.set_t(t)

        # Create new BC control functions at time t
        uin = [Function(spaces.U) for i in range(d)]
        controls.uin.append(uin)

        # Update (initial guess for) controls at time t
        bc = bcs.bcu.inflow
        for i in range(d):
            self.womersley[i].set_t(t)

            dbc = DirichletBC(spaces.U, self.womersley[i], self.facet_domains, bc.region)
            dbc.apply(uin[i].vector())

            # Assign current control function values to the BC functions
            bc.functions[i].assign(uin[i])

        # TODO: Update observations
        #observations.z ...
        # TODO: Add contribution to cost functionals at time t
        #cost_functionals.J += (u-z**2)*dx
        #... + regularization alpha*g**2*dsr + alpha*grad(g)**2*dsr

def main():
    #parameters["adjoint"]["stop_annotating"] = True

    set_log_level(100)

    # Setup problem
    problem = Womersley2D(
        ParamDict(
            dt=1e-3,
            T=0.01,#8,
            num_periods=None,
            refinement_level=2,
            )
        )

    # Setup scheme
    scheme = CoupledPicard(
        ParamDict(
            # BC params
            nietche=ParamDict(
                enable=True,
                formulation=1,
                stabilize=True,
                gamma=100.0,
                ),

            # Variational formulation params
            scale_by_dt=True,
            enable_convection=True, # False = Stokes

            # Nonlinear solver params
            picard_newton_fraction=1.0, # 0.0 = Picard, 1.0 = Newton, (0.0,1.0) = mix
            nonlinear_solver=ParamDict(
                newton_solver=ParamDict(
                    report=True,
                    ),
                ),

            # Form compiler params
            form_compiler_parameters=ParamDict(
                quadrature_degree="auto",
                ),
            )
        )

    # Create casedir name
    params_string = "__".join("{}_{}".format(k, scheme.params.nietche[k]) for k in scheme.params.nietche)
    equation = "NavierStokes" if scheme.params.enable_convection else "Stokes"
    casedir = "results_demo_%s_%s_%s_%s" % (problem.shortname(), scheme.shortname(), equation, params_string)

    # Setup postprocessor
    postprocessor = PostProcessor(
        ParamDict(
            casedir=casedir,
            )
        )
    plot_and_save = dict(plot=True, save=True)
    postprocessor.add_fields([
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        LocalCfl(plot_and_save),
        ])

    # Setup and run solver
    solver = NSSolver(problem, scheme, postprocessor,
        ParamDict(
            timer_frequency=1,
            check_memory_frequency=1,
            enable_annotation=False,
            )
        )

    # Step through forward simulation
    for data in solver.isolve():
        print "TIMESTEP:", data.timestep

    # Try replaying with dolfin-adjoint
    success = da.replay_dolfin()
    print success


if __name__ == "__main__":
    main()
