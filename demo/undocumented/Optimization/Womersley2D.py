#!/usr/bin/env python

from cbcpost import ParamDict, PostProcessor, Parameterized

import dolfin_adjoint as da

from cbcflow import *
from cbcflow.dol import *

from collections import namedtuple
import numpy as np
import time

# TODO: Add types like these to cbcflow to document BC interface?
InitialConditions = namedtuple("InitialConditions", ["icu", "icp"])
VelocityBC = namedtuple("VelocityBC", ["functions", "region", "method"])
PressureBC = namedtuple("PressureBC", ["function", "region"])
BoundaryConditions = namedtuple("BoundaryConditions", ["bcu", "bcp"])

# Create some types to make tuple usage safer and more readable
Observations = namedtuple("Observations", [])
Controls = namedtuple("Controls", ["icu", "uin"])
CostFunctionals = namedtuple("CostFunctionals", ["J"])
VelocityBoundaryConditions = namedtuple("VelocityBoundaryConditions", ["inflow", "noslip"])
PressureBoundaryConditions = namedtuple("PressureBoundaryConditions", ["outflow"])


class Geometry(Parameterized):
    def __init__(self, params=None):
        Parameterized.__init__(self, params)

        # Get parameters
        l = self.params.length
        r = self.params.radius
        nr = self.params.refinement_level

        # Compute mesh size
        refinements = [4, 8, 16, 32, 64]
        n = refinements[nr]
        m = int(n * l / (2.0 * r) + 0.5)

        # Create mesh
        mesh = UnitSquareMesh(m, n)
        mesh.coordinates()[:,0] = l * mesh.coordinates()[:,0]
        mesh.coordinates()[:,1] = 2.0 * r * (mesh.coordinates()[:,1] - 0.5)

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
                return on_boundary and x[0] > l*(1.0 - DOLFIN_EPS)

        # Create boundary markers
        facet_domains = FacetFunction("size_t", mesh)
        facet_domains.set_all(boundary_ids.undefined)
        DomainBoundary().mark(facet_domains, boundary_ids.wall)
        Left().mark(facet_domains, boundary_ids.left)
        Right().mark(facet_domains, boundary_ids.right)

        # Store as member data
        self.mesh = mesh
        self.facet_domains = facet_domains
        self.boundary_ids = boundary_ids

    @classmethod
    def default_params(cls):
        params = ParamDict(
            refinement_level=2,
            length=10.0,
            radius=0.5,
            )
        return params


class ProblemBase(NSProblem):
    "2D pipe test problem with known stationary analytical solution."

    def __init__(self, params):
        NSProblem.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSProblem.default_params()

        # Override default time parameters from NSProblem
        params.replace(
            T=None,
            dt=1e-2,
            period=0.8,
            num_periods=1.0,
            )

        # Override default physical parameters from NSProblem
        params.replace(
            rho = 1.0,
            mu=1.0/30.0,
            )

        return params

    def boundary_conditions(self, spaces, u, p, t, controls):
        d = len(u)
        bids = self.geometry.boundary_ids

        # Create inflow BC control functions to be returned and used in scheme forms
        uin = [Function(spaces.U) for i in range(d)]
        inflow = VelocityBC(uin, bids.left, "nietche")

        # Create no-slip bcs
        u0 = [Constant(0.0) for i in range(d)]
        noslip = VelocityBC(u0, bids.wall, "strong")

        # Create outflow bcs for pressure
        outflow = PressureBC(Constant(0.0), bids.right)

        # Return bcs in two lists
        bcu = VelocityBoundaryConditions(inflow, noslip)
        bcp = PressureBoundaryConditions(outflow)
        return BoundaryConditions(bcu, bcp)


class AnalyticProblem(ProblemBase):
    "2D pipe test problem with known stationary analytical solution."

    def __init__(self, params=None):
        ProblemBase.__init__(self, params)
        self.geometry = Geometry(self.params.geometry)
        self.initialize_geometry(self.geometry.mesh, self.geometry.facet_domains)
        self.init_womersley()

    @classmethod
    def default_params(cls):
        params = ProblemBase.default_params()

        # Analytical solution parameters
        womersley = ParamDict(
            Q=1.0,
            Qfloor=0.25,
            Qpeak=1.0,
            transient_Q=True,
            num_womersley_coefficients=25,
            )

        # Geometry parameters
        geometry = Geometry.default_params()

        # Add parameter groups
        params.update(
            geometry=geometry,
            womersley=womersley,
            )

        return params

    def init_womersley(self):
        mesh = self.geometry.mesh
        facet_domains = self.geometry.facet_domains
        bids = self.geometry.boundary_ids

        # Get parameters
        nu = self.params.mu / self.params.rho
        P = self.params.period
        wp = self.params.womersley

        if wp.transient_Q:
            # Setup coefficients for a transient flow rate
            time_values = np.linspace(0.0, P)
            t = np.mod((P - time_values) / P, P)
            time_profile = wp.Qfloor + (wp.Qpeak - wp.Qfloor) * np.sin(pi*t**3)**2
            Q_coeffs = zip(time_values, wp.Q * time_profile)
        else:
            # Setup coefficients for a constant flow rate
            Q_coeffs = [(0.0, wp.Q), (1.0, wp.Q)]

        # Create womersley objects
        self.womersley = make_womersley_bcs(Q_coeffs, mesh, bids.left, nu, None, facet_domains,
                                            "Q", num_fourier_coefficients=wp.num_womersley_coefficients)

    def observations(self, spaces, t):
        # Can be omitted, keeping here for transparency
        return ()

    def controls(self, spaces):
        # Can be omitted, keeping here for transparency
        return ()

    def cost_functionals(self, spaces, t, observations, controls):
        # Can be omitted, keeping here for transparency
        return ()

    def initial_conditions(self, spaces, controls):
        d = spaces.d

        # Create IC functions
        icu = [Function(spaces.U) for i in range(d)]
        icp = Function(spaces.Q) # Is this actually used?

        # Interpolate womersley solution into entire initial condition (only possible for this simple pipe case!)
        for u, w in zip(icu, self.womersley):
            w.set_t(0.0)
            u.interpolate(w)

        return InitialConditions(icu, icp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls, cost_functionals):
        d = len(u)
        facet_domains = self.geometry.facet_domains

        # Update boundary control functions at time t
        bc = bcs.bcu.inflow
        for u, w in zip(bc.functions, self.womersley):
            w.set_t(t)
            dbc = DirichletBC(spaces.U, w, facet_domains, bc.region)
            u.vector().zero()
            dbc.apply(u.vector())


class AssimilationProblem(ProblemBase):
    "2D pipe test problem with known stationary analytical solution."

    def __init__(self, params, geometry, observations):
        ProblemBase.__init__(self, params)
        self.geometry = geometry
        self.initialize_geometry(self.geometry.mesh, self.geometry.facet_domains)

        # Store observations for later
        self._observations = observations

        # Control id tuples to make the below code more generic (deliberately not including wall here)
        self.controlled_boundary_ids = (self.geometry.boundary_ids.left,)
        self.uncontrolled_boundary_ids = (self.geometry.boundary_ids.right,)

    @classmethod
    def default_params(cls):
        params = ProblemBase.default_params()
        J = ParamDict(
            # Initial control regularization
            alpha_u0=0.0,
            alpha_u0_div=0.0,
            alpha_u0_grad=0.0,
            alpha_u0_grad_controlled=0.0,
            alpha_u0_grad_uncontrolled=0.0,
            alpha_u0_wall=0.0,
            # Boundary control regularization
            alpha_g_volume=0.0,
            alpha_g_left=0.0,
            alpha_g_right=0.0,
            alpha_g_grad_volume=0.0,
            alpha_g_grad_left=0.0,
            alpha_g_grad_right=0.0,
            )
        params.update(J=J)
        return params

    def observations(self, spaces, t):
        return self._observations

    def controls(self, spaces):
        d = spaces.d

        # Initial condition control (may choose to use this or not when computing gradient later)
        icu = [Function(spaces.U) for i in range(d)]

        # Boundary control uin is extended in update() for each timestep
        uin = []

        return Controls(icu, uin)

    def cost_functionals(self, spaces, t, observations, controls):
        ds = self.ds
        dx = self.dx

        # Define static terms in cost functional (transient terms are added in update())
        J = 0

        jp = self.params.J

        # Add initial condition control terms
        if controls.icu:
            u0 = as_vector(controls.icu)
            J += Constant(jp.alpha_u0)              * u0**2       * dx
            J += Constant(jp.alpha_u0_div)          * div(u0)**2  * dx
            J += Constant(jp.alpha_u0_grad)         * grad(u0)**2 * dx
            J += Constant(jp.alpha_u0_grad_controlled)   * grad(u0)**2 * ds(self.controlled_boundary_ids)
            J += Constant(jp.alpha_u0_grad_uncontrolled) * grad(u0)**2 * ds(self.uncontrolled_boundary_ids)
            J += Constant(jp.alpha_u0_wall)              * u0**2       * ds(self.geometry.boundary_ids.wall)

        return CostFunctionals(J)

    def initial_conditions(self, spaces, controls):
        d = spaces.d

        # Get/create IC functions
        icu = controls.icu
        icp = Function(spaces.Q) # Is this actually used?

        # Copy initial condition from observations
        t0, z0 = self._observations[0]
        assert abs(t0) < self.params.dt*0.1, "Expecting t = 0!"
        for i in range(d):
            icu[i].assign(z0[i])

        return InitialConditions(icu, icp)

    def update(self, spaces, u, p, t, timestep, bcs, observations, controls, cost_functionals):
        d = len(u)

        # Get observations from this timestep
        tz, z = observations[timestep]
        z = as_vector(z)
        assert abs(t - tz) < self.problem.dt * 0.01, "Expecting matching time!"

        # ... Handle transient boundary control functions

        # Create new BC control functions at time t
        uin = [Function(spaces.U) for i in range(d)]

        # Make a record of BC controls in list
        controls.uin.append(uin)

        # Assign current control function values to the BC functions used in scheme forms
        bc = bcs.bcu.inflow
        for i in range(d):
            bc.functions[i].assign(uin[i])

        # ... Update cost functional with transient contributions
        ds = self.ds
        dx = self.dx
        jp = self.params.J
        J = cost_functionals.J

        # Add distance to observations at time t to cost functional
        J += (u-z)**2*dx

        # Add regularization of boundary control function to cost functional at time t
        J += Constant(jp.alpha_g)        * g**2       * ds(self.controlled_boundary_ids)
        J += Constant(jp.alpha_g_grad)   * grad(g)**2 * ds(self.controlled_boundary_ids)
        J += Constant(jp.alpha_g_volume)      * g**2       * dx
        J += Constant(jp.alpha_g_grad_volume) * grad(g)**2 * dx


def main():
    #parameters["adjoint"]["stop_annotating"] = True

    set_log_level(100)

    # Setup parameters shared between analytic and assimilation problem
    shared_problem_params = ParamDict(
            # Time
            dt=1e-3,
            T=0.01,#8,
            num_periods=None,
            )

    # Setup parameters shared between analytic and assimilation scheme
    shared_scheme_params = ParamDict(
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

    # Setup analytic problem to produce observations
    problem = AnalyticProblem(shared_problem_params)

    # Setup scheme
    scheme = CoupledPicard(shared_scheme_params)

    # Create unique casedir name
    params_string = str(hash(str(problem.params) + str(scheme.params)))[:8]
    date = "{t.tm_year}_{t.tm_mon}_{t.tm_mday}_{t.tm_hour}_{t.tm_min}".format(t=time.localtime())
    casedir = "results_demo_{}_{}_{}_{}".format(problem.shortname(), scheme.shortname(), date, params_string)


    # Setup postprocessor
    postprocessor = PostProcessor(ParamDict(
            casedir=casedir,
        ))
    plot_and_save = dict(plot=True, save=True)
    postprocessor.add_fields([
        Pressure(plot_and_save),
        Velocity(plot_and_save),
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
    observations = []
    for data in solver.isolve():
        z = [project(u, data.spaces.U) for u in data.u]
        observations.append((data.t, z))



    # Setup assimilation problem to reproduce observations
    daproblem = AssimilationProblem(shared_problem_params, problem.geometry, observations)

    # Setup scheme
    dascheme = CoupledPicard(shared_scheme_params)

    # Setup postprocessor
    dacasedir = casedir + "_da"
    dapostprocessor = PostProcessor(ParamDict(
            casedir=dacasedir,
        ))
    plot_and_save = dict(plot=True, save=True)
    dapostprocessor.add_fields([
        Pressure(plot_and_save),
        Velocity(plot_and_save),
        ])

    # Setup and run solver
    dasolver = NSSolver(daproblem, dascheme, dapostprocessor,
        ParamDict(
            timer_frequency=1,
            check_memory_frequency=1,
            enable_annotation=False,
            )
        )

    # Step through forward simulation
    diffs = []
    #daobservations = []
    for i, data in enumerate(dasolver.isolve()):
        # Compare observations with 'observations' from assimilation problem
        dz = [project(u, data.spaces.U) for u in data.u]
        t, z = observations[i]
        diff = assemble((as_vector(dz)-as_vector(z))**2*dx)
        diffs.append(diff)
    print "Diff in observations:", sqrt(sum(di for di in diffs))


    # TODO: Control annotation (off for problem, on for daproblem) and try replay


    # Try replaying with dolfin-adjoint
    success = da.replay_dolfin()
    print success


if __name__ == "__main__":
    main()
