#!/usr/bin/env python

import dolfin
import dolfin_adjoint as da

from cbcpost import ParamDict, PostProcessor, Parameterized

from cbcflow import *
from cbcflow.dol import *

from collections import namedtuple
import numpy as np
import time


# TODO: Improving scheme:
# - BDF-2 time scheme
# - Coarser time discretization of g(t) control with piecwise time UFL functions
# TODO: Testing:
# - Setting u0 initial control values to (beta*z(0) + eta(0, x))
# - Setting g initial control values to (beta*z(t) + eta(t, x))


# TODO: Move to cbcpost utils
def timestamp():
    return "{t.tm_year}_{t.tm_mon}_{t.tm_mday}_{t.tm_hour}_{t.tm_min}_{t.tm_sec}".format(t=time.localtime())

def bsum(seq, zero=0):
    "Add items in sequence in a binary tree structure."
    n = len(seq)
    if n <= 2:
        if n == 0:
            return zero
        elif n == 1:
            return seq[0]
        elif n == 2:
            return seq[0] + seq[1]
    m = n // 2
    return bsum(seq[:m]) + bsum(seq[m:])

import ufl
def rebalance_expr(expr): # FIXME: Validate this
    output_terms = []
    input_terms = [expr]
    while input_terms:
        term = input_terms.pop()
        if isinstance(term, ufl.classes.Sum):
            input_terms.extend(term.operands())
        else:
            output_terms.append(term)
    return bsum(output_terms)

def rebalance_form(form): # FIXME: Validate this
    integrals = [itg.reconstruct(integrand=rebalance_expr(itg.integrand())) for itg in form.integrals()]
    return ufl.Form(integrals)

class NoOp(object):
    def __mul__(self, other):
        return other
    def __rmul__(self, other):
        return other


# TODO: Add types like these to cbcflow to document BC interface?
InitialConditions = namedtuple("InitialConditions", ["icu", "icp"])
VelocityBC = namedtuple("VelocityBC", ["functions", "region", "method"])
PressureBC = namedtuple("PressureBC", ["function", "region"])
BoundaryConditions = namedtuple("BoundaryConditions", ["bcu", "bcp"])


# Create some problem specific types to make tuple usage safer and more readable
Controls = namedtuple("Controls", ["initial_velocity", "g_timesteps"])
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
        mesh.coordinates()[:, 0] = l * mesh.coordinates()[:, 0]
        mesh.coordinates()[:, 1] = 2.0 * r * (mesh.coordinates()[:, 1] - 0.5)

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
        uin = [Function(spaces.U, name="bc_g%d" % i) for i in range(d)]
        inflow = VelocityBC(uin, bids.left, "nietche")

        # Create no-slip bcs
        u_noslip = [Constant(0.0, name="bc_noslip%d" % i) for i in range(d)]
        noslip = VelocityBC(u_noslip, bids.wall, "strong")

        # Create outflow bcs for pressure
        outflow = PressureBC(Constant(0.0, name="bc_p_out"), bids.right)

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
            time_values = np.linspace(0.0, P, 4*wp.num_womersley_coefficients)
            t = np.mod((P - time_values) / P, 1.0)
            time_profile = wp.Qfloor + (wp.Qpeak - wp.Qfloor) * np.sin(pi*t**3)**2
            Q_coeffs = zip(time_values, wp.Q * time_profile)
        else:
            # Setup coefficients for a constant flow rate
            Q_coeffs = [(0.0, wp.Q), (1.0, wp.Q)]

        # Create womersley objects
        self.womersley = make_womersley_bcs(Q_coeffs, mesh, bids.left, nu, None, facet_domains,
                                            "Q", num_fourier_coefficients=wp.num_womersley_coefficients)

    def initial_conditions(self, spaces, controls):
        d = spaces.d

        # Create IC functions
        icu = [Function(spaces.U, name="ic_u%d" % i) for i in range(d)]
        icp = Function(spaces.Q, name="ic_p")

        # Interpolate womersley solution into entire initial condition (only possible for this simple pipe case!)
        for i in range(d):
            self.womersley[i].set_t(0.0)
            icu[i].interpolate(self.womersley[i])

        return InitialConditions(icu, icp)

    def update(self, spaces,
               # Input from time t:
               u, p, t, timestep,
               # Output to be updated to time t:
               bcs, observations, controls, cost_functionals):
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
    """Problem setup for dolfin-adjoint annotation.

    Observations and initial condition for velocity are given.

    Controls are both boundary inflow and velocity initial condition.

    The cost functional is parameterized with a bunch of
    regularization parameters, by default most of these are zero.
    """

    def __init__(self, params, geometry, initial_velocity, observations):
        ProblemBase.__init__(self, params)
        self.geometry = geometry
        self.initialize_geometry(self.geometry.mesh, self.geometry.facet_domains)

        # Store observations for later
        self._observations = observations
        self.initial_velocity = initial_velocity

        # Control id tuples to make the below code more generic (deliberately not including wall here)
        self.controlled_boundary_ids = (self.geometry.boundary_ids.left,)
        self.uncontrolled_boundary_ids = (self.geometry.boundary_ids.right,)

        # Wrap regularization parameters in Constants
        for k in self.params.J.keys():
            alpha = self.params.J[k]
            if abs(alpha) < 1e-13:
                # This is to make the cost functional terms disappear
                # (unfortunately increasing recompilation)
                alpha = 0
            else:
                # This is to avoid recompiling the cost functional for variations in alpha values
                alpha = Constant(alpha, name=k)
            setattr(self, k, alpha)

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
            alpha_g=0.0,
            alpha_g_grad=0.0,
            alpha_g_volume=0.0,
            alpha_g_grad_volume=0.0,
            )
        params.update(J=J)
        return params

    def observations(self, spaces, t):
        return self._observations

    def controls(self, spaces):
        d = spaces.d

        # Initial condition control (may choose to use this or not when computing gradient later)
        initial_velocity = [Function(spaces.U, name="ic_u%d" % i) for i in range(d)]

        # Boundary control list is extended in update() for each timestep
        g_timesteps = []

        return Controls(initial_velocity, g_timesteps)

    def cost_functionals(self, spaces, t, observations, controls):
        ds = self.ds
        dx = self.dx

        # Define static terms in cost functional (transient terms are added in update())
        cost_functionals = []

        # Choose time measure to apply to static cost functional terms
        dtt = NoOp()
        #dtt = dt[START_TIME]
        #dtt = dt[FINISH_TIME]

        # Add initial condition control terms
        if controls.initial_velocity:
            u0 = as_vector(controls.initial_velocity)
            J_terms = [
                self.alpha_u0              * u0**2       * dx * dtt,
                self.alpha_u0_div          * div(u0)**2  * dx * dtt,
                self.alpha_u0_grad         * grad(u0)**2 * dx * dtt,
                self.alpha_u0_grad_controlled   * grad(u0)**2 * ds(self.controlled_boundary_ids) * dtt,
                self.alpha_u0_grad_uncontrolled * grad(u0)**2 * ds(self.uncontrolled_boundary_ids) * dtt,
                self.alpha_u0_wall              * u0**2       * ds(self.geometry.boundary_ids.wall) * dtt,
                ]
            J = sum(J_terms)
            cost_functionals.append(J)

        return cost_functionals

    def initial_conditions(self, spaces, controls):
        d = spaces.d

        # Get/create IC functions
        icu = controls.initial_velocity
        icp = Function(spaces.Q, name="ic_p")

        # Copy initial condition from observations
        for i in range(d):
            icu[i].assign(self.initial_velocity[i])

        return InitialConditions(icu, icp)

    def update(self, spaces,
               # Input from time t:
               u, p, t, timestep,
               # Output to be updated to time t:
               bcs, observations, controls, cost_functionals):
        d = len(u)
        facet_domains = self.geometry.facet_domains

        # Get observations from this timestep
        tz, z = observations[timestep]
        z = as_vector(z)
        assert abs(float(t) - tz) < self.params.dt * 0.01, "Expecting matching time!"

        # ... Handle transient boundary control functions

        bc = bcs.bcu.inflow

        # Create new BC control functions at time t
        g = [Function(spaces.U, name="g%d_%d" % (i,timestep)) for i in range(d)]

        # Make a record of BC controls in list
        controls.g_timesteps.append(g)

        # Update 'initial guess' of the boundary control functions at time t by copying from initial condition
        for i in range(d):
            #dbc = DirichletBC(spaces.U, self.initial_velocity[i], facet_domains, bc.region)
            #dbc.apply(g[i].vector())
            g[i].assign(self.initial_velocity[i])

        # Assign 'initial guess' control function values to the BC functions used in scheme forms
        for i in range(d):
            bc.functions[i].assign(g[i], annotate=True)

        # ... Update cost functional with transient contributions
        ds = self.ds
        dx = self.dx
        g = as_vector(g)

        # Choose time measure to apply to per-timestep terms in cost functional
        dtt = NoOp()
        #dtt = dt[timestep]
        #dtt = dt[float(t)]
        #dtt = dt[FINISH_TIME]

        # Project the velocity state into the observation function space
        u_t = project(u, spaces.V, name="u_at_ts%d" % timestep)
        J_terms = [
            # Distance between state and observations at time t
            (u_t - z)**2 * dx * dtt,
            # Add regularization of boundary control function to cost functional at time t
            self.alpha_g             * g**2       * ds(self.controlled_boundary_ids) * dtt,
            self.alpha_g_grad        * grad(g)**2 * ds(self.controlled_boundary_ids) * dtt,
            self.alpha_g_volume      * g**2       * dx * dtt,
            self.alpha_g_grad_volume * grad(g)**2 * dx * dtt,
            ]
        J = sum(J_terms)
        cost_functionals.append(J)


class FinalReplayProblem(ProblemBase):
    "Problem setup to run with the final control functions."

    def __init__(self, params, geometry, controls):
        ProblemBase.__init__(self, params)
        self.geometry = geometry
        self.initialize_geometry(self.geometry.mesh, self.geometry.facet_domains)
        self._controls = controls

    def initial_conditions(self, spaces, controls):
        icu = self._controls.initial_velocity
        icp = Function(spaces.Q, name="ic_p")
        return InitialConditions(icu, icp)

    def update(self, spaces,
               # Input from time t:
               u, p, t, timestep,
               # Output to be updated to time t:
               bcs, observations, controls, cost_functionals):
        bc = bcs.bcu.inflow
        us = bc.functions
        gs = self._controls.g_timesteps[timestep-1]
        for u, g in zip(us, gs):
            u.assign(g)


def main():
    # Don't annotate the initial observation production
    parameters["adjoint"]["stop_annotating"] = True

    set_log_level(100)

    # Configure geometry
    geometry_params = ParamDict(
            refinement_level=2,
            length=2.0,
            radius=0.5,
            )

    # Setup parameters shared between analytic and assimilation problem
    shared_problem_params = ParamDict(
            # Time
            dt=1e-2,
            T=0.4,#8,
            num_periods=None,
            )

    # Setup parameters for cost functional
    J_params = ParamDict(
        # Initial control regularization
        alpha_u0=0.0,
        alpha_u0_div=0.0,
        alpha_u0_grad=0.0,
        alpha_u0_grad_controlled=0.0,
        alpha_u0_grad_uncontrolled=0.0,
        alpha_u0_wall=0.0,
        # Boundary control regularization
        alpha_g=1e-8,
        alpha_g_grad=1e-8,
        alpha_g_volume=1e-12,
        alpha_g_grad_volume=0.0,
        )

    # Setup problem params
    problem_params = ParamDict(shared_problem_params,
                               geometry=geometry_params)

    # Setup daproblem params
    daproblem_params = ParamDict(shared_problem_params,
                                 J=J_params)

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
            picard_newton_fraction=1.0, # 0.0 = Picard, 1.0 = Newton, in (0.0, 1.0) = mix
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
    problem = AnalyticProblem(problem_params)

    # Setup scheme
    scheme = CoupledScheme(shared_scheme_params)

    # Create unique casedir name
    params_string = str(hash(str(problem.params) + str(scheme.params)))[:8]
    date = timestamp()
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
    solver_params = ParamDict(
        enable_annotation=False,
        #timer_frequency=1,
        #check_memory_frequency=1,
        )
    solver = NSSolver(problem, scheme, postprocessor, solver_params)

    # Step through forward simulation
    observations = []
    for data in solver.isolve():
        z = [project(u, data.spaces.U, name="z_at_ts%d" % data.timestep) for u in data.u]

        # Observations look good:
        #plot(as_vector(z), title="z")
        #plot(data.u, title="u")
        #interactive()

        observations.append((data.t, z))

    # Extract initial velocity from observations,
    # for another problem this may be given separately,
    # e.g. if the observations are not everywhere, noisy,
    # or in a different space, then the initial velocity
    # could be e.g. the solution to an artificially constructed
    # forward problem with a reasonable flow rate
    assert abs(observations[0][0]) < problem.params.dt * 0.01, "Expecting t = 0!"
    initial_velocity = observations[0][1]


    # DO annotate the forward model
    parameters["adjoint"]["stop_annotating"] = False
    adj_reset()
    parameters["adjoint"]["test_derivative"] = True

    # Setup assimilation problem to reproduce observations
    daproblem = AssimilationProblem(daproblem_params, problem.geometry, initial_velocity, observations)

    # Setup scheme
    dascheme = CoupledScheme(shared_scheme_params)

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
    dasolver_params = ParamDict(
        enable_annotation=True,
        #timer_frequency=1,
        #check_memory_frequency=1,
        )
    dasolver = NSSolver(daproblem, dascheme, dapostprocessor, dasolver_params)

    # Step through forward simulation
    diffs = []
    for i, data in enumerate(dasolver.isolve()):

        # Avoid annotating these steps
        parameters["adjoint"]["stop_annotating"] = True

        # Compare observations with 'observations' from assimilation problem
        dz = [project(u, data.spaces.U) for u in data.u]
        t, z = observations[i]
        diff = assemble((as_vector(dz)-as_vector(z))**2*dx) / assemble((as_vector(z))**2*dx)
        diffs.append(diff)

        parameters["adjoint"]["stop_annotating"] = False

    # Stop annotating after forward simulation
    parameters["adjoint"]["stop_annotating"] = True

    # Always write tape to files for eventual debugging
    adj_html("forward.html", "forward")
    adj_html("adjoint.html", "adjoint")

    # Try replaying with dolfin-adjoint (works within a small tolerance around 1e-12 to 1e-14)
    try_replay = False
    if try_replay:
        success = da.replay_dolfin()
        print "replay success:", success


    # Accumulate terms from cost functional.
    # This gives us a ridiculously long form...
    # Trying to use a hierarchic binary tree sum structure to avoid
    # recursion limit problems in ufl, but this doesn't really work...
    J = rebalance_form(sum(data.cost_functionals))

    # Compute norms of boundary control functions
    dsi = problem.ds(problem.geometry.boundary_ids.left)
    print "|m| =", [sqrt(assemble(as_vector(g)**2*dsi)) for g in data.controls.g_timesteps]
    #print "|u-obs| on inflod =", [sqrt(assemble((u-z)**2*dsi)) for g in data.controls.g_timesteps]
    #print "Diff in observations:", sqrt(sum(di for di in diffs)/len(diffs))

    # Setup controls and functionals for dolfin-adjoint
    gcs = [gc for g in data.controls.g_timesteps for gc in g]
    m = [da.Control(gc) for gc in gcs]
    J = da.Functional(J)
    RJ = da.ReducedFunctional(J, m)


    # Plot gradients
    try_plotting_gradients = False
    if try_plotting_gradients:
        # Compute dJ/dm at m and plot animated gradients on boundary:
        dJdm = da.compute_gradient(J, m, forget=False)
        Vdj = dJdm[0].function_space()
        djp = Function(Vdj)
        for i, dj in enumerate(dJdm):
            # Using dbc trick to set only boundary values of djp, ignoring the rest of the domain
            dbc = DirichletBC(Vdj, dj, problem.geometry.facet_domains, data.bcs.bcu.inflow.region)
            djp.vector().zero()
            dbc.apply(djp.vector())
            plot(djp, title="dj %d" % i)
            time.sleep(0.1)
        interactive()
        return


    # Not capable of running taylor test...
    try_taylor_test = False
    if try_taylor_test:
        # Compute J(m)
        Jm = RJ(gcs)
        print "J(m) =", Jm

        # Compute dJ/dm at m
        dJdm = da.compute_gradient(J, m, forget=False)

        # TODO: Run taylor test!
        conv_rate = da.taylor_test(RJ, m, Jm, dJdm) # This fails, bug reported in dolfin-adjoint
        print "conv_rate =", conv_rate
        return


    run_minimize = True
    if run_minimize:
        # Try optimizing
        for gc in gcs:
            gc.vector().zero()
        bc = data.bcs.bcu.inflow
        for fc in bc.functions:
            fc.vector().zero()
        m_opt = minimize(RJ, tol=1e-8, method="L-BFGS-B", options={"disp": True, "maxiter": 100})
        #m_opt = minimize(RJ)

        #print m_opt
        #for i, mo in enumerate(m_opt):
        #    if i % 2 == 0:
        #        plot(mo, title="m%d" % (i//2))
        #interactive()
        #print "|m_opt| =", [sqrt(assemble(as_vector(g)**2*dsi)) for g in m_opt]


    run_final = True
    if run_final:
        # Run final problem
        u0 = data.controls.initial_velocity
        d = len(u0)
        glist = [[m_opt[d*i+j] for j in range(d)] for i in range(len(m_opt)//d)]
        controls = Controls(u0, glist)
        finproblem = FinalReplayProblem(shared_problem_params, problem.geometry, controls)

        # Setup scheme
        finscheme = CoupledScheme(shared_scheme_params)

        # Setup postprocessor
        fincasedir = casedir + "_fin"
        finpostprocessor = PostProcessor(ParamDict(
                casedir=fincasedir,
            ))
        plot_and_save = dict(plot=True, save=True)
        finpostprocessor.add_fields([
            Pressure(plot_and_save),
            Velocity(plot_and_save),
            ])

        # Setup and run solver
        finsolver_params = ParamDict(
            enable_annotation=True,
            #timer_frequency=1,
            #check_memory_frequency=1,
            )
        finsolver = NSSolver(finproblem, finscheme, finpostprocessor, finsolver_params)

        # Step through forward simulation
        diffs = []
        for i, data in enumerate(finsolver.isolve()):
            pass


if __name__ == "__main__":
    main()
