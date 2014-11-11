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
# - g in CG1
# - Coarser time discretization of g(t) control with piecwise time UFL functions
# - BDF-2 time scheme
# - (Improve stability in scheme with N unannotated Picard iterations plus annotated Newton iterations)
# - Consider moving cost_functionals and observations out of problem and into external loop?
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

    def advance(self, t0, t, timestep, spaces, state, boundary_conditions, body_force, controls):
        pass # TODO: Move to problem base and document


class AnalyticProblem(ProblemBase):
    "2D pipe test problem with known stationary analytical solution."

    def __init__(self, params=None):
        ProblemBase.__init__(self, params)
        self.geometry = Geometry(self.params.geometry)
        self.initialize_geometry(self.geometry.mesh, self.geometry.facet_domains)
        self._init_womersley()

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

    def _init_womersley(self):
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

    def advance(self, t0, t, timestep, spaces, state, boundary_conditions, body_force, controls):
        "Advance boundary condition functions to time t."
        d = self.geometry.mesh.geometry().dim()
        facet_domains = self.geometry.facet_domains
        bc = boundary_conditions.bcu.inflow
        for u, w in zip(bc.functions, self.womersley):
            # Update Womersley Expression state
            w.set_t(t)
            # Trick to efficiently set u = w only at boundary dofs
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
        self._initial_velocity = initial_velocity

    @classmethod
    def default_params(cls):
        params = ProblemBase.default_params()
        return params

    def observations(self, spaces, t):
        return self._observations

    def controls(self, spaces):
        d = spaces.d

        # Initial condition control (may choose to use this or not when computing gradient later)
        initial_velocity = [Function(spaces.U, name="ic_u%d" % i) for i in range(d)]

        # Boundary control list is extended in advance() for each timestep
        g_timesteps = []

        return Controls(initial_velocity, g_timesteps)

    def initial_conditions(self, spaces, controls):
        d = spaces.d

        # Get/create IC functions
        icu = controls.initial_velocity
        icp = Function(spaces.Q, name="ic_p")

        # Copy given initial condition
        for i in range(d):
            icu[i].assign(self._initial_velocity[i])

        return InitialConditions(icu, icp)

    def advance(self, t0, t, timestep, spaces, state, boundary_conditions, body_force, controls):
        "Advance boundary condition functions to control values at time t."

        bc = boundary_conditions.bcu.inflow

        # Create new BC control functions at time t
        d = len(state[0])
        Ug = spaces.U # TODO: Configure control space
        g_at_t = [Function(Ug, name="g%d_%d" % (i,timestep)) for i in range(d)]


        # Update 'initial guess' of the boundary control functions at time t by copying from initial condition
        for g_control, initial_u in zip(g_at_t, self._initial_velocity):
            #dbc = DirichletBC(spaces.U, initial_u, self.geometry.facet_domains, bc.region)
            #dbc.apply(g_control.vector())
            g_control.assign(initial_u) # TODO: Implement alternative initial guesses (e.g. beta*u0 + eta)


        # Make a record of BC controls in list
        controls.g_timesteps.append((float(t), g_at_t))
        assert len(controls.g_timesteps) == timestep

        # Assign 'initial guess' control function values to the BC functions used in scheme forms
        for g_bc, g_control in zip(bc.functions, g_at_t):
            g_bc.assign(g_control, annotate=True) # Force annotation


class CostFunctional(Parameterized):
    def __init__(self, params, geometry, initial_velocity, observations):
        Parameterized.__init__(self, params)

        self._geometry = geometry
        self._observations = observations
        self._initial_velocity = initial_velocity

        self.J_timesteps = []

        # TODO: Configure with these id tuples
        # Control id tuples to make the below code more generic
        self.controlled_boundary_ids = (geometry.boundary_ids.left,)
        self.uncontrolled_boundary_ids = (geometry.boundary_ids.right,)
        self.wall_boundary_ids = (geometry.boundary_ids.wall,)

        # Setup integration measures
        self.dx = Measure("dx", domain=geometry.mesh)
        self.ds = Measure("ds", domain=geometry.mesh, subdomain_data=geometry.facet_domains)
        self.dsw = ds(self.wall_boundary_ids)
        self.dsc = ds(self.controlled_boundary_ids)
        self.dsu = ds(self.uncontrolled_boundary_ids)

        # Wrap regularization parameters in Constants
        for k in self.params.keys():
            alpha = self.params[k]
            if abs(alpha) < 1e-13:
                # This is to make the cost functional terms disappear
                # (unfortunately increasing recompilation)
                # TODO: This didn't work out exactly as planned
                alpha = 0
            else:
                # This is to avoid recompiling the cost functional for variations in alpha values
                alpha = Constant(alpha, name=k)
            setattr(self, k, alpha)

    @classmethod
    def default_params(cls):
        params = ParamDict(
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
        return params

    def update(self, timestep, t, spaces, state, controls):
        "Update cost functionals."
        u, p = state
        d = len(u)

        # Setup integration measures
        dx = self.dx
        dsw = self.dsw
        dsc = self.dsc
        dsu = self.dsu

        # Choose time measure to apply to static cost functional terms
        dts = NoOp()
        #dts = dt[START_TIME]
        #dts = dt[FINISH_TIME]

        # Choose time measure to apply to transient terms in cost functional
        dtt = NoOp()
        #dtt = dt[timestep]
        #dtt = dt[float(t)]
        #dtt = dt[FINISH_TIME]


        # ... Update cost functional with distance between state and observations at time t
        # Get observations from this timestep
        tz, z = self._observations[timestep]
        assert abs(float(t) - tz) < 1e-6, "Expecting matching time!"

        # Project the velocity state into the observation function space
        Uz = z[0].function_space()
        u_t = [project(u[i], Uz, name="u%d_at_ts%d" % (i, timestep)) for i in range(d)]

        J_terms = [
            (as_vector(u_t) - as_vector(z))**2 * dx * dtt
            ]

        # ... Update cost functional with regularization of initial velocity control
        if timestep == 0 and controls.initial_velocity:
            u0 = as_vector(controls.initial_velocity)
            J_terms += [
                self.alpha_u0              * u0**2       * dx * dts,
                self.alpha_u0_div          * div(u0)**2  * dx * dts,
                self.alpha_u0_grad         * grad(u0)**2 * dx * dts,
                self.alpha_u0_grad_controlled   * grad(u0)**2 * dsc * dts,
                self.alpha_u0_grad_uncontrolled * grad(u0)**2 * dsu * dts,
                self.alpha_u0_wall              * u0**2       * dsw * dts,
                ]

        # ... Update cost functional with boundary control regularization
        if timestep >= 1:
            # Get last boundary control function (there's none at timestep 0)
            assert len(controls.g_timesteps) == timestep, "Expecting control functions at every timestep >= 1."
            tg, g = controls.g_timesteps[timestep-1]
            assert abs(float(t) - tg) < 1e-6, "Expecting matching time!"
            g = as_vector(g)

            # Add regularization of boundary control function to cost functional at time t
            J_terms += [
                self.alpha_g             * g**2       * dsc * dtt,
                self.alpha_g_grad        * grad(g)**2 * dsc * dtt,
                self.alpha_g_volume      * g**2       * dx * dtt,
                self.alpha_g_grad_volume * grad(g)**2 * dx * dtt,
                ]

        # Accumulate all contributions and add to list
        self.J_timesteps += [rebalance_form(sum(J_terms))]


    def accumulate(self):
        "Collect terms gathered up during timestep updates."
        # This gives us a ridiculously long form...
        # Using a hierarchic binary tree sum structure to
        # avoid python recursion limit problems in ufl.
        return rebalance_form(sum(self.J_timesteps))


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

    def advance(self, t0, t, timestep, spaces, state, boundary_conditions, body_force, controls):
        "Advance boundary condition functions to control values at time t."
        bc = boundary_conditions.bcu.inflow
        tg, g_at_t = self._controls.g_timesteps[timestep-1]
        assert abs(float(t) - tg) < 1e-6, "Expecting matching time!"
        for g_bc, g_control in zip(bc.functions, g_at_t):
            g_bc.assign(g_control)


def main():
    set_log_level(100)


    # Configure geometry
    geometry_params = ParamDict(
            refinement_level=2,
            length=2.0,
            radius=0.5,
            )

    # Configure time
    time_params = ParamDict(
            dt=1e-2,
            T=0.04,#8,
            num_periods=None,
            )

    # Configure cost functional
    J_params = ParamDict(
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

    # Configure problems for the three stages
    analytic_problem_params = ParamDict(time_params, geometry=geometry_params)
    forward_problem_params = ParamDict(time_params)
    final_problem_params = ParamDict(time_params)


    # Configure scheme (reused in all stages)
    Scheme = CoupledScheme
    scheme_params = ParamDict(
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


    # Create unique casedir name
    all_params = [scheme_params, analytic_problem_params, forward_problem_params, final_problem_params]
    params_string = str(hash(' '.join(map(str, all_params))))[:8]
    date = timestamp()
    casedir = "results_{}_{}_{}".format(Scheme.__name__, date, params_string)


    # TODO: Make this a function
    run_analytic_problem = True
    if run_analytic_problem:
        # DON'T annotate the initial observation production
        parameters["adjoint"]["stop_annotating"] = True

        # Setup analytic problem to produce observations
        analytic_problem = AnalyticProblem(analytic_problem_params)

        # Setup postprocessor
        fp = dict(save=True, plot=False)
        analytic_postprocessor = PostProcessor(dict(casedir=casedir+"_analytic"))
        analytic_postprocessor.add_fields([
            Pressure(fp),
            Velocity(fp),
            ])

        # Setup analytic_solver
        analytic_solver_params = ParamDict(
            enable_annotation=False,
            #timer_frequency=1,
            #check_memory_frequency=1,
            )
        analytic_solver = NSSolver(analytic_problem, CoupledScheme(scheme_params), analytic_postprocessor, analytic_solver_params)

        # Step through analytic simulation
        observations = []
        for data in analytic_solver.isolve():
            Uz = data.spaces.U # TODO: Configure different observation space
            z = [project(u, Uz, name="z_at_ts%d" % data.timestep) for u in data.u]
            observations.append((data.t, z))

            # Debugging evaluations:
            if 0:
                # Should have that z(t) = u(t) = g(t) on boundary bc.region.
                # Compute and print. If this fails, the weak boundary conditions have failed.
                bc = data.bcs.bcu.inflow
                dsc = analytic_problem.ds(bc.region)

                u = as_vector(data.u) # u(t)
                g = as_vector(bc.functions) # g(t)
                z = as_vector(z) # g(t)

                u2 = assemble(u**2*dsc)
                g2 = assemble(g**2*dsc)
                z2 = assemble(z**2*dsc)

                diffg = assemble((u-g)**2*dsc)
                diffz = assemble((u-z)**2*dsc)

                print "ts ={}, t={}, u(t)^2={}, g(t)^2={}, z(t)^2={}, (u(t)-g(t))^2={}, (u(t)-z(t))^2={}".format(
                    data.timestep, float(data.t), u2, g2, z2, diffg, diffz)
                #assert diffg < 1e-10
                #assert diffz < 1e-10

        # Consistency check of observation times
        dt = time_params.dt
        T = time_params.T
        assert abs(observations[0][0]) < dt * 0.01, "Expecting t = 0 at first observation!"
        assert T - dt*0.01 < abs(observations[-1][0]) < T + dt * 0.51, "Expecting t = T at last observation!"

        # Extract initial velocity from observations,
        # for another problem this may be given separately,
        # e.g. if the observations are not everywhere, noisy,
        # or in a different space, then the initial velocity
        # could be e.g. the solution to an artificially constructed
        # forward problem with a reasonable flow rate
        initial_velocity = observations[0][1]

        # Extract geometry from analytic problem,
        # for another problem this may be
        geometry = analytic_problem.geometry


    # TODO: Make this a function
    run_forward_problem = True
    if run_forward_problem:
        # DO annotate the forward model
        parameters["adjoint"]["stop_annotating"] = False

        # Enable derivative testing (TODO: Does this actually do anything?)
        parameters["adjoint"]["test_derivative"] = True

        # Setup assimilation problem to reproduce observations
        forward_problem = AssimilationProblem(forward_problem_params, geometry, initial_velocity, observations)

        # Setup postprocessor
        fp = dict(save=True, plot=False)
        forward_postprocessor = PostProcessor(dict(casedir=casedir+"_forward"))
        forward_postprocessor.add_fields([
            Pressure(fp),
            Velocity(fp),
            ])

        # Setup and run solver
        forward_solver_params = ParamDict(
            enable_annotation=True,
            #timer_frequency=1,
            #check_memory_frequency=1,
            )
        forward_solver = NSSolver(forward_problem, CoupledScheme(scheme_params), forward_postprocessor, forward_solver_params)

        # Setup cost functional
        CF = CostFunctional(J_params, geometry, initial_velocity, observations)
        dx = CF.dx
        ds = CF.ds
        dsc = CF.dsc

        # Step through forward simulation and build cost functional
        diffs = []
        for i, data in enumerate(forward_solver.isolve()):
            CF.update(data.timestep, data.t, data.spaces, data.state, data.controls)

            # Compare observations with 'observations' from assimilation problem
            if 1:
                # Avoid annotating these steps
                parameters["adjoint"]["stop_annotating"] = True
                dx = CF.dx
                dz = [project(u, data.spaces.U, name="dz%d"%j) for j,u in enumerate(data.u)]
                t, z = observations[i]
                diff = assemble((as_vector(dz)-as_vector(z))**2*dx) / assemble((as_vector(z))**2*dx)
                diffs.append(diff)
                parameters["adjoint"]["stop_annotating"] = False

        # Stop annotating after forward simulation
        parameters["adjoint"]["stop_annotating"] = True

        # HACK drop last control
        # FIXME XXX TODO: Need this still? TODO: Add timestamp to g_timesteps to assert for timeshifting bugs.
        #data.controls.g_timesteps.pop()

        # Extract flat list of boundary control components
        all_g_components = [gc for (tg, g) in data.controls.g_timesteps for gc in g]

        # Accumulate terms from cost functional.
        J = CF.accumulate()
        #print str(J)
        J = da.Functional(J)

        # TODO: Clean up and formalize norm computations here and elsewhere, maybe
        #       add to postprocessor so we get them stored in the result directories?
        # Compute norms of boundary control functions
        #dsc = CF.dsc
        #print "|g| =", [sqrt(assemble(as_vector(g)**2*dsc)) for (tg, g) in data.controls.g_timesteps]
        #print "|u-obs| on inflod =", [sqrt(assemble((u-z)**2*dsc)) for (tg, g) in data.controls.g_timesteps]
        #print "Diff in observations:", sqrt(sum(di for di in diffs)/len(diffs))


    # Always write tape to files for eventual debugging
    adj_html("forward.html", "forward")
    adj_html("adjoint.html", "adjoint")


    # Try replaying with dolfin-adjoint (works within a small tolerance around 1e-12 to 1e-14)
    run_replay = False
    if run_replay:
        success = da.replay_dolfin()
        print "replay success:", success
        return


    # TODO: Make this a function
    # Plot gradients
    run_plotting_gradients = False
    if run_plotting_gradients:
        # Select controls
        m_values = list(all_g_components)
        #m_values = [all_g_components[0]] # DEBUGGING
        #m_values = [all_g_components[-1]] # DEBUGGING
        m = [da.Control(mv) for mv in m_values]

        # Compute dJ/dm at m and plot animated gradients on boundary:
        dJdm = da.compute_gradient(J, m, forget=False)
        Vdj = dJdm[0].function_space()
        djp = Function(Vdj)
        bc = data.bcs.bcu.inflow
        for i, dj in enumerate(dJdm):
            # Using dbc trick to set only boundary values of djp, ignoring the rest of the domain
            dbc = DirichletBC(Vdj, dj, geometry.facet_domains, bc.region)
            djp.vector().zero()
            dbc.apply(djp.vector())
            plot(djp, title="dj %d" % i)
            time.sleep(0.1)
        interactive()
        return


    # TODO: Make this a function
    # Test gradient computation with convergence order tests
    run_taylor_test = True
    if run_taylor_test:
        # This is to remove rounding errors caused by tensor representation when computing norms
        parameters["form_compiler"]["representation"] = "quadrature"

        dsc = CF.dsc
        if 0:
            print "Compute z0-g"
            tz, z = observations[0]
            for i, (tg, g) in enumerate(data.controls.g_timesteps):
                assert abs(tz - tg) < 1e-6, "Expecting matching time!"
                for j in range(len(g)):
                    print j, assemble((g[j]-z[j])**2*dsc)

            # Copy observations into boundary controls
            for i, (tg, g) in enumerate(data.controls.g_timesteps):
                tz, z = observations[i+1]
                assert abs(tz - tg) < 1e-6, "Expecting matching time!"
                for j in range(len(g)):
                    g[j].assign(z[j], annotate=False)

        # Update boundary control functions at time t
        for i, (tg, g) in enumerate(data.controls.g_timesteps):
            for u, w in zip(g, analytic_problem.womersley):
                w.set_t(tg)
                u.vector().zero()
                for bid in CF.controlled_boundary_ids:
                    dbc = DirichletBC(data.spaces.U, w,
                                      geometry.facet_domains, bid)
                    dbc.apply(u.vector())

        print "Compute zi-g"
        for i, (tg, g) in enumerate(data.controls.g_timesteps):
            tz, z = observations[i+1]
            assert abs(tz - tg) < 1e-6, "Expecting matching time!"
            for j in range(len(g)):
                print j, assemble((g[j]-z[j])**2*dsc)

        #timesteps 2,3,4:  (u-z)**2*dx != 0

        # Select controls for componentwise taylor test
        m_values = list(all_g_components)
        #m_values = [all_g_components[0]] # DEBUGGING
        #m_values = [all_g_components[-1]] # DEBUGGING

        for i, mv in enumerate(m_values):
            m = da.Control(mv)
            RJ = da.ReducedFunctional(J, m)
            Jm = RJ(mv)
            dJdm = da.compute_gradient(J, m, forget=False)
            conv_rate = da.taylor_test(RJ, m, Jm, dJdm)
            print
            print "m = g_%d(t_%d)" % ((i%2), (i//2))
            print "  J(m)      =", Jm
            print "  dJdm(m)   =", assemble(dJdm**2*dx)
            print "  conv_rate =", conv_rate
        return


    # TODO: Make this a function
    # Optimize! This is where the main work goes!
    run_minimize = True
    if run_minimize:
        # Clear problem functions so we're sure dolfin-adjoint gets its values from the tape!
        for gc in all_g_components:
            gc.vector().zero()
        bc = data.bcs.bcu.inflow
        for fc in bc.functions:
            fc.vector().zero()

        # Select controls TODO: Test adding initial_velocity here
        m_values = list(all_g_components)
        #m_values = [all_g_components[0]] # DEBUGGING
        #m_values = [all_g_components[-1]] # DEBUGGING
        m = [da.Control(mv) for mv in m_values]

        # TODO: Expose parameters
        RJ = da.ReducedFunctional(J, m)
        m_opt = minimize(RJ, tol=1e-8, method="L-BFGS-B", options={"disp": True, "maxiter": 100})
        #m_opt = minimize(RJ)

        #print m_opt
        #for i, mo in enumerate(m_opt):
        #    if i % 2 == 0:
        #        plot(mo, title="m%d" % (i//2))
        #interactive()
        #print "|m_opt| =", [sqrt(assemble(as_vector(g)**2*dsc)) for g in m_opt]


    # TODO: Make this a function
    # Run final problem with controls input from optimization run
    run_final_problem = True
    if run_final_problem:
        # Setup controls on the right format
        u0 = data.controls.initial_velocity
        d = len(u0)
        glist = [[m_opt[d*i+j] for j in range(d)] for i in range(len(m_opt)//d)] # NB! Assumes m_values above use all g components!
        controls = Controls(u0, glist)

        # Setup problem
        final_problem = FinalReplayProblem(final_problem_params, geometry, controls)

        # Setup postprocessor
        fp = dict(save=True, plot=False)
        final_postprocessor = PostProcessor(dict(casedir=casedir + "_final"))
        final_postprocessor.add_fields([
            Pressure(fp),
            Velocity(fp),
            ])

        # Setup and run solver
        final_solver_params = ParamDict(
            enable_annotation=False,
            #timer_frequency=1,
            #check_memory_frequency=1,
            )
        final_solver = NSSolver(final_problem, CoupledScheme(scheme_params), final_postprocessor, final_solver_params)

        # Step through forward simulation
        diffs = []
        for i, data in enumerate(final_solver.isolve()):
            pass


if __name__ == "__main__":
    main()
