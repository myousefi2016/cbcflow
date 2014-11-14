#!/usr/bin/env python

import dolfin
import dolfin_adjoint as da

from cbcpost import ParamDict, PostProcessor, Parameterized

from cbcflow import *
from cbcflow.dol import *

from collections import namedtuple
import numpy as np
import time


# TODO: coupled scheme:
# - Add time discretization schemes: BDF-2 and theta-rule
# - Improve stability in scheme with N unannotated Picard iterations plus annotated Newton iterations?
# TODO: optimization setup:
# - Allow z in CG1/DG0
# - Allow g in CG1 (will be projected into spaces.V in scheme?)
# - Coarser time discretization of g(t) control with piecwise time UFL functions
# TODO: Testing:
# - Setting u0 initial control values to (beta*z(0) + eta(0, x))
# - Setting g initial control values to (beta*z(t) + eta(t, x))
# TODO: Cbcflow structure:
# - Want to get u in V and p in Q as the physically scaled pressure in advance() and from yield?
#   Get rid of converters to achieve this? Only problem is performance of subfunction assignment.

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

        # Geometry parameters
        geometry = Geometry.default_params()

        # Analytical solution parameters
        womersley = ParamDict(
            Q=1.0,
            Qfloor=0.25,
            Qpeak=1.0,
            transient_Q=True,
            num_womersley_coefficients=25,
            )

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

        # Create initial condition functions for velocity
        icu = [Function(spaces.U, name="ic_u%d" % i) for i in range(d)]

        # Interpolate womersley solution into entire initial condition (only possible for this simple pipe case!)
        for i in range(d):
            self.womersley[i].set_t(0.0)
            icu[i].interpolate(self.womersley[i])

        return InitialConditions(icu, 0.0)

    def advance(self, t0, t, timestep, spaces, state, boundary_conditions, body_force, controls):
        "Advance boundary condition functions to time t."
        facet_domains = self.geometry.facet_domains
        bc = boundary_conditions.bcu.inflow
        d = len(bc.functions)
        for i in range(d):
            # Update internal Womersley Expression state
            self.womersley[i].set_t(t)
            # Trick to efficiently set g = w only at boundary dofs,
            # avoiding full interpolation of the expensive Womersley expression.
            # This is not recognized by dolfin-adjoint.
            bc.functions[i].vector().zero()
            dbc = DirichletBC(spaces.U, self.womersley[i], facet_domains, bc.region)
            dbc.apply(bc.functions[i].vector())


class AssimilationProblem(ProblemBase):
    """Problem setup for dolfin-adjoint annotation.

    Observations and initial condition for velocity are given.

    Controls are both boundary inflow and velocity initial condition.
    """

    def __init__(self, params, geometry, initial_velocity, observations):
        ProblemBase.__init__(self, params)
        self.geometry = geometry
        self.initialize_geometry(self.geometry.mesh, self.geometry.facet_domains)

        # Store observations for use in setting initial values of boundary control functions
        self._observations = observations
        self._initial_velocity_value = initial_velocity

    @classmethod
    def default_params(cls):
        params = ProblemBase.default_params()
        params.update(
            initial_g="z",
            )
        return params

    def controls(self, spaces):
        d = len(self._initial_velocity_value)

        # Create initial condition control functions
        initial_velocity = [Function(spaces.U, name="ic_u%d" % i) for i in range(d)]

        # Set initial condition control function values
        for i in range(d):
            initial_velocity[i].assign(self._initial_velocity_value[i])
            #initial_velocity[i].interpolate(self._initial_velocity_value[i]) # TODO: Handle nonmatching spaces

        # Boundary control list is extended in advance() for each timestep
        g_timesteps = []

        return Controls(initial_velocity, g_timesteps)

    def initial_conditions(self, spaces, controls):
        return InitialConditions(controls.initial_velocity, 0.0)

    def advance(self, t0, t, timestep, spaces, state, boundary_conditions, body_force, controls):
        "Advance boundary condition functions to control values at time t."
        Ug = spaces.U # TODO: Configure control space

        bc = boundary_conditions.bcu.inflow

        # Create new BC control functions at time t
        d = len(state[0])
        g_at_t = [Function(Ug, name="g%d_%d" % (i,timestep)) for i in range(d)]

        # Make a record of BC controls in list
        controls.g_timesteps.append((float(t), g_at_t))
        assert len(controls.g_timesteps) == timestep

        # Assign initial values to BC control functions at time t
        if self.params.initial_g == "u0":
            # Get 'initial guess' from initial condition, i.e. set g(t) = u(t=0)
            for i in range(d):
                g_at_t[i].assign(self._initial_velocity_value[i])
                #g_at_t[i].interpolate(self._initial_velocity_value[i]) # TODO: Handle nonmatching spaces

        elif self.params.initial_g == "z":
            # Get observation at time t
            tz, z = self._observations[timestep]
            assert abs(float(t) - tz) < 1e-6, "Expecting matching times!"
            # Get 'initial guess' from observation, i.e. set g(t) = z(t)
            for i in range(d):
                g_at_t[i].assign(z[i])
                #g_at_t[i].interpolate(z[i]) # TODO: Handle nonmatching spaces

        # Assign 'initial guess' control function values to the BC functions used in scheme forms
        for i in range(d):
            bc.functions[i].assign(g_at_t[i], annotate=True) # Force annotation
            #bc.functions[i].interpolate(g_at_t[i], annotate=True) # Force annotation # TODO: Handle nonmatching spaces


class CostFunctional(Parameterized):
    def __init__(self, params, geometry, observations):
        Parameterized.__init__(self, params)

        self._geometry = geometry
        self._observations = observations

        self.J_timesteps = []

        # TODO: Configure with these id tuples
        # Control id tuples to make the below code more generic
        self.controlled_boundary_ids = (geometry.boundary_ids.left,)
        self.uncontrolled_boundary_ids = (geometry.boundary_ids.right,)
        self.wall_boundary_ids = (geometry.boundary_ids.wall,)

        # Setup integration measures
        self.dx = Measure("dx", domain=geometry.mesh)
        self.ds = Measure("ds", domain=geometry.mesh, subdomain_data=geometry.facet_domains)
        self.dsw = self.ds(self.wall_boundary_ids)
        self.dsc = self.ds(self.controlled_boundary_ids)
        self.dsu = self.ds(self.uncontrolled_boundary_ids)

        # Wrap regularization parameters in Constants
        for k in self.params.keys():
            alpha = self.params[k]
            if abs(alpha) < 1e-15:
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
            alpha_u0_dn_controlled=0.0,
            alpha_u0_dn_uncontrolled=0.0,
            alpha_u0_wall=0.0,
            # Boundary control regularization
            alpha_g=0.0,
            alpha_g_t=0.0,
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
        if timestep == 0:
            u0 = as_vector(controls.initial_velocity)
            J_terms += [
                self.alpha_u0              * u0**2       * dx * dts,
                self.alpha_u0_div          * div(u0)**2  * dx * dts,
                self.alpha_u0_grad         * grad(u0)**2 * dx * dts,
                self.alpha_u0_grad_controlled   * grad(u0)**2 * dsc * dts,
                self.alpha_u0_grad_uncontrolled * grad(u0)**2 * dsu * dts,
                self.alpha_u0_dn_controlled     * Dn(u0)**2   * dsc * dts,
                self.alpha_u0_dn_uncontrolled   * Dn(u0)**2   * dsu * dts,
                self.alpha_u0_wall              * u0**2       * dsw * dts,
                ]

        # ... Update cost functional with boundary control regularization
        if timestep >= 1:
            # Get last boundary control function (there's none at timestep 0)
            assert len(controls.g_timesteps) == timestep, "Expecting control functions at every timestep >= 1."
            tg, g = controls.g_timesteps[timestep-1]
            assert abs(float(t) - tg) < 1e-6, "Expecting matching time!"
            g = as_vector(g)

            # Get dt
            #tz_prev = self._observations[timestep-1][0]
            #dt = Constant(tz - tz_prev)
            dt = 1.0 # TODO: Add dt to params

            # Get g from previous timestep, use u0 for first timestep
            if timestep == 1:
                g_prev = as_vector(controls.initial_velocity)
            else:
                g_prev = as_vector(controls.g_timesteps[timestep-2])

            # Add regularization of boundary control function to cost functional at time t
            J_terms += [
                self.alpha_g             * g**2       * dsc * dtt,
                self.alpha_g_t           * (1.0/dt)*(g-g_prev)**2 * dsc * dtt, # TODO: This messed up the initial condition control! Why?
                self.alpha_g_grad        * grad(g)**2 * dsc * dtt, # TODO: Use trace grad on boundary
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
        return InitialConditions(self._controls.initial_velocity, 0.0)

    def advance(self, t0, t, timestep, spaces, state, boundary_conditions, body_force, controls):
        "Advance boundary condition functions to control values at time t."
        bc = boundary_conditions.bcu.inflow
        d = len(bc.functions)
        tg, g_at_t = self._controls.g_timesteps[timestep-1]
        assert abs(float(t) - tg) < 1e-6, "Expecting matching time!"
        for i in range(d):
            bc.functions[i].assign(g_at_t[i])


def main():
    set_log_level(100)

    # This is to remove rounding errors caused by tensor representation when computing norms
    parameters["form_compiler"]["representation"] = "quadrature"


    # Configure geometry
    geometry_params = ParamDict(
            refinement_level=2,
            length=2.0,
            radius=0.5,
            )

    # Configure time
    time_params = ParamDict(
            dt=1e-2,
            T=0.01,#4,#8,
            num_periods=None,
            )

    # Configure cost functional
    J_params = ParamDict(
        # Initial control regularization
        alpha_u0=0.0,
        alpha_u0_div=1.0,
        alpha_u0_grad=0.0,
        alpha_u0_grad_controlled=0.0,
        alpha_u0_grad_uncontrolled=0.0,
        alpha_u0_dn_controlled=0.0,
        alpha_u0_dn_uncontrolled=0.0,
        alpha_u0_wall=1.0,
        # Boundary control regularization
        alpha_g=0.0,
        alpha_g_t=0.0,
        alpha_g_grad=0.0,
        alpha_g_volume=0.0,
        alpha_g_grad_volume=0.0,
        )

    # Configure problems for the three stages
    analytic_problem_params = ParamDict(time_params, geometry=geometry_params)
    forward_problem_params = ParamDict(time_params)
    final_problem_params = ParamDict(time_params)


    # Configure scheme (reused in all stages)
    scheme_params = ParamDict(
            # BC params
            nietche=ParamDict(
                enable=True,
                formulation=1,
                stabilize=True,
                gamma=1000.0,
                ),

            # Variational formulation params
            scale_by_dt=True,
            enable_convection=True, # False = Stokes

            # Annotation
            annotate=True,

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
    casedir = "results_{}_{}_{}".format(CoupledScheme.__name__, date, params_string)


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

        # Setup scheme
        scheme_params.annotate = False
        scheme = CoupledScheme(scheme_params)

        # Setup analytic_solver
        analytic_solver_params = ParamDict(
            enable_annotation=False,
            #timer_frequency=1,
            #check_memory_frequency=1,
            )
        analytic_solver = NSSolver(analytic_problem, scheme, analytic_postprocessor, analytic_solver_params)

        # Step through simulation
        observations = []
        for data in analytic_solver.isolve():
            # Store projection of u into observation space
            Uz = data.spaces.U # TODO: Configure different observation space
            u, p = data.state
            z = [project(u_comp, Uz, name="z_at_ts%d" % data.timestep) for u_comp in u]
            observations.append((data.t, z))

            # Debugging evaluations:
            compute_boundary_diffs = False
            if compute_boundary_diffs:
                # Should have that z(t) = u(t) = g(t) on boundary bc.region.
                # If this fails, the weak boundary conditions have failed.
                # Increase and decrease the gamma parameter to see these
                # comparisons tighten and loosen.
                bc = data.boundary_conditions.bcu.inflow
                dsc = analytic_problem.ds(bc.region)

                u, p = data.state
                u = as_vector(u) # u(t)
                u2 = assemble(u**2*dsc)

                z = as_vector(z) # g(t)
                z2 = assemble(z**2*dsc)
                diffz = sqrt(assemble((u-z)**2*dsc) / u2)

                if data.timestep >= 1:
                    g = as_vector(bc.functions) # g(t)
                    g2 = assemble(g**2*dsc)
                    diffg = sqrt(assemble((u-g)**2*dsc) / u2)
                    diffgz = sqrt(assemble((z-g)**2*dsc) / u2)

                print "Comparing control boundary values:"
                print "ts={}; t={}".format(data.timestep, float(data.t))
                print "    u(t)^2={}".format(u2)
                print "    z(t)^2={}".format(u2)
                if data.timestep >= 1:
                    print "    g(t)^2={}".format(g2)
                    print "    |u(t)-g(t)| / |u(t)|  = {}".format(diffg)
                    print "    |g(t)-z(t)| / |u(t)|  = {}".format(diffgz)
                print "    |u(t)-z(t)| / |u(t)| = {}".format(diffz)
                print

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
        z0 = observations[0][1]
        d = len(z0)
        Uz = z0[0].function_space()
        initial_velocity = [Function(Uz) for i in range(d)]

        if 1: #enable_u0_control: # TODO: Make parameter
            if 0:
                # Set initial velocity to 0
                pass
                #for j in range(d):
                #    initial_velocity[j].vector().zero()
            elif 1:
                # Set initial velocity to scaling of exact initial
                for j in range(d):
                    initial_velocity[j].vector().axpy(0.1, z0[j].vector())
            else:
                # Keep exact initial condition
                pass


        # DEBUGGING: Plot suspect initial condition directly
        if 0:
            print z0[0].id()
            print initial_velocity[0].id()
            plot(as_vector(z0), title="observed initial velocity")
            plot(as_vector(initial_velocity), title="applied initial velocity")
            interactive()
            return

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

        # Setup scheme
        scheme_params.annotate = True
        scheme = CoupledScheme(scheme_params)

        # Setup and run solver
        forward_solver_params = ParamDict(
            enable_annotation=True,
            #timer_frequency=1,
            #check_memory_frequency=1,
            )
        forward_solver = NSSolver(forward_problem, scheme, forward_postprocessor, forward_solver_params)

        # Setup cost functional
        CF = CostFunctional(J_params, geometry, observations)
        dx = CF.dx
        ds = CF.ds
        dsc = CF.dsc

        # Step through forward simulation and build cost functional
        diffs = []
        for data in forward_solver.isolve():
            CF.update(data.timestep, data.t, data.spaces, data.state, data.controls)

            # Compare observations with 'observations' from assimilation problem
            compute_observation_diffs = False
            if compute_observation_diffs:
                # Avoid annotating these steps
                parameters["adjoint"]["stop_annotating"] = True
                dx = CF.dx
                u, p = data.state
                t, z = observations[data.timestep]
                z2 = assemble(as_vector(z)**2*dx)
                error = sqrt(assemble((as_vector(u) - as_vector(z))**2*dx) / z2)
                print data.timestep, "ERROR", error, z2
                parameters["adjoint"]["stop_annotating"] = False

        # Stop annotating after forward simulation
        parameters["adjoint"]["stop_annotating"] = True

        # Extract flat list of boundary control components
        all_g_components = [gc for (tg, g) in data.controls.g_timesteps for gc in g]
        all_u0_components = data.controls.initial_velocity

        # Accumulate terms from cost functional.
        J = CF.accumulate()
        #print str(J)
        J = da.Functional(J)


    # Always write tape to files for eventual debugging
    adj_html("forward.html", "forward")
    adj_html("adjoint.html", "adjoint")


    # Try replaying with dolfin-adjoint (works within a small tolerance)
    run_replay = True
    if run_replay:
        success = da.replay_dolfin(forget=False, tol=1e-13)
        print
        print "Replay success (within tiny tolerance):", success
        print


    # TODO: Make this a function
    # Plot gradients
    run_plot_gradients = False
    if run_plot_gradients:
        # Select controls
        m_values = list(all_g_components)
        m = [da.Control(mv) for mv in m_values]

        # Compute dJ/dm at m and plot animated gradients on boundary:
        dJdm = da.compute_gradient(J, m, forget=False)
        Vdj = dJdm[0].function_space()
        djp = Function(Vdj)
        bc = data.boundary_conditions.bcu.inflow
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
    # Compare various functions by computing some norms
    run_consistency_computations = True
    if run_consistency_computations:
        dsc = CF.dsc

        # Check that boundary control functions and observations have matching timesteps:
        print "Checking timestamps of boundary controls vs observations"
        assert abs(observations[0][0]) < 1e-6, "Expecting first observation at time 0."
        for timestep in range(1, len(observations)):
            tg = data.controls.g_timesteps[timestep-1][0]
            tz = observations[timestep][0]
            assert abs(tz - tg) < 1e-6, "Expecting matching time!"

        # Compute some norms for reference values
        print "Computing |z(t)|, |g(t)| on control boundary for each timestep:"
        tz0, z0 = observations[0]
        z02 = assemble(as_vector(z0)**2*dsc)
        print "timestep %d;  |z(t)|^2 = %g;" % (0, z02)
        for timestep in range(1, len(observations)):
            tg, g = data.controls.g_timesteps[timestep-1]
            tz, z = observations[timestep]
            assert abs(tz - tg) < 1e-6, "Expecting matching time!"
            g2 = sqrt(assemble(as_vector(g)**2*dsc) / z02)
            z2 = sqrt(assemble(as_vector(z)**2*dsc) / z02)
            print "timestep %d;  |g(t)| / |z(0)| = %g" % (timestep, g2)
            print "              |z(t)| / |z(0)| = %g" % (z2,)

        # This should be zero if this g(t) = z(t0) is the initial value used in the forward problem
        print "Computing |z(t0)-g(t)| on control boundary for each timestep:"
        tz0, z0 = observations[0]
        for timestep in range(1, len(observations)):
            tg, g = data.controls.g_timesteps[timestep-1]
            tz, z = observations[timestep]
            assert abs(tz - tg) < 1e-6, "Expecting matching time!"
            gz0diff = sqrt(assemble((as_vector(g) - as_vector(z0))**2*dsc) / z02)
            gzdiff = sqrt(assemble((as_vector(g) - as_vector(z))**2*dsc) / z02)
            zz0diff = sqrt(assemble((as_vector(z) - as_vector(z0))**2*dsc) / z02)
            print "timestep %d;  |g(t)-z(0)| / |z(0)| = %g" % (timestep, gz0diff)
            print "              |g(t)-z(t)| / |z(0)| = %g" % (gzdiff,)
            print "              |g(t)-z(t)| / |z(0)| = %g" % (zz0diff,)


    # TODO: Make this a function
    # Test gradient computation with convergence order tests
    run_taylor_test = False
    if run_taylor_test:
        set_u0_choice = 1 # TODO: Make parameter?
        set_g_choice = 1 # TODO: Make parameter?

        dsc = CF.dsc
        Ug = data.spaces.U # TODO: Configure control space
        d = data.spaces.d

        # Set value of boundary control functions
        if set_g_choice == 0:
            print "Setting g(t) = 0"
            for i, (tg, g) in enumerate(data.controls.g_timesteps):
                for j in range(d):
                    g[j].vector().zero()
        if set_g_choice == 1:
            print "Setting g(t) = z(t)"
            for i, (tg, g) in enumerate(data.controls.g_timesteps):
                for j in range(d):
                    g[j].assign(observations[i+1][1][j])
                    #g[j].interpolate(observations[i+1][1][j]) # TODO: Handle nonmatching spaces

        # Set value of initial velocity control functions
        if set_u0_choice == 0:
            print "Setting u(0) = 0"
            for j in range(d):
                data.controls.initial_velocity[j].vector().zero()
        elif set_u0_choice == 1:
            print "Setting u(0) = z(0)"
            for j in range(d):
                data.controls.initial_velocity[j].assign(observations[0][1][j])
                #data.controls.initial_velocity[j].interpolate(observations[0][1][j]) # TODO: Handle nonmatching spaces

        if set_u0_choice == 1 and set_g_choice == 1:
            print "Since u0,g==analytic solution, J(m) should be zero."
        else:
            print "Since u0,g!=analytic solution, J(m) should be nonzero."

        # Select controls for componentwise taylor test
        m_values = list(all_u0_components)
        #m_values = list(all_u0_components) + list(all_g_components)
        #m_values = list(all_g_components)
        #m_values = [all_g_components[0]] # DEBUGGING
        #m_values = [all_g_components[-1]] # DEBUGGING

        conv_rates = []
        for i, mv in enumerate(m_values):
            m = da.Control(mv)
            RJ = da.ReducedFunctional(J, m)
            Jm = RJ(mv)
            dJdm = da.compute_gradient(J, m, forget=False)
            conv_rate = da.taylor_test(RJ, m, Jm, dJdm)
            conv_rates.append(conv_rate)
            print
            if i < 2:
                print "m = u_%d(t_%d)" % ((i%2), (i//2))
            else:
                print "m = g_%d(t_%d)" % ((i%2), (i//2))
            print "  J(m)      =", Jm
            print "  dJdm(m)   =", assemble(dJdm**2*dx)
            print "  conv_rate =", conv_rate
        print
        print "min_conv_rate =", min(conv_rates)
        return


    # TODO: Make this a function
    # Optimize! This is where the main work goes!
    run_minimize = True
    if run_minimize:
        enable_u0_control = True # TODO: Make a parameter
        enable_g_control = False # TODO: Make a parameter

        # Select controls
        m_values = []
        if enable_u0_control:
            m_values += list(all_u0_components)
        if enable_g_control:
            m_values += list(all_g_components)
        m = [da.Control(mv) for mv in m_values]

        # Clear control functions
        for i in range(len(m_values)):
            m_values[i].vector().zero()


        # TODO: Expose optimization parameters
        # TODO: Scale tolerance to problem somehow?
        # TODO: Solver method choice?
        RJ = da.ReducedFunctional(J, m)
        m_opt = minimize(RJ, tol=1e-14, method="L-BFGS-B", options={"disp": True, "maxiter": 100})


        # DEBUGGING: Plot suspect initial condition directly
        if 0:
            plot(as_vector(m_opt))
            interactive()


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
        d = data.spaces.d

        # Setup controls on the right format
        if enable_u0_control:
            u0_opt = m_opt[:d]
            m_opt = m_opt[d:]
        else:
            #u0_opt = data.controls.initial_velocity
            u0_opt = observations[0][1]

        if enable_g_control:
            g_opt = []
            for i, (tg, g) in enumerate(data.controls.g_timesteps):
                gv = m_opt[d*i: d*(i+1)]
                g_opt.append((tg, gv))
        else:
            assert len(m_opt) == 0
            g_opt = []
            for i, (tg, g) in enumerate(data.controls.g_timesteps):
                gv = observations[i+1][1]
                g_opt.append((tg, gv))

        controls = Controls(u0_opt, g_opt)


        # DEBUGGING: Plot suspect initial condition directly
        if 0:
            plot(as_vector(u0_opt))
            interactive()
            # TODO: This causes a floating point exception?!?!

        # Setup problem
        final_problem = FinalReplayProblem(final_problem_params, geometry, controls)

        # Setup postprocessor
        fp = dict(save=True, plot=False)
        final_postprocessor = PostProcessor(dict(casedir=casedir + "_final"))
        final_postprocessor.add_fields([
            Pressure(fp),
            Velocity(fp),
            ])

        # Setup scheme
        scheme_params.annotate = False
        scheme = CoupledScheme(scheme_params)

        # Setup and run solver
        final_solver_params = ParamDict(
            enable_annotation=False,
            #timer_frequency=1,
            #check_memory_frequency=1,
            )
        final_solver = NSSolver(final_problem, scheme, final_postprocessor, final_solver_params)

        # Step through forward simulation
        diffs = []
        for data in final_solver.isolve():
            pass


if __name__ == "__main__":
    main()
