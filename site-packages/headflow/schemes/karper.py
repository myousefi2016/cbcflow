from __future__ import division

__author__ = "Martin Sandve Alnaes <martinal@simula.no>"
__date__ = "2013-08-13"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ..core.nsscheme import *
from ..core.timesteps import compute_regular_timesteps
from ..core.schemeutils import assign_ics_mixed, make_velocity_bcs, make_rhs_pressure_bcs
from ..core.spaces import NSSpacePoolMixed


class Karper(NSScheme):
    "TODO: Describe..."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to CR1-P0
            u_degree = 1,
            p_degree = 0,
            u_family = "CR",

            reuse_lu_data=False,
            verbose=False,
            )
        return params

    def solve(self, problem, update):
        # Spatial parameters
        mesh = problem.mesh
        n  = FacetNormal(mesh)
        dx = problem.dx
        ds = problem.ds

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[0])

        # Function spaces
        spaces = NSSpacePoolMixed(mesh, self.params.u_degree, self.params.p_degree, u_family=self.params.u_family)
        W = spaces.W

        # Test and trial functions
        u, p = TrialFunctions(W)
        v, q = TestFunctions(W)

        # Solution functions
        up0 = Function(W, name="up0") # Previous timestep
        up1 = Function(W, name="up1") # Current timestep
        u0, p0 = split(up0)
        u1, p1 = split(up1)

        # Get problem specific functions
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)
        ics = problem.initial_conditions(spaces, controls)
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)

        # Apply initial conditions and use it as initial guess
        assign_ics_mixed(up0, spaces, ics)
        up1.assign(up0)

        # Make scheme-specific representation of bcs
        bcu = make_velocity_bcs(problem, spaces, bcs)
        Lbc = make_rhs_pressure_bcs(problem, spaces, bcs, v)

        # Problem parameters
        nu = Constant(problem.params.mu/problem.params.rho, name="nu")
        k  = Constant(dt, name="dt")
        f  = as_vector(problem.body_force(spaces, t))

        # Variational forms
        def neg(f):
            return Min(f, 0.0)
        def pos(f):
            return Max(f, 0.0)

        #Up = (
        #      cell_avg(u('+')) * dot(facet_avg(u0('-')), nn)
        #    + cell_avg(u('-')) * dot(facet_avg(u0('+')), nn)
        #    )

        # Outward pointing normal from E- towards E+
        nn = n('-')
        # Flux from E- to E+
        Q = dot(facet_avg(u0('-')), nn)
        # Upwind term, pick avg velocity from cell where it flows from
        Up = cell_avg(u('-')) * pos(Q) + cell_avg(u('+')) * neg(Q)

        # Debugging, these are printed each timestep:
        Up1 = cell_avg(u1('-')) * pos(Q) + cell_avg(u1('+')) * neg(Q)
        Ups = [Up1]
        Ups = []

        # Build upwind terms from boundary conditions
        Qu = dot(facet_avg(u0), n) # Outward flux from E
        abcUp = sum(
            dot(cell_avg(u) * pos(Qu), cell_avg(v))*ds(region)
            for (ubc, region) in bcs[0])
        LbcUp = sum(
            dot(as_vector(map(facet_avg,ubc)) * neg(Qu), cell_avg(v))*ds(region)
            for (ubc, region) in bcs[0])
        abcUp2 = sum(
            dot(cell_avg(u) * pos(Qu) + facet_avg(u) * neg(Qu), cell_avg(v))*ds(region)
            for (pbc, region) in bcs[1])

        a = (
            dot(cell_avg(u), v)*dx()
            + k('+')*dot(Up, cell_avg(v('+')) - cell_avg(v('-')))*dS()
            + k*nu*inner(grad(u), grad(v))*dx()
            - k*p*div(v)*dx()
            - k*div(u)*q*dx() # TODO: Does the sign here matter? Closer to symmetric.
            + k*abcUp
            + k*abcUp2
            )
        L = (
            dot(cell_avg(u0), v)*dx()
            + k*dot(f,v)*dx()
            - k*LbcUp
            + k*Lbc
            )

        # Create solver
        linear_problem = LinearVariationalProblem(a, L, up1, bcu)
        solver = LinearVariationalSolver(linear_problem)

        # Turn on the noise
        if self.params.verbose:
            solver.parameters["lu_solver"]["report"] = True
            solver.parameters["lu_solver"]["verbose"] = True
            #solver.parameters["newton_solver"]["report"] = True

        # Speed up solvers with reuse
        if self.params.reuse_lu_data:
            solver.parameters["lu_solver"]["reuse_factorization"] = True
            solver.parameters["lu_solver"]["same_nonzero_pattern"] = True
            #solver.parameters["reset_jacobian"] = False

        # Define stricter stopping criteria for "newton" (actually fixed point) solver
        #solver.parameters["newton_solver"]["absolute_tolerance"] = self.params.picard_tolerance
        #solver.parameters["newton_solver"]["relative_tolerance"] = self.params.picard_tolerance
        #solver.parameters["newton_solver"]["maximum_iterations"] = self.params.picard_max_iterations
        #solver.parameters["newton_solver"]["error_on_nonconvergence"] = self.params.picard_error_on_nonconvergence

        # Call update() with initial conditions
        update(u0, p0, float(t), start_timestep, spaces)

        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)

            for i,Upi in enumerate(Ups):
                print "TOTAL UP %s: %g" % (i, assemble(Upi**2*dS()))

            # Solve for up1
            #solver.solve()
            solve(a == L, up1, bcu)

            # Rotate functions for next timestep
            up0.assign(up1)

            # Update postprocessing
            # TODO: Pass controls and observations here?
            update(u0, p0, float(t), timestep, spaces)

        # Make sure annotation gets that the timeloop is over
        finalize_time(t)

        # Return some quantities from the local namespace
        states = (u0, p0)
        namespace = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return namespace
