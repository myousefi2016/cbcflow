from __future__ import division

__author__ = "Oyvind Evju <oyvinev@simula.no>"
__date__ = "2013-05-13"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from ..core.nsscheme import *
from ..core.rhsgenerator import *
from ..core.timesteps import compute_regular_timesteps
from ..core.utils import Timer, is_periodic, cbcflow_log
#from ..core.adaptivetimestepping import AdaptiveTimestepping
#from ..core.constanttimestepping import ConstantTimestepping
from ..core.schemeutils import (assign_ics_segregated,
                                make_segregated_velocity_bcs,
                                make_pressure_bcs,
                                make_penalty_pressure_bcs)
from ..core.spaces import NSSpacePoolSegregated


class IPCS_Stable(NSScheme):
    "Incremental pressure-correction scheme, stable version."

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.update(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            theta = 0.5,
            )
        return params

    def solve(self, problem, update, restart=None):
        timer = Timer(self.params.enable_timer)

        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)
        dims = range(mesh.topology().dim())
        theta = self.params.theta

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[start_timestep])

        # Define function spaces
        spaces = NSSpacePoolSegregated(mesh, self.params.u_degree, self.params.p_degree)
        V = spaces.V
        U = spaces.U
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(U)
        q = TestFunction(Q)
        u = TrialFunction(U)
        p = TrialFunction(Q)

        # Functions
        u0 = as_vector([Function(U, name="u0_%d"%d) for d in dims]) # u^{n-1}
        u1 = as_vector([Function(U, name="u1_%d"%d) for d in dims]) # u^n
        u2 = as_vector([Function(U, name="u2_%d"%d) for d in dims]) # u^{n+1}
        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")
        p2 = Function(Q, name="p2")

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Apply initial conditions and use it as initial guess
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_segregated(u0, p0, spaces, ics)
        for d in dims: u1[d].assign(u0[d])
        for d in dims: u2[d].assign(u1[d])
        p1.assign(p0)
        p2.assign(p1)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u1, p1, t, controls)
        bcu = make_segregated_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Remove boundary stress term is problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))

        # Tentative velocity step. Crank-Nicholson time-stepping is used for diffusion and convection.
        a1 = (1/k) * inner(v, u) * dx()
        a2 = inner(grad(v), nu*grad(u)) * dx()

        # Convection linearized as in Simo/Armero (1994)
        a_conv = v * sum((1.5*u1[r] - 0.5*u0[r]) * u.dx(r) for r in dims) * dx()
        if theta < 1.0:
            # Set a_conv to match rhs theta-weighting for RHSGenerator
            a_conv = Constant(1-theta)*a_conv
            Kconv_axpy_factor = theta/(1-theta)
        else:
            Kconv_axpy_factor = 1.0
        Kconv = Matrix() # assembled from a_conv in the time loop

        # Create the static part of the coefficient matrix for the tentative
        # velocity step. The convection is added in the time loop. We use a
        # single matrix for all dimensions, which means that the BCs must match
        # (they must apply to the full vector space at the vertex, not a
        # subspace.
        A_u_tent = assemble(a1+theta*a2)

        # Create matrices for generating the RHS
        B = assemble(a1-(1-theta)*a2)
        M = assemble(v*u*dx())

        # Define how to create the RHS for the tentative velocity. The RHS is
        # (of course) different for each dimension.
        rhs_u_tent = [None]*len(dims)
        for d in dims:
            C = assemble(-v*p*n[d]*ds() + v.dx(d)*p*dx())
            rhs_u_tent[d] = RhsGenerator(U)
            rhs_u_tent[d] += B, u1[d]
            rhs_u_tent[d] += C, p1
            rhs_u_tent[d] += M, f[d]
            if theta < 1.0:
                rhs_u_tent[d] -= Kconv, u1[d]

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx())
        rhs_p_corr = RhsGenerator(Q)
        rhs_p_corr += A_p_corr, p1
        Ku = [None]*len(dims)
        for d in dims:
            Ku[d] = assemble(-(1/k)*q*u.dx(d)*dx()) # TODO: Store forms in list, this is copied below
            rhs_p_corr += Ku[d], u2[d]

        # Velocity correction. Like for the tentative velocity, a single LHS is used.
        A_u_corr = M
        rhs_u_corr = [None]*len(dims)
        Kp = [None]*len(dims)
        for d in dims:
            Kp[d] = assemble(-k*inner(v, grad(p)[d])*dx())
            rhs_u_corr[d] = RhsGenerator(U)
            rhs_u_corr[d] += M, u2[d]
            rhs_u_corr[d] += Kp[d], p2
            rhs_u_corr[d] -= Kp[d], p1

        # Apply BCs to LHS
        for bc in bcu:
            bc[0].apply(A_u_tent)
            bc[0].apply(A_u_corr)

        for bc in bcp:
            bc.apply(A_p_corr)

        # Create solvers
        if self.params.solver_p:
            solver_p_params = self.params.solver_p
        elif len(bcp) == 0 or is_periodic(bcp):
            solver_p_params = self.params.solver_p_neumann
        else:
            solver_p_params = self.params.solver_p_dirichlet

        solver_u_tent = LinearSolver(*self.params.solver_u_tent)
        solver_p_corr = LinearSolver(*solver_p_params)
        solver_u_corr = LinearSolver(*self.params.solver_u_corr)

        for A,S in [(A_u_tent, solver_u_tent), (A_p_corr, solver_p_corr), (A_u_corr, solver_u_corr)]:
            S.set_operator(A)
            if 'preconditioner' in S.parameters:
                S.parameters['preconditioner']['structure'] = 'same'

        timer.completed("problem initialization")

        # Call update() with initial conditions
        update(u1, p1, float(t), start_timestep, spaces)
        timer.completed("initial postprocessor update")

        # Time loop
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u1, p1, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Assemble the u-dependent convection matrix. It is important that
            # it is assembled into the same tensor, because the tensor is
            # also stored in rhs. (And it's faster).
            if Kconv.size(0) == 0:
                # First time, just assemble normally
                assemble(a_conv, tensor=Kconv, reset_sparsity=True)
            else:
                # Subtract the convection for previous time step before re-assembling Kconv
                A_u_tent.axpy(-Kconv_axpy_factor, Kconv, True)
                assemble(a_conv, tensor=Kconv, reset_sparsity=False)

            # Either zero BC rows in Kconv, or re-apply BCs to A_u_tent after
            # the axpy (it doesn't matter which)
            #for bc in bcu:
            #    bc[0].zero(Kconv)

            A_u_tent.axpy(Kconv_axpy_factor, Kconv, True)
            for bc in bcu:
                bc[0].apply(A_u_tent)
            timer.completed("u_tent assemble convection & construct lhs")

            # Compute tentative velocity step
            for d in dims:
                b = rhs_u_tent[d]()

                for bc in bcu: bc[d].apply(b)
                timer.completed("u_tent construct rhs")
                iter = solver_u_tent.solve(u2[d].vector(), b)

                timer.completed("u_tent solve (%s, %d dofs, %d iter)"%(', '.join(self.params.solver_u_tent), b.size(), iter))

            # Pressure correction
            b = rhs_p_corr()
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            for bc in bcp: bc.apply(b)
            timer.completed("p_corr construct rhs")

            iter = solver_p_corr.solve(p2.vector(), b)
            if len(bcp) == 0 or is_periodic(bcp): normalize(p2.vector())
            timer.completed("p_corr solve (%s, %d dofs, %d iter)"%(', '.join(solver_p_params), b.size(), iter))

            # Velocity correction
            for d in dims:
                b = rhs_u_corr[d]()
                for bc in bcu: bc[d].apply(b)
                timer.completed("u_corr construct rhs")

                iter = solver_u_corr.solve(u2[d].vector(), b)
                timer.completed("u_corr solve (%s, %d dofs, %d iter)"%(', '.join(self.params.solver_u_corr), b.size(), iter))

            # TODO: Reimplement this in some way
            '''
            # Check and adapt timestepping
            if self.params.adaptive_timestepping:
                avg_num_iterations = sum(tentative_vel_iterations) / len(dims)
                timestep_modified, dt = timestepper.adjusted_timestep(avg_num_iterations)

                if not timestep_modified:
                    cbcflow_log(INFO, "Timestep not modified. (dt=%.3e)" %dt)

                if timestep_modified:
                    old_dt = float(k)
                    k.assign(dt)

                    cbcflow_log(INFO, "Modified timestep: %.3e -> %.3e" % (old_dt, dt))

                    # Reassemble dt-dependent matrices
                    assemble(a1+a2, tensor=A_u_tent)
                    assemble(a1-a2, tensor=B)

                    # Reset Kconv to avoid subtraction from new A_u_tent
                    Kconv.assign(Matrix()) # TODO: Better to zero, to avoid resetting sparsity?

                    for d in dims:
                        assemble(-(1/k)*q*u.dx(d)*dx(), tensor=Ku[d]) # TODO: Store forms in list, this is copied from above
                        assemble(-k*inner(v, grad(p)[d])*dx(), tensor=Kp[d])

                    # Restart timestep if dt has decreased
                    if dt < old_dt:
                        continue
            '''
            # Rotate functions for next timestep
            for d in dims: u0[d].assign(u1[d])
            for d in dims: u1[d].assign(u2[d])
            p1.assign(p2)

            # Update postprocessing
            update(u1, p1, float(t), timestep, spaces)
            #timesteps.append(float(t))

        # Return some quantities from the local namespace
        states = (u1, p1)
        namespace = {
            "spaces": spaces,
            "observations": observations,
            "controls": controls,
            "states": states,
            "t": t,
            "timesteps": timesteps,
            }
        return namespace
