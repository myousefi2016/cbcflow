# Copyright (C) 2010-2014 Simula Research Laboratory
#
# This file is part of CBCFLOW.
#
# CBCFLOW is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CBCFLOW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CBCFLOW. If not, see <http://www.gnu.org/licenses/>.
from __future__ import division



from cbcflow.core.nsscheme import *
from cbcflow.utils.common import epsilon, sigma, is_periodic
from cbcflow.utils.schemes import (compute_regular_timesteps,
                                      assign_ics_split,
                                      make_velocity_bcs,
                                      make_pressure_bcs,
                                      make_penalty_pressure_bcs)
from cbcflow.utils.core import NSSpacePoolSplit


class PISO(NSScheme):
    "PISO, Issa 1985, implemented according to algorithm in Versteeg and Malalasekera"

    def __init__(self, params=None):
        NSScheme.__init__(self, params)

    @classmethod
    def default_params(cls):
        params = NSScheme.default_params()
        params.replace(
            # Default to P1-P1
            u_degree = 1,
            p_degree = 1,
            )
        params.update(
            # Solver and preconditioners
            solver_p_periodic=("gmres", "hypre_euclid"),
            solver_p_dirichlet=("gmres", "ml_amg"),
            solver_u_corr=("bicgstab", "hypre_euclid"),
            )
        return params

    def solve(self, problem, update, timer):
        parameters["linear_algebra_backend"] = "PETSc"
        parameters["form_compiler"]["optimize"]     = True
        parameters["form_compiler"]["cpp_optimize"] = True

        # Get problem parameters
        mesh = problem.mesh
        dx = problem.dx
        ds = problem.ds
        n  = FacetNormal(mesh)

        # Timestepping
        dt, timesteps, start_timestep = compute_regular_timesteps(problem)
        t = Time(t0=timesteps[0])

        # Define function spaces
        spaces = NSSpacePoolSplit(mesh, self.params.u_degree, self.params.p_degree)
        V = spaces.V
        Q = spaces.Q

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions # TODO: Describe what these are
        u0 = Function(V, name="u0")
        #u1 = Function(V, name="u1")
        u2 = Function(V, name="u2")
        un = Function(V, name="un")
        unc = Function(V, name="unc")

        p0 = Function(Q, name="p0")
        p1 = Function(Q, name="p1")

        # Debug print
        if 0:
            num_dofs = u0.vector().size()
            print 'num u dofs is ', num_dofs
            num_dofs = u0.vector().size() + p0.vector().size()
            print 'tot num dofs is ', num_dofs

        # Get functions for data assimilation
        observations = problem.observations(spaces, t)
        controls = problem.controls(spaces)

        # Get initial conditions
        ics = problem.initial_conditions(spaces, controls)
        assign_ics_split(u0, p0, spaces, ics)
        #u1.assign(u0)
        p1.assign(p0)

        # Make scheme-specific representation of bcs
        bcs = problem.boundary_conditions(spaces, u0, p0, t, controls)
        bcu = make_velocity_bcs(problem, spaces, bcs)
        bcp = make_pressure_bcs(problem, spaces, bcs)

        # Remove boundary stress term if problem is periodic
        beta = 0 if is_periodic(bcp) else 1

        # Problem coefficients
        nu = Constant(problem.params.mu/problem.params.rho)
        k  = Constant(dt)
        f  = as_vector(problem.body_force(spaces, t))



        if False: # Cebral implementation of SUPG parameter
            dl = ((12/sqrt(2))*(tetrahedron.volume))**.33333333
            supg =   (1/k)*(dl**2/(2*sqrt(inner(u0,u0))*dl+4*nu))
            supgII = (1/k)*(dl**2/(2*sqrt(inner(u2,u2))*dl+4*nu))
            print 'Using supg with dl**2/(2*sqrt(inner(u,u))*dl+4*nu)/k'
        else:
            supg = Constant(1)
            supgII = Constant(1)
            print 'Using supg with Constant(1)'

        # Predictor
        F_u_tent = ( (1/k)*inner(v, un - u0)*dx()
                     + inner(v,grad(p0))*dx()
                     + inner(v, grad(un)*un)*dx()
                     + nu*inner(grad(v), grad(un))*dx()
                     - inner(v, f)*dx()
                     + supg*k*inner(grad(v)*un, grad(un)*(un))*dx()
                     )

        # C1, Pressure correction I
        a_p_corr = inner(grad(q), grad(p))*dx()
        L_p_corr = inner(grad(q), grad(p0))*dx() - (1/k)*q*div(un)*dx()

        # C1, Velocity correction I
        a_u_corr = inner(v, u)*dx()
        L_u_corr = inner(v, un)*dx() - k*inner(v, grad(p1-p0))*dx()

        # C2, Pressure correction II
        a_p_corrII = inner(grad(q), grad(p))*dx()
        L_p_corrII = inner(grad(q), grad(p1))*dx() - (1/k)*q*div(u2)*dx()

        # C2, Velocity correction I
        F_u_corrm = ( (1/k)*inner(v, unc - u0)*dx()
                      + inner(v, grad(p1))*dx()
                      + inner(v, grad(unc)*unc)*dx()
                      + nu*inner(grad(v), grad(unc))*dx()
                      - inner(v, f)*dx()
                      + supgII*k*inner(grad(v)*unc, grad(unc)*(unc))*dx()
            )

        # Assemble matrices
        A_p_corr = assemble(a_p_corr)
        A_u_corr = assemble(a_u_corr)
        A_p_corrII = assemble(a_p_corrII)
        print 'assemble  ok'

        # Setup nonlinear problems
        J = derivative(F_u_tent, un, u)
        uproblem = NonlinearVariationalProblem(F_u_tent, un, bcs=bcu, J=J)
        usolver = NonlinearVariationalSolver(uproblem)
        prm = usolver.parameters
        prm["newton_solver"]["absolute_tolerance"] = 1E-5
        prm["newton_solver"]["relative_tolerance"] = 1E-5
        prm["newton_solver"]["maximum_iterations"] = 20
        prm["newton_solver"]["relaxation_parameter"] = .99
        prm["newton_solver"]["error_on_nonconvergence"] = True
        prm['linear_solver'] = 'gmres'
        prm['preconditioner'] = 'hypre_euclid'
        prm['krylov_solver']['absolute_tolerance'] = 1E-7
        prm['krylov_solver']['relative_tolerance'] = 1E-6
        prm['krylov_solver']['maximum_iterations'] = 20000
        prm['krylov_solver']['monitor_convergence'] = True
        prm['krylov_solver']['nonzero_initial_guess'] = False
        prm['krylov_solver']['gmres']['restart'] = 40

        Jc = derivative(F_u_corrm, unc, u)
        ucproblem = NonlinearVariationalProblem(F_u_corrm, unc, bcs=bcu, J=Jc)
        ucsolver = NonlinearVariationalSolver(ucproblem)
        prmc = ucsolver.parameters
        prmc["newton_solver"]["absolute_tolerance"] = 1E-5
        prmc["newton_solver"]["relative_tolerance"] = 1E-5
        prmc["newton_solver"]["maximum_iterations"] = 50
        prmc["newton_solver"]["relaxation_parameter"] = .99
        prmc["newton_solver"]["error_on_nonconvergence"] = True
        prmc['linear_solver'] = 'gmres'
        prmc['preconditioner'] = 'hypre_euclid'
        prmc['krylov_solver']['absolute_tolerance'] = 1E-7
        prmc['krylov_solver']['relative_tolerance'] = 1E-6
        prmc['krylov_solver']['maximum_iterations'] = 20000
        prmc['krylov_solver']['monitor_convergence'] = True
        prmc['krylov_solver']['nonzero_initial_guess'] = False
        prmc['krylov_solver']['gmres']['restart'] = 40
        print 'Second Newton solver ok'


        if len(bcp) == 0 or is_periodic(bcp):
            solver_p_params = self.params.solver_p_periodic
        else:
            solver_p_params = self.params.solver_p_dirichlet


        # Call update() with initial conditions
        update(u0, p0, float(t), start_timestep, spaces)


        # Loop over fixed timesteps
        for timestep in xrange(start_timestep+1,len(timesteps)):
            assign_time(t, timesteps[timestep])

            # Update various functions
            problem.update(spaces, u0, p0, t, timestep, bcs, observations, controls)
            timer.completed("problem update")

            # Compute tentative velocity step
            usolver.solve()

            # Pressure correction I
            b = assemble(L_p_corr)
            if len(bcp) == 0 or is_periodic(bcp):
                normalize(b)
            for bc in bcp:
                bc.apply(A_p_corr, b)
            iter = solve(A_p_corr, p1.vector(), b, *solver_p_params)
            if len(bcp) == 0 or is_periodic(bcp):
                normalize(p1.vector())

            # Velocity correction I
            b = assemble(L_u_corr)
            for bc in bcu:
                bc.apply(A_u_corr, b)
            iter = solve(A_u_corr, u2.vector(), b, *self.params.solver_u_corr)

            # Pressure correction II
            b = assemble(L_p_corrII)
            if len(bcp) == 0 or is_periodic(bcp):
                normalize(b)
            for bc in bcp:
                bc.apply(A_p_corrII, b)
            iter = solve(A_p_corrII, p1.vector(), b, *solver_p_params)
            if len(bcp) == 0 or is_periodic(bcp):
                normalize(p1.vector())

            # Velocity correction momentum eq
            unc.assign(u2)
            ucsolver.solve()

            # Rotate functions for next timestep
            u0.assign(unc)
            p0.assign(p1)

            # Update postprocessing
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
